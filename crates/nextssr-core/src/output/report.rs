//! Self-contained HTML report (`.report.html`).
//!
//! Unlike MISA's `.statistics` file (still available as a format), the report
//! gives the standard figures used in SSR surveys, i.e. relative abundance
//! (SSRs/Mb) and relative density (SSR bp/Mb) overall, per unit size and per
//! motif class, length and repeat-number distributions, per-sequence figures,
//! compound/imperfect summaries, primer design and in-silico PCR results,
//! complete parameters and provenance, and a ready-to-use methods paragraph.
//! It has no external dependencies (works offline, prints cleanly, follows the
//! light/dark preference of the viewer).

use std::collections::BTreeMap;
use std::fmt::Write;

use crate::epcr::EpcrParams;
use crate::motif::{canonical, misa_class, unit_name};
use crate::params::{Config, Mode};
use crate::stats::Stats;

/// In-silico PCR results to include.
#[derive(Debug, Clone)]
pub struct EpcrInfo<'a> {
    pub params: &'a EpcrParams,
    pub pairs: usize,
    pub by_status: BTreeMap<&'static str, u64>,
    pub loci_with_specific_pair: u64,
}

/// Provenance of the run.
#[derive(Debug, Clone, Default)]
pub struct ReportInfo<'a> {
    /// Input file name.
    pub source: &'a str,
    pub sha256: Option<&'a str>,
    pub bytes: Option<u64>,
    pub command: Option<&'a str>,
    /// e.g. "command line" or "web browser (WebAssembly)".
    pub runtime: &'a str,
    pub threads: Option<usize>,
    pub seconds: Option<f64>,
    /// Date/time of the analysis, already formatted.
    pub generated: Option<&'a str>,
    pub epcr: Option<EpcrInfo<'a>>,
}

fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => o.push_str("&amp;"),
            '<' => o.push_str("&lt;"),
            '>' => o.push_str("&gt;"),
            '"' => o.push_str("&quot;"),
            '\'' => o.push_str("&#39;"),
            _ => o.push(c),
        }
    }
    o
}

/// `1234567` -> `1,234,567`.
fn int(n: u64) -> String {
    let s = n.to_string();
    let mut o = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            o.push(',');
        }
        o.push(c);
    }
    o
}

fn dec(x: f64, d: usize) -> String {
    if !x.is_finite() {
        return "–".into();
    }
    let s = format!("{x:.d$}");
    // thousands separator on the integer part
    let (ip, fp) = s.split_once('.').map_or((s.as_str(), ""), |(a, b)| (a, b));
    let neg = ip.starts_with('-');
    let digits = ip.trim_start_matches('-');
    let mut o = String::new();
    if neg {
        o.push('-');
    }
    o.push_str(&int(digits.parse().unwrap_or(0)));
    if !fp.is_empty() {
        o.push('.');
        o.push_str(fp);
    }
    o
}

fn bp(n: u64) -> String {
    let x = n as f64;
    if n < 1_000 {
        format!("{n} bp")
    } else if n < 1_000_000 {
        format!("{} kb", dec(x / 1e3, 2))
    } else if n < 1_000_000_000 {
        format!("{} Mb", dec(x / 1e6, 2))
    } else {
        format!("{} Gb", dec(x / 1e9, 2))
    }
}

fn pct(a: u64, b: u64) -> String {
    if b == 0 {
        "–".into()
    } else {
        format!("{}%", dec(100.0 * a as f64 / b as f64, 1))
    }
}

/// Median of a value -> count histogram.
fn median(h: &BTreeMap<usize, u64>) -> Option<f64> {
    let n: u64 = h.values().sum();
    if n == 0 {
        return None;
    }
    let at = |k: u64| {
        let mut acc = 0;
        for (v, c) in h {
            acc += c;
            if acc > k {
                return *v as f64;
            }
        }
        0.0
    };
    Some(if n % 2 == 1 {
        at(n / 2)
    } else {
        (at(n / 2 - 1) + at(n / 2)) / 2.0
    })
}

fn mean(h: &BTreeMap<usize, u64>) -> Option<f64> {
    let n: u64 = h.values().sum();
    (n > 0).then(|| h.iter().map(|(v, c)| *v as f64 * *c as f64).sum::<f64>() / n as f64)
}

/// Bin a histogram into at most `max_bins` unit-width bins starting at its
/// minimum; the remainder is folded into a final `≥ x` bin.
fn fold(h: &BTreeMap<usize, u64>, max_bins: usize) -> Vec<(String, u64)> {
    let (Some((&lo, _)), Some((&hi, _))) = (h.first_key_value(), h.last_key_value()) else {
        return Vec::new();
    };
    let last = if hi - lo + 1 > max_bins {
        lo + max_bins - 1
    } else {
        hi
    };
    let mut out = Vec::new();
    for v in lo..last {
        out.push((v.to_string(), h.get(&v).copied().unwrap_or(0)));
    }
    let tail: u64 = h.range(last..).map(|(_, c)| c).sum();
    out.push((
        if last < hi {
            format!("≥{last}")
        } else {
            last.to_string()
        },
        tail,
    ));
    out
}

/// Column chart (single series) as inline SVG; every column has a tooltip.
fn columns(bins: &[(String, u64)], x_title: &str, unit: &str) -> String {
    let (w, h) = (640.0, 200.0);
    let (ml, mr, mt, mb) = (48.0, 8.0, 10.0, 40.0);
    let pw = w - ml - mr;
    let ph = h - mt - mb;
    let max = bins.iter().map(|b| b.1).max().unwrap_or(0).max(1);
    // round the axis maximum to 1, 2 or 5 × 10^k
    let raw = max as f64 / 4.0;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|m| m * mag)
        .find(|s| *s >= raw)
        .unwrap_or(raw)
        .max(1.0)
        .ceil();
    let top = step * (max as f64 / step).ceil();
    let n = bins.len().max(1) as f64;
    let band = pw / n;
    let bw = (band - 2.0).clamp(1.0, 24.0);
    let mut s = format!(
        r#"<svg class="chart" viewBox="0 0 {w} {h}" role="img" aria-label="{}"><g class="grid">"#,
        esc(x_title)
    );
    let mut v = 0.0;
    while v <= top + 1e-9 {
        let y = mt + ph - ph * v / top;
        let _ = write!(
            s,
            r#"<line x1="{ml}" x2="{}" y1="{y:.1}" y2="{y:.1}"/><text x="{}" y="{:.1}" class="yl">{}</text>"#,
            w - mr,
            ml - 6.0,
            y + 4.0,
            int(v as u64)
        );
        v += step;
    }
    s.push_str("</g><g class=\"cols\">");
    let label_every = ((n / 16.0).ceil() as usize).max(1);
    for (i, (label, c)) in bins.iter().enumerate() {
        let x = ml + band * i as f64 + (band - bw) / 2.0;
        let bh = ph * *c as f64 / top;
        let y0 = mt + ph;
        let _ = write!(
            s,
            r#"<g class="col"><rect class="hit" x="{:.1}" y="{mt}" width="{band:.1}" height="{ph}"/>"#,
            ml + band * i as f64
        );
        if bh > 0.0 {
            let r = 4f64.min(bw / 2.0).min(bh);
            let _ = write!(
                s,
                r#"<path d="M{x:.1},{y0:.1}V{:.1}Q{x:.1},{:.1} {:.1},{:.1}H{:.1}Q{:.1},{:.1} {:.1},{:.1}V{y0:.1}Z"/>"#,
                y0 - bh + r,
                y0 - bh,
                x + r,
                y0 - bh,
                x + bw - r,
                x + bw,
                y0 - bh,
                x + bw,
                y0 - bh + r
            );
        }
        let _ = write!(
            s,
            "<title>{}: {} {}</title></g>",
            esc(label),
            int(*c),
            esc(unit)
        );
        if i % label_every == 0 || i + 1 == bins.len() {
            let _ = write!(
                s,
                r#"<text x="{:.1}" y="{:.1}" class="xl">{}</text>"#,
                x + bw / 2.0,
                y0 + 15.0,
                esc(label)
            );
        }
    }
    let _ = write!(
        s,
        r#"</g><text x="{:.1}" y="{:.1}" class="xt">{}</text></svg>"#,
        ml + pw / 2.0,
        h - 4.0,
        esc(x_title)
    );
    s
}

/// Table cell with a proportional bar.
fn bar_cell(v: u64, max: u64) -> String {
    let w = if max == 0 {
        0.0
    } else {
        100.0 * v as f64 / max as f64
    };
    format!(r#"<td class="barcell"><span class="bar" style="width:{w:.1}%"></span></td>"#)
}

fn tile(o: &mut String, value: &str, label: &str, note: &str) {
    let _ = write!(
        o,
        r#"<div class="tile"><b>{value}</b><span>{}</span>{}</div>"#,
        esc(label),
        if note.is_empty() {
            String::new()
        } else {
            format!("<small>{}</small>", esc(note))
        }
    );
}

const CSS: &str = r#"
:root{color-scheme:light;--surface:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#77756f;--line:#e4e3df;--grid:#ecebe7;--series:#2a78d6;--series-soft:#2a78d61f;--accent:#0f766e;--good:#067647;--warn:#b54708;--bad:#b42318;--code:#f1f0ec}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])){color-scheme:dark;--surface:#1a1a19;--panel:#222220;--ink:#f0efec;--ink2:#c3c2b7;--muted:#9a988f;--line:#3a3935;--grid:#2e2d2a;--series:#3987e5;--series-soft:#3987e526;--accent:#2dd4bf;--good:#47cd89;--warn:#fdb022;--bad:#f97066;--code:#2b2a27}}
:root[data-theme="dark"]{color-scheme:dark;--surface:#1a1a19;--panel:#222220;--ink:#f0efec;--ink2:#c3c2b7;--muted:#9a988f;--line:#3a3935;--grid:#2e2d2a;--series:#3987e5;--series-soft:#3987e526;--accent:#2dd4bf;--good:#47cd89;--warn:#fdb022;--bad:#f97066;--code:#2b2a27}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1080px;margin:0 auto;padding:28px 16px 48px}
header.top{border-bottom:1px solid var(--line);padding-bottom:14px;margin-bottom:8px}
.kicker{color:var(--accent);font-weight:600;font-size:13px;letter-spacing:.04em;text-transform:uppercase;margin:0}
h1{font-size:26px;margin:4px 0 6px;overflow-wrap:anywhere}
h2{font-size:19px;margin:34px 0 4px}
h3{font-size:15px;margin:18px 0 6px;color:var(--ink2)}
.sub{color:var(--ink2);margin:0}
p.lead{color:var(--ink2);margin:2px 0 12px;max-width:80ch}
nav.toc{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:13.5px;margin:12px 0 0}
nav.toc a{color:var(--ink2)}
a{color:var(--accent)}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.9em}
code{background:var(--code);padding:1px 5px;border-radius:4px}
.tiles{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px;margin:14px 0}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.tile b{display:block;font-size:21px;font-variant-numeric:tabular-nums;line-height:1.25}
.tile span{color:var(--ink2);font-size:13px}
.tile small{display:block;color:var(--muted);font-size:12px}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13.5px;font-variant-numeric:tabular-nums}
th,td{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{color:var(--ink2);font-weight:600;font-size:12.5px}
td.n,th.n{text-align:right}
tbody tr:hover{background:var(--series-soft)}
td.barcell{width:28%;min-width:90px}
.bar{display:block;height:10px;background:var(--series);border-radius:0 4px 4px 0;min-width:1px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin:10px 0}
.methods{white-space:normal;max-width:90ch}
.copy{float:right;font:inherit;font-size:12.5px;border:1px solid var(--line);background:var(--surface);color:var(--ink2);border-radius:6px;padding:2px 9px;cursor:pointer}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}
svg.chart{width:100%;height:auto;display:block}
svg .grid line{stroke:var(--grid);stroke-width:1}
svg text{fill:var(--muted);font-size:11px;font-family:inherit}
svg .yl{text-anchor:end}
svg .xl,svg .xt{text-anchor:middle}
svg .xt{fill:var(--ink2)}
svg .cols path{fill:var(--series)}
svg .hit{fill:transparent}
svg .col:hover .hit{fill:var(--series-soft)}
.st{display:inline-block;padding:0 8px;border-radius:99px;font-size:12.5px;font-weight:600}
.st.good{color:var(--good);background:color-mix(in srgb,var(--good) 14%,transparent)}
.st.warn{color:var(--warn);background:color-mix(in srgb,var(--warn) 14%,transparent)}
.st.bad{color:var(--bad);background:color-mix(in srgb,var(--bad) 14%,transparent)}
dl.kv{display:grid;grid-template-columns:minmax(150px,max-content) 1fr;gap:4px 16px;margin:0;font-size:13.5px}
dl.kv dt{color:var(--ink2)}
dl.kv dd{margin:0;overflow-wrap:anywhere}
details summary{cursor:pointer;color:var(--ink2);font-size:13.5px;margin:6px 0}
.note{color:var(--muted);font-size:12.5px}
footer{margin-top:40px;color:var(--muted);font-size:12.5px;border-top:1px solid var(--line);padding-top:12px}
@media print{body{background:#fff}main{max-width:none;padding:0}.copy,nav.toc{display:none}.panel,.tile{break-inside:avoid}h2{break-after:avoid}details{display:block}details>*{display:block}}
"#;

fn units_list(cfg: &Config) -> String {
    let parts: Vec<String> = cfg
        .definition
        .min_repeats
        .iter()
        .map(|(u, r)| format!("{r} for {}", unit_name(*u)))
        .collect();
    match parts.len() {
        0 => String::new(),
        1 => parts[0].clone(),
        _ => format!(
            "{} and {}",
            parts[..parts.len() - 1].join(", "),
            parts[parts.len() - 1]
        ),
    }
}

fn methods(cfg: &Config, info: &ReportInfo) -> String {
    let mut m = format!(
        "Microsatellites (SSRs) were identified with nextSSR v{} (https://github.com/LaBiOmicS/nextSSR) ",
        crate::VERSION
    );
    match cfg.mode {
        Mode::Native => m.push_str(
            "in native mode, which reports every maximal perfect repeat of a primitive motif",
        ),
        Mode::Misa => m.push_str(
            "in MISA-compatible mode, which reproduces MISA v2.0 (Thiel et al. 2003; Beier et al. 2017)",
        ),
    }
    let _ = write!(
        m,
        ". The minimum number of repeats was {} motifs.",
        units_list(cfg)
    );
    match cfg.max_interruption {
        Some(d) if d <= 0 => m.push_str(
            " Each SSR was reported separately, except SSRs that were adjacent or overlapping, which were reported together as contiguous compound SSRs.",
        ),
        Some(d) => {
            let _ = write!(
                m,
                " SSRs separated by at most {d} bp were reported as compound SSRs."
            );
        }
        None => m.push_str(" SSRs were not grouped into compound SSRs."),
    }
    if let Some(ip) = &cfg.imperfect {
        let _ = write!(
            m,
            " Imperfect SSRs were assembled from perfect seeds of the same motif separated by at most {} bp, with a minimum purity of {}.",
            ip.max_gap, ip.min_purity
        );
    }
    if let Some(p) = &cfg.primers {
        let t = &p.thermo;
        let _ = write!(
            m,
            " PCR primers flanking each locus were designed with nextSSR (product size {}–{} bp; primer length {}–{} nt, optimum {}; Tm {}–{} °C, optimum {} °C, maximum difference {} °C; GC {}–{}%), with melting temperatures computed with the SantaLucia (1998) nearest-neighbour parameters for {} mM monovalent cations, {} mM Mg²⁺, {} mM dNTP and {} nM primer.",
            p.product_min, p.product_max, p.min_size, p.max_size, p.opt_size, p.min_tm, p.max_tm, p.opt_tm, p.max_tm_diff, p.min_gc, p.max_gc, t.mv_conc, t.dv_conc, t.dntp_conc, t.dna_conc
        );
        if let Some(e) = &info.epcr {
            let ep = e.params;
            let _ = write!(
                m,
                " Primer specificity was assessed by in-silico PCR against the input sequences on both strands, allowing up to {} mismatches per primer and {} within the {} 3′-terminal bases, and products up to {} bp; a pair was considered specific when its designed product was the only one.",
                ep.max_mismatches, ep.max_three_prime_mismatches, ep.three_prime_len, ep.max_product
            );
        }
    }
    m
}

/// Render the report.
pub fn html(st: &Stats, cfg: &Config, info: &ReportInfo) -> String {
    let mb = st.total_bp as f64 / 1e6;
    let per_mb = |x: u64| {
        if mb > 0.0 {
            dec(x as f64 / mb, 2)
        } else {
            "–".into()
        }
    };
    let ssr_bp: u64 = st.per_unit_bp.values().sum();
    let mut o = String::with_capacity(64 << 10);
    let _ = write!(
        o,
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>nextSSR report — {}</title><style>{CSS}</style></head><body><main>",
        esc(info.source)
    );
    let _ = write!(
        o,
        r##"<header class="top"><p class="kicker">nextSSR {} · microsatellite report</p><h1>{}</h1><p class="sub">{} mode · {}{}</p>
<nav class="toc"><a href="#overview">Overview</a><a href="#methods">Methods</a><a href="#units">Unit sizes</a><a href="#motifs">Motif classes</a><a href="#dist">Distributions</a><a href="#seqs">Sequences</a>{}<a href="#params">Parameters</a></nav></header>"##,
        crate::VERSION,
        esc(info
            .source
            .rsplit(['/', '\\'])
            .next()
            .unwrap_or(info.source)),
        match cfg.mode {
            Mode::Native => "Native",
            Mode::Misa => "MISA-compatible",
        },
        esc(info.runtime),
        info.generated
            .map(|g| format!(" · {}", esc(g)))
            .unwrap_or_default(),
        if cfg.primers.is_some() {
            r##"<a href="#primers">Primers</a>"##
        } else {
            ""
        }
    );

    // ---------- overview ----------
    o.push_str(r#"<section id="overview"><h2>Overview</h2><div class="tiles">"#);
    tile(&mut o, &int(st.sequences), "sequences", "");
    tile(
        &mut o,
        &bp(st.total_bp),
        "total size",
        &format!("{} bp", int(st.total_bp)),
    );
    tile(
        &mut o,
        &format!("{}%", dec(st.gc_percent(), 1)),
        "GC content",
        &{
            let n = st.composition[4];
            if n > 0 {
                format!("{} non-ACGT bases ({})", int(n), pct(n, st.total_bp))
            } else {
                String::new()
            }
        },
    );
    tile(
        &mut o,
        &int(st.ssrs),
        "SSRs",
        &format!("in {} loci", int(st.loci)),
    );
    tile(
        &mut o,
        &per_mb(st.ssrs),
        "SSRs per Mb",
        "relative abundance",
    );
    tile(&mut o, &per_mb(ssr_bp), "SSR bp per Mb", "relative density");
    tile(
        &mut o,
        &pct(st.loci_bp, st.total_bp),
        "of the sequence in SSR loci",
        &bp(st.loci_bp),
    );
    tile(
        &mut o,
        &int(st.seqs_with_ssr),
        "sequences with SSRs",
        &pct(st.seqs_with_ssr, st.sequences),
    );
    if cfg.max_interruption.is_some() {
        tile(
            &mut o,
            &int(st.compounds),
            "compound SSRs",
            &format!(
                "{} of SSRs are in compounds",
                pct(st.ssrs_in_compounds, st.ssrs)
            ),
        );
    }
    if cfg.imperfect.is_some() {
        tile(
            &mut o,
            &int(st.imperfect),
            "imperfect SSRs",
            &if st.imperfect > 0 {
                format!(
                    "mean purity {}",
                    dec(st.imperfect_purity_sum / st.imperfect as f64, 3)
                )
            } else {
                String::new()
            },
        );
    }
    if cfg.primers.is_some() {
        tile(
            &mut o,
            &int(st.primers_designed),
            "loci with primers",
            &pct(st.primers_designed, st.loci),
        );
        if let Some(e) = &info.epcr {
            tile(
                &mut o,
                &int(e.loci_with_specific_pair),
                "loci with a specific pair",
                &format!(
                    "{} of loci with primers",
                    pct(e.loci_with_specific_pair, st.primers_designed)
                ),
            );
        }
    }
    o.push_str("</div></section>");

    // ---------- methods ----------
    let _ = write!(
        o,
        r#"<section id="methods"><h2>Methods text</h2><p class="lead">A description of this analysis for a manuscript; edit as needed.</p><div class="panel"><button class="copy" type="button" onclick="navigator.clipboard&&navigator.clipboard.writeText(document.getElementById('mtext').innerText)">Copy</button><p id="mtext" class="methods">{}</p></div></section>"#,
        esc(&methods(cfg, info))
    );

    // ---------- unit sizes ----------
    let mut per_unit_reps: BTreeMap<usize, BTreeMap<usize, u64>> = BTreeMap::new();
    for (m, reps) in &st.per_motif {
        let e = per_unit_reps.entry(m.len()).or_default();
        for (r, c) in reps {
            *e.entry(*r).or_default() += c;
        }
    }
    o.push_str(r#"<section id="units"><h2>SSRs by unit size</h2><p class="lead">Relative abundance is the number of SSRs per Mb of sequence; relative density is the number of SSR base pairs per Mb.</p><div class="panel scroll"><table><thead><tr><th>Unit</th><th class="n">SSRs</th><th class="n">%</th><th></th><th class="n">per Mb</th><th class="n">bp per Mb</th><th class="n">Mean length</th><th class="n">Median length</th><th class="n">Length range</th><th class="n">Mean repeats</th></tr></thead><tbody>"#);
    let max_u = st.per_unit.values().copied().max().unwrap_or(0);
    for (u, c) in &st.per_unit {
        let lens = st.per_unit_len.get(u).cloned().unwrap_or_default();
        let reps = per_unit_reps.get(u).cloned().unwrap_or_default();
        let _ = write!(
            o,
            r#"<tr title="{}: {} SSRs"><td>{}</td><td class="n">{}</td><td class="n">{}</td>{}<td class="n">{}</td><td class="n">{}</td><td class="n">{}</td><td class="n">{}</td><td class="n">{}–{}</td><td class="n">{}</td></tr>"#,
            unit_name(*u),
            int(*c),
            unit_name(*u),
            int(*c),
            pct(*c, st.ssrs),
            bar_cell(*c, max_u),
            per_mb(*c),
            per_mb(st.per_unit_bp.get(u).copied().unwrap_or(0)),
            mean(&lens).map_or("–".into(), |x| dec(x, 1)),
            median(&lens).map_or("–".into(), |x| dec(x, 1)),
            lens.keys().next().copied().unwrap_or(0),
            lens.keys().last().copied().unwrap_or(0),
            mean(&reps).map_or("–".into(), |x| dec(x, 1)),
        );
    }
    if st.per_unit.is_empty() {
        o.push_str(r#"<tr><td colspan="10" class="note">No SSRs found.</td></tr>"#);
    }
    o.push_str("</tbody></table></div></section>");

    // ---------- motif classes ----------
    let classes = st.per_class();
    let mut class_bp: BTreeMap<String, u64> = BTreeMap::new();
    for (m, b) in &st.per_motif_bp {
        *class_bp.entry(misa_class(m)).or_default() += b;
    }
    let mut rows: Vec<(String, u64, u64, String)> = classes
        .iter()
        .map(|(k, reps)| {
            let n: u64 = reps.values().sum();
            let first = k.split('/').next().unwrap_or(k);
            (
                k.clone(),
                n,
                class_bp.get(k).copied().unwrap_or(0),
                canonical(first),
            )
        })
        .collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let shown = rows.len().min(25);
    let max_c = rows.first().map_or(0, |r| r.1);
    let _ = write!(
        o,
        r#"<section id="motifs"><h2>Motif classes</h2><p class="lead">Motifs are grouped with their rotations and reverse complements (e.g. AG, GA, CT and TC form the class AG/CT). {}</p><div class="panel scroll"><table><thead><tr><th>Class</th><th class="n">SSRs</th><th class="n">% of SSRs</th><th></th><th class="n">per Mb</th><th class="n">bp per Mb</th><th class="n">Mean length</th><th class="n">Max repeats</th></tr></thead><tbody>"#,
        if rows.len() > shown {
            format!("The {shown} most frequent of {} classes are shown; the TSV output has every locus.", rows.len())
        } else {
            String::new()
        }
    );
    for (k, n, b, _) in rows.iter().take(shown) {
        let maxr = classes[k].keys().last().copied().unwrap_or(0);
        let _ = write!(
            o,
            r#"<tr title="{}: {} SSRs"><td class="mono">{}</td><td class="n">{}</td><td class="n">{}</td>{}<td class="n">{}</td><td class="n">{}</td><td class="n">{}</td><td class="n">{}</td></tr>"#,
            esc(k),
            int(*n),
            esc(k),
            int(*n),
            pct(*n, st.ssrs),
            bar_cell(*n, max_c),
            per_mb(*n),
            per_mb(*b),
            if *n > 0 {
                dec(*b as f64 / *n as f64, 1)
            } else {
                "–".into()
            },
            maxr
        );
    }
    o.push_str("</tbody></table></div></section>");

    // ---------- distributions ----------
    let mut all_len: BTreeMap<usize, u64> = BTreeMap::new();
    for lens in st.per_unit_len.values() {
        for (l, c) in lens {
            *all_len.entry(*l).or_default() += c;
        }
    }
    o.push_str(r#"<section id="dist"><h2>Distributions</h2>"#);
    if !all_len.is_empty() {
        let _ = write!(
            o,
            r#"<h3>SSR length (all unit sizes)</h3><div class="panel">{}</div>"#,
            columns(&fold(&all_len, 48), "SSR length (bp)", "SSRs")
        );
        o.push_str(r#"<h3>Number of repeats by unit size</h3><div class="grid2">"#);
        for (u, reps) in &per_unit_reps {
            let _ = write!(
                o,
                r#"<div class="panel"><b>{}</b> <span class="note">({} SSRs)</span>{}</div>"#,
                unit_name(*u),
                int(reps.values().sum()),
                columns(&fold(reps, 24), "repeats", "SSRs")
            );
        }
        o.push_str("</div>");
        // full motif x repeats table, as in MISA but per class
        let min_r = classes
            .values()
            .filter_map(|r| r.keys().next().copied())
            .min()
            .unwrap_or(0);
        let max_r = classes
            .values()
            .filter_map(|r| r.keys().last().copied())
            .max()
            .unwrap_or(0);
        let cap = max_r.min(min_r + 29);
        let _ = write!(
            o,
            r#"<details><summary>Table: motif class × number of repeats ({} classes)</summary><div class="panel scroll"><table><thead><tr><th>Class</th>"#,
            rows.len()
        );
        for r in min_r..=cap {
            let _ = write!(
                o,
                r#"<th class="n">{}{r}</th>"#,
                if r == cap && cap < max_r { "≥" } else { "" }
            );
        }
        o.push_str(r#"<th class="n">total</th></tr></thead><tbody>"#);
        for (k, n, _, _) in &rows {
            let _ = write!(o, r#"<tr><td class="mono">{}</td>"#, esc(k));
            for r in min_r..=cap {
                let c: u64 = if r == cap {
                    classes[k].range(r..).map(|(_, c)| c).sum()
                } else {
                    classes[k].get(&r).copied().unwrap_or(0)
                };
                let _ = write!(
                    o,
                    r#"<td class="n">{}</td>"#,
                    if c > 0 { int(c) } else { String::new() }
                );
            }
            let _ = write!(o, r#"<td class="n">{}</td></tr>"#, int(*n));
        }
        o.push_str("</tbody></table></div></details>");
    } else {
        o.push_str(r#"<p class="note">No SSRs found.</p>"#);
    }
    if cfg.max_interruption.is_some() || cfg.imperfect.is_some() {
        o.push_str(r#"<h3>Compound and imperfect SSRs</h3><div class="panel"><dl class="kv">"#);
        if cfg.max_interruption.is_some() {
            let _ = write!(
                o,
                "<dt>Compound loci</dt><dd>{} ({} with members separated by an interruption, <code>c</code>; {} with adjacent or overlapping members, <code>c*</code>)</dd><dt>SSRs in compounds</dt><dd>{} ({} of SSRs)</dd>",
                int(st.compounds),
                int(st.compounds - st.compounds_contiguous),
                int(st.compounds_contiguous),
                int(st.ssrs_in_compounds),
                pct(st.ssrs_in_compounds, st.ssrs)
            );
        }
        if cfg.imperfect.is_some() {
            let _ = write!(
                o,
                "<dt>Imperfect SSRs</dt><dd>{} ({} of SSRs){}</dd>",
                int(st.imperfect),
                pct(st.imperfect, st.ssrs),
                if st.imperfect > 0 {
                    format!(
                        ", mean purity {}",
                        dec(st.imperfect_purity_sum / st.imperfect as f64, 3)
                    )
                } else {
                    String::new()
                }
            );
        }
        o.push_str("</dl></div>");
    }
    o.push_str("</section>");

    // ---------- sequences ----------
    let mut seqs: Vec<&crate::stats::SeqStat> = st.per_sequence.iter().collect();
    let limit = 40;
    let truncated = seqs.len() > limit;
    if truncated {
        seqs.sort_by(|a, b| b.length.cmp(&a.length).then_with(|| a.id.cmp(&b.id)));
        seqs.truncate(limit);
    }
    let max_d = seqs
        .iter()
        .map(|s| {
            if s.length > 0 {
                (s.ssrs as f64 / s.length as f64 * 1e9) as u64
            } else {
                0
            }
        })
        .max()
        .unwrap_or(0);
    let _ = write!(
        o,
        r#"<section id="seqs"><h2>Sequences</h2><p class="lead">{}</p><div class="panel scroll"><table><thead><tr><th>Sequence</th><th class="n">Length</th><th class="n">GC</th><th class="n">SSRs</th><th class="n">Loci</th><th class="n">SSRs per Mb</th><th></th><th class="n">bp per Mb</th></tr></thead><tbody>"#,
        if truncated {
            format!(
                "The {limit} longest of {} sequences ({} of them contain SSRs).",
                int(st.sequences),
                int(st.seqs_with_ssr)
            )
        } else {
            "Every sequence, in input order.".to_string()
        }
    );
    for s in &seqs {
        let d = |x: u64| {
            if s.length > 0 {
                dec(x as f64 * 1e6 / s.length as f64, 1)
            } else {
                "–".into()
            }
        };
        let dv = if s.length > 0 {
            (s.ssrs as f64 / s.length as f64 * 1e9) as u64
        } else {
            0
        };
        let _ = write!(
            o,
            r#"<tr title="{}: {} SSRs per Mb"><td class="mono">{}</td><td class="n">{}</td><td class="n">{}%</td><td class="n">{}</td><td class="n">{}</td><td class="n">{}</td>{}<td class="n">{}</td></tr>"#,
            esc(&s.id),
            d(s.ssrs),
            esc(&s.id),
            int(s.length),
            dec(s.gc, 1),
            int(s.ssrs),
            int(s.loci),
            d(s.ssrs),
            bar_cell(dv, max_d),
            d(s.ssr_bp)
        );
    }
    o.push_str("</tbody></table></div></section>");

    // ---------- primers ----------
    if let Some(p) = &cfg.primers {
        let tried = st.primers_designed + st.primer_failures.values().sum::<u64>();
        let _ = write!(
            o,
            r#"<section id="primers"><h2>Primers</h2><div class="tiles">"#
        );
        tile(
            &mut o,
            &format!("{} / {}", int(st.primers_designed), int(tried)),
            "loci with a primer pair",
            &pct(st.primers_designed, tried),
        );
        if p.num_return > 1 {
            tile(
                &mut o,
                &p.num_return.to_string(),
                "pairs requested per locus",
                "",
            );
        }
        if st.primers_designed > 0 {
            tile(
                &mut o,
                &format!(
                    "{} °C",
                    dec(st.primer_tm_sum / (2 * st.primers_designed) as f64, 1)
                ),
                "mean primer Tm",
                "best pairs",
            );
            let prod_mean = mean(&st.primer_products).unwrap_or(0.0);
            tile(
                &mut o,
                &format!("{} bp", dec(prod_mean, 0)),
                "mean product size",
                "best pairs",
            );
        }
        if let Some(e) = &info.epcr {
            tile(
                &mut o,
                &int(e.loci_with_specific_pair),
                "loci with a specific pair",
                &pct(e.loci_with_specific_pair, st.primers_designed),
            );
        }
        o.push_str("</div>");
        if !st.primer_products.is_empty() {
            // 10 bp bins over the product range
            let mut bins: Vec<(String, u64)> = Vec::new();
            let w = 10usize;
            let lo = p.product_min / w * w;
            let mut x = lo;
            while x <= p.product_max {
                let c: u64 = st.primer_products.range(x..x + w).map(|(_, c)| c).sum();
                bins.push((format!("{x}–{}", x + w - 1), c));
                x += w;
            }
            let _ = write!(
                o,
                r#"<h3>Product size of the best pair</h3><div class="panel">{}</div>"#,
                columns(&bins, "product size (bp)", "loci")
            );
        }
        if !st.primer_failures.is_empty() {
            o.push_str(r#"<h3>Loci without primers</h3><div class="panel scroll"><table><thead><tr><th>Reason</th><th class="n">Loci</th><th>Meaning</th></tr></thead><tbody>"#);
            for (r, c) in &st.primer_failures {
                let meaning = match r.as_str() {
                    "no_left_primer" => "no acceptable primer upstream of the SSR (often a sequence end, Ns or another SSR)",
                    "no_right_primer" => "no acceptable primer downstream of the SSR",
                    "no_compatible_pair" => "primers exist on both sides but no pair meets the product size, Tm difference or complementarity limits",
                    "target_too_long" => "the locus is longer than the maximum product allows",
                    _ => "",
                };
                let _ = write!(
                    o,
                    r#"<tr><td class="mono">{}</td><td class="n">{}</td><td>{}</td></tr>"#,
                    esc(r),
                    int(*c),
                    meaning
                );
            }
            o.push_str("</tbody></table></div>");
        }
        if let Some(e) = &info.epcr {
            o.push_str(r#"<h3>In-silico PCR against the input</h3><div class="panel scroll"><table><thead><tr><th>Status</th><th class="n">Pairs</th><th class="n">%</th><th>Meaning</th></tr></thead><tbody>"#);
            for (k, cls, meaning) in [
                ("specific", "good", "the designed product is the only product"),
                ("nonspecific", "warn", "more than one product (repeated region, paralogous locus or permissive primers)"),
                ("unexpected", "warn", "a single product that is not the designed one"),
                ("no_product", "bad", "no product within the search limits"),
            ] {
                let n = e.by_status.get(k).copied().unwrap_or(0);
                if n == 0 && k != "specific" && k != "nonspecific" {
                    continue;
                }
                let _ = write!(
                    o,
                    r#"<tr><td><span class="st {cls}">{}</span></td><td class="n">{}</td><td class="n">{}</td><td>{meaning}</td></tr>"#,
                    k.replace('_', " "),
                    int(n),
                    pct(n, e.pairs as u64)
                );
            }
            o.push_str(r#"</tbody></table><p class="note">The recommended pair of each locus (the best-ranked specific one) is marked in <code>.primers.tsv</code>; every product is listed in <code>.epcr.tsv</code>.</p></div>"#);
        }
        o.push_str("</section>");
    }

    // ---------- parameters & provenance ----------
    o.push_str(r#"<section id="params"><h2>Parameters and provenance</h2><div class="panel"><dl class="kv">"#);
    let mut kv = |k: &str, v: String| {
        let _ = write!(o, "<dt>{}</dt><dd>{}</dd>", esc(k), v);
    };
    kv("nextSSR version", esc(crate::VERSION));
    kv("Runtime", esc(info.runtime));
    if let Some(g) = info.generated {
        kv("Date", esc(g));
    }
    if let Some(c) = info.command {
        kv("Command", format!("<code>{}</code>", esc(c)));
    }
    kv(
        "Input",
        format!("<span class=\"mono\">{}</span>", esc(info.source)),
    );
    if let Some(b) = info.bytes {
        kv("Input bytes", int(b));
    }
    if let Some(h) = info.sha256 {
        kv(
            "Input SHA-256",
            format!("<span class=\"mono\">{}</span>", esc(h)),
        );
    }
    kv(
        "Detection mode",
        match cfg.mode {
            Mode::Native => "native (every maximal perfect repeat)".into(),
            Mode::Misa => "MISA-compatible (MISA v2.0 rules)".into(),
        },
    );
    kv(
        "Minimum repeats",
        format!(
            "<code>{}</code> (unit-repeats)",
            esc(&cfg.definition.to_misa_string())
        ),
    );
    kv(
        "Compound SSRs",
        cfg.max_interruption
            .map_or("not grouped".into(), |d| format!("members ≤ {d} bp apart")),
    );
    if let Some(ip) = &cfg.imperfect {
        kv(
            "Imperfect SSRs",
            format!(
                "max. gap {} bp, min. purity {}, seeds <code>{}</code>",
                ip.max_gap,
                ip.min_purity,
                cfg.definition
                    .units()
                    .map(|u| format!("{u}-{}", ip.seed_for(u)))
                    .collect::<Vec<_>>()
                    .join(" ")
            ),
        );
    }
    if let Some(p) = &cfg.primers {
        kv(
            "Primer design",
            format!(
                "product {}–{} bp; size {}/{}/{} nt; Tm {}/{}/{} °C (Δ ≤ {} °C); GC {}–{}%; poly-X ≤ {}; GC clamp {}; {} pair(s) per locus; {}",
                p.product_min, p.product_max, p.min_size, p.opt_size, p.max_size, p.min_tm, p.opt_tm, p.max_tm, p.max_tm_diff, p.min_gc, p.max_gc, p.max_poly_x, p.gc_clamp, p.num_return,
                if p.avoid_other_ssrs { "primers avoid other SSRs" } else { "primers may overlap other SSRs" }
            ),
        );
        let t = &p.thermo;
        kv(
            "Thermodynamics",
            format!(
                "SantaLucia (1998); {} mM monovalent, {} mM Mg²⁺, {} mM dNTP, {} nM oligo",
                t.mv_conc, t.dv_conc, t.dntp_conc, t.dna_conc
            ),
        );
    }
    if let Some(e) = &info.epcr {
        let p = e.params;
        kv(
            "In-silico PCR",
            format!(
                "≤ {} mismatches per primer, ≤ {} in the {} 3′ bases; products {}–{} bp; {} pairs",
                p.max_mismatches,
                p.max_three_prime_mismatches,
                p.three_prime_len,
                p.min_product,
                p.max_product,
                int(e.pairs as u64)
            ),
        );
    }
    if let Some(t) = info.threads {
        kv("Threads", t.to_string());
    }
    if let Some(s) = info.seconds {
        kv("Run time", format!("{} s", dec(s, 2)));
    }
    o.push_str("</dl></div></section>");

    o.push_str(r#"<footer>Generated by <a href="https://github.com/LaBiOmicS/nextSSR">nextSSR</a>. MISA-format statistics are available with <code>-f statistics</code>. If you use nextSSR, please cite the repository; when using the MISA-compatible mode, please also cite Thiel et al. (2003) <i>Theor Appl Genet</i> 106:411–422 and Beier et al. (2017) <i>Bioinformatics</i> 33:2583–2585.</footer></main></body></html>"#);
    o.push('\n');
    o
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_all;
    use crate::pipeline::analyze_record;
    use crate::primer::PrimerParams;

    #[test]
    fn helpers() {
        assert_eq!(int(0), "0");
        assert_eq!(int(1234567), "1,234,567");
        assert_eq!(dec(12345.678, 1), "12,345.7");
        assert_eq!(bp(12_071_326), "12.07 Mb");
        let h: BTreeMap<usize, u64> = [(10, 3), (11, 1), (20, 1)].into_iter().collect();
        assert_eq!(median(&h), Some(10.0));
        assert_eq!(mean(&h), Some(12.2));
        let f = fold(&h, 5);
        assert_eq!(f.len(), 5);
        assert_eq!(f[4], ("≥14".to_string(), 1));
        assert_eq!(esc("<a&\"b\">"), "&lt;a&amp;&quot;b&quot;&gt;");
    }

    #[test]
    fn report_contents() {
        let fa = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/NC_000932.1_Arabidopsis_cp.fasta"
        ))
        .unwrap();
        let cfg = Config {
            primers: Some(PrimerParams::default()),
            ..Config::default()
        };
        let mut st = Stats::new();
        for rec in parse_all(&fa) {
            let r = analyze_record(&rec, &cfg);
            st.add_record(&r.id, &rec.seq, &r.loci);
        }
        let html = html(
            &st,
            &cfg,
            &ReportInfo {
                source: "cp<1>.fa",
                runtime: "test",
                ..ReportInfo::default()
            },
        );
        assert!(html.starts_with("<!doctype html>"));
        assert!(html.contains("cp&lt;1&gt;.fa"));
        assert!(!html.contains("cp<1>"));
        // GC content of the Arabidopsis chloroplast genome is 36.3%
        assert!(html.contains("36.3%"), "GC");
        assert!(html.contains("A/T"));
        assert!(html.contains("Microsatellites (SSRs) were identified with nextSSR"));
        assert!(html.contains(r#"id="primers""#));
        assert!(!html.contains("NaN"));
        assert_eq!(html.matches("<svg").count(), html.matches("</svg>").count());
        // no external resources
        assert!(!html.contains("src=\"http") && !html.contains("<link"));
    }
}
