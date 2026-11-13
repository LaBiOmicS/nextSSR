//! nextSSR command-line interface.

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use rayon::prelude::*;
use sha2::{Digest, Sha256};

use nextssr_core::epcr::{self, Epcr, EpcrParams, EpcrReport, PairSpec};
use nextssr_core::output::{bed, gff3, misa, p3in, primers, report, tsv, Format};
use nextssr_core::primer::PrimerParams;
use nextssr_core::stats::misa_statistics;
use nextssr_core::thermo::{tm, ThermoParams};
use nextssr_core::{
    analyze_record, Config, Definition, FastaParser, ImperfectParams, MisaIni, Mode, Record,
    RecordResult, Stats,
};

const DEFAULT_INI: &str =
    "definition(unit_size,min_repeats):                   1-10 2-6 3-5 4-5 5-5 6-5\n\
interruptions(max_difference_between_2_SSRs):        100\n\
GFF:                                                 false\n";

#[derive(Parser)]
#[command(
    name = "nextssr",
    version,
    about = "Fast, MISA-compatible microsatellite (SSR) finder with imperfect/compound SSRs and primer design",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Find SSRs (nextSSR semantics by default) and write TSV/GFF3/BED/MISA/statistics/JSON outputs
    Run(Box<RunArgs>),
    /// Drop-in replacement for `misa.pl <FASTA> [misa.ini]`: writes <FASTA>.misa and <FASTA>.statistics
    Misa(MisaArgs),
    /// Print a default misa.ini
    InitIni,
    /// Compute primer melting temperatures (SantaLucia 1998, Primer3 default conditions)
    Tm(TmArgs),
    /// In-silico PCR: find every product of a set of primer pairs in FASTA sequences
    Epcr(Box<EpcrCmdArgs>),
}

#[derive(Clone, Copy, ValueEnum)]
enum ModeArg {
    /// All maximal perfect repeats, correct compound notation
    Native,
    /// Exact misa.pl v2.0 semantics
    Misa,
}

#[derive(Args)]
struct RunArgs {
    /// Input FASTA (plain or gzip/bgzip); `-` reads standard input
    input: String,

    /// Output prefix (default: the input path, as misa.pl does)
    #[arg(short, long)]
    output: Option<String>,

    /// Comma-separated output formats: tsv,gff3,bed,misa,statistics,json,p3in,primers,report
    /// (`report` = HTML report; `statistics` = MISA-format statistics;
    /// `primers` is added automatically with --primers)
    #[arg(short, long, default_value = "tsv,gff3,misa,report,json")]
    formats: String,

    /// SSR definition "unit-min_repeats" pairs
    #[arg(short, long, default_value = "1-10 2-6 3-5 4-5 5-5 6-5")]
    definition: String,

    /// Read definition and interruptions from a misa.ini file (overrides -d and -i)
    #[arg(long)]
    ini: Option<PathBuf>,

    /// Maximum number of bases between two SSRs reported as one compound SSR
    /// [default: 0 in native mode (only touching or overlapping SSRs),
    /// 100 in MISA mode]
    #[arg(short = 'i', long)]
    max_interruption: Option<i64>,

    /// Do not group SSRs into compound SSRs
    #[arg(long)]
    no_compound: bool,

    /// Detection semantics
    #[arg(long, value_enum, default_value = "native")]
    mode: ModeArg,

    /// Also detect imperfect (interrupted) SSRs
    #[arg(long)]
    imperfect: bool,

    /// Imperfect: minimum copies of perfect seeds, "unit-copies" pairs
    /// (default: 2 copies and >= 8 bp, i.e. 1-8 2-4 3-3 4-2 5-2 6-2)
    #[arg(long)]
    seed_repeats: Option<String>,

    /// Imperfect: maximum interruption (bp) between seeds
    #[arg(long, default_value_t = 5)]
    max_gap: usize,

    /// Imperfect: minimum purity (fraction of the span covered by seeds)
    #[arg(long, default_value_t = 0.85)]
    min_purity: f64,

    /// Report this many flanking bases on each side (TSV)
    #[arg(long, default_value_t = 0)]
    flank: usize,

    #[command(flatten)]
    primer: PrimerArgs,

    /// Check the specificity of every designed pair by in-silico PCR against the
    /// input sequences (requires --primers and a file input); writes <prefix>.epcr.tsv
    #[arg(long, requires = "primers")]
    epcr: bool,

    #[command(flatten)]
    epcr_opts: EpcrArgs,

    /// Worker threads (default: all logical CPUs)
    #[arg(short, long)]
    threads: Option<usize>,

    /// Do not compute the SHA-256 of the input
    #[arg(long)]
    no_hash: bool,

    /// Suppress the summary on stderr
    #[arg(short, long)]
    quiet: bool,
}

#[derive(Args)]
struct PrimerArgs {
    /// Design PCR primers for every locus
    #[arg(long)]
    primers: bool,
    /// Product size range "min-max"
    #[arg(long, default_value = "100-300")]
    product_size: String,
    /// Primer size "min,opt,max"
    #[arg(long, default_value = "18,20,27")]
    primer_size: String,
    /// Primer Tm (°C) "min,opt,max"
    #[arg(long, default_value = "57,60,63")]
    primer_tm: String,
    /// Primer GC% "min,max"
    #[arg(long, default_value = "30,70")]
    primer_gc: String,
    /// Maximum Tm difference between primers (°C)
    #[arg(long, default_value_t = 3.0)]
    max_tm_diff: f64,
    /// G/C bases required among the last five 3' bases
    #[arg(long, default_value_t = 0)]
    gc_clamp: usize,
    /// Longest mononucleotide run allowed in a primer
    #[arg(long, default_value_t = 5)]
    max_poly_x: usize,
    /// Allow primers to overlap other SSR loci
    #[arg(long)]
    allow_primers_on_ssrs: bool,
    /// Alternative primer pairs reported per locus (best first)
    #[arg(long, default_value_t = 1)]
    num_pairs: usize,
    #[command(flatten)]
    thermo: ThermoArgs,
}

#[derive(Args, Clone)]
struct ThermoArgs {
    /// Monovalent cation concentration (mM)
    #[arg(long, default_value_t = 50.0)]
    mv: f64,
    /// Divalent cation (Mg2+) concentration (mM)
    #[arg(long, default_value_t = 1.5)]
    dv: f64,
    /// dNTP concentration (mM)
    #[arg(long, default_value_t = 0.6)]
    dntp: f64,
    /// Oligonucleotide concentration (nM)
    #[arg(long, default_value_t = 50.0)]
    dna: f64,
}

impl ThermoArgs {
    fn params(&self) -> ThermoParams {
        ThermoParams {
            mv_conc: self.mv,
            dv_conc: self.dv,
            dntp_conc: self.dntp,
            dna_conc: self.dna,
        }
    }
}

#[derive(Args, Clone)]
struct EpcrArgs {
    /// e-PCR: maximum mismatches per primer binding site
    #[arg(long, default_value_t = 2)]
    mismatches: usize,
    /// e-PCR: length of the primer 3'-end region
    #[arg(long, default_value_t = 5)]
    three_prime_len: usize,
    /// e-PCR: maximum mismatches inside the 3'-end region
    #[arg(long, default_value_t = 0)]
    three_prime_mismatches: usize,
    /// e-PCR: product size range "min-max" (min 0 = longer primer)
    #[arg(long, default_value = "0-2000")]
    epcr_product: String,
    /// e-PCR: amplicons listed per pair (all are counted)
    #[arg(long, default_value_t = 100)]
    max_amplicons: usize,
}

impl EpcrArgs {
    fn params(&self, with_sequence: bool) -> Result<EpcrParams> {
        let ps: Vec<usize> = parse_list(&self.epcr_product, 2, "--epcr-product")?;
        let p = EpcrParams {
            max_mismatches: self.mismatches,
            three_prime_len: self.three_prime_len,
            max_three_prime_mismatches: self.three_prime_mismatches,
            min_product: ps[0],
            max_product: ps[1],
            max_amplicons: self.max_amplicons,
            with_sequence,
            ..EpcrParams::default()
        };
        p.validate()?;
        Ok(p)
    }
}

#[derive(Args)]
struct EpcrCmdArgs {
    /// Target FASTA files (plain or gzip/bgzip); `-` reads standard input
    #[arg(required = true)]
    targets: Vec<String>,

    /// Primer list: a table with forward/reverse (or left_primer/right_primer)
    /// columns such as nextSSR .primers.tsv, or lines "ID FORWARD REVERSE"
    #[arg(short, long, conflicts_with_all = ["forward", "reverse"])]
    primers: Option<PathBuf>,

    /// Forward primer of a single pair (with --reverse)
    #[arg(short = 'F', long, requires = "reverse")]
    forward: Option<String>,

    /// Reverse primer of a single pair (with --forward)
    #[arg(short = 'R', long, requires = "forward")]
    reverse: Option<String>,

    /// Output prefix: writes <prefix>.epcr.tsv (amplicons) and
    /// <prefix>.epcr_pairs.tsv (one line per pair). Default: amplicons to stdout
    #[arg(short, long)]
    output: Option<String>,

    /// Include the amplicon sequences
    #[arg(long)]
    sequence: bool,

    #[command(flatten)]
    epcr: EpcrArgs,

    /// Worker threads (default: all logical CPUs)
    #[arg(short, long)]
    threads: Option<usize>,

    /// Suppress the summary on stderr
    #[arg(short, long)]
    quiet: bool,
}

#[derive(Args)]
struct MisaArgs {
    /// Input FASTA file
    fasta: String,
    /// misa.ini (default: ./misa.ini, else MISA defaults)
    ini: Option<PathBuf>,
    /// Worker threads (default: all logical CPUs)
    #[arg(short, long)]
    threads: Option<usize>,
}

#[derive(Args)]
struct TmArgs {
    /// Oligonucleotide sequences
    #[arg(required = true)]
    sequences: Vec<String>,
    #[command(flatten)]
    thermo: ThermoArgs,
}

fn parse_list<T: std::str::FromStr>(s: &str, n: usize, what: &str) -> Result<Vec<T>> {
    let v: Vec<T> = s
        .split([',', '-', ' '])
        .filter(|x| !x.is_empty())
        .map(|x| {
            x.parse::<T>()
                .map_err(|_| anyhow::anyhow!("invalid {what}: {s:?}"))
        })
        .collect::<Result<_>>()?;
    if v.len() != n {
        bail!("{what} expects {n} values, got {s:?}");
    }
    Ok(v)
}

fn primer_params(a: &PrimerArgs) -> Result<PrimerParams> {
    let ps: Vec<usize> = parse_list(&a.product_size, 2, "--product-size")?;
    let sz: Vec<usize> = parse_list(&a.primer_size, 3, "--primer-size")?;
    let t: Vec<f64> = parse_list(&a.primer_tm, 3, "--primer-tm")?;
    let gc: Vec<f64> = parse_list(&a.primer_gc, 2, "--primer-gc")?;
    let p = PrimerParams {
        min_size: sz[0],
        opt_size: sz[1],
        max_size: sz[2],
        min_tm: t[0],
        opt_tm: t[1],
        max_tm: t[2],
        min_gc: gc[0],
        max_gc: gc[1],
        max_poly_x: a.max_poly_x,
        gc_clamp: a.gc_clamp,
        max_tm_diff: a.max_tm_diff,
        product_min: ps[0],
        product_max: ps[1],
        avoid_other_ssrs: !a.allow_primers_on_ssrs,
        num_return: a.num_pairs,
        thermo: a.thermo.params(),
        ..PrimerParams::default()
    };
    p.validate()?;
    Ok(p)
}

fn build_config(a: &RunArgs) -> Result<Config> {
    let (definition, max_int) = match &a.ini {
        Some(path) => {
            let text = std::fs::read_to_string(path)
                .with_context(|| format!("reading {}", path.display()))?;
            let ini = MisaIni::parse(&text)?;
            (ini.definition, ini.max_interruption)
        }
        None => (
            Definition::parse(&a.definition)?,
            a.max_interruption.unwrap_or(match a.mode {
                ModeArg::Native => nextssr_core::params::NATIVE_MAX_INTERRUPTION,
                ModeArg::Misa => nextssr_core::params::MISA_MAX_INTERRUPTION,
            }),
        ),
    };
    let imperfect = if a.imperfect {
        let mut ip = ImperfectParams {
            max_gap: a.max_gap,
            min_purity: a.min_purity,
            ..ImperfectParams::default()
        };
        if let Some(sr) = &a.seed_repeats {
            for (u, r) in Definition::parse(sr)?.min_repeats {
                ip.seed_repeats.insert(u, r);
            }
        }
        Some(ip)
    } else {
        None
    };
    let cfg = Config {
        definition,
        mode: match a.mode {
            ModeArg::Native => Mode::Native,
            ModeArg::Misa => Mode::Misa,
        },
        max_interruption: if a.no_compound { None } else { Some(max_int) },
        imperfect,
        primers: if a.primer.primers {
            Some(primer_params(&a.primer)?)
        } else {
            None
        },
        flank: a.flank,
    };
    cfg.validate()?;
    Ok(cfg)
}

/// Reader that hashes everything read through it.
struct HashingReader<R> {
    inner: R,
    hasher: Option<Sha256>,
    bytes: u64,
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if let Some(h) = &mut self.hasher {
            h.update(&buf[..n]);
        }
        self.bytes += n as u64;
        Ok(n)
    }
}

struct InputInfo {
    bytes: u64,
    sha256: Option<String>,
}

type Input = Box<dyn Read + Send>;

/// Open the input (a path or `-` for standard input). Called before any output
/// file is created, so that a missing input leaves nothing behind.
fn open_input(input: &str) -> Result<Input> {
    Ok(if input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(input).with_context(|| format!("cannot open {input}"))?)
    })
}

/// Stream FASTA records in batches to `on_batch`, reading on a separate thread.
fn stream_records(
    raw: Input,
    hash: bool,
    mut on_batch: impl FnMut(Vec<Record>) -> Result<()>,
) -> Result<InputInfo> {
    let mut hr = HashingReader {
        inner: raw,
        hasher: hash.then(Sha256::new),
        bytes: 0,
    };
    let mut magic = [0u8; 2];
    let mut got = 0;
    while got < 2 {
        let n = hr.read(&mut magic[got..])?;
        if n == 0 {
            break;
        }
        got += n;
    }
    let prefix = magic[..got].to_vec();
    let is_gz = got == 2 && magic == [0x1f, 0x8b];

    let (tx, rx) = crossbeam_channel::bounded::<Vec<Record>>(2);
    let reader = std::thread::spawn(move || -> Result<HashingReader<Box<dyn Read + Send>>> {
        const BATCH_BP: usize = 16 << 20;
        const BATCH_RECORDS: usize = 8192;
        let mut parser = FastaParser::new();
        let mut batch: Vec<Record> = Vec::new();
        let mut batch_bp = 0usize;
        let flush = |records: Vec<Record>, batch: &mut Vec<Record>, batch_bp: &mut usize| -> bool {
            for r in records {
                *batch_bp += r.seq.len();
                batch.push(r);
            }
            if *batch_bp >= BATCH_BP || batch.len() >= BATCH_RECORDS {
                *batch_bp = 0;
                return tx.send(std::mem::take(batch)).is_ok();
            }
            true
        };
        let mut buf = vec![0u8; 1 << 20];
        let mut out = Vec::new();
        let mut aborted = false;
        {
            let mut src: Box<dyn Read + '_> = if is_gz {
                Box::new(flate2::read::MultiGzDecoder::new(
                    io::Cursor::new(prefix.clone()).chain(&mut hr),
                ))
            } else {
                parser.push(&prefix, &mut out);
                Box::new(&mut hr)
            };
            loop {
                let n = src.read(&mut buf).context("reading input")?;
                if n == 0 {
                    break;
                }
                parser.push(&buf[..n], &mut out);
                if !out.is_empty() && !flush(std::mem::take(&mut out), &mut batch, &mut batch_bp) {
                    aborted = true;
                    break;
                }
            }
        }
        if aborted {
            return Ok(hr);
        }
        parser.finish(&mut out);
        flush(std::mem::take(&mut out), &mut batch, &mut batch_bp);
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        Ok(hr)
    });

    let mut result = Ok(());
    for batch in rx.iter() {
        if let Err(e) = on_batch(batch) {
            result = Err(e);
            break;
        }
    }
    drop(rx);
    let hr = reader
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;
    result?;
    Ok(InputInfo {
        bytes: hr.bytes,
        sha256: hr.hasher.map(|h| format!("{:x}", h.finalize())),
    })
}

fn create(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        1 << 20,
        File::create(path).with_context(|| format!("cannot create {}", path.display()))?,
    ))
}

fn out_path(prefix: &str, f: Format) -> PathBuf {
    PathBuf::from(format!("{prefix}.{}", f.extension()))
}

fn setup_threads(threads: Option<usize>) {
    if let Some(t) = threads {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(t.max(1))
            .build_global();
    }
}

fn analyze_batch(batch: &[Record], cfg: &Config) -> Vec<RecordResult> {
    batch.par_iter().map(|r| analyze_record(r, cfg)).collect()
}

/// Scan every record of `input` with `e`, merging results in input order.
fn epcr_scan(input: Input, e: &Epcr, report: &mut EpcrReport) -> Result<()> {
    let max = e.params().max_amplicons;
    stream_records(input, false, |batch| {
        let hits: Vec<_> = batch.par_iter().map(|r| e.scan(r.id(), &r.seq)).collect();
        for (r, h) in batch.iter().zip(hits) {
            report.add(r.seq.len(), h, max);
        }
        Ok(())
    })?;
    Ok(())
}

fn run(a: RunArgs) -> Result<()> {
    let t0 = Instant::now();
    setup_threads(a.threads);
    let cfg = build_config(&a)?;
    let mut formats = Vec::new();
    for f in a.formats.split(',').filter(|s| !s.trim().is_empty()) {
        match Format::parse(f) {
            Some(fm) => formats.push(fm),
            None => bail!(
                "unknown format {f:?} (expected tsv,gff3,bed,misa,statistics,json,p3in,primers,report)"
            ),
        }
    }
    if cfg.primers.is_some() {
        formats.push(Format::Primers);
    } else if formats.contains(&Format::Primers) {
        bail!("the primers format requires --primers");
    }
    formats.sort();
    formats.dedup();
    let epcr_params = if a.epcr {
        if a.input == "-" {
            bail!("--epcr reads the input twice and cannot be used with standard input");
        }
        let mut p = a.epcr_opts.params(false)?;
        // the designed products must always be within the searched range
        let pmax = cfg.primers.as_ref().map_or(0, |p| p.product_max);
        p.max_product = p.max_product.max(pmax);
        Some(p)
    } else {
        None
    };
    let prefix = a.output.clone().unwrap_or_else(|| {
        if a.input == "-" {
            "nextssr".into()
        } else {
            a.input.clone()
        }
    });
    let p3_params = cfg.primers.clone().unwrap_or_default();
    let tsv_cols = tsv::TsvColumns {
        primers: cfg.primers.is_some(),
        flanks: cfg.flank > 0,
    };

    let source = open_input(&a.input)?;
    let mut writers: Vec<(Format, BufWriter<File>)> = Vec::new();
    for &f in &formats {
        if matches!(
            f,
            Format::Statistics | Format::Json | Format::Primers | Format::Report
        ) {
            continue;
        }
        let mut w = create(&out_path(&prefix, f))?;
        match f {
            Format::Tsv => w.write_all(tsv::header(tsv_cols).as_bytes())?,
            Format::Gff3 => w.write_all(gff3::header().as_bytes())?,
            Format::Misa => w.write_all(misa::HEADER.as_bytes())?,
            _ => {}
        }
        writers.push((f, w));
    }

    let mut stats = Stats::new();
    let mut primer_rows: Vec<primers::PrimerRow> = Vec::new();
    let info = stream_records(source, !a.no_hash, |batch| {
        let results = analyze_batch(&batch, &cfg);
        for (rec, r) in batch.iter().zip(&results) {
            stats.add_record(&r.id, &rec.seq, &r.loci);
            if cfg.primers.is_some() {
                primer_rows.extend(primers::rows(r));
            }
            for (f, w) in writers.iter_mut() {
                let mut s = String::new();
                match f {
                    Format::Tsv => tsv::write_record(&mut s, r, tsv_cols),
                    Format::Gff3 => gff3::write_record(&mut s, r),
                    Format::Bed => bed::write_record(&mut s, r),
                    Format::Misa => misa::write_record(&mut s, r),
                    Format::P3in => p3in::write_record(&mut s, r, &rec.seq, &p3_params),
                    _ => {}
                }
                w.write_all(s.as_bytes())?;
            }
        }
        Ok(())
    })?;
    for (_, w) in writers.iter_mut() {
        w.flush()?;
    }

    // e-PCR specificity of the designed pairs (second pass over the input)
    let mut epcr_out: Option<(Epcr, EpcrReport)> = None;
    if let Some(ep) = epcr_params {
        let specs: Vec<PairSpec> = primer_rows.iter().map(|r| r.spec()).collect();
        let e = Epcr::new(specs, ep)?;
        let mut report = EpcrReport::new(e.pairs().len());
        if !e.pairs().is_empty() {
            epcr_scan(open_input(&a.input)?, &e, &mut report)?;
        }
        let mut s = epcr::amplicons_header(false);
        epcr::write_amplicons(&mut s, &report, e.pairs(), false);
        std::fs::write(format!("{prefix}.epcr.tsv"), s)?;
        epcr_out = Some((e, report));
    }
    if formats.contains(&Format::Primers) {
        let mut s = primers::header(epcr_out.is_some());
        primers::write(&mut s, &primer_rows, epcr_out.as_ref().map(|(_, r)| r));
        std::fs::write(out_path(&prefix, Format::Primers), s)?;
    }
    let loci_specific = epcr_out.as_ref().map(|(_, rep)| {
        let mut n = 0u64;
        let mut last: Option<(&str, usize)> = None;
        for (row, res) in primer_rows.iter().zip(&rep.results) {
            let key = (row.seq_id.as_str(), row.nr);
            if res.status(true) == epcr::Status::Specific && last != Some(key) {
                n += 1;
                last = Some(key);
            }
        }
        n
    });

    if formats.contains(&Format::Statistics) {
        let text = misa_statistics(
            &stats,
            &cfg.definition,
            cfg.max_interruption,
            &a.input,
            cfg.mode,
        );
        std::fs::write(out_path(&prefix, Format::Statistics), text)?;
    }
    let elapsed = t0.elapsed().as_secs_f64();
    if formats.contains(&Format::Report) {
        let command = std::env::args().collect::<Vec<_>>().join(" ");
        let generated = utc_now();
        let info = report::ReportInfo {
            source: &a.input,
            sha256: info.sha256.as_deref(),
            bytes: Some(info.bytes),
            command: Some(&command),
            runtime: "command line",
            threads: Some(rayon::current_num_threads()),
            seconds: Some(elapsed),
            generated: Some(&generated),
            epcr: epcr_out.as_ref().map(|(e, rep)| report::EpcrInfo {
                params: e.params(),
                pairs: e.pairs().len(),
                by_status: rep.status_counts(e.pairs()),
                loci_with_specific_pair: loci_specific.unwrap_or(0),
            }),
        };
        std::fs::write(
            out_path(&prefix, Format::Report),
            report::html(&stats, &cfg, &info),
        )?;
    }
    if formats.contains(&Format::Json) {
        let mut summary = serde_json::json!({
            "tool": "nextSSR",
            "version": nextssr_core::VERSION,
            "command": std::env::args().collect::<Vec<_>>(),
            "input": { "path": a.input, "bytes_read": info.bytes, "sha256": info.sha256 },
            "config": cfg,
            "threads": rayon::current_num_threads(),
            "runtime_seconds": (elapsed * 1000.0).round() / 1000.0,
            "stats": stats,
            "per_class": stats.per_class(),
        });
        if let Some((e, rep)) = &epcr_out {
            summary["epcr"] = serde_json::json!({
                "params": e.params(),
                "pairs": e.pairs().len(),
                "pairs_by_status": rep.status_counts(e.pairs()),
                "loci_with_specific_pair": loci_specific,
            });
        }
        std::fs::write(
            out_path(&prefix, Format::Json),
            serde_json::to_string_pretty(&summary)?,
        )?;
    }
    if !a.quiet {
        eprintln!(
            "nextSSR {}: {} sequences, {} bp, {} SSRs in {} loci ({} compound, {} imperfect){} in {:.2}s with {} threads",
            nextssr_core::VERSION,
            stats.sequences,
            stats.total_bp,
            stats.ssrs,
            stats.loci,
            stats.compounds,
            stats.imperfect,
            if cfg.primers.is_some() {
                format!(", primers for {}/{} loci", stats.primers_designed, stats.loci)
            } else {
                String::new()
            },
            elapsed,
            rayon::current_num_threads()
        );
        if let (Some((e, rep)), Some(n)) = (&epcr_out, loci_specific) {
            let c = rep.status_counts(e.pairs());
            eprintln!(
                "  e-PCR: {} pairs: {} specific, {} nonspecific; {}/{} loci have a specific pair",
                e.pairs().len(),
                c.get("specific").copied().unwrap_or(0),
                c.get("nonspecific").copied().unwrap_or(0),
                n,
                stats.primers_designed
            );
        }
        for f in &formats {
            eprintln!("  {}", out_path(&prefix, *f).display());
        }
        if epcr_out.is_some() {
            eprintln!("  {prefix}.epcr.tsv");
        }
    }
    Ok(())
}

/// Current UTC time as `YYYY-MM-DD HH:MM UTC`.
fn utc_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // civil-from-days (H. Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02} {:02}:{:02} UTC",
        rem / 3600,
        rem % 3600 / 60
    )
}

fn run_epcr(a: EpcrCmdArgs) -> Result<()> {
    let t0 = Instant::now();
    setup_threads(a.threads);
    let params = a.epcr.params(a.sequence)?;
    let specs = match (&a.primers, &a.forward, &a.reverse) {
        (Some(path), _, _) => {
            let text = std::fs::read_to_string(path)
                .with_context(|| format!("reading {}", path.display()))?;
            epcr::parse_primer_list(&text)?
        }
        (None, Some(f), Some(r)) => vec![PairSpec {
            id: "pair1".into(),
            forward: f.clone(),
            reverse: r.clone(),
            expected: None,
        }],
        _ => bail!("give a primer list with --primers, or one pair with --forward/--reverse"),
    };
    if specs.is_empty() {
        bail!("no primer pairs found");
    }
    if a.targets.iter().filter(|t| *t == "-").count() > 1 {
        bail!("standard input can be read only once");
    }
    let e = Epcr::new(specs, params)?;
    let inputs: Vec<Input> = a
        .targets
        .iter()
        .map(|t| open_input(t))
        .collect::<Result<_>>()?;
    let mut report = EpcrReport::new(e.pairs().len());
    for input in inputs {
        epcr_scan(input, &e, &mut report)?;
    }
    let mut amp = epcr::amplicons_header(a.sequence);
    epcr::write_amplicons(&mut amp, &report, e.pairs(), a.sequence);
    match &a.output {
        Some(prefix) => {
            std::fs::write(format!("{prefix}.epcr.tsv"), amp)?;
            let mut s = String::new();
            epcr::write_pairs(&mut s, &report, e.pairs());
            std::fs::write(format!("{prefix}.epcr_pairs.tsv"), s)?;
        }
        None => io::stdout().lock().write_all(amp.as_bytes())?,
    }
    if !a.quiet {
        let c = report.status_counts(e.pairs());
        let total: u64 = report.results.iter().map(|r| r.amplicons).sum();
        eprintln!(
            "nextSSR e-PCR: {} pairs vs {} sequences ({} bp): {} amplicons; {} pairs with one product, {} with several, {} with none; {:.2}s with {} threads",
            e.pairs().len(),
            report.sequences,
            report.bp,
            total,
            c.get("specific").copied().unwrap_or(0),
            c.get("nonspecific").copied().unwrap_or(0),
            c.get("no_product").copied().unwrap_or(0),
            t0.elapsed().as_secs_f64(),
            rayon::current_num_threads()
        );
        if report.results.iter().any(|r| r.saturated) {
            eprintln!("  warning: some primers have too many binding sites; their counts are lower bounds");
        }
        if let Some(prefix) = &a.output {
            eprintln!("  {prefix}.epcr.tsv\n  {prefix}.epcr_pairs.tsv");
        }
    }
    Ok(())
}

fn run_misa(a: MisaArgs) -> Result<()> {
    setup_threads(a.threads);
    let ini_text = match &a.ini {
        Some(p) => std::fs::read_to_string(p).with_context(|| {
            format!(
                "Userdefined specifications file doesn't exist: {}",
                p.display()
            )
        })?,
        None => match std::fs::read_to_string("misa.ini") {
            Ok(t) => t,
            Err(_) => {
                eprintln!("nextssr misa: ./misa.ini not found, using MISA defaults (1-10 2-6 3-5 4-5 5-5 6-5, interruptions 100)");
                DEFAULT_INI.to_string()
            }
        },
    };
    let ini = MisaIni::parse(&ini_text)?;
    let cfg = Config::misa(&ini);
    cfg.validate()?;
    if !Path::new(&a.fasta).exists() {
        bail!("FASTA file doesn't exist: {}", a.fasta);
    }
    let source = open_input(&a.fasta)?;
    let mut w = create(Path::new(&format!("{}.misa", a.fasta)))?;
    w.write_all(misa::HEADER.as_bytes())?;
    let mut gff = if ini.gff {
        let mut g = create(Path::new(&format!("{}.gff3", a.fasta)))?;
        g.write_all(gff3::header().as_bytes())?;
        Some(g)
    } else {
        None
    };
    let mut stats = Stats::new();
    stream_records(source, false, |batch| {
        for (rec, r) in batch.iter().zip(analyze_batch(&batch, &cfg)) {
            stats.add_record(&r.id, &rec.seq, &r.loci);
            let mut s = String::new();
            misa::write_record(&mut s, &r);
            w.write_all(s.as_bytes())?;
            if let Some(g) = gff.as_mut() {
                let mut s = String::new();
                gff3::write_record(&mut s, &r);
                g.write_all(s.as_bytes())?;
            }
        }
        Ok(())
    })?;
    w.flush()?;
    if let Some(g) = gff.as_mut() {
        g.flush()?;
    }
    let text = misa_statistics(
        &stats,
        &cfg.definition,
        cfg.max_interruption,
        &a.fasta,
        Mode::Misa,
    );
    std::fs::write(format!("{}.statistics", a.fasta), text)?;
    Ok(())
}

fn run_tm(a: TmArgs) -> Result<()> {
    let p = a.thermo.params();
    for s in &a.sequences {
        match tm(s.as_bytes(), &p) {
            Some(t) => println!("{s}\t{t:.2}"),
            None => println!("{s}\tNA"),
        }
    }
    Ok(())
}

fn main() {
    let cli = Cli::parse();
    let res = match cli.command {
        Command::Run(a) => run(*a),
        Command::Misa(a) => run_misa(a),
        Command::InitIni => {
            print!("{DEFAULT_INI}");
            Ok(())
        }
        Command::Tm(a) => run_tm(a),
        Command::Epcr(a) => run_epcr(*a),
    };
    if let Err(e) = res {
        eprintln!("nextssr: error: {e:#}");
        std::process::exit(1);
    }
}
