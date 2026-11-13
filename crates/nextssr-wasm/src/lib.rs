//! WebAssembly bindings: streaming SSR analysis inside the browser.
//!
//! ```js
//! const a = new Analyzer(JSON.stringify(options), "genome.fa");
//! for await (const chunk of file.stream()) a.push(chunk);
//! const summary = JSON.parse(a.finish());
//! const tsv = a.output("tsv");
//!
//! // specificity of the designed primers (second pass over the same file)
//! const e = a.epcr_scanner();
//! for await (const chunk of file.stream()) e.push(chunk);
//! e.finish();
//! a.set_epcr(e);
//! const primers = a.output("primers");
//!
//! // or any primer list against any FASTA
//! const s = EpcrScanner.from_list(primerText, JSON.stringify(options));
//! ```

use std::io::Write;

use flate2::write::MultiGzDecoder;
use serde::Deserialize;
use wasm_bindgen::prelude::*;

use nextssr_core::epcr::{self, Epcr, EpcrParams, EpcrReport, PairSpec, Status};
use nextssr_core::output::{bed, gff3, misa, p3in, primers, tsv};
use nextssr_core::primer::PrimerParams;
use nextssr_core::stats::misa_statistics;
use nextssr_core::thermo::ThermoParams;
use nextssr_core::{
    analyze_record, Config, Definition, FastaParser, ImperfectParams, Mode, Record, RecordResult,
    Stats,
};

/// Options accepted from JavaScript (all optional).
#[derive(Deserialize, Default)]
#[serde(default)]
struct Options {
    definition: Option<String>,
    mode: Option<String>,
    max_interruption: Option<i64>,
    compound: Option<bool>,
    imperfect: bool,
    max_gap: Option<usize>,
    min_purity: Option<f64>,
    primers: bool,
    product_min: Option<usize>,
    product_max: Option<usize>,
    tm_min: Option<f64>,
    tm_opt: Option<f64>,
    tm_max: Option<f64>,
    size_min: Option<usize>,
    size_opt: Option<usize>,
    size_max: Option<usize>,
    gc_min: Option<f64>,
    gc_max: Option<f64>,
    max_tm_diff: Option<f64>,
    num_pairs: Option<usize>,
    flank: Option<usize>,
    preview_limit: Option<usize>,
    // e-PCR
    mismatches: Option<usize>,
    three_prime_len: Option<usize>,
    three_prime_mismatches: Option<usize>,
    epcr_min_product: Option<usize>,
    epcr_max_product: Option<usize>,
    max_amplicons: Option<usize>,
    amplicon_sequence: bool,
}

fn epcr_params(o: &Options, min_max_product: usize) -> Result<EpcrParams, String> {
    let d = EpcrParams::default();
    let p = EpcrParams {
        max_mismatches: o.mismatches.unwrap_or(d.max_mismatches),
        three_prime_len: o.three_prime_len.unwrap_or(d.three_prime_len),
        max_three_prime_mismatches: o
            .three_prime_mismatches
            .unwrap_or(d.max_three_prime_mismatches),
        min_product: o.epcr_min_product.unwrap_or(d.min_product),
        max_product: o
            .epcr_max_product
            .unwrap_or(d.max_product)
            .max(min_max_product),
        max_amplicons: o.max_amplicons.unwrap_or(d.max_amplicons),
        with_sequence: o.amplicon_sequence,
        ..d
    };
    p.validate().map_err(|e| e.to_string())?;
    Ok(p)
}

fn parse_options(json: &str) -> Result<Options, JsError> {
    if json.trim().is_empty() {
        Ok(Options::default())
    } else {
        serde_json::from_str(json).map_err(|e| JsError::new(&format!("invalid options: {e}")))
    }
}

fn build_config(o: &Options) -> Result<Config, String> {
    let definition = match &o.definition {
        Some(d) if !d.trim().is_empty() => Definition::parse(d).map_err(|e| e.to_string())?,
        _ => Definition::misa_default(),
    };
    let mode = match o.mode.as_deref() {
        Some("misa") => Mode::Misa,
        _ => Mode::Native,
    };
    let imperfect = if o.imperfect && mode == Mode::Native {
        let mut ip = ImperfectParams::default();
        if let Some(g) = o.max_gap {
            ip.max_gap = g;
        }
        if let Some(p) = o.min_purity {
            ip.min_purity = p;
        }
        Some(ip)
    } else {
        None
    };
    let primers = if o.primers {
        let d = PrimerParams::default();
        Some(PrimerParams {
            product_min: o.product_min.unwrap_or(d.product_min),
            product_max: o.product_max.unwrap_or(d.product_max),
            min_tm: o.tm_min.unwrap_or(d.min_tm),
            opt_tm: o.tm_opt.unwrap_or(d.opt_tm),
            max_tm: o.tm_max.unwrap_or(d.max_tm),
            min_size: o.size_min.unwrap_or(d.min_size),
            opt_size: o.size_opt.unwrap_or(d.opt_size),
            max_size: o.size_max.unwrap_or(d.max_size),
            min_gc: o.gc_min.unwrap_or(d.min_gc),
            max_gc: o.gc_max.unwrap_or(d.max_gc),
            max_tm_diff: o.max_tm_diff.unwrap_or(d.max_tm_diff),
            num_return: o.num_pairs.unwrap_or(d.num_return),
            ..d
        })
    } else {
        None
    };
    let cfg = Config {
        definition,
        mode,
        max_interruption: if o.compound == Some(false) {
            None
        } else {
            Some(o.max_interruption.unwrap_or(100))
        },
        imperfect,
        primers,
        flank: o.flank.unwrap_or(0),
    };
    cfg.validate().map_err(|e| e.to_string())?;
    Ok(cfg)
}

/// Incremental FASTA input: plain or gzip/bgzip (detected from the first bytes).
struct Input {
    parser: FastaParser,
    /// `None` until the first bytes reveal whether the input is gzip.
    gzip: Option<Option<MultiGzDecoder<Vec<u8>>>>,
    head: Vec<u8>,
    bytes_in: u64,
    finished: bool,
}

impl Input {
    fn new() -> Self {
        Self {
            parser: FastaParser::new(),
            gzip: None,
            head: Vec::new(),
            bytes_in: 0,
            finished: false,
        }
    }

    /// Feed raw bytes; complete records are appended to `out`.
    fn push(&mut self, chunk: &[u8], out: &mut Vec<Record>) -> Result<(), JsError> {
        self.bytes_in += chunk.len() as u64;
        if self.gzip.is_none() {
            self.head.extend_from_slice(chunk);
            if self.head.len() < 2 {
                return Ok(());
            }
            let head = std::mem::take(&mut self.head);
            let is_gz = head[0] == 0x1f && head[1] == 0x8b;
            self.gzip = Some(is_gz.then(|| MultiGzDecoder::new(Vec::new())));
            return self.decode(&head, out);
        }
        self.decode(chunk, out)
    }

    fn decode(&mut self, data: &[u8], out: &mut Vec<Record>) -> Result<(), JsError> {
        match self.gzip.as_mut() {
            Some(Some(dec)) => {
                dec.write_all(data)
                    .map_err(|e| JsError::new(&format!("gzip: {e}")))?;
                let buf = std::mem::take(dec.get_mut());
                self.parser.push(&buf, out);
            }
            _ => self.parser.push(data, out),
        }
        Ok(())
    }

    /// Flush the last record (idempotent).
    fn finish(&mut self, out: &mut Vec<Record>) -> Result<(), JsError> {
        if self.finished {
            return Ok(());
        }
        if self.gzip.is_none() && !self.head.is_empty() {
            let head = std::mem::take(&mut self.head);
            self.gzip = Some(None);
            self.parser.push(&head, out);
        }
        if let Some(Some(dec)) = self.gzip.as_mut() {
            dec.try_finish()
                .map_err(|e| JsError::new(&format!("gzip: {e}")))?;
            let buf = std::mem::take(dec.get_mut());
            self.parser.push(&buf, out);
        }
        self.parser.finish(out);
        self.finished = true;
        Ok(())
    }
}

/// Streaming analyzer.
#[wasm_bindgen]
pub struct Analyzer {
    cfg: Config,
    opts: Options,
    source: String,
    input: Input,
    records: Vec<Record>,
    stats: Stats,
    tsv: String,
    gff3: String,
    bed: String,
    misa: String,
    p3in: String,
    primer_rows: Vec<primers::PrimerRow>,
    epcr: Option<(Vec<PairSpec>, EpcrParams, EpcrReport)>,
    preview: Vec<serde_json::Value>,
    preview_limit: usize,
}

#[wasm_bindgen]
impl Analyzer {
    /// Create an analyzer from a JSON options string.
    #[wasm_bindgen(constructor)]
    pub fn new(options_json: &str, source_name: &str) -> Result<Analyzer, JsError> {
        let o = parse_options(options_json)?;
        let cfg = build_config(&o).map_err(|e| JsError::new(&e))?;
        let cols = tsv::TsvColumns {
            primers: cfg.primers.is_some(),
            flanks: cfg.flank > 0,
        };
        Ok(Analyzer {
            source: source_name.to_string(),
            input: Input::new(),
            records: Vec::new(),
            stats: Stats::new(),
            tsv: tsv::header(cols),
            gff3: gff3::header(),
            bed: String::new(),
            misa: misa::HEADER.to_string(),
            p3in: String::new(),
            primer_rows: Vec::new(),
            epcr: None,
            preview: Vec::new(),
            preview_limit: o.preview_limit.unwrap_or(2000),
            opts: o,
            cfg,
        })
    }

    fn process(&mut self) {
        let cols = tsv::TsvColumns {
            primers: self.cfg.primers.is_some(),
            flanks: self.cfg.flank > 0,
        };
        let p3 = self.cfg.primers.clone().unwrap_or_default();
        for rec in std::mem::take(&mut self.records) {
            let r: RecordResult = analyze_record(&rec, &self.cfg);
            self.stats.add_record(r.length, &r.loci);
            tsv::write_record(&mut self.tsv, &r, cols);
            gff3::write_record(&mut self.gff3, &r);
            bed::write_record(&mut self.bed, &r);
            misa::write_record(&mut self.misa, &r);
            p3in::write_record(&mut self.p3in, &r, &rec.seq, &p3);
            if self.cfg.primers.is_some() {
                self.primer_rows.extend(primers::rows(&r));
            }
            for l in &r.loci {
                if self.preview.len() >= self.preview_limit {
                    break;
                }
                let mut v = serde_json::to_value(l).unwrap_or_default();
                if let Some(obj) = v.as_object_mut() {
                    obj.insert("seq_id".into(), r.id.clone().into());
                    obj.remove("members");
                }
                self.preview.push(v);
            }
        }
    }

    /// Feed a chunk of FASTA bytes (plain, gzip or bgzip; detected automatically).
    pub fn push(&mut self, chunk: &[u8]) -> Result<(), JsError> {
        self.input.push(chunk, &mut self.records)?;
        if !self.records.is_empty() {
            self.process();
        }
        Ok(())
    }

    /// Compressed (raw) bytes received so far.
    pub fn bytes_in(&self) -> f64 {
        self.input.bytes_in as f64
    }

    /// Bytes buffered for the sequence currently being read.
    pub fn buffered(&self) -> usize {
        self.input.parser.buffered_bytes()
    }

    /// Progress as JSON: `{"sequences":..,"bp":..,"ssrs":..}`.
    pub fn progress(&self) -> String {
        serde_json::json!({
            "sequences": self.stats.sequences,
            "bp": self.stats.total_bp,
            "ssrs": self.stats.ssrs,
        })
        .to_string()
    }

    /// Flush the last record and return the summary as JSON.
    pub fn finish(&mut self) -> Result<String, JsError> {
        self.input.finish(&mut self.records)?;
        self.process();
        Ok(self.summary_json())
    }

    fn loci_with_specific_pair(&self, rep: &EpcrReport) -> u64 {
        let mut n = 0;
        let mut last: Option<(&str, usize)> = None;
        for (row, res) in self.primer_rows.iter().zip(&rep.results) {
            let key = (row.seq_id.as_str(), row.nr);
            if res.status(true) == Status::Specific && last != Some(key) {
                n += 1;
                last = Some(key);
            }
        }
        n
    }

    fn summary_json(&self) -> String {
        let mut v = serde_json::json!({
            "tool": "nextSSR",
            "version": nextssr_core::VERSION,
            "runtime": "webassembly",
            "input": { "name": self.source },
            "config": self.cfg,
            "stats": self.stats,
            "per_class": self.stats.per_class(),
        });
        if let Some((specs, params, rep)) = &self.epcr {
            v["epcr"] = serde_json::json!({
                "params": params,
                "pairs": specs.len(),
                "pairs_by_status": rep.status_counts(specs),
                "loci_with_specific_pair": self.loci_with_specific_pair(rep),
            });
        }
        serde_json::to_string_pretty(&v).unwrap_or_default()
    }

    /// Summary JSON (after `finish`, and again after `set_epcr`).
    pub fn summary(&self) -> String {
        self.summary_json()
    }

    /// e-PCR scanner for the designed pairs (call after `finish`); stream the
    /// same input through it, then pass it to `set_epcr`.
    pub fn epcr_scanner(&self) -> Result<EpcrScanner, JsError> {
        let pmax = self.cfg.primers.as_ref().map_or(0, |p| p.product_max);
        let params = epcr_params(&self.opts, pmax).map_err(|e| JsError::new(&e))?;
        let specs = self.primer_rows.iter().map(|r| r.spec()).collect();
        EpcrScanner::build(specs, params)
    }

    /// Attach the results of a finished scanner from `epcr_scanner`.
    pub fn set_epcr(&mut self, scanner: &EpcrScanner) -> Result<(), JsError> {
        if scanner.report.results.len() != self.primer_rows.len() {
            return Err(JsError::new(
                "e-PCR scanner does not belong to this analysis",
            ));
        }
        self.epcr = Some((
            scanner.engine.pairs().to_vec(),
            scanner.engine.params().clone(),
            scanner.report.clone(),
        ));
        Ok(())
    }

    /// Output text for `tsv`, `gff3`, `bed`, `misa`, `statistics`, `json`,
    /// `p3in`, `primers` or `epcr` (amplicons of the designed pairs).
    pub fn output(&self, format: &str) -> Result<String, JsError> {
        Ok(match format {
            "tsv" => self.tsv.clone(),
            "gff3" | "gff" => self.gff3.clone(),
            "bed" => self.bed.clone(),
            "misa" => self.misa.clone(),
            "p3in" => self.p3in.clone(),
            "primers" => {
                let rep = self.epcr.as_ref().map(|(_, _, r)| r);
                let mut s = primers::header(rep.is_some());
                primers::write(&mut s, &self.primer_rows, rep);
                s
            }
            "epcr" => match &self.epcr {
                Some((specs, params, rep)) => {
                    let mut s = epcr::amplicons_header(params.with_sequence);
                    epcr::write_amplicons(&mut s, rep, specs, params.with_sequence);
                    s
                }
                None => return Err(JsError::new("e-PCR has not been run")),
            },
            "statistics" | "stats" => misa_statistics(
                &self.stats,
                &self.cfg.definition,
                self.cfg.max_interruption,
                &self.source,
                self.cfg.mode,
            ),
            "json" => self.summary_json(),
            other => return Err(JsError::new(&format!("unknown format {other}"))),
        })
    }

    /// First loci (up to `preview_limit`) as a JSON array for display.
    pub fn preview(&self) -> String {
        serde_json::to_string(&self.preview).unwrap_or_else(|_| "[]".into())
    }
}

/// Streaming in-silico PCR.
#[wasm_bindgen]
pub struct EpcrScanner {
    engine: Epcr,
    report: EpcrReport,
    input: Input,
    records: Vec<Record>,
}

impl EpcrScanner {
    fn build(specs: Vec<PairSpec>, params: EpcrParams) -> Result<EpcrScanner, JsError> {
        let engine = Epcr::new(specs, params).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(EpcrScanner {
            report: EpcrReport::new(engine.pairs().len()),
            engine,
            input: Input::new(),
            records: Vec::new(),
        })
    }

    fn process(&mut self) {
        let max = self.engine.params().max_amplicons;
        for rec in std::mem::take(&mut self.records) {
            let hits = self.engine.scan(rec.id(), &rec.seq);
            self.report.add(rec.seq.len(), hits, max);
        }
    }
}

#[wasm_bindgen]
impl EpcrScanner {
    /// Scanner for a primer list (see `nextssr epcr --primers` for the formats)
    /// with e-PCR options (`mismatches`, `three_prime_len`,
    /// `three_prime_mismatches`, `epcr_min_product`, `epcr_max_product`,
    /// `max_amplicons`, `amplicon_sequence`).
    pub fn from_list(primers_text: &str, options_json: &str) -> Result<EpcrScanner, JsError> {
        let o = parse_options(options_json)?;
        let params = epcr_params(&o, 0).map_err(|e| JsError::new(&e))?;
        let specs =
            epcr::parse_primer_list(primers_text).map_err(|e| JsError::new(&e.to_string()))?;
        if specs.is_empty() {
            return Err(JsError::new("no primer pairs found"));
        }
        Self::build(specs, params)
    }

    /// Number of primer pairs.
    pub fn pairs(&self) -> usize {
        self.engine.pairs().len()
    }

    /// Feed a chunk of FASTA bytes (plain, gzip or bgzip).
    pub fn push(&mut self, chunk: &[u8]) -> Result<(), JsError> {
        self.input.push(chunk, &mut self.records)?;
        if !self.records.is_empty() {
            self.process();
        }
        Ok(())
    }

    /// Progress as JSON: `{"sequences":..,"bp":..}`.
    pub fn progress(&self) -> String {
        serde_json::json!({ "sequences": self.report.sequences, "bp": self.report.bp }).to_string()
    }

    /// Flush the last record and return the summary as JSON, with one entry per
    /// pair: `{id, forward, reverse, status, amplicons, perfect, saturated}`.
    pub fn finish(&mut self) -> Result<String, JsError> {
        self.input.finish(&mut self.records)?;
        self.process();
        Ok(self.summary())
    }

    pub fn summary(&self) -> String {
        let pairs: Vec<serde_json::Value> = self
            .engine
            .pairs()
            .iter()
            .zip(&self.report.results)
            .map(|(p, r)| {
                serde_json::json!({
                    "id": p.id,
                    "forward": p.forward,
                    "reverse": p.reverse,
                    "status": r.status(p.expected.is_some()),
                    "amplicons": r.amplicons,
                    "perfect": r.perfect,
                    "saturated": r.saturated,
                })
            })
            .collect();
        serde_json::json!({
            "sequences": self.report.sequences,
            "bp": self.report.bp,
            "params": self.engine.params(),
            "pairs_by_status": self.report.status_counts(self.engine.pairs()),
            "pairs": pairs,
        })
        .to_string()
    }

    /// `amplicons` (one line per product) or `pairs` (one line per pair).
    pub fn output(&self, format: &str) -> Result<String, JsError> {
        let ws = self.engine.params().with_sequence;
        Ok(match format {
            "amplicons" | "epcr" => {
                let mut s = epcr::amplicons_header(ws);
                epcr::write_amplicons(&mut s, &self.report, self.engine.pairs(), ws);
                s
            }
            "pairs" => {
                let mut s = String::new();
                epcr::write_pairs(&mut s, &self.report, self.engine.pairs());
                s
            }
            other => return Err(JsError::new(&format!("unknown e-PCR output {other}"))),
        })
    }
}

/// Engine version.
#[wasm_bindgen]
pub fn version() -> String {
    nextssr_core::VERSION.to_string()
}

/// Melting temperature (°C) under Primer3 default conditions; NaN if invalid.
#[wasm_bindgen]
pub fn melting_temperature(seq: &str) -> f64 {
    nextssr_core::thermo::tm(seq.as_bytes(), &ThermoParams::default()).unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analyzer_end_to_end() {
        let mut a = Analyzer::new(r#"{"primers": false, "imperfect": true}"#, "x.fa")
            .ok()
            .unwrap();
        let fa = b">s1\nGGGACACACACACACACTTTGCAAGAAGAAGAAGAAG\n>s2\nAAAAAAAAAAAA\n";
        for c in fa.chunks(7) {
            a.push(c).ok().unwrap();
        }
        let summary: serde_json::Value = serde_json::from_str(&a.finish().ok().unwrap()).unwrap();
        assert_eq!(summary["stats"]["sequences"], 2);
        assert_eq!(summary["stats"]["ssrs"], 3);
        assert!(a.output("misa").ok().unwrap().contains("(AC)7"));
        assert_eq!(
            serde_json::from_str::<Vec<serde_json::Value>>(&a.preview())
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn primers_and_epcr() {
        let mut a = Analyzer::new(r#"{"primers": true, "num_pairs": 2}"#, "cp.fa")
            .ok()
            .unwrap();
        let fa = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/NC_000932.1_Arabidopsis_cp.fasta"
        ))
        .unwrap();
        for c in fa.chunks(1 << 16) {
            a.push(c).ok().unwrap();
        }
        a.finish().ok().unwrap();
        let mut e = a.epcr_scanner().ok().unwrap();
        for c in fa.chunks(1 << 16) {
            e.push(c).ok().unwrap();
        }
        let s: serde_json::Value = serde_json::from_str(&e.finish().ok().unwrap()).unwrap();
        assert!(s["pairs_by_status"]["specific"].as_u64().unwrap() > 50);
        a.set_epcr(&e).ok().unwrap();
        let t = a.output("primers").ok().unwrap();
        assert!(t.lines().next().unwrap().ends_with("\trecommended"));
        assert!(t.contains("\tspecific\t1\t1\tyes"));
        assert!(a.output("epcr").ok().unwrap().lines().count() > 50);
        let summary: serde_json::Value = serde_json::from_str(&a.summary()).unwrap();
        assert!(summary["epcr"]["loci_with_specific_pair"].as_u64().unwrap() > 40);

        // standalone list against the same sequence
        let mut e = EpcrScanner::from_list(&t, r#"{"mismatches": 1}"#)
            .ok()
            .unwrap();
        e.push(&fa).ok().unwrap();
        let s: serde_json::Value = serde_json::from_str(&e.finish().ok().unwrap()).unwrap();
        assert_eq!(s["pairs"].as_array().unwrap().len(), e.pairs());
    }

    #[test]
    fn gzip_multi_member_input() {
        use flate2::{write::GzEncoder, Compression};
        let gz = |d: &[u8]| {
            let mut e = GzEncoder::new(Vec::new(), Compression::default());
            e.write_all(d).unwrap();
            e.finish().unwrap()
        };
        // two concatenated members, like bgzip
        let mut data = gz(b">s1\nGGGACACACACACACAC");
        data.extend(gz(b"ACTTT\n>s2\nCCCAGAGAGAGAGAGTT\n"));
        let mut a = Analyzer::new("{}", "x.fa.gz").ok().unwrap();
        for c in data.chunks(5) {
            a.push(c).ok().unwrap();
        }
        let summary: serde_json::Value = serde_json::from_str(&a.finish().ok().unwrap()).unwrap();
        assert_eq!(summary["stats"]["sequences"], 2);
        let m = a.output("misa").ok().unwrap();
        assert!(m.contains("(AC)8") && m.contains("(AG)6"), "{m}");
    }
}
