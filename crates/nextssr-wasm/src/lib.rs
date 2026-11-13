//! WebAssembly bindings: streaming SSR analysis inside the browser.
//!
//! ```js
//! const a = new Analyzer(JSON.stringify(options), "genome.fa");
//! for await (const chunk of file.stream()) a.push(chunk);
//! const summary = JSON.parse(a.finish());
//! const tsv = a.output("tsv");
//! ```

use std::io::Write;

use flate2::write::MultiGzDecoder;
use serde::Deserialize;
use wasm_bindgen::prelude::*;

use nextssr_core::output::{bed, gff3, misa, p3in, tsv};
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
    flank: Option<usize>,
    preview_limit: Option<usize>,
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

/// Streaming analyzer.
#[wasm_bindgen]
pub struct Analyzer {
    cfg: Config,
    source: String,
    parser: FastaParser,
    records: Vec<Record>,
    stats: Stats,
    tsv: String,
    gff3: String,
    bed: String,
    misa: String,
    p3in: String,
    preview: Vec<serde_json::Value>,
    preview_limit: usize,
    finished: bool,
    /// `None` until the first bytes reveal whether the input is gzip.
    gzip: Option<Option<MultiGzDecoder<Vec<u8>>>>,
    head: Vec<u8>,
    bytes_in: u64,
}

#[wasm_bindgen]
impl Analyzer {
    /// Create an analyzer from a JSON options string.
    #[wasm_bindgen(constructor)]
    pub fn new(options_json: &str, source_name: &str) -> Result<Analyzer, JsError> {
        let o: Options = if options_json.trim().is_empty() {
            Options::default()
        } else {
            serde_json::from_str(options_json)
                .map_err(|e| JsError::new(&format!("invalid options: {e}")))?
        };
        let cfg = build_config(&o).map_err(|e| JsError::new(&e))?;
        let cols = tsv::TsvColumns {
            primers: cfg.primers.is_some(),
            flanks: cfg.flank > 0,
        };
        Ok(Analyzer {
            source: source_name.to_string(),
            parser: FastaParser::new(),
            records: Vec::new(),
            stats: Stats::new(),
            tsv: tsv::header(cols),
            gff3: gff3::header(),
            bed: String::new(),
            misa: misa::HEADER.to_string(),
            p3in: String::new(),
            preview: Vec::new(),
            preview_limit: o.preview_limit.unwrap_or(2000),
            finished: false,
            gzip: None,
            head: Vec::new(),
            bytes_in: 0,
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

    fn feed(&mut self, data: &[u8]) {
        self.parser.push(data, &mut self.records);
        if !self.records.is_empty() {
            self.process();
        }
    }

    /// Feed a chunk of FASTA bytes (plain, gzip or bgzip; detected automatically).
    pub fn push(&mut self, chunk: &[u8]) -> Result<(), JsError> {
        self.bytes_in += chunk.len() as u64;
        if self.gzip.is_none() {
            self.head.extend_from_slice(chunk);
            if self.head.len() < 2 {
                return Ok(());
            }
            let head = std::mem::take(&mut self.head);
            let is_gz = head[0] == 0x1f && head[1] == 0x8b;
            self.gzip = Some(is_gz.then(|| MultiGzDecoder::new(Vec::new())));
            return self.push_inner(&head);
        }
        self.push_inner(chunk)
    }

    fn push_inner(&mut self, data: &[u8]) -> Result<(), JsError> {
        let mut decoded = None;
        if let Some(Some(dec)) = self.gzip.as_mut() {
            dec.write_all(data)
                .map_err(|e| JsError::new(&format!("gzip: {e}")))?;
            decoded = Some(std::mem::take(dec.get_mut()));
        }
        match decoded {
            Some(buf) => self.feed(&buf),
            None => self.feed(data),
        }
        Ok(())
    }

    /// Compressed (raw) bytes received so far.
    pub fn bytes_in(&self) -> f64 {
        self.bytes_in as f64
    }

    /// Bytes buffered for the sequence currently being read.
    pub fn buffered(&self) -> usize {
        self.parser.buffered_bytes()
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
        if !self.finished {
            if self.gzip.is_none() && !self.head.is_empty() {
                let head = std::mem::take(&mut self.head);
                self.gzip = Some(None);
                self.feed(&head);
            }
            if let Some(Some(dec)) = self.gzip.as_mut() {
                dec.try_finish()
                    .map_err(|e| JsError::new(&format!("gzip: {e}")))?;
                let buf = std::mem::take(dec.get_mut());
                self.feed(&buf);
            }
            self.parser.finish(&mut self.records);
            self.process();
            self.finished = true;
        }
        Ok(self.summary_json())
    }

    fn summary_json(&self) -> String {
        serde_json::to_string_pretty(&serde_json::json!({
            "tool": "nextSSR",
            "version": nextssr_core::VERSION,
            "runtime": "webassembly",
            "input": { "name": self.source },
            "config": self.cfg,
            "stats": self.stats,
            "per_class": self.stats.per_class(),
        }))
        .unwrap_or_default()
    }

    /// Output text for `tsv`, `gff3`, `bed`, `misa`, `statistics`, `json` or `p3in`.
    pub fn output(&self, format: &str) -> Result<String, JsError> {
        Ok(match format {
            "tsv" => self.tsv.clone(),
            "gff3" | "gff" => self.gff3.clone(),
            "bed" => self.bed.clone(),
            "misa" => self.misa.clone(),
            "p3in" => self.p3in.clone(),
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
