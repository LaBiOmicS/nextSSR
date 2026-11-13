//! Search parameters, `misa.ini` parsing and run configuration.

use std::collections::BTreeMap;
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::primer::PrimerParams;

/// Largest supported repeat unit (motif length).
pub const MAX_UNIT: usize = 12;

/// Error raised for invalid parameters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParamError(pub String);

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ParamError {}

fn err<T>(msg: impl Into<String>) -> Result<T, ParamError> {
    Err(ParamError(msg.into()))
}

/// Detection semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    /// Exact `misa.pl` v2.0 behaviour (regex masking, MISA compound notation).
    Misa,
    /// nextSSR semantics: all maximal perfect repeats, correct compound notation.
    #[default]
    Native,
}

/// SSR definition: minimum number of repeats for each unit size
/// (MISA `definition(unit_size,min_repeats)`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Definition {
    /// unit size -> minimum number of repeats
    pub min_repeats: BTreeMap<usize, usize>,
}

impl Default for Definition {
    fn default() -> Self {
        Self::misa_default()
    }
}

impl Definition {
    /// MISA default: `1-10 2-6 3-5 4-5 5-5 6-5`.
    pub fn misa_default() -> Self {
        Self::from_pairs(&[(1, 10), (2, 6), (3, 5), (4, 5), (5, 5), (6, 5)]).expect("valid default")
    }

    /// Build from `(unit, min_repeats)` pairs (later pairs override earlier ones).
    pub fn from_pairs(pairs: &[(usize, usize)]) -> Result<Self, ParamError> {
        let mut min_repeats = BTreeMap::new();
        for &(u, r) in pairs {
            min_repeats.insert(u, r);
        }
        let d = Self { min_repeats };
        d.validate()?;
        Ok(d)
    }

    /// Parse `"1-10 2-6 3-5"` (separators: whitespace, commas or semicolons).
    pub fn parse(s: &str) -> Result<Self, ParamError> {
        let nums = extract_numbers(s);
        if nums.is_empty() {
            return err(format!("empty SSR definition: {s:?}"));
        }
        if nums.len() % 2 != 0 {
            return err(format!(
                "SSR definition must contain unit-min_repeats pairs, got {s:?}"
            ));
        }
        let pairs: Vec<(usize, usize)> = nums.chunks(2).map(|c| (c[0], c[1])).collect();
        Self::from_pairs(&pairs)
    }

    /// Check unit and repeat ranges.
    pub fn validate(&self) -> Result<(), ParamError> {
        if self.min_repeats.is_empty() {
            return err("SSR definition is empty");
        }
        for (&u, &r) in &self.min_repeats {
            if u == 0 || u > MAX_UNIT {
                return err(format!("unit size {u} out of range 1..={MAX_UNIT}"));
            }
            if r < 2 {
                return err(format!(
                    "minimum repeats for unit {u} must be >= 2 (got {r})"
                ));
            }
        }
        Ok(())
    }

    /// Unit sizes in ascending order.
    pub fn units(&self) -> impl Iterator<Item = usize> + '_ {
        self.min_repeats.keys().copied()
    }

    /// Minimum repeats for a unit size (if defined).
    pub fn get(&self, unit: usize) -> Option<usize> {
        self.min_repeats.get(&unit).copied()
    }

    /// Smallest minimum-repeat threshold (MISA `$min_repeats`).
    pub fn min_threshold(&self) -> usize {
        self.min_repeats.values().copied().min().unwrap_or(1000)
    }

    /// MISA textual form, e.g. `1-10 2-6 3-5 4-5 5-5 6-5`.
    pub fn to_misa_string(&self) -> String {
        self.min_repeats
            .iter()
            .map(|(u, r)| format!("{u}-{r}"))
            .collect::<Vec<_>>()
            .join(" ")
    }
}

/// All unsigned integers appearing in a string (Perl `/(\d+)/g`).
fn extract_numbers(s: &str) -> Vec<usize> {
    let mut out = Vec::new();
    let mut cur: Option<usize> = None;
    for c in s.chars() {
        if let Some(d) = c.to_digit(10) {
            cur = Some(
                cur.unwrap_or(0)
                    .saturating_mul(10)
                    .saturating_add(d as usize),
            );
        } else if let Some(v) = cur.take() {
            out.push(v);
        }
    }
    if let Some(v) = cur {
        out.push(v);
    }
    out
}

/// Parsed `misa.ini` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MisaIni {
    pub definition: Definition,
    /// `interruptions(max_difference_between_2_SSRs)`; MISA uses 0 when absent.
    pub max_interruption: i64,
    pub gff: bool,
}

impl MisaIni {
    /// Parse with the same rules as `misa.pl`: lines beginning (case-insensitively)
    /// with `def`, `int` or `GFF`, followed by a non-space word and whitespace.
    pub fn parse(text: &str) -> Result<Self, ParamError> {
        let mut definition: Option<Definition> = None;
        let mut max_interruption = 0i64;
        let mut gff = false;
        for line in text.lines() {
            let lower = line.to_ascii_lowercase();
            // `^def\S*\s+(.*)`
            if let Some(rest) = keyword_value(line, &lower, "def") {
                let nums = extract_numbers(rest);
                let mut map = BTreeMap::new();
                for c in nums.chunks(2) {
                    // Perl hash from an odd list leaves the last key undefined; ignore it.
                    if c.len() == 2 {
                        map.insert(c[0], c[1]);
                    }
                }
                let d = Definition { min_repeats: map };
                d.validate()?;
                definition = Some(d);
            }
            // `^int\S*\s+(\d+)`
            if let Some(rest) = keyword_value(line, &lower, "int") {
                let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
                if !digits.is_empty() {
                    max_interruption = digits.parse().unwrap_or(0);
                }
            }
            // `^GFF\S*\s+true`
            if let Some(rest) = keyword_value(line, &lower, "gff") {
                if rest.len() >= 4 && rest[..4].eq_ignore_ascii_case("true") {
                    gff = true;
                }
            }
        }
        let Some(definition) = definition else {
            return err("misa.ini has no 'definition' line");
        };
        Ok(Self {
            definition,
            max_interruption,
            gff,
        })
    }
}

/// Emulates `^kw\S*\s+(.*)` (case-insensitive) and returns the captured value.
fn keyword_value<'a>(line: &'a str, lower: &str, kw: &str) -> Option<&'a str> {
    if !lower.starts_with(kw) {
        return None;
    }
    let bytes = line.as_bytes();
    let mut i = kw.len();
    while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    let ws_start = i;
    while i < bytes.len() && bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    if i == ws_start {
        return None;
    }
    Some(&line[i..])
}

/// Parameters for imperfect (interrupted) SSR detection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ImperfectParams {
    /// Minimum copies of a perfect seed, per unit size. Default: at least two
    /// copies and at least 8 bp (1-8 2-4 3-3 4-2 5-2 6-2 ...).
    pub seed_repeats: BTreeMap<usize, usize>,
    /// Maximum interruption (bp) between two seeds of the same motif.
    pub max_gap: usize,
    /// Minimum fraction of the span covered by seeds.
    pub min_purity: f64,
}

impl Default for ImperfectParams {
    fn default() -> Self {
        let seed_repeats = (1..=MAX_UNIT).map(|u| (u, Self::default_seed(u))).collect();
        Self {
            seed_repeats,
            max_gap: 5,
            min_purity: 0.85,
        }
    }
}

impl ImperfectParams {
    /// Default seed copies: `max(2, ceil(8 / unit))`.
    pub fn default_seed(unit: usize) -> usize {
        8usize.div_ceil(unit).max(2)
    }

    pub fn seed_for(&self, unit: usize) -> usize {
        self.seed_repeats
            .get(&unit)
            .copied()
            .unwrap_or_else(|| Self::default_seed(unit))
            .max(2)
    }

    pub fn validate(&self) -> Result<(), ParamError> {
        if !(0.0..=1.0).contains(&self.min_purity) {
            return err("imperfect min_purity must be within [0, 1]");
        }
        Ok(())
    }
}

/// Default compound distance in native mode: join only adjacent or
/// overlapping SSRs (misa.ini uses 100).
pub const NATIVE_MAX_INTERRUPTION: i64 = 0;
/// Default compound distance of misa.pl (`interruptions` in misa.ini).
pub const MISA_MAX_INTERRUPTION: i64 = 100;

/// Complete configuration of an analysis run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Config {
    pub definition: Definition,
    pub mode: Mode,
    /// Maximum distance (bp) between two SSRs forming a compound SSR;
    /// `None` disables compound formation. MISA uses 100; the native default
    /// is 0, so only SSRs that touch or overlap are reported together.
    pub max_interruption: Option<i64>,
    /// Imperfect SSR detection (native mode only).
    pub imperfect: Option<ImperfectParams>,
    /// Primer design (None = disabled).
    pub primers: Option<PrimerParams>,
    /// Length of flanking sequence to report on each side (0 = none).
    pub flank: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            definition: Definition::misa_default(),
            mode: Mode::Native,
            max_interruption: Some(NATIVE_MAX_INTERRUPTION),
            imperfect: None,
            primers: None,
            flank: 0,
        }
    }
}

impl Config {
    /// Configuration reproducing `misa.pl` with a parsed `misa.ini`.
    pub fn misa(ini: &MisaIni) -> Self {
        Self {
            definition: ini.definition.clone(),
            mode: Mode::Misa,
            max_interruption: Some(ini.max_interruption),
            imperfect: None,
            primers: None,
            flank: 0,
        }
    }

    pub fn validate(&self) -> Result<(), ParamError> {
        self.definition.validate()?;
        if let Some(d) = self.max_interruption {
            if d < 0 {
                return err("max interruption must be >= 0");
            }
        }
        if self.mode == Mode::Misa && self.imperfect.is_some() {
            return err("imperfect SSRs are not available in MISA-compatible mode");
        }
        if let Some(i) = &self.imperfect {
            i.validate()?;
        }
        if let Some(p) = &self.primers {
            p.validate()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_definition() {
        let d = Definition::parse("1-10 2-6 3-5 4-5 5-5 6-5").unwrap();
        assert_eq!(d, Definition::misa_default());
        assert_eq!(d.min_threshold(), 5);
        assert!(Definition::parse("1-10 2").is_err());
        assert!(Definition::parse("2-1").is_err());
        assert!(Definition::parse("0-5").is_err());
    }

    #[test]
    fn parse_misa_ini() {
        let ini = "definition(unit_size,min_repeats):                   1-10 2-6 3-5 4-5 5-5 6-5\n\
                   interruptions(max_difference_between_2_SSRs):        100\n\
                   GFF:\t\t\t\t\t\t\t    true\n";
        let m = MisaIni::parse(ini).unwrap();
        assert_eq!(m.definition, Definition::misa_default());
        assert_eq!(m.max_interruption, 100);
        assert!(m.gff);
    }

    #[test]
    fn misa_ini_without_interruptions_defaults_to_zero() {
        let m = MisaIni::parse("def 2-6 3-5\n").unwrap();
        assert_eq!(m.max_interruption, 0);
        assert!(!m.gff);
        assert_eq!(m.definition.get(2), Some(6));
    }
}
