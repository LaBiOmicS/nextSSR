//! End-to-end tests of the `nextssr` binary.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_nextssr"))
}

fn run_ok(cmd: &mut Command) {
    let out = cmd.output().expect("run nextssr");
    assert!(
        out.status.success(),
        "nextssr failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// `nextssr misa` must reproduce misa.pl v2.0 byte for byte. The expected files
/// were produced by `scripts/misa_port.py`, a line-by-line port of misa.pl;
/// `scripts/compare_misa.sh` checks nextSSR against the original Perl script.
#[test]
fn misa_mode_is_byte_identical() {
    for fasta in ["tricky.fasta", "NC_000932.1_Arabidopsis_cp.fasta"] {
        let dir = tempfile::tempdir().unwrap();
        fs::copy(data(fasta), dir.path().join(fasta)).unwrap();
        fs::copy(data("misa.ini"), dir.path().join("misa.ini")).unwrap();
        run_ok(
            bin()
                .current_dir(dir.path())
                .args(["misa", fasta, "misa.ini"]),
        );
        for ext in ["misa", "statistics"] {
            let got = fs::read(dir.path().join(format!("{fasta}.{ext}"))).unwrap();
            let exp = fs::read(data(&format!("{fasta}.{ext}.expected"))).unwrap();
            assert!(
                got == exp,
                "{fasta}.{ext} differs from the misa.pl reference"
            );
        }
    }
}

fn run_native(dir: &Path, input: &str, extra: &[&str]) -> String {
    let prefix = dir.join("out").to_string_lossy().into_owned();
    let mut cmd = bin();
    cmd.args([
        "run",
        input,
        "-o",
        &prefix,
        "-q",
        "-f",
        "tsv,gff3,bed,misa,statistics,json,p3in",
    ]);
    cmd.args(extra);
    run_ok(&mut cmd);
    prefix
}

#[test]
fn native_outputs_are_well_formed() {
    let dir = tempfile::tempdir().unwrap();
    let input = data("tricky.fasta");
    let prefix = run_native(
        dir.path(),
        input.to_str().unwrap(),
        &["--imperfect", "--primers", "--flank", "30"],
    );

    // TSV: constant column count
    let tsv = fs::read_to_string(format!("{prefix}.tsv")).unwrap();
    let mut lines = tsv.lines();
    let ncol = lines.next().unwrap().split('\t').count();
    let mut n = 0;
    for l in lines {
        assert_eq!(l.split('\t').count(), ncol, "bad TSV row: {l}");
        n += 1;
    }
    assert!(n > 100);

    // GFF3: 9 columns, valid coordinates, no unescaped reserved characters in attributes
    let gff = fs::read_to_string(format!("{prefix}.gff3")).unwrap();
    assert!(gff.starts_with("##gff-version 3\n"));
    for l in gff.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f.len(), 9, "{l}");
        let (s, e): (usize, usize) = (f[3].parse().unwrap(), f[4].parse().unwrap());
        assert!(s >= 1 && s <= e, "{l}");
        for attr in f[8].split(';') {
            let kv: Vec<&str> = attr.splitn(2, '=').collect();
            assert_eq!(kv.len(), 2, "{l}");
            assert!(!kv[1].contains('='), "unescaped '=' in {l}");
        }
    }

    // BED: 0-based half-open
    let bed = fs::read_to_string(format!("{prefix}.bed")).unwrap();
    for l in bed.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f.len(), 6);
        assert!(f[1].parse::<usize>().unwrap() < f[2].parse::<usize>().unwrap());
    }

    // JSON summary
    let json: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(format!("{prefix}.summary.json")).unwrap())
            .unwrap();
    assert_eq!(json["stats"]["sequences"], 60);
    assert_eq!(json["input"]["sha256"].as_str().unwrap().len(), 64);
    assert!(json["stats"]["imperfect"].as_u64().unwrap() > 0);
    assert!(json["stats"]["primers_designed"].as_u64().unwrap() > 0);

    // Primer3 input records end with '='
    let p3 = fs::read_to_string(format!("{prefix}.p3in")).unwrap();
    assert_eq!(p3.matches("\n=\n").count(), n);
}

#[test]
fn deterministic_across_threads_gzip_and_stdin() {
    let dir = tempfile::tempdir().unwrap();
    let input = data("tricky.fasta");
    let raw = fs::read(&input).unwrap();

    let d1 = dir.path().join("t1");
    let d4 = dir.path().join("t4");
    let dz = dir.path().join("gz");
    let ds = dir.path().join("stdin");
    for d in [&d1, &d4, &dz, &ds] {
        fs::create_dir(d).unwrap();
    }
    let a = run_native(
        &d1,
        input.to_str().unwrap(),
        &["-t", "1", "--imperfect", "--primers"],
    );
    let b = run_native(
        &d4,
        input.to_str().unwrap(),
        &["-t", "4", "--imperfect", "--primers"],
    );

    let gz_path = dir.path().join("tricky.fasta.gz");
    {
        let mut e = flate2::write::GzEncoder::new(
            fs::File::create(&gz_path).unwrap(),
            flate2::Compression::default(),
        );
        e.write_all(&raw).unwrap();
        e.finish().unwrap();
    }
    let c = run_native(
        &dz,
        gz_path.to_str().unwrap(),
        &["--imperfect", "--primers"],
    );

    let prefix = ds.join("out").to_string_lossy().into_owned();
    let mut child = bin()
        .args([
            "run",
            "-",
            "-o",
            &prefix,
            "-q",
            "-f",
            "tsv,gff3",
            "--imperfect",
            "--primers",
        ])
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&raw).unwrap();
    assert!(child.wait().unwrap().success());

    for ext in ["tsv", "gff3", "bed", "misa"] {
        let x = fs::read(format!("{a}.{ext}")).unwrap();
        assert_eq!(
            x,
            fs::read(format!("{b}.{ext}")).unwrap(),
            "{ext}: 1 vs 4 threads"
        );
        assert_eq!(
            x,
            fs::read(format!("{c}.{ext}")).unwrap(),
            "{ext}: plain vs gzip"
        );
    }
    for ext in ["tsv", "gff3"] {
        assert_eq!(
            fs::read(format!("{a}.{ext}")).unwrap(),
            fs::read(format!("{prefix}.{ext}")).unwrap(),
            "{ext}: file vs stdin"
        );
    }
}

#[test]
fn rejects_invalid_parameters() {
    let input = data("tricky.fasta");
    for args in [
        vec!["run", input.to_str().unwrap(), "-d", "2-1"],
        vec!["run", input.to_str().unwrap(), "-f", "xml"],
        vec![
            "run",
            input.to_str().unwrap(),
            "--mode",
            "misa",
            "--imperfect",
        ],
        vec!["run", "does-not-exist.fa"],
    ] {
        let out = bin().args(&args).output().unwrap();
        assert!(!out.status.success(), "{args:?} should fail");
        assert!(String::from_utf8_lossy(&out.stderr).contains("error"));
    }
    // a missing input must not leave empty output files behind
    let dir = tempfile::tempdir().unwrap();
    let out = bin()
        .current_dir(dir.path())
        .args(["run", "does-not-exist.fa"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn tm_subcommand() {
    let out = bin()
        .args(["tm", "AGCTGCATGCAGTCAGTCAG", "ACGN"])
        .output()
        .unwrap();
    let s = String::from_utf8(out.stdout).unwrap();
    let mut lines = s.lines();
    let t: f64 = lines
        .next()
        .unwrap()
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!((55.0..65.0).contains(&t));
    assert!(lines.next().unwrap().ends_with("NA"));
}

#[test]
fn primers_with_epcr_self_check() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("cp");
    let p = prefix.to_str().unwrap();
    run_ok(bin().args([
        "run",
        data("NC_000932.1_Arabidopsis_cp.fasta").to_str().unwrap(),
        "-o",
        p,
        "--primers",
        "--num-pairs",
        "3",
        "--epcr",
        "-q",
    ]));
    let table = fs::read_to_string(format!("{p}.primers.tsv")).unwrap();
    let mut lines = table.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert!(rows.len() > 100);
    // every designed pair amplifies at least its own product, which is found
    for r in &rows {
        assert_ne!(r[col("epcr_status")], "no_product");
        assert_ne!(r[col("epcr_amplicons")], "0");
    }
    // at most one recommended pair per locus, and it is specific
    let mut seen = std::collections::HashSet::new();
    for r in rows.iter().filter(|r| r[col("recommended")] == "yes") {
        assert_eq!(r[col("epcr_status")], "specific");
        assert!(seen.insert((r[col("seq_id")], r[col("ssr_nr")])));
    }
    // the chloroplast inverted repeats make some pairs amplify twice
    assert!(rows.iter().any(|r| r[col("epcr_status")] == "nonspecific"));
    let amplicons = fs::read_to_string(format!("{p}.epcr.tsv")).unwrap();
    assert!(amplicons.lines().any(|l| l.ends_with("\tyes")));
    let summary: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(format!("{p}.summary.json")).unwrap()).unwrap();
    assert_eq!(
        summary["epcr"]["pairs"].as_u64().unwrap(),
        rows.len() as u64
    );
    assert_eq!(
        summary["epcr"]["loci_with_specific_pair"].as_u64().unwrap(),
        seen.len() as u64
    );

    // standalone e-PCR of the same table reproduces the counts
    run_ok(bin().args([
        "epcr",
        "-p",
        &format!("{p}.primers.tsv"),
        data("NC_000932.1_Arabidopsis_cp.fasta").to_str().unwrap(),
        "-o",
        p,
        "-q",
    ]));
    let pairs = fs::read_to_string(format!("{p}.epcr_pairs.tsv")).unwrap();
    let counts: Vec<&str> = pairs
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(4).unwrap())
        .collect();
    let expected: Vec<&str> = rows.iter().map(|r| r[col("epcr_amplicons")]).collect();
    assert_eq!(counts, expected);
}

#[test]
fn epcr_single_pair_to_stdout() {
    // a pair from the chloroplast table, also given as a gzip target
    let fa = fs::read(data("NC_000932.1_Arabidopsis_cp.fasta")).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let gz = dir.path().join("cp.fa.gz");
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    e.write_all(&fa).unwrap();
    fs::write(&gz, e.finish().unwrap()).unwrap();
    let out = bin()
        .args([
            "epcr",
            "-F",
            "ATCCGCCCCTACGCTACTAT",
            "-R",
            "AAGGTGGAATTTGCTACCTTTTT",
            "--sequence",
            gz.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = s.lines().collect();
    assert_eq!(lines.len(), 2, "{s}");
    let f: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(&f[..5], ["pair1", "NC_000932.1", "71", "171", "101"]);
    assert_eq!(f.last().unwrap().len(), 101);
    assert!(f.last().unwrap().starts_with("ATCCGCCCCTACGCTACTAT"));
}

#[test]
fn epcr_rejects_bad_usage() {
    let input = data("tricky.fasta");
    let input = input.to_str().unwrap();
    for args in [
        vec!["run", "-", "--primers", "--epcr"],
        vec!["run", input, "-f", "primers"],
        vec!["run", input, "--epcr"],
        vec!["epcr", input],
        vec![
            "epcr",
            input,
            "-F",
            "ACGTACGTACGTACGT",
            "-R",
            "ACGTXCGTACGTACGT",
        ],
        vec![
            "epcr",
            input,
            "-F",
            "ACGTACGTACGTACGT",
            "-R",
            "ACGTACGTACGTACGT",
            "--mismatches",
            "1",
            "--three-prime-mismatches",
            "2",
        ],
    ] {
        let out = bin().args(&args).stdin(Stdio::null()).output().unwrap();
        assert!(!out.status.success(), "{args:?} should fail");
    }
}

#[test]
fn default_outputs_include_html_report() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t");
    let p = prefix.to_str().unwrap();
    run_ok(bin().args(["run", data("tricky.fasta").to_str().unwrap(), "-o", p, "-q"]));
    let html = fs::read_to_string(format!("{p}.report.html")).unwrap();
    assert!(html.starts_with("<!doctype html>"));
    assert!(html.contains("SSRs per Mb") && html.contains("Methods text"));
    assert!(html.contains("Input SHA-256"));
    // the MISA statistics file is opt-in in `run`
    assert!(!Path::new(&format!("{p}.statistics")).exists());
    run_ok(bin().args([
        "run",
        data("tricky.fasta").to_str().unwrap(),
        "-o",
        p,
        "-q",
        "-f",
        "statistics",
    ]));
    assert!(Path::new(&format!("{p}.statistics")).exists());
}
