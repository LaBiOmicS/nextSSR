//! Nucleotide and motif utilities.

/// Case-insensitive nucleotide code: A=0, C=1, G=2, T=3, anything else=4.
pub static CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

/// `true` for A, C, G, T (any case).
#[inline]
pub fn is_acgt(b: u8) -> bool {
    CODE[b as usize] < 4
}

/// Complement of an (uppercase or lowercase) base; non-ACGT bases map to `N`.
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

/// Reverse complement.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Reverse complement of a string.
pub fn revcomp_str(s: &str) -> String {
    String::from_utf8(revcomp(s.as_bytes())).expect("ASCII")
}

/// `true` if the motif is not a repetition of a shorter motif
/// (rejects e.g. `TT` or `ACAC`, as `misa.pl` does). Case-insensitive.
pub fn is_primitive(motif: &[u8]) -> bool {
    let k = motif.len();
    'outer: for d in 1..k {
        if k % d != 0 {
            continue;
        }
        for i in d..k {
            if CODE[motif[i] as usize] != CODE[motif[i % d] as usize] {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// Lexicographically smallest rotation.
pub fn min_rotation(m: &str) -> String {
    let n = m.len();
    if n == 0 {
        return String::new();
    }
    let doubled = format!("{m}{m}");
    (0..n)
        .map(|i| &doubled[i..i + n])
        .min()
        .expect("non-empty")
        .to_string()
}

/// MISA motif class, combining rotations and the reverse complement,
/// e.g. `AG` -> `AG/CT`, `AT` -> `AT/AT`, `C` -> `C/G`.
pub fn misa_class(motif: &str) -> String {
    let a = min_rotation(motif);
    let r = min_rotation(&revcomp_str(motif));
    if a < r {
        format!("{a}/{r}")
    } else {
        format!("{r}/{a}")
    }
}

/// Strand-independent canonical motif: smallest rotation of the motif or of its
/// reverse complement (e.g. `GA`, `TC`, `CT` -> `AG`).
pub fn canonical(motif: &str) -> String {
    let a = min_rotation(motif);
    let r = min_rotation(&revcomp_str(motif));
    a.min(r)
}

/// Standard name of a motif class by unit size.
pub fn unit_name(unit: usize) -> &'static str {
    match unit {
        1 => "mononucleotide",
        2 => "dinucleotide",
        3 => "trinucleotide",
        4 => "tetranucleotide",
        5 => "pentanucleotide",
        6 => "hexanucleotide",
        7 => "heptanucleotide",
        8 => "octanucleotide",
        9 => "nonanucleotide",
        10 => "decanucleotide",
        _ => "oligonucleotide",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primitive() {
        assert!(is_primitive(b"A"));
        assert!(is_primitive(b"AG"));
        assert!(!is_primitive(b"AA"));
        assert!(!is_primitive(b"ACAC"));
        assert!(!is_primitive(b"AAGAAG"));
        assert!(is_primitive(b"AAGAAC"));
        assert!(!is_primitive(b"acAC"));
    }

    #[test]
    fn classes() {
        assert_eq!(misa_class("AG"), "AG/CT");
        assert_eq!(misa_class("GA"), "AG/CT");
        assert_eq!(misa_class("TC"), "AG/CT");
        assert_eq!(misa_class("AT"), "AT/AT");
        assert_eq!(misa_class("T"), "A/T");
        assert_eq!(misa_class("AAG"), "AAG/CTT");
        assert_eq!(canonical("CTT"), "AAG");
        assert_eq!(revcomp_str("AACG"), "CGTT");
    }
}
