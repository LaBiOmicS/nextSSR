#!/usr/bin/env bash
# Compare `nextssr misa` with the original misa.pl v2.0 on one or more FASTA files.
#
#   cargo build --release
#   scripts/compare_misa.sh tests/data/tricky.fasta path/to/genome.fa
#
# misa.pl is downloaded from IPK Gatersleben (it is not redistributed here).
# Environment: NEXTSSR (binary, default target/release/nextssr), INI (default tests/data/misa.ini).
set -euo pipefail
nextssr=$(realpath "${NEXTSSR:-target/release/nextssr}")
ini=$(realpath "${INI:-tests/data/misa.ini}")
url=https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
curl -sSfL -o "$work/misa.zip" "$url"
python3 -I -m zipfile -e "$work/misa.zip" "$work/misa_src"
status=0
for fa in "$@"; do
  name=$(basename "$fa")
  for side in misa nextssr; do
    mkdir -p "$work/$side"
    cp "$fa" "$work/$side/$name"
    cp "$ini" "$work/$side/misa.ini"
  done
  t0=$(date +%s.%N)
  (cd "$work/misa" && perl "$work/misa_src/misa.pl" "$name" misa.ini)
  t1=$(date +%s.%N)
  (cd "$work/nextssr" && "$nextssr" misa "$name" misa.ini)
  t2=$(date +%s.%N)
  printf '%s: misa.pl %.2fs, nextssr %.2fs\n' "$name" "$(echo "$t1 - $t0" | bc)" "$(echo "$t2 - $t1" | bc)"
  for ext in misa statistics; do
    if cmp -s "$work/misa/$name.$ext" "$work/nextssr/$name.$ext"; then
      echo "  IDENTICAL  $name.$ext"
    else
      echo "  DIFFERENT  $name.$ext"
      diff "$work/misa/$name.$ext" "$work/nextssr/$name.$ext" | head -10 || true
      status=1
    fi
  done
  rm -rf "$work/misa" "$work/nextssr"
done
exit $status
