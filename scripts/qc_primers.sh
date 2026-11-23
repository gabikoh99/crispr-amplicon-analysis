#!/usr/bin/env bash
# Primer QC: for each R1 FASTQ, count reads containing the RAB11A on-target forward primer
# and the forward primers of the four off-target amplicons (Nam et al. 2025, Suppl. Table 3).
# Usage: bash scripts/qc_primers.sh data/fastq results/primer_qc.csv
set -euo pipefail
dir="${1:-data/fastq}"
out="${2:-results/primer_qc.csv}"
mkdir -p "$(dirname "$out")"
echo "file,total_reads,on_target,offtarget1,offtarget2,offtarget3,offtarget4,on_target_fraction" > "$out"
for f in "$dir"/*_1.fastq.gz; do
  gunzip -c "$f" | awk -v n="$(basename "$f")" 'BEGIN{OFS=","}
    NR%4==2 { t++
      if (index($0,"ACCCAGTCCGGCAGTTGAA")) a++
      if (index($0,"TTCCCGGCTGCCTGAC")) b++
      if (index($0,"TCATGACCTCATGCCTGTCTC")) c++
      if (index($0,"CATGCTAGGGAGCAATGGTAGA")) d++
      if (index($0,"CACCCGAGCTGTTGTAGGTC")) e++ }
    END { printf "%s,%d,%d,%d,%d,%d,%d,%.3f\n", n, t, a, b, c, d, e, (t ? a/t : 0) }' >> "$out"
done
echo "Wrote $out"
