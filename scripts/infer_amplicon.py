#!/usr/bin/env python3
"""Infer the reference amplicon from the negative-control (NC) reads.

Merges R1/R2 with fastp (installed together with CRISPResso2), keeps merged reads
that start with the forward primer and end with the reverse-complemented reverse
primer, and reports the most frequent one as the reference amplicon.
Then checks that the guide (either strand) and its PAM are present.

Usage:
  python scripts/infer_amplicon.py --r1 data/fastq/SRR32811088_1.fastq.gz \
      --r2 data/fastq/SRR32811088_2.fastq.gz --out data/metadata/amplicon.txt
"""
import argparse, gzip, subprocess, sys, tempfile, os
from collections import Counter

FWD = "acccagtccggcagttgaa"          # Supplementary Table 3, RAB11A on-target forward
REV = "aaccgcacaaaaaggaaggg"         # Supplementary Table 3, RAB11A on-target reverse
GUIDE = "GGTAGTCGTACTCGTCGTCG"       # Supplementary Table 2, RAB11A

def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]

def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def merged_reads(path):
    with open_any(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1"); ap.add_argument("--r2")
    ap.add_argument("--merged", help="already merged FASTQ (skips fastp)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fwd-only", action="store_true",
                    help="keep merged reads that start with the forward primer, ignoring the reverse primer "
                         "(use when a library seems to have a different reverse primer)")
    ap.add_argument("--ref", help="reference amplicon file (e.g. data/metadata/amplicon.txt) to compare against")
    a = ap.parse_args()

    fwd, rev_rc = FWD.upper(), revcomp(REV).upper()
    tmp = None
    merged = a.merged
    if not merged:
        if not (a.r1 and a.r2):
            sys.exit("Give --r1 and --r2, or --merged")
        tmp = tempfile.mkdtemp()
        merged = os.path.join(tmp, "merged.fastq.gz")
        subprocess.run(["fastp", "-i", a.r1, "-I", a.r2, "--merge",
                        "--merged_out", merged,
                        "-j", os.path.join(tmp, "f.json"),
                        "-h", os.path.join(tmp, "f.html")], check=True)

    total = 0
    full = Counter()
    for seq in merged_reads(merged):
        total += 1
        if seq.startswith(fwd) and (a.fwd_only or seq.endswith(rev_rc)):
            full[seq] += 1
    if not full:
        sys.exit("No merged read starts with the forward primer and ends with the reverse primer. "
                 "Check primer orientation or inspect the reads manually.")

    kept = sum(full.values())
    label = "with forward primer" if a.fwd_only else "with both primers"
    print(f"merged reads: {total}; {label}: {kept} ({100*kept/total:.1f}%)")
    print("Top sequences (count, length):")
    for seq, n in full.most_common(5):
        print(f"  {n}\t{len(seq)}")
    amplicon, n = full.most_common(1)[0]
    print(f"\nChosen amplicon: {len(amplicon)} bp, {n} reads ({100*n/kept:.1f}% of reads {label})")

    print(f"Last 30 nt of chosen amplicon: ...{amplicon[-30:]}")
    if a.ref:
        ref = open(a.ref).read().strip().upper()
        shared = 0
        for x, y in zip(amplicon, ref):
            if x != y:
                break
            shared += 1
        print(f"Reference amplicon: {len(ref)} bp; identical to chosen one for the first {shared} nt")
        print(f"Last 30 nt of reference:       ...{ref[-30:]}")

    g = GUIDE.upper()
    pos_f, pos_r = amplicon.find(g), amplicon.find(revcomp(g))
    if pos_f >= 0:
        pam = amplicon[pos_f + 20: pos_f + 23]
        print(f"Guide found on + strand at {pos_f}; PAM = {pam} ({'OK, NGG' if pam[1:] == 'GG' else 'NOT NGG'})")
    elif pos_r >= 0:
        pam_rc = amplicon[pos_r - 3: pos_r]
        print(f"Guide found on - strand at {pos_r}; PAM (revcomp) = {revcomp(pam_rc)} "
              f"({'OK, NGG' if pam_rc[:2] == 'CC' else 'NOT NGG'})")
    else:
        print("WARNING: guide not found in the amplicon. Do not use this amplicon without checking.")

    with open(a.out, "w") as fh:
        fh.write(amplicon + "\n")
    print(f"Saved to {a.out}")

if __name__ == "__main__":
    main()
