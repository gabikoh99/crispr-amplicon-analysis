#!/usr/bin/env python3
"""Build the CRISPRessoBatch settings file (tab-separated) from samples.csv.

Usage:
  python scripts/make_batch.py --samples data/metadata/samples.csv \
      --amplicon data/metadata/amplicon.txt --out data/metadata/batch.tsv
"""
import argparse, csv

GUIDE = "GGTAGTCGTACTCGTCGTCG"   # RAB11A, Supplementary Table 2

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--amplicon", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

amplicon = open(a.amplicon).read().strip()
with open(a.samples) as fh, open(a.out, "w") as out:
    out.write("n\tr1\tr2\ta\tg\n")
    for row in csv.DictReader(fh):
        name = row["sample_name"].replace(" ", "_")
        out.write(f"{name}\t{row['fastq_r1']}\t{row['fastq_r2']}\t{amplicon}\t{GUIDE}\n")
print(f"Wrote {a.out}")
