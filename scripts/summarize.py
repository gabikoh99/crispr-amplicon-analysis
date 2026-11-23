#!/usr/bin/env python3
"""Summarise CRISPRessoBatch output per condition and plot it.

Indel % = (Modified reads - reads with ONLY substitutions) / aligned reads,
so that substitutions (mostly sequencing/PCR errors) are not counted as editing.

Usage:
  python scripts/summarize.py \
      --quant results/CRISPRessoBatch_on_rab11a/CRISPRessoBatch_quantification_of_editing_frequency.txt \
      --samples data/metadata/samples.csv --outdir results
"""
import argparse, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN_ALIGNED = 1000      # below this a sample is excluded from means
MIN_ALIGNED_FRAC = 0.5  # below this a sample is flagged

# Ranges reported by Nam et al. 2025 (Suppl. Fig. S1): ~80% and 60-70% indels
PAPER_HIGH = {"Cas9_WT", "FECO", "Brex", "DrRecA_20AA"}
PAPER_LOW = {"RECA", "RAD51", "DrRecA"}
ORDER = ["negative_control", "Cas9_WT", "Cas9_WT_RNP", "FECO", "Brex", "DrRecA_20AA", "RECA", "RAD51", "DrRecA"]

ap = argparse.ArgumentParser()
ap.add_argument("--quant", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--outdir", default="results")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

q = pd.read_csv(a.quant, sep="\t")
s = pd.read_csv(a.samples)
s["Batch"] = s["sample_name"].str.replace(" ", "_")
d = q.merge(s[["Batch", "condition", "replicate", "delivery"]], on="Batch", how="left")
if d["condition"].isna().any():
    raise SystemExit("Some batch names are not in samples.csv: " + ", ".join(d.loc[d.condition.isna(), "Batch"]))

d["indel_pct"] = 100 * (d["Modified"] - d["Only Substitutions"]) / d["Reads_aligned"]
d["aligned_frac"] = d["Reads_aligned"] / d["Reads_in_input"]
d["flag"] = ""
d.loc[d["aligned_frac"] < MIN_ALIGNED_FRAC, "flag"] = "low alignment rate"
d.loc[d["Reads_aligned"] < MIN_ALIGNED, "flag"] = "excluded: too few aligned reads"

cols = ["Batch", "condition", "replicate", "delivery", "Reads_in_input", "Reads_aligned",
        "aligned_frac", "indel_pct", "flag"]
d[cols].to_csv(os.path.join(a.outdir, "indel_summary_per_sample.csv"), index=False)

ok = d[~d["flag"].str.startswith("excluded")]
m = ok.groupby("condition")["indel_pct"].agg(["mean", "min", "max", "count"]).reindex(ORDER)
m.to_csv(os.path.join(a.outdir, "indel_summary_per_condition.csv"))
print(d[cols].round(2).to_string(index=False))
print("\nPer condition (excluded samples left out):")
print(m.round(1).to_string())

fig, ax = plt.subplots(figsize=(9, 4.8))
ax.axhspan(75, 85, color="#2a9d8f", alpha=0.12, label="Paper: ~80%")
ax.axhspan(60, 70, color="#e9c46a", alpha=0.25, label="Paper: 60-70%")
for i, c in enumerate(ORDER):
    sub = d[d["condition"] == c]
    good = sub[~sub["flag"].str.startswith("excluded")]
    if len(good):
        ax.bar(i, good["indel_pct"].mean(), color="#9bb1c7", width=0.6, zorder=2)
    for _, r in sub.iterrows():
        if r["flag"].startswith("excluded"):
            ax.scatter(i, 0, marker="x", color="red", zorder=4)
            ax.annotate("excluded", (i, 2), ha="center", fontsize=7, color="red")
        else:
            ax.scatter(i + (0.12 if r["replicate"] == 2 else -0.12), r["indel_pct"],
                       color="#1d3557", zorder=3, s=22,
                       facecolors="none" if r["flag"] else "#1d3557")
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=35, ha="right")
ax.set_ylabel("Indel reads (% of aligned reads)")
ax.set_ylim(0, 100)
ax.set_title("RAB11A editing in K562: reanalysis of PRJNA1240319 with CRISPResso2")
ax.legend(loc="upper right", fontsize=8, frameon=False)
fig.text(0.01, 0.01, "Dots = replicates (hollow = low alignment rate). Bars = mean of valid replicates.",
         fontsize=7, color="gray")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(a.outdir, "indel_by_condition.png"), dpi=200)
print("Saved table and figure in", a.outdir)
