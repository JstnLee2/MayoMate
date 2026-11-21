"""Raw-call concordance, CLC vs own caller, for every paired isolate (pairs.json), before any
MayoMate filtering. Per isolate, then summarized per genotype and overall.

Writes to comparison_outputs/:
  raw_calls_per_isolate.csv        one row per isolate and file type (SNPs / SNVs)
  raw_calls_by_genotype.csv        medians per genotype
  raw_calls_tests.csv              paired Wilcoxon (own vs CLC) overall and per genotype, and
                                   Kruskal-Wallis: does the own-minus-CLC difference differ between genotypes?
  raw_calls_by_genotype.png
"""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from mayo.settings import config as cfg
from mayo.settings.genotype_colors import genotype_color

OUT = "comparison_outputs"
os.makedirs(OUT, exist_ok=True)
KEY = ["Chromosome", "Region", "Reference", "Allele"]
PATHS = {"SNPs": cfg["parental_SNPs_path"], "SNVs": cfg["SNVs_path"]}
PAIRS = json.load(open("pairs.json", encoding="utf-8"))["pairs"]
GENO_ORDER = list(dict.fromkeys(p["genotype"] for p in sorted(PAIRS.values(), key=lambda p: p["genotype"])))
FILES = {k: os.listdir(v) for k, v in PATHS.items()}


def load(kind, key):
    f = [x for x in FILES[kind] if x.startswith(key) and x.endswith(".csv")]
    assert len(f) == 1, (kind, key, f)
    return pd.read_csv(os.path.join(PATHS[kind], f[0]), dtype={"Region": str},
                       usecols=lambda c: c in KEY + ["Type", "Zygosity", "Coverage", "Frequency"])


rows = []
for iso, p in PAIRS.items():
    for kind in ["SNPs", "SNVs"]:
        a, b = load(kind, p["CLC"]), load(kind, p["own"])
        m = a.merge(b, on=KEY, how="outer", suffixes=("_CLC", "_own"), indicator=True)
        both = m[m["_merge"] == "both"]
        r = dict(isolate=iso, genotype=p["genotype"], file=kind, calls_CLC=len(a), calls_own=len(b), shared=len(both),
                 only_CLC=int((m["_merge"] == "left_only").sum()), only_own=int((m["_merge"] == "right_only").sum()),
                 jaccard=len(both) / max(1, len(m)),
                 zygosity_agreement_pct=100 * (both.Zygosity_CLC == both.Zygosity_own).mean() if len(both) else np.nan,
                 freq_pearson_r=stats.pearsonr(both.Frequency_CLC, both.Frequency_own)[0] if len(both) > 2 else np.nan,
                 cov_pearson_r=stats.pearsonr(both.Coverage_CLC, both.Coverage_own)[0] if len(both) > 2 else np.nan,
                 cov_median_CLC=a.Coverage.median(), cov_median_own=b.Coverage.median(),
                 freq_median_CLC=a.Frequency.median(), freq_median_own=b.Frequency.median())
        r["cov_ratio_own_over_CLC"] = r["cov_median_own"] / r["cov_median_CLC"] if r["cov_median_CLC"] else np.nan
        for t in ["SNV", "Insertion", "Deletion", "MNV", "Replacement"]:
            r[f"{t}_CLC"], r[f"{t}_own"] = int((a.Type == t).sum()), int((b.Type == t).sum())
        rows.append(r)
    print(f"{iso} done", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/raw_calls_per_isolate.csv", index=False)
num = [c for c in df.columns if c not in ("isolate", "genotype", "file")]
df.groupby(["file", "genotype"])[num].median().reset_index().to_csv(f"{OUT}/raw_calls_by_genotype.csv", index=False)

tests = []
for kind in ["SNPs", "SNVs"]:
    d = df[df.file == kind]
    for metric in ["calls", "cov_median", "freq_median", "SNV", "Insertion", "Deletion"]:
        for g, dd in [("ALL", d)] + [(g, d[d.genotype == g]) for g in GENO_ORDER]:
            x, y = dd[f"{metric}_CLC"].values, dd[f"{metric}_own"].values
            diff = y - x
            tests.append(dict(file=kind, metric=metric, genotype=g, n=len(dd), median_CLC=np.median(x), median_own=np.median(y),
                              median_own_minus_CLC=np.median(diff), pct_isolates_own_higher=100 * (diff > 0).mean(),
                              pct_isolates_own_lower=100 * (diff < 0).mean(),
                              wilcoxon_p=stats.wilcoxon(x, y).pvalue if (diff != 0).any() else 1.0))
        diffs = [d[d.genotype == g][f"{metric}_own"].values - d[d.genotype == g][f"{metric}_CLC"].values for g in GENO_ORDER]
        tests.append(dict(file=kind, metric=metric, genotype="between genotypes (Kruskal-Wallis on own-CLC)", n=len(d),
                          wilcoxon_p=np.nan, kruskal_p=stats.kruskal(*diffs).pvalue))
    for metric in ["jaccard", "cov_ratio_own_over_CLC"]:
        groups = [d[d.genotype == g][metric].dropna().values for g in GENO_ORDER]
        tests.append(dict(file=kind, metric=metric, genotype="between genotypes (Kruskal-Wallis)", n=len(d), kruskal_p=stats.kruskal(*groups).pvalue))
tests = pd.DataFrame(tests)
tests.to_csv(f"{OUT}/raw_calls_tests.csv", index=False)

fig, axes = plt.subplots(2, 3, figsize=(20, 11))
panels = [("SNPs", "jaccard", "SNPs: CLC/own overlap (Jaccard)"), ("SNPs", "cov_ratio_own_over_CLC", "SNPs: coverage, own / CLC"),
          ("SNPs", "zygosity_agreement_pct", "SNPs: zygosity agreement (%)"), ("SNVs", "jaccard", "SNVs: CLC/own overlap (Jaccard)"),
          ("SNVs", "calls_ratio", "SNVs: calls, own / CLC"), ("SNVs", "cov_ratio_own_over_CLC", "SNVs: coverage, own / CLC")]
df["calls_ratio"] = df.calls_own / df.calls_CLC
for ax, (kind, col, title) in zip(axes.flat, panels):
    d = df[df.file == kind]
    data = [d[d.genotype == g][col].dropna() for g in GENO_ORDER]
    bp = ax.boxplot(data, patch_artist=True, showfliers=False)
    for patch, g in zip(bp["boxes"], GENO_ORDER):
        patch.set_facecolor(genotype_color(g)); patch.set_alpha(0.6)
    for i, v in enumerate(data, start=1):
        ax.scatter(np.random.normal(i, 0.06, len(v)), v, s=8, color="k", alpha=0.5, zorder=3)
    if "ratio" in col:
        ax.axhline(1, color="gray", ls="--", lw=1)
    ax.set_xticks(range(1, len(GENO_ORDER) + 1))
    ax.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(GENO_ORDER, data)], rotation=45, ha="right", fontsize=8)
    ax.set_title(title)
plt.suptitle(f"Raw calls: CLC vs own caller, {len(PAIRS)} paired isolates, by genotype", fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUT}/raw_calls_by_genotype.png", dpi=200)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(tests[tests.genotype.isin(["ALL"]) | tests.genotype.str.startswith("between")].to_string())
