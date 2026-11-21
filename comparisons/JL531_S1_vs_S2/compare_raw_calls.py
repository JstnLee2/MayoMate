"""Raw-call concordance between JL531_S1 (CLC) and JL531_S2 (own script), before any MayoMate filtering.
Writes comparison_outputs/raw_calls_summary.csv and figures to comparison_outputs/."""
import glob, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

OUT = "comparison_outputs"
os.makedirs(OUT, exist_ok=True)
BASE = "data/sequencing_calls/CombinedSet"
KEY = ["Chromosome", "Region", "Reference", "Allele"]


def load(kind, tag):
    pat = "JL531_S105*" if tag == "S1" else "JL531_S2.csv"
    f = glob.glob(f"{BASE}/{kind}/{pat}")[0]
    df = pd.read_csv(f, dtype={"Region": str})
    return df


rows = []
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for r, kind in enumerate(["SNPs", "SNVs"]):
    s1, s2 = load(kind, "S1"), load(kind, "S2")
    m = s1.merge(s2, on=KEY, how="outer", suffixes=("_S1", "_S2"), indicator=True)
    both = m[m["_merge"] == "both"]
    n1, n2, nb = len(s1), len(s2), len(both)
    row = {
        "file": kind, "calls_S1": n1, "calls_S2": n2, "shared": nb,
        "only_S1": (m["_merge"] == "left_only").sum(), "only_S2": (m["_merge"] == "right_only").sum(),
        "jaccard": nb / len(m), "pct_S1_found_in_S2": 100 * nb / n1, "pct_S2_found_in_S1": 100 * nb / n2,
        "zygosity_agreement_pct": 100 * (both["Zygosity_S1"] == both["Zygosity_S2"]).mean(),
        "freq_pearson_r": stats.pearsonr(both["Frequency_S1"], both["Frequency_S2"])[0],
        "freq_median_diff_S2_minus_S1": (both["Frequency_S2"] - both["Frequency_S1"]).median(),
        "cov_pearson_r": stats.pearsonr(both["Coverage_S1"], both["Coverage_S2"])[0],
        "cov_median_S1": s1["Coverage"].median(), "cov_median_S2": s2["Coverage"].median(),
        "freq_median_S1": s1["Frequency"].median(), "freq_median_S2": s2["Frequency"].median(),
        "wilcoxon_paired_freq_p": stats.wilcoxon(both["Frequency_S1"], both["Frequency_S2"]).pvalue,
        "wilcoxon_paired_cov_p": stats.wilcoxon(both["Coverage_S1"], both["Coverage_S2"]).pvalue,
    }
    for t in sorted(set(s1["Type"]) | set(s2["Type"])):
        row[f"type_{t}_S1"] = (s1["Type"] == t).sum()
        row[f"type_{t}_S2"] = (s2["Type"] == t).sum()
    rows.append(row)

    ax = axes[r, 0]
    ax.bar(["only S1", "shared", "only S2"], [row["only_S1"], nb, row["only_S2"]], color=["tab:blue", "gray", "tab:orange"])
    ax.set_title(f"{kind}: call overlap (Jaccard {row['jaccard']:.2f})")
    ax = axes[r, 1]
    ax.scatter(both["Frequency_S1"], both["Frequency_S2"], s=3, alpha=0.3)
    ax.plot([0, 100], [0, 100], "k--", lw=1)
    ax.set_xlabel("S1 frequency (%)"); ax.set_ylabel("S2 frequency (%)")
    ax.set_title(f"{kind}: allele frequency, shared calls (r={row['freq_pearson_r']:.2f})")
    ax = axes[r, 2]
    lim = np.percentile(np.r_[both["Coverage_S1"], both["Coverage_S2"]], 99.5)
    ax.scatter(both["Coverage_S1"], both["Coverage_S2"], s=3, alpha=0.3)
    ax.plot([0, lim], [0, lim], "k--", lw=1); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("S1 coverage"); ax.set_ylabel("S2 coverage")
    ax.set_title(f"{kind}: coverage, shared calls (r={row['cov_pearson_r']:.2f})")
    m.to_csv(f"{OUT}/raw_{kind}_merged.csv", index=False)

plt.tight_layout()
plt.savefig(f"{OUT}/raw_calls_concordance.png", dpi=200)
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/raw_calls_summary.csv", index=False)
pd.set_option("display.width", 250)
print(df.T)
