"""Raw-call concordance, CLC vs non-CLC (own caller), for 4 isolates (JL1123-JL1126), before any
MayoMate filtering. Writes comparison_outputs/raw_calls_summary.csv, raw_calls_paired_tests.csv and
comparison_outputs/raw_calls_concordance.png."""
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

OUT = "comparison_outputs"
os.makedirs(OUT, exist_ok=True)
BASE = "data/sequencing_calls/CombinedSet"
KEY = ["Chromosome", "Region", "Reference", "Allele"]
ISOLATES = [1123, 1124, 1125, 1126]
ISO_COLORS = dict(zip(ISOLATES, ["tab:green", "tab:red", "tab:purple", "tab:brown"]))
C_CLC, C_OWN = "tab:blue", "tab:orange"


def load(kind, n, clc):
    pat = f"JL0{n}_*" if clc else f"JL{n}_S.csv"
    files = glob.glob(f"{BASE}/{kind}/{pat}")
    assert len(files) == 1, (kind, pat, files)
    return pd.read_csv(files[0], dtype={"Region": str})


rows, merged = [], {}
for kind in ["SNPs", "SNVs"]:
    for n in ISOLATES:
        a, b = load(kind, n, True), load(kind, n, False)
        m = a.merge(b, on=KEY, how="outer", suffixes=("_CLC", "_own"), indicator=True)
        both = m[m["_merge"] == "both"]
        row = {
            "file": kind, "isolate": f"JL{n}", "calls_CLC": len(a), "calls_own": len(b), "shared": len(both),
            "only_CLC": int((m["_merge"] == "left_only").sum()), "only_own": int((m["_merge"] == "right_only").sum()),
            "jaccard": len(both) / len(m),
            "pct_CLC_found_in_own": 100 * len(both) / len(a), "pct_own_found_in_CLC": 100 * len(both) / len(b),
            "zygosity_agreement_pct": 100 * (both["Zygosity_CLC"] == both["Zygosity_own"]).mean(),
            "freq_pearson_r": stats.pearsonr(both["Frequency_CLC"], both["Frequency_own"])[0],
            "freq_median_diff_own_minus_CLC": (both["Frequency_own"] - both["Frequency_CLC"]).median(),
            "cov_pearson_r": stats.pearsonr(both["Coverage_CLC"], both["Coverage_own"])[0],
            "cov_median_CLC": a["Coverage"].median(), "cov_median_own": b["Coverage"].median(),
            "freq_median_CLC": a["Frequency"].median(), "freq_median_own": b["Frequency"].median(),
        }
        for t in sorted(set(a["Type"]) | set(b["Type"])):
            row[f"type_{t}_CLC"] = int((a["Type"] == t).sum())
            row[f"type_{t}_own"] = int((b["Type"] == t).sum())
        rows.append(row)
        merged[(kind, n)] = m
        m.to_csv(f"{OUT}/raw_{kind}_merged_JL{n}.csv", index=False)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/raw_calls_summary.csv", index=False)

# paired tests across the 4 isolates (CLC vs own, same isolate). With n=4 the smallest possible
# two-sided Wilcoxon p is 0.125, so these show direction/consistency, not significance; the
# paired t-test is given alongside.
tests = []
for kind in ["SNPs", "SNVs"]:
    d = df[df.file == kind]
    for metric in ["calls", "cov_median", "freq_median"]:
        x, y = d[f"{metric}_CLC"].values, d[f"{metric}_own"].values
        tests.append(dict(
            file=kind, metric=metric, CLC_values=", ".join(f"{v:g}" for v in x), own_values=", ".join(f"{v:g}" for v in y),
            mean_CLC=x.mean(), mean_own=y.mean(), own_minus_CLC_all_same_sign=bool(np.all(y > x) or np.all(y < x)),
            wilcoxon_p=stats.wilcoxon(x, y).pvalue if np.any(x != y) else 1.0,
            paired_t_p=stats.ttest_rel(x, y).pvalue if np.any(x != y) else 1.0))
pd.DataFrame(tests).to_csv(f"{OUT}/raw_calls_paired_tests.csv", index=False)

# figure: overlap per isolate, allele frequency and coverage of shared calls (all isolates)
fig, axes = plt.subplots(2, 3, figsize=(17, 10))
for r, kind in enumerate(["SNPs", "SNVs"]):
    d = df[df.file == kind].reset_index(drop=True)
    ax = axes[r, 0]
    x = np.arange(len(d))
    ax.bar(x, d.only_CLC, color=C_CLC, label="only CLC")
    ax.bar(x, d.shared, bottom=d.only_CLC, color="gray", label="shared")
    ax.bar(x, d.only_own, bottom=d.only_CLC + d.shared, color=C_OWN, label="only own caller")
    for i, j in enumerate(d.jaccard):
        ax.text(i, (d.only_CLC + d.shared + d.only_own)[i], f"J={j:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, (d.only_CLC + d.shared + d.only_own).max() * 1.35)
    ax.set_xticks(x); ax.set_xticklabels(d.isolate); ax.set_ylabel("calls"); ax.legend(fontsize=8, loc="upper left", ncol=3)
    ax.set_title(f"{kind}: raw call overlap per isolate")
    for c, (col, lab) in enumerate([("Frequency", "allele frequency (%)"), ("Coverage", "coverage")], start=1):
        ax = axes[r, c]
        allv = []
        for n in ISOLATES:
            both = merged[(kind, n)]
            both = both[both["_merge"] == "both"]
            ax.scatter(both[f"{col}_CLC"], both[f"{col}_own"], s=3, alpha=0.25, color=ISO_COLORS[n], label=f"JL{n}")
            allv.append(np.r_[both[f"{col}_CLC"], both[f"{col}_own"]])
        lim = 100 if col == "Frequency" else np.percentile(np.concatenate(allv), 99.5)
        ax.plot([0, lim], [0, lim], "k--", lw=1); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel(f"CLC {lab}"); ax.set_ylabel(f"own caller {lab}")
        ax.set_title(f"{kind}: {lab}, shared calls"); ax.legend(markerscale=4, fontsize=8)
plt.suptitle("Raw calls: CLC vs own caller, 4 ung1∆ isolates (JL1123-JL1126)", fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUT}/raw_calls_concordance.png", dpi=200)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(df.T.to_string())
print(pd.DataFrame(tests).to_string())
print("saved", f"{OUT}/raw_calls_concordance.png")
