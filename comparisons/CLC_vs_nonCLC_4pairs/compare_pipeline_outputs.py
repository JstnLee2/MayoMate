"""CLC vs own caller (non-CLC) comparison of MayoMate outputs for 4 ung1∆ isolates (JL1123-JL1126).
Run after main_playground_clean.ipynb, cluster_analysis.py and cluster_coordination.py in this folder.

For every isolate, the CLC sample (JL0####_) is compared with the own-caller sample (JL####_S) of
the same isolate. Across the 4 isolates, each metric gets a paired test (Wilcoxon signed-rank and
paired t-test). With 4 pairs the smallest possible two-sided Wilcoxon p is 0.125, so "all 4 isolates
move the same way" is reported as well. Counts are also pooled over isolates (binomial test).

Writes to comparison_outputs/:
  pipeline_comparison_per_isolate.csv   every metric, per isolate, CLC vs own
  pipeline_comparison_paired_tests.csv  paired tests across the 4 isolates
  mutations_CLC_vs_own_JL####.csv       final mutations, merged per isolate
  clusters_by_pval_per_isolate.csv      JT clusters per p-value threshold
  own_only_clusters.csv                 own-caller clusters with no overlapping CLC cluster (artifact check)
  pipeline_comparison_overview.png, pipeline_comparison_clusters.png
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from mayo.base_functions import load_payload_from_pickle

OUT = "comparison_outputs"
os.makedirs(OUT, exist_ok=True)
ISOLATES = [1123, 1124, 1125, 1126]
PAIRS = {n: (f"JL0{n}_", f"JL{n}_S") for n in ISOLATES}          # isolate -> (CLC sample, own sample)
G_CLC, G_OWN = "ung1∆CLC", "ung1∆"
C_CLC, C_OWN = "tab:blue", "tab:orange"
PVALS = ["0.01", "0.001", "0.0001", "1e-05", "1e-06"]
KEY = ["Chromosome", "Region", "Reference", "Allele"]
rows = []  # per isolate: section, metric, isolate, CLC, own


def add(section, metric, n, v_clc, v_own, note=""):
    rows.append(dict(section=section, metric=metric, isolate=f"JL{n}", CLC=v_clc, own=v_own, note=note))


def count_test(k1, k2):
    """Exact test that two counts come from equal Poisson rates (same genome, same exposure)."""
    return stats.binomtest(int(k1), int(k1 + k2), 0.5).pvalue if (k1 + k2) else np.nan


# ---------- 1. Parental SNPs / haplotype assignment ----------
par = {s.name: s for s in load_payload_from_pickle("outputs/partial_pickles/samples_parental_SNPs_1snp.pkl")}
for n, (a, b) in PAIRS.items():
    pa, pb = par[a].output_df[["Chromosome", "Region", "Parent"]], par[b].output_df[["Chromosome", "Region", "Parent"]]
    pm = pa.merge(pb, on=["Chromosome", "Region"], suffixes=("_CLC", "_own"))
    add("parental SNPs", "parental SNPs genotyped", n, len(pa), len(pb))
    add("parental SNPs", "same parent assigned at shared positions (%)", n, 100 * (pm.Parent_CLC == pm.Parent_own).mean(), np.nan,
        f"{len(pm)} shared positions, {(pm.Parent_CLC != pm.Parent_own).sum()} discordant")
    add("parental SNPs", "APSC (avg per-SNP coverage)", n, par[a].apsc, par[b].apsc)
    add("parental SNPs", "% parental (largest parent)", n, par[a].percent, par[b].percent)
    add("parental SNPs", "aneuploid chromosomes", n, ",".join(par[a].aneuploid_chromosomes) or "none", ",".join(par[b].aneuploid_chromosomes) or "none")
    add("parental SNPs", "heterozygous chromosomes", n, ",".join(par[a].heterozygous_chromosomes) or "none", ",".join(par[b].heterozygous_chromosomes) or "none")

# ---------- 2. Final filtered mutations ----------
mut = pd.read_csv("outputs/all_mutations.tsv", sep="\t")
mut["sub"] = mut.Reference + ">" + mut.Allele
a3a = lambda d: d["sub"].isin(["C>T", "G>A"])
own_only_cov, shared_cov = [], []
for n, (a, b) in PAIRS.items():
    m1, m2 = mut[mut.sample_name == a], mut[mut.sample_name == b]
    allm = m1.merge(m2, on=KEY, how="outer", suffixes=("_CLC", "_own"), indicator=True)
    allm.to_csv(f"{OUT}/mutations_CLC_vs_own_JL{n}.csv", index=False)
    sh = (allm._merge == "both").sum()
    add("mutations (after all filters)", "mutations", n, len(m1), len(m2))
    add("mutations (after all filters)", "shared", n, sh, sh)
    add("mutations (after all filters)", "only in this caller", n, (allm._merge == "left_only").sum(), (allm._merge == "right_only").sum())
    add("mutations (after all filters)", "% of own mutations found by the other caller", n, 100 * sh / max(1, len(m1)), 100 * sh / max(1, len(m2)))
    add("mutations (after all filters)", "Jaccard", n, sh / max(1, len(allm)), sh / max(1, len(allm)))
    add("mutations (after all filters)", "% C>T or G>A (A3A-type)", n, 100 * a3a(m1).mean(), 100 * a3a(m2).mean())
    add("mutations (after all filters)", "median coverage", n, m1.Coverage.median(), m2.Coverage.median())
    add("mutations (after all filters)", "median allele frequency", n, m1.Frequency.median(), m2.Frequency.median())
    own_only_cov.append(allm.loc[allm._merge == "right_only", "Coverage_own"])
    shared_cov.append(allm.loc[allm._merge == "both", "Coverage_own"])
spec = pd.crosstab(mut.genotype, mut["sub"])
spec.to_csv(f"{OUT}/mutation_spectrum_pooled.csv")

# ---------- 3. JT clusters at each significance threshold ----------
def overlapping(a, b):
    """Rows of a that overlap any cluster of b on the same chromosome."""
    return np.array([((b.Chromosome == r.Chromosome) & (b.Start <= r.End) & (b.End >= r.Start)).any() for r in a.itertuples()], dtype=bool)


clust_rows, own_only_clusters = [], []
clusters = {}
for pv in PVALS:
    c = pd.read_csv(f"data/new_clusters/mutation_clusters/mutation_clusters_pval_{pv}.tsv", sep="\t")
    c = c[c.Mutations >= 3]
    clusters[pv] = c
    for n, (a, b) in PAIRS.items():
        c1, c2 = c[c["sample"] == a], c[c["sample"] == b]
        o1, o2 = overlapping(c1, c2), overlapping(c2, c1)
        clust_rows.append(dict(
            cluster_pval=pv, isolate=f"JL{n}", clusters_CLC=len(c1), clusters_own=len(c2),
            CLC_clusters_overlapping_own=int(o1.sum()), own_clusters_overlapping_CLC=int(o2.sum()),
            median_length_CLC=c1.Length.median(), median_length_own=c2.Length.median(),
            median_mutations_CLC=c1.Mutations.median(), median_mutations_own=c2.Mutations.median(),
            clustered_mutations_CLC=c1.Mutations.sum(), clustered_mutations_own=c2.Mutations.sum(),
            ssDNA_kb_CLC=c1.Length.sum() / 1000, ssDNA_kb_own=c2.Length.sum() / 1000))
        if pv == "0.0001":
            add("JT clusters (p=0.0001)", "clusters (>=3 mutations)", n, len(c1), len(c2))
            add("JT clusters (p=0.0001)", "clusters overlapping a cluster from the other caller", n, f"{o1.sum()}/{len(c1)}", f"{o2.sum()}/{len(c2)}")
            add("JT clusters (p=0.0001)", "clustered mutations", n, c1.Mutations.sum(), c2.Mutations.sum())
            add("JT clusters (p=0.0001)", "median cluster length (bp)", n, c1.Length.median(), c2.Length.median())
            add("JT clusters (p=0.0001)", "total clustered ssDNA (kb)", n, c1.Length.sum() / 1000, c2.Length.sum() / 1000)
            own_only_clusters.append(c2[~o2].assign(isolate=f"JL{n}"))
dfc = pd.DataFrame(clust_rows)
dfc.to_csv(f"{OUT}/clusters_by_pval_per_isolate.csv", index=False)

# own-caller clusters with no CLC counterpart: depth of their mutations (JL531 showed low-depth artifact patches)
ooc = pd.concat(own_only_clusters, ignore_index=True) if own_only_clusters else pd.DataFrame()
if len(ooc):
    covs = []
    for r in ooc.itertuples():
        b = PAIRS[int(r.isolate[2:])][1]
        mm = mut[(mut.sample_name == b) & (mut.Chromosome == r.Chromosome) & (mut.Region.between(r.Start, r.End))]
        covs.append(mm.Coverage.median())
    ooc["median_coverage_of_mutations"] = covs
ooc.to_csv(f"{OUT}/own_only_clusters.csv", index=False)

# ---------- 4. Switches and recombination events ----------
sw = pd.read_csv("outputs/all_switches_1snp.tsv", sep="\t")
ev = pd.read_csv("outputs/all_recombination_events_haploid_chromosomes_1snp.csv", sep="\t", encoding="utf-16")


def near(a, b, d=1000):
    return sum(((b.Chromosome == x.Chromosome) & ((b.Switch_Center - x.Switch_Center).abs() <= d)).any() for x in a.itertuples())


def near_ev(a, b, d=2000):
    return sum(((b.Chromosome == x.Chromosome) & (b.Event == x.Event) & ((b.Region - x.Region).abs() <= d)).any() for x in a.itertuples())


for n, (a, b) in PAIRS.items():
    w1, w2 = sw[sw.Sample_Name == a], sw[sw.Sample_Name == b]
    e1, e2 = ev[ev.Sample == a], ev[ev.Sample == b]
    add("recombination", "switches (1-SNP support)", n, len(w1), len(w2))
    add("recombination", "switches with a switch from the other caller within 1 kb", n, f"{near(w1, w2)}/{len(w1)}", f"{near(w2, w1)}/{len(w2)}")
    for e in ["GC", "CO"]:
        add("recombination", f"{e} events", n, int((e1.Event == e).sum()), int((e2.Event == e).sum()))
    add("recombination", "events matched by the other caller (same type, within 2 kb)", n, f"{near_ev(e1, e2)}/{len(e1)}", f"{near_ev(e2, e1)}/{len(e2)}")

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/pipeline_comparison_per_isolate.csv", index=False)

# ---------- 5. Paired tests across the 4 isolates (numeric metrics) ----------
tests = []
for (section, metric), d in res.groupby(["section", "metric"], sort=False):
    x, y = pd.to_numeric(d.CLC, errors="coerce").values, pd.to_numeric(d.own, errors="coerce").values
    if np.isnan(x).any() or np.isnan(y).any() or metric in ("shared", "Jaccard"):
        continue
    diff = y - x
    row = dict(section=section, metric=metric,
               CLC_values=", ".join(f"{v:.4g}" for v in x), own_values=", ".join(f"{v:.4g}" for v in y),
               mean_CLC=x.mean(), mean_own=y.mean(), mean_own_minus_CLC=diff.mean(),
               direction=("own > CLC in all 4" if (diff > 0).all() else "own < CLC in all 4" if (diff < 0).all() else "mixed"),
               wilcoxon_p=stats.wilcoxon(x, y).pvalue if (diff != 0).any() else 1.0,
               paired_t_p=stats.ttest_rel(x, y).pvalue if (diff != 0).any() and np.std(diff) > 0 else np.nan)
    # pooled count test for count metrics (sum over isolates; same genomes, so equal rates expected)
    if metric in ("mutations", "only in this caller", "switches (1-SNP support)", "GC events", "CO events",
                  "clusters (>=3 mutations)", "clustered mutations", "parental SNPs genotyped"):
        row["pooled_counts"] = f"{int(x.sum())} vs {int(y.sum())}"
        row["pooled_binomial_p"] = count_test(x.sum(), y.sum())
    tests.append(row)

# pooled distribution / proportion tests (all 4 isolates together)
pooled = []
chi = stats.chi2_contingency(spec.loc[:, (spec > 0).sum() > 0])
pooled.append(dict(test="substitution spectrum, CLC vs own (chi-square)", p_value=chi[1], note=f"chi2={chi[0]:.2f}, dof={chi[2]}"))
m_clc, m_own = mut[mut.genotype == G_CLC], mut[mut.genotype == G_OWN]
pooled.append(dict(test="A3A-type fraction, CLC vs own (Fisher exact)",
                   p_value=stats.fisher_exact([[a3a(m_clc).sum(), (~a3a(m_clc)).sum()], [a3a(m_own).sum(), (~a3a(m_own)).sum()]])[1],
                   note=f"CLC {100 * a3a(m_clc).mean():.1f}%, own {100 * a3a(m_own).mean():.1f}%"))
c = clusters["0.0001"]
cg = c.assign(group=np.where(c["sample"].isin([p[0] for p in PAIRS.values()]), "CLC", "own"))
for col in ["Length", "Mutations"]:
    x, y = cg[cg.group == "CLC"][col], cg[cg.group == "own"][col]
    if len(x) and len(y):
        pooled.append(dict(test=f"cluster {col.lower()} (p=0.0001), CLC vs own (Mann-Whitney U)", p_value=stats.mannwhitneyu(x, y).pvalue,
                           note=f"median CLC {x.median():g} (n={len(x)}), own {y.median():g} (n={len(y)})"))
oc = pd.concat(own_only_cov); shc = pd.concat(shared_cov)
if len(oc) and len(shc):
    pooled.append(dict(test="coverage: own-only vs shared mutations, own caller (Mann-Whitney U)", p_value=stats.mannwhitneyu(oc, shc).pvalue,
                       note=f"median own-only {oc.median():g} (n={len(oc)}), shared {shc.median():g} (n={len(shc)})"))
meta = pd.read_csv("outputs/cluster_length_meta_df_JT.csv", sep="\t", encoding="utf-16")
meta = meta[meta.cluster_p_val == 0.0001] if "cluster_p_val" in meta else meta
for assoc in ["all_clusters", "rec_associated", "not_associated"]:
    t = pd.crosstab(meta[meta.cluster_association == assoc].genotype, meta[meta.cluster_association == assoc].cluster_type)
    t.to_csv(f"{OUT}/cluster_types_counts_{assoc}.csv")
    if t.shape[0] == 2 and t.shape[1] >= 2:
        p = stats.fisher_exact(t.values)[1] if t.shape[1] == 2 else stats.chi2_contingency(t.values)[1]
        pooled.append(dict(test=f"cluster coordination types ({assoc}, p=0.0001), CLC vs own", p_value=p,
                           note="; ".join(f"{g}: " + ", ".join(f"{k} {int(v)}" for k, v in t.loc[g].items()) for g in t.index)))
evg = ev.assign(group=np.where(ev.Sample.isin([p[0] for p in PAIRS.values()]), "CLC", "own"))
t = pd.crosstab(evg.group, evg.Event)
if t.shape == (2, 2):
    pooled.append(dict(test="GC:CO proportion, CLC vs own (Fisher exact)", p_value=stats.fisher_exact(t.values)[1],
                       note="; ".join(f"{g}: GC {t.loc[g, 'GC']}, CO {t.loc[g, 'CO']}" for g in t.index)))
tests_df = pd.concat([pd.DataFrame(tests), pd.DataFrame(pooled).rename(columns={"test": "metric"}).assign(section="pooled over 4 isolates")], ignore_index=True)
tests_df.to_csv(f"{OUT}/pipeline_comparison_paired_tests.csv", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80); pd.set_option("display.max_columns", 20)
print(res.to_string())
print(tests_df.to_string())
print(dfc.to_string())
if len(ooc):
    print(ooc[["isolate", "Chromosome", "Start", "End", "Length", "Mutations", "median_coverage_of_mutations"]].to_string())

# ---------- Figures ----------
def paired_panel(ax, metric, title=None, section=None):
    d = res[(res.metric == metric) & ((res.section == section) if section else True)]
    x, y = pd.to_numeric(d.CLC).values, pd.to_numeric(d.own).values
    for i, (a, b) in enumerate(zip(x, y)):
        ax.plot([0, 1], [a, b], color="gray", lw=1, zorder=1)
        ax.text(1.05, b, d.isolate.iloc[i], fontsize=7, va="center")
    ax.scatter(np.zeros(len(x)), x, color=C_CLC, s=40, zorder=2)
    ax.scatter(np.ones(len(y)), y, color=C_OWN, s=40, zorder=2)
    t = tests_df[tests_df.metric == metric]
    p = t.paired_t_p.iloc[0] if len(t) else np.nan
    ax.set_xticks([0, 1]); ax.set_xticklabels(["CLC", "own caller"]); ax.set_xlim(-0.3, 1.45)
    ptxt = "identical in all 4" if np.allclose(x, y) else f"paired t p={p:.2g}"
    ax.set_title(f"{title or metric}\n({ptxt})", fontsize=10)


panels = [("parental SNPs genotyped", None), ("APSC (avg per-SNP coverage)", "APSC"),
          ("mutations", "final mutations"), ("% C>T or G>A (A3A-type)", "% A3A-type (C>T, G>A)"),
          ("clusters (>=3 mutations)", "JT clusters, p=0.0001"), ("clustered mutations", "clustered mutations, p=0.0001"),
          ("switches (1-SNP support)", "switches"), ("GC events", None), ("CO events", None)]
fig, axes = plt.subplots(3, 4, figsize=(19, 14))
for ax, (metric, title) in zip(axes.flat, panels):
    paired_panel(ax, metric, title)
# overlap of final mutations per isolate
ax = axes.flat[9]
d = res[res.metric == "only in this caller"].reset_index(drop=True)
sh = pd.to_numeric(res[res.metric == "shared"].CLC).values
x = np.arange(len(d))
ax.bar(x, pd.to_numeric(d.CLC), color=C_CLC, label="only CLC")
ax.bar(x, sh, bottom=pd.to_numeric(d.CLC), color="gray", label="shared")
ax.bar(x, pd.to_numeric(d.own), bottom=pd.to_numeric(d.CLC) + sh, color=C_OWN, label="only own caller")
ax.set_xticks(x); ax.set_xticklabels(d.isolate); ax.legend(fontsize=8); ax.set_title("final mutations: overlap per isolate", fontsize=10)
# pooled spectrum
ax = axes.flat[10]
sp = spec.div(spec.sum(1), axis=0) * 100
sp = sp.reindex([G_CLC, G_OWN])
sp.T.plot.bar(ax=ax, color=[C_CLC, C_OWN]); ax.set_ylabel("% of mutations"); ax.legend(["CLC", "own caller"], fontsize=8)
ax.set_title("substitution spectrum (pooled)", fontsize=10)
# coverage of own-only vs shared mutations
ax = axes.flat[11]
ax.boxplot([shc, oc], labels=["shared with CLC", "own caller only"], showfliers=False)
for j, v in enumerate([shc, oc]):
    ax.scatter(np.random.normal(j + 1, 0.05, len(v)), v, s=6, alpha=0.4, color=["gray", C_OWN][j], zorder=3)
ax.set_ylabel("coverage"); ax.set_title("own caller: coverage of its final mutations", fontsize=10)
plt.suptitle("MayoMate outputs: CLC vs own caller, same 4 ung1∆ isolates (JL1123-JL1126)", fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUT}/pipeline_comparison_overview.png", dpi=200)

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
ax = axes[0]
agg = dfc.groupby("cluster_pval", sort=False)[["clusters_CLC", "clusters_own"]].sum()
x = np.arange(len(agg)); wd = 0.38
ax.bar(x - wd / 2, agg.clusters_CLC, wd, color=C_CLC, label="CLC")
ax.bar(x + wd / 2, agg.clusters_own, wd, color=C_OWN, label="own caller")
ax.set_xticks(x); ax.set_xticklabels(agg.index); ax.set_xlabel("JT cluster p-value threshold"); ax.set_ylabel("clusters (4 isolates)")
ax.legend(); ax.set_title("JT clusters per threshold (sum of 4 isolates)")
for i, (col, lab) in enumerate([("Length", "cluster length (bp)"), ("Mutations", "mutations per cluster")], start=1):
    ax = axes[i]
    d = [cg[cg.group == "CLC"][col], cg[cg.group == "own"][col]]
    ax.boxplot(d, labels=["CLC", "own caller"], showfliers=False)
    for j, (vals, col_) in enumerate(zip(d, [C_CLC, C_OWN])):
        ax.scatter(np.random.normal(j + 1, 0.05, len(vals)), vals, color=col_, s=20, zorder=3)
    p = stats.mannwhitneyu(*d).pvalue if len(d[0]) and len(d[1]) else np.nan
    ax.set_title(f"{lab}, p=0.0001 clusters (MWU p={p:.2g})"); ax.set_ylabel(lab)
plt.tight_layout()
plt.savefig(f"{OUT}/pipeline_comparison_clusters.png", dpi=200)
print("saved", f"{OUT}/pipeline_comparison_overview.png", f"{OUT}/pipeline_comparison_clusters.png")
