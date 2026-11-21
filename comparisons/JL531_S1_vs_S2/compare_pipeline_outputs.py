"""S1 (CLC) vs S2 (own script) comparison of MayoMate notebook outputs for JL531.
Run after main_playground_clean.ipynb, cluster_analysis.py and cluster_coordination.py in this folder.
Writes comparison_outputs/pipeline_comparison_stats.csv, per-level tables and figures."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from mayo.base_functions import load_payload_from_pickle

OUT = "comparison_outputs"
os.makedirs(OUT, exist_ok=True)
S1, S2 = "JL531_S1", "JL531_S2"
C1, C2 = "tab:blue", "tab:orange"
PVALS = ["0.01", "0.001", "0.0001", "1e-05", "1e-06"]
rows = []


def add(section, metric, v1, v2, test="", p=np.nan, note=""):
    rows.append(dict(section=section, metric=metric, S1=v1, S2=v2, test=test, p_value=p, note=note))


def count_test(k1, k2):
    """Exact test that two counts come from equal Poisson rates (same genome, same exposure)."""
    return stats.binomtest(int(k1), int(k1 + k2), 0.5).pvalue if (k1 + k2) else np.nan


def overlap(a, b, key):
    sa, sb = set(map(tuple, a[key].values)), set(map(tuple, b[key].values))
    return len(sa & sb), len(sa - sb), len(sb - sa), len(sa & sb) / max(1, len(sa | sb))


# ---------- 1. Parental SNPs / haplotype assignment ----------
par = {s.name: s for s in load_payload_from_pickle("outputs/partial_pickles/samples_parental_SNPs_1snp.pkl")}
p1 = par[S1].output_df[["Chromosome", "Region", "Parent"]]
p2 = par[S2].output_df[["Chromosome", "Region", "Parent"]]
pm = p1.merge(p2, on=["Chromosome", "Region"], suffixes=("_S1", "_S2"))
add("parental SNPs", "parental SNPs genotyped", len(p1), len(p2), "binomial (equal counts)", count_test(len(p1), len(p2)))
add("parental SNPs", "shared positions", len(pm), len(pm))
add("parental SNPs", "same parent assigned (% of shared)", 100 * (pm.Parent_S1 == pm.Parent_S2).mean(), "",
    note=f"{(pm.Parent_S1 != pm.Parent_S2).sum()} discordant positions")
add("parental SNPs", "APSC (avg per-SNP coverage)", par[S1].apsc, par[S2].apsc)
add("parental SNPs", "% parental (largest parent)", par[S1].percent, par[S2].percent)
add("parental SNPs", "aneuploid chromosomes", ",".join(par[S1].aneuploid_chromosomes) or "none",
    ",".join(par[S2].aneuploid_chromosomes) or "none")
add("parental SNPs", "heterozygous chromosomes", ",".join(par[S1].heterozygous_chromosomes) or "none",
    ",".join(par[S2].heterozygous_chromosomes) or "none")
pm[pm.Parent_S1 != pm.Parent_S2].to_csv(f"{OUT}/parental_SNPs_discordant.csv", index=False)

# ---------- 2. Final filtered mutations ----------
mut = pd.read_csv("outputs/all_mutations.tsv", sep="\t")
mut["sub"] = mut.Reference + ">" + mut.Allele
m1, m2 = mut[mut.sample_name == S1], mut[mut.sample_name == S2]
key = ["Chromosome", "Region", "Reference", "Allele"]
sh, o1, o2, jac = overlap(m1, m2, key)
add("mutations (after all filters)", "mutations", len(m1), len(m2), "binomial (equal counts)", count_test(len(m1), len(m2)))
add("mutations (after all filters)", "shared / only S1 / only S2", f"{sh} / {o1} / -", f"{sh} / - / {o2}", note=f"Jaccard = {jac:.3f}")
add("mutations (after all filters)", "% of own mutations found in other", 100 * sh / len(m1), 100 * sh / len(m2))
a3a = lambda d: d["sub"].isin(["C>T", "G>A"])
add("mutations (after all filters)", "% C>T or G>A (A3A-type)", 100 * a3a(m1).mean(), 100 * a3a(m2).mean(), "Fisher exact",
    stats.fisher_exact([[a3a(m1).sum(), (~a3a(m1)).sum()], [a3a(m2).sum(), (~a3a(m2)).sum()]])[1])
spec = pd.crosstab(mut.sample_name, mut["sub"])
spec.to_csv(f"{OUT}/mutation_spectrum.csv")
chi = stats.chi2_contingency(spec.loc[:, (spec > 0).sum() > 0])
add("mutations (after all filters)", "substitution spectrum (12 classes)", "", "", "chi-square", chi[1], f"chi2={chi[0]:.2f}, dof={chi[2]}")
mm = m1.merge(m2, on=key, suffixes=("_S1", "_S2"))
add("mutations (after all filters)", "allele frequency, shared mutations (median)", mm.Frequency_S1.median(), mm.Frequency_S2.median(),
    "Wilcoxon signed-rank (paired)", stats.wilcoxon(mm.Frequency_S1, mm.Frequency_S2).pvalue if (mm.Frequency_S1 != mm.Frequency_S2).any() else 1.0)
add("mutations (after all filters)", "coverage, shared mutations (median)", mm.Coverage_S1.median(), mm.Coverage_S2.median(),
    "Wilcoxon signed-rank (paired)", stats.wilcoxon(mm.Coverage_S1, mm.Coverage_S2).pvalue)
allm = m1.merge(m2, on=key, how="outer", suffixes=("_S1", "_S2"), indicator=True)
allm.to_csv(f"{OUT}/mutations_S1_vs_S2.csv", index=False)

# ---------- 3. JT clusters at each significance threshold ----------
clust_rows = []
for pv in PVALS:
    c = pd.read_csv(f"data/new_clusters/mutation_clusters/mutation_clusters_pval_{pv}.tsv", sep="\t")
    c = c[c.Mutations >= 3]
    c1, c2 = c[c["sample"] == S1], c[c["sample"] == S2]
    # a cluster is "matched" if it overlaps any cluster of the other sample on the same chromosome
    def matched(a, b):
        return sum(((b.Chromosome == r.Chromosome) & (b.Start <= r.End) & (b.End >= r.Start)).any() for r in a.itertuples())
    mt1, mt2 = matched(c1, c2), matched(c2, c1)
    def mw(col):
        return stats.mannwhitneyu(c1[col], c2[col]).pvalue if len(c1) and len(c2) else np.nan
    def ks(col):
        return stats.ks_2samp(c1[col], c2[col]).pvalue if len(c1) and len(c2) else np.nan
    clust_rows.append(dict(
        cluster_pval=pv, clusters_S1=len(c1), clusters_S2=len(c2), p_count_binomial=count_test(len(c1), len(c2)),
        S1_clusters_overlapping_S2=mt1, S2_clusters_overlapping_S1=mt2,
        median_length_S1=c1.Length.median(), median_length_S2=c2.Length.median(),
        p_length_mannwhitney=mw("Length"), p_length_ks=ks("Length"),
        median_mutations_S1=c1.Mutations.median(), median_mutations_S2=c2.Mutations.median(),
        p_mutations_mannwhitney=mw("Mutations"),
        clustered_mutations_S1=c1.Mutations.sum(), clustered_mutations_S2=c2.Mutations.sum(),
        ssDNA_kb_S1=c1.Length.sum() / 1000, ssDNA_kb_S2=c2.Length.sum() / 1000))
dfc = pd.DataFrame(clust_rows)
dfc.to_csv(f"{OUT}/clusters_by_pval_S1_vs_S2.csv", index=False)
r = dfc[dfc.cluster_pval == "0.0001"].iloc[0]
add("JT clusters (p=0.0001)", "clusters (>=3 mutations)", r.clusters_S1, r.clusters_S2, "binomial (equal counts)", r.p_count_binomial)
add("JT clusters (p=0.0001)", "clusters overlapping a cluster in the other sample", f"{r.S1_clusters_overlapping_S2}/{r.clusters_S1}", f"{r.S2_clusters_overlapping_S1}/{r.clusters_S2}")
add("JT clusters (p=0.0001)", "cluster length, median bp", r.median_length_S1, r.median_length_S2, "Mann-Whitney U", r.p_length_mannwhitney, f"KS p={r.p_length_ks:.3g}")
add("JT clusters (p=0.0001)", "mutations per cluster, median", r.median_mutations_S1, r.median_mutations_S2, "Mann-Whitney U", r.p_mutations_mannwhitney)
add("JT clusters (p=0.0001)", "total clustered ssDNA (kb)", r.ssDNA_kb_S1, r.ssDNA_kb_S2)

# cluster coordination (switching pattern) categories: Fisher/chi-square on counts
meta = pd.read_csv("outputs/cluster_length_meta_df_JT.csv", sep="\t", encoding="utf-16")
for assoc in ["all_clusters", "rec_associated", "not_associated"]:
    t = pd.crosstab(meta[meta.cluster_association == assoc].genotype, meta[meta.cluster_association == assoc].cluster_type)
    t.to_csv(f"{OUT}/cluster_types_counts_{assoc}.csv")
    if t.shape == (2, 2):
        p = stats.fisher_exact(t.values)[1]; test = "Fisher exact"
    elif t.shape[0] == 2 and t.shape[1] > 2:
        p = stats.chi2_contingency(t.values)[1]; test = "chi-square (low counts!)"
    else:
        p, test = np.nan, "n/a"
    fmt = lambda g: ", ".join(f"{k}: {int(v)}" for k, v in t.loc[g].items()) if g in t.index else "none"
    add("cluster coordination types", f"cluster type counts ({assoc})", fmt("JL531-S1"), fmt("JL531-S2"), test, p)

# ---------- 4. Switches and recombination events ----------
sw = pd.read_csv("outputs/all_switches_1snp.tsv", sep="\t")
w1, w2 = sw[sw.Sample_Name == S1], sw[sw.Sample_Name == S2]
add("recombination", "switches (1-SNP support)", len(w1), len(w2), "binomial (equal counts)", count_test(len(w1), len(w2)))
def near(a, b, d=1000):
    return sum(((b.Chromosome == x.Chromosome) & ((b.Switch_Center - x.Switch_Center).abs() <= d)).any() for x in a.itertuples())
add("recombination", "switches with a switch in the other sample within 1 kb", f"{near(w1, w2)}/{len(w1)}", f"{near(w2, w1)}/{len(w2)}")
ev = pd.read_csv("outputs/all_recombination_events_haploid_chromosomes_1snp.csv", sep="\t", encoding="utf-16")
for e in ["GC", "CO"]:
    k1, k2 = ((ev.Sample == S1) & (ev.Event == e)).sum(), ((ev.Sample == S2) & (ev.Event == e)).sum()
    add("recombination", f"{e} events", k1, k2, "binomial (equal counts)", count_test(k1, k2))
t = pd.crosstab(ev.Sample, ev.Event)
add("recombination", "GC:CO proportion", "", "", "Fisher exact" if t.shape == (2, 2) else "chi-square",
    stats.fisher_exact(t.values)[1] if t.shape == (2, 2) else stats.chi2_contingency(t.values)[1])
e1, e2 = ev[ev.Sample == S1], ev[ev.Sample == S2]
def near_ev(a, b, d=2000):
    return sum(((b.Chromosome == x.Chromosome) & (b.Event == x.Event) & ((b.Region - x.Region).abs() <= d)).any() for x in a.itertuples())
add("recombination", "events matched in other sample (same type, within 2 kb)", f"{near_ev(e1, e2)}/{len(e1)}", f"{near_ev(e2, e1)}/{len(e2)}")

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/pipeline_comparison_stats.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70); pd.set_option("display.max_columns", 20)
print(res.to_string())
print(dfc.to_string())

# ---------- Figures ----------
fig, ax = plt.subplots(2, 3, figsize=(17, 10))
a = ax[0, 0]
a.bar(["only S1", "shared", "only S2"], [o1, sh, o2], color=[C1, "gray", C2])
a.set_title(f"Final mutations (Jaccard {jac:.2f})"); a.set_ylabel("mutations")
a = ax[0, 1]
spec_pct = spec.div(spec.sum(1), axis=0) * 100
spec_pct.T.plot.bar(ax=a, color=[C1, C2]); a.set_ylabel("% of mutations"); a.set_title("Substitution spectrum")
a = ax[0, 2]
x = np.arange(len(dfc)); wd = 0.38
a.bar(x - wd / 2, dfc.clusters_S1, wd, color=C1, label="S1 (CLC)")
a.bar(x + wd / 2, dfc.clusters_S2, wd, color=C2, label="S2 (own)")
a.set_xticks(x); a.set_xticklabels(dfc.cluster_pval); a.set_xlabel("JT cluster p-value threshold"); a.set_ylabel("clusters"); a.legend()
a.set_title("JT clusters per threshold")
c = pd.read_csv("data/new_clusters/mutation_clusters/mutation_clusters_pval_0.0001.tsv", sep="\t"); c = c[c.Mutations >= 3]
for i, (col, lab) in enumerate([("Length", "cluster length (bp)"), ("Mutations", "mutations per cluster")]):
    a = ax[1, i]
    d = [c[c["sample"] == S1][col], c[c["sample"] == S2][col]]
    a.boxplot(d, labels=["S1 (CLC)", "S2 (own)"], showfliers=False)
    for j, (vals, col_) in enumerate(zip(d, [C1, C2])):
        a.scatter(np.random.normal(j + 1, 0.05, len(vals)), vals, color=col_, s=25, zorder=3)
    p = stats.mannwhitneyu(*d).pvalue
    a.set_title(f"{lab}, p=0.0001 clusters (MWU p={p:.2f})"); a.set_ylabel(lab)
a = ax[1, 2]
cats = ["switches", "GC", "CO", "mutations/10", "parental SNPs/1000"]
v1 = [len(w1), (e1.Event == "GC").sum(), (e1.Event == "CO").sum(), len(m1) / 10, len(p1) / 1000]
v2 = [len(w2), (e2.Event == "GC").sum(), (e2.Event == "CO").sum(), len(m2) / 10, len(p2) / 1000]
x = np.arange(len(cats))
a.bar(x - wd / 2, v1, wd, color=C1, label="S1 (CLC)"); a.bar(x + wd / 2, v2, wd, color=C2, label="S2 (own)")
a.set_xticks(x); a.set_xticklabels(cats, rotation=20); a.legend(); a.set_title("Recombination and totals")
plt.suptitle("JL531: S1 (CLC) vs S2 (own variant caller), MayoMate outputs", fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUT}/pipeline_comparison_overview.png", dpi=200)
print("saved", f"{OUT}/pipeline_comparison_overview.png")
