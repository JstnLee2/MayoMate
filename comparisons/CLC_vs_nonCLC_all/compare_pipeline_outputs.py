"""CLC vs own caller comparison of MayoMate outputs for every paired isolate (pairs.json), kept
separate by genotype. Run after the notebook and run_scripts.py in this folder.

Writes to comparison_outputs/:
  pipeline_per_isolate.csv         one row per isolate: every metric for CLC and own caller
  genotype_summary.csv             per genotype and caller: mean / median of every metric
  paired_tests.csv                 own vs CLC (paired Wilcoxon), overall and per genotype, plus
                                   Kruskal-Wallis: does the caller difference depend on the genotype?
  pooled_tests.csv                 per genotype and overall: spectrum, A3A fraction, cluster length /
                                   size, GC:CO, cluster types, coverage of caller-specific mutations
  caller_specific_clusters.csv     clusters found by only one caller, with the depth of their mutations
  genotype_effects_by_caller.csv   every "genotype vs ung1∆" test from the scripts, CLC set vs own set,
                                   and whether the conclusion (p<0.05) is the same
  scripts_CLC_vs_own_by_genotype.csv  the scripts' own CLC-vs-own tests, one run per genotype
  figures: paired_differences_by_genotype.png, CLC_vs_own_scatter.png,
           genotype_effects_by_caller.png, caller_specific_mutations.png
"""
import glob
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from mayo.base_functions import load_payload_from_pickle
from mayo.settings.genotype_colors import genotype_color

OUT = "comparison_outputs"
os.makedirs(OUT, exist_ok=True)
PAIRS = json.load(open("pairs.json", encoding="utf-8"))["pairs"]
GENO_ORDER = ["ung1∆", "ung1∆msh6∆", "ung1∆_NonSelected", "ung1∆msh6∆_NonSelected", "ung1∆ 37C_Selected",
              "ung1∆ 37C_NonSelected", "ung1∆_zymolyaseheat", "ung1∆SNPless", "ung1∆SNPlessmsh6∆", "ung1∆SNPlessmsh3∆"]
GENO_ORDER = [g for g in GENO_ORDER if g in {p["genotype"] for p in PAIRS.values()}]
SIG = "0.0001"
KEY = ["Chromosome", "Region", "Reference", "Allele"]
C_CLC, C_OWN = "tab:blue", "tab:orange"
CLC_NAMES = {p["CLC"] for p in PAIRS.values()}


def safe(f, default=np.nan):
    try:
        return f()
    except Exception:
        return default


def count_test(k1, k2):
    return stats.binomtest(int(k1), int(k1 + k2), 0.5).pvalue if (k1 + k2) else np.nan


def overlapping(a, b):
    return np.array([((b.Chromosome == r.Chromosome) & (b.Start <= r.End) & (b.End >= r.Start)).any() for r in a.itertuples()], dtype=bool)


# ---------- inputs ----------
par = {s.name: s for s in load_payload_from_pickle("outputs/partial_pickles/samples_parental_SNPs_1snp.pkl")}
mut = pd.read_csv("outputs/all_mutations.tsv", sep="\t")
mut["sub"] = mut.Reference + ">" + mut.Allele
clu = pd.read_csv(f"data/new_clusters/mutation_clusters/mutation_clusters_pval_{SIG}.tsv", sep="\t")
clu = clu[clu.Mutations >= 3]
sw = pd.read_csv("outputs/all_switches_1snp.tsv", sep="\t")
ev = pd.read_csv("outputs/all_recombination_events_haploid_chromosomes_1snp.csv", sep="\t", encoding="utf-16")
a3a = lambda d: d["sub"].isin(["C>T", "G>A"])

# ---------- per isolate ----------
rows, caller_specific_mut, caller_specific_clu = [], [], []
for iso, p in PAIRS.items():
    a, b, g = p["CLC"], p["own"], p["genotype"]
    r = dict(isolate=iso, genotype=g, CLC_sample=a, own_sample=b,
             in_parental_data=(a in par) and (b in par),
             in_mutation_data=(a in set(mut.sample_name)) and (b in set(mut.sample_name)))
    # parental SNPs
    for tag, s in [("CLC", a), ("own", b)]:
        o = par.get(s)
        r[f"parental_SNPs_{tag}"] = safe(lambda: len(o.output_df))
        r[f"APSC_{tag}"] = safe(lambda: o.apsc)
        r[f"pct_parental_{tag}"] = safe(lambda: o.percent)
        r[f"aneuploid_chr_{tag}"] = safe(lambda: ",".join(o.aneuploid_chromosomes) or "none", "")
        r[f"heterozygous_chr_{tag}"] = safe(lambda: ",".join(o.heterozygous_chromosomes) or "none", "")
    def parent_concordance():
        pa, pb = par[a].output_df[["Chromosome", "Region", "Parent"]], par[b].output_df[["Chromosome", "Region", "Parent"]]
        pm = pa.merge(pb, on=["Chromosome", "Region"])
        return 100 * (pm.Parent_x == pm.Parent_y).mean() if len(pm) else np.nan
    r["parent_concordance_pct"] = safe(parent_concordance)
    r["same_ploidy_calls"] = (r["aneuploid_chr_CLC"] == r["aneuploid_chr_own"]) and (r["heterozygous_chr_CLC"] == r["heterozygous_chr_own"])
    # final mutations
    m1, m2 = mut[mut.sample_name == a], mut[mut.sample_name == b]
    allm = m1.merge(m2, on=KEY, how="outer", suffixes=("_CLC", "_own"), indicator=True)
    sh = int((allm._merge == "both").sum())
    r.update(mutations_CLC=len(m1), mutations_own=len(m2), mutations_shared=sh,
             only_CLC=int((allm._merge == "left_only").sum()), only_own=int((allm._merge == "right_only").sum()),
             mutation_jaccard=sh / len(allm) if len(allm) else np.nan,
             pct_A3A_CLC=100 * a3a(m1).mean() if len(m1) else np.nan, pct_A3A_own=100 * a3a(m2).mean() if len(m2) else np.nan,
             mut_coverage_CLC=m1.Coverage.median(), mut_coverage_own=m2.Coverage.median())
    for kind, tag, cov in [("left_only", "CLC only", "Coverage_CLC"), ("right_only", "own only", "Coverage_own"), ("both", "shared", "Coverage_own")]:
        x = allm[allm._merge == kind]
        caller_specific_mut.append(pd.DataFrame(dict(isolate=iso, genotype=g, set=tag, coverage=x[cov].values,
                                                      sub=(x[f"sub_CLC" if kind == "left_only" else "sub_own"]).values)))
    # clusters
    c1, c2 = clu[clu["sample"] == a], clu[clu["sample"] == b]
    o1, o2 = overlapping(c1, c2), overlapping(c2, c1)
    r.update(clusters_CLC=len(c1), clusters_own=len(c2), clusters_matched_CLC=int(o1.sum()), clusters_matched_own=int(o2.sum()),
             clustered_mut_CLC=int(c1.Mutations.sum()), clustered_mut_own=int(c2.Mutations.sum()),
             ssDNA_kb_CLC=c1.Length.sum() / 1000, ssDNA_kb_own=c2.Length.sum() / 1000)
    for tag, cc, oo, s in [("CLC only", c1, o1, a), ("own only", c2, o2, b)]:
        for cr in cc[~oo].itertuples():
            mm = mut[(mut.sample_name == s) & (mut.Chromosome == cr.Chromosome) & mut.Region.between(cr.Start, cr.End)]
            caller_specific_clu.append(dict(isolate=iso, genotype=g, found_by=tag, Chromosome=cr.Chromosome, Start=cr.Start, End=cr.End,
                                            Length=cr.Length, Mutations=cr.Mutations, median_coverage=mm.Coverage.median(),
                                            pct_A3A=100 * a3a(mm).mean() if len(mm) else np.nan))
    # recombination
    w1, w2 = sw[sw.Sample_Name == a], sw[sw.Sample_Name == b]
    e1, e2 = ev[ev.Sample == a], ev[ev.Sample == b]
    r.update(switches_CLC=len(w1), switches_own=len(w2),
             GC_CLC=int((e1.Event == "GC").sum()), GC_own=int((e2.Event == "GC").sum()),
             CO_CLC=int((e1.Event == "CO").sum()), CO_own=int((e2.Event == "CO").sum()),
             events_matched_CLC=sum(((e2.Chromosome == x.Chromosome) & (e2.Event == x.Event) & ((e2.Region - x.Region).abs() <= 2000)).any() for x in e1.itertuples()),
             events_matched_own=sum(((e1.Chromosome == x.Chromosome) & (e1.Event == x.Event) & ((e1.Region - x.Region).abs() <= 2000)).any() for x in e2.itertuples()))
    rows.append(r)

df = pd.DataFrame(rows)
df["genotype"] = pd.Categorical(df.genotype, GENO_ORDER, ordered=True)
df = df.sort_values(["genotype", "isolate"])
# Isolates with a heterozygous chromosome (by either caller): CLC switches between parents far more often there,
# so switch/GC/CO differences are reported separately for these isolates.
df["has_het_chr"] = df.heterozygous_chr_CLC.fillna("none").ne("none") | df.heterozygous_chr_own.fillna("none").ne("none")


def switches_on(sample, chroms, on):
    s = sw[sw.Sample_Name == sample]
    return int(s.Chromosome.isin(chroms).sum() if on else (~s.Chromosome.isin(chroms)).sum())


for tag in ("CLC", "own"):
    het_sets = [{c.strip() for c in f"{a},{b}".split(",") if c.strip() not in ("", "none", "nan")}
                for a, b in zip(df.heterozygous_chr_CLC.astype(str), df.heterozygous_chr_own.astype(str))]
    df[f"switches_het_chr_{tag}"] = [switches_on(s, h, True) for s, h in zip(df[f"{tag}_sample"], het_sets)]
    df[f"switches_other_chr_{tag}"] = [switches_on(s, h, False) for s, h in zip(df[f"{tag}_sample"], het_sets)]
df.to_csv(f"{OUT}/pipeline_per_isolate.csv", index=False)
csm = pd.concat(caller_specific_mut, ignore_index=True)
csc = pd.DataFrame(caller_specific_clu)
csc.to_csv(f"{OUT}/caller_specific_clusters.csv", index=False)

METRICS = ["parental_SNPs", "APSC", "pct_parental", "mutations", "pct_A3A", "mut_coverage", "clusters", "clustered_mut",
           "ssDNA_kb", "switches", "GC", "CO"]

# ---------- genotype summary (genotypes kept separate) ----------
summ = []
for g in GENO_ORDER:
    d = df[df.genotype == g]
    for tag in ["CLC", "own"]:
        row = dict(genotype=g + ("CLC" if tag == "CLC" else ""), caller=tag, isolates=len(d))
        for m in METRICS:
            row[f"{m}_mean"] = d[f"{m}_{tag}"].mean()
            row[f"{m}_median"] = d[f"{m}_{tag}"].median()
        summ.append(row)
    summ.append(dict(genotype=g + " (own - CLC)", caller="difference", isolates=len(d),
                     **{f"{m}_mean": (d[f"{m}_own"] - d[f"{m}_CLC"]).mean() for m in METRICS},
                     **{f"{m}_median": (d[f"{m}_own"] - d[f"{m}_CLC"]).median() for m in METRICS}))
pd.DataFrame(summ).to_csv(f"{OUT}/genotype_summary.csv", index=False)

# ---------- paired tests ----------
tests = []
STRATA = [("ALL", df), ("ALL, no heterozygous chr", df[~df.has_het_chr]), ("ALL, with heterozygous chr", df[df.has_het_chr])]
for m in METRICS + ["switches_het_chr", "switches_other_chr"]:
    for g, d in STRATA + [(g, df[df.genotype == g]) for g in GENO_ORDER]:
        dd = d[[f"{m}_CLC", f"{m}_own"]].dropna()
        x, y = dd[f"{m}_CLC"].values, dd[f"{m}_own"].values
        diff = y - x
        tests.append(dict(metric=m, genotype=g, n=len(dd), mean_CLC=x.mean() if len(x) else np.nan, mean_own=y.mean() if len(y) else np.nan,
                          mean_own_minus_CLC=diff.mean() if len(diff) else np.nan,
                          pct_own_higher=100 * (diff > 0).mean() if len(diff) else np.nan,
                          pct_equal=100 * (diff == 0).mean() if len(diff) else np.nan,
                          pct_own_lower=100 * (diff < 0).mean() if len(diff) else np.nan,
                          wilcoxon_p=safe(lambda: stats.wilcoxon(x, y).pvalue) if (diff != 0).any() else 1.0,
                          pooled_binomial_p=count_test(x.sum(), y.sum()) if m in ("mutations", "clusters", "clustered_mut", "switches", "GC", "CO") else np.nan))
    diffs = [(df[df.genotype == g][f"{m}_own"] - df[df.genotype == g][f"{m}_CLC"]).dropna().values for g in GENO_ORDER]
    diffs = [x for x in diffs if len(x)]
    tests.append(dict(metric=m, genotype="between genotypes (Kruskal-Wallis on own-CLC)", n=sum(map(len, diffs)),
                      kruskal_p=safe(lambda: stats.kruskal(*diffs).pvalue)))
tests = pd.DataFrame(tests)
ok = tests.genotype != "between genotypes (Kruskal-Wallis on own-CLC)"
# Benjamini-Hochberg over all paired Wilcoxon tests
pv = tests.loc[ok, "wilcoxon_p"].values
order = np.argsort(pv); q = np.empty_like(pv); n = len(pv)
q[order] = np.minimum.accumulate((pv[order] * n / np.arange(1, n + 1))[::-1])[::-1]
tests.loc[ok, "wilcoxon_q_BH"] = np.minimum(q, 1)
tests.to_csv(f"{OUT}/paired_tests.csv", index=False)

# ---------- pooled tests per genotype ----------
pooled = []
mut["caller"] = np.where(mut.sample_name.isin(CLC_NAMES), "CLC", "own")
mut["base_genotype"] = mut.genotype.str.replace(r"CLC$", "", regex=True)
clu["caller"] = np.where(clu["sample"].isin(CLC_NAMES), "CLC", "own")
sample_geno = {s: p["genotype"] for p in PAIRS.values() for s in (p["CLC"], p["own"])}
clu["base_genotype"] = clu["sample"].map(sample_geno)
ev["caller"] = np.where(ev.Sample.isin(CLC_NAMES), "CLC", "own")
ev["base_genotype"] = ev.Sample.map(sample_geno)
for g in ["ALL"] + GENO_ORDER:
    mg = mut if g == "ALL" else mut[mut.base_genotype == g]
    mg = mg[mg.sample_name.isin(sample_geno)]
    sp = pd.crosstab(mg.caller, mg["sub"])
    if sp.shape[0] == 2:
        sp = sp.loc[:, (sp > 0).sum() > 0]
        pooled.append(dict(genotype=g, test="substitution spectrum (chi-square)", p_value=safe(lambda: stats.chi2_contingency(sp)[1])))
        k = [[a3a(mg[mg.caller == c]).sum(), (~a3a(mg[mg.caller == c])).sum()] for c in ["CLC", "own"]]
        pooled.append(dict(genotype=g, test="A3A-type fraction (Fisher exact)", p_value=stats.fisher_exact(k)[1],
                           note=f"CLC {100 * k[0][0] / max(1, sum(k[0])):.1f}%, own {100 * k[1][0] / max(1, sum(k[1])):.1f}%"))
    cg = clu if g == "ALL" else clu[clu.base_genotype == g]
    cg = cg[cg.base_genotype.notna()]
    for col in ["Length", "Mutations"]:
        x, y = cg[cg.caller == "CLC"][col], cg[cg.caller == "own"][col]
        if len(x) and len(y):
            pooled.append(dict(genotype=g, test=f"cluster {col.lower()} (Mann-Whitney U)", p_value=stats.mannwhitneyu(x, y).pvalue,
                               note=f"median CLC {x.median():g} (n={len(x)}), own {y.median():g} (n={len(y)})"))
    eg = ev if g == "ALL" else ev[ev.base_genotype == g]
    t = pd.crosstab(eg.caller, eg.Event)
    if t.shape == (2, 2):
        pooled.append(dict(genotype=g, test="GC:CO proportion (Fisher exact)", p_value=stats.fisher_exact(t.values)[1],
                           note="; ".join(f"{c}: GC {t.loc[c, 'GC']}, CO {t.loc[c, 'CO']}" for c in t.index)))
    cm = csm if g == "ALL" else csm[csm.genotype == g]
    for tag in ["own only", "CLC only"]:
        x, y = cm[cm.set == tag].coverage.dropna(), cm[cm.set == "shared"].coverage.dropna()
        if len(x) and len(y):
            pooled.append(dict(genotype=g, test=f"coverage of {tag} vs shared mutations (Mann-Whitney U)", p_value=stats.mannwhitneyu(x, y).pvalue,
                               note=f"median {tag} {x.median():g} (n={len(x)}), shared {y.median():g} (n={len(y)})"))
meta_files = glob.glob("script_runs/pair_*/outputs/cluster_length_meta_df_JT.csv")
for f in meta_files:
    meta = pd.read_csv(f, sep="\t", encoding="utf-16")
    meta = meta[(meta.cluster_p_val == float(SIG)) & (meta.cluster_association == "all_clusters")]
    t = pd.crosstab(meta.genotype, meta.cluster_type)
    if t.shape[0] == 2 and t.shape[1] >= 2:
        g = [x for x in t.index if not x.endswith("CLC")][0]
        pooled.append(dict(genotype=g, test="cluster coordination types, all clusters (chi-square)", p_value=safe(lambda: stats.chi2_contingency(t.values)[1]),
                           note="; ".join(f"{x}: " + ", ".join(f"{k} {int(v)}" for k, v in t.loc[x].items()) for x in t.index)))
pooled = pd.DataFrame(pooled)
pooled.to_csv(f"{OUT}/pooled_tests.csv", index=False)

# ---------- scripts: CLC vs own per genotype, and genotype effects per caller ----------
def read_any(f):
    for kw in [dict(sep="\t", encoding="utf-16"), dict(sep=","), dict(sep="\t")]:
        try:
            d = pd.read_csv(f, **kw)
            if d.shape[1] > 2:
                return d
        except Exception:
            pass
    return pd.DataFrame()


def script_tests(run, reference):
    # Only the files for this run's reference: a run folder can also hold the previous run's last output file.
    out = []
    for f in glob.glob(f"script_runs/{run}/outputs/clusters/stats_vs_{reference}_*.csv"):
        d = read_any(f)
        for r in d.itertuples():
            out.append(dict(script="cluster_analysis", test=r.metric, genotype=r.genotype, reference=r.reference,
                            median_genotype=r.median_genotype, median_reference=r.median_reference, p_value=r.p_mannwhitney_all))
    for f in glob.glob(f"script_runs/{run}/outputs/cluster_length_meta_test_df_JT_vs_{reference}.csv"):
        d = read_any(f)
        d = d[d.cluster_p_val == float(SIG)]
        for r in d.itertuples():
            out.append(dict(script="cluster_coordination", test=f"cluster length, {r.cluster_type} ({r.cluster_association})",
                            genotype=r.genotype, reference=r.reference, p_value=r.p_val))
    for f in glob.glob(f"script_runs/{run}/outputs/cluster_association_chi2_*_vs_{reference}.csv"):
        d = read_any(f)
        for r in d.itertuples():
            if r.genotype != r.reference:
                out.append(dict(script="cluster_association_chi2", test="recombination-associated clusters (chi-square)",
                                genotype=r.genotype, reference=r.reference, median_genotype=r.percent_associated, p_value=r.p_value))
    return pd.DataFrame(out)


pair_runs = pd.concat([script_tests(os.path.basename(d), os.path.basename(d)[len("pair_"):] + "CLC").assign(run=os.path.basename(d)) for d in sorted(glob.glob("script_runs/pair_*"))], ignore_index=True) \
    if glob.glob("script_runs/pair_*") else pd.DataFrame()
pair_runs.to_csv(f"{OUT}/scripts_CLC_vs_own_by_genotype.csv", index=False)

eff = pd.DataFrame()
if os.path.isdir("script_runs/CLC_set") and os.path.isdir("script_runs/own_set"):
    a = script_tests("CLC_set", "ung1∆CLC"); b = script_tests("own_set", "ung1∆")
    a["base_genotype"] = a.genotype.str.replace(r"CLC$", "", regex=True)
    b["base_genotype"] = b.genotype
    eff = a.merge(b, on=["script", "test", "base_genotype"], suffixes=("_CLC", "_own"), how="outer")
    eff["significant_CLC"] = eff.p_value_CLC < 0.05
    eff["significant_own"] = eff.p_value_own < 0.05
    eff["same_conclusion"] = eff.significant_CLC == eff.significant_own
    eff = eff[["script", "test", "base_genotype", "reference_CLC", "reference_own", "median_genotype_CLC", "median_genotype_own",
               "p_value_CLC", "p_value_own", "significant_CLC", "significant_own", "same_conclusion"]]
    eff.to_csv(f"{OUT}/genotype_effects_by_caller.csv", index=False)

# ---------- print ----------
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 90); pd.set_option("display.max_columns", 30)
print(tests[tests.genotype.isin(["ALL"]) | tests.genotype.str.startswith("between")].round(4).to_string())
print(pooled[pooled.genotype == "ALL"].to_string())
print("caller-specific clusters:", csc.groupby(["found_by"]).size().to_dict() if len(csc) else 0)
if len(eff):
    print("genotype-effect tests:", len(eff), "| same conclusion with both callers:", int(eff.same_conclusion.sum()))
    print(eff[~eff.same_conclusion].to_string())

# ---------- figures ----------
fig_metrics = [("mutations", "final mutations"), ("pct_A3A", "% A3A-type"), ("mut_coverage", "median coverage of mutations"),
               ("clusters", f"JT clusters (p={SIG})"), ("clustered_mut", "clustered mutations"), ("ssDNA_kb", "clustered ssDNA (kb)"),
               ("switches", "switches"), ("GC", "GC events"), ("CO", "CO events")]
fig, axes = plt.subplots(3, 3, figsize=(21, 16))
for ax, (m, title) in zip(axes.flat, fig_metrics):
    data = [(df[df.genotype == g][f"{m}_own"] - df[df.genotype == g][f"{m}_CLC"]).dropna() for g in GENO_ORDER]
    bp = ax.boxplot(data, patch_artist=True, showfliers=False)
    for patch, g in zip(bp["boxes"], GENO_ORDER):
        patch.set_facecolor(genotype_color(g)); patch.set_alpha(0.6)
    for i, v in enumerate(data, start=1):
        ax.scatter(np.random.normal(i, 0.06, len(v)), v, s=8, color="k", alpha=0.5, zorder=3)
    ax.axhline(0, color="gray", ls="--", lw=1)
    t = tests[(tests.metric == m)]
    pall = t[t.genotype == "ALL"].wilcoxon_p.iloc[0]
    pk = t[t.genotype.str.startswith("between")].kruskal_p.iloc[0]
    ax.set_title(f"{title}: own caller minus CLC\nall isolates Wilcoxon p={pall:.2g}; differs by genotype? KW p={pk:.2g}", fontsize=10)
    ax.set_xticks(range(1, len(GENO_ORDER) + 1))
    ax.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(GENO_ORDER, data)], rotation=45, ha="right", fontsize=8)
plt.suptitle(f"Paired differences (own caller - CLC) per isolate, by genotype ({len(df)} isolates)", fontweight="bold")
plt.tight_layout(); plt.savefig(f"{OUT}/paired_differences_by_genotype.png", dpi=180)

fig, axes = plt.subplots(2, 3, figsize=(18, 11))
for ax, (m, title) in zip(axes.flat, [("mutations", "final mutations"), ("clusters", "JT clusters"), ("clustered_mut", "clustered mutations"),
                                      ("GC", "GC events"), ("CO", "CO events"), ("APSC", "APSC")]):
    for g in GENO_ORDER:
        d = df[df.genotype == g]
        ax.scatter(d[f"{m}_CLC"], d[f"{m}_own"], s=18, color=genotype_color(g), label=g, alpha=0.8, edgecolor="k", lw=0.3)
    lim = np.nanmax(df[[f"{m}_CLC", f"{m}_own"]].values) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=1); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("CLC"); ax.set_ylabel("own caller"); ax.set_title(title)
axes.flat[0].legend(fontsize=7)
plt.suptitle("Per isolate: CLC vs own caller", fontweight="bold")
plt.tight_layout(); plt.savefig(f"{OUT}/CLC_vs_own_scatter.png", dpi=180)

fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
ax = axes[0]
data, labels = [], []
for g in GENO_ORDER:
    for tag, col in [("shared", "gray"), ("CLC only", C_CLC), ("own only", C_OWN)]:
        data.append(csm[(csm.genotype == g) & (csm.set == tag)].coverage.dropna()); labels.append((g, tag, col))
bp = ax.boxplot(data, patch_artist=True, showfliers=False, positions=np.arange(len(data)) + np.repeat(np.arange(len(GENO_ORDER)), 3))
for patch, (_, _, col) in zip(bp["boxes"], labels):
    patch.set_facecolor(col); patch.set_alpha(0.7)
ax.set_xticks(np.arange(len(GENO_ORDER)) * 4 + 1); ax.set_xticklabels(GENO_ORDER, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("coverage"); ax.set_title("Coverage of final mutations: shared (gray), CLC only (blue), own only (orange)")
ax = axes[1]
x = np.arange(len(GENO_ORDER)); wd = 0.38
n_iso = [len(df[df.genotype == g]) for g in GENO_ORDER]
for k, (tag, col) in enumerate([("CLC only", C_CLC), ("own only", C_OWN)]):
    cnt = [len(csc[(csc.genotype == g) & (csc.found_by == tag)]) / n if len(csc) else 0 for g, n in zip(GENO_ORDER, n_iso)]
    ax.bar(x + (k - 0.5) * wd, cnt, wd, color=col, label=tag)
ax.set_xticks(x); ax.set_xticklabels(GENO_ORDER, rotation=45, ha="right", fontsize=8); ax.legend()
ax.set_ylabel("clusters per isolate"); ax.set_title(f"Clusters (p={SIG}) found by only one caller")
plt.tight_layout(); plt.savefig(f"{OUT}/caller_specific_mutations.png", dpi=180)

if len(eff):
    fig, ax = plt.subplots(figsize=(8.5, 8))
    e = eff.dropna(subset=["p_value_CLC", "p_value_own"])
    for scr, mk in [("cluster_analysis", "o"), ("cluster_coordination", "s"), ("cluster_association_chi2", "^")]:
        d = e[e.script == scr]
        ax.scatter(-np.log10(d.p_value_CLC.clip(lower=1e-12)), -np.log10(d.p_value_own.clip(lower=1e-12)), marker=mk, s=30,
                   c=[genotype_color(g) for g in d.base_genotype], edgecolor="k", lw=0.3, label=scr)
    lim = max(2, np.nanmax(-np.log10(e[["p_value_CLC", "p_value_own"]].clip(lower=1e-12).values)) * 1.05)
    ax.plot([0, lim], [0, lim], "k--", lw=1)
    ax.axvline(-np.log10(0.05), color="r", lw=0.8); ax.axhline(-np.log10(0.05), color="r", lw=0.8)
    ax.set_xlabel("-log10 p, CLC data (genotype vs ung1∆CLC)"); ax.set_ylabel("-log10 p, own-caller data (genotype vs ung1∆)")
    ax.set_title(f"Genotype effects: same conclusion in {int(e.same_conclusion.sum())}/{len(e)} tests\n(off-diagonal quadrants = conclusion depends on the caller)")
    ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(f"{OUT}/genotype_effects_by_caller.png", dpi=180)
print("figures saved to", OUT)
