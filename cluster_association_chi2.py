"""
Chi-squared test: are mutation clusters associated with recombination events more (or less)
often in each genotype than in the reference genotype?

Replaces notebook cells "change the sample.genotype to the strain you want..." (cluster counts)
and "from scipy.stats import chi2_contingency" (manual 2x2 test). Instead of typing counts by
hand, it counts them for every genotype and tests each one against REFERENCE_GENOTYPE.

Input: the pickle the main notebook saves after the association step
       (outputs/partial_pickles/samples_parental_SNPs_JT_<SIG>_<SNP_NUM>snp.pkl), so it can be
       run any time later without rerunning the notebook.

For each genotype, a 2x2 table is tested against the reference:
                     not associated   associated
    reference              a               b
    genotype               c               d

Usage (from the MayoMate folder):
    python cluster_association_chi2.py
    python cluster_association_chi2.py --reference ung1∆SNPless
    python cluster_association_chi2.py --feature GC      # GC-associated only (or CO)
    python cluster_association_chi2.py --include-het-chr # also count clusters on heterozygous/aneuploid chromosomes

Clusters on a sample's heterozygous or aneuploid chromosomes are left out by default (as in the notebook's
association table), because switches there are not reliable. The number left out is in the
excluded_het_aneuploid_chr column.

Outputs:
    outputs/cluster_association_chi2_<feature>_vs_<reference>.csv    one row per genotype
    outputs/logs/cluster_association_chi2.log                       full log incl. 2x2 tables
"""
import argparse
import logging
import os
import pickle

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

from mayo.settings.genotype_colors import ANALYSIS_GENOTYPES, REFERENCE_GENOTYPE as DEFAULT_REFERENCE, genotype_order

# ============================== SETTINGS ==============================
REFERENCE_GENOTYPE = DEFAULT_REFERENCE   # genotype all others are compared against (default "ung1∆")
GENOTYPES = ANALYSIS_GENOTYPES           # genotypes to test; all by default
FEATURE = "both"                         # "both" (any recombination event), "GC" or "CO"
CLUSTER_CALLER, SIG, SNP_NUM = "JT", 0.0001, 1
REMOVE_RTG = False                       # True = drop samples in rtg_list (the notebook cell did not)
INCLUDE_HET_CHR = False                  # True = also count clusters on heterozygous/aneuploid chromosomes
LOG_FILE = "outputs/logs/cluster_association_chi2.log"
# =====================================================================

FEATURE_COLUMNS = {"both": "Clustered_Proximal", "GC": "Clustered_GC_Proximal", "CO": "Clustered_CO_Proximal"}


def count_associated_clusters(samples, feature_column: str, include_het_chr: bool = False) -> pd.DataFrame:
    """One row per sample: number of clusters not associated / associated with recombination.

    By default, clusters on the sample's heterozygous or aneuploid chromosomes are left out, as in the
    notebook's association table (overlapClustersWithSwitches). Switches on those chromosomes are not
    reliable, so association there is noise. include_het_chr=True counts them too (the old behavior).
    """
    rows = []
    for sample in samples:
        df = getattr(sample, "df_SNPs", None)
        if df is None or "Cluster_ID" not in df or feature_column not in df:
            continue
        per_cluster = df.dropna(subset=["Cluster_ID"]).groupby("Cluster_ID").agg(
            chromosome=("Chromosome", "first"), feature=(feature_column, "first"))
        excluded = (set(getattr(sample, "heterozygous_chromosomes", None) or [])
                    | set(getattr(sample, "aneuploid_chromosomes", None) or []))
        on_excluded = per_cluster["chromosome"].isin(excluded)
        kept = per_cluster if include_het_chr else per_cluster[~on_excluded]
        rows.append({
            "sample": sample.name,
            "genotype": sample.genotype,
            "not_associated": int((kept["feature"] == False).sum()),  # noqa: E712 (values can be numpy bools)
            "associated": int((kept["feature"] == True).sum()),      # noqa: E712
            "excluded_het_aneuploid_chr": 0 if include_het_chr else int(on_excluded.sum()),
        })
    per_sample = pd.DataFrame(rows)
    return per_sample


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", default=REFERENCE_GENOTYPE)
    p.add_argument("--genotypes", nargs="+", default=None)
    p.add_argument("--feature", default=FEATURE, choices=list(FEATURE_COLUMNS))
    p.add_argument("--pickle", default=f"outputs/partial_pickles/samples_parental_SNPs_{CLUSTER_CALLER}_{SIG}_{SNP_NUM}snp.pkl")
    p.add_argument("--remove-rtg", action="store_true", default=REMOVE_RTG)
    p.add_argument("--include-het-chr", action="store_true", default=INCLUDE_HET_CHR,
                   help="also count clusters on heterozygous/aneuploid chromosomes (old behavior)")
    p.add_argument("--log", default=LOG_FILE)
    args = p.parse_args()

    os.makedirs(os.path.dirname(args.log), exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s:%(levelname)s: %(message)s",
                        handlers=[logging.FileHandler(args.log, mode="w", encoding="utf-8"), logging.StreamHandler()],
                        force=True)
    log = logging.getLogger(__name__)

    reference = args.reference
    genotypes = genotype_order(set(args.genotypes or GENOTYPES) | {reference})
    feature_column = FEATURE_COLUMNS[args.feature]
    log.info(f"Reference: {reference} | feature: {args.feature} ({feature_column}) | pickle: {args.pickle} | "
             f"clusters on heterozygous/aneuploid chromosomes: {'included' if args.include_het_chr else 'excluded'}")

    with open(args.pickle, "rb") as fh:
        samples_SNVs, _switches_list, _chromosomes = pickle.load(fh)

    per_sample = count_associated_clusters(samples_SNVs, feature_column, include_het_chr=args.include_het_chr)
    if args.remove_rtg:
        from mayo.settings import rtg_list
        per_sample = per_sample[~per_sample["sample"].isin(rtg_list)]
    per_sample = per_sample[per_sample["genotype"].isin(genotypes)]
    counts = per_sample.groupby("genotype")[["not_associated", "associated", "excluded_het_aneuploid_chr"]].sum()
    counts["samples"] = per_sample.groupby("genotype").size()
    counts = counts.reindex([g for g in genotypes if g in counts.index])

    if reference not in counts.index:
        raise SystemExit(f"Reference genotype '{reference}' has no clusters in {args.pickle}. Genotypes found: {list(counts.index)}")

    ref = counts.loc[reference]
    results = []
    for genotype, row in counts.iterrows():
        table = np.array([[ref["not_associated"], ref["associated"]],
                          [row["not_associated"], row["associated"]]])
        total = row["not_associated"] + row["associated"]
        res = {
            "genotype": genotype,
            "reference": reference,
            "feature": args.feature,
            "samples": int(row["samples"]),
            "not_associated": int(row["not_associated"]),
            "associated": int(row["associated"]),
            "percent_associated": round(100 * row["associated"] / total, 2) if total else np.nan,
            "excluded_het_aneuploid_chr": int(row["excluded_het_aneuploid_chr"]),
        }
        if genotype != reference:
            try:
                chi2, p_val, dof, expected = chi2_contingency(table)
                res.update(chi2=chi2, p_value=p_val, dof=dof, min_expected=float(expected.min()))
            except ValueError as e:  # e.g. a row or column of zeros
                res.update(chi2=np.nan, p_value=np.nan, dof=np.nan, min_expected=np.nan, note=str(e))
            # Fisher's exact test as a check when expected counts are small (< 5)
            res["p_value_fisher"] = fisher_exact(table)[1]
            log.info(f"\n{genotype} vs {reference}\n{pd.DataFrame(table, index=[reference, genotype], columns=['not associated', 'associated'])}\n"
                     f"chi2 = {res.get('chi2')}, p = {res.get('p_value')}, Fisher p = {res['p_value_fisher']}")
        results.append(res)

    out = pd.DataFrame(results)
    out_path = f"outputs/cluster_association_chi2_{args.feature}_vs_{reference}.csv"
    out.to_csv(out_path, index=False, encoding="utf-8")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        log.info("\n" + out.to_string(index=False))
    log.info(f"Saved: {out_path}")
    from mayo.log_archive import keep_log_copy
    log.info(f"Copy of this run's log kept at {keep_log_copy(args.log)}")


if __name__ == "__main__":
    main()
