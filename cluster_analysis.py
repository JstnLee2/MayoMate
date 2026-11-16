# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
Cluster analysis: compares mutation-cluster statistics between genotypes.

1. Create combined cluster table (aggregates the per-sample cluster data into one table).
2. Conduct cluster data analysis: every genotype is compared against REFERENCE_GENOTYPE
   (t-test, Mann-Whitney U, KS, Levene) and violin plots are drawn.

Settings: edit the block below, or override from the command line, e.g.
    python cluster_analysis.py --reference ung1∆SNPless --genotypes ung1∆SNPless ung1∆SNPlessmsh2∆

Genotype names, colors and plot order come from mayo/settings/genotype_colors.py.

Outputs:
    outputs/logs/cluster_analysis.log                       full text log (all test results)
    outputs/clusters/stats_vs_<reference>_<settings>.csv    one row per genotype x metric with p-values
    figures/clusters/...                                    violin plots
'''
import argparse
import logging
import sys
import os
from mayo.settings import genotype_dict as genotype_dict
from mayo.settings.genotype_colors import ANALYSIS_GENOTYPES, REFERENCE_GENOTYPE as DEFAULT_REFERENCE, genotype_palette, genotype_order
from mayo.clusters import create_combined_cluster_table
from mayo.clusters import conduct_cluster_data_analysis

# ============================== SETTINGS ==============================
REFERENCE_GENOTYPE = DEFAULT_REFERENCE   # genotype all others are compared against (default "ung1∆"); e.g. "ung1∆SNPless"
GENOTYPES = ANALYSIS_GENOTYPES           # genotypes to include; all by default. Or a list, e.g.
# GENOTYPES = ['ung1∆', 'ung1∆msh2∆', 'ung1∆msh3∆', 'ung1∆msh6∆', 'ung1∆SNPless', 'ung1∆SNPlessmsh2∆']
FIGSIZE = None                           # None = scale with the number of genotypes, or e.g. (5, 5)

SNP_NUM = 1
ASSOCIATION = None # None, True or False
ASSOCIATION_FEATURE = "both" # "GC", "CO", "both"
CLUSTER_CALLER = "JT" # JT or PMACD (deprecated)
SIG = 0.0001 #JT; 0.01, 0.001, 0.0001, 0.00001
PLOIDY = "all" # all, haploid, aneuploid
min_cluster_mutations = 3
LOG_FILE = "outputs/logs/cluster_analysis.log"
# =====================================================================


def parse_args():
    p = argparse.ArgumentParser(description="Compare mutation-cluster statistics between genotypes.")
    p.add_argument("--reference", default=REFERENCE_GENOTYPE, help=f"reference genotype (default: {REFERENCE_GENOTYPE})")
    p.add_argument("--genotypes", nargs="+", default=None, help="genotypes to include (default: GENOTYPES in this file)")
    p.add_argument("--sig", type=float, default=SIG)
    p.add_argument("--snp-num", type=int, default=SNP_NUM)
    p.add_argument("--ploidy", default=PLOIDY, choices=["all", "haploid", "aneuploid"])
    p.add_argument("--log", default=LOG_FILE)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    reference = args.reference
    genotypes = list(args.genotypes or GENOTYPES)
    if reference not in genotypes:
        genotypes = [reference] + genotypes   # the reference must be in the data to compare against it
    order = genotype_order(genotypes)
    palette = genotype_palette(order)
    figsize = FIGSIZE or (max(5, 1.1 * len(order)), 5)

    os.makedirs(os.path.dirname(args.log), exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:.%(funcName)s: %(message)s',
        handlers=[logging.FileHandler(args.log, mode='w', encoding='utf-8'), logging.StreamHandler()],
        force=True)
    logger = logging.getLogger(__name__)
    logger.info(f"Reference genotype: {reference}")
    logger.info(f"Genotypes: {order}")
    logger.info(f"SIG={args.sig}, SNP_NUM={args.snp_num}, ploidy={args.ploidy}, association={ASSOCIATION}, feature={ASSOCIATION_FEATURE}")

    ### Part 1: Create combined cluster table
    # make sure the right cluster association is read in if you are using GC or CO
    create_combined_cluster_table(
        cluster_caller = CLUSTER_CALLER, #JT or PMACD
        snp_num=args.snp_num,
        sig=args.sig,
        files = ['1_master_para_bed_sorted_anz5.txt'],
        genotype_dict_master = genotype_dict,
        association = ASSOCIATION, # None, True or False
        association_feature = ASSOCIATION_FEATURE,
        min_cluster_mutations = min_cluster_mutations
        )

    ### Part 2: Conduct cluster data analysis (printed test results go to the same log file)
    with open(args.log, 'a', encoding='utf-8') as f:
        original_stdout = sys.stdout
        sys.stdout = f
        try:
            conduct_cluster_data_analysis(
                cluster_caller = CLUSTER_CALLER, #JT or PMACD
                snp_num = args.snp_num,
                sig = args.sig,
                association = ASSOCIATION, # None, True or False
                association_feature = ASSOCIATION_FEATURE,
                ploidy = args.ploidy, #all, haploid, aneuploid
                genotype1 = reference,   # reference: every other genotype is compared against it
                genotype2 = '',          # optional: genotype for the first detailed comparison ('' = first other genotype)
                save_plt = True, # save the plot as a png file
                save_plotly = False, # save the plot as an interactive html file,
                genotypes_include = order,
                figsize = figsize,
                order = order,
                palette = palette,
                )
        finally:
            sys.stdout = original_stdout
    from mayo.log_archive import keep_log_copy
    logger.info(f"Done. Full results in {args.log}")
    logger.info(f"Copy of this run's log kept at {keep_log_copy(args.log)}")
