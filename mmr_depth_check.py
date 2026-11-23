"""Summarise read depth over the mismatch-repair genes and check it against the genotype labels.

Reads the samtools bedcov files written by mmr_depth_check.sge (OUT_DIR/per_sample/*.bedcov)
and, for every sample and gene, calls:

    DELETED    >= 90% of the ORF's 100 bp bins have (almost) no reads
    PARTIAL    less than that, but >= 2 consecutive no-read bins touch the ORF
               (part of the ORF is gone, e.g. only the 5' end / start codon)
    PRESENT    otherwise
    LOW_DEPTH  the chromosome's mean depth is too low to tell

A bin "has no reads" when its depth is below 10% of the chromosome's mean depth.
no_read_interval in the output gives where the missing stretch starts and ends.
The call is then compared with the genotype in mayo/settings/sample_genotypes.py
(a gene is expected deleted when the genotype contains e.g. "msh2∆").

Usage (run from the MayoMate folder; the .sge job does this for you):
    python mmr_depth_check.py outputs/mmr_depth_check/<run>
    python mmr_depth_check.py <run> --name-chars 8 --lost-fraction 0.1 --min-depth 10
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

# 1-based, inclusive (S288C R64-4-1, data/orf_coding.fasta). Keep in sync with GENES in mmr_depth_check.sge.
GENES = {
    "UNG1": (229734, 230813),
    "MSH2": (147382, 150276),
    "MSH3": (276764, 279820),
    "MSH6": (640109, 643837),
    "MLH1": (594886, 597195),
}
DELTA = "∆"  # the ∆ used in sample_genotypes.py


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_dir", type=Path, help="OUT_DIR of mmr_depth_check.sge (holds per_sample/)")
    p.add_argument("--name-chars", type=int, default=8,
                   help="sample name = first N characters of the BAM name, as in the genotype dictionary (default 8)")
    p.add_argument("--genotype-set", default="justin_set", help="which set in genotype_dict to use (default justin_set)")
    p.add_argument("--lost-fraction", type=float, default=0.10,
                   help="a bin has no reads when its depth < this x chromosome mean (default 0.10)")
    p.add_argument("--min-depth", type=float, default=10, help="chromosome mean depth below this -> LOW_DEPTH (default 10)")
    p.add_argument("--min-partial-bins", type=int, default=2,
                   help="consecutive no-read bins touching the ORF needed to call PARTIAL (default 2 = ~200 bp)")
    return p.parse_args()


def load_bedcov(per_sample: Path) -> pd.DataFrame:
    files = sorted(per_sample.glob("*.bedcov"))
    if not files:
        sys.exit(f"!!! no .bedcov files in {per_sample}")
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                        names=["contig", "start0", "end", "label", "sum"])
        d["bam"] = f.name[: -len(".bedcov")]
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    parts = d["label"].str.split("|", expand=True)
    d["gene"], d["kind"] = parts[0], parts[1]
    d["bin_start"] = pd.to_numeric(parts[2], errors="coerce") if parts.shape[1] > 2 else float("nan")
    d["depth"] = d["sum"] / (d["end"] - d["start0"])
    return d


def longest_lost_run(bins: pd.DataFrame, orf_start: int, orf_end: int, bin_size: int):
    """Run of consecutive no-read bins with the most bins touching the ORF.

    Returns (interval 'start-end' 1-based, number of its bins that touch the ORF), or ('', 0).
    The interval can reach outside the ORF (e.g. a deletion that also takes the promoter).
    """
    runs, run = [], []
    for r in bins.itertuples():
        if r.lost:
            run.append(r.bin_start)
        elif run:
            runs.append(run)
            run = []
    if run:
        runs.append(run)
    best, best_n = "", 0
    for run in runs:
        n_in_orf = sum(1 for b in run if b <= orf_end and b + bin_size - 1 >= orf_start)
        if n_in_orf > best_n:
            best, best_n = f"{int(run[0])}-{int(run[-1]) + bin_size - 1}", n_in_orf
    return best, best_n


def call_genes(d: pd.DataFrame, lost_fraction: float, min_depth: float,
               min_partial_bins: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, bin_rows = [], []
    for (bam, gene), g in d.groupby(["bam", "gene"], sort=True):
        s, e = GENES[gene]
        chrom_depth = g.loc[g.kind == "chrom", "depth"].iloc[0]
        orf_depth = g.loc[g.kind == "orf", "depth"].iloc[0]
        bins = g[g.kind == "bin"].sort_values("bin_start").copy()
        bin_size = int((bins.end - bins.start0).iloc[0])
        bins["rel_depth"] = bins.depth / chrom_depth if chrom_depth > 0 else float("nan")
        bins["lost"] = bins.rel_depth < lost_fraction
        in_orf = bins[(bins.bin_start >= s) & (bins.bin_start + bin_size - 1 <= e)]
        frac_lost = in_orf.lost.mean()
        interval, run_bins = longest_lost_run(bins, s, e, bin_size)
        if chrom_depth < min_depth:
            call = "LOW_DEPTH"
        elif frac_lost >= 0.9:
            call = "DELETED"
        elif run_bins >= min_partial_bins:
            call = "PARTIAL"
        else:
            call = "PRESENT"
        rows.append(dict(bam=bam, gene=gene, call=call, chrom_depth=round(chrom_depth, 1),
                         orf_depth=round(orf_depth, 1), orf_rel_depth=round(orf_depth / chrom_depth, 3) if chrom_depth else None,
                         frac_orf_bins_no_reads=round(frac_lost, 3), no_read_interval=interval))
        bins["bam"], bins["gene"] = bam, gene
        bin_rows.append(bins[["bam", "gene", "contig", "bin_start", "depth", "rel_depth"]])
    return pd.DataFrame(rows), pd.concat(bin_rows, ignore_index=True)


def mismatch(row) -> str:
    if row.genotype == "NOT_IN_DICT" or row.call == "LOW_DEPTH":
        return ""
    if row.expected_deleted and row.call == "PRESENT":
        return f"labelled {row.gene.lower()}{DELTA} but gene PRESENT"
    if row.expected_deleted and row.call == "PARTIAL":
        return f"labelled {row.gene.lower()}{DELTA}: only PARTIAL deletion"
    if not row.expected_deleted and row.call in ("DELETED", "PARTIAL"):
        return f"{row.call} but genotype has no {row.gene.lower()}{DELTA}"
    return ""


def plot_profiles(bins: pd.DataFrame, calls: pd.DataFrame, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out.mkdir(exist_ok=True)
    info = calls.set_index(["bam", "gene"])
    for gene, (s, e) in GENES.items():
        b = bins[bins.gene == gene]
        if b.empty:
            continue
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True, sharey=True)
        groups = [("labelled " + gene.lower() + DELTA, lambda r: r.expected_deleted and r.genotype != "NOT_IN_DICT"),
                  ("other genotypes", lambda r: (not r.expected_deleted) and r.genotype != "NOT_IN_DICT"),
                  ("not in genotype dictionary", lambda r: r.genotype == "NOT_IN_DICT")]
        for ax, (title, keep) in zip(axes, groups):
            n = n_flag = 0
            for bam, bb in b.groupby("bam"):
                r = info.loc[(bam, gene)]
                if not keep(r):
                    continue
                n += 1
                flagged = bool(r["mismatch"])
                n_flag += flagged
                ax.plot(bb.bin_start + 50, bb.rel_depth.clip(upper=2.5), lw=1.2 if flagged else 0.6,
                        color="tab:red" if flagged else "0.4", alpha=0.9 if flagged else 0.25, zorder=3 if flagged else 2)
            ax.axvspan(s, e, color="tab:blue", alpha=0.08, zorder=1)
            ax.axhline(1, color="0.7", lw=0.5, ls="--")
            ax.set_title(f"{title}: {n} samples" + (f", {n_flag} disagree with label (red)" if n_flag else ""), fontsize=10)
            ax.set_ylabel("depth / chromosome mean")
        axes[-1].set_xlabel(f"position (bp); shaded = {gene} ORF {s:,}-{e:,}")
        axes[0].set_ylim(0, 2.5)
        fig.suptitle(f"{gene}: read depth in 100 bp bins", fontsize=12)
        fig.tight_layout()
        fig.savefig(out / f"{gene}_depth_profile.png", dpi=150)
        plt.close(fig)


def main():
    args = parse_args()
    run = args.run_dir
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from mayo.settings.sample_genotypes import genotype_dict

    gdict = genotype_dict[args.genotype_set]
    d = load_bedcov(run / "per_sample")
    calls, bins = call_genes(d, args.lost_fraction, args.min_depth, args.min_partial_bins)

    calls.insert(1, "sample", calls.bam.str[: args.name_chars])
    calls.insert(2, "genotype", calls["sample"].map(gdict).fillna("NOT_IN_DICT"))
    calls["expected_deleted"] = [f"{g.lower()}{DELTA}" in gt for g, gt in zip(calls.gene, calls.genotype)]
    calls["mismatch"] = calls.apply(mismatch, axis=1)

    calls.to_csv(run / "mmr_calls_long.tsv", sep="\t", index=False)
    bins.to_csv(run / "mmr_bins.tsv", sep="\t", index=False)

    wide = calls.pivot(index=["bam", "sample", "genotype"], columns="gene", values="call")[list(GENES)]
    wide["mismatch"] = calls.groupby(["bam", "sample", "genotype"])["mismatch"].apply(lambda x: "; ".join(f for f in x if f))
    wide = wide.reset_index().sort_values(["genotype", "sample"])
    wide.to_csv(run / "mmr_calls_by_sample.tsv", sep="\t", index=False)
    wide[wide["mismatch"] != ""].to_csv(run / "mmr_mismatches.tsv", sep="\t", index=False)

    counts = calls.groupby(["gene", "genotype", "call"]).size().unstack("call", fill_value=0).reset_index()
    counts.to_csv(run / "mmr_genotype_vs_call.tsv", sep="\t", index=False)

    plot_profiles(bins, calls, run / "plots")

    # --- screen report ---
    pd.set_option("display.width", 200)
    pd.set_option("display.max_rows", 500)
    print(f"\n{wide.shape[0]} BAMs; {(wide.genotype == 'NOT_IN_DICT').sum()} not in genotype_dict['{args.genotype_set}'] "
          f"(sample name = first {args.name_chars} characters of the BAM name)")
    for gene in GENES:
        c = counts[(counts.gene == gene)]
        c = c[c.genotype.str.contains(gene.lower() + DELTA, regex=False) | (c.drop(columns=["gene", "genotype"]).get("DELETED", 0) > 0)
              | (c.drop(columns=["gene", "genotype"]).get("PARTIAL", 0) > 0)]
        print(f"\n== {gene}: genotypes labelled {gene.lower()}{DELTA}, or with any DELETED/PARTIAL call")
        print(c.drop(columns="gene").to_string(index=False) if len(c) else "   (none)")
    mism = wide[wide["mismatch"] != ""]
    print(f"\n== {len(mism)} samples disagree with their genotype label (full list: mmr_mismatches.tsv)")
    if len(mism):
        print(mism[["sample", "genotype", "mismatch"]].to_string(index=False))
    print(f"\nResults in {run}: mmr_calls_by_sample.tsv, mmr_mismatches.tsv, mmr_calls_long.tsv, plots/")


if __name__ == "__main__":
    main()
