# Running MayoMate end-to-end (locally or on an HPC)

This guide explains how to run the whole analysis without clicking through the notebook,
how to submit it to an HPC (qsub / Grid Engine), and where to change genotypes, colors and the
reference genotype.

## What's new

| File | What it does |
|---|---|
| `run_notebook.py` | Runs `main_playground_clean.ipynb` top to bottom with no user input. Saves an executed copy of the notebook (all outputs) after every cell and writes a report of each cell's status, runtime and memory. |
| `run_mayomate.sge` | HPC job script for `qsub` (Open Grid Scheduler / Grid Engine): sets up Python, runs the notebook with `run_notebook.py`, then runs the three analysis scripts below. |
| `mayo/sample_cache.py` | Per-sample cache: the first run processes every sample, later runs only process samples added since. |
| `mayo/settings/genotype_colors.py` | **One place** for genotype names, colors, plotting order and the reference genotype. Used by the notebook and every script. |
| `cluster_association_chi2.py` | Chi-squared test of recombination-associated clusters, every genotype vs. the reference (replaces the old notebook cells that needed counts typed in by hand). |
| `cluster_analysis.py` | Updated: shared colors/genotypes, `REFERENCE_GENOTYPE` setting, log file, CSV of all p-values. |
| `cluster_coordination.py` | Updated: shared colors/genotypes, `REFERENCE_GENOTYPE` setting (was hard-coded to ung1∆). |
| `report_rtg_risk()` (in `mayo`) | Lists samples at risk of being Return-To-Growth (APSC < 30) with their coverage plot. Runs in notebook cell 7. |

## 1. Run everything on your own computer

From the MayoMate folder, with `myenv` set up (see README):

```
myenv\Scripts\activate                      (Windows)   or   source myenv/bin/activate   (Mac/Linux)
python run_notebook.py --continue-on-error
```

- `--continue-on-error` keeps going if a cell fails, so one run shows every problem.
  Without it, the run stops at the first error.
- Progress is printed one line per cell. Plots are saved to `figures/` as usual but not shown on screen.
- Each run gets its own folder, `outputs/run_logs/<date_time>/` (see [Logs](#logs-every-run-is-kept)):
  - `run_report.md`: table of every cell: ok / error / skipped, time, peak memory, error message.
  - `main_playground_clean_executed.ipynb`: the notebook with all outputs. Open it in VS Code to see any cell's printout or error.
  - `logs/`: copies of the logs written during the run.

Other options:

```
python run_notebook.py --list                 # show cell numbers
python run_notebook.py --skip 5 9             # skip cells (e.g. the long ones, if their pickles already exist)
python run_notebook.py --only 2 20-29         # run only some cells
python run_notebook.py --rerun-all            # reprocess every sample (ignore the sample cache)
python run_notebook.py --rerun JL102_R1       # reprocess only these samples
```

**Only new samples are processed.** The slow per-sample steps (import, filters, ploidy, recombination
maps, JT clusters) are saved per sample in `outputs/sample_cache/` the first time a sample is run.
Later runs load those and only process samples that were added (or whose input file changed).
Steps that compare samples (common-mutation filtering, statistics, figures) always use all samples.
Use `--rerun-all` after changing how samples are processed in the `mayo` code. On the HPC:
`qsub -v RUN_NOTEBOOK_ARGS="--rerun-all" run_mayomate.sge`.

Cell numbers count from 1 at the top of the notebook, **including text cells**, the same way
you count them in VS Code.

## 2. Run it on an HPC (qsub / Open Grid Scheduler)

Job script: **`run_mayomate.sge`** (for Open Grid Scheduler / Grid Engine, submitted with `qsub`).

1. **Copy the project to the cluster**, including `data/` (about 9 GB of sequencing calls) and
   `S288C_reference_sequencewoUra3wUra329wNAT.fa`. Do not copy `myenv/` (Windows-only); the job builds
   its own Linux environment (`venv_hpc/`) on the first run. For example, from your computer:
   ```
   scp -r MayoMate  <username>@<cluster-login-address>:~/
   ```
   (Or use Globus / WinSCP if your HPC recommends it for large transfers.)
2. **Check the lines marked `<-- EDIT` in `run_mayomate.sge`**. On the cluster's login node:
   - `#$ -pe smp 4`: the parallel environment name. `qconf -spl` lists them (often `smp`, `shm`, `threaded` or `openmp`).
   - `module load python/3.10`: the Python 3.10 module name. `module avail python` lists them.
     If your cluster has no modules but has `python3.10` installed, the line is skipped automatically.
   - `#$ -q`: only if your cluster requires a queue name (`qconf -sql` lists them).
   - Memory (`h_vmem`, per CPU) and time (`h_rt`): defaults are 4 x 8G = 32G and 8 h, from the measured numbers below.
3. **Submit** from the MayoMate folder on the cluster:
   ```
   qsub run_mayomate.sge
   ```
4. **Watch progress**:
   ```
   qstat -u $USER            # qw = waiting, r = running, Eqw = error
   tail -f mayomate.o<jobid>
   ```
5. **Results**: figures and tables in `figures/` and `outputs/` as usual. Everything about the run
   (report, executed notebook, all logs, the job's screen output) is in
   `outputs/run_logs/<date_time>_job<jobid>/`.

Options when submitting (no need to edit the file):

```
qsub -v REFERENCE=ung1∆SNPless run_mayomate.sge           # different reference genotype
qsub -v RUN_SCRIPTS=0 run_mayomate.sge                    # notebook only
qsub -v RUN_NOTEBOOK_ARGS="--skip 5 9" run_mayomate.sge   # skip cells
qsub -l h_vmem=12G run_mayomate.sge                      # more memory (per CPU) without editing the file
```

The first job takes a few extra minutes to install packages; later jobs reuse `venv_hpc/`. If the
first job fails during setup, the half-built `venv_hpc/` is deleted, so just fix the problem and resubmit.

## 3. Genotypes, colors and the reference genotype

Everything is in **`mayo/settings/genotype_colors.py`**:

```python
GENOTYPE_COLORS = {
    "ung1∆":       "tab:blue",
    "ung1∆fun30∆": "tab:olive",
    ...
}
EXCLUDE_FROM_ANALYSIS = ["ung1∆ 37C_A3A", "ung1∆ 37C_EV"]
REFERENCE_GENOTYPE = "ung1∆"
```

- **Add a genotype**: add its samples to `mayo/settings/sample_genotypes.py`, then add one line to
  `GENOTYPE_COLORS`. It's automatically included in every "all genotypes" list (notebook and scripts)
  and gets the same color in every plot. Its position in the dictionary is its position in plots.
- **Change a color**: edit it in `GENOTYPE_COLORS`. Any matplotlib color name or hex code works.
- **Leave a genotype out of the main analysis** but keep its color: add it to `EXCLUDE_FROM_ANALYSIS`.
- **Reference genotype**: change `REFERENCE_GENOTYPE` there (default for everything), or per script:
  - at the top of `cluster_analysis.py`, `cluster_coordination.py`, `cluster_association_chi2.py`
    (`REFERENCE_GENOTYPE = ...`), or
  - on the command line: `python cluster_analysis.py --reference ung1∆SNPless`
  The reference is added automatically if it isn't in the genotype list.

In the notebook, `GENOTYPES` (set in cell 2) is the "all genotypes" list, and `GENOTYPE_QUERY` is the
same list as the `"ung1∆+ung1∆fun30∆+..."` string used by `get_samples_with_genotype()`.

Note: the colors are the ones from `cluster_analysis.py`, as chosen. A few genotypes share a color
(msh3∆/rnh201∆ cyan; msh6∆/37C orange; fun30∆/pol3_01 olive; mph1∆/pol3_5DV brown;
mlh2∆/SNPlessmsh2∆/rnh1∆rnh201∆ green). `ung1∆rnh1∆rnh201∆` had no color in `cluster_analysis.py`
and was given `tab:green`. Change any of them in the one file.

## 4. Return-To-Growth (RTG) check

Notebook cell 7 runs `report_rtg_risk()` after the parental SNP data are loaded. It:

1. Takes each sample's APSC (average per-SNP coverage, the number in the coverage plot names
   `figures/coverage/<sample>_(APSC-<value>)_SNP_coverage_map.png`).
2. Flags samples with APSC < 30 (`apsc_threshold=30` to change).
3. Writes `outputs/rtg_risk_report.csv` (lowest APSC first) with genotype, APSC, the coverage plot path,
   whether the sample is already in `rtg_list`, and an `action` column (`CHECK coverage plot`, `already in rtg_list`, `ok`).
4. Prints the samples to check as ready-to-paste lines, e.g. `"JL587_S3",  # APSC 10.84, ung1∆SNPlessmsh2∆`.

Open the coverage plot of each flagged sample, then paste its line into `mayo/settings/chromosome_lengths.py`:

- **RTG** → `rtg_list`: the sample is removed from the analysis (notebook cell 14).
- **Not RTG** → `rtg_checked_ok` (just below `rtg_list`): the sample stays in the analysis and is
  **not flagged again** in later runs. Adding the date to the comment helps, e.g.
  `"JL587_S3",  # APSC 10.84, ung1∆SNPlessmsh2∆ - checked 2026-11-11, OK`.

Each run then only lists new low-APSC samples. The report's `action` column shows
`CHECK coverage plot`, `already in rtg_list`, `checked OK (rtg_checked_ok)` or `ok`. To re-check a sample,
delete its line from `rtg_checked_ok`.

## 5. Chi-squared test of cluster association (old notebook cells 42-43)

```
python cluster_association_chi2.py                      # all genotypes vs ung1∆
python cluster_association_chi2.py --reference ung1∆SNPless
python cluster_association_chi2.py --feature GC         # GC-associated clusters only (or CO)
python cluster_association_chi2.py --include-het-chr    # also count clusters on heterozygous/aneuploid chromosomes
```

It reads `outputs/partial_pickles/samples_parental_SNPs_JT_0.0001_1snp.pkl` (saved by the notebook),
counts clusters associated / not associated with recombination for every genotype, and tests each
genotype against the reference with a chi-squared test (and Fisher's exact test as a check).
Clusters on a sample's heterozygous or aneuploid chromosomes are left out by default, as in the notebook's
association table, because switches there are not reliable. The `excluded_het_aneuploid_chr` column gives
how many were left out per genotype. Before 2026-11-20 they were counted.
Results: `outputs/cluster_association_chi2_<feature>_vs_<reference>.csv`; details incl. every 2x2
table: `outputs/logs/cluster_association_chi2.log` (latest; every run is also kept, see [Logs](#logs-every-run-is-kept)).

## 6. cluster_analysis.py and cluster_coordination.py

Both now take their genotypes and colors from `genotype_colors.py` and have a `REFERENCE_GENOTYPE`
setting at the top (or `--reference` / `--genotypes` on the command line).

| Script | Log | Main results |
|---|---|---|
| `cluster_analysis.py` | `outputs/logs/cluster_analysis.log` | `outputs/clusters/stats_vs_<reference>_<settings>.csv` (p-values of every genotype vs. reference, for each metric), `figures/clusters/` |
| `cluster_coordination.py` | `outputs/logs/cluster_coordination.log` | `outputs/cluster_length_meta_test_df_JT_vs_<reference>.csv` (Mann-Whitney p-values per cluster type), `figures/cluster_types_by_pval/`, `figures/cluster_length_violin_by_pval/` |

## 6b. Check the mismatch-repair deletions in the BAMs (`mmr_depth_check.sge`)

Confirms each mutant carries its deletion by measuring read depth (`samtools bedcov`) over
MSH2, MSH3, MSH6 and MLH1, with UNG1 as a control (deleted in every strain). It needs the
BAM files; CLC read mappings must first be exported as BAM (CLC: *Export → BAM*).

```
qsub -v BAM_DIR=/path/to/bams mmr_depth_check.sge
qsub -v BAM_DIR=/path/to/bams,MIN_MAPQ=30 mmr_depth_check.sge
```

- Check the `<-- EDIT` lines first (`-pe`, `module load samtools`), as for `run_mayomate.sge`.
- Sample names are the first 8 characters of the BAM name (`NAME_CHARS=8`), matched against
  `sample_genotypes.py`. The screen report says how many BAMs were not found there. If that's most of
  them, the BAM names don't start like the CSV names; set `NAME_CHARS` or rename.
- Any contig naming works (`ref|NC_001147|`, `NC_001147.6`, `chrXV`, `XV`, `chr15`).
- Each gene gets **DELETED** (≥90% of the ORF has no reads), **PARTIAL** (≥200 bp of the ORF has no
  reads, e.g. only the 5′ end), **PRESENT** or **LOW_DEPTH**. "No reads" = below 10% of the chromosome's mean depth.
- Results go to `outputs/mmr_depth_check/<date_time>_job<id>/`:
  `mmr_mismatches.tsv` (samples whose depth disagrees with their label), `mmr_calls_by_sample.tsv`,
  `mmr_calls_long.tsv` (depths and the no-read interval, ±100 bp), and `plots/<GENE>_depth_profile.png`.
- The summary can be re-run with other thresholds without re-reading the BAMs:
  `python mmr_depth_check.py outputs/mmr_depth_check/<run> --min-partial-bins 3`.

## 7. What changed in the notebook

A backup of the notebook before these changes is in `outputs/backups/main_playground_clean_before_cleanup_20261110.ipynb`.

- Removed: sector code (`specific_genotypes_2/3/4`, sector filters and tables) from the SNV import
  cells (it only affected sector samples, none of which are in the data); the "I don't think I need
  this (JL)" sector counts cell; the rtg True/False check; the duplicate plotting cell (old 27); the
  `df_filter` cell (old 35); the cluster-count and chi-squared cells (old 42-43, now
  `cluster_association_chi2.py`); the notebook copy of `drawSNPMap` (its changes were moved into the
  package version); the null-hypothesis plot that read a missing file.
- Commented out (code kept): the OptionalJL indel analysis (cell 3).
- Moved: the null-hypothesis plot now comes after the simulation cell that creates its input.
- Added: RTG risk check (cell 7); note pointing to the chi-squared script.
- All "all genotypes" lists now use `GENOTYPES` / `GENOTYPE_QUERY`; Windows-only `\` paths fixed so
  the notebook runs on Linux.

## Logs: every run is kept

The logs in `outputs/logs/` always hold the **latest** run (they are overwritten each time).
Copies of every run are kept in `outputs/run_logs/`, which is never overwritten:

```
outputs/run_logs/
├── 20261111_082323/                       # a local run of run_notebook.py
│   ├── run_report.md / run_report.json    # per-cell status, time, memory, errors
│   ├── main_playground_clean_executed.ipynb
│   └── logs/
│       └── main_mayo_0.0001_1snp_JT.log   # the notebook's own log
├── 20261112_090000_job123456/             # an HPC job (run_mayomate.sge)
│   ├── run_report.md / .json
│   ├── main_playground_clean_executed.ipynb
│   ├── logs/                              # notebook log + cluster_analysis.log,
│   │                                      #   cluster_coordination.log, cluster_association_chi2.log
│   └── sge_123456.out                     # the job's screen output
└── scripts/                               # dated copies each time a script is run by itself
    ├── cluster_analysis_20261111_101500.log
    ├── cluster_coordination_20261111_102000.log
    └── cluster_association_chi2_20261111_082240.log
```

- On the HPC the job's screen output is also written to `mayomate.o<jobid>` in the MayoMate folder
  while it runs (for `tail -f`); a copy goes into the run folder when the job ends, even if it fails part-way.
- Running a cell by hand in VS Code still overwrites `outputs/logs/main_mayo_...log`; only runs
  through `run_notebook.py` / the HPC job are archived.
- The notebook's log is saved as UTF-16. To read it on the cluster:
  `iconv -f UTF-16 -t UTF-8 <file> | less` (VS Code opens it directly).
- Old run folders are never deleted automatically; remove ones you no longer need.

## Test results

Full test on 2026-11-11 (Windows laptop, `myenv`, a fresh copy of the project, all 620 samples):

| Step | Result | Time |
|---|---|---|
| Notebook, all 38 code cells | 37 ok, 1 error (cell 33, fixed since, see below) | 2 h 26 min |
| `cluster_analysis.py` | ok, wrote `stats_vs_ung1∆_...csv` | ~1 min |
| `cluster_coordination.py` | ok, wrote `cluster_length_meta_test_df_JT_vs_ung1∆.csv` | ~2 min |
| `cluster_association_chi2.py` | ok, wrote `cluster_association_chi2_both_vs_ung1∆.csv` | ~2 min |

- **Peak memory: 12.4 GB** (the notebook's Python process, measured every 2 s).
- Long cells: 5 (load SNPs, 15 min), 9 (maps, 9 min), 23 (38 min), 34 (cluster maps, 8 min),
  36 (6 min), 37 (27 min), 42 (combined maps, 26 min). Everything else takes under 3 min.
- Cell 33 (expression heatmap from Table S1) failed because it reads `df_features_summary.csv`,
  which no notebook cell creates. This was the same before the cleanup. The cell now prints
  "Skipped: ..." when the file is missing. To use it, create the file with
  `findTranscriptAndDirection()` in `mayo/prototypes.py`.
- Chi-squared note: `ung1∆SNPless` has 0 associated clusters (no parental SNPs means no
  recombination events can be detected), so its p-value against ung1∆ is not meaningful.

The job scripts ask for 8 hours (about 3x the laptop runtime, since cluster cores can be slower) and
32G (`run_mayomate.sge`: 4 CPUs x `h_vmem=8G`; Grid Engine limits *virtual* memory, which for Python
is well above real use). Raise them if the job is killed.

## Troubleshooting

- **A cell fails**: open the executed notebook in `outputs/run_logs/<run>/` and scroll to that cell to see the full error.
- **Job killed / "DeadKernelError" / "h_vmem exceeded"**: usually out of memory; resubmit with
  `qsub -l h_vmem=12G run_mayomate.sge`.
- **Job hits the time limit**: raise `h_rt`, or skip finished long cells with
  `RUN_NOTEBOOK_ARGS="--skip ..."` (the "reload from pickle" cells right after them load their saved results).
- **Job stuck in `Eqw`**: `qstat -j <jobid> | grep error` shows why (often a wrong `-pe` or `-q` name);
  fix the line, then `qdel <jobid>` and resubmit.
- **"no Python found" / `module load python/3.10` fails**: run `module avail python` on the cluster and use the listed name.
- **`/bin/bash^M: bad interpreter`**: the script got Windows line endings while being copied. Fix with
  `sed -i 's/\r$//' run_mayomate.sge` (the repository's `.gitattributes` prevents this for git checkouts).
