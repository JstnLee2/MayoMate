# Resuming the all-pairs CLC vs own-caller run

**Status: complete (2026-11-20, 19:30).** The run was resumed and finished. Results are in `README.md`. The steps
below are only needed to re-run it.

## What is saved

All 312 samples (156 isolates × CLC / own caller) are cached, so these steps will not be redone:

| Saved | Where |
|---|---|
| Parental SNP import, filters, ploidy | `outputs/sample_cache/parental_SNPs/` (312) |
| Recombination maps | `outputs/sample_cache/recombination_1snp/` (312) |
| SNV import and filters | `outputs/sample_cache/SNVs/` (312) |
| JT clusters (5 p-values) | `outputs/sample_cache/jt_clusters/` (1560) |
| Parental reference | `outputs/sample_cache/parental_origin/` |
| Notebook checkpoints from cells 5–10 | `outputs/partial_pickles/` |

Already finished, and they don't need the notebook:
- `comparison_outputs/raw_calls_*`: the raw-call comparison (`compare_raw_calls.py`).

## Resume

From this folder (`comparisons/CLC_vs_nonCLC_all`), with the laptop plugged in and set not to sleep:

```
..\..\myenv\Scripts\python.exe run_notebook.py --continue-on-error > nb_run_stdout.txt 2>&1
..\..\myenv\Scripts\python.exe run_scripts.py
..\..\myenv\Scripts\python.exe compare_pipeline_outputs.py > comparison_stdout.txt 2>&1
```

1. **Notebook.** Re-run the whole notebook, not only the remaining cells, so that every variable is set again.
   Cells 5, 9 and 10 load the cached samples instead of reprocessing them, which saves about 25 minutes. Expect
   roughly 1–1.5 hours, mostly the long cells after cell 11 (23, 37, 42).
2. **`run_scripts.py`.** Runs `cluster_analysis.py`, `cluster_coordination.py` and `cluster_association_chi2.py`
   12 times: once per genotype (CLC vs own caller), on the CLC set alone, and on the own-caller set alone.
   Each run's outputs are kept in `script_runs/<run>/`. Expect about 15–30 minutes.
3. **`compare_pipeline_outputs.py`.** Runs the full comparison and writes its tables and figures to
   `comparison_outputs/`. Expect about 5 minutes.

Do not use `--rerun-all`, because it throws away the cache. If a run is interrupted again, run the same commands:
samples that were already processed are always loaded from the cache.

## What is in this workspace

- **Samples:** `mayo/settings/sample_genotypes.py`. At the end, `justin_set` is replaced by the 156 paired
  isolates (also listed in `pairs.json`). CLC files are `<genotype>CLC`, own-caller files are `<genotype>`.
- **Genotype colors:** `mayo/settings/genotype_colors.py` has the 20 genotypes (10 × 2), with the main project's
  colors.
- **Sample files:** `mayo/settings/config.py` reads them from the main project's `data/` folder (`../../data/...`).
  They are not copied here.
- **Common-mutation filter:** in the notebook, `max_duplicates` is 20 instead of 10. Each isolate is in the data
  twice, so 20 samples is still 10 isolates.
