# CLC vs own caller: 4 ung1∆ isolates (JL1123–JL1126)

A separate copy of MayoMate that runs only 8 samples: the same 4 isolates, each called two ways.

| Group (pseudo-genotype) | Files | Variant caller |
|---|---|---|
| `ung1∆CLC` (reference) | `JL01123_` … `JL01126_` | CLC Genomics |
| `ung1∆` | `JL1123_S` … `JL1126_S` | own caller |

How it differs from the main project:
- `mayo/settings/sample_genotypes.py`: `justin_set` is replaced at the end of the file by these 8 samples.
- `mayo/settings/genotype_colors.py`: only the two groups above, with `REFERENCE_GENOTYPE = "ung1∆CLC"`.
- `main_playground_clean.ipynb`: copied from `comparisons/JL531_S1_vs_S2`. Its genotype lists use
  `GENOTYPES`, and the simulation subset uses `genotype_dict` instead of `haploid_samples_ung1d`.
- The parent files are renamed to short names (`Mal_JK69_parent.csv`, `Mal_JK70_parent.csv`) to stay
  under Windows' 260-character path limit.

## Re-run (from this folder)

```
..\..\myenv\Scripts\python.exe run_notebook.py --continue-on-error
..\..\myenv\Scripts\python.exe cluster_analysis.py
..\..\myenv\Scripts\python.exe cluster_coordination.py
..\..\myenv\Scripts\python.exe cluster_association_chi2.py
..\..\myenv\Scripts\python.exe compare_raw_calls.py
..\..\myenv\Scripts\python.exe compare_pipeline_outputs.py
```

The notebook gives 37 ok and 1 error. The error is cell 13 (optional unreliable-SNP finder,
`KeyError: 'Chromosome'`), which also fails in the JL531 comparison.

## Outputs

- `comparison_outputs/`: the CLC vs own-caller comparison.
  - `pipeline_comparison_per_isolate.csv`: every metric, per isolate.
  - `pipeline_comparison_paired_tests.csv`: paired tests across the 4 isolates, plus tests on the
    pooled data.
  - `raw_calls_summary.csv`, `raw_calls_paired_tests.csv`: the raw call files, before any filtering.
  - `own_only_clusters.csv`: clusters that only the own caller reports (artifact check).
  - Figures: `pipeline_comparison_overview.png`, `pipeline_comparison_clusters.png`,
    `raw_calls_concordance.png`.
- `outputs/`, `figures/`, `data/new_clusters/`: the normal MayoMate outputs for these 8 samples.
- `outputs/clusters/stats_vs_ung1∆CLC_*.csv`, `outputs/cluster_length_meta_test_df_JT_vs_ung1∆CLC.csv`,
  `outputs/cluster_association_chi2_both_vs_ung1∆CLC.csv`: the three scripts' statistical tests
  (own caller vs CLC).

With 4 pairs, the smallest possible two-sided Wilcoxon p is 0.125. Paired t-tests and the
"all 4 isolates move the same way" column carry the per-isolate comparison; the pooled tests use all
mutations and clusters.
