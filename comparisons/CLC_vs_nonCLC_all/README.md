# CLC vs own caller: all 156 paired isolates

Each isolate was called twice: by CLC Genomics (`<genotype>CLC`) and by the own caller (`<genotype>`). Both
versions went through the full MayoMate notebook and the three cluster scripts together. Genotypes stay
separate throughout. How to re-run is in `RESUME.md`. The pairs are listed in `pairs.json`.

| Genotype | Isolates | Genotype | Isolates |
|---|---:|---|---:|
| ung1∆SNPless | 30 | ung1∆ 37C_Selected | 14 |
| ung1∆SNPlessmsh6∆ | 27 | ung1∆SNPlessmsh3∆ | 13 |
| ung1∆msh6∆_NonSelected | 20 | ung1∆ | 11 |
| ung1∆msh6∆ | 18 | ung1∆ 37C_NonSelected | 10 |
| ung1∆_zymolyaseheat | 9 | ung1∆_NonSelected | 4 |

## Results

**1. Coverage.** The own caller reports about 15% lower coverage at every site. This holds in 100% of isolates
and every genotype: average per-SNP coverage is 70 vs 83, and coverage at mutations is 72 vs 84. It is a
systematic offset, probably different read or duplicate filtering, not a genotype effect.

**2. Mutations agree closely.** Final mutation counts are 282 (CLC) vs 285 (own) per isolate, and 41,136
mutations are shared. Mutations found by only one caller have low coverage (median 56 own-only, 45 CLC-only,
74 shared), so they are mostly borderline calls. The A3A fraction is 96.5% vs 97.0%.

**3. Clusters agree.** JT clusters are 11.9 vs 11.9 per isolate, and clustered mutations and ssDNA length are not
significantly different overall. Cluster length and size distributions are the same. There are 155 CLC-only and
160 own-only clusters, mostly in ung1∆msh6∆ and ung1∆msh6∆_NonSelected.

**4. Recombination agrees, once heterozygous chromosomes are excluded.** The first run showed far fewer crossovers
with the own caller (25.9 vs 42.6). That came from a bug: `get_all_rec_events_df` left heterozygous and aneuploid
chromosomes out for gene conversions but not for crossovers. Fixed on 2026-11-20.

| Isolates | n | CO, CLC | CO, own | CO identical |
|---|---:|---:|---:|---:|
| no heterozygous chromosome | 89 | 17.1 | 17.2 | 97% of isolates |
| with heterozygous chromosome(s), haploid chromosomes only | 67 | 14.0 | 13.2 | 75% |
| all | 156 | 15.8 | 15.5 | 87% (Wilcoxon q = 0.23) |

No genotype differs between callers (all q ≥ 0.13). Gene conversions are 8.1 vs 8.0. Raw switches
(`all_switches_1snp.tsv`, no filter) still differ on heterozygous chromosomes, where CLC flips between parental
haplotypes (96 vs 26 per isolate), but they no longer feed into GC/CO counts. Crossovers per isolate now
count haploid chromosomes only, so genotypes with many heterozygous isolates have lower counts. Compare
them per haploid chromosome or per Mb, not as raw counts.

**5. Genotype conclusions hold with either caller.** Comparing each genotype with ung1∆ gives the same
significant or not-significant result in 178 of 183 tests. All 5 that differ are borderline (p between 0.02 and
0.04): average mutations per cluster (3) and cluster_coordination cluster lengths (2).
- **Recombination-associated clusters (chi-square)** now agree for every genotype. `cluster_association_chi2.py`
  used to count clusters on heterozygous and aneuploid chromosomes, unlike the notebook's association table. It now
  leaves them out (fixed 2026-11-20; `--include-het-chr` gives the old behavior). Percent associated, CLC vs own:
  ung1∆ 46 / 46, msh6∆ 56 / 51, msh6∆_NonSelected 39 / 50, 37C_Selected 58 / 59, 37C_NonSelected 55 / 55,
  zymolyaseheat 53 / 55, SNPlessmsh3∆ 30 / 31. None differ from ung1∆ except SNPless and SNPlessmsh6∆ (about 1%
  and 0%, p < 1e-10 with both callers).
- **Fewer clusters left.** Most clusters in some genotypes are on heterozygous chromosomes: msh6∆_NonSelected
  keeps 18 of 215, 37C_NonSelected 55 of 379, SNPlessmsh3∆ 56 of 233, msh6∆ 54 of 217. ung1∆_NonSelected has none
  left (all 88 are on heterozygous chromosomes), so it can't be tested.
- Between callers within a genotype (the pair runs), only 2 of 126 tests have p < 0.05, both cluster_coordination
  cluster lengths, with no correction for multiple testing.

**6. Samples to check.**
- **JL1022 (ung1∆msh6∆):** the callers share almost no mutations (452 CLC-only, 655 own-only). It has very low
  coverage (APSC 13) and only 2,659 to 5,772 parental SNPs. It is a low-quality sample, or the two files may not
  be the same library.
- **JL1009 (ung1∆msh6∆):** 211 own-only mutations, and 7 vs 24 clusters.

## Outputs (`comparison_outputs/`)

| File | Contents |
|---|---|
| `pipeline_per_isolate.csv` | every metric per isolate, both callers, and `has_het_chr` |
| `paired_tests.csv` | own vs CLC paired Wilcoxon (BH q) for all isolates, with and without heterozygous chromosomes, and per genotype; Kruskal-Wallis between genotypes |
| `pooled_tests.csv` | spectrum, A3A, cluster length and size, GC:CO, coverage of caller-specific mutations |
| `genotype_effects_by_caller.csv/.png` | each genotype vs ung1∆, with CLC data and with own-caller data |
| `scripts_CLC_vs_own_by_genotype.csv` | the three scripts run on gCLC vs g |
| `caller_specific_clusters.csv` | clusters found by only one caller |
| `raw_calls_*` | raw call files, before any MayoMate filtering |
| `paired_differences_by_genotype.png`, `CLC_vs_own_scatter.png`, `caller_specific_mutations.png` | figures |

## Known failures (not caused by the callers)

- Notebook cell 13 (unreliable-SNP finder) fails with `KeyError: 'Chromosome'`, as in the main project.
- Cell 24: the figure filename joins all 20 genotypes and goes over Windows' 260-character path limit.
- Cells 25–26: `KeyError: 'Reference'` when one genotype has no mutations near recombination events.
- `cluster_coordination.py` fails for pair_SNPless, pair_SNPlessmsh6∆ and pair_ung1∆_NonSelected because those
  genotypes have almost no switches or gene conversions. It runs normally on CLC_set and own_set.
