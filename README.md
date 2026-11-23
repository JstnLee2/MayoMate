# MayoMate

`MayoMate` is a Python toolkit for analyzing and interpreting meiotic mutagenesis data in yeast. It finds and visualizes recombination events, mutations, and mutation clusters in Whole Genome Sequencing (WGS) data from *Saccharomyces cerevisiae* meiotic outcomes.

`MayoMate` processes Illumina read alignments and variant calls from CLC Genomics Workbench 20 to identify high-confidence de novo SNVs and mutation clusters. It can build a Parental Reference Genome so that parental SNPs can be called against a single modified reference. Its classes store per-sample information (SNVs, mutation clusters, parental SNPs, recombination events, sample names, and genotypes), which supports event classification, genomic feature correlation, hypermutable ssDNA estimation, mutation rate normalization, statistical testing, and custom plots, summaries, and tables.

### Features
- Identify and visualize recombination events (gene conversions and crossovers) from a panel of parental SNPs.
- Identify, visualize, and classify mutation clusters.
- Associate mutation clusters with recombination events.
- Infer the length of hypermutable single-stranded DNA (ssDNA).
- Estimate mutation rates around promoters and other genomic features.
- Compare statistics between mutant genotypes.
- Measure parental SNP density and test whether mutation clusters are enriched in SNP-dense regions.
- Plot pre- and post-meiotic mutation frequencies.

### Pipeline overview
![Meiotic Pipeline_new_bg](https://github.com/user-attachments/assets/62bee6ce-a6d3-4793-94d5-b0d6fe9d154b)

## Getting started

### Requirements
MayoMate was developed and tested with Python 3.10 on a standard desktop computer with 16 GB of RAM. Its dependencies are pinned in `requirements.txt`:

| Package | Version |
|---|---|
| betterbeeswarm | 0.2.0 |
| bio-aid | 0.3.1 |
| matplotlib | 3.5.2 |
| natsort | 8.3.1 |
| numpy | 1.22.3 |
| pandas | 1.4.2 |
| plotly | 5.11.0 |
| scipy | 1.8.1 |
| seaborn | 0.13.2 |

You also need Jupyter to run the main notebook.

### Installation
Clone the repository and install the dependencies into a virtual environment. This should take a few minutes.

```bash
git clone https://github.com/malkovalab/MayoMate.git
cd MayoMate
python -m venv myenv

# Windows
myenv\Scripts\activate
# macOS / Linux
source myenv/bin/activate

pip install -r requirements.txt
pip install jupyter
```

Run every script and the notebook from the repository root. All paths in the settings files are relative to it.

### Configuration
Before running anything, check the settings in [`mayo/settings/`](mayo/settings/):

| File | What it controls |
|---|---|
| `config.py` | `config`: paths to the parental, parental-SNP, and SNV call folders and the reference FASTA, plus analysis parameters (`SNP_NUM`, `SIG`, `CLUST_INTER_MUT_MAX`, `CLUSTER_TYPE_ANALYSIS`). `config_cluster`: inputs and outputs for `cluster_analysis.py`. |
| `sample_genotypes.py` | Maps each sample name to its genotype and defines sample sets (for example, haploid samples). Add new samples here. |
| `chromosome_lengths.py` | S288C chromosome names and lengths, sample sets, and the RTG lists `rtg_list` / `rtg_checked_ok` (see [RTG samples](#rtg-samples-and-heterozygous--aneuploid-chromosomes)). |
| `problematic_positions.py` | Genomic positions to exclude from analysis. |
| `non_essential_genes.py` | List of non-essential genes used in feature analyses. |

The defaults point to the example data under `data/sequencing_calls/CombinedSet/` and the reference `S288C_reference_sequencewoUra3wUra329wNAT.fa`.

## Usage

### Main workflow
Start with [`main_playground_clean.ipynb`](main_playground_clean.ipynb). It walks through MayoMate's functions step by step: loading variant calls, filtering SNVs, calling parental SNPs and recombination events, calling and classifying mutation clusters, and making figures.

The repository includes a small subset of the full study's data. Use it to confirm your environment is set up correctly. The notebook and its log files show example outputs and compute times for this subset. The notebook writes tables to `outputs/` and figures to `figures/`, and the scripts below take those outputs as input.

### Running everything without user input (locally or on an HPC)
`python run_notebook.py --continue-on-error` runs the whole notebook top to bottom and writes a per-cell report; `qsub run_mayomate.sge` does the same on an HPC (Open Grid Scheduler / Grid Engine) and then runs the analysis scripts. The first run processes every sample; later runs only process samples added since (`--rerun-all` to redo everything, see `mayo/sample_cache.py`). Genotype names, colors and the reference genotype are set in one file, `mayo/settings/genotype_colors.py`. See [HPC_RUN_GUIDE.md](HPC_RUN_GUIDE.md).

### Analysis scripts
These run after the notebook (`run_mayomate.sge` runs them automatically) and read its outputs. Run them from the repository root with `python <script>.py`; settings are at the top of each file or given on the command line (e.g. `--reference`).

| Script | Purpose |
|---|---|
| `cluster_analysis.py` | Analyzes called mutation clusters and compares cluster statistics between mutants. Uses `config_cluster`. |
| `cluster_coordination.py` | Classifies mutation clusters by their relation to recombination switches (for example "single-switch", "non-switching", "multiple-switch") and plots their ssDNA contributions. |
| `cluster_association_chi2.py` | Chi-squared test of recombination-associated clusters, each genotype vs. the reference genotype. Reads the notebook's saved pickle. |

### RTG samples and heterozygous / aneuploid chromosomes
Return-to-growth (RTG) samples are diploids that entered meiosis and went back to mitotic growth. Their parental-SNP switches are not meiotic gene conversions or crossovers, so they are left out of the analysis.

**Finding them.** Notebook cell 7 (`report_rtg_risk`) lists samples with an average per-SNP coverage (APSC) below 30 in `outputs/rtg_risk_report.csv`, with the path to each sample's coverage plot. Check each plot, then add the sample to one of two lists in `mayo/settings/chromosome_lengths.py`:
- `rtg_list`: RTG, removed from the analysis.
- `rtg_checked_ok`: not RTG, kept, and not flagged again.

**Where `rtg_list` applies.**
- Notebook: RTG samples are removed before switch calling (cell 9) and from the SNV samples before they are paired with their switches (cell 11). The switch hotspot plot (cell 27) also leaves them out.
- Scripts: `cluster_analysis.py`, `cluster_association_chi2.py` and `cluster_coordination.py` each drop `rtg_list` samples themselves. After adding a sample to `rtg_list`, rerunning the scripts updates their results. Switch and GC/CO results need the notebook rerun from cell 9.

**Heterozygous and aneuploid chromosomes.** Cell 5 finds them for each sample from coverage and SNP frequencies. Switches are called on every chromosome, and each analysis then leaves out the affected chromosomes:

| Analysis (notebook cell) | Heterozygous / aneuploid chromosomes |
|---|---|
| Recombination events file and GC/CO location plot (16, 28) | Excluded |
| Cluster–switch association (36) and chi-squared script | Excluded |
| Switch hotspot plot (27) | Excluded |
| Mutations around recombination events (25, 26) | Heterozygous excluded; aneuploid kept |
| GC/CO distribution (24) | Whole sample skipped if it has either |
| Haplotype block lengths (37, 38) | Whole sample skipped if heterozygous; aneuploid kept |

## Repository layout

```
MayoMate/
├── mayo/                      # Core package
│   ├── settings/              # Configuration: paths, parameters, sample genotypes
│   ├── clusters/              # Mutation cluster tables and summaries
│   ├── nbinom.py              # Negative binomial cluster calling (JT method)
│   ├── find_GC_CO.py          # Gene conversion / crossover detection
│   ├── *Table*.py             # Sample, SNP, parental and switch table classes
│   └── sample_cache.py        # Per-sample cache: later runs only process new samples
├── data/                      # Input data: CLC variant calls, genomic features, cluster tables
├── outputs/                   # Generated tables, pickles, logs, sample cache
├── figures/                   # Generated figures
├── main_playground_clean.ipynb
├── run_notebook.py            # Runs the notebook without user input
├── run_mayomate.sge           # HPC job (qsub): notebook + analysis scripts
├── cluster_*.py               # Analysis scripts (see above)
├── S288C_reference_sequencewoUra3wUra329wNAT.fa   # Modified S288C reference
└── requirements.txt
```

## Data
This repository includes only a small subset of the data for testing. The full dataset comes from WGS of *S. cerevisiae* meiotic outcomes. Reads were processed with CLC Genomics Workbench 20, and the resulting variant calls are MayoMate's input.

Raw reads are available from NCBI BioProject / SRA under accession **[PRJNA1225307](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1225307)**.

## Contributing
Contributions are welcome. To contribute:

1. Fork the repository.
2. Create a branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -am 'Add new feature'`).
4. Push the branch (`git push origin feature-branch`).
5. Open a pull request.

## Authors
MayoMate was developed by Jerzy Twarowski in the Malkova Lab (University of Iowa). Later analysis scripts were added by lab contributors.

## License
MayoMate is released under the [MIT License](LICENSE).
