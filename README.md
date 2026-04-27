# Bioinformatics Internship: NicheNet Analysis

This repository contains the assignment workflow for a NicheNet ligand activity analysis of matched donor gene count files.

## Comparison

The analysis uses five matched donors:

| Donor | Myotubes sample | SkMVECs sample |
|---|---|---|
| 788 | GC122806 | GC122786 |
| 790 | GC122807 | GC122787 |
| 792 | GC122808 | GC122788 |
| 794 | GC122809 | GC124881 |
| 796 | GC122810 | GC122790 |

Default biological direction:

```text
Sender:   SkMVECs
Receiver: Myotubes
```

This asks which ligands expressed by SkMVECs could explain genes upregulated in Myotubes.

## Repository Layout

```text
data/
  counts/                 Raw .count files
  sample_metadata.tsv     Sample-to-donor mapping
scripts/
  NicheNet_assignment_Myotubes_SkMVECs.R
results/                  Output folder, ignored by Git
```

The NicheNet prior model files are downloaded automatically into `nichenet_networks/`. They are large and should not be committed.

## Install R Packages

Run once in R:

```r
install.packages(c("tidyverse", "BiocManager", "devtools", "cowplot", "ggpubr"))
BiocManager::install("edgeR")
devtools::install_github("saeyslab/nichenetr")
```

## Run

From the repository root:

```powershell
Rscript scripts/NicheNet_assignment_Myotubes_SkMVECs.R
```

Or open `scripts/NicheNet_assignment_Myotubes_SkMVECs.R` in RStudio and run the script.

## Outputs

The script writes results to:

```text
results/nichenet_assignment_results/
```

Important output files:

- `combined_counts_matrix.csv`
- `sample_metadata.csv`
- `DE_receiver_vs_sender_edgeR.csv`
- `geneset_of_interest_receiver_upregulated.csv`
- `potential_ligands.csv`
- `ligand_activities.csv`
- `active_ligand_target_links.csv`
- `active_ligand_receptor_links.csv`
- `04_nichenet_summary_figure.png`

## Method Summary

Raw count files were combined into one count matrix. Differential expression was performed with `edgeR` using a paired design:

```r
~ donor + cell_type
```

The receiver gene set of interest was defined as genes significantly upregulated in the receiver cell type. NicheNet then ranked sender-cell ligands by how well their predicted target genes matched this receiver-cell gene set.

The results should be interpreted as candidate ligand-receptor-target hypotheses, not direct proof of causal signaling.
