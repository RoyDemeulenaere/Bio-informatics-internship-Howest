# Bioinformatics Internship: NicheNet Analysis

This repository contains the assignment workflow for a NicheNet ligand activity analysis of matched donor gene count files.
To find out how myoblasts and Skeletal Microvascular Endothelial cells (SkMVECs) interact with each other.


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


## Method Summary

Raw count files were combined into one count matrix. Differential expression was performed with `edgeR` using a paired design:

```r
~ donor + cell_type
```

The receiver gene set of interest is defined as genes significantly upregulated in the receiver cell type. NicheNet ranks sender-cell ligands by how well they match with the predicted target genes to the receiver-cell.
The results should be interpreted as candidate ligand-receptor-target hypotheses. 
Afterwards a Gene Ontology enrichment is performed on the predicted target genes. Using Cytoscape, networks are created consisting of ligand-receptro-target connections.
Also Gene Ontology terms are performed to understand the biological meaning of the interaction network.
