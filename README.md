# Bioinformatics Internship: NicheNet Analysis

Title research project: In vitro interactomics of donor-derived SkMVECs and myotubes using NicheNet.

This repository contains the code used during my bioinformatics internship at KU Leuven Kukal. The project focuses on predicted cell-cell communication between skeletal muscle microvascular endothelial cells (SkMVECs), muscle cells (myotubes) and human umbilical vein endothelial cells (HUVECs) using RNA-seq data, differential expression analysis, NicheNet analysis, Gene Ontology and network visualisation.

The main biological question of the internship is:
**Which predicted cell-cell interactions occur between SkMVECs and myotubes, as identified using NicheNet analysis?**

## Project overview

The repository contains scripts for: 
1. Processing matched donor RNA-seq count files.
2. Performing differntial expression analysis between SkMVECs and myotubes.
3. Running NicheNet ligand activity analysis in the forward direction.
4. Running NicheNet ligand activity analysis in the reverse direction.
5. Creating ligand activity, ligand-target and ligand-receptor visualisations.
6. Creating a Circos plot of predicted ligand-target relationships.
7. Performing Gene Ontology enrichment on NicheNet results.
8. Comparing selected QuickGO gene lists between confluent SkMVECs and confluent HUVECs.

## Biological comparisons

**1. Matched donor SkMVECs and myotubes**
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

This analysis asks which ligands expressed by SkMVECs could explain genes upregulated in Myotubes.

**2. Confluent SkMVECs versus confluent HUVECs**
The analysis uses five donors (2D confluent):

| HUCECs sample | SkMVECs sample |
|---|---|
| GC122791 | GC122786 |
| GC122792 | GC122787 |
| GC122793 | GC122788 |
| GC122794 | GC124881 |
| GC122795 | GC122790 |

A separate DESeq2 analysis is performed to compare confluent SkMVECs with confluent HUVECs.

Positive log2 fold change = higher expression in SkMEVCs
Negative log2 fold change = higher expression in HUVECs

This analysis is also used to compare selected QuickGO-associated gene lists related to:
- Myoblast fusion GO:0007520
- Myoblast proliferation GO:0051450
- Sprouting angiogenesis GO:0002040

## Repository Layout

```text
data/
  counts/                 Raw .count files
  sample_metadata.tsv     Sample-to-donor mapping
scripts/
  NicheNet_assignment_SkMVECs_Myotubes.R
  NicheNet_reverse_Myotubes_sender_SkMVECs_receiver.R
  
  
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

The receiver gene set of interest is defined as genes significantly upregulated in the receiver cell type. NicheNet ranks sender-cell ligands by how well they match with the predicted target genes to the receiver-cell.
The results should be interpreted as candidate ligand-receptor-target hypotheses. 
Afterwards a Gene Ontology enrichment is performed on the predicted target genes. Using Cytoscape, networks are created consisting of ligand-receptro-target connections.
Also Gene Ontology terms are performed to understand the biological meaning of the interaction network.
