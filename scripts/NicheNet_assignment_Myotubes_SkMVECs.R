# NicheNet ligand activity analysis for the internship assignment
# Comparison: matched donors, Myotubes vs SkMVECs 2D confluent
#
# This script reads the provided raw gene count files, performs paired
# differential expression with edgeR, and runs NicheNet on the genes
# upregulated in the receiver cell type.

# Load required packages
#install.packages("BiocManager")
#BiocManager::install("edgeR")

library(tidyverse)  # data manipulation and plotting
library(edgeR)      # differential expression analysis
library(nichenetr)  # NicheNet ligand activity analysis
library(cowplot)    # combining plots
library(ggpubr)     # publication-style plots

### Settings/parameters -------------------------------------------------------

# Main biological question:
# Which ligands expressed by SkMVECs could explain the changing gene activity observed in Myotubes?

# To run the opposite direction, change these two values:
sender_celltype <- "SkMVECs"
receiver_celltype <- "Myotubes"

# NicheNet analysis parameters
# Organism-specific ligand-receptor and ligand-target network.
organism <- "human"
# NicheNet ranks ligands by predicted activity, keep top 30 ligands for downstream analysis.
top_n_ligands <- 30
# Counts per million: filtering threshold for gene expression.
min_cpm <- 1
# This threshold removes genes with very low or unreliable expression.
min_samples_expressed <- 2
# False discovery rate: adjusted p-value threshold
fdr_cutoff <- 0.05
# log2 fold change
logfc_cutoff <- 1

setwd("C:/Users/royde/OneDrive/Documenten/BIT11 Internship/repo_upload")

# Set project folders
project_dir <- getwd()

count_dir <- file.path(project_dir, "data", "counts")
output_dir <- file.path(project_dir, "results", "nichenet_assignment_results")
network_dir <- file.path(project_dir, "nichenet_networks")

# Create output folders if they do not exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)

# Create sample metadata
sample_info <- tribble(
  ~sample_id, ~donor, ~cell_type, ~file_name,
  "GC122806", "788", "Myotubes", "GC122806_KG.count",
  "GC122807", "790", "Myotubes", "GC122807_KG.count",
  "GC122808", "792", "Myotubes", "GC122808_KG.count",
  "GC122809", "794", "Myotubes", "GC122809_KG.count",
  "GC122810", "796", "Myotubes", "GC122810_KG.count",
  "GC122786", "788", "SkMVECs", "GC122786_KG.count",
  "GC122787", "790", "SkMVECs", "GC122787_KG.count",
  "GC122788", "792", "SkMVECs", "GC122788_KG.count",
  "GC124881", "794", "SkMVECs", "GC124881_KG.count",
  "GC122790", "796", "SkMVECs", "GC122790_KG.count"
)

# Add full file paths and convert variables
sample_info <- sample_info %>%
  mutate(
    file_path = file.path(count_dir, file_name),
    donor = factor(donor),
    cell_type = factor(cell_type)
  )

# Set SkMVECs as the reference cell type
sample_info$cell_type <- relevel(sample_info$cell_type, ref = "SkMVECs")

# Check if all count files are present
missing_files <- sample_info$file_path[!file.exists(sample_info$file_path)]

if (length(missing_files) > 0) {
  stop("Some count files are missing: ", paste(missing_files, collapse = ", "))
}


### 1. Read count files -------------------------------------------------------

# Function to read one count file
read_count_file <- function(file_path, sample_id) {
  
  counts_one_sample <- read_tsv(
    file_path,
    col_names = c("gene", sample_id),
    show_col_types = FALSE
  )
# if the gene appears more than once, this adds the counts together  
  counts_one_sample <- counts_one_sample %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    summarise(
      count = sum(.data[[sample_id]], na.rm = TRUE),
      .groups = "drop"
    )
  
  colnames(counts_one_sample) <- c("gene", sample_id)
  
  return(counts_one_sample)
}

# Read all count files into a list
count_tables <- list()

for (i in 1:nrow(sample_info)) {
  file_path <- sample_info$file_path[i]
  sample_id <- sample_info$sample_id[i]
  
  count_tables[[i]] <- read_count_file(file_path, sample_id)
}

# Combine all count tables into one table
counts_tbl <- count_tables[[1]]

for (i in 2:length(count_tables)) {
  counts_tbl <- full_join(counts_tbl, count_tables[[i]], by = "gene")
}

# Replace missing values with zero
counts_tbl[is.na(counts_tbl)] <- 0

# Convert the table into a count matrix
counts_df <- as.data.frame(counts_tbl)

rownames(counts_df) <- counts_df$gene
counts_df$gene <- NULL

counts_df <- counts_df[, sample_info$sample_id]

counts <- as.matrix(counts_df)
storage.mode(counts) <- "integer"

head(rownames(counts))

# Save the combined count matrix and metadata
write_csv(counts_tbl, file.path(output_dir, "combined_counts_matrix.csv"))
write_csv(sample_info, file.path(output_dir, "sample_metadata.csv"))

### 2. Paired differential expression ----------------------------------------

# Differential expression analysis using edgeR.
# Which genes are more highly expressed in the receiver cell type (Myotubes), 
# compared with the sender cell type (sKMVECs), while accounting for matched donors?

# Create an edgeR object from the count matrix
y <- DGEList(
  counts = counts,
  samples = sample_info
)

# Remove genes with very low expression
keep_genes <- filterByExpr(
  y,
  group = sample_info$cell_type
)

y <- y[keep_genes, , keep.lib.sizes = FALSE]

# Normalize the count data
y <- calcNormFactors(y)

# Create the design matrix
# This compares cell types while correcting for donor differences
design <- model.matrix(
  ~ donor + cell_type,
  data = sample_info
)

# Show the design matrix columns
colnames(design)

# In this analysis, SkMVECs is the reference group
# So this coefficient tests Myotubes compared with SkMVECs
receiver_coef <- "cell_typeMyotubes"

# Check if the coefficient exists
if (!(receiver_coef %in% colnames(design))) {
  stop("The coefficient cell_typeMyotubes was not found in the design matrix.")
}

# Estimate variation between samples
y <- estimateDisp(y, design)

# Fit the statistical model: glmQLFit is the step where edgeR learns the expression pattern and variability of each gene 
# so we can test which genes are upregulated or downregulated between cell types.
fit <- glmQLFit(y, design)

# Test for genes upregulated in Myotubes compared with SkMVECs
qlf <- glmQLFTest(
  fit,
  coef = receiver_coef
)

# Get the differential expression results
de_results <- topTags(qlf, n = Inf)$table

# Add gene names as a normal column
de_results$gene <- rownames(de_results)

# Convert to tibble and sort by FDR
de_results <- as_tibble(de_results)
de_results <- de_results %>%
  select(gene, everything()) %>%
  arrange(FDR)

# Save the differential expression results
write_csv(
  de_results,
  file.path(output_dir, "DE_receiver_vs_sender_edgeR.csv")
)

# Select significantly upregulated genes in Myotubes
geneset_oi_raw <- de_results %>%
  filter(FDR <= 0.05, logFC >= 1) %>%
  pull(gene)

# Check if enough genes were selected for NicheNet
if (length(geneset_oi_raw) < 20) {
  warning("Less than 20 genes were selected. Consider using FDR <= 0.10 or logFC >= 0.5.")
}


### 3. Load NicheNet prior model ---------------------------------------------

# Give R more time to download the large NicheNet files
options(timeout = 1200)
# NicheNet network files for human data
lr_network_url <- "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
ligand_target_url <- "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
weighted_networks_url <- "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"

# File names where the downloaded files will be saved
lr_network_file <- file.path(network_dir, "lr_network_human_21122021.rds")
ligand_target_file <- file.path(network_dir, "ligand_target_matrix_nsga2r_final.rds")
weighted_networks_file <- file.path(network_dir, "weighted_networks_nsga2r_final.rds")

# Download files only if they are not already present
if (!file.exists(lr_network_file)) {
  download.file(lr_network_url, lr_network_file, mode = "wb")
}

if (!file.exists(ligand_target_file)) {
  download.file(ligand_target_url, ligand_target_file, mode = "wb")
}

if (!file.exists(weighted_networks_file)) {
  download.file(weighted_networks_url, weighted_networks_file, mode = "wb")
}

# Load the NicheNet networks into R
lr_network <- readRDS(lr_network_file)
ligand_target_matrix <- readRDS(ligand_target_file)
weighted_networks <- readRDS(weighted_networks_file)

# Keep only unique ligand-receptor combinations
lr_network <- lr_network %>%
  distinct(from, to)

lr_network <- lr_network %>%
  mutate(
    from = toupper(from),
    to = toupper(to)
  )

rownames(ligand_target_matrix) <- toupper(rownames(ligand_target_matrix))
colnames(ligand_target_matrix) <- toupper(colnames(ligand_target_matrix))



### 4. Define sender ligands and receiver receptors --------------------------

# Calculate CPM values for all genes and samples
plain_cpm <- cpm(y)

# Get sender samples
sender_samples <- sample_info %>%
  filter(cell_type == sender_celltype) %>%
  pull(sample_id)

# Get receiver samples
receiver_samples <- sample_info %>%
  filter(cell_type == receiver_celltype) %>%
  pull(sample_id)

# Find genes expressed in sender samples
expressed_genes_sender <- rownames(plain_cpm)[
  rowSums(plain_cpm[, sender_samples] >= min_cpm) >= min_samples_expressed
]

# Find genes expressed in receiver samples
expressed_genes_receiver <- rownames(plain_cpm)[
  rowSums(plain_cpm[, receiver_samples] >= min_cpm) >= min_samples_expressed
]

message("Expressed sender genes: ", length(expressed_genes_sender))
message("Expressed receiver genes: ", length(expressed_genes_receiver))

head(expressed_genes_sender)
head(expressed_genes_receiver)
