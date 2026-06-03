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
# Which ligands expressed by Myotubes could explain the changing gene activity observed in SkMVECs?

# To run the opposite direction, change these two values:
sender_celltype <- "Myotubes"
receiver_celltype <- "SkMVECs"

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
output_dir <- file.path(project_dir, "results", "nichenet_assignment_results_reverse")
network_dir <- file.path(project_dir, "nichenet_networks")

# Create output folders if they do not exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)

# Create sample metadata    (tribble: create small table mannually)
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

# Set the sender cell type as the reference cell type
sample_info$cell_type <- relevel(sample_info$cell_type, ref = sender_celltype)

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
# if the gene appears more than once, this adds the counts together (isoform!)
# isoform-information is lost, NicheNet works with gene symbols, not transcript isoforms
# NicheNet expects one value per gene. 
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
# Which genes are more highly expressed in the receiver cell type,
# compared with the sender cell type, while accounting for matched donors?

# Create an edgeR object from the count matrix
# DGEList puts the raw count data and sample metadata into a format that edgeR can use for further analysis.
y <- DGEList(
  counts = counts,
  samples = sample_info
)

# Remove genes with very low expression (how is this filterByExpr filtered?!: A gene is kept if it has enough expression in enough samples from at least one group. Filter keeps genes with sufficient counts per million.)
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

# The sender cell type is the reference group.
# This coefficient tests the receiver cell type compared with the sender cell type.
receiver_coef <- paste0("cell_type", receiver_celltype)

# Check if the coefficient exists
if (!(receiver_coef %in% colnames(design))) {
  stop("The receiver coefficient was not found in the design matrix: ", receiver_coef)
}

# Estimate variation between samples
y <- estimateDisp(y, design)

# Fit the statistical model: glmQLFit is the step where edgeR learns the expression pattern and variability of each gene 
# so we can test which genes are upregulated or downregulated between cell types.
fit <- glmQLFit(y, design)

# Test for genes upregulated in the receiver cell type compared with the sender cell type
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
  dplyr::select(gene, everything()) %>%
  dplyr::arrange(FDR)

# Save the differential expression results
write_csv(
  de_results,
  file.path(output_dir, "DE_receiver_vs_sender_edgeR.csv")
)

# Select significantly upregulated genes in the receiver cell type
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
# nsga2r: Non-dominated sorting genetic algorithm 2 --> genetic algorithm to find the best solutions when you have multiple goals to optimize at once.
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

# Get all ligands and receptors from the NicheNet ligand-receptor network
all_ligands <- unique(lr_network$from)
all_receptors <- unique(lr_network$to)

# Keep only ligands expressed in the sender cells
expressed_ligands <- intersect(all_ligands, expressed_genes_sender)

# Keep only receptors expressed in the receiver cells
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

# Keep ligand-receptor pairs where both genes are expressed
possible_lr_pairs <- lr_network %>%
  filter(from %in% expressed_ligands) %>%
  filter(to %in% expressed_receptors)

# Get the possible ligands from these ligand-receptor pairs
potential_ligands <- unique(possible_lr_pairs$from)

# Keep only ligands that are present in the ligand-target matrix
potential_ligands <- intersect(
  potential_ligands,
  colnames(ligand_target_matrix)
)

# Keep genes of interest that are present in the NicheNet model
geneset_oi <- intersect(
  geneset_oi_raw,
  rownames(ligand_target_matrix)
)

# Define background genes:
# expressed receiver genes that are present in the NicheNet model
background_expressed_genes <- intersect(
  expressed_genes_receiver,
  rownames(ligand_target_matrix)
)

# Save the NicheNet input files
write_csv(
  tibble(gene = geneset_oi),
  file.path(output_dir, "geneset_of_interest_receiver_upregulated.csv")
)

write_csv(
  tibble(gene = background_expressed_genes),
  file.path(output_dir, "background_receiver_expressed_genes.csv")
)

write_csv(
  tibble(ligand = potential_ligands),
  file.path(output_dir, "potential_ligands.csv")
)

# Print summary information
message("Sender cell type: ", sender_celltype)
message("Receiver cell type: ", receiver_celltype)
message("Number of expressed sender ligands: ", length(expressed_ligands))
message("Number of expressed receiver receptors: ", length(expressed_receptors))
message("Number of potential ligands: ", length(potential_ligands))
message("Number of genes of interest: ", length(geneset_oi))
message("Number of background genes: ", length(background_expressed_genes))

# Stop if too few genes of interest are available for NicheNet
if (length(geneset_oi) < 10) {
  stop("Too few genes of interest overlap with the NicheNet model.")
}


### 5. NicheNet ligand activity analysis -------------------------------------
# Which sender ligands best explain the genes that are upregulated in the receiver cells?
# Ranks Myotube ligands based on how well they could explain the SkMVEC gene expression changes.

# Run NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Sort ligands from best to worst
ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected)) #Sorts ligands on corrected Area Under the Precision-Recall curve (AUPR).
# A higher AUPR means the ligand is better at predicting your receiver genes of interest.

# Add a rank column
ligand_activities <- ligand_activities %>%
  mutate(rank = row_number())

# Save all ligand activity results
write_csv(
  ligand_activities,
  file.path(output_dir, "ligand_activities.csv")
)
# top_n_ligands = 30
# Select the top ligands
best_upstream_ligands <- ligand_activities %>%
  slice_head(n = top_n_ligands)

# Keep only the ligand names
best_upstream_ligands <- best_upstream_ligands$test_ligand

# Save the top ligand names
write_lines(
  best_upstream_ligands,
  file.path(output_dir, "best_upstream_ligands.txt")
)


### 6. Infer ligand-target and ligand-receptor links -------------------------
#1. Which target genes are predicted for these ligands?
#2. Which receptors could these ligands bind to?
#3. Can we visualize these links as heatmaps and/or Circos plot?

# Get predicted target genes for each top ligand: Ligand-Target Links
ligand_target_list <- lapply(
  best_upstream_ligands,
  get_weighted_ligand_target_links,
  geneset = geneset_oi,
  ligand_target_matrix = ligand_target_matrix,
  n = 200
)

# Combine all ligand-target results into one table
active_ligand_target_links_df <- bind_rows(ligand_target_list)

# Save ligand-target links
write_csv(
  active_ligand_target_links_df,
  file.path(output_dir, "active_ligand_target_links.csv")
)

# Prepare ligand-target data for heatmap
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25
)

# Choose ligands and target genes for the heatmap
order_ligands <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
)

order_ligands <- rev(order_ligands)

# Select only the strongest target genes to make the heatmap readable
target_genes <- active_ligand_target_links_df %>%
  arrange(desc(weight)) %>%
  pull(target)

target_genes <- unique(target_genes)

target_genes <- head(target_genes, 100)

order_targets <- intersect(
  target_genes,
  rownames(active_ligand_target_links)
)

# Create the matrix for the ligand-target heatmap
vis_ligand_target <- active_ligand_target_links[
  order_targets,
  order_ligands
]

vis_ligand_target <- t(vis_ligand_target)

# Make ligand-target heatmap
p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  y_name = paste("Prioritized", sender_celltype, "ligands"),
  x_name = paste(receiver_celltype, "upregulated genes"),
  color = "purple",
  legend_title = "Regulatory potential"
)

p_ligand_target <- p_ligand_target +
  scale_fill_gradient2(
    low = "whitesmoke",
    high = "purple"
  ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  )

# Get ligand-receptor links for the top ligands: Ligand-Receptor links
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands,
  expressed_receptors,
  lr_network,
  weighted_networks$lr_sig
)

# Save ligand-receptor links
write_csv(
  ligand_receptor_links_df,
  file.path(output_dir, "active_ligand_receptor_links.csv")
)

# Prepare ligand-receptor data for heatmap
vis_ligand_receptor <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both"
)

# Make ligand-receptor heatmap
p_ligand_receptor <- make_heatmap_ggplot(
  t(vis_ligand_receptor),
  y_name = paste("Prioritized", sender_celltype, "ligands"),
  x_name = paste("Receptors expressed by", receiver_celltype),
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
)


### 7. Summary figures --------------------------------------------------------

# Make a heatmap matrix for ligand activity
vis_ligand_aupr <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  dplyr::select(test_ligand, aupr_corrected)

# Use ligand names as row names
vis_ligand_aupr <- as.data.frame(vis_ligand_aupr)
rownames(vis_ligand_aupr) <- vis_ligand_aupr$test_ligand
vis_ligand_aupr$test_ligand <- NULL

# Sort ligands by AUPR score
vis_ligand_aupr <- vis_ligand_aupr %>%
  arrange(aupr_corrected)

# Convert to matrix for NicheNet heatmap function
vis_ligand_aupr <- as.matrix(vis_ligand_aupr)

# Create ligand activity heatmap
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  y_name = paste("Prioritized", sender_celltype, "ligands"),
  x_name = "Ligand activity",
  color = "darkorange",
  legend_title = "AUPR"
)

# Show plots in RStudio
print(p_ligand_aupr)
print(p_ligand_target)
print(p_ligand_receptor)

# Save individual plots
ggsave(file.path(output_dir, "01_ligand_activity_aupr.png"), p_ligand_aupr, width = 5, height = 7, dpi = 300)
ggsave(file.path(output_dir, "02_ligand_target_heatmap.png"), p_ligand_target, width = 14, height = 7, dpi = 300)
ggsave(file.path(output_dir, "03_ligand_receptor_heatmap.png"), p_ligand_receptor, width = 10, height = 7, dpi = 300)

# Combine the three plots into one figure
combined_without_legends <- plot_grid(
  p_ligand_aupr,
  p_ligand_target,
  p_ligand_receptor,
  ncol = 1
)

# Show the combined figure
print(combined_without_legends)

# Save the combined figure
ggsave(
  file.path(output_dir, "04_nichenet_summary_figure.png"),
  combined_without_legends,
  width = 12,
  height = 13,
  dpi = 300
)


