# ============================================================
# Differential expression analysis:
# confluent SkMVECs versus confluent HUVECs
# ============================================================
#
# This script uses only the 10 selected samples:
# - 5 confluent SkMVEC samples
# - 5 confluent HUVEC samples
#
# HUVECs are used as the reference group.
# This means:
# positive log2FoldChange = higher in SkMVECs
# negative log2FoldChange = higher in HUVECs

# Install once if needed:
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages(c("ggplot2", "pheatmap"))

# Load packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)


# ------------------------------------------------------------
# 1. Define folders
# ------------------------------------------------------------

project_dir <- "C:/Users/royde/OneDrive/Documenten/BIT11 Internship"

count_directory <- file.path(
  project_dir,
  "data",
  "Confluent_SkMVECs_vs_Confluent_HUVECs"
)

output_directory <- file.path(
  project_dir,
  "results",
  "SkMVEC_vs_HUVEC_DESeq2_results"
)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Create the sample-information table
# ------------------------------------------------------------

sample_table <- data.frame(
  sampleName = c(
    "GC122786_KG",
    "GC122787_KG",
    "GC122788_KG",
    "GC124881_KG",
    "GC122790_KG",
    "GC122791_KG",
    "GC122792_KG",
    "GC122793_KG",
    "GC122794_KG",
    "GC122795_KG"
  ),

  fileName = c(
    "GC122786_KG.count",
    "GC122787_KG.count",
    "GC122788_KG.count",
    "GC124881_KG.count",
    "GC122790_KG.count",
    "GC122791_KG.count",
    "GC122792_KG.count",
    "GC122793_KG.count",
    "GC122794_KG.count",
    "GC122795_KG.count"
  ),

  condition = c(
    rep("confluent_SkMVECs", 5),
    rep("confluent_HUVECs", 5)
  )
)

rownames(sample_table) <- sample_table$sampleName

# Set HUVECs as the reference group
sample_table$condition <- factor(
  sample_table$condition,
  levels = c("confluent_HUVECs", "confluent_SkMVECs")
)

write.csv(
  sample_table,
  file.path(output_directory, "sample_metadata.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Check whether all input files exist
# ------------------------------------------------------------

file_paths <- file.path(count_directory, sample_table$fileName)

if (!all(file.exists(file_paths))) {
  missing_files <- file_paths[!file.exists(file_paths)]
  stop(
    "One or more count files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 4. Import the HTSeq-count files
# ------------------------------------------------------------

# DESeqDataSetFromHTSeqCount reads the separate .count files
# and combines them into one DESeq2 object.
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_table,
  directory = count_directory,
  design = ~ condition
)


# ------------------------------------------------------------
# 5. Remove genes with very low counts
# ------------------------------------------------------------

# Keep genes that have at least 10 reads in at least 2 samples.
# This removes genes that are probably too low to analyse reliably.
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]


# ------------------------------------------------------------
# 6. Run the differential-expression analysis
# ------------------------------------------------------------

dds <- DESeq(dds)

# Compare confluent SkMVECs with confluent HUVECs.
results_object <- results(
  dds,
  contrast = c(
    "condition",
    "confluent_SkMVECs",
    "confluent_HUVECs"
  ),
  alpha = 0.05
)

# Convert the results to a normal table and sort by adjusted p-value.
results_table <- as.data.frame(results_object)
results_table$gene <- rownames(results_table)

results_table <- results_table[
  order(results_table$padj),
]

results_table <- results_table[
  c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
]


# ------------------------------------------------------------
# 7. Save all differential-expression results
# ------------------------------------------------------------

write.csv(
  results_table,
  file.path(
    output_directory,
    "confluent_SkMVECs_vs_confluent_HUVECs_all_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Save significant genes
# ------------------------------------------------------------

significant_genes <- results_table[
  !is.na(results_table$padj) &
    results_table$padj < 0.05 &
    abs(results_table$log2FoldChange) >= 1,
]

write.csv(
  significant_genes,
  file.path(
    output_directory,
    "confluent_SkMVECs_vs_confluent_HUVECs_significant_genes.csv"
  ),
  row.names = FALSE
)

# Save SkMVEC-upregulated genes separately
skmvec_upregulated_genes <- significant_genes[
  significant_genes$log2FoldChange >= 1,
]

write.csv(
  skmvec_upregulated_genes,
  file.path(
    output_directory,
    "confluent_SkMVECs_upregulated_vs_HUVECs.csv"
  ),
  row.names = FALSE
)

# Save HUVEC-upregulated genes separately
huvec_upregulated_genes <- significant_genes[
  significant_genes$log2FoldChange <= -1,
]

write.csv(
  huvec_upregulated_genes,
  file.path(
    output_directory,
    "confluent_HUVECs_upregulated_vs_SkMVECs.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Save normalized counts
# ------------------------------------------------------------

normalized_counts <- as.data.frame(
  counts(dds, normalized = TRUE)
)

normalized_counts$gene <- rownames(normalized_counts)

normalized_counts <- normalized_counts[
  c("gene", sample_table$sampleName)
]

write.csv(
  normalized_counts,
  file.path(output_directory, "normalized_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Compare QuickGO myoblast fusion genes with DESeq2 results
# ------------------------------------------------------------

# This section compares genes from the GO terms:
# GO:0007520 - myoblast fusion
# 
# The goal is to check whether these myoblast-fusion genes:
# 1. are present in the SkMVEC vs HUVEC dataset
# 2. are expressed in SkMVECs and/or HUVECs
# 3. are higher in SkMVECs or higher in HUVECs
# 4. are statistically significantly different

# Important:
# In this script, the contrast is:
# confluent_SkMVECs versus confluent_HUVECs
#
# Therefore:
# positive log2FoldChange = higher in SkMVECs
# negative log2FoldChange = higher in HUVECs

positive_lfc_group <- "SkMVECs"
negative_lfc_group <- "HUVECs"


# ------------------------------------------------------------
# 10.1 Add the QuickGO gene list
# ------------------------------------------------------------

## GO:0007520 - Myoblast fusion
## Homo Sapiens annotations from QuickGO

go_myoblast_fusion_genes <- c(
  "ADAM12",
  "ADAMTS15",
  "ADAMTS5",
  "ADGRB3",
  "CACNA1H",
  "CD81",
  "CD9",
  "DOCK1",
  "DOCK2",
  "DOCK5",
  "ERVFRD-1",
  "ERVW-1",
  "FER1L5",
  "ITGB1",
  "KCNH1",
  "MYMK",
  "MYMX",
  "NOS1",
  "NPHS1",
  "PTGFRN"
)

# ------------------------------------------------------------
# 10.2 Prepare the DESeq2 results table
# ------------------------------------------------------------

# Make a copy of the full DESeq2 results table.
# We use the full table, not only the significant genes,
# because some QuickGO genes may be present but not significant.

results_for_go <- as.data.frame(results_table)

# If gene names are not stored in a gene column,
# add them from the row names.
if (!"gene" %in% colnames(results_for_go)) {
  results_for_go <- rownames_to_column(results_for_go, var = "gene")
}

# Inspect the first genes and the availalbe columns
head(results_for_go)
colnames(results_for_go)
head(results_for_go$gene,20)

# Clean the gene names by removing extra spaces.
# This makes it easier to match the QuickGO genes with the DESeq2 results.

results_for_go <- results_for_go %>%
  mutate(
    gene_match = trimws(as.character(gene))
  )

# ------------------------------------------------------------
# 10.3 Create a table from the GO gene list
# ------------------------------------------------------------

go_myoblast_fusion_table <- tibble(
  go_term = "GO:0007520",
  go_description = "myoblast fusion",
  go_gene = go_myoblast_fusion_genes,
  gene_match = trimws(go_myoblast_fusion_genes)
)

# ------------------------------------------------------------
# 10.4 Compare GO genes with the DESeq2 results
# ------------------------------------------------------------

go_myoblast_fusion_comparison <- go_myoblast_fusion_table %>%
  left_join(
    results_for_go %>%
      dplyr::select(
        gene_match,
        matched_gene = gene,
        baseMean,
        log2FoldChange,
        lfcSE,
        pvalue,
        padj
      ),
    by = "gene_match"
  ) %>%
  mutate(
    found_in_results_table = !is.na(matched_gene),
    statistically_significant = !is.na(padj) & padj <= 0.05
  )


# ------------------------------------------------------------
# 10.5 Add mean normalized expression per cell type
# ------------------------------------------------------------

# DESeq2 normalized counts correct for differences in sequencing depth.
# This allows us to compare average expression between samples more fairly.

normalized_count_matrix <- counts(dds, normalized = TRUE)

sample_information <- as.data.frame(colData(dds))

# Select the SkMVEC and HUVEC samples based on the condition column.
skmvec_samples <- rownames(sample_information)[
  sample_information$condition == "confluent_SkMVECs"
]

huvec_samples <- rownames(sample_information)[
  sample_information$condition == "confluent_HUVECs"
]

expression_summary <- tibble(
  gene = rownames(normalized_count_matrix),
  
  mean_normalized_count_SkMVECs = rowMeans(
    normalized_count_matrix[, skmvec_samples, drop = FALSE]
  ),
  
  mean_normalized_count_HUVECs = rowMeans(
    normalized_count_matrix[, huvec_samples, drop = FALSE]
  )
) %>%
  mutate(
    gene_match = trimws(as.character(gene))
  )

go_myoblast_fusion_comparison <- go_myoblast_fusion_comparison %>%
  left_join(
    expression_summary %>%
      dplyr::select(
        gene_match,
        mean_normalized_count_SkMVECs,
        mean_normalized_count_HUVECs
      ),
    by = "gene_match"
  )


# ------------------------------------------------------------
# 10.6 Add expressed / not expressed classification
# ------------------------------------------------------------

# Here, a gene is called expressed if it has at least 10 raw reads
# in at least 2 samples of that cell type.
#
# This is similar to the low-count filtering used earlier.

raw_count_matrix <- counts(dds, normalized = FALSE)

min_raw_count <- 10
min_number_of_samples <- 2

expression_detection_summary <- tibble(
  gene = rownames(raw_count_matrix),
  
  number_of_SkMVEC_samples_expressed = rowSums(
    raw_count_matrix[, skmvec_samples, drop = FALSE] >= min_raw_count
  ),
  
  number_of_HUVEC_samples_expressed = rowSums(
    raw_count_matrix[, huvec_samples, drop = FALSE] >= min_raw_count
  )
) %>%
  mutate(
    gene_match = trimws(as.character(gene)),
    
    expressed_in_SkMVECs =
      number_of_SkMVEC_samples_expressed >= min_number_of_samples,
    
    expressed_in_HUVECs =
      number_of_HUVEC_samples_expressed >= min_number_of_samples
  )

go_myoblast_fusion_comparison <- go_myoblast_fusion_comparison %>%
  left_join(
    expression_detection_summary %>%
      dplyr::select(
        gene_match,
        number_of_SkMVEC_samples_expressed,
        number_of_HUVEC_samples_expressed,
        expressed_in_SkMVECs,
        expressed_in_HUVECs
      ),
    by = "gene_match"
  )


# ------------------------------------------------------------
# 10.7 Add an automatic interpretation
# ------------------------------------------------------------

go_myoblast_fusion_comparison <- go_myoblast_fusion_comparison %>%
  mutate(
    interpretation = case_when(
      !found_in_results_table ~
        "Not found in the filtered DESeq2 dataset",
      
      is.na(padj) ~
        "Detected, but no adjusted p-value available",
      
      padj <= 0.05 & log2FoldChange >= 1 ~
        paste0("Significantly higher in ", positive_lfc_group),
      
      padj <= 0.05 & log2FoldChange <= -1 ~
        paste0("Significantly higher in ", negative_lfc_group),
      
      padj <= 0.05 & abs(log2FoldChange) < 1 ~
        "Statistically significant, but absolute log2FoldChange is smaller than 1",
      
      padj > 0.05 & log2FoldChange > 0 ~
        paste0("Higher in ", positive_lfc_group, ", but not statistically significant"),
      
      padj > 0.05 & log2FoldChange < 0 ~
        paste0("Higher in ", negative_lfc_group, ", but not statistically significant"),
      
      TRUE ~
        "No clear difference"
    )
  )

################## Compact table
go_myoblast_fusion_word_table <- go_myoblast_fusion_comparison %>%
  select(
    go_gene,
    found_in_results_table,
    baseMean,
    log2FoldChange,
    padj,
    statistically_significant,
    interpretation
  ) %>%
  arrange(desc(log2FoldChange))

go_myoblast_fusion_word_table <- go_myoblast_fusion_word_table %>%
  rename(
    Gene = go_gene,
    `Found in dataset` = found_in_results_table,
    `Mean normalized expression` = baseMean,
    `Log2 fold change` = log2FoldChange,
    `Adjusted p-value` = padj,
    `Statistically significant` = statistically_significant,
    Interpretation = interpretation
  )

go_myoblast_fusion_word_table <- go_myoblast_fusion_word_table %>%
  mutate(
    `Mean normalized expression` =
      round(`Mean normalized expression`, 2),
    
    `Log2 fold change` =
      round(`Log2 fold change`, 2),
    
    `Adjusted p-value` =
      signif(`Adjusted p-value`, 3)
  )

view(go_myoblast_fusion_word_table)
##################

# ------------------------------------------------------------
# 10.8 Save the final comparison table
# ------------------------------------------------------------

write_csv(
  go_myoblast_fusion_comparison,
  file.path(
    output_directory,
    "GO_0007520_myoblast_fusion_complete_comparison.csv"
  )
)

# Open the table in RStudio
View(go_myoblast_fusion_comparison)


# ------------------------------------------------------------
# 11. Inspect normalized counts for myoblast-fusion genes
# ------------------------------------------------------------

# The normalized_counts table contains:
# one gene column
# ten sample columns
#
# These values are useful to compare expression levels between samples,
# because DESeq2 corrects for differences in sequencing depth.

normalized_counts_for_go <- normalized_counts %>%
  filter(gene %in% go_myoblast_fusion_genes)

# Save this table so it can be opened in Excel.
write_csv(
  normalized_counts_for_go,
  file.path(
    output_directory,
    "GO_0007520_myoblast_fusion_normalized_counts_per_sample.csv"
  )
)

View(normalized_counts_for_go)

# Convert the normalized count table to long format.
# This makes it easier to plot in ggplot2.

normalized_counts_long <- normalized_counts_for_go %>%
  pivot_longer(
    cols = -gene,
    names_to = "sampleName",
    values_to = "normalized_count"
  ) %>%
  left_join(
    sample_table,
    by = "sampleName"
  )

# Calculate the mean normalized expression per gene per cell type.

go_expression_summary <- normalized_counts_long %>%
  group_by(gene, condition) %>%
  summarise(
    mean_normalized_count = mean(normalized_count, na.rm = TRUE),
    sd_normalized_count = sd(normalized_count, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  go_expression_summary,
  file.path(
    output_directory,
    "GO_0007520_myoblast_fusion_mean_expression_per_cell_type.csv"
  )
)

View(go_expression_summary)

go_myoblast_fusion_word_table <- go_myoblast_fusion_comparison %>%
  dplyr::select(
    go_gene,
    mean_normalized_count_SkMVECs,
    mean_normalized_count_HUVECs,
    log2FoldChange,
    padj
  ) %>%
  dplyr::rename(
    `Mean SkMVECs` = mean_normalized_count_SkMVECs,
    `Mean HUVECs` = mean_normalized_count_HUVECs
  ) %>%
  dplyr::mutate(
    `Mean SkMVECs` = round(`Mean SkMVECs`, 1),
    `Mean HUVECs` = round(`Mean HUVECs`, 1),
    log2FoldChange = round(log2FoldChange, 2),
    padj = signif(padj, 3)
  )

# Show the clean table in RStudio
View(go_myoblast_fusion_word_table)

# Save the table as a CSV file
write_csv(
  go_myoblast_fusion_word_table,
  file.path(
    output_directory,
    "GO_0007520_myoblast_fusion_word_table.csv"
  )
)

# ------------------------------------------------------------
# 11. Create a PCA plot
# ------------------------------------------------------------

# vst transforms the count data so it can be used better for PCA.
vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
) +
  ggtitle("PCA of confluent SkMVECs and HUVECs")

print(pca_plot)

ggsave(
  file.path(output_directory, "PCA_plot.pdf"),
  plot = pca_plot,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_directory, "PCA_plot.png"),
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 12. Create a sample-distance heatmap
# ------------------------------------------------------------

sample_distances <- dist(t(assay(vsd)))
sample_distance_matrix <- as.matrix(sample_distances)

annotation_table <- data.frame(
  condition = sample_table$condition
)

rownames(annotation_table) <- sample_table$sampleName

pdf(
  file.path(output_directory, "sample_distance_heatmap.pdf"),
  width = 8,
  height = 7
)

pheatmap(
  sample_distance_matrix,
  annotation_col = annotation_table,
  annotation_row = annotation_table
)

dev.off()


# ------------------------------------------------------------
# 13. Create an MA plot
# ------------------------------------------------------------

pdf(
  file.path(output_directory, "MA_plot.pdf"),
  width = 7,
  height = 5
)

plotMA(results_object, ylim = c(-5, 5))

dev.off()




