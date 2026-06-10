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
# 10. Compare selected QuickGO gene lists with DESeq2 results
# ------------------------------------------------------------

# This section compares genes from three GO terms:
# GO:0007520 - myoblast fusion
# GO:0051450 - myoblast proliferation
# GO:0002040 - sprouting angiogenesis
#
# The goal is to check whether these GO genes:
# 1. are present in the SkMVEC vs HUVEC dataset
# 2. are expressed in SkMVECs and/or HUVECs
# 3. are higher in SkMVECs or higher in HUVECs
# 4. are statistically significantly different
#
# This helps compare biological themes that are relevant for the project:
# muscle-cell fusion, muscle-cell proliferation and endothelial sprouting.

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
# 10.1 Define GO gene lists
# ------------------------------------------------------------

# GO:0007520 - myoblast fusion
# Homo sapiens annotations from QuickGO

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

# GO:0051450 - myoblast proliferation
# Homo sapiens annotations from QuickGO

go_myoblast_proliferation_genes <- c(
  "ABL1",
  "ATOH8",
  "FES",
  "FOS",
  "IGF1",
  "KRAS"
)

# GO:0002040 - sprouting angiogenesis
# Homo sapiens annotations from QuickGO

go_sprouting_angiogenesis_genes <- c(
  "ADGRA2",
  "ANGPT1",
  "CCBE1",
  "CDH13",
  "E2F7",
  "E2F8",
  "ENG",
  "ESM1",
  "FLT1",
  "FLT4",
  "JMJD6",
  "LOXL2",
  "NAXE",
  "NRP1",
  "OTULIN",
  "PARVA",
  "PGF",
  "PTK2B",
  "RAMP2",
  "RECK",
  "RNF213",
  "RSPO3",
  "SEMA3E",
  "TEK",
  "TGFB1",
  "THBS1",
  "TMEM215",
  "VEGFA",
  "VEGFB",
  "VEGFC",
  "VEGFD",
  "YJEFN3"
)

# ------------------------------------------------------------
# 10.2 Function to compare a GO gene list with SkMVEC vs HUVEC results
# ------------------------------------------------------------

# This function compares one GO gene list with the DESeq2 results.
# This way, we do not have to copy the same code three times.
#
# For every GO gene, the function checks:
# 1. Is the gene present in the DESeq2 result table?
# 2. Is the gene expressed in SkMVECs?
# 3. Is the gene expressed in HUVECs?
# 4. Is the gene higher in SkMVECs or HUVECs?
# 5. Is the difference statistically significant?

compare_go_genes <- function(go_genes, go_id, go_description) {
  
  # ------------------------------------------------------------
  # Prepare the full DESeq2 results table
  # ------------------------------------------------------------
  
  # We use the full results_table, not only significant genes.
  # This is important because a GO gene can be present in the dataset
  # even if it is not significantly different between SkMVECs and HUVECs.
  
  results_for_go <- as.data.frame(results_table)
  
  # If gene names are not already in a column, add them from the row names.
  if (!"gene" %in% colnames(results_for_go)) {
    results_for_go <- rownames_to_column(results_for_go, var = "gene")
  }
  
  # Clean gene names by removing extra spaces.
  results_for_go <- results_for_go %>%
    mutate(
      gene_match = trimws(as.character(gene))
    )
  
  # ------------------------------------------------------------
  # Create a table from the GO gene list
  # ------------------------------------------------------------
  
  go_table <- tibble(
    go_term = go_id,
    go_description = go_description,
    go_gene = go_genes,
    gene_match = trimws(go_genes)
  )
  
  # ------------------------------------------------------------
  # Match GO genes with the DESeq2 results
  # ------------------------------------------------------------
  
  go_comparison <- go_table %>%
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
  # Add mean normalized expression per cell type
  # ------------------------------------------------------------
  
  # Normalized counts are corrected for sequencing depth.
  # This allows us to compare expression values between samples more fairly.
  
  normalized_count_matrix <- counts(dds, normalized = TRUE)
  sample_information <- as.data.frame(colData(dds))
  
  # Select SkMVEC and HUVEC samples.
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
  
  go_comparison <- go_comparison %>%
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
  # Add expressed / not expressed information
  # ------------------------------------------------------------
  
  # A gene is called expressed if it has at least 10 raw reads
  # in at least 2 samples of that cell type.
  
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
  
  go_comparison <- go_comparison %>%
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
  # Add automatic interpretation
  # ------------------------------------------------------------
  
  # Because the DESeq2 contrast is SkMVECs vs HUVECs:
  # positive log2FoldChange = higher in SkMVECs
  # negative log2FoldChange = higher in HUVECs
  
  go_comparison <- go_comparison %>%
    mutate(
      interpretation = case_when(
        !found_in_results_table ~
          "Not found in the filtered DESeq2 dataset",
        
        is.na(padj) ~
          "Detected, but no adjusted p-value available",
        
        padj <= 0.05 & log2FoldChange >= 1 ~
          "Significantly higher in SkMVECs",
        
        padj <= 0.05 & log2FoldChange <= -1 ~
          "Significantly higher in HUVECs",
        
        padj <= 0.05 & abs(log2FoldChange) < 1 ~
          "Statistically significant, but absolute log2FoldChange is smaller than 1",
        
        padj > 0.05 & log2FoldChange > 0 ~
          "Higher in SkMVECs, but not statistically significant",
        
        padj > 0.05 & log2FoldChange < 0 ~
          "Higher in HUVECs, but not statistically significant",
        
        TRUE ~
          "No clear difference"
      )
    )
  
  # ------------------------------------------------------------
  # Save full table
  # ------------------------------------------------------------
  
  # Make a clean file name.
  file_name <- paste0(go_id, "_", go_description)
  file_name <- gsub(":", "_", file_name)
  file_name <- gsub(" ", "_", file_name)
  
  write_csv(
    go_comparison,
    file.path(
      output_directory,
      paste0(file_name, "_complete_comparison.csv")
    )
  )
  
  # ------------------------------------------------------------
  # Save smaller Word-friendly table
  # ------------------------------------------------------------
  
  # This table is easier to copy into Word or use in your EIN.
  word_table <- go_comparison %>%
    dplyr::select(
      go_gene,
      mean_normalized_count_SkMVECs,
      mean_normalized_count_HUVECs,
      log2FoldChange,
      padj,
      interpretation
    ) %>%
    dplyr::rename(
      `Mean SkMVECs` = mean_normalized_count_SkMVECs,
      `Mean HUVECs` = mean_normalized_count_HUVECs
    ) %>%
    mutate(
      `Mean SkMVECs` = round(`Mean SkMVECs`, 1),
      `Mean HUVECs` = round(`Mean HUVECs`, 1),
      log2FoldChange = round(log2FoldChange, 2),
      padj = signif(padj, 3)
    )
  
  write_csv(
    word_table,
    file.path(
      output_directory,
      paste0(file_name, "_word_table.csv")
    )
  )
  
  return(go_comparison)
}

# ------------------------------------------------------------
# 10.3 Run GO gene comparisons
# ------------------------------------------------------------

# First GO term: myoblast fusion
myoblast_fusion_comparison <- compare_go_genes(
  go_genes = go_myoblast_fusion_genes,
  go_id = "GO:0007520",
  go_description = "myoblast fusion"
)

# Second GO term: myoblast proliferation
myoblast_proliferation_comparison <- compare_go_genes(
  go_genes = go_myoblast_proliferation_genes,
  go_id = "GO:0051450",
  go_description = "myoblast proliferation"
)

# Third GO term: sprouting angiogenesis
sprouting_angiogenesis_comparison <- compare_go_genes(
  go_genes = go_sprouting_angiogenesis_genes,
  go_id = "GO:0002040",
  go_description = "sprouting angiogenesis"
)


# ------------------------------------------------------------
# 10.4 Create summary table for the three biological themes
# ------------------------------------------------------------

# This table summarizes the three GO themes.
# It is useful for your EIN because it gives a quick overview.

make_go_summary <- function(go_comparison, go_id, go_description) {
  
  summary_table <- tibble(
    go_term = go_id,
    go_description = go_description,
    
    number_of_genes_in_GO_list = nrow(go_comparison),
    
    found_in_DESeq2_results = sum(
      go_comparison$found_in_results_table,
      na.rm = TRUE
    ),
    
    expressed_in_SkMVECs = sum(
      go_comparison$expressed_in_SkMVECs,
      na.rm = TRUE
    ),
    
    expressed_in_HUVECs = sum(
      go_comparison$expressed_in_HUVECs,
      na.rm = TRUE
    ),
    
    significantly_higher_in_SkMVECs = sum(
      go_comparison$padj <= 0.05 &
        go_comparison$log2FoldChange >= 1,
      na.rm = TRUE
    ),
    
    significantly_higher_in_HUVECs = sum(
      go_comparison$padj <= 0.05 &
        go_comparison$log2FoldChange <= -1,
      na.rm = TRUE
    )
  )
  
  return(summary_table)
}

go_theme_summary <- bind_rows(
  make_go_summary(
    myoblast_fusion_comparison,
    "GO:0007520",
    "myoblast fusion"
  ),
  
  make_go_summary(
    myoblast_proliferation_comparison,
    "GO:0051450",
    "myoblast proliferation"
  ),
  
  make_go_summary(
    sprouting_angiogenesis_comparison,
    "GO:0002040",
    "sprouting angiogenesis"
  )
)

write_csv(
  go_theme_summary,
  file.path(
    output_directory,
    "GO_theme_summary_SkMVEC_vs_HUVEC.csv"
  )
)

View(go_theme_summary)


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




