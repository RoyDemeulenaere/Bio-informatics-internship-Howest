# Circos / chord plot for NicheNet ligand-target results
#
# Run this script after running NicheNet_assignment_Myotubes_SkMVECs.R.
# It uses the NicheNet output file active_ligand_target_links.csv.

library(tidyverse)
library(circlize)

# Set folders
project_dir <- "C:/Users/royde/OneDrive/Documenten/BIT11 Internship"
output_dir <- file.path(project_dir, "results", "nichenet_assignment_results")

# Read ligand-target links from NicheNet
ligand_target_links <- read_csv(
  file.path(output_dir, "active_ligand_target_links.csv"),
  show_col_types = FALSE
)

# Select the top ligands based on total interaction weight
top_ligands <- ligand_target_links %>%
  group_by(ligand) %>%
  summarise(total_weight = sum(weight, na.rm = TRUE)) %>%
  arrange(desc(total_weight)) %>%
  slice_head(n = 10) %>%
  pull(ligand)

# Select the top 5 target genes for each top ligand
chord_links <- ligand_target_links %>%
  filter(ligand %in% top_ligands) %>%
  arrange(desc(weight)) %>%
  group_by(ligand) %>%
  slice_head(n = 5) %>%
  ungroup()

# Keep only the columns needed for the chord diagram
chord_links <- chord_links %>%
  select(
    from = ligand,
    to = target,
    value = weight
  )

# Define colors
ligands <- unique(chord_links$from)
targets <- unique(chord_links$to)

ligand_colors <- setNames(rep("dodgerblue3", length(ligands)), ligands)
target_colors <- setNames(rep("tomato2", length(targets)), targets)

node_colors <- c(ligand_colors, target_colors)

circos.clear()

chordDiagram(
  chord_links,
  grid.col = node_colors,
  transparency = 0.35,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(
      x = mean(xlim),
      y = ylim[1],
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0.5, 0.5),
      cex = 0.65
    )
  },
  bg.border = NA
)

title("Predicted ligand-target interactions")

legend(
  "bottomleft",
  legend = c("SkMVEC ligands", "Myotube target genes"),
  fill = c("dodgerblue3", "tomato2"),
  border = NA,
  bty = "n",
  cex = 1.1
)

circos.clear()

