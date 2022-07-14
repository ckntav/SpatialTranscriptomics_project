setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(tidyverse)
source("scripts/utils/custom_seurat_functions.R")
source("scripts/utils/ckn_utils_savePlot.R")

#
today <- get_today()

#
method_nb <- "method2"
spt_data_filepath <- file.path("output/ABCD_object", paste0("20220626_ABCD_object_", method_nb, ".rds"))
spt_data <- readRDS(spt_data_filepath)

# Let’s plot the distribution among clusters using our custom function:
plot_integrated_clusters(spt_data) 

clusters_index <- c(0, 2, 4, 6, 7, 9, 11, 13, 15)

sliceA1 <- 
SpatialDimPlot(spt_data, ncol = 3, images = "sliceA1",
               cells.highlight = CellsByIdentities(object = spt_data, idents = clusters_index),
               facet.highlight = TRUE) +
  plot_annotation(title = "sliceA1 | Harmony")

sliceB1 <- 
SpatialDimPlot(spt_data, ncol = 3, images = "sliceB1",
               cells.highlight = CellsByIdentities(object = spt_data, idents = clusters_index),
               facet.highlight = TRUE) +
  plot_annotation(title = "sliceB1 | Harmony")

#
savePlotPNG(plot = sliceA1, output_dir = "output/viz", output_file = paste(today, "sliceA1", "commonClustersAB", sep = "_"),
            width_val = 7, height_val = 7, res_val = 500)

#
savePlotPNG(plot = sliceB1, output_dir = "output/viz", output_file = paste(today, "sliceB1", "commonClustersAB", sep = "_"),
            width_val = 7, height_val = 7, res_val = 500)
