setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(SingleCellExperiment)
library(SingleR)
library(Seurat)
library(tidyverse)
source("scripts/utils/ckn_utils_savePlot.R")

#
today <- get_today()

#
method_nb <- "method2"
spt_data_filepath <- file.path("output/ABCD_object", paste0("20220626_ABCD_object_", method_nb, ".rds"))
spt_data <- readRDS(spt_data_filepath)

# ref
# blueprintEncode.ref <- celldex::BlueprintEncodeData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()
# monaco.ref <- celldex::MonacoImmuneData()

# Let’s convert our Seurat object to single cell experiment (SCE) for convenience.
# After this, using SingleR becomes very easy:
# sce <- as.SingleCellExperiment(DietSeurat(spt_data))
sce <- as.SingleCellExperiment(spt_data)
sce

# cell type nnotations
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
hpca.main <- SingleR(test = sce, ref = hpca.ref, labels = hpca.ref$label.main)


# hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
# saveRDS(hpca.main, file = file.path("output/annotation_clusters_celltype", paste0(today, "_", method_nb, "_annotation_clusters_celltype_hpca_main.rds")))
# saveRDS(hpca.fine, file = file.path("output/annotation_clusters_celltype", paste0(today, "_", method_nb, "_annotation_clusters_celltype_hpca_fine.rds")))
hpca.main <- readRDS(file = file.path("output/annotation_clusters_celltype", paste0("20220628_", method_nb, "_annotation_clusters_celltype_hpca_main.rds")))

# summary of cell type annotations
table(hpca.main$pruned.labels)
table(hpca.fine$pruned.labels)

# Let’s add the annotations to the Seurat object metadata so we can use them:
spt_data@meta.data$hpca.main <- hpca.main$pruned.labels
spt_data@meta.data$hpca.fine <- hpca.fine$pruned.labels

# Finally, let’s visualize the annotations.
plot1 <- DimPlot(spt_data, label = T , repel = T, label.size = 3, group.by = "hpca.main")
plot2 <- DimPlot(spt_data, label = T , repel = T, label.size = 3, group.by = "ident")
plot1 + plot2

plot3 <- DimPlot(spt_data, label = T , repel = T, label.size = 3, group.by = "hpca.fine")
plot4 <- DimPlot(spt_data, label = T , repel = T, label.size = 3, group.by = "ident") + NoLegend()
plot3 + plot4
# plotScoreHeatmap(hpca.main)




library(tidyverse)
source("scripts/utils/custom_seurat_functions.R")






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
