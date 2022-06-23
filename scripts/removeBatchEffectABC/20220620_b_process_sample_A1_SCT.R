setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium

# Load dataset
sampleA1 <- Load10X_Spatial(data.dir = "output/sampleA1/outs", slice = "sliceA1")
sampleA1
names(sampleA1)

# Data preprocessing
plot1_A1 <- VlnPlot(sampleA1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2_A1 <- SpatialFeaturePlot(sampleA1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1_A1, plot2_A1)

sampleA1_sct <- SCTransform(sampleA1, assay = "Spatial", verbose = TRUE)

isg15 <- SpatialFeaturePlot(sampleA1, features = "ISG15")
isg15_sct <- SpatialFeaturePlot(sampleA1_sct, features = "ISG15")
wrap_plots(isg15, isg15_sct)

# Dimensionality reduction, clustering, and visualization
sampleA1 <- SCTransform(sampleA1, assay = "Spatial", verbose = TRUE)
sampleA1 <- RunPCA(sampleA1, assay = "SCT", verbose = FALSE)
DimPlot(sampleA1)
sampleA1 <- FindNeighbors(sampleA1, reduction = "pca", dims = 1:30)
sampleA1 <- FindClusters(sampleA1, verbose = FALSE)
sampleA1 <- RunUMAP(sampleA1, reduction = "pca", dims = 1:30)

p1 <- DimPlot(sampleA1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(sampleA1, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(sampleA1, cells.highlight = CellsByIdentities(object = sampleA1, idents = c(1, 2, 3, 4, 5, 6, 7 ,8)), facet.highlight = TRUE, ncol = 3)

# Identification of Spatially Variable Features
sampleA1 <- FindSpatiallyVariableFeatures(sampleA1, assay = "SCT", features = VariableFeatures(sampleA1)[1:1000], selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(sampleA1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(sampleA1, features = top.features, ncol = 3, alpha = c(0.1, 1))