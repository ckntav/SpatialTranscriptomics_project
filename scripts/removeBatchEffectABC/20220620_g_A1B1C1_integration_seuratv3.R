setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(Seurat)
# library(SeuratDisk)
# library(SeuratWrappers)

library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)

source("scripts/utils/custom_seurat_functions.R")

# https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#seurat-v3-3-vs-5-10k-pbmc

##### 1. Load data
sampleA1_data <- Read10X(data.dir = "output/sampleA1/outs/filtered_feature_bc_matrix")
sampleA1 <- CreateSeuratObject(counts = sampleA1_data, project = "A1")
sampleA1

sampleB1_data <- Read10X(data.dir = "output/sampleB1/outs/filtered_feature_bc_matrix")
sampleB1 <- CreateSeuratObject(counts = sampleB1_data, project = "B1")
sampleB1

sampleC1_data <- Read10X(data.dir = "output/sampleC1/outs/filtered_feature_bc_matrix")
sampleC1 <- CreateSeuratObject(counts = sampleC1_data, project = "C1")
sampleC1

#####
sampleA1[["percent.mt"]]  <- PercentageFeatureSet(sampleA1, pattern = "^MT-")
sampleA1[["percent.rbp"]] <- PercentageFeatureSet(sampleA1, pattern = "^RP[SL]")
VlnPlot(sampleA1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

sampleB1[["percent.mt"]]  <- PercentageFeatureSet(sampleB1, pattern = "^MT-")
sampleB1[["percent.rbp"]] <- PercentageFeatureSet(sampleB1, pattern = "^RP[SL]")
VlnPlot(sampleB1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

sampleC1[["percent.mt"]]  <- PercentageFeatureSet(sampleC1, pattern = "^MT-")
sampleC1[["percent.rbp"]] <- PercentageFeatureSet(sampleC1, pattern = "^RP[SL]")
VlnPlot(sampleC1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

#####
table(rownames(sampleA1) %in% rownames(sampleB1)) 
table(rownames(sampleA1) %in% rownames(sampleC1)) 

##### Quick filtering of the datasets removes dying cells and putative doublets
sampleA1 <- subset(sampleA1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
sampleB1 <- subset(sampleB1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
sampleC1 <- subset(sampleC1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)

#####
sample_list <- list()
sample_list[["sampleA1"]] <- sampleA1
sample_list[["sampleB1"]] <- sampleB1
sample_list[["sampleC1"]] <- sampleC1

for (i in 1:length(sample_list)) {
  sample_list[[i]] <- NormalizeData(sample_list[[i]], verbose = F)
  sample_list[[i]] <- FindVariableFeatures(sample_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}

#####
sample_anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:30)

#####
sample_seurat <- IntegrateData(anchorset = sample_anchors, dims = 1:30)

##### Let’s remove all the datastructures we’re not using to save the RAM:
rm(sample_list)
rm(sample_anchors)

##### Seurat integration creates a unified object that contains both original data (‘RNA’ assay)
# as well as integrated data (‘integrated’ assay). Let’s set the assay to RNA and visualize
# the datasets before integration.
DefaultAssay(sample_seurat) <- "RNA"

#####
sample_seurat <- NormalizeData(sample_seurat, verbose = F)
sample_seurat <- FindVariableFeatures(sample_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
sample_seurat <- ScaleData(sample_seurat, verbose = F)
sample_seurat <- RunPCA(sample_seurat, npcs = 30, verbose = F)
sample_seurat <- RunUMAP(sample_seurat, reduction = "pca", dims = 1:30, verbose = F)

#####
DimPlot(sample_seurat,reduction = "umap") + plot_annotation(title = "A1B1C1, before integration")

##### Now let’s change the assay to integrated and do the same do the same thing
# in the integrated assay (it’s already normalized and HVGs are selected)
DefaultAssay(sample_seurat) <- "integrated"
sample_seurat <- ScaleData(sample_seurat, verbose = F)
sample_seurat <- RunPCA(sample_seurat, npcs = 30, verbose = F)
sample_seurat <- RunUMAP(sample_seurat, reduction = "pca", dims = 1:30, verbose = F)

#####
DimPlot(sample_seurat, reduction = "umap") + plot_annotation(title = "A1B1C1, after integration (Seurat 3)")

#####
DimPlot(sample_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

##### Now let’s cluster the integrated matrix and look how clusters are distributed between the sets:
sample_seurat <- FindNeighbors(sample_seurat, dims = 1:30, k.param = 10, verbose = F)
sample_seurat <- FindClusters(sample_seurat, verbose = F)
DimPlot(sample_seurat,label = T) + NoLegend()

##### We can now calculate the number of cells in each cluster that came for each samplet:
count_table <- table(sample_seurat@meta.data$seurat_clusters, sample_seurat@meta.data$orig.ident)
count_table

##### Let’s plot the distribution among clusters using our custom function:
plot_integrated_clusters(sample_seurat) 

#### Puis sur les lames?