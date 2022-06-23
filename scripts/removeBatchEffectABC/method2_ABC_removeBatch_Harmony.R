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

# https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#harmony-3-vs-5-10k-pbmc

##### 1. Load data
# sampleA1_data <- Read10X(data.dir = "output/sampleA1/outs/filtered_feature_bc_matrix")
sampleA1 <- Load10X_Spatial(data.dir = "output/sampleA1/outs", slice = "sliceA1")
sampleA1$orig.ident <- "sampleA1"
sampleA1

# sampleB1_data <- Read10X(data.dir = "output/sampleB1/outs/filtered_feature_bc_matrix")
sampleB1 <- Load10X_Spatial(data.dir = "output/sampleB1/outs", slice = "sliceB1")
sampleB1$orig.ident <- "sampleB1"
sampleB1

# sampleC1_data <- Read10X(data.dir = "output/sampleC1/outs/filtered_feature_bc_matrix")
sampleC1 <- Load10X_Spatial(data.dir = "output/sampleC1/outs", slice = "sliceC1")
sampleC1$orig.ident <- "sampleC1"
sampleC1

#####
sampleA1[["percent.mt"]]  <- PercentageFeatureSet(sampleA1, pattern = "^MT-")
sampleA1[["percent.rbp"]] <- PercentageFeatureSet(sampleA1, pattern = "^RP[SL]")
VlnPlot(sampleA1, features = c("nFeature_Spatial","nCount_Spatial","percent.mt","percent.rbp"), group.by = "orig.ident", ncol = 4)

sampleB1[["percent.mt"]]  <- PercentageFeatureSet(sampleB1, pattern = "^MT-")
sampleB1[["percent.rbp"]] <- PercentageFeatureSet(sampleB1, pattern = "^RP[SL]")
VlnPlot(sampleB1, features = c("nFeature_Spatial","nCount_Spatial","percent.mt","percent.rbp"), group.by = "orig.ident", ncol = 4)

sampleC1[["percent.mt"]]  <- PercentageFeatureSet(sampleC1, pattern = "^MT-")
sampleC1[["percent.rbp"]] <- PercentageFeatureSet(sampleC1, pattern = "^RP[SL]")
VlnPlot(sampleC1, features = c("nFeature_Spatial","nCount_Spatial","percent.mt","percent.rbp"), group.by = "orig.ident", ncol = 4)

#####
table(rownames(sampleA1) %in% rownames(sampleB1)) 
table(rownames(sampleA1) %in% rownames(sampleC1)) 

##### Quick filtering of the datasets removes dying cells and putative doublets
# sampleA1 <- subset(sampleA1, subset = nFeature_Spatial > 500 & nFeature_Spatial < 5000 & percent.mt < 15)
# sampleB1 <- subset(sampleB1, subset = nFeature_Spatial > 500 & nFeature_Spatial < 5000 & percent.mt < 15)
# sampleC1 <- subset(sampleC1, subset = nFeature_Spatial > 500 & nFeature_Spatial < 5000 & percent.mt < 15)
sampleA1 <- subset(sampleA1, subset = nFeature_Spatial > 500 & percent.mt < 15)
sampleB1 <- subset(sampleB1, subset = nFeature_Spatial > 500 & percent.mt < 15)
sampleC1 <- subset(sampleC1, subset = nFeature_Spatial > 500 & percent.mt < 15)

##### 
ABC_harmony <- merge(sampleA1, y = c(sampleB1, sampleC1))

#####
ABC_harmony <- NormalizeData(ABC_harmony, verbose = F)
ABC_harmony <- FindVariableFeatures(ABC_harmony, selection.method = "vst", nfeatures = 2000, verbose = F)
ABC_harmony <- ScaleData(ABC_harmony, verbose = F)
ABC_harmony <- RunPCA(ABC_harmony, npcs = 30, verbose = F)
ABC_harmony <- RunUMAP(ABC_harmony, reduction = "pca", dims = 1:30, verbose = F)

#####
DimPlot(ABC_harmony, reduction = "umap", group.by = "orig.ident") + plot_annotation(title = "A1B1C1, before integration, before removing batch effect")

#####
ABC_harmony <- ABC_harmony %>% RunHarmony("orig.ident", plot_convergence = T, assay.use = "Spatial")

##### Check the generated embeddings:
harmony_embeddings <- Embeddings(ABC_harmony, 'harmony')
harmony_embeddings[1:5, 1:5]

##### Check the PCA plot after:
p1 <- DimPlot(object = ABC_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = ABC_harmony, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
plot_grid(p1,p2)

##### Do UMAP and clustering:
ABC_harmony <- ABC_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()

# Finally, so same UMAP plots of integrated datasets as above:
# pbmc_harmony <- SetIdent(pbmc_harmony,value = "orig.ident")
DimPlot(ABC_harmony,reduction = "umap", group.by = "orig.ident") + plot_annotation(title = "A1B1C1, after integration (Harmony)")
DimPlot(ABC_harmony,reduction = "umap", group.by = "orig.ident", split.by = "orig.ident") + plot_annotation(title = "A1B1C1, after integration (Harmony)")
DimPlot(ABC_harmony,reduction = "umap", group.by = "ident", split.by = "orig.ident", label = TRUE, label.box = FALSE) + plot_annotation(title = "A1B1C1, after integration (Harmony)")

##### Finally, let’s take a look at the cluster content:
# We can now calculate the number of cells in each cluster that came for each samplet:
count_table <- table(ABC_harmony@meta.data$seurat_clusters, ABC_harmony@meta.data$orig.ident)
cluster_size <- rowSums(count_table)
count_df <- cbind(test, cluster_size = cluster_size)
rownames(count_df) <- paste0("cluster", rownames(count_df))
knitr::kable(count_df)

# Let’s plot the distribution among clusters using our custom function:
plot_integrated_clusters(ABC_harmony) 
# 
# #### Puis sur les lames?
# SpatialDimPlot(ABC_harmony, ncol = 2)
SpatialDimPlot(ABC_harmony, ncol = 3)
SpatialPlot(ABC_harmony, group.by = "ident")
SpatialDimPlot(ABC_harmony, ncol = 3, images = "sliceA1",
               cells.highlight = CellsByIdentities(object = ABC_harmony, idents = c(0, 2, 3, 4, 6, 7, 9, 11, 12)),
               facet.highlight = TRUE)
SpatialDimPlot(ABC_harmony, ncol = 3, images = "sliceB1",
               cells.highlight = CellsByIdentities(object = ABC_harmony, idents = c(0, 1, 2, 5, 6, 7, 11, 15, 16, 19)),
               facet.highlight = TRUE)

SpatialDimPlot(ABC_harmony, ncol = 3,
               cells.highlight = CellsByIdentities(object = ABC_harmony, idents = c(12)))
# SpatialDimPlot(ABC_harmony, ncol = 3,
#                cells.highlight = CellsByIdentities(object = ABC_harmony, idents = c(1, 2)), facet.highlight = TRUE)

# saveRDS
saveRDS(ABC_harmony, file = "output/ABC_object/20220622_ABC_object_method2.rds")


# 
# # splitData <- SplitObject(sample_seurat, split.by = "orig.ident")
# # sampleA1b <- splitData$sampleA1
# # SpatialDimPlot(sampleA1b, cells.highlight = CellsByIdentities(object = sampleA1b, idents = c(1, 2, 3, 4, 5, 6, 7 ,8)), facet.highlight = TRUE, ncol = 3)
# 
###### Export Corrected Projections
# split the object
split.data <- SplitObject(ABC_harmony, split.by = "orig.ident")

# edit the barcodes into a format that is compatible with the Loupe Browser
A1.barcode <- rownames(Embeddings(object = split.data$sampleA1, reduction = "umap"))
# A1.barcode <- gsub("A1_", "", A1.barcode)
A1.barcode <- gsub('.{2}$', '', A1.barcode)
A1.barcode <- gsub('.{2}$', '', A1.barcode)
A1.barcode <- paste(A1.barcode,"-1", sep="")
A1.proj <- Embeddings(object = split.data$sampleA1, reduction = "umap")
UMAP.A1 <- cbind("Barcode" = A1.barcode, A1.proj)

B1.barcode <- rownames(Embeddings(object = split.data$sampleB1, reduction = "umap"))
# B1.barcode <- gsub("B1_","",B1.barcode)
B1.barcode <- gsub('.{2}$','',B1.barcode)
B1.barcode <- gsub('.{2}$','',B1.barcode)
B1.barcode <- paste(B1.barcode,"-2", sep="")
B1.proj<- Embeddings(object = split.data$sampleB1, reduction = "umap")
UMAP.B1 <- cbind("Barcode" = B1.barcode, B1.proj)

C1.barcode <- rownames(Embeddings(object = split.data$sampleC1, reduction = "umap"))
# C1.barcode <- gsub("C1_","",C1.barcode)
C1.barcode <- gsub('.{2}$','',C1.barcode)
C1.barcode <- gsub('.{2}$','',C1.barcode)
C1.barcode <- paste(C1.barcode,"-3", sep="")
C1.proj<- Embeddings(object = split.data$sampleC1, reduction = "umap")
UMAP.C1 <- cbind("Barcode" = C1.barcode, C1.proj)

# merge the two samples back into the same object, and export into a CSV
split.umap <- rbind(UMAP.A1, UMAP.B1, UMAP.C1)
write.table(split.umap, file = "output/correctedABC/20220622_method2_corrected_umap_ABC.csv", sep = ",", quote = F, row.names = F, col.names = T)

###### 7. Export Clusters
clusters.A1 = Idents(split.data$sampleA1)
clusters.A1.data <- cbind("Barcode" = A1.barcode, data.frame("clusters" = clusters.A1))

clusters.B1 = Idents(split.data$sampleB1)
clusters.B1.data <- cbind("Barcode" = B1.barcode, data.frame("clusters" = clusters.B1))

clusters.C1 = Idents(split.data$sampleC1)
clusters.C1.data <- cbind("Barcode" = C1.barcode, data.frame("clusters" = clusters.C1))

split.cluster <- rbind(clusters.A1.data, clusters.B1.data, clusters.C1.data)
write.table(split.cluster, file="output/correctedABC/20220622_method2_corrected_clusters_ABC.csv", sep = ",", quote = F, row.names = F, col.names = T)