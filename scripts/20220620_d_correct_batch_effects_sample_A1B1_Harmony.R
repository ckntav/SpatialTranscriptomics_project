setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(harmony)

# https://www.10xgenomics.com/resources/analysishttps://www.singlecellcourse.org/scrna-seq-dataset-integration.html#harmony-3-vs-5-10k-pbmc
##### 1. Create Directories, Download data

##### 2. Load and Combine Data Sets 
sampleA1_data <- Read10X(data.dir = "output/sampleA1/outs/filtered_feature_bc_matrix")
sampleA1 <- CreateSeuratObject(counts = sampleA1_data, project = "A1")
# sampleA1 <-  SCTransform(sampleA1, do.scale = TRUE, do.center = TRUE)
sampleA1

sampleB1_data <- Read10X(data.dir = "output/sampleB1/outs/filtered_feature_bc_matrix")
sampleB1 <- CreateSeuratObject(counts = sampleB1_data, project = "B1")
# sampleB1 <- SCTransform(sampleB1, do.scale = TRUE, do.center = TRUE)
sampleB1

##### 3. Merge Objects
A1B1 <- merge(sampleA1, y = sampleB1, add.cell.ids = c("A1", "B1"), project = "A1B1")
A1B1
# colnames(A1B1)
head(colnames(A1B1))
tail(colnames(A1B1))

#
A1B1 <- NormalizeData(A1B1, verbose = F)
A1B1 <- FindVariableFeatures(A1B1, selection.method = "vst", nfeatures = 2000, verbose = F)
A1B1 <- ScaleData(A1B1, verbose = F)
A1B1 <- RunPCA(A1B1, npcs = 30, verbose = F)
A1B1 <- RunUMAP(A1B1, reduction = "pca", dims = 1:30, verbose = F)
DimPlot(A1B1, reduction = "umap") + plot_annotation(title = "A1B1, before integration")

#s
A1B1 <- A1B1 %>% RunHarmony("orig.ident", plot_convergence = T)

#
A1B1_embeddings <- Embeddings(A1B1, 'harmony')
A1B1_embeddings[1:5, 1:5]

p1 <- DimPlot(object = A1B1, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = A1B1, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
plot_grid(p1,p2)

#
A1B1 <- A1B1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()

# A1B1 <- SetIdent(A1B1, value = "orig.ident")
DimPlot(A1B1, reduction = "umap") + plot_annotation(title = "A1B1, after integration (Harmony)")
DimPlot(A1B1, reduction = "umap", split.by = "orig.ident") + plot_annotation(title = "A1B1, after integration (Harmony)")

DimPlot(A1B1, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()
DimPlot(A1B1, reduction = "umap", group.by = "orig.ident", pt.size = .1)


# A1B1 <- SetIdent(A1B1,value = "seurat_clusters")
DimPlot(A1B1, label = T) + NoLegend()

# A1B1 <- SetIdent(A1B1, value = "orig.ident")
DimPlot(A1B1) + NoLegend()






#
VariableFeatures(A1B1) <- c(VariableFeatures(sampleA1), VariableFeatures(sampleB1))
A1B1 <- RunPCA(A1B1, verbose = FALSE)
DimPlot(A1B1, group.by = "orig.ident")
# A1B1 <- FindNeighbors(A1B1, dims = 1:30)
# A1B1 <- FindClusters(A1B1, verbose = FALSE)
# A1B1 <- RunUMAP(A1B1, dims = 1:30)
# DimPlot(A1B1, reduction = "umap", group.by = c("ident", "orig.ident"))

#
A1B1 <- RunHarmony(A1B1, group.by.vars = "orig.ident", assay.use = "SCT")
A1B1 <- RunUMAP(A1B1, reduction = "harmony", dims = 1:30)
A1B1 <- FindNeighbors(A1B1, reduction = "harmony", dims = 1:30) %>% FindClusters() %>% identity()
# "orig.ident" = original identity
DimPlot(A1B1, group.by = "orig.ident")
# "ident" = identity, which are clusters
DimPlot(A1B1, group.by = "ident", split.by = 'orig.ident')
DimPlot(A1B1, group.by = "ident")
DimPlot(A1B1, group.by = c("ident", "orig.ident"))

plot_integrated_clusters(pbmc_harmony)


A1B1 <- FindVariableFeatures(A1B1)

A1B1 <- RunPCA(A1B1, assay= "SCT", verbose = TRUE, npcs = 50)
DimPlot(A1B1, group.by = "orig.ident")



##### 4. Visualize
# First, visualize the data before running batch effect correction.
# Within this command, we will also be normalizing, scaling the data,
# and calculating gene and feature variance which will be used to run a PCA and UMAP.
A1B1 <- NormalizeData(A1B1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = TRUE)
DimPlot(A1B1, group.by = "orig.ident")
A1B1 <- RunUMAP(A1B1, dims = 1:30)
DimPlot(A1B1, group.by = "orig.ident")

##### 5. Run Harmony
A1B1 <- RunHarmony(A1B1, group.by.vars = "orig.ident")
A1B1 <- RunUMAP(A1B1, reduction = "harmony", dims = 1:30)
A1B1 <- FindNeighbors(A1B1, reduction = "harmony", dims = 1:30) %>% FindClusters()
# "orig.ident" = original identity
DimPlot(A1B1, group.by = "orig.ident")
# "ident" = identity, which are clusters
DimPlot(A1B1, group.by = "ident", split.by = 'orig.ident')
DimPlot(A1B1, group.by = "ident")

###### 6. Export Corrected Projections
# split the object
corrected.data <- SplitObject(A1B1, split.by = "orig.ident")

# edit the barcodes into a format that is compatible with the Loupe Browser
A1.barcode <- rownames(Embeddings(object = corrected.data$A1, reduction = "umap"))
A1.barcode <- gsub("A1_", "", A1.barcode)
A1.barcode <- gsub('.{2}$', '', A1.barcode)
A1.barcode <- paste(A1.barcode,"-1", sep="")
A1.proj <- Embeddings(object = corrected.data$A1, reduction = "umap")
UMAP.A1 <- cbind("Barcode" = A1.barcode, A1.proj)

B1.barcode <- rownames(Embeddings(object = corrected.data$B1, reduction = "umap"))
B1.barcode <- gsub("B1_","",B1.barcode)
B1.barcode <- gsub('.{2}$','',B1.barcode)
B1.barcode <- paste(B1.barcode,"-2", sep="")
B1.proj<- Embeddings(object = corrected.data$B1, reduction = "umap")
UMAP.B1 <- cbind("Barcode" = B1.barcode, B1.proj)

# merge the two samples back into the same object, and export into a CSV
corrected.umap <- rbind(UMAP.A1, UMAP.B1)
write.table(corrected.umap, file = "output/correctedA1B1/corrected_umap_A1B1.csv", sep = ",", quote = F, row.names = F, col.names = T)

###### 7. Export Clusters
clusters.A1 = Idents(corrected.data$A1)
clusters.A1.data <- cbind("Barcode" = A1.barcode, data.frame("clusters" = clusters.A1))

clusters.B1 = Idents(corrected.data$B1)
clusters.B1.data <- cbind("Barcode" = B1.barcode, data.frame("clusters" = clusters.B1))

corrected.cluster <- rbind(clusters.A1.data, clusters.B1.data)
write.table(corrected.cluster, file="output/correctedA1B1/corrected_clusters_A1B1.csv", sep = ",", quote = F, row.names = F, col.names = T)
