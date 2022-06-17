setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(harmony)

# https://www.10xgenomics.com/resources/analysis-guides/correcting-batch-effects-in-visium-data

##### 1. Create Directories, Download data

##### 2. Load and Combine Data Sets 
sampleA1_data <- Read10X(data.dir = "outputBN/sampleA1/outs/filtered_feature_bc_matrix")
sampleA1 <- CreateSeuratObject(counts = sampleA1_data, project = "A1")
sampleA1

sampleB1_data <- Read10X(data.dir = "outputBN/sampleB1/outs/filtered_feature_bc_matrix")
sampleB1 <- CreateSeuratObject(counts = sampleB1_data, project = "B1")
sampleB1

##### 3. Merge Objects
A1B1 <- merge(sampleA1, y = sampleB1, add.cell.ids = c("A1", "B1"), project = "A1B1")
A1B1
# colnames(A1B1)
head(colnames(A1B1))
tail(colnames(A1B1))

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
write.table(corrected.umap, file = "outputBN/corrected/corrected_umap_A1B1.csv", sep = ",", quote = F, row.names = F, col.names = T)

###### 7. Export Clusters
clusters.A1 = Idents(corrected.data$A1)
clusters.A1.data <- cbind("Barcode" = A1.barcode, data.frame("clusters" = clusters.A1))

clusters.B1 = Idents(corrected.data$B1)
clusters.B1.data <- cbind("Barcode" = B1.barcode, data.frame("clusters" = clusters.B1))

corrected.cluster <- rbind(clusters.A1.data, clusters.B1.data)
write.table(corrected.cluster, file="outputBN/corrected/corrected_clusters_A1B1.csv", sep = ",", quote = F, row.names = F, col.names = T)
