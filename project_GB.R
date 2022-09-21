library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)

#library(R.utils)
#gunzip("/Users/garancehaefliger/Desktop/GB_project/data/GSE122466_Merged5347cells_RAW.csv.gz", remove = FALSE)

############# PRE-PROCESS THE 10X DATA #################################

data <- read.csv(file ="/Users/garancehaefliger/Desktop/GB_project/data/GSE122466_Merged5347cells_RAW.csv", sep=';')
rownames <- data[,1]
data <- data[,2:ncol(data)]
row.names(data) <- rownames

#data@meta.data

data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

data <- NormalizeData(data)
data <- ScaleData(data)
data <- FindVariableFeatures(data, selection.method = 'vst')

data <- RunPCA(data, features = VariableFeatures(data))
#DimPlot(data, reduction = "pca")

data <- RunHarmony(data, group.by.vars = 'orig.ident')
#DimPlot(data, reduction='harmony')

data <- FindNeighbors(data, dims = 1:10, reduction='harmony')
data <- FindClusters(data, graph.name = 'RNA_nn')

data <- RunTSNE(data, dims = 1:5, reduction="harmony")
DimPlot(data, reduction = "tsne")

data <- RunUMAP(data, dims = 1:10, min.dist = 0.3, n.neighbors = 30L, reduction = "harmony") # 3 embbeded dim ??
DimPlot(data, reduction = "umap")

# Rename cluster
#new.cluster.ids <- c("Progenitors 1", "Progenitors 2", "RGC 1", " Neuroblast-G2M", "Neuroblast-G1", "Amacrine", "Photoreceptors", "RGC 2", "U/RGC", "Progenitors 4", "Amacrine/Horizontal", "Progenitors 3", "12", "13")
#names(new.cluster.ids) <- levels(data)
#data <- RenameIdents(data, new.cluster.ids)
#DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



######### EXTRACT THE SUBSET #####################
subset_data <- subset(data, idents = c(0,1,9,11), invert = TRUE)

subset_data <- NormalizeData(subset_data)
subset_data <- ScaleData(subset_data)
subset_data <- FindVariableFeatures(subset_data, selection.method = 'vst')

subset_data <- RunPCA(subset_data, features = VariableFeatures(subset_data))
#DimPlot(subset_data, reduction = "pca")

subset_data <- RunHarmony(subset_data, group.by.vars = 'orig.ident')
#DimPlot(subset_data, reduction='harmony')

subset_data <- RunTSNE(subset_data, dims = 1:5, reduction="harmony")
DimPlot(subset_data, reduction = "tsne")

new.cluster.ids_ <- c( "RGC 1", "Neuroblast-G2M", "Neuroblast-G1", "Amacrine", "Photoreceptors", "RGC 2", "U/RGC", "Amacrine/Horizontal", "12", "13")
names(new.cluster.ids_) <- levels(subset_data)
subset_data <- RenameIdents(subset_data, new.cluster.ids_)

subset_data <- RunUMAP(subset_data, dims = 1:10, min.dist = 0.3, n.neighbors = 30L, reduction = "harmony") # 3 embbeded dim ??
DimPlot(subset_data, reduction = "umap", label = TRUE)

# find marker of yellow cluster 
cluster12.markers <- FindMarkers(subset_data, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 30)


cluster2.markers <- FindMarkers(subset_data, ident.1 ="RGC 1", min.pct = 0.25)
head(cluster2.markers, n = 30)

# Quality control of cluster 12 vs the rest
VlnPlot(subset_data, features = c("nFeature_RNA", "nCount_RNA"))
#cluster_12 <- subset(subset_data, idents = c(12))
#VlnPlot(cluster_12, features = c("nFeature_RNA", "nCount_RNA"))
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

######## SAVING DATA FOR PYTHON ##################

# assuming that you have some Seurat object called seurat_obj:
seurat_obj <- subset_data
# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='/Users/garancehaefliger/Desktop/GB_project/data/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
out_data_dir = ''
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(out_data_dir, '/Users/garancehaefliger/Desktop/GB_project/data/counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='/Users/garancehaefliger/Desktop/GB_project/data/pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='/Users/garancehaefliger/Desktop/GB_project/data/gene_names.csv',
  quote=F,row.names=F,col.names=F
)



