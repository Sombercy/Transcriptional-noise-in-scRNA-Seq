library(dplyr)
library(Seurat)
library(Matrix)


pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- SCTransform(pbmc)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))


all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")



genes.ordered<- order(pbmc@assays$SCT@meta.features$sct.gmean, decreasing = TRUE)
bins <-lapply(0:9, function(x){
  start=floor(x/10*length(genes.ordered))+1
  end=floor((x+1)/10*length(genes.ordered))
  genes.ordered[start:end]
})
bins<- bins[2:9]
bins<-lapply(bins, function(bin){
  bin[order(pbmc@assays$SCT@meta.features$sct.variance[bin])]
})


nv.genes<-lapply(bins, function(bin){
  bin[1:(length(bin)%/%10)]
})
nv.genes_unlisted<- unlist(nv.genes)
pbmc.nv<-pbmc[nv.genes_unlisted,]

pbmc.nv[,pbmc.nv@active.ident==5]@assays$RNA


noise<-lapply(levels(pbmc.nv@active.ident), function(l){
  center<-apply(sqrt(pbmc.nv[,pbmc.nv@active.ident==l]@assays$RNA@counts),1, mean)
  apply(sqrt(pbmc.nv[,pbmc.nv@active.ident==l]@assays$RNA@counts),2, function(x){
    sqrt(sum((x-center)^2))
  })
})

noise_data <- lapply(seq_along(noise), function(n){
  data.frame(noise=noise[[n]], cluster=n)
})

# rbind(noise_data[[1]], noise_data[[2]], noise_data[[3]])[1200,]
noise_data <- do.call(rbind, noise_data)
unique(noise_data$cluster)
noise_data

boxplot(noise ~ cluster, noise_data)

# noise_data$age <- rbinom(nrow(noise_data), 1, .5)
# library(ggplot2)
# ggplot(noise_data, aes(factor(cluster), noise, fill = factor(age))) +
#   geom_boxplot()
