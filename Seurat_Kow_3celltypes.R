#single experiment script to find DE genes within
#another example from satija lab
# adapted by L.V.Bystrykh
#from https://satijalab.org/seurat/immune_alignment.html
#later tips in https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
install.packages("Seurat")
#library(Matrix)
library(dplyr)
library(Seurat)
#try the same raw data
#setwd("/home/leonid/Documents/courses/GED_Milan_2021/3.RNAseq/SingleCell_R/scater")
setwd("~/Documents/review_SC_python")
f='GSE59114_C57BL6_Y_all.csv'
pdt.data <- read.csv(f,  header=TRUE, row.names=1, sep="\t")
#adjust column names, it will be used for definition of cells
colnames(pdt.data)<-sub("young_","",colnames(pdt.data))
#Kowalchyk data are already log-transformed, restore read counts
pdt.data<-round(2^(pdt.data)-1,0)
pdt <- CreateSeuratObject(counts = pdt.data, 
                          project = "Young", min.cells = 5, min.features=3000)
pdt<-AddMetaData(object=pdt, metadata=pdt@meta.data$orig.ident, col.name="cell_types")
#you can check genes and reads per cell as scatter plot
plot(pdt@meta.data$nCount_RNA, pdt@meta.data$nFeature_RNA,
     xlab="counts",
     ylab="features", 
     col=pdt@meta.data$orig.ident,
     pch=3, cex=0.3,
     main=paste(f, "\ncounts and features"))
# check how data are distributed
hist(as.matrix(log2(pdt.data+1)), main="Data")
#pdt@meta.data$cell_types<-c(rep("LT",167), rep("ST",330-167),rep("MPP", 498-330))
#scRNA-seq QC metric.
mito.genes <- grep(pattern = "^mt-", x = rownames(pdt@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(pdt@assays$RNA@data[mito.genes, ])/Matrix::colSums(pdt@assays$RNA@data)
ERCC.genes<-grep(pattern="^ERCC-",x=rownames(pdt@assays$RNA@data), value=TRUE)
percent.ERCC <- Matrix::colSums(pdt@assays$RNA@data[ERCC.genes, ])/Matrix::colSums(pdt@assays$RNA@data)
pdt <- AddMetaData(object = pdt, metadata = percent.mito, col.name = "percent.mito")
pdt <- AddMetaData(object = pdt, metadata = percent.ERCC, col.name = "percent.ERCC")
#show data
VlnPlot(pdt, features = c("nFeature_RNA", "nCount_RNA", "percent.ERCC","percent.mito"), ncol = 2)
#the same as self-made plot above, now as a package function
FeatureScatter(pdt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#normalize if not yet normalized SCTransform does normalization to the spiking genes
#this option will not work with Jackstraw method below, skip if you want more of DE and clouds
#pdt <- SCTransform(pdt, 
 #                 # vars.to.regress = "percent.ERCC", 
  #                 verbose = FALSE)
#alternatively all sets can be uniformly normalised as follows:
pdt <- NormalizeData(pdt, normalization.method = "LogNormalize", scale.factor = 10000)
#check transformed data
as.matrix(pdt@assays$RNA@data)[1:5,1:20]
plot1 <- FeatureScatter(pdt, feature1 = "nCount_RNA", feature2 = "percent.mito")
#if SCTransform was used, probably not used
#plot2 <- FeatureScatter(pdt, feature1 = "nCount_SCT", feature2 = "nFeature_SCT")
#if NormalizeData was used
plot2 <- FeatureScatter(pdt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#find variable genes, #you can try parameters, but it is not always effective
pdt <- FindVariableFeatures(object = pdt,
                            mean.cutoff = c(3,Inf),
                            dispersion.cutoff = c(3,Inf),
                            nfeatures=3000) 
length(VariableFeatures(pdt))
top10 <- head(VariableFeatures(pdt), 10)
plot1<-VariableFeaturePlot(pdt)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
#check if the plot really shows what you want by cutoffs, probably not

#scale data, it scales all genes to average=0, sd=1
pdt <- ScaleData(object = pdt)
#run PCA on variable genes
pdt <- RunPCA(object = pdt, 
              pc.genes = pdt@var.genes, 
              do.print = TRUE, pcs.print = 1:3)
#results in pdt@reductions$pca
# Examine and vine dots as 
#visualize PCA results
PCAPlot(object = pdt)#, dim.1 = 1, dim.2 = 2)
#you can define dots as PCAPlot and extract dots coords
#plot(dots$data[,1],dots$data[,2], xlab="comp.1", ylab="comp.2", main="PCA coordinates")

#optional, shows most influential genes for selected component
PCHeatmap(object = pdt, dims = 1, cells = 90, #adjust cells use to real numbers!
          balanced = TRUE)#, label.columns = FALSE)

#find significant components. It seems to work only with normalization by NormalizeData()
pdt <- JackStraw(object = pdt, num.replicate = 100, dims=10)#, display.progress = FALSE)
pdt <- ScoreJackStraw(object= pdt, dims=1:10)
#enrichment for low Pval genes plot
JackStrawPlot(object = pdt, dims=1:10)

#refine cells clustering, play with resolution and perplexity
pdt <- FindNeighbors(object = pdt)
pdt <- FindClusters(object = pdt, reduction.type = "pca", dims.use = 1:5, 
                     resolution = 0.3, print.output = 0, algorithm= 3, #change for big sets
                    save.SNN = TRUE)
#t-SNE plot
pdt <- RunTSNE(object = pdt, dims.use = 1:5, do.fast = TRUE, perplexity=15, ncomponents=2)
#coords are in reducedDims
clustrs<-TSNEPlot(object = pdt)#, group_by="seurat_clusters")
clustrs$data
clustrs
dots<-TSNEPlot(object = pdt, group.by="cell_types")
dots$data
dots
dots+clustrs
pdt<-RunUMAP(object=pdt, dims=1:5)
UMAPPlot(object=pdt, reduction="umap")
dots<-UMAPPlot(object=pdt, reduction="umap", group.by="cell_types")
write.table(dots$data, "Kow_Seurat_UMAP.tsv", sep="\t")
#cluster IDs in pdt$seurat_clusters or pdt@meta.data$seurat_clusters
head(pdt@meta.data,10)
#finding DE genes for selected cluster (you need at least two)
cluster1.markers <- FindMarkers(object = pdt, ident.1 = 0, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
# find markers distinguishing one cluster from another or group
cluster0.markers <- FindMarkers(object = pdt, 
                                ident.1 = 1, #ST
                                ident.2 = c(2,3), #MPP
                                min.pct = 0.25)
print(x = head(x = cluster0.markers, n = 5))
# all positive markers for all groups
pdt.markers <- FindAllMarkers(object = pdt, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pdt.markers %>% group_by(cluster) #%>% top_n(10, avg_logFC)
#other test for significance=roc
#cluster01.markers <- FindMarkers(object = pdt, ident.1 = 0, ident.2 = 1,
 #                               test.use = "roc", only.pos = TRUE)
#top01<-cluster01.markers  %>% top_n(20, avg_logFC)
#top01<-subset(top01, p_val_adj<0.05)

#plot specific gene expressions in clouds
VlnPlot(object = pdt, features = c("Abca5", "Vwf","Selp"))
FeaturePlot(object = pdt, features = c("Abca5", "Selp","Mpo","Cd48"), cols = c("grey", "blue"), 
            reduction = "tsne")
#show heatmap for markers
#top10 <- pdt.markers %>% group_by(cluster) %>% top_n(10,  avg_logFC)#-p_val_adj)#
#DoHeatmap(object = pdt, features = top10$gene)
#change Idents parameters from clusters to "young" and "old"
Idents(pdt)<- pdt@meta.data$cell_types
LT_ST_markers <- FindMarkers(pdt, ident.1 = "LT", ident.2 = "ST", 
                            # group.by="age",
                            # subset.ident=set.combined@meta.data$age, 
                            verbose = FALSE)
print(x = head(x = LT_ST_markers, n = 5))
FeaturePlot(object = pdt, features = c("Gnb1", "Ly6a","Mpo","Cd48"), cols = c("grey", "blue"), 
            reduction = "umap")
All_markers<- FindAllMarkers(object = pdt, 
                             only.pos = TRUE, 
                             logfc.threshold = 1, #log2 fold change
                             min.pct = 0.2) #minimal fraction of cells with this gene
write.table(All_markers,"Kowalczyk_all_markers_F3000.tsv", sep="\t")
