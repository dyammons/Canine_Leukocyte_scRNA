#!/usr/bin/Rscript

library(Seurat)
library(ggplot2)

file <- "./output/seu_hVoWadj_rmPal_reg_S2.rds"
seu.integrated.obj <- readRDS(file)

###MODIFY based on clustree DATA###
final.dims <- 50
final.res <- 1.3
outName <- "221004_hVoWadj_rmPal_reg"
stashID <- "clusterID"

#perform clustering: Neighbor finding and resolution parameters
seu.integrated.obj <- FindNeighbors(object = seu.integrated.obj, dims = 1:final.dims)
seu.integrated.obj <- FindClusters(object = seu.integrated.obj, algorithm = 3, resolution = final.res)

#choose appropriate clustering resolution
res <- paste("integrated_snn_res.", final.res,sep = "") 
seu.integrated.obj <- SetIdent(object = seu.integrated.obj, value = res)

# Run UMAP embedding
seu.integrated.obj <- RunUMAP(seu.integrated.obj,
  dims = 1:final.dims,
  min.dist = 0.6,
  n.neighbors = 75
)

#save cluster number as "clusterID"
seu.integrated.obj[[stashID]] <- seu.integrated.obj@active.ident

#build hierarchical clustering based on clustering results
DefaultAssay(seu.integrated.obj) <- "RNA"
seu.integrated.obj <- NormalizeData(seu.integrated.obj)
seu.integrated.obj <- BuildClusterTree(seu.integrated.obj, assay = "RNA", dims = 1:final.dims)

#export UMAP colored by sample ID and cluster
outfile <- paste("./output/", outName, "_res", final.res, "_dims", final.dims, "_cluster_S3.png", sep = "") 
p <- DimPlot(seu.integrated.obj, label = TRUE, reduction = "umap", group.by = stashID)
ggsave(outfile)

outfile <- paste("./output/", outName, "_res", final.res, "_dims", final.dims, "_sample_S3.png", sep = "") 
p <- DimPlot(seu.integrated.obj, reduction = "umap", group.by = "orig.ident")
ggsave(outfile)

features <- c("PTPRC", "CD3E", "CD8A", "GZMA", 
                "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                "CD4", "MS4A1", "PPBP","HBM")

#visulize UMAP featPlots
outfile <- paste("./output/", outName,"_featPlotDefault_S3.png", sep = "")
p <- FeaturePlot(seu.integrated.obj,features = features)
ggsave(outfile, width = 12, height = 8)

#stash identy by healthy vs OS & store as cellSource
Idents(seu.integrated.obj) <- "orig.ident"
seu.integrated.obj$cellSource <- ifelse(grepl("ealthy", seu.integrated.obj@meta.data$orig.ident), "Healthy", "Osteosarcoma")

#save seurat object as rds
outfile <- paste("./output/", outName, "_res", final.res, "_dims", final.dims, "_S3.rds", sep = "") 
saveRDS(seu.integrated.obj, file = outfile)

