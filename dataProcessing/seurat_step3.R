#!/usr/bin/Rscript

library(Seurat)
library(clustree)
library(ggplot2)
library(stringr)

file <- "./output/seu_hVoWadj_rmPal_reg_S2.rds"
outName <- "221002_hVoWadj_rmPal_Wreg"
seu.integrated.obj <- readRDS(file)

#set finalized number of dimensions 
test_dims <- c(50,45,40,35,30)

for (dimz in test_dims) {

    seu.test <- FindNeighbors(object = seu.integrated.obj, dims = 1:dimz)

    seu.test <- FindClusters(object = seu.test, algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1))) ###can MODIFY resolution vector as desired###

    p <- clustree::clustree(seu.test, prefix = "integrated_snn_res.") + ggtitle(paste0("The number of dims used:", dimz, sep = ""))

    outfile <- paste("./output/",dimz ,"_clustree_", outName,".png", sep = "") 
    png(file = outfile, height = 1100, width = 2000)

    print(p)
    dev.off()

}

