#!/usr/bin/Rscript

##################################################################################################################
#       ### Step 1: process 10X data, subset data based on QC and identify putative cell dublets] ###            #
#                                                                                                                #
# Required file structure: create a directory that has 2 sub-directories 1) "input" 2) "output"                  #
#  >> Put this script in the in the base directory -- this sets working directory as needed for proper run       #
#  >> Put files to be analyzed into "input" folder - each folder in the input directory should have              #
#     unique sample name that contains the three 10X output files                                                #
#  >> The unique file name will be used to name the files that are generated throughout the pipeline             #
#                                                                                                                #
# Step 1 will:                                                                                                   #
#  >> loop through all input directories, convert to seurat object, then subset based on QC values               #
#  >> various figures will be generated include QC plots                                                         #
#     >> be sure to look at QC figures to ensure the parameters are set appropriately (set permissively!)        #
#     >> may need to modify the QC cutoffs and re-run                                                            #
#  >> putative doublets will be identified using doubletFinder then removed                                      #
##################################################################################################################

library(Seurat)
library(stringr)
library(DoubletFinder)
library(tidyverse)
######################## set pwds ####################

din <- "inputIntrons"
dout <- "output_S1_wPalpct_introns"
outName <- "221005_introns_wPalScore"


######################## set qc cut offs #######################

#if testQC value == TRUE then the script will only generate the violin plots for each input
testQC <- FALSE

#run with testQC == TRUE to get an idea of how to set the cut offs
nFeature_RNA_high <- 3500 
nFeature_RNA_low <- 200
nCount_RNA_high <- 15000 
nCount_RNA_low <- 100
percent.mt_high <- 20

#after identifying putative doublets you can either choose to remove them (removeDubs == TRUE) or save the metadata tag and leave them in the dataset (removeDubs == FALSE)
removeDubs <- TRUE

#after identifying putative doublets you can either choose to remove them (removeDubs == TRUE) or save the metadata tag and leave them in the dataset (removeDubs == FALSE)
removeRBC_pal <- TRUE
#################################################################

if (testQC == FALSE){
    print(paste0("The QC parameters are: nFeature_RNA < ", nFeature_RNA_high, " & nFeature_RNA > ", nFeature_RNA_low, " & percent.mt < ", percent.mt_high, " & nCount_RNA < ", nCount_RNA_high," & nCount_RNA > ", nCount_RNA_low, sep = ""))
    }

fpath <-  paste("./", din,"/", sep = "") 
files <- list.files(path = fpath, pattern=NULL, all.files=FALSE,
    full.names=FALSE)

#load in pal gene set
#pal_1 <- FindMarkers(seu.integrated.obj, ident.1 = 19, min.pct = 0.50)
#pal_2 <- FindMarkers(seu.integrated.obj, ident.1 = 24, min.pct = 0.50)

pal_1 <- read.csv("./output/pal_1.csv", row.names = 1)
pal_2 <- read.csv("./output/pal_2.csv", row.names = 1)
pal_1 <- pal_1[pal_1$avg_log2FC > 0, ]
pal_2 <- pal_2[pal_2$avg_log2FC > 0, ]

pal_feats <- rownames(pal_1)[rownames(pal_1) %in% rownames(pal_2)[!grepl("^MT-|^RPS|^RPL",rownames(pal_2))]]

df.list <- list()
for (infile in files) {
  
  #set up df for export
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df) <- c("Initial","Filtered","Platelet_rbc_rm","Singlets")
    
  #set import path
  pwd <- paste("./", din,"/", infile, "/", sep = "") 
  
  #read 10X data
  indata <- Read10X(pwd)
    
  #create seurat object
  seu_obj <- CreateSeuratObject(counts = indata,
                     project = infile,
                     min.cells = 3,
                     min.features = 200)
  
  
  #Add mitochondrial QC data to seurat metadata
  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
  seu_obj[["percent.hbm"]] <- PercentageFeatureSet(seu_obj, pattern = "HBM")
  seu_obj[["percent.ppbp"]] <- PercentageFeatureSet(seu_obj, pattern = "PPBP")
    
  #visualize QC metrics as a violin plot
  outfile <- paste("./",dout,"/", infile,"_QC_S1.png", sep = "") 
  p <- VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  ggsave(outfile)

  if (testQC == TRUE){
      next
  }
    
  df[1,1] <- dim(seu_obj)[2]

  #set QC cutoffs based on above plots ### MODIFY VALUES based on your data!! ###
  seu_obj <- subset(seu_obj,
                    subset = nFeature_RNA < nFeature_RNA_high & nFeature_RNA > nFeature_RNA_low &
                    percent.mt < percent.mt_high & nCount_RNA < nCount_RNA_high & nCount_RNA > nCount_RNA_low
  )
  
  #next steps normalize, scale, and run UMAP
  seu_obj <- NormalizeData(seu_obj,
                    normalization.method = "LogNormalize", 
                    scale.factor = 10000) ###change method to scran???

  seu_obj <- FindVariableFeatures(seu_obj,
                    selection.method = "vst", 
                    nfeatures = 2500) #can change number of feats used
  
  all.genes <- rownames(seu_obj)
  seu_obj <- ScaleData(seu_obj, features = all.genes) ###this step is optional! & possibly not recommended
  
  seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
  
  seu_obj <- FindNeighbors(seu_obj, 
                    dims = 1:10
                    ) #can change dims
  
  seu_obj <- FindClusters(seu_obj, 
                    resolution = 0.1
                    ) #can change resolution
  
  seu_obj <- RunUMAP(seu_obj, 
                    dims = 1:15
                    ) #can change dims
  
  features <- c("PTPRC", "CD3E", "CD8A", "GZMA", 
                "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                "CD4", "MS4A1", "PPBP","HBM")

  #visulize UMAP featPlots
  outfile <- paste("./",dout,"/", infile,"_featPlotDefault_S1.png", sep = "")
  p <- FeaturePlot(seu_obj,features = features)
  ggsave(outfile, width = 12, height = 8)

  #visulize UMAP
  outfile <- paste("./",dout,"/", infile,"_uMAP_S1.png", sep = "") 
  p <- DimPlot(seu_obj, reduction = "umap")
  ggsave(outfile)
    
  #stash initial cluster IDs
  seu_obj[["int.clusID"]] <- Idents(object = seu_obj)
    
  df[1,2] <- dim(seu_obj)[2]

  if (removeRBC_pal == TRUE){
  #find the rbc cluster
  if(length(AverageExpression(seu_obj, features = "HBM")) != 0){
  rbc.df <- as.data.frame(AverageExpression(seu_obj, features = "HBM"), header = TRUE)
  rbc.clus <- str_split(as.character(colnames(rbc.df)[max.col(rbc.df,ties.method="first")]),"[.]")[[1]][2]
    
  sus.rbc.size <- as.numeric(length(WhichCells(seu_obj, expression = HBM > 0)))
  sus.rbc.clus.size <- as.numeric(length(WhichCells(seu_obj, idents = rbc.clus )))

  sus.rbc.comp.pct <- sus.rbc.size/sus.rbc.clus.size
  rbc.clus <- ifelse(sus.rbc.comp.pct > 0.8, rbc.clus, "NULL")
  print(rbc.clus)
  } else{rbc.clus <- "NULL"}
    
    
  #find the palelet cluster
  if(length(AverageExpression(seu_obj, features = "PPBP")) != 0){
  pal.df <- as.data.frame(AverageExpression(seu_obj, features = "PPBP"), header = TRUE)
  pal.clus <- str_split(as.character(colnames(pal.df)[max.col(pal.df,ties.method="first")]),"[.]")[[1]][2]

  sus.pal.size <- as.numeric(length(WhichCells(seu_obj, expression = PPBP > 0)))
  sus.pal.clus.size <- as.numeric(length(WhichCells(seu_obj, idents = pal.clus )))

  sus.pal.comp.pct <- sus.pal.size/sus.pal.clus.size
  pal.clus <- ifelse(sus.pal.comp.pct > 0.8, pal.clus, "NULL")
  print(pal.clus)
  } else{pal.clus <- "NULL"}

  if(rbc.clus != "NULL" | pal.clus != "NULL"){
  seu_obj <- subset(seu_obj,
                  subset = int.clusID != pal.clus & int.clusID != rbc.clus)
  } else{print("No platelets of rbcs detected in sample!")}
  }
  df[1,3] <- dim(seu_obj)[2]
  
  data <- lapply(pal_feats, function(x){PercentageFeatureSet(seu_obj, pattern = x)})
  data <- as.data.frame(bind_cols(data))
  seu_obj[["percent.pal"]] <- rowSums(data)
    
  #next steps complete doublet identification using DoubletFinder
  sweep.res.list_seu_obj <- paramSweep_v3(seu_obj, PCs = 1:10, sct = FALSE)
  sweep.stats_seu_obj <- summarizeSweep(sweep.res.list_seu_obj, GT = FALSE)
  bcmvn_seu_obj <- find.pK(sweep.stats_seu_obj)
  pk_val <- as.numeric(as.character(bcmvn_seu_obj$pK[bcmvn_seu_obj$BCmetric == max(bcmvn_seu_obj$BCmetric)])) 
  annotations <- seu_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  #determine expected doublet rate -- this formula assumes 0.5% doublet per 1000 cells (lower than 10x recommendation to account for homoypic doublets and cells removed during QC)
  expceted_dub_rate <- dim(seu_obj)[2]/1000*0.5/100
  
  nExp_poi <- round(expceted_dub_rate*length(seu_obj@meta.data$orig.ident))  ### MODIFY VALUE; should be percentage expected to be doublets based on 10X expected values ###
  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  pN_value <- 0.25 ### MODIFY VALUE as needed ###
  
  pANN_value <- paste0("pANN_",pN_value,"_",pk_val,'_',nExp_poi, sep = "")
  seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:10, pN = pN_value, pK = pk_val, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  dfClass <- paste0("DF.classifications_",pN_value,"_",pk_val,'_',nExp_poi, sep = "")
  dfmeta <- paste0("seu_obj@meta.data$",dfClass, sep = "")
  seu_obj[["doublet"]] <- eval(parse(text = dfmeta))
  
  #export UMAP highlighting doublet vs singlet cells
  outfile <- paste("./",dout,"/", infile,"_DF_S1.png", sep = "") 
  DimPlot(seu_obj,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = dfClass )
  ggsave(outfile)
    
  if (removeDubs == TRUE){
  #remove putative doublets
  seu_obj <- subset(seu_obj,
                  subset = doublet == "Singlet"
  )}    

  #update user
  df[1,4] <- dim(seu_obj)[2]
  rownames(df) <- infile
  df.list[[which(infile == files)]] <- df
    
  #recluster the data with all unwanted cells removed
  seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
  
  seu_obj <- FindNeighbors(seu_obj, 
                    dims = 1:10
                    ) #can change dims
  
  seu_obj <- FindClusters(seu_obj, 
                    resolution = 0.1
                    ) #can change resolution
  
  seu_obj <- RunUMAP(seu_obj, 
                    dims = 1:15
                    ) #can change dims
  
  features <- c("PTPRC", "CD3E", "CD8A", "GZMA", 
                "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                "CD4", "MS4A1", "PPBP","HBM")

  outfile <- paste("./",dout,"/", infile,"_featPlotDefault_postRBC_pal_rm_S1.png", sep = "")
  p <- FeaturePlot(seu_obj,features = features)
  ggsave(outfile, width = 12, height = 8)
    
  #export final umap
  outfile <- paste("./",dout,"/", infile,"_uMAP_postRBC_pal_rm_S1.png", sep = "") 
  DimPlot(seu_obj, reduction = "umap")
  ggsave(outfile)
    
  #identify cycling cells
  seu_obj <- CellCycleScoring(
  object = seu_obj,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes, set.ident =
    FALSE
  )
  
  #export UMAP with cycling data
  outfile <- paste("./",dout,"/", infile,"_uMAP_cellCylce_postRBC_pal_rm_S1.png", sep = "") 
  DimPlot(seu_obj, reduction = "umap", group.by = "Phase")
  ggsave(outfile)

  #store cell cycle state
  seu_obj[["clusters"]] <- Idents(object = seu_obj)
  
  #export processed seurat object as an .RDS file
  outfile <- paste("./",dout,"/", infile,"_S1.rds", sep = "") 
  saveRDS(seu_obj, file = outfile)
}

cellCounts <- do.call(rbind, df.list)
outfile <- paste("./",dout,"/", outName, "_cell_counts_S1.csv", sep = "")
write.csv(cellCounts, file = outfile)
