#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

# #QC metrics <<- math to calculate qc (i.e. sequencing depth)
# qc.df <- read.csv("compiledQCmetrics.csv")
# mean(as.numeric(gsub(",","",qc.df[1:7,]$Mean.Reads.per.Cell)))

##############################
### BEGIN HEALTHY ANALYSIS ###
##############################

###Analysis note: 
#our preliminary analysis indicated that platlet-cell doublets were confounding clustering, so we attempted to reduce the impact of platelets (pals) by calculating a gene signature that was determined by compelitng an intitial analysis with pals included - note the approach used is similar to just using percent PPBP instead of an entire list of pal assocaited features
#similar gene signatures can be obtained by running the custom 'load10x' funciton with 'removeRBC_pal' set to FALSE then running:

    #pal_1 <- FindMarkers(seu.integrated.obj, ident.1 = 19, min.pct = 0.50) # where ident.1 = a pal cluster
    #pal_2 <- FindMarkers(seu.integrated.obj, ident.1 = 24, min.pct = 0.50) # where ident.1 = a pal cluster (if there is a second)
    
    #pal_1 <- read.csv("./output/pal_1.csv", row.names = 1)
    #pal_2 <- read.csv("./output/pal_2.csv", row.names = 1)
    #pal_1 <- pal_1[pal_1$avg_log2FC > 0, ]
    #pal_2 <- pal_2[pal_2$avg_log2FC > 0, ]

    #pal_feats <- rownames(pal_1)[rownames(pal_1) %in% rownames(pal_2)[!grepl("^MT-|^RPS|^RPL",rownames(pal_2))]]

#to avoid unnesscissary work during a reproducible run, the pal gene signatire is provided as pal_feats and is incooperated into the analysis
pal_feats = c('TIMP1', 'NAA10', 'ENSCAFG00000037735', 'GP6', 'SEC11C', 'FTL', 'NRGN', 'ACOT7', 'VCL', 'RSU1', 'ITGB1', 'H3-3A', 'RABGAP1L', 'SELP', 'SH3GLB1', 'ACTB', 'ENSCAFG00000008221', 'TLN1', 'GSN', 'AMD1', 'TREM2', 'SH3BGRL2', 'MYH9', 'PLEK', 'ENSCAFG00000042554', 'RAP1B', 'ENSCAFG00000004260', 'NAP1L1', 'PPBP', 'RASA3', 'ITGA2B', 'EIF1', 'ACTG1', 'C9H17orf64', 'JMJD6', 'CCL14', 'GNG11', 'IGF2BP3', 'TBXAS1', 'VDAC3', 'MARCHF2', 'TPM4', 'TKT', 'FTH1.1', 'FERMT3', 'RTN3', 'PRKAR2B', 'SVIP', 'ENSCAFG00000030286', 'ADA', 'MYL9', 'TUBB1', 'TUBA1B', 'METTL7A', 'THBS1', 'SERF2', 'PIF1', 'B2M', 'GAS2L1', 'YWHAH', 'HPSE', 'ATG3', 'ENSCAFG00000015217', 'ITGA6','RGS18', 'SUB1', 'LGALS1', 'CFL1', 'BIN2', 'CAT', 'RGS10', 'MGST3', 'TMBIM6', 'PFN1', 'CD63', 'RALBP1', 'GNAS', 'SEPTIN7', 'TPT1', 'UBB', 'ATF4', 'BBLN', 'MTDH', 'ENSCAFG00000017655','FYB1', 'ENO1', 'GABARAP', 'SSR4', 'MSN', 'ENSCAFG00000011134', 'ENSCAFG00000046637', 'COX8A', 'DLA-64', 'CD47', 'VASP', 'DYNLRB1', 'DLA88', 'SMDT1', 'ATP5PF','ELOB', 'ENSCAFG00000029155', 'ARPC3', 'VPS28', 'LRRFIP1', 'SRP14', 'ABRACL', 'ENSCAFG00000043577', 'ENSCAFG00000042598')



### Process and analyze healthy samples first
load10x(din = "./healthyIntrons/", dout = "./output/s1/healthy/", outName = "221005_h7_introns", testQC = F, nFeature_RNA_high = 4500, nFeature_RNA_low = 200, percent.mt_high = 10, nCount_RNA_high = 20000, nCount_RNA_low = 500, pal_feats = pal_feats)

#integrate the data into one object
sctIntegrate(din = "./output/s1/healthy/", dout = "./output/s2/", outName = "221005_h7_regPal_wPalpct_introns", vars.to.regress = c("percent.mt", "percent.pal"), nfeatures = 2000)

#load in integrated obj
seu.obj <- readRDS(file = "./output/s2/221005_h7_regPal_wPalpct_introns_seu.integrated.obj_S2.rds")

#determine optimal cluster parameters
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "221005_h7_regPal_wPalpct_introns", test_dims = c(50,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#reduce and visulize the dataset
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "221005_h7_regPal_wPalpct_introns", final.dims = 45, final.res = 1.1, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#remove low quality clusters
seu.obj <- subset(seu.obj,
                  subset = 
                  clusterID !=  "18" & clusterID ==  "19" & clusterID ==  "28" & clusterID ==  "31") 

#complete independent reclustering on subset data
indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "hVoWadj_CLEAN_introns", preSub = T,
                      vars.to.regress = c("percent.mt","percent.pal")
                       )

#determine optimal cluster parameters
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "hVoWadj_CLEAN_introns", test_dims = c(50,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#reduce and visulize the final dataset
final.dims = 45
final.res = 1.6
min.dist = 0.6
n.neighbors = 75

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "hVoWadj_CLEAN_introns", final.dims = final.dims, final.res = final.res, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = min.dist, n.neighbors = n.neighbors, assay = "integrated", saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


######################
### BEGIN FIGURE 1 ###
######################


### Complete healthy analysis
seu.obj <- readRDS("./output/s3/h7_CLEAN_introns_res1.6_dims45_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")
colArray <- read.csv("./h7_majorID_wIntrons.csv")
outName <- "fig1"


### Extra: Supplemental data: Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj.healthy, groupBy = "clusterID", numOfFeats = 24, outName = "h7",
                      outDir = "./output/h7/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)
                     
                            
### Fig 1a: plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              cols = colArray$color,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol, nudge_x = c(rep(0,length(colArray$labCol)-1),0.75)) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_a_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig supp: transfer TRDC data from ROS data
seu.obj.ros.cyto <- readRDS(file = "./ros_output/HvO_mt_r7_res1.5_dims50_S3.rds")

seu.obj <- AddMetaData(seu.obj, metadata = seu.obj.ros.cyto@assays$RNA@data[rownames(seu.obj.ros.cyto@assays$RNA@data) == "TRDC",], col.name = "TRDC")
seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 

features = c("TRDC")
p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 1, ncol = 1, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp_TRDC_UMAP.png", sep = ""), width = 7, height = 7)


### Fig 1b: key feature plots
features <- c("CD3G","CD8A", "GZMA",
              "IL7R","CD4", "S100A12",
              "DLA-DRA","FLT3", "ANPEP",
              "MS4A1","JCHAIN", "IL5RA",
              "TOP2A","GATA3", "CD34")
colorz <- c("black", "#3267AD", "#3267AD",
            "#90D293", "#933C81", "#933C81",
            "black", "#CC1A36", "#F39300",
            "#645A9F", "#DBB9EC", "#AF615B",
            "#13B9AD", "#00366C", "#7E7E7E"
           )
fig1b <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = colorz, order = F, legJust = "top") 
ggsave(paste("./output/", outName, "/", outName, "_b_featPlots.png", sep = ""), width =9, height = 15)



### Fig 1c: key dot plot features
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./h7_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "colorID")

fig1c <- majorDot(seu.obj = seu.obj, groupBy = "colorID",
                  yAxis = c("Monocyte","DC","B cell","Neutrophil","CD4 T cell","CD8/NK cell","Eosinophil","Basophil","gd T cell","Cycling T cell","CD34+ unk"),
                  features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
                               "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
                               "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
                               "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_c_majorDot.png", sep = ""), width =8, height = 6)


### Fig 1d: umap by sample
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name",
              cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = FALSE,
              shuffle = TRUE
)
fig1d <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "_d_umap_bySample.png", sep = ""), width =7, height = 7)


### Fig 1e: stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "colorID") +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12)) + NoLegend() + scale_x_discrete(limits=rev(c("Monocyte","DC","B cell","Neutrophil","CD4 T cell","CD8/NK cell","Eosinophil","Basophil","gd T cell","Cycling T cell","CD34+ unk")),expand = c(0, 0))
ggsave(paste("./output/", outName, "/", outName, "_e_stackedBar.png", sep = ""), width =7, height = 5)


### Fig supp: singleR to use human reference for cell classification
singleR(seu.obj = seu.obj, outName = "h7_CLEAN_introns", clusters = "clusterID", outDir = "./output/fig1/singleR/")


### Fig supp: stacked bar graph to supplement umap by sample
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "clusterID") +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1, byrow =T))
ggsave(paste("./output/", outName, "/", outName, "_supp_stackedBar.png", sep = ""), width =8, height = 12)

#################################################
### END FIGURE 1  | BEGIN HvO DATA PROCESSING ###
#################################################

### Prepare healthy vs OS dataset
#load in 10x data and qc filter eeach sample
load10x(din = "./inputIntrons/", dout = "./output/s1/", outName = "221005_introns", testQC = F, nFeature_RNA_high = 4500, nFeature_RNA_low = 200, percent.mt_high = 10, nCount_RNA_high = 20000, nCount_RNA_low = 500, pal_feats = pal_feats)

#integrate the data into one object
sctIntegrate(din = "./output/s1/", dout = "./output/s2/", outName = "221005_hVoWadj_regPal_wPalpct_introns", vars.to.regress = c("percent.mt", "percent.pal"), nfeatures = 2000)

#load in integrated obj
seu.obj <- readRDS(file = "./output/s2/221005_hVoWadj_regPal_wPalpct_introns_seu.integrated.obj_S2.rds")

#determine optimal cluster parameters
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "221005_hVoWadj_regPal_wPalpct_introns", test_dims = c(50,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#reduce and visulize the dataset
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "hVoWadj_regPal_wPalpct", final.dims = 45, final.res = 1.1, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#remove low quality clusters
seu.obj <- subset(seu.obj,
                  subset = 
                  clusterID !=  "18" & clusterID ==  "19" & clusterID ==  "28" & clusterID ==  "31") 

#complete independent reclustering on subset data
indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "hVoWadj_CLEAN_introns", preSub = T,
                      vars.to.regress = c("percent.mt","percent.pal")
                       )

#determine optimal cluster parameters
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "hVoWadj_CLEAN_introns", test_dims = c(50,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#reduce and visulize the final dataset
final.dims = 45
final.res = 1.9
min.dist = 0.6
n.neighbors = 75

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "hVoWadj_CLEAN_introns", final.dims = final.dims, final.res = final.res, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = min.dist, n.neighbors = n.neighbors, assay = "integrated", saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#stash addtional metadata
Idents(seu.obj) <- "orig.ident"
seu.obj$cellSource <- ifelse(grepl("ealthy", seu.obj@meta.data$orig.ident), "Healthy", "Osteosarcoma")

#save the RDS for downstream analysis
outfile <- paste("./output/s3/hVoWadj_CLEAN_introns_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors,"_S3.rds", sep = "")
saveRDS(seu.obj, file = outfile)


######################
### BEGIN FIGURE 2 ###
######################

### Laod in processed data
seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "ageGroup")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "prettyName")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")
#reorder majorID factor levels
seu.obj$majorID <- factor(seu.obj$majorID, levels = c("CD4 T cell", "CD8/NK cell","gd T cell", "DN T cell",
                                                      "Cycling T cell", "Monocyte", "DC", "Neutrophil",
                                                     "B cell","Plasma cell","Granulocyte", "CD34+ unk"))

colArray <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
outName <- "fig2"


# for testing purposes
# seu.obj$ageGroup <- factor(seu.obj$ageGroup, levels = c("Healthy", "Middle-aged_OS", "Old_OS"))
# seu.obj <- subset(seu.obj,
#                   subset = 
#                   ageGroup !=  "Healthy"
#                  ) 
# seu.obj$ageGroup <- droplevels(seu.obj$ageGroup)
# seu.obj$cellSource <- factor(seu.obj$cellSource)
# seu.obj$name <- droplevels(seu.obj$name)

### Fig supp: transfer TRDC data from ROS data
seu.obj.ros <- readRDS(file = "./ros_output/HvO_mt_r7_res1.5_dims50_S3.rds")

#plot inital cluster umap
pi <- DimPlot(seu.obj.ros,
        reduction = "umap", 
        group.by = "clusterID",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_supp_ros_byClus.png", sep = ""), width =7, height = 7)

# for testing purposes
# get barcodes for Clusters 12 & 14
# seu.obj.ros.trdcClus <- subset(seu.obj.ros.cyto, subset = 
#                                clusterID ==  12 | clusterID == 14
#                               ) 

highlight <- colnames(seu.obj)[colnames(seu.obj) %in% WhichCells(seu.obj.ros, expression = TRDC > 0)]

pi <- DimPlot(seu.obj,#seu.integrated.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = "grey",
        pt.size = 0.5,
        label = F,
        cells.highlight= highlight,
        cols.highlight = "orchid",
        order = F
        #label.box = TRUE
 )

p <- formatUMAP(pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_supp_highlight_trdcPos.png", sep = ""), width =7, height = 7)


### Fig 2a: plot final colorized cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              cols = colArray$color,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_a_umap_final.png", sep = ""), width =7, height = 7)


### Fig 2b: umap by sample
# NOTE: down sample data to get equal number of cells from each sample
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name",
              cols = levels(seu.obj.ds$colz),
              pt.size = 0.25,
              label = FALSE,
              shuffle = TRUE
)
p <- formatUMAP(pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_b_UMAPbySample.png", sep = ""), width =7, height = 7)


### Fig supp: singleR classifier
singleR(seu.obj = seu.obj, outName = "hVoWadj_CLEAN_introns", clusters = "clusterID", outDir = "./output/fig2/singleR/")



### Fig supp: feature plots
features <- c("CD3G","CD8A", "GZMA",
              "IL7R","CD4", "S100A12",
              "DLA-DRA","FLT3", "CD79A",
              "MS4A1","JCHAIN", "IL5RA",
              "TOP2A","GATA3", "CD34")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp_keyFeatures.png", sep = ""), width =7, height = 7)


### Supplemetal data: Generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = "all_hVoWadj_CLEAN_introns",
                     outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )


### Fig 2c: frequency plots
fig2c <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source",
               namez = "name", 
               colz = "colz"
              ) + theme(axis.text.x = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_c_freqPlots_majorID.png", sep = ""), width =7, height = 6)


### Fig supp: stacked bar hvo
supp <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "orig.ident", 
                   clusters = "clusterID") + scale_fill_manual(labels = levels(seu.obj$name),
                                                               values = levels(seu.obj$colz)
                                                              ) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow =T))
ggsave(paste("./output/", outName, "/", outName, "_supp_stackedBar_hvo.png", sep = ""), width =8, height = 12)


### Fig 1d: Complete linDEG in pseudobulk-type format
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = F, groupBy = "allCells", comparision = "cellSource", outDir = "./output/fig2/", outName = "fig2_d_all", colUp = "red", colDwn = "blue",subtitle = F, labCutoff = 10, pValCutoff = 0.01, flipLFC = T, saveGeneList = T)


### Fig 1d: Complete linDEG in each major cell type
linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = F, groupBy = "prettyName", comparision = "cellSource", outDir = "./output/fig2/", outName = "fig2_d_major", colUp = "red", colDwn = "blue",subtitle = F, labCutoff = 10, flipLFC = T, pValCutoff = 0.01, saveGeneList = T)


### Fig supp: explore age effects
#reorder majorID factor levels
seu.obj$ageGroup <- factor(seu.obj$ageGroup, levels = c("Healthy", "Middle-aged_OS", "Old_OS"))

#remove OS_1,2,4,5,and 10 from analysis
seu.obj.ageSub <- subset(seu.obj,
                  subset = 
                  name !=  "OS_1" & name !=  "OS_2" & name !=  "OS_4" & name !=  "OS_5" & name !=  "OS_10"
                 ) 

p <- freqPlots(seu.obj = seu.obj.ageSub, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source", comp = "ageGroup",
               namez = levels(seu.obj.ageSub$name)[c(1:7,10,13,14,15,16)], 
               colz = levels(seu.obj.ageSub$colz)[c(1:7,10,13,14,15,16)]
                                                  )+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste("./output/", outName, "/", outName, "_supp_ageMatched_freqPlots_hvo.png", sep = ""), width =12, height = 6)


p <- freqPlots(seu.obj = seu.obj, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source", comp = "ageGroup",
               namez = "name", 
               colz = "colz"
                                                  ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste("./output/", outName, "/", outName, "_supp_ageMatched_freqPlots2_hvo.png", sep = ""), width =12, height = 6)


seu.obj.ageSub <- subset(seu.obj,
                  subset = 
                  ageGroup !=  "Healthy"
                 ) 

p <- freqPlots(seu.obj = seu.obj.ageSub, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source", comp = "ageGroup",
               namez = levels(seu.obj.ageSub$name)[8:17], 
               colz = levels(seu.obj.ageSub$colz)[8:17]
                                                  )+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste("./output/", outName, "/", outName, "_supp_ageMatched_freqPlots3_hvo.png", sep = ""), width =12, height = 6)


### Fig Extra: during preeration of uploading data to UCSC cell browser I noticed clus38 was largely coming from OS dogs -- using the following code it is apprent that the cluster has pal contaminatoin and is not of interest -- as suggested through independednt reclustering
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID", idents.1 = "38", idents.2 = c("3","4","7","23"), bioRep = "name",
                      padj_cutoff = 0.05, lfcCut = 0.58, minCells = 5, outDir = paste0("./output/", outName, "/"), 
                      title = "c38_vs_cyto", idents.1_NAME = "c38", idents.2_NAME = "cyto", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = "",lowFilter = T, dwnSam = F
                     )


###########################################
### END FIGURE 2  | BEGIN CYTO ANALYSIS ###
###########################################


#load in parent  
seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")

highlight <- c("CD8/NK cell","DN T cell")
colArray <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
colArray <- colArray %>% mutate(high_col = ifelse(majorID %in% highlight, color,"grey"))

umapHighLight <- pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank())  + NoLegend()

seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "Cytotoxic cell" | majorID ==  "DN T cell"
                 ) 

indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "cytotoxic_noGD", preSub = T,
         vars.to.regress = c("percent.mt", "percent.pal")
         )

seu.obj <- readRDS(file = "./output/s2/cytotoxic_noGD_S2.rds")

clusTree(seu.obj = seu.obj, outName = "cytotoxic_noGD", dout = "./output/clustree/", test_dims = c(45,40,35,30,25), algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "cytotoxic_noGD_hVoWadj_CLEAN_introns", 
           final.dims = 35, final.res = 0.8, returnFeats = F, saveRDS = T, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.3,
            n.neighbors = 30, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )


######################
### BEGIN FIGURE 3 ###
######################

#load in processed data
seu.obj <- readRDS(file = "./output/s3/cytotoxic_noGD_hVoWadj_CLEAN_introns_res0.8_dims35_dist0.3_neigh30_S3.rds")

#load in metadata
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "ageGroup")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "classification")
#reorder majorID factor levels
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("Naive CD8", "Effector CD8","Memory CD8",
                                                      "NK/NKT cell", "DN T cell", "CD8 gd T cell"))

#load in color data for plotting beauty
colArray <- read.csv("./cyto_hVoWadj.csv", header = T)
colArray$clusterID_sub <- as.factor(sort(as.numeric(as.character(colArray$clusterID_sub))))
outName <- "fig3"

# for testing purposes
# seu.obj$ageGroup <- factor(seu.obj$ageGroup, levels = c("Healthy", "Middle-aged_OS", "Old_OS"))
# seu.obj <- subset(seu.obj,
#                   subset = 
#                   ageGroup !=  "Old_OS"
#                  ) 

### Fig supp: transfer TRDC data from ROS data
seu.obj.ros.cyto <- readRDS(file = "./ros_output/HvO_mt_r7_res1.5_dims50_S3.rds")

seu.obj <- AddMetaData(seu.obj, metadata = seu.obj.ros.cyto@assays$RNA@data[rownames(seu.obj.ros.cyto@assays$RNA@data) == "TRDC",], col.name = "TRDC")
seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 

features = c("TRDC")
p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 1, ncol = 1, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_/supp_trdc_cyto_featPlots.png", sep = ""), width =3, height = 3)

pi <- DimPlot(seu.obj,#seu.integrated.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = "grey",
        pt.size = 0.5,
        label = F,
        cells.highlight= highlight,
        cols.highlight = "orchid",
        order = F
        #label.box = TRUE
 )

p <- formatUMAP(pi)
ggsave("./output/umap_cyto_highlight_trdcPos.png")
ggsave(paste("./output/", outName, "/", outName, "_supp_trdc_cyto_Highlight.png", sep = ""), width =3, height = 3)


### Extra: Supplemental data: Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "cyto_hVoWadj_CLEAN_introns",
                     outDir = "./output/fig3/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )


### Fig Extra: Dot plot by feature type
namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub
p <- dotPlotBY_TYPE(seu_obj = seu.obj, pwdTOvilnCSVoutput = "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/fig3/viln/cyto_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave(paste("./output/", outName, "/", outName, "_supp_bigDots.png", sep = ""), width =20, height = 12)


### Fig Extra: SingleR classifications
singleR(seu.obj = seu.obj, outName = "cyto", clusters = "clusterID_sub", outDir = "./output/fig3/singleR/")


### Fig 3a: plot clustering results - req umapHighLight fig generated from parent
pi <- DimPlot(seu.obj,
              reduction = "umap", 
              group.by = "clusterID",
              cols = colArray$colour,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

fig3a <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol) + NoLegend()
fig3a <- fig3a + inset_element(umapHighLight, left=0, bottom=0, right=0.25, top=0.25)
ggsave(paste("./output/", outName, "/", outName, "_a_umap_cyto_raw_clusterID.png", sep = ""), width = 7, height = 7)


### Fig 3a: create custom legend for UMAP
colArray$majorID <- "CD8, DN T, & NK cell subtypes"
leg <- cusLeg(legend = colArray, colz = 6, rowz = 2, clusLabel = "clusterID_sub", legLabel = "classification", colorz = "colour",
                   groupLabel = "majorID", groupBy = "majorID", sortBy = "clusterID_sub", labCol = "labCol", headerSize = 6,
                   cellHeader = T, bump = 0, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 2.5, compress_x = 0.6, spaceBtwnCols = c(0.5,0.5,0.5,0.5,0.4)
             )
ggsave(paste("./output/", outName, "/", outName, "_a_legend.png", sep = ""), width = 12, height = 3)


### Fig supp: plot sankey plot -- colorization differernt (but data the same) due to package update :(
colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
colArray_all <- colArray_all %>% filter(majorID %in% highlight)

p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_sub", old.ident = "clusterID", old.colorz = colArray_all$color,
                       new.colorz = colArray$colour, old.labCol = colArray_all$labCol, new.labCol = colArray$labCol, flowCol = "grey"
                    )

ggsave(paste("./output/", outName, "/", outName, "_supp_cyto_sankey.png", sep = ""), width = 7, height = 10)


### Fig 3b: create violin plots
features <- c("CD3E", "CD4", "CD8A", "IL2RB",
              "CCL4", "CCL5", "CCR5", "CCR7",
              "NCR3", "KLRB1", "CD96","IL7R",
              "NFKBID","GZMA", "GZMB", "GZMK",
              "ITGAM","RARRES2","KLRF1", "LGALS3", 
              "PI3","ID3", "GATA3", "IKZF2")

pi <- VlnPlot(
    object = seu.obj,
    cols = colArray$colour,
    pt.size = 0,
    same.y.lims = T,
    combine = FALSE,
    features = features
)
    

p <- lapply(pi, function(x) x + theme(axis.title = element_blank(),
                                      axis.text = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks = element_blank(),
                                      title = element_text(size = 20),
                                      legend.position = "none",
                                      plot.margin = margin(5.5, 0, 0, 0, "pt")
                                     ) + coord_cartesian(expand = TRUE)
            )

g <- ggplot(colArray, aes(x = clusterID_sub, y = 1, fill = clusterID_sub, label = clusterID_sub)) + 
    geom_tile(fill = "transparent") +
    geom_point(shape = 21, stroke = 0.75, size = 11) +
    geom_text(fontface = "bold", size = 6, color = colArray$labCol) + theme_bw(base_size = 12) +
    scale_fill_manual(values = colArray$colour) + scale_y_discrete(expand = c(0, 0)) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_blank(), 
              panel.border = element_blank(),
              plot.background = element_blank(), 
              plot.margin = margin(0, 0, 0, 0, "pt"),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()) #+ xlab("Cluster")

      
nrow = 7
ncol = 4
patch <- area()
counter=0
for (i in 1:nrow) {
    for (x in 1:ncol) {
        counter = counter+1
        if (counter <= ncol*nrow) {
            patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }
}

pi <- Reduce( `+`, p ) + g + g + g + g + plot_layout(design = patch, widths = c(rep.int(1, ncol)), heights = c(rep.int(1, nrow-1),0.3)) 

ggsave(paste("./output/", outName, "/", outName, "_b_vilnPlots_cyto_clusterID_sub.png", sep = ""), width = 22.45, height = 13.4)


### Fig supp: plot same features as feat plot
fig4b <- prettyFeats(seu.obj = seu.obj, nrow = 6, ncol = 4, features = features, color = "black", order = F,bottomLeg = T) 
ggsave(paste("./output/", outName, "/", outName, "_supp_featPlots.png", sep = ""), width = 12, height = 16)


### Fig supp: plot enrichment scores
ecLists <- read.csv("ecScores.csv", header = F)
modulez <- split(ecLists$V2, ecLists$V1)
names(modulez) <- paste0(names(modulez),"_SIG")

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

fig_supp <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 30, "pt")) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", outName, "_supp_cyto_modScores.png", sep = ""), width = 10, height = 6)


### Fig 3c: DEGs for DN T cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "classification", idents.1 = "DN T cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/fig3/btwnClusDEG/", title = "DN T cell vs Clusters 0-7,9-11", idents.1_NAME= "DN_T_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

pi <- p[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_c_dnTs_vs_otherCD8tCells_volc.png", sep = ""), width = 7, height = 7)


### Fig supp: DEGs for NK T cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "classification", idents.1 = "NK T cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/fig3/btwnClusDEG/", title = "NK T cell vs Clusters 0-9,11", idents.1_NAME= "NK_T_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

pi <- p[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_supp_nkTs_vs_otherCD8tCells_volc.png", sep = ""), width = 7, height = 7)


### Fig 3d: DEGs for NK cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "classification", idents.1 = "NK cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/fig3/btwnClusDEG/", title = "NK cell vs Clusters 1-8,10,11", idents.1_NAME= "NK_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

pi <- p[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_d_nks_vs_otherCD8tCells_volc.png", sep = ""), width = 7, height = 7)


### Fig 3e: DEGs for cd8 gd T cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "majorID_sub", idents.1 = "CD8 gd T cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/fig3/btwnClusDEG/", title = "CD8+ gd T cell vs Clusters 1-10", idents.1_NAME= "CD8_gd_T_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

pi <- p[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_e_cd8gdT_vs_otherCD8tCells_volc.png", sep = ""), width = 7, height = 7)


### Fig 3f: freqPlots by majorID_sub
pi <- freqPlots(seu.obj, method = 1, nrow= 2, groupBy = "majorID_sub", legTitle = "Cell source",
               namez = "name", 
               colz = "colz",
              ) + scale_x_discrete(labels=c('Healthy', 'OS'))

ggsave(paste("./output/", outName, "/", outName, "_f_freqPlots_cyto_majorID_sub.png", sep = ""), width = 6, height = 6)


### Fig Extra: linDEG by cluster
linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = F, groupBy = "clusterID_sub", comparision = "cellSource", outDir = "./output/fig3/linDEG/", outName = "cyto", colUp = "red", colDwn = "blue",subtitle = F)


### Fig 3g: plot key feats altered in disease - as determined by linDEG/FindMarkers
colArray$clusterID_sub <- as.factor(sort(as.numeric(as.character(colArray$clusterID_sub))))
p <- vilnSplitCompxGene(seu.obj = seu.obj, groupBy = "clusterID_sub", comp = "cellSource", metaAdd = "majorID", features = c("CCL5","FCER1G", "IL2RB"), 
                               cols = c("mediumseagreen","mediumpurple1"), save = FALSE, outName = "", outDir = "", 
                               height = 4, width = 8, labelData = colArray
                              )
ggsave(paste("./output/", outName, "/", outName, "_g_viln_cyto_genezviln_cyto_genez.png", sep = ""), width = 7, height = 4)


#############################################
### END FIGURE 3  | BEGIN HELPER ANALYSIS ###
#############################################

#load in parent  
seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")

highlight <- c("CD4 T cell")
colArray <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
colArray <- colArray %>% mutate(high_col = ifelse(majorID %in% highlight, color,"grey"))

umapHighLight <- pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank())

seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "CD4 T cell"
                 ) 

indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "helper_hVoWadj_CLEAN_introns", preSub = T,
         vars.to.regress = c("percent.mt", "percent.pal")
         )

seu.obj <- readRDS(file = "./output/s2/helper_hVoWadj_CLEAN_introns_S2.rds")

clusTree(seu.obj = seu.obj, outName = "helper_hVoWadj_CLEAN_introns", dout = "./output/clustree/", test_dims = c(45,40,35,30,25), algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "helper_hVoWadj_CLEAN_introns", 
           final.dims = 40, final.res = 0.9, returnFeats = F, saveRDS = T, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.2,
            n.neighbors = 20, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )


######################
### BEGIN FIGURE 4 ###
######################

### Load in processed data
seu.obj <- readRDS(file = "./output/s3/helper_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.2_neigh20_S3.rds")

#load in meta data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "clusterID_sub_clean")
seu.obj$clusterID_sub_clean <- factor(x = Idents(seu.obj), levels = sort(as.numeric(as.character(levels(seu.obj)))))
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")

#create subset object to remove por qulity cells from downstream analysis
seu.obj.sub <- subset(seu.obj, subset = majorID_sub != "exclude")

#reorder majorID factor levels
seu.obj.sub$majorID_sub <- droplevels(seu.obj.sub$majorID_sub)
seu.obj.sub$majorID_sub_forPLOT <- gsub("CD4\\+ ","",seu.obj.sub$majorID_sub)
seu.obj.sub$majorID_sub_forPLOT <- gsub("CD4\\+, ","",seu.obj.sub$majorID_sub_forPLOT)
seu.obj.sub$majorID_sub_forPLOT <- factor(seu.obj.sub$majorID_sub_forPLOT, levels = rev(c("Naive", "TCM","TEM", "TEM, Th1-like", 
                                                              "TEM, Th2-like","TEM, Th17-like","T reg",
                                                              "IFN signature")))

#load in color data for pretty plots
colArray <- read.csv("./helper_hVoWadj_forLeg_020623.csv", header = T)
colArray <- colArray[-14,] # remove the first poor quality pop from list

namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub_clean
outName <- "fig4"

### Fig supp: singleR classifier
singleR(seu.obj = seu.obj.sub, outName = "cd4", clusters = "clusterID_sub", outDir = "./output/fig4/singleR/")


### Supplmental data: Generate violing plots of defining features
vilnPlots(seu.obj = seu.obj.sub, groupBy = "clusterID_sub_clean", numOfFeats = 24, outName = "cd4_hVoWadj_CLEAN_introns",
                     outDir = "./output/fig4/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )


### Fie Extra: create DotPlot by feature type
p <- dotPlotBY_TYPE(seu_obj = seu.obj.sub, pwdTOvilnCSVoutput = "./output/fig4/viln/cd4_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub_clean",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave(paste("./output/", outName, "/", outName, "_helper_bigDots.png", sep = ""), width = 20, height = 12)


### Fig 4a: plot colorized umap
p <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub_clean",
              cols = colArray$colour,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

p <- cusLabels(plot = p, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol)
p <- p + inset_element(umapHighLight, left=0, bottom=0, right=0.25, top=0.25)
ggsave(paste("./output/", outName, "/", outName, "_a_UMAP_helper.png", sep = ""), width = 7, height = 7)


### Fig 4a: create custom legend
leg <- cusLeg(legend = colArray, colz = 4, rowz = 4, clusLabel = "clusterID_sub_clean", legLabel = "classification", colorz = "colour",
                   groupLabel = "title", groupBy = "title", sortBy = "clusterID_sub", labCol = "labCol", headerSize = 6,
                   cellHeader = T, bump = 0, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 2, compress_x = 0.8, spaceBtwnCols = c(0.4,0.45,0.3)
             )
ggsave(paste("./output/", outName, "/", outName, "_a_UMAP_legend_helper.png", sep = ""), width = 12, height = 3)



### Fig supp: create sankey plot
colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
saveCells <- c("CD4 T cell")
colArray_all <- colArray_all %>% filter(majorID %in% saveCells)

p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_sub_clean", old.ident = "clusterID", old.colorz = colArray_all$color,
                       new.colorz = colArray$colour, old.labCol = colArray_all$labCol, new.labCol = colArray$labCol, flowCol = "grey"
                    )
ggsave(paste("./output/", outName, "/", outName, "_supp_helper_sankey.png", sep = ""), width = 7, height = 10)


### Fig supp: create stacked bar graph
stackedBar(seu.obj = seu.obj.sub, downSampleBy = "cellSource", groupBy = "orig.ident", clusters = "clusterID_sub_clean") +
scale_fill_manual(labels = c("H_1", "H_2","H_3","H_4", "H_5","H_6", "H_7", 
                         "OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6", "OS_7", "OS_8", "OS_9", "OS_10"), 
               values = c(sequential_hcl(palette = "Greens", n = 10)[1:7], 
                        sequential_hcl(palette = "Purp", n = 12)[1:10])) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow =T))

ggsave(paste("./output/", outName, "/", outName, "_supp_stacked.png", sep = ""), width = 10, height = 7.5)


### Fig 4b: freqPlots by majorID_sub
features <- c("CD3G","CD4", "CD8A","IL7R",
              "SELL","CCR7", "CD40LG", "CD28",
              "IL18R1","CCR5","CTLA4","IKZF2",
              "S100A5", "GATA3", "SYTL3", "RUNX2"
              )
fig4b <- prettyFeats(seu.obj = seu.obj, nrow = 4, ncol = 4, features = features, color = "black", order = F,bottomLeg = T) 
ggsave(paste("./output/", outName, "/", outName, "_b_helper_feats.png", sep = ""), width = 12, height = 12)


### Fig supp: pathway analysis on Cluster 14
clus.markers <- read.csv(file = "./output/fig4/viln/cd4_hVoWadj_CLEAN_introns_gene_list.csv")
clus14.markers <- clus.markers[clus.markers$cluster == 14,]

can_gene_sets <- as.data.frame(msigdbr(species = "dog", category = "H"))
msigdbr_list <- split(x = can_gene_sets$gene_symbol, f = can_gene_sets$gs_name)
datas <- can_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

enriched <- as.data.frame(enricher(gene = clus14.markers$gene, TERM2GENE = datas, pvalueCutoff = 1))

enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.05,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Hallmarks pathway enrichment")

ggsave(paste("./output/", outName, "/", outName, "_supp_gsea_c14.png", sep = ""), width = 8, height = 7)


### Fig 4c: dotPlot by keyFeats
seu.obj.sub <- ScaleData(seu.obj.sub)
p <- majorDot(seu.obj = seu.obj.sub, groupBy = "majorID_sub",
                     features = c("ZNF536","LEF1","SELL",
                                  "CCR7", "TSHZ2",
                                  "CD38","IL7R","CD44",
                                  "TBX21", "EOMES","GZMK",
                                  "LGALS3","CCDC3","GATA3",
                                  "RORC","RORA","SMAD2","TGFB1","TNFRSF8",
                                  "FOXP3", "CTLA4",
                                  "CCL5","IFNG","BCL6","IRF1"
                                 )
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 24, "pt")) + scale_y_discrete(position = "right")
ggsave(paste("./output/", outName, "/", outName, "_c_helper_dots.png", sep = ""), width = 8, height = 4)


### Fig supp: plot enrichment scores
ecLists <- read.csv("ecScores.csv", header = F)
modulez <- split(ecLists$V2, ecLists$V1)

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                          features = modulez,
                         name = "_score")

names(seu.obj.sub@meta.data)[grep("_score", names(seu.obj.sub@meta.data))] <- paste0(names(modulez),"_SIG")

p <- majorDot(seu.obj = seu.obj.sub, groupBy = "majorID_sub_forPLOT",
                     features = c("Naïve_SIG", "TH1_SIG","TH2_SIG", "TH17_SIG",
                                  "TREG_SIG")
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 24, "pt")) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", outName, "_supp_helper_ecScores.png", sep = ""), width = 8, height = 4)


### Fig Extra: investigated gene expression changes in os -- not much going on....
vilnSplitComp(seu.obj = seu.obj.sub, groupBy = "clusterID_sub_clean", refVal = "cellSource", 
              outName = "cd4", nPlots = 20, outDir = "./output/fig4/spltViln/"
                       )


### Fig Extra: plot cell frequency plots -- again not much going on....
p <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "majorID_sub", legTitle = "Cell source",
               namez = c("H_1", "H_2","H_3","H_4", "H_5","H_6", "H_7", 
                         "OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6", "OS_7", "OS_8", "OS_9", "OS_10"), 
               colz = c(sequential_hcl(palette = "Greens", n = 10)[1:7], 
                        sequential_hcl(palette = "Purp", n = 12)[1:10])
              ) + theme(legend.position="none") 
ggsave(paste("./output/", outName, "/", outName, "_supp_freqPlots.png", sep = ""), width = 7, height = 5)


### Fig 4d: volcano plots for cd4 subtypes
#clean matedata solt
seu.obj.sub$majorID_sub_forPLOT2 <- gsub(",","_",seu.obj.sub$majorID_sub_forPLOT)
seu.obj.sub$majorID_sub_forPLOT2 <- gsub(" ","_",seu.obj.sub$majorID_sub_forPLOT2)
seu.obj.sub$majorID_sub_forPLOT2 <- gsub("-","_",seu.obj.sub$majorID_sub_forPLOT2)
seu.obj.sub$majorID_sub_forPLOT2 <- factor(seu.obj.sub$majorID_sub_forPLOT2, levels = rev(c("Naive", "TCM","TEM", "TEM__Th1_like", 
                                                              "TEM__Th2_like","TEM__Th17_like","T_reg",
                                                              "IFN_signature")))

#set up for pathway analysis
can_gene_sets <- as.data.frame(msigdbr(species = "dog", category = "C7"))
msigdbr_list <- split(x = can_gene_sets$gene_symbol, f = can_gene_sets$gs_name)
datas <- can_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

#narrow the search....
datas <- datas[grepl("CD4|TREG|TH1|TH17|TH2",datas$gs_name),]


### volc for Th1s
fig4d1 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub_forPLOT2", idents.1 = "TEM__Th1_like", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/fig4/btwnClusDEG/", title = "TEM,Th1-like vs other CD4 clusters", idents.1_NAME= "TEM__Th1_like", idents.2_NAME = "Other_T_cells", returnVolc = T, doLinDEG = F,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

fig4d1 <- fig4d1[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0) 
ggsave(paste("./output/", outName, "/", outName, "_d_th1_volc.png", sep = ""), width = 7, height = 7)


#use features up Th1 for pathway analysis
th1.df <- read.csv("./output/fig4/btwnClusDEG/TEM__Th1_like_vs_Other_T_cells_all_genes.csv")
th1.df <- th1.df[th1.df$log2FoldChange > 0,]

enriched <- as.data.frame(enricher(gene = th1.df$gene, TERM2GENE = datas, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.01,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Hallmarks pathway enrichment")

ggsave(paste("./output/", outName, "/", outName, "_supp_th1_gsea.png", sep = ""), width = 16, height = 10)


### volc for Th2s
fig4d2 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub_forPLOT2", idents.1 = "TEM__Th2_like", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/fig4/btwnClusDEG/", title = "TEM, Th2-like vs other CD4 clusters", idents.1_NAME= "Th2", idents.2_NAME = "Other_T_cells", returnVolc = T,doLinDEG = F, lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

fig4d2 <- fig4d2[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_d_th2_volc.png", sep = ""), width = 7, height = 7)

#run gsea
th2.df <- read.csv("./output/fig4/btwnClusDEG/Th2_vs_Other_T_cells_all_genes.csv")
th2.df <- th2.df[th2.df$log2FoldChange > 0,]

enriched <- as.data.frame(enricher(gene = th2.df$gene, TERM2GENE = datas, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.01,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Hallmarks pathway enrichment")

ggsave(paste("./output/", outName, "/", outName, "_supp_th2_gsea.png", sep = ""), width = 16, height = 10)



### volc for Th17s
fig4d3 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub_forPLOT2", idents.1 = "TEM__Th17_like", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/fig4/btwnClusDEG/", title = "TEM, Th17-like vs other CD4 clusters", idents.1_NAME= "Th17", idents.2_NAME = "Other_T_cells", returnVolc = T, doLinDEG = F,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

fig4d3 <- fig4d3[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_d_th17_volc.png", sep = ""), width = 7, height = 7)

#run gsea
th17.df <- read.csv("./output/fig4/btwnClusDEG/Th17_vs_Other_T_cells_all_genes.csv")
th17.df <- th17.df[th17.df$log2FoldChange > 0,]

enriched <- as.data.frame(enricher(gene = th17.df$gene, TERM2GENE = datas, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.01,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Hallmarks pathway enrichment")
ggsave(paste("./output/", outName, "/", outName, "_supp_th17_gsea.png", sep = ""), width = 16, height = 10)


### volc for Tregs
fig4d4 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub_forPLOT2", idents.1 = "T_reg", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/fig4/btwnClusDEG/", title = "T reg vs other CD4 clusters", idents.1_NAME= "T_reg", idents.2_NAME = "Other_T_cells", returnVolc = T,, doLinDEG = F,lowFilter = T, dwnSam = F, setSeed = 12
                    ) 

fig4d4 <- fig4d4[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_d_treg_volc.png", sep = ""), width = 7, height = 7)

#run gsea
tregs.df <- read.csv("./output/fig4/btwnClusDEG/T_reg_vs_Other_T_cells_all_genes.csv")
tregs.df <- tregs.df[tregs.df$log2FoldChange < 0,]

enriched <- as.data.frame(enricher(gene = tregs.df$gene, TERM2GENE = datas, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust > 0.01,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Hallmarks pathway enrichment")

ggsave(paste("./output/", outName, "/", outName, "_supp_treg_gsea.png", sep = ""), width = 16, height = 10)


### Fig supp: run SlingShot
sce.obj <- as.SingleCellExperiment(seu.obj.sub)
rd1 <- Embeddings(seu.obj.sub, reduction = "pca")[,1:2]
rd2 <- Embeddings(seu.obj.sub, reduction = "umap")
colnames(rd2) <- c('UMAP1', 'UMAP2')

#assign origin -- largest naive cd4 t cell clus
start.clus <- '0'

reducedDims(sce.obj) <- SimpleList(PCA = rd1, UMAP = rd2)

sce.obj <- slingshot(sce.obj, clusterLabels = 'clusterID_sub_clean', reducedDim = 'UMAP', start.clus = start.clus)

#identify lineages
lin1 <- getLineages(Embeddings(seu.obj.sub, reduction = "umap"), sce.obj$clusterID_sub_clean, start.clus = start.clus)
branchData <- SlingshotDataSet(lin1)@lineages

#plot the lineages
colArray <- read.csv("./helper_hVoWadj_forLeg.csv")
plot <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "clusterID_sub_clean",
              cols = colArray$colour[-c(14)],
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

p <- cleanSling(plot = plot, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = branchData)
ggsave(paste("./output/", outName, "/", outName, "_supp_cd4Bracnch.png", sep = ""), width = 7, height = 7)

#get the meta data over - to the seurat object
metaData <- as.data.frame(sce.obj@colData@rownames)
colnames(metaData) <- "barcode"
metaData$slingPseudotime_1 <- sce.obj$slingPseudotime_1
metaData$slingPseudotime_2 <- sce.obj$slingPseudotime_2
metaData$slingPseudotime_3 <- sce.obj$slingPseudotime_3
metaData$slingPseudotime_4 <- sce.obj$slingPseudotime_4
metaData$slingPseudotime_5 <- sce.obj$slingPseudotime_5
metaData$slingPseudotime_6 <- sce.obj$slingPseudotime_6
metaData$slingPseudotime_7 <- sce.obj$slingPseudotime_7

metaData <- metaData %>% rowwise %>% mutate(pseudoTime = mean(c(slingPseudotime_1,slingPseudotime_2,slingPseudotime_3,slingPseudotime_4,slingPseudotime_5,slingPseudotime_6,slingPseudotime_7),na.rm=TRUE))

seuMeta <- seu.obj.sub@meta.data %>% mutate(barcode = rownames(.))

newMeta <- seuMeta %>% left_join(metaData, by = 'barcode')
rownames(newMeta) <- newMeta$barcode
seu.obj.sub@meta.data <- newMeta

#plot the lineages on original umap rep
features <- c("slingPseudotime_5","slingPseudotime_2","slingPseudotime_3","slingPseudotime_1")
titles <- c("Lineage 1","Lineage 2","Lineage 3","Lineage 4")

p <- prettyFeats(seu.obj = seu.obj.sub, nrow = 2, ncol = 2, features = features, color = "black", order = F, titles = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

ggsave(paste("./output/", outName, "/", outName, "_supp_pseudoTime_helper.png", sep = ""), width = 8, height = 8)


##############################################
### END FIGURE 4  | BEGIN MYELOID ANALYSIS ###
##############################################

seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")

highlight <- c("Monocyte","Neutrophil","Granulocyte","DC")
colArray <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
colArray <- colArray %>% mutate(high_col = ifelse(majorID %in% highlight, color,"grey"))

umapHighLight <- pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank()) + NoLegend()


seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "Monocyte" | majorID ==  "Neutrophil" | majorID ==  "Granulocyte" | majorID ==  "DC"
                 ) 

indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "myeloid", preSub = T,
         vars.to.regress = c("percent.mt", "percent.pal")
         )

seu.obj <- readRDS(file = "./output/s2/myeloid_S2.rds")

clusTree(seu.obj = seu.obj, outName = "myeloid", dout = "./output/clustree/", test_dims = c(45,40,35,30,25), algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "myeloid_hVoWadj_CLEAN_introns", 
           final.dims = 40, final.res = 0.9, returnFeats = F, saveRDS = F, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.35,
            n.neighbors = 20, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )


######################
### BEGIN FIGURE 5 ###
######################

### Load in processed data
seu.obj <- readRDS(file = "./output/s3/myeloid_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.35_neigh20_S3.rds")

#load in key meta data and remove poor quality cells for downstream analysis
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "clusterID_sub_clean")
seu.obj$clusterID_sub_clean <- factor(x = Idents(seu.obj), levels = sort(as.numeric(as.character(levels(seu.obj)))))

seu.obj.sub <- subset(seu.obj, subset = majorID_sub != "exclude")
seu.obj.sub <- loadMeta(seu.obj = seu.obj.sub, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "ageGroup2")

# for testing purposes
# seu.obj.sub$ageGroup2 <- factor(seu.obj.sub$ageGroup2, levels = c("Healthy", "Middle-aged_OS", "Old_OS"))
# seu.obj.sub <- subset(seu.obj.sub,
#                   subset = 
#                   ageGroup2 !=  "Old_OS"
#                  ) 
# seu.obj.sub$ageGroup2 <- droplevels(seu.obj.sub$ageGroup2)
# seu.obj.sub$cellSource <- factor(seu.obj.sub$cellSource)
# seu.obj.sub$name <- droplevels(seu.obj.sub$name)

#read in color meta data for ploting uniformity
colArray <- read.csv("./coldaf_myeloid_new.csv")
namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub_clean
outName <- "fig5"

### Fig supp: singleR classifier
singleR(seu.obj = seu.obj.sub, outName = "myeloid", clusters = "clusterID_sub_clean", outDir = "./output/fig5/singleR/")


### Supplmental data: Generate violing plots of defining features
vilnPlots(seu.obj = seu.obj.sub, groupBy = "clusterID_sub_clean", numOfFeats = 24, outName = "myeloid_hVoWadj_CLEAN_introns",
                     outDir = "./output/fig5/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )


### Fig Extra: create DotPlot by feature type
p <- dotPlotBY_TYPE(seu_obj = seu.obj.sub, pwdTOvilnCSVoutput = "./output/fig5/viln/myeloid_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub_clean",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave(paste("./output/", outName, "/", outName, "_myeloid_bigDots.png", sep = ""), width = 20, height = 12)


### Fig 5a: plot colorized umap
p <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub_clean",
              cols = c(colArray$colour[-c(16,19)],"grey"),
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )
p <- cusLabels(plot = p, shape = 21, size = 8, alpha = 0.8, labCol = c(colArray$labCol[-c(16,19)],NA)) + NoLegend()
p <- p + inset_element(umapHighLight, left=0, bottom=0.75, right=0.25, top=1)
ggsave(paste("./output/", outName, "/", outName, "_a_myeloid_umap.png", sep = ""), width = 20, height = 12)


### Fig 5a: create custom legend - with manual hack to force legend into position
colArray_forLeg <- read.csv("./coldaf_myeloid_new_forLeg.csv", header = T)
colArray_forLeg$clusterID_sub <- as.factor(sort(as.numeric(as.character(colArray_forLeg$clusterID_sub))))

leg <- cusLeg(legend = colArray_forLeg, colz = 5, rowz = 5, clusLabel = "clusterID_sub", 
              legLabel = "majorID_sub", colorz = "colour",
              groupLabel = "title", groupBy = "title", sortBy = "clusterID_sub", 
              labCol = "labCol", headerSize = 6, cellHeader = T, bump = 2, 
              nudge_left = 0, nudge_right = 0, topBuffer = 1.1, ymin = 0, 
              compress_y = 12, compress_x = 0.75, 
              titleOrder = c("Monocyte","Neutrophil","Granulocyte", "Low quality","Dendritic cell"),
              spaceBtwnCols = c(0.5,0.5,0.5), breakGroups = T, returnData = T 
             )
# write.csv(leg[1], "./output/legend.csv")
# write.csv(leg[2], "./output/header.csv")

leg <- read.csv("./output/legend.csv")
header <- read.csv("./output/header.csv")
overrideData <- list(leg,header)

leg <- cusLeg(legend = colArray_forLeg, colz = 5, rowz = 5, clusLabel = "clusterID_sub", 
              legLabel = "majorID_sub", colorz = "colour",
              groupLabel = "title", groupBy = "title", sortBy = "clusterID_sub", 
              labCol = "labCol", headerSize = 6, cellHeader = T, bump = 2, 
              nudge_left = 0, nudge_right = 0, topBuffer = 1.1, ymin = 0, 
              compress_y = 12, compress_x = 0.75, 
              titleOrder = c("Monocyte","Neutrophil","Granulocyte", "Low quality","Dendritic cell"),
              spaceBtwnCols = c(0.5,0.5,0.5), breakGroups = T, returnData = F, overrideData = overrideData
             )

ggsave(paste("./output/", outName, "/", outName, "_a_myeloid_legend_umap.png", sep = ""), width = 14, height = 10)


### Fig supp: sankey plot to show how clustering changed
colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
saveCells <- c("Monocyte","Neutrophil","Granulocyte","DC")
colArray_all <- colArray_all %>% filter(majorID %in% saveCells)

p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_sub", old.ident = "clusterID", 
                old.colorz = colArray_all$color, new.colorz = colArray$colour, 
                old.labCol = colArray_all$labCol, new.labCol = colArray$labCol, 
                flowCol = "#D87FAB"
                    )
ggsave(paste("./output/", outName, "/", outName, "_supp_sankey.png", sep = ""), width = 10, height = 14)


### Fig 5b: create feature plots
features <- c("DLA-DRA","ANPEP", "FLT3",
              "S100A12","CD3G", "CD4",
              "LTF","IL5RA", "BATF3",
              "IRF8","CADM1", "TCF4")

# #some additional features used to try to ID DCs
# features <- c("FLT3","IRF8", "IRF4",
#               "BATF3","FCER1A", "KLF4",
#               "RUNX3","AXL", "CADM1",
#               "BTLA","KIT", "CD1C",
#               "IL3RA","CDH1", "SDC2")

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_b_featPlots.png", sep = ""), width = 10, height = 15)


### Fig 5c: create frequency plots
#order myeloid cell populations
seu.obj.sub$majorID_sub2 <- factor(seu.obj.sub$majorID_sub2, levels = c("DC","Monocyte, CD4-", "Monocyte, IFN signature",
                                                                      "Monocyte, CD4+", "M-MDSC", "PMN-MDSC",
                                                                      "Neutrophil","Eosinophil","Basophil"))

p <- freqPlots(seu.obj = seu.obj.sub, method = 1, nrow = 3, 
                   groupBy = "majorID_sub2", legTitle = "Cell source",
                   refVal = "name",comp = "cellSource", no_legend = T,
                   namez = "name", colz = "colz"
                  ) + scale_x_discrete(labels=c('Healthy', 'OS'))#+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(paste("./output/", outName, "/", outName, "_c_freqPlots.png", sep = ""), width = 8, height = 6)


### Fig supp: linDEG on c10 (IFN sig)
geneList <- linDEG(seu.obj = seu.obj.sub, threshold = 0.5, thresLine = F, groupBy = "majorID_sub", comparision = "cellSource", outDir = "./output/fig5/", outName = "c10", colUp = "red", colDwn = "blue",subtitle = T, cluster = "Monocyte, IFN signature", returnUpList = F, returnDwnList = T, useLineThreshold = F, pValCutoff = 0.01, flipLFC = F, saveGeneList = F, addLabs = "")

### Fig supp: gsea on healthy linDEG gene sig
can_gene_sets <- as.data.frame(msigdbr(species = "dog", category = "H"))
msigdbr_list <- split(x = can_gene_sets$gene_symbol, f = can_gene_sets$gs_name)
datas <- can_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

enriched <- as.data.frame(enricher(gene = unlist(geneList), TERM2GENE = datas, pvalueCutoff = 1))

enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.05,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Hallmarks pathway enrichment")

ggsave(paste("./output/", outName, "/", outName, "_supp_clus10UpinHealth.png", sep = ""), width = 8, height = 7)


### Fig supp: additional feature plots
features <- c("DLA-DRA", "FCGR1A", "CD86", "CCR2")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 4, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp_featPlots.png", sep = ""), width = 13, height = 3)


### Extra analysis: compare to ros data -- must request ros data from authors; sorry for the inconvenience
ros_myeloid <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/ros_output/sub_H vO_Myeloid_mt_r7_res1.1_dims40_S3.rds")
features <- c("DLA-DRA", "FCGR1A", "CD86", "CCR2", "S100A12", 
              "CD14", "CD4", "ANPEP", "CD68")

p <- prettyFeats(seu.obj = ros_myeloid, nrow = 3, ncol = 3, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp_featPlots.png", sep = ""), width = 9, height = 9)

#transfer TRDC data from ROS data
highlight <- colnames(seu.obj.sub)[colnames(seu.obj.sub) %in% WhichCells(ros_myeloid, expression = CD14 > 0)]

pi <- DimPlot(seu.obj.sub,#seu.integrated.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = "grey",
        pt.size = 0.5,
        label = F,
        cells.highlight= highlight,
        cols.highlight = "orchid",
        order = F
        #label.box = TRUE
 )

p <- formatUMAP(pi)
ggsave(paste("./output/", outName, "/", outName, "_supp_CD14_highlight.png", sep = ""), width = 7, height = 7)


### Fig 5e and Fig 5g: complete linDEG for each major myeloid cell group
linDEG(seu.obj = seu.obj.sub, threshold = 0.5, thresLine = F, groupBy = "majorID_sub", 
       comparision = "cellSource", outDir = "./output/fig5/linDEG/", outName = "majorID_sub_degs", 
       cluster = NULL, labCutoff = 20, colUp = "red", colDwn = "blue", 
       subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, flipLFC = T, saveGeneList = F, addLabs = ""
                  )


### Fig Extra - another way to plot the data
vilnSplitComp(seu.obj = seu.obj.sub, groupBy = "majorID_sub", refVal = "cellSource", 
              outName = "myeloid", nPlots = 20, outDir = "./output/fig5/spltCompViln/"
                       )
                             

### Fig 5d: volc for PMN-MDSC
p <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", 
                      idents.1 = "PMN-MDSC", idents.2 = "Neutrophil", bioRep = "orig.ident",
                      padj_cutoff = 0.01, lfcCut = 0.58, minCells = 25, 
                      outDir = "./output/fig5/btwnClusDEG/", title = "PMN-MDSC vs Clusters 1, 7, 9 & 20", 
                      idents.1_NAME= "PMN-MDSC", idents.2_NAME = "Neutrophil", 
                      returnVolc = T,doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 12
                    ) 
p <- p[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_d_pmnMDSCvNeuts.png", sep = ""), width = 7, height = 7)

                             
### Fig 5f: volc for M-MDSC
p <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", 
                     idents.1 = "M-MDSC", idents.2 = c("Monocyte, CD4-","Monocyte, CD4+","Monocyte, IFN signature"), 
                     bioRep = "orig.ident", padj_cutoff = 0.01, lfcCut = 0.58, minCells = 25, 
                     outDir = "./output/fig5/btwnClusDEG/", title = "M-MDSC vs Clusters 0, 2, 3, 4, 6, 10 & 11", 
                     idents.1_NAME= "M-MDSC", idents.2_NAME = "Monocyte", 
                     returnVolc = T, doLinDEG = F, paired = T, addLabs = c("S100A5", "BPT", "TFPT", "RACK1", "DLA-DRA", "DLA-DQA1", "HLA-DQB2", "NRG","F13A1", "MT2A", "FN1"),
                 lowFilter = T, dwnSam = F, setSeed = 12, topn=c(20,10)
                    ) 
p <- p[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave(paste("./output/", outName, "/", outName, "_f_mMDSCvMono.png", sep = ""), width = 7, height = 7)


#################################################
### END FIGURE 5  | BEGIN SUPP BCELL ANALYSIS ###
#################################################

seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")

highlight <- c("B cell", "Plasma cell")
colArray <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
colArray <- colArray %>% mutate(high_col = ifelse(majorID %in% highlight, color,"grey"))

umapHighLight <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank())


seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "B cell" | majorID ==  "Plasma cell" | majorID ==  "Granulocyte" | majorID ==  "DC"
                 ) 

indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "bcell_hVoWadj_CLEAN_introns", preSub = T,
         vars.to.regress = c("percent.mt", "percent.pal")
         )

seu.obj <- readRDS(file = "./output/s2/bcell_hVoWadj_CLEAN_introns_S2.rds")

clusTree(seu.obj = seu.obj, outName = "bcell_hVoWadj_CLEAN_introns", dout = "./output/clustree/", test_dims = c(45,40,35,30,25), algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "bcell_hVoWadj_CLEAN_introns", 
           final.dims = 30, final.res = 0.6, returnFeats = F, saveRDS = F, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.6,
            n.neighbors = 75, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )


##################################
### BEGIN SUPP B CELL ANALYSIS ###
##################################

### Load in processed data
seu.obj <- readRDS(file = "./output/s3/bcell_hVoWadj_CLEAN_introns_res0.6_dims30_S3.rds")

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")

#load in colors
colArray <- read.csv("./bcell_hVoWadj.csv")
outName <- "bcell"


vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "bcell_hVoWadj_CLEAN_introns",
                     outDir = "./output/bcell/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub

p <- dotPlotBY_TYPE(seu_obj = seu.obj, pwdTOvilnCSVoutput = "./output/bcell/viln/bcell_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave(paste("./output/", outName, "/", outName, "_supp_bigDots.png", sep = ""), width = 20, height = 12)


### Fig supp: plot inital cluster umap
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./bcell_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "pctID")
seu.obj.sub <- subset(seu.obj, subset = pctID != "exclude")
seu.obj.sub$pctID <- droplevels(seu.obj.sub$pctID)

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              cols = colArray$colour,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8,labCol = colArray$labCol)
ggsave(paste("./output/", outName, "/", outName, "_supp_bcell_umap.png", sep = ""), width = 7, height = 7)


### Fig supp: create custom legend                             
colArray$title <- "B cell subtypes"
leg <- cusLeg(legend = colArray, colz = 3, rowz = 3, clusLabel = "clusterID_sub", legLabel = "pctID", colorz = "colour",
                   groupLabel = "title", groupBy = "title", sortBy = "clusterID_sub", labCol = "labCol", headerSize = 6,
                   cellHeader = T, bump = 0, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 2, compress_x = 0.8, spaceBtwnCols =c(0.5,0.4)
             )
ggsave(paste("./output/", outName, "/", outName, "_supp_bcell_legend_umap.png", sep = ""), width = 9, height = 3)


### Fig supp: create sankey plot to show how things changed                           
colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
saveCells <- c("Bcell")
colArray_all <- colArray_all %>% filter(colorID %in% saveCells)
                             
p <- sankeyPlot(seu_obj = seu.obj.sub, new.ident = "clusterID_sub", old.ident = "clusterID", old.colorz = colArray_all$color,
                       new.colorz = colArray$colour, old.labCol = colArray_all$labCol, new.labCol = colArray$labCol, flowCol = "grey"
                    )

ggsave(paste("./output/", outName, "/", outName, "_supp_bcell_sankey.png", sep = ""), width = 7, height = 10)


### Fig supp: b cell feat plots
features <- c("MS4A1","VPREB3","IGF1R", "PAX5","IGHM",
              "CXCR4","IRF5", "CD80", "CD40",
              "NAPSA","RAG1","DNTT"
              )
fig4b <- prettyFeats(seu.obj = seu.obj, nrow = 4, ncol = 4, features = features, color = "black", order = F,bottomLeg = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp_bcell_feats.png", sep = ""), width = 12, height = 9)


### Fig supp: b cell freq plots
p <- freqPlots(seu.obj.sub, method = 1, nrow= 2, groupBy = "pctID", legTitle = "Cell source",refVal = "name",comp = "cellSource", no_legend = T,
               namez = "name", 
               colz = "colz"
              ) + theme(axis.text.x = element_blank()) 

ggsave(paste("./output/", outName, "/", outName, "_supp_bcell_freqPlots.png", sep = ""), width = 6, height = 4)


### Fig supp: b cell linDEG
linDEG(seu.obj = seu.obj.sub, threshold = 1, thresLine = F, groupBy = "pctID", comparision = "cellSource", outDir = "./output/bcell/linDEG/", outName = "bcell", cluster = NULL, labCutoff = 20,
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, pValCutoff = 0.01, flipLFC = F, saveGeneList = F, addLabs = ""
                  )


################################################################
### END SUPP BCELL ANALYSIS | BEGIN HUAMN HOMOLOGY ANALSYSIS ###
################################################################

#load in data downloaded from https://zenodo.org/record/4021967
#then subset on the heaty humans in datset
seu.obj.human <- readRDS(file = "./blish_covid.seu.rds")
seu.obj.human <- subset(seu.obj.human,
                        subset = 
                        Status ==  "Healthy") 

#load in final canine healthy dataset -- avaliable as "GSE225599_final_dataSet_H.rds.gz" at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225599
seu.obj.dog <- readRDS(file = "./output/s3/final_dataSet_H.rds")

seu.obj.crossSpecies <- merge(seu.obj.dog,seu.obj.human)
indReClus(seu.obj = seu.obj.crossSpecies, outDir = "./output/s2/", subName = "cross_Species_2000", preSub = T, nfeatures = 2000,
                      vars.to.regress = "percent.mt"
                       )

seu.obj <- readRDS(file = "./output/s2/cross_Species_2000_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "cross_Species_2000-5", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "cross_Species_2000", final.dims = 45, final.res = 0.8, stashID = "clusterID2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


### Load in integrated dataset
seu.obj <- readRDS("./output/s3/cross_Species_2000_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj$cellSource2 <- ifelse(is.na(seu.obj$cell.type), "Canine", "Human")    
seu.obj$type <- ifelse(is.na(seu.obj$cell.type), paste("can_",seu.obj$celltype.l3,sep=""), paste("hu_",seu.obj$cell.type.fine,sep=""))    
seu.obj$mergedAnno <- ifelse(is.na(seu.obj$cell.type), paste0(seu.obj$celltype.l3), paste0(seu.obj$cell.type.fine))
outName <- "human_dog_naive"


### Fig supp 9a: show integrated umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "mergedAnno",
              split.by = "cellSource2",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_umap_humanCellType.png", sep = ""), width = 14, height = 7)


### Fig supp 9b: Create dendrogram
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$integrated@data)) #use sct count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

# avg_expression <- expression %>% group_by(anno_merge) %>% summarise_each(mean) %>% as.data.frame()
avg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()

rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL
M <- (1- cor(t(avg_expression),method="pearson"))/2

outfile <- paste("./output/", outName, "/", outName, "_phylo.png", sep = "")
png(file = outfile, width=4000, height=4000, res=400)
par(mfcol=c(1,1))
## dd <- dist(M)
hc <- hclust(as.dist(M),method="complete")
p <- plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE,tip.col = "black")
dev.off()

node.df <- as.data.frame(list("node" = c(73,72,68,64,78,79,67,66,60),
                              "colour" = rep(c("lightblue","lightgrey"),5)[1:9],
                              "cellType" = c("Mono/Neut","Gran","DC","B cell","Plasma cell", "Cycling", "CD4 T & other", "CD4","CD8 T & NK")))


ggtree(as.phylo(hc),layout = "fan") + geom_tiplab() + geom_hilight(data=node.df, mapping=aes(node=node, fill=colour)) + scale_fill_manual(values=c("lightblue","lightgrey")) +  NoLegend()+ xlim(NA, 0.75)

ggsave(paste("./output/", outName, "/", outName, "_ggTree.png", sep = ""), width = 9, height = 9) 


###############################################
### END HOMOLGY ANALYSIS | MAKE FINAL TABLE ###
###############################################

#Note: this code is not very efficient, but it is how I transfered labels and made the tables
#Note: if you wish to explore the cell percentatges, it wil likely be far more efficient if you download "GSE225599_final_dataSet_HvO.rds.gz" from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225599 and use the process file directly -- i.e. skip the noted out lines and load the .rds file in to recreate the table
#Note: feel free to submit an issue or email me directly with questions - dylan

### SKIP ###

# seu.obj.all <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "majorID", metaAdd = "color")
# colorData.all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)

# Idents(seu.obj.all) <- "majorID"

# seu.obj.all_laggers <- subset(seu.obj.all, idents =  c("gd T cell", "Cycling T cell","CD34+ Unclassified"))
# seu.obj.all_laggers$majorID <- droplevels(seu.obj.all_laggers$majorID)
# seu.obj.all_laggers$color <- droplevels(seu.obj.all_laggers$color)
# colorData.all <- colorData.all[colorData.all$majorID %in% c("gd T cell", "Cycling T cell","CD34+ Unclassified"),]

# seu.obj.cyto <- readRDS(file = "./output/s3/cytotoxic_noGD_hVoWadj_CLEAN_introns_res0.8_dims35_dist0.3_neigh30_S3.rds")
# seu.obj.cyto <- loadMeta(seu.obj = seu.obj.cyto, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "classification")
# seu.obj.cyto <- loadMeta(seu.obj = seu.obj.cyto, metaFile = "./cyto_hVoWadj.csv", groupBy = "classification", metaAdd = "colour")
# colorData.cyto <- read.csv("./cyto_hVoWadj.csv", header = T)

# seu.obj.helper <- readRDS(file = "./output/s3/helper_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.2_neigh20_S3.rds")
# seu.obj.helper <- loadMeta(seu.obj = seu.obj.helper, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
# seu.obj.helper <- loadMeta(seu.obj = seu.obj.helper, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "majorID_sub", metaAdd = "colour")
# colorData.helper <- read.csv("./helper_hVoWadj_forLeg_020623.csv", header = T)

# seu.obj.myeloid <- readRDS(file = "./output/s3/myeloid_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.35_neigh20_S3.rds")
# seu.obj.myeloid <- loadMeta(seu.obj = seu.obj.myeloid, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")    
# seu.obj.myeloid <- loadMeta(seu.obj = seu.obj.myeloid, metaFile = "./coldaf_myeloid_new.csv", groupBy = "majorID_sub", metaAdd = "colour")    
# colorData.myeloid <- read.csv("./coldaf_myeloid_new.csv", header = T)

# seu.obj.bcell <- readRDS(file = "./output/s3/bcell_hVoWadj_CLEAN_introns_res0.6_dims30_S3.rds")
# seu.obj.bcell <- loadMeta(seu.obj = seu.obj.bcell, metaFile = "./bcell_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "pctID")
# seu.obj.bcell <- loadMeta(seu.obj = seu.obj.bcell, metaFile = "./bcell_hVoWadj.csv", groupBy = "pctID", metaAdd = "colour")
# colorData.bcell <- read.csv("./bcell_hVoWadj.csv", header = T)

# classifications <- c(seu.obj.helper$majorID_sub,
#                      seu.obj.all_laggers$majorID,
#                      seu.obj.myeloid$majorID_sub,
#                      seu.obj.cyto$classification,
#                      seu.obj.bcell$pctID)
        
# levels(classifications) <- ifelse(levels(classifications) == "exclude", "Low quality", levels(classifications))
# levels(classifications) <- ifelse(levels(classifications) == "Poor quality", "Low quality", levels(classifications))

# seu.obj.all <- AddMetaData(seu.obj.all, classifications, col.name = "finalID")
        
# colourz <- as.data.frame(c(levels(seu.obj.helper$colour), 
#                            levels(seu.obj.all_laggers$color),
#                            levels(seu.obj.myeloid$colour),
#                            levels(seu.obj.cyto$colour),
#                            levels(seu.obj.bcell$colour)
#                           )
#                         )
# colnames(colourz) <- "colour"
        
# colourz$classification <- c(levels(seu.obj.helper$majorID_sub), 
#                             levels(seu.obj.all_laggers$majorID), 
#                             levels(seu.obj.myeloid$majorID_sub), 
#                             levels(seu.obj.cyto$classification), 
#                             levels(seu.obj.bcell$pctID)
#                            )
# colourz$classification <- ifelse(colourz$classification == "exclude", "Low quality", colourz$classification)
# colourz$classification <- ifelse(colourz$classification == "Poor quality", "Low quality", colourz$classification)

# colourz$classificationFinal <- colourz$classification

# colourz$majorGroup <- c(rep("CD4 T cell", length(levels(seu.obj.helper$majorID_sub))), 
#                         rep("Miscellaneous", length(levels(seu.obj.all_laggers$majorID))), 
#                         rep("Myeloid", length(levels(seu.obj.myeloid$majorID_sub))), 
#                         rep("CD8/NK cell", length(levels(seu.obj.cyto$classification))), 
#                         rep("B cell", length(levels(seu.obj.bcell$pctID)))
#                        )

# rm(list = c("seu.obj.helper", "seu.obj.myeloid", "seu.obj.all_laggers", "seu.obj.cyto", "seu.obj.bcell"))

# #hack to fix annoying error where labCol cannot be brought over
# colnames(colorData.all)[4] <- "colour"

# labColData <- as.data.frame(mapply(c, colorData.helper[,c("majorID_sub","labCol")],
#                     colorData.all[,c("majorID","labCol")],
#                     colorData.myeloid[,c("majorID_sub","labCol")],
#                     colorData.cyto[,c("classification","labCol")],
#                     colorData.bcell[,c("pctID","labCol")]))
        
# labColData <- labColData[!duplicated(labColData$majorID_sub),]
# colnames(labColData)[1] <- "classification"
# colourz <- colourz %>% left_join(labColData, by = "classification", keep = F)
        
# #clean things up for plotting and downstream table
# getSomeOrder <- as.data.frame(levels(seu.obj.all$finalID)) %>% left_join(colourz, by = c("levels(seu.obj.all$finalID)" = "classificationFinal"), keep = F)
# getSomeOrder <- getSomeOrder[!duplicated(getSomeOrder$classification),]

# #hardcode to remove low quality cells
# clusterID_final <- table(seu.obj.all$finalID) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
#         mutate(clusterID_final=case_when(row_number() < 15 ~ row_number()-1,
#                                          row_number() == 15 ~ 100,
#                                          row_number() > 15 ~ row_number()-2)) %>% arrange(clusterID_final) 

# newID <- clusterID_final$clusterID_final
# names(newID) <- clusterID_final$Var1

# Idents(seu.obj.all) <- "finalID"
# seu.obj.all <- RenameIdents(seu.obj.all, newID)
# table(Idents(seu.obj.all))
# seu.obj.all$clusterID_final <- Idents(seu.obj.all)

# getSomeOrder <- clusterID_final %>% left_join(getSomeOrder, by = c("Var1" = "levels(seu.obj.all$finalID)"), keep = F) %>% arrange(clusterID_final)

# getSomeOrder[getSomeOrder$clusterID_final == 100,][7] <- NA

### RESUME ###

########## >>>>>>>>>> IF SKIP RUN THE FOLLOWING!!!!  <<<<<<<<< ########


seu.obj.all <- readRDS("./output/s3/final_dataSet_HvO.rds") #path to GSE225599_final_dataSet_HvO.rds.gz downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225599
seu.obj.all$finalID <- seu.obj.all$celltype.l3 
getSomeOrder <- read.csv("majorGroups.csv", header = T)
getSomeOrder <- getSomeOrder[-37,]


#note pretty colors do not follow if you skip
###################################################
### END LABEL/COLOR TRANSFER | MAKE FINAL TABLE ###
###################################################
outName <- "tablesANDsupp8"

### Fig supp 8a: final umap
plot <- DimPlot(seu.obj.all, 
              reduction = "umap", 
              group.by = "clusterID_final",
              cols = getSomeOrder$colour,
              pt.size = 0.25,
              label = T,
              label.box = T
)

fig7c <- cusLabels(plot = plot, shape = 21, labCol = getSomeOrder$labCol, size = 8, 
                   alpha = 0.8, rm.na = T, nudge_x = c(rep(0,27),-0.45,rep(0,8)), nudge_y = c(rep(0,27),1.15,rep(0,8))) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_supp_8a_finalUMAP.png", sep = ""), width = 7, height = 7) 


#remove Poor quality cells from table        
seu.obj.all <- subset(seu.obj.all, subset = finalID != "Low quality")

#BONUS: do age matched stats
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
# seu.obj.all <- subset(seu.obj.all,
#                   subset = 
#                   name !=  "OS_1" & name !=  "OS_2" & name !=  "OS_4" & name !=  "OS_5" & name !=  "OS_10"
#                  ) 
# getSomeOrder <- read.csv("majorGroups_3.csv", header = T)

seu.obj.all$finalID <- droplevels(seu.obj.all$finalID)

tableData <- as.data.frame(table(seu.obj.all$finalID, seu.obj.all$orig.ident))

statData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   pct = Freq/samSize*100,
                                                   cellSource = ifelse(grepl("ealthy",Var2),"Healthy","OS")) %>% group_by(Var1) %>% as.data.frame()

getSomeOrder <- getSomeOrder %>% mutate(newName = paste(clusterID_final,Var1,sep="_")) 
statData <- statData %>% left_join(getSomeOrder[,c("Var1","newName")], by = "Var1")
        
statData$newName <- factor(statData$newName,levels = getSomeOrder$newName)

tstRes <- ggpubr::compare_means(pct ~ cellSource, statData, group.by = "Var1") %>% as.data.frame()

p <- ggplot(statData, aes_string(y = "pct", x = "cellSource")) + 
        labs(x = NULL, y = "Proportion (%)") + 
        theme_bw() + 
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(), 
              strip.background = element_rect(fill = NA, color = NA), 
              strip.text = element_text(face = "bold"), 
              axis.ticks.x = element_blank(), 
              axis.text = element_text(color = "black")         
             )
        
pi <- p + facet_wrap("newName", scales = "free_y", nrow = 4) + 
        guides(fill = "none") + 
        geom_boxplot(aes_string(x = "cellSource"), alpha = 0.25, outlier.color = NA) + 
        geom_point(size = 2, position = position_jitter(width = 0.25),
                   aes_string(x = "cellSource", y = "pct", color = "newName")) +
         ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3) + 
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        theme(panel.grid.major = element_line(color = "grey", size = 0.25),
              text = element_text(size = 12)
             ) + NoLegend() + scale_colour_manual(values = getSomeOrder$colour)

ggsave(paste("./output/", outName, "/", outName, "_supp_8b_stats.png", sep = ""), width = 24, height = 12) 


coorData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   pct = Freq/samSize*100,
                                                   cellSource = ifelse(grepl("ealthy",Var2),"Healthy","OS")) %>% as.data.frame()
# #BONUS: more age-assocaited analysis
# ageData <- read.csv("refColz.csv")
# ageData$orig.ident <- factor(ageData$orig.ident)
# ageData$cellSource <- ifelse(grepl("ealthy", ageData$orig.ident), "Healthy", "Osteosarcoma")

# ageRes <- ggpubr::compare_means(age ~ cellSource, ageData) %>% as.data.frame()

# p <- ggplot(ageData, aes_string(y = "age", x = "cellSource")) + 
#         labs(x = NULL, y = "Age (years)") + 
#          geom_boxplot(aes_string(x = "cellSource"), alpha = 0.25, outlier.color = NA) + 
#         geom_point(size = 2, position = position_jitter(width = 0.25),
#                    aes_string(x = "cellSource", y = "age")) +
#         theme_bw() + 
#         theme(panel.grid.minor = element_blank(),
#               panel.grid.major = element_blank(), 
#               strip.background = element_rect(fill = NA, color = NA), 
#               strip.text = element_text(face = "bold"), 
#               axis.ticks.x = element_blank(), 
#               axis.text = element_text(color = "black")         
#              )
# ggsave(paste("./output/", outName, "/", outName, "_supp_age_Stats.png", sep = ""), width = 7, height = 7) 


# coorData <- coorData %>% left_join(ageData[,-5], by = c("Var2" = "orig.ident"))

# p <- ggplot(coorData, aes(x=age, y=pct)) + 
#                 stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
#                 #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
#                 geom_point(aes(color = cellSource.x)) +
#                 #scale_x_continuous(limits = c(0, 100),expand = c(0, 0)) +
#                 #scale_y_continuous(limits = c(0, 100),expand = c(0, 0)) +
#                 labs(x = "Age (years)", y = "Percentage") +
#                 guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
#                 #geom_text(x = 75, y = 5, label = as.character(as.expression(eq)), parse = TRUE) +
#                 theme(panel.background = element_rect(fill = "transparent",colour = NA),
#                       plot.background = element_rect(fill = "transparent",colour = NA),
#                       legend.background = element_rect(fill = "transparent",colour = NA),
#                       legend.key = element_rect(fill = "transparent",colour = NA),
#                       panel.grid.major = element_line(color = "gray"), 
#                       panel.grid.minor = element_line(color = "gray"),#axis.text = element_blank(), 
#                       axis.ticks = element_blank(),
#                       axis.title = element_text(size= 14),
#                       #plot.title = element_text(hjust = 0.02),
#                       #plot.title.position = "plot",
#                       plot.title = element_blank(),
#                       title = element_text(size= 16),
#                       axis.line = element_blank(),
#                       panel.border = element_rect(color = "black",
#                                                   fill = NA,
#                                                   size = 2)
#                       ) 
                
# pi <- p + facet_wrap("Var1", scales = "free_y", nrow = 8) + stat_cor() + theme(plot.background = element_rect(fill = "white"))

# ggsave(paste("./output/", outName, "/", outName, "_supp_cellByAgeCoor.png", sep = ""), width = 12, height = 18)


allData <- statData %>% group_by(Var1) %>% summarize(avg = mean(pct),
                                  std = sd(pct)
                                 ) %>% mutate(cellSource = "All") %>% relocate(cellSource, .after = Var1) %>% as.data.frame()

# OPTIONAL: manually add major groups
# write.csv(getSomeOrder, "majorGroups.csv", row.names= F) # If you did not skip the ealrier lines, you will have to save the getSomeOrder obj and modify it to add major groups
# groupData <- read.csv("majorGroups.csv", header = T) # then load it back in
        

#################################################
### MAKE FINAL TABLE | EXTRACT DATA FOR TABLE ###
#################################################

headerData_all_stat <- tableData %>% left_join(groupData[,c("Var1","majorGroup")], by = "Var1") %>% group_by(majorGroup,Var2) %>%  summarise(Freq = sum(Freq)) %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   pct = Freq/samSize*100)

headerData_all <- headerData_all_stat %>% group_by(majorGroup) %>% summarize(avg = mean(pct),
                                                                                                            std = sd(pct)
                                                                                                           ) %>% mutate(cellSource = "All") %>% relocate(cellSource, .after = majorGroup) %>% as.data.frame()
                
tableData$cellSource <- ifelse(grepl("ealthy", tableData$Var2), "Healthy", "Osteosarcoma")
subData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   pct = Freq/samSize*100) %>% group_by(Var1, cellSource) %>% summarize(avg = mean(pct),
                                                                                                            std = sd(pct)
                                                                                                           ) %>% as.data.frame()
        
        
headerData_sub_stat <- tableData %>% left_join(groupData[,c("Var1","majorGroup")], by = "Var1") %>% group_by(majorGroup,Var2) %>%  summarise(Freq = sum(Freq)) %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   pct = Freq/samSize*100,
                                                   cellSource = ifelse(grepl("ealthy", Var2), "Healthy", "Osteosarcoma"))
        
headerData_sub <- headerData_sub_stat %>% group_by(majorGroup, cellSource) %>% summarize(avg = mean(pct),
                                                                                                            std = sd(pct)
                                                                                                           ) 

compiledData <- rbind(subData,allData) %>% arrange(Var1)
        
compiledData <- reshape(compiledData, idvar = "Var1", timevar = "cellSource", direction = "wide") %>% filter(avg.Healthy != 0) %>% 
        column_to_rownames(var="Var1") %>% round(.,1) %>% 
        transmute(All = paste(avg.All," ± ",std.All,sep=""),
                  Healthy = paste(avg.Healthy," ± ",std.Healthy,sep=""),
                  Osteosarcoma = paste(avg.Osteosarcoma," ± ",std.Osteosarcoma,sep="")
              )

compiledData <- compiledData %>% rownames_to_column() %>% left_join(tstRes[,c("Var1","p.format")], by = c("rowname" = "Var1")) %>% remove_rownames %>% column_to_rownames(var="rowname")

write.csv(compiledData, file = paste("./output/", outName, "/", outName, "compiledTableData.csv", sep = ""))


######################################################
### MAKE FINAL TABLE | EXTRACT MORE DATA FOR TABLE ###
######################################################

compiledData <- rbind(headerData_sub,headerData_all) %>% arrange(majorGroup) %>% as.data.frame()

statData <- rbind(headerData_sub_stat,headerData_all_stat) %>% arrange(majorGroup) %>% as.data.frame()

tstRes <- ggpubr::compare_means(pct ~ cellSource, statData, group.by = "majorGroup") %>% as.data.frame()

compiledData <- reshape(compiledData, idvar = "majorGroup", timevar = "cellSource", direction = "wide") %>% remove_rownames() %>%
        column_to_rownames(var="majorGroup") %>% round(.,1) %>% 
        transmute(All = paste(avg.All," ± ",std.All,sep=""),
                  Healthy = paste(avg.Healthy," ± ",std.Healthy,sep=""),
                  Osteosarcoma = paste(avg.Osteosarcoma," ± ",std.Osteosarcoma,sep="")
              )

compiledData <- compiledData %>% rownames_to_column() %>% left_join(tstRes[,c("majorGroup","p.format")], by = c("rowname" = "majorGroup")) %>% remove_rownames %>% column_to_rownames(var="rowname")
        
write.csv(compiledData, file = paste("./output/", outName, "/", outName, "header_compiledTableData.csv", sep = ""))


############################################################
### NOW YOU SHOULD HAVE 2 .csv FILES TO CREATE THE TABLE ###
############################################################

#####################################################
### NEXT UP -- EVALUTE COORELATION WITH FLOW DATA ###
#####################################################

#prepare scRNA data
groupData <- read.csv("majorGroupsForCoor.csv", header = T)
sclong <- tableData %>% left_join(groupData[,c("Var1","majorGroup")], by = "Var1") %>% group_by(majorGroup,Var2) %>%  summarise(Freq = sum(Freq)) %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   value = Freq/samSize*100) %>% as.data.frame()
sclong <- sclong[sclong$majorGroup != "Other",]
sclong <- sclong[,-c(3,4)]
colnames(sclong) <- c("cellType","Dog","Percent")

#read and clean flow data
decode <- read.csv("randID.csv")
decode$randID <- as.integer(decode$randID)     

flow_df <- read.csv("scRNA_pbmc_flow.csv")
flow_df <- flow_df %>% left_join(decode, by = c("randID"))
flow_df <- flow_df[,-c(1,5)]
flowlong <- melt(flow_df, id.vars = c("dog"))

colnames(flowlong) <- c("Dog","cellType","Percent")

flowlong$cellType <- gsub("tcell", "T cell", flowlong$cellType)
flowlong$cellType <- gsub("myeloid", "Myeloid", flowlong$cellType)
flowlong$cellType <- gsub("bcell", "B cell", flowlong$cellType)

plotData <- inner_join(flowlong, sclong, by = c("Dog", "cellType"))

lmModel <- lm(data = plotData, Percent.y ~ Percent.x)

eq <- substitute(italic(y) == a + b ~italic(x)*","~~italic(R)^2~"="~r2, 
         list(a = format(unname(coef(lmModel)[1]), digits = 2),
              b = format(unname(coef(lmModel)[2]), digits = 2),
             r2 = format(summary(lmModel)$r.squared, digits = 3)))

### Fig supp 8c: make the coor plot with flow -- NOTE: warning is b/c there is one NA in the df
ggplot(plotData, aes(x=Percent.x, y=Percent.y)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
                #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point(aes(color = cellType)) +
                scale_x_continuous(limits = c(0, 100),expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, 100),expand = c(0, 0)) +
                labs(x = "Percentage by flow cytometry", y = "Percentage by scRNA") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                geom_text(x = 75, y = 5, label = as.character(as.expression(eq)), parse = TRUE) +
                theme(panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.background = element_rect(fill = "transparent",colour = NA),
                      legend.key = element_rect(fill = "transparent",colour = NA),
                      panel.grid.major = element_line(color = "gray"), 
                      panel.grid.minor = element_line(color = "gray"),#axis.text = element_blank(), 
                      axis.ticks = element_blank(),
                      axis.title = element_text(size= 14),
                      #plot.title = element_text(hjust = 0.02),
                      #plot.title.position = "plot",
                      plot.title = element_blank(),
                      title = element_text(size= 16),
                      axis.line = element_blank(),
                      panel.border = element_rect(color = "black",
                                                  fill = NA,
                                                  size = 2)
                      ) + scale_colour_manual(values=c("#645A9F", "#FF755F", "#009DA5"))
                
ggsave(paste("./output/", outName, "/", outName, "_supp_8c_FlowCoorPlot.png.png", sep = ""), width = 7, height = 7)


### IF YOU SKIPPED THE LINES AND ARE DOING A REPRODUCIBLE RUN -- RECCOMEND SAVING THE RDS OBJ NOW
#saveRDS(seu.obj.all, file = "./output/processedData.rds")
#seu.obj.all2 <- readRDS(file = "./output/processedData.rds")


###### generate feat list #####
seu.obj.all <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds") # can just load in from GEO

vilnPlots(seu.obj = seu.obj.all, groupBy = "celltype.l3", numOfFeats = 24, outName = "HvO_cell.l3",
                      outDir = "./output/viln_finalID_HvO_cell.l3/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

vilnPlots(seu.obj = seu.obj.all, groupBy = "celltype.l2", numOfFeats = 24, outName = "HvO_cell.l2",
                      outDir = "./output/viln_finalID_HvO_cell.l2/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

vilnPlots(seu.obj = seu.obj.all, groupBy = "celltype.l1", numOfFeats = 24, outName = "HvO_cell.l1",
                      outDir = "./output/viln_finalID_HvO_cell.l1/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )



seu.obj_h <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/s3/final_dataSet_H.rds") # can just load in from GEO
vilnPlots(seu.obj = seu.obj_h, groupBy = "celltype.l1", 
          numOfFeats = 24, outName = "H_cell.l1",
                      outDir = "./output/viln_finalID_H_cell.l1/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

vilnPlots(seu.obj = seu.obj_h, groupBy = "celltype.l2", 
          numOfFeats = 24, outName = "H_cell.l2",
                      outDir = "./output/viln_finalID_H_cell.l2/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

vilnPlots(seu.obj = seu.obj_h, groupBy = "celltype.l3", 
          numOfFeats = 24, outName = "H_cell.l3",
                      outDir = "./output/viln_finalID_H_cell.l3/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

####################################################################
### EXTRACT DATA TO MAKE SUPP TABLE 3 and TABLES 1 & 2 ON GITHUB ###
####################################################################
featOut.df <- read.csv("./output/viln_finalID_H_cell.l1/H_cell.l1_gene_list.csv")
featOutPROCESSED.df <- as.data.frame(lapply(split(featOut.df$gene,featOut.df$cluster), head, n = 50))
write.csv(featOutPROCESSED.df, file = paste("./output/", outName, "/", outName, "H_cell.l1_cleanGeneList.csv", sep = ""))

featOut.df <- read.csv("./output/viln_finalID_H_cell.l2/H_cell.l2_gene_list.csv")
featOutPROCESSED.df <- as.data.frame(lapply(split(featOut.df$gene,featOut.df$cluster), head, n = 50))
write.csv(featOutPROCESSED.df, file = paste("./output/", outName, "/", outName, "H_cell.l2_cleanGeneList.csv", sep = ""))

featOut.df <- read.csv("./output/viln_finalID_H_cell.l3/H_cell.l3_gene_list.csv")
featOutPROCESSED.df <- as.data.frame(lapply(split(featOut.df$gene,featOut.df$cluster), head, n = 50))
write.csv(featOutPROCESSED.df, file = paste("./output/", outName, "/", outName, "H_cell.l3_cleanGeneList.csv", sep = ""))



###IF SKIPPED PREVIOUS LINES THIS CODE WILL CLEAN UP THE OBJ TO GET IT TO THE STATE UPLOADED TO GEO

##################################################
### LAST -- CLEAN UP OBJECTS FOR UPLOAD TO GEO ###
##################################################

seu.obj.all@meta.data <- seu.obj.all@meta.data[,!grepl("DF|pANN", colnames(seu.obj.all@meta.data))]
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./majorGroups.csv", groupBy = "finalID", metaAdd = "middleID")
colnames(seu.obj.all@meta.data)[23] <- "celltype.l1"
colnames(seu.obj.all@meta.data)[25] <- "celltype.l3"
colnames(seu.obj.all@meta.data)[27] <- "celltype.l2"
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj.all@meta.data <- seu.obj.all@meta.data[,-c(7,8,9,11,15)]

#reorder for readability
seu.obj.all@meta.data <- seu.obj.all@meta.data[,c(seq(1,12),16,13,17,23,15,14,21,18,22,20)]

saveRDS(seu.obj.all, file = "./output/s3/final_dataSet_HvO.rds")

seu.final.H <- subset(seu.obj.all,
                  subset = 
                  cellSource ==  "Healthy") 

saveRDS(seu.final.H, file = "./output/s3/final_dataSet_H.rds")

#remove additional meta slots for cell browser
seu.obj < readRDS("./output/s3/final_dataSet_H.rds")
seu.obj@meta.data <- seu.obj@meta.data[,-c(5,6,13,14,15)]
saveRDS(seu.obj, file = "./output/s3/final_dataSet_H.rds")

#clean helpers
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/helper.rds")
seu.obj@meta.data <- seu.obj@meta.data[,-c(5,6,13,14,15,16)]

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "clusterID_sub_clean")
seu.obj$clusterID_sub_clean <- factor(x = Idents(seu.obj), levels = sort(as.numeric(as.character(levels(seu.obj)))))
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")

#create subset object to remove por qulity cells from downstream analysis
seu.obj <- subset(seu.obj, subset = majorID_sub != "exclude")
seu.obj$majorID_sub <- droplevels(seu.obj$majorID_sub)
colnames(seu.obj@meta.data)[17] <- "celltype.l3"
saveRDS(seu.obj, "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/helper.rds")

#clean myeloid
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/myeloid.rds")
seu.obj@meta.data <- seu.obj@meta.data[,-c(5,6,13,14,15,16,17)]

#load in key meta data and remove poor quality cells for downstream analysis
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "clusterID_sub_clean")
seu.obj$clusterID_sub_clean <- factor(x = Idents(seu.obj), levels = sort(as.numeric(as.character(levels(seu.obj)))))

seu.obj <- subset(seu.obj, subset = majorID_sub != "exclude")
seu.obj@meta.data <- seu.obj@meta.data[,-c(17,18,20)]
colnames(seu.obj@meta.data)[16] <- "celltype.l3"
saveRDS(seu.obj, "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/myeloid.rds")

#clean bcell
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/bcell.rds")
seu.obj@meta.data <- seu.obj@meta.data[,-c(5,6,13,14,15,16)]

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./bcell_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "pctID")
seu.obj.sub <- subset(seu.obj, subset = pctID != "exclude")
seu.obj.sub$pctID <- droplevels(seu.obj.sub$pctID)
colnames(seu.obj@meta.data)[16] <- "celltype.l3"
saveRDS(seu.obj, "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/bcell.rds")


#clean cyto
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/cytotoxic.rds")
seu.obj@meta.data <- seu.obj@meta.data[,-c(5,6,13,14,15,16)]


#load in metadata
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
#reorder majorID factor levels
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("Naive CD8", "Effector CD8","Memory CD8",
                                                      "NK/NKT cell", "DN T cell", "CD8 gd T cell"))

colnames(seu.obj@meta.data)[16] <- "celltype.l3"
saveRDS(seu.obj, "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/inputIntrons/forGEO/cytotoxic.rds")


#clean colors
df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/sheets/majorGroups.csv")
df <- df[-37,]
df <- df[,c("colour","finalID")]
write.csv(df, "./output/cellbrowser/celltype.l3_colours.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/sheets/majorGroups.csv")
df <- df[-37,]
df <- df[,c("colour","middleID")]
write.csv(df, "./output/cellbrowser/celltype.l2_colours.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/sheets/majorGroups.csv")
df <- df[-37,]
df <- df[,c("colour","majorGroup")]
write.csv(df, "./output/cellbrowser/celltype.l1_colours.csv", row.names=F)


##########################################
### COMPILE -- SUPP DATA 3 (DEG SHEET) ###
##########################################
files <- list.files(path = c("./output/fig3/btwnClusDEG/","./output/fig4/btwnClusDEG/","./output/fig5/btwnClusDEG/"), pattern=".csv", all.files=FALSE,
                        full.names=T)

df.list <- lapply(files, read.csv)

compiled.df <- do.call(rbind, df.list) %>% as.data.frame()
compiled.df$gs_base <- gsub("OTHER_T_CELLS","OTHER_CD4T_CELLS",compiled.df$gs_base)
write.csv(compiled.df, file = "./output/supplemental_data_3.csv")


##############################################
### CLEAN GENE LISTS FOR UCSC CELL BROWSER ###
##############################################
df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_HvO_cell.l1/HvO_cell.l1_gene_list.csv")[,c("cluster","gene","p_val_adj","avg_log2FC")]
write.csv(df, "./output/cellbrowser/HvO_cell.l1_gene_list.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_HvO_cell.l2/HvO_cell.l2_gene_list.csv")[,c("cluster","gene","p_val_adj","avg_log2FC")]
write.csv(df, "./output/cellbrowser/HvO_cell.l2_gene_list.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_HvO_cell.l3/HvO_cell.l3_gene_list.csv")[,c("cluster","gene","p_val_adj","avg_log2FC")]
write.csv(df, "./output/cellbrowser/HvO_cell.l3_gene_list.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_H_cell.l1/H_cell.l1_gene_list.csv")[,c("cluster","gene","p_val_adj","avg_log2FC")]
write.csv(df, "./output/cellbrowser/H_cell.l1_gene_list.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_H_cell.l2/H_cell.l2_gene_list.csv")[,c("cluster","gene","p_val_adj","avg_log2FC")]
write.csv(df, "./output/cellbrowser/H_cell.l2_gene_list.csv", row.names=F)

df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_H_cell.l3/H_cell.l3_gene_list.csv")[,c("cluster","gene","p_val_adj","avg_log2FC")]
write.csv(df, "./output/cellbrowser/H_cell.l3_gene_list.csv", row.names=F)



### BONUS: final evlautation of age-driven differneces
seu.obj <- readRDS("./output/s3/final_dataSet_HvO.rds")

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "ageGroup")
seu.obj$ageGroup <- factor(seu.obj$ageGroup, levels = c("Healthy", "Middle-aged_OS", "Old_OS"))
seu.obj$ageGroup <- droplevels(seu.obj$ageGroup)
seu.obj$cellSource <- factor(seu.obj$cellSource)
seu.obj$name <- droplevels(seu.obj$name)

seu.obj$ageGroup <- factor(seu.obj$ageGroup, levels = c("Healthy", "Middle-aged_OS", "Old_OS"))
seu.obj.os <- subset(seu.obj,
                  subset = 
                  ageGroup !=  "Healthy"
                 ) 


#complete linDEG in each major cell type
linDEG(seu.obj = seu.obj.os, threshold = 0.5, thresLine = T, groupBy = "celltype.l3", comparision = "ageGroup", outDir = "./output/ageTest/", outName = "ageTest", colUp = "red", colDwn = "blue",subtitle = F)


seu.obj.oldvH <- subset(seu.obj,
                  subset = 
                  ageGroup !=  "Middle-aged_OS"
                 ) 


#complete linDEG in each major cell type
linDEG(seu.obj = seu.obj.oldvH, threshold = 0.5, thresLine = T, groupBy = "celltype.l3", comparision = "cellSource", outDir = "./output/ageTest/", outName = "HvOLD", colUp = "red", colDwn = "blue",subtitle = F)


seu.obj.MatchedvH <- subset(seu.obj,
                  subset = 
                  ageGroup !=  "Old_OS"
                 ) 


#complete linDEG in each major cell type
linDEG(seu.obj = seu.obj.MatchedvH, threshold = 0.5, thresLine = T, groupBy = "celltype.l3", comparision = "cellSource", outDir = "./output/ageTest/", outName = "HvMATCH", colUp = "red", colDwn = "blue",subtitle = F)


#logging data
system("cat /pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
sessionInfo()
