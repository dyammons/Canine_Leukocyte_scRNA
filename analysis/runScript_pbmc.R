#!/usr/bin/Rscript

#FOR PBMCs

#QC metrics
qc.df <- read.csv("compiledQCmetrics.csv")
mean(as.numeric(gsub(",","",qc.df[1:7,]$Mean.Reads.per.Cell)))

### hVoWadj

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### process and analyze healthy samples first ###
load10x(din = "./healthyIntrons/", dout = "./output/s1/healthy/", outName = "221005_h7_introns", testQC = F, nFeature_RNA_high = 4500, nFeature_RNA_low = 200, percent.mt_high = 25, nCount_RNA_high = 20000, nCount_RNA_low = 100) #update to include pal.pct calc (see seurat_1.R)

#integrate the data into one object
sctIntegrate(din = "./output/s1/healthy/", dout = "./output/s2/", outName = "221005_h7_regPal_wPalpct_introns", vars.to.regress = c("percent.mt", "percent.pal"), nfeatures = 2000, isolatePalRBC = T)

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

seu.obj <- readRDS("./output/s3/h7_CLEAN_introns_res1.6_dims45_S3.rds")

#load in meta data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")

colArray <- read.csv("./h7_majorID_wIntrons.csv")

#plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              cols = colArray$color,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol, nudge_x = c(rep(0,length(colArray$labCol)-1),0.75))
ggsave("./output/fig1a_umap_h7_CLEAN_introns.png")

#transfer TRDC data from ROS data
seu.obj.ros.cyto <- readRDS(file = "./ros_output/HvO_mt_r7_res1.5_dims50_S3.rds")

seu.obj <- AddMetaData(seu.obj, metadata = seu.obj.ros.cyto@assays$RNA@data[rownames(seu.obj.ros.cyto@assays$RNA@data) == "TRDC",], col.name = "TRDC")
seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 

features = c("TRDC")
p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 1, ncol = 1, features = features, color = "black", order = F) 
ggsave(plot = p, "./output/supp_trdc_allCells_featPlots_introns.png", width=3, height=3)


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
fig1b <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = colorz, order = F) 
ggsave("./output/fig1b_featPlots_h7_CLEAN_introns.png", width=9, height=15)

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./h7_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "colorID")

fig1c <- majorDot(seu.obj = seu.obj, groupBy = "colorID",
                  yAxis = c("Monocyte","DC","B cell","Neutrophil","CD4 T cell","CD8/NK cell","Eosinophil","Basophil","gd T cell","Cycling T cell","CD34+ Unclassified"),
                  features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
                               "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
                               "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
                               "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
                 )

ggsave("./output/fig1c_majorDot_h7_CLEAN_introns.png",width = 10)

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
ggsave("./output/fig1d_umap_bySample_h7_CLEAN_introns_introns.png")

#use singleR to use human reference for cell classification
singleR(seu.obj = seu.obj, outName = "h7_CLEAN_introns", clusters = "clusterID", outDir = "./output/singleR/")

#use conocial markers to ID cells
features <- c("CD3G","CD8A", "GZMA",
              "IL7R","CD4", "S100A12",
              "DLA-DRA","FLT3", "ANPEP",
              "MS4A1","JCHAIN", "IL5RA",
              "TOP2A","GATA3", "CD34")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = "black", order = F) 
ggsave(plot = p, "./output/featPlots_h7_CLEAN_introns.png", width=9, height=15)

p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "clusterID") +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1, byrow =T))

ggsave(file = './output/supp_stackedBar_h7.png', width = 8, height = 12)

markers <- FindMarkers(seu.obj, ident.1 = "CMP", lfc = 0.25, only.pos = T)

##### prepare healthy vs OS dataset #####

#load in 10x data and qc filter eeach sample
load10x(din = "./inputIntrons/", dout = "./output/s1/", outName = "221005_introns", testQC = F, nFeature_RNA_high = 4500, nFeature_RNA_low = 200, percent.mt_high = 25, nCount_RNA_high = 20000, nCount_RNA_low = 100) #update to include pal.pct calc (see seurat_1.R)

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


##### overview analysis of healthy vs OS dataset #####
#laod in processed data
seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")

#load in meta data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")

#transfer TRDC data from ROS data
seu.obj.ros.cyto <- readRDS(file = "./ros_output/sub_HvO_Cyto_mt_r7_res1.5_dims30_S3.rds")

features = c("CD3E", "CD4", "CD8A", "IL2RB",
             "CCL4", "CCL5", "CCR5", "CCR7",
             "NCR3", "KLRB1", "CD96","IL7R",
             "NFKBID","GZMA", "GZMB", "GZMK",
             "TRDC","RARRES2","KLRF1", "LGALS3", 
             "PI3","ID3", "GATA3", "IKZF2")
p <- prettyFeats(seu.obj = seu.obj.ros.cyto, nrow = 6, ncol = 4, features = features, color = "black", order = F) 
ggsave(plot = p, "./output/supp_featPlots_ros_introns.png", width=12, height=18)

#plot inital cluster umap
pi <- DimPlot(seu.obj.ros.cyto, 
        reduction = "umap", 
        group.by = "clusterID",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8)
ggsave("./output/umap_raw_ros_introns.png")

#get barcodes for Clusters 12 & 14
# seu.obj.ros.trdcClus <- subset(seu.obj.ros.cyto, subset = 
#                                clusterID ==  12 | clusterID == 14
#                               ) 

highlight <- colnames(seu.obj)[colnames(seu.obj) %in% WhichCells(seu.obj.ros.cyto, expression = TRDC > 0)]

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
ggsave("./output/umap_allCells_highlight_trdcPos.png")


#plot umap to evaluate sample distibution (Figure xx)
#NOTE: down sample data to get equal number of cells from each sample
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
fig2c <- formatUMAP(pi)
ggsave("./output/fig2b_umap_bySample_hVoWadj_CLEAN_introns_introns.png")

#use singleR to use human reference for cell classification
singleR(seu.obj = seu.obj, outName = "hVoWadj_CLEAN_introns", clusters = "clusterID", outDir = "./output/singleR/")



#use conocial markers to ID cells
features <- c("CD3G","CD8A", "GZMA",
              "IL7R","CD4", "S100A12",
              "DLA-DRA","FLT3", "ANPEP",
              "MS4A1","JCHAIN", "IL5RA",
              "TOP2A","GATA3", "CD34")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = "black", order = F) 
ggsave(plot = p, "./output/supp_featPlots_hVoWadj_CLEAN_introns.png", width=9, height=15)

#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = "all_hVoWadj_CLEAN_introns",
                     outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

#load in cell metadata
colArray <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
#plot funal colorized cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              cols = colArray$color,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )
fig2b <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol)
ggsave("./output/umap_final_hVoWadj_CLEAN_introns_introns.png")

#stash cell identies by majorID
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")

#reorder majorID factor levels
seu.obj$majorID <- factor(seu.obj$majorID, levels = c("CD4 T cell", "CD8/NK cell","gd T cell", "DN T cell",
                                                      "Cycling T cell", "Monocyte", "DC", "Neutrophil",
                                                     "B cell","Plasma cell","Granulocyte", "CD34+ Unclassified"))

levels(seu.obj$majorID)[12] <- "CD34+ Unk"
#plot cell frequency plots
fig2c <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source",
               namez = "name", 
               colz = "colz"
              ) + theme(axis.text.x = element_blank())
ggsave("./output/fig2c_freqPlots_majorID_hVoWadj_CLEAN_introns.png", width = 7, height = 6) 

supp <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "orig.ident", clusters = "clusterID") +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow =T))

ggsave(file = './output/stackedBar_hvo.png', width = 8, height = 12)

#####
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "clusterID") +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(legend.position = "top") + guides(fill = guide_legend(nrow = 1, byrow =T))

ggsave(file = './output/supp_stackedBar_h7.png', width = 10, height = 10)

#stash cell identies by majorID
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "prettyName")

#complete linDEG in pseudobulk-type format
linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = T, groupBy = "All cells", comparision = "cellSource", outDir = "./output/", outName = "fig2d_all", colUp = "red", colDwn = "blue",subtitle = F)

#complete linDEG in each major cell type
linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = T, groupBy = "prettyName", comparision = "cellSource", outDir = "./output/", outName = "fig2d_all", colUp = "red", colDwn = "blue",subtitle = F)

#remove OS_1,2,4,5,and 10 from analysis
seu.obj.ageSub <- subset(seu.obj,
                  subset = 
                  name !=  "OS_1" & name !=  "OS_2" & name !=  "OS_4" & name !=  "OS_5" & name !=  "OS_10"
                 ) 

p <- freqPlots(seu.obj = seu.obj.ageSub, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source", comp = "cellSource",
               namez = levels(seu.obj.ageSub$name)[c(1:7,10,13,14,15,16)], 
               colz = levels(seu.obj.ageSub$colz)[c(1:7,10,13,14,15,16)]
                                                  )
ggsave("./output/freqPlots_majorID_hVoWadj_CLEAN_introns_ageMatched.png", width = 12, height = 6) 

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "age")

#reorder majorID factor levels
seu.obj$age <- factor(seu.obj$age, levels = c("Healthy", "Age-matched", "Old"))


p <- freqPlots(seu.obj = seu.obj, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source", comp = "age",
               namez = "name", 
               colz = "colz"
                                                  ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("./output/freqPlots_majorID_hVoWadj_CLEAN_introns_ageMatched_2.png", width = 12, height = 6) 

seu.obj.ageSub <- subset(seu.obj,
                  subset = 
                  age !=  "Healthy"
                 ) 

p <- freqPlots(seu.obj = seu.obj.ageSub, method = 1, nrow= 3, groupBy = "majorID", legTitle = "Cell source", comp = "age",
               namez = levels(seu.obj.ageSub$name)[8:17], 
               colz = levels(seu.obj.ageSub$colz)[8:17]
                                                  )
ggsave("./output/freqPlots_majorID_hVoWadj_CLEAN_introns_ageMatched_3.png", width = 12, height = 6) 


##########################################
##              cytotoxic               ##
##########################################
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
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank())

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

seu.obj <- readRDS(file = "./output/s3/cytotoxic_noGD_hVoWadj_CLEAN_introns_res0.8_dims35_dist0.3_neigh30_S3.rds")

#transfer TRDC data from ROS data
seu.obj.ros.cyto <- readRDS(file = "./ros_output/HvO_mt_r7_res1.5_dims50_S3.rds")

seu.obj <- AddMetaData(seu.obj, metadata = seu.obj.ros.cyto@assays$RNA@data[rownames(seu.obj.ros.cyto@assays$RNA@data) == "TRDC",], col.name = "TRDC")
seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 

features = c("TRDC")
p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 1, ncol = 1, features = features, color = "black", order = F) 
ggsave(plot = p, "./output/supp_trdc_cyto_featPlots_introns.png", width=3, height=3)


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



#plot inital cluster umap
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_sub",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8)
ggsave("./output/umap_cyto_raw_clusterID_sub_hVoWadj_CLEAN_introns_introns.png")

#plot cell frequency plots
pi <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "clusterID_sub", legTitle = "Cell source",
               namez = c("H_1", "H_2","H_3","H_4", "H_5","H_6", "H_7", 
                         "OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6", "OS_7", "OS_8", "OS_9", "OS_10"), 
               colz = c(sequential_hcl(palette = "Greens", n = 10)[1:7], 
                        sequential_hcl(palette = "Purp", n = 12)[1:10])
              ) + scale_x_discrete(labels=c('Healthy', 'OS'))

pi <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "clusterID_sub", legTitle = "Cell source",
               namez = c("H_1", "H_2","H_3","H_4", "H_5","H_6", "H_7", 
                         "OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6", "OS_7", "OS_8", "OS_9", "OS_10"), 
               colz = c(sequential_hcl(palette = "Greens", n = 10)[1:7], 
                        sequential_hcl(palette = "Purp", n = 12)[1:10])
              ) + scale_x_discrete(labels=c('Healthy', 'OS'))

ggsave("./output/freqPlots_cyto_clusterID_sub_hVoWadj_CLEAN_introns.png", width = 7, height = 6)


#generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "cyto_hVoWadj_CLEAN_introns",
                     outDir = "./output/cyto/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

#load in color data for plotting beauty
colArray <- read.csv("./cyto_hVoWadj.csv", header = T)
colArray$clusterID_sub <- as.factor(sort(as.numeric(as.character(colArray$clusterID_sub))))

namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub

p <- dotPlotBY_TYPE(seu_obj = seu.obj, pwdTOvilnCSVoutput = "./output/cyto/cyto_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave("./output/cyto/cyto_bigDots.png", width = 20, height = 12)

singleR(seu.obj = seu.obj, outName = "cyto", clusters = "clusterID_sub", outDir = "./output/")

#incredibly ineffcient way to do this... cionvert to list of list and run in 1 cmd.... duh
ecLists <- read.csv("ecScores.csv", header = F)
naive_list <- ecLists[ecLists$V1 == "Naïve",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(naive_list),
                  name="NAIVE_SIG")

eff_list <- ecLists[ecLists$V1 == "Effector",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(eff_list),
                  name="EFFECTOR_SIG")

cd8.tem_list <- ecLists[ecLists$V1 == "CD8_TEM",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(cd8.tem_list),
                  name="CD8_TEM_SIG")

cd8.tcm_list <- ecLists[ecLists$V1 == "CD8_TCM",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(cd8.tcm_list),
                  name="CD8_TCM_SIG")

mait_list <- ecLists[ecLists$V1 == "MAIT",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(mait_list),
                  name="MAIT_SIG")

cd8.reg_list <- ecLists[ecLists$V1 == "CD8_REG",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(cd8.reg_list),
                  name="CD8_REG_SIG")

cd8.reg.mem_list <- ecLists[ecLists$V1 == "CD8_REG",2]

seu.obj <- AddModuleScore(seu.obj,
                  features = list(cd8.reg.mem_list),
                  name="CD8_REG_SIG")


features <- c("NAIVE_SIG1","EFFECTOR_SIG1", "CD8_TEM_SIG1",
              "CD8_TCM_SIG1", "MAIT_SIG1","CD8_REG_SIG1", "CD8_REG_SIG", "NCR1", "NCR3"
              )
test <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 4, features = features, color = "black", order = F,bottomLeg = T,  min.cutoff = "q20", pt.size = 0.001)
ggsave("./output/eScores_helper_hVoWadj_CLEAN_introns.png", width=12, height=6)

pi <- VlnPlot(
    object = seu.obj,
    cols = colArray$colour,
    pt.size = 0,
    same.y.lims = T,
    combine = FALSE,
            features = c("CD3E", "CD4", "CD8A", "IL2RB",
                         "CCL4", "CCL5", "CCR5", "CCR7",
                         "NCR3", "KLRB1", "CD96","IL7R",
                         "NFKBID","GZMA", "GZMB", "GZMK",
                         "ITGAM","RARRES2","KLRF1", "LGALS3", 
                         "PI3","ID3", "GATA3", "IKZF2")
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
                  

ggsave("./output/fig_3b_vilnPlots_cyto_clusterID_sub_hVoWadj_CLEAN_introns.png", width = 22.45, height = 13.4)

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              cols = colArray$colour,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

fig3a <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol)

fig3a <- fig3a + inset_element(umapHighLight, left=0, bottom=0, right=0.25, top=0.25)
ggsave("./output/fig3a_umap_cyto_final_hVoWadj_CLEAN_introns.png")

colArray$majorID <- "CD8, DN T, & NK cell subtypes"

#create custom legend
leg <- cusLeg(legend = colArray, colz = 4, rowz = 3, clusLabel = "clusterID_sub", legLabel = "classification", colorz = "colour",
                   groupLabel = "majorID", groupBy = "majorID", sortBy = "clusterID_sub", labCol = "labCol", headerSize = 6,
                   cellHeader = T, bump = 0, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 2.5, compress_x = 0.70
             )
ggsave("./output/fig3a_legend_cyto_hVoWadj_CLEAN_introns.png", width = 9, height = 3) 

#supplenntal to show how clustering changed
colArrayParent <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
colArrayParent <- colArrayParent %>% mutate(high_col = ifelse(majorID %in% highlight, color,"grey"))

sankeyPlot(seu.obj = seu.obj, new.ident = "clusterID_sub", old.ident = "clusterID", old.colorz = NULL,
                       new.colorz = NULL, old.labCol = NULL, new.labCol = NULL
                    )

#stash cell identies by majorID
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")

#reorder majorID factor levels
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("Naive CD8", "Effector CD8","Memory CD8",
                                                      "NK/NKT cell", "DN T cell", "CD8 gd T cell"))

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "classification")

#volc for DN T cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "classification", idents.1 = "DN T cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/", title = "DN T cell vs Clusters 0-7,9-11", idents.1_NAME= "DN_T_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T
                    ) 

pi <- p[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave("./output/fig3c_dnTs_vs_otherCD8tCells_volc.png")


#volc for NK cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "classification", idents.1 = "NK cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/", title = "NK cell vs Clusters 1-8,10,11", idents.1_NAME= "NK_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T
                    ) 

pi <- p[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave("./output/fig3d_NKs_vs_otherCD8tCells_volc.png")

#volc for Cluster 11 cells
p <- btwnClusDEG(seu.obj = seu.obj, groupBy = "majorID_sub", idents.1 = "CD8 gd T cell", idents.2 = NULL, bioRep = "orig.ident",
                        minCells = 25, outDir = "./output/", title = "CD8+ gd T cell vs Clusters 1-10", idents.1_NAME= "CD8_gd_T_cell", idents.2_NAME = "Other_cyto_cells", returnVolc = T
                    ) 

pi <- p[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave("./output/fig3e_cd8gdT_vs_otherCD8tCells_volc.png")

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")

#plot cell frequency plots
fig3 <- freqPlots(seu.obj, method = 1, nrow= 2, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name",comp = "cellSource", no_legend = T,
               namez = "name", 
               colz = "colz"
              ) + scale_x_discrete(labels=c('Healthy', 'OS'))
ggsave("./output/freqPlots_cyto_majorID_sub_hVoWadj_CLEAN_introns.png", width = 5, height = 5) 

#test linDEGs
linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = T, groupBy = "clusterID_sub", comparision = "cellSource", outDir = "./output/cyto/", outName = "cyto", colUp = "red", colDwn = "blue",subtitle = F)



vilnSplitComp(seu.obj = seu.obj, groupBy = "clusterID_sub", refVal = "cellSource", outName = "cyto", nPlots = 9,
                          outDir = "./output/cyto/"
                       )


colArray$clusterID_sub <- as.factor(sort(as.numeric(as.character(colArray$clusterID_sub))))
p <- vilnSplitCompxGene(seu.obj = seu.obj, groupBy = "clusterID_sub", comp = "cellSource", metaAdd = "majorID", features = c("CCL5","FCER1G", "IL2RB"), 
                               cols = c("mediumseagreen","mediumpurple1"), save = FALSE, outName = "", outDir = "", 
                               height = 4, width = 8, labelData = colArray
                              )
ggsave("./output/viln_cyto_genez_hVoWadj_CLEAN_introns.png", width = 7, height = 4) 

namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub

#ensembl connection is intemetendtly bad, try later
p <- dotPlotBY_TYPE(seu_obj = seu.obj, pwdTOvilnCSVoutput = "./output/cyto/cyto_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub", namedCols = namedCols, database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30"
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave("bigDots.png", width = 20, height = 12)
##########################################
##               helpers                ##
##########################################
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


seu.obj <- readRDS(file = "./output/s3/helper_hVoWadj_CLEAN_introns_res0.9_dims40_S3.rds")


seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "helper_hVoWadj_CLEAN_introns", 
           final.dims = 40, final.res = 0.9, returnFeats = F, saveRDS = T, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.2,
            n.neighbors = 20, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )

### load in processed data ###
seu.obj <- readRDS(file = "./output/s3/helper_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.2_neigh20_S3.rds")

#load in meta data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "clusterID_sub_clean")
seu.obj$clusterID_sub_clean <- factor(x = Idents(seu.obj), levels = sort(as.numeric(as.character(levels(seu.obj)))))
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")

#create subset object to remove por qulity cells from downstream analysis
seu.obj.sub <- subset(seu.obj, subset = majorID_sub != "exclude")

#load in color data for pretty plots
colArray <- read.csv("./helper_hVoWadj_forLeg_020623.csv", header = T)
colArray <- colArray[-14,] # remove the first poor quality pop from list

namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub_clean
singleR(seu.obj = seu.obj.sub, outName = "cd4", clusters = "clusterID_sub", outDir = "./output/")

#generate violing plots of defining features
vilnPlots(seu.obj = seu.obj.sub, groupBy = "clusterID_sub_clean", numOfFeats = 24, outName = "cd4_hVoWadj_CLEAN_introns",
                     outDir = "./output/helper/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

#create DotPlot by feature type
p <- dotPlotBY_TYPE(seu_obj = seu.obj.sub, pwdTOvilnCSVoutput = "./output/helper/cd4_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub_clean",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave("bigDots.png", width = 20, height = 12)

#Pathway analysis on Cluster 14
clus.markers <- read.csv(file = "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/helper/cd4_hVoWadj_CLEAN_introns_gene_list.csv")
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

ggsave("./output/supp_clus14_gsea.png", width = 8)

ecLists <- read.csv("ecScores.csv", header = F)
naive_list <- ecLists[ecLists$V1 == "Naïve",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(naive_list),
                  name="NAIVE_SIG")

th1_list <- ecLists[ecLists$V1 == "TH1",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(th1_list),
                  name="TH1_SIG")

th2_list <- ecLists[ecLists$V1 == "TH2",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(th2_list),
                  name="TH2_SIG")

th17_list <- ecLists[ecLists$V1 == "TH17",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(th17_list),
                  name="TH17_SIG")

cd4.ctl_list <- ecLists[ecLists$V1 == "CD4_CTL",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(cd4.ctl_list),
                  name="CD4_CTL_SIG")

treg_list <- ecLists[ecLists$V1 == "TREG",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(treg_list),
                  name="TREG_SIG")

tfh_list <- ecLists[ecLists$V1 == "TFH",2]

seu.obj.sub <- AddModuleScore(seu.obj.sub,
                  features = list(tfh_list),
                  name="TFH_SIG")


features <- c("TH1_SIG1","TH2_SIG1", "TH17_SIG1",
              "CD4_CTL_SIG1", "TREG_SIG1","TFH_SIG1"
              )
test <- prettyFeats(seu.obj = seu.obj.sub, nrow = 2, ncol = 3, features = features, color = "black", order = F,bottomLeg = T,  min.cutoff = "q20", pt.size = 0.001)
ggsave("./output/eScores_helper_hVoWadj_CLEAN_introns.png", width=9, height=6)




# Plot scores
FeaturePlot(seu.obj.sub,
            features = "TH17_SIG1", label = TRUE, repel = TRUE, min.cutoff = "q20", pt.size = 0.1)+
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  

ggsave("./output/test.png")

#plot fig4a
fig4a <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub_clean",
              cols = colArray$colour,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

fig4a <- cusLabels(plot = fig4a, shape = 21, size = 8, alpha = 0.8, labCol = colArray$labCol)
fig4a <- fig4a + inset_element(umapHighLight, left=0, bottom=0, right=0.25, top=0.25)
ggsave("./output/fig4a_umap_helper_final_hVoWadj_CLEAN_introns.png")

leg <- cusLeg(legend = colArray, colz = 4, rowz = 4, clusLabel = "clusterID_sub_clean", legLabel = "classification", colorz = "colour",
                   groupLabel = "title", groupBy = "title", sortBy = "clusterID_sub", labCol = "labCol", headerSize = 6,
                   cellHeader = T, bump = 0, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 2, compress_x = 0.8, spaceBtwnCols = c(0.4,0.45,0.3)
             )
ggsave("./output/fig4a_legend_helper_hVoWadj_CLEAN_introns.png", width = 12, height = 3) 


#supplenntal to show how clustering changed
colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
saveCells <- c("CD4 T cell")
colArray_all <- colArray_all %>% filter(majorID %in% saveCells)

p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_sub", old.ident = "clusterID", old.colorz = colArray_all$color,
                       new.colorz = colArray$colour, old.labCol = colArray_all$labCol, new.labCol = colArray$labCol, flowCol = "grey"
                    )
ggsave("./output/supp5a_sankey_helper_hVoWadj_CLEAN_introns.png", width=7, height=10)


features <- c("CD3G","CD4", "CD8A","IL7R",
              "SELL","CCR7", "CD40LG", "CD28",
              "IL18R1","CCR5","CTLA4","IKZF2",
              "S100A5", "GATA3", "SYTL3", "RUNX2"
              )
fig4b <- prettyFeats(seu.obj = seu.obj, nrow = 4, ncol = 4, features = features, color = "black", order = F,bottomLeg = T) 
ggsave("./output/fig4b_featPlots_helper_hVoWadj_CLEAN_introns.png", width=12, height=12)


#reorder majorID factor levels
seu.obj.sub$majorID_sub <- factor(seu.obj.sub$majorID_sub, levels = rev(c("Naive", "TCM","TEM", "TEM, Th1-like", 
                                                              "TEM, Th2-like","TEM, Th17-like","T reg",
                                                              "IFN signature")))
seu.obj.sub <- ScaleData(seu.obj.sub)

fig4c <- majorDot(seu.obj = seu.obj.sub, groupBy = "majorID_sub",
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

ggsave("./output/fig4c_dots_helper_hVoWadj_CLEAN_introns.png",width = 8,height=4)

fig_supp <- majorDot(seu.obj = seu.obj.sub, groupBy = "majorID_sub",
                     features = c("NAIVE_SIG1", "TH1_SIG1","TH2_SIG1", "TH17_SIG1",
                                  "TREG_SIG1","TFH_SIG1"
                                 )
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 24, "pt")) + scale_y_discrete(position = "right")
ggsave("./output/sup_helper_hVoWadj_CLEAN_introns.png",width = 8,height=4)



#then run dot plot by feat type
stackedBar(seu.obj = seu.obj.sub, downSampleBy = "cellSource", groupBy = "orig.ident", clusters = "clusterID_sub_clean") +
scale_fill_manual(labels = c("H_1", "H_2","H_3","H_4", "H_5","H_6", "H_7", 
                         "OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6", "OS_7", "OS_8", "OS_9", "OS_10"), 
               values = c(sequential_hcl(palette = "Greens", n = 10)[1:7], 
                        sequential_hcl(palette = "Purp", n = 12)[1:10])) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow =T))

ggsave(file = './output/stackedBar_helpers.png', width = 10, height = 7.5)

#investigated gene expression changes in os
vilnSplitComp(seu.obj = seu.obj.sub, groupBy = "clusterID_sub_clean", refVal = "cellSource", 
              outName = "cd4", nPlots = 20, outDir = "./output/helper/"
                       )

#plot cell frequency plots
fig3 <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "majorID_sub", legTitle = "Cell source",
               namez = c("H_1", "H_2","H_3","H_4", "H_5","H_6", "H_7", 
                         "OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6", "OS_7", "OS_8", "OS_9", "OS_10"), 
               colz = c(sequential_hcl(palette = "Greens", n = 10)[1:7], 
                        sequential_hcl(palette = "Purp", n = 12)[1:10])
              ) + theme(legend.position="none") 
ggsave("./output/freqPlots_helper_majorID_sub_hVoWadj_CLEAN_introns.png", width = 7, height = 5) 

#set up for pathway analysis
can_gene_sets <- as.data.frame(msigdbr(species = "dog", category = "C7"))
msigdbr_list <- split(x = can_gene_sets$gene_symbol, f = can_gene_sets$gs_name)
datas <- can_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

#narrow the search....
datas <- datas[grepl("CD4|TREG|TH1|TH17|TH2",datas$gs_name),]

#input ref data
datas_2 <- read.csv("./ecScores.csv")
colnames(datas_2) <- c("gs_name", "gene_symbol")

#volc for Th1s
fig4d1 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", idents.1 = "TEM, Th1-like", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/", title = "TEM,Th1-like vs other CD4 clusters", idents.1_NAME= "TEM__Th1_like", idents.2_NAME = "Other_T_cells", returnVolc = T, doLinDEG = F
                    ) 

fig4d1 <- fig4d1[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0) 
ggsave("./output/fig4d1_Th1_vs_otherCD4tCells_volc.png")


#use features up Th1 for pathway analysis
th1.df <- read.csv("./output/Th1_vs_Other_T_cells_all_genes.csv")
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

ggsave("./output/supp_th1_path.png", width = 16)

#use short gene lists
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.05,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Th subtype comparision")

ggsave("./output/supp_th1_vs_ref.png", width = 8)

#volc for Th2s
fig4d2 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", idents.1 = "TEM, Th2-like", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/", title = "TEM, Th2-like vs other CD4 clusters", idents.1_NAME= "Th2", idents.2_NAME = "Other_T_cells", returnVolc = T,doLinDEG = F
                    ) 

fig4d2 <- fig4d2[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave("./output/fig4d2_Th2_vs_otherCD4tCells_volc.png")

th2.df <- read.csv("./output/Th2_vs_Other_T_cells_all_genes.csv")
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

ggsave("./output/supp_th2_path.png", width = 16)

#use short gene lists
enriched <- as.data.frame(enricher(gene = th2.df$gene, TERM2GENE = datas_2, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.05,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Th subtype comparision")

ggsave("./output/supp_th2_vs_ref.png", width = 8)

#volc for Th17s
fig4d3 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", idents.1 = "TEM, Th17-like", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/", title = "TEM, Th17-like vs other CD4 clusters", idents.1_NAME= "Th17", idents.2_NAME = "Other_T_cells", returnVolc = T,doLinDEG = F
                    ) 

fig4d3 <- fig4d3[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave("./output/fig4d3_Th17_vs_otherCD4tCells_volc.png")

th17.df <- read.csv("./output/Th17_vs_Other_T_cells_all_genes.csv")
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

ggsave("./output/supp_th17_path.png", width = 16)

#use short gene lists
enriched <- as.data.frame(enricher(gene = th17.df$gene, TERM2GENE = datas_2, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.05,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Th subtype comparision")

ggsave("./output/supp_th17_vs_ref.png", width = 8)


#volc for Tregs
fig4d4 <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", idents.1 = "T reg", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.58, 
                        minCells = 25, outDir = "./output/", title = "T reg vs other CD4 clusters", idents.1_NAME= "T reg", idents.2_NAME = "Other_T_cells", returnVolc = T,doLinDEG = F
                    ) 

fig4d4 <- fig4d4[[1]] + theme(legend.position = "none", 
                     axis.title.y = element_blank()) + scale_x_symmetric(mid = 0)
ggsave("./output/fig4d4_Tregs_vs_otherCD4tCells_volc.png")

tregs.df <- read.csv("./output/T_reg_vs_Other_T_cells_all_genes.csv")
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

ggsave("./output/supp_tregs_path.png", width = 16)

#use short gene lists
enriched <- as.data.frame(enricher(gene = tregs.df$gene, TERM2GENE = datas_2, pvalueCutoff = 1))
enriched$GeneRatio <- sapply(enriched$GeneRatio, function(x){eval(parse(text=x))})
enriched <- enriched[enriched$p.adjust < 0.05,]
p <- ggplot(enriched, aes(x = GeneRatio, y = fct_reorder(ID,GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = p.adjust)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        ylab(NULL) +
        ggtitle("Th subtype comparision")

ggsave("./output/supp_tregs_vs_ref.png", width = 8)


#run slingshot
sce.obj <- as.SingleCellExperiment(seu.obj.sub)
rd1 <- Embeddings(seu.obj.sub, reduction = "pca")[,1:2]
rd2 <- Embeddings(seu.obj.sub, reduction = "umap")
colnames(rd2) <- c('UMAP1', 'UMAP2')

#assign origin
start.clus <- '0'

reducedDims(sce.obj) <- SimpleList(PCA = rd1, UMAP = rd2)

sce.obj <- slingshot(sce.obj, clusterLabels = 'clusterID_sub_clean', reducedDim = 'UMAP', start.clus = start.clus) #need to get it so i can specify tart pooint here - its defulting to a weird sopt

#identify lineages
lin1 <- getLineages(Embeddings(seu.obj.sub, reduction = "umap"), sce.obj$clusterID_sub_clean, start.clus = start.clus)
branchData <- SlingshotDataSet(lin1)@lineages

#plot the lineages
colArray <- read.csv("helper_hVoWadj_forLeg.csv")
plot <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "clusterID_sub_clean",
              cols = colArray$colour[-c(14)],
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

p <- cleanSling(plot = plot, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = branchData)
ggsave("./output/supp_branch_cd4.png")

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

# features <- "SYTL3"
# colorBy <- "clusterID_sub"
# plotBy <- "slingPseudotime_2"
# rmZeros <- F

# counts <- as.data.frame(seu.obj.sub@assays$RNA@data[rownames(seu.obj.sub@assays$RNA@data) %in% features,])
# meta <- as.data.frame(seu.obj.sub@meta.data)
# meta <- meta[,c(colorBy,plotBy)]

# plotData <- cbind(counts,meta)
# colnames(plotData) <- c("countz","colorByz","plotByz")

# if(rmZeros){
#     plotData <- plotData[plotData$countz !=0,]
# }

# lineData <- plotData %>% na.omit() %>% arrange(plotByz) %>% mutate(bin = ntile(plotByz, nrow(plotData)/50)) %>% group_by(bin) %>% mutate(avg = mean(countz)) %>% arrange(plotByz) %>% as.data.frame()

# p <- ggplot()  + geom_point(data = plotData, aes(x = plotByz, y = countz, colour = colorByz), alpha = 0.5) + geom_smooth(data = lineData, aes(x = plotByz, y = avg)) + labs(title = features, x = "Pseudotime", y = "Normalized count") + guides(colour = guide_legend(ncol = 3))
# ggsave("test.png", width = 12, height = 2)

features <- c("slingPseudotime_5","slingPseudotime_2","slingPseudotime_3","slingPseudotime_1")
titles <- c("Lineage 1","Lineage 2","Lineage 3","Lineage 4")

p <- prettyFeats(seu.obj = seu.obj.sub, nrow = 2, ncol = 2, features = features, color = "black", order = F, title = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")
ggsave(plot = p, "./output/supp_pseudoTime_helper.png", width=6, height=6)

##########################################
##               myeloid                ##
##########################################
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
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank())

seu.obj <- readRDS(file = "./output/s3/myeloid_hVoWadj_CLEAN_introns_res0.8_dims40_S3.rds")


seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "myeloid_hVoWadj_CLEAN_introns", 
           final.dims = 40, final.res = 0.9, returnFeats = F, saveRDS = F, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.35,
            n.neighbors = 20, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )

#load in processed data
seu.obj <- readRDS(file = "./output/s3/myeloid_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.35_neigh20_S3.rds")

#load in key meta data and remove poor quality cells for downstream analysis
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub2")

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "clusterID_sub_clean")
seu.obj$clusterID_sub_clean <- factor(x = Idents(seu.obj), levels = sort(as.numeric(as.character(levels(seu.obj)))))

seu.obj.sub <- subset(seu.obj, subset = majorID_sub != "exclude")


#read in color meta data for ploting uniformity
colArray <- read.csv("./coldaf_myeloid_new.csv")

#exploratory anakysis to define each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, 
          outName = "myeloid_hVoWadj_CLEAN_introns",
          outDir = "./output/myeloid/", outputGeneList = T, 
          filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )
namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub_clean

p <- dotPlotBY_TYPE(seu_obj = seu.obj.sub, pwdTOvilnCSVoutput = "./output/myeloid/myeloid_hVoWadj_CLEAN_introns_gene_list.csv", groupBy = "clusterID_sub_clean",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave("bigDots.png", width = 20, height = 12)

singleR(seu.obj = seu.obj.sub, outName = "myeloid", clusters = "clusterID_sub", outDir = "./output/")

#plot figure 5a
fig5a <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub_clean",
              cols = c(colArray$colour[-c(16,19)],"grey"),
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )
fig5a <- cusLabels(plot = fig5a, shape = 21, size = 8, alpha = 0.8, labCol = c(colArray$labCol[-c(16,19)],NA))
fig5a <- fig5a + inset_element(umapHighLight, left=0, bottom=0.75, right=0.25, top=1)
ggsave("./output/fig5a_umap_myeloid_final_hVoWadj_CLEAN_introns.png")

#load in data for legend
colArray_forLeg <- read.csv("./coldaf_myeloid_new_forLeg.csv", header = T)
colArray_forLeg$clusterID_sub <- as.factor(sort(as.numeric(as.character(colArray_forLeg$clusterID_sub))))

#create custom legend - with manual hack to force legend into position
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

ggsave("./output/fig5a_legend_myeloid_hVoWadj_CLEAN_introns.png", width = 14, height = 10) 

#supplenntal to show how clustering changed
colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
saveCells <- c("Monocyte","Neutrophil","Granulocyte","DC")
colArray_all <- colArray_all %>% filter(majorID %in% saveCells)
                             
p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_sub", old.ident = "clusterID", 
                old.colorz = colArray_all$color, new.colorz = colArray$colour, 
                old.labCol = colArray_all$labCol, new.labCol = colArray$labCol, 
                flowCol = "#D87FAB"
                    )
ggsave("./output/supp_sankey_myeloid_hVoWadj_CLEAN_introns.png", width=10, height=14)

#create feature plots for fig5b
features <- c("DLA-DRA","ANPEP", "FLT3",
              "S100A12","CD3G", "CD4",
              "LTF","IL5RA", "BATF3",
              "IRF8","CADM1", "TCF4")

#some additional features used to try to ID DCs
features <- c("FLT3","IRF8", "IRF4",
              "BATF3","FCER1A", "KLF4",
              "RUNX3","AXL", "CADM1",
              "BTLA","KIT", "CD1C",
              "IL3RA","CDH1", "SDC2")

fig5b <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = "black", order = F) 
#ggsave(plot = fig5b, "./output/fig5b_featPlots_myeloid_hVoWadj_CLEAN_introns.png", width=13, height=9)
ggsave(plot = fig5b, "./output/supp_dc_featPlots_myeloid_hVoWadj_CLEAN_introns.png", width=10, height=15)

#order myeloid cell populations for frequency plots and plot fig5c
seu.obj.sub$majorID_sub2 <- factor(seu.obj.sub$majorID_sub2, levels = c("DC","Monocyte, CD4-", "Monocyte, IFN signature",
                                                                      "Monocyte, CD4+", "M-MDSC", "PMN-MDSC",
                                                                      "Neutrophil","Eosinophil","Basophil"))

fig5c <- freqPlots(seu.obj = seu.obj.sub, method = 1, nrow = 3, 
                   groupBy = "majorID_sub2", legTitle = "Cell source",
                   refVal = "name",comp = "cellSource", no_legend = T,
                   namez = "name", colz = "colz"
                  ) + scale_x_discrete(labels=c('Healthy', 'OS'))#+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave("./output/fig5c_freqPlots_meyloid_majorID_sub_hVoWadj_CLEAN_introns.png", width = 8, height = 6)

#complete linDEG on cluster 10 (IFN responsive)
geneList <- linDEG(seu.obj = seu.obj, threshold = 0.5, thresLine = T, groupBy = "majorID_sub", comparision = "cellSource", outDir = "./output/myeloid/", outName = "myeloid", colUp = "red", colDwn = "blue",subtitle = T, cluster = "Monocyte, IFN responsive", returnUpList = F, returnDwnList = T)

#use features up in healthys of cluster 10 for pathway analysis
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

ggsave("./output/supp_clus10Upinealth.png", width = 8)


features <- c("DLA-DRA", "FCGR1A", "CD86", "CCR2")


fig5b <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 4, features = features, color = "black", order = F) 
ggsave(plot = fig5b, "./output/supp_featPlots_myeloid_hVoWadj_CLEAN_introns.png", width=13, height=3)

#compare to ros data
ros_myeloid <- readRDS("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/ros_output/sub_HvO_Myeloid_mt_r7_res1.1_dims40_S3.rds")
features <- c("DLA-DRA", "FCGR1A", "CD86", "CCR2", "S100A12", 
              "CD14", "CD4", "ANPEP", "CD68")

fig5b <- prettyFeats(seu.obj = ros_myeloid, nrow = 3, ncol = 3, features = features, color = "black", order = F) 
ggsave(plot = fig5b, "./output/supp_featPlots_myeloid_hVoWadj_CLEAN_introns.png", width=9, height=9)


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
ggsave("./output/umap_myeloid_CD14_highlight_trdcPos.png")

#complete linDEG for each major myeloid cell group -- makes fig5E and fig5G
linDEG(seu.obj = seu.obj.sub, threshold = 0.5, thresLine = T, groupBy = "majorID_sub", 
       comparision = "cellSource", outDir = "./output/", outName = "majorID_sub_degs", 
       cluster = NULL, labCutoff = 20, colUp = "red", colDwn = "blue", 
       subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F
                  )

#addtional statsitical valudation
vilnSplitComp(seu.obj = seu.obj.sub, groupBy = "majorID_sub", refVal = "cellSource", 
              outName = "myeloid", nPlots = 20, outDir = "./output/myeloid/"
                       )
                             

#volc for PMN-MDSC - fig5d
fig5d <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", 
                      idents.1 = "PMN-MDSC", idents.2 = "Neutrophil", bioRep = "orig.ident",
                      padj_cutoff = 0.01, lfcCut = 0.58, minCells = 25, 
                      outDir = "./output/", title = "PMN-MDSC vs Clusters 1, 7, 9 & 20", 
                      idents.1_NAME= "PMN-MDSC", idents.2_NAME = "Neutrophil", 
                      returnVolc = T,doLinDEG = F
                    ) 
fig5d <- fig5d[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave("./output/fig5d1_fig4d_PMN-MDSC_vs_otherCD4tCells_volc.png")                    
                             
#volc for M-MDSC - fig5f
fig5f <- btwnClusDEG(seu.obj = seu.obj.sub, groupBy = "majorID_sub", 
                     idents.1 = "M-MDSC", idents.2 = c("Monocyte, CD4-","Monocyte, CD4+","Monocyte, IFN signature"), 
                     bioRep = "orig.ident", padj_cutoff = 0.01, lfcCut = 0.58, minCells = 25, 
                     outDir = "./output/", title = "M-MDSC vs Clusters 0, 2, 3, 4, 6, 10 & 11", 
                     idents.1_NAME= "M-MDSC", idents.2_NAME = "Monocyte", 
                     returnVolc = T, doLinDEG = F, paired = T, addLabs = c("S100A5", "BPT", "TFPT", "RACK1", "DLA-DRA", "DLA-DQA1", "HLA-DQB2")
                    ) 
fig5f <- fig5f[[1]] + theme(legend.position = "none") + scale_x_symmetric(mid = 0)
ggsave("./output/fig5f_M-MDSCvMonos_volc.png")                    

highlight <- c("M-MDSC")
colArray <- colArray %>% mutate(high_col = ifelse(majorID_sub %in% highlight, colour,"grey"))

umapHighLight <- pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_sub",
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(umapHighLight) + theme(axis.title = element_blank())
ggsave(file = "./output/mmdscHighlight.png")

#convert seurat object to sce for Slingshot
seu.obj.sub@meta.data <- seu.obj.sub@meta.data[,1:63]
sce.obj <- as.SingleCellExperiment(seu.obj.sub)
rd1 <- Embeddings(seu.obj.sub, reduction = "pca")[,1:2]
rd2 <- Embeddings(seu.obj.sub, reduction = "umap")
colnames(rd2) <- c('UMAP1', 'UMAP2')

reducedDims(sce.obj) <- SimpleList(PCA = rd1, UMAP = rd2)

#assign origin
start.clus <- '11'

#OR run in unsupervised manner
#start.clus <- NULL 

sce.obj <- slingshot(sce.obj, clusterLabels = 'clusterID_sub', reducedDim = 'UMAP', start.clus = start.clus)

#identify lineages
lin1 <- getLineages(Embeddings(seu.obj.sub, reduction = "umap"), sce.obj$clusterID_sub, start.clus = start.clus)
branchData <- SlingshotDataSet(lin1)@lineages

#plot the lineages
# colArray <- colArray[-c(16,19),]
plot <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              cols = colArray$colour,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

p <- cleanSling(plot = plot, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = branchData)
ggsave("./output/test.png")

metaData <- as.data.frame(sce.obj@colData@rownames)
colnames(metaData) <- "barcode"
metaData$slingPseudotime_1 <- sce.obj$slingPseudotime_1
metaData$slingPseudotime_2 <- sce.obj$slingPseudotime_2
metaData$slingPseudotime_3 <- sce.obj$slingPseudotime_3
metaData$slingPseudotime_4 <- sce.obj$slingPseudotime_4
metaData$slingPseudotime_5 <- sce.obj$slingPseudotime_5
metaData$slingPseudotime_6 <- sce.obj$slingPseudotime_6
metaData$slingPseudotime_7 <- sce.obj$slingPseudotime_7
metaData$slingPseudotime_8 <- sce.obj$slingPseudotime_8
metaData$slingPseudotime_9 <- sce.obj$slingPseudotime_9
metaData$slingPseudotime_10 <- sce.obj$slingPseudotime_10
metaData$slingPseudotime_11 <- sce.obj$slingPseudotime_11

metaData <- metaData %>% rowwise %>% mutate(pseudoTime = mean(c(slingPseudotime_1,slingPseudotime_2,slingPseudotime_3,slingPseudotime_4,slingPseudotime_5,slingPseudotime_6,slingPseudotime_7,slingPseudotime_8,slingPseudotime_9,slingPseudotime_10,slingPseudotime_11),na.rm=TRUE)) 

seuMeta <- seu.obj.sub@meta.data[,1:63] %>% mutate(barcode = rownames(.))

newMeta <- seuMeta %>% left_join(metaData, by = 'barcode')
rownames(newMeta) <- newMeta$barcode
seu.obj.sub@meta.data <- newMeta

features <- c("slingPseudotime_1","slingPseudotime_2","slingPseudotime_3","slingPseudotime_4","slingPseudotime_5",
              "slingPseudotime_6","slingPseudotime_7","slingPseudotime_8","slingPseudotime_9","slingPseudotime_10",
              "slingPseudotime_11", "pseudoTime")
titles <- c("Lineage 1","Lineage 2","Lineage 3","Lineage 4","Lineage 5",
           "Lineage 6","Lineage 7","Lineage 8","Lineage 9","Lineage 10", "Lineage 11", "PseudoTime")

p <- prettyFeats(seu.obj = seu.obj.sub, nrow = 4, ncol = 4, features = features, color = "black", order = F, title = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")
ggsave(plot = p, "./output/pseudoTime_myeloid.png", width=12, height=9)


features <- "IL1B"
colorBy <- "clusterID_sub"
plotBy <- "slingPseudotime_5"
rmZeros <- F

counts <- as.data.frame(seu.obj.sub@assays$RNA@data[rownames(seu.obj.sub@assays$RNA@data) %in% features,])
meta <- as.data.frame(seu.obj.sub@meta.data)
meta <- meta[,c(colorBy,plotBy)]

plotData <- cbind(counts,meta)
colnames(plotData) <- c("countz","colorByz","plotByz")

if(rmZeros){
    plotData <- plotData[plotData$countz !=0,]
}

lineData <- plotData %>% na.omit() %>% arrange(plotByz) %>% mutate(bin = ntile(plotByz, nrow(plotData)/100)) %>% group_by(bin) %>% mutate(avg = mean(countz)) %>% arrange(plotByz) %>% as.data.frame()

p <- ggplot()  + geom_point(data = plotData, aes(x = plotByz, y = countz, colour = colorByz), alpha = 0.5) + geom_smooth(data = lineData, aes(x = plotByz, y = avg)) + labs(title = features, x = "Pseudotime", y = "Normalized count") + guides(colour = guide_legend(ncol = 3)) + scale_colour_manual(values = colArray$colour)
ggsave("./output/test.png", width = 12, height = 2)

#focus on monocytes
seu.obj <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")

seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "Monocyte"
                 ) 

indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "monocyte", preSub = T,
         vars.to.regress = c("percent.mt", "percent.pal")
         )

seu.obj <- readRDS(file = "./output/s2/monocyte_S2.rds")

clusTree(seu.obj = seu.obj, outName = "monos", dout = "./output/clustree/", test_dims = c(45,40,35,30,25), algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj,
           outDir = "./output/s3/", outName = "monos", 
           final.dims = 25, final.res = 1.0, returnFeats = F, saveRDS = F, return_obj = T,
           stashID = "clusterID_sub", algorithm = 3, 
            prefix = "integrated_snn_res.", min.dist = 0.35,
            n.neighbors = 20, features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                           "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                           "CD4", "MS4A1", "PPBP","HBM")
           )


#extra exploratory analysis
vilnSplitComp(seu.obj = seu.obj, groupBy = "clusterID_sub", refVal = "cellSource", 
              outName = "myeloid", nPlots = 20, outDir = "./output/myeloid/"
                       )
                             
##########################################
##               b cell                 ##
##########################################

seu.obj <- readRDS(file = "./output/s3/bcell_hVoWadj_CLEAN_introns_res0.6_dims30_S3.rds")

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")

#load in colors
colArray <- read.csv("./bcell_hVoWadj.csv")

vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "bcell_hVoWadj_CLEAN_introns",
                     outDir = "./output/bcell/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

namedCols <- colArray$colour
names(namedCols) <- colArray$clusterID_sub

p <- dotPlotBY_TYPE(seu_obj = seu.obj, pwdTOvilnCSVoutput = "./output/bcell/bcell_hVoWadj_CLEAN_introns_res0.6_dims30_gene_list.csv", groupBy = "clusterID_sub",  database = "clfamiliaris_gene_ensembl", exlcude = "^ENSCAF", boxColor = "grey30", namedCols = namedCols
                          ) + theme(panel.background = element_rect(fill='white'))
ggsave("bigDots.png", width = 20, height = 12)


#plot inital cluster umap
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_sub",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8)
ggsave("./output/umap_bcell_new_raw_hVoWadj_CLEAN_introns_introns.png")

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
ggsave("./output/umap_bcell_new_raw_hVoWadj_CLEAN_introns_introns.png")


features <- c("MS4A1","VPREB3","IGF1R", "PAX5","IGHM",
              "CXCR4","IRF5", "CD80", "CD40",
              "NAPSA","RAG1","DNTT"
              )
fig4b <- prettyFeats(seu.obj = seu.obj, nrow = 4, ncol = 4, features = features, color = "black", order = F,bottomLeg = F) 
ggsave("./output/supp_featPlots_bcell_hVoWadj_CLEAN_introns.png", width=12, height=9)


figBcell <- freqPlots(seu.obj.sub, method = 1, nrow= 2, groupBy = "pctID", legTitle = "Cell source",refVal = "name",comp = "cellSource", no_legend = T,
               namez = "name", 
               colz = "colz"
              ) + theme(axis.text.x = element_blank()) 

ggsave("./output/freqPlots_bcell_majorID_sub_hVoWadj_CLEAN_introns.png", width = 6, height = 4)
                             
colArray$title <- "B cell subtypes"
leg <- cusLeg(legend = colArray, colz = 3, rowz = 3, clusLabel = "clusterID_sub", legLabel = "pctID", colorz = "colour",
                   groupLabel = "title", groupBy = "title", sortBy = "clusterID_sub", labCol = "labCol", headerSize = 6,
                   cellHeader = T, bump = 0, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 2, compress_x = 0.8, spaceBtwnCols =c(0.5,0.4)
             )
ggsave("./output/supp_legend_bcell_hVoWadj_CLEAN_introns.png", width = 9, height = 3) 

linDEG(seu.obj = seu.obj.sub, threshold = 1, thresLine = T, groupBy = "pctID", comparision = "cellSource", outDir = "./output/", outName = "bcell", cluster = NULL, labCutoff = 20,
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F
                  )

colArray_all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)
saveCells <- c("Bcell")
colArray_all <- colArray_all %>% filter(colorID %in% saveCells)
                             
                             
p <- sankeyPlot(seu_obj = seu.obj.sub, new.ident = "clusterID_sub", old.ident = "clusterID", old.colorz = colArray_all$color,
                       new.colorz = colArray$colour[-9], old.labCol = colArray_all$labCol, new.labCol = colArray$labCol[-9], flowCol = "grey"
                    )

ggsave("./output/supp_sankey_bcell_hVoWadj_CLEAN_introns.png", width=7, height=10)

                            
########### healthy dataset ############
seu.obj.healthy <- readRDS(file = "./output/s3/h7_CLEAN_introns_res1.6_dims45_S3.rds")

vilnPlots(seu.obj = seu.obj.healthy, groupBy = "clusterID", numOfFeats = 24, outName = "h7",
                      outDir = "./output/h7/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)
                     
                             
singleR(seu.obj = seu.obj.healthy, outName = "h7_singleR", clusters = "clusterID", outDir = "./output/")
                             
Idents(seu.obj.healthy) <- "clusterID"
high_cells <- WhichCells(seu.obj.healthy, ident = "38")
                             
pi <- DimPlot(seu.obj.healthy, 
              reduction = "umap", 
              group.by = "clusterID",
              cols.highlight = "hotpink",
              cells.highlight = high_cells,
              pt.size = 0.25,
              label = F,
              label.box = F
 )
supp <- formatUMAP(pi)
ggsave("./output/umap_showCD4negs_1.png")
                             
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              cols.highlight = "hotpink",
              cells.highlight = high_cells,
              pt.size = 0.25,
              label = F,
              label.box = F
 )
supp <- formatUMAP(pi)
ggsave("./output/umap_showCD4negs_2.png")
                             



####### Make final table ##########

seu.obj.all <- readRDS(file = "./output/hVoWadj_CLEAN_introns_res1.9_dims45_S3.rds")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./HvOvadj_majorID_wIntrons.csv", groupBy = "majorID", metaAdd = "color")
colorData.all <- read.csv("./HvOvadj_majorID_wIntrons.csv", header = T)

Idents(seu.obj.all) <- "majorID"

seu.obj.all_laggers <- subset(seu.obj.all, idents =  c("gd T cell", "Cycling T cell","CD34+ Unclassified"))
seu.obj.all_laggers$majorID <- droplevels(seu.obj.all_laggers$majorID)
seu.obj.all_laggers$color <- droplevels(seu.obj.all_laggers$color)
colorData.all <- colorData.all[colorData.all$majorID %in% c("gd T cell", "Cycling T cell","CD34+ Unclassified"),]

seu.obj.cyto <- readRDS(file = "./output/s3/cytotoxic_noGD_hVoWadj_CLEAN_introns_res0.8_dims35_dist0.3_neigh30_S3.rds")
seu.obj.cyto <- loadMeta(seu.obj = seu.obj.cyto, metaFile = "./cyto_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "classification")
seu.obj.cyto <- loadMeta(seu.obj = seu.obj.cyto, metaFile = "./cyto_hVoWadj.csv", groupBy = "classification", metaAdd = "colour")
colorData.cyto <- read.csv("./cyto_hVoWadj.csv", header = T)

seu.obj.helper <- readRDS(file = "./output/s3/helper_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.2_neigh20_S3.rds")
seu.obj.helper <- loadMeta(seu.obj = seu.obj.helper, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")
seu.obj.helper <- loadMeta(seu.obj = seu.obj.helper, metaFile = "./helper_hVoWadj_forLeg_020623.csv", groupBy = "majorID_sub", metaAdd = "colour")
colorData.helper <- read.csv("./helper_hVoWadj_forLeg_020623.csv", header = T)

seu.obj.myeloid <- readRDS(file = "./output/s3/myeloid_hVoWadj_CLEAN_introns_res0.9_dims40_dist0.35_neigh20_S3.rds")
seu.obj.myeloid <- loadMeta(seu.obj = seu.obj.myeloid, metaFile = "./coldaf_myeloid_new.csv", groupBy = "clusterID_sub", metaAdd = "majorID_sub")    
seu.obj.myeloid <- loadMeta(seu.obj = seu.obj.myeloid, metaFile = "./coldaf_myeloid_new.csv", groupBy = "majorID_sub", metaAdd = "colour")    
colorData.myeloid <- read.csv("./coldaf_myeloid_new.csv", header = T)

seu.obj.bcell <- readRDS(file = "./output/s3/bcell_hVoWadj_CLEAN_introns_res0.6_dims30_S3.rds")
seu.obj.bcell <- loadMeta(seu.obj = seu.obj.bcell, metaFile = "./bcell_hVoWadj.csv", groupBy = "clusterID_sub", metaAdd = "pctID")
seu.obj.bcell <- loadMeta(seu.obj = seu.obj.bcell, metaFile = "./bcell_hVoWadj.csv", groupBy = "pctID", metaAdd = "colour")
colorData.bcell <- read.csv("./bcell_hVoWadj.csv", header = T)

classifications <- c(seu.obj.helper$majorID_sub,
                     seu.obj.all_laggers$majorID,
                     seu.obj.myeloid$majorID_sub,
                     seu.obj.cyto$classification,
                     seu.obj.bcell$pctID)
        
levels(classifications) <- ifelse(levels(classifications) == "exclude", "Low quality", levels(classifications))
levels(classifications) <- ifelse(levels(classifications) == "Poor quality", "Low quality", levels(classifications))

seu.obj.all <- AddMetaData(seu.obj.all, classifications, col.name = "finalID")
        
colourz <- as.data.frame(c(levels(seu.obj.helper$colour), 
                           levels(seu.obj.all_laggers$color),
                           levels(seu.obj.myeloid$colour),
                           levels(seu.obj.cyto$colour),
                           levels(seu.obj.bcell$colour)
                          )
                        )
colnames(colourz) <- "colour"
        
colourz$classification <- c(levels(seu.obj.helper$majorID_sub), 
                            levels(seu.obj.all_laggers$majorID), 
                            levels(seu.obj.myeloid$majorID_sub), 
                            levels(seu.obj.cyto$classification), 
                            levels(seu.obj.bcell$pctID)
                           )
colourz$classification <- ifelse(colourz$classification == "exclude", "Low quality", colourz$classification)
colourz$classification <- ifelse(colourz$classification == "Poor quality", "Low quality", colourz$classification)

colourz$classificationFinal <- colourz$classification

colourz$majorGroup <- c(rep("CD4 T cell", length(levels(seu.obj.helper$majorID_sub))), 
                        rep("Miscellaneous", length(levels(seu.obj.all_laggers$majorID))), 
                        rep("Myeloid", length(levels(seu.obj.myeloid$majorID_sub))), 
                        rep("CD8/NK cell", length(levels(seu.obj.cyto$classification))), 
                        rep("B cell", length(levels(seu.obj.bcell$pctID)))
                       )

rm(list = c("seu.obj.helper", "seu.obj.myeloid", "seu.obj.all_laggers", "seu.obj.cyto", "seu.obj.bcell"))

#some annoying error where labCol cannont be brought over
colnames(colorData.all)[4] <- "colour"

labColData <- as.data.frame(mapply(c, colorData.helper[,c("majorID_sub","labCol")],
                    colorData.all[,c("majorID","labCol")],
                    colorData.myeloid[,c("majorID_sub","labCol")],
                    colorData.cyto[,c("classification","labCol")],
                    colorData.bcell[,c("pctID","labCol")]))
        
labColData <- labColData[!duplicated(labColData$majorID_sub),]
colnames(labColData)[1] <- "classification"
colourz <- colourz %>% left_join(labColData, by = "classification", keep = F)
        
#clean things up for plotting and downstream table
getSomeOrder <- as.data.frame(levels(seu.obj.all$finalID)) %>% left_join(colourz, by = c("levels(seu.obj.all$finalID)" = "classificationFinal"), keep = F)
getSomeOrder <- getSomeOrder[!duplicated(getSomeOrder$classification),]

#hardcode to set up to remove low quality cells
clusterID_final <- table(seu.obj.all$finalID) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
        mutate(clusterID_final=case_when(row_number() < 15 ~ row_number()-1,
                                         row_number() == 15 ~ 100,
                                         row_number() > 15 ~ row_number()-2)) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1

Idents(seu.obj.all) <- "finalID"
seu.obj.all <- RenameIdents(seu.obj.all, newID)
table(Idents(seu.obj.all))
seu.obj.all$clusterID_final <- Idents(seu.obj.all)

getSomeOrder <- clusterID_final %>% left_join(getSomeOrder, by = c("Var1" = "levels(seu.obj.all$finalID)"), keep = F) %>% arrange(clusterID_final)

getSomeOrder[getSomeOrder$clusterID_final == 100,][7] <- NA

plot <- DimPlot(seu.obj.all, 
              reduction = "umap", 
              group.by = "clusterID_final",
              cols = getSomeOrder$colour,
              pt.size = 0.25,
              label = T,
              label.box = T
)

fig7c <- cusLabels(plot = plot, shape = 21, labCol = getSomeOrder$labCol, size = 8, 
                   alpha = 0.8, rm.na = T, nudge_x = c(rep(0,27),-0.45,rep(0,8)), nudge_y = c(rep(0,27),1.15,rep(0,8)))
ggsave("./output/umap_final_hVoWadj_CLEAN_introns.png")

#remove Poor quality cells from table        
seu.obj.all <- subset(seu.obj.all, subset = finalID != "Low quality")

#do age matched stats
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

ggsave("./output/supp_tableStats.png",width = 24, height = 12)

coorData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
                                                   pct = Freq/samSize*100,
                                                   cellSource = ifelse(grepl("ealthy",Var2),"Healthy","OS")) %>% as.data.frame()
ageData <- read.csv("refColz.csv")
ageData$orig.ident <- factor(ageData$orig.ident)
ageData$cellSource <- ifelse(grepl("ealthy", ageData$orig.ident), "Healthy", "Osteosarcoma")

ageRes <- ggpubr::compare_means(age ~ cellSource, ageData) %>% as.data.frame()

p <- ggplot(ageData, aes_string(y = "age", x = "cellSource")) + 
        labs(x = NULL, y = "Age (years)") + 
         geom_boxplot(aes_string(x = "cellSource"), alpha = 0.25, outlier.color = NA) + 
        geom_point(size = 2, position = position_jitter(width = 0.25),
                   aes_string(x = "cellSource", y = "age")) +
        theme_bw() + 
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(), 
              strip.background = element_rect(fill = NA, color = NA), 
              strip.text = element_text(face = "bold"), 
              axis.ticks.x = element_blank(), 
              axis.text = element_text(color = "black")         
             )
ggsave("./output/supp_age_Stats.png")


coorData <- coorData %>% left_join(ageData[,-5], by = c("Var2" = "orig.ident"))

p <- ggplot(coorData, aes(x=age, y=pct)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
                #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point(aes(color = cellSource)) +
                #scale_x_continuous(limits = c(0, 100),expand = c(0, 0)) +
                #scale_y_continuous(limits = c(0, 100),expand = c(0, 0)) +
                labs(x = "Age (years)", y = "Percentage") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                #geom_text(x = 75, y = 5, label = as.character(as.expression(eq)), parse = TRUE) +
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
                      ) 
                
pi <- p + facet_wrap("Var1", scales = "free_y", nrow = 8) + stat_cor() + theme(plot.background = element_rect(fill = "white"))

ggsave("./output/cellByAgeCoor.png", height = 12, width = 18)

allData <- statData %>% group_by(Var1) %>% summarize(avg = mean(pct),
                                  std = sd(pct)
                                 ) %>% mutate(cellSource = "All") %>% relocate(cellSource, .after = Var1) %>% as.data.frame()

#manually add major groups
#write.csv(getSomeOrder, "majorGroups.csv", row.names= F)

groupData <- read.csv("majorGroups.csv", header = T)
        
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

write.csv(compiledData, file = "./output/compiledTableData.csv")
        
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
        
write.csv(compiledData, file = "./output/header_compiledTableData.csv")


#### coor plot ####

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
                
ggsave("./output/suppFlowCoorPlot.png", height = 5, width = 7)

##### create input ref for cybersort #####
        
createCYBRsort(seu.obj = seu.obj.all, groupBy = "classificationFinal", downSample = F, 
         outName = "cybr_no_ds_NORM", outDir = "./output/cyberSort/"
                    )
        
createCYBRsort(seu.obj = seu.obj.all, groupBy = "classificationFinal", downSample = T, 
         outName = "cybr_with_ds_NORM", outDir = "./output/cyberSort/"
                    )

vilnPlots(seu.obj = seu.obj.all, groupBy = "classificationFinal", numOfFeats = 24, outName = "classificationFinal",
                      outDir = "./output/vilnFinalClass/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )
        
#load viln output to filter for defining features
classy.df <- read.csv("./output/vilnFinalClass/classificationFinal_gene_list.csv", row.names = 1, header = T)
classy.df <- classy.df[classy.df$avg_log2FC > 0.85 & classy.df$pct.1 > 0.75, ]
classy.df <- classy.df %>% filter(duplicated(gene) == FALSE)
        
#load a filter .csv files for cyber sort
cyberDS.df <- read.csv("./output/cyberSort/cybr_with_ds_NORM_cyberSort_matrix.csv", row.names = 1, header = T)
cyberDS.df.filtered <- cyberDS.df[row.names(cyberDS.df) %in% classy.df$gene,]
write.csv(cyberDS.df.filtered, file = "./output/cyberSort/cyberDS_NORM_df_filtered_cyberSort_matrix.csv")

        
cyberNoDS.df <- read.csv("./output/cyberSort/cybr_no_ds_NORM_cyberSort_matrix.csv", row.names = 1, header = T)
cyberNoDS.df.filtered <- cyberNoDS.df[row.names(cyberDS.df) %in% classy.df$gene,]
write.csv(cyberNoDS.df.filtered, file = "./output/cyberSort/cyberNoDS_NORM_df_filtered_cyberSort_matrix.csv")

        
#load in new class
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./finalUMAP_hVoWadj_final.csv", groupBy = "classificationFinal", metaAdd = "cyberGroup")
        
vilnPlots(seu.obj = seu.obj.all, groupBy = "cyberGroup", numOfFeats = 24, outName = "cyberGroup",
                      outDir = "./output/vilnFinalClass/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )
        
        
createCYBRsort(seu.obj = seu.obj.all, groupBy = "cyberGroup", downSample = T, 
         outName = "cyberGroup_ds_NORM", outDir = "./output/cyberSort/"
                    )

createCYBRsort(seu.obj = seu.obj.all.sub, groupBy = "cyberGroup", downSample = T, 
         outName = "cyberGroup_NO_OTHERS_ds_NORM", outDir = "./output/cyberSort/"
                    )
        
classy.df <- read.csv("./output/vilnFinalClass/cyberGroup_gene_list.csv", row.names = 1, header = T)
classy.df <- classy.df[classy.df$avg_log2FC > 0.85 & classy.df$pct.1 > 0.75, ]
classy.df <- classy.df %>% filter(duplicated(gene) == FALSE)
        
        
cyberDS.df <- read.csv("./output/cyberSort/cyberGroup_ds_NORM_cyberSort_matrix.csv", row.names = 1, header = T)
cyberDS.df.filtered <- cyberDS.df[row.names(cyberDS.df) %in% classy.df$gene,]
write.csv(cyberDS.df.filtered, file = "./output/cyberSort/cyberGroup_ds_NORM_df_filtered_cyberSort_matrix.csv")

classy.df <- read.csv("./output/vilnFinalClass/cyberGroup_gene_list.csv", row.names = 1, header = T)
classy.df <- classy.df[classy.df$avg_log2FC > 0.85 & classy.df$pct.1 > 0.75 & classy.df$cluster != "Other", ]
classy.df <- classy.df %>% filter(duplicated(gene) == FALSE)
        
cyberDS.df <- read.csv("./output/cyberSort/cyberGroup_NO_OTHERS_ds_NORM_cyberSort_matrix.csv", row.names = 1, header = T)
cyberDS.df.filtered <- cyberDS.df[row.names(cyberDS.df) %in% classy.df$gene,]
write.csv(cyberDS.df.filtered, file = "./output/cyberSort/cyberGroup_NO_OTHERS_ds_NORM_df_filtered_cyberSort_matrix.csv")
        
        
#saveRDS(seu.obj.all, file = "./output/processedData.rds")
#seu.obj.all <- readRDS(file = "./output/processedData.rds")


###### generate feat list #####
vilnPlots(seu.obj = seu.obj.all, groupBy = "finalID", numOfFeats = 24, outName = "HvO_cell.l2",
                      outDir = "./output/viln_finalID_HvO_cell.l2/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

vilnPlots(seu.obj = seu.obj.all, groupBy = "majorID", numOfFeats = 24, outName = "HvO_cell.l1",
                      outDir = "./output/viln_finalID_HvO_cell.l1/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


seu.final.H <- subset(seu.obj.all,
                  subset = 
                  cellSource ==  "Healthy") 

vilnPlots(seu.obj = seu.final.H, groupBy = "majorID", numOfFeats = 24, outName = "H_cell.l1",
                      outDir = "./output/viln_finalID_H_cell.l1/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

vilnPlots(seu.obj = seu.final.H, groupBy = "finalID", numOfFeats = 24, outName = "H_cell.l2",
                      outDir = "./output/viln_finalID_H_cell.l2/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )



#Extract data for key feat table (supp table 3)
featOut.df <- read.csv("./output/viln_finalID_H_cell.l2/H_cell.l2_gene_list.csv")
featOutPROCESSED.df <- as.data.frame(lapply(split(featOut.df$gene,featOut.df$cluster), head, n = 50))
write.csv(featOutPROCESSED.df, file = "./output/viln_finalID_H_cell.l2/cleanGeneList.csv")

featOut.df <- read.csv("./output/viln_finalID_HvO_cell.l2/HvO_cell.l2_gene_list.csv")
featOutPROCESSED.df <- as.data.frame(lapply(split(featOut.df$gene,featOut.df$cluster), head, n = 50))
write.csv(featOutPROCESSED.df, file = "./output/viln_finalID_HvO_cell.l2/cleanGeneList.csv")


#evalutate clssification accuracy
#use custom classifier
#geneLists <- read.csv("./output/seurat_h5/seurath5_l3_gene_list.csv")
geneLists <- read.csv(file = "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_HvO_cell.l2/HvO_cell.l2_gene_list.csv")

datas <- geneLists[,c("cluster","gene")]
colnames(datas) <- c("gs_name", "gene_symbol")
datas <- datas %>% group_by(gs_name) %>% top_n(50) %>% dplyr::distinct(gene_symbol) %>% as.data.frame()

#clus.markers <- read.csv(file = "/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_HvO_cell.l2/HvO_cell.l2_gene_list.csv")
#clus.markers <- read.csv("./output/seurat_h5/seurath5_l3_gene_list.csv")
clus.markers <- read.csv(file = "/pl/active/dow_lab/dylan/k9_atlas_scRNA/Adam_BM_CD34_rds/BM_gene_list.csv")
clusters <- unique(clus.markers$cluster)

df.list <- list()
for (cluster in clusters) {
    clus_sub <- clus.markers[clus.markers$cluster == cluster, ]

    #run enrichR
    enriched <- as.data.frame(enricher(gene = clus_sub$gene, TERM2GENE = datas, pvalueCutoff = 1))
    if(nrow(enriched) > 0){
        enriched$cluster <- cluster
        enriched <- head(enriched) #only takes the top 6 terms - can modify if desired
        df.list[[which(cluster == clusters)]] <- enriched
    }
}

cellCalls <- do.call(rbind, df.list)
outfile <- paste("./output/cell_classification_azTest.csv", sep = "")
write.csv(cellCalls, file = outfile)

###comment### this out if you are using "goodTerms"
cellCalls_clean <- cellCalls

######### FINAL TWEAK TO df ##########


plot <- ggplot(data = cellCalls_clean, mapping = aes_string(x = 'cluster', y = 'ID')) +
    geom_point(mapping = aes_string(size = 'Count', color = -log2(cellCalls_clean$p.adjust))) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t = 150,
                               b = 5),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "transparent",colour = NA),
          panel.grid.major = element_line(color = "gray"), 
          panel.grid.minor = element_line(color = "gray")
          ) + 
    scale_colour_viridis(option="magma", name='-log2(padj)') +
    guides(size=guide_legend(title="Gene count", override.aes=list(fill=NA))) +
    #geom_tile(aes(fill = cluster, y = 0), size = 1, show.legend = FALSE) + #add this 
    #geom_tile(aes(fill = cluster, y = as.numeric(length(unique(PanglaoDB_Augmented_2021.Term)))+1), size = 1, show.legend = FALSE) + #and this back in if you like tiles
    #geom_point(aes(y = 0, x = cluster, size = max(Count), fill = cluster), shape = 21, stroke = 0.75, show.legend = FALSE) +
    #geom_point(aes(y = as.numeric(length(unique(ID)))+1, x=cluster, size = max(Count), fill = cluster), shape = 21, stroke = 0.75, show.legend = FALSE) +
    #geom_text(aes(y = 0, label = cluster), size = 3.5) +
    geom_text(aes(y = as.numeric(length(unique(ID))), label = cluster), size = 3.5,vjust = 0, angle=90, hjust = 0) +
    #scale_y_discrete(limits = rev) #+
    coord_cartesian(expand = TRUE, clip = "off")

#check path is correct
ggsave("./output/supp_gsea_bm.png", width = 14, height = 10)



#create annotated reference files
seu.obj.all@meta.data <- seu.obj.all@meta.data[,!grepl("DF|pANN", colnames(seu.obj.all@meta.data))]
colnames(seu.obj.all@meta.data)[23] <- "celltype.l1"
colnames(seu.obj.all@meta.data)[25] <- "celltype.l2"
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj.all@meta.data <- seu.obj.all@meta.data[,-c(7,8,9,11,15)]

#reorder for readability
seu.obj.all@meta.data <- seu.obj.all@meta.data[,c(seq(2,12),16,13,17,1,22,15,14,21,19,18,20)]

saveRDS(seu.obj.all, file = "./output/s3/final_dataSet_HvO.rds")

seu.final.H <- subset(seu.obj.all,
                  subset = 
                  cellSource ==  "Healthy") 

saveRDS(seu.final.H, file = "./output/s3/final_dataSet_H.rds")

### ref map ### testing
reference <- LoadH5Seurat("./pbmc_multimodal.h5seurat")
#reference <- LoadH5Seurat("./output/seu_final_healthy.h5seurat")
seu.all <- readRDS(file = "./output/s3/seu.obj.allFinal.rds")
# seu.all <- subset(seu.all,
#                   subset = 
#                   cellSource ==  "Osteosarcoma") 

# seu.all <- readRDS(file = "/pl/active/dow_lab/dylan/k9_atlas_scRNA/Adam_BM_CD34_rds/outputcombined_BM_res0.8_dims45_S3.rds")

anchors <- FindTransferAnchors(reference = reference, query = seu.all,
    dims = 1:50, reference.reduction = "pca", normalization.method = "SCT")
predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l2,
    dims = 1:50)
seu.all <- AddMetaData(seu.all, metadata = predictions)

p <- DimPlot(seu.all, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave("./output/cd8_azimuthMap_cellTypel3.png")
# seu.all@meta.data <- seu.all@meta.data[,!grepl("prediction.score", colnames(seu.all@meta.data))]

features <- colnames(seu.all@meta.data)[grepl("prediction.score", colnames(seu.all@meta.data))]
titles <- unlist(lapply(features, substring, 18))
supp <- prettyFeats(seu.obj = seu.all, nrow = 7, ncol = 10, features = features, order = F,titles = titles) 
ggsave("./output/azimuth_predScore_cellTypel3.png", width=30, height=15)

reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")

DefaultAssay(reference) <- "integrated"

anchors <- FindTransferAnchors(
  reference = reference,
  query = seu.all,
  normalization.method = "SCT",
    reference.reduction = "pca", 
   
    dims= 1:50
)

# ,
#     reference.reduction = "umap",
#     dims= 1:2

seu.all <- MapQuery(
    anchorset = anchors,
    query = seu.all,
    reference = reference,
    refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2"
    )#,
   # reference.reduction = "pca", #i dont know why these dont work....
   # reduction.model = "umap"  #i dont know why these dont work....
)

p1 = DimPlot(seu.all, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(seu.all, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2


Idents(reference) <- "orig.ident"
reference_0 <- subset(reference, idents = c("P1_0", "P2_0", "P3_0", "P4_0", "P5_0", "P6_0", "P7_0", "P8_0"))
reference_0[["percent.mt"]] <- PercentageFeatureSet(reference_0, pattern = "^MT-")
seu.reference_0.list <- SplitObject(reference_0, split.by = "orig.ident")

seu.obj <- readRDS("./output/s3/h7_CLEAN_introns_res1.6_dims45_S3.rds")
seu.h7.list <- SplitObject(seu.obj, split.by = "orig.ident")

seu.list <- c(seu.h7.list,seu.reference_0.list)
        
indReClus(seu.list = seu.list, outDir = "./output/", subName = "h7WpbmcMulti", preSub = T,
                      vars.to.regress = "percent.mt"
                       )
        
#logging data
system("cat /pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
sessionInfo()
