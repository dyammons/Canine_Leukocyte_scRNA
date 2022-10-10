#!/usr/bin/Rscript

library(Seurat)
library(ggplot2)

###change###
ident <- "clusterID"
subName <- "hVoWadj_CLEAN_introns"
inFile <- "./output/221005_hVoWadj_regPal_wPalpct_introns_res1.4_dims45_S3.rds"

seu.integrated.obj <- readRDS(file = inFile)

############

#read in data
Idents(seu.integrated.obj) <- ident

seu.sub <- subset(seu.integrated.obj,
                  subset = 
                    clusterID !=  "1" & clusterID != "5" & clusterID !=  "10" & clusterID != "12" 
)   

z <- table(seu.sub@meta.data$orig.ident)
print(z)
k <- ifelse(min(z)>100, 100, min(z))
print(k)

seu.sub.list <- SplitObject(seu.sub, split.by = "orig.ident")

seu.obj <- lapply(seu.sub.list, 
                     SCTransform, 
                     vars.to.regress = c("percent.mt", "percent.pal"),
                     verbose = TRUE,
                     conserve.memory=TRUE)

SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                       nfeatures = 2000)

seu.integrated <- PrepSCTIntegration(
  object.list = seu.obj,
  anchor.features = SelectedFeatures,
  verbose = TRUE
)

gc()

seu.integrated.anchors <- FindIntegrationAnchors(
  object.list = seu.integrated,
  normalization.method = "SCT",
  anchor.features = SelectedFeatures,
  verbose = TRUE
)

#clean up environment a bit
rm(seu.integrated.obj)
rm(seu.sub)
rm(seu.sub.list)
rm(seu.obj)
rm(seu.integrated)
gc()

#create list of genes to keep - still might not be retaining all genes for feat plots
#to_integrate <- Reduce(intersect, lapply(seu.integrated.anchors@object.list, rownames)) #skipped for hvo helpers

#integrate data and keep full gene set - still might not be retaining all genes
seu.integrated.obj <- IntegrateData(
  anchorset = seu.integrated.anchors,
  normalization.method = "SCT",
  #features.to.integrate = to_integrate,
  k.weight = k,
  verbose = TRUE
)

#clean up environment a bit
rm(seu.integrated.anchors)
gc()
seu.integrated.obj <- RunPCA(seu.integrated.obj )

outfile <- paste("./output/", subName,"_S2_elbow.png", sep = "")
p <- ElbowPlot(seu.integrated.obj, ndims = 50)
ggsave(outfile)

DefaultAssay(seu.integrated.obj) <- "integrated"

outfile <- paste("./output/", subName,"_S2.rds", sep = "")
saveRDS(seu.integrated.obj, file = outfile)
