#!/usr/bin/Rscript
##################################################################################################################
#       ### Step 2: read in cleaned up individual files then integreated them into one seurat object ###         #
#                                                                                                                #
# Required files: after running seurat_step1_v1.R the required files should be located in the "output" directory #
#  >> Put this script in the in the base directory -- this sets working directory as needed for proper run       #
#                                                                                                                #
# Step 2 will:                                                                                                   #
#  >> loop through all the preprocessed seurat ojects created in step 1                                          #
#  >> all files will be integrated by finding anchors then compelting SCT normalization/intergration             #
#  >> the output will be located in the "output" directory called seu.integrated.obj_S2.rds                      #
#  >> WARNING: this step is very memory intensive frequently over 100GB of ram required!                         #
#     >> it is recommened to us a high memory compute node to complete this step (i.e. smem node)               #
##################################################################################################################

#load required packages
library(Seurat)
library(tidyverse)

######################## set pwds ####################

din <- "output_S1_wPalpct"
dout <- "output"
outName <- "221005_hVoWadj_regPal_usingPalpct"

#################################################################

#get seurat objects to process and integrate
fpath <- paste("./", din,"/", sep = "") 
files <- list.files(path = fpath, pattern= "S1.rds", all.files=FALSE,
    full.names=TRUE)

create.seu.call <- function(x) {
  readRDS(x)
}

#load all seurat objects and store in a large list
seu.obj <- mapply(create.seu.call, files)

#use SCT inegration to merge the samples
#regress out percent.mt; can regress other things (like HBM) if desired
seu.obj <- lapply(seu.obj, 
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
rm(seu.integrated)
gc()

#integrate data and keep full gene set - still might not be retaining all genes
seu.integrated.obj <- IntegrateData(
  anchorset = seu.integrated.anchors,
  normalization.method = "SCT",
  verbose = TRUE
)

#clean up environment a bit
rm(seu.integrated.anchors)
gc()
seu.integrated.obj <- RunPCA(seu.integrated.obj )

p <- ElbowPlot(seu.integrated.obj, ndims = 50)
outfile <- paste("./",dout,"/", outName, "_seu.integrated.obj_S2_elbow.png", sep = "")
ggsave(outfile)

DefaultAssay(seu.integrated.obj) <- "integrated"

outfile <- paste("./",dout,"/", outName, "_seu.integrated.obj_S2.rds", sep = "")
saveRDS(seu.integrated.obj, file = outfile)
