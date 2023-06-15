# Canine_Leukocyte_scRNA
[![DOI](https://zenodo.org/badge/546752777.svg)](https://zenodo.org/badge/latestdoi/546752777)

This GitHub repository contains all the analysis code used in, "A single-cell RNA sequencing atlas of circulating leukocytes from healthy and osteosarcoma affected dogs."

If you use our raw/processed data, extract data using the UCSC Cell Browser portal, or use portions of our code in your analysis, please cite:
> Ammons DT, Harris RA, Hopkins LS, Kurihara J, Weishaar K and Dow S (2023) A single-cell RNA sequencing atlas of circulating leukocytes from healthy and
osteosarcoma affected dogs. Front. Immunol. 14:1162700. doi: 10.3389/fimmu.2023.1162700

## Repository goals: 
- provide a resource to make the data generated from this project accessible
- enable reproducible/transparent data reporting
- provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue, contact the corresponding author(s), and/or contact Dylan Ammons at dylan.ammons @ colostate dot edu.

## File structure:
- [:file\_folder: input](/input) contains relevant metadata files and instructions for obtaining data associated with this study
- [:file\_folder: analysis](/analysis) contains the analysis code and source file used to complete the data analysis

## Supplemental data and potential uses:
1. [Browse the data](#1-browse-the-complete-annotated-dataset)
2. [Cell type annotations](#2-cell-type-annotations-with-defining-markers)
3. [Reference Mapping](#3-using-the-data-to-complete-reference-mapping)
4. [GSEA using dataset](#4-gene-set-enrichment-analysis)
5. [Module scoring](#5-module-scoring)
6. [CIBERSORTx](#6-cibersortx)

### 1. Browse the complete annotated dataset

The proccessed dataset is avaliable for browsing via the UCSC Cell Browser portal.
Using the portal you can explore feature expression throughout the dataset as well as obtain the transcriptomic signatures of each cell type though an interactive webpage.

Note: the cell type gene lists on UCSC Cell Browser are ordered by P value in descending order by default, you can toggle it to ascending to get the enriched markers at the top of the list. I will update if/when the default setting is modified. 

Link to the dataset: https://canine-leukocyte-atlas.cells.ucsc.edu

Link to UCSC Cell Browser documentation: https://cellbrowser.readthedocs.io/en/master/

### 2. Cell type annotations with defining markers

Cell markers lists were curated using 7 healthy canine leukocyte samples. The top 50 defining features (identified using `FindMarkers` for each cell type were considered, with the top 24 features evaluated for specificity using violin plots and preference given to unique features only found in the top 50 of one cell type.

<details><summary>cellType.l1</summary>
<p>
  
|Cell Type         |Marker                                                             |
|------------------|-------------------------------------------------------------------|
|CD8/NK cell       |CCL5, GZMB, KLRB1, NCR3, GZMA, IL12RB2, TRPM3, KLRK1, KLRD1, CD96  |
|CD4 T cell        |IL7R, LEF1, DOCK3, ATP10A, CD52, CSTA, ICOS, CD3E                  |
|Monocyte          |LYZ, MAFB, BPI, FN1, F13A1, TCF7L2, S100P, CD83, NRG1, PLXDC2, MITF|
|Dendritic Cell    |IGF1, NCAM2, FLT3, RAB3C, FCER1A, FGF12, NRXN1, TLR3               |
|Neutrophil        |S100A12, CYP4F22, RGS2, CD4, ALDH1A2, SGK1, IL1R2                  |
|Granulocyte       |DACH1, TGM2, ADGRE2, CA8, SMPD1, IL5RA                             |
|B cell            |TNFRSF13C, PAX5, EBF1, BTLA, HTR1F, NRIP1, MS4A1, DLA-DRA          |
|Plasma cell       |JCHAIN, TXNDC5, IRF4, TNFRSF17, DERL3, CCR10, MPP6, TNFRSF13B      |
|DN T cell         |KANK1, TMEM132D, NMB, SYNJ2, GZMK, LEF1                            |
|gd T cell         |PDE11A, PSD3, RHEX, IL17RB, CDH4, GATA3, FAT1, ZNF683              |
|Cycling T cell    |TOP2A, MKI67, RRM2, H1-5, DIAPH3, TK1, KIF11, TPX2                 |
|CD34+ Unclassified|ZNF521, CD34, KIT, CD109, IGF2BP2                                  |
  
</p>
</details>

<details><summary>cellType.l2</summary>
<p>

|Cell type         |                                                                    |Marker                                                              |
|------------------|--------------------------------------------------------------------|--------------------------------------------------------------------|
|CD8 T cell        |                                                                    |                                                                    |
|                  |CD8 T cell                                                          |CCL5, GZMB, TRPM3, KLRK1, KLRB1, CCL4, IL2RB, KLRD1, CD96           |
|                  |CD8+ Memory                                                         |GZMK, KLRK1, GZMB, KLRB1, KLRD1, KLRG1, CD8A, FASLG, CCR5           |
|CD4 T cell        |                                                                    |                                                                    |
|                  |CD4+ Naive                                                          |ATP10A, RGS10, LEF1, CTPS1, ZNF536, SATB1, CSTA, ITGA1, COL6A5      |
|                  |CD4+ T reg                                                          |IKZF2, S100A5, CTLA4, GATA3, IL2RA, TOX, CD80, ZNF831               |
|                  |CD4+ TCM                                                            |LEF1, CCR7, TSHZ2, LTB, CD52, CSTA, IL7R, TCF7                      |
|                  |CD4+ TEM                                                            |IL7R, RORA, PLCL1, MAF, SLC9A9, ICOS, GALNT17, PTPN13, DOCK3        |
|NK cell           |                                                                    |                                                                    |
|                  |NK T cell                                                           |TGFBR3, GPA33, RARRES1, IL2RB, IKZF2, KLRK1, CD96                   |
|                  |NK cell                                                             |GZMA, PI3, IL2RB, KLRF1, PTPRM, PAX4, CD96, IGSF3, STMN2, F2RL3     |
|Monocyte/DC       |                                                                    |                                                                    |
|                  |DC                                                                  |PKIB, TCF4, FLT3, HAVCR1, SDC2, NCAM2, MRC1, IGF1, FCER1A           |
|                  |M-MDSC                                                              |IL18, IL1B, MEFV, CPXM2, LTF, CYP4F22, STEAP4, KCNJ2, S100A12       |
|                  |Monocyte                                                            |ARHGAP45, F13A1, CD83, FN1, LYZ, NRG1, RETN, SLAMF9, TCF7L2         |
|Granulocyte       |                                                                    |                                                                    |
|                  |Eosinophil                                                          |CSTB, SLCO4C1, TGM2, ADGRE2, C30H15orf48, PTPRN2, SMPD1, CA8        |
|                  |PMN-MDSC                                                            |PGLYRP1, CAMP, MMP9, CRISP2, TCN1, MMP8, CD177, FADS1, LTF, S100A12 |
|                  |Neutrophil                                                          |CYP4F22, S100A12, CD4, SOD2, SERPINA1, IL1R2, ALDH1A2, S100A8       |
|                  |Basophil                                                            |SLCO4C1, IL5RA, DAPK2, DACH1, CA8, ADGRE2, EEPD1, ANKRD33B, HK2     |
|B cell            |                                                                    |                                                                    |
|                  |B cell                                                              |TNFRSF13C, PAX5, BANK1, EBF1, BTLA, PLEKHO1, HTR1F, NRIP1           |
|                  |Plasma cell                                                         |MZB1, JCHAIN, LMAN1, TXNDC5, LAP3, IRF4, RARRES2, TNFRSF17          |
|Miscellaneous     |                                                                    |                                                                    |
|                  |DN T cell                                                           |KIAA0825, TMEM132D, NMB, SYNJ2, GZMK, MYB, KANK1, PLCL1, SLF1, CTLA4|
|                  |gd T cell                                                           |PDE11A, PSD3, CRLF2, IL17RB, VSTM4, RHEX, FAT1, TOX2                |
|                  |Other T cell                                                        |RRM2, MKI67, SPC24, TK1, CENPF, TOP2A, CLSPN, NCAPG                 |
|                  |CD34+ Unclassified                                                  |CD34, NDST3, TFPI, CLEC3B, KIT, NAV3, CD109                         |

</p>
</details>

<details open><summary>cellType.l3 (Default)</summary>
<p>

|Cell Type     |                       |Marker                                                                   |
|--------------|-----------------------|-------------------------------------------------------------------------|
|CD8 T cell    |                       |                                                                         |
|              |CD8+ Naive             |ITGA1, LEF1, PTGDR, IL2RB, ADGRG1, NBEA                                  |
|              |CD8+ Effector          |CCL5, TRPM3, IL12RB2, GZMB, KLRB1, GZMA, NCR3, IL2RB, KLRD1, CD96        |
|              |CD8+ Memory            |GZMK, GZMB, PI3, BTBD11, CTSW, CCR5, CCL4, KLRG1, FASLG                  |
|              |CD8+ gd T cell         |PTHLH, IGF2BP2, ABTB2, AKAP12, SOX4, CTSW, SLC16A10, PXT1, ZNRF3, SULT2B1|
|CD4 T cell    |                       |                                                                         |
|              |CD4+ Naive             |LEF1, CSTA, RGS10, ZNF536, CCR7, COL6A5, LTB, TNFSF8                     |
|              |CD4+ TCM               |LEF1, TSHZ2, CD52, CCR7, IL7R, CTPS1, EFHC2, CARMIL1                     |
|              |CD4+ TEM               |IL7R, SLC9A9, ICOS, MAF, CD28, SKAP1, CD40LG                             |
|              |CD4+ TEM, Th1-like     |IL7R, PTPN13, IL18R1, CD28, RCAN2, CCR9, CCR5, IL12RB2, CD52, PRUNE2     |
|              |CD4+ TEM, Th2-like     |RNF220, ITGA2, GATA3, CCDC3, LGALS3, PTPN13, S100A2, PPEF1, CMA1         |
|              |CD4+ TEM, Th17-like    |NTRK2, PTPN13, ADAM12, NRG2, RGS17, DNAH8, CCR6, NPAS2, RORA, LTBP1      |
|              |CD4+ T reg             |IKZF2, CTLA4, RGS1, ICOS, IL2RA, CD28, ZNF831                            |
|              |CD4+, IFN signature    |CXCL10, IFI44, OAS1, ISG15, IFI44L, IFGGB2, CTLA4, STAT1, DDX58, XAF1    |
|Monocyte      |                       |                                                                         |
|              |Monocyte, CD4-         |LYZ, BPI, LRMDA, MT2A, F13A1, FN1, NRG1, CCDC88A, CD83, RETN             |
|              |Monocyte, CD4+         |IL1B, MAFB, NFKBIA, CXCL8, FN1, BLOC1S6, CD83, S100P, BPI, NRG1          |
|              |M-MDSC                 |IL18, IL1B, LTF, MEFV, KCNJ2, CPXM2, S100A12, STEAP4, CSF3R, IL31RA      |
|              |Monocyte, IFN signature|RSAD2, OAS1, OAS2, DDX58, HERC6, OAS3, RTP4, EIF2AK2, IFIT2              |
|Dendritic cell|                       |                                                                         |
|              |Pre-DC                 |FGF12, GPHA2, MTUS2, FCER1A, PLCE1, PTPRS, IGF1, NECTIN1, IL3RA, AK8     |
|              |Myeloid cDC1           |ZNF366, SDC2, DISC1, ECRG4, TMEM163, RIMS2, KIT, OTOF, RTKN2, RAB7B      |
|              |Myeloid cDC2           |PKIB, CD300H, SDC2, CD1C, NCAM2, CD86, BATF3, ZNF366, PID1, ECM1         |
|              |Plasmacytoid DC        |COBLL1, RAB3C, IGF1, FCER1A, RYR1, PRKG1, CCND1, STYXL2, ANK1, OCIAD2    |
|              |Unclassified DC        |PLCB4, ZNF366, KCNK13, STRIP2, SDC2, OTOF, HACD1, C5, SLC8A1, CNTLN      |
|Granulocyte   |                       |                                                                         |
|              |Neutrophil             |S100A12, CD4, SERPINA1, SGK1, S100A8, ALDH1A2, FNDC3B, GGH, SRGN, IL1R2  |
|              |PMN-MDSC               |CAMP, PGLYRP1, CRISP2, MMP9, MMP8, TCN1, CD177, LTF, FADS1, S100A12      |
|              |Eosinophil             |C30H15orf48, TGM2, DACH1, PADI3, SMPD1, CA8, IL5RA                       |
|              |Basophil               |DACH1, CA8, IL5RA, DAPK2, TGFA, ANKRD33B, HK2, PRR5L                     |
|B cell        |                       |                                                                         |
|              |Immature B cell        |SYT1, PAX5, VPREB3, ERC2, TMTC2, KLHL14, F8, TEX9, TDRP, ADGRF1          |
|              |Naive B cell           |TNFRSF13C, BANK1, HTR1F, PAX5, EBF1, BTLA, NRIP1, ADAM9                  |
|              |Class switched B cell  |TNFRSF13C, GOLM1, BANK1, BTLA, EBF1, DYNC1I1, MTMR2, PAX5                |
|              |Activated B cell       |IGKC, CACNB2, PAX5, TNFRSF13C, IGHM, RASGRF2, AOX2, BCAR3, ADAM32        |
|              |Plasma cell            |JCHAIN, MZB1, TXNDC5, LMAN1, FKBP11, LAP3, DERL3, CCR10, MKI67, TNFRSF13B|
|Miscellaneous |                       |                                                                         |
|              |DN T cell              |KIAA0825, TMEM132D, KANK1, NMB, CTLA4, SYNJ2, BICDL1, SLF1, ID3, KIAA1549|
|              |gd T cell              |PARD3B, RHEX, IL17RB, CDH4, GATA3, FAT1, TOX2, ADARB1, ZNF683, TGFBR3    |
|              |NK cell                |KLRF1, STMN2, PAX4, NCR3, F2RL3, CD96, IL2RB, IGSF3, FREM1, FASLG        |
|              |Cycling T cell         |TOP2A, MKI67, RRM2, H1-5, DIAPH3, TK1, KIF11, TPX2, ASPM                 |
|              |NK T cell              |GPA33, TGFBR3, KLRK1, CD96, SYTL2, MOV10L1, SLA2, DSTN, RARRES1          |
|              |CD34+ Unclassified     |TFPI, ZNF521, CD34, NDST3, GUCY1A1, HPGD, CLEC3B, KIT, CD109, DNTT       |

</p>
</details>

### 3. Using the data to complete reference mapping
Reference mapping is useful tool to facilitate the identification of cell types in single cell datasets. The approach described here uses Seurat functions to identify anchors between a query dataset (external/personal data) and the reference datasets generated in this study. The default approach describes how to use the healthy only dataset, but it will also work with the combined dataset if you load that file in as the reference.

Before running the reference mapping code, a Seurat object need to be preprocessed and stored as an object named `seu.obj`.
```r
#set the path to the location in which the reference file is saved
reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds")

#prepare the reference
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

#find conserved anchors with query and reference
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims= 1:50
)

#select meta.data slot to use for label transfer -- change refdata value to use alternate labels (i.e., refdata = reference$celltype.l1)
predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l3,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)

#generate and save the image
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
ggsave("./output/referenceMap.png", width = 7, height = 7)
```

### 4. Gene set enrichment analysis

The data generated from this work have the potential to provide supporting evidence to evaluate/confirm the cell identity of sorted bulk RNA sequencing dataset. One approach to do this is to use gene set enrichment analysis (GSEA) with the terms representing the cell type identified in our dataset.

Required input: a list of gene symbols that you wish to query. In this case the genelists ate stored in a dataframe called `clus.markers`

These gene lists could be generated by simply using the features with the highest level of expression after normalizing your dataset, comparing the transcriptome of a cell population of interest (i.e., blood-derived macrophage) verses a reference (i.e., total PBMCs), or any other relevant approach to identify genes of interest.

Example data frame format:

```r
> str(clus.markers)
'data.frame':   400 obs. of  2 variables:
 $ gene   : chr  "B2M" "CD74" "DLA-64" "PPBP" ...
 $ cluster: chr  "sample_1" "sample_2" "sample_3" "sample_4" ...
```

```r
#read in the one of the supplemntal data files provided with the publication
geneLists <- read.csv(file = "./input/supplementalData_1.csv")

#clean the reference
datas <- geneLists[,c("cluster","gene")]
colnames(datas) <- c("gs_name", "gene_symbol")
datas <- datas %>% group_by(gs_name) %>% top_n(50) %>% dplyr::distinct(gene_symbol) %>% as.data.frame()

#run GSEA using clusterProfiler
clusters <- unique(clus.markers$cluster)
df.list <- list()
for (cluster in clusters) {
    clus_sub <- clus.markers[clus.markers$cluster == cluster, ]

    #run enricher
    enriched <- as.data.frame(clusterProfiler::enricher(gene = clus_sub$gene, TERM2GENE = datas, pvalueCutoff = 1))
    if(nrow(enriched) > 0){
        enriched$cluster <- cluster
        enriched <- head(enriched) #only takes the top 6 terms - can modify if desired
        df.list[[which(cluster == clusters)]] <- enriched
    }
}

cellCalls <- do.call(rbind, df.list)
outfile <- paste("./output/cell_classification.csv", sep = "")
write.csv(cellCalls, file = outfile)

#plot the data
plot <- ggplot(data = cellCalls, mapping = aes_string(x = 'cluster', y = 'ID')) +
    geom_point(mapping = aes_string(size = 'Count', color = -log2(cellCalls$p.adjust))) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white",colour = NA),
          plot.background = element_rect(fill = "white",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "transparent",colour = NA),
          panel.grid.major = element_line(color = "gray"), 
          panel.grid.minor = element_line(color = "gray")
          ) + 
    scale_colour_viridis(option="magma", name='-log2(padj)') +
    guides(size=guide_legend(title="Gene count", override.aes=list(fill=NA))) +
    geom_text(aes(y = as.numeric(length(unique(ID))), label = cluster), size = 3.5,vjust = -.3, angle=45, hjust = -0.1) +
    coord_cartesian(expand = TRUE, clip = "off") +
    xlab("Sample") + ylab("GSEA term")

#check path is correct
ggsave("./output/gsea_scRNA_terms.png", width = 6, height = 4)
```

### 5. Module scoring
Module scoring is a supplemental approach that can be applied to single cell datasets with the goal of providing further insights into cell identities. The approach described below uses the Seurat function `AddModuleScore` and the gene lists presented in Table 3 (also found in supplemental data 4) of our associated manuscript. 

The concept of the AddModuleScore() function is similar to GSEA, but also distinct in many ways. Read the [Seurat documentation](https://satijalab.org/seurat/reference/addmodulescore) and/or check out [this webpage](https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/) for more details.

```r
#load in the reference file from supplemental data
ref.df <- read.csv("supplementalData_4", row.names = 1, header = T)

#organize the data
modulez <- split(datas$gene_symbol, datas$cellType_l2)

#complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFeunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              #axis.ticks = element_blank(),
                              #legend.justification = "left",
                              #plot.margin = margin(7, 21, 7, 7, "pt")
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_dots_celltypes.png", sep = ""),width = 10,height=6)
```

### 6. CIBERSORTx

Under development

The reference dataset generated from this study provides the data required to make a comprehensive CIBERSORTx reference to be used for deconvolution of bulk RNA sequencing of peripheral blood. We do not feel comfortable providing a reference that can be directly used, as this study did not benchmark CIBERSORTx references.

We have plans to complete futher work up and release references at a later date. They will be linked here when available.

In the meantime, feel free to generate your own cibersort references using this dataset. Please cite this publication if using the dataset to generate a CIBERSORTx reference. Also, I would recommend completing validation studies to evaluate the performance of the reference before applying it to your system.

With the being said, here is the code I would use to prepare our dataset for CIBERSORTx reference generation:

```
##### create input ref for cibersort #####
#load in processed data from GEO -- see input dir for instructions (healthy only will do, but could used combined dataset if desired
#NOTE: it may work better if only using 1 dog to generate the reference, but the code here does not do that)
seu.obj.h <- readRDS(file = "./output/s3/final_dataSet_H.rds")

#randomly downsample the subset data to obtain equal number of cells in each celltype
set.seed(12)
seu.obj.h <- subset(x = seu.obj.h, downsample = min(table(seu.obj.h$celltype.l1))) #can pull any of the 3 celltype metadata slots
seu.obj.h <- NormalizeData(seu.obj.h)

#extract necessary data
rownameCode <- as.data.frame(seu.obj.h$celltype.l1)
cntMatrix <- seu.obj.h@assays$RNA@data
colnames(rownameCode)[1] <- "name"
colnames(cntMatrix) <- rownameCode$name[match(colnames(cntMatrix), rownames(rownameCode))]
cntMatrix <- as.data.frame(cntMatrix)
cntMatrix$gene <- rownames(cntMatrix)
cntMatrix <- cntMatrix %>% relocate(gene)

#save the matrix in the format required for cibersort -- follow instructions at https://cibersortx.stanford.edu/
write.table(cntMatrix, file="./output/ciberSort/cibersort_forRef_ds.txt", row.names = F, col.names = T, quote=FALSE, sep='\t')
```




