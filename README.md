# Canine_Leukocyte_scRNA

This GitHub repository contains all the analysis code used in, "A single cell RNA sequencing atlas of circulating leukocytes from healthy and osteosarcoma affected dogs." [add citation]

NOTE: to be finalized by time of publication (excuse current typos)
## Repository goals: 
- provide a resource to make the data generated from this project accessible
- enable reproducible/transparent data reporting
- provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue, contact the correspoining author(s), and/or contact Dylan Ammons at dyammons95@gmail.com.

## File structure:
- [:file\_folder: input](/input) contains relevent metadata files and Cell Ranger output count matrices required for a repoducible run
- [:file\_folder: analysis](/analysis) contains the analysis code and source file used to compelte the data analysis

## Supplemental data and potential uses:

1. [Cell type annotations](#1-cell-type-annotations-with-defining-markers)
2. [Reference Mapping](#2-using-the-data-to-complete-reference-mapping)
3. [GSEA using dataset](#3-gene-set-enrichment-analysis)
4. [Module scoring](#4-module-scoring)
5. [CIBERSORT](#4-cibersort)

### 1. Cell type annotations with defining markers

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

|Cell Type         |Marker                                                             |
|------------------|-------------------------------------------------------------------|
|TBD      |TBD  |

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

### 2. Using the data to complete reference mapping
Reference mapping is useful tool to faciliate the identification of cell types in single cell datasets. The apporach described here uses Seurat functions to identify anchors between a query dataset (external/personal data) and the refernce datasets generate in this study. The default approach descrubes how to use the healthy only dataset, but it will also work with the combined dataset if you load that file in as the refernce.

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

### 3. Gene set enrichment analysis

These data have the potential to provide supporting evdience to evaluate cell idniey of sorted bulk RNA sequencing dataset. One apporach to do this is to use GSEA with the terms repsenting the cell type idnified in this dataset.

Required input would be a list a gene names that you wish to query to determine what cell type they most resmeble. These lists could be generated by simply using the features with the highest level of expression after normalizing your dataset, comparing the transcriptome of a cell population of interest (i.e., blood-derrived macrophage) verses a refecne (i.e., total PBMCs), or any other relvent approach.

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

#run GSEA using enrichr
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

### 4. Module scoring

ref.df <- read.csv("/pl/active/dow_lab/dylan/k9_PBMC_scRNA/analysis/output/viln_finalID_H_cell.l3/H_cell.l3_gene_list.csv", row.names = 1, header = T)

datas <- ref.df[,c("cluster","gene")]
colnames(datas) <- c("gs_name", "gene_symbol")
datas <- datas %>% group_by(gs_name) %>% top_n(10) %>% dplyr::distinct(gene_symbol) %>% as.data.frame()

modulez <- split(datas$gene_symbol, datas$gs_name)

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

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

ggsave(paste("./output/", outName, "/", outName, "_dots_drugTargs.png", sep = ""),width = 10,height=6)

### 5. CIBERSORT

Under development
