# Canine_Leukocyte_scRNA

This GitHub repository contains all the analysis code used in, "A single cell RNA sequencing atlas of circulating leukocytes from healthy and osteosarcoma affected dogs." [add citation]

## Repository goals: 
- provide a resource to make the data generated from this project accessible
- enable reproducible/transparent data reporting
- provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue, contact the correspoining author(s), and/or contact Dylan Ammons at dyammons95@gmail.com.

## File structure:
- [:file\_folder: input](/input) contains relevent metadata files and Cell Ranger output count matrices required for a repoducible run
- [:file\_folder: analysis](/analysis) contains the analysis code and source file used to compelte the data analysis

## Supplemental data and potential uses:
### Cell type annotations with defining markers

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




## Potential uses:
- [:file\_folder: cibersort](/cibersort)
- [:file\_folder: scRNA_referenceMap](/scRNA_referenceMap)
- [:file\_folder: loupeBrowser](/loupeBrowser)
- [:file\_folder: gsea](/gsea)
