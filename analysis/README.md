# Analysis code overview
### To complete a reproducible run from raw data:
Retrieve the fastq files from SRA and then align to the canine genome. 
Instructions to download the fastq files can be found in [:file\_folder: input](/input). 
The alignment code is currently not provided, but can be shared if desired. If you want it, you can request through GitHub.

### To complete a reproducible run from cellranger output count matricies:
Download the supplementary zip folder on the NCBI GEO project page as decribed in [:file\_folder: input](/input) directory.
In addition to the count matrices you will also need the .csv files located in the [:file\_folder: sheets](/input/sheets) directory.
Once all files are in place, you will need the [runScript_pbmc.R](/analysis/runScript_pbmc.R) which contains the analysis code and you will have to source the [customFunctions.R](/analysis/customFunctions.R).

NOTE: the file paths in `runScript_pbmc.R` will likely need to be modified when reading/writing data to run on your system.

### To reproduce/explore data using processed data:
Download the supplementary .rds file on the NCBI GEO project page as described in [:file\_folder: input](/input) directory.
There is an .rds file for just the healthy dogs, combined dataset, and a file for each major immune cell subset.
These files can be loaded into [runScript_pbmc.R](/analysis/runScript_pbmc.R) to bypass some of the more computationally intense steps, or they can be used for your own exploration!

### Envrionement to reproduce:
This work was completed on the UC boulder Alpine supercomputer on a Linux operating system. All code was run in a conda environment and I am happy to share the .yml file upon request.

Otherwise, see below for packages loaded into the system.

```r
> sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.4 (Ootpa)

Matrix products: default
BLAS/LAPACK: /projects/dyammons@colostate.edu/software/anaconda/envs/r_env/lib/libopenblasp-r0.3.18.so

locale:
 [1] LC_CTYPE=C.UTF-8    LC_NUMERIC=C        LC_TIME=C          
 [4] LC_COLLATE=C        LC_MONETARY=C       LC_MESSAGES=C      
 [7] LC_PAPER=C          LC_NAME=C           LC_ADDRESS=C       
[10] LC_TELEPHONE=C      LC_MEASUREMENT=C    LC_IDENTIFICATION=C

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggtree_3.2.1                ape_5.6-2                  
 [3] scuttle_1.4.0               scRNAseq_2.8.0             
 [5] ggpubr_0.4.0                slingshot_2.2.1            
 [7] TrajectoryUtils_1.2.0       SingleCellExperiment_1.16.0
 [9] princurve_2.1.6             clusterProfiler_4.2.2      
[11] msigdbr_7.5.1               ggsankey_0.0.99999         
[13] lemon_0.4.5                 reshape_0.8.9              
[15] viridis_0.6.2               viridisLite_0.4.1          
[17] SingleR_1.8.1               SeuratDisk_0.0.0.9019      
[19] RColorBrewer_1.1-3          pheatmap_1.0.12            
[21] DESeq2_1.34.0               SummarizedExperiment_1.24.0
[23] Biobase_2.54.0              MatrixGenerics_1.6.0       
[25] matrixStats_0.62.0          GenomicRanges_1.46.1       
[27] GenomeInfoDb_1.30.1         IRanges_2.28.0             
[29] S4Vectors_0.32.2            BiocGenerics_0.40.0        
[31] colorspace_2.0-3            ggrepel_0.9.1              
[33] cowplot_1.1.1               scales_1.2.1               
[35] patchwork_1.1.2             DoubletFinder_2.0.3        
[37] clustree_0.4.4              ggraph_2.0.5               
[39] forcats_0.5.2               stringr_1.4.1              
[41] dplyr_1.0.10                purrr_0.3.5                
[43] readr_2.1.2                 tidyr_1.2.1                
[45] tibble_3.1.8                ggplot2_3.3.6              
[47] tidyverse_1.3.1             SeuratObject_4.1.3         
[49] Seurat_4.3.0               

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                rtracklayer_1.54.0           
  [3] scattermore_0.8               bit64_4.0.5                  
  [5] knitr_1.40                    irlba_2.3.5                  
  [7] DelayedArray_0.20.0           data.table_1.14.4            
  [9] AnnotationFilter_1.18.0       KEGGREST_1.34.0              
 [11] RCurl_1.98-1.6                generics_0.1.3               
 [13] GenomicFeatures_1.46.5        ScaledMatrix_1.2.0           
 [15] RSQLite_2.2.18                shadowtext_0.1.1             
 [17] RANN_2.6.1                    future_1.29.0                
 [19] bit_4.0.4                     tzdb_0.3.0                   
 [21] enrichplot_1.14.2             spatstat.data_3.0-0          
 [23] xml2_1.3.3                    lubridate_1.9.0              
 [25] httpuv_1.6.5                  assertthat_0.2.1             
 [27] xfun_0.34                     hms_1.1.2                    
 [29] babelgene_22.9                promises_1.2.0.1             
 [31] restfulr_0.0.15               progress_1.2.2               
 [33] fansi_1.0.3                   dbplyr_2.1.1                 
 [35] readxl_1.4.1                  igraph_1.2.11                
 [37] DBI_1.1.3                     geneplotter_1.72.0           
 [39] htmlwidgets_1.5.4             spatstat.geom_3.0-5          
 [41] ellipsis_0.3.2                backports_1.4.1              
 [43] annotate_1.72.0               biomaRt_2.50.3               
 [45] deldir_1.0-6                  sparseMatrixStats_1.6.0      
 [47] vctrs_0.5.1                   ensembldb_2.18.3             
 [49] ROCR_1.0-11                   abind_1.4-5                  
 [51] cachem_1.0.6                  withr_2.5.0                  
 [53] ggforce_0.4.1                 progressr_0.11.0             
 [55] sctransform_0.3.5             GenomicAlignments_1.30.0     
 [57] treeio_1.18.1                 prettyunits_1.1.1            
 [59] goftest_1.2-3                 cluster_2.1.4                
 [61] DOSE_3.20.1                   ExperimentHub_2.2.1          
 [63] lazyeval_0.2.2                crayon_1.5.2                 
 [65] genefilter_1.76.0             hdf5r_1.3.7                  
 [67] spatstat.explore_3.0-5        labeling_0.4.2               
 [69] pkgconfig_2.0.3               tweenr_2.0.2                 
 [71] ProtGenerics_1.26.0           nlme_3.1-160                 
 [73] rlang_1.0.6                   globals_0.16.1               
 [75] lifecycle_1.0.3               miniUI_0.1.1.1               
 [77] downloader_0.4                filelock_1.0.2               
 [79] BiocFileCache_2.2.1           modelr_0.1.8                 
 [81] rsvd_1.0.5                    AnnotationHub_3.2.2          
 [83] cellranger_1.1.0              polyclip_1.10-4              
 [85] lmtest_0.9-40                 Matrix_1.5-1                 
 [87] aplot_0.1.2                   carData_3.0-5                
 [89] zoo_1.8-9                     reprex_2.0.1                 
 [91] ggridges_0.5.4                rjson_0.2.21                 
 [93] png_0.1-7                     bitops_1.0-7                 
 [95] KernSmooth_2.23-20            Biostrings_2.62.0            
 [97] blob_1.2.3                    DelayedMatrixStats_1.16.0    
 [99] qvalue_2.26.0                 parallelly_1.32.1            
[101] spatstat.random_3.1-3         rstatix_0.7.0                
[103] gridGraphics_0.5-1            ggsignif_0.6.4               
[105] beachmat_2.10.0               memoise_2.0.1                
[107] magrittr_2.0.2                plyr_1.8.7                   
[109] ica_1.0-3                     zlibbioc_1.40.0              
[111] scatterpie_0.1.7              compiler_4.1.1               
[113] BiocIO_1.4.0                  fitdistrplus_1.1-8           
[115] Rsamtools_2.10.0              cli_3.6.0                    
[117] XVector_0.34.0                listenv_0.8.0                
[119] pbapply_1.5-0                 MASS_7.3-58.1                
[121] tidyselect_1.2.0              stringi_1.7.8                
[123] yaml_2.3.6                    GOSemSim_2.20.0              
[125] BiocSingular_1.10.0           locfit_1.5-9.6               
[127] grid_4.1.1                    fastmatch_1.1-3              
[129] tools_4.1.1                   timechange_0.1.1             
[131] future.apply_1.10.0           parallel_4.1.1               
[133] rstudioapi_0.14               gridExtra_2.3                
[135] farver_2.1.1                  Rtsne_0.16                   
[137] BiocManager_1.30.19           digest_0.6.30                
[139] shiny_1.7.1                   Rcpp_1.0.9                   
[141] car_3.1-0                     broom_1.0.1                  
[143] BiocVersion_3.14.0            later_1.3.0                  
[145] RcppAnnoy_0.0.20              httr_1.4.4                   
[147] AnnotationDbi_1.56.2          rvest_1.0.3                  
[149] XML_3.99-0.12                 fs_1.5.2                     
[151] tensor_1.5                    reticulate_1.26              
[153] splines_4.1.1                 uwot_0.1.14                  
[155] yulab.utils_0.0.5             tidytree_0.3.9               
[157] spatstat.utils_3.0-1          graphlayouts_0.8.3           
[159] sp_1.5-1                      ggplotify_0.1.0              
[161] plotly_4.10.1                 xtable_1.8-4                 
[163] jsonlite_1.8.3                tidygraph_1.2.2              
[165] ggfun_0.0.8                   R6_2.5.1                     
[167] pillar_1.8.1                  htmltools_0.5.3              
[169] mime_0.12                     glue_1.6.2                   
[171] fastmap_1.1.0                 BiocParallel_1.28.3          
[173] BiocNeighbors_1.12.0          interactiveDisplayBase_1.32.0
[175] codetools_0.2-18              fgsea_1.20.0                 
[177] utf8_1.2.2                    lattice_0.20-45              
[179] spatstat.sparse_3.0-0         curl_4.3.3                   
[181] leiden_0.4.3                  GO.db_3.14.0                 
[183] limma_3.50.3                  survival_3.4-0               
[185] munsell_0.5.0                 DO.db_2.9                    
[187] GenomeInfoDbData_1.2.7        haven_2.4.3                  
[189] reshape2_1.4.4                gtable_0.3.1           
```
