# Analysis code overview
To complete a reproducible run, clone this repositiory, then follow the runScript file in R or Rstuido. The analysis for this paper was completed on the UC Boulder Summit supercomputer (cite), but the code should be able to be run on a computer with 64 Gb RAM.


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
 [1] scuttle_1.4.0               scRNAseq_2.8.0             
 [3] ggpubr_0.4.0                slingshot_2.2.1            
 [5] TrajectoryUtils_1.2.0       SingleCellExperiment_1.16.0
 [7] princurve_2.1.6             clusterProfiler_4.2.2      
 [9] msigdbr_7.5.1               ggsankey_0.0.99999         
[11] lemon_0.4.5                 reshape_0.8.9              
[13] viridis_0.6.2               viridisLite_0.4.1          
[15] SingleR_1.8.1               SeuratDisk_0.0.0.9019      
[17] RColorBrewer_1.1-3          pheatmap_1.0.12            
[19] DESeq2_1.34.0               SummarizedExperiment_1.24.0
[21] Biobase_2.54.0              MatrixGenerics_1.6.0       
[23] matrixStats_0.62.0          GenomicRanges_1.46.1       
[25] GenomeInfoDb_1.30.1         IRanges_2.28.0             
[27] S4Vectors_0.32.2            BiocGenerics_0.40.0        
[29] colorspace_2.0-3            ggrepel_0.9.1              
[31] cowplot_1.1.1               scales_1.2.1               
[33] patchwork_1.1.2             DoubletFinder_2.0.3        
[35] clustree_0.4.4              ggraph_2.0.5               
[37] forcats_0.5.2               stringr_1.4.1              
[39] dplyr_1.0.10                purrr_0.3.5                
[41] readr_2.1.2                 tidyr_1.2.1                
[43] tibble_3.1.8                ggplot2_3.3.6              
[45] tidyverse_1.3.1             SeuratObject_4.1.3         
[47] Seurat_4.0.6               

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                rtracklayer_1.54.0           
  [3] scattermore_0.8               bit64_4.0.5                  
  [5] knitr_1.40                    irlba_2.3.5                  
  [7] DelayedArray_0.20.0           data.table_1.14.4            
  [9] rpart_4.1.19                  AnnotationFilter_1.18.0      
 [11] KEGGREST_1.34.0               RCurl_1.98-1.6               
 [13] generics_0.1.3                GenomicFeatures_1.46.5       
 [15] ScaledMatrix_1.2.0            RSQLite_2.2.18               
 [17] shadowtext_0.1.1              RANN_2.6.1                   
 [19] future_1.29.0                 bit_4.0.4                    
 [21] tzdb_0.3.0                    enrichplot_1.14.2            
 [23] spatstat.data_3.0-0           xml2_1.3.3                   
 [25] lubridate_1.9.0               httpuv_1.6.5                 
 [27] assertthat_0.2.1              xfun_0.34                    
 [29] hms_1.1.2                     babelgene_22.9               
 [31] promises_1.2.0.1              restfulr_0.0.15              
 [33] progress_1.2.2                fansi_1.0.3                  
 [35] dbplyr_2.1.1                  readxl_1.4.1                 
 [37] igraph_1.2.11                 DBI_1.1.3                    
 [39] geneplotter_1.72.0            htmlwidgets_1.5.4            
 [41] spatstat.geom_3.0-5           ellipsis_0.3.2               
 [43] backports_1.4.1               annotate_1.72.0              
 [45] biomaRt_2.50.3                deldir_1.0-6                 
 [47] sparseMatrixStats_1.6.0       vctrs_0.5.1                  
 [49] ensembldb_2.18.3              ROCR_1.0-11                  
 [51] abind_1.4-5                   cachem_1.0.6                 
 [53] withr_2.5.0                   ggforce_0.4.1                
 [55] progressr_0.11.0              sctransform_0.3.3            
 [57] GenomicAlignments_1.30.0      treeio_1.18.1                
 [59] prettyunits_1.1.1             goftest_1.2-3                
 [61] cluster_2.1.4                 DOSE_3.20.1                  
 [63] ExperimentHub_2.2.1           ape_5.6-2                    
 [65] lazyeval_0.2.2                crayon_1.5.2                 
 [67] genefilter_1.76.0             hdf5r_1.3.7                  
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
[121] mgcv_1.8-41                   tidyselect_1.2.0             
[123] stringi_1.7.8                 yaml_2.3.6                   
[125] GOSemSim_2.20.0               BiocSingular_1.10.0          
[127] locfit_1.5-9.6                grid_4.1.1                   
[129] fastmatch_1.1-3               tools_4.1.1                  
[131] timechange_0.1.1              future.apply_1.10.0          
[133] parallel_4.1.1                rstudioapi_0.14              
[135] gridExtra_2.3                 farver_2.1.1                 
[137] Rtsne_0.16                    BiocManager_1.30.19          
[139] digest_0.6.30                 shiny_1.7.1                  
[141] Rcpp_1.0.9                    car_3.1-0                    
[143] broom_1.0.1                   BiocVersion_3.14.0           
[145] later_1.3.0                   RcppAnnoy_0.0.20             
[147] httr_1.4.4                    AnnotationDbi_1.56.2         
[149] rvest_1.0.3                   XML_3.99-0.12                
[151] fs_1.5.2                      tensor_1.5                   
[153] reticulate_1.26               splines_4.1.1                
[155] yulab.utils_0.0.5             uwot_0.1.14                  
[157] tidytree_0.3.9                spatstat.utils_3.0-1         
[159] graphlayouts_0.8.3            sp_1.5-1                     
[161] ggplotify_0.1.0               plotly_4.10.1                
[163] xtable_1.8-4                  ggtree_3.2.1                 
[165] jsonlite_1.8.3                tidygraph_1.2.2              
[167] ggfun_0.0.8                   R6_2.5.1                     
[169] pillar_1.8.1                  htmltools_0.5.3              
[171] mime_0.12                     glue_1.6.2                   
[173] fastmap_1.1.0                 BiocParallel_1.28.3          
[175] BiocNeighbors_1.12.0          interactiveDisplayBase_1.32.0
[177] codetools_0.2-18              fgsea_1.20.0                 
[179] utf8_1.2.2                    lattice_0.20-45              
[181] spatstat.sparse_3.0-0         curl_4.3.3                   
[183] leiden_0.4.3                  GO.db_3.14.0                 
[185] survival_3.4-0                munsell_0.5.0                
[187] DO.db_2.9                     GenomeInfoDbData_1.2.7       
[189] haven_2.4.3                   reshape2_1.4.4               
[191] gtable_0.3.1                  spatstat.core_2.4-0          
```
