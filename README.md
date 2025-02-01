# Canine tumor-infiltrating immune cells versus circulating immune cells

This GitHub repository contains all the analysis code used in, "Characterization of canine tumor-infiltrating leukocyte transcriptomic signatures reveals conserved expression patterns with human osteosarcoma."

If you use our raw/processed data, extract data using the UCSC Cell Browser portal, or use portions of our code in your analysis, please cite:
> accepted, publication pending

If you have any questions or concerns, please submit an issue, contact the corresponding author(s), and/or contact Dylan Ammons at dylan.ammons @ colostate dot edu.

## File structure:
- [:file\_folder: scripts](/scripts) contains the analysis code, source file, and metadata used to complete analysis

## Browse the annotated dataset
The proccessed dataset is avaliable for browsing via the UCSC Cell Browser portal.
Using the portal you can explore feature expression throughout the dataset as well as obtain the transcriptomic signatures of each cell type though an interactive webpage.

Link to the dataset: https://canine-tils-blood.cells.ucsc.edu

Link to UCSC Cell Browser documentation: https://cellbrowser.readthedocs.io/en/master/

## Software versions
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
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] ComplexHeatmap_2.13.2       ggtree_3.2.1               
 [3] ape_5.7-1                   scuttle_1.4.0              
 [5] scRNAseq_2.8.0              ggpubr_0.4.0               
 [7] slingshot_2.7.0             TrajectoryUtils_1.2.0      
 [9] SingleCellExperiment_1.16.0 princurve_2.1.6            
[11] clusterProfiler_4.2.2       msigdbr_7.5.1              
[13] ggsankey_0.0.99999          lemon_0.4.5                
[15] reshape_0.8.9               viridis_0.6.2              
[17] viridisLite_0.4.1           SingleR_1.8.1              
[19] SeuratDisk_0.0.0.9019       RColorBrewer_1.1-3         
[21] pheatmap_1.0.12             DESeq2_1.34.0              
[23] SummarizedExperiment_1.24.0 Biobase_2.54.0             
[25] MatrixGenerics_1.6.0        matrixStats_1.0.0          
[27] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
[29] IRanges_2.28.0              S4Vectors_0.32.4           
[31] BiocGenerics_0.40.0         colorspace_2.0-3           
[33] ggrepel_0.9.1               cowplot_1.1.1              
[35] scales_1.3.0                patchwork_1.3.0.9000       
[37] DoubletFinder_2.0.3         clustree_0.4.4             
[39] ggraph_2.0.5                forcats_0.5.2              
[41] stringr_1.4.1               dplyr_1.0.10               
[43] purrr_0.3.5                 readr_2.1.2                
[45] tidyr_1.2.1                 tibble_3.1.8               
[47] ggplot2_3.3.6               tidyverse_1.3.1            
[49] SeuratObject_4.1.3          Seurat_4.3.0               

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                    ica_1.0-3                    
  [3] Rsamtools_2.10.0              foreach_1.5.2                
  [5] lmtest_0.9-40                 crayon_1.5.2                 
  [7] MASS_7.3-58.1                 nlme_3.1-160                 
  [9] backports_1.4.1               reprex_2.0.1                 
 [11] GOSemSim_2.20.0               rlang_1.0.6                  
 [13] XVector_0.34.0                ROCR_1.0-11                  
 [15] readxl_1.4.1                  irlba_2.3.5                  
 [17] filelock_1.0.2                BiocParallel_1.37.1          
 [19] rjson_0.2.21                  bit64_4.0.5                  
 [21] glue_1.6.2                    sctransform_0.3.5            
 [23] parallel_4.1.1                spatstat.sparse_3.0-0        
 [25] AnnotationDbi_1.56.2          DOSE_3.20.1                  
 [27] spatstat.geom_3.0-5           haven_2.4.3                  
 [29] tidyselect_1.2.0              fitdistrplus_1.1-8           
 [31] XML_3.99-0.12                 zoo_1.8-9                    
 [33] GenomicAlignments_1.30.0      xtable_1.8-4                 
 [35] magrittr_2.0.3                cli_3.6.0                    
 [37] zlibbioc_1.40.0               rstudioapi_0.14              
 [39] miniUI_0.1.1.1                sp_1.5-1                     
 [41] fastmatch_1.1-3               ensembldb_2.18.3             
 [43] treeio_1.18.1                 shiny_1.7.1                  
 [45] BiocSingular_1.10.0           xfun_0.39                    
 [47] clue_0.3-62                   cluster_2.1.4                
 [49] tidygraph_1.2.2               KEGGREST_1.34.0              
 [51] interactiveDisplayBase_1.32.0 listenv_0.8.0                
 [53] Biostrings_2.62.0             png_0.1-7                    
 [55] future_1.29.0                 withr_2.5.0                  
 [57] bitops_1.0-7                  ggforce_0.4.1                
 [59] plyr_1.8.7                    cellranger_1.1.0             
 [61] AnnotationFilter_1.18.0       pillar_1.8.1                 
 [63] GlobalOptions_0.1.2           cachem_1.0.6                 
 [65] GenomicFeatures_1.46.5        fs_1.5.2                     
 [67] hdf5r_1.3.8                   GetoptLong_1.0.5             
 [69] DelayedMatrixStats_1.16.0     vctrs_0.5.1                  
 [71] ellipsis_0.3.2                generics_0.1.3               
 [73] tools_4.1.1                   munsell_0.5.0                
 [75] tweenr_2.0.2                  fgsea_1.20.0                 
 [77] DelayedArray_0.20.0           fastmap_1.1.0                
 [79] compiler_4.1.1                abind_1.4-5                  
 [81] httpuv_1.6.5                  rtracklayer_1.54.0           
 [83] ExperimentHub_2.2.1           plotly_4.10.1                
 [85] GenomeInfoDbData_1.2.7        gridExtra_2.3                
 [87] lattice_0.20-45               deldir_1.0-6                 
 [89] utf8_1.2.2                    later_1.3.0                  
 [91] BiocFileCache_2.2.1           jsonlite_1.8.3               
 [93] ScaledMatrix_1.2.0            tidytree_0.3.9               
 [95] pbapply_1.5-0                 carData_3.0-5                
 [97] sparseMatrixStats_1.6.0       genefilter_1.76.0            
 [99] lazyeval_0.2.2                promises_1.2.0.1             
[101] car_3.1-0                     doParallel_1.0.17            
[103] R.utils_2.12.1                goftest_1.2-3                
[105] spatstat.utils_3.0-1          reticulate_1.34.0            
[107] Rtsne_0.16                    downloader_0.4               
[109] uwot_0.1.14                   igraph_1.5.1                 
[111] survival_3.4-0                yaml_2.3.6                   
[113] htmltools_0.5.3               memoise_2.0.1                
[115] BiocIO_1.4.0                  locfit_1.5-9.6               
[117] graphlayouts_0.8.3            digest_0.6.30                
[119] assertthat_0.2.1              mime_0.12                    
[121] rappdirs_0.3.3                RSQLite_2.2.18               
[123] yulab.utils_0.0.5             future.apply_1.10.0          
[125] data.table_1.14.4             blob_1.2.3                   
[127] R.oo_1.25.0                   splines_4.1.1                
[129] AnnotationHub_3.2.2           ProtGenerics_1.26.0          
[131] RCurl_1.98-1.12               broom_1.0.1                  
[133] hms_1.1.2                     modelr_0.1.8                 
[135] BiocManager_1.30.19           shape_1.4.6                  
[137] aplot_0.1.2                   Rcpp_1.0.11                  
[139] RANN_2.6.1                    circlize_0.4.16              
[141] enrichplot_1.14.2             fansi_1.0.3                  
[143] tzdb_0.3.0                    parallelly_1.32.1            
[145] R6_2.5.1                      ggridges_0.5.4               
[147] lifecycle_1.0.3               curl_4.3.3                   
[149] ggsignif_0.6.4                leiden_0.4.3                 
[151] DO.db_2.9                     Matrix_1.5-1                 
[153] qvalue_2.26.0                 RcppAnnoy_0.0.20             
[155] iterators_1.0.14              spatstat.explore_3.0-5       
[157] htmlwidgets_1.5.4             beachmat_2.10.0              
[159] polyclip_1.10-4               biomaRt_2.50.3               
[161] shadowtext_0.1.1              timechange_0.1.1             
[163] gridGraphics_0.5-1            rvest_1.0.3                  
[165] globals_0.16.1                spatstat.random_3.1-3        
[167] progressr_0.11.0              codetools_0.2-18             
[169] lubridate_1.9.0               GO.db_3.14.0                 
[171] prettyunits_1.1.1             dbplyr_2.1.1                 
[173] R.methodsS3_1.8.2             gtable_0.3.1                 
[175] DBI_1.1.3                     ggfun_0.0.8                  
[177] tensor_1.5                    httr_1.4.4                   
[179] KernSmooth_2.23-20            stringi_1.7.8                
[181] progress_1.2.2                reshape2_1.4.4               
[183] farver_2.1.1                  annotate_1.72.0              
[185] xml2_1.3.3                    BiocNeighbors_1.12.0         
[187] restfulr_0.0.15               geneplotter_1.72.0           
[189] ggplotify_0.1.0               scattermore_0.8              
[191] BiocVersion_3.14.0            bit_4.0.4                    
[193] scatterpie_0.1.7              spatstat.data_3.0-0          
[195] pkgconfig_2.0.3               babelgene_22.9               
[197] rstatix_0.7.0                 knitr_1.40    
```
