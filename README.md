# Introduction



This repository contains all scripts, custom functions, source files, and data used to generate the plots for the manuscript: **"RNA tomography reveals spatial gene expression maps of *Arabidopsis thaliana* roots infected with *Heterodera schachtii*"**



## R environment

```R
> sessionInfo()

R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Amsterdam
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] here_1.0.1            future_1.40.0         tibble_3.2.1          openxlsx_4.2.8        dplyr_1.1.4           writexl_1.5.4         RColorBrewer_1.1-3   
 [8] gridtext_0.1.5        tidyr_1.3.1           ComplexHeatmap_2.18.0 glmGamPoi_1.14.3      Seurat_5.2.1          SeuratObject_5.1.0    sp_2.2-0             
[15] Matrix_1.6-5          stringr_1.5.1         readxl_1.4.5          pheatmap_1.0.12       heatmaply_1.5.0       viridis_0.6.5         viridisLite_0.4.2    
[22] plotly_4.10.4         ggplot2_3.5.2         tomoda_1.12.0        

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.3.2               later_1.4.2                 bitops_1.0-9                cellranger_1.1.0           
  [6] polyclip_1.10-7             fastDummies_1.7.5           lifecycle_1.0.4             rprojroot_2.0.4             doParallel_1.0.17          
 [11] globals_0.17.0              lattice_0.22-7              MASS_7.3-60.0.1             dendextend_1.19.0           magrittr_2.0.3             
 [16] rmarkdown_2.29              yaml_2.3.10                 httpuv_1.6.16               sctransform_0.4.1           zip_2.3.2                  
 [21] spam_2.11-1                 askpass_1.2.1               spatstat.sparse_3.1-0       reticulate_1.42.0           cowplot_1.1.3.9000         
 [26] pbapply_1.7-2               abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1       
 [31] purrr_1.0.4                 BiocGenerics_0.48.1         RCurl_1.98-1.17             circlize_0.4.16             seriation_1.5.7            
 [36] GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6               irlba_2.3.5.1              
 [41] listenv_0.9.1               spatstat.utils_3.1-3        umap_0.2.10.0               goftest_1.2-3               RSpectra_0.16-2            
 [46] spatstat.random_3.3-3       fitdistrplus_1.2-2          parallelly_1.43.0           DelayedMatrixStats_1.24.0   codetools_0.2-20           
 [51] DelayedArray_0.28.0         xml2_1.3.8                  tidyselect_1.2.1            shape_1.4.6.1               farver_2.1.2               
 [56] TSP_1.2-4                   matrixStats_1.5.0           stats4_4.3.2                spatstat.explore_3.4-2      webshot_0.5.5              
 [61] jsonlite_2.0.0              GetoptLong_1.0.5            progressr_0.15.1            ggridges_0.5.6              survival_3.8-3             
 [66] iterators_1.0.14            foreach_1.5.2               tools_4.3.2                 ica_1.0-3                   Rcpp_1.0.14                
 [71] glue_1.8.0                  gridExtra_2.3               SparseArray_1.2.4           xfun_0.52                   MatrixGenerics_1.14.0      
 [76] GenomeInfoDb_1.38.8         ca_0.71.1                   withr_3.0.2                 fastmap_1.2.0               openssl_2.3.2              
 [81] digest_0.6.37               R6_2.6.1                    mime_0.13                   colorspace_2.1-1            Cairo_1.6-2                
 [86] scattermore_1.2             tensor_1.5                  spatstat.data_3.1-6         utf8_1.2.4                  generics_0.1.3             
 [91] data.table_1.17.0           httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.2.1              uwot_0.2.3                 
 [96] pkgconfig_2.0.3             gtable_0.3.6                registry_0.5-1              lmtest_0.9-40               XVector_0.42.0             
[101] htmltools_0.5.8.1           dotCall64_1.2               clue_0.3-66                 scales_1.4.0                Biobase_2.62.0             
[106] png_0.1-8                   spatstat.univar_3.1-2       knitr_1.50                  rstudioapi_0.17.1           reshape2_1.4.4             
[111] rjson_0.2.23                nlme_3.1-168                zoo_1.8-14                  GlobalOptions_0.1.2         KernSmooth_2.23-26         
[116] parallel_4.3.2              miniUI_0.1.2                pillar_1.10.2               vctrs_0.6.5                 RANN_2.6.2                 
[121] promises_1.3.2              xtable_1.8-4                cluster_2.1.8.1             evaluate_1.0.3              magick_2.8.6               
[126] cli_3.6.5                   compiler_4.3.2              rlang_1.1.6                 crayon_1.5.3                future.apply_1.11.3        
[131] labeling_0.4.3              plyr_1.8.9                  stringi_1.8.7               deldir_2.0-4                assertthat_0.2.1           
[136] lazyeval_0.2.2              spatstat.geom_3.3-6         RcppHNSW_0.6.0              patchwork_1.3.0             sparseMatrixStats_1.14.0   
[141] shiny_1.10.0                SummarizedExperiment_1.32.0 ROCR_1.0-11                 igraph_2.1.4         
```



## Citation

If you use or adapt code/data from this repository, please cite the manuscript above.