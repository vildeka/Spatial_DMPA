Integrate ST Data
================
3/13/23

## Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(harmony)
source("../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../results/01_QC_st_data/"
#input_dir <- "../results/00_load_st_data/"
result_dir <- "../results/02_integrate_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
DATA <- readRDS(paste0(input_dir,"seuratObj_filtered.RDS"))
```

## Identify Highly Variable Genes (HVG) across samples

``` r
################################
# SPLIT INTO SEPERATE DATASETS #
################################
DATA_nested <- DATA %>%
  mutate(batch = orig.ident) %>%
  nest(data = -batch) %>%
  mutate(data = imap(
    data, ~ .x %>%
      NormalizeData(., normalization.method = "LogNormalize", 
                    verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", 
                           nfeatures = 2000, 
                           verbose = FALSE) )) %>%
  mutate(data = setNames(.[["data"]], .$batch))

#########################################
# FIND HIGLY VARIABLE GENES PER DATASET #
#########################################
hvgs_heat <- DATA_nested %>%
  .$data %>%
  map(., ~ .x@assays$RNA@var.features) %>%
  ( function(x){unique(unlist(x)) ->> hvgs_all; return(x)} ) %>%
  # intersect across all samples:
  ( function(x){Reduce(intersect, x) ->> hvgs; return(x)} ) %>% 
  imap_dfc(., ~hvgs_all %in% .x, .id=.y) %>%
  mutate(rownames = hvgs_all) %>%
  column_to_rownames(var = "rownames")

# choose the hvg present in at least two samples:
hig_var <- rownames(hvgs_heat)[rowSums(hvgs_heat)>2]

# remove all VDJ-genes from list of HVG
remove <- str_subset(hig_var, "^IGH|^IGK|^IGL|^TRA|^TRB|^TRD|^TRG")
hig_var <- setdiff(hig_var, remove)
```

## Heatmap of HVG in all samples

``` r
pheatmap::pheatmap(t(hvgs_heat * 1), cluster_rows = F, color = c("grey90", "grey20"))
```

<img src=".../Figures/02/02a_HVG_heatmap.png" data-fig-align="center" />

## Integration

``` r
############
# HARMONY #
###########
DATA <- DATA %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                      nfeatures = 4000,
                      verbose = FALSE) %>%
  ScaleData(verbose = FALSE, features = hig_var ) %>%
  RunPCA(verbose = FALSE, npcs = 50) %>%
  RunUMAP(dims = 1:20, n.components = 2L, n.neighbors = 10) 

DATA <- DATA %>%
  RunHarmony(group.by.vars = "orig.ident", 
             reduction = "pca", 
             dims.use = 1:20, 
             assay.use = "RNA") #%>%

DATA <-   DATA %>%
  RunUMAP(dims = 1:10, 
          n.neighbors = 10,
          reduction = "harmony",
          reduction.name = "umap_harmony")
```

## Plot before and after integration

``` r
#  dev.new(height=6, width=6.6929133858, noRStudioGD = TRUE)
res <- c("PC", "harmony", "UMAP", "umap_harmony")
title <- c("PCA raw data", "PCA Harmony integrated", "UMAP raw data", "UMAP Harmony integrated")
p <- map2(res, title, 
          ~plot_clusters.fun(DATA, 
                             cluster="orig.ident", txt_size = 9,
                             red=.x, lable=FALSE, title=.y))
plot_grid(ncol = 2, 
         plotlist = p)
```

<img src=".../Figures/02/02b_Plot_dim_reduction.png"
data-fig-align="center" />

## Plot marker genes

``` r
#  dev.new(height=3, width=8, noRStudioGD = TRUE)
################################
# VISUALIZE EXPR. OF KEY GENES #
################################
# col <- c("grey90","grey80","grey60","navy","black")
col <- c("lightgray", "mistyrose", "red", "dark red", "black")
genes <- c("KRT1", "KRT15", "CDH1")
# genes <- c("CD8A", "SFRP2", "CD3E")
# genes <- c("CD8A", "MYOZ2", "CD3E", "EPCAM", "COL6A1", "CD4")

p <- map(genes, ~plot_genes.fun(DATA, .x, col = col, lable = FALSE))
plot_grid(ncol = 3, 
          plotlist = p)
```

<img src=".../Figures/02/02c_plot_marker_genes.png"
data-fig-align="center" />

# Paulos base R code

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
#saveRDS(DATA, paste0(result_dir,"SeuratObj_harmony.RDS"))
saveRDS(DATA, paste0(result_dir,"SeuratObj_harmony_filt.RDS"))
#DATA <- readRDS(paste0(result_dir,"SeuratObj_harmony.RDS"))
```

### Session info

``` r
sessionInfo()
```

    R version 4.1.2 (2021-11-01)
    Platform: x86_64-apple-darwin13.4.0 (64-bit)
    Running under: macOS Big Sur 10.16

    Matrix products: default
    BLAS/LAPACK: /Users/vilkal/Applications/miniconda3/envs/Spatial_DMPA/lib/libopenblasp-r0.3.21.dylib

    locale:
    [1] sv_SE.UTF-8/sv_SE.UTF-8/sv_SE.UTF-8/C/sv_SE.UTF-8/sv_SE.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] harmony_0.1.1      Rcpp_1.0.9         cowplot_1.1.1      tidyseurat_0.5.3  
     [5] ttservice_0.2.2    SeuratObject_4.1.3 Seurat_4.3.0       forcats_0.5.2     
     [9] stringr_1.5.0      dplyr_1.0.10       purrr_1.0.1        readr_2.1.3       
    [13] tidyr_1.2.1        tibble_3.1.8       ggplot2_3.4.0      tidyverse_1.3.2   

    loaded via a namespace (and not attached):
      [1] readxl_1.4.1           backports_1.4.1        plyr_1.8.8            
      [4] igraph_1.3.5           lazyeval_0.2.2         sp_1.5-1              
      [7] splines_4.1.2          listenv_0.9.0          scattermore_0.8       
     [10] digest_0.6.31          htmltools_0.5.4        fansi_1.0.3           
     [13] magrittr_2.0.3         tensor_1.5             googlesheets4_1.0.1   
     [16] cluster_2.1.4          ROCR_1.0-11            tzdb_0.3.0            
     [19] globals_0.16.2         modelr_0.1.10          matrixStats_0.63.0    
     [22] timechange_0.2.0       spatstat.sparse_3.0-0  colorspace_2.0-3      
     [25] rvest_1.0.3            ggrepel_0.9.2          haven_2.5.1           
     [28] xfun_0.36              crayon_1.5.2           jsonlite_1.8.4        
     [31] progressr_0.13.0       spatstat.data_3.0-0    survival_3.5-0        
     [34] zoo_1.8-11             glue_1.6.2             polyclip_1.10-4       
     [37] gtable_0.3.1           gargle_1.2.1           leiden_0.4.3          
     [40] future.apply_1.10.0    abind_1.4-5            scales_1.2.1          
     [43] pheatmap_1.0.12        DBI_1.1.3              spatstat.random_3.0-1 
     [46] miniUI_0.1.1.1         viridisLite_0.4.1      xtable_1.8-4          
     [49] reticulate_1.27        htmlwidgets_1.6.1      httr_1.4.4            
     [52] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
     [55] farver_2.1.1           pkgconfig_2.0.3        uwot_0.1.14           
     [58] dbplyr_2.2.1           deldir_1.0-6           utf8_1.2.2            
     [61] labeling_0.4.2         tidyselect_1.2.0       rlang_1.0.6           
     [64] reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
     [67] cellranger_1.1.0       tools_4.1.2            cli_3.6.0             
     [70] generics_0.1.3         broom_1.0.2            ggridges_0.5.4        
     [73] evaluate_0.19          fastmap_1.1.0          yaml_2.3.6            
     [76] goftest_1.2-3          knitr_1.41             fs_1.5.2              
     [79] fitdistrplus_1.1-8     RANN_2.6.1             pbapply_1.6-0         
     [82] future_1.30.0          nlme_3.1-161           mime_0.12             
     [85] xml2_1.3.3             compiler_4.1.2         rstudioapi_0.14       
     [88] plotly_4.10.1          png_0.1-8              spatstat.utils_3.0-1  
     [91] reprex_2.0.2           stringi_1.7.12         lattice_0.20-45       
     [94] Matrix_1.5-3           vctrs_0.5.1            pillar_1.8.1          
     [97] lifecycle_1.0.3        spatstat.geom_3.0-3    lmtest_0.9-40         
    [100] RcppAnnoy_0.0.20       data.table_1.14.6      irlba_2.3.5.1         
    [103] httpuv_1.6.8           patchwork_1.1.2        R6_2.5.1              
    [106] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
    [109] parallelly_1.33.0      codetools_0.2-18       MASS_7.3-58.1         
    [112] assertthat_0.2.1       withr_2.5.0            sctransform_0.3.5     
    [115] parallel_4.1.2         hms_1.1.2              grid_4.1.2            
    [118] rmarkdown_2.20         googledrive_2.0.0      Rtsne_0.16            
    [121] spatstat.explore_3.0-5 shiny_1.7.4            lubridate_1.9.0       
