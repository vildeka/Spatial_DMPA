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

#########
# PATHS #
#########
input_dir <- "../../results/02_integrate_ref_data/"
result_dir <- "../../results/03_integrate_ref_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
seuratObj <- readRDS(paste0(input_dir,"seuratObj_harmony.RDS"))
```

``` r
my_theme <-
  list(
    #scale_fill_manual(values = friendly_cols),
    #scale_color_manual(values = friendly_cols),
    theme_bw() +
      #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        #legend.position = "bottom",
        #aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )
```

# Clustering

``` r
##################################
# EVALUATE CLUSTERING RESOLUTION #
##################################
seuratObj <- FindNeighbors(seuratObj, reduction = "harmony", dims = 1:20, k.param = 20, prune.SNN = 1/15) 

# Clustering with louvain (algorithm 1) or leiden (algorithm 4)
for (res in c(0.1, 0.25, 0.5, 1, 1.5, 2)) {
    seuratObj <- FindClusters(seuratObj, resolution = res, algorithm = 1)
}
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20539
    ## Number of edges: 696725
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9669
    ## Number of communities: 10
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20539
    ## Number of edges: 696725
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9353
    ## Number of communities: 12
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20539
    ## Number of edges: 696725
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8918
    ## Number of communities: 17
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20539
    ## Number of edges: 696725
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8433
    ## Number of communities: 22
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20539
    ## Number of edges: 696725
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8116
    ## Number of communities: 27
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20539
    ## Number of edges: 696725
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7863
    ## Number of communities: 31
    ## Elapsed time: 4 seconds

``` r
# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.

seuratObj <- SetIdent(seuratObj, value = "RNA_snn_res.1.5")
```

``` r
# obj <- seuratObj
# cluster <- sym("RNA_snn_res.0.1")
plot_clusters.fun <- function(obj, cluster, red = "umap_harmony", color = "Brew_all", lable = TRUE){
  if(color == "Brew_all"){
    pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2") )}
  
  cluster <- sym(cluster)
  
  if(lable == TRUE){ lab <- cluster
    text <- NoLegend() #+ labs(color= "Clusters")
  }else{lab <- sym(lable)
    text <- guides(color = "none")}
  
  feat <- obj %>%
    select(.cell, !!(cluster), !!(lab), sample_name, nCount_RNA, nFeature_RNA) %>%
    group_by(!!(cluster)) %>%
    add_tally() %>%
    arrange(nFeature_RNA) %>%
    arrange(desc(n))
  
  lable_df <- feat %>%
    ungroup() %>%
    group_by(!!(lab)) %>%
    select(!!(lab), contains(red)) %>% 
    summarize_all(mean)
  
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  p <- ggplot(feat, aes(!!(red_1), !!(red_2), 
                        color = !!cluster), label=TRUE) + 
    geom_point(alpha = 0.5, size=.5) + ggtitle(as_label(cluster)) +
    geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5) +
    #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
    scale_color_manual(values = pal)  +
    my_theme + text
  return(p)
}

res <- c("RNA_snn_res.1", "RNA_snn_res.1.5")
p <- map(res, ~plot_clusters.fun(seuratObj, cluster=.x))
plot_grid(ncol = 2, 
         plotlist = p)
```

<img src="./Figures/03a_plot_resolution.png" style="display: block; margin: auto;" />

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(seuratObj, paste0(result_dir,"seuratObj_clustered.RDS"))
# seuratObj<- readRDS(paste0(result_dir,"seuratObj_clustered.RDS"))
```

## Session info

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /Users/vilkal/Applications/miniconda3/envs/Spatial_DMPA/lib/libopenblasp-r0.3.18.dylib
    ## 
    ## locale:
    ## [1] sv_SE.UTF-8/sv_SE.UTF-8/sv_SE.UTF-8/C/sv_SE.UTF-8/sv_SE.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] harmony_0.1.0      Rcpp_1.0.8.3       cowplot_1.1.1      tidyseurat_0.5.1  
    ##  [5] ttservice_0.1.2    SeuratObject_4.0.4 Seurat_4.1.0       forcats_0.5.1     
    ##  [9] stringr_1.4.0      dplyr_1.0.8        purrr_0.3.4        readr_2.1.2       
    ## [13] tidyr_1.2.0        tibble_3.1.6       ggplot2_3.3.5      tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1          backports_1.4.1       plyr_1.8.7           
    ##   [4] igraph_1.3.0          lazyeval_0.2.2        splines_4.1.2        
    ##   [7] listenv_0.8.0         scattermore_0.8       digest_0.6.29        
    ##  [10] htmltools_0.5.2       fansi_1.0.3           magrittr_2.0.3       
    ##  [13] tensor_1.5            cluster_2.1.2         ROCR_1.0-11          
    ##  [16] tzdb_0.2.0            globals_0.14.0        modelr_0.1.8         
    ##  [19] matrixStats_0.61.0    spatstat.sparse_2.1-0 colorspace_2.0-3     
    ##  [22] rvest_1.0.2           ggrepel_0.9.1         haven_2.4.3          
    ##  [25] xfun_0.30             crayon_1.5.1          jsonlite_1.8.0       
    ##  [28] spatstat.data_2.1-4   survival_3.2-13       zoo_1.8-9            
    ##  [31] glue_1.6.2            polyclip_1.10-0       gtable_0.3.0         
    ##  [34] leiden_0.3.9          future.apply_1.8.1    abind_1.4-5          
    ##  [37] scales_1.1.1          DBI_1.1.2             spatstat.random_2.2-0
    ##  [40] miniUI_0.1.1.1        viridisLite_0.4.0     xtable_1.8-4         
    ##  [43] reticulate_1.24       spatstat.core_2.4-2   htmlwidgets_1.5.4    
    ##  [46] httr_1.4.2            RColorBrewer_1.1-3    ellipsis_0.3.2       
    ##  [49] ica_1.0-2             pkgconfig_2.0.3       farver_2.1.0         
    ##  [52] sass_0.4.1            uwot_0.1.11           dbplyr_2.1.1         
    ##  [55] deldir_1.0-6          utf8_1.2.2            tidyselect_1.1.2     
    ##  [58] labeling_0.4.2        rlang_1.0.2           reshape2_1.4.4       
    ##  [61] later_1.3.0           munsell_0.5.0         cellranger_1.1.0     
    ##  [64] tools_4.1.2           cli_3.2.0             generics_0.1.2       
    ##  [67] broom_0.7.12          ggridges_0.5.3        evaluate_0.15        
    ##  [70] fastmap_1.1.0         yaml_2.3.5            goftest_1.2-3        
    ##  [73] knitr_1.38            fs_1.5.2              fitdistrplus_1.1-8   
    ##  [76] RANN_2.6.1            pbapply_1.5-0         future_1.24.0        
    ##  [79] nlme_3.1-157          mime_0.12             xml2_1.3.3           
    ##  [82] compiler_4.1.2        rstudioapi_0.13       plotly_4.10.0        
    ##  [85] png_0.1-7             spatstat.utils_2.3-0  reprex_2.0.1         
    ##  [88] bslib_0.3.1           stringi_1.7.6         highr_0.9            
    ##  [91] lattice_0.20-45       Matrix_1.4-1          vctrs_0.4.0          
    ##  [94] pillar_1.7.0          lifecycle_1.0.1       spatstat.geom_2.4-0  
    ##  [97] lmtest_0.9-40         jquerylib_0.1.4       RcppAnnoy_0.0.19     
    ## [100] data.table_1.14.2     irlba_2.3.5           httpuv_1.6.5         
    ## [103] patchwork_1.1.1       R6_2.5.1              promises_1.2.0.1     
    ## [106] KernSmooth_2.23-20    gridExtra_2.3         parallelly_1.31.0    
    ## [109] codetools_0.2-18      MASS_7.3-56           assertthat_0.2.1     
    ## [112] withr_2.5.0           sctransform_0.3.3     mgcv_1.8-40          
    ## [115] parallel_4.1.2        hms_1.1.1             grid_4.1.2           
    ## [118] rpart_4.1.16          rmarkdown_2.11        Rtsne_0.15           
    ## [121] shiny_1.7.1           lubridate_1.8.0
