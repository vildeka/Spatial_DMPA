# Clustering filtered ST data

3/1/23

## Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
source("../bin/spatial_visualization.R")
source("../bin/plotting_functions.R")
#library(harmony)

#########
# PATHS #
#########
input_dir <- "../results/02_integrate_st_data/"
result_dir <- "../results/03_clustering_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
DATA <- readRDS(paste0(input_dir,"SeuratObj_harmony_filt.RDS"))
#DATA <- readRDS(paste0(input_dir,"seuratObj_harmony.RDS"))

#################
# COLOUR PALLET #
#################
clus <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2") )
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

## Clustering

``` r
##################################
# EVALUATE CLUSTERING RESOLUTION #
##################################
DATA <- FindNeighbors(DATA, reduction = "harmony", dims = 1:20, k.param = 20, prune.SNN = 1/15) 

# Clustering with louvain (algorithm 1) or leiden (algorithm 4)
for (res in c(0.1, 0.25, 0.5, 1, 1.5, 2)) {
    DATA <- FindClusters(DATA, resolution = res, algorithm = 1)
}
```

    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6508
    Number of edges: 244213

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9412
    Number of communities: 3
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6508
    Number of edges: 244213

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8947
    Number of communities: 4
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6508
    Number of edges: 244213

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8470
    Number of communities: 8
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6508
    Number of edges: 244213

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7844
    Number of communities: 13
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6508
    Number of edges: 244213

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7414
    Number of communities: 15
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6508
    Number of edges: 244213

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7044
    Number of communities: 19
    Elapsed time: 0 seconds

``` r
# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.
```

``` r
res <- c("RNA_snn_res.1", "RNA_snn_res.1.5")
p <- map(res, ~plot_clusters.fun(DATA, cluster=.x))
plot_grid(ncol = 2, 
         plotlist = p)
```

```{=html}
<img src="../Figures/03/03a_plot_resolution.png"
data-fig-align="center" />
```
``` r
DATA <- DATA %>%
  rename(Clusters="RNA_snn_res.1.5") %>%
  SetIdent(., value = "Clusters")
```

### Plot clusters on tissue section:

``` r
#is.character(pull(DATA, RNA_snn_res.1.5))

# Old single image version of function:
# plots <- DATA %>%
#   mutate(group = orig.ident) %>%
#   nest(., data = -group) %>%
#   pmap(., 
#     ~plot_spatial.fun(..2,
#       sampleid = ..1,
#       colors = clus,
#       geneid = "Clusters",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
#       zoom = "zoom",
#       img_alpha = 0,
#       point_size = 1)
#     )
# 
# 
# 
# legend <- get_legend(plots[[1]] + theme(legend.position="right"))
# combined <- wrap_plots(plots, ncol=2) & theme(legend.position="none")
# combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .2)) 
# combined

############################
# PLOT FACET WRAP CLUSTERS #
############################
plot_st_meta.fun( DATA,  # filter(spe, decon_columns[[..1]]=="1"), # removes spots with no % of given cell type
        feat = "Clusters",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
        zoom = "zoom",
        colors = clus,
        #annot_col = "#dbd9d9",
        annot_line = .1,
        img_alpha = 0,
        point_size = 0.5)
```

<img src="../Figures/03/03c_clust_plot.png" data-fig-align="center"/>

# Paulos Code

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
 saveRDS(DATA, paste0(result_dir,"seuratObj_clustered_filt.RDS"))
# saveRDS(DATA, paste0(result_dir,"seuratObj_clustered.RDS"))
# DATA<- readRDS(paste0(result_dir,"seuratObj_clustered.RDS"))
```

## Session info

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
     [1] patchwork_1.1.2    cowplot_1.1.1      tidyseurat_0.5.3   ttservice_0.2.2   
     [5] SeuratObject_4.1.3 Seurat_4.3.0       forcats_0.5.2      stringr_1.5.0     
     [9] dplyr_1.0.10       purrr_1.0.1        readr_2.1.3        tidyr_1.2.1       
    [13] tibble_3.1.8       ggplot2_3.4.0      tidyverse_1.3.2   

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
     [43] DBI_1.1.3              spatstat.random_3.0-1  miniUI_0.1.1.1        
     [46] Rcpp_1.0.9             viridisLite_0.4.1      xtable_1.8-4          
     [49] reticulate_1.27        htmlwidgets_1.6.1      httr_1.4.4            
     [52] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
     [55] pkgconfig_2.0.3        farver_2.1.1           uwot_0.1.14           
     [58] dbplyr_2.2.1           deldir_1.0-6           utf8_1.2.2            
     [61] tidyselect_1.2.0       labeling_0.4.2         rlang_1.0.6           
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
    [103] httpuv_1.6.8           R6_2.5.1               promises_1.2.0.1      
    [106] KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.33.0     
    [109] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
    [112] withr_2.5.0            sctransform_0.3.5      parallel_4.1.2        
    [115] hms_1.1.2              grid_4.1.2             rmarkdown_2.20        
    [118] googledrive_2.0.0      Rtsne_0.16             spatstat.explore_3.0-5
    [121] shiny_1.7.4            lubridate_1.9.0       

### Paulo section
