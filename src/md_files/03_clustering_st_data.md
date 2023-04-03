Clustering filtered spatial data
================
3/21/23

### Load data and libraries

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
DATA <- readRDS(paste0(input_dir,"seuratObj_harmony_filt.RDS"))
#DATA <- readRDS(paste0(input_dir,"seuratObj_harmony.RDS"))

#################
# COLOUR PALLET #
#################
clus_1.5 <- c("#F8766D", "#7CAE00", "#CD9600", "#00A9FF", "#984EA3","#C77CFF", "#00BFC4", "#FF7F00",
              "#FFFF33", "#E41A1C","#377EB8","#FF61CC", "#4DAF4A",  "#A65628", "#F781BF", "#999999")
clus_1 <- c( "#F8766D","#7CAE00", "#CD9600", "#00A9FF", "#377EB8","#984EA3", "#E41A1C", "#C77CFF",
             "#00BFC4", "#FF7F00","#FFFF33", "#FF61CC", "#4DAF4A",  "#A65628", "#F781BF", "#999999")

# clus <- c(scales::hue_pal()(8),
#              RColorBrewer::brewer.pal(9,"Set1"),
#              RColorBrewer::brewer.pal(8,"Set2"),
#              RColorBrewer::brewer.pal(8,"Accent"),
#              RColorBrewer::brewer.pal(9,"Pastel1"),
#              RColorBrewer::brewer.pal(8,"Pastel2") )
# scales::show_col(clus_1.5)
# scales::show_col(RColorBrewer::brewer.pal(9,"Set1"))
```

## Clustering

``` r
##################################
# EVALUATE CLUSTERING RESOLUTION #
##################################
DATA <- FindNeighbors(DATA, reduction = "harmony", dims = 1:30, k.param = 15, prune.SNN = 1/15) 

# Clustering with louvain (algorithm 1) or leiden (algorithm 4)
for (res in c(0.1, 0.5, 1, 1.5, 2)) {
    DATA <- FindClusters(DATA, resolution = res, algorithm = 1)
}
```

    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6598
    Number of edges: 299026

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9375
    Number of communities: 2
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6598
    Number of edges: 299026

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8316
    Number of communities: 7
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6598
    Number of edges: 299026

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7608
    Number of communities: 11
    Elapsed time: 0 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6598
    Number of edges: 299026

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7115
    Number of communities: 13
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

    Number of nodes: 6598
    Number of edges: 299026

    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.6669
    Number of communities: 16
    Elapsed time: 0 seconds

``` r
# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.
```

### UMAP of cluster resolutions

``` r
res <- c("RNA_snn_res.1", "RNA_snn_res.1.5")
p <- map2(res, list(clus_1,clus_1.5),
          ~plot_clusters.fun(DATA, cluster=.x, txt_size = 7, 
                             color = .y, red = "umap_harmony"))
plot_grid(ncol = 2, 
          plotlist = p)
```

<img src="../Figures/03/03a_plot_resolution.png"
data-fig-align="center" />

### Cluster resolutions on tissue

``` r
# dev.new(width=6.6929133858, height=14, noRStudioGD = TRUE)
plots <- DATA %>%
  mutate(group = orig.ident) %>%
  nest(., data = -group) %>%
  mutate( "res_1" = pmap(., 
    ~plot_spatial.fun(..2, sampleid=..1, geneid="RNA_snn_res.1", 
                      point_size = 0.8, zoom="zoom", colors = clus_1))) %>%
  mutate( "res_2" = pmap(., 
    ~plot_spatial.fun(..2, sampleid=..1, geneid="RNA_snn_res.1.5",
                      point_size = 0.8, zoom="zoom", colors = clus_1.5)))

legend_1 <- get_legend(plots$res_1[[2]] + theme(legend.position="right"))
legend_2 <- get_legend(plots$res_2[[1]] + theme(legend.position="right"))
legend <- plot_grid( legend_1, legend_2, ncol = 1)
combined <- wrap_plots(plotlist=c(plots$res_1, plots$res_2), nrow = 8, byrow = F) & theme(legend.position="none")
combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .3)) 
combined
```

<img src="../Figures/03/03b_plot_resolutions_on_tissue.png"
data-fig-align="center" />

### Set cluster resolution

``` r
DATA <- DATA %>%
  rename(Clusters="RNA_snn_res.1") %>%
  SetIdent(., value = "Clusters") %>%
  select(-any_of(contains("RNA_snn_res")))
```

### Spot distribution by clusters

``` r
DATA_sub <- DATA %>%
  mutate(gr = .$groups) %>%
  mutate(ID = .$orig.ident) %>%
  nest(., data=-c(gr, orig.ident)) %>%
  mutate(epi =  map(data, ~filter(.x, !(sp_annot == "SubMuc"))),
         subMuc =  map(data, ~filter(.x, sp_annot == "SubMuc"))) %>%
  mutate(across(c("epi", "subMuc"), ~map(., ~table(.x$Clusters)), .names = "{.col}_n_before")) %>%
  mutate(across(contains("_n_"), ~set_names(.x, paste0(.data[["gr"]],"_",.data[["orig.ident"]]))))
      
table(DATA$Clusters)

bind_cols(DATA_sub$epi_n_before, "Clus" = paste0("**",names(table(DATA$Clusters)),"**")) %>%
  rowwise() %>% 
  mutate(DMPA_sum = sum(c_across(starts_with("DMPA_"))),
         ctrl_sum = sum(c_across(starts_with("ctrl_")))) %>%
  select(sort(colnames(.)[1:8]), everything()) %>%
  knitr::kable(., caption = "Distribution of epithelial spots per cluster per subject")

bind_cols(DATA_sub$subMuc_n_before, "Clus" = paste0("**",names(table(DATA$Clusters)),"**")) %>%
  rowwise() %>% 
  mutate(DMPA_sum = sum(c_across(starts_with("DMPA_"))),
         ctrl_sum = sum(c_across(starts_with("ctrl_")))) %>%
  select(sort(colnames(.)[1:8]), everything()) %>%
  knitr::kable(., caption = "Distribution of submucosal spots per cluster per subject")
```


       0    1    2    3    4    5    6    7    8    9   10 
    1148  884  747  656  646  558  511  507  471  432   38 

| ctrl_P031 | ctrl_P080 | ctrl_P105 | ctrl_P118 | DMPA_P097 | DMPA_P107 | DMPA_P108 | DMPA_P114 | Clus   | DMPA_sum | ctrl_sum |
|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|:-------|---------:|---------:|
|         0 |         0 |         1 |         0 |         0 |         1 |         0 |         0 | **0**  |        1 |        1 |
|         1 |         0 |         0 |         0 |         0 |         1 |         0 |         1 | **1**  |        2 |        1 |
|         0 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | **2**  |        0 |        0 |
|         3 |         0 |         3 |         0 |         0 |         0 |         0 |         1 | **3**  |        1 |        6 |
|         0 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | **4**  |        0 |        0 |
|        79 |        67 |        47 |        77 |        67 |        31 |        37 |        67 | **5**  |      202 |      270 |
|        60 |        66 |        43 |        69 |        30 |        76 |        29 |       112 | **6**  |      247 |      238 |
|        89 |        67 |        34 |        77 |        49 |        51 |        38 |        61 | **7**  |      199 |      267 |
|        12 |         9 |         1 |         6 |         4 |         7 |         3 |         1 | **8**  |       15 |       28 |
|        92 |        69 |        33 |        53 |        51 |        51 |        26 |        56 | **9**  |      184 |      247 |
|         0 |         0 |         1 |         0 |         0 |         0 |         0 |         1 | **10** |        1 |        1 |

Distribution of epithelial spots per cluster per subject

| ctrl_P031 | ctrl_P080 | ctrl_P105 | ctrl_P118 | DMPA_P097 | DMPA_P107 | DMPA_P108 | DMPA_P114 | Clus   | DMPA_sum | ctrl_sum |
|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|:-------|---------:|---------:|
|       150 |       208 |       124 |       201 |       227 |        20 |        34 |       182 | **0**  |      463 |      683 |
|       131 |       192 |       118 |       121 |       182 |        10 |        36 |        91 | **1**  |      319 |      562 |
|        37 |       155 |        72 |        91 |       231 |        35 |         9 |       117 | **2**  |      392 |      355 |
|        46 |       100 |        76 |       130 |       141 |        25 |        43 |        88 | **3**  |      297 |      352 |
|        24 |        73 |        76 |       192 |       123 |        37 |        21 |       100 | **4**  |      281 |      365 |
|         1 |         5 |         8 |        20 |        41 |         1 |         1 |         9 | **5**  |       52 |       34 |
|         0 |         3 |         0 |         0 |         4 |        19 |         0 |         0 | **6**  |       23 |        3 |
|         1 |        33 |         0 |         2 |         4 |         0 |         0 |         1 | **7**  |        5 |       36 |
|        24 |        46 |        59 |        85 |       102 |        14 |        20 |        78 | **8**  |      214 |      214 |
|         0 |         0 |         0 |         0 |         1 |         0 |         0 |         0 | **9**  |        1 |        0 |
|         1 |         1 |         0 |         0 |         5 |        29 |         0 |         0 | **10** |       34 |        2 |

Distribution of submucosal spots per cluster per subject

### Plot final clusters on tissue section:

``` r
# Horizontal 
# dev.new(width=7, height=3.5, noRStudioGD = TRUE)
############################
# PLOT FACET WRAP CLUSTERS #
############################
(p <- plot_st_meta.fun( DATA,
        feat =  "Clusters",
        zoom = "zoom",
        colors = clus_1,
        alpha = .9,
        #annot_col = "#dbd9d9",
        annot_line = .1,
        img_alpha = 0,
        point_size = .35))
```

<img src="../Figures/03/03c_clust_plot.png" data-fig-align="center" />

``` r
###################
# ADD ANNOTATION #
##################
ord1 <- c("Sup_1","Sup_2","Basal_2","Basal_1","0","1","2","3","4","8","10")
ord2 <- c("6","9","7","5","0","1","2","3","4","8","10")
epi_layers <- set_names(ord1, ord2)

DATA <- DATA %>%
mutate(layers = factor(epi_layers[as.character(.$Clusters)], levels = ord1))
```

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(DATA, paste0(result_dir,"seuratObj_clustered_filt.RDS"))
# saveRDS(DATA, paste0(result_dir,"seuratObj_clustered.RDS"))
# DATA<- readRDS(paste0(result_dir,"seuratObj_clustered.RDS"))
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
     [91] reprex_2.0.2           stringi_1.7.12         highr_0.10            
     [94] lattice_0.20-45        Matrix_1.5-3           vctrs_0.5.1           
     [97] pillar_1.8.1           lifecycle_1.0.3        spatstat.geom_3.0-3   
    [100] lmtest_0.9-40          RcppAnnoy_0.0.20       data.table_1.14.6     
    [103] irlba_2.3.5.1          httpuv_1.6.8           R6_2.5.1              
    [106] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
    [109] parallelly_1.33.0      codetools_0.2-18       MASS_7.3-58.1         
    [112] assertthat_0.2.1       withr_2.5.0            sctransform_0.3.5     
    [115] parallel_4.1.2         hms_1.1.2              grid_4.1.2            
    [118] rmarkdown_2.20         googledrive_2.0.0      Rtsne_0.16            
    [121] spatstat.explore_3.0-5 shiny_1.7.4            lubridate_1.9.0       

## Paulo’s Code
