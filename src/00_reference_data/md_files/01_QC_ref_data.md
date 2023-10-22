Quality Control of refrence single-cell data
================
4/22/23

### Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)

source("../../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../../results/00_load_ref_data/"
result_dir <- "../../results/01_QC_ref_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
seuratObj <- readRDS(paste0(input_dir,"seuratObj_merged.RDS"))
```

### Features and counts histogram QC plots

``` r
###########################################
# QUALITY CONTROLL FEATURE & COUNTS PLOTS #
###########################################
p1 <- ggplot() +
  geom_histogram(data = seuratObj@meta.data, aes(nFeature_RNA), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per cell") 

p2 <- ggplot() +
  geom_histogram(data = seuratObj@meta.data, aes(nCount_RNA), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per cell")

gene_attr <- data.frame(nUMI = Matrix::rowSums(seuratObj@assays$RNA@counts), 
                        nCells = Matrix::rowSums(seuratObj@assays$RNA@counts > 0))
p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  ggtitle("Total counts per gene (log10 scale)")

p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nCells), fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Total cells per gene")

plot_grid(ncol = 1, 
          p1 + p2 + p3 + p4)
```

<img src="../Figures/01/01a_Feature_and_counts.png"
data-fig-align="center" />

### QC violin plots

``` r
################################
# CALC. % MITO/RIBO/HEMO-GENES #
################################
seuratObj <- PercentageFeatureSet(seuratObj, "^MT-", col.name = "percent_mito")
seuratObj <- PercentageFeatureSet(seuratObj, "^HB[^(P)]", col.name = "percent_hb")
seuratObj <- PercentageFeatureSet(seuratObj, "^RP[SL]", col.name = "percent_ribo")


################
# SEURAT PLOT #
################
# VlnPlot(seuratObj, features = c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo"), 
#         pt.size = 0.1, ncol = 1, y.max =100000) + NoLegend()
# 
# FeatureScatter(seuratObj, "nCount_RNA", "nFeature_RNA", group.by = "sample_name", pt.size = 0.5)

################
# GGPLOT PLOT #
################
friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#F8D0A4", "#E3E6AD", "#be6a7d", "#f1a6b1", "#A8EDFC", "#7fe2e9", "#c4ce96", "#9aacce", "#e7b993", "#ffc8d9")

feature <-  c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo")
gg_violin_fun <- function(obj, feature, fill="sample_name", col_pal=friendly_cols, n=2){
  m <- max(obj[[feature]])/n
obj %>%
tidyseurat::ggplot(aes(sample_name, .data[[feature]], fill=.data[[fill]])) +
  geom_violin() + ylim(c(0, m)) + ggtitle(feature) +
  geom_jitter(width = 0.3, alpha = 0.3, size=.1) +
  scale_fill_manual(values = col_pal) +
  my_theme + NoLegend() +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank()) 
}

 p <-  map(feature, ~gg_violin_fun(seuratObj, .x, ))
 plot_grid(plotlist=p, ncol = 1)
```

<img src="../Figures/01/01b_QC_plots.png" data-fig-align="center" />

### Filtering

``` r
##########################
# FILTER GENES AND CELLS #
##########################
# filter genes present in less than 3 cells:
filt_low_genes <- function(x, n_cell = 3) x[rowSums(x) >= n_cell]
# remove specified genes:
remove_genes <- function(x, gene_name) x[!(grepl(gene_name, rownames(x[["RNA"]]))), ]
# identify transcripts within the 0.005 percentile:
percentile <- function(x, nF) between(nF,quantile(nF,probs = c(0.005)), quantile(nF,probs = c(0.995)))

seuratObj_ <- seuratObj

seuratObj <- seuratObj_ %>%
  filter(., percentile(., .$nFeature_RNA)) %>%
  # cells with less than 200 and more than 15% mt genes and less than 5% hb:
  filter(., nFeature_RNA > 200 & percent_mito < 15 & percent_ribo > 0.05 & percent_hb < 5) %>% 
  filt_low_genes(., n_cell = 3) %>%
  remove_genes(., "^MT-|MALAT1") 

table(seuratObj$sample_name)
```


     CX1  CX2  CX3  CX4  CX5 CX6A CX6B CX7A CX7B  CX8 hg19   N1   N2   N3   N4   N5 
    1800 1587 3129 4197 2889 1141  660 2154 1333 1645 2698 3257 4120 3287 2393 4344 

### Replotting QC after filtering

``` r
############################
# GGPLOT PLOT FILTERED OBJ #
############################
 p_ <-  map(feature, ~gg_violin_fun(seuratObj, .x, n=1))
 plot_grid(plotlist=p_, ncol = 1)
```

<img src="../Figures/01/01c_plot_filtered.png" data-fig-align="center" />

### Plot top abundant genes

``` r
# C = seuratObj@assays$RNA@counts
# C@x = C@x/rep.int(colSums(C), diff(C@p))
# most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
# top_genes <- rownames(C[most_expressed, ])
# 
# 
# boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per spot", 
#     col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

#############################
# GET TOP 20 ABUNDANT GENES #
#############################
top_genes <- seuratObj@assays$RNA@counts %>%
  Matrix::rowSums(.) %>%
  sort(., decreasing = T) %>%
  .[1:20]

percent.fun <- function(df, sample_name, gene, count){
  sample_name <- enquo(sample_name)
  gene <- enquo(gene)
  count <- enquo(count)

  percent <- df %>%
    select(!!sample_name, !!gene, !!count) %>%
    group_by(!!sample_name) %>%
    mutate(Percent = (!!count/sum(!!count)*100)) %>%
    select(-!!count) %>%
    ungroup() 
  
  return(percent$Percent)
}
col = (scales::hue_pal())(20)[20:1]

##################
# PLOT TOP GENES #
##################
(genes_plot <- seuratObj %>%
   join_features(features = names(top_genes) ) %>%
   mutate(.feature = factor(.feature, levels = rev(names(top_genes)))) %>%
   mutate("% total count per cell" = percent.fun(., .cell, .feature, .abundance_RNA),
          .after=.abundance_RNA) %>%
   ggplot(aes(y=`% total count per cell`, x=.feature, fill=.feature)) +
   stat_boxplot(geom = "errorbar", width = 0.2) +
   geom_boxplot(outlier.alpha = 0.1, outlier.size = .5) +
   scale_fill_manual(values = col) + my_theme +
   theme(plot.title = element_text(hjust = 0.5),
         axis.title.y = element_blank()) +
   NoLegend() + coord_flip() )
```

<img src="../Figures/01/01d_top_abundante_genes.png"
data-fig-align="center" />

### Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(seuratObj, paste0(result_dir,"seuratObj_filtered.RDS"))
# seurat_object_list <- readRDS(paste0(result_dir,"seurat_object_list_not_modified.RDS"))
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
     [1] cowplot_1.1.1      tidyseurat_0.5.3   ttservice_0.2.2    SeuratObject_4.1.3
     [5] Seurat_4.3.0       forcats_1.0.0      stringr_1.5.0      dplyr_1.1.1       
     [9] purrr_1.0.1        readr_2.1.3        tidyr_1.3.0        tibble_3.2.1      
    [13] ggplot2_3.4.2      tidyverse_1.3.2   

    loaded via a namespace (and not attached):
      [1] readxl_1.4.1           backports_1.4.1        plyr_1.8.8            
      [4] igraph_1.4.1           lazyeval_0.2.2         sp_1.5-1              
      [7] splines_4.1.2          listenv_0.9.0          scattermore_0.8       
     [10] digest_0.6.31          htmltools_0.5.5        fansi_1.0.4           
     [13] magrittr_2.0.3         tensor_1.5             googlesheets4_1.0.1   
     [16] cluster_2.1.4          ROCR_1.0-11            tzdb_0.3.0            
     [19] globals_0.16.2         modelr_0.1.10          matrixStats_0.63.0    
     [22] timechange_0.2.0       spatstat.sparse_3.0-0  colorspace_2.1-0      
     [25] rvest_1.0.3            ggrepel_0.9.3          haven_2.5.1           
     [28] xfun_0.38              crayon_1.5.2           jsonlite_1.8.4        
     [31] progressr_0.13.0       spatstat.data_3.0-0    survival_3.5-5        
     [34] zoo_1.8-11             glue_1.6.2             polyclip_1.10-4       
     [37] gtable_0.3.3           gargle_1.2.1           leiden_0.4.3          
     [40] future.apply_1.10.0    abind_1.4-5            scales_1.2.1          
     [43] DBI_1.1.3              spatstat.random_3.0-1  miniUI_0.1.1.1        
     [46] Rcpp_1.0.10            viridisLite_0.4.1      xtable_1.8-4          
     [49] reticulate_1.28        htmlwidgets_1.6.2      httr_1.4.5            
     [52] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
     [55] pkgconfig_2.0.3        farver_2.1.1           uwot_0.1.14           
     [58] dbplyr_2.2.1           deldir_1.0-6           utf8_1.2.3            
     [61] tidyselect_1.2.0       labeling_0.4.2         rlang_1.1.0           
     [64] reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
     [67] cellranger_1.1.0       tools_4.1.2            cli_3.6.1             
     [70] generics_0.1.3         broom_1.0.4            ggridges_0.5.4        
     [73] evaluate_0.20          fastmap_1.1.1          yaml_2.3.7            
     [76] goftest_1.2-3          knitr_1.42             fs_1.6.1              
     [79] fitdistrplus_1.1-8     RANN_2.6.1             pbapply_1.7-0         
     [82] future_1.32.0          nlme_3.1-162           mime_0.12             
     [85] xml2_1.3.3             compiler_4.1.2         rstudioapi_0.14       
     [88] plotly_4.10.1          png_0.1-8              spatstat.utils_3.0-1  
     [91] reprex_2.0.2           stringi_1.7.12         lattice_0.20-45       
     [94] Matrix_1.5-3           vctrs_0.6.1            pillar_1.9.0          
     [97] lifecycle_1.0.3        spatstat.geom_3.0-3    lmtest_0.9-40         
    [100] RcppAnnoy_0.0.20       data.table_1.14.6      irlba_2.3.5.1         
    [103] httpuv_1.6.9           patchwork_1.1.2        R6_2.5.1              
    [106] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
    [109] parallelly_1.35.0      codetools_0.2-19       MASS_7.3-58.3         
    [112] assertthat_0.2.1       withr_2.5.0            sctransform_0.3.5     
    [115] parallel_4.1.2         hms_1.1.2              grid_4.1.2            
    [118] rmarkdown_2.21         googledrive_2.0.0      Rtsne_0.16            
    [121] spatstat.explore_3.0-5 shiny_1.7.4            lubridate_1.9.0       
