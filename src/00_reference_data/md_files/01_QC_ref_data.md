## Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)

#########
# PATHS #
#########
input_dir <- "../../results/00_load_ref_data/"
result_dir <- "../../results/02_QC_ref_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
seuratObj <- readRDS(paste0(input_dir,"seuratObj_merged.RDS"))
```

``` r
my_theme <-
  list(
    #scale_fill_manual(values = friendly_cols),
    #scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        #aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )
```

## Features and counts histogram QC plots

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

<img src="./Figures/01a_Feature_and_counts.png" style="display: block; margin: auto;" />

## QC violin plots

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

<img src="./Figures/01b_QC_plots.png" style="display: block; margin: auto;" />

## Filtering

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

    ## 
    ##  CX1  CX2  CX3  CX4  CX5 CX6A CX6B CX7A CX7B  CX8   N1   N2   N3   N4   N5 
    ## 1800 1587 3129 4197 2889 1141  663 2154 1334 1645 3257 4122 3290 2394 4346

## Replotting QC after filtering

``` r
############################
# GGPLOT PLOT FILTERED OBJ #
############################
 p_ <-  map(feature, ~gg_violin_fun(seuratObj, .x, n=1))
 plot_grid(plotlist=p_, ncol = 1)
```

<img src="./Figures/01c_plot_filtered.png" style="display: block; margin: auto;" />

## Plot top abundant genes

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

<img src="./Figures/01d_top_abundante_genes.png" style="display: block; margin: auto;" />

## Save seurat object

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
    ##  [1] cowplot_1.1.1      tidyseurat_0.5.1   ttservice_0.1.2    SeuratObject_4.0.4
    ##  [5] Seurat_4.1.0       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.8       
    ##  [9] purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.6      
    ## [13] ggplot2_3.3.5      tidyverse_1.3.1   
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
    ##  [40] miniUI_0.1.1.1        Rcpp_1.0.8.3          viridisLite_0.4.0    
    ##  [43] xtable_1.8-4          reticulate_1.24       spatstat.core_2.4-2  
    ##  [46] htmlwidgets_1.5.4     httr_1.4.2            RColorBrewer_1.1-3   
    ##  [49] ellipsis_0.3.2        ica_1.0-2             pkgconfig_2.0.3      
    ##  [52] farver_2.1.0          sass_0.4.1            uwot_0.1.11          
    ##  [55] dbplyr_2.1.1          deldir_1.0-6          utf8_1.2.2           
    ##  [58] tidyselect_1.1.2      labeling_0.4.2        rlang_1.0.2          
    ##  [61] reshape2_1.4.4        later_1.3.0           munsell_0.5.0        
    ##  [64] cellranger_1.1.0      tools_4.1.2           cli_3.2.0            
    ##  [67] generics_0.1.2        broom_0.7.12          ggridges_0.5.3       
    ##  [70] evaluate_0.15         fastmap_1.1.0         yaml_2.3.5           
    ##  [73] goftest_1.2-3         knitr_1.38            fs_1.5.2             
    ##  [76] fitdistrplus_1.1-8    RANN_2.6.1            pbapply_1.5-0        
    ##  [79] future_1.24.0         nlme_3.1-157          mime_0.12            
    ##  [82] xml2_1.3.3            compiler_4.1.2        rstudioapi_0.13      
    ##  [85] plotly_4.10.0         png_0.1-7             spatstat.utils_2.3-0 
    ##  [88] reprex_2.0.1          bslib_0.3.1           stringi_1.7.6        
    ##  [91] highr_0.9             lattice_0.20-45       Matrix_1.4-1         
    ##  [94] vctrs_0.4.0           pillar_1.7.0          lifecycle_1.0.1      
    ##  [97] spatstat.geom_2.4-0   lmtest_0.9-40         jquerylib_0.1.4      
    ## [100] RcppAnnoy_0.0.19      data.table_1.14.2     irlba_2.3.5          
    ## [103] httpuv_1.6.5          patchwork_1.1.1       R6_2.5.1             
    ## [106] promises_1.2.0.1      KernSmooth_2.23-20    gridExtra_2.3        
    ## [109] parallelly_1.31.0     codetools_0.2-18      MASS_7.3-56          
    ## [112] assertthat_0.2.1      withr_2.5.0           sctransform_0.3.3    
    ## [115] mgcv_1.8-40           parallel_4.1.2        hms_1.1.1            
    ## [118] grid_4.1.2            rpart_4.1.16          rmarkdown_2.11       
    ## [121] Rtsne_0.15            shiny_1.7.1           lubridate_1.8.0
