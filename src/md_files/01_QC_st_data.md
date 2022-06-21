# Spatial transcriptomics

------------------------------------------------------------------------

### Load packages

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(tidyseurat)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(patchwork)
# remotes::install_github("czarnewski/niceRplots",force=T)
library(niceRplots)
source("../bin/spatial_visualization.R")
```

### Load ST data

``` r
#########
# PATHS #
#########
input_dir <- "../results/00_load_st_data/"
result_dir <- "../results/01_QC_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
#metadata <- read_csv("../data/Clinincal_data_Spatial_DMPA.csv")
DATA <- readRDS(paste0(input_dir,"seuratObj_merged.RDS"))
#sample_names <- unique(DATA$orig.ident)

#################
# COLOUR PALLET #
#################
friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#F8D0A4", "#E3E6AD", "#be6a7d", "#f1a6b1", "#A8EDFC", "#7fe2e9", "#c4ce96", "#9aacce", "#e7b993", "#ffc8d9")
```

``` r
###########################################
# QUALITY CONTROLL FEATURE & COUNTS PLOTS #
###########################################
p1 <- ggplot() +
  geom_histogram(data = DATA@meta.data, aes(nFeature_RNA), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per spot") 

p2 <- ggplot() +
  geom_histogram(data = DATA@meta.data, aes(nCount_RNA), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per spots")

gene_attr <- data.frame(nUMI = Matrix::rowSums(DATA@assays$RNA@counts), 
                        nSpots = Matrix::rowSums(DATA@assays$RNA@counts > 0))
p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  ggtitle("Total counts per gene (log10 scale)")

p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nSpots), fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Total spots per gene")

plot_grid(ncol = 1, 
          p1 + p2 + p3 + p4)
```

<img src="./Figures/01a_Feature_and_counts.png" style="display: block; margin: auto;" />

# Quality control

------------------------------------------------------------------------

## QC violin plots

``` r
################################
# CALC. % MITO/RIBO/HEMO-GENES #
################################
DATA <- PercentageFeatureSet(DATA, "^MT-", col.name = "percent_mito")
DATA <- PercentageFeatureSet(DATA, "^HB[^(P)]", col.name = "percent_hb")
DATA <- PercentageFeatureSet(DATA, "^RP[SL]", col.name = "percent_ribo")


################
# SEURAT PLOT #
################
# VlnPlot(DATA, features = c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo"), 
#         pt.size = 0.1, ncol = 1, y.max =100000) + NoLegend()
# 
# FeatureScatter(DATA, "nCount_RNA", "nFeature_RNA", group.by = "sample_name", pt.size = 0.5)

################
# VIOLIN PLOT #
################
feature <-  c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo")
# gg_violin_fun <- function(obj, feature, fill="sample_name", col_pal=friendly_cols, n=2){
#   m <- max(obj[[feature]])/n
# obj %>%
# tidyseurat::ggplot(aes(orig.ident, .data[[feature]], fill=.data[[fill]])) +
#   geom_violin() + ylim(c(0, m)) + ggtitle(feature) +
#   geom_jitter(width = 0.3, alpha = 0.3, size=.1) +
#   scale_fill_manual(values = col_pal) +
#   my_theme + NoLegend() +
#   theme(axis.text.x = element_text(angle = 30, hjust=1),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.y = element_blank()) 
# }

 p <-  map(feature, ~violin.fun(DATA, .x, fill="orig.ident", col_pal=friendly_cols))
 plot_grid(plotlist=p, ncol = 1)
```

<img src="./Figures/01b_QC_plots.png" style="display: block; margin: auto;" />

We can also plot the same data onto the tissue section.

``` r
# percentage of mitochondria
plots <- DATA %>%
  mutate(group = orig.ident) %>%
  nest(., data = -group) %>%
  pmap(., 
    ~plot_spatial.fun(..2,
      sampleid = ..1,
      geneid = "percent_mito",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      img_alpha = 0,
      point_size = 1)
    )
legend <- get_legend(plots[[1]] + theme(legend.position="right"))
combined <- wrap_plots(plots, ncol=2) & theme(legend.position="none")
combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .2)) 
combined
```

<img src="./Figures/01c_sp_mt_plot.png" style="display: block; margin: auto;" />

``` r
# number of genes per spot
plots_f <- DATA %>%
  mutate(group = orig.ident) %>%
  nest(., data = -group) %>%
  pmap(., 
    ~plot_spatial.fun(..2,
      sampleid = ..1,
      geneid = "nFeature_RNA",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      img_alpha = 0,
      point_size = 1)
    )
legend <- get_legend(plots_f[[1]] + theme(legend.position="right"))
combined <- wrap_plots(plots_f, ncol=2) & theme(legend.position="none")
combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .2)) 
combined
```

<img src="./Figures/01d_sp_feat_plot.png" style="display: block; margin: auto;" />

``` r
# number of reads per spot
plots_c <- DATA %>%
  mutate(group = orig.ident) %>%
  nest(., data = -group) %>%
  pmap(., 
    ~plot_spatial.fun(..2,
      sampleid = ..1,
      geneid = "nCount_RNA",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      img_alpha = 0,
      point_size = 1)
    )
legend <- get_legend(plots_c[[1]] + theme(legend.position="right"))
combined <- wrap_plots(plots_c, ncol=2) & theme(legend.position="none")
combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .2)) 
combined
```

<img src="./Figures/01e_sp_count_plot.png" style="display: block; margin: auto;" />

As you can see, the spots with low number of counts/features and high
mitochondrial content is mainly towards the edges of the tissue. It is
quite likely that these regions are damaged tissue. You may also see
regions within a tissue with low quality if you have tears or folds in
your section.

But remember, for some tissue types, the amount of genes expressed and
proportion mitochondria may also be a biological features, so bear in
mind what tissue you are working on and what these features mean.

### Filter

Select all spots with less than 25% mitocondrial reads, less than 20%
hb-reads and 1000 detected genes. You must judge for yourself based on
your knowledge of the tissue what are appropriate filtering criteria for
your dataset.

## Filtering

``` r
##########################
# FILTER GENES AND CELLS #
##########################
# filter genes present in less than 2 spots:
filt_low_genes <- function(x, n_cell = 2) x[rowSums(x) >= n_cell]
# remove specified genes:
remove_genes <- function(x, gene_name) x[!(grepl(gene_name, rownames(x[["RNA"]]))), ]
# identify transcripts within the 0.005 percentile:
percentile <- function(x, nF) between(nF,quantile(nF,probs = c(0.005)), quantile(nF,probs = c(0.995)))

DATA_ <- DATA

DATA <- DATA_ %>%
  filter(., percentile(., .$nFeature_RNA)) %>%
  # filter out spots with less than 500 genes and 0.05% ribo and more than 25% mt and 20% hb:
  mutate(filt = ifelse(nFeature_RNA > 200 & percent_ribo > 0.05 & percent_mito < 25 & percent_hb < 20,
                       "keep", "filt")) %>%
  {. ->> temp } %>%
  filter(., nFeature_RNA > 300 & percent_ribo > 0.05 & percent_mito < 25 & percent_hb < 20) %>% 
  filt_low_genes(., n_cell = 2) %>%
  remove_genes(., "MALAT1|^HB[^(P)]") # "^MT-|MALAT1|^HB[^(P)]"

#################
# SUMMARY STATS #
#################
dim(DATA_)
```

    ## [1] 36601  6700

``` r
dim(DATA)
```

    ## [1] 21950  6527

``` r
#rm(DATA_)

DATA %>%
  group_by(orig.ident) %>%
  summarise(across(c(nCount_RNA:percent_hb), list(min = min, max = max))) %>%
  rename_with(., ~str_replace(., "_RNA", "")) # |percent_
```

    ## # A tibble: 8 × 9
    ##   orig.ident nCount_min nCount_max nFeature_min nFeature_max percent_mito_min
    ##   <chr>           <dbl>      <dbl>        <int>        <int>            <dbl>
    ## 1 P031              355      32525          301         5880            0.832
    ## 2 P080              364      18944          297         5487            0.748
    ## 3 P097              383      35452          313         6353            0.921
    ## 4 P105              366      34872          300         6484            1.17 
    ## 5 P107              601      32542          316         5997            0.737
    ## 6 P108              415      44508          313         6419            1.31 
    ## 7 P114              496      17684          304         5184            1.08 
    ## 8 P118              406      30163          304         6518            1.58 
    ## # … with 3 more variables: percent_mito_max <dbl>, percent_hb_min <dbl>,
    ## #   percent_hb_max <dbl>

### Replotting QC after filtering

``` r
############################
# GGPLOT PLOT FILTERED OBJ #
############################
 p_ <-  map(feature, ~violin.fun(DATA, .x, fill="orig.ident", col_pal=friendly_cols, n=1))
 plot_grid(plotlist=p_, ncol = 1)
```

<img src="./Figures/01f_QC_plot_filtered.png" style="display: block; margin: auto;" />

### And replot onto tissue section:

``` r
plots <- temp %>%
  mutate(group = orig.ident) %>%
  nest(., data = -group) %>%
  pmap(., 
    ~plot_spatial.fun(..2,
      sampleid = ..1,
      geneid = "filt",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      img_alpha = 0,
      point_size = 1)
    )
legend <- get_legend(plots[[1]] + theme(legend.position="right"))
combined <- wrap_plots(plots, ncol=2) & theme(legend.position="none")
combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .2)) 
combined
```

<img src="./Figures/01g_filtered_spots.png" style="display: block; margin: auto;" />

### Plot top expressed genes

``` r
#############################
# GET TOP 20 ABUNDANT GENES #
#############################
top_genes <- DATA@assays$RNA@counts %>%
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
(genes_plot <- DATA %>%
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

<img src="./Figures/01h_top_abundante_genes.png" style="display: block; margin: auto;" />

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(DATA_, paste0(result_dir,"seuratObj_merged.RDS"))
saveRDS(DATA, paste0(result_dir,"seuratObj_filtered.RDS"))
#DATA <- readRDS(paste0(result_dir,"seuratObj_filtered.RDS"))
```

### Session info

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
    ##  [1] niceRplots_0.1.0   patchwork_1.1.1    cowplot_1.1.1      RColorBrewer_1.1-3
    ##  [5] Seurat_4.1.0       tidyseurat_0.5.3   SeuratObject_4.0.4 ttservice_0.1.2   
    ##  [9] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.8        purrr_0.3.4       
    ## [13] readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
    ## [17] tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15            colorspace_2.0-3      deldir_1.0-6         
    ##   [4] ellipsis_0.3.2        ggridges_0.5.3        fs_1.5.2             
    ##   [7] spatstat.data_2.1-4   rstudioapi_0.13       farver_2.1.0         
    ##  [10] leiden_0.3.9          listenv_0.8.0         ggrepel_0.9.1        
    ##  [13] fansi_1.0.3           lubridate_1.8.0       xml2_1.3.3           
    ##  [16] codetools_0.2-18      splines_4.1.2         knitr_1.38           
    ##  [19] polyclip_1.10-0       jsonlite_1.8.0        broom_0.7.12         
    ##  [22] ica_1.0-2             cluster_2.1.2         dbplyr_2.1.1         
    ##  [25] png_0.1-7             uwot_0.1.11           spatstat.sparse_2.1-0
    ##  [28] sctransform_0.3.3     shiny_1.7.1           compiler_4.1.2       
    ##  [31] httr_1.4.2            backports_1.4.1       lazyeval_0.2.2       
    ##  [34] assertthat_0.2.1      Matrix_1.4-1          fastmap_1.1.0        
    ##  [37] cli_3.3.0             later_1.3.0           htmltools_0.5.2      
    ##  [40] tools_4.1.2           igraph_1.3.0          gtable_0.3.0         
    ##  [43] glue_1.6.2            reshape2_1.4.4        RANN_2.6.1           
    ##  [46] Rcpp_1.0.8.3          scattermore_0.8       cellranger_1.1.0     
    ##  [49] jquerylib_0.1.4       vctrs_0.4.1           nlme_3.1-157         
    ##  [52] lmtest_0.9-40         spatstat.random_2.2-0 xfun_0.30            
    ##  [55] globals_0.14.0        rvest_1.0.2           mime_0.12            
    ##  [58] miniUI_0.1.1.1        lifecycle_1.0.1       irlba_2.3.5          
    ##  [61] goftest_1.2-3         future_1.24.0         MASS_7.3-57          
    ##  [64] zoo_1.8-9             scales_1.2.0          spatstat.core_2.4-2  
    ##  [67] spatstat.utils_2.3-0  hms_1.1.1             promises_1.2.0.1     
    ##  [70] parallel_4.1.2        yaml_2.3.5            gridExtra_2.3        
    ##  [73] reticulate_1.24       pbapply_1.5-0         sass_0.4.1           
    ##  [76] rpart_4.1.16          stringi_1.7.6         highr_0.9            
    ##  [79] rlang_1.0.2           pkgconfig_2.0.3       matrixStats_0.61.0   
    ##  [82] evaluate_0.15         lattice_0.20-45       tensor_1.5           
    ##  [85] ROCR_1.0-11           labeling_0.4.2        htmlwidgets_1.5.4    
    ##  [88] tidyselect_1.1.2      parallelly_1.31.0     RcppAnnoy_0.0.19     
    ##  [91] plyr_1.8.7            magrittr_2.0.3        R6_2.5.1             
    ##  [94] generics_0.1.2        DBI_1.1.2             mgcv_1.8-40          
    ##  [97] pillar_1.7.0          haven_2.4.3           withr_2.5.0          
    ## [100] fitdistrplus_1.1-8    abind_1.4-5           survival_3.2-13      
    ## [103] future.apply_1.8.1    modelr_0.1.8          crayon_1.5.1         
    ## [106] KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_2.4-0  
    ## [109] plotly_4.10.0         tzdb_0.2.0            rmarkdown_2.11       
    ## [112] grid_4.1.2            readxl_1.3.1          data.table_1.14.2    
    ## [115] reprex_2.0.1          digest_0.6.29         xtable_1.8-4         
    ## [118] httpuv_1.6.5          munsell_0.5.0         viridisLite_0.4.0    
    ## [121] bslib_0.3.1
