Quality Control Spatial data
================
10/22/23

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
friendly_cols <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC") 
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

<img src="../Figures/01/01a_Feature_and_counts.png"
data-fig-align="center" />

### Add QC features to DATA

``` r
################################
# CALC. % MITO/RIBO/HEMO-GENES #
################################
DATA <- PercentageFeatureSet(DATA, "^MT-", col.name = "percent_mito")
DATA <- PercentageFeatureSet(DATA, "^HB[^(P)]", col.name = "percent_hb")
DATA <- PercentageFeatureSet(DATA, "^RP[SL]", col.name = "percent_ribo")
```

### Summary stats before filtering

``` r
#################
# SUMMARY STATS #
#################
feature <-  c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo")
sapply(DATA@meta.data[feature], summary) %>% 
  as_tibble(rownames = "stat") %>% 
  knitr::kable(digits = 1)
```

| stat    | nCount_RNA | nFeature_RNA | percent_mito | percent_hb | percent_ribo |
|:--------|-----------:|-------------:|-------------:|-----------:|-------------:|
| Min.    |        7.0 |          7.0 |          0.0 |        0.0 |          0.0 |
| 1st Qu. |     1267.5 |        838.0 |          2.9 |        0.0 |          8.4 |
| Median  |     2465.5 |       1355.0 |          3.6 |        0.0 |         11.4 |
| Mean    |     4866.6 |       1797.7 |          3.9 |        0.1 |         11.2 |
| 3rd Qu. |     6199.0 |       2430.2 |          4.6 |        0.0 |         13.9 |
| Max.    |    51014.0 |       7865.0 |         28.6 |        8.3 |         31.2 |

### Plot feature data on the tissue sections

``` r
# dev.new(width=7, height=3.5, noRStudioGD = TRUE)
# percentage of mitochondria
(plots_m <- DATA %>%
  plot_spatial.fun(., 
      sampleid = unique(.$orig.ident),
      geneid = "percent_mito",
      zoom = "zoom",
      ncol = 4,
      img_alpha = 0,
      point_size = .5)
    )
```

<img src="../Figures/01/01b_sp_mt_plot.png" data-fig-align="center" />

``` r
# number of genes per spot
(plots_f <- DATA %>%
  plot_spatial.fun(., 
      sampleid = unique(.$orig.ident),
      geneid = "nFeature_RNA",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      ncol = 4,
      img_alpha = 0,
      point_size = .5)
    )
```

<img src="../Figures/01/01c_sp_feat_plot.png" data-fig-align="center" />

``` r
# number of reads per spot
(plots_c <- DATA %>%
  plot_spatial.fun(., 
      sampleid = unique(.$orig.ident),
      geneid = "nCount_RNA",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      ncol = 4,
      img_alpha = 0,
      point_size = .5)
    )
```

<img src="../Figures/01/01d_sp_count_plot.png" data-fig-align="center" />

## Filtering

Select all spots with less than 15% mitochondrial reads, less than 10%
hb-reads and 100 detected genes.<br/> Filter genes present in less than
2 spots and remove hemoglobin and MALAT1 genes.

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

DATA <- DATA %>%
  # filter out spots with less than 100 genes and more than 15% mt and 10% hb:
  mutate(filt = case_when(nFeature_RNA < 100 ~ 'filt',
                          percent_mito > 15 ~ 'filt',
                          percent_hb > 10 ~ 'filt',
                              TRUE ~ "keep")) %>%
  {. ->> temp } %>%
  #filter(., percentile(., .$nFeature_RNA)) %>%
  filter(filt == "keep") %>%
  filt_low_genes(., n_cell = 2) %>%
  remove_genes(., "MALAT1|^HB[^(P)]") %>% # "^MT-|MALAT1|^HB[^(P)]"
  select(-filt)
```

### Summary stats

Dimension of DATA before filtering, genes: 36601, spots: 6612<br/>
Dimension of DATA after filtering, genes: 22026, spots: 6598

``` r
###########################
# SUMMARY STATS ALL SPOTS #
###########################
feature <-  c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo")
sapply(DATA@meta.data[feature], summary) %>% 
  as_tibble(rownames = "stat") %>% 
  knitr::kable(digits = 1)
```

| stat    | nCount_RNA | nFeature_RNA | percent_mito | percent_hb | percent_ribo |
|:--------|-----------:|-------------:|-------------:|-----------:|-------------:|
| Min.    |      111.0 |        104.0 |          0.7 |        0.0 |          2.1 |
| 1st Qu. |     1221.2 |        840.2 |          2.9 |        0.0 |          8.4 |
| Median  |     2423.5 |       1355.0 |          3.6 |        0.0 |         11.4 |
| Mean    |     4836.2 |       1799.4 |          3.9 |        0.1 |         11.2 |
| 3rd Qu. |     6181.0 |       2429.0 |          4.6 |        0.0 |         13.9 |
| Max.    |    50999.0 |       7860.0 |         14.1 |        8.3 |         31.2 |

``` r
#########################
# SUMMARY STATS GROUPS #
#########################
DATA@meta.data %>%
  split(.$groups, drop = T) %>% 
  map(., ~.x %>%
        select(., any_of(feature[1:2])) %>%
        map(summary) %>%
        bind_rows(.id = "stat")
        #tibble(.x, .name_repair="unique")
      ) %>%
  bind_rows(., .id = "groups") %>%
  arrange(stat) %>%
  knitr::kable(digits = 1)
```

| groups | stat         | Min. | 1st Qu. | Median |     Mean | 3rd Qu. |  Max. |
|:-------|:-------------|-----:|--------:|-------:|---------:|--------:|------:|
| ctrl   | nCount_RNA   |  111 |     974 |   2279 | 4743.045 |    6458 | 46060 |
| DMPA   | nCount_RNA   |  183 |    1497 |   2585 | 4952.667 |    5749 | 50999 |
| ctrl   | nFeature_RNA |  104 |     708 |   1276 | 1735.319 |    2420 |  7860 |
| DMPA   | nFeature_RNA |  158 |     993 |   1422 | 1879.477 |    2453 |  7363 |

### Plotting QC before and after filtering

``` r
# dev.new(width=6, height=3, noRStudioGD = TRUE)
################################
# VIOLIN PLOT BEFORE FILTERING #
################################
feature <-  c("nCount_RNA", "nFeature_RNA","percent_mito","percent_hb", "percent_ribo")
p_ <-  map(feature, ~violin.fun(temp, feature="nCount_RNA", fill="orig.ident", col_pal=friendly_cols, txt_size=20))
# plot_grid(plotlist=p_, ncol = 1)

################################
# VIOLIN PLOT AFTER FILTERING #
################################
p <-  map(feature, ~violin.fun(DATA, feature=.x, fill="orig.ident", col_pal=friendly_cols, n=1, txt_size=20))
plot_grid(plotlist=c(p_, p), nrow = 5, byrow = F)
```

<img src="../Figures/01/01e_QC_plot_filtered.png"
data-fig-align="center" />

### Filtered spots

``` r
temp %>%
  filter(filt == "filt") %>%
  arrange(nFeature_RNA) %>%
  as_tibble() %>%
  select(-sp_annot2) %>%
  knitr::kable(digits = 1)
```

| .cell                 | groups | sp_annot | orig.ident | nCount_RNA | nFeature_RNA | percent_mito | percent_hb | percent_ribo | filt |
|:----------------------|:-------|:---------|:-----------|-----------:|-------------:|-------------:|-----------:|-------------:|:-----|
| P097_TCGCGTAGCAGTGTCC | DMPA   | SubMuc   | P097       |          7 |            7 |         28.6 |        0.0 |         14.3 | filt |
| P118_CATGGTAAGTAGCGTT | ctrl   | SubMuc   | P118       |         14 |           13 |          0.0 |        0.0 |          0.0 | filt |
| P080_TCGTGTACTATGGATG | ctrl   | SubMuc   | P080       |         17 |           17 |          0.0 |        0.0 |          5.9 | filt |
| P031_AGAAGAGCGCCGTTCC | ctrl   | SubMuc   | P031       |         30 |           29 |          3.3 |        0.0 |         10.0 | filt |
| P097_GCATCGGCCGTGTAGG | DMPA   | SubMuc   | P097       |         30 |           29 |          3.3 |        0.0 |         10.0 | filt |
| P114_TCGCGTAGCAGTGTCC | DMPA   | SubMuc   | P114       |         41 |           38 |          2.4 |        0.0 |         19.5 | filt |
| P080_AGAAGAGCGCCGTTCC | ctrl   | SubMuc   | P080       |         48 |           39 |          2.1 |        0.0 |          8.3 | filt |
| P114_TCGTGTACTATGGATG | DMPA   | SubMuc   | P114       |         58 |           56 |          6.9 |        0.0 |         13.8 | filt |
| P107_CATGGTAAGTAGCGTT | DMPA   | epi      | P107       |         63 |           59 |          4.8 |        0.0 |         15.9 | filt |
| P105_AGAAGAGCGCCGTTCC | ctrl   | SubMuc   | P105       |         63 |           62 |         11.1 |        0.0 |         19.0 | filt |
| P097_TCGTGTACTATGGATG | DMPA   | SubMuc   | P097       |         88 |           82 |          3.4 |        0.0 |          8.0 | filt |
| P108_CATGGTAAGTAGCGTT | DMPA   | SubMuc   | P108       |         88 |           82 |          2.3 |        0.0 |         17.0 | filt |
| P031_AGATACCGGTGTTCAC | ctrl   | SubMuc   | P031       |        790 |          439 |         15.6 |        0.0 |         13.0 | filt |
| P031_GTGCCATCACACGGTG | ctrl   | SubMuc   | P031       |       1136 |          566 |         16.6 |        0.1 |         12.3 | filt |

### Plot filtered spots

``` r
# dev.new(width=7, height=3.5, noRStudioGD = TRUE)
(plots <- temp %>%
  plot_spatial.fun(., 
      sampleid = unique(.$orig.ident),
      geneid = "filt",#"KRT15", #"PTPRC",#"sp_annot",#"CDH1",
      zoom = "zoom",
      ncol = 4,
      img_alpha = 0,
      point_size = 0.5)
    )
```

<img src="../Figures/01/01f_filtered_spots.png"
data-fig-align="center" />

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
   mutate("% total count per spot" = percent.fun(., .cell, .feature, .abundance_RNA),
          .after=.abundance_RNA) %>%
   ggplot(aes(y=`% total count per spot`, x=.feature, fill=.feature)) +
   stat_boxplot(geom = "errorbar", width = 0.2) +
   geom_boxplot(outlier.alpha = 0.1, outlier.size = .5) +
   scale_fill_manual(values = col) + my_theme +
   theme(plot.title = element_text(hjust = 0.5),
         axis.title.y = element_blank()) +
   NoLegend() + coord_flip() )
```

<img src="../Figures/01/01g_top_abundante_genes.png"
data-fig-align="center" />

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(DATA, paste0(result_dir,"seuratObj_filtered.RDS"))
#DATA <- readRDS(paste0(result_dir,"seuratObj_filtered.RDS"))
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
     [1] niceRplots_0.1.0   patchwork_1.1.2    cowplot_1.1.1      RColorBrewer_1.1-3
     [5] Seurat_4.3.0       tidyseurat_0.5.3   SeuratObject_4.1.3 sp_1.5-1          
     [9] ttservice_0.2.2    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.2       
    [13] purrr_1.0.1        readr_2.1.3        tidyr_1.3.0        tibble_3.2.1      
    [17] ggplot2_3.4.3      tidyverse_1.3.2   

    loaded via a namespace (and not attached):
      [1] googledrive_2.0.0      Rtsne_0.16             colorspace_2.1-0      
      [4] deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4        
      [7] fs_1.6.2               spatstat.data_3.0-0    rstudioapi_0.14       
     [10] farver_2.1.1           leiden_0.4.3           listenv_0.9.0         
     [13] ggrepel_0.9.3          fansi_1.0.4            lubridate_1.9.0       
     [16] xml2_1.3.3             codetools_0.2-19       splines_4.1.2         
     [19] knitr_1.42             polyclip_1.10-4        jsonlite_1.8.5        
     [22] broom_1.0.4            ica_1.0-3              cluster_2.1.4         
     [25] dbplyr_2.2.1           png_0.1-8              uwot_0.1.14           
     [28] spatstat.sparse_3.0-0  sctransform_0.3.5      shiny_1.7.4           
     [31] compiler_4.1.2         httr_1.4.5             backports_1.4.1       
     [34] lazyeval_0.2.2         assertthat_0.2.1       Matrix_1.6-1          
     [37] fastmap_1.1.1          gargle_1.2.1           cli_3.6.1             
     [40] later_1.3.0            htmltools_0.5.5        tools_4.1.2           
     [43] igraph_1.4.1           gtable_0.3.4           glue_1.6.2            
     [46] reshape2_1.4.4         RANN_2.6.1             Rcpp_1.0.10           
     [49] scattermore_0.8        cellranger_1.1.0       vctrs_0.6.3           
     [52] nlme_3.1-163           spatstat.explore_3.0-5 progressr_0.13.0      
     [55] lmtest_0.9-40          spatstat.random_3.0-1  xfun_0.38             
     [58] globals_0.16.2         rvest_1.0.3            timechange_0.2.0      
     [61] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
     [64] irlba_2.3.5.1          goftest_1.2-3          googlesheets4_1.0.1   
     [67] future_1.32.0          MASS_7.3-60            zoo_1.8-11            
     [70] scales_1.2.1           spatstat.utils_3.0-1   hms_1.1.2             
     [73] promises_1.2.0.1       parallel_4.1.2         yaml_2.3.7            
     [76] gridExtra_2.3          pbapply_1.7-0          reticulate_1.28       
     [79] stringi_1.7.12         rlang_1.1.1            pkgconfig_2.0.3       
     [82] matrixStats_0.63.0     evaluate_0.21          lattice_0.21-8        
     [85] tensor_1.5             ROCR_1.0-11            labeling_0.4.3        
     [88] htmlwidgets_1.6.2      tidyselect_1.2.0       parallelly_1.36.0     
     [91] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
     [94] R6_2.5.1               generics_0.1.3         DBI_1.1.3             
     [97] pillar_1.9.0           haven_2.5.1            withr_2.5.0           
    [100] fitdistrplus_1.1-8     abind_1.4-5            survival_3.5-5        
    [103] future.apply_1.10.0    modelr_0.1.10          crayon_1.5.2          
    [106] KernSmooth_2.23-20     utf8_1.2.3             spatstat.geom_3.0-3   
    [109] plotly_4.10.1          tzdb_0.3.0             rmarkdown_2.21        
    [112] grid_4.1.2             readxl_1.4.1           data.table_1.14.6     
    [115] reprex_2.0.2           digest_0.6.31          xtable_1.8-4          
    [118] httpuv_1.6.9           munsell_0.5.0          viridisLite_0.4.2     
