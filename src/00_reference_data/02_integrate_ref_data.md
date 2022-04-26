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
input_dir <- "../../results/01_QC_ref_data/"
result_dir <- "../../results/02_integrate_ref_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
seuratObj <- readRDS(paste0(input_dir,"seuratObj_filtered.RDS"))

seuratObj <- seuratObj %>% 
  #filter(p == "Li_Y")
  filter(p == "Peng")
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

## Identify Highly Variable Genes (HVG) across samples

``` r
################################
# SPLIT INTO SEPERATE DATASETS #
################################
seuratObj_nested <- seuratObj %>%
  mutate(batch = sample_name) %>%
  nest(data = -batch) %>%
  #filter(batch == "CX3" | batch == "CX1") %>%
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
hvgs_heat <- seuratObj_nested %>%
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

<img src="./Figures/02a_HVG_heatmap.png" style="display: block; margin: auto;" />

## Integration

``` r
############################
# CCA REDUCTION INTERATION #
############################
# seuratObj_CCA <- seuratObj_split %>%
#   .$data %>%
#   FindIntegrationAnchors(normalization.method = "LogNormalize",
#                          anchor.features = hig_var,
#                          dims = 1:30,
#                          reduction = "cca") %>%
#   IntegrateData(., dims = 1:30, new.assay.name = "CCA")
# 
# saveRDS(seuratObj_CCA, paste0(result_dir,"seuratObj_CCA.RDS"))

############
# HARMONY #
###########
seuratObj <- seuratObj %>%
  #filter(sample_name == "CX3" | sample_name == "CX1") %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                      nfeatures = 4000,
                      verbose = FALSE) %>%
  ScaleData(verbose = FALSE, features = hig_var ) %>%
  RunPCA(verbose = FALSE, npcs = 50) %>%
  RunUMAP(dims = 1:50, n.components = 2L) 

seuratObj <- seuratObj %>%
  RunHarmony(group.by.vars = "sample_name", 
             reduction = "pca", 
             dims.use = 1:50, 
             assay.use = "RNA") %>%
  RunUMAP(dims = 1:30, 
          reduction = "harmony",
          reduction.name = "umap_harmony")

######################
# SCT NORMALIZATION #
#####################
# seuratObj_sct <- seuratObj_nested %>%
#   mutate(data = map(data, ~ SCTransform(.x, verbose = FALSE))) %>%
#   .$data %>%
#   ( function(x){SelectIntegrationFeatures(x) ->> af; return(x)} ) %>%
#   PrepSCTIntegration(anchor.features = af) %>%
#   FindIntegrationAnchors(normalization.method = "SCT",
#                          anchor.features = af) %>%
#   IntegrateData(normalization.method = "SCT") 
# 
# seuratObj_sct <- readRDS(paste0(result_dir,"seuratObj_sct.RDS"))
```

## Plot before and after integration

``` r
plot_clusters.fun <- function(obj, cluster="sample_name", red = "umap_harmony", color = "Brew_all", lable = TRUE, title = NULL){
  if(color == "Brew_all"){
    pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2") )}
  
  cluster <- sym(cluster)
  feat <- obj %>%
    select(.cell, !!(cluster), sample_name, nCount_RNA, nFeature_RNA) %>%
    group_by(!!(cluster)) %>%
    add_tally() %>%
    arrange(nFeature_RNA) %>%
    arrange(desc(n))
    
  lable_df <- feat %>%
    select(!!(cluster), contains(red)) %>% 
    summarize_all(mean)
  
  if(lable == TRUE){
    text <- geom_text(data = lable_df, aes(label = !!cluster), col="black") +
      NoLegend() + labs(color= "Clusters")}
  else{text <- labs(color= "")}
  
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  p <- ggplot(feat, aes(!!(red_1), !!(red_2), 
                        color = !!cluster), label = lable) + 
    geom_point(alpha = 0.5, size=.5) + ggtitle(title) + text +
    #ggrepel::geom_label_repel(data = lable, aes(label = !!cluster)) +
    scale_colour_manual(values = pal)  +
    my_theme 
  return(p)
}

res <- c("PC", "harmony", "UMAP", "umap_harmony")
title <- c("PCA raw data", "PCA Harmony integrated", "UMAP raw data", "UMAP Harmony integrated")
p <- map2(res, title, ~plot_clusters.fun(seuratObj, red=.x, lable=FALSE, title=.y))
plot_grid(ncol = 2, 
         plotlist = p)
```

<img src="./Figures/02b_Plot_dim_reduction.png" style="display: block; margin: auto;" />

## Plot marker genes

``` r
################################
# VISUALIZE EXPR. OF KEY GENES #
################################
library(viridis)
col=c("grey90","grey80","grey60","navy","black")
genes <- c("CD8A", "SFRP2", "CD3E")
# genes <- c("CD8A", "MYOZ2", "CD3E", "EPCAM", "COL6A1", "CD4")
seuratObj <- seuratObj %>%
  #join_features(features = genes ) %>%
  mutate(., FetchData(., vars = c("CD8A", "SFRP2", "CD3E")) )

# obj <- seuratObj       
# gene <- sym("SFRP2")
plot_genes.fun <- function(obj, gene, mins=NULL, maxs=NULL, red = "umap_harmony"){
  gene <- sym(gene)
  feat <- pull(obj, gene)
  
  # Colour pal:
  if(is.null(mins)){mins <- min(c(feat,0),na.rm = T)}
  if(is.null(maxs)){
    maxs <- quantile(feat,0.99,na.rm = T)
    if(maxs==0){maxs <- max(feat,na.rm = T)}
  }
  if(max(feat,na.rm = T) != 0){
    feat <- (feat - mins) / ( maxs - mins)
    feat[feat > 1] <- 1}
  
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  myPalette <-  colorRampPalette(col[-1])
  obj <- arrange(obj,nCount_RNA)
  
  p <- ggplot(obj, aes(!!(red_1), !!(red_2), color = !!(gene)) ) + 
    
    geom_point(alpha = 0.5, size=.5) + ggtitle(as_label(gene)) +
    #scale_color_viridis(option = "D", na.value="#EBECF0") +
    scale_colour_gradientn(colours = c( col[1], myPalette(99)), limits=c(mins, maxs))  +
    my_theme + theme_void() + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) #+ NoLegend()
  return(p)
}

# get approximate max value for among all features (genes)
red <- seuratObj@reductions[["umap_harmony"]]@cell.embeddings
max <- quantile(red,0.99,na.rm = T)

#p <- map(genes, ~plot_genes.fun(seuratObj, .x, maxs = max))
p <- map(genes, ~plot_genes.fun(seuratObj, .x))
plot_grid(ncol = 3, 
          plotlist = p)
```

<img src="./Figures/02c_plot_marker_genes.png" style="display: block; margin: auto;" />

# Paulos base R code

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(seuratObj, paste0(result_dir,"seuratObj_harmony.RDS"))
#seuratObj <- readRDS(paste0(result_dir,"seuratObj_harmony.RDS"))
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
    ##  [1] viridis_0.6.2      viridisLite_0.4.0  harmony_0.1.0      Rcpp_1.0.8.3      
    ##  [5] cowplot_1.1.1      tidyseurat_0.5.1   ttservice_0.1.2    SeuratObject_4.0.4
    ##  [9] Seurat_4.1.0       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.8       
    ## [13] purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.6      
    ## [17] ggplot2_3.3.5      tidyverse_1.3.1   
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
    ##  [37] scales_1.1.1          pheatmap_1.0.12       DBI_1.1.2            
    ##  [40] spatstat.random_2.2-0 miniUI_0.1.1.1        xtable_1.8-4         
    ##  [43] reticulate_1.24       spatstat.core_2.4-2   htmlwidgets_1.5.4    
    ##  [46] httr_1.4.2            RColorBrewer_1.1-3    ellipsis_0.3.2       
    ##  [49] ica_1.0-2             farver_2.1.0          pkgconfig_2.0.3      
    ##  [52] sass_0.4.1            uwot_0.1.11           dbplyr_2.1.1         
    ##  [55] deldir_1.0-6          utf8_1.2.2            labeling_0.4.2       
    ##  [58] tidyselect_1.1.2      rlang_1.0.2           reshape2_1.4.4       
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
    ##  [91] RSpectra_0.16-0       lattice_0.20-45       Matrix_1.4-1         
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
