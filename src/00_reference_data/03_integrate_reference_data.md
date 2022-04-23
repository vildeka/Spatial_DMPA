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
input_dir <- "../../results/02_QC_reference_data/"
result_dir <- "../../results/03_integrate_reference_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
seuratObj <- readRDS(paste0(input_dir,"seurat_object_filtered.RDS"))

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

## Integration

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

pheatmap::pheatmap(t(hvgs_heat * 1), cluster_rows = F, color = c("grey90", "grey20"))
```

<img src="./Figures/Find_HVG-1.png" style="display: block; margin: auto;" />

``` r
# choose the hvg present in at least two samples:
hig_var <- rownames(hvgs_heat)[rowSums(hvgs_heat)>2]
```

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

saveRDS(seuratObj, paste0(result_dir,"seuratObj_harmony.RDS"))
#seuratObj <- readRDS(paste0(result_dir,"seuratObj_harmony.RDS"))

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

<img src="./Figures/Plot_dim_reduction-1.png" style="display: block; margin: auto;" />

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
  my_theme + theme_void() + theme(legend.position = "bottom")
  return(p)
}

p <- map(genes, ~plot_genes.fun(seuratObj, .x, maxs = 10))
plot_grid(ncol = 3, 
         plotlist = p )
```

<img src="./Figures/visualize_marker_genes-1.png" style="display: block; margin: auto;" />

# Paulos base R code

# Clustering
