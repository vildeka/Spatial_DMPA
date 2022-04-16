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
        legend.position = "bottom",
        #aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )
```

## Integration

``` r
seuratObj.list <- SplitObject(seuratObj, split.by = "sample_name")

for (i in 1:length(seuratObj.list)) {
    seuratObj.list[[i]] <- NormalizeData(seuratObj.list[[i]], verbose = FALSE)
    seuratObj.list[[i]] <- FindVariableFeatures(seuratObj.list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}

hvgs_per_dataset <- lapply(seuratObj.list, function(x) {
    x@assays$RNA@var.features
})
# venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn
# = 1,cexil = 1,lwd=1,col='white',frame=F,borders = NA)

temp <- unique(unlist(hvgs_per_dataset))
overlap <- sapply(hvgs_per_dataset, function(x) {
    temp %in% x
})
rownames(overlap) <- temp
pheatmap::pheatmap(t(overlap * 1), cluster_rows = F, color = c("grey90", "grey20"))
```

<img src="./Figures/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

``` r
hig_var <- rownames(overlap)[rowSums(overlap)>2]
seuratObj <- seuratObj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                      nfeatures = 4000,
                      verbose = FALSE) %>%
  ScaleData(verbose = FALSE, features = hig_var ) %>%
  RunPCA(verbose = FALSE, npcs = 50) %>%
  RunUMAP(dims = 1:50, n.components = 2L) 

seuratObj <- RunHarmony(seuratObj, group.by.vars = "sample_name", reduction = "pca",
    dims.use = 1:50, assay.use = "RNA")

seuratObj <- RunUMAP(seuratObj, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony")
```

## Plot before and after integration

``` r
# fig.height=10,
seuratObj_red <- readRDS(paste0(result_dir,"seurat_object_red.RDS"))
#seuratObj_int_sct <- readRDS(paste0(result_dir,"seurat_object_int_sct_red.RDS"))
seuratObj_int_log <- readRDS(paste0(result_dir,"seurat_object_int_log_red.RDS"))

library(cowplot)
plot_grid(ncol = 2,
  DimPlot(seuratObj, reduction = "pca", group.by = "sample_name")+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(seuratObj, reduction = "umap", group.by = "sample_name")+NoAxes()+ggtitle("UMAP raw_data"),
  
  DimPlot(seuratObj, reduction = "harmony", group.by = "sample_name")+NoAxes()+ggtitle("harmony integrated"),
  DimPlot(seuratObj, reduction = "umap_harmony", group.by = "sample_name")+NoAxes()+ggtitle("UMAP harmony integrated")
  
  # DimPlot(seuratObj_int_log, reduction = "pca", group.by = "sample_name")+NoAxes()+ggtitle("PCA integrated"),
  # DimPlot(seuratObj_int_log, reduction = "tsne", group.by = "sample_name")+NoAxes()+ggtitle("tSNE integrated"),
  # DimPlot(seuratObj_int_log, reduction = "umap", group.by = "sample_name")+NoAxes()+ggtitle("UMAP integrated")
)
```

<img src="./Figures/Plot_dimreduction-1.png" style="display: block; margin: auto;" />

## Dimentionality reduction

## Plot before and after integration

``` r
##################################
# EVALUATE CLUSTERING RESOLUTION #
##################################
seuratObj <- FindNeighbors(seuratObj, reduction = "harmony", dims = 1:20, k.param = 20, prune.SNN = 1/15) 

for (res in c(0.1, 0.25, 0.5, 1, 1.5, 2)) {
    seuratObj <- FindClusters(seuratObj, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20696
    ## Number of edges: 710503
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9669
    ## Number of communities: 10
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20696
    ## Number of edges: 710503
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9352
    ## Number of communities: 13
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20696
    ## Number of edges: 710503
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8911
    ## Number of communities: 16
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20696
    ## Number of edges: 710503
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8376
    ## Number of communities: 22
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20696
    ## Number of edges: 710503
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8066
    ## Number of communities: 28
    ## Elapsed time: 4 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 20696
    ## Number of edges: 710503
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7817
    ## Number of communities: 35
    ## Elapsed time: 3 seconds

``` r
library(niceRplots)
par(mfrow = c(1,2))
plot_meta(seuratObj, red = "umap_harmony", feat = "RNA_snn_res.0.25", cex = 0.5)
plot_meta(seuratObj, red = "umap_harmony", feat = "RNA_snn_res.1", cex = 0.5, label = T)
```

<img src="./Figures/Clustering_resolution-1.png" style="display: block; margin: auto;" />

``` r
# plot_grid(ncol = 3, DimPlot(seuratObj, reduction = "umap", group.by = "RNA_snn_res.0.5") +
#     ggtitle("leiden_0.5"), DimPlot(seuratObj, reduction = "umap", group.by = "RNA_snn_res.1") +
#     ggtitle("leiden_1"), DimPlot(seuratObj, reduction = "umap", group.by = "RNA_snn_res.2") +
#     ggtitle("leiden_2"))
```
