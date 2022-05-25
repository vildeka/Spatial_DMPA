## Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(tidyseurat)
library(cowplot)

#########
# PATHS #
#########
input_dir <- "../../results/03_clustering_ref_data/"
result_dir <- "../../results/04_DGE_ref_data/"
marker_dir <- "./Marker_genes/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
if( isFALSE(dir.exists(marker_dir)) ) { dir.create(marker_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
seuratObj <- readRDS(paste0(input_dir,"seuratObj_clustered.RDS"))
```

``` r
#######################
# PLOT GENES FUNCTION #
#######################
# obj <- seuratObj
 # gene <- sym("CTSK")
plot_genes.fun <- function(obj, gene, mins=NULL, maxs=NULL, red = "umap_harmony", lable = TRUE){
  gene <- sym(gene)
  obj <- obj %>%
    mutate(lab = obj@active.ident) %>%
    mutate(., FetchData(., vars = c(as_label(gene))) ) %>%
    mutate(feat = !!(gene))
  feat_vec <- pull(obj, as_label(gene))
  
  # Colour pal:
  if(is.null(mins)){
    mins <- min(c(feat_vec, 0),na.rm = T)} # get 0 or negative value
  if(is.null(maxs)){maxs <- quantile(feat_vec,0.99,na.rm = T) # get percentile
    if(maxs==0){maxs <- max(feat_vec,na.rm = T)}
  }
  if(max(feat_vec, na.rm=T) != 0){
    # Calculate percentage:
    obj <- obj %>%
      mutate(feat = (!!(gene) - mins) / ( maxs - mins) ) %>%
      mutate(feat = ifelse(.$feat > 1, 1, .$feat))
  }
  
  obj <- obj %>%
      #select(1:3, !!(gene)) %>%
      mutate(feat = round(.$feat*98)+1) %>%
      mutate(pal = c( col[1],colorRampPalette(col[-1])(99))[.$feat] ) %>%
      arrange(!!(gene))
  
  # reduction method:
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  # txt lable:
  if(lable == FALSE){l = FALSE
    text <- NoLegend() #+ labs(color= "Clusters")
  }else{l = TRUE
    if(lable != TRUE){obj <- mutate(obj, lab = pull(obj, lable))}
    
    lable_df <- obj %>%
      group_by(lab) %>%
      select(lab, contains(red)) %>% 
      summarize_all(mean) 
    
    text <- geom_text(data = lable_df, aes(label = lab), col="black", size=2.5) }
  
  p <- ggplot(obj, aes(!!(red_1), !!(red_2), label=l , color = pal) ) +
    geom_point(alpha = 0.5, size=.5) + ggtitle(as_label(gene)) +
    text + #scale_color_viridis(option = "D", na.value="#EBECF0") +
    scale_colour_identity() +
    my_theme + theme_void() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) 
  return(p)
}

##########################
# PLOT CLUSTERS FUNCTION #
##########################
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

## Identify marker genes to seperate clusters

``` r
# Compute differential expression one clusters against rest 
markers_genes <- seuratObj %>%
  FindAllMarkers(.,
    log2FC.threshold = 0.1, 
    test.use = "wilcox",
    min.pct = 0.1, 
    min.diff.pct = 0.1, 
    only.pos = TRUE, 
    max.cells.per.ident = 50,
    assay = "RNA", slot = "data")


#markers_genes$pct.diff <- -markers_genes$pct.2-markers_genes$pct.1
#markers_genes$log.pct.diff <- -log2(markers_genes$pct.2/markers_genes$pct.1)
```

## Save seurat object

``` r
###########################
# SAVE INTERMEDIATE OJECT #
###########################
# saveRDS(markers_genes, paste0(result_dir,"markers_genes.RDS"))
# markers_genes <- readRDS(paste0(result_dir,"markers_genes.RDS"))
```

``` r
# Get difference in expression between one cluster and the rest
markers_genes <- markers_genes %>%
  mutate(pct.diff = -markers_genes$pct.2-markers_genes$pct.1) %>%
  mutate(log.pct.diff = -log2(markers_genes$pct.2/markers_genes$pct.1))

# Identify the top genes that have a high difference in expression
top25 <- markers_genes %>%
  group_by(cluster) %>%
  top_n(-40, p_val_adj) %>%
  top_n(20, pct.diff) %>%
  top_n(10, log.pct.diff) 

# remove all VDJ-genes from list of HVG
remove <- str_subset(top25$gene, "^IGH|^IGK|^IGL|^TRA|^TRB|^TRD|^TRG")
top25_gr <- top25 %>%
  ungroup() %>%
  filter(., !(.$gene %in% remove)) %>%
  group_by(cluster) %>%
  group_split() %>%
  set_names(., seq(0, length(.)-1) )
  

top25 <- top25 %>%
  arrange(cluster)
# top25[grep("NCR1",top25$gene),]
```

``` r
pdf(paste0(marker_dir,"top_DEG.pdf"), #width = 4*10*300, height = 30*4*300#, res = 300
    )
par(mfrow=c(2, 5), mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
    barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
        horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}
dev.off()
```

    ## pdf 
    ##   2

``` r
par(mfrow=c(2, 5), mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)[1:10]) {
    barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
        horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}
```

<img src="./Figures/04a_barplot_top25-1.png" style="display: block; margin: auto;" />

``` r
col=c("grey90","grey80","grey60","navy","black")

clus_plot <- plot_clusters.fun(seuratObj, cluster="RNA_snn_res.1.5") + theme_void() + NoLegend()

grid_genes <- function(plot_list, title){
  title <- ggdraw() + draw_label(paste0("Top 10 Markers for Cluster ", title), fontface='bold')
  g <- plot_grid(plotlist = plot_list,
            ncol = 4)
  g_ <- plot_grid(title, g, ncol=1, rel_heights=c(0.1, 1))
  return(g_)
}

plots <- imap(top25_gr, "gene") %>%
  map(., map, ~plot_genes.fun(seuratObj, .x, lable = "RNA_snn_res.1.5"))

cluster_markers <- plots %>%
  map(., ~c("Clusters"=list(clus_plot), .x)) %>%
  imap(., ~grid_genes(.x, .y ))

cluster_markers[[1]]
```

<img src="./Figures/graph_plot_top25-1.png" style="display: block; margin: auto;" />

``` r
imap(cluster_markers, ~ggsave(paste0(marker_dir,"Marker_genes_cluster_", .y, ".jpg"), plot=.x))
```

    ## $`0`
    ## [1] "./Marker_genes/Marker_genes_cluster_0.jpg"
    ## 
    ## $`1`
    ## [1] "./Marker_genes/Marker_genes_cluster_1.jpg"
    ## 
    ## $`2`
    ## [1] "./Marker_genes/Marker_genes_cluster_2.jpg"
    ## 
    ## $`3`
    ## [1] "./Marker_genes/Marker_genes_cluster_3.jpg"
    ## 
    ## $`4`
    ## [1] "./Marker_genes/Marker_genes_cluster_4.jpg"
    ## 
    ## $`5`
    ## [1] "./Marker_genes/Marker_genes_cluster_5.jpg"
    ## 
    ## $`6`
    ## [1] "./Marker_genes/Marker_genes_cluster_6.jpg"
    ## 
    ## $`7`
    ## [1] "./Marker_genes/Marker_genes_cluster_7.jpg"
    ## 
    ## $`8`
    ## [1] "./Marker_genes/Marker_genes_cluster_8.jpg"
    ## 
    ## $`9`
    ## [1] "./Marker_genes/Marker_genes_cluster_9.jpg"
    ## 
    ## $`10`
    ## [1] "./Marker_genes/Marker_genes_cluster_10.jpg"
    ## 
    ## $`11`
    ## [1] "./Marker_genes/Marker_genes_cluster_11.jpg"
    ## 
    ## $`12`
    ## [1] "./Marker_genes/Marker_genes_cluster_12.jpg"
    ## 
    ## $`13`
    ## [1] "./Marker_genes/Marker_genes_cluster_13.jpg"
    ## 
    ## $`14`
    ## [1] "./Marker_genes/Marker_genes_cluster_14.jpg"
    ## 
    ## $`15`
    ## [1] "./Marker_genes/Marker_genes_cluster_15.jpg"
    ## 
    ## $`16`
    ## [1] "./Marker_genes/Marker_genes_cluster_16.jpg"
    ## 
    ## $`17`
    ## [1] "./Marker_genes/Marker_genes_cluster_17.jpg"
    ## 
    ## $`18`
    ## [1] "./Marker_genes/Marker_genes_cluster_18.jpg"
    ## 
    ## $`19`
    ## [1] "./Marker_genes/Marker_genes_cluster_19.jpg"
    ## 
    ## $`20`
    ## [1] "./Marker_genes/Marker_genes_cluster_20.jpg"
    ## 
    ## $`21`
    ## [1] "./Marker_genes/Marker_genes_cluster_21.jpg"
    ## 
    ## $`22`
    ## [1] "./Marker_genes/Marker_genes_cluster_22.jpg"
    ## 
    ## $`23`
    ## [1] "./Marker_genes/Marker_genes_cluster_23.jpg"
    ## 
    ## $`24`
    ## [1] "./Marker_genes/Marker_genes_cluster_24.jpg"
    ## 
    ## $`25`
    ## [1] "./Marker_genes/Marker_genes_cluster_25.jpg"
    ## 
    ## $`26`
    ## [1] "./Marker_genes/Marker_genes_cluster_26.jpg"
    ## 
    ## $`27`
    ## [1] "./Marker_genes/Marker_genes_cluster_27.jpg"

``` r
# ggsave("multipage.pdf", gridExtra::marrangeGrob(grobs=cluster_markers, ncol=1, nrow=1))
```

# Paulos Code

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
    ##  [16] tzdb_0.2.0            limma_3.50.1          globals_0.14.0       
    ##  [19] modelr_0.1.8          matrixStats_0.61.0    spatstat.sparse_2.1-0
    ##  [22] colorspace_2.0-3      rvest_1.0.2           ggrepel_0.9.1        
    ##  [25] haven_2.4.3           xfun_0.30             crayon_1.5.1         
    ##  [28] jsonlite_1.8.0        spatstat.data_2.1-4   survival_3.2-13      
    ##  [31] zoo_1.8-9             glue_1.6.2            polyclip_1.10-0      
    ##  [34] gtable_0.3.0          leiden_0.3.9          future.apply_1.8.1   
    ##  [37] abind_1.4-5           scales_1.1.1          DBI_1.1.2            
    ##  [40] spatstat.random_2.2-0 miniUI_0.1.1.1        Rcpp_1.0.8.3         
    ##  [43] viridisLite_0.4.0     xtable_1.8-4          reticulate_1.24      
    ##  [46] spatstat.core_2.4-2   htmlwidgets_1.5.4     httr_1.4.2           
    ##  [49] RColorBrewer_1.1-3    ellipsis_0.3.2        ica_1.0-2            
    ##  [52] pkgconfig_2.0.3       farver_2.1.0          sass_0.4.1           
    ##  [55] uwot_0.1.11           dbplyr_2.1.1          deldir_1.0-6         
    ##  [58] utf8_1.2.2            labeling_0.4.2        tidyselect_1.1.2     
    ##  [61] rlang_1.0.2           reshape2_1.4.4        later_1.3.0          
    ##  [64] munsell_0.5.0         cellranger_1.1.0      tools_4.1.2          
    ##  [67] cli_3.2.0             generics_0.1.2        broom_0.7.12         
    ##  [70] ggridges_0.5.3        evaluate_0.15         fastmap_1.1.0        
    ##  [73] yaml_2.3.5            goftest_1.2-3         knitr_1.38           
    ##  [76] fs_1.5.2              fitdistrplus_1.1-8    RANN_2.6.1           
    ##  [79] pbapply_1.5-0         future_1.24.0         nlme_3.1-157         
    ##  [82] mime_0.12             xml2_1.3.3            compiler_4.1.2       
    ##  [85] rstudioapi_0.13       plotly_4.10.0         png_0.1-7            
    ##  [88] spatstat.utils_2.3-0  reprex_2.0.1          bslib_0.3.1          
    ##  [91] stringi_1.7.6         highr_0.9             lattice_0.20-45      
    ##  [94] Matrix_1.4-1          vctrs_0.4.0           pillar_1.7.0         
    ##  [97] lifecycle_1.0.1       spatstat.geom_2.4-0   lmtest_0.9-40        
    ## [100] jquerylib_0.1.4       RcppAnnoy_0.0.19      data.table_1.14.2    
    ## [103] irlba_2.3.5           httpuv_1.6.5          patchwork_1.1.1      
    ## [106] R6_2.5.1              promises_1.2.0.1      KernSmooth_2.23-20   
    ## [109] gridExtra_2.3         parallelly_1.31.0     codetools_0.2-18     
    ## [112] MASS_7.3-56           assertthat_0.2.1      withr_2.5.0          
    ## [115] sctransform_0.3.3     mgcv_1.8-40           parallel_4.1.2       
    ## [118] hms_1.1.1             grid_4.1.2            rpart_4.1.16         
    ## [121] rmarkdown_2.11        Rtsne_0.15            shiny_1.7.1          
    ## [124] lubridate_1.8.0
