Differential gene expression
================
3/21/23

### Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(tidyseurat)
library(cowplot)
library(ggrepel)
library(niceRplots)
library(MAST)
library(scran)
library(openxlsx)

source("../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../results/04_deconvolution_st_data/"
result_dir <- "../results/06_DGE_condition_st_data/"
marker_dir <- "./marker_genes_condition/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
if( isFALSE(dir.exists(marker_dir)) ) { dir.create(marker_dir,recursive = TRUE) }
PlosPath_DGEs <- "/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/data/Bulk_data/DMPA_PlosPath_DEGs.csv"

#############
# LODA DATA #
#############
DATA <- readRDS(paste0(input_dir,"seuratObj_deconvolution_scdc.RDS"))
# DATA <- readRDS(paste0("../results/03_clustering_st_data/","seuratObj_clustered.RDS"))
DGEs_PP <- read_csv(PlosPath_DGEs)
epi_clus <- "^5$|^6$|^7|^9" # non-filt

#################
# COLOUR PALLET #
#################
clus <- c( "#F8766D","#7CAE00", "#CD9600", "#00A9FF", "#377EB8","#984EA3", "#E41A1C", "#C77CFF",
             "#00BFC4", "#FF7F00","#FFFF33", "#FF61CC", "#4DAF4A",  "#A65628", "#F781BF", "#999999")
```

``` r
########################
# SEURAT DGEs FUNCTION #
########################
# obj <- DATA_sub$data[[1]]
# condition <- 'groups'
# subgroup <- "Clusters"
DEGs_fun <- function(obj, condition, subgroup){
  obj <- SetIdent(obj, value = condition)
  
  DATA_degs <- obj %>%
  mutate(clus = .data[[subgroup]]) %>%
  nest(data = -c("clus")) %>%
  mutate( DEGs = map(data, 
        ~FindAllMarkers(.x,
                       test.use = "wilcox",
                       only.pos = F,
                       return.thresh = 1,
                       max.cells.per.ident = Inf,
                       logfc.threshold = -Inf,
                       assay = "RNA",
                       min.pct = -Inf)
                    ))

  DEGs_table <- DATA_degs %>%
    mutate(DEGs = setNames( .[["DEGs"]], .$clus)) %>%
    .$DEGs %>%
    map(., ~as_tibble(.x)) %>%
    #map2(., comb$comb, ~mutate(.x, Combination = .y)) %>%
    bind_rows(., .id = "subgroup") %>%
    mutate(pct.diff = -.$pct.2-.$pct.1) %>%
    mutate(log.pct.diff = -log2(.$pct.2/.$pct.1))
  return(DEGs_table)
}
```

### Spot distribution after subsetting

``` r
####################################
# SUBSET SEURAT OBJECT PER CLUSTER #
####################################
# creating a subseted object with 25 spots per sampleID for each cluster
DATA_sub <- DATA %>%
  mutate(gr = .$groups) %>%
  mutate(ID = .$orig.ident) %>%
  nest(., data=-c(gr, orig.ident)) %>%
  mutate(epi =  map(data, ~filter(.x, !(sp_annot == "SubMuc"))),
         subMuc =  map(data, ~filter(.x, sp_annot == "SubMuc"))) %>%
  select(-data) %>%
  mutate(across(c("epi", "subMuc"), ~map(., ~table(.x$Clusters)), .names = "{.col}_n_before")) %>%
  mutate(epi = map(epi, ~filter(.x, .cell %in% WhichCells(., downsample = 25)))) %>%
  mutate(subMuc = map(subMuc, ~filter(.x, .cell %in% WhichCells(., downsample = 50)))) %>%
  mutate(across(c("epi", "subMuc"), ~map(., ~table(.x$Clusters)), .names = "{.col}_n_after")) %>%
  mutate(across(contains("_n_"), ~set_names(.x, paste0(.data[["gr"]],"_",.data[["orig.ident"]]))))

bind_cols(DATA_sub$epi_n_after, "Clus" = paste0("**",names(table(DATA$Clusters)),"**")) %>%
  rowwise() %>% 
  mutate(DMPA_sum = sum(c_across(starts_with("DMPA_"))),
         ctrl_sum = sum(c_across(starts_with("ctrl_")))) %>%
  select(sort(colnames(.)[1:8]), everything()) %>%
  knitr::kable(., caption = "Distribution of epithelial spots per cluster per subject")
```

| ctrl_P031 | ctrl_P080 | ctrl_P105 | ctrl_P118 | DMPA_P097 | DMPA_P107 | DMPA_P108 | DMPA_P114 | Clus   | DMPA_sum | ctrl_sum |
|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|:-------|---------:|---------:|
|         0 |         0 |         1 |         0 |         0 |         1 |         0 |         0 | **0**  |        1 |        1 |
|         1 |         0 |         0 |         0 |         0 |         1 |         0 |         1 | **1**  |        2 |        1 |
|         0 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | **2**  |        0 |        0 |
|         3 |         0 |         3 |         0 |         0 |         0 |         0 |         1 | **3**  |        1 |        6 |
|         0 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | **4**  |        0 |        0 |
|        25 |        25 |        25 |        25 |        25 |        25 |        25 |        25 | **5**  |      100 |      100 |
|        25 |        25 |        25 |        25 |        25 |        25 |        25 |        25 | **6**  |      100 |      100 |
|        25 |        25 |        25 |        25 |        25 |        25 |        25 |        25 | **7**  |      100 |      100 |
|        12 |         9 |         1 |         6 |         4 |         7 |         3 |         1 | **8**  |       15 |       28 |
|        25 |        25 |        25 |        25 |        25 |        25 |        25 |        25 | **9**  |      100 |      100 |
|         0 |         0 |         1 |         0 |         0 |         0 |         0 |         1 | **10** |        1 |        1 |

Distribution of epithelial spots per cluster per subject

``` r
bind_cols(DATA_sub$subMuc_n_after, "Clus" = paste0("**",names(table(DATA$Clusters)),"**")) %>%
  rowwise() %>% 
  mutate(DMPA_sum = sum(c_across(starts_with("DMPA_"))),
         ctrl_sum = sum(c_across(starts_with("ctrl_")))) %>%
  select(sort(colnames(.)[1:8]), everything()) %>%
  knitr::kable(., caption = "Distribution of submucosal spots per cluster per subject")
```

| ctrl_P031 | ctrl_P080 | ctrl_P105 | ctrl_P118 | DMPA_P097 | DMPA_P107 | DMPA_P108 | DMPA_P114 | Clus   | DMPA_sum | ctrl_sum |
|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|:-------|---------:|---------:|
|        50 |        50 |        50 |        50 |        50 |        20 |        34 |        50 | **0**  |      154 |      200 |
|        50 |        50 |        50 |        50 |        50 |        10 |        36 |        50 | **1**  |      146 |      200 |
|        37 |        50 |        50 |        50 |        50 |        35 |         9 |        50 | **2**  |      144 |      187 |
|        46 |        50 |        50 |        50 |        50 |        25 |        43 |        50 | **3**  |      168 |      196 |
|        24 |        50 |        50 |        50 |        50 |        37 |        21 |        50 | **4**  |      158 |      174 |
|         1 |         5 |         8 |        20 |        41 |         1 |         1 |         9 | **5**  |       52 |       34 |
|         0 |         3 |         0 |         0 |         4 |        19 |         0 |         0 | **6**  |       23 |        3 |
|         1 |        33 |         0 |         2 |         4 |         0 |         0 |         1 | **7**  |        5 |       36 |
|        24 |        46 |        50 |        50 |        50 |        14 |        20 |        50 | **8**  |      134 |      170 |
|         0 |         0 |         0 |         0 |         1 |         0 |         0 |         0 | **9**  |        1 |        0 |
|         1 |         1 |         0 |         0 |         5 |        29 |         0 |         0 | **10** |       34 |        2 |

Distribution of submucosal spots per cluster per subject

### Run differential gene expression analysis

``` r
##############################################
# DGEs BETWEEN CONDITION FOR EACH SUBGROUPS #
##############################################
DEGs_table_epi <- DATA_sub %>%
  unnest(epi) %>%
  filter(grepl(epi_clus, .$Clusters)) %>%
  DEGs_fun(., "groups", "Clusters") %>%
  filter(cluster == "DMPA") 

DEGs_table_subMuc <- DATA_sub %>%
  unnest(subMuc) %>%
  filter(!(grepl(epi_clus, .$Clusters))) %>%
  DEGs_fun(., "groups", "Clusters") %>%
  filter(cluster == "DMPA") 

DEGs_table <- bind_rows(DEGs_table_epi, DEGs_table_subMuc)

###################
# ADD ANNOTATION #
##################
ord1 <- c("Sup_1","Sup_2","Basal_2","Basal_1","0","1","2","3","4","8","10")
ord2 <- c("6","9","7","5","0","1","2","3","4","8","10")
epi_layers <- set_names(ord1, ord2)

DEGs_table <- DEGs_table %>%
  mutate(subgroup = factor(.$subgroup, levels = ord2)) %>%
  mutate(layers = factor(epi_layers[as.character(.$subgroup)], levels = ord1)) 
```

### Save deconvolution objects

``` r
# saveRDS(list(DATA_sub=DATA_sub, DEGs_table=DEGs_table), paste0(result_dir, "DATA_sub.RDS"))
# write_csv(DEGs_table, paste0(result_dir, "DGEs_condition_wilcox", ".csv"))

DEGs_list <- DEGs_table %>% 
  mutate(morf = ifelse(grepl(epi_clus, .$subgroup), "epi", "SubMuc")) %>%
  nest(data=-c(morf)) %>%
  mutate(data = setNames( .[["data"]], .$morf)) %>%
  .$data %>%
  map(., ~split(.x, .x$layers)) 

imap(DEGs_list, ~write.xlsx(.x, keepNA=TRUE, na.string="NA", overwrite=TRUE,
           file=paste0(result_dir,"DGEs_condition_wilcox_",.y,"_25.xlsx")))
```

### Load allready saved data

``` r
# Load data
DATA_sub <- readRDS(paste0(result_dir, "DATA_sub.RDS"))
DEGs_table <- DATA_sub$DEGs_table %>% filter(cluster == "DMPA")
DATA_sub <- DATA_sub$DATA_sub 
```

### Volcano plots of differentially expressed genes per cluster

``` r
###########################
# VOLCANO PLOT EPITHELIUM #
###########################
DEGs_filt <- DEGs_table %>% 
  filter(grepl(epi_clus, .$subgroup)) %>%
  filter(p_val < 0.099) 

Volcano.fun_logFC(DEGs_filt, "layers", y.axis="p-value", 
                  up=c(.2, 0.05), down = c(-.2, 0.05)) # labeling: (logFC, p-value)
```

<img src="../Figures/06/06a_DEGs_volcano_plot_epi.png"
data-fig-align="center" />

``` r
###########################
# VOLCANO PLOT SUBMUCOSA #
###########################
DEGs_filt <- DEGs_table %>% 
  filter(!(grepl(epi_clus, .$subgroup))) %>%
  filter(p_val < 0.099) 

Volcano.fun_logFC(DEGs_filt, "layers", y.axis="p-value", 
                  up=c(.2, 0.05), down = c(-.2, 0.05)) # labeling: (logFC, p-value)
```

<img src="../Figures/06/06b_DEGs_volcano_plot_subMuc.png"
data-fig-align="center" />

### Compare DEGs in spatial data with DEGs on bulk data

``` r
# Significant genes from the PlosPath study
# Genes with FDR-adjustedP-values<0.05 were considered significant
PP_sig <- DGEs_PP %>%
  filter(DE != 0) %>%
  {. ->> sig_PP } %>%
  group_by(DE) %>%
  top_n(15, abs(logFC))

# top 15 genes sorted by logFC for each epithelial layer
sig_nest <- DEGs_table %>%
  mutate(Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
  filter(p_val_adj < 0.05) %>%
  {. ->> sig_table } %>%
  group_by(layers, Direction) %>%
  top_n(15, abs(avg_log2FC)) %>% 
  #top_n(-30, p_val) %>% 
  nest() %>% mutate(n = map_dbl(data, nrow)) %>%
  arrange(layers, Direction)

sig_nest_epi <- sig_nest %>% ungroup() %>% filter(grepl("_", .$layers))
sig_genes_epi <- unnest(sig_nest_epi, c(layers, Direction, data))

sig_nest
```

    # A tibble: 20 × 4
    # Groups:   layers, Direction [20]
       layers  Direction data                   n
       <fct>   <chr>     <list>             <dbl>
     1 Sup_1   DOWN      <tibble [12 × 10]>    12
     2 Sup_1   UP        <tibble [7 × 10]>      7
     3 Sup_2   DOWN      <tibble [6 × 10]>      6
     4 Sup_2   UP        <tibble [11 × 10]>    11
     5 Basal_2 DOWN      <tibble [8 × 10]>      8
     6 Basal_2 UP        <tibble [15 × 10]>    15
     7 Basal_1 DOWN      <tibble [5 × 10]>      5
     8 Basal_1 UP        <tibble [11 × 10]>    11
     9 0       DOWN      <tibble [8 × 10]>      8
    10 0       UP        <tibble [12 × 10]>    12
    11 1       DOWN      <tibble [3 × 10]>      3
    12 1       UP        <tibble [13 × 10]>    13
    13 2       DOWN      <tibble [12 × 10]>    12
    14 2       UP        <tibble [15 × 10]>    15
    15 3       DOWN      <tibble [12 × 10]>    12
    16 3       UP        <tibble [15 × 10]>    15
    17 4       DOWN      <tibble [10 × 10]>    10
    18 4       UP        <tibble [7 × 10]>      7
    19 8       DOWN      <tibble [15 × 10]>    15
    20 8       UP        <tibble [15 × 10]>    15

``` r
#knitr::kable(sig_nest)
```

### Check overlapping genes

``` r
######################
# CHECK GENE OVERLAP #
######################
intersect(PP_sig$symbol, sig_genes_epi$gene)
```

    [1] "SPRR2B" "SPRR2G" "OLFM4"  "TGM3"  

``` r
intersect(sig_PP$symbol, sig_genes_epi$gene)
```

     [1] "PKP1"    "SPRR2D"  "SPRR2A"  "SPRR2B"  "SPRR2E"  "SPRR2G"  "CSTA"   
     [8] "SAMD9"   "GPRC5A"  "KRT6C"   "KRT6A"   "OLFM4"   "KRT13"   "KRT19"  
    [15] "KRT14"   "CD177"   "KRTDAP"  "DMKN"    "CEACAM1" "TGM3"   

``` r
# common DGEs between bulk and individual spatial clusters
t <- sig_table %>% 
  group_by(layers) %>%
  nest() %>%
  ungroup() %>%
  mutate(genes = map( data, ~pluck(.x, "gene"))) %>% 
  mutate(genes = setNames(.[["genes"]], .$layers)) %>%
  arrange(layers)

# all significant (FDR <0.05) genes overlap for each epi layer 
overlap <- map(t$genes, ~intersect(sig_PP$symbol, .x) )

#######################
# PLOTS VENN DIAGRAM #
######################
library(nVennR)
d <- plotVenn(t$genes[1:4], 
         sNames = names(t$genes[1:4]),
         outFile=paste0(params$fig.path, "Venn_epi_sig.svg")
         )

# showSVG(d)
knitr::include_graphics(paste0(params$fig.path, "Venn_epi_sig.svg"))
```

<img src="../Figures/06/Venn_epi_sig.svg" data-fig-align="center" />

### Identify marker genes to seperate clusters

``` r
#######################################
# FILTER BY P-VAL logFC AND pct.diff #
######################################
# Identify the top genes that have a high difference in expression between the clusters
top20 <- DEGs_table %>%
  mutate(rank = sign(avg_log2FC) * log10(p_val)) %>%
  mutate(Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
  group_by(subgroup, layers, Direction) %>%
  top_n(-30, p_val_adj) %>%
  top_n(20, abs(avg_log2FC))
  
top20 %>% nest() %>% mutate(n = map_dbl(data, nrow)) %>% arrange(layers)
```

    # A tibble: 20 × 5
    # Groups:   subgroup, layers, Direction [20]
       subgroup layers  Direction data                   n
       <fct>    <fct>   <chr>     <list>             <dbl>
     1 6        Sup_1   DOWN      <tibble [20 × 10]>    20
     2 6        Sup_1   UP        <tibble [20 × 10]>    20
     3 9        Sup_2   UP        <tibble [20 × 10]>    20
     4 9        Sup_2   DOWN      <tibble [20 × 10]>    20
     5 7        Basal_2 UP        <tibble [20 × 10]>    20
     6 7        Basal_2 DOWN      <tibble [20 × 10]>    20
     7 5        Basal_1 UP        <tibble [20 × 10]>    20
     8 5        Basal_1 DOWN      <tibble [20 × 10]>    20
     9 0        0       DOWN      <tibble [20 × 10]>    20
    10 0        0       UP        <tibble [20 × 10]>    20
    11 1        1       DOWN      <tibble [20 × 10]>    20
    12 1        1       UP        <tibble [20 × 10]>    20
    13 2        2       UP        <tibble [20 × 10]>    20
    14 2        2       DOWN      <tibble [20 × 10]>    20
    15 3        3       UP        <tibble [20 × 10]>    20
    16 3        3       DOWN      <tibble [20 × 10]>    20
    17 4        4       UP        <tibble [20 × 10]>    20
    18 4        4       DOWN      <tibble [20 × 10]>    20
    19 8        8       UP        <tibble [20 × 10]>    20
    20 8        8       DOWN      <tibble [20 × 10]>    20

``` r
# DEGs ordered by rank (logFC*p-value)
toprank <- DEGs_table %>%
  mutate(rank = sign(avg_log2FC) * log10(p_val)) %>%
  arrange(rank) 

# Identify marker genes to separate condition
# epithelium:
genes_e <- top20 %>% ungroup %>% filter(grepl(epi_clus, .$subgroup)) %>% .$gene %>% unique()
# Submucosa:
genes_s <- top20 %>% ungroup %>% filter(!(grepl(epi_clus, .$subgroup))) %>% .$gene %>% unique()
```

### Dotplot of DEGs between ctr and DMPA for each cluster

They are selected by first choosing the top 30 p-value and then
selecting the top 20 logFC among them

``` r
# dev.new(width=8.27, height=length(genes)/6+2, noRStudioGD = TRUE)
#######################
# DOT PLOT EPITHELIUM #
#######################
DATA_top <-  DATA %>%
  filter(grepl(epi_clus, .$Clusters)) %>%
  mutate(., feature = paste0(.$layers, "_", .$groups))
ord <- niceRplots::getcluster(DATA_top, genes_e, "feature")

#pdf( paste0("./top_10_DEG_subset_epi.pdf"),width=8.27,height=length(ord)/6+2)
par(mfrow=c(1,1), mar=c(7,6,1,5))
plot_dots(DATA_top, names(sort(ord)), clustering = "feature", show_grid = T,main = "top cluster markers",cex.main=1,font.main=1,cex.col = 1,srt = 90,cex.row = 1.1)
abline(v=cumsum(c(table(sub( "_.*","",names(table(DATA_top$feature))))))+0.5)
```

<img src="../Figures/06/06d_dot_plott_top_DEGs_epi.png"
data-fig-align="center" />

``` r
#dev.off()
```

``` r
# dev.new(width=8.27, height=length(genes)/6+2, noRStudioGD = TRUE) 
#######################
# DOT PLOT SUBMUCOSA #
#######################
DATA_top <-  DATA %>%
  filter(!(grepl(epi_clus, .$Clusters))) %>%
  mutate(., feature = paste0(.$Clusters, "_", .$groups))
ord <- niceRplots::getcluster(DATA_top, genes_s, "feature")

#pdf( paste0("./top_10_DEG_subset_epi.pdf"),width=8.27,height=length(ord)/6+2)
par(mfrow=c(1,1), mar=c(7,6,1,5))
plot_dots(DATA_top, names(sort(ord)), clustering = "feature", show_grid = T,main = "top cluster markers",cex.main=1,font.main=1,cex.col = 1,srt = 90,cex.row = 1.1)
abline(v=cumsum(c(table(sub( "_.*","",names(table(DATA_top$feature))))))+0.5)
```

<img src="../Figures/06/06e_dot_plott_top_DEGs_subMuc.png"
data-fig-align="center" />

``` r
#dev.off()
```

### Violin plots of top DEGs

``` r
# dev.new(width=6.6929133858, height=(round((length(unique(toprank$gene[1:3])))/3))*3, noRStudioGD = TRUE) 
################
# VIOLIN PLOTS #
################
library(patchwork)
violin_split.fun <- function(obj, gene_list, split.by, group.by){
  plots <- VlnPlot(obj, features = gene_list, adjust = 1/2,
                   split.by = split.by, group.by = group.by, 
                   pt.size = 0, combine = FALSE, split.plot = TRUE) 
  plots <- map(plots, ~.x + geom_point(alpha = .3, size=.1, show.legend = F,
                                       position = position_jitterdodge(jitter.width=.5)) +
                 theme(axis.title.x = element_blank()))
  p <- wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect")
  h <- round(length(plots)/3)*3
return(tibble(plot=list(p), height = list(h)))
}
plots <- DATA %>%
  filter(grepl(epi_clus, .$Clusters)) %>%
  violin_split.fun(., unique(toprank$gene[1:3]), "groups", "layers")
#pmap(plots, ~ggsave(paste0(result_dir,"top_DEGs_violin_top_rank",".jpeg"), plot=..1, width = 10, height = ..2))
plots$plot
```

    [[1]]

<img src="../Figures/06/06f_split_violin_plots.png"
data-fig-align="center" />

### Barplot of top 10 down and up regulated DEGs

The top DEGs are selected by rank (abs(avg_log2FC) \* log10(p_val))

``` r
# dev.new(height=6, width=7, noRStudioGD = TRUE)
#####################
# BARPLOTS TOP DEGs #
#####################8
top <- top20 %>% 
  slice_min(., n = 10, order_by = p_val_adj, with_ties=F) #%>% nest() %>%
n_gr <- length(unique(top20$subgroup))

#pdf(paste0(result_dir,"top_DEG.pdf"), height=(n_gr/4)*5, width=n_gr*2.5)
par(mfrow=c(ceiling(n_gr/4), 4), mar = c(4, 6, 3, 1))
for (i in unique(top$layers)) {
    barplot(sort(setNames(top$avg_log2FC, top$gene)[top$layers == i], F),
        horiz = T, las = 1, main = paste0( i), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}
#dev.off()
```

<img src="../Figures/06/06g_barplot_top20.png"
data-fig-align="center" />

### Top 10 Marker genes for each cluster

Selected only by adjusted p-value

``` r
# dev.new(height=7, width=7, noRStudioGD = TRUE)
top <- top20 %>% 
  slice_min(., n = 10, order_by = p_val_adj, with_ties=F) %>%
  slice_max(., n = 1, order_by = abs(log.pct.diff), with_ties=F) #%>% nest()

# remove all VDJ-genes from list of HVG
remove <- str_subset(top$gene, "^IGH|^IGK|^IGL|^TRA|^TRB|^TRD|^TRG")
top_gr <- top %>%
  ungroup() %>%
  #filter(., !(.$gene %in% remove)) %>%
  group_by(subgroup) %>%
  group_split() %>%
  set_names(., unique(top$layers))

##################################
# UMAP OF TOP 5 DOWN AND UP DEGs #
##################################
# col=c("grey90","grey80","grey60","navy","black")
col <- c("lightgray", "mistyrose", "red", "dark red", "black")

grid_genes <- function(plot_list, title){
  title <- ggdraw() + draw_label(paste0("Top 10 Markers for Cluster ", title), fontface='bold')
  g <- plot_grid(plotlist = plot_list,
            ncol = 3)
  g_ <- plot_grid(title, g, ncol=1, rel_heights=c(0.1, 1))
  return(g_)
}

clus_plot <- plot_clusters.fun(DATA, red = "umap_harmony", 
                               cluster="Clusters", color = clus) + theme_void() + NoLegend()
plots <- imap(top_gr, "gene") %>%
  map(., map, ~plot_genes.fun(DATA, red = "umap_harmony", .x, lable = "Clusters", col = col, point_size = .2,))

cluster_markers <- plots %>%
  map(., ~c("Clusters"=list(clus_plot), .x)) %>%
  imap(., ~grid_genes(.x, .y ))

cluster_markers[[1]]

imap(cluster_markers, ~ggsave(paste0(result_dir,"Marker_genes_condition_", .y, ".jpg"), 
                              plot=.x, height = 8 ))
ggsave(paste0("./Figures/06/","06_marker_genes_condition_pct.diff", ".pdf"), 
       gridExtra::marrangeGrob(grobs=cluster_markers, ncol=1, nrow=1, height = 7))
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
    [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] patchwork_1.1.2             nVennR_0.2.3               
     [3] openxlsx_4.2.5.1            scran_1.22.1               
     [5] scuttle_1.4.0               MAST_1.20.0                
     [7] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
     [9] Biobase_2.54.0              GenomicRanges_1.46.1       
    [11] GenomeInfoDb_1.30.1         IRanges_2.28.0             
    [13] S4Vectors_0.32.4            BiocGenerics_0.40.0        
    [15] MatrixGenerics_1.6.0        matrixStats_0.63.0         
    [17] niceRplots_0.1.0            ggrepel_0.9.2              
    [19] cowplot_1.1.1               tidyseurat_0.5.3           
    [21] ttservice_0.2.2             SeuratObject_4.1.3         
    [23] Seurat_4.3.0                forcats_0.5.2              
    [25] stringr_1.5.0               dplyr_1.0.10               
    [27] purrr_1.0.1                 readr_2.1.3                
    [29] tidyr_1.2.1                 tibble_3.1.8               
    [31] ggplot2_3.4.0               tidyverse_1.3.2            

    loaded via a namespace (and not attached):
      [1] utf8_1.2.2                spatstat.explore_3.0-5   
      [3] reticulate_1.27           tidyselect_1.2.0         
      [5] htmlwidgets_1.6.1         grid_4.1.2               
      [7] BiocParallel_1.28.3       Rtsne_0.16               
      [9] munsell_0.5.0             ScaledMatrix_1.2.0       
     [11] codetools_0.2-18          ica_1.0-3                
     [13] statmod_1.5.0             future_1.30.0            
     [15] miniUI_0.1.1.1            withr_2.5.0              
     [17] spatstat.random_3.0-1     colorspace_2.0-3         
     [19] progressr_0.13.0          highr_0.10               
     [21] knitr_1.41                rstudioapi_0.14          
     [23] ROCR_1.0-11               tensor_1.5               
     [25] listenv_0.9.0             labeling_0.4.2           
     [27] GenomeInfoDbData_1.2.7    polyclip_1.10-4          
     [29] farver_2.1.1              bit64_4.0.5              
     [31] parallelly_1.33.0         vctrs_0.5.1              
     [33] generics_0.1.3            xfun_0.36                
     [35] timechange_0.2.0          R6_2.5.1                 
     [37] rsvd_1.0.5                locfit_1.5-9.7           
     [39] bitops_1.0-7              spatstat.utils_3.0-1     
     [41] DelayedArray_0.20.0       assertthat_0.2.1         
     [43] vroom_1.6.0               promises_1.2.0.1         
     [45] scales_1.2.1              googlesheets4_1.0.1      
     [47] gtable_0.3.1              beachmat_2.10.0          
     [49] globals_0.16.2            goftest_1.2-3            
     [51] rlang_1.0.6               splines_4.1.2            
     [53] lazyeval_0.2.2            gargle_1.2.1             
     [55] spatstat.geom_3.0-3       broom_1.0.2              
     [57] yaml_2.3.6                reshape2_1.4.4           
     [59] abind_1.4-5               modelr_0.1.10            
     [61] backports_1.4.1           httpuv_1.6.8             
     [63] tools_4.1.2               ellipsis_0.3.2           
     [65] RColorBrewer_1.1-3        ggridges_0.5.4           
     [67] Rcpp_1.0.9                plyr_1.8.8               
     [69] sparseMatrixStats_1.6.0   zlibbioc_1.40.0          
     [71] RCurl_1.98-1.9            deldir_1.0-6             
     [73] pbapply_1.6-0             zoo_1.8-11               
     [75] haven_2.5.1               cluster_2.1.4            
     [77] fs_1.5.2                  magrittr_2.0.3           
     [79] data.table_1.14.6         scattermore_0.8          
     [81] lmtest_0.9-40             reprex_2.0.2             
     [83] RANN_2.6.1                googledrive_2.0.0        
     [85] fitdistrplus_1.1-8        hms_1.1.2                
     [87] mime_0.12                 evaluate_0.19            
     [89] xtable_1.8-4              readxl_1.4.1             
     [91] gridExtra_2.3             compiler_4.1.2           
     [93] KernSmooth_2.23-20        crayon_1.5.2             
     [95] htmltools_0.5.4           later_1.3.0              
     [97] tzdb_0.3.0                lubridate_1.9.0          
     [99] DBI_1.1.3                 dbplyr_2.2.1             
    [101] MASS_7.3-58.1             Matrix_1.5-3             
    [103] cli_3.6.0                 metapod_1.2.0            
    [105] parallel_4.1.2            igraph_1.3.5             
    [107] pkgconfig_2.0.3           sp_1.5-1                 
    [109] plotly_4.10.1             spatstat.sparse_3.0-0    
    [111] xml2_1.3.3                dqrng_0.3.0              
    [113] XVector_0.34.0            rvest_1.0.3              
    [115] digest_0.6.31             sctransform_0.3.5        
    [117] RcppAnnoy_0.0.20          spatstat.data_3.0-0      
    [119] rmarkdown_2.20            cellranger_1.1.0         
    [121] leiden_0.4.3              edgeR_3.36.0             
    [123] uwot_0.1.14               DelayedMatrixStats_1.16.0
    [125] shiny_1.7.4               lifecycle_1.0.3          
    [127] nlme_3.1-161              jsonlite_1.8.4           
    [129] BiocNeighbors_1.12.0      limma_3.50.3             
    [131] viridisLite_0.4.1         fansi_1.0.3              
    [133] pillar_1.8.1              lattice_0.20-45          
    [135] fastmap_1.1.0             httr_1.4.4               
    [137] survival_3.5-0            glue_1.6.2               
    [139] zip_2.2.2                 png_0.1-8                
    [141] bit_4.0.5                 bluster_1.4.0            
    [143] stringi_1.7.12            BiocSingular_1.10.0      
    [145] irlba_2.3.5.1             future.apply_1.10.0      
