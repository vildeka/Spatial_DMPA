Figure 7
================
11/15/24

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
library(gridtext)
library(png)
library(grid)
library(scatterpie)
library(patchwork)
library(openxlsx)
library(readxl)

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

source("../../bin/spatial_visualization.R")
source("../../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../../results/08_spatial_dist/SM/"
input_dir <- "../../results/04_deconvolution_st_data/"
result_dir <- "../../results/09_figures/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
epi_clus <- "^5$|^6$|^7|^9" # non-filt

#############
# LOAD DATA #
#############
# DEGs_table <- read_csv(paste0(input_dir,"subset_100/DGEs_condition_wilcox_epi_100.csv")) %>% filter(groups == "DMPA")
DEGs_table <- read_csv(paste0("../../results/05_DGE_clusters_st_data/","DGEs_clusters_wilcox.csv"))
DATA <- readRDS(paste0(input_dir,"seuratObj_deconvolution_scdc.RDS"))
#DATA <- readRDS(paste0(input_dir,"seuratObj_spatial_dist_SM.RDS"))

# Zalenskya sig DEGs
Zal_DEGs <- read_xlsx("/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/data/Bulk_data/Zalenskaya_DEGs.xlsx", 
                         sheet = 'DMPA-vs-BL', col_names = TRUE, skip = 1 )
PlosPath_DGEs <- "/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/data/Bulk_data/DMPA_PlosPath_DEGs.csv"
DGEs_PP <- read_csv(PlosPath_DGEs)


DEGs_table <- read_csv(paste0("../../results/06_DGE_condition_st_data/","DGEs_condition_wilcox.csv"))
DEGs_table_E2lvl <- read_csv(paste0("../../results/06_DGE_condition_st_data/","DGEs_condition_E2lvl.csv")) 
```

``` r
# top 15 genes sorted by logFC for each epithelial layer
sig_nest <- DEGs_table %>%
  mutate(Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
  filter(p_val_adj < 0.05) %>%
  {. ->> sig_table } %>%
  group_by(layers, Direction) %>%
  #top_n(15, abs(avg_log2FC)) %>% 
  #top_n(-30, p_val) %>% 
  nest() %>% mutate(n = map_dbl(data, nrow)) %>%
  arrange(layers, Direction) %>%
  ungroup()

sig_nest_epi <- sig_nest %>% ungroup() %>% filter(grepl("_", .$layers))
sig_genes_epi <- unnest(sig_nest_epi, c(layers, Direction, data))
sig_genes <- unnest(sig_nest, c(layers, Direction, data))


# Significant genes from the PlosPath study
# Genes with FDR-adjusted P-values<0.05 were considered significant
PP_sig <- DGEs_PP %>%
  filter(DE != 0) %>%
  {. ->> sig_PP } %>%
  group_by(DE) %>%
  top_n(15, abs(logFC)) %>%
  arrange(logFC) %>%
  ungroup() %>%
  mutate(symbol = factor(.$symbol, levels=unique(.$symbol)))
```

### Check overlapping genes (bulk data)

``` r
######################
# CHECK GENE OVERLAP #
######################
E2 <- DEGs_table_E2lvl %>% filter(p_val_adj < 0.05) #%>% filter((grepl(epi_clus, .$subgroup)))
Z <- intersect(Zal_DEGs$symbol, sig_genes$gene)
B <- intersect(sig_PP$symbol, sig_genes$gene)
# intersect(Z, B)
# intersect(E2$gene, sig_genes$gene)
# intersect(Zal_DEGs$symbol, E2$gene)
# intersect(sig_PP$symbol, E2$gene)
# setdiff(sig_genes$gene, E2$gene)
# setdiff(E2$gene, sig_genes$gene)
# intersect(sig_genes$gene, E2$gene)
# intersect(E2$gene, sig_PP$symbol)

# Estrogen idenpendent genes
non_E2 <- setdiff(sig_genes$gene, E2$gene)

# Estrogen dependent genes
E2_dep <- setdiff(E2$gene, sig_genes$gene)

g <- intersect(E2_dep, Zal_DEGs$symbol)
# intersect(E2_dep, sig_PP$symbol)


#####################################
# TABLES WITH DEGs FROM ALL STUDIES #
#####################################
t <- list(Brad=sig_PP$symbol, Zal=Zal_DEGs$symbol, ST=sig_genes$gene, E2=E2$gene) %>%
  tibble(sig_DEGs=., study=names(.))

# Create table with DEGs from all clinical DMPA studies
# NB! Zalenskaya use fold change not -logFC
All_DEGs <- list(ST=sig_genes, Zal=Zal_DEGs, Brad=sig_PP, E2=E2) %>%
  map(., ~dplyr::select(.x, symbol=matches("symbol|gene"),
                            FDR=matches("p_val_adj|value"),
                            logFC=matches("FC|Fold"), cluster=matches("subgroup") )) %>%
  bind_rows(., .id = "study") %>%
  mutate(logFC = ifelse(.$study == "Zal"& .$logFC < 0, -log2(abs(.$logFC)),
                      ifelse(.$study == "Zal"& .$logFC > 0, log2(.$logFC), .$logFC)) )
write_csv(All_DEGs, "../../results/06_DGE_condition_st_data/DEGs_Zal_Brad_ST.csv")
All_DEGs <- read_csv("../../results/06_DGE_condition_st_data/DEGs_Zal_Brad_ST.csv")
```

``` r
#######################
# PLOTS VENN DIAGRAM #
######################
# "./Figures/07/"
library(nVennR)
# remotes::install_github("vqf/nVennR")
t <- list(Zal=Zal_DEGs$symbol, ST=unique(sig_genes$gene), E2=unique(E2$gene), Brad=sig_PP$symbol) %>%
  tibble(sig_DEGs=., study=names(.))

getVennOverlap <- function(lsvenn) {
  ItemsList <- gplots::venn(lsvenn, show.plot = FALSE)
  print(lengths(attributes(ItemsList)$intersections))
  t <- attributes(ItemsList)$intersections %>%
    tibble(name = names(.), "Overlapping genes" = .) %>%
    mutate(n = map_dbl(.$`Overlapping genes`, ~length(.x)), .after="name")
  return(t)
}
col <- set_names(c('#ca0020','#f4a582','#0571b0',"#018571"), t$study)
# all four comparisons
# plotVenn(t$sig_DEGs, 
#          labelRegions=FALSE,
#          #sNames = names(t$sig_DEGs), 
#          systemShow=F) %>%
#   {. ->> d} %>%
#   showSVG(., 
#          setColors = col[names(d$orig)],
#          labelRegions=FALSE,
#          showNumbers=FALSE,
#          fontScale = 1,
#          systemShow=T,
#          outFile=paste0("./Figures/07/", "Venn_E2_all.svg")
#          )
# int_genes <- getVennOverlap(t$sig_DEGs) 

# Zalenskaya, ST and E2
d <- plotVenn(t$sig_DEGs[1:3], 
         labelRegions=FALSE,
         sNames = names(t$sig_DEGs), 
         systemShow=F) 

d <- showSVG(d, 
         setColors = col[names(d$orig)],
         #labelRegions=FALSE,
         #showNumbers=FALSE,
         fontScale = 2,
         systemShow=F,
         outFile=paste0("./Figures/07/", "Venn_E2_Zal.svg")
         )

int_genes <- getVennOverlap(t$sig_DEGs[1:3]) 
```

           E2        ST       Zal     ST:E2    Zal:E2    Zal:ST Zal:ST:E2 
          589        70       163        61        49         6        17 

``` r
write.xlsx(int_genes, paste0(result_dir, "Venn_table_genes.xlsx"))
knitr::include_graphics(paste0("./Figures/07/", "Venn_E2_Zal.svg"))
```

<img src="../Figures/07/Venn_E2_Zal.svg" data-fig-align="center" />

``` r
# get genes for the overlaps:
# getVennRegion(d, c("E2", "Brad"))

#plot(d)

#p <- grImport2::readPicture(paste0("./Figures/07/", "Venn_E2_Zal.svg"))
```

``` r
#### Tassos r object ####
umap_obj_Tassos <- readRDS("/Users/vilkal/work/Brolidens_work/Projects/DMPA/data/umap_obj_Tassos.RDS")

ptc_br <- c('#FB5273', '#4FCEEF')

#### meta data ####
meta_dt <- umap_obj_Tassos %>% 
  mutate(group = ifelse(.$group == "non", "Ctrl", .$group)) %>%
  #df %>% select(ID) %>% unique() %>%
  mutate(st_group = case_match(ID,
        c("P114", "P107") ~ "DMPA high",
        c("P097", "P108") ~ "DMPA low",
        c("P031", "P105","P118", "P080") ~ "ST ctrl",
        .default = group)) %>%
  mutate(txt = ifelse(grepl("high|low|ST", .$st_group), .$ID, NA)) %>%
  mutate(st_group = factor(.$st_group, levels= c("DMPA", "Ctrl", "DMPA high", "DMPA low", "ST ctrl" )))

#### use Tassos UMAP ####
cpm1_3_umap_dt <- meta_dt %>% 
  #left_join(., select(meta_dt, ID, txt, st_group="group"), by="ID") %>% 
  dplyr::rename(UMAP1 = "X1", UMAP2 = "X2") %>%
  add_row(., UMAP1 = 0.370858994263568, UMAP2 = -1.15998967007587, txt="P118",
    ML = "ML8310", ID = "P118", group = "Ctrl", st_group = "ST ctrl") 
  #filter(!(.$ID == "P118")) %>%
  #filter(!(is.na(.$txt))) 

#### UMAP plotting ####
txt_df <- cpm1_3_umap_dt[!(is.na(cpm1_3_umap_dt$txt)),] %>%
  mutate(UMAP2 = UMAP2+c(0, .1, rep(0,6)))

(C <- ggplot(cpm1_3_umap_dt, aes(x=UMAP1, y=UMAP2, fill=st_group))+ 
  #geom_jitter( shape=21, size=3, color="white", width=.5, height=.5) +  # Tassos used jitter
  geom_point( shape=21, size=2, color="white", alpha = .7) +  
  scale_fill_manual(values=c('#4FCEEF','#FB5273', "#6B51A3","#9D9AC8", "#e68633"), name="Group")+ 
  guides(fill = guide_legend( keywidth = .1, keyheight = 1) ) + #, override.aes = list(size=10)
  geom_text(data=txt_df, aes(x=UMAP1, y=UMAP2, label=txt), size=3,
            hjust = 0, nudge_x = 0.07, color="gray51") +
  theme_bw()+ # base_size=14
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 8),
        legend.text = element_text(size = 9),
        legend.spacing.x = unit(4, 'pt'),
        legend.justification = "right",
        legend.margin=margin(0,2,-7,2), # moves the legend box
        legend.box.margin=margin(-5,1,-3,-30), # moves the legend
        line=element_blank(),
        axis.text = element_blank(),
        axis.title.x.bottom = element_text(margin = margin(t = 1)),
        axis.title.y.left = element_text(margin = margin(r = 0)),
        axis.ticks=element_blank(), 
        aspect.ratio=.7, 
        )   )
```

<img src="../Figures/07/07c_bulk_UMAP.png" data-fig-align="center" />

``` r
# dev.new(height=16.11, width=16.11, units="cm")
# ggsave(paste0("./Figures/07/", "Bulk_UMAP.png"), p)
```

``` r
library(ComplexHeatmap)
library(circlize)

#### GENES THAT ARE UNIQELY OVERLAPING ####
getVennOverlap <- function(lsvenn) {
  
  ItemsList <- gplots::venn(lsvenn, show.plot = FALSE)
  print(lengths(attributes(ItemsList)$intersections))
  return(attributes(ItemsList)$intersections)
}
################
# SELECT GENES #
################
int_genes <- t$sig_DEGs[2:4] %>%
  getVennOverlap() %>% 
  tibble(name = names(.), Uniqe_int = .) %>%
  mutate(n = map_dbl(.$Uniqe_int, ~length(.x)), .after="name")
```

    ST:E2:Brad         E2         ST       Brad      ST:E2    E2:Brad    ST:Brad 
            28        479         61       2115         50        159         15 

``` r
E2_split <- DEGs_table_E2lvl %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(morf = ifelse(grepl(epi_clus, .$subgroup), "epi", "submuc")) %>%
  split(., ~morf) 


E2_lvl_ord <- c("P114", "P107", "P108", "P097", "P118", "P080", "P105","P031")
# lvl_layers <- c("Sup_1","Sup_2","Basal_2","Basal_1","8","3","0","4","1","2" ,"10")
lvl_layers <- c("Superficial","Upper IM","Lower IM","Basal","8","3","0","4","1","2","10")

################
# PLOT HEATMAP #
################
heatmap.fun <- function(genes){
meta_df <- DATA %>% as_tibble() %>% dplyr::select(1:4, layers) 
###############
# GET MATRIX #
###############
# expression ctrl vs DMPA
tbl <- DATA %>% 
  FetchData(vars = genes ) %>% 
  cbind(., "id"=meta_df$orig.ident,"groups"=meta_df$groups, "layers"=meta_df$layers) %>%
  mutate(id=factor(.$id, levels=E2_lvl_ord)) %>%
  {if(morf=="epi") filter(., !(grepl("\\d+", .$layers))) else .} %>%
  {if(morf=="subMuc") filter(., grepl("8|3|^0", .$layers)) else .} %>% # select epithelium only |8|4|^0
  #filter(., grepl("8|3|^0", .$layers)) %>% 
  as_tibble(., rownames = "cell_id") 

# select column grouping
matx <- tbl %>%
  # group_by(groups, layers) %>% # summed by group + layer
  mutate(id = factor(.$id, levels = E2_lvl_ord)) %>%
  group_by(id, groups, layers) %>% # summed by id + group +layer
  summarise(across(-cell_id, ~sum(.)), .groups = "keep") %>%
  mutate(gr = pmap_chr(cur_group(), paste, sep = "_"), .before = "id") %>%
  #mutate(groups = paste0(.$groups,"_", .$layers)) %>%
  arrange(id) %>%
  arrange(groups) %>%
  arrange(layers) %>%
  ungroup() %>%
  {. ->> meta_df} %>%
  dplyr::select(gr, any_of(genes)) %>%
  column_to_rownames(var = "gr") %>%
  as.matrix() %>% scale() %>% t()

##############
# ANNOTATION #
##############
annot_row <- E2 %>%
  select(gene, layers) %>%
  #filter(layers %in% tbl$layers) %>%
  mutate(row = .$gene, rown= 1) %>%
  #mutate(Z = ifelse(.$genes %in% g_Z, "Z", NA)) %>%
  mutate(B = ifelse(.$gene %in% sig_PP$symbol, "B", NA)) %>%
  mutate(ST = ifelse(.$gene %in% sig_genes$gene, "ST", NA)) %>%
  pivot_wider(., names_from = "layers", values_from = "rown", values_fill = NA) %>%
  column_to_rownames(var = "row") %>%
  .[genes, ] # arranges the rows in same order as matx

left_anno_row <- rowAnnotation(df = annot_row %>% select(all_of(sort(unique(tbl$layers)))),
  simple_anno_size = unit(.1, "cm"),
  na_col = "white",
  show_legend = FALSE, show_annotation_name = FALSE,
  col = list("Superficial" = c("1"="#E41A1C"),"Upper IM" = c("1"="#FF7F00"),"Lower IM" = c("1"="#C77CFF"),"Basal" = c("1"="#984EA3"),
            "8" = c("1"="#00BFC4"),"3" = c("1"="#00A9FF"),"0" = c("1"="#CD9600")) )

right_anno_row <- rowAnnotation(
  show_annotation_name = F,
  simple_anno_size = unit(.15, "cm"),
  "B" = annot_row$B,
  "ST" = annot_row$ST,
  na_col = "white",
  show_legend = FALSE,
  col = list("ST" = c("ST"="#FFDACC"), "B" = c("B"="#99CEC6")) )

annot_col <- meta_df %>%
  dplyr::select(1:4) %>%
  filter(layers %in% tbl$layers) %>%
  mutate(gr = paste0(.$layers, "_",.$groups)) %>%
  mutate(E2_lvl = ifelse(grepl("P114|P107", .$id), "high", 
                         ifelse(grepl("P097|P108", .$id), "low", "ctrl"))) %>%
  mutate(E2_lvl = ifelse(grepl("P031|P105", .$id), "ctrl_high", .$E2_lvl))

column_labels = structure(str_extract(colnames(matx), "P\\d\\d\\d"), names = colnames(matx))

col_E2 <- c("ctrl_high"="#FF7F00", "ctrl"="#FED9A6", "high"="#6A51A3", "low"="#9E9AC8" )
lgd = Legend(
    title = "E2", 
    title_gp = gpar(fontsize = 7),
    labels = gt_render(c("ctrl_high", "ctrl", "high", "low")),
    legend_gp = gpar(fill =col_E2)
)

################
# DRAW HEATMAP #
################

quantile(matx, c(0.1, 0.95))
# average expression:
set.seed(123)
H <- Heatmap(matx, name = " ",
             col = circlize::colorRamp2(c(min(matx), 0, 4), c("#440154FF", "#21908CFF", "#FDE725FF")),
             #row_km = 6, #column_km = 2, # kmeans change every time you run it
             column_split =  factor(as.character(annot_col$layers), levels = unique(annot_col$layers)),
             row_split = 5,  # hierarchical static
             show_row_dend = FALSE,
             cluster_columns = F,
             
             # text
             row_title = gt_render("", padding = unit(c(0, 0, 0, 0), "pt")),
             
             column_labels = gt_render(column_labels, padding = unit(c(0, 0, 0, 0), "pt")),
             column_title_gp = grid::gpar(fontsize = 8),
             #column_labels_gp =  grid::gpar(fontsize = 8),
             
             column_names_gp = grid::gpar(fontsize = 6),
             row_names_gp = grid::gpar(fontsize = 6),
            
             
             # annotation
             right_annotation = right_anno_row, left_annotation = left_anno_row,
             top_annotation = 
               columnAnnotation(E2=annot_col$E2_lvl, 
                                #show_legend = FALSE,
                                show_annotation_name = F,
                                annotation_legend_param = list(grid_height = unit(.1, "mm"), grid_width = unit(2, "mm"), title = "",
                                                               labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8)),
                                simple_anno_size = unit(.1, "cm"),
                                #gap = unit(1, "cm"),
                                col=list( E2= set_names(c("#6A51A3","#9E9AC8","#FED9A6","#FF7F00"), unique(annot_col$E2_lvl)) ) ),
             
             # legend
             heatmap_legend_param = list(legend_height = unit(20, "mm"), grid_width = unit(2, "mm"), 
                                         labels_gp = gpar(fontsize = 7))
             #annotation_legend_param = list(size = unit(2, "mm"))

             ) #colorRampPalette(c(col))(10) 

# set.seed(1) Tassos seed

return(H)
}
# Heatmap global options:
ht_opt$COLUMN_ANNO_PADDING = unit(.05, "cm")
ht_opt$HEATMAP_LEGEND_PADDING = unit(-0, "cm") #unit(c(0, 0, 0, -1), "cm") #b,l,t,r
ht_opt$TITLE_PADDING = unit(.05, "cm")
ht_opt$DIMNAME_PADDING = unit(.05, "cm")


# Epi
genes <- intersect(E2_split$epi$gene, Zal_DEGs$symbol)
morf <- "epi"
H_epi <- heatmap.fun(genes)
H_epi <- grid.grabExpr(draw(H_epi, show_heatmap_legend = FALSE, merge_legend = TRUE))

# SubMuc
genes <- intersect(E2_split$submuc$gene, Zal_DEGs$symbol)
morf <- "subMuc"
H_sub <- heatmap.fun(genes)
H_sub <- grid.grabExpr(draw(H_sub, heatmap_legend_side = "left", merge_legend = F, background = NA, 
             show_annotation_legend = FALSE, padding =  unit(c(3, 9, -2, 3), "mm"))) # b,l,t,r


################
# COMBINE PLOTS #
################
# dev.new(width=3.5, height=7.48, noRStudioGD = TRUE)
(B <- plot_grid(H_epi, H_sub, labels = c('B'), ncol = 1, rel_heights = c(1, .7)))
```

<img src="../Figures/07/07b_E2-lvl-gene-heatmap.png"
data-fig-align="center" />

``` r
################
# SAVE RESULTS #
################
# dev.new(width=8.5, height=8, noRStudioGD = TRUE,  res = 300)
# png(file=paste0("./Figures/07/", "E2lvl_heatmap_Sub_Zal.png"), 
#     units = "in", res = 300,
#     #width = 8.5, height = 10 # Epi
#     width = 8.5, height = 8 # Sub
# )
# H_epi 
# dev.off()
```

``` r
# get higest logFC values from top sig genes list
filter_hig_pval.fun <- function(df, n){
  top_u <- df %>%
  mutate(Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
  #mutate(Direction = ifelse(p_val_adj > 0.05, "NOT SIG.", .$Direction)) %>%
  #filter(Regulation != "NOT SIG.") %>%
  mutate(Direction = factor(.$Direction, levels = c("UP", "DOWN"))) %>%
  slice_max(n=n, order_by = tibble(abs(avg_log2FC), p_val), by="Direction", with_ties = F) %>%
  filter(n()==1 | n()>1 & p_val==min(p_val), .by="gene") # filter duplicate genes

  g <- set_names(top_u$gene, top_u$Direction) %>% .[!duplicated(.)]

  DEGs <- top_u %>%  
    filter(., gene %in% c(g[names(g)=="UP"][1:15], g[names(g)=="DOWN"][1:15]))
return(DEGs)
}

DEGs_E2lvl <- DEGs_table_E2lvl %>%
  filter(Morphology == "epi") %>%
  filter_hig_pval.fun(., n=30)

# top_df <- DEGs_E2high
# genes <- DEGs_E2high$gene
top_df <- DEGs_E2lvl
genes <- DEGs_E2lvl$gene

df <- DATA %>%
  mutate(E2_lvl = ifelse(grepl("P114|P107", .$orig.ident), "high", 
                         ifelse(grepl("P097|P108", .$orig.ident), "low", "ctrl"))) %>%
  mutate(E2_lvl = ifelse(grepl("P031|P105", .$orig.ident), "ctrl high", .$E2_lvl)) %>%
  mutate(., FetchData(., vars = c(genes), slot = "counts")) %>%
  filter(sp_annot == "epi") %>%
  as_tibble() %>%
  select(1:5, E2_lvl, layers, all_of(genes)) %>%
  pivot_longer(cols = any_of(genes), names_to = "gene", values_to = "values") %>%
  group_by( E2_lvl, gene) %>%
  summarize(sum_counts = sum(values), .groups="drop") %>%
  mutate("Sum counts(log10)" = log10(.$sum_counts)) %>%
  left_join(., select(top_df, Direction, gene), by="gene") %>%
  arrange(`Sum counts(log10)`) %>%
  mutate(gene = factor(.$gene, levels=unique(.$gene)))

col <- c("ctrl high"="#FF7F00", "ctrl"="#FED9A6", "high"="#6A51A3", "low"="#9E9AC8" )

(D <- ggplot2::ggplot(data=df, aes(x=gene, y=`Sum counts(log10)`)) +
  geom_point(aes(col=E2_lvl), size = 1) +
  scale_colour_manual(values = col) +
  facet_wrap("Direction", ncol = 2, strip.position = "top", scales = "free_x", shrink = F) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0,-9,5,1),units = "pt"),
        legend.margin = margin(-15,2,-8,2), # moves the legend box
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(4, 'pt'),
        axis.title = element_text(size=6),
        axis.text.x = element_text(size=6.5, angle=45, hjust=1, color="black"),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
  ) )
```

<img src="../Figures/07/07d_E2-lvl-top_DEGs.png"
data-fig-align="center" />

``` r
#############################
# COMBINE ALL FIGURE PANELS #
#############################
# dev.new(width=3.74, height=7.48, noRStudioGD = TRUE)
A_C_D <- plot_grid( NULL, C, D, ncol=1, rel_heights = c(.75,1,1), labels = c('A','C','D'),  hjust = -.2) 

# dev.new(width=7.48, height=7.48, noRStudioGD = TRUE)
Figure7 <- plot_grid( A_C_D, B, ncol=2, rel_widths = c(1,1)) 
# ggsave("./Figures/Figure-7.pdf", Figure7, width = 7.48, height = 7.48, bg = "white", dpi = 1000)

Figure7
```

<img src="../Figures/07/Figure-7.png" data-fig-align="center" />
