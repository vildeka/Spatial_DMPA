Figure 1
================
4/12/24

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
library(png)
library(grid)
library(scatterpie)
library(patchwork)
library(ggnewscale)
library(openxlsx)

source("../../bin/spatial_visualization.R")
source("../../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../../results/04_deconvolution_st_data/"
result_dir <- "../../results/09_figures/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LOAD DATA #
#############
# DEGs_table <- read_csv(paste0(input_dir,"subset_100/DGEs_condition_wilcox_epi_100.csv")) %>% filter(groups == "DMPA")
DEGs_table <- read_csv(paste0("../../results/06_DGE_condition_st_data/","DGEs_condition_wilcox.csv"))
DEGs_table_morf <- read_csv(paste0("../../results/05_DGE_clusters_st_data/","DGEs_clusters_morf.csv"))
DATA <- readRDS(paste0(input_dir,"seuratObj_deconvolution_scdc.RDS"))
DATA_r <- readRDS(paste0("../../results/06_plot_annotation_ref_data/","seuratObj_annotated_celltypes.RDS"))

#img <- readPNG("/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/data/Spatial_data/P097/spatial/tissue_hires_image.png")
#img <- readPNG("/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/data/Spatial_data/P118/spatial/tissue_lowres_image.png")
#img <- readPNG("../../resources/Schematic figure ST.png")
#g <- rasterGrob(img, interpolate=TRUE) #, height = 1, width = 1
#g <- grid::rasterGrob(img_[[1]], interpolate=TRUE)

# g <- ggplot() +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   #xlim(260, 800)
#   #annotation_custom(g, xmin=1, xmax=10, ymin=1, ymax=10) +
#   #coord_fixed(ratio = nrow(img)/ncol(img)) +
#   theme(panel.background = element_rect(fill = "transparent"))
#     #plot.margin = unit(c(-2,0,-2,0), "cm"))
```

``` r
######################
# CLUSTERS ON TISSUE #
######################
colors <- c( "#CD9600", "#7CAE00", "#e0e067", "#00A9FF", "#377EB8","#984EA3", "#E41A1C", "#C77CFF",
             "#00BFC4", "#FF7F00","#FFFF33")
clus_n <- 0:10
colors <- set_names(colors, clus_n)

clus <- DATA %>%
  #mutate(Clusters = as.character(.$Clusters)) %>%
  plot_spatial.fun( .,
                 geneid = "Clusters",
                 sampleid = c("P118", "P097"),
                 save_space = F,
                 lab = F,
                 zoom = "zoom",
                 col = colors,
                 alpha = .9,
                 ncol = 1, 
                 #annot_col = "#dbd9d9",
                 annot_line = .1,
                 img_alpha = 0,
                 point_size = .6) + 
  theme(legend.position = "none") # moves the legend)
  
clus
```

<img src="../Figures/01/01c_clusters_on_tissue.png"
data-fig-align="center" />

``` r
# dev.new(height=3.3, width=2, noRStudioGD = TRUE)
# ggsave("./Figures/S1/clusters_on_tissue.pdf", clus, width = 7, height = 3.3, dpi = 500, bg = "white")
# ggsave("./Figures/01/clusters_on_tissue.png", clus, width = 2, height = 3.3, dpi = 1000)
```

``` r
#########################
# UMAP CLUSTERS ST DATA #
#########################
clus_col <- c("#E41A1C","#FF7F00", "#C77CFF","#984EA3","#00BFC4", "#00A9FF","#377EB8","#CD9600","#7CAE00", "#e0e067","#FFFF33")
clus_lvl <- c("6", "9", "7", "5","8","3","4","0","1","2","10")
clus_col <- set_names(clus_col, clus_lvl)

UMAP_ST <- plot_clusters.fun(DATA, cluster="Clusters", txt_size = 10, dot_size = 0.2,
                        color = clus_col, red = "umap_harmony", lable_size = 4) + theme_void() +
  #xlab("UMAP 1") + ylab("UMAP 2") +
  #coord_equal() +
  theme(legend.position = "none",
        plot.title = element_blank(),
        text = element_text(size = 7),
        plot.margin = unit(c(5,2,1,7), "pt"), #t,r,b,l c(.1,.1,-.4,-.1)
        #axis.title.x.bottom = element_blank(),
        #axis.title.y.left = element_text(margin = margin(r = 3), angle = 90)
      )

# dev.new(height=2.8, width=3, noRStudioGD = TRUE)
# ggsave("./Figures/01/st_UMAP_clusters.png", UMAP_ST, width = 3, height = 2.8, bg = "white", dpi = 1000)
UMAP_ST
```

<img src="../Figures/01/01c_st_UMAP_clusters.png"
data-fig-align="center" />

``` r
############################
# MARKER GENES FOR DOTPLOT #
############################
marker_genes_morf <- DEGs_table_morf %>%
  mutate(Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
  group_by(Direction) %>%
  top_n(-70, p_val_adj) %>%
  top_n(40, abs(log.pct.diff)) %>%
  {.->> top_40} %>%
  top_n(10, abs(avg_log2FC)) %>%
  arrange(pct.diff, Direction)

################################
# DOTPLOT MARKER GENE FUNCTION #
################################
get_order <- function(df, markers, markers_id, value){
  m <- set_names(markers, markers_id)
  value <- enquo(value) 
  
  df <- df %>%
    mutate(marker_id =  map_chr(.$marker, ~names(m)[which(m == as.character(.x))][1])) %>%
    nest(., .by=c(groups)) %>%
    mutate(., data = pmap(., ~filter(..2, marker_id == .[[1]]))) %>%
    mutate(., data = pmap(., ~arrange(..2, desc(!!(value))))) %>%
    unnest(cols = c("data")) %>%
    #mutate(marker = fct_reorder2(marker, Avg, Pct))
    mutate(row_id = cur_group_rows(), .by=c(groups))
  
  return(df[["marker"]]) # use rev() to reorder the dots depending on horizontal/vertical 
}
get_df_for_dotplot.fun <- function(obj, markers, groups="Clusters", gr_col, rect=TRUE, alpha=.2, gr_lvl=NULL){
  col <- c("#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D") # 
  markers_lvl <- set_names(seq_along(markers), markers)
  gr <- unique(obj[[groups]])[[1]]
  gr_lvl <- set_names(seq_along(gr), if(is.null(gr_lvl)){gr}else{gr_lvl})
  
  df <- obj %>%
    mutate(., FetchData(., vars = c(markers)) ) %>%
    as_tibble() %>%
    select(., .cell, any_of(c(groups, markers))) %>%
    pivot_longer(., cols = -c(".cell", groups), 
                 names_to = "marker", values_to = "values") %>%
    #mutate(Clusters = factor(.$Clusters, levels=names(clus_lvl))) %>%
    group_by(!!sym(groups), marker) %>%
    summarise(Avg = mean(values),
              Avg_ = mean(values)/max(values), # gives the median as a percentage of max value among spots if .5 you have a normal dist.
              Pct = sum(values > 0) / length(values) * 100, .groups="drop") %>%
    #group_by(marker) %>%
    #nest()
    mutate(Avg_p = Avg/max(Avg), .by = "marker") %>%
    mutate(marker = factor(.$marker, levels=get_order(., markers, markers_id, Avg))) %>% 
    mutate(., ymin = gr_lvl[as.character(.[[groups]])]-.5,
              ymax = gr_lvl[as.character(.[[groups]])]+.5) %>%
    mutate(., xmin = markers_lvl[as.character(.$marker)]-0.5,
              xmax = markers_lvl[as.character(.$marker)]+0.5) %>%
    mutate(!!sym(groups) := factor(unlist(.[,1]), levels=names(gr_lvl)))
  
  df_rect <- df %>% select(1, ymin,  ymax) %>% slice_max(., n=1, order_by=ymin, by=groups, with_ties=F)
  
  p <- ggplot(df, aes(x=marker, y=.data[[groups]])) +
    {if(rect == TRUE)
      list(new_scale_fill(),
      
      annotate(geom = "rect", ymin=df_rect$ymin, ymax=df_rect$ymax,
                  xmin=rep(0, length(gr)), xmax=rep(length(markers_lvl)+1, length(gr)),
                  fill = gr_col[1:length(gr)], colour = NA, alpha = alpha),
      coord_cartesian(expand = F) )} +
    
    new_scale_fill() +
    geom_point(aes(size = Pct, fill = Avg_), col = "gray", shape = 21) + #"#6A2595"
    scale_fill_gradientn(colours = col,
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression") +
    
    
    theme_bw() + labs(size="% expressed") +
    guides(size = guide_legend(override.aes = list(color = "black") )) +
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
          axis.text.y = element_text(size=12, color="black",hjust=.9 ),
          legend.text = element_text(size=10),
          legend.margin=margin(0,0,0,0), #t,r,b,l
          legend.box.margin=margin(0,-6,-10,-8),
          axis.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(.2,.2,.2,.2),units = "cm") #trbl
          ) 
  return(p)
}

###################################
# DOTPLOT TOP 10 DOWN AND UP DEGs #
###################################
# dev.new(width=6.5, height=1.9, noRStudioGD = TRUE)
gr_col <- c( "#B3B6FF","#F0B6AD") #"#BACCE5"
markers <-  marker_genes_morf$gene
markers_id <- rep("epi", length(markers))
groups="sp_annot"
(plot <- get_df_for_dotplot.fun(DATA, markers, groups="sp_annot", gr_col, rect=T, alpha=.27 ) + 
    guides(fill = guide_colourbar(barwidth = 5, barheight = .5),
           size = guide_legend(keyheight = .5, keywidth = .1 )) +
    theme(legend.position="top",
          
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.margin=margin(0,0,0,0), # moves the legend box
          legend.box.margin=margin(-3,1,-11,-30) )) # moves the legend 
```

<img src="../Figures/01/01d_top5_marker_genes_morf.png"
data-fig-align="center" />

``` r
# ggsave(paste0("./Figures/01/","top5_marker_genes_morf", ".pdf"), plot,  width = 6.5, height = 1.9)
# ggsave(paste0("./Figures/01/","top5_marker_genes_morf", ".png"), plot,  width = 6.5, height = 1.9, dpi = 1000)
```

``` r
# DATA@misc$cell_annot <- DATA@assays[["misc"]][["cell_annot"]]
DATA@misc$cell_annot %>% 
  mutate(cell_annot_1 = ifelse(.$cell_annot_1 == "Epithelial", "Keratinocyte", .$cell_annot_1)) %>%
  {. ->> DATA@misc$cell_annot }
```

``` r
####################################
# COLOUR PALLET SUBGROUPS FUNCTION #
####################################
# maximum levels in the subgroup is 16
s_col <- c("#FFD92F","#8DA0CB","#FC8D62","#E78AC3","#66C2A5","#E6F5C9","#eb6062","#377EB8","#984EA3","#4DAF4A","#E5C494","#B3B3B3","#FF7F00","#FFFF33","#A65628","#F781BF")
e_col <- c("#FFF2AE","#CBD5E8","#FDCDAC","#F4CAE4","#B3E2CD","#A6D854","#FBB4AE","#B3CDE3","#DECBE4","#CCEBC5","#F1E2CC","#CCCCCC","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC")

#### function ####
ColourPalleteMulti <- function(df, group, subgroup){

  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  #categories <- arrange(categories, group)
  s <- s_col # c(RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(8,"Set1"))
  e <- e_col # c(RColorBrewer::brewer.pal(8,"Pastel2"), RColorBrewer::brewer.pal(8,"Pastel1"))
  category.start <- (s[1:nrow(categories)]) # Set the top of the colour pallete
  category.end  <- (e[1:nrow(categories)]) # set the bottom

  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                          function(i){
                            colorRampPalette(colors = c(category.start[i],
                                                        category.end[i]))(categories[i,2])}))
  return(colours)
  #return(categories)
}

# choose higer and lower level:
# colours <- ColourPalleteMulti(candida_df, "rank", "taxName") %>% set_names(levels(candida_df$taxName))
#scales::show_col(scales::hue_pal()(8))
```

### Barplot of cell composition for st dataset

``` r
##################
# PLOT FUNCTIONS #
##################
get_order <- function(df, id, value){
  id <- enquo(id)
  id_name <- as_label(id)
  value <- enquo(value) 
  
  df %>%
    group_by(!!({id})) %>%
    add_tally(wt = !!(value), sort = F) %>%
    ungroup() %>%
    mutate(!!(id_name) := fct_reorder(!!(id), n))}

get_legend.fun <- function(df, title, id, colours, txt_size = 10){
  # https://stackoverflow.com/questions/73711918/how-to-remove-majority-of-blank-space-from-get-legend
  p <- ggplot() + geom_col(data = df, aes(fill=annot, y=values, x=.data[[id]])) +
    scale_fill_manual(title, values = colours[unique(pull(df, "annot"))]) +
    guides(fill = guide_legend(override.aes = list(size=2), keyheight = .7, keywidth = .7)) +
    theme(text = element_text(size = txt_size),)
  l <- cowplot::get_legend(p)
  return(l)
}

###################
# GET DECON ASSAY #
###################
# dev.new(height=6, width=6.6, noRStudioGD = TRUE)
decon_column_st <- DATA@misc$cell_annot 

####################
# ARRANGE ID ORDER #
####################
meta_group <- "cell_annot_1"
subgroup <- "cell_annot_2"
meta_group_lvl <- c("B cell", "T cell", "Keratinocyte", "fibroblast", "endothelial") #, "granulocytes",
meta_group_lvl <- c("B cell", "T cell", "Keratinocyte", "Fibroblast", "Endothelial","Myeloid",  "NK cells",  "Granulocyte", "ILC", "Lymphatics") # ,"NA" ,"Smooth muscle"

df <- decon_column_st %>%
  filter(grepl(paste0(meta_group_lvl,collapse = "|"),.$cell_annot_1)) %>%
  mutate(ID = str_extract(.cell, "^P\\d+")) %>%
  mutate(Percent_morf = 100 *values/sum(values), .by = c("ID", "sp_annot")) %>%
  mutate(cell_annot_1 = replace_na(as.character(.$cell_annot_1), "NA")) %>%

  # arrange the leves and order of colours:
  mutate(meta_group = factor(.[[meta_group]], levels = meta_group_lvl)) %>%
  mutate(annot = .[[subgroup]]) %>%
  mutate("fibroblast" = ifelse(grepl("14|1|22|0|9", .$clus), .$values, NA)) %>%
  get_order(., id = `ID`,`fibroblast`) %>%
  arrange(meta_group) %>%
  mutate(annot = factor(.[[subgroup]], levels = unique(.[[subgroup]])))

  # check that the percent adds upt to 100:
  df %>%
  nest(data = -c(sp_annot, ID)) %>%
  mutate(p = map_dbl(data, ~sum(.x[["Percent_morf"]])) )
```

    # A tibble: 16 × 4
       ID    sp_annot data                      p
       <fct> <chr>    <list>                <dbl>
     1 P031  SubMuc   <tibble [1,547 × 15]>   100
     2 P080  SubMuc   <tibble [1,849 × 15]>   100
     3 P097  SubMuc   <tibble [4,211 × 15]>   100
     4 P105  SubMuc   <tibble [2,168 × 15]>   100
     5 P105  epi      <tibble [371 × 15]>     100
     6 P107  epi      <tibble [419 × 15]>     100
     7 P108  SubMuc   <tibble [714 × 15]>     100
     8 P114  SubMuc   <tibble [1,855 × 15]>   100
     9 P114  epi      <tibble [490 × 15]>     100
    10 P118  SubMuc   <tibble [3,016 × 15]>   100
    11 P031  epi      <tibble [521 × 15]>     100
    12 P097  epi      <tibble [428 × 15]>     100
    13 P107  SubMuc   <tibble [421 × 15]>     100
    14 P108  epi      <tibble [359 × 15]>     100
    15 P118  epi      <tibble [710 × 15]>     100
    16 P080  epi      <tibble [443 × 15]>     100

``` r
colours <- ColourPalleteMulti(df, "meta_group", "annot") %>% 
  set_names(unique(df$annot))
# scales::show_col(colours)
##########
# LEGEND #
##########
c <- df %>%
  nest(data = -meta_group) %>%
  mutate(legend = pmap(., ~get_legend.fun(..2, ..1, ".cell", colours, txt_size = 15)))

# plot_grid(plotlist = c$legend)
##########
# TABLE #
##########
round.fun <- function(x){
    df <- x %>%
    mutate_if(
    is.numeric,
    ~ round(.x, digits = 1))
  }

decon_percent <- df %>%
  #nest(data = -meta_group) %>%
  summarise(sum = sum(Percent_morf), .by = c("sp_annot","meta_group")) %>%
  arrange(sp_annot) %>%
  pivot_wider(names_from = c("sp_annot"), values_from = "sum", values_fill= 0 ) %>%
  mutate("All spots" = SubMuc+epi)
  
t <- list()
decon_percent <- df %>%
  mutate( ID = str_extract(.cell, "^P\\d\\d\\d")) %>%
  summarise(sum = sum(Percent_morf), .by = c("meta_group", "ID")) %>%
  arrange(ID) %>%
  pivot_wider(names_from = c("ID"), values_from = "sum", values_fill= 0) %>% 
  group_by(meta_group) %>%
  group_walk(., ~ {t <<- append(t, wilcox.test(data=.x, 
                                               x=as.double(.x[1:4]), 
                                               y=as.double(.x[5:8]))$p.value)}) %>%
  ungroup() %>%
  mutate(p.value = unlist(t)) %>%
  left_join(., decon_percent, by="meta_group") %>%
  mutate(across(-any_of(c("meta_group", "p.value")), ~./sum(.)*100)) %>%
  mutate(DMPA = rowMeans(.[,2:5]), .before = "p.value" ) %>%
  mutate(Control = rowMeans(.[,6:9]), .before = "p.value" ) %>%
  mutate(Average = rowMeans(.[,2:9]), .before = "DMPA" ) %>%
  round.fun(.) %>%
  add_column(" " = "", .before = "DMPA") %>%
  add_column("  " = "", .before = "SubMuc") %>%
  dplyr::rename("   "="meta_group", "Submucosal spots"="SubMuc", "Epithelial spots"="epi")

write.xlsx(decon_percent, paste0(result_dir, "deconvolution_scdc_percent_.xlsx"))
#########################################
# ST CELL TYPE DISTRIBUTION PER SUBJECT #
#########################################
bar <- ggplot(df, aes(fill=annot, y=Percent_morf, x=ID)) + 
    geom_bar(position="stack", stat="identity") +
    theme_minimal() + 
    scale_fill_manual(values = colours) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.4, hjust=1),
          axis.title = element_blank(),
          text = element_text(size = 15),
          legend.position = "none") +
    #facet_wrap(~group ) +
    facet_grid(cols = vars(groups), rows = vars(sp_annot), 
               labeller = as_labeller(c(epi='Epithelial spots', SubMuc='Submucosal spots', ctrl="Controls", DMPA="DMPA")),
               scales = "free_x", space = "free_x", switch = "x")

l <- plot_grid(plotlist = c$legend, nrow =4 , byrow = T , 
               rel_heights = c(.3,.3,.3,.3), rel_widths = c(.5,0.2), align = "hv", axis = "t")
(p <- plot_grid(bar, NULL,rel_widths = c(.7,1)))
```

<img src="../Figures/01/01e_decon_bar_plot.png"
data-fig-align="center" />

``` r
# ggsave(paste0("./Figures/01/", "decon_bar_plot.png"), p, height=6, width=6.6, bg = "white")
```

``` r
col <- c("#FFD92F","#8DA0CB","#FC8D62","#66C2A5","#E78AC3","#A6D854","#377EB8","#eb6062","#4DAF4A","#984EA3","#B3B3B3","#E5C494","#FF7F00","#FFFF33","#A65628","#F781BF") 
cell_type <- c("B cell", "T cell","Epithelial","Endothelial", "Fibroblast", "Myeloid", "Granulocyte", "NK cells","Lymphatics", "ILC", "Unknown")
colours <- set_names(col[1:length(cell_type)],cell_type)

# dev.new(height=5, width=5, noRStudioGD = TRUE)
UMAP <- DATA_r %>%
  filter(!(cell_annot_1 == "Smooth muscle")) %>%
  plot_clusters.fun(., seed = 127, cluster="cell_annot_2", title = "",
                  txt_size = 15, lable = "cell_annot_1",
                  red = "UMAP", color = colours, lable_size = 5) + 
  xlab("") + ylab("UMAP 2") + theme_nothing() + theme(axis.title.x = element_text(vjust=1))

ggsave("./Figures/01/sc_reference_UMAP.png", UMAP, width = 4.4, height = 4.4, bg = "white", dpi = 400)
UMAP
```

<img src="../Figures/01/01e_sc_reference_UMAP.png"
data-fig-align="center" />

### Plot deconvolution cell proportion on tissue

``` r
#################################
# CELL TYPE PROPORTION PER SPOT #
#################################
# dev.new(width=4, height=6, noRStudioGD = TRUE)
# dev.new(width=7, height=6, noRStudioGD = TRUE)
col <- c("#FFD92F","#8DA0CB","#FC8D62","#66C2A5","#E78AC3","#A6D854","#377EB8","#eb6062","#4DAF4A","#984EA3","#B3B3B3","#E5C494","#FF7F00","#FFFF33","#A65628","#F781BF") # "#E5C494"
cell_type <- c("B cell", "T cell", "Keratinocyte","Endothelial","Fibroblast", "Myeloid", "Granulocyte", "NK cells","Lymphatics", "ILC", "Unknown") # "Smooth muscle"
cell_col <- set_names(col[1:length(cell_type)], cell_type)

DefaultAssay(DATA) <- "RNA"
# DATA@misc$cell_annot <- DATA@assays[["misc"]][["cell_annot"]]
DAT <- DATA %>%  filter(orig.ident == "P118")

(plot_ <- DAT %>%
  plot_cell_pie.fun(.,
                   assay = "RNA",
                   ct.res = "cell_annot_1",
                   ct.select = cell_type, 
                   radius_adj = -2,
                   zoom = "zoom",
                   col = cell_col,
                   alpha = 1,
                   ncol = 1,
                   #annot_col = "#dbd9d9",
                   annot_line = .5,
                   img_alpha = 0) + 
    guides(fill = guide_legend(keyheight = .7, keywidth = .7)) +
    theme(plot.margin = unit(c(-3, -1, -3, -3), "lines"),
          legend.title = element_blank(),
          legend.box.margin=margin(0,10,-0,-30), # moves the legend
          legend.text = element_text(size = 8)))
```

<img src="../Figures/01/01f_cell_prop_tissue.png"
data-fig-align="center" />

``` r
# dev.new(height=2.3, width=3.3, noRStudioGD = TRUE)
# ggsave("./Figures/04/cell_prop_tissue.png", plot_, width = 3.6, height = 2.6, dpi = 500, bg = "white")
# ggsave("./Figures/04/cell_prop_tissue.png", plot_, width = 9, height = 7, dpi = 500, bg = "white")
```

``` r
sampleid <- c("P107", "P108", "P114", "P097","P118", "P105", "P080", "P031")
sample_id <- c("P118", "P105", "P080", "P031", "P097", "P108", "P114", "P107")
# dev.new(height=7, width=6, noRStudioGD = TRUE)
# dev.new(height=3.5, width=8, noRStudioGD = TRUE)
plot <- plot_spatial.fun(DATA, 
      sampleid= c("P118","P097"), #sample_id,
      txt = F,
      geneid = "nFeature_RNA",
      zoom = "zoom",
      ncol = 1,
      img_alpha = 0,
      point_size = .9) + 
  #ylim(600,l$min_row) +
  theme(panel.spacing.y = unit(-9, "lines"), plot.margin = unit(c(-9, -1, -5, -1), "lines"))


plots <- map(c("P118","P097"), 
    ~plot_spatial.fun(DATA, 
      sampleid=.x,
      txt = F,
      geneid = "sp_annot",
      zoom = "zoom",
      ncol = 1,
      img_alpha = 1,
      point_size = 0) + theme_nothing() + coord_equal())

plots_ <- list(plots[[2]] + theme(plot.margin = unit(c(5,-1,-1,-1.5), "cm")), #t,r,b,l
               plots[[1]] + theme(plot.margin = unit(c(-1,-1,-1,-1.5), "cm")) ) #l,b,r,t

(p <- plots_[[1]] + inset_element(plots_[[2]], 0, .4, .7, 1, align_to = 'full', on_top=F) )

p+theme(plot.margin = unit(c(-3,-1,-1,-1.5), "cm"))


legend_1 <- get_legend(plots$res_1[[2]] + theme(legend.position="right"))
legend_2 <- get_legend(plots$res_2[[1]] + theme(legend.position="right"))
legend <- plot_grid( legend_1, legend_2, ncol = 1)
combined <- wrap_plots(plotlist=c(plots$res_1, plots$res_2), nrow = 8, byrow = F) & theme(legend.position="none")
combined <- plot_grid( combined, legend, ncol = 2, rel_widths = c(1, .3)) 
combined
```
