Figure 3
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
library(ggrepel)
library(png)
library(grid)
library(patchwork)
library(openxlsx)

source("../../bin/spatial_visualization.R")
source("../../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../../results/08_spatial_dist/"
result_dir <- "../../results/09_figures/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
epi_clus <- "^5$|^6$|^7|^9" # non-filt

#############
# LOAD DATA #
#############
# DEGs_table <- read_csv(paste0(input_dir,"subset_100/DGEs_condition_wilcox_epi_100.csv")) %>% filter(groups == "DMPA")
DEGs_table <- read_csv(paste0("../../results/05_DGE_clusters_st_data/","DGEs_clusters_wilcox.csv"))
DATA <- readRDS(paste0(input_dir,"seuratObj_spatial_dist.RDS"))
```

### UMAP of cluster resolutions

``` r
# fig.asp=5/10,
###############
# FILTER DATA #
###############
keep <- "Basal|Upper IM|Lower IM|Superficial"
names <-  c("Superficial","Upper IM","Lower IM","Basal")
DAT <- DATA %>%
  arrange(spatial_dist) %>%
  filter( sp_annot == "epi" ) %>%
  filter(!(grepl("\\d", .$layers))) %>%
  mutate(groups = if_else(groups =="ctrl", "Control", .$groups))

################
# COLOR PALLET #
################
# centroids2d <- as.matrix(t(t(DAT@reductions$umap_harmony@cell.embeddings) %*% mm)/Matrix::colSums(mm))
lable_df <- DAT %>%
    group_by(Clusters) %>%
    select(Clusters, contains("umapharmony")) %>% 
    summarize_all(mean)

col1 <- set_names(c("#E41A1C","#FF7F00","#C77CFF","#984EA3"), c("6","9","7","5"))
col2 <- c("grey90","grey70", "orange3", "firebrick", "purple4")

############
# PLOTTING #
############
p1 <- plot_clusters.fun(DAT, cluster="Clusters", txt_size = 10, dot_size = 0.05,
                        color = col1, red = "umapharmony", lable_size = 3) + theme_void() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  #coord_equal() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA),
        plot.title = element_blank(),
        text = element_text(size = 7),
        plot.margin = unit(c(5,2,1,7), "pt"), #t,r,b,l c(.1,.1,-.4,-.1)
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(margin = margin(r = 3), angle = 90)
      )
  
p2 <- plot_genes.fun(DAT, gene="spatial_dist", point_size = .05,
                     col = col2, red = "umapharmony", lable_size = 3) + 
  xlab("UMAP 1") + ylab("UMAP 2") +
  #coord_equal() +
  theme(panel.border = element_rect(colour = "black", fill=NA),
        plot.title = element_blank(),
        text = element_text(size = 7),
        plot.margin = unit(c(1,2,1,7), "pt"), #t,r,b,l c(.1,.1,-.4,-.1)
        axis.title.x.bottom = element_text(margin = margin(t = 3)),
        axis.title.y.left = element_text(margin = margin(r = 3), angle = 90)
      )

# dev.new(width=1.3, height=1.9, noRStudioGD = TRUE)
(A <- plot_grid(ncol = 1, #labels = c('A'),
          plotlist = list(p1, p2)) )
```

<img src="../Figures/03/03a_plot_resolution.png"
data-fig-align="center" />

### Plot spatial distance on tissue

``` r
####################
# SPATIAL DISTANCE #
####################
# dev.new(width=2.3, height=1.9, noRStudioGD = TRUE)
col1 <- set_names(c("#E41A1C","#FF7F00","#C77CFF","#984EA3"), c("6","9","7","5"))
col2 <- c("grey70", "orange3", "firebrick", "purple4")
scale_color_pseudotime <- function(palette = col2, direction = 1, ...) {
  # get normalized values of spatial dist
  pseudotime <- as.data.frame(list("sp"=DAT$spatial_dist))
  pseudotime <-  arrange(pseudotime, sp)
  x <- rowMeans(pseudotime)
  x <- x/max(na.omit(x))
  
  # generate colors
  pal <- palette
  my_pal <- colorRampPalette(pal)(99)[x * 98 + 1]
  # map them as ggplot pallet
  scale_color_gradientn(colors = my_pal, na.value = "grey90", ...)
}

feat <- list("Clusters", "spatial_dist")
p <- map(feat,
          ~plot_spatial.fun(DATA,
                   geneid = .x,
                   sampleid = c("P118", "P097"),
                   zoom = "zoom",
                   lab = F,
                   save_space = F,
                   col = col1,
                   alpha = 1,
                   ncol = 1, 
                   #annot_col = "#dbd9d9",
                   annot_line = .1,
                   img_alpha = 0,
                   point_size = .35) + #.3
    #annotate("text", x = 500, y = 450, size = 4, label = "Control") + # col = "#CB251C",
    theme(legend.position = "none",
          panel.spacing.y = unit(-3, "lines"),
          plot.margin = unit(c(-.3, -.4, -.1, -.9), "lines")) )

p[[2]] <- p[[2]] + scale_color_pseudotime(palette = col2) 

# dev.new(width=2.6, height=1.9, noRStudioGD = TRUE)
B <- plot_grid(ncol = 2, #labels = c('B'),  hjust = -.2,
          plotlist = p) 

(B <- B + draw_label("DMPA", x = c(.05), y = c(.48), hjust = 0, size = 9) +
    draw_label("Control", x = c(.1), y = c( .94), hjust = 0, size = 9) +
    draw_label("Gene exp. similarity", x = c(.53), y = c( .94), hjust = 0, size = 8) )
```

<img src="../Figures/03/03b_spatial_trajectory_tissue.png"
data-fig-align="center" />

``` r
# dev.new(width=3.35, height=1.9, noRStudioGD = TRUE)

(A_B <- plot_grid( A, B, ncol=2, rel_widths = c(.55,1), labels = c('A','B'),  hjust = -.2))
```

<img src="../Figures/03/03ab_UMAP_tissue.png" data-fig-align="center" />

``` r
# ggsave("./Figures/03/UMAP_tissue.png", A_B, width = 3.35, height = 1.9, bg = "white", dpi = 300)
```

### Plot dotplot and tissue

``` r
####################
# DOTPLOT FUNCTION #
####################
ord <-  c("Superficial","Upper IM","Lower IM","Basal", "8","3","4","0","2","1")
clus_1 <- c("#F40014","#FF8A1E","#C77CFF","#A750B2") %>% set_names(., ord[1:4])
ID <- c("P118", "P105", "P080", "P031", "P097", "P108", "P114", "P107")
col_gr <- set_names(c("#813340", "#AC5A64", "#DA84A4", "#FFB1BA","#125E6B","#3A7D8B", "#80C2D0", "#C0DFEE"),ID)

layer_dotplot.fun <- function(DATA, feat, 
                              facet = TRUE, 
                              txt_size = 5,
                              density = FALSE,  
                              x_max=NULL, 
                              col){
  
  rects <- DATA %>%
    group_by(layers) %>%
    summarise(., ystart=min(spatial_dist, na.rm=T), yend=max(spatial_dist, na.rm=T),
              Q1=quantile(spatial_dist, probs = 0.009, na.rm=T),
              Q3=quantile(spatial_dist, probs = 0.90, na.rm=T)) %>%
    filter(!(is.infinite(.$ystart))) %>%
    mutate(Q1 = ifelse(.$Q1 == min(.$Q1), 0,.$Q1)) %>%
    mutate(Q3 = ifelse(.$Q3 == max(.$Q3), max(.$yend),.$Q3)) %>%
    mutate(Q3 = ifelse(.$layers == "Lower IM", .$Q3+.3,.$Q3)) %>%
    #mutate(Q1 = ifelse(.$layers == "Basal", .$Q1+.3,.$Q1)) %>%
    arrange(layers) %>% ungroup()
  
  # rects <- DATA %>%
  #   group_by(layers) %>%
  #   summarise(., ystart=min(spatial_dist, na.rm=T), yend=max(spatial_dist, na.rm=T),
  #             Q1=quantile(spatial_dist, probs = 0.1, na.rm=T),
  #             Q3=quantile(spatial_dist, probs = .9, na.rm=T)) %>%
  #   filter(!(is.infinite(.$ystart))) %>%
  #   mutate(Q1 = ifelse(.$Q1 == min(.$Q1), 0,.$Q1)) %>%
  #   mutate(Q3 = ifelse(.$Q3 == max(.$Q3), max(.$yend),.$Q3)) %>%
  #   #mutate(Q3 = ifelse(.$layers == "Basal_2", .$Q3+.3,.$Q3)) %>%
  #   arrange(layers) %>% ungroup() %>%
  #   filter(., grepl(keep, .$layers))
  
  DAT <- DATA %>%
    #filter(., grepl("SubMuc", DATA$sp_annot) & !(grepl(epi_clus, DATA$Clusters))) %>%
    #filter(., (grepl(keep, DATA$layers))) %>%
    mutate(., FetchData(., vars = c(feat), slot = "data") ) %>%
    select(groups, ID="orig.ident", layers, all_of(c(feat)), spatial_dist) %>%
    mutate(ID = factor(.$ID, levels = names(col_gr)))
   
  dens <- ggplot() + 
    geom_density(data = DAT %>% filter(., .[[feat]] != 0 ), aes(y=spatial_dist, color=ID)) + 
    scale_y_continuous(expand = c(0, 0)) + ggpubr::clean_theme()+
    scale_color_manual(values = col_gr) +
    theme(panel.border = element_blank(), 
          plot.margin = unit(c(0,0,0,0), "inches"))

  if(facet == TRUE){facets <- facet_wrap(~groups, ncol = 2)
                    w <- 3}
  else{facets <- NULL
       w <- 1}
  
  dot <- ggplot() +
    #ggtitle(feature) +
    geom_rect(data = rects, alpha = 0.1, show.legend=FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3, fill = layers)) +
    geom_jitter(data = DAT, aes(x=.data[[feat]], y=spatial_dist, col=layers), 
                width = 0.1, alpha = 0.7, size=.3) + #
    scale_fill_manual(values = col) + 
    scale_colour_manual(values = col) +
    coord_cartesian(clip = "off") +
    guides(col = guide_legend(override.aes = list(size=2), keyheight = .7, keywidth = .7, title = "Cluster")) +
    scale_y_reverse(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0)) +
    #scale_x_continuous(expand = c(0, 0)) +
    {if(!(is.null(x_max))){xlim(-.5, x_max)}} +
    facets +
    my_theme + ylab("Similarity in gene expression") +
    theme(legend.position="top", legend.justification = "left",
          #legend.position = c(.1, 1.15), legend.direction = "horizontal",
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.margin = unit(c(2,.0,-2,.5), "lines"),
          panel.border = element_blank(),
          text = element_text(size = txt_size),
          axis.line = element_line(linewidth = .3),
          #axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x.bottom = element_text(margin = margin(t = .5)),
          axis.title.y.left = element_text(margin = margin(r = 5)),
          panel.grid.major = element_line(linewidth = 0.2),
          panel.grid.minor = element_line(linewidth = 0.1))
  
  if(density == TRUE){dot <- dot + dens + plot_layout(ncol=2, nrow=1, widths=c(w, 1), guides = "collect")}
  
  return(dot)
}
```

``` r
###############
# FILTER DATA #
###############
keep <- "Basal|Upper IM|Lower IM|Superficial"
DAT <- DATA %>%
  filter(!(grepl("\\d", .$layers))) %>%
  mutate(groups = if_else(groups =="ctrl", "Control", .$groups))

##############################
# CLUS DOTPLOT PER EPI LAYER #
##############################
# dev.new(height=1.9, width=6.7, noRStudioGD = TRUE)
# dev.new(height=3.2, width=3.54, noRStudioGD = TRUE) #with legend
feat <- c("KRT78","TGM1","KRT6C","MIR205HG" )
dot_fig <- map(feat, 
               ~layer_dotplot.fun(DAT, .x, 
                                  facet = F, 
                                  txt_size = 9,
                                  col = clus_1)) # , x_max = 6
# dot_fig[[1]]
C_1 <- wrap_plots(dot_fig, ncol = 4, guides = "collect" ) + 
  plot_layout(axis_titles = "collect") & 
  theme(#legend.position = c(.1, 1.15), legend.direction = "horizontal",
        plot.margin = unit(c(.2,0,0,.2), "lines"),
        legend.position = 'top', 
        legend.justification = "right", 
        legend.background = element_rect(colour ="gray", linewidth=.2),
        legend.margin=margin(1,2,1,1), # moves the legend box
        legend.box.margin=margin(1,1,-4,0), # moves the legend
        #legend.box.margin=margin(-30,0,-15,0), 
        # axis.title.y.left = element_text(margin = margin(r = 5))
        ) 

C_1 <- plot_grid(C_1, NULL, rel_widths = c(1,.01))
ggsave("./Figures/03/Fig_03C1.png", C_1, width = 6.7, height = 2.4) #, dpi = 300

############################
# CLUS EXPRESION ON TISSUE #
############################
DAT <- DATA %>%
    filter(grepl("P118", .$orig.ident))

col_feat <- c("#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D") # Purples
p <- map(feat, 
        ~plot_st_feat.fun( DAT,
                           geneid = .x,
                           txt = F,
                           save_space = F,
                           zoom = "zoom",
                           col = col_feat,
                           alpha = .9,
                           ncol = 2, 
                           #annot_col = "#dbd9d9",
                           annot_line = .1,
                           img_alpha = 0,
                           point_size = .6)) 

# dev.new(height=1, width=6.7/4, noRStudioGD = TRUE)
# p_list[[1]]
# dev.new(width=6.7, height=1, noRStudioGD = TRUE) 

p_list <- map(p, ~ .x + theme(plot.margin = unit(c(-35, -2, -30, -15), "pt"), #t,r,b,l
                #panel.background = element_rect(fill = "white"),
                panel.spacing.y = unit(-1.5, "lines"),
                panel.spacing.x = unit(-2, "lines"),
                legend.key.height = unit(0.2, "cm"),
                legend.box.margin=margin(-20,0,0,0),
                legend.title = element_blank(),
                legend.position = "none"
                ) )

#########################
# COMBINE PANEL A AND C #
#########################
# patchwork does not respect margins when nesting plots, use cowplot instead!
# dev.new(width=3.54, height=3.2, noRStudioGD = TRUE) 
c <- plot_grid( plotlist = p_list, ncol = 4, rel_heights = c(1), rel_widths = c(1))
C_2 <- plot_grid(NULL, c, rel_widths = c(.05,1))
(C <- plot_grid( C_1, C_2, ncol=1, rel_heights = c(3,.9), labels = c('C'),  hjust = -.2))
```

<img src="../Figures/03/03c_dot_and_tissue.png"
data-fig-align="center" />

``` r
# ggsave("./Figures/03/dot_and_tissue.png", C, width = 3.4, height = 3.2, bg = "white", dpi = 300)
```

``` r
# dev.new(width=3.4, height=1.9+3, noRStudioGD = TRUE)
A_B_C <- plot_grid( A_B, C, ncol=1, rel_heights = c(.6,1))

 plot_grid( A_B_C, ncol=2, rel_widths = c(1,.03))
```

<img src="../Figures/03/03abc_combined.png" data-fig-align="center" />

``` r
################################
# DOTPLOT MARKER GENE FUNCTION #
################################
library(ggnewscale)
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
    #mutate(Pct = .$Pct*) %>% # scaling
    #group_by(marker) %>%
    #nest()
    mutate(Avg_p = Avg/max(Avg), .by = "marker") %>%
    mutate(marker = factor(.$marker, levels=get_order(., markers, markers_id, Avg))) %>% 
    mutate(., ymin = gr_lvl[as.character(.[[groups]])]-.5,
              ymax = gr_lvl[as.character(.[[groups]])]+.5) %>%
    mutate(., xmin = markers_lvl[as.character(.$marker)]-0.5,
              xmax = markers_lvl[as.character(.$marker)]+0.5) %>%
    mutate(!!sym(groups) := factor(unlist(.[,1]), levels=names(gr_lvl)))
  
  # df <<- df
  df_rect <- df %>% select(1, ymin,  ymax) %>% slice_max(., n=1, order_by=ymin, by=groups, with_ties=F)
  
  p <- ggplot(df, aes(x=marker, y=.data[[groups]])) +
    {if(rect == TRUE)
      list(new_scale_fill(),
      
      annotate(geom = "rect", ymin=df_rect$ymin, ymax=df_rect$ymax,
                  xmin=rep(0, length(gr)), xmax=rep(length(markers_lvl)+1, length(gr)),
                  fill = gr_col[1:length(gr)], colour = NA, alpha = alpha),
      coord_cartesian(expand = F) )} +
    
    new_scale_fill() +
    geom_point(aes(size = Pct, fill = Avg_, col = Avg_), shape = 21) + #"#6A2595"
    scale_fill_gradientn(colours = col, aesthetics = c("fill","colour"),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression") +
    
    #scale_size(max_size = 3) +
    scale_size_area(max_size = 3) +
    theme_bw() + labs(size="% expressed") +
    guides(size = guide_legend(override.aes = list(color = "black") )) +
    theme(text =  element_text(size=10),
          axis.text.x = element_text( angle=45, hjust=1, color="black"),
          axis.text.y = element_text( color="black",hjust=.9 ),
          legend.text = element_text(size=4),
          legend.margin=margin(0,0,0,0), #t,r,b,l
          legend.box.margin=margin(0,-6,-10,-8),
          axis.title = element_blank(),
          panel.grid.major = element_line(linewidth = .1), 
          #panel.grid.minor = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(.2,0,.2,0),units = "cm") #trbl
          ) 
  return(p)
}
```

### Dotplot of epithelial marker genes

``` r
###################################
# DOTPLOT TOP 10 MARKER GENES EPI #
###################################
library(ggnewscale)
# dev.new(width=6.5, height=3, noRStudioGD = TRUE)
gr_col <-c("#EB5253","#FA9F50","#D7A1FF","#984EA3")
markers <- c("RNASE7", "HMOX1", "DUOXA2", "DHRS9", "CHAC1", "PLA2G4D", "CEACAM7", "VSIG10L", "FAM3D", "LINC02487", "ENDOU", "LTB4R", "TMEM184A", "FAM83G", "ALOX15B", "COL7A1", "STMN1", "THSD4", "FGFR2", "LAMA5")
markers_id <- c(rep("Superficial",5), rep("Upper IM",5),rep("Lower IM",5),rep("Basal",5))
groups="layers"
gr_lvl <- rev(c("Superficial", "Upper IM", "Lower IM", "Basal"))
D <- DATA %>%
    filter(grepl(epi_clus, .$Clusters) & sp_annot == "epi") %>%
    get_df_for_dotplot.fun(., markers, groups="layers", gr_col, rect=T, gr_lvl=gr_lvl)  +
    theme(legend.position="bottom") + guides(fill = guide_colourbar(barwidth = 5, barheight = .5),
                                             size = guide_legend(keyheight = 1, keywidth = .1 )
                                             ) +# horizontal
  theme(#legend.box.margin=margin(-10,0,0,0),
        panel.border = element_rect(colour = "gray", fill=NA),
        legend.spacing.x = unit(1.5, 'pt'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.justification = "right",
        legend.margin=margin(0,0,0,0), # moves the legend box
        legend.box.margin=margin(-13,1,-5,-30), # moves the legend
        axis.text.y = element_text(margin = margin(r = -1)) )  
# dev.new(width=8, height=3, noRStudioGD = TRUE) # horizontal
# dev.new(width=4, height=5.5, noRStudioGD = TRUE) # vertical
# ggsave(paste0("./Figures/05/","Marker_genes_epi_h", ".pdf"), plot,  width = 6, height = 3)
# ggsave(paste0("./Figures/05/","Marker_genes_epi_v", ".pdf"), plot,  width = 4, height = 5.5)

# dev.new(width=3.4, height=2, noRStudioGD = TRUE) # horizontal
(D <- plot_grid(D, labels = c('D'), vjust = .9,  hjust = -.2))
```

<img src="../Figures/03/03d_marker_genes_epi.png"
data-fig-align="center" />

``` r
#############################
# COMBINE ALL FIGURE PANELS #
#############################
# dev.new(width=3.54, height=7, noRStudioGD = TRUE)
Figure3 <- plot_grid( A_B, C, D, ncol=1, rel_heights = c(.6,.85,.55)) 
ggsave("./Figures/Figure-3.pdf", Figure3, width = 3.54, height = 6.7, bg = "white", dpi = 1000)

Figure3
```

<img src="../Figures/03/Figure-3.png" data-fig-align="center" />
