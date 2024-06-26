Figure 6
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
input_dir <- "../../results/08_spatial_dist/SM/"
result_dir <- "../../results/09_figures/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
epi_clus <- "^5$|^6$|^7|^9" # non-filt

#############
# LOAD DATA #
#############
# DEGs_table <- read_csv(paste0(input_dir,"subset_100/DGEs_condition_wilcox_epi_100.csv")) %>% filter(groups == "DMPA")
DEGs_table <- read_csv(paste0("../../results/06_DGE_condition_st_data/","DGEs_condition_wilcox.csv"))
DATA <- readRDS(paste0(input_dir,"seuratObj_spatial_dist_SM.RDS"))
ihc <- readPNG("/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/resources/CD20_IHC.png")
```

### Volcanoplot and IHC image (Figure 6A-C)

``` r
# width 170 mm, 6.7 inches
###########################
# VOLCANO PLOT SUBMUCOSA #
###########################
ord <- c("Superficial","Upper IM","Lower IM","Basal","8","3","4","0","2","1")
DEGs_filt <- DEGs_table %>% 
  filter((grepl(paste0(ord, collapse = "|^"), .$subgroup))) %>%
  mutate(layers = factor(.$layers, levels = ord)) %>%
  filter(p_val < 0.099) #%>%
  #filter(p_val_adj < 0.05)

set.seed(1);A <- DEGs_filt %>%
  Volcano.fun_logFC(., "layers", 
        y.axis="p-value", 
        lab_size = 2.1, dot_size = .3,
        up=c(.2, 0.05), down = c(-.2, 0.05)) + # labeling: (logFC, p-value)
  ylab("avg. Log2 fold change")  +
  theme(legend.position = c(.92,.89), #y,x
        plot.margin = unit(c(1, .5, .2, 1), "lines")) #t,r,b,l

##################
# CD20 IHC IMAGE #
##################
# dev.new(width=1.48, height=3.1, noRStudioGD = TRUE) 
g <- grid::rasterGrob(ihc, interpolate=TRUE) #, height = 1, width = 1

C <- ggplot() +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_nothing() +
  theme(rect = element_blank(), # removes the box around the plot)
        plot.margin = unit(c(10,0,0,0), "pt")) #t,r,b,l  

#########################
# COMBINE PANEL A AND C #
#########################
# dev.new(width=6.7, height=3.2, noRStudioGD = TRUE) 
C <- plot_grid(C, rel_widths = c(1), labels = c('C'), hjust = .5)
(A_C <- plot_grid(A, C, rel_widths = c(3,.85), labels = c('A'), hjust = 0))
```

<img src="../Figures/06/05ac_volcano_plot_IHC.png"
data-fig-align="center" />

### Plot dotplot and tissue

``` r
####################
# DOTPLOT FUNCTION #
####################
names <- c("Basal", "0", "1", "2", "3", "4", "8", "10")
clus_1 <- c("#984EA3","#CD9600","#7CAE00","#e0e067","#00A9FF","#377EB8","#00BFC4","#FFFF33") %>% set_names(., names)
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
              Q1=quantile(spatial_dist, probs = 0.1, na.rm=T),
              Q3=quantile(spatial_dist, probs = .9, na.rm=T)) %>%
    filter(!(is.infinite(.$ystart))) %>%
    mutate(Q1 = ifelse(.$Q1 == min(.$Q1), 0,.$Q1)) %>%
    mutate(Q3 = ifelse(.$Q3 == max(.$Q3), max(.$yend),.$Q3)) %>%
    #mutate(Q3 = ifelse(.$layers == "Basal_2", .$Q3+.3,.$Q3)) %>%
    arrange(layers) %>% ungroup() #%>%
    #filter(., grepl(keep, .$layers))
  
  DAT <- DATA %>%
    #filter(., grepl("SubMuc", DATA$sp_annot) & !(grepl(epi_clus, DATA$Clusters))) %>%
    #filter(., (grepl(keep, DATA$layers))) %>%
    mutate(., FetchData(., vars = c(feat), slot = "data") ) %>%
    select(groups, ID="orig.ident", layers, all_of(c(feat)), spatial_dist) %>%
    mutate(ID = factor(.$ID, levels = names(col_gr)))
   
  dens <- ggplot() + 
    geom_density(data = DAT %>% filter(., .[[feat]] != 0 ), aes(y=spatial_dist, color=ID)) + 
    scale_y_reverse(expand = c(0, 0)) + ggpubr::clean_theme()+
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
    #scale_y_reverse(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    #scale_x_continuous(expand = c(0, 0)) +
    {if(!(is.null(x_max))){xlim(-.5, x_max)}} +
    facets +
    my_theme + ylab("Similarity in gene expression") +
    theme(legend.position="top", legend.justification = "left",
          #legend.position = c(.1, 1.15), legend.direction = "horizontal",
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.margin = unit(c(2,.2,-2,.3), "lines"),
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
keep <- "^5|8|3|^0"
DAT <- DATA %>%
  filter(grepl(keep, .$Clusters)) %>%
  mutate(groups = if_else(groups =="ctrl", "Control", .$groups))

##############################
# CLUS DOTPLOT PER EPI LAYER #
##############################
# dev.new(height=1.9, width=6.7, noRStudioGD = TRUE)
# dev.new(height=2.4, width=6.7, noRStudioGD = TRUE)  # with legend
feat <- c("IGHA2","IGLC1","JCHAIN","OLFM4" )
dot_fig <- map(feat, 
               ~layer_dotplot.fun(DAT, .x, 
                                  facet = T, 
                                  txt_size = 9,
                                  col = clus_1)) # , x_max = 6
# dot_fig[[1]]
B_1  <- wrap_plots(dot_fig, nrow = 1, guides = "collect" ) + 
  plot_layout(axis_titles = "collect") & 
  theme(#legend.position = c(.1, 1.15), legend.direction = "horizontal",
        plot.margin = unit(c(1,.2,0,.3), "lines"),
        legend.position = 'top', 
        legend.justification = "right", 
        legend.background = element_rect(colour ="gray", linewidth=.2),
        legend.margin=margin(1,2,1,1),
        legend.box.margin=margin(-25,0,-10,0), 
        #legend.box.margin=margin(-30,0,-15,0), 
        # axis.title.y.left = element_text(margin = margin(r = 5))
        )

# ggsave("./Figures/06/gene_dotplot.png", B_1, width = 6.7, height = 2.4) #, dpi = 300

############################
# CLUS EXPRESION ON TISSUE #
############################
DAT <- DATA %>%
    filter(grepl("P118|P097", .$orig.ident))

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

p_list <- map(p, ~ .x + theme(plot.margin = unit(c(-25, -10, -10, -15), "pt"), #t,r,b,l
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
# dev.new(width=6.7, height=3.4, noRStudioGD = TRUE) 
c <- plot_grid( plotlist = p_list, ncol = 4, rel_heights = c(1), rel_widths = c(1))
B_2 <- plot_grid(NULL, c, rel_widths = c(.04,1))
(B <- plot_grid( B_1, B_2, ncol=1, rel_heights = c(3,1.1), labels = c('B'), hjust = 0))
```

<img src="../Figures/06/05b_dotplot_and_tissue_plot.png"
data-fig-align="center" />

``` r
# ggsave("./Figures/06/dotplot_and_tissue_plot.png", B, width = 6.7, height = 3.2, bg = "white", dpi = 300)
```

``` r
#############################
# COMBINE ALL FIGURE PANELS #
#############################
# dev.new(width=6.7, height=6.7, noRStudioGD = TRUE)
Figure6 <- plot_grid( A_C, B, ncol=1, rel_heights = c(1,1)) 
ggsave("./Figures/Figure-6.pdf", Figure6, width = 6.7, height = 6.7, bg = "white", dpi = 1000)

Figure6
```

<img src="../Figures/06/Figure-6.png" data-fig-align="center" />
