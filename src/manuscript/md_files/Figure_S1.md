Supplemental Figure 2
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
library(scatterpie)
library(patchwork)
library(openxlsx)

source("../../bin/spatial_visualization.R")
source("../../bin/plotting_functions.R")

#########
# PATHS #
#########
input_dir <- "../../results/04_deconvolution_st_data/"
result_dir <- "../../results/09_figures/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
epi_clus <- "^5$|^6$|^7|^9" # non-filt

#############
# LOAD DATA #
#############
# DEGs_table <- read_csv(paste0(input_dir,"subset_100/DGEs_condition_wilcox_epi_100.csv")) %>% filter(groups == "DMPA")
DEGs_table <- read_csv(paste0("../../results/05_DGE_clusters_st_data/","DGEs_clusters_wilcox.csv"))
DEGs_table_EvsS <- read_csv(paste0("../../results/05_DGE_clusters_st_data/","DGEs_clusters_morf.csv"))
DATA <- readRDS(paste0(input_dir,"seuratObj_deconvolution_scdc.RDS"))
DATA_r <- readRDS(paste0("../../results/06_plot_annotation_ref_data/","seuratObj_annotated_celltypes.RDS"))
```

``` r
DATA@misc$cell_annot %>% 
  mutate(cell_annot_1 = ifelse(.$cell_annot_1 == "Epithelial", "Keratinocyte", .$cell_annot_1)) %>%
  {. ->> DATA@misc$cell_annot }
```

``` r
colors <- c( "#CD9600", "#7CAE00", "#e0e067", "#00A9FF", "#377EB8","#984EA3", "#E41A1C", "#C77CFF",
             "#00BFC4", "#FF7F00","#FFFF33")
clus_n <- 0:10
colors <- set_names(colors, clus_n)

clus <- DATA %>%
  #mutate(Clusters = as.character(.$Clusters)) %>%
  plot_spatial.fun( .,
                 geneid = "Clusters",
                 sampleid = c("P118", "P080","P031", "P105", "P097","P108", "P114", "P107"),
                 save_space = T,
                 lab = F,
                 zoom = "zoom",
                 col = colors,
                 alpha = .9,
                 ncol = 4, 
                 #annot_col = "#dbd9d9",
                 annot_line = .1,
                 img_alpha = 0,
                 point_size = .6) + 
  theme(legend.box.margin=margin(0,30,-0,-30)) # moves the legend)
  
clus
```

<img src="../Figures/S1/S1c_clusters-on-tissue.png"
data-fig-align="center" />

``` r
# dev.new(height=3.3, width=7, noRStudioGD = TRUE)
# ggsave("./Figures/S1/clusters_on_tissue.pdf", clus, width = 7, height = 3.3, dpi = 500, bg = "white")
```

### Plot deconvolution cell proportion on tissue

``` r
#################################
# CELL TYPE PROPORTION PER SPOT #
#################################
col <- c("#FFD92F","#8DA0CB","#FC8D62","#66C2A5","#E78AC3","#A6D854","#377EB8","#eb6062","#4DAF4A","#984EA3","#B3B3B3","#E5C494","#FF7F00","#FFFF33","#A65628","#F781BF") # "#E5C494"
cell_type <- c("B cell", "T cell", "Keratinocyte","Endothelial","Fibroblast", "Myeloid", "Granulocyte", "NK cells","Lymphatics", "ILC", "Unknown") # "Smooth muscle"
cell_col <- set_names(col[1:length(cell_type)], cell_type)

DefaultAssay(DATA) <- "RNA"

(plot <- DATA %>%
  plot_cell_pie.fun(.,
                   assay = "RNA",
                   ct.res = "cell_annot_1",
                   ct.select = cell_type, 
                   radius_adj = 3,
                   zoom = "zoom",
                   col = cell_col,
                   alpha = 1,
                   ncol = 4,
                   #annot_col = "#dbd9d9",
                   annot_line = .2,
                   img_alpha = 0)  + 
   guides(fill = guide_legend(keyheight = .7, keywidth = .7)) +
   theme(legend.box.margin=margin(0,25,-0,-20),
         legend.title = element_blank()) )# moves the legend)
```

<img src="../Figures/S1/S1d_cell_prop_tissue.png"
data-fig-align="center" />

``` r
# dev.new(height=3.3, width=7, noRStudioGD = TRUE)
# ggsave("./Figures/S1/cell_prop_tissue.pdf", plot, width = 7, height = 3.3, dpi = 500, bg = "white")
```
