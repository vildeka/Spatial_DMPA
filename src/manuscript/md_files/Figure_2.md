Figure 2
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

# remotes::install_github('linxihui/NNLM')
library(NNLM)
library(parallel)

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
#######################
# DEFINE MARKER GENES #
#######################
Cell_marker <- c("MS4A1","CD79A", "POU2AF1", "MZB1", "FCRL5", "CD38","CD247","CD8A","CD3D","CD3G","CD4","CD28",  
                  "NCAM1", "GZMA", "GNLY",  "FCER1A", "CD1A", "CD1C",  
                 "COL3A1","COL1A1", "MYH11","LAMA2", "CDH5","EMCN","FLT1","SELE",   
                  "ATG9B", "KPRP","ALOX12", "PRSS3", "BICDL2", "MT1X")
names <- c("B cell", "B cell","B cell", "Plasma cell","Plasma cell", "Plasma cell", "T cell", "T cell",
           "T cell", "T cell", "T cell", "T cell", "NK cells", "NK cells", "NK cells", 
           "Myeloid", "Myeloid", "Myeloid", "Fibroblast", "Fibroblast", "Fibroblast","Fibroblast", 
           "Endothelial","Endothelial", "Endothelial", "Endothelial", "Keratinocyte supra", "Keratinocyte supra", 
           "Keratinocyte supra","Keratinocyte supra", "Keratinocyte basal", "Keratinocyte basal")
cell_type <- set_names(names, Cell_marker)
clus_lvl <- rev(c("6", "9", "7", "5","8","3","4","0","1","2","10")) 

clus_lvl <- set_names(seq_along(clus_lvl), clus_lvl)
gene_lvl <- set_names(seq_along(Cell_marker), Cell_marker)
```

``` r
################
# FETCH GENES #
################
cell_type <- set_names(names, Cell_marker)
df <- DATA %>%
  mutate(., FetchData(., vars = c(Cell_marker)) ) %>%
  as_tibble() %>%
  select(., .cell, Clusters, any_of(Cell_marker)) %>%
  pivot_longer(., cols = -c(".cell", "Clusters"), 
               names_to = "marker", values_to = "values") %>%
  mutate(., cell = cell_type[as.character(.$marker)]) %>%
  mutate( marker_id = paste0(.$cell," (",.$marker,")")) %>%
  mutate(Clusters = factor(.$Clusters, levels=names(clus_lvl))) %>%
  mutate(marker = factor(.$marker, levels=Cell_marker)) %>%
  group_by(Clusters, marker, cell) %>%
  summarise(Avg = mean(values),
            Pct = sum(values > 0) / length(values) * 100, .groups="drop") %>%
  mutate(., ymin = clus_lvl[as.character(.$Clusters)]-0.5,
            ymax = clus_lvl[as.character(.$Clusters)]+0.5) %>%
  mutate(., xmin = gene_lvl[as.character(.$marker)]-0.5,
            xmax = gene_lvl[as.character(.$marker)]+0.5)

########################
# MARKER GENES DOTPLOT #
########################
library(ggnewscale)
cell_col <- c("#FFD92F","#FFFFCC","#8DA0CB","#eb6062","#A6D854","#E78AC3","#66C2A5","#FC8D62","#FED9A6",
         "#377EB8","#4DAF4A","#B3B3B3","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")

clus_col <- rev(c("#FBAAB1","#FFDAB8","#F1D2FF","#E6B3E9","#8BCFCE","#92DBFF","#ABC9E1", "#E5C264","#BBD99B", "#E1E2A4", "#FDFFD3"))
lab <- c("B cell", "Plasma cell", "T cell", "NK cells", "Myeloid", "Fibroblast", 
"Endothelial", "Keratinocytes")

ym <- max(df$ymax)
xm <- length(cell_type)+3.5

clus_lvl <- rev(c("6", "9", "7", "5","8","3","4","0","1","2","10")) 
# DATA %>%
#   mutate(Clusters = factor(.$Clusters, levels = clus_lvl)) %>%
# ggplot(., aes(x=Clusters, y=nCount_RNA, fill=Clusters)) + geom_violin() + scale_fill_manual(values = clus_col)

(A <- ggplot(df, aes(x=marker, y=Clusters)) +
  geom_point(aes(size = Pct, fill = Avg), color="white", shape=21) +
  scale_fill_gradientn(colours = viridisLite::magma(100),
                       guide = guide_colorbar(ticks.colour = "white",
                                              frame.colour = "white",
                                              barwidth = .5, barheight = 4),
                       name = "Average\nexpression") +
  #facet_grid(~ cell, scales = "free_x") +
  # Cell type colour bar
  new_scale_fill() +
  geom_rect(aes(ymin=max(ymax), ymax=max(ymax)+.4,
                xmin=xmin, xmax=xmax,
                fill = cell),data=df,alpha = 0.1,show.legend=F) +
  annotate("text", x = c(1.8, 5, 10, 14, 17, 20.5, 24.6, 29.5 ), y = 12.2, label = lab) +
  geom_rect(aes(ymin=min(ymin), ymax=min(ymin)+.3,
                xmin=xmin, xmax=xmax,
                fill = cell),data=df,alpha = 0.1,show.legend=F) +
  scale_fill_manual(values = set_names(cell_col[1:length(unique(cell_type))], unique(cell_type))) +
    
  # add extra borders for the colour bars
  # annotate(x = c(.5,xm,-.5,-.5,.5,-.5), xend=c(.5,xm,xm,xm,xm,-.5),
  #            y=c(1,1,.5,ym,ym+.5,.5), yend=c(ym+.5,ym+.5,.5,ym,ym+.5,ym),
  #            geom="segment",colour = "black", linewidth = .6, alpha = 1) +
  # Cluster bar annnotation
  new_scale_fill() +
  geom_rect(aes(ymin=ymin, ymax=ymax,
                xmin=0.4, xmax=-.5,
                fill = Clusters),data=df,alpha = 1,show.legend=F) +  
  scale_fill_manual(values = clus_col) +
  coord_cartesian(clip="off", xlim=c(.5,xm),ylim=c(.5,ym),expand = F) +
  scale_size("% detected", range = c(0,6)) +

  ylab("Cluster") + xlab("") +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(color = "black"), keywidth = .6, keyheight = .6)) +
  theme(axis.text.x = element_text(size=8, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=9, color="black",hjust=.9 ),
        axis.title = element_text(size=9),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 9),
        legend.margin=margin(0,0,-0,0),
        plot.margin = unit(c(1,0,-.5,0),units = "cm") #trbl
        ) 
)
```

<img src="../Figures/02/02a_marker-gene-dotplot.png"
data-fig-align="center" />

``` r
# dev.new(width=8, height=5, noRStudioGD = TRUE)
# ggsave("./Figures/04/marker-gene-dotplot.pdf", A, width = 8, height = 5) # cell_marker
```

``` r
####################################
# COLOUR PALLET SUBGROUPS FUNCTION #
####################################
# maximum levels in the subgroup is 16
s_col <- c("#FFD92F","#8DA0CB","#A6D854","#eb6062","#FC8D62","#E78AC3","#66C2A5","#B3B3B3","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#E5C494")
e_col <- c("#FFF2AE","#CBD5E8","#E6F5C9","#FBB4AE","#FDCDAC","#F4CAE4","#B3E2CD","#CCCCCC","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F1E2CC")

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

``` r
################
# FETCH GENES #
################
cell_type <- set_names(names, Cell_marker)
df <- DATA %>%
  mutate(., FetchData(., vars = c(Cell_marker)) ) %>%
  as_tibble() %>%
  select(., .cell, Clusters, any_of(Cell_marker)) %>%
  pivot_longer(., cols = -c(".cell", "Clusters"), 
               names_to = "marker", values_to = "values") %>%
  mutate(., cell_type = cell_type[as.character(.$marker)]) %>%
  mutate( marker_id = paste0(.$cell_type," (",.$marker,")")) %>%
  mutate(Clusters = factor(.$Clusters, levels=names(clus_lvl))) %>%
  mutate(marker = factor(.$marker, levels=Cell_marker)) 

###########################
# COLOUR PALLET SUBGROUPS #
###########################
#s_col <- c("#FFD92F","#66C2A5","#A6D854","#eb6062","#377EB8","#4DAF4A","#B3B3B3","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")
#e_col <- c("#FFF2AE","#B3E2CD","#E6F5C9","#FBB4AE","#B3CDE3","#CCEBC5","#CCCCCC","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC")

meta_group <- "cell_type"
subgroup <- "marker_id"
meta_group_lvl <- c("B cell","T cell","Myeloid","NK cells","Keratinocyte","Fibroblast","Endothelial") # all

col <- df %>%
  filter(grepl(paste0(meta_group_lvl,collapse = "|"),.$cell_type)) %>%
  mutate(cell_type = ifelse( grepl("Keratinocyte", .$cell_type), "Keratinocyte", .$cell_type)) %>%
  mutate(meta_group = factor(.[[meta_group]], levels = meta_group_lvl)) %>%
  mutate(annot = .[[subgroup]]) %>%
  arrange(meta_group) %>%
  mutate(annot = factor(.[[subgroup]], levels = unique(.[[subgroup]]))) %>%
  {. ->> df_} %>%
  ColourPalleteMulti(., "meta_group", "annot")

cell_type <- levels(df_$annot)
cell_col <- set_names(col[1:length(cell_type)], cell_type)
# DATA@misc$markers <- list(cell_annot = df_)  

###################################
# CELL MARKER PROPORTION PER SPOT #
###################################
# dev.new(width=4, height=6, noRStudioGD = TRUE)
immune <- cell_type[grepl(paste0(c("B cell","T cell","Myeloid","NK cells"),collapse = "|"), cell_type)] # immune
select <- list(NULL, immune)

DAT <- DATA %>% filter(orig.ident == "P097")
DAT@misc$cell_annot <- df_ 

plot <- map(select,
          ~plot_cell_pie.fun(DAT,
                           ct.res = "marker_id",
                           #ct.res = cell_type,
                           ct.select = .x, 
                           radius_adj = -2,
                           zoom = "zoom",
                           colors = cell_col,
                           alpha = 1,
                           ncol = 2,
                           #annot_col = "#dbd9d9",
                           annot_line = .1,
                           img_alpha = 0) +
            coord_equal() +
            theme(legend.position = "none",
                  panel.spacing.x = unit(-1, "lines"),
                  panel.spacing.y = unit(-3, "lines"),
                  plot.margin = unit(c(-2.5, -3.5, -2.5, -4.2), "lines") # t,r,b,l
                  ) )

p_1 <- plot[[1]] + theme(plot.margin = unit(c(-1.3, -3.5, .5, -4), "lines"))
p_2 <- plot[[2]] + theme(plot.margin = unit(c(.5, -3.5, -1.3, -4.2), "lines"))

#p <- plot_grid(plot[[1]], plot[[2]], ncol = 2)
# dev.new(width=5.5, height=2.4, noRStudioGD = TRUE)
(B_ <- plot_grid(NULL, p_1, NULL, p_2, nrow = 1, rel_widths = c(.15,1.1,.1,1.1)))
```

<img src="../Figures/02/02b_cell-marker-prop-on-tissue.png"
data-fig-align="center" />

``` r
col <- c("#FFD92F","#8DA0CB","#FC8D62","#66C2A5","#E78AC3","#A6D854","#eb6062","#B3B3B3","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")
cell_type <- c("B cell", "T cell", "Keratinocyte","Endothelial","Fibroblast", "Myeloid", "NK cells")
cell_col <- set_names(col[1:length(cell_type)], cell_type)

# create freestanding legend
dat <- data.frame( X = rnorm(length(cell_type)), c = col[1:length(cell_type)], b=factor(cell_type, levels = cell_type)  )
p <- ggplot(dat, aes(X,  fill = b)) + scale_fill_manual(values = cell_col) +  geom_bar()
(p_ <- p + guides(fill = guide_legend(keywidth = .6, keyheight = .6) ) + 
    theme(legend.title = element_blank(),
          legend.margin=margin(0,0,-0,0),
          legend.background = element_rect(fill='transparent')
          ) )
```

<img src="../Figures/02/02b_legend.png" data-fig-align="center" />

``` r
# combine plots and legend
legend <- get_legend(p_)
B <- B_ + draw_grob(legend, hjust = -.09, vjust = -.28)

# dev.new(width=7.4, height=2.4, noRStudioGD = TRUE)
B_C <- plot_grid(B, NULL, rel_widths = c(1,.6),  labels = c('B','C'))
# ggsave("./Figures/02/Figure02B.png", B_C, width = 7.4, height = 3, bg = "white", dpi = 400)
```

## Semla deconvolution

``` r
###################
# SEMLA FUNCTIONS #
###################
FactorGeneLoadingPlot <- function (
  object,
  factor = 1,
  topn = 20,
  dark.theme = FALSE
) {

  # require NNLM
  if (!requireNamespace("NNLM")) {
    #devtools::install_github('linxihui/NNLM')
  }
  
  # Check if NMF has been computed
  if (!"NMF" %in% names(object@reductions)) stop("NMF has not been computed ... \n", call. = FALSE)

  ftr <- paste0("factor_", factor)
  nmf <- object@reductions$NMF@feature.loadings[, ftr]
  gene <- names(nmf)
  df <- data.frame(gene, val = nmf, stringsAsFactors = F)
  df <- df[order(df$val, decreasing = T), ]
  df <- df[1:topn, ]
  df$gene <- factor(df$gene, levels = df$gene)
  p <- ggplot(df[1:topn, ], aes(reorder(gene, val), val)) +
    geom_bar(stat = "identity", fill = ifelse(dark.theme, "dark gray", "lightgray"), color = ifelse(dark.theme, "lightgray", "black"), width = 0.7) +
    coord_flip() +
    labs(x = "gene", y = "value")

  if (dark.theme) p <- p + DarkTheme() else p <- p + theme_minimal()

  return(p)
}


#' Run Non-negative Matrix Factorization
#'
#' Decompose an expression matrix A with non-negative elements into matrices WxH, also with
#' non-negative elements. W is the feature loading matrix (features x factors) and H is the
#' low dimensional embedding of the spots (factors x spots).
#'
#' @param object Seurat object
#' @param assay Assay Name of Assay NMF is being run on
#' @param slot Slot to pull data from.
#' @param features Features to compute the NMF for. Note that these features must be present in the
#' slot used to compute the NMF. By default, the `features` is set to `VariableFeatures(object)`
#' to include the most variable features selected in the normalization step.
#' @param nfactors Total Number of factors to compute and store (20 by default)
#' @param rescale Rescale data to make sure that values of the input matrix are non-n
#' @param reduction.name Dimensional reduction name, "NMF" by default
#' @param reduction.key Dimensional reduction key, specifies the prefix of the factor ids, e.g.
#' "factor_1", "factor_2", etc.
#' @param n.cores Number of threads to use in computation
#' @param order.by.spcor Order factors by spatial correlation
#' @param sort.spcor.by.var Sort factors by decreasing variance
#' @param ... Additional parameters
#'
#' @importFrom parallel detectCores
#' @importFrom Seurat CreateDimReducObject DefaultAssay VariableFeatures GetAssayData
#'
#' @export
#'
RunNMF <- function (
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  nfactors = 20,
  rescale = TRUE,
  reduction.name = "NMF",
  reduction.key = "factor_",
  n.cores = NULL,
  order.by.spcor = FALSE,
  sort.spcor.by.var = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  var.genes <- features %||% VariableFeatures(object)
  norm.counts <- GetAssayData(object, slot = slot, assay = assay)
  if (rescale) {
    norm.counts <- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }
  if (min(norm.counts) < 0) stop("Negative values are not allowed")
  nmf.results <- rnmf(A = norm.counts[var.genes, ], k = nfactors)
  #nmf.results$W <- swne::ProjectFeatures(norm.counts, nmf.results$H, n.cores = n.cores)
  feature.loadings <- nmf.results$W
  cell.embeddings <- t(nmf.results$H)

  # Set cores
  n.cores <- n.cores %||% {detectCores() - 1}

  # Order factors based on spatial correlation
  if (order.by.spcor) {
    CN <- do.call(rbind, GetSpatNet(object = object, nNeighbours = NULL, maxdist = NULL))
    resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
    resCN[resCN > 0] <- 1
    empty.CN <- matrix(0, nrow = nrow(cell.embeddings), ncol = nrow(cell.embeddings), dimnames = list(rownames(cell.embeddings), rownames(cell.embeddings)))
    colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
    colnames(resCN) <- gsub(pattern = "^X", replacement = "", x = colnames(resCN))
    empty.CN[rownames(resCN), colnames(resCN)] <- resCN
    listw <- spdep::mat2listw(empty.CN)
    fun <- function (x) spdep::lag.listw(listw, x, TRUE)

    # Calculate the lag matrix from the network
    tablag <- apply(cell.embeddings, 2, fun)

    # Split sp.cor by sample
    if (sort.spcor.by.var) {
      sp.cor.split <- do.call(rbind, lapply(unique(GetStaffli(object)@meta.data$sample), function(s) {
        tablag.split <- tablag[GetStaffli(object)@meta.data$sample == s, ]
        cell.embeddings.split <- cell.embeddings[GetStaffli(object)@meta.data$sample == s, ]
        unlist(lapply(1:ncol(cell.embeddings.split), function(i) {
          cor(tablag.split[, i], cell.embeddings.split[, i])
        }))
      }))
      order.vec <- order(apply(sp.cor.split, 2, var))
    } else {
      sp.cor <- unlist(lapply(1:ncol(cell.embeddings), function(i) {
        cor(cell.embeddings[, i], tablag[, i])
      }))
      order.vec <- order(sp.cor, decreasing = TRUE)
    }

    cell.embeddings <- cell.embeddings[, order.vec]
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
  }

  rownames(x = feature.loadings) <- var.genes
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nfactors)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject (
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}


#' Run NMF with ica init
#' 
#' @param A Expression matrix
#' @param k = Number of factors
#' @param alpha Alpha parameter value
#' @param init Method to initiate NMF by
#' @param n.cores Number of threads to run
#' @param loss Loss method
#' @param max.iter Maximum number of iterations for NMF
#' @param ica.fast Should fast ica be run?
#' 
#' @importFrom stats runif
#' 
rnmf <- function (
  A,
  k,
  alpha = 0,
  init = "ica",
  n.cores = 1,
  loss = "mse",
  max.iter = 500,
  ica.fast = F
) {
  if (any(A < 0))
    stop("The input matrix contains negative elements !")
  if (k < 3)
    stop("k must be greater than or equal to 3 to create a viable SWNE plot")
  if (!init %in% c("ica", "nnsvd", "random")) {
    stop("Invalid initialization method")
  }
  A <- as.matrix(A)
  if (any(A < 0)) {
    stop("Input matrix has negative values")
  }
  if (init == "ica") {
    nmf.init <- ica_init(A, k, ica.fast = ica.fast)
  }
  else if (init == "nnsvd") {
    nmf.init <- nnsvd_init(A, k, LINPACK = T)
  }
  else {
    nmf.init <- NULL
  }
  if (is.null(nmf.init)) {
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, n.threads = n.cores,
                          loss = loss, max.iter = max.iter)
  }
  else {
    A.mean <- mean(A)
    zero.eps <- 1e-06
    nmf.init$W[nmf.init$W < zero.eps] <- 0
    nmf.init$H[nmf.init$H < zero.eps] <- 0
    zero.idx.w <- which(nmf.init$W == 0)
    zero.idx.h <- which(nmf.init$H == 0)
    nmf.init$W[zero.idx.w] <- runif(length(zero.idx.w), 0,
                                    A.mean/100)
    nmf.init$H[zero.idx.h] <- runif(length(zero.idx.h), 0,
                                    A.mean/100)
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, init = nmf.init,
                          n.threads = n.cores, loss = loss, max.iter = max.iter)
  }
  colnames(nmf.res$W) <- rownames(nmf.res$H) <- sapply(1:ncol(nmf.res$W),
                                                       function(i) paste("factor", i, sep = "_"))
  return(nmf.res)
}



#' Summarize features associated with cselected factors
#'
#' Extracts the top driving features per factor and returns
#'
#' @param object Seurat object
#' @param dims Factors to use
#' @param features.return Number of features to return per factor
#' @param features.use Select features (genes) to subset the data on
#'
#' @export
#'
SummarizeAssocFeatures <- function (
  object,
  dims = NULL,
  features.return = 10,
  features.use = NULL
) {

  if (!"NMF" %in% names(object@reductions)) stop(paste0("No factors available in Seurat object. Run RunNMF() first "), call. = FALSE)

  feature.factor.assoc <- object@reductions[["NMF"]]@feature.loadings
  if (!is.null(features.use)) {
    feature.factor.assoc <- feature.factor.assoc[features.use, ]
  }
  if (!is.null(dims)) {
    feature.factor.assoc <- feature.factor.assoc[, dims]
  }
  factor.features.df <- do.call("rbind", lapply(1:ncol(feature.factor.assoc),
                                                function(i) {
                                                  features.df <- data.frame(assoc_score = feature.factor.assoc[, i])
                                                  features.df$feature <- rownames(feature.factor.assoc)
                                                  features.df$factor <- colnames(feature.factor.assoc)[[i]]
                                                  features.df <- features.df[order(features.df$assoc_score, decreasing = T), ]
                                                  head(features.df, n = features.return)
                                                }))
  rownames(factor.features.df) <- NULL
  gene.loadings.selected <- feature.factor.assoc[unique(factor.features.df$feature), ]
  return(list(factor.features.df, gene.loadings.selected))
}


# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Initiate NMF using ICA
#'
#' @param A input matrix
#' @param k number of components to compute
#' @param ica.fast Should a fast implementation of ICA be used?
#'

ica_init <- function (
    A, 
    k, 
    ica.fast = F
) {
  
  if (!requireNamespace("irlba")) install.packages("irlba")
  if (!requireNamespace("ica")) install.packages("ica")
  
  if (ica.fast) {
    pc.res.h <- irlba::irlba(t(A), nv = 50, maxit = 100,
                             center = rowMeans(A))
    ica.res.h <- ica::icafast(pc.res.h$u, nc = k, maxit = 25,
                              tol = 1e-04)
    return(list(W = (A - Matrix::rowMeans(A)) %*% ica.res.h$S,
                H = t(ica.res.h$S)))
  }
  else {
    ica.res <- ica::icafast(t(A), nc = k, maxit = 25, tol = 1e-04)
    return(list(W = ica.res$M, H = t(ica.res$S)))
  }
}

nnsvd_init <- function (A, k, LINPACK)
{
  size <- dim(A)
  m <- size[1]
  n <- size[2]
  W <- matrix(0, m, k)
  H <- matrix(0, k, n)
  s = svd(A, k, k, LINPACK = LINPACK)
  U <- s$u
  S <- s$d
  V <- s$v
  W[, 1] = sqrt(S[1]) * abs(U[, 1])
  H[1, ] = sqrt(S[1]) * abs(t(V[, 1]))
  for (i in seq(2, k)) {
    uu = U[, i]
    vv = V[, i]
    uup = .pos(uu)
    uun = .neg(uu)
    vvp = .pos(vv)
    vvn = .neg(vv)
    n_uup = .norm(uup)
    n_vvp = .norm(vvp)
    n_uun = .norm(uun)
    n_vvn = .norm(vvn)
    termp = n_uup %*% n_vvp
    termn = n_uun %*% n_vvn
    if (termp >= termn) {
      W[, i] = sqrt(S[i] * termp) * uup/n_uup
      H[i, ] = sqrt(S[i] * termp) * vvp/n_vvp
    }
    else {
      W[, i] = sqrt(S[i] * termn) * uun/n_uun
      H[i, ] = sqrt(S[i] * termn) * vvn/n_vvn
    }
  }
  return(list(W = W, H = H))
}



  
# plot_genes.fun(DATA_st, gene="factor_1" ,red = "UMAP")
# 
# ID <- c("P107", "P108", "P114", "P097","P118", "P105", "P080", "P031")
# ID <- c("P107")
# plot_spatial.fun(DATA_st, geneid = c("factor_1") , 
#                  zoom = "zoom", 
#                  assay = "celltypeprops", 
#                  sampleid = ID,
#                  alpha = 1, img_alpha = 0)
# 
# plot_st_feat.fun(DATA_st, geneid = c("factor_1") , zoom = "zoom", assay = "celltypeprops" )
# 
# plot_clusters.fun(DATA_nest$data[[2]], cluster=.x, txt_size = 10, dot_size = 0.2,
#                              color = .y, red = "harmony") + xlab("UMAP 1") + ylab("UMAP 2")
# 
# DimPlot(., group.by = "factor_1", reduction = "umap")
```

``` r
#################################
# SEMLA DECONVOLUTION FUNCTIONS #
#################################
RunNNLS <- function (
  object,
  singlecell_matrix,
  groups,
  nCells_per_group = 50L,
  minCells_per_celltype = 10L,
  min_prop = 0.01,
  L1 = 0.01,
  seed = 1337L,
  return_expression_profiles = FALSE,
  verbose = TRUE,
  ...
) {

  # Check input objects
  if (!inherits(object, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    rlang::abort()(glue::glue("Invalid format of spatial expression matrix: '{class(object)}'"))
  if (any(dim(object) == 0))
    rlang::abort()(glue::glue("Invalid dimensions of spatial expression matrix: '{paste(dim(object), collapse = 'x')}'"))
  if (!inherits(singlecell_matrix, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    rlang::abort()(glue::glue("Invalid format of single-cell expression matrix: '{class(singlecell_matrix)}'"))
  if (any(dim(singlecell_matrix) == 0))
    rlang::abort()(glue::glue("Invalid dimensions of single-cell expression matrix matrix: '{paste(dim(singlecell_matrix), collapse = 'x')}'"))


  # Check that groups are valid
  if (!inherits(groups, what = "character"))
    rlang::abort()(glue::glue("Invalid class '{class(groups)}' for groups, expected a 'character' vector"))
  if (length(groups) == 0)
    rlang::abort()("'groups' is empty")
  if (length(groups) != ncol(singlecell_matrix))
    rlang::abort()(glue::glue("'groups' does not match the single-cell expression matrix. ",
               "Expected {ncol(singlecell_matrix)}, got {length(groups)}"))


  # Make sure that dev version of RcppML is installed
  if (!requireNamespace("RcppML"))
    rlang::abort()("Package 'RcppML' is required.") # compile dev version
  if ((packageVersion("RcppML") |> as.character()) != "0.3.7")
    rlang::warn(c("The NNLS function might break if using a dev version of RcppML is used. ",
           "If RcppML::project(...) fails, try installing CRAN version 0.3.7 of RcppML."))

  # Prepare data
  if (verbose) cli::cli_alert_info("Preparing data for NNLS")
  W <- .get_w_matrix(object = object,
                     singlecell_matrix = singlecell_matrix,
                     nCells_per_group = nCells_per_group,
                     minCells_per_celltype = minCells_per_celltype,
                     groups = groups,
                     seed = seed,
                     verbose = verbose)

  # Run NNLS
  if (verbose) cli::cli_alert_info("Predicting cell type proportions with NNLS for {ncol(W)} cell types")
  proj_expr <- try({RcppML::project(object, W, L1 = L1, ...)}, silent = TRUE)
  
  if (inherits(proj_expr, what = "try-error")) {
    proj_expr <- RcppML::project(data = object, w = W, L1 = L1, ...)
  }

  # Convert predicted values to proportions
  prop <- apply(proj_expr, 2, function(x) {prop.table(x)})
  prop[prop < min_prop] <- 0
  prop <- apply(prop, 2, function(x) {prop.table(x)})
  rownames(prop) <- colnames(W)
  colnames(prop) <- colnames(object)

  if (return_expression_profiles) {
    return(list(prop = prop, W = W))
  } else {
    return(prop)
  }
}
RunNNLS.Seurat <- function (
    object,
    singlecell_object,
    groups = NULL,
    features = NULL,
    singlecell_assay = "RNA",
    spatial_assay = "Spatial",
    slot = "data",
    min_prop = 0.01,
    nCells_per_group = 50L,
    minCells_per_celltype = 10L,
    assay_name = "celltypeprops",
    dimred_name = "nnls",
    return_as_dimred = FALSE,
    L1 = 0.01,
    seed = 1337L,
    verbose = TRUE,
    ...
) {

  # Validate input objects
  stopifnot(inherits(object, what = "Seurat"),
            inherits(singlecell_object, what = "Seurat"))

  if (verbose) cli::cli_h2("Predicting cell type proportions")
  if (verbose) cli::cli_alert_info("Fetching data from Seurat objects")

  # Get features
  features <- .prep_features(object, singlecell_object, features, verbose)

  # Get groups
  groups <- .prep_groups(singlecell_object, groups, verbose)

  # Get expression matrix for single cell data and Visium data
  x_singlecell <- GetAssayData(singlecell_object, slot = slot, assay = singlecell_assay)[features, ]
  x_spatial <- GetAssayData(object, slot = slot, assay = spatial_assay)[features, ]
  # TODO: check if RNA data slot is normalized if slot = "data"

  results <- RunNNLS(
    object = x_spatial,
    singlecell_matrix = x_singlecell,
    groups = groups,
    nCells_per_group = nCells_per_group,
    minCells_per_celltype = minCells_per_celltype,
    min_prop = min_prop,
    L1 = L1,
    seed = seed,
    verbose = verbose,
    return_expression_profiles = TRUE,
    ... = ...
  )

  # Return results as a DimReduc object or an Assay object
  object <- .return_as(
    st_object = object,
    prop = results$prop,
    W = results$W,
    return_as_dimred = return_as_dimred,
    st_assay = spatial_assay,
    assay_name = assay_name,
    dimred_name = dimred_name,
    verbose = verbose
  )

  if (verbose) cli::cli_alert_success("Finished")

  return(object)
}


#' Utility function to format returned data
#'
#' @param st_object A \code{Seurat} object woth Visium data
#' @param prop A matrix with cell type proportion estimates
#' @param W A matrix with cell type expression profiles
#' @param return_as_dimred Logical specifying if the data should be returned
#' as a \code{DimReduc} object or an \code{Assay} object
#' @param st_assay Assay used for spatial data
#' @param assay_name Assay name for returned data. Only used if \code{return_as_dimred=FALSE}
#' @param dimred_name Name for \code{DimReduc} object. Only used if \code{return_as_dimred=TRUE}
#' @param dimred_prefix Prefix for \code{DimReduc} object vectors.
#' @param verbose Print messages
#'
#' @import glue
#' @import cli
#' @importFrom Seurat CreateDimReducObject CreateAssayObject DefaultAssay `DefaultAssay<-`
#'
#' @noRd
.return_as <- function (
    st_object,
    prop,
    W,
    return_as_dimred,
    st_assay,
    assay_name,
    dimred_name,
    dimred_prefix = "cFactor_",
    verbose
) {
  # Return results as a DimReduc object or an Assay object
  if (return_as_dimred) {
    if (verbose) cli::cli_alert_info("Returning results as a 'DimReduc' object")
    prop <- t(prop)
    colnames(prop) <- paste0(dimred_prefix, 1:ncol(prop))
    feature_loadings <- W
    colnames(feature_loadings) <- paste0(dimred_prefix, 1:ncol(prop))
    props_dimreduc <- CreateDimReducObject(embeddings = prop,
                                           loadings = feature_loadings,
                                           key = dimred_prefix,
                                           assay = st_assay)
    st_object[[dimred_name]] <- props_dimreduc
  } else {
    if (verbose) cli::cli_alert_info("Returning results in a new 'Assay' named '{assay_name}'")
    props_assay <- CreateAssayObject(data = as(prop, "dgCMatrix"))
    st_object[[assay_name]] <- props_assay
    if (verbose) cli::cli_alert_info("Setting default assay to '{assay_name}'")
    #DefaultAssay(st_object) <- assay_name
  }

  return(st_object)
}


#' Prepare features to keep for NNLS
#'
#' @import cli
#' @import rlang
#' @import glue
#' @importFrom Seurat VariableFeatures
#'
#' @noRd
.prep_features <- function (
  st_object,
  sc_object,
  features,
  verbose
) {
  if (is.null(features)) {
    if (length(VariableFeatures(sc_object)) == 0) rlang::abort()("No variable features found in single cell Seurat object")
    if (length(VariableFeatures(st_object)) == 0) rlang::abort()("No variable features found in spatial Seurat object")
    features <- features %||% intersect(VariableFeatures(sc_object), VariableFeatures(st_object))
    if (length(features) < 1e3) rlang::warn(glue::glue("Only {length(features)} shared features detected"))
    if (length(features) < 100) rlang::abort()(glue::glue("At least 100 shared features required, found {length(features)}"))
  }

  # Filter features
  if (verbose) cli::cli_alert("  Filtering out features that are only present in one data set")
  keep_features <- intersect(intersect(rownames(sc_object), rownames(st_object)), features)
  if (verbose) cli::cli_alert("  Kept {length(keep_features)} features for deconvolution")

  return(keep_features)
}


#' @param sc_object A \code{Seurat} object with single-cell data
#' @param groups A character vector with group labels
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import cli
#' @import dplyr
#' @importFrom Seurat Idents
#'
#' @noRd
.prep_groups <- function (
  sc_object,
  groups,
  verbose
) {
  # Select groups
  if (is.null(groups)) {
    if (verbose) cli::cli_alert_info("Using current identity as groups.")
    groups <- as.character(Idents(sc_object))
  } else if (is.character(groups) & length(groups) == 1) {
    if (!groups %in% colnames(sc_object[[]])) rlang::abort()(glue::glue("'{groups}' is not a valid column in Seurat object meta.data slot"))
    groups <- sc_object[[]] |>  pull(all_of(groups))
  } else {
    rlang::abort()("Invalid value for 'groups'")
  }
  return(groups |> as.character())
}


#' Calculates expression profiles for groups of cells
#'
#' @param object An expression matrix with Visium data
#' @param singlecell_matrix An expression matrix with single-cell data
#' @param nCells_per_group Number of cells per cell type to keep
#' @param minCells_per_celltype Minimum number of cells accepted for a cell type
#' @param groups A character vector with group (cell type) labels
#' @param seed A seed for reproducibility
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import cli
#' @importFrom tibble tibble
#'
#' @noRd
.get_w_matrix <- function (
    object,
    singlecell_matrix,
    nCells_per_group,
    minCells_per_celltype,
    groups,
    seed,
    verbose
) {

  # Set global variables to NULL
  group <- barcode <- NULL

  # Rescale data to ensure positive values
  if (any(object < 0)) {
    if (verbose) rlang::warn("Found negative values in input matrix")
    if (verbose) cli::cli_alert_info("Rescaling data to ensure positive values")
    object <- t(apply(object, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }

  # Sample barcodes
  if (verbose) cli::cli_alert(glue::glue("  Downsampling scRNA-seq data to include a maximum of ",
                                   "{nCells_per_group} cells per cell type"))
  set.seed(seed)
  barcodes <- tibble(barcode = colnames(singlecell_matrix), group = groups) |>
    group_by(group) |>
    dplyr::slice(sample(min(nCells_per_group, n())))

  # Check that all celltypes have enough cells
  cells_n <- barcodes |>
    summarize(n = n())
  if (any(cells_n$n < minCells_per_celltype)) {
    too_small <- cells_n |>
      filter(n < minCells_per_celltype)
    if (verbose) cli::cli_alert(cli::col_br_magenta("  Cell type(s) {paste(too_small$group, collapse = ',')} ",
                           "have too few cells (<{minCells_per_celltype}) and will be excluded"))
    barcodes <- barcodes |>
      filter(!group %in% too_small$group)
  }
  if (verbose) cli::cli_alert("  Kept {length(unique(barcodes$group))} cell types after filtering")
  barcodes <- barcodes |>
    pull(barcode)
  new_groups <<- setNames(groups, nm = colnames(singlecell_matrix))[barcodes]
  new_groups <<- replace_na(new_groups, replace = "named_NA")

  # Compute cell type expression profiles
  x <<- singlecell_matrix[, barcodes]
  if (verbose) cli::cli_alert("  Calculating cell type expression profiles")
  W <- .get_expression_profiles(x = singlecell_matrix[, barcodes], groups =  new_groups)

  return(W)
}


#' Get single cell data expression profiles
#'
#' This function is used to obtain cell type expression profiles represented by
#' enrichment scores calculated for each cell.
#'
#' @param x Expression matrix to calculate expression profiles from
#' @param groups Character vector with group (cell type) labels
#'
#' @importFrom Matrix rowMeans
#' @importFrom tibble tibble
#' @import dplyr
#'
#' @noRd
.get_expression_profiles <- function (
    x,
    groups
){

  # Calculate means
  row_means <- do.call(cbind, lapply(unique(groups), function(grp) {
    Matrix::rowMeans(x[, groups == grp])
  }))
  colnames(row_means) <- unique(groups)

  # Calculate enrichment scores
  W <- do.call(bind_cols, lapply(colnames(row_means), function(grp) {
    x1 <- row_means[, grp]
    x2 <- row_means[, -which(colnames(row_means) == grp)]
    x2 <- rowMeans(x2) + 1
    y <- tibble(x1/x2) |> setNames(nm = grp)
    return(y)
  })) |>
    as.matrix()
  rownames(W) <- rownames(x)

  W <- apply(W, 2, function(x) {
    x/max(x)
  })

  return(W)
}
```

``` r
################
# MARKER GENES #
################
Cell_list <- list(
B_cell = c("MS4A1","LCN10","BANK1","BLK","TNFRSF13C","CD79A","MZB1"),# "CD19",
Plasma_cell = c("TNFRSF13B","FCRL5","TNFRSF17","CD38"), # "TLR10",
T_cell = c("CD247","CD8A","CD8B","CD3E","CD3D","CD3G","CD4","CD28"), #,"CD7"
NK_cell = c("KLRC1","NCAM1", "NCAM2", "GZMA","GNLY"), #ITGA1
DC_cell = c("FCER1A","CD1A","CD1C"),
Fibroblast = c("MYH11","LAMA2","COL3A1","COL1A1","FBN1","LUM"), # "OLFML3",
Endothelial = c("CDH5","EMCN","FLT1","SELE"), # "SELE","ACKR1","FAM110D", "EGFL7"
Epithelial_supra = c("KPRP","ALOX12", "ATG9B", "PRSS3"), # "KRT78", "LCE3E","CEACAM7","SLURP1", "SLURP2"
Epithelial_basal = c("MT1X", "BICDL2") #"COL17A1",
)
#############################
# PREPARE FOR DECONVOLUTION #
#############################
# identify common marker gnes between ST and SC datasets
markers_genes_annot1 <- readRDS(paste0("../../results/06_plot_annotation_ref_data/","markers_genes_annot1_filt.RDS")) 
markers_genes_annot1 <- readRDS(paste0("../../results/06_plot_annotation_ref_data/","markers_genes_annot1.RDS"))
st_genes <- rownames(DATA@assays$RNA@data)

markers_filt <- markers_genes_annot1 %>%
  bind_rows(., DEGs_table_EvsS) %>%
  filter(gene %in% st_genes) %>%
  filter(!(grepl("^IGK|IGH|JCHAIN", .$gene))) %>%
  filter(., !(ifelse(.$cluster == "Lymphatics" & pct.2 > .001, TRUE, FALSE))) %>%
  filter(., !(ifelse(.$cluster == "Endothelial" & pct.2 > .001, TRUE, FALSE))) %>%
  group_by(cluster) %>%
  filter(pct.1 >= 0.1 & pct.2 <= 0.1) %>%
  filter(log.pct.diff > 3) %>%
  nest() %>%
  filter(!(cluster == "Epithelial")) %>%
  ungroup() %>%
  mutate(cluster = ifelse(.$cluster == "epi_all", "Epithelial", .$cluster)) %>%
  unnest(cols = c(data)) 
  
var_features <- unique(c( unlist(Cell_list), markers_filt$gene ))

# Normalize data and find variable features for Visium data
DATA <- DATA %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 10000)

# Normalize data and run vanilla analysis to create UMAP embedding
DATA_r <- DATA_r %>%
  #filter(!(grepl("Smooth muscle|Lymphatics|NA", .$cell_annot_1))) %>%
  filter(!(is.na(DATA_r$cell_annot_1) | DATA_r$cell_annot_1=="Smooth muscle" | DATA_r$cell_annot_1=="Lymphatics")) %>% #| DATA_r$cell_annot_1=="Endothelial" )) |>
  #filter(!(grepl("Smooth muscle|Lymphatics|Endothelial", .$cell_annot_1))) |> 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# Alternative 1: Rerun FindVariableFeatures to increase the number before cell type deconvolution
# DATA_r <- DATA_r %>% 
#   FindVariableFeatures(nfeatures = 10000)

# Alternative 2: manually choose variable features
DefaultAssay(DATA) <- "RNA"
VariableFeatures(DATA_r) <- var_features
VariableFeatures(DATA) <- var_features

DimPlot(DATA_r, group.by = "cell_annot_1")
```

<img src="../Figures/02/prepare-for-semla-decon.png"
data-fig-align="center" />

``` r
##########################
# RUN NNLS DECONVOLUTION #
##########################
ti <- Sys.time()
# uses specified genes as variable features
# DATA_cell <- RunNNLS.Seurat(object = DATA,  spatial_assay = "RNA",
#                             features = int,
#                             nCells_per_group = 200,
#                             singlecell_object = DATA_r, 
#                             groups = "cell_annot_1")

# extracts variable features from the objects and uses intersect
DATA_cell <- RunNNLS.Seurat(object = DATA,  spatial_assay = "RNA",
                            # return_as_dimred = TRUE, dimred_name = "cell",
                            nCells_per_group = 200,
                            singlecell_object = DATA_r, 
                            groups = "cell_annot_1")

rownames(DATA_cell@assays$celltypeprops@data)
```

    [1] "B cell"      "Endothelial" "Epithelial"  "Fibroblast"  "Granulocyte"
    [6] "ILC"         "Myeloid"     "NK cells"    "T cell"     

## Plot cell type contribution

``` r
cell_type <- c("T cell", "B cell", "Myeloid", "NK cells", "Epithelial", "Fibroblast", "Endothelial")
clus_lvl <- c("6", "9", "7", "5","8","3","4","0","1","2","10")

DefaultAssay(DATA_cell) <- "celltypeprops"
sum_clus <- FetchData(DATA_cell, cell_type) %>% 
  #na.omit() %>%
  bind_cols(., select(DATA_cell,Clusters)) %>%
  mutate(Clusters = factor(.$Clusters, levels = clus_lvl)) %>%
  summarise(across(all_of(cell_type), ~sum(.)),.by="Clusters") %>%
  mutate(., sum_c = rowSums(across(all_of(cell_type)))) %>%
  #mutate(b=.$`B cell`/.$sum_c)
  mutate(., across(all_of(cell_type), ~(./sum_c))) |> 
  tidyr::pivot_longer(cols = all_of(cell_type), 
                      names_to = "Cell", 
                      values_to = "Weight") %>%
  mutate(Cell = factor(.$Cell, levels = rev(cell_type)))

(C <- ggplot(sum_clus, aes(x=Clusters, y=Cell, size=Weight, color=Weight)) +
  geom_point() +
  labs( title="Cell type contribution", x="Clusters", y = "Cell type", 
       color = "", size = "Scaled weight") +
  scale_color_viridis_c(direction = -1, option = "mako") +
  scale_size("% expressd") +
  theme_bw() +
  
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, margin = margin(t=-1)),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust = .1, margin = margin(b = 1)), 
        legend.text = element_text(size = 6),
        legend.margin=margin(-5,4,2,2), # moves the legend box
        legend.box.margin=margin(-5,15,-3,-30), # moves the legend
        legend.spacing.x = unit(2, 'pt'),
        plot.margin = unit(c(5,1,1,0),units = "pt"),
        legend.title = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9, margin = margin(b = -2)),
        panel.grid = element_blank()) +
    guides(col = guide_colorbar(barwidth = 4, barheight = .5),
           size = guide_legend( keywidth = .5),)
)
```

<img src="../Figures/02/02c_cell-type-contribution-plot.png"
data-fig-align="center" />

## Combine all plots

``` r
#############################
# COMBINE ALL FIGURE PANELS #
#############################
# dev.new(width=7.4, height=2.4, noRStudioGD = TRUE)
B_C <- plot_grid(B, C, rel_widths = c(1,.6),  labels = c('B','C'))

# dev.new(width=7.4, height=7.4, noRStudioGD = TRUE)
Figure2 <- plot_grid( A, B_C, ncol=1, rel_heights = c(1,.85), hjust = -.2, labels = c('A')) 
ggsave("./Figures/Figure-2.pdf", Figure2, width = 7.4, height = 7.4, bg = "white", dpi = 1000)

Figure2
```

<img src="../Figures/02/Figure-2.png" data-fig-align="center" />
