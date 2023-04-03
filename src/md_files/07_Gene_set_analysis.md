Gene Set Analysis (GSEA)
================
4/3/23

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
library(ggridges)
library(scales)
library(niceRplots)
library(fgsea)
library(openxlsx)

source("../bin/plotting_functions.R")
source("../../broliden_5325/code/enrichment_function.R")

#########
# PATHS #
#########
input_dir <- "../results/06_DGE_condition_st_data/"
result_dir <- "../results/07_GSEA_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

GO_path <- "../../broliden_5325/resources/KEGG_GO_database/c5.bp.v6.2.symbols.gmt.txt"
KEGG_path <- "../../broliden_5325/resources/KEGG_GO_database/c2.cp.kegg.v6.2.symbols.gmt.txt"

#############
# LOAD DATA #
#############
# DEGs_table <- read_csv(paste0(input_dir,"subset_100/DGEs_condition_wilcox_epi_100.csv")) %>% filter(groups == "DMPA")
DEGs_table <- read_csv(paste0(input_dir,"DGEs_condition_wilcox.csv"))
DATA <- readRDS(paste0("../results/04_deconvolution_st_data/","seuratObj_deconvolution_scdc.RDS"))

GO_database <- fgsea::gmtPathways(GO_path)
KEGG_database <- fgsea::gmtPathways(KEGG_path)
```

``` r
ord1 <- c("Sup_1","Sup_2","Basal_2","Basal_1","0","1","2","3","4","8","10")
ord2 <- c("6","9","7","5","0","1","2","3","4","8","10")

DEGs_table <- DEGs_table %>%
  mutate(subgroup = factor(.$subgroup, levels = ord2)) %>%
  mutate(layers = factor(.$layers, levels = ord1)) %>%
  dplyr::rename(Clusters="layers") 
```

``` r
####################
# CREATE GENE RANK #
####################
# log.pct.diff , p_val_adj
DEGs_list <- DEGs_table %>% 
  group_split(., subgroup) %>% 
  set_names(., sort(unique(DEGs_table$Clusters)))

get_gene_rank.fun <- function(DEGs_list, stat = "avg_log2FC"){
  stat <- sym(stat)
  gene_rank_df <- DEGs_list %>%
    map(., ~.x %>%
      select(gene, avg_log2FC) %>%
        mutate(.x, "logFC*p-val" = sign(avg_log2FC) * -log10(p_val))  %>%
        mutate(.x, stat_ord = !!(stat)) %>%
        mutate(.x, rank = rank(stat_ord))
        #mutate(.x, statsAdj = sign(stat_ord) * (abs(stat_ord) ^ 1) ) %>%
        #mutate(.x, statsAdj = statsAdj / max(abs(statsAdj))) %>%
  ) %>% map(., ~arrange(.x, stat_ord))
}
gene_rank_df <- get_gene_rank.fun(DEGs_list, "avg_log2FC")
#gene_rank_df <- get_gene_rank.fun(DEGs_list, "logFC*p-val")

gene_rank_list <- gene_rank_df %>%
  map(., ~pull(.x, stat_ord)) %>% 
  map2(., DEGs_list, ~set_names(.x, pull(.y, "gene")))
```

Here is a documentation for GSEA, pointing to Enrichment Scores in
particular:
https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#*Enrichment_Score*(ES)

And here is one of the clearest and most thorough explanations I’ve
found: https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/

Answering your specific questions:

Larger magnitudes are “better” The sign on the score simply indicates
when end of your ranked gene list is enriched. You provide the rank list
of genes, so the biiological interpretation is up to you. If, for
example, you provide a gene list ranked by a combination of fold change
and p-value (e.g., sign(FC) \* log10(pvalue)), then the positive scores
are associated with upregulated genes and negative scores are associated
with downregulated genes. Caution: some tools reverse this, so manually
check a few to see which convention they are using. ES values range from
-1 to 1. Normalized ES values will go a bit beyond these bounds. If you
are seeing ES values \> 1, then I would suspect something is wrong.

Additional comment: You probably want to look at NES (Normalized ES) if
that is provided by the tool you are using. The first two points above
apply just the same.

``` r
# http://www.baderlab.org/ChangjiangXu/GSEA
# Create a gene rank based on the parameter chosen in the previous chunck
#####################
# RUN GESA ANALYSIS #
#####################
fgseaRes_GO <- gene_rank_list %>%
  map(., ~fgsea(pathways=GO_database, stats=.x)) %>%
  map(., ~as_tibble(.x)) %>%
  map(., ~arrange(.x, pval))

#fgseaRes_KEGG <- NULL
#fgseaRes_GO <- NULL
fgseaRes_KEGG <- gene_rank_list %>%
  map(., ~fgsea(pathways=KEGG_database, stats=.x)) %>%
  map(., ~as_tibble(.x)) %>%
  map(., ~arrange(.x, pval))

fgseaRes <- append(fgseaRes_GO, fgseaRes_KEGG)

# Show in a nice table:
# fgseaResTidy <- fgseaRes$`8` %>% 
#   top_n(15, abs(NES)) %>% 
#   dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#   arrange(padj) #%>% 
#   #DT::datatable()
```

### save GESA object

``` r
saveRDS(fgseaRes, paste0(result_dir,"fgseaRes_tables.RDS"))
```

``` r
fgseaRes <- readRDS(paste0(result_dir,"fgseaRes_tables.RDS"))
```

``` r
# gene_rank_df <- gene_rank_df[[6]]
# fgseaRes <- fgseaRes[[6]]
# top_n <- 20
db <- c(rep("GO", length(fgseaRes)/2),rep("KEGG", length(fgseaRes)/2))
GSEA_table.fun <- function(fgseaRes, gene_rank_df, top_n, Cluster, n = 40){
  f <- fgseaRes  %>% # fgseaRes[["Basal_2"]]
    mutate(Direction = ifelse(ES > 0, "UP", "DOWN")) %>%
    group_by(Direction) %>%
    top_n(-pval, n = top_n) %>%
    ungroup() %>%
    #filter(padj < 0.05) %>%
    mutate(overlap = map_dbl(leadingEdge, ~length(.x)), .after = "size") %>%
    unnest(leadingEdge) %>%
    left_join(., select(gene_rank_df, gene, avg_log2FC), by=c("leadingEdge"="gene")) %>% # gene_rank_df$Basal_2
    mutate(., dir = ifelse(.$avg_log2FC < 0, "down", "up")) %>%
    group_by(pathway) %>% 
    mutate(min_logFC = min(avg_log2FC)) %>%
    mutate(max_logFC = max(avg_log2FC)) %>% 
    add_count(., dir) %>%
    ungroup() %>%
    select(pathway, overlap, dir, n, min_logFC, max_logFC) %>%
    unique() %>% #filter(pathway == "GO_PEPTIDE_CROSS_LINKING")
    pivot_wider(.,  names_from = dir, values_from = n)
  
  f_ <- fgseaRes %>% # fgseaRes[["Basal_2"]]
    left_join(., f, by="pathway") %>%
    mutate(dir = ifelse(is.na(.$down) & is.na(.$up), NA,
                        ifelse(is.na(.$up), "down",
                          ifelse(is.na(.$down) | .$down < .$up, "up", "down")))) %>%
    mutate(Direction = ifelse(ES > 0, "UP", "DOWN")) %>%
    select(-leadingEdge, leadingEdge)
  return(f_)
}

ridge_plots <- tibble(Clusters=factor(names(fgseaRes), levels = levels(DEGs_table$Clusters)),
                      db = db,
                      fgseaRes = fgseaRes) %>% arrange(Clusters)

epi_clus <- "Sup_0|Sup_1|Sup_2|Intrmed|Basal_2|Basal_1"
tab <- ridge_plots %>%
  mutate(tab = pmap(., ~GSEA_table.fun(fgseaRes=..3, gene_rank_df[[..1]], top_n=25, Cluster=..1))) %>%
  mutate(tab = map(tab, ~split(.x, .x$dir))) %>%
  unnest(c(tab)) %>%
  mutate(dir = names(tab)) %>%
  mutate(morf = ifelse(grepl(epi_clus, .$Clusters), "epi", "SubMuc")) %>%
  mutate(tab = set_names(.$tab, paste0(.data[["Clusters"]],"_",.data[["dir"]]))) %>%
  select(db, tab, morf) %>%
  group_by(morf, db) %>%
  summarise(tab = list(tab), .groups = "drop")

pmap(tab, ~write.xlsx(..3, keepNA=TRUE, na.string="NA", overwrite=TRUE,
                             file=paste0(result_dir,"GESA_",..1,"_",..2,".xlsx"))) 
```

### Split violin plots of genes in GO terms of intrest

``` r
# dev.new(width=6.6929133858, height=5*3, noRStudioGD = TRUE)
######################
# GENES IN GESA TERM #
#####################
terms <- c("GO_KERATINOCYTE_DIFFERENTIATION", "GO_PEPTIDE_CROSS_LINKING", "GO_SKIN_DEVELOPMENT", "GO_B_CELL_MEDIATED_IMMUNITY")
terms.fun <- function(tab, term){
  term_genes <- tab %>%
    unnest(c(tab)) %>%
    mutate(clus = str_replace(names(.$tab), "_down|_up", "") ) %>%
    unnest(c(tab)) %>%
    filter(pathway == term & morf == "epi") %>%
    mutate(clus_all = paste(unique(.$clus), collapse = "$|^")) %>%
    arrange(desc(overlap)) %>%
    #filter(overlap == max(overlap)) %>%
    summarise(genes = unique(unlist(.$leadingEdge)), clus = .$clus[1], clus_all = .$clus_all[1]) %>%
    mutate( genes = factor(.$genes, levels = gene_rank_df[[.$clus[1]]]$gene)) %>%
    nest(data = c(genes)) %>%
    
  return(term_genes)
}

term_genes <- terms %>%
  map(., ~terms.fun(tab, .x)) %>% list_rbind()

################
# VIOLIN PLOTS #
################
library(patchwork)
violin_split.fun <- function(obj, gene_list, split.by, group.by, filt_pattern){
  obj <- filter(obj, grepl(filt_pattern, obj$layers))
  plots <- VlnPlot(obj, features = gene_list, adjust = 1/2, cols =  c("#88CCEE", "#CC6677"),
                   split.by = split.by, group.by = group.by, 
                   pt.size = 0, combine = FALSE, split.plot = TRUE) 
  plots <- map(plots, ~.x + geom_point(alpha = .3, size=.1, show.legend = F,
                                       position = position_jitterdodge(jitter.width=.5)) +
                 #scale_fill_manual(values = c("#88CCEE", "#CC6677")) +
                 #scale_color_manual(values = c("#88CCEE", "#CC6677")) +
                 theme(axis.title.x = element_blank()))
  p <- wrap_plots(plots, ncol = 10) + plot_layout(guides = "collect")
  h <- round(length(plots)/10)*4
return(tibble(plot=list(p), height = list(h), p = list(plots)))
}

plots <- term_genes %>%
  mutate(data = map(data, ~arrange(.x, genes))) %>%
  mutate(data = map(data, ~unlist(.x))) %>%
  pmap(., ~violin_split.fun(DATA, ..3, "groups", "layers", ..2)) %>%
  list_rbind() %>%
  mutate(term = terms)

pmap(plots, ~ggsave(paste0("./Figures/07/","DEGs_",..4,".png"), plot=..1, width = 40, height = ..2))
```

    [[1]]
    [1] "./Figures/07/DEGs_GO_KERATINOCYTE_DIFFERENTIATION.png"

    [[2]]
    [1] "./Figures/07/DEGs_GO_PEPTIDE_CROSS_LINKING.png"

    [[3]]
    [1] "./Figures/07/DEGs_GO_SKIN_DEVELOPMENT.png"

    [[4]]
    [1] "./Figures/07/DEGs_GO_B_CELL_MEDIATED_IMMUNITY.png"

``` r
wrap_plots(plots[["p"]][[2]], ncol = 3) + plot_layout(guides = "collect")
```

<img src="./Figures/07/07a_violin_plots.png" data-fig-align="center" />

``` r
# dev.new(width=12, height=3, noRStudioGD = TRUE)
# width = aspect * height 4.857143
# height = width/aspect 10/4.857143
##################
# GESA RANK PLOT #
##################
db <- c(rep("GO", length(fgseaRes)/2),rep("KEGG", length(fgseaRes)/2))
plot_GseaTable.fun <- function(db, fgseaRes, gene_rank_list, n = 10, padj = 0.05){
    if(db == "GO"){database <- GO_database}else{database <- KEGG_database}
    
    #fgseaRes <- mutate(fgseaRes, pathway = formating_str.fun(fgseaRes$pathway))
  
    topPathwaysUp <- fgseaRes[fgseaRes$ES > 0 & fgseaRes$padj < padj,]
    topPathwaysDown <- fgseaRes[fgseaRes$ES < 0 & fgseaRes$padj < padj,]
    
    topPathwaysUp <- unlist(topPathwaysUp[head(order(topPathwaysUp$pval), n=n), "pathway"])
    topPathwaysDown <- unlist(topPathwaysDown[head(order(topPathwaysDown$pval), n=n), "pathway"])
  
    p <- plotGseaTable(database[c(topPathwaysUp, rev(topPathwaysDown))], 
                  gene_rank_list, fgseaRes, gseaParam = 0.5, render=F)
    return(tibble(rank_plots=list(p), n_pathways = length(topPathwaysUp) +length(topPathwaysDown)))  # tibble(plots=list(plot), df=list(f))
    }

# fgseaRes <- fgseaRes %>%
#   map(., ~mutate(.x, pathway = formating_str.fun(.x$pathway)))

ridge_plots <- tibble(Clusters=names(fgseaRes),db = db,fgseaRes = fgseaRes) 
plots_df <- pmap_dfr(ridge_plots,
              ~plot_GseaTable.fun(db=..2, fgseaRes=..3, n=10,
                                  gene_rank_list[[..1]]) ) %>%
  bind_cols(ridge_plots, .) %>%
  mutate(rank_plots = setNames(.[["rank_plots"]], .$Clusters)) %>%
  filter(n_pathways > 0)

if( isFALSE(dir.exists("./Figures/07/GO/")) ) { dir.create("./Figures/07/GO/",recursive = TRUE) }
if( isFALSE(dir.exists("./Figures/07/KEGG/")) ) { dir.create("./Figures/07/KEGG/",recursive = TRUE) }

pmap(plots_df, ~ggsave(paste0("./Figures/07/",..2,"/GSEA_rank_plot_", ..1,"_",..2, ".png"), #"GO_filt/",
                          plot=..4, width = 13, height = (..5*.3)+2 ))
```

    [[1]]
    [1] "./Figures/07/GO/GSEA_rank_plot_0_GO.png"

    [[2]]
    [1] "./Figures/07/GO/GSEA_rank_plot_1_GO.png"

    [[3]]
    [1] "./Figures/07/GO/GSEA_rank_plot_2_GO.png"

    [[4]]
    [1] "./Figures/07/GO/GSEA_rank_plot_3_GO.png"

    [[5]]
    [1] "./Figures/07/GO/GSEA_rank_plot_4_GO.png"

    [[6]]
    [1] "./Figures/07/GO/GSEA_rank_plot_8_GO.png"

    [[7]]
    [1] "./Figures/07/GO/GSEA_rank_plot_Basal_1_GO.png"

    [[8]]
    [1] "./Figures/07/GO/GSEA_rank_plot_Basal_2_GO.png"

    [[9]]
    [1] "./Figures/07/GO/GSEA_rank_plot_Sup_1_GO.png"

    [[10]]
    [1] "./Figures/07/GO/GSEA_rank_plot_Sup_2_GO.png"

    [[11]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_0_KEGG.png"

    [[12]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_1_KEGG.png"

    [[13]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_2_KEGG.png"

    [[14]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_3_KEGG.png"

    [[15]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_4_KEGG.png"

    [[16]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_8_KEGG.png"

    [[17]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_Sup_1_KEGG.png"

    [[18]]
    [1] "./Figures/07/KEGG/GSEA_rank_plot_Sup_2_KEGG.png"

``` r
grid::grid.draw(plots_df$rank_plots$Basal_2)
```

<img src="./Figures/07/07b_GSEA_rank_plot.png"
data-fig-align="center" />
