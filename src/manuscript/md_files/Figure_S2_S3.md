Figure S2 & S3
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
library(patchwork)
library(ggrepel)
library(ggpubr)
library(scales)
library(niceRplots)
#library(fgsea)
library(openxlsx)
library(jsonlite)
library(ggraph)
library(tidygraph)
library(igraph)
library(oaqc)
library(ggnewscale)

# BiocManager::install("DOSE")
# BiocManager::install("GSEABase")
# BiocManager::install("GO.db")
# BiocManager::install("org.Hs.eg.db")
library(DOSE)
library(GSEABase)
library(GO.db)
library(org.Hs.eg.db)
select <- dplyr::select

source("../../bin/plotting_functions.R")
#source("../../broliden_5325/code/enrichment_function.R")

#########
# PATHS #
#########
input_dir <- "../results/06_DGE_condition_st_data/"
result_dir <- "../../results/07_GSEA_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

KEGG_path <- "../../resources/c2.cp.kegg_legacy.v2023.2.Hs.json"
GO_path <- "../../resources/c5.all.v2023.2.Hs.json"
epi_clus <- "^5$|^6$|^7|^9"
submuc_remove <- "1|2"

#############
# LOAD DATA #
#############
DEGs_table <- read_csv(paste0("../../results/06_DGE_condition_st_data/","DGEs_condition_wilcox.csv"))
DATA <- readRDS(paste0("../../results/04_deconvolution_st_data/","seuratObj_deconvolution_scdc.RDS"))
fgseaRes <- readRDS(paste0(result_dir,"fgseaRes_tables_updated_DB.RDS"))

read_json.fun <- function(db){
  db <- db %>%
  map(., ~as_tibble(.x[4])) %>%
  bind_rows(., .id = "pathway") %>%
  mutate(db = map(db, ~pluck(.x, "geneSymbols")))
  return(db)}
GO_database <- read_json.fun(fromJSON(GO_path))
KEGG_database <- read_json.fun(fromJSON(KEGG_path))

# set seed for all operations
c <- addTaskCallback(function(...) {set.seed(123);TRUE})
# removeTaskCallback(c)
# sample(1:500, 3)
```

``` r
####################
# CREATE GENE RANK #
####################
# log.pct.diff , p_val_adj
DEGs_list <- DEGs_table %>% 
  group_split(., subgroup) %>% 
  set_names(., sort(unique(DEGs_table$layers)))

get_gene_rank.fun <- function(DEGs_list, stat = "avg_log2FC"){
  stat <- sym(stat)
  gene_rank_df <- DEGs_list %>%
    map(., ~.x %>%
      dplyr::select(gene, avg_log2FC) %>%
        mutate(.x, "logFC*p-val" = sign(avg_log2FC) * -log10(p_val))  %>%
        mutate(.x, stat_ord = !!(stat)) %>%
        mutate(.x, spatial_dist = rank(stat_ord)) %>%
        mutate(.x, rank = rank(stat_ord))
        #mutate(.x, statsAdj = sign(stat_ord) * (abs(stat_ord) ^ 1) ) %>%
        #mutate(.x, statsAdj = statsAdj / max(abs(statsAdj))) %>%
  ) %>% map(., ~arrange(.x, desc(stat_ord)))
}
gene_rank_df <- get_gene_rank.fun(DEGs_list, "avg_log2FC")
#gene_rank_df <- get_gene_rank.fun(DEGs_list, "logFC*p-val")

gene_rank_list <- enframe(gene_rank_df) %>%
  mutate(gene_list = pmap(., ~set_names(pull(..2, "stat_ord"), pull(..2, "gene") )) ) %>%
  mutate(gene_list = setNames(.[["gene_list"]], .$name)) %>%
  .$gene_list

# seems faulty, so now I don't trust the GSEA results
# gene_rank_list <- gene_rank_df %>%
#   map(., ~pull(.x, stat_ord)) %>% 
#   map2(., DEGs_list, ~set_names(.x, pull(.y, "gene")))
```

``` r
# Calculate leading edge analysis from a fgesa results table
gseaScores <- getFromNamespace("gseaScores", "DOSE")
# tmp_res <- fgseaRes$Basal_1
# geneSets <- GO_database
# geneList <- gene_rank_list$Basal

# NB! this function assumes that the genes rank is ordered by increasing values (pos logFC to neg logFC)
# this affects the order of the core genes, but not the 
leading_edge <- function(res, geneSets, geneList){
  
    # Observed info (gseaScore):
    observed_info <- lapply(geneSets[res$pathway], function(gs)
        gseaScores(geneSet=gs,
                   geneList=geneList,
                   exponent=1)
        )
    
    observed_info <<- observed_info
    # Leading edge analysis:
    core_enrichment <- lapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            leading_gene <- runningES$gene[1:i]
        } else {
            i <- which.min(runningES$runningScore)
            leading_gene <- runningES$gene[-c(1:i)]
        }
        return(leading_gene)
    })

    rank <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            rr <- which.max(runningES$runningScore)
        } else {
            i <- which.min(runningES$runningScore)
            rr <- nrow(runningES) - i + 1
        }
        return(rr)
    })

    tags <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    ll <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) { # if true assumes up regulation
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else { # assumes down regulation
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    
    N <- length(geneList)
    signal <- tags * (1-ll) * (N / (N - res$size))

    tags <- paste0(round(tags * 100), "%")
    ll <- paste0(round(ll * 100), "%")
    signal <- paste0(round(signal * 100), "%")
    LE_analysis <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)
    
    res$rank_at_max <- rank
    res$leading_edge <- LE_analysis
    res$core_enrichment <- sapply(core_enrichment, paste0, collapse=',')
    return(res)

}
db_ids <- list("GO" = set_names(GO_database$exactSource, GO_database$pathway),
               "KEGG" = set_names(KEGG_database$exactSource, KEGG_database$pathway))

db <- c(rep("GO", length(fgseaRes)/2),rep("KEGG", length(fgseaRes)/2))
fgseaRes_df <-  tibble(Clusters=names(fgseaRes),
                      db = db, fgseaRes = fgseaRes) %>%
  mutate(fgseaRes = pmap(., ~mutate(..3, id = db_ids[[..2]][.data[["pathway"]] ], .before="pathway" ) )) %>% # 
  mutate(fgseaRes = pmap(., ~filter(..3, padj < 0.05) )) %>%
  filter(unlist(map(.$fgseaRes, ~nrow(.x) > 0)))

fgseaRes_LE <- fgseaRes_df %>% 
  mutate(fgseaRes = pmap(., ~leading_edge(..3, get(paste0(..2,"_database"))$db, gene_rank_list[[..1]]) )) %>%
  mutate(fgseaRes = set_names(.$fgseaRes, paste0(.data[["Clusters"]]) )) %>%
  {. ->> fgseaRes_filt} %>%
  group_by(db) %>%
  summarise(fgseaRes = list(fgseaRes), .groups = "drop")

# save object:
saveRDS(fgseaRes_filt, paste0(result_dir,"fgseaRes_filt_tables.RDS")) 
# fgseaRes_LE <- readRDS(paste0(result_dir,"fgseaRes_filt_tables.RDS"))
```

``` r
path_slim <- "/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/resources/goslim_agr.obo"
slim <- getOBOCollection(path_slim)

mappedIds <- function(df, collection, OFFSPRING){
    map <- as.list(GO.db::GOBPOFFSPRING[rownames(df)])
    mapped <- lapply(map, intersect, ids(collection))
    #df[["go_ids"]] <- vapply(unname(mapped), paste, collapse = ",", character(1L))
    df[["go_ids"]] <- mapped
    df[["go_terms"]] <- map(unname(mapped), ~paste(go[.x]), collapse = ";") # paste(go[unlist(mapped)], collapse = ";")
    df
  }

goterms <- tibble("Term"=Term(GOTERM), "id"=names(Term(GOTERM))) %>%
  mutate(t = paste0("GOBP_", toupper(.$Term)) ) %>%
  mutate(t = gsub(x = .$t, " |-|, |/","_" ))
go <- set_names(goterms$Term, goterms$id)

slim_df <- fgseaRes_df %>%
  filter(db == "GO") %>%
  mutate(morf = ifelse(grepl("\\d", .$Clusters), "SubMuc", "epi")) %>%
  filter(!(grepl(submuc_remove, .$Clusters))) %>%
  # the old go ids from the GOTERM db that was missing some ids:
    mutate(fgseaRes = pmap(., ~left_join(..3, select(goterms, t, GOTERM_ids="id"), by=c("pathway"="t")) )) %>%
    mutate(GOTERM_ids = pmap(., ~..3[["GOTERM_ids"]])) %>%
  mutate(go_ids = pmap(., ~..3[["id"]])) %>%
  dplyr::select(morf, go_ids) %>%
  
  # GO slim data:
  summarise(go_ids = list(unique(unlist(.$go_ids[cur_group_rows()]))), .by="morf") %>%
  mutate(myCollection = pmap(., ~GOCollection(na.omit(..2), ontology="BP") )) %>%
  mutate(slimdf = pmap(., ~goSlim(..3, slim, "BP") )) %>%
  mutate(slim = pmap(., ~as_tibble(mappedIds(..4, ..3, GOMFOFFSPRING), rownames = "SlimID") )) %>%
  
  mutate(slim = setNames(.[["slim"]], .$morf))

list("Epithelium"=slim_df$slim[[1]], "Submucosa"=slim_df$slim[[2]]) %>%
  {. ->> slim_lst} #%>%
  #write.xlsx(., paste0(result_dir,"GOslim_summary_New.xlsx"))
```

``` r
############################
# SAVE SUPPLEMENTARY TABLE #
############################
append(slim_lst, fgseaRes_LE$fgseaRes[[1]]) %>%
write.xlsx(., keepNA=TRUE, na.string="NA", overwrite=TRUE,
                             file=paste0(result_dir,"Supplemental Table 5.xlsx"))
```

``` r
pal <- c(RColorBrewer::brewer.pal(10,"Paired"),
         RColorBrewer::brewer.pal(8,"Set2"),
         RColorBrewer::brewer.pal(8,"Dark2") )
         #RColorBrewer::brewer.pal(9,"Pastel1"),
         #RColorBrewer::brewer.pal(8,"Pastel2")

df <- map(slim_df$slim, ~ .x %>%
  filter(!(Percent == 0)) %>%
  arrange(desc(Percent)) %>%
  mutate(Term = factor(.$Term, levels = .$Term)) %>%
  mutate(ystart = c(0, cumsum(Percent)[1:nrow(.)-1])) %>%
  mutate(yend = lead(ystart, default = 100)) 
)

ggplot(df$epi) + 
  # ajust the width of the circle by changing max and min x values
  # to get a rageular pice chart set xmin=0
  geom_rect(aes(fill=Term, ymax=yend, ymin=ystart, xmax=4, xmin=3)) +
  # geom_rect(aes(fill=SlimID, ymax=yend, ymin=ystart, xmax=3, xmin=0)) +
  xlim(c(0, 4)) + scale_fill_manual(values = pal) +
  theme(aspect.ratio=1) + theme_void() +
  coord_polar(theta="y")  
```

<img src="../Figures/S2&amp;3/piechart-GOslim.png"
data-fig-align="center" />

``` r
sig_epi_GO <- fgseaRes_filt$fgseaRes[1:4] %>% bind_rows(., .id = "clus") %>% filter(!(grepl("\\d", .$clus))) %>% filter(padj<=0.05) #%>% .$pathway %>% unique()
sig_sub_GO <- fgseaRes_filt$fgseaRes[5:8] %>% bind_rows(., .id = "clus") %>% filter(grepl("\\d", .$clus)) %>% filter(padj<=0.05) #%>% .$pathway %>% unique()
#sig_sub_GO <- fgseaRes[5:8] %>% bind_rows(., .id = "clus") %>% filter(grepl("\\d", .$clus)) %>% filter(padj<=0.05)

# because i filtered away the terms (n=38) with no match in the GOTERM db, the result is different 
# when not removing these, much more messy plot, here are the missing ids:
go_filt <- c("GO:0034341", "GO:0140888", "GO:0048871", "GO:1903829", "GO:0034341", "GO:0051603", "GO:0071346", "GO:0043603", "GO:0048871", "GO:0007249", "GO:0010608", "GO:0140888", "GO:0071346", "GO:0034341", "GO:0140888", "GO:0051603", "GO:0009063", "GO:0007249", "GO:0048871", "GO:0034341", "GO:0071346", "GO:0048871", "GO:0071824", "GO:0006974", "GO:0034341", "GO:0071346", "GO:0043603", "GO:0034248", "GO:0051603", "GO:0071824", "GO:0010608", "GO:0006458", "GO:0006278", "GO:1903050", "GO:0071826", "GO:0040029", "GO:0034249", "GO:0043603", "GO:0010608", "GO:0034248", "GO:0051603", "GO:0048871", "GO:0034250", "GO:0034766", "GO:0071826", "GO:0034341", "GO:1903050", "GO:1903829", "GO:0071824", "GO:0051897", "GO:1903828", "GO:1903051", "GO:0071346", "GO:0015986", "GO:0043271", "GO:1902895", "GO:0006278", "GO:0015793", "GO:0032781")
```

``` r
# Modified code by YuLab from enrichplot package:
# https://github.com/YuLab-SMU/enrichplot/blob/devel/R/emapplot_utilities.R

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    # gene overlap between nodes are used as weights:
    length(intersect(x, y))/length(unique(c(x,y))) }

##' Get the similarity matrix
##'
##' @param y A data.frame of enrichment result
##' @param geneSets A list, the names of geneSets are term ids,
##' and every object is a vertor of genes.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @param semData GOSemSimDATA object
##' @noRd
d <- GOSemSim::godata('org.Hs.eg.db', ont="BP")

# gene overlap between nodes are used as weights:
get_similarity_matrix <- function(y, geneSets, method="JC", semData = NULL) {
    id <- y$pathway
    geneSets <- set_names(y$leadingEdge, id)
    
    n <- nrow(y)
    y_id <- unlist(strsplit(y$pathway, "_"))[1]
    ## Choose the method to calculate the similarity
    if (method == "JC") {
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$pathway
        for (i in seq_len(n-1)) {
            for (j in (i+1):n) {
                w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        return(w)
    }

    if (y_id == "GO" | y_id == "GOBP") {
        if(is.null(semData)) {
            stop("The semData parameter is missing,
                and it can be obtained through godata function in GOSemSim package.")
        }
        w <- GOSemSim::mgoSim(id, id, semData=semData, measure="Wang",
                              combine=NULL)
    }

    if (y_id == "DOID") w <- DOSE::doSim(id, id, measure="Wang")
    rownames(y) <- y$pathway
    rownames(w) <- colnames(w) <- y[colnames(w), "pathways"]
    return(w)
}
build_emap_graph <- function(sim_m, fgesa_df, method = "JC", cex_line = 1, min_edge = 0.2){
  # Use similarity as the weight(length) of an edge
  wd <- pivot_longer(cols = -1, as_tibble(sim_m, rownames = "from"), names_to = "to", values_to = "weight")
      wd <- wd[wd[,1] != wd[,2],]
      # remove NA
      wd <- wd[!is.na(wd[,3]),]
      
      # Create network tidy:
      nodes <- fgesa_df %>% dplyr::select(pathway)
      edges <- wd %>% mutate(width = sqrt(.$weight * 5)* cex_line)
      
      graph <- list()
      g <- tbl_graph(nodes, edges, directed = FALSE)
      # filter edges
      g <- g %>% activate(edges) %>% filter(weight > min_edge) %>% activate(nodes)
      graph$tidy <- g 
      
      # Create network igraph:
      # g <- graph.data.frame(wd[, -3], directed=FALSE)
      # E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
      # # Use similarity as the weight(length) of an edge
      # E(g)$weight <- wd[, 3]
      # # filter edges
      # g <- delete.edges(g, E(g)[wd[, 3] > min_edge])
      # graph$igraph <- g
  
  return(graph)
}
```

``` r
dge <- DEGs_table %>% 
  filter(grepl("[^\\d]", .$layers)) %>% 
  filter(p_val<=0.01) %>%
  #filter(pct.1 > .4 & p_val < 0.01)  
  split(., ~layers) %>%  map(., ~.x$gene)
# temp_ <- fgseaRes$Basal_1 %>% filter(padj<=0.05)
# m <- get_similarity_matrix(temp_, geneSets=GO_database, semData=d) 
# graph <- build_emap_graph(sim_m=m, fgesa_df=temp_, cex_line=1, min_edge=0)

ord1 <- c("Superficial","Upper IM","Lower IM","Basal","8","3","4","0","1","2","10")
ord2 <- c("Sup_1", "Sup_2", "Basal_2", "Basal_2","8","3","0","4","1","2","10")
l <- set_names(ord1, ord2)

GO_comb_epi <- sig_epi_GO %>% 
  #mutate(., clus = factor(l[as.character(.$clus)], levels = ord1)) %>%
  add_count(pathway) %>% 
  dplyr::rename(GOid="id") %>%
  left_join( dplyr::select(goterms, Term, t), by=c("pathway"="t")) %>%
  #mutate(Term = str_replace_all(tolower(.$pathway), "_", " ")) %>%
  #mutate(Term = str_replace(.$Term, "^gobp ", "")) %>%
  # filter leading edge by parameters set above:
  #mutate("p<0.01" = pmap(., ~intersect(..10, dge[[..1]]))) %>%
  # combine same terms from different clusters:
  mutate("leadingEdge_comb" = list(unique(unlist(.$leadingEdge[cur_group_rows()]))), .by="pathway") %>%
  mutate(pathway = paste0(.$clus, ":", .$GOid)) %>%
  mutate("Sig. in" = paste0(.$clus[cur_group_rows()], collapse = "|"), .by="GOid") %>%
  #filter((is.na(.$id))) %>%
  filter(!(.$GOid %in% go_filt)) %>%
  # removes the idential terms:
  slice_max(n=1, order_by = tibble(abs(NES), size), by="Term") %>%
  ungroup()
  
#### get similarity matrix and create graph ####
# weight is the percentage overlap between the two terms
m <- get_similarity_matrix(GO_comb_epi, geneSets=GO_database, semData=d) 
set.seed(123); graph <- build_emap_graph(sim_m=m, fgesa_df=GO_comb_epi, cex_line=1, min_edge=0.1)
# G <- as_tbl_graph(graph$igraph)

#### find number of edges between nodes: ####
m %>% as_tibble(., rownames="from")
```

    # A tibble: 324 × 325
       from     Superficial:GO:00022…¹ Superficial:GO:00197…² Superficial:GO:00024…³
       <chr>                     <dbl>                  <dbl>                  <dbl>
     1 Superfi…                     NA                  0.276                  0.537
     2 Superfi…                     NA                 NA                      0.42 
     3 Superfi…                     NA                 NA                     NA    
     4 Superfi…                     NA                 NA                     NA    
     5 Basal:G…                     NA                 NA                     NA    
     6 Superfi…                     NA                 NA                     NA    
     7 Superfi…                     NA                 NA                     NA    
     8 Basal:G…                     NA                 NA                     NA    
     9 Superfi…                     NA                 NA                     NA    
    10 Superfi…                     NA                 NA                     NA    
    # ℹ 314 more rows
    # ℹ abbreviated names: ¹​`Superficial:GO:0002250`, ²​`Superficial:GO:0019724`,
    #   ³​`Superficial:GO:0002449`
    # ℹ 321 more variables: `Superficial:GO:0002460` <dbl>,
    #   `Basal:GO:0044419` <dbl>, `Superficial:GO:0002443` <dbl>,
    #   `Superficial:GO:0006958` <dbl>, `Basal:GO:0006956` <dbl>,
    #   `Superficial:GO:0002455` <dbl>, `Superficial:GO:0002252` <dbl>, …

``` r
summary <- graph$tidy %>% activate(edges) %>% 
  as_tibble() %>%
  summarise(l = sum(weight), # n is number of shared genes with other terms, l is the sum of weights
            n_genes = sum(weight > 0), .by = "from" )  %>% 
  arrange(n_genes)
s <- set_names(summary$n_genes, as.character(summary$from))

# graph$tidy %>% activate(edges) %>% as_tibble() %>% filter(weight < 0)
# https://tidygraph.data-imaginist.com/reference/index.html

#### Add info to Graph ####
# https://bookdown.org/markhoff/social_network_analysis/centrality.html
set.seed(123);G <- graph$tidy %>% 
  mutate(gene_module = str_extract(as_tibble(.)$pathway, "[^:]+") ) %>%
  left_join(., dplyr::select(GO_comb_epi, size, pathway, "Sig. in"), by="pathway") %>%
  mutate(pathway = str_extract(as_tibble(.)$pathway, "GO:\\d+") ) %>%
  mutate(degree = centrality_degree(),
         btwn = centrality_betweenness(),
         closness = centrality_closeness(),
         community = as.character(group_louvain())) %>%
  mutate(community = factor( as_tibble(.)$community, 
                             levels = as.character(1:max(as.integer(as_tibble(.)$community))) )) %>%
  left_join(goterms, by=c("pathway"="id")) %>%
  
  #filter(community == "1" | community == "2" | community == "3" )  %>% #  | community == "3" | community == "4" | community == "5"
  {. ->> G_nodes} %>%
  activate(edges) %>% 
  left_join(summary, by="from") %>%
  mutate(term_from = .N()$Term[from]) %>%
  mutate(term_to = .N()$Term[to]) %>%
  mutate(comm = .N()$community[from]) %>%
  {. ->> G_edges} %>%
  activate(nodes)

#### get table version of graph ###
G_nodes_e <- G_nodes %>% as_tibble() # %>% split(., ~community)
Reduce(intersect, as.list(GOBPANCESTOR[na.omit(G_nodes_e$pathway)]))
```

    [1] "GO:0008150" "all"       

``` r
G_edges %>%
  #filter(weight > .6) 
  filter(!(comm == "1" & weight > .2))
```

    # A tbl_graph: 324 nodes and 2107 edges
    #
    # An undirected simple graph with 31 components
    #
    # A tibble: 2,107 × 9
       from    to weight width     l n_genes term_from                term_to  comm 
      <int> <int>  <dbl> <dbl> <dbl>   <int> <chr>                    <chr>    <fct>
    1     1     7  0.132 0.811  8.64      46 adaptive immune response complem… 1    
    2     1     8  0.108 0.736  8.64      46 adaptive immune response complem… 1    
    3     1     9  0.145 0.851  8.64      46 adaptive immune response humoral… 1    
    4     1    22  0.152 0.873  8.64      46 adaptive immune response regulat… 1    
    5     1    25  0.138 0.830  8.64      46 adaptive immune response regulat… 1    
    6     1    26  0.135 0.823  8.64      46 adaptive immune response positiv… 1    
    # ℹ 2,101 more rows
    #
    # A tibble: 324 × 10
      pathway    gene_module  size `Sig. in`   degree  btwn closness community Term 
      <chr>      <chr>       <int> <chr>        <dbl> <dbl>    <dbl> <fct>     <chr>
    1 GO:0002250 Superficial   339 Superficia…     46  785.  0.00138 1         adap…
    2 GO:0019724 Superficial    99 Superficia…     21  570.  0.00122 6         B ce…
    3 GO:0002449 Superficial   215 Superficia…     38  325.  0.00131 1         lymp…
    # ℹ 321 more rows
    # ℹ 1 more variable: t <chr>

``` r
#### filter graph ####
# https://bookdown.org/markhoff/social_network_analysis/centrality.html
t <- c("phagocytosis, engulfment","membrane invagination", "positive regulation of guanylate cyclase activity", "multicellular organismal water homeostasis","humoral immune response", "leukocyte mediated immunity", "adaptive immune response", "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", "positive regulation of B cell activation","retina homeostasis", "leukocyte mediated immunity")
set.seed(123);G_epi <- G %>% 
  #filter(!(grepl(paste(t,collapse="|"), Term))) %>%
  activate(edges) %>%
  #filter(grepl("|^1$|^2$|", comm) & weight < .2) %>%
  #filter(!(comm == "1" & weight > .8)) %>%
  filter(!(comm == "1" & weight < .3 )) %>%
  filter(!(comm == "3" & weight < .3 )) %>%
  filter(!(comm == "2" & weight < .3 )) %>%
  filter(!(comm == "4" & weight < .3 )) %>%
  #filter(weight > .3) %>%
  #filter(grepl("|^1$|", comm) & weight > .6) %>%
  arrange(desc(weight)) %>%
  #filter(n > 10) %>%
  # mutate(width = .["width"]*.1) %>%
  activate(nodes) %>%
  # re-calculate degree after filtering and remove nodes with no connections
  # before re-clustering
  mutate(degree = centrality_degree()) %>%
  filter(degree > 0) %>%
  mutate(com = as.character(group_louvain()))

G_nodes_epi <- G_epi %>% as_tibble() 

g_table_epi <- G_epi %>% as_tibble() %>%
  #mutate(col = pal[.$com]) %>%
  dplyr::select(community, com, pathway, Term, degree) %>% unique() %>% arrange(community) %>%
  mutate(community = factor( .$community, levels = as.character(1:max(as.integer(.$community))) )) %>%
  mutate(count = n_distinct(Term), .by = c("com")) %>%
  # determines the number of rows are included in the table:
  filter(count > 4) %>%
  left_join(., dplyr::select(GO_comb_epi, Term, ES), by="Term") %>%
  slice_max(n=2, order_by = tibble(degree, abs(ES)), by="com", with_ties = F) %>%
  slice_max(n=1, order_by = tibble(abs(ES)), by="com", with_ties = F) %>%
  dplyr::mutate(n = row_number(), .by="community") %>%
  #mutate(Term = str_replace(.$Term, "(.{10,}?)\\s", "\\1\n")) %>%
  arrange(community) %>%
  mutate(Cluster = paste0(.$community,"^",.$n), .after="community" ) %>%
  dplyr::select(-com, -count) 
  

#### Plot graph ####
pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2"), c( "#364B9A", "#4A7BB7", "#6EA6CD","#98CAE1"), col )

#pal <- set_names(pal[1:length(unique(G_nodes_e$com))], sort(unique(G_nodes_e$com)))
#pal <- set_names(pal[1:length(unique(G_nodes_e$community))], sort(unique(G_nodes_e$community)))

c <- c(  "1",     "10",     "11",     "12",     "13",     "14",     "15",     "16",     "17",     "2",      "3",      "4",      "5",      "6",      "7",      "8",      "9")
col <- c("#FFD92F","#7CAE00","#E41A1C","#377EB8","#4DAF4A","#F781BF","#B3CDE3","#BEAED4","#BF5B17","#CCEBC5","#FBB4AE","#DECBE4","#984EA3","#F8766D","#FED9A6","#66C2A5","#00A9FF")
pal_epi <- set_names(col, c)

# 'graphopt', 'stress', 'kk', 'mds', 'lgl', 'backbone', 
set.seed(1);p <- G_epi %>%
  left_join(., dplyr::select(g_table_epi, pathway, lab="Term", n), by="pathway") %>%
  
  ggraph(layout = 'backbone') + # 'graphopt'
  geom_edge_fan(width = .3, aes(colour=width), show.legend = FALSE) +
  scale_edge_colour_gradient(low="#DFDFDF", high = "#616161") + # colors=c("#DFDFDF" ,"#616161")
  new_scale_color() + theme_graph() +
  #geom_node_point(aes(color="red")) +
  geom_node_point(aes(color=community), size=2) + # community 
  #geom_node_point(aes(color=btwn), size=2) + 
  scale_colour_manual(values = pal_epi) +
  guides(colour = guide_legend(byrow = TRUE)) +
  #geom_node_text(aes(label = n, colour="white"), show.legend = F, nudge_y = .4) + #repel = TRUE, 
  geom_node_text(aes(label = n), show.legend = F,  nudge_x = .4) +
  theme(plot.margin = unit(c(-.1,0,-.2,-.2), "inches"),
        legend.title = element_blank(),
        legend.key.size = unit(0.008, "cm"),
        legend.spacing.y = unit(.1, 'cm'),
        legend.position = c(.96, .70),
        legend.justification = c("right", "top") ) 
  #facet_nodes(~gene_module)
  
# dev.new(width=5.5, height=5, noRStudioGD = TRUE)
# ggsave(paste0("./Figures/S2&3/","./enrichment-graph-epi.png"), plot = p, width = 6, height = 5, dpi = 1000 )
p
```

<img src="../Figures/S2&amp;3/S1a_enrichment-graph-epi.png"
data-fig-align="center" />

``` r
dge <- DEGs_table %>% 
  filter(grepl("[\\d]", .$layers)) %>% 
  filter(p_val<=0.01) %>%
  #filter(pct.1 > .4 & p_val < 0.01)  
  split(., ~layers) %>%  map(., ~.x$gene)

# because i filtered away the terms (n=38) with no match in the GOTERM db, the result is different 
# when not removing these, much more messy plot
GO_comb_sub <- sig_sub_GO %>% 
  #mutate(., clus = factor(l[as.character(.$clus)], levels = ord1)) %>%
  add_count(pathway) %>% 
  #mutate(Term = str_replace_all(tolower(.$pathway), "_", " ")) %>%
  #mutate(Term = str_replace(.$Term, "^gobp ", "")) %>%
  dplyr::rename(GOid = "id") %>%
  left_join( dplyr::select(goterms, Term, t, id), by=c("pathway"="t")) %>%
  # filter leading edge by parameters set above:
  #mutate("p<0.01" = pmap(., ~intersect(..10, dge[[..1]]))) %>%
  # combine same terms from different clusters:
  mutate("leadingEdge_comb" = list(unique(unlist(.$leadingEdge[cur_group_rows()]))), .by="pathway") %>%
  mutate(pathway = paste0(.$clus, ":", .$id)) %>%
  mutate("Sig. in" = paste0(.$clus[cur_group_rows()], collapse = "|"), .by="Term") %>%
  # filter((is.na(.$id))) %>%
  filter(!(.$GOid %in% go_filt)) %>%
  # removes the idential terms:
  slice_max(n=1, order_by = tibble(abs(NES), size), by="Term") %>%
  ungroup()
  
#### get similarity matrix and create graph ####
# weight is the percentage overlap between the two terms
m <- get_similarity_matrix(GO_comb_sub, geneSets=GO_database, semData=d) 
set.seed(123); graph <- build_emap_graph(sim_m=m, fgesa_df=GO_comb_sub, cex_line=1, min_edge=0.1)
# G <- as_tbl_graph(graph$igraph)

#### find number of edges between nodes: ####
m %>% as_tibble(., rownames="from")
```

    # A tibble: 730 × 731
       from         `8:GO:0045229` `8:GO:0030199` `8:GO:0006091` `8:GO:0006119`
       <chr>                 <dbl>          <dbl>          <dbl>          <dbl>
     1 8:GO:0045229             NA          0.319              0          0    
     2 8:GO:0030199             NA         NA                  0          0    
     3 8:GO:0006091             NA         NA                 NA          0.417
     4 8:GO:0006119             NA         NA                 NA         NA    
     5 8:GO:0050853             NA         NA                 NA         NA    
     6 8:GO:0009060             NA         NA                 NA         NA    
     7 8:GO:0032963             NA         NA                 NA         NA    
     8 3:GO:0019724             NA         NA                 NA         NA    
     9 8:GO:0045333             NA         NA                 NA         NA    
    10 8:GO:0002460             NA         NA                 NA         NA    
    # ℹ 720 more rows
    # ℹ 726 more variables: `8:GO:0050853` <dbl>, `8:GO:0009060` <dbl>,
    #   `8:GO:0032963` <dbl>, `3:GO:0019724` <dbl>, `8:GO:0045333` <dbl>,
    #   `8:GO:0002460` <dbl>, `8:GO:0002449` <dbl>, `8:GO:0044281` <dbl>,
    #   `8:GO:0002443` <dbl>, `8:GO:0042742` <dbl>, `8:GO:0015980` <dbl>,
    #   `8:GO:0098542` <dbl>, `8:GO:0019731` <dbl>, `3:GO:0051254` <dbl>,
    #   `8:GO:0019637` <dbl>, `8:GO:0002250` <dbl>, `8:GO:0002252` <dbl>, …

``` r
summary <- graph$tidy %>% activate(edges) %>% 
  as_tibble() %>%
  summarise(l = sum(weight), # n is number of shared genes with other terms, l is the sum of weights
            n_genes = sum(weight > 0), .by = "from" )  %>% 
  arrange(n_genes)
s <- set_names(summary$n_genes, as.character(summary$from))

# graph$tidy %>% activate(edges) %>% as_tibble() %>% filter(weight < 0)
# https://tidygraph.data-imaginist.com/reference/index.html

#### Add info to Graph ####
# https://bookdown.org/markhoff/social_network_analysis/centrality.html
set.seed(123);G <- graph$tidy %>% 
  mutate(gene_module = str_extract(as_tibble(.)$pathway, "[^:]+") ) %>%
  left_join(., dplyr::select(GO_comb_sub, size, pathway, "Sig. in"), by="pathway") %>%
  mutate(pathway = str_extract(as_tibble(.)$pathway, "GO:\\d+") ) %>%
  mutate(degree = centrality_degree(),
         btwn = centrality_betweenness(),
         closness = centrality_closeness(),
         community = as.character(group_louvain())) %>%
  mutate(community = factor( as_tibble(.)$community, 
                             levels = as.character(1:max(as.integer(as_tibble(.)$community))) )) %>%
  left_join(goterms, by=c("pathway"="id")) %>%
  #mutate(GOslim = GOslim_l[as_tibble(.)$pathway]) %>%
  
  #filter(community == "1" | community == "2" | community == "3" )  %>% #  | community == "3" | community == "4" | community == "5"
  {. ->> G_nodes} %>%
  activate(edges) %>% 
  left_join(summary, by="from") %>%
  mutate(term_from = .N()$Term[from]) %>%
  mutate(term_to = .N()$Term[to]) %>%
  mutate(comm = .N()$community[from]) %>%
  {. ->> G_edges} %>%
  activate(nodes)

#### get table version of graph ###
G_nodes_s <- G_nodes %>% as_tibble() # %>% split(., ~community)
G_edges <- G_edges %>% as_tibble() 
Reduce(intersect, as.list(GOBPANCESTOR[na.omit(G_nodes$pathway)]))
```

    NULL

``` r
G_edges %>%
  #filter(weight > .6) 
  filter(!(comm == "1" & weight > .2))
```

    # A tibble: 6,801 × 9
        from    to weight width     l n_genes term_from                term_to comm 
       <int> <int>  <dbl> <dbl> <dbl>   <int> <chr>                    <chr>   <fct>
     1     1     2  0.319 1.26   1.54      10 external encapsulating … collag… 4    
     2     1     7  0.173 0.931  1.54      10 external encapsulating … collag… 4    
     3     1    29  0.122 0.780  1.54      10 external encapsulating … skelet… 4    
     4     1    32  0.125 0.791  1.54      10 external encapsulating … connec… 4    
     5     1    35  0.111 0.745  1.54      10 external encapsulating … cartil… 4    
     6     1    40  0.236 1.09   1.54      10 external encapsulating … extrac… 4    
     7     1    90  0.106 0.727  1.54      10 external encapsulating … ossifi… 4    
     8     1    97  0.128 0.801  1.54      10 external encapsulating … extrac… 4    
     9     1   126  0.101 0.711  1.54      10 external encapsulating … cellul… 4    
    10     1   601  0.115 0.758  1.54      10 external encapsulating … regula… 4    
    # ℹ 6,791 more rows

``` r
#### filter graph ####
# https://bookdown.org/markhoff/social_network_analysis/centrality.html
t <- c("phagocytosis")
set.seed(123);G_sub <- G %>% 
  #filter(!(grepl(paste(t,collapse="|"), Term))) %>%
  activate(edges) %>%
  #filter(grepl("|^1$|^2$|", comm) & weight < .2) %>%
  #filter(!(comm == "1" & weight > .8)) %>%
  filter(!(comm == "1" & weight < .3 )) %>%
  filter(!(comm == "3" & weight < .3 )) %>%
  filter(!(comm == "2" & weight < .3 )) %>%
  filter(!(comm == "4" & weight < .3 )) %>%
  filter(!(comm == "5" & weight < .3 )) %>%
  #filter(!(comm == "6" & weight < .3 )) %>%
  filter(!(comm == "7" & weight < .3 )) %>%
  #filter(weight > .3) %>%
  #filter(grepl("|^1$|", comm) & weight > .6) %>%
  arrange(desc(weight)) %>%
  #filter(n > 10) %>%
  # mutate(width = .["width"]*.1) %>%
  activate(nodes) %>%
  # re-calculate degree after filtering and remove nodes with no connections
  # before re-clustering
  mutate(degree = centrality_degree()) %>%
  filter(degree > 0) %>%
  mutate(com = as.character(group_louvain()))

G_nodes_sub <- G_sub %>% as_tibble() 

# select terms to label in graph and diplay in table:
g_table_sub <- G_sub %>% as_tibble() %>%
  # mutate(col = pal[.$com]) %>%
  dplyr::select(community, com, pathway, Term, degree) %>% unique() %>% arrange(community) %>%
  mutate(count = n_distinct(Term), .by = c("com")) %>%
  # determines the number of rows are included in the table:
  filter(count > 4) %>%
  left_join(., dplyr::select(GO_comb_sub, id, ES), by=c("pathway"="id")) %>%
  # mutate(ES = ) %>%
  slice_max(n=2, order_by = tibble(degree, abs(ES)), by="com", with_ties = F) %>%
  slice_max(n=1, order_by = tibble(abs(ES)), by="com", with_ties = F) %>%
  dplyr::mutate(n = row_number(), .by="community") %>%
  mutate(community = factor( .$community, levels = as.character(1:17)) ) %>%
  #mutate(Term = str_replace(.$Term, "(.{10,}?)\\s", "\\1\n")) %>%
  arrange(community) %>%
  mutate(Cluster = paste0(.$community,"^",.$n), .after="community" ) 

#### Plot graph ####
library(ggnewscale)
pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2"), c( "#364B9A", "#4A7BB7", "#6EA6CD","#98CAE1") )

#pal <- set_names(pal[1:length(unique(G_nodes_sub$com))], sort(unique(G_nodes_sub$com)))
#pal <- set_names(pal[1:length(unique(G_nodes_sub$community))], sort(unique(G_nodes_sub$community)))
c <- c(  "1",     "10",     "11",     "12",     "13",     "14",     "15",     "16",     "17",     "2",      "3",      "4",      "5",      "6",      "7",      "8",      "9",
         "18",    "19",     "20",     "21",     "22",     "23" )
col <- c("#FFD92F","#7CAE00","#E41A1C","#377EB8","#4DAF4A","#F781BF","#B3CDE3","#BEAED4","#BF5B17","#CCEBC5","#FBB4AE","#DECBE4","#984EA3","#F8766D","#FED9A6","#66C2A5","#00A9FF",  "#FC8D62","#8DA0CB","#364B9A","#A6D854","#FFF2AE","#E5C494")
pal_sub <- set_names(col, c)

# 'graphopt', 'stress', 'kk', 'mds', 'lgl', 'backbone', 
set.seed(1);p <- G_sub %>%
  left_join(., dplyr::select(g_table_sub, pathway, lab="Term", n), by="pathway") %>%
  #filter(!is.na(n)) %>% as_tibble() -> l
  
  ggraph(layout = 'backbone') + # 'graphopt'
  geom_edge_fan(width = .3, aes(colour=width), show.legend = FALSE) +
  scale_edge_colour_gradient(low="#DFDFDF", high = "#616161") + # colors=c("#DFDFDF" ,"#616161")
  new_scale_color() + theme_graph() +
  #geom_node_point(aes(color="red")) +
  geom_node_point(aes(color=community), size=2) + # community 
  #geom_node_point(aes(color=btwn), size=2) + 
  scale_colour_manual(values = pal_sub) +
  guides(colour = guide_legend(byrow = TRUE)) +
  #geom_node_text(aes(label = n, colour="white"), show.legend = F, nudge_y = .4) + #repel = TRUE, 
  geom_node_text(aes(label = n), show.legend = F,  nudge_x = .4) +
  theme(plot.margin = unit(c(-.1,0,-.2,-.2), "inches"),
        legend.title = element_blank(),
        legend.key.size = unit(0.008, "cm"),
        legend.spacing.y = unit(.1, 'cm'),
        legend.position = c(.96, .70),
        legend.justification = c("right", "top") ) 
  #facet_nodes(~gene_module)
  
# dev.new(width=5.5, height=5, noRStudioGD = TRUE)
# ggsave(paste0("./Figures/S2&3/","./enrichment-graph-submuc.png"), plot = p, width = 6, height = 5, dpi = 1000 )
p
```

<img src="../Figures/S2&amp;3/S1b_enrichment-graph-submuc.png"
data-fig-align="center" />

``` r
m <- set_names(c("#E41A1C","#FF7F00", "#C77CFF","#984EA3"), c("Sup", "Upper", "Lower", "Basal"))

set.seed(1);p <- imap(m, 
  ~mutate(G_epi, col = factor(c(.x, "gray"))[ifelse(grepl(.y, G_nodes_epi$`Sig. in`), 1, 2)]) ) %>%
  
  imap(., ~ggraph(.x, layout = 'backbone') + # 'graphopt'
  geom_edge_fan(width = .1, aes(colour=width), show.legend = FALSE) +
  scale_edge_colour_gradient(low="#DFDFDF", high = "#616161") + # colors=c("#DFDFDF" ,"#616161")
  new_scale_color() +
  geom_node_point(aes(color=col), size=.05) + # community `Sig. in`
  scale_colour_identity() +
  theme_graph() + theme(plot.margin = unit(c(0,2,0,0), "pt")) ) 
  #facet_nodes(~gene_module)
p <- wrap_plots(p, ncol = 1)

# ggsave(filename = "./Figures/S2/graph_per_layer_epi.pdf", plot = p, width =1 , height =3.3, dpi = 1000 )

m <- set_names(c("#00BFC4","#00A9FF","#377EB8","#CD9600"), c("8", "3","4", "0"))

set.seed(1);p <- imap(m, 
  ~mutate(G_sub, col = factor(c(.x, "gray"))[ifelse(grepl(.y, G_nodes_sub$`Sig. in`), 1, 2)]) ) %>%
  
  imap(., ~ggraph(.x, layout = 'backbone') + # 'graphopt'
  geom_edge_fan(width = .1, aes(colour=width), show.legend = FALSE) +
  scale_edge_colour_gradient(low="#DFDFDF", high = "#616161") + # colors=c("#DFDFDF" ,"#616161")
  new_scale_color() +
  geom_node_point(aes(color=col), size=.05) + # community `Sig. in`
  scale_colour_identity() +
  theme_void() + theme(plot.margin = unit(c(0,-10,0,-10), "pt")) )
  #facet_nodes(~gene_module)

# dev.new(width=5, height=15, noRStudioGD = TRUE)
(p <- wrap_plots(p, ncol = 1) )#
```

<img src="../Figures/S2&amp;3/S2a_graph-per-layer.png"
data-fig-align="center" />

``` r
# ggsave(filename = "./Figures/S2&3/graph_per_layer_sub.pdf", plot = p, width =1 , height =3.9, dpi = 1000, bg = "white" )
```

``` r
pal2 <- c(RColorBrewer::brewer.pal(10,"Paired"),
         RColorBrewer::brewer.pal(8,"Set2"),
         RColorBrewer::brewer.pal(8,"Dark2") )
         #RColorBrewer::brewer.pal(9,"Pastel1"),
         #RColorBrewer::brewer.pal(8,"Pastel2")
  
G_nodes <- list("epi"=G_nodes_e, "SubMuc"=G_nodes_s)
pal_list <- list("epi"=pal_epi, "SubMuc"=pal_sub)

df <- imap(slim_df$slim, ~ .x %>% # slim_df$slim[[2]] %>%
          filter(!(Percent == 0)) %>%
          arrange(desc(Percent)) %>%
          mutate(Term = factor(.$Term, levels = .$Term)) %>%
          mutate(ystart = c(0, cumsum(Percent)[1:nrow(.)-1])) %>%
          mutate(yend = lead(ystart, default = 100)) %>% 
          {. ->> temp}  %>%
          unnest(go_terms) %>%
          left_join(., select(G_nodes[[.y]], community, go_terms="Term"), by=c("go_terms")) %>%
          #left_join(., select(G_nodes[[.y]], community, go_terms="Term"), by=c("go_terms")) %>%
          summarise(count = n_distinct(go_terms), .by = c("Term", "Count", "Percent","community", "ystart", "yend")) %>%
          #summarise(count = n_distinct(go_ids), .by = c("Term", "Count", "Percent","community", "ystart", "yend")) %>%
          mutate(Percent = count/sum(temp$Count)*100) %>%
          mutate(ystart2 = c(0, cumsum(Percent)[1:nrow(.)-1]/sum(Percent))*100) %>%
          mutate(yend2 = lead(ystart2, default = 100))  %>%
          mutate(community = ifelse(as.integer(.$community) > length(pal_list[[.y]]), NA, as.character(.$community))) %>%
          mutate(community = factor( .$community, levels = as.character(1:length(pal_list[[.y]]))) ) )


pie_plots <- map2(df, pal_list,
                 ~ggplot(.x) + 
                # ajust the width of the circle by changing max and min x values
                # to get a rageular pice chart set xmin=0
                geom_rect(aes(fill=Term, ymax=yend, ymin=ystart, xmax=2.7, xmin=1)) +
                scale_fill_manual(values = pal2) + new_scale_fill() +
                geom_rect(aes(fill=community, ymax=yend2, ymin=ystart2, xmax=3.3, xmin=2.7)) +
                scale_fill_manual(values = .y, na.value = "gray", guide = "none") + # , guide = "none"
                xlim(c(0, 4)) + 
                theme_void() +
                theme(aspect.ratio=1,
                      legend.title = element_blank(),
                      legend.box.margin=margin(0,10,-0,-30)) + 
                coord_polar(theta="y")  )
pie_plots[[1]]
```

<img src="../Figures/S2&amp;3/S2a_piechart-GOslim-graph-clus.png"
data-fig-align="center" />

``` r
pie_plots[[2]]
```

<img src="../Figures/S2&amp;3/S2a_piechart-GOslim-graph-clus-2.png"
data-fig-align="center" />

``` r
# ggsave(filename = "./Figures/S2&3/piechart-GOslim-graph-clus_epi.pdf", plot = pie_plots[[1]], width =10 , height =10, dpi = 1000 )
# ggsave(filename = "./Figures/S2&3/piechart-GOslim-graph-clus_submuc.pdf", plot = pie_plots[[2]], width =10 , height =10, dpi = 1000 )
```

``` r
# Coloring the table conditionally using `ggpubr::table_cell_bg()`
tabel_col.fun<- function(ggtab, g, pal){
  for(i in 1:nrow(g)){
    #print(i)
    col <- g$community[i]
    #print(col)
    ggtab <- table_cell_bg(ggtab, 
                           row = i+1, column = 1, 
                           fill = pal[as.character(col)], color = pal[as.character(col)] )
    # ggtab <- table_cell_bg(ggtab, i+1, column = 4, linewidth = 1,
    #                        color = pal_epi[as.character(col)])
  }
  print(ggtab)}

base_size <- 8
g_table <- list("epi"=g_table_epi, "SubMuc"=g_table_sub)
pal_list <- list("epi"=pal_epi, "SubMuc"=pal_sub)
options(scipen=999)

ggtab <- imap(g_table, ~ .x[,-1] %>%
  select(-any_of(c("com", "count","n"))) %>%
  mutate(ES = format(.$ES, digits = 2 ,scientific = FALSE) ) %>%
  #mutate(ES = factor(.$ES) ) %>%
  ggtexttable(., rows = NULL, 
            theme = # ttheme("light"))
            ttheme(base_style = "light",
                  base_colour = "#616161",
                  padding = unit(c(2, 2), "mm"),
                  #core = list(fg_params = list(parse=TRUE)),
                  colnames.style = colnames_style(size = base_size+2, fill = "white", color = "#616161"),
                  rownames.style = rownames_style(size = base_size+1),
                  tbody.style = 
                    tbody_style(
                      color = "#616161", size = base_size, fill = "transparent", #parse = T,
                      hjust = as.vector(matrix(c(.5, 0, .5, .5, 1), ncol = 5, nrow = nrow(.), byrow = TRUE)),
                      x = as.vector(matrix(c(0.5, .1, .5,.5, .9), ncol = 5, nrow = nrow(.), byrow = TRUE)) 
                  ) ))  %>%
  tabel_col.fun(., g_table[[.y]], pal_list[[.y]]) )
```

<img src="../Figures/S2&amp;3/S2b_enrichment_table.png"
data-fig-align="center" />

<img src="../Figures/S2&amp;3/S2b_enrichment_table-2.png"
data-fig-align="center" />

``` r
#ggtab[[1]]
#ggtab[[2]]

# dev.new(width = 5.5, height = 5)
# ggsave(paste0("./Figures/S2&3/","enrichment_table_epi.png"), plot = ggtab[[1]], width = 5, height = 3, dpi = 1000 )
# ggsave(paste0("./Figures/S2&3/","enrichment_table_submuc.png"), plot = ggtab[[2]], width = 6, height = 5.5, dpi = 1000 )
```
