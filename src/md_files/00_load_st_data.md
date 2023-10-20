Load Spatial data
================
10/20/23

### Load libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(tidyseurat)
library(Seurat)
library(hdf5r)
# remotes::install_github("czarnewski/niceRplots",force=T)
library(niceRplots)

library(xml2) # loads the image 
library(sp)
source("../bin/help_functions.R")
source("../bin/spatial_visualization.R")

#################
# COLOR PALETTS #
#################
pal <- rep(c(RColorBrewer::brewer.pal(9,"Set1"),
         RColorBrewer::brewer.pal(9,"Pastel1"),
         RColorBrewer::brewer.pal(8,"Accent"),
         RColorBrewer::brewer.pal(8,"Set2"),
         RColorBrewer::brewer.pal(8,"Pastel2") ,
         scales::hue_pal()(8)),99)
```

### Load Visium data

``` r
#########
# PATHS #
#########
input_dir <- "../data/spatial_data"
result_dir <- "../results/00_load_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
meta <- read_csv("../data/Clinical_data/Clinical_data_Spatial_DMPA.csv")
h5_files <- list.dirs(path = input_dir,
                      full.names = T, recursive = T) %>%
            grep("P\\d\\d\\d$", ., value = TRUE) %>%
            set_names(., str_extract(., "P\\d\\d\\d$"))

#sample_names <- names(h5_files)
image <- map(h5_files, 
             ~Read10X_Image(
               filter.matrix = T,
               image.dir = paste0(.x, "/spatial"),
               image.name = "tissue_hires_image.png"))

# Read in h5 files and create Seurat Object
seuratObj_list <- pmap(list(h5_files, image, names(h5_files)),
                       ~Load10X_Spatial(
                         filename = "filtered_feature_bc_matrix.h5",
                         filter.matrix = T,
                         assay = "RNA",
                         data.dir = ..1,
                         image =  ..2,
                         slice = ..3)) 

sample_id <- c("P031", "P080", "P097", "P105", "P107", "P108", "P114", 
"P118") %>% set_names()
```

### Tidy up the seurat object

``` r
####################################
# RENAME SAMPLES, SPOTS AND IMAGES #
####################################
seuratObj_list <- seuratObj_list %>%
  imap(., ~AddMetaData(object = .x, 
                       metadata = rep(.y, length(Idents(.x))), 
                       col.name = "orig.ident")) %>%
  map(.,  ~SetIdent(., value = .@meta.data$orig.ident)) %>%
  imap(., ~RenameCells(.x, 
                      new.names = paste0(.y,"_", gsub("-.*","",colnames(.x[["RNA"]])))) ) %>%
  imap(., ~ {.x@images <- set_names(.@images,.y); .x})

##################
# MERGE SAMPLES #
#################
# Merge datasets into one single seurat object
DATA  <- merge(seuratObj_list[[1]], y = seuratObj_list[2:length(seuratObj_list)])
```

### Load morphology annotation

``` r
##########################
# ADD MANUAL ANNOTATION #
#########################
# scale.factors <- map(DATA@images, ~pluck(.x, "scale.factors")$hires)
# img_coord <- map(DATA@images, ~pluck(.x, "coordinates"))
#   
# img_coord <- img_coord %>%
#   map2(., scale.factors, 
#        ~mutate(.x, imagecol = .$imagecol * .y,
#                    imagerow = .$imagerow * .y) )
  
a2 <- map(sample_id, ~xml2::as_list( read_xml( paste0( input_dir,"/",.x,"/",.x,".svg")) ) )
dd <- map(a2, ~as.numeric( strsplit( attributes(.x$svg)$viewBox , " ")[[1]] ))
dd2 <- map(sample_id, ~dim(DATA@images[[.x]]@image))

# sample_id <- "P031"
# a2 <- a2[["P031"]]
# dd <- dd[["P031"]]
# dd2 <- dd2[["P031"]]

# add image coordinates to the seurat object
get_img_coord <- function(DATA, sample_id){
  img_coord <- DATA@images[[sample_id]]@coordinates
  DATA@images[[sample_id]]@coordinates <<- img_coord 
}

# add manual spatial annotation
get_sp_annot <- function(a2, dd, dd2, sample_id){
  scale.factor <- DATA@images[[sample_id]]@scale.factors$hires
  img_coord <- DATA@images[[sample_id]]@coordinates
  
  id <- a2$svg %>%
    map_chr(., ~attr(.x,"id")) %>% 
    set_names(seq_along(.), .)
  
  annot_coord <- id %>%
    map(., ~get_shape(a2$svg[[.x]]) ) %>%
    map(., ~as_tibble(.x, .name_repair="unique"))  %>%
    map(., ~mutate(.x, x = .x[[1]]*dd2[2]/dd[3],
                       y = .x[[2]]*dd2[1]/dd[4] )) %>%
    map(., ~rowid_to_column(., var = "path_idx")) %>%
    {. ->>  temp} %>%
    bind_rows(., .id="name") %>%
    group_by(., name) %>%
    mutate(elem_idx = cur_group_id()) %>% # group_indices(., name)
    ungroup() %>%
    mutate(colour = ifelse(grepl("fov|zoom|full_image",.$name),
                           "transparent", "black")) %>%
    list(.) %>%
    set_names(., sample_id[1])
  
  DATA@tools <<- append(DATA@tools, annot_coord )
  
  img_coord <- temp %>%
    list_modify("fov" = NULL, "full_image" = NULL, "zoom" = NULL) %>%
    compact() %>%
    imap(., ~mutate(img_coord, !!.y := sp::point.in.polygon(
         point.x = img_coord$imagecol*scale.factor,
         point.y = img_coord$imagerow*scale.factor,
         pol.x = .x$x,
         pol.y = .x$y )) ) %>%
    map(., ~select(.x, last_col())) %>%
    cbind(img_coord, .)
  
  DATA@images[[sample_id]]@coordinates <<- img_coord 
  
  sp_annot <- img_coord %>%
    rownames_to_column(., var = "barcodes") %>%
    mutate(across(7:ncol(.), ~ifelse(. == 0, NA, .)) ) %>%
    pivot_longer(., cols = 7:ncol(.), names_to ="sp_annot", values_to = "count") %>%
    filter(!(is.na(count))) %>%
    group_by(barcodes) %>%
    mutate(dupp = row_number()) %>%
    ungroup() %>%
    filter(., .$dupp == 1) %>%
    select(., barcodes, sp_annot)
  
  #DATA <<- left_join( DATA, sp_annot, by=c(".cell"="barcodes", "sp_annot"="sp_annot")) 
  
  return(list(coord=annot_coord, annot=sp_annot))
}
annot <- list(a2, dd, dd2, names(a2)) %>%
  pmap(., ~get_sp_annot(..1, ..2, ..3, ..4)) %>%
  set_names(., names(a2))

sp <- map(annot, 2) %>% bind_rows()
DATA <- left_join( DATA, sp, by=c(".cell"="barcodes")) 

DATA 
```

    # A Seurat-tibble abstraction: 6,700 × 5
    # [90mFeatures=36601 | Cells=6700 | Active assay=RNA | Assays=RNA[0m
       .cell                 orig.ident nCount_RNA nFeature_RNA sp_annot
       <chr>                 <chr>           <dbl>        <int> <chr>   
     1 P031_AAACGAGACGGTTGAT P031              463          356 SubMuc  
     2 P031_AAACTGCTGGCTCCAA P031             4081         2232 SubMuc  
     3 P031_AAAGTAGCATTGCTCA P031             7595         2588 epi     
     4 P031_AAAGTGTGATTTATCT P031            11394         3873 epi     
     5 P031_AAAGTTGACTCCCGTA P031             4617         2270 epi     
     6 P031_AAATACCTATAAGCAT P031             5538         2507 epi     
     7 P031_AAATCGTGTACCACAA P031             1765         1084 SubMuc  
     8 P031_AAATGGCCCGTGCCCT P031              789          477 <NA>    
     9 P031_AAATTAACGGGTAGCT P031              639          465 SubMuc  
    10 P031_AAATTTGCGGGTGTGG P031             3570         1948 SubMuc  
    # ℹ 6,690 more rows

### Add meta data

``` r
meta <- meta %>%
  select(orig.ident="ID", groups=Contraception) %>%
  mutate(groups = ifelse(.$groups=="no HC", "ctrl", .$groups))

DATA <-  DATA %>%
  mutate(sp_annot2 = .$sp_annot) %>%
  mutate(sp_annot = ifelse(grepl("epi", .$sp_annot2), "epi", 
                           ifelse(grepl("SubMuc", .$sp_annot2), "SubMuc", .$sp_annot2 ))) %>%
  left_join(., meta, by="orig.ident") %>%
  select(groups, sp_annot, everything())
```

``` r
DATA %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "sp_annot",
          zoom = "zoom",
          ncol = 2,
          annot_line = .1,
          img_alpha = 0,
          point_size = 0.8
        )
```

<img src="../Figures/00/00a_plot_sp_annot.png" data-fig-align="center" />

### Identify spots with missing morphology annotation

``` r
##########################
# SP MANUAL ANNOTATION #
#########################
keep_epi <- c("P097_CTATTTGCTTGGAGGA", "P031_ATAGAGTTATCAACTT", "P031_GCGGAGAGGGAGAACG", "P031_CCTATCTATATCGGAA", "P031_GGATCTTGACTCAACC","P031_GCCCTAGCCGTCGCGA", "P031_AGCTCTTTACTCAGTT", "P031_AGTCAACACCACCATC", "P031_CGAACCCGCATGCGTC", "P080_TCGAGCCAGGCAGGCC", "P080_AACCCGACAACCCGTG","P080_GGCCGTTTGGGTTTCA","P080_CATTTGAGTGGTACGT", "P080_TTGGTCACACTCGTAA", "P080_CGAACCCGCATGCGTC", "P080_TTGCTGATCATGTTCG", "P080_CGAGACCCTAGAGTGT", "P080_CCAGCCTGGACCAATA", "P080_CGCATGGTGCGATGCT", "P080_ATAGACAACGGGACCT", "P080_CTCATTAACGTTGCCC", "P097_CGAGTTCTGTCCCACC", "P097_GTATGAAATTTCACTC", "P097_AGCAACCGAAAGTAAT", "P097_GAAGCCACTGATTATG", "P097_GCACTGCCTACCTTTA", "P097_CACAGCACCCACGGCA", "P097_AGTTTGGCCAGACCTA", "P097_ACAGAACTGAGAACAA", "P097_GCTTTCAGAGGAGGTG", "P107_ACATCCCGGCCATACG", "P107_ACGCAAACTAATAGAT","P107_ACTTGACTCCCTCTTT", "P107_GATCTTGGAGGGCATA", "P105_AAGACTGCAAGCTACT", "P105_ACACGGGAACTTAGGG", "P105_CACATTCTTTCGATGG", "P105_CCCGACCATAGTCCGC", "P105_CCTATGGGTTACCGTC", "P105_CCTCTAATCTGCCAAG", "P105_CGAGGCTAAATATGGC", "P105_CGTTTCACTTCGGGCG", "P105_GCGCTAATTGAATAGA", "P105_GATATGCGGTAGCCAA", "P105_GGCTCTGCTCCAACGC", "P105_GGGCTGCCTAGGGCGA", "P105_GGCTCTGCTCCAACGC", "P105_GGGCTGCCTAGGGCGA", "P105_GGTTTACAATCTCAAT", "P105_TGTGGCGGGCTTCTGG", "P105_TTAATCAGTACGTCAG", "P105_TTCATGGCGCAACAGG", "P105_TTGGGACACTGCCCGC", "P114_ACCTGCGTGTCATGTT", "P114_CAAGGATCGCATGTTC", "P114_AATGTTGTCGTGAGAC", "P108_ATACGCCGGCGAAACC", "P108_CACGTCGGCAACCTCT", "P108_TCGCCGAAGTTGCGTC", "P108_TTGATTAGCTGTTTCT", "P118_AACCCGACAACCCGTG", "P118_TGAGCCATACAGTCTC", "P118_ATGGGACCTGCTGAAC", "P118_AGCTAACAAGCAATGT", "P118_TTGATTAGCTGTTTCT", "P118_CGTGTCTCGTTACGAC", "P118_GATCAACATAAAGGGA", "P118_TAACAGCGTTTGTGCT", "P118_TGAGCCATACAGTCTC", "P118_GATTACTGAATTTGGG", "P118_TCCAACTTTAAATTCT", "P118_CGATCCTCGCAACATA", "P118_AACCCGACAACCCGTG", "P118_TCTCTTACCGCGAACC", "P118_TATGTAAAGTGCTTAA", "P118_GAAGTTTCCACTCAAT", "P118_TAGTCCCGGAGACCAC", "P118_CGGCCAGAGCGACCAT", "P118_ATAAAGGCTCGGTCGT", "P118_TATTCGTGCCAGAATA", "P118_GGGAGTTAATGAGGCG", "P118_CCGGGCGGTCTCGTCA") 
keep_SubMuc <- c("P105_TGACATCGAGCGGACC", "P118_GTATCAAACGTTAGCT", "P097_AGGTGGTGACCTTCGC", "P097_AGCCCGGCATTAGAGG", "P031_CTAGTTGGGCCCGGTA","P031_TCTACCGTCCACAAGC", "P031_AAATGGCCCGTGCCCT", "P031_AATCTGGCTTTCTAGT", "P031_AAATGGCCCGTGCCCT", "P031_AATCTGGCTTTCTAGT", "P080_AGGCATTGTCGTAGGG", "P097_GGCAGCAAACCTATGC", "P097_ACAAAGCATGACCTAG", "P097_CCTCCTGTTGTGTCGT", "P097_CGTCGGATAGTGTTGA", "P107_TCTTGATGCGTAGCGA", "P105_ACGATCATCTTGTAAA", "P105_AGATATAATACGACTA", "P105_CTAGGTCTGAAGGAAT", "P114_AAGGGTTTGATTTCAG", "P114_GCACGTGGTTTACTTA", "P114_TAGGCCTATATAGTCT", "P107_TCTTCGATACCAATAA", "P108_GTGCGAAATCGAACAC", "P118_TCGAGACCAACACCGT", "P118_TCGGAGTACATGAGTA", "P118_TATTCAATTCTAATCC", "P118_AGGAGGCCTTCGCGCG") 

###############################
# IDENTIFY MISSING ANNOTATION #
###############################
df <- map(sample_id, ~pluck(DATA@images, .x, "coordinates")) %>%
  bind_rows() %>%
  #cbind(.,as_tibble(select(DATA, filt))) %>%
  cbind(.,as_tibble(select(DATA, orig.ident))) %>%
  cbind(.,as_tibble(select(DATA, nCount_RNA))) %>%
  cbind(.,as_tibble(select(DATA, nFeature_RNA))) %>%
  cbind(.,as_tibble(select(DATA, sp_annot))) %>%
  rownames_to_column(var = "barcode") %>%
  as_tibble() %>%
  select(-starts_with(c("epi_", "SubMuc_")))

f <- df %>% 
  filter(orig.ident == "P097" & is.na(.$sp_annot)) #%>%
  #filter(orig.ident == "P118" & filt == "keep")

# dev.new(width=10, height=10, noRStudioGD = TRUE)
# dev.new(width=5, height=3, noRStudioGD = TRUE)
# DATA %>%
#   mutate(sp_annot = case_when(colnames(DATA) %in% keep_epi ~ 'epi',
#                               colnames(DATA) %in% keep_SubMuc ~ 'SubMuc',
#                               TRUE ~ .$sp_annot)) %>%
#   mutate(filt = case_when(is.na(.$sp_annot) & nCount_RNA < 500 ~ 'filt',
#                           colnames(DATA) %in% filt ~ 'filt',
#                           TRUE ~ 'keep')) %>%
#   filter(grepl("P097_TAGAATAGCCGATGAA", `.cell`)) %>%
#   plot_st_meta.fun(.,  
#           assay="RNA",
#           feat = "filt",
#           zoom = "full_image",
#           ncol = 2,
#           annot_line = .1,
#           img_alpha = 0,
#           point_size = 1
#         )
```

### Add missing morphology annotation

``` r
# dev.new(width=10, height=10, noRStudioGD = TRUE)
DATA <-  DATA %>%
  mutate(sp_annot = case_when(colnames(DATA) %in% keep_epi ~ 'epi',
                              colnames(DATA) %in% keep_SubMuc ~ 'SubMuc',
                              TRUE ~ .$sp_annot))

DATA %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "sp_annot",
          zoom = "zoom",
          ncol = 2,
          annot_line = .1,
          img_alpha = 0,
          point_size = 0.8
        )
```

<img src="../Figures/00/00b_plot_new_sp_annot.png"
data-fig-align="center" />

## Plot spots to be removed

``` r
# dev.new(width=10, height=10, noRStudioGD = TRUE)

DATA <-  DATA %>%
  mutate(filt = case_when(is.na(.$sp_annot)  ~ 'filt',
                          TRUE ~ 'keep')) 
DATA %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "filt",
          zoom = "full_image",
          ncol = 2,
          annot_line = .1,
          img_alpha = 0,
          point_size = 1
        )
```

<img src="../Figures/00/00c_plot_spots_to_remove.png"
data-fig-align="center" />

### Remove NA spots

``` r
dim(DATA)
```

    [1] 36601  6700

``` r
# Filter spots outside manual annotation
DATA <- DATA[, !(is.na(DATA$sp_annot))]
DATA$filt <- NULL

dim(DATA)
```

    [1] 36601  6612

``` r
DATA
```

    # A Seurat-tibble abstraction: 6,612 × 7
    # [90mFeatures=36601 | Cells=6612 | Active assay=RNA | Assays=RNA[0m
       .cell            groups sp_annot orig.ident nCount_RNA nFeature_RNA sp_annot2
       <chr>            <chr>  <chr>    <chr>           <dbl>        <int> <chr>    
     1 P031_AAACGAGACG… ctrl   SubMuc   P031              463          356 SubMuc   
     2 P031_AAACTGCTGG… ctrl   SubMuc   P031             4081         2232 SubMuc   
     3 P031_AAAGTAGCAT… ctrl   epi      P031             7595         2588 epi      
     4 P031_AAAGTGTGAT… ctrl   epi      P031            11394         3873 epi      
     5 P031_AAAGTTGACT… ctrl   epi      P031             4617         2270 epi      
     6 P031_AAATACCTAT… ctrl   epi      P031             5538         2507 epi      
     7 P031_AAATCGTGTA… ctrl   SubMuc   P031             1765         1084 SubMuc   
     8 P031_AAATGGCCCG… ctrl   SubMuc   P031              789          477 <NA>     
     9 P031_AAATTAACGG… ctrl   SubMuc   P031              639          465 SubMuc   
    10 P031_AAATTTGCGG… ctrl   SubMuc   P031             3570         1948 SubMuc   
    # ℹ 6,602 more rows

``` r
# dev.new(width=5, height=2, noRStudioGD = TRUE)
DATA %>%
  #filter(orig.ident == "P118" | orig.ident == "P097") %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "sp_annot",
          zoom = "zoom",
          ncol = 2,
          annot_line = .1,
          annot_col = "black",
          img_alpha = 0,
          point_size = 0.8
        )
```

<img src="../Figures/00/final-sp_annot.png" data-fig-align="center" />

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(DATA, paste0(result_dir,"seuratObj_merged.RDS"))
#DATA <- readRDS(paste0(result_dir,"seuratObj_merged.RDS"))
```

### Session info

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
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] xml2_1.3.3         niceRplots_0.1.0   hdf5r_1.3.8        Seurat_4.3.0      
     [5] tidyseurat_0.5.3   SeuratObject_4.1.3 sp_1.5-1           ttservice_0.2.2   
     [9] forcats_1.0.0      stringr_1.5.0      dplyr_1.1.2        purrr_1.0.1       
    [13] readr_2.1.3        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.4.3     
    [17] tidyverse_1.3.2   

    loaded via a namespace (and not attached):
      [1] readxl_1.4.1           backports_1.4.1        plyr_1.8.8            
      [4] igraph_1.4.1           lazyeval_0.2.2         splines_4.1.2         
      [7] listenv_0.9.0          scattermore_0.8        digest_0.6.31         
     [10] htmltools_0.5.5        fansi_1.0.4            magrittr_2.0.3        
     [13] tensor_1.5             googlesheets4_1.0.1    cluster_2.1.4         
     [16] ROCR_1.0-11            tzdb_0.3.0             globals_0.16.2        
     [19] modelr_0.1.10          matrixStats_0.63.0     vroom_1.6.0           
     [22] timechange_0.2.0       spatstat.sparse_3.0-0  colorspace_2.1-0      
     [25] rvest_1.0.3            ggrepel_0.9.3          haven_2.5.1           
     [28] xfun_0.38              crayon_1.5.2           jsonlite_1.8.5        
     [31] progressr_0.13.0       spatstat.data_3.0-0    survival_3.5-5        
     [34] zoo_1.8-11             glue_1.6.2             polyclip_1.10-4       
     [37] gtable_0.3.4           gargle_1.2.1           leiden_0.4.3          
     [40] future.apply_1.10.0    abind_1.4-5            scales_1.2.1          
     [43] DBI_1.1.3              spatstat.random_3.0-1  miniUI_0.1.1.1        
     [46] Rcpp_1.0.10            viridisLite_0.4.2      xtable_1.8-4          
     [49] reticulate_1.28        bit_4.0.5              htmlwidgets_1.6.2     
     [52] httr_1.4.5             RColorBrewer_1.1-3     ellipsis_0.3.2        
     [55] ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3       
     [58] uwot_0.1.14            dbplyr_2.2.1           deldir_1.0-6          
     [61] utf8_1.2.3             labeling_0.4.3         tidyselect_1.2.0      
     [64] rlang_1.1.1            reshape2_1.4.4         later_1.3.0           
     [67] munsell_0.5.0          cellranger_1.1.0       tools_4.1.2           
     [70] cli_3.6.1              generics_0.1.3         broom_1.0.4           
     [73] ggridges_0.5.4         evaluate_0.21          fastmap_1.1.1         
     [76] yaml_2.3.7             goftest_1.2-3          knitr_1.42            
     [79] bit64_4.0.5            fs_1.6.2               fitdistrplus_1.1-8    
     [82] RANN_2.6.1             pbapply_1.7-0          future_1.32.0         
     [85] nlme_3.1-163           mime_0.12              compiler_4.1.2        
     [88] rstudioapi_0.14        plotly_4.10.1          png_0.1-8             
     [91] spatstat.utils_3.0-1   reprex_2.0.2           stringi_1.7.12        
     [94] lattice_0.21-8         Matrix_1.6-1           vctrs_0.6.3           
     [97] pillar_1.9.0           lifecycle_1.0.3        spatstat.geom_3.0-3   
    [100] lmtest_0.9-40          RcppAnnoy_0.0.20       data.table_1.14.6     
    [103] cowplot_1.1.1          irlba_2.3.5.1          httpuv_1.6.9          
    [106] patchwork_1.1.2        R6_2.5.1               promises_1.2.0.1      
    [109] KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.36.0     
    [112] codetools_0.2-19       MASS_7.3-60            assertthat_0.2.1      
    [115] withr_2.5.0            sctransform_0.3.5      parallel_4.1.2        
    [118] hms_1.1.2              grid_4.1.2             rmarkdown_2.21        
    [121] googledrive_2.0.0      Rtsne_0.16             spatstat.explore_3.0-5
    [124] shiny_1.7.4            lubridate_1.9.0       
