Load Spatial data
================
12/13/22

``` r
source("../bin/render_with_jobs.R")
file_name <- "./00_plotting_st_data.md"
lab_dir <- "../lab_book/00_load_st_data/"

file <- paste0(basename(xfun::sans_ext(file_name)), '_', Sys.Date(), '.html')

# quarto
# render_html_with_job(out_dir = lab_dir)
# fs::file_move(path = file, new_path = paste0(lab_dir, file))

# currently using quarto for github and kniter for html du to source code option 
render_git_with_job()

# kniter
knit_html_with_job(out_dir = lab_dir)
```

### Load packages

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

#########
# PATHS #
#########
input_dir <- "../data/spatial_data"
result_dir <- "../results/00_load_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }
```

### Load ST data

``` r
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
    # … with 6,690 more rows

## add meta info

``` r
meta <- meta %>%
  select(orig.ident="ID", groups=Contraception) %>%
  mutate(groups = ifelse(.$groups=="no HC", "ctrl", .$groups))

DATA <-  DATA %>%
  rename(sp_annot2="sp_annot") %>%
  mutate(sp_annot = ifelse(grepl("epi_1|epi_2|epi_3", .$sp_annot2), "epi", .$sp_annot2 )) %>%
  mutate(sp_annot = ifelse(grepl("SubMuc_1|SubMuc_2|SubMuc_3", .$sp_annot2), "SubMuc", .$sp_annot2 )) %>%
  left_join(., meta, by="orig.ident") %>%
  select(groups, sp_annot, everything())
```

``` r
DATA %>%
  mutate(filt = ifelse(is.na(.$sp_annot),"filt","keep")) %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "filt",
          zoom = "zoom",
          annot_line = .1,
          img_alpha = 0,
          point_size = 0.5
        )
```

<img src="../Figures/plot_spots_to_remove.jpeg"
data-fig-align="center" />

``` r
# Filter spots outside manual annotation
DATA <- DATA[, !(is.na(DATA$sp_annot))]
DATA
```

    # A Seurat-tibble abstraction: 6,508 × 7
    # [90mFeatures=36601 | Cells=6508 | Active assay=RNA | Assays=RNA[0m
       .cell                 groups sp_annot orig.ident nCount_RNA nFeatur…¹ sp_an…²
       <chr>                 <chr>  <chr>    <chr>           <dbl>     <int> <chr>  
     1 P031_AAACGAGACGGTTGAT ctrl   SubMuc   P031              463       356 SubMuc 
     2 P031_AAACTGCTGGCTCCAA ctrl   SubMuc   P031             4081      2232 SubMuc 
     3 P031_AAAGTAGCATTGCTCA ctrl   epi      P031             7595      2588 epi    
     4 P031_AAAGTGTGATTTATCT ctrl   epi      P031            11394      3873 epi    
     5 P031_AAAGTTGACTCCCGTA ctrl   epi      P031             4617      2270 epi    
     6 P031_AAATACCTATAAGCAT ctrl   epi      P031             5538      2507 epi    
     7 P031_AAATCGTGTACCACAA ctrl   SubMuc   P031             1765      1084 SubMuc 
     8 P031_AAATTAACGGGTAGCT ctrl   SubMuc   P031              639       465 SubMuc 
     9 P031_AAATTTGCGGGTGTGG ctrl   SubMuc   P031             3570      1948 SubMuc 
    10 P031_AACCCTACTGTCAATA ctrl   SubMuc   P031             1438       940 SubMuc 
    # … with 6,498 more rows, and abbreviated variable names ¹​nFeature_RNA,
    #   ²​sp_annot2

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
    BLAS/LAPACK: /Users/vilkal/Applications/miniconda3/envs/Spatial_DMPA/lib/libopenblasp-r0.3.18.dylib

    locale:
    [1] sv_SE.UTF-8/sv_SE.UTF-8/sv_SE.UTF-8/C/sv_SE.UTF-8/sv_SE.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] sp_1.5-0           xml2_1.3.3         niceRplots_0.1.0   hdf5r_1.3.5       
     [5] Seurat_4.1.0       tidyseurat_0.5.3   SeuratObject_4.0.4 ttservice_0.1.2   
     [9] forcats_0.5.1      stringr_1.4.1      dplyr_1.0.7        purrr_0.3.4       
    [13] readr_2.1.2        tidyr_1.2.0        tibble_3.1.8       ggplot2_3.3.6     
    [17] tidyverse_1.3.1   

    loaded via a namespace (and not attached):
      [1] readxl_1.3.1          backports_1.4.1       plyr_1.8.7           
      [4] igraph_1.3.0          lazyeval_0.2.2        splines_4.1.2        
      [7] listenv_0.8.0         scattermore_0.8       digest_0.6.30        
     [10] htmltools_0.5.3       fansi_1.0.3           magrittr_2.0.3       
     [13] tensor_1.5            cluster_2.1.4         ROCR_1.0-11          
     [16] tzdb_0.2.0            globals_0.14.0        modelr_0.1.8         
     [19] matrixStats_0.61.0    vroom_1.5.7           spatstat.sparse_2.1-0
     [22] colorspace_2.0-3      rvest_1.0.2           ggrepel_0.9.1        
     [25] haven_2.4.3           xfun_0.33             crayon_1.5.2         
     [28] jsonlite_1.8.2        spatstat.data_2.1-4   survival_3.4-0       
     [31] zoo_1.8-9             glue_1.6.2            polyclip_1.10-0      
     [34] gtable_0.3.1          leiden_0.3.9          future.apply_1.8.1   
     [37] abind_1.4-5           scales_1.2.1          DBI_1.1.2            
     [40] spatstat.random_2.2-0 miniUI_0.1.1.1        Rcpp_1.0.9           
     [43] viridisLite_0.4.1     xtable_1.8-4          reticulate_1.24      
     [46] spatstat.core_2.4-2   bit_4.0.4             htmlwidgets_1.5.4    
     [49] httr_1.4.4            RColorBrewer_1.1-3    ellipsis_0.3.2       
     [52] ica_1.0-2             pkgconfig_2.0.3       farver_2.1.1         
     [55] uwot_0.1.11           dbplyr_2.1.1          deldir_1.0-6         
     [58] utf8_1.2.2            labeling_0.4.2        tidyselect_1.2.0     
     [61] rlang_1.0.6           reshape2_1.4.4        later_1.3.0          
     [64] munsell_0.5.0         cellranger_1.1.0      tools_4.1.2          
     [67] cli_3.4.1             generics_0.1.3        broom_0.7.12         
     [70] ggridges_0.5.3        evaluate_0.18         fastmap_1.1.0        
     [73] yaml_2.3.5            goftest_1.2-3         knitr_1.40           
     [76] bit64_4.0.5           fs_1.5.2              fitdistrplus_1.1-8   
     [79] RANN_2.6.1            pbapply_1.5-0         future_1.24.0        
     [82] nlme_3.1-160          mime_0.12             compiler_4.1.2       
     [85] rstudioapi_0.13       plotly_4.10.0         png_0.1-7            
     [88] spatstat.utils_2.3-0  reprex_2.0.1          stringi_1.7.8        
     [91] lattice_0.20-45       Matrix_1.5-3          vctrs_0.4.2          
     [94] pillar_1.8.1          lifecycle_1.0.3       spatstat.geom_2.4-0  
     [97] lmtest_0.9-40         RcppAnnoy_0.0.19      data.table_1.14.2    
    [100] cowplot_1.1.1         irlba_2.3.5           httpuv_1.6.5         
    [103] patchwork_1.1.1       R6_2.5.1              promises_1.2.0.1     
    [106] KernSmooth_2.23-20    gridExtra_2.3         parallelly_1.31.0    
    [109] codetools_0.2-18      MASS_7.3-58.1         assertthat_0.2.1     
    [112] withr_2.5.0           sctransform_0.3.3     mgcv_1.8-40          
    [115] parallel_4.1.2        hms_1.1.1             grid_4.1.2           
    [118] rpart_4.1.16          rmarkdown_2.18        Rtsne_0.15           
    [121] shiny_1.7.1           lubridate_1.8.0      
