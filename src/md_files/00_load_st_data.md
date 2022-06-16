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
scale.factors <- map(DATA@images, ~pluck(.x, "scale.factors")$hires)
img_coord <- map(DATA@images, ~pluck(.x, "coordinates"))
  
img_coord <- img_coord %>%
  map2(., scale.factors, 
       ~mutate(.x, imagecol = .$imagecol * .y,
                   imagerow = .$imagerow * .y) )
  
a2 <- map(names(img_coord), ~xml2::as_list( read_xml( paste0( input_dir,"/",.x,"/",.x,".svg")) ) )
dd <- map(a2, ~as.numeric( strsplit( attributes(.x$svg)$viewBox , " ")[[1]] ))
dd2 <- map(sample_id, ~dim(DATA@images[[.x]]@image))

# red <- "P080"
# sample_id <- "P080"
# a2 <- xml2::as_list( read_xml( paste0( input_dir,"/",red,"/",red,".svg")) )
# dd <- as.numeric( strsplit( attributes(a2$svg)$viewBox , " ")[[1]] )
# dd2 <- dim(DATA@images[[red]]@image)
# img_coord <- img_coord[[2]][1:5]
get_sp_annot <- function(a2, dd, dd2, img_coord, sample_id){
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
    ungroup() 
  
  img_coord <- temp %>%
    list_modify("fov" = NULL, full_image = NULL) %>%
    compact() %>%
    imap(., ~mutate(img_coord, !!.y := sp::point.in.polygon(
         point.x = img_coord$imagecol,
         point.y = img_coord$imagerow,
         pol.x = .x$x,
         pol.y = .x$y )) ) %>%
    map(., ~rownames_to_column(., var = "barcodes")) %>%
    Reduce(dplyr::full_join, .)
  
  sp_annot <- img_coord %>%
    mutate(across(7:ncol(.), ~ifelse(. == 0, NA, .)) ) %>%
    pivot_longer(., cols = 7:ncol(.), names_to ="sp_annot", values_to = "count") %>%
    filter(!(is.na(count))) %>%
    group_by(barcodes) %>%
    mutate(dupp = row_number()) %>%
    ungroup() %>%
    filter(., .$dupp == 1) %>%
    select(., barcodes, sp_annot)
  
  return(list(coord=annot_coord, annot=sp_annot))
}

annot <- list(a2, dd, dd2, img_coord, names(img_coord)) %>%
  pmap(., ~get_sp_annot(..1, ..2, ..3, ..4, ..5)) %>%
  set_names(., names(img_coord))

DATA@tools <- map(annot, 1)

sp <- map(annot, 2) %>% bind_rows()
DATA <- left_join( DATA, sp, by=c(".cell"="barcodes")) 

DATA 
```

    ## # A Seurat-tibble abstraction: 6,700 × 5
    ## # [90mFeatures=36601 | Cells=6700 | Active assay=RNA | Assays=RNA[0m
    ##    .cell                 orig.ident nCount_RNA nFeature_RNA sp_annot
    ##    <chr>                 <chr>           <dbl>        <int> <chr>   
    ##  1 P031_AAACGAGACGGTTGAT P031              463          356 SubMuc  
    ##  2 P031_AAACTGCTGGCTCCAA P031             4081         2232 SubMuc  
    ##  3 P031_AAAGTAGCATTGCTCA P031             7595         2588 epi     
    ##  4 P031_AAAGTGTGATTTATCT P031            11394         3873 epi     
    ##  5 P031_AAAGTTGACTCCCGTA P031             4617         2270 epi     
    ##  6 P031_AAATACCTATAAGCAT P031             5538         2507 epi     
    ##  7 P031_AAATCGTGTACCACAA P031             1765         1084 SubMuc  
    ##  8 P031_AAATGGCCCGTGCCCT P031              789          477 <NA>    
    ##  9 P031_AAATTAACGGGTAGCT P031              639          465 SubMuc  
    ## 10 P031_AAATTTGCGGGTGTGG P031             3570         1948 SubMuc  
    ## # … with 6,690 more rows

## Save seurat object

``` r
DATA <- DATA %>%
  select(sp_annot, everything())
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

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /Users/vilkal/Applications/miniconda3/envs/Spatial_DMPA/lib/libopenblasp-r0.3.18.dylib
    ## 
    ## locale:
    ## [1] sv_SE.UTF-8/sv_SE.UTF-8/sv_SE.UTF-8/C/sv_SE.UTF-8/sv_SE.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] sp_1.5-0           xml2_1.3.3         niceRplots_0.1.0   hdf5r_1.3.5       
    ##  [5] Seurat_4.1.0       tidyseurat_0.5.3   SeuratObject_4.0.4 ttservice_0.1.2   
    ##  [9] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.8        purrr_0.3.4       
    ## [13] readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
    ## [17] tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1          backports_1.4.1       plyr_1.8.7           
    ##   [4] igraph_1.3.0          lazyeval_0.2.2        splines_4.1.2        
    ##   [7] listenv_0.8.0         scattermore_0.8       digest_0.6.29        
    ##  [10] htmltools_0.5.2       fansi_1.0.3           magrittr_2.0.3       
    ##  [13] tensor_1.5            cluster_2.1.2         ROCR_1.0-11          
    ##  [16] tzdb_0.2.0            globals_0.14.0        modelr_0.1.8         
    ##  [19] matrixStats_0.61.0    spatstat.sparse_2.1-0 colorspace_2.0-3     
    ##  [22] rvest_1.0.2           ggrepel_0.9.1         haven_2.4.3          
    ##  [25] xfun_0.30             crayon_1.5.1          jsonlite_1.8.0       
    ##  [28] spatstat.data_2.1-4   survival_3.2-13       zoo_1.8-9            
    ##  [31] glue_1.6.2            polyclip_1.10-0       gtable_0.3.0         
    ##  [34] leiden_0.3.9          future.apply_1.8.1    abind_1.4-5          
    ##  [37] scales_1.2.0          DBI_1.1.2             spatstat.random_2.2-0
    ##  [40] miniUI_0.1.1.1        Rcpp_1.0.8.3          viridisLite_0.4.0    
    ##  [43] xtable_1.8-4          reticulate_1.24       spatstat.core_2.4-2  
    ##  [46] bit_4.0.4             htmlwidgets_1.5.4     httr_1.4.2           
    ##  [49] RColorBrewer_1.1-3    ellipsis_0.3.2        ica_1.0-2            
    ##  [52] farver_2.1.0          pkgconfig_2.0.3       sass_0.4.1           
    ##  [55] uwot_0.1.11           dbplyr_2.1.1          deldir_1.0-6         
    ##  [58] utf8_1.2.2            tidyselect_1.1.2      rlang_1.0.2          
    ##  [61] reshape2_1.4.4        later_1.3.0           munsell_0.5.0        
    ##  [64] cellranger_1.1.0      tools_4.1.2           cli_3.3.0            
    ##  [67] generics_0.1.2        broom_0.7.12          ggridges_0.5.3       
    ##  [70] evaluate_0.15         fastmap_1.1.0         yaml_2.3.5           
    ##  [73] goftest_1.2-3         knitr_1.38            bit64_4.0.5          
    ##  [76] fs_1.5.2              fitdistrplus_1.1-8    RANN_2.6.1           
    ##  [79] pbapply_1.5-0         future_1.24.0         nlme_3.1-157         
    ##  [82] mime_0.12             compiler_4.1.2        rstudioapi_0.13      
    ##  [85] plotly_4.10.0         png_0.1-7             spatstat.utils_2.3-0 
    ##  [88] reprex_2.0.1          bslib_0.3.1           stringi_1.7.6        
    ##  [91] lattice_0.20-45       Matrix_1.4-1          vctrs_0.4.1          
    ##  [94] pillar_1.7.0          lifecycle_1.0.1       spatstat.geom_2.4-0  
    ##  [97] lmtest_0.9-40         jquerylib_0.1.4       RcppAnnoy_0.0.19     
    ## [100] data.table_1.14.2     cowplot_1.1.1         irlba_2.3.5          
    ## [103] httpuv_1.6.5          patchwork_1.1.1       R6_2.5.1             
    ## [106] promises_1.2.0.1      KernSmooth_2.23-20    gridExtra_2.3        
    ## [109] parallelly_1.31.0     codetools_0.2-18      MASS_7.3-57          
    ## [112] assertthat_0.2.1      withr_2.5.0           sctransform_0.3.3    
    ## [115] mgcv_1.8-40           parallel_4.1.2        hms_1.1.1            
    ## [118] grid_4.1.2            rpart_4.1.16          rmarkdown_2.11       
    ## [121] Rtsne_0.15            shiny_1.7.1           lubridate_1.8.0
