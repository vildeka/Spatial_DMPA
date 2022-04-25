### Load data and libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(readxl)

#########
# PATHS #
#########
cell_types <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-20358-y/MediaObjects/41467_2020_20358_MOESM5_ESM.xlsx"
ref_info <- "../../data/Reference_datasets.xlsx"
data_dir <- "../../data/Reference_data"
result_dir <- "../../results/00_load_ref_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
reference_info <- read_xlsx(ref_info, na = "NA", sheet = "selection")

ref <- reference_info %>%
  separate_rows(., Files, sep=", ") %>%
  mutate(G = str_extract(.$Files, "^G\\w\\w\\d+"), .after="accesion_number") %>%
  mutate(type = ifelse(grepl("^GSE", .$Files), "series", "samples"), .after="G") %>%
  select(type, G, Files, accesion_number, dataset)
```

### Download data from GEO

``` r
#########################
# DOWNLOAD SUPPL. DATA #
#########################
download.file(cell_types, destfile=paste0(
  data_dir, "/GSE151202/","Cell_types.xlsx"), method = "auto")

##########################
# DOWNLOAD DATA FROM GEO #
##########################
ftp <- "https://ftp.ncbi.nlm.nih.gov/geo"


# Create folder for each accession
map(reference_info$accesion_number, 
    ~if( isFALSE(dir.exists(paste0(data_dir,"/",.x))) ) {
      dir.create(paste0(data_dir,"/",.x),recursive = TRUE) })

# download files from GEO
pwalk(ref, 
      ~download.file(paste(ftp, ..1, str_replace(..2, "\\d{3}$", "nnn"),
                           ..2, "suppl", ..3, sep="/"), method="auto",
                     destfile=paste0(data_dir,"/",..4,"/",..3)))

# Untar and unzip files
tar <- list.files(data_dir, 
                  pattern = ".+\\.tar.gz$|.+\\.tar$", recursive = T, full.names = T)
map(tar, ~untar(.x, exdir = dirname(.x)))
zip <- list.files(data_dir, 
                  pattern = ".+\\.gz$", recursive = T, full.names = T)
#map(zip, ~unzip(.x, exdir = dirname(.x), list=T))
system(paste0("gunzip ",data_dir,"/*/*.gz"))
```

``` r
######################
# LOAD SPARSE MATRIX #
######################
f <- c("matrix.mtx$", "genes.tsv$|features.tsv$", "barcodes.tsv$") %>% set_names()

f_list <- map_df(f, ~list.files(path = data_dir, pattern = .x,
                      full.names = T, recursive = T) ) %>%
  mutate(accesion_number = str_match(.$"matrix.mtx$", paste0(data_dir,"\\/([^\\/]+)\\/"))[,2]) %>%
  mutate(sample_name = str_match(.$"matrix.mtx$", ".*\\/([^\\/]+).matrix\\.mtx$")[,2]) %>%
  mutate(sample_name = ifelse(.$accesion_number == "GSE173231", 
                              str_extract(.$sample_name, "(?<=_).+(?=_)"), .$sample_name)) %>%
  group_by(accesion_number) %>%
  mutate(sample_id = paste0(accesion_number, "_", row_number())) %>%
  ungroup()

matrix_list <- pmap(f_list, ~ReadMtx(mtx = ..1, features = ..2, cells = ..3 ) ) %>%
                     set_names(., f_list$sample_name)

#### other formats: ####
# rds <- readRDS("../data/Reference_data/GSE111976/GSE111976_ct_endo_10x.rds")
# tbl <- read.table(file = "../data/Reference_data/GSE134355/GSM3980130_Adult-Cervix1_dge.txt",
#                   header = TRUE, row.names = 1)
# 
# # matrix_list <- append(expression_matrix, list("GSE111976"=rds, atlas=tbl))
# matrix_list <- append(expression_matrix, list("GSE134355_1"=tbl))
```

### Create seurat object

``` r
###################
# SAMPLE METADATA #
###################
meta_list <- f_list %>%
  select(sample_name, accesion_number) %>%
  #add_row(sample_name="atlas", accesion_number="GSE134355", name="GSE134355_1") %>%
  left_join(., select(ref, dataset, accesion_number), by="accesion_number") %>% 
  unique() %>%
  mutate(p = str_replace(str_match(.$dataset, "(.+) (et\\.\\sal\\.)")[,2], "\\s", "_")) %>%
  {. ->> meta } %>% 
  split(., .$sample_name) #%>%
  #set_names(., sort(unique(meta$sample_name)))

#######################
# CREATE SEURAT OJECT #
#######################
seuratObj_list <- map(matrix_list, ~CreateSeuratObject(counts = .x )) %>%
  imap(., ~AddMetaData(object = .x, 
                       metadata = cbind(.@meta.data, 
                                        slice(meta_list[[.y]], rep(1:n(), each=nrow(.@meta.data))))) ) %>%
  map(., ~SetIdent(., value = .@meta.data$sample_name)) %>%
  imap(., ~RenameCells(.x, 
                      new.names = paste0(str_extract(colnames(.x[["RNA"]]), "^[ATCG]+"),"_",.y)) )

# Endometrium data:
# meta <- read.csv("../data/Reference_data/GSE111976/GSE111976_summary_10x_day_donor_ctype.csv", row.names = 1)
# seuratObj_list[["endo"]] <- AddMetaData(seuratObj_list[["endo"]], metadata = meta)
# AddMetaData(seuratObj_list[["GSE134355_1"]], metadata = meta)
```

### Add cell annotation

``` r
########################
# ADD CELL ANNOTATION  #
########################
path <- paste0(data_dir, "/GSE151202/","Cell_types.xlsx")

cell_types <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_xlsx, path = path) %>%
  bind_rows() %>%
  filter(grepl("^N.+", .$sample)) %>%
  mutate(barcode = paste0(str_extract(.$barcode, "^[ATCG]+"),"_",.$sample)) %>%
  select(.cell="barcode", cell_id="cell types")
  

seuratObj_list <- seuratObj_list %>%
  map(., ~left_join(.x, cell_types, by = ".cell")) %>%
  map(., ~tidyseurat::select(., .cell, cell_id, sample_name, everything(), -orig.ident))
```

### Merge datasets into one seurat object

``` r
merged_Seurat <- merge(seuratObj_list[[1]], y = seuratObj_list[2:length(seuratObj_list)])
```

### Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(merged_Seurat, paste0(result_dir,"seurat_object_merged.RDS"))
# seuratObj_list <- readRDS(paste0(result_dir,"seuratObj_list_not_modified.RDS"))
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
    ##  [1] readxl_1.3.1       tidyseurat_0.5.1   ttservice_0.1.2    SeuratObject_4.0.4
    ##  [5] Seurat_4.1.0       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.8       
    ##  [9] purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.6      
    ## [13] ggplot2_3.3.5      tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15            colorspace_2.0-3      deldir_1.0-6         
    ##   [4] ellipsis_0.3.2        ggridges_0.5.3        fs_1.5.2             
    ##   [7] spatstat.data_2.1-4   rstudioapi_0.13       leiden_0.3.9         
    ##  [10] listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.3          
    ##  [13] lubridate_1.8.0       xml2_1.3.3            codetools_0.2-18     
    ##  [16] splines_4.1.2         knitr_1.38            polyclip_1.10-0      
    ##  [19] jsonlite_1.8.0        broom_0.7.12          ica_1.0-2            
    ##  [22] cluster_2.1.2         dbplyr_2.1.1          png_0.1-7            
    ##  [25] uwot_0.1.11           spatstat.sparse_2.1-0 sctransform_0.3.3    
    ##  [28] shiny_1.7.1           compiler_4.1.2        httr_1.4.2           
    ##  [31] backports_1.4.1       lazyeval_0.2.2        assertthat_0.2.1     
    ##  [34] Matrix_1.4-1          fastmap_1.1.0         cli_3.2.0            
    ##  [37] later_1.3.0           htmltools_0.5.2       tools_4.1.2          
    ##  [40] igraph_1.3.0          gtable_0.3.0          glue_1.6.2           
    ##  [43] reshape2_1.4.4        RANN_2.6.1            Rcpp_1.0.8.3         
    ##  [46] scattermore_0.8       cellranger_1.1.0      jquerylib_0.1.4      
    ##  [49] vctrs_0.4.0           nlme_3.1-157          lmtest_0.9-40        
    ##  [52] spatstat.random_2.2-0 xfun_0.30             globals_0.14.0       
    ##  [55] rvest_1.0.2           mime_0.12             miniUI_0.1.1.1       
    ##  [58] lifecycle_1.0.1       irlba_2.3.5           goftest_1.2-3        
    ##  [61] future_1.24.0         MASS_7.3-56           zoo_1.8-9            
    ##  [64] scales_1.1.1          spatstat.core_2.4-2   spatstat.utils_2.3-0 
    ##  [67] hms_1.1.1             promises_1.2.0.1      parallel_4.1.2       
    ##  [70] RColorBrewer_1.1-3    yaml_2.3.5            gridExtra_2.3        
    ##  [73] reticulate_1.24       pbapply_1.5-0         sass_0.4.1           
    ##  [76] rpart_4.1.16          stringi_1.7.6         rlang_1.0.2          
    ##  [79] pkgconfig_2.0.3       matrixStats_0.61.0    evaluate_0.15        
    ##  [82] lattice_0.20-45       tensor_1.5            ROCR_1.0-11          
    ##  [85] htmlwidgets_1.5.4     patchwork_1.1.1       cowplot_1.1.1        
    ##  [88] tidyselect_1.1.2      parallelly_1.31.0     RcppAnnoy_0.0.19     
    ##  [91] plyr_1.8.7            magrittr_2.0.3        R6_2.5.1             
    ##  [94] generics_0.1.2        DBI_1.1.2             mgcv_1.8-40          
    ##  [97] pillar_1.7.0          haven_2.4.3           withr_2.5.0          
    ## [100] fitdistrplus_1.1-8    abind_1.4-5           survival_3.2-13      
    ## [103] future.apply_1.8.1    modelr_0.1.8          crayon_1.5.1         
    ## [106] KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_2.4-0  
    ## [109] plotly_4.10.0         tzdb_0.2.0            rmarkdown_2.11       
    ## [112] grid_4.1.2            data.table_1.14.2     reprex_2.0.1         
    ## [115] digest_0.6.29         xtable_1.8-4          httpuv_1.6.5         
    ## [118] munsell_0.5.0         viridisLite_0.4.0     bslib_0.3.1
