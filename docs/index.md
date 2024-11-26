## Workflow

### Analysis scripts

1.  [00_load_st_data.Rmd](https://vildeka.github.io/Spatial_DMPA/00_load_st_data){.uri}
2.  [01_QC_st_data](https://vildeka.github.io/Spatial_DMPA/01_QC_st_data){.uri}
3.  

### Manuscript figures

``` r
ord1 <- c("Sup_1","Sup_2","Basal_2","Basal_1","0","1","2","3","4","8","10")
ord2 <- c("6","9","7","5","0","1","2","3","4","8","10")

DEGs_table <- DEGs_table %>%
  mutate(subgroup = factor(.$subgroup, levels = ord2)) %>%
  mutate(layers = factor(.$layers, levels = ord1)) %>%
  dplyr::rename(Clusters="layers") 
```