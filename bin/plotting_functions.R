#################
# GGPLOT THEME #
#################
my_theme <-
  list(
    #scale_fill_manual(values = friendly_cols),
    #scale_color_manual(values = friendly_cols),
    theme_bw() +
      #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        #legend.position = "bottom",
        #aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
        #axis.title.y = element_text(size = txt_size, margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )

##########################
# PLOT CLUSTERS FUNCTION #
##########################
# obj <- seuratObj
# cluster <- sym("RNA_snn_res.0.1")
plot_clusters.fun <- function(obj, cluster, 
                              red = "umap_harmony", 
                              color = "Brew_all", 
                              lable = TRUE, 
                              txt_size = 5,
                              dot_size = 0.5,
                              assay="RNA"){
  if(color[[1]] == "Brew_all"){
    pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2") )}else{pal=color}
  
  cluster <- sym(cluster)
  
  if(lable == TRUE){ lab <- cluster
  l=T
  t <- NoLegend() #+ labs(color= "Clusters")
  }else if(lable == FALSE){ lab <- cluster
  l=F
  t <- NoLegend()
  text <- geom_blank() #+ labs(color= "Clusters"){
  }else{lab <- sym(lable)
  l=T
  t <- guides(color = "none")}
  
  DefaultAssay(obj) <- assay
  
  feat <- obj %>%
    select(.cell, !!(cluster), !!(lab), nCount_RNA, nFeature_RNA) %>%
    group_by(!!(cluster)) %>%
    add_tally() %>%
    arrange(nFeature_RNA) %>%
    arrange(desc(n))
  
  lable_df <- feat %>%
    ungroup() %>%
    group_by(!!(lab)) %>%
    select(!!(lab), contains(red)) %>% 
    summarize_all(mean)
  
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  if(!(lable== FALSE)){text <- geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5) } 
  
  p <- ggplot(feat, aes(!!(red_1), !!(red_2), 
                        color = !!cluster), label=l) + 
    geom_point(alpha=.5, size=dot_size) + ggtitle(as_label(cluster), ) +
    #if(!(lable == FALSE)){geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5)} +
    #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
    text +
    scale_color_manual(values = pal)  +
    my_theme + t +
    theme(axis.title = element_text(size=txt_size, margin=margin(t=10, r=10, b=10, l=10)),
          axis.text=element_text(size=txt_size),
          plot.title = element_text(size=txt_size))
  return(p)
}

#######################
# PLOT GENES FUNCTION #
#######################
# obj <- seuratObj
# gene <- sym("CTSK")
plot_genes.fun <- function(obj, gene, mins=NULL, maxs=NULL, red = "umap_harmony", lable = TRUE){
  gene <- sym(gene)
  obj <- obj %>%
    mutate(lab = obj@active.ident) %>%
    mutate(., FetchData(., vars = c(as_label(gene))) ) %>%
    mutate(feat = !!(gene))
  feat_vec <- pull(obj, as_label(gene))
  
  # Colour pal:
  if(is.null(mins)){
    mins <- min(c(feat_vec, 0),na.rm = T)} # get 0 or negative value
  if(is.null(maxs)){maxs <- quantile(feat_vec,0.99,na.rm = T) # get percentile
  if(maxs==0){maxs <- max(feat_vec,na.rm = T)}
  }
  if(max(feat_vec, na.rm=T) != 0){
    # Calculate percentage:
    obj <- obj %>%
      mutate(feat = (!!(gene) - mins) / ( maxs - mins) ) %>%
      mutate(feat = ifelse(.$feat > 1, 1, .$feat))
  }
  
  obj <- obj %>%
    #select(1:3, !!(gene)) %>%
    mutate(feat = round(.$feat*98)+1) %>%
    mutate(pal = c( col[1],colorRampPalette(col[-1])(99))[.$feat] ) %>%
    arrange(!!(gene))
  
  # reduction method:
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  # txt lable:
  if(lable == FALSE){l = FALSE
  text <- NoLegend() #+ labs(color= "Clusters")
  }else{l = TRUE
  if(lable != TRUE){obj <- mutate(obj, lab = pull(obj, lable))}
  
  lable_df <- obj %>%
    group_by(lab) %>%
    select(lab, contains(red)) %>% 
    summarize_all(mean) 
  
  text <- geom_text(data = lable_df, aes(label = lab), col="black", size=2.5) }
  
  p <- ggplot(obj, aes(!!(red_1), !!(red_2), label=l , color = pal) ) +
    geom_point(alpha = 0.5, size=.5) + ggtitle(as_label(gene)) +
    text + #scale_color_viridis(option = "D", na.value="#EBECF0") +
    scale_colour_identity() +
    my_theme + theme_void() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) 
  return(p)
}
