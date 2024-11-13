#################
# GGPLOT THEME #
#################
txt_size = 7
my_theme <-
  list(
    #scale_fill_manual(values = friendly_cols),
    #scale_color_manual(values = friendly_cols),
    theme_classic() +
      #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
      theme(
        panel.border = element_rect(fill=NA, linewidth = .7),
        axis.line = element_line(),
        panel.grid.major = element_blank(), #element_line(linewidth = 0.2),
        panel.grid.minor = element_blank(), #element_line(linewidth = 0.1),
        text = element_text(size = txt_size),
        plot.title = element_text(hjust = 0.5),
        #legend.position = "bottom",
        #aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title = element_text(margin = margin(t = 12, r = 10, b = 10, l = 10))
        #axis.title.y = element_text(size = txt_size, margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )

##########################
# PLOT CLUSTERS FUNCTION #
##########################
# obj <- seuratObj
# cluster <- sym("RNA_snn_res.0.1")
plot_clusters.fun <- function(obj, cluster, 
                              red = "umapharmony", 
                              color = "Brew_all", 
                              lable = TRUE,
                              nudge = F,
                              lable_size = 4,
                              txt_size = 7,
                              dot_size = .5,
                              title = "colname",
                              seed = 127,
                              assay="RNA"){
  if(color[[1]] == "Brew_all" ){
    pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2") )}else{pal <- color}
  
  cluster <- sym(cluster)
  if(title == "colname"){title <- as_label(cluster)}
    else{title <- title}
  
  if(lable == TRUE){ lab <- cluster
  l=T
  t <- NoLegend() #+ labs(color= "Clusters")
    }
    else if(lable == FALSE){ 
    lab <- cluster
    l=F
    t <- NoLegend()
    text <- geom_blank() 
    #+ labs(color= "Clusters"){
      }
      else{lab <- sym(lable)
      l=T
      t <- guides(color = "none")}
  
  DefaultAssay(obj) <- assay
  
  feat <- obj %>%
    select(.cell, !!(cluster), !!(lab), nCount_RNA, nFeature_RNA) %>%
    group_by(!!(cluster)) %>%
    add_tally() %>%
    arrange(nFeature_RNA) %>%
    arrange(desc(n))

  # reduction method:
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  set.seed(seed)
  lable_df <- feat %>%
    ungroup() %>%
    group_by(!!(lab)) %>%
    select(!!(lab), contains(red)) %>% 
    summarize_all(mean) %>%
    mutate(vnudge = rnorm(nrow(.),-.01,0.5)) %>%
    mutate(hnudge = rnorm(nrow(.),-.01,0.5)) %>%
    {if(nudge) mutate(!!(red_1) := !!(red_1)+hnudge) else .} %>%
    {if(nudge) mutate(!!(red_2) := !!(red_2)+vnudge) else .}
    #mutate(!!(red_1) := !!(red_1)+hnudge) %>%
    #mutate(!!(red_2) := !!(red_2)+vnudge)
  
  lable_df<<-lable_df
  
  if(!(lable == FALSE)){text <- geom_text(data = lable_df, aes(label = !!(lab)),
                                          col="black", size=lable_size,  hjust=.7, vjust=-.3,) }
  else{text <- NULL}
  
  p <- ggplot(feat, aes(!!(red_1), !!(red_2), 
                        color = !!cluster), label=l) + 
    geom_point(alpha=.5, size=dot_size) + ggtitle(title) +
    #if(!(lable == FALSE)){geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5)} +
    #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
    text +
    scale_color_manual(values = pal, na.value = "transparent")  +
    my_theme + t +
    theme(
      title = element_blank(),
      plot.margin = unit(c(.05,.1,0,0), "cm"), #t,r,b,l c(.1,.1,-.4,-.1)
      axis.title.x = element_text(size=txt_size, vjust=3.5), # margin=margin(t=10, r=10, b=10, l=10,)
      axis.title.y = element_text(size=txt_size, vjust=-1),
      axis.text = element_text(size=txt_size),
      plot.title = element_text(size=txt_size, vjust=-9, hjust = 0.05))
  return(p)
}

#######################
# PLOT GENES FUNCTION #
#######################
# obj <- seuratObj
# gene <- sym("CTSK")
plot_genes.fun <- function(obj, 
                           gene, 
                           point_size = .5,
                           lable_size = 4,
                           mins=NULL, maxs=NULL, 
                           red="umapharmony", 
                           col=c("grey90","grey80","grey60","navy","black") ,
                           lable = TRUE){
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
    mutate(feat = round(.$feat*98)+1) %>%
    mutate(pal = c( col[1],colorRampPalette(col[-1])(99))[.$feat] ) %>%
    arrange(!!(gene))
  
  # reduction method:
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  # txt lable:
  if(lable == FALSE){l = FALSE
  text <- NoLegend()
  }else{l = TRUE
  if(lable != TRUE){obj <- mutate(obj, lab = pull(obj, lable))}
  
  lable_df <- obj %>%
    group_by(lab) %>%
    select(lab, contains(red)) %>% 
    summarize_all(mean) 
  
  text <- list(
    geom_point(data = lable_df, aes(x=!!(red_1), y=!!(red_2), ), col="#FFFFFF99", size=lable_size, alpha=.6),
    geom_text(data = lable_df, aes(label = lab), col="black", size=lable_size)) }
  
  p <- ggplot(obj, aes(!!(red_1), !!(red_2), label=l , color = pal) ) +
    geom_point(alpha = 0.5, size=point_size) + ggtitle(as_label(gene)) +
    text +
    scale_colour_identity() +
    theme_void() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) 
  return(p)
}

#################
# VOLCANO PLOT #
#################
revlog_trans <- function(base = exp(1)){
  ## Define the desired transformation.
  trans <- function(x){
    -log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(-x)
  }
  ## Creates the transformation
  scales::trans_new(paste("revlog-", base, sep = ""),
                    trans, ## The transformation function (can be defined using anonymous functions)
                    inv,  ## The reverse of the transformation
                    scales::log_breaks(base = base), ## default way to define the scale breaks
                    domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}


Volcano.fun_logFC <- function(DEGs_table, group, y.axis, 
                              up=c(1, 0.001), down = c(-1, 0.001),
                              lab_size = 4, dot_size = .3){
  tt <- DEGs_table %>% 
    mutate('p-value treshold' = ifelse(avg_log2FC >= 0 & p_val_adj <= 0.05 ,"Up", 
                                       ifelse(avg_log2FC <= -0 & p_val_adj  <= 0.05, "Down", 'NotSig'))) %>%
    mutate('Lable' = ifelse(avg_log2FC >= up[1] & p_val_adj <= up[2] | avg_log2FC <= down[1] & p_val_adj <= down[2],.$gene,NA)) %>%
    arrange(desc(p_val))
  
  if(y.axis == 'p-value') {
    plot <- ggplot(tt, aes(x = .data[[group]], y = avg_log2FC, colour = `p-value treshold`)) +
      #scale_x_continuous(trans = revlog_trans(), expand = c(0.005, 0.05)) +
      #expand_limits(x = c(0.001, 1)) +
      #geom_point(data = tt, alpha = 0.5, lable = tt$Lable) +
      geom_jitter(width = 0.3, alpha = 0.3, size=dot_size) +
      geom_text_repel(data = tt, label= tt$Lable, colour = "black", size=lab_size, #vjust = -0.6,
                      show.legend=FALSE, segment.color = NA,
                      box.padding = 0.1,       # Reduce padding around text
                      point.padding = 0.1 ) +  # Reduce space between text and points
      geom_hline(yintercept = 0, linetype = "solid", linewidth =.1) +
      # geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      guides(col = guide_legend(override.aes = list(size=1), keyheight = .6)) +
      theme_minimal() + 
      theme(axis.title = element_blank(),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 8, vjust = -1.5),
            axis.text.y = element_text(size = 7, hjust=1),
            axis.text.x = element_text(vjust = .5, size = 9)) +   # hjust=1, angle = 45, 
      scale_colour_manual(values = c("Up"= "red", "NotSig"= "grey90", "Down"="blue"),
                          name = paste0("FDR < ",up[2])) #+
    #facet_wrap(~cluster, nrow = 1)
    
  } else {
    pos <- which(abs(TopTable$FDR-0.05)==min(abs(TopTable$FDR-0.05)))
    plot <- ggplot(tt, aes(x = logFC, y = -log10(P.Value), colour = `p-value treshold`)) +
      geom_point(data = tt, alpha = 0.5, size=3) +
      geom_text(data = tt, label= tt$Lable, colour = "black", size=lab_size, vjust = -0.6,
                show.legend=FALSE, 
                #check_overlap = TRUE 
                #,point.padding = NA, segment.color = NA,
      ) + 
      # aes(label=taxa$Genus, color='Taxa')
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.001), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(TopTable$P.Value[pos]), linetype = "dashed", alpha = 0.5) +
      #geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      theme_linedraw() + theme_classic() + # Set the theme
      xlab(bquote('                                             '~log[2]~ "(fold change)")) + ylab(bquote('-'~log[10]~ "(p-value)")) + # Relabel the axes
      guides(col = guide_legend(override.aes = list(size=2), keyheight = .7)) +
      theme(#legend.position="none", 
        axis.title.x = element_text(hjust=0.001),
        legend.title = element_blank(),) + # Hide the legend
      scale_colour_manual(values = c("Up"= "red", "NotSig"= "black", "Down"="blue"),
                          breaks = c("Up", "Down"),
                          labels = c("Upregualated\nin DMPA\n", #"Non sig. diff.\nexpressed", 
                                     "Downregulated\nin DMPA")) +
      annotate(geom="text", x = Inf, y = -log10(TopTable$P.Value[pos]), size = 3, label = "FDR 0.05",
               color = "black", hjust = -.1) + 
      annotate(geom="text", x = Inf, y = -log10(0.001), size = 3, label = "p 0.001",
               color = "black", hjust = -.1) + 
      annotate(geom="text", x = Inf, y = -log10(0.05), size = 3, label = "p 0.05",
               color = "black", hjust = -.1) + coord_cartesian(clip = 'off') +
      annotate(geom="text", x = -3.5, y = 0, size = 8, label = "down in DMPA",
               color = "black", hjust = -.1, fontface =2) +
      annotate(geom="text", x = 1.5, y = 0, size = 8, label = "up in DMPA",
               color = "black", hjust = -.1, fontface =2)
  }
  return(plot)
}

###################################
# PLOT CLUSTERS PROP PIE FUNCTION #
###################################
# obj <- seuratObj
# cluster <- sym("RNA_snn_res.0.1")
plot_cell_pie_c.fun <- function(obj, markers, 
                              red = "umap_harmony", 
                              color = "Brew_all", 
                              radius_adj = 0,
                              #lable = TRUE, 
                              #lable_size = 4,
                              txt_size = 7,
                              dot_size = 0.5,
                              title = "colname",
                              assay="RNA"){
  if(color[[1]] == "Brew_all"){
    pal <- c(scales::hue_pal()(8),
             RColorBrewer::brewer.pal(9,"Set1"),
             RColorBrewer::brewer.pal(8,"Set2"),
             RColorBrewer::brewer.pal(8,"Accent"),
             RColorBrewer::brewer.pal(9,"Pastel1"),
             RColorBrewer::brewer.pal(8,"Pastel2") )}else{pal <- color}
  
  #cluster <- sym(cluster)
  
  
  DefaultAssay(obj) <- assay
  
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  feat <- obj %>%
    mutate(., FetchData(., vars = markers) ) %>%
    select(.cell, !!(red_1), !!(red_2), any_of(markers), nCount_RNA, nFeature_RNA) %>%
    as_tibble()
  
  radius = (max(feat[[red_1]]) - min(feat[[red_2]])) * (max(feat[[red_1]]) - min(feat[[red_2]]))
  radius = radius / nrow(feat)
  radius = radius / pi
  radius = sqrt(radius) * 0.85
  
  p <- ggplot(feat, aes(!!(red_1), !!(red_2))) + 
    #geom_point(alpha=.5, size=dot_size) + ggtitle(title) +
    geom_scatterpie(data=feat,aes(!!(red_1), !!(red_2),r=radius+radius_adj), cols=markers, color=NA) + 
    #if(!(lable == FALSE)){geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5)} +
    #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +

    scale_fill_manual(values = pal, na.value = "transparent")  +
    my_theme +
    theme(
      plot.margin = unit(c(.05,.1,0,0), "cm"), #t,r,b,l c(.1,.1,-.4,-.1)
      axis.title.x = element_text(size=txt_size, vjust=3.5), # margin=margin(t=10, r=10, b=10, l=10,)
      axis.title.y = element_text(size=txt_size, vjust=-1),
      axis.text = element_text(size=txt_size),
      plot.title = element_text(size=txt_size, vjust=-9, hjust = 0.05))
  return(p)
}
