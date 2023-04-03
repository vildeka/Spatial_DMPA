#################
# GGPLOT THEME #
#################
txt_size = 7
my_theme <-
  list(
    #scale_fill_manual(values = friendly_cols),
    #scale_color_manual(values = friendly_cols),
    theme_bw() +
      #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1),
        text = element_text(size = txt_size),
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
  
  lable_df <- feat %>%
    ungroup() %>%
    group_by(!!(lab)) %>%
    select(!!(lab), contains(red)) %>% 
    summarize_all(mean)
  
  red_1 <- sym(paste0(red, "_1"))
  red_2 <- sym(paste0(red, "_2"))
  
  if(!(lable == FALSE)){text <- geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5) }else{text <- NULL}
  
  p <- ggplot(feat, aes(!!(red_1), !!(red_2), 
                        color = !!cluster), label=l) + 
    geom_point(alpha=.5, size=dot_size) + ggtitle(title) +
    #if(!(lable == FALSE)){geom_text(data = lable_df, aes(label = !!(lab)), col="black", size=2.5)} +
    #guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
    text +
    scale_color_manual(values = pal)  +
    my_theme + t +
    theme(
      plot.margin = unit(c(.05,.1,0,0), "cm"), #t,r,b,l c(.1,.1,-.4,-.1)
      axis.title.x = element_text(size=txt_size, vjust=3.5), # margin=margin(t=10, r=10, b=10, l=10,)
      axis.title.y = element_text(size=txt_size, vjust=-1),
      axis.text = element_text(size=txt_size),
      plot.title = element_text(size=txt_size, vjust=-5, hjust = 0.05))
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
                           mins=NULL, maxs=NULL, 
                           red="umap_harmony", 
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
    geom_point(alpha = 0.5, size=point_size) + ggtitle(as_label(gene)) +
    text + #scale_color_viridis(option = "D", na.value="#EBECF0") +
    scale_colour_identity() +
    my_theme + theme_void() +
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


Volcano.fun_logFC <- function(DEGs_table, group, y.axis, up=c(1, 0.001), down = c(-1, 0.001)){
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
      geom_jitter(width = 0.3, alpha = 0.3, size=.1) +
      geom_text_repel(data = tt, label= tt$Lable, colour = "black", size=2, #vjust = -0.6,
                      show.legend=FALSE, segment.color = NA,
                      #check_overlap = TRUE 
                      #,point.padding = NA, segment.color = NA,
      )+
      geom_hline(yintercept = 0, linetype = "solid") +
      #geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
      scale_colour_manual(values = c("Up"= "red", "NotSig"= "grey90", "Down"="blue"),
                          name = paste0("FDR < ",up[2])) #+
    #facet_wrap(~cluster, nrow = 1)
    
  } else {
    pos <- which(abs(TopTable$FDR-0.05)==min(abs(TopTable$FDR-0.05)))
    plot <- ggplot(tt, aes(x = logFC, y = -log10(P.Value), colour = `p-value treshold`)) +
      geom_point(data = tt, alpha = 0.5, size=3) +
      geom_text(data = tt, label= tt$Lable, colour = "black", size=4, vjust = -0.6,
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
      theme(#legend.position="none", 
        axis.title.x = element_text(hjust=0.001),
        legend.title = element_blank()) + # Hide the legend
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
