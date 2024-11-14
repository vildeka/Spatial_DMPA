library(RColorBrewer)
library(scatterpie)
#################
# GEOM SPATIAL #
################
geom_spatial <-  function(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE,
    ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

########################
# PLOT RAW DATA VALUES #
########################
# can also plot single images and meta data
# spe <- DATA_st
# geneid <- sym("CDH1")
# zoom <- sym("zoom")
# assay="SCDC"
# img_alpha = 0
# sampleid = c("P105")

plot_spatial.fun <- function(
    spe,
    assay="RNA",
    sp_annot = TRUE,
    sampleid = c("P080"),
    geneid = "CD3E",
    title = " ",
    lab = TRUE,
    image_id = "hires",
    alpha = 1,
    ncol = 2,
    save_space = T,
    spectral = TRUE,
    annot_line = .3,
    colors = NULL, # lightgray
    point_size = 1.75,
    img_alpha = 0,
    zoom = NULL ) {
  
  # select samples to plot:
  sampleid <- set_names(sampleid)
  spe <- spe %>% filter(., orig.ident %in% sampleid)
  
  # Set default assay
  DefaultAssay(spe) <- assay
  
  ## get feature to plot:
  if (!(c(geneid) %in% colnames(spe@meta.data))) {
    spe <- spe %>%
      mutate(., FetchData(., vars = c(geneid)) ) 
  }
  
  ## Colour pallets:
  if (is.numeric(pull(spe, geneid))){
    if (is.null(colors)){
      #cont_colors <- c("grey90", "mistyrose", "red", "dark red", "black")
      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
      cont_colors <- myPalette(100)
    }else{cont_colors <- colors}
    colour_pallet <- scale_color_gradientn(geneid, colours = cont_colors,
                                           na.value = "#FFFFFF",
                                           guide = guide_colourbar(barwidth = .5, barheight = 5 ))
    guides <- NULL 
  }else{
    if (is.null(colors)){
      # scales::show_col(disc_colors)
      disc_colors <- c(RColorBrewer::brewer.pal(9,"Pastel1"),
                       RColorBrewer::brewer.pal(9,"Set1"),
                       scales::hue_pal()(8),
                       RColorBrewer::brewer.pal(8,"Set2"),
                       RColorBrewer::brewer.pal(8,"Accent"),
                       
                       RColorBrewer::brewer.pal(8,"Pastel2") )
    }else{disc_colors <- colors}
    colour_pallet <- scale_fill_manual(geneid, values = disc_colors, aesthetics = c("colour"), na.value = "grey90")
    guides <- guides(colour = guide_legend(override.aes = list(size=2), keyheight = .5))
  }
  
  # get scale factor:
  scale_fact <- map_dbl(sampleid, ~pluck(spe@images, .x, "scale.factors", image_id)) #spe@images[[sampleid]]@scale.factors[[image_id]]
  geneid <- enquo(geneid)
  
  # get spot coordinates:
  df <- map(sampleid, ~pluck(spe@images, .x, "coordinates")) %>%
    map2(., scale_fact, ~mutate(.x, scale_fact = .y)) %>%
    bind_rows(., .id = "orig.ident") %>%
    mutate(imagecol = .$imagecol * scale_fact) %>%
    mutate(imagerow = .$imagerow * scale_fact) %>%
    {. ->> tools} %>%
    rownames_to_column(var = "barcode") %>%
    left_join(.,select(spe, barcode=".cell",!!(geneid)), by="barcode") %>%
    as_tibble() 
  
  if(lab){
    gr <- spe@meta.data %>% group_by(groups, orig.ident ) %>% nest() %>% arrange(match(orig.ident, sampleid)) %>% pull(., "groups")
    text_annot <- tibble(sample_id = sampleid, x=500, y=800, orig.ident = sampleid, gr = gr) 
    txt <- list(geom_text(aes(label = sample_id, x=x, y=y), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt), # sample ID
                geom_text(aes(label = gr, x=x, y=y+90), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt))
  }else{txt <- NULL}
  
  # select viewframe:
  if(length(spe@tools)!=0){
    if(is.null(zoom)){zoom <- "zoom"}
    tools <- sampleid %>%
      map(., ~pluck(spe@tools, .x)) %>% 
      bind_rows(.id = "ID") %>% 
      filter(.data[["name"]] == zoom) }
  
  l <- tools %>% 
    mutate("_row"=y, "_col"=x) %>%
    summarise(across("_row":"_col", 
                     list(min=min, max=max), 
                     .names = "{.fn}{.col}")) 
  lims <- list(xlim(l$min_col,l$max_col), ylim(l$max_row,l$min_row))
  
  width <- l$max_col-l$min_col
  height <- l$max_row-l$min_row
  aspect <- width/height
  
  ## Spatial image:
  if (!(img_alpha == 0)){
    
    im <- map(sampleid, ~pluck(spe@images, .x, "image")) 
    # im is a 3D matrix, each dimension representing red, green or blue
    # a 4th dimension can be added to specify the opacity

    # im <- map(im, ~ matrix(
    #   rgb(.x[,,1],.x[,,2],.x[,,3], .x[4,,]* img_alpha), nrow=dim(.x)[1]))
    # 
    # img <- map(im, ~as.raster(.x))
    
    # add alpha dimension to the rgb matrix
    alpha_m <- map(sampleid, ~spe@misc$alpha[[.x]]*img_alpha )
    im <- imap(im, ~array(c(.x, alpha_m[[.y]]), dim = c(dim(.x)[1:2],4)) )
    # select the size of the viewbox
    img_ <- map(im, ~.x[l$min_row:l$max_row,l$min_col:l$max_col,])
    
    # img <- map(im, ~as.raster(.x))
    # plot(img[[1]])
    
    # get grob and save as list
    grob <- map(img_, ~grid::rasterGrob(.x, width=unit(1,"npc"), height=unit(1,"npc")))
    images_tibble <- tibble(sample=factor(sampleid), grob=grob, orig.ident = sampleid)
    
    spatial_image <- geom_spatial(data=images_tibble, aes(grob=grob), x=0.5, y=0.5)}
  else{spatial_image <- NULL}
  
  ## Spatial annotation:
  if(sp_annot){
    tools <- map(sampleid, ~pluck(spe@tools, .x)) %>% bind_rows(., .id = "orig.ident") %>%
      filter(!(grepl("fov|zoom|full_image", .$name))) # removes the black frame around image
    spatial_annotation <- geom_path(
      data=tools, 
      show.legend = FALSE, linewidth = annot_line,
      aes(x=x, y=y, group=interaction(elem_idx)), colour="#808080")
  }
  else{spatial_annotation <- NULL}
  
  p <- ggplot()+
    spatial_image +
    geom_point(data=df, aes(x=imagecol,y=imagerow, colour=.data[[geneid]]),
               shape = 16, #colour = "transparent", stroke = 0.5,
               size = point_size, alpha = alpha) +
    spatial_annotation + 
    colour_pallet +
    coord_fixed(ratio = aspect, expand=FALSE) +

    lims +
    txt +
    facet_wrap(~factor(orig.ident, levels = sampleid), ncol = ncol, dir = if(ncol > 2){"h"}else{"v"})
  
  # Hexagon shape:             
  # p <- p +
  #   ggstar::geom_star(
  #     starshape = "hexagon",
  #     size = point_size,
  #     #stroke = 0,
  #     alpha = alpha
  #   )
  
  # Define theme settings for tight space
  if(save_space){theme_tight <- theme(
    panel.spacing.x = unit(-1, "lines"),
    panel.spacing.y = unit(-3, "lines"),
    legend.box.margin = margin(0,10,0,-10), # moves the legend
    legend.title = element_text(size = 8),
    plot.margin = if (ncol == 2) unit(c(-2, -1, -2, -1), "lines") # t,r,b,l
    else unit(c(-2.5, -1, -2, -1), "lines") )}else{theme_tight <-NULL}
  
  p <- p +
    xlab("") +
    ylab("") +

    theme_void() +
    guides +
    theme(strip.text.x = element_blank(), # removes facet title
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.spacing.y = unit(-2, "lines"),
    ) + theme_tight

  return(p)
}

#################
# SP ANNOTATION #
#################
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

################
# GGPLOT THEME #
################
my_theme <-
  list(
    #scale_fill_manual(values = friendly_cols),
    #scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        #aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )
###############
# VIOLIN PLOT #
###############
# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
violin.fun <- function(
    obj, 
    feature, 
    facet="orig.ident", 
    fill="sample_name", 
    col_pal=NULL, 
    txt_size=7,
    dot_size=.1,
    n=1){
  if(is.null(col_pal)){col_pal <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC") }
  m <- max(obj[[feature]])/n # try e.g 2
  obj %>%
    ggplot(aes(.data[[facet]], .data[[feature]], fill=.data[[fill]])) +
    geom_violin(linewidth = .2) + ggtitle(feature) +
    geom_jitter(width = 0.3, alpha = 0.2, size=dot_size, col="black", stroke=0, shape=20) +
    scale_fill_manual(values = col_pal) +
    theme_minimal(base_size = 6) + NoLegend() + ylim(c(0, m)) +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 30, hjust=1),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) 
}

##############################
# FACET_WRAP DISCRETE GROUPS #
##############################
# orig.ident = sym("orig.ident")
# spe <- DATA_st
# geneid <- sym("CDH1")
# zoom <- sym("zoom")
# assay="SCDC"
# img_alpha = 0
# colors=c("grey90","grey80","grey60","navy","black")

plot_st_meta.fun <- function(
    spe,
    assay="RNA",
    sp_annot = TRUE,
    feat = "groups",
    save_space = TRUE,
    orig.ident = "orig.ident",
    lvls = c("P118", "P080", "P031", "P105", "P097", "P108", "P114", "P107"), 
    title = " ",
    image_id = "hires",
    alpha = 1,
    ncol = 4,
    spectral = TRUE,
    colors = NULL,
    annot_col = "#808080",
    annot_line = .3,
    point_size = 1.75,
    img_alpha = .5,
    zoom = "zoom" # one of "zoom" or "fov"
    ){
  
  spe <- mutate(spe, orig.ident = as.character(.data[[orig.ident]]))
  
  feature <- enquo(feat)
  feat <- sym(feat)
  orig.ident <- enquo(orig.ident)
  
  ID <- unique(pull(spe, orig.ident)) %>% set_names(.)
  # Set default assay
  DefaultAssay(spe) <- assay
  
  ## get feature to plot:
  if (!(as_label(feat) %in% colnames(spe@meta.data))) {
    spe <- spe %>%
      mutate(., FetchData(., vars = c(feat)) ) 
  }
  
  ## Colour pallets:
  if (is.null(colors)){
    # scales::show_col(disc_colors)
    disc_colors <- c(RColorBrewer::brewer.pal(9,"Pastel1"),
                     RColorBrewer::brewer.pal(9,"Set1"),
                     scales::hue_pal()(8),
                     RColorBrewer::brewer.pal(8,"Set2"),
                     RColorBrewer::brewer.pal(8,"Accent"),
                     
                     RColorBrewer::brewer.pal(8,"Pastel2") )
  }else{disc_colors <- colors}
  colour_pallet <- scale_fill_manual(values = disc_colors, aesthetics = c("colour"))
  guides <- guides(#fill = guide_legend(override.aes = list(size=2), keyheight = .5),
                   colour = guide_legend(override.aes = list(size=2), keyheight = .5))
  
  # get all spot coordinates:
  scale_fact <- map_dbl(ID, ~pluck(spe@images, .x, "scale.factors", "hires"))
  df <- map(ID, ~pluck(spe@images, .x, "coordinates")) %>%
    map2(., scale_fact, ~mutate(.x, scale_fact = .y)) %>%
    bind_rows() %>%
    mutate(imagecol = .$imagecol * .$scale_fact) %>%
    mutate(imagerow = .$imagerow * .$scale_fact) %>%
    cbind(.,as_tibble(select(spe, !!(feature)))) %>%
    cbind(.,as_tibble(select(spe, "orig.ident"=!!(orig.ident)))) %>%
    rownames_to_column(var = "barcode") %>%
    as_tibble() 
  
  gr <- unique(spe@meta.data[,c("orig.ident", "groups")])[,"groups"]
  text_annot <- tibble(sample_id = ID, x=450, y=800, orig.ident = ID, gr = gr) 
  
  # select viewframe:
  if (is.null(zoom)){zoom <- "zoom"}
  tools <- map(ID, ~pluck(spe@tools, .x)) %>% bind_rows(., .id = "orig.ident")
  l <- tools %>% 
    filter(.data[["name"]] == zoom) %>% 
    dplyr::rename("_row"=y, "_col"=x) %>%
    summarise(across("_row":"_col", 
                     list(min=min, max=max), 
                     .names = "{.fn}{.col}"))
  
  ## Spatial image:
  if (!(img_alpha == 0)){
    
    im <- map(ID, ~pluck(spe@images, .x, "image")) 
    # im is a 3D matrix, each dimension representing red, green or blue
    # a 4th dimension can be added to specify the opacity
    
    # add alpha dimension to the rgb matrix
    alpha_m <- map(ID, ~spe@misc$alpha[[.x]]*img_alpha )
    im <- imap(im, ~array(c(.x, alpha_m[[.y]]), dim = c(dim(.x)[1:2],4)) )
    
    # select the size of the viewbox
    # if using "full_img" you have to modify the code to use smallest max_col value of the images
    # otherwise you get a out of bounds error. I have not attempted to make this work
    img_ <- map(im, ~.x[l$min_row:l$max_row,l$min_col:l$max_col,])
    # img <- map(im_, ~as.raster(.x)) 
    # base::plot(img[[1]])
    
    # get grob and save as list
    grob <- map(img_, ~grid::rasterGrob(.x, width=unit(1,"npc"), height=unit(1,"npc")))
    images_tibble <- tibble(sample=factor(ID), grob=grob, orig.ident = ID)
    
    spatial_image <- geom_spatial(data=images_tibble, aes(grob=grob), x=0.5, y=0.5)}
  else{spatial_image <- NULL}
  
  ## Spatial annotation:
  if(sp_annot){
    tools <- map(ID, ~pluck(spe@tools, .x)) %>% bind_rows(., .id = "orig.ident") %>%
      filter(!(grepl("fov|zoom|full_image", .$name))) # removes the black frame around image
    spatial_annotation <- geom_path(
      data=tools, 
      show.legend = FALSE, linewidth = annot_line,
      aes(x=x, y=y, group=interaction(elem_idx)), colour=annot_col)
  }
  else{spatial_annotation <- NULL}
  #id_lab <- tibble(id=ID, x=rep(-Inf, length(ID)), y=rep(Inf, length(ID)))

  p <- ggplot() +
    spatial_image + 
    geom_point(data=df, aes(x=imagecol,y=imagerow, colour = .data[[feat]]),
               shape = 16, size = point_size, alpha = alpha) +
    
    colour_pallet +
    spatial_annotation + 
    #scale_color_manual(values=c(annot_col, "transparent")) +
    coord_equal(expand=FALSE) + 
    xlim(l$min_col,l$max_col) +
    ylim(l$max_row,l$min_row) +
    geom_text(aes(label = sample_id, x=x, y=y), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt) + # sample ID
    geom_text(aes(label = gr, x=x, y=y+90), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt) + # condition
    facet_wrap(~factor(orig.ident, levels = lvls), ncol = ncol, dir = if(ncol > 2){"h"}else{"v"} )
  
  #Hexagon shape:
  # p <- p +
  #   ggstar::geom_star(data=df, aes(x=imagecol,y=imagerow, fill=.data[[feat]], colour=.data[[feat]]),
  #     starshape = "hexagon",
  #     size = point_size,
  #     #stroke = 0,
  #     alpha = alpha
  #   )
  
  # Define theme settings for tight space
  if(save_space){theme_tight <- theme(
    axis.ticks.length = unit(0, "cm"),
    panel.spacing.x = unit(-1, "lines"),
    panel.spacing.y = unit(-3, "lines"),
    legend.box.margin = margin(0,15,0,-25), # moves the legend
    plot.margin = if (ncol == 2) unit(c(-2, -1, -2, -1), "lines") # t,r,b,l
    else unit(c(-2.5, -1, -2, -1), "lines") )}else{theme_tight <-NULL}
  
  p <- p +
    theme_void() +
    guides +
    theme(strip.text.x = element_blank(), # reoves facet title
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.spacing.y = unit(-2, "lines"),
    ) + theme_tight
  
  return(p)
}

##################################
# FACET_WRAP CONTINEOUS FEATURES #
##################################
# orig.ident = sym("orig.ident")
# spe <- DATA_st
# geneid <- sym("IGKC")
# zoom <- sym("zoom")
# assay="RNA"
# img_alpha = 0
# colors=c("grey90","grey80","grey60","navy","black")

plot_st_feat.fun <- function(
    spe,
    assay="RNA",
    sp_annot = TRUE,
    txt = TRUE,
    save_space = TRUE,
    geneid = "nFeature_RNA",
    orig.ident = "orig.ident",
    lvls = c("P118", "P080", "P031", "P105", "P097", "P108", "P114", "P107"),
    title = " ",
    image_id = "hires",
    alpha = 1,
    ncol = 4,
    col = c("grey95","grey70","navy","black"),
    scale = TRUE,
    mins = NULL, maxs = NULL,
    annot_col = "#808080",
    annot_line = .3,
    point_size = 1.75,
    img_alpha = .5,
    zoom = "zoom" # one of "zoom" or "fov"
    ){
  
  spe <- mutate(spe, orig.ident = as.character(.data[[orig.ident]]))
  
  gene <- sym(geneid)
  orig.ident <- enquo(orig.ident) 
  ID <- unique(pull(spe, orig.ident)) %>% set_names(.)
  
  # Set default assay
  DefaultAssay(spe) <- assay
  
  ## get feature to plot:
  if (!(as_label(gene) %in% colnames(spe@meta.data))) {
    spe <- spe %>%
      mutate(., FetchData(., vars = c(as_label(gene))) ) %>%
      mutate(feat = !!(gene))
  }
  feat_vec <- pull(spe, as_label(gene))
  
  # Colour pal:
  if(is.null(mins)){
    mins <- min(c(feat_vec, 0),na.rm = T)} # get 0 or negative value
  if(is.null(maxs)){maxs <- quantile(feat_vec,0.99,na.rm = T) 
  if(maxs==0){maxs <- max(feat_vec,na.rm = T)}} # get percentile
  
  if(max(feat_vec, na.rm=T) != 0){
    spe <- spe %>%
      mutate(feat = (!!(gene) - mins) / ( maxs - mins) ) %>%
      mutate(feat = ifelse(.$feat > 1, 1, .$feat))
  } # Calculate percentage
  
  spe <- spe %>%
    #select(1:3, feat, !!(gene)) %>%
    mutate(feat_val = !!(gene)) %>%
    #mutate(!!(gene) := round(.$feat*98)+1) #%>%
    mutate(feat = round(.$feat*98)+1) #%>%
    #mutate(pal = c( col[1],colorRampPalette(col[-1])(99))[.$feat] ) 
  
  # get scale factor:
  scale_fact <- map_dbl(ID, ~pluck(spe@images, .x, "scale.factors", "hires"))
  # get all spot coordinates:
  df <- map(ID, ~pluck(spe@images, .x, "coordinates")) %>%
    map2(., scale_fact, ~mutate(.x, scale_fact = .y)) %>%
    bind_rows() %>%
    mutate(imagecol = .$imagecol * .$scale_fact) %>%
    mutate(imagerow = .$imagerow * .$scale_fact) %>%
    cbind(.,as_tibble(select(spe, feat_val))) %>%
    cbind(.,as_tibble(select(spe, feat))) %>%
    cbind(.,as_tibble(select(spe, !!(gene)))) %>%
    cbind(.,as_tibble(select(spe, "orig.ident"=!!(orig.ident)))) %>%
    rownames_to_column(var = "barcode") %>%
    as_tibble() %>%
    arrange(feat_val) 
  
  if(scale == TRUE){df <- df %>% mutate(!!(gene) := feat)}
  
  # select viewframe:
  if (is.null(zoom)){zoom <- "zoom"}
  tools <- map(ID, ~pluck(spe@tools, .x)) %>% bind_rows(., .id = "orig.ident")
  l <- tools %>% 
    filter(.data[["name"]] == zoom) %>% 
    dplyr::rename("_row"=y, "_col"=x) %>%
    summarise(across("_row":"_col", 
                     list(min=min, max=max), 
                      .names = "{.fn}{.col}")) 

  ## Spatial image:
  if (!(img_alpha == 0)){
    
    im <- map(ID, ~pluck(spe@images, .x, "image")) 
    # im is a 3D matrix, each dimension representing red, green or blue
    # a 4th dimension can be added to specify the opacity
    
    # add alpha dimension to the rgb matrix
    alpha_m <- map(ID, ~spe@misc$alpha[[.x]]*img_alpha )
    im <- imap(im, ~array(c(.x, alpha_m[[.y]]), dim = c(dim(.x)[1:2],4)) )
    
    # select the size of the viewbox
    img_ <- map(im, ~.x[l$min_row:l$max_row,l$min_col:l$max_col,])
    # img <- map(im_, ~as.raster(.x)) 
    # base::plot(img[[1]])
    
    # get grob and save as list
    grob <- map(img_, ~grid::rasterGrob(.x, width=unit(1,"npc"), height=unit(1,"npc")))
    images_tibble <- tibble(sample=factor(ID), grob=grob, orig.ident = ID)
    
    spatial_image <- geom_spatial(data=images_tibble, aes(grob=grob), x=0.5, y=0.5)}
  else{spatial_image <- NULL}
  
  ## Spatial annotation:
  if(sp_annot){
    tools <- map(ID, ~pluck(spe@tools, .x)) %>% 
      bind_rows(., .id = "orig.ident") %>%
      filter(!(grepl("fov|zoom|full_image", .$name)))
    spatial_annotation <- geom_path(
      data=tools, 
      show.legend = FALSE, linewidth = annot_line,
      aes(x=x, y=y, group=interaction(elem_idx)), colour=annot_col)
  }
  else{spatial_annotation <- NULL}
  
  if(txt){
    gr <- unique(spe@meta.data[,c("orig.ident", "groups")])[,"groups"]
    text_annot <- tibble(sample_id = ID, x=450, y=700, orig.ident = ID, gr = gr) 
    txt_spacing <- if(save_space){c(130,60)}else{c(0,120)}
    
    txt <- list(geom_text(aes(label = sample_id, x=x, y=y+txt_spacing[2]), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt), # sample ID
                geom_text(aes(label = gr, x=x, y=y+txt_spacing[1]), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt) ) # condition
    }else{txt <- NULL}
  
  p <- ggplot() +
    spatial_image +
    geom_point(data=df, aes(x=imagecol,y=imagerow, colour=.data[[gene]]),
               stroke = 0, size = point_size, alpha = alpha) +
    #scale_colour_identity() +
    scale_color_gradientn(colours = col,
                          values = seq(from=0, to=1, along.with=col),
                          breaks = if(scale== TRUE){c(1,25,50,75,99)}else{waiver()},
                          na.value = "grey90",
                          guide = guide_colourbar(barwidth = .5, barheight = 5 )) + 
    txt +
    spatial_annotation + 
    coord_equal(expand=FALSE) +
    xlim(l$min_col,l$max_col) +
    ylim(l$max_row,l$min_row) +
    facet_wrap(~factor(orig.ident, levels = lvls), ncol = ncol, dir = if(ncol > 2){"h"}else{"v"})
  
  #Hexagon shape:
  # p <- p +
  #   ggstar::geom_star(data=df, aes(x=imagecol,y=imagerow, fill=.data[[gene]], colour=.data[[gene]]),
  #     starshape = "hexagon",
  #     size = point_size,
  #     #stroke = 0,
  #     alpha = alpha
  #   )
  
  # Define theme settings for tight space
  if(save_space){theme_tight <- theme(
    axis.ticks.length = unit(0, "cm"),
    panel.spacing.x = unit(-1, "lines"),
    panel.spacing.y = unit(-3, "lines"),
    legend.box.margin = margin(0,10,0,-10), # moves the legend
    plot.margin = if (ncol == 2) unit(c(-2, -1, -2, -1), "lines") # t,r,b,l
    else unit(c(-2.5, -1, -2, -1), "lines") )}else{theme_tight <-NULL}
  
  p <- p +
    theme_void() +
    theme(strip.text.x = element_blank(), # removes facet title
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.spacing.y = unit(-2, "lines"),
    ) + theme_tight
  
  return(p)
}

# spe <- DATA_st
# ID = c("P097")
# ct.res <- sym("cell_annot_1")
########################
# TISSUE PROP PIE PLOT #
########################
plot_cell_pie.fun <- function(
    spe,
    assay="RNA",
    sp_annot = TRUE,
    orig.ident = "orig.ident",
    ct.res = NULL,
    ct.select = NULL,
    radius_adj = -1,
    save_space = T,
    lvls = c("P118", "P080","P031", "P105", "P097","P108", "P114", "P107"),
    title = " ",
    image_id = "hires",
    alpha = 1,
    ncol = 4,
    spectral = TRUE,
    colors = NULL,
    annot_col = "#808080",
    annot_line = .3,
    img_alpha = 0,
    zoom = "zoom" ) {
  
  if(isFALSE(is.null(ct.res))){ct.res <- enquo(ct.res)}
  ID <- unique(pull(spe, orig.ident)) %>% set_names(.)
  ID_ <- ID
  orig.ident <- enquo(orig.ident)
  # Set default assay
  DefaultAssay(spe) <- assay
  
  # get all cell annotations
  if(assay == "celltypeprops"){
    cell_annot <- t(spe@assays$celltypeprops@data) %>%
      as_tibble(., rownames = "barcode") 
  }else{
    cell_annot <- spe@misc$cell_annot %>% 
    filter(grepl(paste0(ID_, collapse="|"), .cell)) %>%
    select(barcode=".cell", values, !!(ct.res)) %>%
    pivot_wider(., names_from = !!(ct.res), values_fill = 0,
                values_from = values, values_fn = function(x) sum(x)) }
  
    if(is.null(ct.select)){ct.select <- colnames(cell_annot)[2:ncol(cell_annot)]}
  
  # get all spot coordinates:
  scale_fact <- map_dbl(ID, ~pluck(spe@images, .x, "scale.factors", "hires"))
  df <- map(ID, ~pluck(spe@images, .x, "coordinates")) %>%
    map2(., scale_fact, ~mutate(.x, scale_fact = .y)) %>%
    bind_rows(., .id = "orig.ident") %>%
    mutate(imagecol = .$imagecol * .$scale_fact) %>%
    mutate(imagerow = .$imagerow * .$scale_fact) %>%
    rownames_to_column(var = "barcode") %>%
    select(-any_of(c("epi", "SubMuc")), -contains("_")) %>%
    left_join(., cell_annot, by="barcode") %>%
    #mutate(orig.ident = !!(facet)) %>%
    #mutate(orig.ident = factor(.data[["orig.ident"]], levels = lvls)) %>%
    #cbind(.,as_tibble(select(spe, groups))) %>%
    as_tibble() 
  
  gr <- unique(spe@meta.data[,c("orig.ident", "groups")])[,"groups"]
  text_annot <- tibble(sample_id = ID, x=500, y=500, orig.ident = ID, gr = gr) 
  
  ## Colour pallets:
  if (is.null(colors)){
    # scales::show_col(disc_colors)
    cell_col <- c(RColorBrewer::brewer.pal(9,"Pastel1"),
                  RColorBrewer::brewer.pal(9,"Set1"),
                  scales::hue_pal()(8),
                  RColorBrewer::brewer.pal(8,"Set2"),
                  RColorBrewer::brewer.pal(8,"Accent"),
                  
                  RColorBrewer::brewer.pal(8,"Pastel2") )
    cell_col <- set_names(cell_col[1:length(ct.select)], ct.select)
  }else{cell_col <- colors}
  colour_pallet <- scale_fill_manual(values = cell_col, na.value = "grey90" )
  guides <- guides(fill=guide_legend(ncol=1,title = ""))
  
  # select viewframe:
  if (!(is.null(zoom))){
    tools <- map(ID, ~pluck(spe@tools, .x)) %>% bind_rows(., .id = "orig.ident")
    l <- tools %>% 
      filter(.data[["name"]] == zoom) %>% # zoom <- "zoom"
      #select(row=imagerow)
      dplyr::rename("_row"=y, "_col"=x) %>%
      summarise(across("_row":"_col", 
                       list(min=min, max=max), 
                       .names = "{.fn}{.col}")) 
  }
  # else{l <- tibble( min_col = 0, max_col = ncol(img),
  #                   min_row = 0, max_row = nrow(img))}
  
  ## Spatial image:
  if (!(img_alpha == 0)){
    
    # set image alpha:
    im <- map(ID, ~pluck(spe@images, .x, "image")) 
    im <- map(im, ~ matrix(
      rgb(.x[,,1],.x[,,2],.x[,,3], .x[4,,]* img_alpha), nrow=dim(.x)[1]))
    
    img <- map(im, ~as.raster(.x))
    img_ <- map(img, ~.x[l$min_row:l$max_row,l$min_col:l$max_col])
    
    # get grob and save as list
    grob <- map(img_, ~grid::rasterGrob(.x, width=unit(1,"npc"), height=unit(1,"npc")))
    images_tibble <- tibble(sample=factor(ID), grob=grob, orig.ident = ID)
    
    spatial_image <- geom_spatial(data=images_tibble, aes(grob=grob), x=0.5, y=0.5)
  }
  else{spatial_image <- NULL}
  
  ## Spatial annotation:
  if(sp_annot){
    tools <- map(ID, ~pluck(spe@tools, .x)) %>% bind_rows(., .id = "orig.ident") %>%
      filter(!(grepl("fov|zoom|full_image", .$name))) # removes the black frame around image
    spatial_annotation <- geom_path(
      data=tools, 
      show.legend = FALSE, linewidth = annot_line,
      aes(x=x, y=y, group=interaction(elem_idx)), colour=annot_col)
  }
  else{spatial_annotation <- NULL}
  
  radius = (max(df$imagecol) - min(df$imagecol)) * (max(df$imagerow) - min(df$imagerow))
  radius = radius / nrow(df)
  radius = radius / pi
  radius = sqrt(radius) * 0.85
  
  p <- ggplot() +
    geom_scatterpie(data=na.omit(df), aes(x=imagecol,y=imagerow, r=radius+radius_adj), cols=ct.select, color=NA) +
    
    colour_pallet +
    spatial_image + 
    spatial_annotation +
    coord_equal(expand=FALSE ) + #theme(l) +
    xlim(l$min_col,l$max_col) +
    ylim(l$max_row,l$min_row) +
    #geom_text(aes(label = sample_id, x=x, y=y), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt) + # sample ID
    #geom_text(aes(label = gr, x=x, y=y+120), data = text_annot, inherit.aes = F, hjust = 0, size = 8/.pt) + # condition
    facet_wrap(~factor(orig.ident, levels = lvls), ncol = ncol,  dir = if(ncol > 2){"h"}else{"v"})
  # facet_wrap(vars(!!(orig.ident)), ncol = ncol)
  
  #Hexagon shape:
  # p <- p +
  #   ggstar::geom_star(data=df, aes(x=imagecol,y=imagerow, fill=.data[[feat]], colour=.data[[feat]]),
  #     starshape = "hexagon",
  #     size = point_size,
  #     #stroke = 0,
  #     alpha = alpha
  #   )
  
  
  # Define theme settings for tight space
  if(save_space){theme_tight <- theme(
    panel.spacing.x = unit(-1, "lines"),
    panel.spacing.y = unit(-3, "lines"),
    legend.box.margin = margin(0,10,0,-10), # moves the legend
    legend.title = element_text(size = 8),
    plot.margin = if (ncol == 2) unit(c(-2, -1, -2, -1), "lines") # t,r,b,l
    else unit(c(-2.5, -1, -2, -1), "lines") )}else{theme_tight <-NULL}
  
  p <- p +
    xlab("") +
    ylab("") +

    theme_void()+
    guides+
    theme(strip.text.x = element_blank(), # reoves facet title
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.spacing.y = unit(-2, "lines"),
    ) + theme_tight
  
  return(p)
}

