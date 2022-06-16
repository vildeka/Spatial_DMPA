#################
# GEOM SPATIAL #
################
geom_spatial <-  function(mapping = NULL,
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

#####################
# PLOTTING FUNCTION #
#####################
# spe <- DATA
# geneid <- sym("CDH1")
# zoom <- sym("tissue")

vis_gene_p <- function(
  spe,
  sampleid = c("P080"),
  spatial=TRUE,
  geneid="CDH1", #"LINC01135",# "nFeature_RNA"
  title=" ",
  spectral= TRUE,
  image_id = "hires",
  alpha = .7,
  cont_colors = c("lightgray", "mistyrose", "red", "dark red", "black"),
  point_size = 1.75,
  img_alpha = .5,
  zoom = NULL ) {
  
  if (!(as_label(geneid) %in% colnames(spe@meta.data))) {
    spe <- spe %>%
      mutate(., FetchData(., vars = c(geneid)) ) 
  }
  
  # Colour pallets:
  if (spectral){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    cont_colors <- myPalette(100)
    colour_pallet <- scale_fill_gradientn(colours = cont_colors)
  }
  if (is.character(pull(spe, geneid))){
    # scales::show_col(disc_colors)
    disc_colors <- c(RColorBrewer::brewer.pal(9,"Pastel1"),
                     RColorBrewer::brewer.pal(9,"Set1"),
                     scales::hue_pal()(8),
                     RColorBrewer::brewer.pal(8,"Set2"),
                     RColorBrewer::brewer.pal(8,"Accent"),
                     RColorBrewer::brewer.pal(8,"Pastel2") )
    colour_pallet <- scale_fill_manual(values = disc_colors)
  }
  
  # filter samples:
  spe <- spe %>% filter(orig.ident %in% sampleid)
  # get scale factor:
  scale_fact <- spe@images[[sampleid]]@scale.factors[[image_id]]
  geneid <- enquo(geneid)
  
  df <- spe@images[[sampleid]]@coordinates %>%
    mutate(imagecol = .$imagecol * scale_fact) %>%
    mutate(imagerow = .$imagerow * scale_fact) %>%
    cbind(.,as_tibble(select(spe, !!(geneid)))) %>%
    rownames_to_column(var = "barcode") %>%
    as_tibble() 
  
  # set image alpha:
  im <- spe@images[[sampleid]]@image
  spe@images[[sampleid]]@image <- matrix(
    rgb(im[,,1],im[,,2],im[,,3], im[4,,]* img_alpha), nrow=dim(im)[1]) 
  
  img <- GetImage(spe, image = sampleid, mode = "raster")
  
  # select viewframe:
  if (!(is.null(zoom))) {
    l <- df %>% 
      filter(.data[[zoom]] == 1) %>% 
      #select(row=imagerow)
      rename("_row"=imagerow, "_col"=imagecol) %>%
      summarise(across("_row":"_col", 
                       list(min=min, max=max), 
                       .names = "{.fn}{.col}")) 
  }else{l <- tibble( min_col = 0, max_col = ncol(img),
                     min_row = 0, max_row = nrow(img))}
  img_ <- img[l$min_row:l$max_row,l$min_col:l$max_col]
  
  # get grob and save as list
  grob <- grid::rasterGrob(img_, width=unit(1,"npc"), height=unit(1,"npc"))
  images_tibble <- tibble(sample=factor(sampleid), grob=list(grob))
  
  p <- ggplot()+
    geom_spatial(data=images_tibble, aes(grob=grob), x=0.5, y=0.5)+
    geom_point(data=df, aes(x=imagecol,y=imagerow,fill=.data[[geneid]]),
               shape = 21, #colour = "black",
               size = point_size, stroke = 0.5, alpha = alpha)+
    geom_path(data = spe@tools[[sampleid]], aes(x=x, y=y, group = interaction(elem_idx)))+
    colour_pallet +
    coord_cartesian(expand=FALSE )+
    
    xlim(l$min_col,l$max_col) +
    ylim(l$max_row,l$min_row)
  
  p <- p +
    xlab("") +
    ylab("") +
    ggtitle(sampleid)+
    #labs(fill = "Total UMI")+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  return(p)
}

#################
# PATH FUNCTIONS #
#################
c_curve<-function(start, c1, c2, end, n=100) {
  c1 <- start+c1
  c2 <- start+c2
  end <- start+end
  sapply((0:n)/n, function(i) {
    q1 <- start+i*(c1-start)
    q2 <- c1+i*(c2-c1)
    q3 <- c2+i*(end-c2)
    r1 <- q1+i*(q2-q1)
    r2 <- q2+i*(q3-q2)
    r1+i*(r2-r1)
  })
}
# https://stackoverflow.com/questions/46221703/converting-svg-to-ggplot
parse_svg_path <- function(path) {
  parts <- regmatches(path, gregexpr("[A-Za-z]|,|-?[0-9.]+", path, perl=T))[[1]] 
  parts <- parts[parts!=","]
  vals <- suppressWarnings(as.numeric(parts))
  i <- 1
  points <- matrix(ncol=0, nrow=2)
  while(i < length(parts)) {
    if (parts[i]=="M") {
      points <- cbind(points, c(vals[i+1], vals[i+2]))
      i <- i+3
    } else if (parts[i]=="l") {
      cpoints <- c_curve(
        points[, ncol(points)],
        c(vals[i+1], vals[i+2]),
        c(vals[i+3], vals[i+4]),
        c(vals[i+5], vals[i+6])
      )
      points <- cbind(points, cpoints)
      i <- i+7
    } else {
      stop(paste("unrecognized command", parts[i]))
    }
  }
  points
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
