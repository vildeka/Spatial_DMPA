#####################
# QUARTO KNIT HTML #
#####################
custom_knit <- function(inputFile, fig_path, ...) {
  #file_name = paste0(basename(xfun::sans_ext(inputFile)), '_', Sys.Date(),'.html')
  quarto::quarto_render(
    inputFile,
    output_file = paste0(basename(xfun::sans_ext(inputFile)), '_', Sys.Date(),'.html'),
    execute_params = list(fig.path = fig_path),
    #output_dir = paste0("../lab_book/",out_dir,"/"),
    output_format = "html",
    as_job = FALSE
  )#;
  # quarto::quarto_render(
  #   inputFile,
  #   output_format = "gfm",
  #   output_file = file_name
  #   )
}

# for script to run in hte background
# cat(
#   'res <- quarto::quarto_render("', basename(rmd_path), '",',
#   '"', paste0(basename(xfun::sans_ext(rmd_path)), '_', Sys.Date(), '.html'), '",',
#   'execute_params = list(fig.path =','"',fig_path, '"),',
#   ' output_format = "html"', ')\n',
#   sep = ""
# )


###################
# QUARTO KNIT GIT #
###################
git_knit <- function(inputFile, fig_path, ...) {
  #file_name = paste0(basename(xfun::sans_ext(inputFile)))
  #dir =  paste0(getwd(),"/md_files/")
  out_file = paste0(basename(xfun::sans_ext(inputFile)), ".md")
  #new_filename = paste0(dirname(inputFile),"/md_files/", basename(inputFile))
  #fs::file_copy(path = inputFile, new_path = new_filename)
  quarto::quarto_render(
    new_filename,
    output_file = out_file,
    execute_params = list(fig.path = fig_path),
    #execute_dir = "file",
    #output_dir = paste0("../lab_book/",out_dir,"/"),
    output_format = "gfm"
  )
  #fs::file_delete(new_filename)
  #fs::file_move(path = out_file, new_path = paste0(dirname(inputFile),"/md_files/", out_file))
}

#####################
# GIT KNIT FUNCTION #
#####################
rmd_custom_knit <- function(inputFile, out_dir, ...) {
  rmarkdown::render(
    inputFile,
    output_file = paste0(
      xfun::sans_ext(inputFile), '_', Sys.Date(),'.html'),
    output_dir = paste0("../lab_book/",out_dir,"/"),
    output_format = rmarkdown::html_document(
      dev = c("jpeg", "tiff", "pdf"),
      dpi = 300,
      echo=TRUE,
      #keep_md = TRUE,
      toc = FALSE,                             
      toc_float = FALSE,
      code_folding = "show",
      code_download = TRUE)
  )#;
  # rmarkdown::render(
  #   inputFile,
  #   output_format = rmarkdown::pdf_document()
  #   )
}

################
# KNIT TO GITHUB #
################
pkg_file <- function(...) {
  system.file(..., package = "rmarkdown")
}
pkg_file_arg <- function(...) {
  rmarkdown::pandoc_path_arg(pkg_file(...))
}
css <- pkg_file_arg(
  "rmarkdown/templates/github_document/resources/github.css")

# inputFile <- "/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/src/03_clustering_st_data.Rmd"
# out_dir <- "/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/src/md_files"
# setwd("/Users/vilkal/work/Brolidens_work/Projects/Spatial_DMPA/src")

knit_github <- function(inputFile, out_dir, ...) {
  name <- xfun::sans_ext(inputFile)
  out <- paste0(basename(name), '_', Sys.Date())
  # args <- c("--highlight-style", "pygments",
  #           "--template", pkg_file_arg(
  #             "rmarkdown/templates/github_document/resources/preview.html"),
  #           "--metadata", paste0("title=", basename(name), "pagetitle=PREVIEW"),
  #           "--variable", paste0("github-markdown-css:", css),
  #           c("--metadata", "pagetitle=PREVIEW")
  # )
  rmarkdown::render(inputFile,
                    output_file = out, #c(paste0(out_dir,"/",out,".html"), paste0("./",basename(name), ".md")), #"out, "
                    output_dir = out_dir,
                    output_format = rmarkdown::html_document(
                      highlight = "default",
                      #options=args,
                      dev= c("png","pdf"),
                      dpi = 300,
                      keep_md = TRUE,
                      echo=TRUE,
                      toc = FALSE,
                      toc_float = FALSE,
                      code_folding = "show",
                      code_download = TRUE
                      )
  )
  #fs::file_move(path = paste0(out_dir,"/",out,".md"), new_path = paste0("./md_files/",basename(name), ".md"))
  format <- setNames(c("gfm"), c(".md"))
  #format <- setNames(c("gfm", "pdf", "html"), c(".md",".pdf", '.html'))
  purrr::imap(format,
              ~rmarkdown::pandoc_convert(input= paste0(getwd(),"/",out_dir,"/",out,".md"), # paste0("./",basename(name), ".md"),
                                         to = .x,
                                         output=paste0(getwd(),"/md_files/",basename(name), .y), # paste0(name, .y),
                                         options=paste0('-M title=', basename(name), " standalone=TRUE ", "date=", Sys.Date(),
                                                        ' --extract-media=', getwd(),"./Figures/")
                                         # options=sub("-FILE-", paste0(out,".md"), "--lua-filter additional-metadata.lua --metadata default_meta_file:-FILE-.yaml")
                                         ))
}


# the code in the chuncks are not evaluated after the rmarkdown document has been knited. 
# However there is an option to keep the YAML header for md_documents 

################
# KNIT FIGURES #
################
knit_figure <- function(inputFile, out_dir, ...) {
  rmarkdown::render(
    inputFile,
    output_file = paste0(
      xfun::sans_ext(inputFile), '_', Sys.Date(),'.html'),
    output_dir = paste0("../../lab_book/",out_dir,"/"),
    output_format = rmarkdown::html_document(
      echo=TRUE,
      toc = FALSE,                             
      toc_float = FALSE,
      code_folding = "show",
      code_download = TRUE)
  );
  rmarkdown::render(
    inputFile,
    output_format = "pdf_document"
  )
}

###############
# KNIT TABLES #
###############
knit_table <- function(inputFile, ...) {
  rmarkdown::render(
    inputFile,
    output_file = paste0(
      xfun::sans_ext(inputFile), '_', Sys.Date(),'.html'),
    output_dir = "../lab_book/Clinical_tables/",
    output_format = rmarkdown::html_document(
      echo = TRUE,
      toc = FALSE,                             
      toc_float = FALSE,
      code_folding = "show",
      code_download = TRUE)
  );
  rmarkdown::render(
    inputFile,
    output_format = rmarkdown::word_document(
      reference_docx = '/Users/vilkal/work/Brolidens_work/Projects/DMPA/data/Klinsik_template_02.dotx',
      df_print = "kable")
  )
}


# rmarkdown::render(
#   inputFile,
#   output_format = rmarkdown::pdf_document(
#     dev = c("pdf", "png"),
#     dpi = 300)
# );