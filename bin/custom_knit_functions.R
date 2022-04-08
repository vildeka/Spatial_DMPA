################
# KNIT SCRIPTS #
################
knit_script <- function(inputFile, out_dir, ...) {
  rmarkdown::render(
    inputFile,
    output_file = paste0(
      xfun::sans_ext(inputFile), '_', Sys.Date(),'.html'),
    output_dir = paste0("../lab_book/",out_dir,"/"),
    output_format = rmarkdown::html_document(
      echo = TRUE,
      toc = FALSE,                             
      toc_float = FALSE,
      code_folding = "show",
      code_download = TRUE)
  )
}

################
# KNIT TO GITHUB #
################
knit_github <- function(inputFile, out_dir, ...) {
  name <- xfun::sans_ext(inputFile)
  out <- paste0(out_dir, basename(name), '_', Sys.Date())
  rmarkdown::render(inputFile,
                    output_file = out,
                    output_dir = out_dir,
                    output_format = rmarkdown::html_document(
                      preserve_yaml = TRUE,
                      dev= c("png","pdf"),
                      keep_md = TRUE,
                      echo=TRUE,
                      toc = FALSE,
                      toc_float = FALSE,
                      code_folding = "show",
                      code_download = TRUE
                      )
  )
  format <- setNames(c("gfm", "pdf"), c(".md",".pdf"))
  #format <- setNames(c("gfm", "pdf", "html"), c(".md",".pdf", '.html'))
  purrr::imap(format,
              ~rmarkdown::pandoc_convert(input=paste0(getwd(),"/",out,".md"), to = .x,
                                         output=paste0(name, .y),
                                         #options="-M echo=FALSE"
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