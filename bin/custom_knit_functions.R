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
pkg_file <- function(...) {
  system.file(..., package = "rmarkdown")
}
pkg_file_arg <- function(...) {
  rmarkdown::pandoc_path_arg(pkg_file(...))
}
css <- pkg_file_arg(
  "rmarkdown/templates/github_document/resources/github.css")



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
                      keep_md = TRUE,
                      echo=TRUE,
                      toc = FALSE,
                      toc_float = FALSE,
                      code_folding = "show",
                      code_download = TRUE
                      )
  )
  fs::file_move(path = paste0(out_dir,"/",out,".md"), new_path = paste0("./",basename(name), ".md"))
  format <- setNames(c("gfm"), c(".md"))
  #format <- setNames(c("gfm", "pdf", "html"), c(".md",".pdf", '.html'))
  purrr::imap(format,
              ~rmarkdown::pandoc_convert(input= paste0("./",basename(name), ".md"),#paste0(getwd(),"/",out_dir,"/",out,".md"),
                                         to = .x,
                                         output=paste0(name, .y),
                                         options=paste0('-M title=', basename(name), " standalone=TRUE", "date=", Sys.Date())
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