#' Generate code files required for shiny app (one dataset)
#'
#' Generate code files required for shiny app containing only one dataset. In 
#' particular, two R scripts will be generated, namely \code{server.R} and 
#' \code{ui.R}. If users want to include multiple dataset in one shiny app, 
#' please use \code{makeShinyCodesMulti()} instead. Note that both 
#' \code{makeShinyFiles} and \code{makeShinyCodes} functions are ran when 
#' running the wrapper function \code{makeShinyApp}.
#'
#' @param shiny.title title for shiny app
#' @param shiny.footnotes text for shiny app footnote
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#'
#' @return server.R and ui.R required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table readr glue
#'
#' @examples
#' makeShinyCodes(shiny.title = "scRNA-seq shiny app", shiny.footnotes = '""',
#'                shiny.prefix = "sc1", shiny.dir = "shinyApp/")
#'
#' @export
makeShinyCodes <- function(shiny.title, shiny.footnotes,
                           shiny.prefix, shiny.dir){
  if(packageVersion("readr") >= "1.4.0"){
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr","ggplot2",
        "hdf5r","ggdendro","gridExtra")), file = fname)
    readr::write_file(wrSVload(shiny.prefix), append = TRUE, file = fname)
    readr::write_file(wrSVfix(), append = TRUE, file = fname)
    readr::write_file(wrSVmain(shiny.prefix), append = TRUE, file = fname)
    readr::write_file(wrSVend(), append = TRUE, file = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr")), file = fname)
    readr::write_file(wrUIload(shiny.prefix), append = TRUE, file = fname)
    readr::write_file(wrUIsingle(shiny.title), append = TRUE, file = fname)
    readr::write_file(wrUImain(shiny.prefix), append = TRUE, file = fname)
    readr::write_file(glue::glue(', \n'), append = TRUE, file = fname)
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, file = fname)
    
  } else {
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr","ggplot2",
        "hdf5r","ggdendro","gridExtra")), path = fname)
    readr::write_file(wrSVload(shiny.prefix), append = TRUE, path = fname)
    readr::write_file(wrSVfix(), append = TRUE, path = fname)
    readr::write_file(wrSVmain(shiny.prefix), append = TRUE, path = fname)
    readr::write_file(wrSVend(), append = TRUE, path = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr")), path = fname)
    readr::write_file(wrUIload(shiny.prefix), append = TRUE, path = fname)
    readr::write_file(wrUIsingle(shiny.title), append = TRUE, path = fname)
    readr::write_file(wrUImain(shiny.prefix), append = TRUE, path = fname)
    readr::write_file(glue::glue(', \n'), append = TRUE, path = fname)
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, path = fname)
  }
  
}


