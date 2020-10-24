#' Generate code files required for shiny app (multi datasets)
#'
#' Generate code files required for shiny app containing multiple datasets. In 
#' particular, two R scripts will be generated, namely \code{server.R} and 
#' \code{ui.R}. Note that \code{makeShinyFiles} has to be ran prior to 
#' generate the necessary data files for each of the dataset to be included. 
#' The prefix used in \code{makeShinyFiles} have to be then supplied in this 
#' function.
#'
#' @param shiny.title specify the overall title for shiny app
#' @param shiny.footnotes text for shiny app footnote
#' @param shiny.prefix specify file prefix for each dataset. Must match the 
#'   prefix used in \code{makeShinyFiles}
#' @param shiny.headers specify the tab header names for each dataset. Length 
#'   must match that of \code{shiny.prefix}
#' @param shiny.dir specify directory to create the shiny app in
#'
#' @return server.R and ui.R required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table readr
#'
#' @examples
#' makeShinyCodes(shiny.title = "scRNA-seq shiny app", shiny.footnotes = '""',
#'                shiny.prefix = c("sc1", "sc2"), 
#'                shiny.headers = c("dataset1", "dataset2"),
#'                shiny.dir = "shinyApp/")
#'                
#' @export
makeShinyCodesMulti <- function(shiny.title, shiny.footnotes,
                                shiny.prefix, shiny.headers, shiny.dir){
  ### Checks
  if(length(shiny.prefix) != length(shiny.headers)){
    stop("length of shiny.prefix and shiny.headers does not match!")
  }
  
  if(packageVersion("readr") >= "1.4.0"){
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr","ggplot2",
        "hdf5r","ggdendro","gridExtra")), file = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVload(i), append = TRUE, file = fname)
    }
    readr::write_file(wrSVfix(), append = TRUE, file = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVmain(i), append = TRUE, file = fname)
    }
    readr::write_file(wrSVend(), append = TRUE, file = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr")), file = fname)
    for(i in shiny.prefix){
      readr::write_file(wrUIload(i), append = TRUE, file = fname)
    }
    readr::write_file(wrUIsingle(shiny.title), append = TRUE, file = fname)
    for(i in seq_along(shiny.prefix)){
      hhh = shiny.headers[i]
      readr::write_file(glue::glue('navbarMenu("{hhh}",'), 
                        append = TRUE, file = fname)
      readr::write_file(wrUImain(shiny.prefix[i]), append = TRUE, file = fname)
      readr::write_file(glue::glue('), \n\n\n'), append = TRUE, file = fname)
    }
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, file = fname)
    
  } else {
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr","ggplot2",
        "hdf5r","ggdendro","gridExtra")), path = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVload(i), append = TRUE, path = fname)
    }
    readr::write_file(wrSVfix(), append = TRUE, path = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVmain(i), append = TRUE, path = fname)
    }
    readr::write_file(wrSVend(), append = TRUE, path = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","magrittr")), path = fname)
    for(i in shiny.prefix){
      readr::write_file(wrUIload(i), append = TRUE, path = fname)
    }
    readr::write_file(wrUIsingle(shiny.title), append = TRUE, path = fname)
    for(i in seq_along(shiny.prefix)){
      hhh = shiny.headers[i]
      readr::write_file(glue::glue('navbarMenu("{hhh}",'), 
                        append = TRUE, path = fname)
      readr::write_file(wrUImain(shiny.prefix[i]), append = TRUE, path = fname)
      readr::write_file(glue::glue('), \n\n\n'), append = TRUE, path = fname)
    }
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, path = fname)
  }
  
}


