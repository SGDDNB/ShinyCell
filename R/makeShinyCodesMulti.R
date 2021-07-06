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
#' @param shiny.footnotes text for shiny app footnote. When given as a list, 
#'   citation can be inserted by specifying author, title, journal, volume, 
#'   page, year, doi, link. See example below. 
#' @param shiny.prefix specify file prefix for each dataset. Must match the 
#'   prefix used in \code{makeShinyFiles}
#' @param shiny.headers specify the tab header names for each dataset. Length 
#'   must match that of \code{shiny.prefix}
#' @param shiny.dir specify directory to create the shiny app in
#' @param enableSubset specify whether to enable "Toggle to subset cells" 
#'   functionality in the shiny app. Default is to enable this functionality
#' @param defPtSiz specify default point size for single cells. For example, a 
#'   smaller size can be used if you have many cells in your dataset. A single 
#'   value can be specified to set the point size for all datasets. Otherwise,
#'   users have to specify one value for each dataset
#' @param ganalytics Google analytics tracking ID (e.g. "UA-123456789-0")
#'
#' @return server.R and ui.R required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table readr
#'
#' @examples
#' # Example citation
#' citation = list(
#'   author  = "Liu X., Ouyang J.F., Rossello F.J. et al.",
#'   title   = "",
#'   journal = "Nature",
#'   volume  = "586",
#'   page    = "101-107",
#'   year    = "2020", 
#'   doi     = "10.1038/s41586-020-2734-6",
#'   link    = "https://www.nature.com/articles/s41586-020-2734-6")
#' makeShinyCodes(shiny.title = "scRNA-seq shiny app", shiny.footnotes = "",
#'                shiny.prefix = c("sc1", "sc2"), defPtSiz = c(1.25, 1.5),
#'                shiny.headers = c("dataset1", "dataset2"),
#'                shiny.dir = "shinyApp/")
#'                
#' @export
makeShinyCodesMulti <- function(shiny.title, shiny.footnotes,
                                shiny.prefix, shiny.headers, shiny.dir, 
                                enableSubset = TRUE, defPtSiz = 1.25,
                                ganalytics = NA){
  ### Checks
  if(length(shiny.prefix) != length(shiny.headers)){
    stop("length of shiny.prefix and shiny.headers does not match!")
  }
  subst = "#"
  if(enableSubset){subst = ""}
  if(length(shiny.prefix) != length(defPtSiz)){
    defPtSiz = rep(defPtSiz[1], length(shiny.prefix))
  }
  defPtSiz = as.character(defPtSiz)
  
  if(packageVersion("readr") >= "1.4.0"){
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr","ggplot2",
        "ggrepel","hdf5r","ggdendro","gridExtra")), file = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVload(i), append = TRUE, file = fname)
    }
    readr::write_file(wrSVfix(), append = TRUE, file = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVmain(i, subst), append = TRUE, file = fname)
    }
    readr::write_file(wrSVend(), append = TRUE, file = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr")), file = fname)
    for(i in shiny.prefix){
      readr::write_file(wrUIload(i), append = TRUE, file = fname)
    }
    readr::write_file(wrUIsingle(shiny.title, ganalytics), append = TRUE, file = fname)
    for(i in seq_along(shiny.prefix)){
      hhh = shiny.headers[i]
      readr::write_file(glue::glue('navbarMenu("{hhh}",'), 
                        append = TRUE, file = fname)
      readr::write_file(wrUImain(shiny.prefix[i], subst, defPtSiz[i]), 
                        append = TRUE, file = fname)
      readr::write_file(glue::glue('), \n\n\n'), append = TRUE, file = fname)
    }
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, file = fname)
    
    
    ### Write code for google-analytics.html
    if(!is.na(ganalytics)){
      fname = paste0(shiny.dir, "/google-analytics.html")
      readr::write_file(wrUIga(ganalytics), file = fname)
    }
    
  } else {
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr","ggplot2",
        "ggrepel","hdf5r","ggdendro","gridExtra")), path = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVload(i), append = TRUE, path = fname)
    }
    readr::write_file(wrSVfix(), append = TRUE, path = fname)
    for(i in shiny.prefix){
      readr::write_file(wrSVmain(i, subst), append = TRUE, path = fname)
    }
    readr::write_file(wrSVend(), append = TRUE, path = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr")), path = fname)
    for(i in shiny.prefix){
      readr::write_file(wrUIload(i), append = TRUE, path = fname)
    }
    readr::write_file(wrUIsingle(shiny.title, ganalytics), append = TRUE, path = fname)
    for(i in seq_along(shiny.prefix)){
      hhh = shiny.headers[i]
      readr::write_file(glue::glue('navbarMenu("{hhh}",'), 
                        append = TRUE, path = fname)
      readr::write_file(wrUImain(shiny.prefix[i], subst, defPtSiz[i]), 
                        append = TRUE, path = fname)
      readr::write_file(glue::glue('), \n\n\n'), append = TRUE, path = fname)
    }
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, path = fname)
    
    
    ### Write code for google-analytics.html
    if(!is.na(ganalytics)){
      fname = paste0(shiny.dir, "/google-analytics.html")
      readr::write_file(wrUIga(ganalytics), path = fname)
    }
  }
  
}


