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
#' @param shiny.footnotes text for shiny app footnote. When given as a list, 
#'   citation can be inserted by specifying author, title, journal, volume, 
#'   page, year, doi, link. See example below. 
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#' @param enableSubset specify whether to enable "Toggle to subset cells" 
#'   functionality in the shiny app. Default is to enable this functionality
#' @param defPtSiz specify default point size for single cells. For example, a 
#'   smaller size can be used if you have many cells in your dataset
#' @param ganalytics Google analytics tracking ID (e.g. "UA-123456789-0")
#'
#' @return server.R and ui.R required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table readr glue
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
#'                shiny.prefix = "sc1", shiny.dir = "shinyApp/")
#'
#' @export
makeShinyCodes <- function(shiny.title, shiny.footnotes,
                           shiny.prefix, shiny.dir, 
                           enableSubset = TRUE, defPtSiz = 1.25,
                           ganalytics = NA, markers.all = FALSE, markers.top20 = FALSE, de.genes = FALSE, gene.ranks=FALSE){
  subst = "#"
  if(enableSubset){subst = ""}
  defPtSiz = as.character(defPtSiz)
  
  if(packageVersion("readr") >= "1.4.0"){
    ### Write code for server.R
    fname = paste0(shiny.dir, "/server.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr","ggplot2",
        "ggrepel","hdf5r","ggdendro","gridExtra","AUCell","rbokeh")), file = fname)
    readr::write_file(wrSVload(shiny.prefix, markers.all, markers.top20, de.genes), append = TRUE, file = fname)
    readr::write_file(wrSVfix(), append = TRUE, path = fname)
    readr::write_file(wrSVmain(shiny.prefix, subst), append = TRUE, file = fname)
    readr::write_file(wrSVmarkersDE(shiny.prefix, markers.all, markers.top20, de.genes), append = TRUE, file = fname)
    readr::write_file(wrSVgeneSig(shiny.prefix, gene.ranks), append = TRUE, file = fname)
    readr::write_file(wrSVend(), append = TRUE, file = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr")), file = fname)
    readr::write_file(wrUIload(shiny.prefix), append = TRUE, file = fname)
    readr::write_file(wrUIsingle(shiny.title, ganalytics), append = TRUE, file = fname)
    readr::write_file(wrUImain(shiny.prefix, subst, defPtSiz), append = TRUE, file = fname)
    readr::write_file(wrUImarkersDE(shiny.prefix, markers.all, markers.top20, de.genes), append = TRUE, file = fname)
    readr::write_file(wrUIgeneSig(shiny.prefix, gene.ranks), append = TRUE, file = fname)
    readr::write_file(glue::glue(', \n'), append = TRUE, file = fname)
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
        "ggrepel","hdf5r","ggdendro","gridExtra","AUCell","rbokeh")), path = fname)
    readr::write_file(wrSVload(shiny.prefix, markers.all, markers.top20, de.genes), append = TRUE, path = fname)
    readr::write_file(wrSVfix(), append = TRUE, path = fname)
    readr::write_file(wrSVmain(shiny.prefix, subst), append = TRUE, path = fname)
    readr::write_file(wrSVmarkersDE(shiny.prefix, markers.all, markers.top20, de.genes), append = TRUE, path = fname)
    readr::write_file(wrSVgeneSig(shiny.prefix, gene.ranks), append = TRUE, file = fname)
    readr::write_file(wrSVend(), append = TRUE, path = fname)
    
    
    ### Write code for ui.R
    fname = paste0(shiny.dir, "/ui.R")
    readr::write_file(wrLib(
      c("shiny","shinyhelper","data.table","Matrix","DT","magrittr")), path = fname)
    readr::write_file(wrUIload(shiny.prefix), append = TRUE, path = fname)
    readr::write_file(wrUIsingle(shiny.title, ganalytics), append = TRUE, path = fname)
    readr::write_file(wrUImain(shiny.prefix, subst, defPtSiz), append = TRUE, path = fname)
    readr::write_file(wrUImarkersDE(shiny.prefix, markers.all, markers.top20, de.genes), append = TRUE, path = fname)
    readr::write_file(wrUIgeneSig(shiny.prefix, gene.ranks), append = TRUE, file = fname)
    readr::write_file(glue::glue(', \n'), append = TRUE, path = fname)
    readr::write_file(wrUIend(shiny.footnotes), append = TRUE, path = fname)
    
    
    ### Write code for google-analytics.html
    if(!is.na(ganalytics)){
      fname = paste0(shiny.dir, "/google-analytics.html")
      readr::write_file(wrUIga(ganalytics), path = fname)
    }
  }
  
}


