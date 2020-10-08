#' Reorder the order in which metadata appear in the shiny app
#'
#' Reorder the order in which metadata appear in the dropdown menu in the 
#' shiny app.
#'
#' @param scConf shinycell config data.table
#' @param new.meta.order character vector containing new order. All metadata 
#'   names must be included, which can be found at \code{scConf$ID}
#' 
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = reorderMeta(scConf, scConf$ID[c(1,3,2,4:length(scConf$ID))])
#'
#' @export
reorderMeta <- function(scConf, new.meta.order){
  # Check if new.meta.order matches scConf$ID
  if(!all.equal(sort(new.meta.order), sort(as.character(scConf$ID)))){
    stop("new.meta.order does not match scConf$ID!")
  }
  
  # Start reordering
  scConf$ID = factor(scConf$ID, levels = new.meta.order)
  scConf = scConf[order(ID)]

  return(scConf)
}


