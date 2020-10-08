#' Shows the order in which metadata will be displayed
#'
#' Shows the order in which metadata will be displayed in the shiny app. This 
#' helps users to decide if the display order is ok. If not, users can use 
#' \code{reorderMeta} to change the order in which metadata will be displayed.
#' 
#' @param scConf shinycell config data.table
#'
#' @return table showing the order in which metadata will be displayed
#'
#' @author John F. Ouyang
#'
#' @import data.table grid gridExtra
#'
#' @examples
#' showOrder(scConf)
#' 
#' @export
showOrder <- function(scConf){
  
  # Start!
  ggOut = scConf[, c("ID", "UI", "fID"), with = FALSE]
  ggOut$nlvl = 0
  ggOut$default = as.character(scConf$default)
  for(i in seq_along(ggOut$ID)){
    ggOut[i]$nlvl = length(strsplit(ggOut[i]$fID, "\\|")[[1]])
  }
  ggOut[is.na(fID)]$nlvl = 0
  ggOut[is.na(fID)]$fID = "cont."
  ggOut[nlvl > 0]$fID = "cat."
  ggOut[default == 0]$default = ""
  colnames(ggOut) = c("actual name", "display name", "type", "nlevels", "default")
  
  # Plot table
  grid.newpage()
  grid.table(ggOut)

  return(ggOut)
}


