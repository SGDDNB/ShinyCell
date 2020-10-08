#' Set the default metadata to display
#'
#' Set the default metadata to display in the shiny app. Default1 is used when 
#' plotting metadata with gene expression or when plotting two metadata 
#' simultaneously. Default2 is used when plotting two metadata simultaneously.
#'
#' @param scConf shinycell config data.table
#' @param default1 metadata to set as the 1st default metadata to display. 
#'   Users can either use the actual metadata column names or display names
#' @param default2 metadata to set as the 2nd default metadata to display. 
#'   Users can either use the actual metadata column names or display names
#' 
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = modDefault(scConf, default1 = "library", default2 = "cluster")
#'
#' @export
modDefault <- function(scConf, default1, default2){
  # Check that only one metadata is provided
  if(length(default1) != 1){stop("Please specify only one default1 metadata!")}
  if(length(default2) != 1){stop("Please specify only one default2 metadata!")}
  
  # Check if default1/default2 exist
  if(default1 %in% scConf$ID){
    useID1 = TRUE   # Use IDs
  } else if(default1 %in% scConf$UI){
    useID1 = FALSE  # Use UIs
  } else {
    stop("default1 not found in shinycell config!")
  }
  if(default2 %in% scConf$ID){
    useID2 = TRUE   # Use IDs
  } else if(default2 %in% scConf$UI){
    useID2 = FALSE  # Use UIs
  } else {
    stop("default2 not found in shinycell config!")
  }
  
  # Start changing the defaults
  scConf$default = 0
  if(useID1){
    scConf[ID == default1]$default = 1
  } else {
    scConf[UI == default1]$default = 1
  }
  if(useID2){
    scConf[ID == default2]$default = 2
  } else {
    scConf[UI == default2]$default = 2
  }
  return(scConf)
}


