#' Modify the legend labels for categorical metadata
#'
#' Modify the legend labels for categorical metadata.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.mod metadata for which to modify the legend labels. Users 
#'   can either use the actual metadata column names or display names. Please 
#'   specify only one metadata
#' @param new.labels character vector of new legend labels
#' 
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = modLabels(scConf, meta.to.mod = "library", 
#'                    new.labels = c("Fib", "Primed", "Naive", "RSeT"))
#'
#' @export
modLabels <- function(scConf, meta.to.mod, new.labels){
  # Check that only one metadata is provided
  if(length(meta.to.mod) != 1){
    stop("Please specify only one metadata to modify legend labels!")
  }
  
  # Check if meta.to.mod exist
  if(meta.to.mod %in% scConf$ID){
    useID = TRUE   # Use IDs
  } else if(meta.to.mod %in% scConf$UI){
    useID = FALSE  # Use UIs
  } else {
    stop("meta.to.mod not found in shinycell config!")
  }
  
  # Check if meta.to.mod is categorical and if length(new.labels) matches 
  if(useID){
    res = strsplit(scConf[ID == meta.to.mod]$fUI, "\\|")[[1]]
  } else {
    res = strsplit(scConf[UI == meta.to.mod]$fUI, "\\|")[[1]]
  }
  if(is.na(res[1])){
    stop("meta.to.mod is not a categorical metadata!")
  }
  if(length(res) != length(new.labels)){
    stop("Length of new.labels does not match!")
  }
  
  # Start changing the colours
  if(useID){
    scConf[ID == meta.to.mod]$fUI = paste0(new.labels, collapse = "|")
  } else {
    scConf[UI == meta.to.mod]$fUI = paste0(new.labels, collapse = "|")
  }
  return(scConf)
}


