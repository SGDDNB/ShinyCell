#' Modify the display name of metadata
#'
#' Modify the display name of metadata. It is possible that the original 
#' metadata name is not so informative e.g "orig.ident" or too long e.g. 
#' "seurat_clusters" and users want to shorten the way they are displayed 
#' on the shiny app. This function allows users to specify display names for 
#' metadata i.e. the names that will be displayed on the shiny app. Note that 
#' \code{showLegend} shows the display name instead of the actual name.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.mod metadata for which to modify the display name. Users can 
#'   either use the actual metadata column names or display names. Multiple 
#'   metadata can be specified. It is reccomended to use the original metadata 
#'   column names to reduce confusion.
#' @param new.name new display names for the corresponding metadata
#'
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = modMetaName(scConf, 
#'                      meta.to.mod = c("orig.ident", "seurat_clusters"), 
#'                      new.name = c("library", "cluster"))
#'
#' @export
modMetaName <- function(scConf, meta.to.mod, new.name){
  # Check if meta.to.mod exist
  if(all(meta.to.mod %in% scConf$ID)){
    useID = TRUE   # Use IDs
  } else if(all(meta.to.mod %in% scConf$UI)){
    useID = FALSE  # Use UIs
  } else {
    stop("meta.to.mod not found in shinycell config!")
  }
  
  # Check replacement length
  if(length(meta.to.mod) != length(new.name)){
    stop("Lengths of meta.to.mod and new.name do not match!")
  }
  
  # Start changing display name
  if(useID){
    for(i in seq_along(meta.to.mod)){
      scConf[ID == meta.to.mod[i]]$UI = new.name[i]
    }
  } else {
    for(i in seq_along(meta.to.mod)){
      scConf[UI == meta.to.mod[i]]$UI = new.name[i]
    }
  }
  
  return(scConf)
}


