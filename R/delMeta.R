#' Remove a metadata from being included in the shiny app
#'
#' Remove a metadata from being included in the shiny app.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.del metadata to delete. Users can either use the original 
#'   metadata column names or display names. For more information regarding 
#'   display name, see \code{?modMetaName}. Multiple metadata can be specified.
#'
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = delMeta(scConf, c("orig.ident"))
#'
#' @export
delMeta <- function(scConf, meta.to.del){
  # Check if meta.to.del exist
  if(all(meta.to.del %in% scConf$ID)){
    useID = TRUE   # Use IDs
  } else if(all(meta.to.del %in% scConf$UI)){
    useID = FALSE  # Use UIs
  } else {
    stop("meta.to.del not found in shinycell config!")
  }
  
  # Start removing meta.data
  if(useID){
    scConf = scConf[!ID %in% meta.to.del]
  } else {
    scConf = scConf[!UI %in% meta.to.del]
  }
  
  # Reassign default if it is removed
  if(!1 %in% scConf$default){
    chkname = paste0(scConf[default != 2]$ID, "_", scConf[default != 2]$UI)
    def1 = scConf[default != 2]$ID[grep("ident|library", chkname, 
                                        ignore.case = TRUE)[1]]
    if(is.na(def1)){def1 = scConf[default != 2]$ID[1]}
    scConf[ID == def1]$default = 1
  }
  if(!2 %in% scConf$default){
    chkname = paste0(scConf[default != 1]$ID, "_", scConf[default != 1]$UI)
    def2 = scConf[default != 1]$ID[grep("clust", chkname, 
                                        ignore.case = TRUE)[1]]
    if(is.na(def2)){def2 = scConf[default != 1]$ID[1]}
    scConf[ID == def2]$default = 2
  }
  
  return(scConf)
}


