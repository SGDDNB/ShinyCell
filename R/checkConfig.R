#' Checks if shinycell config data.table contains any errors
#'
#' Checks if shinycell config data.table contains any errors. It is useful and 
#' reccomended to run this function if users have motified the shinycell 
#' config manually. Errors can include (i) levels in scConf does not match 
#' that in the Seurat/SingleCellExperiment object, (ii) number of levels does 
#' not match number of colours and (iii) specified colours are invalid colours.
#'
#' @param scConf shinycell config data.table
#' @param obj input single-cell data object. Both Seurat objects (v3+) and 
#'   SingleCellExperiment objects are accepted.
#'
#' @return any potential error messages
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' checkConfig(scConf, seu)
#'
#' @export
checkConfig <- function(scConf, obj){
  nErr = 0
  # Extract metadata to check
  if(class(obj)[1] == "Seurat"){
    objMeta = obj@meta.data
  } else if (class(obj)[1] == "SingleCellExperiment"){
    objMeta = obj@colData
  } else {
    stop("Only Seurat or SingleCellExperiment objects are accepted!")
  }
  
  # Loop through metadata and check
  for(iMeta in scConf$ID){
    fID = strsplit(scConf[ID == iMeta]$fID, "\\|")[[1]]
    if(!is.na(fID[1])){
      # Check if levels in scConf matches objMeta
      if(!all.equal(sort(fID), sort(as.character(unique(objMeta[[iMeta]]))))){
        message(paste0("ERROR! factors in scConf/obj do not match for: ", iMeta))
        nErr = nErr + 1
      }
      
      # Check if no. factors match no. colours
      fUI = strsplit(scConf[ID == iMeta]$fUI, "\\|")[[1]]
      if(length(fID) != length(fUI)){
        message(paste0("ERROR! no. factors & display do not match for: ", iMeta))
        nErr = nErr + 1
      }
      
      # Check if no. factors match no. colours
      fCL = strsplit(scConf[ID == iMeta]$fCL, "\\|")[[1]]
      if(length(fID) != length(fCL)){
        message(paste0("ERROR! no. factors & colours do not match for: ", iMeta))
        nErr = nErr + 1
      }
      
      # Check if colours are valid colours
      res = try(col2rgb(fCL), silent = TRUE)
      if("try-error" %in% class(res)){
        message(paste0("ERROR! invalid colours found for metadata: ", iMeta))
        nErr = nErr + 1
      }
    }
  }
  
  # Final output
  if(nErr == 0){
    message("ALL OK! No errors found in shinycell config!")
  } else {
    message(paste0("Total of ", nErr, " errors found in shinycell config!"))
  }
  # return(scConf)
}


