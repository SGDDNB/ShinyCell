#' Add a metadata to be included in the shiny app
#'
#' Add a metadata to be included in the shiny app.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.add metadata to add from the single-cell metadata.  
#'   Must match the column names in 
#'   \code{seu@meta.data} i.e. \code{colnames(seu@meta.data)} for 
#'   Seurat objects or the column names in 
#'   \code{sce@colData} i.e. \code{colnames(sce@colData)} for 
#'   SingleCellExperiment objects. 
#' @param obj original single-cell data object used to create shinycell 
#'   config data.table
#' 
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = addMeta(scConf, c("orig.ident"), seu)
#'
#' @export
addMeta <- function(scConf, meta.to.add, obj){
  # Extract corresponding metadata
  if(class(obj)[1] == "Seurat"){
    objMeta = obj@meta.data
  } else if (class(obj)[1] == "SingleCellExperiment"){
    objMeta = obj@colData
  } else {
    stop("Only Seurat or SingleCellExperiment objects are accepted!")
  }
  
  # Check if meta.to.add exist in obj metadata
  if(!all(meta.to.add %in% colnames(objMeta))){
    stop("meta.to.add not found in single-cell data object!")
  }
  
  # Start adding meta.data
  for(iMeta in meta.to.add){
    tmpConf = data.table(ID = iMeta, UI = iMeta, fID = NA, fUI = NA, 
                         fCL = NA, fRow = NA, default = 0, grp = FALSE)
    
    # Convert to factors if metadata contains characters
    if(is.character(objMeta[[iMeta]])){
      objMeta[[iMeta]] = factor(objMeta[[iMeta]])
    }
    
    # Additional preprocessing for categorical metadata
    nLevels = nlevels(objMeta[[iMeta]])
    if(nLevels >= 2){
      tmpConf$fID = paste0(levels(objMeta[[iMeta]]), collapse = "|")
      tmpConf$fUI = tmpConf$fID
      tmpConf$fCL = paste0(colorRampPalette(brewer.pal(12, "Paired"))(nLevels), 
                           collapse = "|")
      tmpConf$fRow = ceiling(nLevels / 4)
      tmpConf$grp = TRUE
    } else if(nLevels == 1){
      tmpConf$fID = levels(objMeta[[iMeta]])
      tmpConf$fUI = tmpConf$fID
      tmpConf$fCL = "black"
      tmpConf$fRow = 1
    }
    scConf = rbindlist(list(scConf, tmpConf))
  }
  return(scConf)
}


