#' Add a metadata to be included in the shiny app
#'
#' Add a metadata to be included in the shiny app.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.add metadata to add from the single-cell metadata.  
#'   Must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: column names in \code{seu@meta.data} 
#'       i.e. \code{colnames(seu@meta.data)}
#'     \item{SCE objects}: column names in \code{sce@colData} 
#'       i.e. \code{colnames(sce@colData)}
#'     \item{h5ad files}: column names in \code{h5ad.obs} 
#'       i.e. \code{h5ad.obs.columns.values} 
#'     \item{loom files}: column names in \code{loom/col_attrs} 
#'       i.e. \code{loom/col_attrs.names}
#'   }
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment 
#'   data or input file path for h5ad / loom files
#' @param maxLevels maximum number of levels allowed for categorical metadata.
#'   Metadata with nlevels > maxLevels will throw up an error message
#' 
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table reticulate hdf5r
#'
#' @examples
#' scConf = addMeta(scConf, c("orig.ident"), seu)
#'
#' @export
addMeta <- function(scConf, meta.to.add, obj, maxLevels = 50){
  # Extract corresponding metadata
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    objMeta = obj@meta.data
    
  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    objMeta = SingleCellExperiment::colData(obj)
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    ad <- import("anndata", convert = FALSE)
    inpH5 = ad$read_h5ad(obj)
    objMeta = data.frame(py_to_r(inpH5$obs$values))
    rownames(objMeta) = py_to_r(inpH5$obs_names$values)
    colnames(objMeta) = py_to_r(inpH5$obs$columns$values)
    for(i in colnames(objMeta)){
      objMeta[[i]] = unlist(objMeta[[i]])   # unlist and refactor
      if(as.character(inpH5$obs[i]$dtype) == "category"){
        objMeta[[i]] = factor(objMeta[[i]], levels = 
                                py_to_r(inpH5$obs[i]$cat$categories$values))
      }
    } 
    
  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    inpLM = H5File$new(obj, mode = "r+")
    cellIdx = which(inpLM[["col_attrs"]]$names == "CellID")
    if(length(cellIdx) != 1){
      stop("CellID attribute not found in col_attrs in loom file!")
    }
    objMeta = data.frame(row.names = inpLM[["col_attrs"]][["CellID"]]$read())
    for(i in inpLM[["col_attrs"]]$names[-cellIdx]){
      tmp = inpLM[["col_attrs"]][[i]]$read()
      if(length(tmp) == nrow(objMeta)){objMeta[[i]] = tmp}
    }
    inpLM$close_all()
    
  } else {
    stop("Only Seurat/SCE objects or h5ad/loom file paths are accepted!")
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
    if(nLevels > maxLevels){
      stop(paste0(iMeta, " has exceeded the maximum number of levels allowed!"))
    }
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


