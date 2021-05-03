#' Updates shinycell config to recognize a metadata as a discrete one 
#'
#' Updates shinycell config to recognize a metadata as a discrete one. This 
#' function is useful when a discrete metadata only contains integers, e.g.
#' unspervised cluster labels starting from 0 to (n-1) clusters. If these 
#' metadata are not factored, then ShinyCell is unable to recognise these 
#' as discrete metadata automatically. This is especially important when 
#' working with loom files as they do not support factors/levels by default.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.input metadata from the single-cell metadata to be made 
#'   discrete. Must match one of the following:
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
#' scConf = makeMetaCategorical(scConf, "clusterID", loom_obj)
#'
#' @export
makeMetaCategorical <- function(scConf, meta.to.input, obj, maxLevels = 50){
  meta = meta.to.input[1]     # use only the first one
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
  
  # Check if meta exist in obj metadata
  if(!all(meta %in% colnames(objMeta))){
    stop("meta.to.input not found in single-cell data object!")
  }
  if(!all(meta %in% scConf$ID)){
    stop("meta.to.input not found in shinycell config!")
  }
  
  # Start discretizing meta.data
  objMeta[[meta]] = factor(objMeta[[meta]])
  nLevels = nlevels(objMeta[[meta]])
  if(nLevels > maxLevels){
    stop(paste0(meta, 
                " has exceeded the maximum number of levels allowed!"))
  }
  
  if(nLevels >= 2){
    scConf[ID == meta]$fID = paste0(levels(objMeta[[meta]]), collapse = "|")
    scConf[ID == meta]$fUI = scConf[ID == meta]$fID
    scConf[ID == meta]$fCL = paste0(colorRampPalette(brewer.pal(12, "Paired"))(nLevels), 
                         collapse = "|")
    scConf[ID == meta]$fRow = ceiling(nLevels / 4)
    scConf[ID == meta]$grp = TRUE
  } else if(nLevels == 1){
    scConf[ID == meta]$fID = levels(objMeta[[meta]])
    scConf[ID == meta]$fUI = scConf[ID == meta]$fID
    scConf[ID == meta]$fCL = "black"
    scConf[ID == meta]$fRow = 1
  }

  return(scConf)
}


