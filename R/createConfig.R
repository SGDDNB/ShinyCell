#' Create a shinycell config data.table
#'
#' Create a shinycell config data.table containing (i) the single-cell 
#' metadata to display on the Shiny app, (ii) ordering of factors / 
#' categories of categorical metadata and (iii) colour palettes associated 
#' with each metadata.
#'
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment 
#'   data or input file path for h5ad / loom files
#' @param meta.to.include columns to include from the single-cell metadata. 
#'   Default is \code{NA}, which is to use all columns. Users can specify 
#'   the columns to include, which must match one of the following:
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
#' @param legendCols maximum number of columns allowed when displaying the 
#'   legends of categorical metadata
#' @param maxLevels maximum number of levels allowed for categorical metadata.
#'   Metadata with nlevels > maxLevels will be discarded automatically
#'
#' @return shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table reticulate hdf5r
#'
#' @examples
#' scConf = createConfig(obj)
#'
#' @export
createConfig <- function(obj, meta.to.include = NA, legendCols = 4,
                         maxLevels = 50){
  # Extract corresponding metadata
  drExist = TRUE
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    objMeta = obj@meta.data
    if(length(names(obj@reductions)) == 0){drExist = FALSE}
    
  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    objMeta = SingleCellExperiment::colData(obj)
    if(length(SingleCellExperiment::reducedDimNames(obj)) == 0){drExist=FALSE}
    
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
    if(length(py_to_r(inpH5$obsm_keys())) == 0){drExist = FALSE}
    
  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    inpLM = hdf5r::H5File$new(obj, mode = "r+")
    cellIdx = which(inpLM[["col_attrs"]]$names == "CellID")
    if(length(cellIdx) != 1){
      stop("CellID attribute not found in col_attrs in loom file!")
    }
    objMeta = data.frame(row.names = inpLM[["col_attrs"]][["CellID"]]$read())
    for(i in inpLM[["col_attrs"]]$names[-cellIdx]){
      tmp = inpLM[["col_attrs"]][[i]]$read()
      if(length(tmp) == nrow(objMeta)){objMeta[[i]] = tmp}
    }
    nDR = inpLM[["col_attrs"]]$names[
      grep("pca|tsne|umap", inpLM[["col_attrs"]]$names, ignore.case = TRUE)]
    if(length(nDR) == 0){drExist = FALSE}
    inpLM$close_all()

  } else {
    stop("Only Seurat/SCE objects or h5ad/loom file paths are accepted!")
  }
  if(!drExist){
    stop(paste0("ShinyCell did not detect any dimension reduction data \n", 
                "       e.g. umap / tsne. Has any analysis been performed?"))
  }
  
  # Checks and get list of metadata to include
  if(is.na(meta.to.include[1])){meta.to.include = colnames(objMeta)}
  if(length(meta.to.include) < 2){stop("At least 2 metadata is required!")}
  
  # Start making config data.table
  scConf = data.table()
  for(iMeta in meta.to.include){
    tmpConf = data.table(ID = iMeta, UI = iMeta, fID = NA, fUI = NA, 
                         fCL = NA, fRow = NA, default = 0, grp = FALSE)
    
    # Convert to factors if metadata contains characters
    if(is.character(objMeta[[iMeta]])){
      objMeta[[iMeta]] = factor(objMeta[[iMeta]])
    }
    
    # Additional preprocessing for categorical metadata
    nLevels = nlevels(objMeta[[iMeta]])
    if(nLevels <= maxLevels){
      if(nLevels >= 2){
        tmpConf$fID = paste0(levels(objMeta[[iMeta]]), collapse = "|")
        tmpConf$fUI = tmpConf$fID
        tmpConf$fCL = paste0(colorRampPalette(brewer.pal(12, "Paired"))(nLevels), 
                             collapse = "|")
        tmpConf$fRow = ceiling(nLevels / legendCols)
        tmpConf$grp = TRUE
      } else if(nLevels == 1){
        tmpConf$fID = levels(objMeta[[iMeta]])
        tmpConf$fUI = tmpConf$fID
        tmpConf$fCL = "black"
        tmpConf$fRow = 1
      }
      scConf = rbindlist(list(scConf, tmpConf))
    }
  }
  
  # Set defaults
  def1 = grep("ident|library", scConf$ID, ignore.case = TRUE)[1]
  def2 = grep("clust", scConf$ID, ignore.case = TRUE)
  def2 = setdiff(def2, def1)[1]
  if(is.na(def1)){def1 = setdiff(c(1,2), def2)[1]}
  if(is.na(def2)){def2 = setdiff(c(1,2), def1)[1]}
  scConf[def1]$default = 1
  scConf[def2]$default = 2
  
  # STOP if there is no single multi-level covariate
  if(nrow(scConf[grp == TRUE]) == 0){
    stop(paste0("ShinyCell did not detect any multi-group cell metadata \n", 
                "       e.g. library / cluster. Has any analysis been performed?"))
  }
  
  return(scConf)
}


