#' Create a shinycell config data.table
#'
#' Create a shinycell config data.table containing (i) the single-cell 
#' metadata to display on the Shiny app, (ii) ordering of factors / 
#' categories of categorical metadata and (iii) colour palettes associated 
#' with each metadata.
#'
#' @param obj input single-cell data object. Both Seurat objects (v3+) and 
#'   SingleCellExperiment objects are accepted.
#' @param meta.to.include columns to include from the single-cell metadata. 
#'   Default is \code{NA}, which is to use all columns. Users can specify 
#'   columns to include, which must match the column names in 
#'   \code{seu@meta.data} i.e. \code{colnames(seu@meta.data)} for 
#'   Seurat objects or the column names in 
#'   \code{sce@colData} i.e. \code{colnames(sce@colData)} for 
#'   SingleCellExperiment objects. 
#' @param legendCols maximum number of columns allowed when displaying the 
#'   legends of categorical metadata
#'
#' @return shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' scConf = createConfig(obj)
#'
#' @export
createConfig <- function(obj, meta.to.include = NA, legendCols = 3){
  # Extract corresponding metadata
  if(class(obj)[1] == "Seurat"){
    objMeta = obj@meta.data
  } else if (class(obj)[1] == "SingleCellExperiment"){
    objMeta = obj@colData
  } else {
    stop("Only Seurat or SingleCellExperiment objects are accepted!")
  }
  
  # Checks and get list of metadata to include
  if(is.na(meta.to.include)){meta.to.include = colnames(objMeta)}
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
  
  # Set defaults
  def1 = grep("ident|library", meta.to.include, ignore.case = TRUE)[1]
  def2 = grep("clust", meta.to.include, ignore.case = TRUE)
  def2 = setdiff(def2, def1)[1]
  if(is.na(def1)){def1 = setdiff(c(1,2), def2)[1]}
  if(is.na(def2)){def2 = setdiff(c(1,2), def1)[1]}
  scConf[def1]$default = 1
  scConf[def2]$default = 2
  
  return(scConf)
}


