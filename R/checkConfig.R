#' Checks if shinycell config data.table contains any errors
#'
#' Checks if shinycell config data.table contains any errors. It is useful and 
#' reccomended to run this function if users have motified the shinycell 
#' config manually. Errors can include (i) levels in scConf does not match 
#' that in the Seurat/SingleCellExperiment object, (ii) number of levels does 
#' not match number of colours and (iii) specified colours are invalid colours.
#'
#' @param scConf shinycell config data.table
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment 
#'   data or input file path for h5ad / loom files
#'
#' @return any potential error messages
#'
#' @author John F. Ouyang
#'
#' @import data.table reticulate hdf5r
#'
#' @examples
#' checkConfig(scConf, seu)
#'
#' @export
checkConfig <- function(scConf, obj){
  nErr = 0
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
    inpLM$close_all()
    
  } else {
    stop("Only Seurat/SCE objects or h5ad/loom file paths are accepted!")
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
      
      # Check if no. factors match no. display
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


