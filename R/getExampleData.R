#' Download example Seurat objects / single-cell data
#'
#' Download example Seurat objects / single-cell data required for 
#' ShinyCell tutorials.
#'
#' @param type can be either "single" or "multi" or "h5ad" or "loom" 
#'   or "plaintext"
#' 
#' @return downloaded Seurat object
#'
#' @author John F. Ouyang
#'
#' @import data.table
#'
#' @examples
#' getExampleData()
#'
#' @export
getExampleData <- function(type = c("single", "multi", "h5ad", "loom", 
                                    "plaintext")){
  # Setup and checks
  files = c("http://files.ddnetbio.com/hrpiFiles/readySeu_rset.rds",
            "http://files.ddnetbio.com/hrpiFiles/readySeu_d21i.rds",
            "http://files.ddnetbio.com/shinyCell/endocrinogenesis_day15.h5ad",
            "http://files.ddnetbio.com/shinyCell/xr7ne3_dim_reduction_13225_output.loom",
            "http://files.ddnetbio.com/shinyCell/rset.txt.gz")
  names(files) = c("./readySeu_rset.rds", 
                   "./readySeu_d21i.rds",
                   "./endocrinogenesis_day15.h5ad",
                   "./xr7ne3_dim_reduction_13225_output.loom",
                   "./rset.txt.gz")
  
  if(type[1] == "single"){
    files = files[1]
  } else if(type[1] == "multi"){
    files = files[1:2]
  } else if(type[1] == "h5ad"){
    files = files[3]
  } else if(type[1] == "loom"){
    files = files[4]  
  } else if(type[1] == "plaintext"){
    files = files[5]  
  } else {
    stop("argument has to be either 'single' or 'multi' or 'h5ad' or 'loom' or 'plaintext'!")
  }
  
  # Download files
  for(i in seq_along(files)){
    if(!file.exists(names(files)[i])) {
      res <- tryCatch(download.file(files[i],
                                    destfile = names(files)[i],
                                    method = "auto"),
                      error = function(e) 1)
      if(res == 1){ 
        stop(paste0("Following file cannot be downloaded: ", names(files)[i]))
      }
    }
  }
  # return(scConf)
}


