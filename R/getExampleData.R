#' Download example Seurat objects
#'
#' Download example Seurat objects required for ShinyCell tutorial.
#'
#' @param type can be either "single" or "multi" which downloads one or two 
#'   Seurat objects respectively for the Quick Start/Detailed Tutorial or 
#'   Multi-dataset Tutorial respectively
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
getExampleData <- function(type = c("single", "multi")){
  # Setup and checks
  files = c("http://files.ddnetbio.com/hrpiFiles/readySeu_rset.rds",
            "http://files.ddnetbio.com/hrpiFiles/readySeu_d21i.rds")
  names(files) = c("./readySeu_rset.rds", "./readySeu_d21i.rds")
  if(type[1] == "single"){
    files = files[1]
  } else if(type[1] != "multi"){
    stop("argument has to be either 'single' or 'multi'!")
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


