#' Write code for loading libraries
#'
#' @param plugInFile filename containing plugin source code
#' @param shiny.prefixSet character vector containing prefix identities. 
#'
#' @rdname writerPlugin
#' @export writerPlugin
#'


writerPlugin <- function(plugInFile,shiny.prefixSet){
	sierra <- grep("sierra",shiny.prefixSet,value=TRUE,ignore.case=TRUE)
	expression <- grep("expression",shiny.prefixSet,value=TRUE,ignore.case=TRUE)
	scRNA <- grep("scRNA",shiny.prefixSet,value=TRUE,ignore.case=TRUE)
	atac <- grep("atac",shiny.prefixSet,value=TRUE,ignore.case=TRUE)
	
	
	
	code <- readLines(plugInFile,warn = F)
	code <- paste0(code,collapse="\n")
	
	code <- gsub("\\{sierra\\}",sierra,code,ignore.case=TRUE)
	code <- gsub("\\{expression\\}",expression,code,ignore.case=TRUE)
	code <- gsub("\\{scRNA\\}",scRNA,code,ignore.case=TRUE)
	code <- gsub("\\{atac\\}",atac,code,ignore.case=TRUE)
	
	code <- paste0(code,"\n")
	#code <- glue::glue(code)
	
	return (code)
}

