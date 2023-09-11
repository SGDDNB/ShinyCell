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
	
	if("sc1" %in% shiny.prefixSet){
		scRNA <- "sc1"
	}

	code <- readLines(plugInFile,warn = F)
	code <- paste0(code,collapse="\n")

	c <- try(gsub("\\{sierra\\}",sierra,code,ignore.case=TRUE))
	if(!identical(class(c),"try-error")){
		code <- c
	}
	c <- try(gsub("\\{expression\\}",expression,code,ignore.case=TRUE))
	if(!identical(class(c),"try-error")){
		code <- c
	}
	c <- try(gsub("\\{scRNA\\}",scRNA,code,ignore.case=TRUE))
	if(!identical(class(c),"try-error")){
		code <- c
	}
	c <- try(gsub("\\{atac\\}",atac,code,ignore.case=TRUE))
	if(!identical(class(c),"try-error")){
		code <- c
	}

	code <- paste0(code,"\n")
	#code <- glue::glue(code)

	return (code)
}