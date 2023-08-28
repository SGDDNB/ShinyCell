output$SierraPeaks.ui <- renderUI({  ### Added DTH
     p( paste0(grep(pattern = input$PeaksGene, x = names({sierra}gene),value = TRUE), collapse=", "))
  }) 
  
 output$SierraGeneModel.ui <- renderUI({
	   	genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
		annotationLibrary <- org.Mm.eg.db::org.Mm.eg.db
		txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
		geneName <- input$GenetoATACPeaks

		a <- try(AnnotationDbi::select(annotationLibrary, keys = geneName, columns=c("ENTREZID", "SYMBOL"),keytype="SYMBOL"),silent=TRUE)
		b <- AnnotationDbi::select(txdb, keys = a$ENTREZID, columns=c('GENEID', 'TXCHROM', 'EXONSTART',  'EXONEND','TXNAME', 'EXONSTRAND'),keytype="GENEID")
		

		transcript.gr <- GenomicRanges::makeGRangesListFromDataFrame(df = b[which(b$TXCHROM == 'chrX'),],
                                           split.field = "TXNAME",seqnames.field="TXCHROM", start.field="EXONSTART",end.field="EXONEND", strand.field="EXONSTRAND")

		Sierra::PlotCoverage(transcript.gr, geneSymbol=geneName,genome='mm10')
	}) 
  
output$ATACPeaks.ui <- renderUI({ 
		genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
		annotationLibrary <- org.Mm.eg.db::org.Mm.eg.db
		txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
		geneName <- input$GenetoATACPeaks # 

		a <- try(AnnotationDbi::select(annotationLibrary, keys = geneName, columns=c("ENTREZID", "SYMBOL"),keytype="SYMBOL"),silent=TRUE)
		b <- AnnotationDbi::select(txdb, keys = a$ENTREZID, columns=c('GENEID', 'TXCHROM', 'TXSTART',  'TXEND','TXNAME', 'TXSTRAND'),keytype="GENEID")

		peakCoord <- names({atac}gene)  
		peakCoord <- do.call(rbind,strsplit(peakCoord,split = "-"))
		idx <- {}
		blurb <- HTML(paste0("Total number of transcripts: ",length(b$TXSTART),br()))
		
		if (input$PeaktoGeneFeature == "TSS")  # capture peaks within range of TSS 
		{
			toQuery <- b$TXSTART 
			if (b$TXSTRAND[1] == "-")
			{	toQuery <- b$TXEND  }
				
			for(i in 1:length(toQuery))
			{    
				minCoord <- toQuery - input$PeaksGeneDistance
				maxCoord <- toQuery  + input$PeaksGeneDistance
				idx <- c(which(peakCoord[,1] == b$TXCHROM[1] & peakCoord[,2] > minCoord & peakCoord[,2] < maxCoord), idx)
			}
			#browser()
			
		}
		else  # Gene body
		{
			minCoord <- min(b$TXSTART,b$TXEND) - input$PeaksGeneDistance
			maxCoord <- max(b$TXSTART,b$TXEND) + input$PeaksGeneDistance

			idx <- which(peakCoord[,1] == b$TXCHROM[1] & peakCoord[,2] > minCoord & peakCoord[,2] < maxCoord)
		}
		
		
		if(length(idx) > 0)
		{  idx <- unique(idx)
		   blurb_peaks <- paste0(peakCoord[idx,1],"-", peakCoord[idx,2],"-",peakCoord[idx,3],collapse=", ")
		   blurb <- HTML(paste0(blurb,br(),"Total number of peaks: ",length(idx),br(),blurb_peaks))
		}
		else
		{  blurb <- 'No peaks in this gene'  }

        p(blurb)
  }) 


################# CellPhoneDB functions
   output$cellPhoneDB_CellList.ui <- renderUI({  
    sub = strsplit({scrna}conf[UI == input${scrna}d2sub1]$fID, "\\|")[[1]]  
    #checkboxGroupInput("cellPhoneDB_metaData", "Select which cells to show", inline = TRUE,  
	radioButtons("cellPhoneDB_metaData", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub[1])  # Select only first entry 
  })  
 

	output$cellPhoneDBTopCandidate.ui <- renderUI({
		cellPhone_cellID <- input$cellPhoneDB_metaData
	
		blurb <- "No CellPhoneDB directory"
		if (file.exists("cellPhoneDB"))
		{	ld <- list.dirs(path="cellPhoneDB",full.names = FALSE)[-1]
			if (length(ld) == 0)
			{	blurb <- "No data in CellPhoneDB directory"   }
			else
			{	# Loop through each directory and prepare plot
				outputPlot <- list()
				blurb <- ''
				for(i in 1:length(ld))
				{	# Read in data and find best candidates
					fn.mean   <- list.files(path=paste0("cellPhoneDB/",ld[i],"/"), pattern = ".*?means.*?",recursive = FALSE,full.names =  TRUE)
					fn.pvalue <- list.files(path=paste0("cellPhoneDB/",ld[i],"/"),pattern = ".*?pvalues.*?",recursive = TRUE,full.names =  TRUE)
				
					data.means <- read.table(file = fn.mean,header = TRUE, sep = "\t")
					data.pvalue <- read.table(file = fn.pvalue,header = TRUE, sep = "\t")
					# Filter out columns that do not contain cell type
					idx.cellID <- grep(pattern=cellPhone_cellID,x=colnames(data.means))
					#Identify which rows meet criteria
					meanCriteria   <- apply(data.means[,idx.cellID],MARGIN=1, FUN=function(x){  any(x>input$cellPhoneDB_ExpressionFilter)})
					pValueCriteria <- apply(data.pvalue[,idx.cellID],MARGIN=1, FUN=function(x){  any(x<0.05)})
					
					
					blurb <- c(blurb,paste(br(),br(),ld[i],":"),data.means$interacting_pair[meanCriteria & pValueCriteria])
					
				}
			}
		}
		HTML(blurb)
	})

 
 
  output$cellPhoneDB.ui <- renderUI({
	plotOutput("cellPhoneDB",height = paste0(input$cellPhoneDB_plotHeights,"px"))# "2000px")

  })
  
  
  output$cellPhoneDB <- renderPlot({ # Prepare bubbleplot for CellphoneDB
  
	text <- "CellPhoneDB directory does not exist"
	outputPlot <- ggplot() +  annotate("text", x = 4, y = 25, size=8, label = text) + theme_void()	
	cellPhone_geneList <- strsplit(input$cellphoneDB_GeneList, ",")[[1]]
	cellPhone_cellID <- input$cellPhoneDB_metaData
	if (file.exists("cellPhoneDB"))
	{ 	ld <- list.dirs(path="cellPhoneDB",full.names = FALSE)[-1]
		if (length(ld) == 0)
		{	text <- "No data in CellPhoneDB directory"
			outputPlot <- ggplot() +  annotate("text", x = 4, y = 25, size=8, label = text) + theme_void()
		}
		else
		{	# Loop through each directory and prepare plot
			outputPlot <- list()
			for(i in 1:length(ld))
			{
				fn.mean   <- list.files(path=paste0("cellPhoneDB/",ld[i],"/"), pattern = ".*?means.*?",recursive = FALSE,full.names =  TRUE)
				fn.pvalue <- list.files(path=paste0("cellPhoneDB/",ld[i],"/"),pattern = ".*?pvalues.*?",recursive = TRUE,full.names =  TRUE)
				
				data.means <- read.table(file = fn.mean,header = TRUE, sep = "\t")
				data.pvalue <- read.table(file = fn.pvalue,header = TRUE, sep = "\t")
				
				ggData = data.frame()
				idx <- {}
				for(j in cellPhone_geneList)
				{  idx <- c(grep(pattern = j,x = data.means$interacting_pair,ignore.case = TRUE),idx)
				}
				idx.columns <- 12:ncol(data.means)
				ggData <- cbind(cellType = colnames(data.means)[idx.columns])

				# Prepare means data frame (long format)
				toCopy <- t(data.means[idx,idx.columns])
				colnames(toCopy) <- data.means$interacting_pair[idx]
				ggData.Means <- data.frame(cbind(ggData,toCopy))
				ggData.Means <-reshape2::melt(data = ggData.Means, id='cellType')
				colnames(ggData.Means)[3] <- 'means'

				# Prepare pvalue data frame (long format)
				toCopy <- t(data.pvalue[idx,idx.columns])
				colnames(toCopy) <- data.pvalue$interacting_pair[idx]
				ggData.pvalue <- data.frame(cbind(ggData,toCopy))
				ggData.pvalue <-reshape2::melt(data = ggData.pvalue, id='cellType')
				colnames(ggData.pvalue)[3] <- 'pvalue'
                                              
				# make final data.frame for plotting
				ggData <- ggData.Means
				ggData$means <- as.numeric(ggData$means)
				ggData$pvalue <- as.numeric(ggData.pvalue$pvalue)

			#	browser()
				idx <- {}
				for(j in cellPhone_cellID)
				{ idx <- c(grep(pattern = j,x = ggData$cellType),idx)
				}
				ggData<- ggData[idx,] # Subset data for selected cells

				outputPlot[[i]] <- ggplot(ggData, aes(cellType, variable, color = pvalue, size = means))  + ggtitle(ld[i]) +
						geom_point() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12)) +
						scale_size_continuous("Avg expression", range = c(0, 8), 
                             limits = c(0, max(ggData$means))) + #, breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
							scale_color_gradientn("pvalue", limits = c(0,1), colours = cList[["Blue-Yellow-Red"]]) + 
							theme(axis.title = element_blank(), legend.box = "vertical")

				
			} # for
				
		} #if
	} # if (file.exists("cellPhoneDB"))
	
	return(gridExtra::grid.arrange(outputPlot[[1]],outputPlot[[2]],outputPlot[[3]],outputPlot[[4]], ncol=1))
	#return(outputPlot)
	

  })
  
  

  observeEvent(input$cellPhoneDBsub1all, {  
    sub = strsplit({scrna}conf[UI == input$cellPhoneDBsub1all]$fID, "\\|")[[1]]  
    updateCheckboxGroupInput(session, inputId = "cellPhoneDB_metaData", label = "Select which cells to show",  
                             choices = sub, selected = sub, inline = TRUE)  
  })  

  observeEvent(input$cellPhoneDBsub1non, {  
    sub = strsplit({scrna}conf[UI == input$cellPhoneDBsub1non]$fID, "\\|")[[1]]  
    updateCheckboxGroupInput(session, inputId = "cellPhoneDB_metaData", label = "Select which cells to show",  
                             choices = sub, selected = NULL, inline = TRUE)  
  })  
