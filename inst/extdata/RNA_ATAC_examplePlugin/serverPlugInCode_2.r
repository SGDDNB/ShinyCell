  
output$ATACPeaks.ui <- renderUI({ 
		genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
		annotationLibrary <- org.Mm.eg.db::org.Mm.eg.db
		txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
		geneName <- input$GenetoATACPeaks # 
		
		a <- try(AnnotationDbi::select(annotationLibrary, keys = geneName, columns=c("ENTREZID", "SYMBOL"),keytype="SYMBOL"),silent=TRUE)
		b <- AnnotationDbi::select(txdb, keys = a$ENTREZID, columns=c('GENEID', 'TXCHROM', 'TXSTART',  'TXEND','TXNAME', 'TXSTRAND'),keytype="GENEID")

		peakCoord <- names({atac}gene)  
		
		if(as.numeric(gregexpr("chr",peakCoord[1])) == -1){ #no chr prefix
			b$TXCHROM <- gsub("chr","",b$TXCHROM)
		}
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

