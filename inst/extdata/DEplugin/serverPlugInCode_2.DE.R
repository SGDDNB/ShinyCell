require(plotly)

DEgeneList <- reactive ({
	require(gprofiler2)
	require(plotly)
	req(input$DE_comparison)
	DEdata <- NULL
	bckgeneList <- NULL
	if ((length(input$DE_comparison) == 1) && (length(input$DE_Cell_Type) == 1))
	{
		DEdata <- DE_Matrix()$DEdata

		if (input$UseBackgroundGenes == "Yes")
        { bckgeneList <- DEdata$bckgeneList  }
	}
	
	if (is.null(DEdata))  # check there is data
	{	return(NULL)   }
	
	geneList <- list(c(DE_Matrix()$geneListUp, DE_Matrix()$geneListDown))
	if (input$upDown_regulated == "Up")
    {    geneList <- list(DE_Matrix()$geneListUp)  }
    if (input$upDown_regulated == "down")
    {    geneList <- list(DE_Matrix()$geneListDown)  }
    
	

    names(geneList) <- paste0(input$DE_Cell_Type, " from comparison: ", input$DE_comparison)
	
	
	gostres <- gost(query = geneList, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = bckgeneList, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
				
	
	return(gostres)
})

output$DE_selectCell.ui <- renderUI({
	req(input$DE_comparison)
	if(length(input$DE_comparison) == 1)
	{
		DEdata <- readRDS(paste0("DE/",input$DE_comparison))
		allMarkersFormat <- grep(pattern = ".*?.AllMarkers.RDS", x= input$DE_comparison)
		cellnames <- {}
		if(length(allMarkersFormat))
		{   # columns:  p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
			cellnames <- unique(DEdata$cluster) 
		}
		else # Data will be a list
		{ cellnames <- names(DEdata) }
		
		selectInput("DE_Cell_Type", "cell type:", choices=cellnames)
	}
	
})

DE_Matrix <- reactive({  # PRepares everything DE related. Returns a list.
	
	req(input$DE_comparison)
	DEdata <- data.frame(error="No file selected")
	bckgeneList  <- NULL 
	geneListUp   <- NULL
	geneListDown <- NULL
	
	if ((length(input$DE_comparison) == 1) && (length(input$DE_Cell_Type) == 1))
	{
		DEdata <- readRDS(paste0("DE/",input$DE_comparison))
		
		# Check to see what type of data is being loaded. 
		allMarkersFormat <- grep(pattern = ".*?.AllMarkers.RDS", x= input$DE_comparison)
		if(length(allMarkersFormat)) # Seurat AllMarkers format
		{   # columns:  p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
			DEdata <- DEdata[DEdata$cluster == input$DE_Cell_Type,] 
			bckgeneList <- DEdata$gene			
			
			idx.sigGenes <- which(DEdata$p_val_adj < input$FDR_cutoff & abs(DEdata$avg_log2FC) > input$FoldChange_cutoff)
			DEdata   <- DEdata[idx.sigGenes,]
			idx.down <- which(DEdata$avg_log2FC < 0)
			idx.up   <- which(DEdata$avg_log2FC > 0)
			geneListUp   <- DEdata$gene[idx.up]
			geneListDown <- DEdata$gene[idx.down]
		}
		else # assume MAST format data
		{ 
			DEdata <- DEdata[[input$DE_Cell_Type]]
			bckgeneList <- DEdata$names
			
			idx.sigGenes <- which(DEdata$fdr < input$FDR_cutoff & abs(DEdata$logFC) > input$FoldChange_cutoff)
			DEdata 	  	 <- DEdata[idx.sigGenes,]
			idx.down 	 <- which(DEdata$logFC < 0)
			idx.up   	 <- which(DEdata$logFC > 0)
			geneListUp   <- DEdata$names[idx.up]
			geneListDown <- DEdata$names[idx.down]
		}
	}
	list(DEdata = DEdata, bckgeneList = bckgeneList, geneListUp = geneListUp, geneListDown=geneListDown)
})

output$download_DE_data_button <- renderUI({
    downloadButton("download_DE_data", "Download selected data into CSV file")
  })

output$ShowDEtable.ui  <- renderDataTable ({
	DE_Matrix()$DEdata
})

output$download_DE_data <-   downloadHandler(
    filename = function() { paste0(gsub( pattern= "RDS$", replacement="csv", x= input$DE_comparison,))},
    content = function(file) {
      write.csv(DE_Matrix()$DEdata, file, row.names = FALSE, quote=FALSE)
})

output$gprofilerPlot  <- renderPlotly({
	
	gostres <- DEgeneList()
		
	idx <- which(gostres$result$term_size < input$MaxTermSize)
    gostres$result <- gostres$result[idx,]
	
	
	gostplot(gostres, capped = (input$CappedPlot=="On"), interactive = TRUE)

})

output$download_gprofilerPlot <- downloadHandler(
        filename = function() {
            "interactive_plot.html"  # Name of the file to download
        },
        content = function(file) {
            # Save the Plotly plot as an HTML widget
			require(htmlwidgets)
            gostres <- DEgeneList()

            idx <- which(gostres$result$term_size < input$MaxTermSize)
            gostres$result <- gostres$result[idx,]


            plot <- gostplot(gostres, capped = (input$CappedPlot=="On"), interactive = TRUE)
			saveWidget(as_widget(plot), file)
			
        }
)

output$download_DEtable <- downloadHandler(
        filename = function() {
			fileFormat <- "DE"
		    allMarkersFormat <- grep(pattern = ".*?.AllMarkers.RDS", x= input$DE_comparison)
			if(length(allMarkersFormat)) # Seurat AllMarkers format
			{  fileFormat <- "allMarkers"  }
			filename = paste0(fileFormat,"_table.cell_",input$DE_Cell_Type,".FDR_lt_",input$FDR_cutoff, ".FC_gt_",input$FoldChange_cutoff, ".csv")  
			# Filename will describe contents of table
        },
        content = function(file) {
            # Save the DE table as a CSV file
			write.csv(x= DE_Matrix()$DEdata, file=file)
        }
)

manhattenData <- reactive ({
	req(input$DE_comparison)
	gene_tables <-  readRDS(paste0("DE/",input$DE_comparison))
	condition_names <- names(gene_tables)
	
	names(condition_names) <- condition_names

	# Add condition names to each data frame and combine them
	combined_data <- do.call(rbind, lapply(seq_along(gene_tables), function(i) {
	  df <- gene_tables[[i]]
	  df$Condition <- condition_names[i]  # Add the condition name
	  df$GeneIndex <- seq_len(nrow(df))   # Add a gene index for plotting
	  df$ConditionIndex <- i 
	  df$JitteredCondition <- i + runif(nrow(df), -0.2, 0.2) # for plotly output
	  return(df)
	}))
	combined_data$Log10FDR <- -log10(combined_data$fdr)
	combined_data
	
})
        
output$manhattenPlot  <- renderPlotly({
	require(plotly)
	combined_data <- manhattenData()
	condition_names <- unique(combined_data$Condition)
    names(condition_names) <- condition_names
	
	plot <- plot_ly(
			data = combined_data,
			x = ~JitteredCondition,
			#y = ~Log10FDR,
			y= get(input$manhatten_yAxis),
			type = 'scatter',
			 mode = 'markers',
			 color = ~as.factor(ConditionIndex),   # would be nice to have this the same colour as Dimplots
			 text = ~paste("Gene:", names, "<br>Condition:", Condition),
			 marker = list(size = 6, opacity = 0.6)
			) %>%
			layout(
				xaxis = list(
				  title = "Condition",
				  tickvals = seq_along(condition_names),
				  ticktext = condition_names,
				  tickangle = 45
				),
				yaxis = list(title = "-log10(FDR)"),
				title = "Interactive Manhattan-like Plot of Differential Gene Expression"
			)


	# Print the plot
	plot
	
	
})