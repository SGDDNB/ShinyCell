# For this plugin there must be excel spreadsheets pre-generated. The name of the spreadsheet should identify the comparison 
#    and the sheet names with the document should reflect the cell types.

navbarMenu("PLUGIN: Differential expression",
  
   tabPanel( 
		HTML("DE settings"), 
		h4("Differential Expression (DE) settings"), 
		"In this tab, users can select cell information and gene select thresholds that identify DE genes.", 
		br(),br(), 
		fluidRow( 
		  column( 12, h4("Comparison"), 
			fluidRow( 
			  column( 
				4, selectInput("DE_comparison", "DE filename:", choices=list.files(path="DE",pattern="*.RDS")) %>%  
						 helper(type = "inline", size = "m", fade = TRUE, 
						 title = "DE comparison", 
						 content = c("Select pre-generated differential expression gene list file")), 
					
					uiOutput("DE_selectCell.ui"), # display cell/cluster IDs            						   
					sliderInput("FDR_cutoff", "FDR cutoff:", min = 0, max = 1, value = 0.05, step = 0.01),
					sliderInput("FoldChange_cutoff", "Minimum fold change:", min = 0, max = 20, value = 1, step = 0.5)
					
			),
			 column(8, h4("DE data"),
				dataTableOutput("ShowDEtable.ui"),
				uiOutput("download_DE_data_button") 
		   )
			) # FluidRow 
		 
		 ) # column 12	   
		 ),  # End of Fluid Row
	 ),  # End of tab panel  DE settings
	tabPanel(
		HTML("GO analysis"),
		h4("Gene ontology (GO) settings"), 
		"In this tab, users can select cell information and gene select thresholds that identify DE genes.", 
		br(),br(), 
		fluidRow(
			column(12, h4(""),
				fluidRow( 
				  column( 
					3, 	selectInput("upDown_regulated", "Gene expression direction:", choices=c("Up","Down", "All"), selected ="All") %>%  
							 helper(type = "inline", size = "m", fade = TRUE, 
							 title = "Gene direction", 
							 content = c("Can select genes that are either up or down regulated. These are then assessed for GO analysis.")),
						selectInput("CappedPlot", "Capped p-values (>=16) :", choices=c("On","off"), selected = "On" ),
						selectInput("UseBackgroundGenes", "Use background genes", choices=c("Yes","No"), selected = "No" ),
						sliderInput("MaxTermSize", "term size cutoff:", min = 0, max = 20000, value = 20000, step = 500),
					),
					column(9, h4("GO output"),
							plotly::plotlyOutput("gprofilerPlot"),
							downloadButton("download_gprofilerPlot", "Download Interactive Plot")
					)
				)
			)
		
		
		
		)
		
		
	), # tabpanel    GO analysis
	tabPanel(
		HTML("Manhatten plot"),
		h4("DE manatten plot settings"), 
		"In this tab, users can select cell information and gene select thresholds that identify DE genes.", 
		br(),br(), 
		fluidRow(
			column(12, h4(""),
				fluidRow(
				  column( 
					3, 
					selectInput("DE_comparison_MP", "DE filename:", choices=list.files(path="DE",pattern="*.RDS")) %>%  
						 helper(type = "inline", size = "m", fade = TRUE, 
						 title = "DE comparison", 
						 content = c("Select pre-generated differential expression gene list file")), 				
					
					selectInput("manhatten_yAxis", "Y axis:", c("Log10FDR", "logFC"), selected ="Log10FDR") %>%  
							 helper(type = "inline", size = "m", fade = TRUE, 
							 title = "Plot option", 
							 content = c("Select what is to be displayed on Y axis."))
					),
					column(9, h4("Manhatten plot output"),
							plotly::plotlyOutput("manhattenPlot")
					)
				)
			)
		)
	)  # tabpanel    Manhatten plot analysis


)