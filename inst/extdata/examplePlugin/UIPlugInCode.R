navbarMenu("PLUGIN: Peaks",### Tab1.a1: cellInfo vs geneExpr on dimRed  #DTH added
  tabPanel( 
    HTML("Peak lookup"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 12, h4("Gene lookup"), 
        fluidRow( 
          column( 
            3, selectInput("PeaksGene", "Gene name:", choices=NULL) %>%  
                     helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))), 
						   
						   
				sliderInput("PeaksGeneDistance", "Flanking distance:", 
                             min = 0, max = 100000, value = 5000, step = 1000)
            
        ),
		 column(9, 
	         h4("Sierra peaks"),
			 # NEED OUPUT FUNCTION FOR SIERRA
			 uiOutput("SierraPeaks.ui"),  
	           
	         h4("ATAC peaks"),
			 # NEED OUTPUT FUNCTION FOR ATAC
			 uiOutput("ATACPeaks.ui")
	   )
		) # FluidRow 
	 
     ) # column 12	   
	 ),  # End of Fluid Row
	 )  # End of tab panel
	)