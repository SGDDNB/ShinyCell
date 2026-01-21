navbarMenu("PLUGIN: Peaks",
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
        ),
		 column(9, 
	         h4("Sierra peaks"),
			 # NEED OUPUT FUNCTION FOR SIERRA
			 uiOutput("SierraPeaks.ui"), 
			 uiOutput("SierraGeneModel.ui")
	           
	   )
		) # FluidRow 
	 
     ) # column 12	   
	 ),  # End of Fluid Row
	 ),  # End of tab panel
   tabPanel( 
    HTML("ATAC Peak lookup"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 12, h4("Gene lookup"), 
        fluidRow( 
          column( 
            3, selectInput("GenetoATACPeaks", "Gene name:", choices=NULL) %>%  
                     helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))), 
						   
						   
				sliderInput("PeaksGeneDistance", "Flanking distance:", 
                             min = 0, max = 100000, value = 5000, step = 1000),
				radioButtons("PeaktoGeneFeature", "Font size:", 
	                       choices = c("Gene", "TSS"), 
	                       selected = "Gene", inline = TRUE) 
            
        ),
		 column(9, 
 
	           
	         h4("ATAC peaks"),
			 # NEED OUTPUT FUNCTION FOR ATAC
			 uiOutput("ATACPeaks.ui")
	   )
		) # FluidRow 
	 
     ) # column 12	   
	 ),  # End of Fluid Row
	 ),  # End of tab panel

###### CellPhoneDB
	 tabPanel( 
    HTML("CellphoneDB"), 
    h4("Exploring cell - cell ligand receptor interaction"), 
    "In this tab, users can access prebuilt cellphone DB analysis ",  
    br(),br(),
    fluidRow( 
      column( 12, h4("Gene lookup"),
        fluidRow( 
          column( 
            3,textAreaInput("cellphoneDB_GeneList", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0({scrna}def$gene1, collapse = ", ")) %>% 
					helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")),  
						   			   
				
				radioButtons("cellphoneDB_addunderscore", "Stringent name search:", 
	                       choices = c("Yes", "No"), 
	                       selected = "Yes", inline = TRUE),
						   
			    selectInput("cellPhoneDB_CellList", "Cell information to subset:", 
	                      choices = {scrna}conf[grp == TRUE]$UI, 
	                      selected = {scrna}def$grp1, width='100%'), 
	            uiOutput("cellPhoneDB_CellList.ui"),  
	            actionButton("cellPhoneDBsub1all", "Select all groups", class = "btn btn-primary"), 
	            actionButton("cellPhoneDBsub1non", "Deselect all groups", class = "btn btn-primary"), 
				sliderInput("cellPhoneDB_plotHeights", "Plot heights:", 
                             min = 1000, max = 10000, value = 2000, step = 200),br(),
				h4("Top interacting candidates for selected cell type (pvalue < 0.05)"),
				sliderInput("cellPhoneDB_ExpressionFilter", "Avg expression filter (greater than):", 
                             min = 1, max = 100, value = 20, step = 1),
				uiOutput("cellPhoneDBTopCandidate.ui")
            
        ),
		 column(9, 
     
	         h4("Cell phone DB output"),
		
			 uiOutput("cellPhoneDB.ui")
	   )
		) # FluidRow 
	 
     ) # column 12	   
	 ),  # End of Fluid Row
	 )  # End of tab panel
)