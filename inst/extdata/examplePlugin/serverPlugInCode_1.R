    
updateSelectizeInput(session, "PeaksGene", choices = names({scRNA}gene), server = TRUE,  #ew20230301scRNAgene
	selected = {scRNA}def$gene1, options = list(  	 #ew20230301scRNAdef
	maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  
updateSelectizeInput(session, "GenetoATACPeaks", choices = names({scRNA}gene), server = TRUE,
							selected = {scRNA}def$gene1, options = list(  
            				maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))	