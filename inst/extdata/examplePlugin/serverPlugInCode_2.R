output$SierraPeaks.ui <- renderUI({  ### Added DTH
     p( paste0(grep(pattern = input$PeaksGene, x = names({sierra}gene),value = TRUE), collapse=", "))
  }) 