library(randomForest)
library(data.table)

shinyServer(function(input, output, session) {
  
  # Loads the Model to memory
  AChE <- file.path("AChE_RF_int.rds")
  BChE <- file.path("BChE_RF_int.rds")
  AChE.RF <- readRDS(AChE)
  BChE.RF <- readRDS(BChE)
  
  # Retrieving the descriptor names from trained model data
  AChE_internal <- readRDS("AChE_SubstructureFingerprinter_internal.rds")
  BChE_internal <- readRDS("BChE_SubstructureFingerprinter_internal.rds")
  
  #AChE.desc.name <- data.frame(names(AChE_internal))
  AChE.desc.name <- data.frame(rownames(AChE.RF$importance))
  #AChE.desc.name <- data.frame(AChE.desc.name[-nrow(AChE.desc.name),])
  AChE.desc.name <- t(AChE.desc.name)
  names(AChE.desc.name) <- as.character(unlist(AChE.desc.name[1,]))
  
  #BChE.desc.name <- data.frame(names(BChE_internal))
  BChE.desc.name <- data.frame(rownames(BChE.RF$importance))
  #BChE.desc.name <- data.frame(BChE.desc.name[-nrow(BChE.desc.name),])
  BChE.desc.name <- t(BChE.desc.name)
  names(BChE.desc.name) <- as.character(unlist(BChE.desc.name[1,]))
  
  observe({
    
    shinyjs::hide("downloadData") # Hide download button before input submission
    if(input$submitbutton>0)
      shinyjs::show("downloadData") # Show download button after input submission
  })
  
  observe({
    COMPOUNDDATA <- ''
    compoundexample <- 'C(Cc1cccc(CCNc2c3CCCCc3nc4ccccc24)n1)Nc5c6CCCCc6nc7ccccc57	CHEMBL521935
CC1=C[C@@H]2C[C@H](C1)c3c(N)c4ccc(Cl)cc4nc3C2	CHEMBL4112162
CN(CCCOc1ccc2C=CC(=O)Oc2c1)Cc3cccc(OC(=O)NCCCCCCCN4CCOCC4)c3	CHEMBL2047529
Fc1ccc(C[n+]2ccc(cc2)C(=O)NCCc3c[nH]c4ccccc34)cc1	CHEMBL3628058
CN(C)Cc1oc(CSCCNc2cc(NN=C(C)C)c(cc2[N+](=O)[O-])[N+](=O)[O-])cc1	CHEMBL106932
COc1ccc2nc3CCCCc3c(NCCNC4=CC(=O)c5ccccc5C4=O)c2c1	CHEMBL3356528
'
    
    if(input$addlink>0) {
      isolate({
        COMPOUNDDATA <- compoundexample
        updateTextInput(session, inputId = "Sequence", value = COMPOUNDDATA)
      })
    }
  })
  
  datasetInput <- reactive({
    
    inFile <- input$file1 
    inTextbox <- input$Sequence
    
    if (is.null(inTextbox)) {
      return("Please insert/upload molecules in SMILES notation")
    } else {
      if (is.null(inFile)) {
        # Read data from text box
        x <- inTextbox
        write.table(x, sep="\t", file = "text.smi", col.names=FALSE, row.names=FALSE, quote=FALSE)
        #x <- read.table("text.smi")
        
        
        # PADEL descriptors for Testing set
        
        #test <- x
        
        try(system("bash PADEL.sh", intern = TRUE, ignore.stderr = TRUE))
        #desc.df <- read.csv("descriptors_output.csv")
        AChE.desc.df <- read.csv("descriptors_output.csv")
        AChE.desc.df2 <- AChE.desc.df[,( names(AChE.desc.df) %in% AChE.desc.name )]
        AChE.mol.desc = data.frame(AChE.desc.df2)
        
        BChE.desc.df <- read.csv("descriptors_output.csv")
        BChE.desc.df2 <- BChE.desc.df[,( names(BChE.desc.df) %in% BChE.desc.name )]
        BChE.mol.desc = data.frame(BChE.desc.df2)
   
        # Predicting unknown sequences
        AChE.prediction <- data.frame(Prediction= predict(AChE.RF,AChE.mol.desc), round(predict(AChE.RF,AChE.mol.desc,type="prob"),3))
        BChE.prediction <- data.frame(Prediction= predict(BChE.RF,BChE.mol.desc), round(predict(BChE.RF,BChE.mol.desc,type="prob"),3))
        compoundname <- data.frame(AChE.desc.df$Name)
        row.names(compoundname) <- AChE.desc.df$Name
        results <- cbind(compoundname, AChE.prediction, BChE.prediction)
        #names(results)[1] <- "Name"
        names(results) <- c("Name","AChE.prediction","AChE.active","AChE.inactive","BChE.prediction","BChE.active","BChE.inactive")
        results <- data.frame(results, row.names=NULL)
        
        print(results)
      } 
      else {  
        # Read data from uploaded file
        x <- read.table(inFile$datapath)
        
        # PADEL descriptors for Testing set
        
        test <- x
        
        try(system("bash PADEL.sh", intern = TRUE, ignore.stderr = TRUE))
        #desc.df <- read.csv("descriptors_output.csv")
        AChE.desc.df <- read.csv("descriptors_output.csv")
        AChE.desc.df2 <- desc.df[,( names(desc.df) %in% desc.name )]
        
        AChE.mol.desc = data.frame(desc.df2)
        
        # Predicting unknown sequences
        AChE.prediction <- data.frame(Prediction= predict(RFalpha,AChE.mol.desc), round(predict(RFalpha,AChE.mol.desc,type="prob"),3))
        BChE.prediction <- data.frame(Prediction= predict(RFbeta,BChE.mol.desc), round(predict(RFbeta,BChE.mol.desc,type="prob"),3))
        compoundname <- data.frame(AChE.mol.desc$Name)
        row.names(compoundname) <- AChE.mol.desc$Name
        results <- cbind(compoundname, AChE.prediction, BChE.prediction)
        #names(results)[1] <- "Name"
        names(results) <- c("Name","AChE.prediction","AChE.a.active","AChE.inactive","BChE.prediction","BChE.active","BChE.inactive")
        results <- data.frame(results, row.names=NULL)
        
        print(results)
      }
    }
  })
  
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } else {
      return("Server is ready for prediction.")
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste('predicted_results', '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file, row.names=FALSE)
    })
  
})