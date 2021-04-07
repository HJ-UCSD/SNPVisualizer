#handles backend work of app, specifically the translations
server = (function(input, output,session) {
  output$Peptides_Output<-renderDataTable({
    
    eligible_Peptides<-MDP_database %>% #filtering based on user inputted SNP choice
      filter(.$Peptide_End >= as.numeric(input$SNP_Location)) %>%
      filter(.$Peptide_Start <= as.numeric(input$SNP_Location))
    
    #searching for RS database for value
    RS_Output<-(RS_Number[RS_Number$mtsnp==as.numeric(input$SNP_Location),])$rs
    if(length(RS_Output) == 0)
    {
      RS_Output<-NA
    }
    
    for(i in 1:nrow(eligible_Peptides))
    {
      #calculating amino acid length and location
      eligible_Peptides$AA_Length[i]<-nchar(eligible_Peptides$Amino_Acid_Sequence[i])
      eligible_Peptides$AA_Location[i]<-ceiling((as.numeric(input$SNP_Location)-eligible_Peptides$Peptide_Start[i]+1)/3)
      if(eligible_Peptides$Forward_Reverse[i] == "Forward")
      {
        eligible_Peptides$Old_AA[i]<-substr(eligible_Peptides$Amino_Acid_Sequence[i],eligible_Peptides$AA_Location[i],eligible_Peptides$AA_Location[i])
        #changing amino acid based on user
        substr(eligible_Peptides$Nucleotide_Sequence_New[i],as.numeric(input$SNP_Location)-eligible_Peptides$Peptide_Start[i]+1,as.numeric(input$SNP_Location)-eligible_Peptides$Peptide_Start[i]+1)<-as.character(input$Base_Choice)
      } else if(eligible_Peptides$Forward_Reverse[i] == "Reverse")
      {
        eligible_Peptides$Old_AA[i]<-substr(eligible_Peptides$Amino_Acid_Sequence[i],eligible_Peptides$AA_Location[i],eligible_Peptides$AA_Location[i])
        substr(eligible_Peptides$Nucleotide_Sequence_New[i],eligible_Peptides$Peptide_End[i]-as.numeric(input$SNP_Location)+1,eligible_Peptides$Peptide_End[i]-as.numeric(input$SNP_Location)+1)<-as.character(input$Base_Choice)
      }
      
      
      #translating AA sequences based on forward/reverse and translation types
      if(eligible_Peptides$Vert_Hg19[i] == "hg19" && eligible_Peptides$Forward_Reverse[i] == "Forward") 
      {
        eligible_Peptides$Amino_Acid_Sequence_New[i] = toString(translate(DNAString(eligible_Peptides$Nucleotide_Sequence_New[i])))
      } else if(eligible_Peptides$Vert_Hg19[i] == "vert" && eligible_Peptides$Forward_Reverse[i] == "Forward")
      {
        eligible_Peptides$Amino_Acid_Sequence_New[i] = toString(translate(DNAString(eligible_Peptides$Nucleotide_Sequence_New[i]), genetic.code=SGC1))
      } else if(eligible_Peptides$Vert_Hg19[i] == "hg19" && eligible_Peptides$Forward_Reverse[i] == "Reverse")
      {
        eligible_Peptides$Amino_Acid_Sequence_New[i] = toString(translate(reverseComplement(DNAString(eligible_Peptides$Nucleotide_Sequence_New[i]))))
      } else if(eligible_Peptides$Vert_Hg19[i] == "vert" && eligible_Peptides$Forward_Reverse[i] == "Reverse")
      {
        eligible_Peptides$Amino_Acid_Sequence_New[i] = toString(translate(reverseComplement(DNAString(eligible_Peptides$Nucleotide_Sequence_New[i])), genetic.code=SGC1))
      }
      
      #determining if AA changed
      if(eligible_Peptides$Amino_Acid_Sequence_New[i] == eligible_Peptides$Amino_Acid_Sequence[i]) 
      {
        eligible_Peptides$AA_Change[i] = "No"
      } else
      {
        eligible_Peptides$AA_Change[i] = "Yes"
      }
      
      eligible_Peptides$New_AA[i]<-substr(eligible_Peptides$Amino_Acid_Sequence_New[i],eligible_Peptides$AA_Location[i],eligible_Peptides$AA_Location[i])
      
      #Converting database naming back to more standard scientific jargon
      if(eligible_Peptides$Vert_Hg19[i] == "hg19")
      {
        eligible_Peptides$Vert_Hg19[i] <- "Cytoplasmic"
      } else
      {
        eligible_Peptides$Vert_Hg19[i] <- "Mitochondrial"
      }
      
      eligible_Peptides$RS_Number[i]<-RS_Output
    }
    
    colnames(eligible_Peptides)[13]<-"Strand"
    colnames(eligible_Peptides)[14]<-"Translation"
    
    eligible_Peptides_output<-eligible_Peptides %>% #trimming output
      select(Peptide_Name,Peptide_Start,Peptide_End,AA_Length,AA_Change,Old_AA,New_AA,AA_Location,Strand,Translation,Mass_Spec,RS_Number) %>%
      left_join(AA_PTM,by=c("Old_AA" = "Abbreviation")) %>%
      left_join(AA_PTM,by=c("New_AA" = "Abbreviation"))
    
    #renaming columns
    colnames(eligible_Peptides_output)[13]<-"Old_AA_PTM"
    colnames(eligible_Peptides_output)[14]<-"New_AA_PTM"
    colnames(eligible_Peptides_output)[8]<-"AA_SNP_Location"
    colnames(eligible_Peptides_output)[4]<-"Triplet_Length"
    
    
    eligible_Peptides_output
  })
  
  
  #logic below handles the password username entry, only rendering the protected data is proper credneitals are provided
  USER <- reactiveValues(Logged = Logged)
  
  observe({ 
    if (USER$Logged == FALSE) {
      if (!is.null(input$Login)) {
        if (input$Login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          Id.username <- which(my_username == Username)
          Id.password <- which(my_password == Password)
          if (length(Id.username) > 0 & length(Id.password) > 0) {
            if (Id.username == Id.password) {
              USER$Logged <- TRUE
            } 
          }
        } 
      }
    }    
  })
  observe({
    if (USER$Logged == FALSE) {
      
      output$page <- renderUI({
        div(class="outer",do.call(bootstrapPage,c("",ui1())))
      })
    }
    if (USER$Logged == TRUE) 
    {
      output$page <- renderUI({
        ui2()
        
      })
      print(ui)
    }
  })
})



#runApp(list(ui = ui, server = server))
