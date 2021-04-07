options(repos = BiocInstaller::biocinstallRepos())
library(shiny)
library(Biostrings)
library(tidyverse)
library(rsconnect)

#!diagnostics off

SGC1 <- getGeneticCode("SGC1")
MDP_database<-read_csv(REDACTEDDATABASENAME)
AA_PTM<-read_csv(REDACTEDDATABASENAME)
RS_Number<-read_csv(REDACTEDDATABASENAME)
colnames(AA_PTM)[2]<-"Abbreviation"
#replacing erroneous data in amino acid name
AA_PTM$Abbreviation[11]<-"L"
AA_PTM$`Amino Acid`<-NULL


MDP_database$Old_AA_PTM<-NA
MDP_database$New_AA_PTM<-NA
MDP_database$AA_Length <-NA
MDP_database$RS_Number<-NA

for(i in 1:nrow(MDP_database))
{
  substr(MDP_database$Amino_Acid_Sequence[i],1,1)<-as.character("M") #preventing weird start codon errors
  
  temp_list<-strsplit(MDP_database$ORF_name_short[i], split = "_")
  if(temp_list[[1]][6]=="dir") #extracting start and ends of amino acid sequences for forward
  {
    MDP_database$Peptide_Start[i] = temp_list[[1]][3]
    MDP_database$Peptide_End[i] = temp_list[[1]][4]
    MDP_database$Forward_Reverse[i] = "Forward"
    
  } else if (temp_list[[1]][6]=="rev") #extracting start and ends of amino acid sequences for reverse
  {
    MDP_database$Peptide_Start[i] = temp_list[[1]][4]
    MDP_database$Peptide_End[i] = temp_list[[1]][3]
    MDP_database$Forward_Reverse[i] = "Reverse"
  }
  if(substr(MDP_database$Codon[i],1,3) == "ATG") #extracting translation types
  {
    MDP_database$Vert_Hg19[i] = "hg19"
  } else
  {
    MDP_database$Vert_Hg19[i] = "vert"
  }
}

MDP_database$Peptide_End<-as.numeric(MDP_database$Peptide_End)
MDP_database$Peptide_Start<-as.numeric(MDP_database$Peptide_Start)

for(i in 1:nrow(MDP_database))
{
  if(MDP_database$Forward_Reverse[i] == "Reverse")
  {
    #taking reverseComplment in order to display It to the user
    MDP_database$Nucleotide_Sequence[i] = tolower(toString(reverseComplement(DNAString(MDP_database$Nucleotide_Sequence[i]))))
  }
}

MDP_database <- MDP_database %>%
  mutate(Nucleotide_Sequence_New = Nucleotide_Sequence) %>% #copying columns for translation on server side
  mutate(Amino_Acid_Sequence_New = Amino_Acid_Sequence)

Logged = FALSE;

my_username <- REDACTEDUSERNAME
my_password <- REDACTEDPASSWORD


#UI1 is the initial interface the user is presented with
ui1 <- function(){
  tagList(
    div(id = "login",
        wellPanel(textInput("userName", "Username"),
                  passwordInput("passwd", "Password"),
                  br(),actionButton("Login", "Log in"))),
    tags$style(type="text/css", "#login {font-size:10px;   text-align: left;position:absolute;top: 40%;left: 50%;margin-top: -100px;margin-left: -150px;}")
  )}

#UI2 is the main interface where the user can visualize SNPs
ui2 <- function(){tagList(fluidPage(
  titlePanel("Identifying effects of amino acids changes in smORFS"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("SNP_Location", value = "11", h3(
        "Enter SNP Location (11-16567):"
      )),
      
      selectInput("Base_Choice", h3(
        "Select nitrogenous base to change"),
        choices = list("Adenine" = "a",
                       "Cytosine" = "c",
                       "Guanine" = "g",
                       "Thymine" = "t"))
    ),
    #help text on the panel
    mainPanel(
      dataTableOutput("Peptides_Output"),
      helpText("If an error is produced when entering a SNP Location, the SNP does not exist. Please enter a different one."),
      helpText("If a cell is empty, the Mass Spec or the Modification does not exist."),
      helpText("Created by Henry Jiao in 2018 for the lab of Dr. Cohen in the USC School of Gerontology")
    )
  )
  
))}

ui = (htmlOutput("page"))