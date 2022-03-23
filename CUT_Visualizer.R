#CUT_visualization.r
#This is the Codon Usage Visualizer Script. The primary function of this to take Codon Usage output
#from the CUT pipeline and easily visualize and disseminate relevant information.

#This if statement exists if the user runs CUT_visualizer independent of CUT.r and initializes the requisite directory paths needed
# to access input files
library("stringr")
library("pracma")
library("ggplot2")
library("reshape2")

if(!exists("updateflag"))
{
  username <- Sys.info()[7]
  User <- paste("/Users/", username, sep ="")
  HomeDir <- paste("/Users/", username, "/Desktop/CUT", sep = "")
  
  path_FASTA <- paste(HomeDir, "/Genbank_Organisms/FASTA_src", sep = "")
  path_GBK <- paste(HomeDir, "/Genbank_Orgnaisms/GBK_src", sep = "")
  path_Output <- paste(HomeDir, "/Genbank_Organisms/CU_Output", sep = "")
  path_Output_Forward <- paste(path_Output, "/Forward", sep = "")
  path_Output_Reverse <- paste(path_Output, "/Reverse", sep = "")
  path_HelperScripts <- paste(HomeDir, "/HelperScripts", sep = "")
  
  cols <- c("T", "G", "C",	"A")
  rows <- c("T[X]T", "T[X]G","T[X]C","T[X]A","G[X]T","G[X]G","G[X]C","G[X]A",
            "C[X]T","C[X]G","C[X]C","C[X]A", "A[X]T","A[X]G", "A[X]C","A[X]A")
  GC3C <- matrix(0L, nrow =16, ncol = 4) #Matrix to store Genome Coding Triplet count
  rownames(GC3C) <- rows
  colnames(GC3C) <- cols
}

codonArrName <- c("TTT", "TGT", "TCT", "TAT",
              "TTG", "TGG", "TCG", "TAG",
              "TTC", "TGC", "TCC", "TAC",
              "TTA", "TGA", "TCA", "TAA",
              "GTT", "GGT", "GCT", "GAT",
              "GTG", "GGG", "GCG", "GAG",
              "GTC", "GGC", "GCC", "GAC",
              "GTA", "GGA", "GCA", "GAA",
              "CTT", "CGT", "CCT", "CAT",
              "CTG", "CGG", "CCG", "CAG",
              "CTC", "CGC", "CCC", "CAC",
              "CTA", "CGA", "CCA", "CAA",
              "ATT", "AGT", "ACT", "AAT",
              "ATG", "AGG", "ACG", "AAG",
              "ATC", "AGC", "ACC", "AAC",
              "ATA", "AGA", "ACA", "AAA")


path_SingleOrgOutput <- paste(path_Output_Forward, "/SingleOrganismVisualization", sep = "")

if(!dir.exists(path_SingleOrgOutput))
{
  dir.create(path_SingleOrgOutput)
}

path_MultiOrgOutput <- paste(path_Output_Forward, "/MultiOrganismVisualization", sep = "")

if(!dir.exists(path_MultiOrgOutput))
{
  dir.create(path_MultiOrgOutput)
}

print("Performing Forward Strand Analysis.")

sourceDir <- path_Output_Forward
names_array <- c()
matrix_index_array <- c()
orglist <- list.files(sourceDir)
barplotname <- "FowardCU_BarPlot.png"
Titlelabel <- "Total Codon Usage Frequency - Forward Strand"
manTitle <- "Codon Usage Manhattan Plot - Forward Strand"
nameflag <- "Forward"

setwd(path_HelperScripts)
source("VisScript.r")

path_SingleOrgOutput <- paste(path_Output_Reverse, "/SingleOrganismVisualization", sep = "")

if(!dir.exists(path_SingleOrgOutput))
{
  dir.create(path_SingleOrgOutput)
}

path_MultiOrgOutput <- paste(path_Output_Reverse, "/MultiOrganismVisualization", sep = "")

if(!dir.exists(path_MultiOrgOutput))
{
  dir.create(path_MultiOrgOutput)
}


print("Performing Reverse Strand Analysis.")
sourceDir <- path_Output_Reverse
names_array <- c()
matrix_index_array <- c()
orglist <- list.files(sourceDir)
nameflag <- "Reverse"
barplotname <- "ReverseCU_BarPlot.png"
Titlelabel <- "Total Codon Usage Frequency - Reverse Strand"
manTitle <- "Codon Usage Manhattan Plot - Reverse Strand"

setwd(path_HelperScripts)
source("VisScript.r")

