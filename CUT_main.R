#### Job Talk Coding Project ######

#We would like to investigate codon usage in microbes. Your goal is to pull in relevant data
# and organize it in a way that can be shared with internal scientists. 

#1. Obtain >1000 microbial genomes through sources like RefSeq
#2. Tabulate protein-coding codons
#3. Provide the results in a tabular format and preferably, with visualization.
#4. (Bonus)  Are there any interesting trends in codon counts#
# or percentages based on phylogeny, niche, etc.?


########Step 0: Install requisite packages#################

## For Loop to install the requisite packages (also used in CDMAP) ######
packages <- c("seqinr", "BiocManager", "lattice", "tidyverse", "vcfR", "stringr", "here", "beepr", "pracma")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
BiocManager::install("genbankr")
library("genbankr") # BiocManager package for genbank file parsing and manipulation



####Step 1: Importing requiste input Bacterial genobank files######



###########################################################

#For a Single Organism
##################################################
#cat("What is the name of the organism?")
#organism <- readLines(con = "stdin", n = 1)

#cat("What is your reference sequence? (please provide the full Path)")
#Path_RefFile <- readLines(con = "stdin", n = 1)

#cat("What is your Genbank file? (please provide the full Path)")
#Path_GBFile <- readLines(con = "stdin", n = 1)
##################################################



#for iterating through an entire directory
##################################################

#simple setting of the working directory in a semi-dynamic fashion, would expand this to be more generic in future
username <- Sys.info()[7]
User <- paste("/Users/", username, sep ="")
HomeDir <- paste("/Users/", username, "/Desktop/CUT", sep = "")

path_FASTA <- paste(HomeDir, "/Genbank_Organisms/FASTA_src", sep = "")
path_GBK <- paste(HomeDir, "/Genbank_Organisms/GBK_src", sep = "")
path_Output <- paste(HomeDir, "/Genbank_Organisms/CU_Output", sep = "")
path_Output_Forward <- paste(path_Output, "/Forward", sep = "")
path_Output_Reverse <- paste(path_Output, "/Reverse", sep = "")
path_HelperScripts <- paste(HomeDir, "/HelperScripts", sep = "")

if(!dir.exists(path_Output))
{
  dir.create(path_Output)
}
if(!dir.exists(path_Output_Forward))
{
  dir.create(path_Output_Forward)
}
if(!dir.exists(path_Output_Reverse))
{
  dir.create(path_Output_Reverse)
}

directory_FASTA <- list.files(path_FASTA)
directory_GBK <- list.files(path_GBK)


##################################################

#Setup of Requisite Matrices
##################################################

cols <- c("T", "G", "C",	"A")
rows <- c("T[X]T", "T[X]G","T[X]C","T[X]A","G[X]T","G[X]G","G[X]C","G[X]A",
          "C[X]T","C[X]G","C[X]C","C[X]A", "A[X]T","A[X]G", "A[X]C","A[X]A")

GC3C <- matrix(0L, nrow =16, ncol = 4) #Matrix to store Genome Coding Triplet count
rownames(GC3C) <- rows
colnames(GC3C) <- cols

CU_Forward <- GC3C
CU_Reverse <- GC3C


##################################################

##########Question 2: Tabulate protein-coding codons ##########################################


#Existing CDMAP Code#######

k <- 1
for(k in 1:length(directory_FASTA))
  {
  
  FASTA_path <- paste(path_FASTA, "/", directory_FASTA[k], sep = "")
  GBK_path <- paste(path_GBK, "/", directory_GBK[k], sep = "")

  organism <- str_remove(directory_GBK[k], ".gbk")

  FASTAobj <- getSequence(read.fasta(FASTA_path)) #Path to fasta file
  GBKobj <- parseGenBank(GBK_path) #parsing the genbank file and converting to a dataframe
    
  RefSeq_arr <- unlist(FASTAobj)

  #Developer note:
  #When instantiating GBKOBJ the following error occured with Clostridioides_difficile_ATCC9689:
  #Error in .Call2("new_XString_from_CHARACTER", class(x0), string, start,  : 
  #                  key 79 (char 'O') not in lookup table
  # Will need to debug this, it's very possible an issue with the formatting of the individual gbk


  featlength <- length(GBKobj$FEATURES) #number of gene features contained in the genbank file

  startindex <- c() #array of all all gene starting positions 
  endindex <- c() #array of all gene ending positions
  featlengthindex <- c() #array of the length of all gene features



  i <- 1

  for(i in  1:featlength)
  {
 
    if(i == 1) #This is to skip the first entry, it's often a copy of the reference nucleotide sequence
    {
      #print("skipping first entry")
      next()
    }
  
  
    featuretype <- GBKobj$FEATURES[[i]]$type #extract the feature field (gene, CDS, etc)
    print(paste("feature #: ", i, " feature type: ", featuretype, sep = ""))
  

    if(featuretype != "gene") #we used the 'gene' identifier to analyze protein coding codon
    {
      print(paste("Feature", i, "skipped. Is not a protein coding sequence", sep = " "))
      next()
    }
    dataobj <- GBKobj[["FEATURES"]][[i]]
    matobj <- data.matrix(GBKobj[["FEATURES"]][[i]]) #extracting an individual gene feature from the data frame
    #matobj <- testlist[[1]]$type
  
    #Grab the start and end position of the coding region
    featstart <- as.numeric(matobj[1,2]) #gene nucleotide start position
    featend <- as.numeric(matobj[1,3]) #gene nucelotide end position
  
  
    startindex <- c(startindex, featstart) 
    endindex <- c(endindex, featend) 
  
    if(featstart == 1)
    {
      upNuc <- RefSeq_arr[length(RefSeq_arr)] 
      downNuc <- RefSeq_arr[featend+1]
    } else if(featend == length(RefSeq_arr))
    {
      upNuc <- RefSeq_arr[featstart-1] 
      downNuc <- RefSeq_arr[1]
    }else
    {
      upNuc <- RefSeq_arr[featstart-1] #leftmost upstream nucleotide 4mer position in gene coding region
      downNuc <- RefSeq_arr[featend+1] # rightmost downstream nucleotide 4mer position in gene coding region
    }
  
  
    #Checks for null reference, skips if a null reference is induced
    if(identical(upNuc, character(0)) | identical(downNuc, character(0)) ) #is.element(featend, endindex)
    {
      next()
    }
  
    gene_arr <- RefSeq_arr[featstart:featend] #isolation of a specific gene feature
    gene_end <- length(gene_arr)-1 
    featlengthindex <- c(featlengthindex, length(gene_arr)) #length of the gene feature

    gene_arr <- RefSeq_arr[featstart:featend] #isolation of a specific gene feature
    gene_end <- length(gene_arr)-1 
    featlengthindex <- c(featlengthindex, length(gene_arr)) #length of the gene feature
  

  
  
    #checks for null reference, skips if null detected
    if(is.na(upNuc))
    {
      i <- i+1
      next
    }
    if(is.na(downNuc))
    {
      i <- i+1
      next
    }  
  
    #checks if an ambiguous nucleotide 'N' assigned, if true, it skips it
    upNuc <- toupper(upNuc)
    downNuc <- toupper(downNuc)
   
    if(upNuc != 'A' & upNuc != 'T' & upNuc != 'C' & upNuc != 'G')
    {
      i <- i+1
      next
    }
    if(downNuc != 'A' & downNuc != 'T' & downNuc != 'C' & downNuc != 'G')
    {
      i <- i+1
      next
    }
  
    orientation <- GBKobj$FEATURES[[i]]$strand
  
    if(orientation == '+')
    {
      #setwd(Path_to_scripts)
      print("Strand in the forward '+' Orientation")
    }
    if(orientation == '-')
    {
      #setwd(Path_to_scripts)
      print("Strand in the reverse '-' Orientation")
    }

  
    j <- 2 #iterator instantiated with respect to center nucleotide, grabs j-1 and j+1 for left and right nucleotide.
    while(j <= gene_end)
    {
      iter <- j
      setwd(path_HelperScripts)
      source('CodonUsage_Tabulator.r')
    
    
      j <- j+3 #iterates by 3, so it checks with respect to each nucleotide triplet.
    }
    #DiagArr <- sort(DiagArr)
  
  }

  name_forward <- paste("Forward_CU_Analysis_", organism, ".csv", sep = "")
  setwd(path_Output_Forward)
  write.csv(CU_Forward, name_forward)

  name_Reverse <- paste("Reverse_CU_Analysis_", organism, ".csv", sep = "")
  setwd(path_Output_Reverse)
  write.csv(CU_Reverse, name_Reverse)
}

beep(3)
print("Codon Directory Analysis Complete!")

cat("Would you like to update your multi-organism analysis folder?")
updateflag <- readLines(con = "stdin", n = 1)
# end loop here