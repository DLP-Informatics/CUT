#VisScript.r

# This is the Decoupled visualizer script necessary for visualizing each organism on a per-organism and inter-organism basis



for(i in 1:length(orglist))
{
  if(orglist[i] == "SingleOrganismVisualization")
  {
    next
  }
  if(orglist[i] == "MultiOrganismVisualization")
  {
    next
  }
  
  name <- unlist(strsplit(orglist[i], split='.', fixed=TRUE))[1] #Grabs the name of the Output organism file
  names_array <- c(names_array, as.character(name)) #Adds the name to the name indexing array
  
  file_path <- as.character(paste(sourceDir, "/", orglist[i], sep = "")) #grabs file path to read in matrix
  mat <- as.matrix(read.csv(file=file_path, header=TRUE, row.names = 1, sep=",")) #assigns matrix to temp variable
  
  assign(paste("matrix", i, sep = ""), mat) # assigns matrix to the ith temp variable, corresponding to the ith entry in the
  # name indexing array to keep track of each matrix
  
  matrix_index_array <- c(matrix_index_array, paste("matrix", i, sep = "")) #stores the the matrix into an array where it's 
  #indice corresponds to the name of the temp matrix
  
}

for(i in 1:length(names_array))
{
  if(nameflag == "Forward")
  {
  names_array[i] <- str_replace(names_array[i], "Forward_CU_Analysis_", "")
  }
  if(nameflag == "Reverse")
  {
  names_array[i] <- str_replace(names_array[i], "Reverse_CU_Analysis_", "")
  }
}

#generate heatmaps of all organisms analyzed by CUT for individual analysis
outputflag <- "Single"
for(i in 1:length(names_array))
{
  setwd(path_HelperScripts)
  source("lattice_visualizer.r")
}

#Initialize empty objects for storing and data manipulation#
sumCodonTable <- GC3C
codonArr <- c()
manhattanMatrix <- matrix(0L, nrow = 64, ncol = 0)
rownames(manhattanMatrix) <- codonArrName
manhatttanCodon <- c()
##########################################################


#create a per-organism, per-codon data frame of all 64 codons ######################

i <- 1
j <- 1
for(i in 1:length(names_array))
{
  nameobj <- names_array[i]
  matobj <- get(matrix_index_array[i])
  manhatttanCodon <- c()
  
  for(j in j:length(matobj[,1]))
  {
    manhatttanCodon<- c(manhatttanCodon, matobj[j,])
  } 
  
  manhattanMatrix <- cbind(manhattanMatrix, manhatttanCodon)
  colnames(manhattanMatrix)[i] <- nameobj
  
}


manhattanDF <- data.frame(manhattanMatrix)
meltDF <- melt(t(manhattanDF))
colnames(meltDF)[1] <- "organism"
colnames(meltDF)[2] <- "Codon"
colnames(meltDF)[3] <- "Frequency"

# Generate the Organism specific manhattan plot of codon usage#
require(ggplot2)
require(reshape2)
ggobj1 <- ggplot(data = meltDF, aes(x=Codon, y = Frequency) ) + 
  theme(aspect.ratio=5/10) + theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 20)) +
  theme(axis.text.x=element_text(angle=90, size = 15, hjust=1), axis.text.y = element_text(size = 15)) + 
  labs(title = manTitle)+ theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
  theme(legend.text=element_text(size=13)) + geom_point(aes(colour=organism))
if(nameflag == "Forward")
{
  ggsave(file = "Org_ForwardCU_ManhattanPlot.png", plot = ggobj1, path = path_MultiOrgOutput, height = 10, width = 20)
}
if(nameflag == "Reverse")
{
  ggsave(file = "Org_ReverseCU_ManhattanPlot.png", plot = ggobj1, path = path_MultiOrgOutput, height = 10, width = 20)
}
##############################################################

# Generate the Codon specific manhattan plot of codon usage#
ggobj2 <- ggplot(data = meltDF, aes(x=Codon, y = Frequency) ) + 
  theme(aspect.ratio=5/10) + theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 20)) +
  theme(axis.text.x=element_text(angle=90, size = 15, hjust=1), axis.text.y = element_text(size = 15)) + 
  labs(title = manTitle)+ theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
  theme(legend.position = 'none') + geom_point(aes(colour=Codon))
if(nameflag == "Forward")
{
  ggsave(file = "Codon_ForwardCU_ManhattanPlot.png", plot = ggobj2, path = path_MultiOrgOutput, height = 10, width = 20)
}
if(nameflag == "Reverse")
{
  ggsave(file = "Codon_ReverseCU_ManhattanPlot.png", plot = ggobj2, path = path_MultiOrgOutput, height = 10, width = 20)
}

##########################################################


# Generic Summed barplot of all 64 codons ###############
i <- 1
for(i in 1:length(manhattanMatrix[,1]))
{
 codonSum <- 0
 codonSum <- sum(manhattanMatrix[i,])
 codonArr <- c(codonArr, codonSum)
}

setwd(path_MultiOrgOutput)
write.csv(manhattanMatrix, "Codon_Data_Frame.csv")
png(file = barplotname, width = 1280, height = 720, units = "px")
barObj <- barplot(codonArr, main = Titlelabel, xlab = "Codons", names.arg = codonArrName, las=2)
dev.off()
##########################################################







