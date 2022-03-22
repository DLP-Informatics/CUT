#Output Conversion script
#David L Patton / 800728881
#Sung Lab
library("lattice")

if(outputflag == "Single")
{
  setwd(path_SingleOrgOutput)
}

if(outputflag == "Multi")
{
  setwd(path_MultiOrgOutput)
}

matrix_name <- names_array[i]
output_data_matrix <- get(matrix_index_array[i])
HeatmapTitle <- paste("Codon Usage Heatmap: ", matrix_name, sep = "")

output_pic_name <- paste( matrix_name, ".png", sep = "")
#output_pic_name <- "lattice_test_image.png"

d <- t(output_data_matrix)
dlin <- linspace(min(d), max(d), n = 50) #this serves as the key scale for the heatmap label
dlin <- round(dlin, digits = 3)
## FIGURE OUT HOW TO MANUALLY CREATE THE NUMERIC SCALE IN LATTICE##
png(output_pic_name, width = 600, height = 800, units = "px")
rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
heatplot <- levelplot(d, col.regions=rgb.palette, scales=list(x=list(rot=90), cex = 1.5),
                      xlab = "", ylab = "", main = list(label = HeatmapTitle, cex = 1.5),
                      aspect="fill",
                      colorkey=list(labels= list(cex = 1.5, as.character(dlin)),
                                    col=(rgb.palette))
                      )

print(heatplot)
dev.off()



