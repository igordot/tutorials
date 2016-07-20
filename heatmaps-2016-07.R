# creating heatmaps in R using gene expression data


# prerequisites

# install the relevant packages if they are not installed
install.packages("RColorBrewer")
install.packages("pheatmap")

# download the data files (it will download to the current working directory)
fileName="ca-genes-fpkm.csv"
download.file(paste0("https://raw.githubusercontent.com/igordot/tutorials/master/", fileName), destfile=fileName)
fileName="ca-genes-stats-sig.csv"
download.file(paste0("https://raw.githubusercontent.com/igordot/tutorials/master/", fileName), destfile=fileName)

# load the relevant packages
library(RColorBrewer)
library(pheatmap)
