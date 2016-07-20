# creating heatmaps in R using gene expression data


# prerequisites

# install the relevant packages if they are not installed
install.packages("RColorBrewer")
install.packages("pheatmap")

# download the data files (it will download to the current working directory)
download.file("https://raw.githubusercontent.com/igordot/tutorials/master/ca-genes-stats-sig.csv", destfile="ca-genes-stats-sig.csv")
download.file("https://raw.githubusercontent.com/igordot/tutorials/master/ca-genes-fpkm.csv", destfile="ca-genes-fpkm.csv")

# load the relevant packages
library(RColorBrewer)
library(pheatmap)
