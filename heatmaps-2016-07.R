# creating heatmaps in R using gene expression data



# part 1 - set working directory (important for any project)

# select Desktop on a Mac (or any other directory)
setwd("~/Desktop")



# part 2 - get prerequisites

# install the relevant packages if they are not installed
install.packages("RColorBrewer")
install.packages("pheatmap")

# download genes with FPKM values (filtered for detectable genes)
fileName = "ca-genes-fpkm.csv"
download.file(paste0("https://raw.githubusercontent.com/igordot/tutorials/master/", fileName), destfile=fileName)

# download genes with stats (filtered for significant genes)
fileName = "ca-genes-stats-sig.csv"
download.file(paste0("https://raw.githubusercontent.com/igordot/tutorials/master/", fileName), destfile=fileName)



# part 3 - load the relevant packages and example data

# load packages
library(pheatmap)
library(RColorBrewer)

# import gene FPKM values
geneVals = read.csv(file="ca-genes-fpkm.csv", row.names=1, check.names=F, stringsAsFactors=F)

# check what you imported
head(geneVals)
dim(geneVals)

# import gene stats
geneStats = read.csv(file="ca-genes-stats-sig.csv", row.names=1, check.names=F, stringsAsFactors=F)

# check what you imported
head(geneStats)
dim(geneStats)



# part 4 - check and adjust expression values

# check the distribution of expression values ("las=2" means labels are perpendicular)
boxplot(geneVals, las=2)

# take a log of values to reduce variance
geneValsLog = log2(geneVals)
head(geneValsLog)

# take a log of values accounting for zeroes
geneValsLog = log2(geneVals + 0.01)
head(geneValsLog)

# check the distribution of logged values
boxplot(geneValsLog, las=2)



# part 5 - generate a basic heatmap

# learn more about pheatmap
?pheatmap

# default heatmap using original values (use 500 genes to save time)
heatmapVals = geneVals[1:500,]
dim(heatmapVals)
pheatmap(heatmapVals)

# default heatmap using logged values (use 500 genes to save time)
heatmapValsLog = geneValsLog[1:500,]
dim(heatmapValsLog)
pheatmap(heatmapValsLog)

# center and scale in the row direction
pheatmap(heatmapValsLog, scale="row")



# part 6 - color options

# manually create color range
myColors = c("green", "black", "red")
myColors

# expand the color range
myColors = colorRampPalette(myColors)(50)
myColors

# visualize results
pheatmap(heatmapValsLog, scale="row", color=myColors)

# learn more about display.brewer.all
?display.brewer.all

# see RColorBrewer color palettes
display.brewer.all()

# select a color palette
myColors = brewer.pal(n=11, name="RdBu")
myColors = colorRampPalette(myColors)(50)
myColors

# visualize results
pheatmap(heatmapValsLog, scale="row", color=myColors)



# part 7 - subset genes

# examine the values data frame
dim(geneValsLog)
head(geneValsLog)

# examine the stats data frame
dim(geneStats)
head(geneStats)

# subset for significant genes
heatmapAllSig = geneValsLog[rownames(geneStats),]
dim(heatmapAllSig)

# visualize results
pheatmap(heatmapAllSig, scale="row", color=myColors)

# filter significant genes based on q-value less than 0.01
geneStats[,"q_value"] < 0.01
geneStatsSubset = geneStats[geneStats[,"q_value"] < 0.01,]
dim(geneStatsSubset)
head(geneStatsSubset)

# sort the data frame of significant genes by fold change
geneStatsSorted = geneStatsSubset[order(geneStatsSubset[,"log2_fold_change"]),] 
head(geneStatsSorted)
geneStatsSorted = geneStatsSubset[order(abs(geneStatsSubset[,"log2_fold_change"]), decreasing=T),] 
head(geneStatsSorted)

# get 50 significant genes with the highest fold change
heatmapTop50 = geneValsLog[rownames(geneStatsSorted[1:50,]),]
dim(heatmapTop50)
head(heatmapTop50)

# visualize results
pheatmap(heatmapTop50, scale="row", color=myColors)

# remove the border around each cell
pheatmap(heatmapTop50, scale="row", color=myColors, border_color=NA)

# reduce row label font
pheatmap(heatmapTop50, scale="row", color=myColors, border_color=NA, fontsize_row=6)



# part 8 - determine gene order

# save heatmap as an object
p = pheatmap(heatmapTop50, scale="row", color=myColors, border_color=NA, fontsize_row=6)

# gene order
p$tree_row$order

# input data frame sorted by heatmap order
heatmapTop50[p$tree_row$order,]



# part 9 - label groups

# column (sample) labels
colLabels = data.frame(group=colnames(geneValsLog), stringsAsFactors=F)
colLabels
rownames(colLabels) = colLabels[,"group"]
colLabels
colLabels[,"group"] = gsub("AL_05.*", "young", x=colLabels[,"group"])
colLabels[,"group"] = gsub("AL_15.*", "old", x=colLabels[,"group"])
colLabels

# visualize results
pheatmap(heatmapTop50, scale="row", color=myColors, border_color=NA, fontsize_row=6,
annotation_col=colLabels)

# column (sample) label colors
colColors = list(group = c("green", "brown"))
colColors
names(colColors$group) = unique(colLabels[,"group"])
colColors

# visualize results
pheatmap(heatmapTop50, scale="row", color=myColors, border_color=NA, fontsize_row=6,
annotation_col=colLabels, annotation_colors=colColors)



# if all else fails, some graphical alternatives:
# https://github.com/igordot/genomics/blob/master/notes/heatmaps.md
