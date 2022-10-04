#Vibrio genome project 
#Heatmap analysis 

library(ggdendro)
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)

# Read in data
# read.delim uses tab delimiter, so be sure to set it to whichver
#    your data is stored as
# I set check.names to FALSE because R can manipulate them if they're numbers
# or anything that doesn't look right
data <- read.csv("~/OneDrive - Flinders/EDgenomes/phenotypic/carbon_phenotypic.csv", header = TRUE)

###############################################################
# Make distance matrix and perform clustering
###############################################################
# Remove first column of names and convert into matrix
mat <- as.matrix(data[, 2:ncol(data)])
rownames(mat) <- data[, 1]

# Remove any rows with incomplete data
mat <- mat[complete.cases(mat),]

# Cluster rows using default parameters
hcr <- hclust(dist(mat))

# Cluster columns using default parameters
#hcc <- hclust(dist(t(mat)))
hcc= col(mat)

# Create lists of the orders of the rows and columns
row.order <- hcr$order
#col.order <- hcc$order
col.order <- hcc
row.order.names <- rownames(mat[row.order,])

###############################################################
# Generate heatmap
###############################################################
# Create the ggplot data frame with the row and column orders
#df <- data.frame(mat[row.order, col.order], row.names=NULL, check.names=FALSE)
df <- data.frame(mat[row.order, col.order], row.names=NULL, check.names=FALSE)

# Create a column of the row names
df$rowname <- row.order.names

# Make the data in long format for ggplot
df.melt <- melt(df, id.vars=c("rowname"), variable.name="colname", value.name="datavalue")

# Make sure order of rows are still in the same order as before the melt command
df.melt$rowname <- factor(df.melt$rowname, levels=row.order.names)

# Create heatmap object
hm <- ggplot(df.melt, aes(x=colname, y=rowname, fill=datavalue)) +  # x = rows, y = columns, fill = value name
  geom_tile(colour="white") +  # This is color of lines inbetween the squares in the heatmap
  scale_fill_gradientn(colours=c("white",
                                 "blue",
                                 "black"),
                       name="Data\nValue") +  # name is the legend title
  theme(axis.ticks=element_blank(),  # Remove tick marks on the y and x axis
        axis.text=element_text(colour="black"),
        axis.text.x=element_text(hjust=1, vjust=0.5, angle=90, size=8),  # hjust, vjust, and angle together for a 90 degree
        # turn on the x-axis text
        axis.text.y=element_text(hjust=0),  # This left-aligns the y axis text
        panel.background=element_rect(fill="white")) +  # This sets the weird grey background ggplot makes into a white one
  xlab("") + ylab("")  # Remove heatmap labels on the y and x axis

# Create left dendrogram
dendro.l <- ggdendrogram(hcr, colour="black", rotate=T, size=3) +
  theme(axis.text.y=element_blank(),  # Remove row names since they're already there in the heatmap object
        axis.text.x=element_text(colour="black")) +
  scale_y_reverse()

# Create top dendrogram
dendro.t <- ggdendrogram(hcc, colour="black", rotate=F, size=3) +
  theme(axis.text.x=element_blank(),  # Remove column names since they're already there in the heatmap object
        axis.text.y=element_text(colour="black"))

blank <- grid.rect(gp=gpar(col="white"))  # Create blank plot to place in the top left

# Plot together
pl <- grid.arrange(dendro.l, hm, blank, ncol=2, padding=0, widths=c(1,3))

# Save heatmap: uncomment these lines to save the heatmap using the specified parameters
# ggsave("heatmap.svg",
#        plot=pl,
#        width=20,  # The width and height need to be changed depending on the number of rows and columns you have
#        height=20,
#        units="cm",
#        dpi=250)
