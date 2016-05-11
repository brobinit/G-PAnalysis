# barGraph.R
# written by Blaire Robinson
#
# Description: reads in a .txt file summarizing data from the PMAnalyzer pipeline (written
# by Daniel Cuevas) and outputs and clustered bar graph comparing growth level between 
# samples accross tested substrates.
#
# Imput file: a text file with Sample name, well ID, compound tested, and associated 
# growth level. 
#
# After bar graph is created, user can save in format desired.


library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggdendro)
library(gridExtra)
library(grid)
library(mvtnorm)

#file selection
data <- read.table("~/Documents/SDSU_Grad/2/spring2016/Thesis/Figures/Analysis_2016mar23/Vcyclitrophicus_params_average_seq.txt", header=T, sep="\t")

# Grab only the sample name, well ID, compound tested, and growth level and
# create a data frame:
# SampleName | A1 | A2 | ...
# Each cell holds the growth level
gls <- dcast(data[c(1, 3, 5, 10)], sample ~ well + compound, value.var="growthlevel") 

# Convert dataframe to matrix and ignore the first column of sample names
#gls.mat <- scale(as.matrix(gls[,2:ncol(gls)]))  # Use when data scaling is required
gls.mat <- as.matrix(gls[,2:ncol(gls)])
rownames(gls.mat) <- gls$sample
gls.mat <- gls.mat[complete.cases(gls.mat),]
hcr <- hclust(dist(gls.mat))
# hcr <- hclust(dist(abs(cor(t(gls.mat), method="pearson"))))  # Use for Pearson correlation
# sets the row order, comment out to specifically order
#row.order <- hcr$order
# row order starts at 0,0 and works up
row.order <-c("ED144","ED252")

hcc <- hclust(dist(t(gls.mat)))
#hcc <- hclust(dist(abs(cor(gls.mat), method="pearson")))  # Use for Pearson correlation

# sets the column order, comment out to specifically order
col.order <- hcc$order

# Use a specific order for the columns
# col.order <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
#                "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12",
#                "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
#                "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
#                "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
#                "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
#                "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
#                "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")

row.order.names <- rownames(gls.mat[row.order,])
col.order.names <- colnames(gls.mat[,col.order])

# Create dataframe to be used for plotting
df <- data.frame(gls.mat[row.order, col.order])
df$sample <- row.order.names

# Flatten dataframe into 3 columns
df.melt <- melt(df, id.vars=c("sample"), variable.name="well", value.name="growthlevel")
df.melt$sample <- factor(df.melt$sample, levels=row.order.names)

#create associated bar graph. 
bar <- ggplot(data=df.melt, aes(x=factor(well), y=growthlevel, fill=sample)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() + 
  scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values=c("green","blue")) + 
  ylim(-0.25,1.5) +
  theme(axis.ticks=element_blank(),
        axis.text.x=element_text(hjust=1, vjust=0.5, angle=90, size=10),
        axis.text.y=element_text(hjust=0, vjust=0.5, size=10),
        legend.title.align=0.5) + 
  xlab("") + ylab("Growth Level")
plot(bar)
