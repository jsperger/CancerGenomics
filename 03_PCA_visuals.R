# This script is a work in progress
# Its goal is to do add some PCA visualizations to help us understand the data
# This follows after using 00_Libs_Data_DL.R in memory (e.g. requires ddr.mat)
# You also need the barcode in data (run with your own respective barcode loaded in memory)

# these are two visualization packages I've been using
# autoplot is one handy plot, so is pairs()
install.packages('ggfortify')
install.packages('car')

# this is how we can abstract from our data in case we need to redo an analysis
# x <- ddr.mat
x <- subtype.mat

# we can see how many samples & genes
# I saw 379 samples, 808 genes
# For subtype, this is 1357 genes
dim(x)

# load the ggfortify plotting tool
# load the car plotting tool
library(ggfortify)
library(car)

# we're only going to look at the first four components, the data is already
# centered and pretty scaled
dnarep.pca <- prcomp(t(x), rank. = 4, center=TRUE) 
#dnarep.pca <- prcomp(t(x), rank. = 4)

# this shows you the proportion of variance explained by each component
screeplot(dnarep.pca)

# here are the principal components if you want to know the values
# dnarep.pca$rotation

# this is a standard plot of the PCA'd data - diagonals can be replaced
# scatterplotMatrix below just calls this with a better interface
# pairs(dnarep.pca$x)

# scatterplot matrix is a much more complex version of pairs
# complete with kernel density estimates - uncomment to see
# scatterplotMatrix(dnarep.pca$x, reg.line=FALSE, smoother=FALSE)

#### Groupings ####
g <- barcode$NMFC3

#printPretty - same as above but with groups
scatterplotMatrix(dnarep.pca$x, groups=g, reg.line=FALSE, smoother=FALSE)
