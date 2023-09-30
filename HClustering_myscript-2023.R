
########################################################################
# Create a dataframe in your work space that contains data from a file #
########################################################################

## read.table function

accdata <- read.table("6-GSE10927-norm.txt", header=TRUE, sep="\t", dec=".", row.names=1)


########################################################################
# Calcultation of Euclidian distances between tumor samples            #
########################################################################

# Calcultation of Euclidian distances between tumor sample

distanceCalculee <- dist(accdata[,1:1905])

# Create the dendrogram parameters

hclust(distanceCalculee, method = "complete", members = NULL)

# Draw the dendogram

jpeg("cluster.jpg", width = 700, height = 700)
plot(clusters)
dev.off()


