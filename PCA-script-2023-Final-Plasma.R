

###########################
# Packages and workspace  #
###########################

### Packages are collections of functions, data and codes that you can import into your R workspace

## To see packages imported into your session
sessionInfo() #search() might also be used

## To see all packages installed on your computer
library()
# NB : to close this list of librairies, and go back to the command prompt: q

## To Import a package into your workspace

library(FactoMineR)	#(example given with FactoMineR that we use for the PCA)

## You may check that FactoMineR has been imported into your R workspace
search()

#################################
# Choose your working directory #
#################################

### R can import or export data to files
### By default, R looks for files inside your working directory
### -> il you want to read files, or write into files, located in a specific directory, you will need to set the working directory


## To know what your actual working directory is:
getwd() 


## To choose your working directory with an absolute path
#setwd("/my/absolute/path/data/") # on adenine, the absolute path should be like this "/srv/home/login/meg_m2_gac" with your won login.
# alternatively, you can use the symbol "~" which stands for your home:
#setwd("~/meg_m2_gac")


########################################################################
# Create a dataframe in your work space that contains data from a file #
########################################################################

### To import data: 
##		- from a file in a table format
##		- to a data frame


## read.table function

accdata <- read.table("6-GSE10927-norm.txt", header=TRUE, sep="\t", dec=".", row.names=1)

# read.table is a function that create a dataframe from a file in a table format; between brackets are the following arguments:
# "6-GSE10927-norm-alldata.txt" is the file name
# header=TRUE indicates the presence of a headline in the file
# sep="\t" indicates that, in the file, columns are separated by tabulations
# dec="." indicates that decimal symbols are dots
# row.names = 1 indicates that names of rows in the dataframe are written in the first column of the file


#########################################################################
# Perform the Principal Component Analysis using the FactoMineR package #
#########################################################################

## PCA function 
## has been imported into your workspace with the FactoMineR librairy

res.pca <- PCA(accdata[,1:1905], quali.sup = 1901:1904, quanti.sup = 1905, graph = FALSE)

# accdata[,1:1905]: the dataframe "accdata" has 1907 columns. Col 1906 and 1907 are excluded since they contain viability data that would disrupt the PCA analysis
# results are sent in the "res" object
# accdata is the dataframe you have just created
# quali.sup = 1901:1904 indicates that the columns 1901 to 1904 of the dataframe contain supplementary qualitative variables not used for PCA
# quanti.sup = 1905 indicates that column 1905 of the dataframe contains a supplementary quantitative variable not used for PCA
# graph = FALSE: by default, graph = TRUE and generates graphs. these "default" graphs would be unlisible because too much variables are used.

### By default, 2 plots are created by the PCA function: one with individuals and one with variables. 
### Too many individuals and variables are used by PCA function on ACC data -> simplification of plots is required!


#############################################################
################## Create plots of the PCA ##################
#############################################################


# Plot of individuals
######################


## plot function (and arguments between brackets)
# Let's plot the PCA result for individuals !

plot(res.pca, choix = "ind", cex=0.6, invisible = "quali")

# res.pca is the name of the object used by the plot function
# invisible = "quali" indicates that qualitative variables centroid shall not be drawn on the plot

## To save your plot:
jpeg("PCA_ind.jpg", width = 700, height = 700)
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali")
dev.off()


# Plot of individuals with annotations (supplementary variables) and try to find a biological meaning to this PCA
#################################################################################################################

# Add the information of supplementary qualitative variables
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", autoLab = "y")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "clusterforACCs", autoLab = "y")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Side", autoLab = "y")

# res.pca is the name of the object used by the plot function
# invisible = "quali" indicates that qualitative variables centroid shall not be drawn on the plot
# habillage colors individuals among a categorical variable (here "Diagnosis", which is the name of the column)

# Add the information of supplementary qualitative variables, for axes 2 and 3
plot(res.pca, axes=c(2,3), choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", autoLab = "y")
plot(res.pca, axes=c(2,3), choix = "ind", cex=0.6, invisible = "quali", habillage = "clusterforACCs", autoLab = "y")
plot(res.pca, axes=c(2,3), choix = "ind", cex=0.6, invisible = "quali", habillage = "Side", autoLab = "y")

# res.pca is the name of the object used by the plot function
# invisible = "quali" indicates that qualitative variables centroid shall not be drawn on the plot
# habillage colors individuals among a categorical variable (here "Diagnosis", which is the name of the column)
# axes=c(2,3) plot the 2nd and 3rd maximum inertia axes



# Plot of variables, try to find which mRNA levels are correlated to biological information
###########################################################################################

plot(res.pca,choix="var", cex=0.6, select="contrib 15", unselect = 1) 

# select="contrib 15" select the first 15 variables that contribute to PC1 and PC2 (axes illustrated by default)

# Have a look on mRNA contributing to axes 2 and 3
plot(res.pca, axes = c(2,3), choix="var", cex=0.6, select="contrib 15", unselect = 0.96)



# Export to files data important for interpretation

## Export the correlation coefficient of mRNA levels with Principal Components

# write.infile function
write.infile(res.pca$var$cor, file="varcor.txt", sep="\t")



#################################################################################
########### Export to files informations important for interpretation ###########
#################################################################################

# To know the informations contained in the res.pca object
res.pca

# Principle Components (coordinates of individual on axes)
write.infile(res.pca$ind$coord, file="PrincipalComponentsInd.txt", sep="\t")

# Correlation of mRNA levels with axes
write.infile(res.pca$var$cor, file="CorrelationVar.txt", sep="\t")

# Coordinates of mRNA levels with
write.infile(res.pca$var$coord, file="CoordinVar.txt", sep="\t")



##############################################################################################
########### Quality of the interpertation depends on your knowledge in the field ! ###########
##############################################################################################



#######################################################################
########### Contribution of Principal Components to Inertia ###########
#######################################################################

## Eigenvalues are the variances of the principal components.
## The following eigen values matrix contains the inertia value explained by each of the 64 PCs 
eig.val <- res.pca$eig

## Let us draw a plot
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        ylim = c(0,25),
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")

## Zoom on the first 25 components

jpeg("QualityProj_25C.jpg", width = 700, height = 700)
eig.val <- res.pca$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        ylim = c(0,25),
        xlim = c(0,25),
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")
dev.off()


## Export to a file the eigenvalues
write.infile(res.pca$eig, file="QualityProj_PC.txt", sep="\t")



summary(res.pca, ncp=3, nbelements=Inf, file="pca-essai.txt" )

# res is the name of the object used by the summary function
# ncp argument gives the number of dimensions for which the summary is given
# nbelements is the number of elements

plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "contrib 5")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "cos2 0.5")



#############################
# PCA considering only ACCs #
#############################

ACConly <- subset(accdata[,1:1905], accdata$Diagnosis == 'ACC Low Grade' | accdata$Diagnosis == 'ACC High Grade')
res2 <- PCA(ACConly, quali.sup = 1901:1904, quanti.sup = 1905, graph = FALSE)
plot(res2, cex=0.6, invisible = "quali", habillage = "Diagnosis", autoLab = "y")
plot(res2, choix = "var", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "contrib 10", unselect = 1)








