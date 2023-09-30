#!/usr/bin/env Rscript

################################################################################
## M1 MEG option bioiformatique-genomique
## R Session 1 - commandes du diaporama + extra examples
## Claire Vandiedonck
################################################################################

#--------------------
# tips before starting: 
# -------------------
# Open this .Rscript in Rstudio for example
# To have a nice display of the comments wrapped, in R studio, select "Tools" in the menue, then "Options", "Code" and check "Soft-wrap R source files"
# You may comment some lines that do not necesserely need to be done (eg: lines 129, 134, 135 and 147) if sourcing the whole script.


#---------------------------------------------------
# the different kinds of variables and objects in R: slides 12-15
# --------------------------------------------------

## Classes of R objects

# slide 17
x <- c(3,7,1,2)
x<2
x==2
s <- "this is a string of characters"
mode(x)
mode(s)
class(x)
class(s)
as.numeric(x<2)

# extra examples:
is.numeric(as.numeric(x))
is.numeric(FALSE)
is.integer(FALSE)
as.numeric(FALSE)
as.numeric(TRUE)
is.logical(NA)
is.logical(x)
is.logical(x<2)
as.integer(x<2)
x <- c(1, "1")
x
is.numeric(x)
c(FALSE, 1, "3")
is.logical(c(FALSE, 1, "3"))
is.numeric(c(FALSE, 1, "3"))
is.character(c(FALSE, 1, "3"))
mode(c(FALSE, 1, "3"))
typeof(c(FALSE, 1, "3"))

#---------------------------------------------------
# vectors: slides 21-29
# --------------------------------------------------

## Vectors

# slide 21:
a <- c()
a
weight <- c(60, 72, 57, 90, 95, 72)
weight
4:10
seq(4,10)
seq(2,10,2)
rep(4,2)

rep(seq(4,10,2))
c(rep(1,4),rep(2,4))

# slide 22:
c(5,s)
mode(c(5,s))
class(c(5,s))

length(1:10)
length(weight)
str(weight)

# slide 24
size <- c(1.75, 1.8, 1.65, 1.9, 1.74, 1.91)
size^2

weight <- c(60, 72, 57, 90, 95, 72)
bmi <- weight/size^2 

# slide 25
size
sort(size)
mean(size)
sd(size)
median(size)
min(size)
max(size)
range(size)
summary(size)

# slide 28
#1:20 # example depends on the width of your screen -> see index at the beginning of each line
#5:10
1:100
50:100

size
size[1]
size[2]
size[6]
size[c(2,6)]
size[c(6,2)]
min(size[c(6,2)])

# silde 29
names(size)
names(size) <- c("Fabien","Pierre","Sandrine","Claire","Bruno","Delphine")
size
str(size)

my_vector <- c("one"=1, "two"=2, "three"=3)
my_vector
str(my_vector)

#------------------------------------
# interactive tutorials: slides 30-31
# -----------------------------------

# slide 30
#install.packages("swirl") # this command has to be done only once if package not installed in your computer or environment : once a packaged is installed, it is present on the computer/server.

# During the session, it was installed on your home in adenine in this directory:  "/srv/home/login/R/x86_64-conda_cos6-linux-gnu-library/4.0"
# You may check it is well installed by listing the files and folders in this directory:
list.files(.libPaths()[1]) # with '.libPaths()' you  get the folders where R libraries are installed; by adding [1], you get the first of these folders which is your default library path for R and is "/srv/home/login/R/x86_64-conda_cos6-linux-gnu-library/4.0; by including this command inside list.files(), you list all the files in this folder.
library(swirl) # to be done at each session in which you want to use swirl
swirl() # then follow the instructions

# Then try the tutorials 1, 3 to 6.

# slide 31: if you cannot select the tutorial "R programming: the basics of R programming"
# quit swirl
# then in R, install the course R_Programming.swc (alerady present in your megm1_ue5_r environment in adenine) with the following command:
# swirl::install_course()
# and select the course with the browser from your environment directory
# you may then restart swirl and follow the instructions as in slide 30


bye() # after selecting 0 if selection required


#---------------------------------------------------
# saving the data: slide 29
# --------------------------------------------------

# slide 29
save.image(file="RSession1.RData")
save(weight, size, bmi, file="anthropo.RData")

# in a new R session, you can then load them all or only the anthropo data:
#load("anthropo.RData")
#ls()
