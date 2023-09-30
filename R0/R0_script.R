#!/usr/bin/env Rscript

################################################################################
## M1 MEG - UE5
## commandes du diaporama
## Claire Vandiedonck
################################################################################

## Some very simple examples
2 + 2
exp(-2)
log(100,base = 10)
?exp # equivalent to help(exp)

## Assigning data into R objects, using and reading them
x <- 2
y <- x + 3
s <- "this is a string of characters"
x
y
s
x + x
x^y
x <- 4
y
y <- x + 3
y

#Managing objects in your R Session
ls()
rm(y)
ls()
rm(list = ls())# pour tout deleter
ls()


## Managing your session
sessionInfo()
getwd()
#setwd()

# Creating a folder in your working directory
list.files()
dir.create("myfolder")
list.files()
ls()

zz <- file("mynewfile.txt", open = "wt") 
list.files()
close(zz)

sink("myRoutputs.Rout") 
1 + 1
x <- 4
is.numeric(x)
sink()

zz <- file("myRoutputs2.Rout", open = "wt") 
sink(zz)
is.numeric(x)
1:10
sink()
close(zz)

#Saving your session, data and history
save(x, file = "x.RData")
save.image(file = "AllMyData.RData")
savehistory(file = "MyHistory.Rhistory")
#load("x.RData")
#load("AllMyData.RData")

# Running a set of commands = run a script
#source ("myscript.R")


