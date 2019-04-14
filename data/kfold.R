# /*! \file kfold.R
#  *
#  * \brief  k-fold 
#  *
#  * \details  This file is part of the LEAC.\n\n
#  * \version 1.0
#  * \date 2015-2017
#  * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
#  * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
#  */
#  Use:
#       Rscript  kfold.R dataset [kfold]
#         Parameters:
#           dataset with header
#           kfold   an integer, default kfold is 10
#
#       For example: Rscript kfold.R iris_z5.data
#!/usr/bin/env Rscript
#
# Install caret
#
# Install function for packages    
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    if (!require(x,character.only=TRUE))
      stop("\nError:\n\tPackage 'caret' not found\n\n")
  }
}
#
packages(caret)
#
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Use:\n Rscript kfold.R dataset [k]\n\n", call.=FALSE)
}
#
numk <-10
#
if (length(args)==2) {
  numk <- as.numeric(args[2])
}
suppressPackageStartupMessages(library(caret))
filename <- args[1]
vecfilename <- strsplit(args[1], "\\.")[[1]]
dirname <- vecfilename[1]
if (!dir.exists(dirname)) {
dir.create(dirname)
}
# load the CSV file from the local directory
dataset <- read.csv(filename, header=TRUE)
ndataset <- nrow(dataset)
ndatatest <- ndataset / numk
#define train/test split of the dataset
split = 1 - ndatatest/ ndataset
for (i in 1:numk) {
filenamett <- paste(dirname, "/", sep="") 
filenamett <- paste(filenamett,dirname, sep="") 
filenamett <- paste(filenamett, "-", sep="")
filenamett <- paste(filenamett, numk, sep="")
filenamett <- paste(filenamett, "-", sep="")
filenamett <- paste(filenamett,toString(i), sep="")
filenametra <- paste(filenamett,"tra.dat",sep="")
filenametst <- paste(filenamett,"tst.dat",sep="")
trainIndex <- createDataPartition(dataset[,ncol(dataset)], p=split, list=FALSE)
data_train <- dataset[ trainIndex,]
data_test <- dataset[-trainIndex,]
write.csv(data_train,file=filenametra,row.names=FALSE,quote=F)
write.csv(data_test,file=filenametst,row.names=FALSE,quote=F)
}
#end kfols.r
