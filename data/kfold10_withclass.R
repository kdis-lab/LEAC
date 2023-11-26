# /*! \file kfold10_withclass.R
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
#       Rscript  --vanilla kfold10_withclass.R iris_z5.dat 
#
#
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
library(caret)
filename <- args[1]
vecfilename <- strsplit(args[1], "\\.")[[1]]
dirname <- vecfilename[1]
if (!dir.exists(dirname)) {
dir.create(dirname)
}
# load the CSV file from the local directory
dataset <- read.csv(filename, header=TRUE)
#
#
#define an 90%/10% train/test split of the dataset
split=0.90
for (i in 1:10) {
filenamett <- paste(dirname, "/", sep="") 
filenamett <- paste(filenamett,dirname, sep="") 
filenamett <- paste(filenamett, "-10-", sep="")
filenamett <- paste(filenamett,toString(i), sep="")
filenametra <- paste(filenamett,"tra.dat",sep="")
filenametst <- paste(filenamett,"tst.dat",sep="")

trainIndex <- createDataPartition(dataset[,ncol(dataset)], p=split, list=FALSE)
data_train <- dataset[ trainIndex,]
data_test <- dataset[-trainIndex,]
write.csv(data_train,file=filenametra,row.names=FALSE,quote=F)
write.csv(data_test,file=filenametst,row.names=FALSE,quote=F)
}

