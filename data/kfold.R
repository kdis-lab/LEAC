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
#
args = commandArgs(trailingOnly=TRUE)
numk <-10
#
if (length(args)==2) {
  numk <- as.numeric(args[2])
}
filename <- args[1]
vecfilename <- strsplit(args[1], "\\.")[[1]]
dirname <- vecfilename[1]
if (!dir.exists(dirname)) {
  dir.create(dirname)
}
# load the CSV file from the local directory
dataset <- read.csv(filename, header=TRUE)

dataset<-dataset[sample(nrow(dataset)),]
#Create numk equally size folds
folds <- cut(seq(1,nrow(dataset)),breaks=numk,labels=FALSE)
#Perform numk fold cross validation
for(i in 1:numk){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  data_test <- dataset[testIndexes, ]
  data_train <- dataset[-testIndexes, ]
  #Use the test and train data partitions however you desire...
  filenamett <- paste(dirname, "/", sep="") 
  filenamett <- paste(filenamett,dirname, sep="") 
  filenamett <- paste(filenamett, "-", sep="")
  filenamett <- paste(filenamett, numk, sep="")
  filenamett <- paste(filenamett, "-", sep="")
  filenamett <- paste(filenamett,toString(i), sep="")
  filenametra <- paste(filenamett,"tra.dat",sep="")
  filenametst <- paste(filenamett,"tst.dat",sep="")
  write.csv(data_train,file=filenametra,row.names=FALSE,quote=F)
  write.csv(data_test,file=filenametst,row.names=FALSE,quote=F)
}
#end kfold.R
