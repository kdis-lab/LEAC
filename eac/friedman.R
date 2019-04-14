#!/usr/bin/env Rscript
#
# /*! \file friedman.R
#  *
#  * \brief  Friedman test
#  *
#  * \details  This file is part of the LEAC.\n\n
#  * \version 1.0
#  * \date 2015-2017
#  * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
#  * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
#  */
#
#
#  Use:
#       Rscript --vanilla friedman.R rand_index_table.csv > test_metric.log 
#
#
#
# Install scmamp
#
# Install function for packages    
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    if (!require("devtools")) {
       install.packages("devtools")
    }
    devtools::install_github("b0rxa/scmamp")
#    install.packages(pkgs=x,repos="http://cran.r-project.org")
    if (!require(x,character.only=TRUE))
      stop("\nError:\n\tPackage 'scmamp' not found\n\n")
  }
}
#
packages(scmamp)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least two argument must be supplied\n\nUse:\n\tRscript friedman.R rand_index_table.csv metric_name > test_metric.log", call.=FALSE)
}
if (length(args)==1) {
  stop("At least two argument must be supplied\n\nUse:\n\tRscript friedman.R rand_index_table.csv metric_name > test_metric.log", call.=FALSE)
}
library('scmamp')
options(width = 360)
filename  <- args[1]
metricname <-  args[2]
#READ FILE CSV VALUES METRIC
valuemetric <- read.csv(file=filename,header=TRUE,sep=",",row.names = 1)
imanDavenportTest(valuemetric)
nm <- nemenyiTest(valuemetric)
nm
nm$diff.matrix
nemenyplot <- paste(metricname, "_nemeny_plot.pdf", sep="")
pdf(nemenyplot)
plotCD(results.matrix=valuemetric,alpha = 0.05)
#Friedman post-hoc test with Bergmann and Hommel’s correction
#cat("\n\n\nFriedman post-hoc test with Bergmann and Hommel’s correction\n")
#for 9 or less algorithms
#test.res <- postHocTest(data =valuemetric,test = 'friedman',  correct =#'bergmann')
cat("\n\n\nFriedman post-hoc test with Bergmann\n")
test.res <- postHocTest(data =valuemetric,test = 'friedman')
test.res
bold <- test.res$corrected.pval < 0.05
bold[is.na(bold)] <- FALSE
nemenymatrix <- paste(metricname, "_nemeny_matrix.tex", sep="")
sink(nemenymatrix)
writeTabular(table = test.res$corrected.pval, format = 'f', bold=bold,hrule=0, vrule = 0)
cat("\\caption{Friedman post-hoc test with Bergmann for: ", metricname)
cat("}\n")
sink()
#
# On the other hand, the results can be graphically shown in a graph that
# represents the algorithms that show no significant differences as connected nodes.
#
average.ranking <- colMeans(rankMatrix(valuemetric))
rankingfile <- paste(metricname, "_ranking.csv", sep="")
sink(rankingfile)
sort(average.ranking)
cat("Average Rankings of the algorithms: ",metricname)
cat("\n")
sink()
rankinggraph <- paste(metricname, "_ranking_graph_friedman.pdf", sep="")
pdf(rankinggraph)
drawAlgorithmGraph(pvalue.matrix = test.res$corrected.pval,mean.value =average.ranking)
