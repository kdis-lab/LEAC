#!/bin/bash
# /*! \file run_alg_kfold.sh
#  *
#  * \brief  k-fold 
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
#       parallel -j10 ./dataset_proc.sh ::: {1..10} > log.log4 &
#
#
#MODIFIQUE Y ESCRIBIR EN DONDE SE ENCUENTRAN LOS BINARIOS 
PATH_AALGORITHMS="~/leac-update/leac-master/eac"
#MODIFIQUE Y ESCRIBA DONDE SEENECUENTRAN LOS 
PATH_DATASETS="~/leac-update/leac-master/data"
#PATH_DATASETS="~/dataset"
#
#-b  --format-file[=NAME]    uci, or keel, by default uci
ADATASET=(
"iris_z5"        "-b uci -h yes -a \"1-4\"  -c 5  --number-clusters=3"
"wine_z5"        "-b uci -h yes -a \"1-12\" -c 13 --number-clusters=3"
"sonar_z5"       "-b uci -h yes -a \"1-60\" -c 61 --number-clusters=2"
"glass_z5"       "-b uci -h yes -a \"1-9\"  -c 10 --number-clusters=7"
"spect_z5"      "-b uci -h yes -a \"1-22\" -c 23 --number-clusters=2"
"haberman_z5"    "-b uci -h yes -a \"1-3\"  -c 4  --number-clusters=2"
"ecoli_z5"       "-b uci -h yes -a \"1-7\"  -c 8  --number-clusters=8"
"ionosphere_z5"  "-b uci -h yes -a \"1-33\" -c 34 --number-clusters=2"
"hayes_roth_z5"  "-b uci -h yes -a \"1-3\"  -c 4  --number-clusters=3"
"zoo_z5"         "-b uci -h yes -a \"1-16\" -c 17  --number-clusters=7"
)
#
#NAME ALGORITHMS AND PARAMETERS OF THE ALGORITHMS
AALGORITHMS=( 
"kga_fkcentroid"  ""
"fgka_fklabel"    ""
"hka_fkmedoid"    ""
"gka_fklabel"     ""
"gaclustering_fklabel"  ""
"gagr_fkcentroid" "--generations=200"
#"gas_fkcentroid"  ""
#"igka_fklabel"    ""
#"gca_fkmedoid"    ""
#"gaclustering_fkcrispmatrix" "--generations=200" 
)
# get number of dataset
NUM_DATASET=${#ADATASET[@]}
# get number algorithms
NUM_ALGORITHMS=${#AALGORITHMS[@]}
#
#
echo "10-FOLD*************************************************************************"
echo "********************************************************************************"
#<<<PARAMETER IN-------------------------------------------------------------
#
VESION_RUN="v4"
#ID_MACHINE=`cat id_machine.dat`
NUM_RUN=1
IFILE="$1"
KFOLD_MOD=10
IFILEPROC=$(( IFILE % KFOLD_MOD  ))
IFILEPROC=$(( IFILEPROC == 0? KFOLD_MOD:IFILEPROC))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
for (( j=0; j<${NUM_ALGORITHMS}; j++ ));
do
NAME_ALGORITHMS=${AALGORITHMS[$j]}
PARAM_ALGORITHMS=${AALGORITHMS[$j+1]}
DIR_OUT=$NAME_ALGORITHMS$"_"$VESION_RUN
#
if [ ! -d $DIR_OUT ]; then
echo $DIR_OUT "File " $DIR_OUT " not exists"
mkdir $DIR_OUT    
fi
cd $DIR_OUT
#
echo
echo "procesando:" $DIR_OUT 
echo
echo $NAME_ALGORITHMS
#
for (( i=0; i<${NUM_DATASET}; i++ ));
do
NAME_PROCESING=${ADATASET[$i]}
PARAM_DATASET=${ADATASET[$i+1]}
NAME_FILE=${ADATASET[$i]}"-10-"
#DIR_OUT=$NAME_ALGORITHMS$VESION_RUN$NAME_PROCESING
#DIR_OUT=$NAME_ALGORITHMS$"_"$VESION_RUN
DIR_DATASET=$PATH_DATASETS"/"$NAME_PROCESING
echo "DIR_DATASET = " $DIR_DATASET

for IRUN in $(seq 1 $NUM_RUN)
do
NUMHAVE=0
#PREFIXFILE=$NAME_ALGORITHMS"_"$NAME_PROCESING"_"$ID_MACHINE$IFILE$IRUN
PREFIXFILE=$NAME_ALGORITHMS"_"$NAME_PROCESING"_"$IFILE"_"$IRUN
if [ -e $PREFIXFILE"_run.csv" ]; then
NUMHAVE=$(grep "_inout,out" $PREFIXFILE"_run.csv" |  wc -l)
fi
if [ "$NUMHAVE" -eq "0" ]; then
FILE_PROC=$DIR_DATASET"/"$NAME_FILE$IFILEPROC
comando="$PATH_AALGORITHMS/$NAME_ALGORITHMS -i $FILE_PROC\"tra.dat\" -t $FILE_PROC\"tst.dat\" $PARAM_DATASET  $PARAM_ALGORITHMS -r 1 -R $PREFIXFILE\"_run.csv\" -C $PREFIXFILE\"_centroids.dat\" -M $PREFIXFILE\"_membership.dat\" -T $PREFIXFILE\"_partitionstable.dat\"&"
#comando="$PATH_AALGORITHMS/$NAME_ALGORITHMS -i $FILE_PROC\"tra.dat\" -t $FILE_PROC\"tst.dat\" $PARAM_DATASET  $PARAM_ALGORITHMS -b keel -r 1  -q -R $PREFIXFILE\"_run.csv\" -C $PREFIXFILE\"_centroids.dat\" -M $PREFIXFILE\"_membership.dat\" -T $PREFIXFILE\"_partitionstable.dat\" --gnuplot=$PREFIXFILE&"
echo $comando
eval $comando
wait
fi
done
((i++))
done #DATASET
cd ..
echo "process finalzado" $NAME_PROCESING
#cat  *_run.csv >> ../$NAME_ALGORITHMS"_metric_all.csv"
((j++))
done #ALGO
