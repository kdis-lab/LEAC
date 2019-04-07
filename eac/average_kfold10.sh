# #
# use:
# sh ../../eac/average_kfold10.sh directory_run  num_run
#   ej  sh average_kfold10.sh hka_fkmedoid 10
#
if [ "$#" -ne 2 ]; then
  echo
  echo "Usage: $0 directory numberRun" >&2
  exit 1
fi

if [ ! -d "$1" ]; then
  echo
  echo "No existe el directory: $1"
  echo
  echo "Usage: $0 directory numberRun" >&2
fi
echo $1
DIRTOPROC=$1
TMPDIRPROC_OUT="__"$1"_proc"
NUM_RUN=$2
NUM_RUNBYFOLD=$( expr 10 '*' $NUM_RUN )
#a=$( expr 2 '*' "$k" + 1 )
# Absolute path to this script
#
SCRIPT=$(readlink -f $0)
# Absolute path this script
SCRIPTPATH=`dirname $SCRIPT`
#
if [ -d "$TMPDIRPROC_OUT" ]; then
rm -r $TMPDIRPROC_OUT
fi
mkdir $TMPDIRPROC_OUT
#cd $DIRTOPROC
#cd /run/media/hermes/Transcend/ga_run20180515/rawrun/garun_20180809_k6z0; mkdir  ../../proc/fixk_z0k6/kga_fkcentroid_BM02a
find  $DIRTOPROC -name $DIRTOPROC"*_run.csv" -print0 | xargs -0 grep '_inout,out' > $TMPDIRPROC_OUT/_data_all.csv
cd $TMPDIRPROC_OUT
# PASO 1)
awk -F"," -v OFS="," '{for(j = 1; j <= NF; j++) {if ($j == "_k") { ++j; kval = $j}};  for(j = 1; j <= NF; j++) {if ($j == "_dataset") { ++j;  dataset = $j}}; printf("%s/_k%s,%s\n",dataset,kval,$0)}' _data_all.csv | sort > _data_all_sort.csv
#NUMERA.AWK
awk -F"," -v OFS="," -v numrun="$NUM_RUN" 'BEGIN { ant = ""; count = 1; } { if ( $1 == ant ) { count++; } else { count = 1; } if ( count <= numrun ) { print $0 } ant = $1 }' _data_all_sort.csv  > _data_all_sort_filter.csv
#
# #ES DIFERENTE PORQUE SE TIENEN CORRIDAS CON DIFERENTE K PARA EL MISMO DATASET
#PREFIXOUT="_"$DIRTOPROC".csv"
#FILEAVG=$DIRTOPROC"_all_sort_filter_procavg.csv"
perl -X $SCRIPTPATH/partir_kfix.pl _data_all_sort_filter.csv > _data_sort_filter_procavg.csv

# # PASO 2) HACER EL ARCHIVO PARA PARTIR DATOS DEL ALGORTIMO POR DATASET
#awk -F"," -v varfile="$PREFIXOUT"   'BEGIN {gato  = "#"; }{cnt[$2]++}END{for (x in cnt) { if (cnt[x] < 200) printf("%s",gato); gx =  "grep \","x",\"";  print varfile > ./"x"$PREFIXOUT  "gato cnt[x]}}'  | sort > $DIRTOPROC".sh"

#warning: missing runs
#echo "NUM_RUNBYFOLD $NUM_RUNBYFOLD"
awk -v numrun="$NUM_RUNBYFOLD"  -F"," 'BEGIN {gato  = "#"; }{cnt[$2]++}END{for (x in cnt) { if (cnt[x] < numrun) msj = " warning: missing runs"; else msj = ""; gx =  "grep \","x",\"";  print  gx " _data_sort_filter_procavg.csv > ./"x"_split.csv  "gato cnt[x] msj;}}' _data_sort_filter_procavg.csv | sort > splitscript.sh
# # PASO 3) PARTIR POR DATASET
cat splitscript.sh
sh splitscript.sh
# # PASO 4)
#
METRICDIR_OUT=../metric_average
if [ ! -d "$METRICDIR_OUT" ]; then
mkdir $METRICDIR_OUT
fi
FILEOUTMETRIC=$METRICDIR_OUT/$DIRTOPROC"_avg_metric.csv"
if [ -e "$FILEOUTMETRIC" ]; then
rm $FILEOUTMETRIC
fi
perl $SCRIPTPATH/average_fold.pl *_split.csv >> $FILEOUTMETRIC
cd ..
if [ -d "$TMPDIRPROC_OUT" ]; then
rm -r $TMPDIRPROC_OUT
fi
