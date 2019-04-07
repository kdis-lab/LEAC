# #
# use:
#    sh ../../eac/scrip_friedman.sh
#
#
SCRIPT=$(readlink -f $0)
# Absolute path this script
SCRIPTPATH=`dirname $SCRIPT`
#
Rscript $SCRIPTPATH/friedman.R cs_measure_table.csv cs_measure > cs_measure.log
echo "cs_measure OK"
Rscript $SCRIPTPATH/friedman.R db_index_table.csv db_index > db_index.log
echo "db_index OK"
Rscript $SCRIPTPATH/friedman.R db_index_kge2_table.csv db_index_kge2 > db_index_kge2.log
echo "db_index_kge2 OK"
Rscript $SCRIPTPATH/friedman.R distortion_table.csv distortion > distortion.log
echo "distortion OK"
Rscript $SCRIPTPATH/friedman.R dunn_s_index_table.csv dunn_s_index > dunn_s_index.log
echo "dunn_s_index OK"
Rscript $SCRIPTPATH/friedman.R entropy_table.csv entropy > entropy.log
echo "entropy OK"
Rscript $SCRIPTPATH/friedman.R percent_run_geqk2_tst_table.csv percent_run_geqk2_tst > percent_run_geqk2_tst.log
echo "percent_run_geqk2_tst OK"
Rscript $SCRIPTPATH/friedman.R execution_time_seconds_table.csv execution_time_seconds > execution_time_seconds.log
echo "execution_time_seconds OK"
Rscript $SCRIPTPATH/friedman.R f_measure_table.csv f_measure > f_measure.log
echo "f_measure OK"
Rscript $SCRIPTPATH/friedman.R index_i_c_fuzzy_table.csv index_i_c_fuzzy > index_i_c_fuzzy.log
echo "index_i_c_fuzzy OK"
Rscript $SCRIPTPATH/friedman.R index_i_table.csv index_i > index_i.log
echo "index_i OK"
Rscript $SCRIPTPATH/friedman.R iterations_need_table.csv iterations_need > iterations_need.log
echo "iterations_need OK"
Rscript $SCRIPTPATH/friedman.R j1_table.csv j1 > j1.log
echo "j1 OK"
Rscript $SCRIPTPATH/friedman.R jaccard_index_table.csv jaccard_index > jaccard_index.log
echo "jaccard_index OK"
Rscript $SCRIPTPATH/friedman.R jm_table.csv jm > jm.log
echo "jm OK"
Rscript $SCRIPTPATH/friedman.R misclassified_table.csv misclassified > misclassified.log
echo "misclassified OK"
Rscript $SCRIPTPATH/friedman.R precision_table.csv precision > precision.log
echo "precision OK"
Rscript $SCRIPTPATH/friedman.R purity_table.csv purity > purity.log
echo "purity OK"
Rscript $SCRIPTPATH/friedman.R rand_index_table.csv rand_index > rand_index.log
echo "rand_index OK"
Rscript $SCRIPTPATH/friedman.R recall_table.csv recall > recall.log
echo "recall OK"
Rscript $SCRIPTPATH/friedman.R score_function_table.csv score_function > score_function.log
echo "score_function OK"
Rscript $SCRIPTPATH/friedman.R sdunn_s_index_table.csv sdunn_s_index > sdunn_s_index.log
echo "sdunn_s_index OK"
Rscript $SCRIPTPATH/friedman.R sed_table.csv sed > sed.log
echo "sed OK"
Rscript $SCRIPTPATH/friedman.R silhouette_table.csv silhouette > silhouette.log
echo "silhouette OK"
Rscript $SCRIPTPATH/friedman.R simplified_silhouette_table.csv simplified_silhouette > simplified_silhouette.log
echo "simplified_silhouette OK"
Rscript $SCRIPTPATH/friedman.R ssb_table.csv ssb > ssb.log
echo "ssb OK"
Rscript $SCRIPTPATH/friedman.R sse_table.csv sse > sse.log
echo "sse OK"
Rscript $SCRIPTPATH/friedman.R time_seconds_need_best_table.csv time_seconds_need_best > time_seconds_need_best.log
echo "time_seconds_need_best OK"
Rscript $SCRIPTPATH/friedman.R variance_ratio_criterion_table.csv variance_ratio_criterion > variance_ratio_criterion.log
echo "variance_ratio_criterion OK"
Rscript $SCRIPTPATH/friedman.R wb_index_table.csv wb_index > wb_index.log
echo "wb_index OK"
Rscript $SCRIPTPATH/friedman.R xie_beni_index_crisp_table.csv xie_beni_index_crisp > xie_beni_index_crisp.log
echo "xie_beni_index_crisp OK"
Rscript $SCRIPTPATH/friedman.R xie_beni_index_table.csv xie_beni_index > xie_beni_index.log
echo "xie_beni_index OK"
