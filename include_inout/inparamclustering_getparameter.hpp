 /*! \file inparamclustering_getparameter.hpp
 *
 * \brief  inparam clustering get parameter
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_CLUSTERING_GETPARAMETER_HPP__
#define __IN_PARAM_CLUSTERING_GETPARAMETER_HPP__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <string.h>
#include <getopt.h>
#include "vector_utils.hpp"
#include "file_utils.hpp"

extern char  *optarg;

#include "verbose_global.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace  inout {

/*---< inparam_usageClustering>------------------------------------------------------------*/

#ifdef __INPARAM_KMEANS__  
  template <typename T_CLUSTERIDX,
	    typename T_FEATURE,         
	    typename T_FEATURE_SUM,
	    typename T_INSTANCES_CLUSTER_K
	    > 
void
inparamclustering_usage
(char                    *argv0, 
 InParamKmeans
 <T_CLUSTERIDX,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>  &aoipc_inParamClustering
 )
#endif /* __INPARAM_KMEANS__ */

  
#ifdef  __INPARAM_SLINK__
template <typename T_CLUSTERIDX,
	  typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K
	 > 
void
inparamclustering_usage
(char                                   *argv0, 
 InParamSlink<T_CLUSTERIDX,T_FEATURE,T_INSTANCES_CLUSTER_K>  &aoipc_inParamClustering
)
#endif /*__INPARAM_SLINK__*/
  

#ifdef  __INPARAM_PCA__ 
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
  void
  inparamclustering_usage
  (char         *argv0, 
   InParamPCA<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>   &aoipc_inParamClustering
   )
#endif /* __INPARAM_PCA__ */

#ifdef  __INPARAM_PCATRANSMATRIX__
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   > 
  void
  inparamclustering_usage
  (char                   *argv0, 
   InParamPCAtransmatrix<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>  &aoipc_inParamClustering
   )
#endif /* __INPARAM_PCATRANSMATRIX__ */

#ifdef  __INPARAM_STDVAR_MILLIGAN_COOPER1988__
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
  void
  inparamclustering_usage
  (char            *argv0, 
   InParamStdVar<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>   &aoipc_inParamClustering
   )
#endif /* __INPARAM_STDVAR_MILLIGAN_COOPER1988__ */
  
#ifdef  __INPARAM_PLOT_CLUSTERING__ 
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	  > 
  void
  inparamclustering_usage
(char                    *argv0, 
 InParamPlotClustering<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>   &aoipc_inParamClustering
 )
#endif /* __INPARAM_PLOT_CLUSTERING__ */


#ifdef __INPARAM_PAM__
template <typename T_CLUSTERIDX,
          typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K
          > 
  void
  inparamclustering_usage
  (char                         *argv0, 
   InParamPAM
   <T_CLUSTERIDX,
   T_FEATURE,
   T_INSTANCES_CLUSTER_K>                   &aoipc_inParamClustering
   )
#endif /*__INPARAM_PAM__*/

#ifdef  __INPARAM_ISODATA__
  template <typename T_CLUSTERIDX,
            typename T_EPSILON,
            typename T_FEATURE,         
	    typename T_FEATURE_SUM,
	    typename T_INSTANCES_CLUSTER_K
	    > 
void 
inparamclustering_usage
(char                       *argv0, 
 InParamIsoData
 <T_CLUSTERIDX,
 T_EPSILON,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>     &aoipc_inParamClustering
 )
#endif /*__INPARAM_ISODATA__*/

#ifdef  __INPARAM_FCM__
  template <typename T_CLUSTERIDX,
            typename T_U,
            typename T_FEATURE,         
	    typename T_INSTANCES_CLUSTER_K
	    >
  void 
  inparamclustering_usage
  (char                    *argv0, 
   InParamFCM
   <T_CLUSTERIDX,
   T_U,
   T_FEATURE,
   T_INSTANCES_CLUSTER_K>  &aoipc_inParamClustering
   )  
#endif /*__INPARAM_FCM__*/

#ifdef __INPARAM_DBSCAN__
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	  >     
  void 
  inparamclustering_usage
  (char *argv0, 
   InParamDBSCAN<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>  &aoipc_inParamClustering
   )

#endif /*__INPARAM_DBSCAN__*/


#ifdef __INPARAM_PCPMFK__    
template<typename T_CLUSTERIDX,
  typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   > 
void 
inparamclustering_usage
(char *argv0, 
 InParamPcPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K
 >                             &aoipc_inParamClustering
) 
#endif /* __INPARAM_PCPMFK__ */

#ifdef __INPARAM_PCPMVK__
template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_usage
(char                                  *argv0, 
 InParamPcPmVk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                &aoipc_inParamClustering
 )
#endif /*__INPARAM_PCPMVK__
	*/

#ifdef __INPARAM_GAPROTOTYPESFK__
template<typename T_BITSIZE,
         typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
  	 > 
void 
inparamclustering_usage
(char *argv0, 
 InParamGAPrototypesFk
 <T_BITSIZE,
 T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K
 >
 &aoipc_inParamClustering
 )
#endif /*__INPARAM_GAPROTOTYPESFK__
	*/

  
#ifdef __INPARAM_TGCA__
template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
         > 
  void 
  inparamclustering_usage
  (char                      *argv0, 
   InParamTGCA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>     &aoipc_inParamClustering
   )  
#endif /*__INPARAM_TGCA__
	*/

#ifdef __INPARAM_GGA__
  template<typename T_CLUSTERIDX,
            typename T_REAL,
            typename T_FEATURE,
	    typename T_FEATURE_SUM,
	    typename T_INSTANCES_CLUSTER_K
	   > 
  void 
  inparamclustering_usage
  (char                      *argv0, 
   InParamGGA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>  &aoipc_inParamClustering
   )
#endif /* __INPARAM_GGA__ */
    
#ifdef __INPARAM_SUBCLUSTERBINARYVK__
  template<typename T_REAL,
           typename T_BITSIZE,
           typename T_CLUSTERIDX,
           typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
          > 
  void 
  inparamclustering_usage
  (char                              *argv0, 
   InParamSubClusterBinaryVk
   <T_REAL,
   T_BITSIZE,
   T_CLUSTERIDX,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>            &aoipc_inParamClustering
   )
#endif /*__INPARAM_SUBCLUSTERBINARYVK__*/

 
#ifdef __INPARAM_GENWOCHGVK__
  template <typename T_BITSIZE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
  void 
  inparamclustering_usage
  (char *argv0, 
   InParamGenWOChgVk
   <T_BITSIZE,T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
   &aoipc_inParamClustering
   )
#endif /*__INPARAM_GENWOCHGVK__*/
  

#ifdef __INPARAM_FEAC__
template<typename T_FEATURE,
	 typename T_REAL,
	 typename T_CLUSTERIDX,
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_usage
(char *argv0, 
 InParamFEAC
 <T_FEATURE,T_REAL,T_CLUSTERIDX,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
 &aoipc_inParamClustering
 )
#endif /*__INPARAM_FEAC__*/


#ifdef __INPARAM_WITHOUTPCPMFK__
template<typename T_CLUSTERIDX,
	 typename T_BITSIZE,
	 typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_usage
(char *argv0, 
 InParamWithoutPcPmFk
 <T_CLUSTERIDX,T_BITSIZE,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
 &aoipc_inParamClustering
 )
#endif /*__INPARAM_WITHOUTPCPMFK__*/

#ifdef  __INPARAM_GCA__
  template<typename T_CLUSTERIDX,
	   typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   > 
void 
inparamclustering_usage
  (char *argv0, 
   InParamGCA
   <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
   &aoipc_inParamClustering
   )
#endif /*__INPARAM_GCA__*/


#ifdef  __INPARAM_CBGA__
template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
         typename T_INSTANCES_CLUSTER_K,
         typename T_INSTANCE_FREQUENCY> 
  void 
  inparamclustering_usage
  (char                             *argv0, 
   InParamCBGA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K,
   T_INSTANCE_FREQUENCY>           aoipc_inParamClustering
   )
#endif /*__INPARAM_CBGA__*/

#ifdef  __INPARAM_HKA__
template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_usage
(char *argv0, 
InParamHKA
<T_CLUSTERIDX,
T_REAL,
T_FEATURE,
T_FEATURE_SUM,
T_INSTANCES_CLUSTER_K>              &aoipc_inParamClustering
)
#endif /* __INPARAM_HKA__ */


#ifdef __INPARAM_PMFK__

template<typename T_CLUSTERIDX,
	 typename T_REAL,
	 typename T_FEATURE,         
	 typename T_FEATURE_SUM,
         typename T_INSTANCES_CLUSTER_K
	> 
void 
  inparamclustering_usage
  (char *argv0, 
   InParamPmFk
   <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
   &aoipc_inParamClustering
   )
#endif /* __INPARAM_PMFK__
	*/

#ifdef  __INPARAM_ADAPTIVEPCPMFK__
template<typename T_CLUSTERIDX,
	 typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
        > 
  void 
  inparamclustering_usage
  (char *argv0, 
   InParamAdaptivePcPmFk
   <T_CLUSTERIDX,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
   &aoipc_inParamClustering
   )
#endif /*__INPARAM_ADAPTIVEPCPMFK__*/

{ 
  const char   *larray_opFormatFile[] = INPARAMCLUSTERING_FORMATINSTANCEFILE;
  const char   *las_optTokensYesNo[] = {"no", "yes", (char *) NULL };
  
  uintidx li_i;

  std::cout << "\nUsage: "
	    <<  argv0 
	    << " [OPTION]\n" 
	    << "\tAbout groups of instances for a set K as well as statistics\n"
	    << "\tof the algorithm used\n\n";
  std::cout << "  -i, --instances=FILE or DIRECTORY\n"
	    << "                              file or directory containing data of instances\n"
	    << "                                to be clustered\n";

#ifndef __NOT_TEST_FILES__
  
  std::cout << "  -x, --select-instances[=PREFIX]\n"
	    << "                              if instances is directory search files with\n"
	    << "                                prefix for training (eg. iris-10-1tra.dat,\n"
	    << "                                iris-10-2tra.dat,... PREFIX=tra.dat)\n";
  std::cout << "  -t, --test[=FILE or PREFIX] if instances is directory search files with\n"
            << "                                prefix for test (eg. iris-10-1tst.dat,\n"
	    << "                                iris-10-2tst.dat,... PREFIX=tst.dat),\n"
	    << "                                in other case only name file\n";
 
#endif /*__NOT_TEST_FILES__*/
  
  std::cout << "  -b  --format-file[=NAME]    ";
  li_i = 0;
  while (larray_opFormatFile[li_i+1] !=  NULL) {
    std::cout << larray_opFormatFile[li_i] << ",";
    ++li_i;
  } 
  std::cout << " or " 
	    << larray_opFormatFile[li_i] << ", by default " 
	    << larray_opFormatFile[aoipc_inParamClustering.getFormatInstanceFile()]
	    << '\n';
  std::cout << "  -h, --with-header[=yes/no]  file contains names of instances or a header,\n"
	    << "                                by default is "
	    << las_optTokensYesNo[aoipc_inParamClustering.getHaveHeaderFileInstance()]
	    << '\n';
  std::cout << "  -u, --number-instances[=NUMBER]\n"
	    << "                              the number of instances the file contains\n" 
	    << "                                instances, if not specified file is obtained\n";
  std::cout << "  -a, --select-attributes[=ARG]\n"
	    << "                              select the attributes to be processed for\n"
	    << "                                example, \"1-2,4\" by default all. Also\n"
	    << "                                used to specify the  number of dimensions\n"
	    << "                                of the instances, unless specified file is\n"
	    << "                                obtained instances\n";
  std::cout << "  -d, --delimit-attributes=[ARG]\n"
	    << "                              separated file by default \""
	    <<  INPARAM_FILE_SEPARATOR_DEFAULT
	    << "\"\n";
  std::cout << "  -c, --class-column[=NUMBER] input file of instances has a class assigned\n"
	    << "                                in the column [NUMBER=undefined]\n";
  std::cout << "  -e, --cluster-column[=NUMBER]\n"
	    << "                              input file of instances has a cluster assigned\n"
	    << "                                in the column [NUMBER=undefined]\n";
  std::cout << "  -l, --idinstances-column[=NUMBER]\n"
	    << "                              the input file instance is assigned a column\n"
	    << "                                instance identifier [NUMBER=undefined]\n";
  std::cout << "  -f, --freq-instances-column[=NUMBER]\n"
	    << "                              the input file instance is assigned a column\n"
	    << "                                frequency instances [NUMBER=undefined]\n";
  std::cout << "  -r, --number-runs[=NUMBER]  number of runs or repetitions of the algorithm\n"
	    << "                                (by default [NUMBER=1])\n";
  std::cout << "  -R, --runtime-filename=[FILE]\n"
	    << "                              out file of times run\n";
 
  /* ONLY CLUSTERING */
#ifdef __ALG_CLUSTERING__
  const char   *las_opFuncDistance[] = INPARAMCLUSTERING_DISTANCE_TYPE;
  std::cout << "  -n  --distance[=NAME]       ";
  li_i = 0;
  size_t lst_lenline = 0;
  while ( las_opFuncDistance[li_i+1] !=  NULL) {
    lst_lenline += strlen(las_opFuncDistance[li_i]);
    if ( lst_lenline > 40 ) {
       std::cout << "\n                                ";
       lst_lenline = 0;
    }
    std::cout << las_opFuncDistance[li_i] << ", ";
    ++li_i;
  }
  lst_lenline += strlen(las_opFuncDistance[li_i]);
  if ( lst_lenline > 40 ) {
       std::cout << "\n                                ";
       lst_lenline = 0;
  }
  std::cout << "or " 
	    << las_opFuncDistance[li_i]
	    << ",\n                                by default " 
	    << las_opFuncDistance[aoipc_inParamClustering.getOpDistance()]
	    << '\n';
  std::cout << "  -z, --random-seed[=NUMBER]  string with integer number seed by, default\n"
	    << "                                is random\n";
  std::cout << "  -w, --max-execution-time[=NUMBER]\n"
	    << "                              real number for max execution time in seconds\n"
	    << "                                by default is "
	    << aoipc_inParamClustering.getMaxExecutiontime()
	    << "\n";
  
  std::cout << "  -C, --centroids-outfile[=FILE]\n"
	    << "                              print centroids, standard output FILE="
	    <<  OUTFILENAME_STDOUT               
	    << "\n";

  std::cout << "      --centroids-format[=yes/no]\n"
	    << "                              print the matrices by rows and columns,\n"
	    << "                                by default is "
	    << las_optTokensYesNo[aoipc_inParamClustering.getPrintCentroidsFormat()]
	    << '\n';
  
  std::cout << "  -M, --membership-outfile[=FILE]\n"
	    << "                              print membership of the instances,\n"
	    <<"                                 standard output FILE="
	    <<  OUTFILENAME_STDOUT
	    << "\n";
  
  std::cout << "  -T, --partitionstable-outfile[=FILE]\n"
	    << "                              print partitions table of the instances,\n"
	    << "                                standard output FILE="
	    << OUTFILENAME_STDOUT
	    << "\n";

  std::cout << "      --table-format[=yes/no]\n"
	    << "                              print the partitions table by rows and\n"
	    << "                                columns, by default is "
	    << las_optTokensYesNo[aoipc_inParamClustering.getPrintTableFormat()]
	    << '\n';

#ifdef _ALG_GRAPH_BASED_
    
  std::cout << "  -G, --graph-outfile[=FILE] for graph based algorithms, prints the\n"
	    << "                               partition in a graph, standard\n"
	    << "                               output FILE="
	    << OUTFILENAME_STDOUT
	    << "\n";

#endif /*_ALG_GRAPH_BASED_*/

  
#ifndef __WITHOUT_PLOT_STAT
  std::cout << "  -P, --gnuplot=FILE          file of  gnuplot to graphics result\n"
	    << "                                (compiling only with WITHOUT_PLOT_STAT)\n"
	    << "  -y, --gnuplot-styles=WORD   plot graphics with: points, lines,\n"
	    << "                                linespoints, and dot "
	    << "[ARG=" << aoipc_inParamClustering.getGnuPlotCoreStyles() 
	    << "]\n";
#endif /*__WITHOUT_PLOT_STAT*/
 
#endif /*__ALG_CLUSTERING__*/


  std::cout << "\nParticular options of the algorithm "
            << aoipc_inParamClustering.getAlgorithmoName()
	    << '\n'
	    << "  based on "
	    << aoipc_inParamClustering.getAlgorithmoAuthor()
	    << '\n';

#ifdef __INPARAM_PCATRANSMATRIX__
  std::cout << "      --transmatrix-outfile[=FILE]\n"
	    << "                              print transformation matrix,\n"
	    << "                                standard output FILE="
	    <<  OUTFILENAME_STDOUT               
	    << "\n";
#endif /*__INPARAM_PCATRANSMATRIX__*/

#ifdef  __INPARAM_PLOT_CLUSTERING__
  
  const char   *las_opProjection[] = INPARAMCLUSTERING_PLOT_PLOJECTION;
  std::cout << "      --projection[=NAME]    ";
  li_i = 0;
  while (las_opProjection[li_i+1] !=  NULL) {
    std::cout << las_opProjection[li_i] << ",";
    ++li_i;
  } 
  std::cout << " or " 
	    << las_opProjection[li_i] << ", by default " 
	    << las_opProjection[aoipc_inParamClustering.getOpProjection()]
	    << '\n';
  std::cout << "      --centroids-infile[=FILE]\n"
	    << "                              centroids cluster, FILE=" << OUTFILENAME_STDOUT
	    << "\n";
  std::cout << "      --centroids-title[=TEXT]\n"
	    << "                              title for centroids (by default\n"
	    << "                                \"" << aoipc_inParamClustering.getTitleCentroids() << "\")\n";
  
  std::cout << "      --centroids-id=[yes/no] centroids with id, by default is yes\n";
  std::cout << "      --centroids-point[=NUMBER]\n"
	    << "                              point type, by default cicle "
	    << aoipc_inParamClustering.getPointTypeCentroids()
	    << "\n";
  std::cout << "      --centroids-size[=NUMBER]\n"
	    << "                              point size, by default "
	    << aoipc_inParamClustering.getPointSizeCentroids()
	    << "\n";  
  std::cout << "      --centroids-infile2[=FILE]\n"
	    << "                              centroids cluster2, FILE=" << OUTFILENAME_STDOUT 
	    << "\n";
  
  std::cout << "      --centroids-title2[=TEXT]\n"
	    << "                              label for centroids2 (by default "
	    <<  "\"" << aoipc_inParamClustering.getTitleCentroids2() << "\")\n";

  std::cout << "      --centroids-id2=[yes/no]\n"
	    <<  "                             centroids with id, by default is yes\n";
  std::cout << "      --centroids-point2[=NUMBER]\n"
	    << "                              point type, by default cicle "
	    << aoipc_inParamClustering.getPointTypeCentroids()
	    << "\n";
  std::cout << "      --centroids-size2[=NUMBER]\n"
	    << "                              point size, by default "
	    << aoipc_inParamClustering.getPointSizeCentroids()
	    << "\n";
  
  std::cout << "      --member-infile[=FILE]  member cluster, FILE=" << OUTFILENAME_STDOUT
	    << "\n";
  std::cout << "      --member-label[=TEXT]   label for member cluster (by default "
	    <<  "\"" << aoipc_inParamClustering.getLabelMemberCluster() << "\")\n";
  std::cout << "      --member-size[=NUMBER]  point size for member cluster (by default "
	    <<  "\"" << aoipc_inParamClustering.getPointSizeMemberCluster() << "\")\n";

  std::cout << "      --id-increment[=NUMBER] to the identification labels of members\n"
	    << "                                and centroids add a number by default\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getIDincrement()
	    << "]\n";
  
  std::cout << "      --graph-infile[=FILE]   graph, FILE=" << OUTFILENAME_STDOUT // << ",\n"
	    << "\n";
  std::cout << "      --x-coord[=NUMBER]      dimension number as used to coordinate x\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getCoordX()
	    << "]\n";
  std::cout << "      --y-coord[=NUMBER]      dimension number as used to coordinate y,\n"
	    << "                                0 not used [NUMBER="
	    << aoipc_inParamClustering.getCoordY()
	    << "]\n";
  std::cout << "      --z-coord[=NUMBER]      dimension number as used to coordinate z,\n"
	    << "                                0 not used [NUMBER="
	    << aoipc_inParamClustering.getCoordZ()
	    << "]\n";
  std::cout << "      --title[=TEXT]          the chart title and location (eg \"\"my title\"\n"
	    << "                                offset 0,-2\")\n";
  std::cout << "      --legend[=TEXT]         legend options (eg hide off \"{left bottom,\n"
	    << "                                right top }\" samplen 1 spacing 1.2 width -4\n"
	    << "                                maxrows 20 font \",8\"\n";
  std::cout << "      --size-instance[=NUMBER]\n"
	    << "                              instance point size by default [NUMBER="
	    << aoipc_inParamClustering.getSizeInstance()
	    << "]\n";
  std::cout << "      --label-offset[=NUMBER] offset for labels of the instances by default\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getSizeInstance()
	    << "]\n";
  std::cout << "      --x-range[=NUMBER:NUMBER]\n"
	    << "                              range of x to graph, for example:\n"
	    << "                                --x-range=-10:10\n";
  std::cout << "      --y-range[=NUMBER:NUMBER]\n"
	    << "                              range of y to graph, for example:\n"
	    << "                                --y-range=-10:10\n";
  std::cout << "      --color-graphics=[yes/no]\n"
	    << "                              color graphics, by default is yes, Only\n"
	    << "                                for output to a file\n";
  std::cout << "      --sleep[=NUMBER]        number of seconds to wait for the graphic\n"
	    << "                                display [NUMBER="
	    << aoipc_inParamClustering.getSleep()
	    << "]\n";
  std::cout << "      --graphics-outfile[=FILE]\n"
	    << "                              file of out for graphics result, if null\n" 
	    << "                                out terminal x11\n";
#endif /*__INPARAM_PLOT_CLUSTERING__*/

#ifdef __INPARAM_STDVAR_MILLIGAN_COOPER1988__
  const char   *las_opStdVarName[] = STD_VAR_NAME;
  
  std::cout << "      --std-var[=ARG]         standardization of variables: Z0, Z1,\n" 
	    << "                                Z2, Z3, Z4, Z5, Z6, Z7 [ARG=" 
	    << las_opStdVarName[aoipc_inParamClustering.getStandardizationVar()] 
	    << "]\n";
#endif /*__INPARAM_STDVAR_MILLIGAN_COOPER1988__*/
  
#ifdef  __INPARAM_KMEANS__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --iterations[=NUMBER]   maximum number of iterations [NUMBER="
	    << aoipc_inParamClustering.getNumMaxIter()
	    << "]\n";
  std::cout << "      --threshold[=NUMBER]    threshold value [NUMBER="  
	    << aoipc_inParamClustering.getMinThreshold(); 
  std::cout << "]\n";
#endif /* __INPARAM_KMEANS__ */

#ifdef __INPARAM_SLINK__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
#endif /*__INPARAM_SLINK__*/
  
#ifdef  __INPARAM_ISODATA__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --iterations[=NUMBER]   maximum number of iterations\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumMaxIter()
	    << "]\n";
  std::cout << "      --epsilon[=NUMBER]      epsilon value [NUMBER=" 
	    << aoipc_inParamClustering.getEpsilon()
	    << "]\n";
  std::cout << "      --optimal-initializacion\n"
	    << "                              if the instances are assigned class is\n"
	    << "                                initialized with the correct partition\n";
#endif /*__INPARAM_ISODATA__*/    

#ifdef  __INPARAM_FCM__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --iterations[=NUMBER]   maximum number of iterations [NUMBER="
	    << aoipc_inParamClustering.getNumMaxIter()
	    << "]\n";
  std::cout << "      --epsilon[=NUMBER]      epsilon value [NUMBER="
	    << aoipc_inParamClustering.getEpsilon()
	    << "]\n"; 
  std::cout << "      --weighting-exponent[=NUMBER]\n"
	    << "                              m = weighting exponent [NUMBER=" 
	    << aoipc_inParamClustering.getWeightingExponent() 
	    << "], m != 1\n"; 
#endif /*ALG_FCM_BEZDEK_H*/

#ifdef __INPARAM_DBSCAN__
  std::cout << "      --epsilon[=NUMBER]               epsilon neighborhood [NUMBER=" 
	    << aoipc_inParamClustering.getEps()
	    << "]\n";
  std::cout << "      --min-pts[=NUMBER]               minimum number of points [NUMBER=" 
	    << aoipc_inParamClustering.getMinPts() << "]\n";
  std::cout << "\nindependent parameters r * algorithm for tree:\n";
  std::cout << "      --buffer-capacity[=NUMBER]       storage buffer-capacity [NUMBER=" 
	    << aoipc_inParamClustering.getRTreeBufferCapacity() << "]\n";
  std::cout << "      --tree-fill-factor[=NUMBER]      tree fill factor [NUMBER=" 
	    << aoipc_inParamClustering.getRTreeFillFactor() << "]\n";
  std::cout << "      --tree-index-capacity[=NUMBER]   tree index capacity [NUMBER=" 
	    << aoipc_inParamClustering.getRTreeCapacity() << "]\n";
  std::cout << "      --tree-leaf-capacity[=NUMBER]    tree leaf capacity [NUMBER=" 
	    << aoipc_inParamClustering.getRTreeleafCapacity() << "]\n";
  

#endif /*__INPARAM_DBSCAN__*/

#ifdef __INPARAM_PCPMFK__

  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                               [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --crossover-probability[=NUMBER]\n"
	    << "                              real number in the interval [0.25, 1]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                              real number in the interval [0, 0.5]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
#endif /*__INPARAM_PCPMFK__*/

  
#ifdef __INPARAM_TGCA__
  std::cout << "      --k-minimum[=NUMBER]    minimum number of clusters per search\n"
	    << "                                [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterKMinimum() << "]\n";
  std::cout << "      --k-maximum[=NUMBER]    maximum number of clusters per search\n" 
	    << "                                if eq -1, k-maximum = N^1/2 [NUMBER="
	    << aoipc_inParamClustering.getNumClusterKMaximum() << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --num-subpopulations-cross[=NUMBER]\n"
	    << "                              number of subpopulations for parallel\n"
	    << "                                crossover. It must be less than half\n"
	    << "                                the size compared to the population\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumSubpopulationsCross()
	    << "]\n";
  std::cout << "      --crossover-probability[=NUMBER]\n"
	    << "                              real number in the interval [0.25, 1]\n"
	    << "                              [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --kmeans-iterations[=NUMBER]\n"
	    << "                              maximum number of iterations for k-means\n"
	    << "                                algorithm [NUMBER="
	    << aoipc_inParamClustering.getKmeansNumMaxIter()
	    << "]\n";
  std::cout << "      --kmeans-threshold[=NUMBER]\n"
	    << "                              threshold value for k-means algorithm\n"
	    << "                                [NUMBER="  
	    << aoipc_inParamClustering.getKmeansMinThreshold(); 
  std::cout << "]\n";
#endif /*__INPARAM_TGCA__*/

  
#ifdef __INPARAM_GGA__
  std::cout << "    --k-minimum[=NUMBER]      minimum number of clusters per search\n"
	    << "                                [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterKMinimum() << "]\n";
  std::cout << "      --k-maximum[=NUMBER]    maximum number of clusters per search\n" 
	    << "                                if eq -1, k-maximum = N^1/2 [NUMBER="
	    << aoipc_inParamClustering.getNumClusterKMaximum() << "]\n";
  std::cout << "      --sub-population-size[=NUMBER]\n"
	    << "                              size of sub-populations (islands)\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getSubPopulationSize()
	    << "]\n";
  std::cout << "      --number-island[=NUMBER]\n"
	    << "                              number of sub-populations or islands\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumIsland() 
	    << "]\n";
  std::cout << "      --pe[=NUMBER]           probability of migration\n"
	    << "                               good individuals between islands\n"
            << "                               [0,1] [NUMBER="
	    << aoipc_inParamClustering.getPe() 
	    << "]\n";
  std::cout << "      --generations[=NUMBER] number of generations or iterations\n"
	    << "                              [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations() 
	    << "]\n";
  std::cout << "      --pci[=NUMBER]          initial probability crossover, real\n"
	    << "                                number in the interval [0,1] must be\n"
	    << "                                high in the first stages [NUMBER="
	    << aoipc_inParamClustering.getPci()
	    << "]\n";
  std::cout << "      --pcf[=NUMBER]          final probability crossover, real\n"
	    << "                                number in the interval [0,1] must\n"
	    << "                                moderate in the last stages [NUMBER="
	    << aoipc_inParamClustering.getPcf() 
	    << "]\n";
  std::cout << "      --pci[=NUMBER]          initial probability mutation, real number\n"
	    << "                                in the interval [0,1] is smaller in the\n"
	    << "                                first generations [NUMBER="
	    << aoipc_inParamClustering.getPmi() 
	    << "]\n";
  std::cout << "      --pcf[=NUMBER]          final probability mutation, real number in\n"
	    << "                                the interval [0,1] is larger in the last\n"
	    << "                                ones [NUMBER="
	    << aoipc_inParamClustering.getPmf() 
	    << "]\n";  
  std::cout << "      --pbi[=NUMBER]          initial probability local search, real\n"
	    << "                                number in the interval [0,1] must be\n"
	    << "                                high in the first stages [NUMBER="
	    << aoipc_inParamClustering.getPbi()
	    << "]\n";
  std::cout << "      --pbf[=NUMBER]          final probability local search, real\n"
	    << "                                number in the interval [0,1] must\n"
	    << "                                moderate in the last stages [NUMBER="
	    << aoipc_inParamClustering.getPbf() 
	    << "]\n";
#endif 
  
    
#ifdef __INPARAM_SUBCLUSTERBINARYVK__
  std::cout << "      --u-parameter[=NUMBER]    parameter u [NUMBER=" 
	    << aoipc_inParamClustering.getU() << "]\n";
  std::cout << "      --lambda[=NUMBER]         parameter lambda [NUMBER=" 
	    << aoipc_inParamClustering.getLambda() << "]\n";
  std::cout << "      --w1[=NUMBER]             smallest value w1 [NUMBER=" 
	    << aoipc_inParamClustering.getW1() << "]\n";
  std::cout << "      --w2[=NUMBER]             largest value w2 [NUMBER=" 
	    << aoipc_inParamClustering.getW2() << "]\n";
  std::cout << "      --generations[=NUMBER]    number of generations or iterations\n"
	    << "                                  [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
            << "                               size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --crossover-probability[=NUMBER]\n"
	    << "                               real number in the interval [0.25, 1]\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                               real number in the interval [0, 0.5]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
#endif /*__INPARAM_SUBCLUSTERBINARYVK__*/

#ifdef __INPARAM_GENWOCHGVK__

  std::cout << "      --notchangestop[=NUMBER]\n"
	    << "                              after a number x of iteratios,\n"
	    << "                                the best chromosome does\n"
            << "                                not change. (We have fixed x = 3)\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumNotChangeStop() << "].\n";
  
  std::cout << "      --generations[=NUMBER]  maximum number of generations is\n"
	    << "                                reached. We chose this number to\n"
	    << "                                be N if equal to 0, where N the\n"
	    << "                                number of instances. [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              if 0 initial population of N x 10,\n"
	    << "                                where N is the number of instances\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --crossover-probability[=NUMBER]\n"
	    << "                              real number in the interval [0.25, 1]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                              real number in the interval [0, 0.5]\n"
	    << "                              [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
#endif /*__INPARAM_GENWOCHGVK__*/


  
#ifdef __INPARAM_PCPMVK__
  std::cout << "      --k-minimum[=NUMBER]    minimum number of clusters per search\n"
	    << "                                [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterKMinimum() << "]\n";
  std::cout << "      --k-maximum[=NUMBER]    maximum number of clusters per search\n"
            << "                                NUMBER = -1, depends on the algorithm,\n"
	    << "                                for example  N/2 or N^1/2. [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterKMaximum() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --crossover-probability[=NUMBER]\n"
	    << "                              real number in the interval [0.25, 1]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                              for operator 1 and 2 real number in the\n"
	    << "                                interval [0, 0.5] [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
#endif /*__INPARAM_PCPMVK__*/


#ifdef __INPARAM_GAPROTOTYPESFK__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population either 10 or 20\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --crossover-probability[=NUMBER]\n"
	    << "                              real number in the interval [0.25, 1]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                              real number in the interval [0, 0.5]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
  std::cout << "      --probability-ini[=NUMBER]\n"
	    << "                              real number in the interval [0.01, 0.1]\n"
	    << "                                depending n and k,if -1 is eq k/n\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getPini()
	    << "]\n";
  std::cout << "      --alpha[=NUMBER]        real number alpha > 0 to force the algorithm\n"
	    << "                                desired number cluster [NUMBER="
	    << aoipc_inParamClustering.getAlpha()
	    << "]\n";
  
#endif /*__INPARAM_GAPROTOTYPESFK__
	*/

#ifdef __INPARAM_FEAC__
  std::cout << "      --k-minimum[=NUMBER]    minimum number of clusters per search\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getNumClusterKMinimum() << "]\n";
  std::cout << "      --k-maximum[=NUMBER]    maximum number of clusters per search,\n"
	    << "                                 if eq -1, k-maximum = N^1/2 [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterKMaximum() << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --desiable-objfunc[=NUMBER]\n"
	    << "                              value desiable of objetive function\n"
	    << "                                (eg. silhouette[-1,1] rand index[0,1])\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getDesiableObjetiveFunc()
	    << "]\n";
  std::cout << "      --kmeans-iterations[=NUMBER]\n"
	    << "                              maximum number of iterations for k-means\n"
	    << "                                algorithm [NUMBER="
	    << aoipc_inParamClustering.getKmeansNumMaxIter()
	    << "]\n";
  std::cout << "      --kmeans-difference[=NUMBER]\n"
	    << "                              maximum absolute difference between centroids\n"
	    << "                                in two consecutive iterations is less than or\n"
	    << "                                equal to [NUMBER="
	    << aoipc_inParamClustering.getKmeansMaxDiffCent()
	    << "]\n";
#endif /*__INPARAM_FEAC__
       */


#ifdef __INPARAM_GCA__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or\n" 
	    << "                              iterations [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population either 10 or 20\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --mix-recombination-probability[=NUMBER] real number\n" 
	    << "                              in the interval [0.5, 0.9] [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --point-mutation-probability[=NUMBER]  real number\n"  
	    << "                              in the interval [0.2, 0.4] [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
  std::cout << "      --mix-mutation-probability[=NUMBER] real number in the\n" 
	    << "                              interval [0, 0.125] [NUMBER="
	    << aoipc_inParamClustering.getProbMixMutation()
	    << "]\n";
#endif /*__INPARAM_GCA__ 
       */

#ifdef  __INPARAM_HKA__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                               [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                               size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --mix-recombination-probability[=NUMBER]\n"
	    << "                               real number in the interval [0.5, 0.9]\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getProbCrossover() 
	    << "]\n";
  std::cout << "      --point-mutation-probability[=NUMBER]\n"
	    << "                               real number in the interval [0.2, 0.4]\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
  std::cout << "      --mix-mutation-probability[=NUMBER]\n"
	    << "                               real number in the interval [0, 0.125]\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getProbMixMutation()
	    << "]\n";
  std::cout << "      --order-tournament[=NUMBER]\n"
	    << "                               order of tournament [NUMBER="
	    << aoipc_inParamClustering.getOrderTournament()
	    << "]\n";
  std::cout << "      --nearest-neighbors[=NUMBER]\n"
	    << "                               number of the nearest neighbors (p)\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getNearestNeighbors()
	    << "]\n";
  std::cout << "      --search-heuristic-probability[=NUMBER]\n"
	    <<  "                              real number in the interval [0.0, 1.0]\n"
	    << "                                 [NUMBER="
	    << aoipc_inParamClustering.getProbSearchHeuristic()
	    << "]\n";
#endif /*__INPARAM_HKA__*/


#ifdef  __INPARAM_CBGA__
  const char   *las_opSelectMethod[] = INPARAMCLUSTERING_CBGA_SELECMETH;

  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                               [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                               size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                              real number in the interval [0, 1.0]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
 
  std::cout << "      --select-method[=NAME] ";
  li_i = 0;
  while (las_opSelectMethod[li_i+1] !=  NULL) {
    std::cout << las_opSelectMethod[li_i] << ", ";
    ++li_i;
  } 
  std::cout << "\n                                         or " 
	    << las_opSelectMethod[li_i] << ", by default " 
	    << las_opSelectMethod[aoipc_inParamClustering.getOpSelectMethod()]
	    << '\n';
  std::cout << "      --gla-iterations[=NUMBER] number of GLA iterations [NUMBER="
	    << aoipc_inParamClustering.getNumGLAIterations()
	    << "]\n";
#endif  /*__INPARAM_CBGA__*/


#ifdef __INPARAM_WITHOUTPCPMFK__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                               [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                               size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --matingpool-size[=NUMBER]\n"
	    << "                               size of mating pool [NUMBER="
	    << aoipc_inParamClustering.getSizeMatingPool() 
	    << "]\n";
#endif  /*__INPARAM_WITHOUTPCPMFK__*/

#ifdef __INPARAM_PMFK__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                               [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                               size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
  std::cout << "      --mutation-probability[=NUMBER]\n"
	    << "                              real number in the interval [0, 1.0]\n"
	    << "                                [NUMBER="
	    << aoipc_inParamClustering.getProbMutation()
	    << "]\n";
#endif /*__INPARAM_PMFK__*/

#ifdef  __INPARAM_ADAPTIVEPCPMFK__
  std::cout << "      --number-clusters[=NUMBER]\n"
	    << "                              number of clusters [NUMBER=" 
	    << aoipc_inParamClustering.getNumClusterK() << "]\n";
  std::cout << "      --generations[=NUMBER]  number of generations or iterations\n"
	    << "                               [NUMBER="
	    << aoipc_inParamClustering.getNumMaxGenerations()
	    << "]\n";
  std::cout << "      --population-size[=NUMBER]\n"
	    << "                              size of population [NUMBER="
	    << aoipc_inParamClustering.getSizePopulation()
	    << "]\n";
#endif  /*__INPARAM_ADAPTIVEPCPMFK__*/


  std::cout << "\n";
  std::cout << "  -v, --verbose[=NUMBER]      explain what is being done (compiled with\n"
	    << "                                VERBOSE=yes)\n";
#ifdef __VERBOSE_YES
  std::cout << "                              NUMBER=[-1,..,9999] Quiet level -1 not,\n"
	    << "                                verbose, default="
	    <<  geiinparam_verboseMax << "\n";
#endif /*__VERBOSE_YES*/
  
  std::cout << "  -q, --bar-progress          progress bar printing,  default is not\n";

  std::cout << "  -?, --help                  help";
  std::cout <<  std::endl;
  
  exit(-1);
}


#ifdef  __INPARAM_STDVAR_MILLIGAN_COOPER1988__ 
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
void 
inparamclustering_getParameter
(InParamStdVar<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aoipc_inParamClustering,
 int                 argc, 
 char                **argv
 )
#endif /* __INPARAM_STDVAR_MILLIGAN_COOPER1988__ */
  

#ifdef  __INPARAM_PCA__

#include "inparam_pca.hpp"
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   > 
void 
inparamclustering_getParameter
  (InParamPCA<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>   &aoipc_inParamClustering,
   int          argc, 
   char         **argv
   )
#endif /* __INPARAM_PCA__ */

#ifdef  __INPARAM_PCATRANSMATRIX__
  template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   > 
  void 
  inparamclustering_getParameter
  (InParamPCAtransmatrix<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>   &aoipc_inParamClustering,
   int          argc, 
   char         **argv
   )
#endif /*__INPARAM_PCATRANSMATRIX__*/
  
  
#ifdef  __INPARAM_PLOT_CLUSTERING__  
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	  > 
void 
inparamclustering_getParameter
(InParamPlotClustering<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aoipc_inParamClustering,
 int                   argc, 
 char                  **argv
)
#endif /* __INPARAM_PLOT_CLUSTERING__ */


#ifdef __INPARAM_KMEANS__
  template <typename T_CLUSTERIDX,
            typename T_FEATURE,         
	    typename T_FEATURE_SUM,
	    typename T_INSTANCES_CLUSTER_K
	    >
  void 
  inparamclustering_getParameter
  (InParamKmeans
   <T_CLUSTERIDX,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>    &aoipc_inParamClustering,
   int                       argc, 
   char                      **argv
   )
#endif /*__INPARAM_KMEANS__ */

#ifdef __INPARAM_SLINK__
template <typename T_CLUSTERIDX,
	  typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K
	 >
  void 
  inparamclustering_getParameter
  (InParamSlink
   <T_CLUSTERIDX,
   T_FEATURE,
   T_INSTANCES_CLUSTER_K>   &aoipc_inParamClustering,
   int                      argc, 
   char                     **argv
   )
#endif /*__INPARAM_SLINK__*/


#ifdef  __INPARAM_PAM__
template <typename T_CLUSTERIDX,
          typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K
          > 
void 
inparamclustering_getParameter
(InParamPAM
 <T_CLUSTERIDX,
 T_FEATURE,
 T_INSTANCES_CLUSTER_K>        &aoipc_inParamClustering,
  int                          argc, 
  char                         **argv
)
#endif /*__INPARAM_PAM__*/


#ifdef  __INPARAM_ISODATA__
  template <typename T_CLUSTERIDX,
	    typename T_EPSILON,
	    typename T_FEATURE,         
	    typename T_FEATURE_SUM,
	    typename T_INSTANCES_CLUSTER_K
	    > 
void 
inparamclustering_getParameter
  (InParamIsoData
   <T_CLUSTERIDX,
   T_EPSILON,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>         &aoipc_inParamClustering,
   int                            argc, 
   char                           **argv
   )
#endif /*__INPARAM_ISODATA__*/


#ifdef  __INPARAM_FCM__
template<typename T_CLUSTERIDX,
         typename T_U,
         typename T_FEATURE,         
	 typename T_INSTANCES_CLUSTER_K
> 
void 
inparamclustering_getParameter
(InParamFCM
<T_CLUSTERIDX,
 T_U,
 T_FEATURE,
 T_INSTANCES_CLUSTER_K>    &aoipc_inParamClustering,
 int                       argc, 
 char                      **argv
)  
#endif /*__INPARAM_FCM__*/


#ifdef __INPARAM_DBSCAN__
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	  >
  void 
  inparamclustering_getParameter
  (InParamDBSCAN
   <T_FEATURE,
   T_INSTANCES_CLUSTER_K,
   T_CLUSTERIDX>           &aoipc_inParamClustering,
   int                     argc, 
   char                    **argv
   )
#endif /*__INPARAM_DBSCAN__*/


#ifdef  __INPARAM_WITHOUTPCPMFK__
  template<typename T_CLUSTERIDX,
           typename T_BITSIZE,
           typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
          > 
  void 
  inparamclustering_getParameter
  (InParamWithoutPcPmFk
   <T_CLUSTERIDX,T_BITSIZE,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
   &aoipc_inParamClustering,
   int                                     argc, 
   char                                    **argv
   )
#endif /*__INPARAM_WITHOUTPCPMFK__*/

  
#ifdef __INPARAM_PCPMFK__

template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	> 
void 
inparamclustering_getParameter
(InParamPcPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K
 >
 &aoipc_inParamClustering,
 int                                     argc, 
 char                                    **argv
 )  
#endif /*__INPARAM_PCPMFK__*/


#ifdef  __INPARAM_GCA__
template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_getParameter
  (InParamGCA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>      &aoipc_inParamClustering,
   int                         argc, 
   char                        **argv
   )
#endif /*__INPARAM_GCA__*/


#ifdef __INPARAM_GAPROTOTYPESFK__
template<typename T_BITSIZE,
         typename T_CLUSTERIDX,
	 typename T_REAL,
	 typename T_FEATURE,         
	 typename T_FEATURE_SUM,
         typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_getParameter
(InParamGAPrototypesFk
 <T_BITSIZE,
 T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,         
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>    &aoipc_inParamClustering,
 int                       argc, 
 char                      **argv
)
#endif /*__INPARAM_GAPROTOTYPESFK__*/


#ifdef  __INPARAM_HKA__
  template<typename T_CLUSTERIDX,
           typename T_REAL,
           typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
  > 
  void 
  inparamclustering_getParameter
  (InParamHKA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>        &aoipc_inParamClustering,
   int                           argc, 
   char                          **argv
   )
#endif /*__INPARAM_HKA__*/

  
#ifdef  __INPARAM_CBGA__
template<typename T_CLUSTERIDX,
	 typename T_REAL,
	 typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K,
	 typename T_INSTANCE_FREQUENCY> 
void 
inparamclustering_getParameter
  (InParamCBGA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K,
   T_INSTANCE_FREQUENCY>      &aoipc_inParamClustering,
   int                        argc, 
   char                       **argv
   )
#endif /*__INPARAM_CBGA__*/


#ifdef __INPARAM_PMFK__
template<typename T_CLUSTERIDX,
         typename T_REAL,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void
inparamclustering_getParameter
(InParamPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>     &aoipc_inParamClustering,
 int                       argc, 
 char                      **argv
 )
#endif /* __INPARAM_PMFK__ */


#ifdef  __INPARAM_ADAPTIVEPCPMFK__
template<typename T_CLUSTERIDX,
         typename T_FEATURE,         
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
        > 
void 
inparamclustering_getParameter
(InParamAdaptivePcPmFk
 <T_CLUSTERIDX,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aoipc_inParamClustering,
 int                   argc, 
 char                  **argv
 )
#endif /*__INPARAM_ADAPTIVEPCPMFK__*/


#ifdef __INPARAM_TGCA__
  template<typename T_CLUSTERIDX,
	   typename T_REAL,
           typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   > 
  void 
  inparamclustering_getParameter
  (InParamTGCA
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,         
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>                 &aoipc_inParamClustering,
   int                                     argc, 
   char                                    **argv
   )  
#endif /*__INPARAM_TGCA__
	*/

#ifdef  __INPARAM_PCPMVK__
  template<typename T_CLUSTERIDX,
           typename T_REAL,
           typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
          > 
  void 
  inparamclustering_getParameter
  (InParamPcPmVk
   <T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,         
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K
   >                                       &aoipc_inParamClustering,
   int                                     argc, 
   char                                    **argv
   )  
#endif /*__INPARAM_PCPMVK__*/


#ifdef __INPARAM_FEAC__
template<typename T_FEATURE,
	 typename T_REAL,
	 typename T_CLUSTERIDX,
	 typename T_FEATURE_SUM,
	 typename T_INSTANCES_CLUSTER_K
	 > 
void 
inparamclustering_getParameter
(InParamFEAC
 <T_FEATURE,
 T_REAL,
 T_CLUSTERIDX,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                &aoipc_inParamClustering,
 int                                   argc, 
 char                                  **argv
 )
#endif /*__INPARAM_FEAC__*/

#ifdef __INPARAM_GENWOCHGVK__
  template < typename T_BITSIZE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
void 
  inparamclustering_getParameter
  (InParamGenWOChgVk
   <T_BITSIZE,
   T_CLUSTERIDX,
   T_REAL,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>            &aoipc_inParamClustering,
   int                               argc, 
   char                              **argv
   )  
#endif /*__INPARAM_GENWOCHGVK__*/


#ifdef __INPARAM_SUBCLUSTERBINARYVK__

  template<typename T_REAL,
	   typename T_BITSIZE,
           typename T_CLUSTERIDX,
           typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
          > 
  void 
  inparamclustering_getParameter
  (InParamSubClusterBinaryVk
   <T_REAL,
   T_BITSIZE,
   T_CLUSTERIDX,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>       &aoipc_inParamClustering,
   int                          argc, 
   char                         **argv
   )  

#endif /*__INPARAM_SUBCLUSTERBINARYVK__*/


#ifdef __INPARAM_GGA__
  template<typename T_CLUSTERIDX,
           typename T_REAL,
           typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   > 
void 
inparamclustering_getParameter
(InParamGGA
 <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
 &aoipc_inParamClustering,
 int                                     argc, 
 char                                    **argv
 )
#endif /*__INPARAM_GGA__
	*/

{ /*BEGIN inparamclustering_getParameter*/

  const char   *larray_opFormatFile[] = INPARAMCLUSTERING_FORMATINSTANCEFILE;

#ifdef __ALG_CLUSTERING__ /* ONLY CLUSTERING */

  const char   *las_opFuncDistance[] = INPARAMCLUSTERING_DISTANCE_TYPE;
  const char   *las_optGnuplotCoreStyles[] = {"points", "lines", "linespoints", "dot", (char *) NULL };

#endif /* __ALG_CLUSTERING__  */

#ifdef  __INPARAM_STDVAR_MILLIGAN_COOPER1988__
  const char   *las_optStanVar[] =  STD_VAR_NAME;
  const char   *las_optStandardizationVar[] =  {"std-var", (char *) NULL };
#endif /* __INPARAM_STDVAR_MILLIGAN_COOPER1988__ */
 
#ifdef __INPARAM_PCATRANSMATRIX__
  const char   *las_optPCAtransmatrix[] = 
    {"transmatrix-outfile", (char *) NULL };
#endif /*__INPARAM_PCATRANSMATRIX__*/

#ifdef __INPARAM_PLOT_CLUSTERING__
  const char   *las_opProjection[] = INPARAMCLUSTERING_PLOT_PLOJECTION;
  const char   *las_optPlotClustering[] = 
    {"projection",          //0
     "centroids-infile",   //1
     "centroids-title",    //2
     "centroids-id",       //3
     "centroids-point",    //4
     "centroids-size",     //5
     "centroids-infile2",  //6
     "centroids-title2",   //7
     "centroids-id2",      //8
     "centroids-point2",   //9
     "centroids-size2",    //10
     "member-infile",      //11
     "member-label",       //12
     "member-size",        //13
     "id-increment",       //14
     "graph-infile",       //15
     "x-coord",            //16
     "y-coord",            //17
     "z-coord",            //18
     "title",              //19
     "legend",             //20
     "size-instance",      //21
     "label-offset",       //22
     "x-range",            //23
     "y-range",            //24
     "color-graphics",     //25
     "sleep",              //26
     "graphics-outfile",   //27
     (char *) NULL
    };
#endif /*__INPARAM_PLOT_CLUSTERING__*/

#ifdef __INPARAM_KMEANS__
  const char   *las_optKmeans[] = 
    {"number-clusters","iterations", "threshold", (char *) NULL };
#endif /*__INPARAM_KMEANS__*/

#ifdef __INPARAM_SLINK__
  const char   *las_optSLINK[] = 
    {"number-clusters", (char *) NULL };
#endif /*__INPARAM_SLINK__*/

#ifdef __INPARAM_PAM__
  const char   *las_optPAM[] = {"number-clusters","iterations", (char *) NULL };
#endif /*__INPARAM_PAM__*/

#ifdef __INPARAM_ISODATA__ 
  const char   *las_optIsoData[] = 
    {"number-clusters", "iterations", "epsilon", "optimal-initializacion", (char *) NULL };
#endif /*__INPARAM_ISODATA__*/

#ifdef  __INPARAM_DBSCAN__
  const char   *las_optAlgDBSCANEsterKriegelSanderXuClustering1996[] = 
    {"epsilon", 
     "min-pts", 
     "buffer-capacity",
     "tree-fill-factor",
     "tree-index-capacity",
     "tree-leaf-capacity",
     (char *) NULL 
    };
#endif /*__INPARAM_DBSCAN__*/


#ifdef  __INPARAM_FCM__
  const char   *las_optFCM[] = 
    {"number-clusters", "iterations", "epsilon", "weighting-exponent", (char *) NULL };
#endif /*__INPARAM_FCM__*/

#ifdef __INPARAM_PCPMVK__
    const char   *las_optPcPmVk[] = 
    {"k-minimum",
     "k-maximum",
     "population-size", 
     "crossover-probability", 
     "mutation-probability", 
     "generations", 
     (char *) NULL
    };
#endif /*__INPARAM_PCPMVK__*/

  
#ifdef __INPARAM_PCPMFK__

  const char   *lastr_inParamPcPmFk[] = 
    {"number-clusters",
     "population-size", 
     "crossover-probability", 
     "mutation-probability", 
     "generations", 
     (char *) NULL
    };
#endif /*__INPARAM_PCPMFK__*/

#ifdef __INPARAM_GAPROTOTYPESFK__
  const char   *las_optGAPrototypesFk[] = 
    {"number-clusters",
     "population-size", 
     "crossover-probability", 
     "mutation-probability",
     "probability-ini",
     "generations",
     "alpha",
     (char *) NULL
    };
#endif /*__INPARAM_GAPROTOTYPESFK__*/


#ifdef __INPARAM_TGCA__
  const char   *las_optTGCA[] = 
    {"k-minimum",
     "k-maximum",
     "population-size",
     "num-subpopulations-cross",
     "crossover-probability", 
     "generations",
     "kmeans-iterations",
     "kmeans-threshold",
     (char *) NULL
    };
#endif /*__INPARAM_TGCA__*/


#ifdef __INPARAM_GGA__
  const char   *las_optGGA[] = 
    {"k-minimum",
     "k-maximum",
     "sub-population-size",
     "number-island",
     "pe", /*probability of migration*/
     "generations",
     "pci",
     "pcf",
     "pmi",
     "pmf",
     "pbi",
     "pbf",
     (char *) NULL
    };
#endif /*__INPARAM_GGA__
       */  

#ifdef __INPARAM_SUBCLUSTERBINARYVK__
  const char   *las_optCLUSTERINGTsengYang2001[] = 
    {"u-parameter",
     "lambda",
     "w1",
     "w2",
     "population-size", 
     "crossover-probability", 
     "mutation-probability", 
     "generations", 
     (char *) NULL
    };
#endif /*__INPARAM_SUBCLUSTERBINARYVK__*/

#ifdef __INPARAM_GENWOCHGVK__
  const char   *las_optGenWOChgVk[] = 
    {"k-minimum",
     "k-maximum", 
      "population-size", 
      "crossover-probability", 
      "mutation-probability", 
      "generations",
      "notchangestop",
      (char *) NULL
    };
#endif /*__INPARAM_GENWOCHGVK__*/

#ifdef __INPARAM_FEAC__
  const char   *las_optInparamFEAC[] = 
    {"k-minimum",
     "k-maximum",
     "population-size", 
     "generations",
     "desiable-objfunc",
     "kmeans-iterations",
     "kmeans-difference",
     (char *) NULL 
    };
#endif /*__INPARAM_FEAC__
       */

#ifdef  __INPARAM_GCA__
  const char   *las_optGCA[] = 
    {"number-clusters",
     "population-size",               /*Population size                x2*/
     "mix-recombination-probability", /*Mix recombination probability  x7 
					(crossover-probability)*/ 
     "point-mutation-probability",    /*Point mutation probability     x8 
					(mutation-probability)*/
     "generations",                   /*Number of generations          x1*/
     "mix-mutation-probability",      /*Mix mutation probability       x9*/
     (char *) NULL 
    };
#endif /*__INPARAM_GCA__*/

#ifdef  __INPARAM_HKA__
  const char   *las_optHKA[] = 
    {"number-clusters",
     "population-size",               /*Population size                x2*/
     "mix-recombination-probability", /*Mix recombination probability  x7 
					(crossover-probability)*/ 
     "point-mutation-probability",    /*Point mutation probability     x8 
					(mutation-probability)*/
     "generations",                   /*Number of generations          x1*/
     "mix-mutation-probability",      /*Mix mutation probability       x9*/
     "order-tournament",
     "nearest-neighbors",
     "search-heuristic-probability",
     (char *) NULL 
    };
#endif /*__INPARAM_HKA__*/


#ifdef __INPARAM_PMFK__
  const char  *las_optGAPmFk[] = 
    {"number-clusters",
     "population-size", 
     "mutation-probability", 
     "generations", 
     (char *) NULL 
    };
#endif 

#ifdef  __INPARAM_ADAPTIVEPCPMFK__
  const char   *las_optAlgChangEtAl2009[] = 
    {"number-clusters", 
     "population-size", 
     "generations", 
     (char *) NULL 
    };
#endif /*__INPARAM_ADAPTIVEPCPMFK__*/

#ifdef  __INPARAM_CBGA__
  const char   *las_optCBGA[] 
    = {"number-clusters",
       "population-size", 
       "mutation-probability", 
       "generations", 
       "select-method", 
       "gla-iterations", 
       (char *) NULL };
  const char   *las_opSelectMethod[] = INPARAMCLUSTERING_CBGA_SELECMETH; /*?CHECAR FALTA*/
#endif /*__INPARAM_CBGA__*/

  
#ifdef __INPARAM_WITHOUTPCPMFK__
  const char   *las_opWithoutPcPm[] = 
    {"number-clusters", 
     "population-size", 
     "matingpool-size", 
     "generations",  (char *) NULL 
    };
#endif /*__INPARAM_WITHOUTPCPMFK__*/

#ifdef __ALG_CLUSTERING__ /* ONLY CLUSTERING */
  
  const char   *las_opGeneral[] = {"centroids-format", "table-format",  (char *) NULL };

#endif /* __ALG_CLUSTERING__ */
  
  int          li_opt;
  char         *lps_optsubValue; 
  int          li_idxSubOpt;

  char *lptstr_fileInstance        = NULL;
  char *lptstr_fileSelectInstances = NULL;
  char *lptstr_fileInstanceTest    = NULL;

  static struct option long_options[] =
    { {"instances",               required_argument, 0, 'i'},
      {"select-instances",        required_argument, 0, 'x'},
      {"test",                    required_argument, 0, 't'},
      {"format-file",             required_argument, 0, 'b'},
      {"with-header",             required_argument, 0, 'h'},
      {"number-instances",        required_argument, 0, 'u'},
      {"select-attributes",       required_argument, 0, 'a'},
      {"delimit-attributes",      required_argument, 0, 'd'},
      {"class-column",            required_argument, 0, 'c'},
      {"cluster-column",          required_argument, 0, 'e'},
      {"idinstance-column",       required_argument, 0, 'l'},
      {"freq-instance-column",    required_argument, 0, 'f'},

#ifdef __MULTI_INSTANCE
      {"idmultiinstance-column",  required_argument, 0, 'm'},
      {"classmultiinst-column",   required_argument, 0, 'p'},      
#endif /*__MULTI_INSTANCE*/
      
      {"distance",                required_argument, 0, 'n'},
      {"number-runs",             required_argument, 0, 'r'},
      {"random-seed",             required_argument, 0, 'z'},
      {"max-execution-time",      required_argument, 0, 'w'},
      {"centroids-outfile",       required_argument, 0, 'C'},
      {"centroids-format",        required_argument, 0, 0},
      {"membership-outfile",      required_argument, 0, 'M'},
      {"partitionstable-outfile", required_argument, 0, 'T'},
      {"table-format",            required_argument, 0, 0},
      
#ifdef _ALG_GRAPH_BASED_
      {"graph-outfile",           required_argument, 0, 'G'},
#endif /*_ALG_GRAPH_BASED_*/
      
      {"runtime-filename",        required_argument, 0, 'R'},
      {"gnuplot",                 required_argument, 0, 'P'},
      {"verbose",                 no_argument,       0, 'v'},
      {"bar-progress",            no_argument,       0, 'q'},
      {"plot-styles",             required_argument, 0, 'y'},
      {"help",                    no_argument,       0, '?'},

#ifdef  __INPARAM_STDVAR_MILLIGAN_COOPER1988__
      {"std-var",                 required_argument, 0, 0},
#endif /* __INPARAM_STDVAR_MILLIGAN_COOPER1988__ */

#ifdef __INPARAM_PCATRANSMATRIX__
      {"transmatrix-outfile",     required_argument, 0, 0},
#endif /*__INPARAM_PCATRANSMATRIX__*/

#ifdef __INPARAM_PLOT_CLUSTERING__
      {"projection",             required_argument, 0, 0},
      {"transmatrix-infile",     required_argument, 0, 0},
      {"centroids-infile",       required_argument, 0, 0},
      {"centroids-title",        required_argument, 0, 0},
      {"centroids-id",           required_argument, 0, 0},
      {"centroids-point",        required_argument, 0, 0},
      {"centroids-size",         required_argument, 0, 0},
      {"centroids-infile2",      required_argument, 0, 0},
      {"centroids-title2",       required_argument, 0, 0},
      {"centroids-id2",          required_argument, 0, 0},
      {"centroids-point2",       required_argument, 0, 0},
      {"centroids-size2",        required_argument, 0, 0},
      {"member-infile",          required_argument, 0, 0},
      {"member-label",           required_argument, 0, 0},
      {"member-size",            required_argument, 0, 0},
      {"id-increment",           required_argument, 0, 0},
      {"graph-infile",           required_argument, 0, 0},
      {"x-coord",                required_argument, 0, 0},
      {"y-coord",                required_argument, 0, 0},
      {"z-coord",                required_argument, 0, 0},
      {"title",                  required_argument, 0, 0},
      {"legend",                 required_argument, 0, 0},
      {"size-instance",          required_argument, 0, 0},
      {"label-offset",           required_argument, 0, 0},
      {"x-range",                required_argument, 0, 0},
      {"y-range",                required_argument, 0, 0},
      {"color-graphics",         required_argument, 0, 0},
      {"sleep",                  required_argument, 0, 0},
      {"graphics-outfile",       required_argument, 0, 0},
#endif /*__INPARAM_PLOT_CLUSTERING__*/
      
      
#ifdef __INPARAM_KMEANS__
      {"number-clusters",         required_argument, 0, 0},
      {"iterations",              required_argument, 0, 0},
      {"threshold",               required_argument, 0, 0},
#endif /*__INPARAM_KMEANS__ */

#ifdef __INPARAM_SLINK__
      {"number-clusters",         required_argument, 0, 0},
#endif /*__INPARAM_SLINK__*/

      
#ifdef __INPARAM_PAM__
      {"number-clusters",         required_argument, 0, 0},
      {"iterations",              required_argument, 0, 0},
#endif /*__INPARAM_PAM__*/


#ifdef  __INPARAM_ISODATA__
      {"number-clusters",         required_argument, 0, 0},
      {"iterations",              required_argument, 0, 0},
      {"epsilon",                 required_argument, 0, 0},
      {"optimal-initializacion",  no_argument,       0, 0},      
#endif /*__INPARAM_ISODATA__*/

#ifdef  __INPARAM_FCM__
      {"number-clusters",         required_argument, 0, 0},
      {"iterations",              required_argument, 0, 0},
      {"epsilon",                 required_argument, 0, 0},
      {"weighting-exponent",      required_argument, 0, 0},
#endif /*__INPARAM_FCM__*/

#ifdef  __INPARAM_DBSCAN__
      {"epsilon",                 required_argument, 0, 0},
      {"min-pts",                 required_argument, 0, 0},
      {"buffer-capacity",         required_argument, 0, 0},
      {"tree-fill-factor",        required_argument, 0, 0},
      {"tree-index-capacity",     required_argument, 0, 0},
      {"tree-leaf-capacity",      required_argument, 0, 0},
#endif /*__INPARAM_DBSCAN__*/

#ifdef __INPARAM_PCPMFK__
      {"number-clusters",         required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif /*__INPARAM_PCPMFK__*/

#ifdef __INPARAM_GAPROTOTYPESFK__
      {"number-clusters",         required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"probability-ini",         required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
      {"alpha",                   required_argument, 0, 0},
#endif /*__INPARAM_GAPROTOTYPESFK__*/
      
#ifdef __INPARAM_PCPMVK__
      {"k-minimum",               required_argument, 0, 0},
      {"k-maximum",               required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif /*__INPARAM_PCPMVK__*/

#ifdef __INPARAM_TGCA__
      {"k-minimum",               required_argument, 0, 0},
      {"k-maximum",               required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"num-subpopulations-cross",required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
      {"kmeans-iterations",       required_argument, 0, 0},
      {"kmeans-threshold",       required_argument, 0, 0},
#endif /*__INPARAM_TGCA__*/

#ifdef __INPARAM_GGA__
      {"k-minimum",               required_argument, 0, 0}, /*0*/
      {"k-maximum",               required_argument, 0, 0}, /*1*/
      {"sub-population-size",     required_argument, 0, 0}, /*2*/
      {"number-island",           required_argument, 0, 0}, /*3*/
      {"pe",                      required_argument, 0, 0}, /*4. probability of migration*/
      {"generations",             required_argument, 0, 0}, /*5.*/
      {"pci",                     required_argument, 0, 0}, /*6.*/
      {"pcf",                     required_argument, 0, 0}, /*7.*/
      {"pmi",                     required_argument, 0, 0}, /*8.*/
      {"pmf",                     required_argument, 0, 0}, /*9.*/
      {"pbi",                     required_argument, 0, 0}, /*10.*/
      {"pbf",                     required_argument, 0, 0}, /*11.*/
#endif /*__INPARAM_GGA__
       */
      
#ifdef __INPARAM_SUBCLUSTERBINARYVK__
      {"u-parameter",             required_argument, 0, 0},
      {"lambda",                  required_argument, 0, 0},
      {"w1",                      required_argument, 0, 0},
      {"w2",                      required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif /*__INPARAM_SUBCLUSTERBINARYVK__*/

#ifdef __INPARAM_GENWOCHGVK__
      {"k-minimum",               required_argument, 0, 0},
      {"k-maximum",               required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
      {"notchangestop",           required_argument, 0, 0},
#endif /*__INPARAM_GENWOCHGVK__*/

#ifdef __INPARAM_PCPMVK__
      {"k-minimum",               required_argument, 0, 0},
      {"k-maximum",               required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"crossover-probability",   required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif /*__INPARAM_PCPMVK__*/

#ifdef __INPARAM_FEAC__
      {"k-minimum",               required_argument, 0, 0},
      {"k-maximum",               required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
      {"desiable-objfunc",        required_argument, 0, 0},
      {"kmeans-iterations",       required_argument, 0, 0},
      {"kmeans-difference",       required_argument, 0, 0},
#endif /*__INPARAM_FEAC__
       */ 

#ifdef  __INPARAM_GCA__
      {"number-clusters",               required_argument, 0, 0},
      {"population-size",               required_argument, 0, 0},
      {"mix-recombination-probability", required_argument, 0, 0},
      {"point-mutation-probability",    required_argument, 0, 0},
      {"generations",                   required_argument, 0, 0},
      {"mix-mutation-probability",      required_argument, 0, 0},
#endif /*__INPARAM_GCA__*/

#ifdef __INPARAM_HKA__
      {"number-clusters",               required_argument, 0, 0},
      {"population-size",               required_argument, 0, 0},
      {"mix-recombination-probability", required_argument, 0, 0},
      {"point-mutation-probability",    required_argument, 0, 0},
      {"generations",                   required_argument, 0, 0},
      {"mix-mutation-probability",      required_argument, 0, 0},
      {"order-tournament",              required_argument, 0, 0},
      {"nearest-neighbors",             required_argument, 0, 0},
      {"search-heuristic-probability",  required_argument, 0, 0},
#endif /*__INPARAM_HKA__*/

#ifdef  __INPARAM_CBGA__
      {"number-clusters",         required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
      {"select-method",           required_argument, 0, 0},
      {"gla-iterations",          required_argument, 0, 0},
#endif /*__INPARAM_CBGA__*/


#ifdef __INPARAM_PMFK__
      {"number-clusters",         required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"mutation-probability",    required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif

#ifdef  __INPARAM_ADAPTIVEPCPMFK__
      {"number-clusters",         required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif /*__INPARAM_ADAPTIVEPCPMFK__*/

#ifdef __INPARAM_WITHOUTPCPMFK__
      {"number-clusters",         required_argument, 0, 0},
      {"population-size",         required_argument, 0, 0},
      {"matingpool-size",         required_argument, 0, 0},
      {"generations",             required_argument, 0, 0},
#endif /*__INPARAM_WITHOUTPCPMFK__*/

      {0, 0, 0, 0}
    };
  int option_index = 0;

  using namespace std;
  istringstream   liss_stringstream;
  uintidx          luintidx_read;
  std::string     ls_read;
 
  while (1) {


#ifdef __MULTI_INSTANCE

    li_opt = getopt_long
      (argc,argv,
       "i:x:t:b:h:u:a:d:c:e:l:f:m:p:r:R:n:z:w:C:M:T:G:P:y:v:q?",
       long_options, 
       &option_index
       );

#else
 
    li_opt = getopt_long
      (argc,argv,
       "i:x:t:b:h:u:a:d:c:e:l:f:r:R:n:z:w:C:M:T:G:P:y:v:q?",
       long_options, 
       &option_index
       );

#endif /*__MULTI_INSTANCE*/

    if (li_opt == -1)
      break;
    switch (li_opt) {
    case 0:

#ifdef __ALG_CLUSTERING__ /* ONLY CLUSTERING */

      if ( strcmp //centroids-format
	   (long_options[option_index].name,
	    las_opGeneral[0] ) == 0 ) 
	{
	  aoipc_inParamClustering.setPrintCentroidsFormat
	    (aoipc_inParamClustering.isYesNo(optarg, argv[0],long_options[option_index].name));
	}
      else if ( strcmp //table-format
		(long_options[option_index].name,
		 las_opGeneral[1] ) == 0 ) 
	{
	  aoipc_inParamClustering.setPrintTableFormat
	    (aoipc_inParamClustering.isYesNo(optarg, argv[0],long_options[option_index].name));
	}
     
#endif /* __ALG_CLUSTERING__ */
      

#ifdef  __INPARAM_STDVAR_MILLIGAN_COOPER1988__

      if ( strcmp
	   (long_options[option_index].name,
	    las_optStandardizationVar[0] ) == 0 ) 
	{
	  //STANDARDIZATION OF VARIABLES: Z0, Z1, Z2,...
	  if ( ( li_idxSubOpt = 
		 getsubopt_getsubopt
		 (&optarg, 
		  las_optStanVar, 
		  &lps_optsubValue)) != -1
	       ) 
	    {
	      aoipc_inParamClustering.setStandardizationVar(li_idxSubOpt);
	    }
	  else 
	    {
	      aoipc_inParamClustering.errorArgument
		(argv[0],
		 long_options[option_index].name, 
		 las_optStanVar
		 );
	    }
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name, 
	     las_optStandardizationVar
	     );
	}
	  
#endif /* __INPARAM_STDVAR_MILLIGAN_COOPER1988__ */
      
#ifdef __INPARAM_PCATRANSMATRIX__
      
      if ( strcmp
	   (long_options[option_index].name,
	    las_optPCAtransmatrix[0] ) == 0 ) 
	{
	  aoipc_inParamClustering.setOutFileTransMatrix(optarg);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],long_options[option_index].name,las_optPCAtransmatrix);
	}
      
#endif /*__INPARAM_PCATRANSMATRIX__*/

#ifdef __INPARAM_PLOT_CLUSTERING__
      
      if ( strcmp
	   (long_options[option_index].name,
	    las_optPlotClustering[0] ) == 0 ) 
	{
	  if ( (li_idxSubOpt = 
		getsubopt_getsubopt
		(&optarg, 
		 las_opProjection, 
		 &lps_optsubValue)) != -1
	       ) 
	    {
	      aoipc_inParamClustering.setOpProjection(li_idxSubOpt);
	    }
	  else 
	    {
	      aoipc_inParamClustering.errorArgument
		(argv[0],
		 long_options[option_index].name, 
		 las_opProjection
		 );
	    }
	}
      else if ( strcmp //centroids-infile
		(long_options[option_index].name,
		 las_optPlotClustering[1] ) == 0 ) 
	{
	  aoipc_inParamClustering.setInFileCentroids(optarg);
	}

      else if ( strcmp //centroids-title
		(long_options[option_index].name,
		 las_optPlotClustering[2] ) == 0 ) 
	{
	  aoipc_inParamClustering.setTitleCentroids(optarg);
	}
      
      else if ( strcmp //centroids-id
		(long_options[option_index].name,
		 las_optPlotClustering[3] ) == 0 ) 
	{
	  aoipc_inParamClustering.setWithIDCentroids
	    (aoipc_inParamClustering.isYesNo(optarg, argv[0],long_options[option_index].name)); 
	}
      else if ( strcmp //centroids-point
		(long_options[option_index].name,
		 las_optPlotClustering[4] ) == 0 ) 
	{ 
	  int  li_centroidsPoint;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> li_centroidsPoint;
	  aoipc_inParamClustering.setPointTypeCentroids(li_centroidsPoint);
	}
      else if ( strcmp //centroids-size
		(long_options[option_index].name,
		 las_optPlotClustering[5] ) == 0 ) 
	{ 
	  double  ld_centroidsSize;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ld_centroidsSize;
	  aoipc_inParamClustering.setPointSizeCentroids(ld_centroidsSize);
	}
           
      else if ( strcmp //centroids-infile2
		(long_options[option_index].name,
		 las_optPlotClustering[6] ) == 0 ) 
	{
	  aoipc_inParamClustering.setInFileCentroids2(optarg);
	}
      else if ( strcmp //centroids-title2
		(long_options[option_index].name,
		 las_optPlotClustering[7] ) == 0 ) 
	{
	  aoipc_inParamClustering.setTitleCentroids2(optarg);
	}

      else if ( strcmp //centroids-id2
		(long_options[option_index].name,
		 las_optPlotClustering[8] ) == 0 ) 
	{
	  aoipc_inParamClustering.setWithIDCentroids2
	    (aoipc_inParamClustering.isYesNo(optarg, argv[0],long_options[option_index].name)); 
	}
      else if ( strcmp //centroids-point2
		(long_options[option_index].name,
		 las_optPlotClustering[9] ) == 0 ) 
	{ 
	  int  li_centroidsPoint;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> li_centroidsPoint;
	  aoipc_inParamClustering.setPointTypeCentroids2(li_centroidsPoint);
	}
      else if ( strcmp //centroids-size2
		(long_options[option_index].name,
		 las_optPlotClustering[10] ) == 0 ) 
	{ 
	  double  ld_centroidsSize;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ld_centroidsSize;
	  aoipc_inParamClustering.setPointSizeCentroids2(ld_centroidsSize);
	}
      else if ( strcmp //member-infile
		(long_options[option_index].name,
		 las_optPlotClustering[11] ) == 0 ) 
	{
	  aoipc_inParamClustering.setInFileMemberCluster(optarg);
	}
      else if ( strcmp //member-label
		(long_options[option_index].name,
		 las_optPlotClustering[12] ) == 0 ) 
	{
	  aoipc_inParamClustering.setLabelMemberCluster(optarg);
	}
      else if ( strcmp //member-size
		(long_options[option_index].name,
		 las_optPlotClustering[13] ) == 0 ) 
	{ 
	  double  ld_sizeMembeCluster;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ld_sizeMembeCluster;
	  aoipc_inParamClustering.setPointSizeMemberCluster(ld_sizeMembeCluster);
	}

      else if ( strcmp //id-increment
		(long_options[option_index].name,
		 las_optPlotClustering[14] ) == 0 ) 
	{ 
	  int  li_idincrement;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> li_idincrement;
	  aoipc_inParamClustering.setIDincrement(li_idincrement);
	}
      
      else if ( strcmp //graph-infile
		(long_options[option_index].name,
		 las_optPlotClustering[15] ) == 0 ) 
	{
	  aoipc_inParamClustering.setInFileGraph(optarg);
	}
      else if ( strcmp //x-coord
		(long_options[option_index].name,
		 las_optPlotClustering[16] ) == 0 ) 
	{ 
	  uintidx  luintidx_coordX;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_coordX;
	  aoipc_inParamClustering.setCoordX(luintidx_coordX);
	}
      else if ( strcmp //y-coord
		(long_options[option_index].name,
		 las_optPlotClustering[17] ) == 0 ) 
	{ 
	  uintidx  luintidx_coordY;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_coordY;
	  aoipc_inParamClustering.setCoordY(luintidx_coordY);
	}
      else if ( strcmp //z-coord
		(long_options[option_index].name,
		 las_optPlotClustering[18] ) == 0 ) 
	{ 
	  uintidx  luintidx_coordZ;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_coordZ;
	  aoipc_inParamClustering.setCoordZ(luintidx_coordZ);
	}
      else if ( strcmp //title
		(long_options[option_index].name,
		 las_optPlotClustering[19] ) == 0 ) 
	{
	  aoipc_inParamClustering.setTitle(optarg);
	}
      else if ( strcmp //legend
		(long_options[option_index].name,
		 las_optPlotClustering[20] ) == 0 ) 
	{
	  aoipc_inParamClustering.setLegendOption(optarg);
	}

      else if ( strcmp //size-instance
		(long_options[option_index].name,
		 las_optPlotClustering[21] ) == 0 ) 
	{
	  double  ld_sizeInstance;
	  
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ld_sizeInstance;
	  aoipc_inParamClustering.setSizeInstance(ld_sizeInstance);
	}


      else if ( strcmp //label-offset
		(long_options[option_index].name,
		 las_optPlotClustering[22] ) == 0 ) 
	{
	  double  ld_labelOffset;
	  
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ld_labelOffset;
	  aoipc_inParamClustering.setLabelOffset(ld_labelOffset);
	}

      else if ( strcmp //x-range
		(long_options[option_index].name,
		 las_optPlotClustering[23] ) == 0 ) 
	{
	  aoipc_inParamClustering.setXRange(optarg);
	}

      else if ( strcmp //y-range
		(long_options[option_index].name,
		 las_optPlotClustering[24] ) == 0 ) 
	{
	  aoipc_inParamClustering.setYRange(optarg);
	}

      else if ( strcmp //color-graphics
		(long_options[option_index].name,
		 las_optPlotClustering[25] ) == 0 ) 
	{
	  aoipc_inParamClustering.setColorGraphics
	    (aoipc_inParamClustering.isYesNo(optarg, argv[0],long_options[option_index].name)); 
	}
      else if ( strcmp //sleep
		(long_options[option_index].name,
		 las_optPlotClustering[26] ) == 0 ) 
	{
	  unsigned int  luint_sleep;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luint_sleep;
	  aoipc_inParamClustering.setSleep(luint_sleep);
	}
      else if ( strcmp //graphics-outfile
		(long_options[option_index].name,
		 las_optPlotClustering[27] ) == 0 ) 
	{ 
	  aoipc_inParamClustering.setOutFileGraphics(optarg);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],long_options[option_index].name,las_optPlotClustering);
	}
      
#endif /*__INPARAM_PLOT_CLUSTERING__*/
      
      
#ifdef __INPARAM_KMEANS__

      else if ( strcmp
	   (long_options[option_index].name,
	    las_optKmeans[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optKmeans[1] ) == 0 ) 
	{ 
	  COMMON_IDOMAIN lT_readNumMaxIter;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxIter;
	  aoipc_inParamClustering.setNumMaxIter(lT_readNumMaxIter);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optKmeans[2] ) == 0 
		) 
	{ 
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setMinThreshold(luintidx_read);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],long_options[option_index].name,las_optKmeans);
	}

#endif /*__INPARAM_KMEANS__ */

#ifdef __INPARAM_SLINK__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optSLINK[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],long_options[option_index].name,las_optSLINK);
	}
#endif /*__INPARAM_SLINK__*/
      
#ifdef __INPARAM_PAM__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optPAM[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optPAM[1] ) == 0 
		)
	{ 
	  COMMON_IDOMAIN lT_readNumMaxIter;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxIter;
	  aoipc_inParamClustering.setNumMaxIter(lT_readNumMaxIter);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name,
	     las_optPAM
	     );
	}
#endif /*__INPARAM_PAM__*/


#ifdef  __INPARAM_ISODATA__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optIsoData[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].
		 name,las_optIsoData[1]) == 0 
		) 
	{ 
	  aoipc_inParamClustering.setNumMaxIter(atoi(optarg));
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optIsoData[2]) == 0 
		) 
	{
	  T_EPSILON  lT_readEpsilon;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readEpsilon;
	  aoipc_inParamClustering.setEpsilon(lT_readEpsilon);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optIsoData[3]) == 0 
		) 
	{
	  aoipc_inParamClustering.setIsOptimalInitializacion(true);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optIsoData
	   );
      }  
#endif /*__INPARAM_ISODATA__*/

#ifdef  __INPARAM_FCM__

      else if ( strcmp
	   (long_options[option_index].name,
	    las_optFCM[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optFCM[1]) == 0 
		) 
	{ 
	  aoipc_inParamClustering.setNumMaxIter(atoi(optarg));
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optFCM[2]) == 0 
		) 
	{
	  T_U  lT_readEpsilon;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readEpsilon;
	  aoipc_inParamClustering.setEpsilon(lT_readEpsilon);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optFCM[3]) == 0 
		) 
	{
	  T_U lT_readWeightingExponent;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readWeightingExponent;
	  aoipc_inParamClustering.setWeightingExponent
	    (lT_readWeightingExponent);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name, 
	   las_optFCM
	   );
      }

#endif /*ALG_FCM_BEZDEK_H*/


#ifdef  __INPARAM_DBSCAN__

      else if ( strcmp
	   (long_options[option_index].name,
	    las_optAlgDBSCANEsterKriegelSanderXuClustering1996[0] ) == 0 ) 
	{ 
	  double lrt_eps;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_eps;
	  aoipc_inParamClustering.setEps(lrt_eps);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optAlgDBSCANEsterKriegelSanderXuClustering1996[1]) == 0 
		) 
	{ 
	  uint32_t lui32t_minPts;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lui32t_minPts;
	  aoipc_inParamClustering.setMinPts(lui32t_minPts);
	}
      
      else if ( strcmp
		(long_options[option_index].name,
		 las_optAlgDBSCANEsterKriegelSanderXuClustering1996[2]) == 0 
		) 
	{ 
	  uint32_t  lrtreeui32t_bufferCapacity; 
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrtreeui32t_bufferCapacity;
	  aoipc_inParamClustering.setRTreeBufferCapacity(lrtreeui32t_bufferCapacity);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optAlgDBSCANEsterKriegelSanderXuClustering1996[3]) == 0 
		) 
	{
	  double lrtreed_fillFactor;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrtreed_fillFactor;
	  aoipc_inParamClustering.setRTreeFillFactor(lrtreed_fillFactor);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optAlgDBSCANEsterKriegelSanderXuClustering1996[4]) == 0 
		) 
	{
	  uint32_t  lrtreeui32t_indexCapacity;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrtreeui32t_indexCapacity;
	  aoipc_inParamClustering.setRTreeCapacity
	    (lrtreeui32t_indexCapacity);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optAlgDBSCANEsterKriegelSanderXuClustering1996[5]) == 0 
		) 
	{
	  uint32_t lrtreeui32t_leafCapacity; /*tree-leaf-capacity*/
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrtreeui32t_leafCapacity;
	  aoipc_inParamClustering.setRTreeleafCapacity
	    (lrtreeui32t_leafCapacity);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name, 
	   las_optAlgDBSCANEsterKriegelSanderXuClustering1996
	   );
      }

#endif /*__INPARAM_DBSCAN__*/

      
#ifdef __INPARAM_PCPMFK__
      else if ( strcmp /*number-clusters*/
	   (long_options[option_index].name,
	    lastr_inParamPcPmFk[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp /*population-size*/
		(long_options[option_index].name,
		 lastr_inParamPcPmFk[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp /*crossover-probability*/
		(long_options[option_index].name,
		 lastr_inParamPcPmFk[2]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp /*mutation-probability*/
		(long_options[option_index].name,
		 lastr_inParamPcPmFk[3]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}
      else if ( strcmp /*generations*/
		(long_options[option_index].name,
		 lastr_inParamPcPmFk[4]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   lastr_inParamPcPmFk
	   );
      }
#endif /*__INPARAM_PCPMFK__*/

#ifdef __INPARAM_GAPROTOTYPESFK__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optGAPrototypesFk[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPrototypesFk[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPrototypesFk[2]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPrototypesFk[3]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}

      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPrototypesFk[4]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityIni;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityIni;
	  aoipc_inParamClustering.setPini(lT_readProbabilityIni);
	}
      
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPrototypesFk[5]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}

      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPrototypesFk[6]) == 0 
		) 
	{
	  T_REAL  lT_readApha;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readApha;
	  aoipc_inParamClustering.setAlpha(lT_readApha);
	}
      
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optGAPrototypesFk
	   );
      }
#endif /*__INPARAM_GAPROTOTYPESFK__*/


#ifdef __INPARAM_PCPMVK__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optPcPmVk[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMinimum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMinimum;
	  aoipc_inParamClustering.setNumClusterKMinimum(lmcidxT_kMinimum);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optPcPmVk[1] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMaximum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMaximum;
	  aoipc_inParamClustering.setNumClusterKMaximum(lmcidxT_kMaximum);
	}
      else if ( strcmp /*Population-size*/
		(long_options[option_index].name,
		 las_optPcPmVk[2]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp /*crossover-probability*/
		(long_options[option_index].name,
		 las_optPcPmVk[3]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optPcPmVk[4]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optPcPmVk[5]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optPcPmVk
	   );
      }
#endif  /*__INPARAM_PCPMVK__*/

#ifdef __INPARAM_TGCA__
      else if ( strcmp /*k-minimum*/
	   (long_options[option_index].name,
	    las_optTGCA[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMinimum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMinimum;
	  aoipc_inParamClustering.setNumClusterKMinimum(lmcidxT_kMinimum);
	}
      else if ( strcmp /*k-maximum*/
		(long_options[option_index].name,
		 las_optTGCA[1] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMaximum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMaximum;
	  aoipc_inParamClustering.setNumClusterKMaximum(lmcidxT_kMaximum);
	}
      else if ( strcmp /*Population-size*/
		(long_options[option_index].name,
		 las_optTGCA[2]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp /*num-subpopulations-cross*/
		(long_options[option_index].name,
		 las_optTGCA[3]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setNumSubpopulationsCross(luintidx_read);
	}
      else if ( strcmp /*crossover-probability*/
		(long_options[option_index].name,
		 las_optTGCA[4]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp /*generations*/
		(long_options[option_index].name,
		 las_optTGCA[5]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else if ( strcmp /*kmeans-iterations*/
		(long_options[option_index].name,
		 las_optTGCA[6]) == 0 
		) 
	{
	  COMMON_IDOMAIN liT_readKmeansNumMaxIter;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> liT_readKmeansNumMaxIter;
	  aoipc_inParamClustering.setKmeansNumMaxIter
	    (liT_readKmeansNumMaxIter);
	}
      else if ( strcmp /*kmeans-threshold*/
		(long_options[option_index].name,
		 las_optTGCA[7]) == 0 
		) 
	{
	  uintidx lui_kmeansNumMinThreshold;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lui_kmeansNumMinThreshold;
	  aoipc_inParamClustering.setKmeansMinThreshold
	    (lui_kmeansNumMinThreshold);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optTGCA
	   );
      }
#endif /*__INPARAM_TGCA__*/

#ifdef __INPARAM_GGA__
      else if ( strcmp /*k-minimum*/
	   (long_options[option_index].name,
	    las_optGGA[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMinimum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMinimum;
	  aoipc_inParamClustering.setNumClusterKMinimum(lmcidxT_kMinimum);
	}
      else if ( strcmp /*k-maximum*/
		(long_options[option_index].name,
		 las_optGGA[1] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMaximum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMaximum;
	  aoipc_inParamClustering.setNumClusterKMaximum(lmcidxT_kMaximum);
	}
      else if ( strcmp /*sub-population-size*/
		(long_options[option_index].name,
		 las_optGGA[2]) == 0 
		) 
	{
	  uintidx luintidx_readSubPopulationSize;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_readSubPopulationSize;
	  aoipc_inParamClustering.setSubPopulationSize(luintidx_readSubPopulationSize);
	}
      else if ( strcmp /*number-island*/ 
		(long_options[option_index].name,
		 las_optGGA[3]) == 0 
		) 
	{
	  uintidx luintidx_readNumIsland;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_readNumIsland;
	  aoipc_inParamClustering.setNumIsland(luintidx_readNumIsland);
	}
      else if ( strcmp /*probability of migration pe*/
		(long_options[option_index].name,
		 las_optGGA[4]) == 0 
		) 
	{
	  T_REAL  lrt_readPe;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPe;
	  aoipc_inParamClustering.setPe(lrt_readPe);
	}
      else if ( strcmp /*generations*/
		(long_options[option_index].name,
		 las_optGGA[5]) == 0 
		) 
	{
	  COMMON_IDOMAIN lit_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lit_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lit_readNumMaxGenerations);
	}
      else if ( strcmp /*crossover probability initial pci*/
		(long_options[option_index].name,
		 las_optGGA[6]) == 0 
		) 
	{
	  T_REAL  lrt_readPci;
	  
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPci;
	  aoipc_inParamClustering.setPci(lrt_readPci);
	}
      else if ( strcmp /*crossover probability final pcf*/
		(long_options[option_index].name,
		 las_optGGA[7]) == 0 
		) 
	{
	  T_REAL  lrt_readPcf;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPcf;
	  aoipc_inParamClustering.setPcf(lrt_readPcf);
	}
      else if ( strcmp /*mutation probability initial pmi*/
		(long_options[option_index].name,
		 las_optGGA[8]) == 0 
		) 
	{
	  T_REAL  lrt_readPmi;
	  
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPmi;
	  aoipc_inParamClustering.setPmi(lrt_readPmi);
	}
      else if ( strcmp /*mutation probability final pmf*/
		(long_options[option_index].name,
		 las_optGGA[9]) == 0 
		) 
	{
	  T_REAL  lrt_readPmf;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPmf;
	  aoipc_inParamClustering.setPmf(lrt_readPmf);
	}

      else if ( strcmp /*local search probability initial pmi*/
		(long_options[option_index].name,
		 las_optGGA[10]) == 0 
		) 
	{
	  T_REAL  lrt_readPbi;
	  
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPbi;
	  aoipc_inParamClustering.setPbi(lrt_readPbi);
	}
      else if ( strcmp /*local search probability final pmf*/
		(long_options[option_index].name,
		 las_optGGA[11]) == 0 
		) 
	{
	  T_REAL  lrt_readPbf;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lrt_readPbf;
	  aoipc_inParamClustering.setPbf(lrt_readPbf);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optGGA
	   );
      }
#endif /*__INPARAM_GGA__
       */      
      
#ifdef __INPARAM_SUBCLUSTERBINARYVK__
      
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optCLUSTERINGTsengYang2001[0] ) == 0 ) 
	{ 
	  T_REAL  ltr_readUParameter;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ltr_readUParameter;
	  aoipc_inParamClustering.setU(ltr_readUParameter);

	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[1] ) == 0 ) 
	{ 

	  T_REAL  ltr_readLambda;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ltr_readLambda;
	  aoipc_inParamClustering.setLambda(ltr_readLambda);
	
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[2] ) == 0 ) 
	{ 

	  T_REAL  ltr_readW1;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ltr_readW1;
	  aoipc_inParamClustering.setW1(ltr_readW1);
	
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[3] ) == 0 ) 
	{ 

	  T_REAL  ltr_readW2;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> ltr_readW2;
	  aoipc_inParamClustering.setW2(ltr_readW2);
	  
	}
      else if ( strcmp /*population-size*/
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[4]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[5]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[6]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCLUSTERINGTsengYang2001[7]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optCLUSTERINGTsengYang2001
	   );
      }

#endif /*__INPARAM_SUBCLUSTERBINARYVK__*/


#ifdef __INPARAM_GENWOCHGVK__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optGenWOChgVk[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMinimum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMinimum;
	  aoipc_inParamClustering.setNumClusterKMinimum(lmcidxT_kMinimum);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGenWOChgVk[1] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMaximum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMaximum;
	  aoipc_inParamClustering.setNumClusterKMaximum(lmcidxT_kMaximum);
	}
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optGenWOChgVk[2]) == 0 
	   ) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGenWOChgVk[3]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGenWOChgVk[4]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGenWOChgVk[5]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGenWOChgVk[6]) == 0 
		) 
	{
	  COMMON_IDOMAIN lit_numNotChangeStop;

	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lit_numNotChangeStop;
	  aoipc_inParamClustering.setNumNotChangeStop(lit_numNotChangeStop);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optGenWOChgVk
	   );
      }

#endif /*__INPARAM_GENWOCHGVK__*/

      
#ifdef __INPARAM_FEAC__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optInparamFEAC[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMinimum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMinimum;
	  aoipc_inParamClustering.setNumClusterKMinimum(lmcidxT_kMinimum);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optInparamFEAC[1] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_kMaximum;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_kMaximum;
	  aoipc_inParamClustering.setNumClusterKMaximum(lmcidxT_kMaximum);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optInparamFEAC[2]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}

      else if ( strcmp
		(long_options[option_index].name,
		 las_optInparamFEAC[3]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      
      else if ( strcmp
		(long_options[option_index].name,
		 las_optInparamFEAC[4]) == 0 
		) 
	{
	  T_REAL lT_readDesiableObjetiveFunc;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readDesiableObjetiveFunc;
	  aoipc_inParamClustering.setDesiableObjetiveFunc(lT_readDesiableObjetiveFunc);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optInparamFEAC[5]) == 0 
		) 
	{
	  COMMON_IDOMAIN liT_readKmeansNumMaxIter;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> liT_readKmeansNumMaxIter;
	  aoipc_inParamClustering.setKmeansNumMaxIter
	    (liT_readKmeansNumMaxIter);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optInparamFEAC[6]) == 0 
		) 
	{
	  T_FEATURE lT_readKmeansmaxDiffCent;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readKmeansmaxDiffCent;
	  aoipc_inParamClustering.setKmeansMaxDiffCent
	    (lT_readKmeansmaxDiffCent);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_optInparamFEAC
	   );
      }
#endif /*__INPARAM_FEAC__
	*/

#ifdef  __INPARAM_WITHOUTPCPMFK__

      else if ( strcmp
	   (long_options[option_index].name,
	    las_opWithoutPcPm[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_opWithoutPcPm[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_opWithoutPcPm[2]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizeMatingPool(luintidx_read);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_opWithoutPcPm[3]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name,
	   las_opWithoutPcPm
	   );
      }
#endif /*__INPARAM_WITHOUTPCPMFK__*/

#ifdef  __INPARAM_GCA__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optGCA[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp  /*Population size x2*/    
		(long_options[option_index].name,
		 las_optGCA[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp /*Mix recombination probability  x7*/ 
		(long_options[option_index].name,
		 las_optGCA[2]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp /*Point mutation probability     x8*/  
      		(long_options[option_index].name,
		 las_optGCA[3]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}		
      else if ( strcmp /*Number of generations          x1*/
		(long_options[option_index].name,
		 las_optGCA[4]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else if ( strcmp /*Mix mutation probability       x9*/
		(long_options[option_index].name,
		 las_optGCA[5]) == 0 
		) 
	{
	  T_REAL lT_readProbabilityMixMutation;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMixMutation;
	  aoipc_inParamClustering.setProbMixMutation
	    ( lT_readProbabilityMixMutation );
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name,
	     las_optGCA
	     );
	}
#endif /*__INPARAM_GCA__*/


#ifdef __INPARAM_HKA__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optHKA[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp  /*Population size  x2*/             
		(long_options[option_index].name,
		 las_optHKA[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp                    /*Mix recombination probability  x7*/ 
		(long_options[option_index].name,
		 las_optHKA[2]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityCrossover;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityCrossover;
	  aoipc_inParamClustering.setProbCrossover(lT_readProbabilityCrossover);
	}
      else if ( strcmp                   /*Point mutation probability     x8*/  
      		(long_options[option_index].name,
		 las_optHKA[3]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}		
      else if ( strcmp           /*Number of generations          x1*/
		(long_options[option_index].name,
		 las_optHKA[4]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else if ( strcmp           /*Mix mutation probability       x9*/
		(long_options[option_index].name,
		 las_optHKA[5]) == 0 
		) 
	{
	  T_REAL lT_readProbabilityMixMutation;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMixMutation;
	  aoipc_inParamClustering.setProbMixMutation
	    ( lT_readProbabilityMixMutation );
	}
      else if ( strcmp                           /*order-tournament*/
		(long_options[option_index].name,
		 las_optHKA[6]) == 0 ) {
	liss_stringstream.clear();
	liss_stringstream.str(optarg);
	liss_stringstream >> luintidx_read;
	aoipc_inParamClustering.setOrderTournament(luintidx_read);
      }
      else if ( strcmp                           /*order-tournament*/
		(long_options[option_index].name,
		 las_optHKA[7]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setNearestNeighbors(luintidx_read);
	}
      else if ( strcmp           /*Mix mutation probability       x9*/
		(long_options[option_index].name,
		 las_optHKA[8]) == 0 
		) 
	{
	  T_REAL lT_readProbSearchHeuristic;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbSearchHeuristic;
	  aoipc_inParamClustering.setProbSearchHeuristic
	    ( lT_readProbSearchHeuristic );
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name,
	     las_optHKA
	     );
	}
#endif /*__INPARAM_HKA__*/


#ifdef  __INPARAM_CBGA__

      else if ( strcmp
	   (long_options[option_index].name,
	    las_optCBGA[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCBGA[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if 
	( strcmp
	  (long_options[option_index].name,
	   las_optCBGA[2]) == 0 
	  ) 
	{
	  T_REAL  lT_readProbabilityMutation;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}
      else if ( strcmp(long_options[option_index].name,
		       las_optCBGA[3]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCBGA[4]) == 0 
		) 
	{
	  if ( (li_idxSubOpt = 
		getsubopt_getsubopt
		(&optarg, las_opSelectMethod, &lps_optsubValue)) != -1) 
	    {
	      aoipc_inParamClustering.setOpSelectMethod(li_idxSubOpt);
	    }
	  else {
	    aoipc_inParamClustering.errorArgument
	      (argv[0],long_options[option_index].name, las_opSelectMethod);
	  }
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optCBGA[5]) == 0 
		) 
	{
	  int li_numGlaIterations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> li_numGlaIterations;
	  aoipc_inParamClustering.setNumGLAIterations(li_numGlaIterations);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],
	   long_options[option_index].name, 
	   las_optCBGA
	   );
      }
#endif /*__INPARAM_CBGA__*/


#ifdef __INPARAM_PMFK__
      else if ( strcmp
	   (long_options[option_index].name,
	    las_optGAPmFk[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPmFk[1]) == 0 
		) 
	{
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> luintidx_read;
	  aoipc_inParamClustering.setSizePopulation(luintidx_read);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPmFk[2]) == 0 
		) 
	{
	  T_REAL  lT_readProbabilityMutation;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readProbabilityMutation;
	  aoipc_inParamClustering.setProbMutation(lT_readProbabilityMutation);
	}
      else if ( strcmp
		(long_options[option_index].name,
		 las_optGAPmFk[3]) == 0 
		) 
	{
	  COMMON_IDOMAIN lT_readNumMaxGenerations;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lT_readNumMaxGenerations;
	  aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	}
      else {
	aoipc_inParamClustering.errorArgument
	  (argv[0],long_options[option_index].name, 
	   las_optGAPmFk
	   );
      }
#endif /*__INPARAM_PMFK__*/

#ifdef  __INPARAM_ADAPTIVEPCPMFK__

      else if ( strcmp
	   (long_options[option_index].name,
	    las_optAlgChangEtAl2009[0] ) == 0 ) 
	{ 
	  T_CLUSTERIDX lmcidxT_numClusterK;
	  liss_stringstream.clear();
	  liss_stringstream.str(optarg);
	  liss_stringstream >> lmcidxT_numClusterK;
	  aoipc_inParamClustering.setNumClusterK(lmcidxT_numClusterK);
	}
      else
	if ( strcmp
	     (long_options[option_index].name,
	      las_optAlgChangEtAl2009[1]) == 0 
	     ) 
	  {
	    liss_stringstream.clear();
	    liss_stringstream.str(optarg);
	    liss_stringstream >> luintidx_read;
	    aoipc_inParamClustering.setSizePopulation(luintidx_read);
	  }
	else if ( strcmp
		  (long_options[option_index].name,
		   las_optAlgChangEtAl2009[2]) == 0 
		  ) 
	  {
	    COMMON_IDOMAIN lT_readNumMaxGenerations;
	    liss_stringstream.clear();
	    liss_stringstream.str(optarg);
	    liss_stringstream >> lT_readNumMaxGenerations;
	    aoipc_inParamClustering.setNumMaxGenerations(lT_readNumMaxGenerations);
	  }
	else {
	  aoipc_inParamClustering.errorArgument
	    (argv[0],long_options[option_index].name, 
	     las_optAlgChangEtAl2009
	     );
	}

#endif /*__INPARAM_ADAPTIVEPCPMFK__*/

      break;
    case 'i':
      lptstr_fileInstance = optarg;
      break;
    case 'x':
      lptstr_fileSelectInstances = optarg;
      break;
    case 't': 
      lptstr_fileInstanceTest = optarg;
      break;
    case 'b':
      if ( (li_idxSubOpt = 
	    getsubopt_getsubopt
	    (&optarg, 
	     larray_opFormatFile, 
	     &lps_optsubValue)) != -1
	   ) 
	{ 
	  switch (li_idxSubOpt) {
	  case 0:
	    aoipc_inParamClustering.setFormatInstanceFile
	      ( INPARAM_FORMATINSTANCEFILE_UCI );
	    break;
	  case 1:
	    aoipc_inParamClustering.setFormatInstanceFile
	      ( INPARAM_FORMATINSTANCEFILE_KEEL );
	    aoipc_inParamClustering.setSeparateAttributes
	      ( INPARAMCLUSTERING_KEEL_FILE_SEPARATOR );
	    aoipc_inParamClustering.setHaveHeaderFileInstance(false);
	    break;
	  default:
	    break;
	  }
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name, 
	     larray_opFormatFile
	     );
	}
      break;
    case 'h': 
      aoipc_inParamClustering.setHaveHeaderFileInstance
	(aoipc_inParamClustering.isYesNo(optarg, argv[0],long_options[option_index].name));
      break;  
    case 'u':
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setNumInstances(luintidx_read);
      break;
    case 'a':
      aoipc_inParamClustering.setSelectAttributes(optarg);
      break;
    case 'd':
      aoipc_inParamClustering.setSeparateAttributes(optarg);
      break;
    case 'c':
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setClassInstanceColumn(luintidx_read);
      break;
    case 'e':
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setClusterInstanceColumn(luintidx_read);
      break;
    case 'l': 
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setIDInstanceColumn(luintidx_read);
      break;

    case 'f': 
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setInstanceFrequencyColumn(luintidx_read);
      break;

#ifdef __MULTI_INSTANCE

    case 'm':
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setIDMultiInstanceColumn(luintidx_read);
      break;

    case 'p':
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> luintidx_read;
      aoipc_inParamClustering.setClassMultiInstColumn(luintidx_read);
      break;

#endif /*__MULTI_INSTANCE*/
  
    case 'r': 
      aoipc_inParamClustering.setTimesRunAlgorithm(atoi(optarg));
      break;
    case 'R': 
      aoipc_inParamClustering.setFileNameTimesRun(optarg);
      break;

#ifdef __ALG_CLUSTERING__ /* ONLY CLUSTERING */

    case 'n':
      if ( (li_idxSubOpt = 
	    getsubopt_getsubopt
	    (&optarg, 
	     las_opFuncDistance, 
	     &lps_optsubValue)) != -1
	   ) 
	{
	  aoipc_inParamClustering.setOpDistance(li_idxSubOpt);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name, 
	     las_opFuncDistance
	     );
	}
      break;
    case 'z': 
      aoipc_inParamClustering.setRandomSeed(optarg);
      break;
    case 'w': 
      aoipc_inParamClustering.setMaxExecutiontime(atof(optarg));
      break;
    case 'C': 
      aoipc_inParamClustering.setOutFileCentroids(optarg);
      break;
    case 'M': 
      aoipc_inParamClustering.setOutFileMemberShip(optarg);
      break;
    case 'T': 
      aoipc_inParamClustering.setOutFilePartitionsTable(optarg);
      break;
      
#ifdef _ALG_GRAPH_BASED_
    case 'G': 
      aoipc_inParamClustering.setOutFileGraph(optarg);
      break;
#endif /*_ALG_GRAPH_BASED_*/
      
    case 'P':  //PLOT FUNCTION OR OTHER OBJECTS WITH PCA
      aoipc_inParamClustering.setFileNamePlotStatObjetiveFunc(optarg);
      break;
    case 'y':
      if ( (li_idxSubOpt = 
	    getsubopt_getsubopt
	    (&optarg, 
	     las_optGnuplotCoreStyles, 
	     &lps_optsubValue)
	    ) != -1) 
	{
	  aoipc_inParamClustering.setGnuPlotCoreStyles
	    (las_optGnuplotCoreStyles[li_idxSubOpt]);
	}
      else 
	{
	  aoipc_inParamClustering.errorArgument
	    (argv[0],
	     long_options[option_index].name, 
	     las_optGnuplotCoreStyles
	     );
	}
      break;

#endif /* __ALG_CLUSTERING__ */
   
    case 'q': 
      aoipc_inParamClustering.setProgressBarPrinting(true);
      break;

    case 'v': 
#ifdef __VERBOSE_YES
      liss_stringstream.clear();
      liss_stringstream.str(optarg);
      liss_stringstream >> geiinparam_verboseMax;
#endif /*__VERBOSE_YES*/
      break;
    
    case '?': 
      inparamclustering_usage(argv[0], aoipc_inParamClustering);
      break;
      
    default: 
      inparamclustering_usage(argv[0], aoipc_inParamClustering);
      break;
    } /* END switch*/
    } /* END while*/

  if (optind < argc) {
    printf ("non-option ARGV-elements: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
    inparamclustering_usage(argv[0], aoipc_inParamClustering); 
  }
    
  if ( lptstr_fileInstance == NULL ) { 
    std::cout << "missing file:\n "
	      << "  -i, --instances=FILE or DIRECTORY\n"
	      <<"                             file or directory containing data of instances to be clustered\n";
    inparamclustering_usage(argv[0], aoipc_inParamClustering);
  }  
 
  if ( fileutils_isdir(lptstr_fileInstance) ) {

    std::vector<std::string>lvectorstr_fileDir;
    try {
      lvectorstr_fileDir = 
	fileutils_listdir(lptstr_fileInstance);
    }
    catch(const std::invalid_argument&) {
      std::cout << "missing dir:\n "
		<< "  -i, --instances=FILE or DIRECTORY\n"
		<< "                          file or directory containing data of instances to be clustered"
		<< std::endl;
      throw std::invalid_argument("inparamclustering_getparameter: instances directory does not exist");
    }
    

    if ( lptstr_fileSelectInstances != NULL ) {

      std::vector<std::string>&& lvectorstr_fileDirInstance =
	vectorutils::foundStrings
	(lvectorstr_fileDir,
	 std::string(lptstr_fileSelectInstances)
	 );

      aoipc_inParamClustering.setVectorFilesInstance
	( lvectorstr_fileDirInstance );

      if ( lptstr_fileInstanceTest != NULL ) {
	
	std::vector<std::string>&& lvectorstr_fileDirInstanceTest =
	  vectorutils::foundStrings
	  (lvectorstr_fileDir,
	   std::string(lptstr_fileInstanceTest)
	   );
	aoipc_inParamClustering.setVectorFilesInstanceTest
	  ( lvectorstr_fileDirInstanceTest );

	if ( aoipc_inParamClustering.getNumFilesInstance() != 
	     aoipc_inParamClustering.getNumFilesInstanceTest() ) {
	  throw std::invalid_argument
	    ("inparamclustering_getParameter: number files training and test is diferent");
	}
      }
    } else {
      aoipc_inParamClustering.setVectorFilesInstance(lvectorstr_fileDir); 
    }
    } else { /* is dir */
      aoipc_inParamClustering.setFileInstance(lptstr_fileInstance);
      if ( lptstr_fileInstanceTest != NULL ) {
      aoipc_inParamClustering.setFileInstanceTest(lptstr_fileInstanceTest);
    }
    }

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "inout::inparamclustering_getParameter";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

#ifdef __VERBOSE_YES

  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelFilesInstance;
    lostrstream_labelFilesInstance << "<FILEINSTANCES:" << lpc_labelFunc;
    inout::containerprint
      (aoipc_inParamClustering.getVectorFilesInstance().begin(),
       aoipc_inParamClustering.getVectorFilesInstance().end(),
       std::cout,
       lostrstream_labelFilesInstance.str().c_str(),
       ','
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelFilesInstanceTest;
    lostrstream_labelFilesInstanceTest << "<FILEINSTANCESTEST:" << lpc_labelFunc;
    inout::containerprint
      (aoipc_inParamClustering.getVectorFilesInstanceTest().begin(),
       aoipc_inParamClustering.getVectorFilesInstanceTest().end(),
       std::cout,
       lostrstream_labelFilesInstanceTest.str().c_str(),
       ','
       );
    std::cout << std::endl;
 
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END */

} /*END namespace inparam 
   */

#endif  /*__IN_PARAM_CLUSTERING_GETPARAMETER_HPP__*/

