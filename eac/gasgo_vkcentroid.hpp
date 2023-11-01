/*! \file gasgo_vkcentroid.hpp
 *
 * \brief GASGO \cite RoblesBerumen:Zafra:Ventura:GAclusteringVarK:GASGO:2023  
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GASGO algorithm based on the paper:\n

 * Robles-Berumen, H., Zafra, A., Ventura, S. (2023). A Novel 
 * Genetic Algorithm with Specialized Genetic Operators for 
 * Clustering. In: García Bringas, P., et al. Hybrid Artificial 
 * Intelligent Systems. HAIS 2023. Lecture Notes in Computer 
 * Science(), vol 14001. Springer, Cham. 
 * <a href="https://doi.org/10.1007/978-3-031-40725-3_39"</a>.\n
 * \n
 * Library Evolutionary Algorithms for Clustering (LEAC) is a library\n
 * for the implementation of evolutionary algorithms\n
 * focused on the partition type clustering problem. Based on the\n
 * current standards of the <a href="http://en.cppreference.com">C++</a> language, as well as on Standard\n
 * Template Library <a href="http://en.cppreference.com/w/cpp/container">STL</a> 
 * and also  <a href="http://www.openblas.net/">OpenBLAS</a> to have a better performance.\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __GASGO_VKCENTROID_HPP__
#define __GASGO_VKCENTROID_HPP__

#include <leac.hpp>
#include "inparam_pcpmfreqvk.hpp"
#include "outparam_gac.hpp"

#include "plot_runtime_function.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace eac {


/*! \fn uintidx  getMaximumNumKmeansIter(const uintidx aiui_numCluster)
  \brief getMaximumNumKmeansIter get the maximum number of kmeans iterations to optimize the partition
*/
uintidx  getMaximumNumKmeansIter(const uintidx aiui_numCluster)
{
  return (uintidx) ((std::round( 0.5 * std::log( (double) aiui_numCluster )
				 / ( std::sqrt(2.0)* std::log(2.0))) <  1.0 )
		   ? 1
		   : std::round( 0.5 * std::log( (double) aiui_numCluster )
				 / ( std::sqrt(2.0)* std::log(2.0))));  
}

/*! \fn gaencode::ChromCodeBook<T_FEATURE,T_CLUSTERIDX, T_INSTANCE_FREQUENCY, T_INSTANCES_CLUSTER_K, T_FEATURE_SUM, T_REAL> gasgo_vkcentroid(inout::OutParamGAC<T_REAL,
 T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamPcPmFreqVk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K,T_INSTANCE_FREQUENCY> &aiinp_inParamPcPmVk, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist)
  \brief  GASGO  \cite RoblesBerumen:Zafra:Ventura:GAclusteringVarK:GASGO:2023  
  \details Implementation of the GASGO algorithm based on \cite RoblesBerumen:Zafra:Ventura:GAclusteringVarK:GASGO:2023  
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid. 
*/  
template < typename T_FEATURE,
	   typename T_REAL,
	   typename T_INSTANCE_FREQUENCY,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, //INSTANCES FREQUENCY SUM
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromCodeBook
<T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL
 >
gasgo_vkcentroid
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                 &aoop_outParamGAC,
 inout::InParamPcPmFreqVk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K,
 T_INSTANCE_FREQUENCY>         &aiinp_inParamPcPmVk,
 const INPUT_ITERATOR          aiiterator_instfirst,
 const INPUT_ITERATOR          aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist
 )
{  
 
  const uintidx lconstui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  const int lconsti_numGLAIterations = 2;
 
#if defined(__FITNESS_SILHOUETTE__) 
  /*CALCULATE MATRIX DISSIMILARITY
   */
  mat::MatrixTriang<T_REAL>&&
    lmatrixtriagT_dissimilarity = 
    dist::getMatrixDissimilarity
    (aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );
#endif //__FITNESS_SILHOUETTE__
  
  T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];
   
     
  T_FEATURE_SUM *larray_sumFeatureTmp =
    new T_FEATURE_SUM[data::Instance<T_FEATURE>::getNumDimensions()];
  
  stats::sumFeactures
    (larray_sumFeatureTmp,
     aiiterator_instfirst,
     aiiterator_instlast,
     T_FEATURE(0)
     );
  
  stats::meanVector
    (larray_centroid1,
     lconstui_numInstances,
     larray_sumFeatureTmp
     );

#if defined(__FITNESS_INDEX_I__)

  T_REAL lmetric_e1;
  const T_REAL airt_p = 2.0;
  
  lmetric_e1 =
    um::e1
    (larray_centroid1,
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );
    
#endif //__FITNESS_INDEX_I__

  delete [] larray_sumFeatureTmp;
 
  gaencode::ChromCodeBook 
    <T_FEATURE,
     T_CLUSTERIDX,
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,    
     T_FEATURE_SUM,
     T_REAL> 
    lochrom_best
    (3,
     3,
     data::Instance<T_FEATURE>::getNumDimensions(),
     lconstui_numInstances
     );
  lochrom_best.setObjetiveFunc(measuare_undefIndexI(T_REAL));
  
  if ( aiinp_inParamPcPmVk.getNumClusterKMaximum() == 
       INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED)
    aiinp_inParamPcPmVk.setNumClusterKMaximum
      (std::round(std::sqrt((double)lconstui_numInstances)));
  
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   */
  std::vector<
    gaencode::ChromCodeBook
    <T_FEATURE,
     T_CLUSTERIDX, //-1, 0, 1, .., K
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,    
     T_FEATURE_SUM,
     T_REAL>
    >
    lvectorchrom_population;
  
  std::vector<
    gaencode::ChromCodeBook
    <T_FEATURE,
     T_CLUSTERIDX, //-1, 0, 1, .., K
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,    
     T_FEATURE_SUM,
     T_REAL>
    >
    lvectorchrom_matingPool;

  std::uniform_int_distribution<T_CLUSTERIDX>
    uniformdis_mmcidxMinMaxK
    (aiinp_inParamPcPmVk.getNumClusterKMinimum(), 
     aiinp_inParamPcPmVk.getNumClusterKMaximum()
     );

  const T_REAL  aifeact_ai = (T_REAL)  aiinp_inParamPcPmVk.getNumClusterKMinimum();
  const T_REAL  aifeact_bi = (T_REAL)  aiinp_inParamPcPmVk.getNumClusterKMaximum();
   
  std::uniform_int_distribution<uintidx> uniformdis_idxInstances
    (0,lconstui_numInstances-1);
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gasgo_vkcentroid:";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA
      << "  IN(" << geiinparam_verbose << ")\n"
      << "\t(output inout::OutParamGAC&: aoop_outParamGAC[" 
      << &aoop_outParamGAC << "]\n"
      << "\t input  InParamClusteringGaProbFk&: aiinp_inParamPcPmVk[" 
      << &aiinp_inParamPcPmVk << "]\n"
      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      << &aifunc2p_dist << ']'
      << "\n\t\tPopulation size = " 
      << aiinp_inParamPcPmVk.getSizePopulation()
      << "\n\t\tProbCrossover = " 
      << aiinp_inParamPcPmVk.getProbCrossover() 
      << "\n\t\tProbMutation  = " 
      << aiinp_inParamPcPmVk.getProbMutation()
      << "\n\t\tk-minimum  = " 
      << aiinp_inParamPcPmVk.getNumClusterKMinimum()
      << "\n\t\tk-maximum  = " 
      << aiinp_inParamPcPmVk.getNumClusterKMaximum()
      << "\n\t)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 

  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamPcPmVk.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream  lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_DBindex = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL> lvectort_statfuncObjetiveFunc;
  
  if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
    
    lvectort_statfuncObjetiveFunc.reserve
      ( aiinp_inParamPcPmVk.getSizePopulation());
    //DEFINE FUNCTION
    lofh_DBindex  = new runtime::RuntimeFunctionValue<T_REAL>
      ("DB", 
       aiinp_inParamPcPmVk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_DBindex);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamPcPmVk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamPcPmVk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamPcPmVk.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open
      (aoop_outParamGAC.getFileNameOutPlotStatObjetiveFunc().c_str(),
       std::ios::out | std::ios::app
       );

    lfileout_plotStatObjetiveFunc.precision(COMMON_COUT_PRECISION);

    //FUNCTION HEADER

    lfileout_plotStatObjetiveFunc 
      <<  llfh_listFuntionHist.getHeaderFuntions() 
      << "\n";
  }

#endif /*__WITHOUT_PLOT_STAT*/
  
 
  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING
   */
  aoop_outParamGAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION
   */
  runtime::ExecutionTime let_executionTime = runtime::start();

  /*Population initialization:
    In the first generation, a random population representing different 
    search space solutions is created at the initial stage. This random 
    sampling is carried out in two steps:
        -Initialization of number of groups
        -Initialization of centroids
  */
  {/*BEGIN CREATE SPACE FOR STORE POPULATION*/
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "CREATE SPACE FOR STORE POPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep  
	<< ": IN(" << geiinparam_verbose << ')'
	<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    /*CREATE SPACE FOR STORE POPULATION-----------------------------------------
     */
    lvectorchrom_population.reserve
      (aiinp_inParamPcPmVk.getSizePopulation());
    lvectorchrom_matingPool.reserve
      (aiinp_inParamPcPmVk.getSizePopulation());
  
    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinp_inParamPcPmVk.getSizePopulation(); 
	 luintidx_i++) 
      {
	
	/*It is generated random values of k (lui_numClusterFk) to establish 
	  the number of initial groups of each chromosome, this number will 
	  be between \in [2,\sqrt{n}} being n the number of instances in 
	  the dataset
	 */
	uintidx lui_numClusterFk  = uintidx(uniformdis_mmcidxMinMaxK(gmt19937_eng));
	
	lvectorchrom_population.push_back
	  (gaencode::ChromCodeBook
	   <T_FEATURE,
	   T_CLUSTERIDX, //-1, 0, 1, .., K
	   T_INSTANCE_FREQUENCY,
	   T_INSTANCES_CLUSTER_K,    
	   T_FEATURE_SUM,
	   T_REAL> 
	   (lui_numClusterFk,
	    lui_numClusterFk,
	    data::Instance<T_FEATURE>::getNumDimensions(),
	    lconstui_numInstances
	    )
	   );	
      }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep
	<< ": OUT(" << geiinparam_verbose << ')'
	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  } /*END CREATE SPACE FOR STORE POPULATION*/


      
  { /*BEGIN POPULATION INITIALIZATION
     */
      
#ifdef __VERBOSE_YES
    const char *geverbosepc_labelstep = "POPULATION INITIALIZATION:";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep
	<< " IN(" << geiinparam_verbose << ')'
	<< std::endl;
    }
#endif /*__VERBOSE_YES*/
   
    long ll_invalidOffspring = 0;
    for ( auto &&liter_chromCodeBook: lvectorchrom_population ) {      
   
#if defined(__FITNESS_DB_INDEX__)     
      T_REAL lor_objetiveValue = 1.0 / measuare_undefDBindex(T_REAL);
#endif //__FITNESS_DB_INDEX__
      
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)     
      T_REAL lor_objetiveValue = measuare_undefSilhouette(T_REAL);
#endif //__FITNESS_SIMPLIFIED_SILHOUETTE__
      
#if defined(__FITNESS_SILHOUETTE__)
      T_REAL lor_objetiveValue = measuare_undefSilhouette(T_REAL);
#endif //__FITNESS_SILHOUETTE__

#if defined(__FITNESS_INDEX_I__)      
      T_REAL lor_objetiveValue = measuare_undefIndexI(T_REAL);
#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
      T_REAL lor_objetiveValue = 1.0/measuare_undefWBIndex(T_REAL);
#endif //__FITNESS_WB_INDEX__
      
#if defined(__FITNESS_VRC__)
      T_REAL lor_objetiveValue = measuare_undefVRC(T_REAL);
#endif //__FITNESS_VRC__)
      
      liter_chromCodeBook.setObjetiveFunc(lor_objetiveValue);

      /*Initialization of centroids. It is selected random instances to be 
	the initial centroids of each solution.
       */
      clusteringop::randomInitialize
	(liter_chromCodeBook.getCodeBook(), 
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );

      liter_chromCodeBook.setOptimalityCBGA(gaencode::OPT_NONE);
	   
      if ( lconsti_numGLAIterations > 0 )  {
	
	gaclusteringop::iterateGLAAux
	  (liter_chromCodeBook,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   lconsti_numGLAIterations,
	   aifunc2p_dist
	   ); 
      }
	   
      /*Generate optimal Partitioning
       */
      std::pair<bool,T_REAL> lpair_validPartCVIDistortion =
	clusteringop::reassignCluster
	(liter_chromCodeBook.getPartition(),
	 liter_chromCodeBook.getCodeBook(),
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 aifunc2p_dist
	 );
      
      //CONVERT Distortion to SSE
      lpair_validPartCVIDistortion.second *=
	((T_REAL) lconstui_numInstances *  
	 (T_REAL) data::Instance<T_FEATURE>::getNumDimensions());
 
      if ( lpair_validPartCVIDistortion.first &&
	   liter_chromCodeBook.getCodeBook().getNumRows() >= 2 ) {

	mat::MatrixRow<T_FEATURE> aimatrixt_centroids
	  (liter_chromCodeBook.getCodeBook().getNumRows(),
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   liter_chromCodeBook.getCodeBook().toArray()
	   );

	
#if defined(__FITNESS_DB_INDEX__)
	
	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartitionlabel_clusters
	  (liter_chromCodeBook.getPartition().getMembersShip(),
	   lconstui_numInstances,
	   (T_CLUSTERIDX) liter_chromCodeBook.getCodeBook().getNumRows() 
	   );
	
	lor_objetiveValue = 
	  1.0 / um::dbindex
	  (aimatrixt_centroids,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   lpartitionlabel_clusters,
	   aifunc2p_dist
	   );
      
#endif //__FITNESS_DB_INDEX__
      
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)

	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartitionlabel_clusters
	  (liter_chromCodeBook.getPartition().getMembersShip(),
	   lconstui_numInstances,
	   (T_CLUSTERIDX) liter_chromCodeBook.getCodeBook().getNumRows() 
	   );
	
	std::vector<T_REAL>&&  lvectort_partialSilhouette =
	  um::simplifiedSilhouette
	  (aimatrixt_centroids,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   lpartitionlabel_clusters,
	   liter_chromCodeBook.getPartition().getNumInstancesClusterK(),
	   aifunc2p_dist
	   );

	T_REAL lmetrict_sumPartialSilhouette = 
	  interfacesse::sum
	  (lvectort_partialSilhouette.data(),
	   (uintidx) lvectort_partialSilhouette.size()
	   );

	lor_objetiveValue = 
	  (lvectort_partialSilhouette.size() != 0)?
	  lmetrict_sumPartialSilhouette / (T_REAL) lvectort_partialSilhouette.size():
	  measuare_undefSilhouette(T_REAL);
	
#endif // __FITNESS_SIMPLIFIED_SILHOUETTE__
	

#if defined(__FITNESS_SILHOUETTE__)

	lor_objetiveValue = 
	  um::silhouette
	  (lmatrixtriagT_dissimilarity,
	   liter_chromCodeBook.getPartition()
	   );
	
#endif //__FITNESS_SILHOUETTE__

  
#if defined(__FITNESS_INDEX_I__)
	
	T_REAL lmetric_dk =
	  um::maxDistCjCjp
	  (aimatrixt_centroids, 	    
	   aifunc2p_dist
	   );
	lor_objetiveValue = ( lpair_validPartCVIDistortion.second > 0.0 )?
	  (( lmetric_e1 / lpair_validPartCVIDistortion.second ) *  lmetric_dk )
	  / T_REAL(aimatrixt_centroids.getNumRows())
	  : measuare_undefIndexI(T_REAL);
  
	lor_objetiveValue = std::pow(lor_objetiveValue,airt_p);
	
#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
	
	T_REAL lmetrict_SSb =
	  um::ssb
	  (aimatrixt_centroids,
	   larray_centroid1,
	   liter_chromCodeBook.getPartition().getNumInstancesClusterK(),
	   aifunc2p_dist
	   );

	lor_objetiveValue =
	  ( aimatrixt_centroids.getNumRows() > 0.0
	    && lpair_validPartCVIDistortion.second > 0.0 )?
	  lmetrict_SSb /( T_REAL(aimatrixt_centroids.getNumRows())
			  * lpair_validPartCVIDistortion.second)
	  :1.0/measuare_undefWBIndex(T_REAL);

#endif //__FITNESS_WB_INDEX__
      
#if defined(__FITNESS_VRC__)

	T_REAL lmetrict_SSb =
	  um::ssb
	  (aimatrixt_centroids,
	   larray_centroid1,
	   liter_chromCodeBook.getPartition().getNumInstancesClusterK(),
	   aifunc2p_dist
	   );

	lor_objetiveValue =
	  ( lpair_validPartCVIDistortion.second > 0.0 )?
	  (lmetrict_SSb  * (T_REAL(lconstui_numInstances) -
			    T_REAL(aimatrixt_centroids.getNumRows())))
	  / (lpair_validPartCVIDistortion.second *
	     (T_REAL(aimatrixt_centroids.getNumRows()) -1) )
	  : measuare_undefIndexI(T_REAL);
	
#endif //__FITNESS_VRC__
	
      } // END if ( liter_chromCodeBook.getCodeBook().getNumRows() >= 2 ) {
      
      liter_chromCodeBook.setValidString(lpair_validPartCVIDistortion.first);
      liter_chromCodeBook.setObjetiveFunc(lor_objetiveValue);
      
#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << "NUM CLUSTER: " << liter_chromCodeBook.getCodeBook().getNumRows()
	  << " OBJECTIVE: " << liter_chromCodeBook.getObjetiveFunc() << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
      if ( liter_chromCodeBook.getValidString() == false ) 
	++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
      lvectort_statfuncObjetiveFunc.push_back(liter_chromCodeBook.getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	liter_chromCodeBook.print();
	std::cout << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    }

    /*METRIC INVALID SOLUTION
     */
    aoop_outParamGAC.sumTotalInvalidOffspring
      (ll_invalidOffspring);
    

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< " OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  } /*END POPULATION INITIALIZATION*/


   { /*BEGIN PRESERVING THE BEST STRING
       */
#ifdef __VERBOSE_YES
      const char *geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/


      auto lit_chromMax = std::max_element
	(lvectorchrom_population.begin(), 
	 lvectorchrom_population.end(), 
	 [](gaencode::ChromCodeBook
	    <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_INSTANCE_FREQUENCY,
	    T_INSTANCES_CLUSTER_K,    
	    T_FEATURE_SUM,
	    T_REAL> & x, 
            gaencode::ChromCodeBook
	    <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_INSTANCE_FREQUENCY,
	    T_INSTANCES_CLUSTER_K,    
	    T_FEATURE_SUM,
	    T_REAL>& y
	    ) 
	 {

	   return x.getObjetiveFunc() < y.getObjetiveFunc();
	   
	 }
	 );
        
      if ( lochrom_best.getObjetiveFunc() < (*lit_chromMax).getObjetiveFunc() ) {
	lochrom_best = *lit_chromMax;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION 
	 */
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  lochrom_best.print();
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      } /*END IF*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END PRESERVING THE BEST STRING*/
  
    /*MEASUREMENT BEST: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT
    if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
      lofh_DBindex->setValue(lochrom_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectort_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectort_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/
    

  while ( 1 ) {

    /*Fitness computation:
      Fusion fitness is evaluated when an individual is obtained 
      through crossing and mutation.
    */

    {/*BEGIN FITNESS COMPUTATION*/
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "FITNESS COMPUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ": IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << " OUT(" << geiinparam_verbose << ')'
	  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    } //END FITNESS COMPUTATION
   
    
    /*TERMINATION CRITERION 
     */
    { /*BEGIN TERMINATION CRITERION*/
#ifdef __VERBOSE_YES

      /*ID PROC
       */
      ++geverboseui_idproc;
      const char *geverbosepc_labelstep = "TERMINATION CRITERION ATTAINED:";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << "generation "
	  << llfh_listFuntionHist.getDomainUpperBound()
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
      if ( !(llfh_listFuntionHist.getDomainUpperBound() 
	     < aiinp_inParamPcPmVk.getNumMaxGenerations() ) 
	   )
	break;
      
      llfh_listFuntionHist.increaseDomainUpperBound();
    } /*END TERMINATION CRITERION
       */

    
    /*GENETIC OPERATIONS
     */
    { /*BEGIN GENETIC OPERATIONS*/

      /*Selection: Conventional proportional selection is applied 
	to the population of individuals for the roulette process.
      */
      { /*BEGIN SELECTION*/

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "SELECTION";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << geverbosepc_labelstep  
	    << ":  IN(" << geiinparam_verbose << ')'
	    << std::endl;
	}
#endif /*__VERBOSE_YES*/
    
	const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	  prob::makeDistRouletteWheel
	  (lvectorchrom_population.begin(),
	   lvectorchrom_population.end(),
	   [](gaencode::ChromCodeBook
	      <T_FEATURE,
	      T_CLUSTERIDX, //-1, 0, 1, .., K
	      T_INSTANCE_FREQUENCY,
	      T_INSTANCES_CLUSTER_K,    
	      T_FEATURE_SUM,
	      T_REAL>& 
	      liter_chromCodeBook
	      )
	   -> T_REAL
	   {
	     
#if defined(__FITNESS_SILHOUETTE__)
	     return liter_chromCodeBook.getObjetiveFunc() + 2.0 ;
#elif defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)
	     return liter_chromCodeBook.getObjetiveFunc() + 2.0 ;
#elif defined(__FITNESS_DB_INDEX__)
	     return liter_chromCodeBook.getObjetiveFunc();
#else //__FITNESS_OTHER__
	     return liter_chromCodeBook.getObjetiveFunc(); 
#endif	     
	   }	   
	   );
	
	/*COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL--------------------------
	 */
	lvectorchrom_matingPool.clear();
	
	for (uintidx luintidx_i = 0; 
	     luintidx_i < aiinp_inParamPcPmVk.getSizePopulation(); 
	     luintidx_i++) 
	  {      
	    uintidx luintidx_chrom = 
	      gaselect::getIdxRouletteWheel
	      (lvectorT_probDistRouletteWheel,
	       uintidx(0)
	       );
	    
	    lvectorchrom_matingPool.push_back
	      (gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>
	       (lvectorchrom_population.at(luintidx_chrom))
	       );
	  }

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << " OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

      } /*END SELECTION*/


      /*Crossover
       */
      { /*BEGIN CROSSOVER*/

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "CROSSOVER";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep  
		    << ": IN(" << geiinparam_verbose << ')'
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	long ll_invalidOffspring = 0;

	T_REAL lrt_sumFitness = 
	  std::accumulate
	  (lvectorchrom_matingPool.begin(),
	   lvectorchrom_matingPool.end(),
	   T_REAL(0.0), //VALUE INITIAL
	   [&](T_REAL lT_sumPartial,
	       gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>& lchromfixleng_iter)
	   {
	     return lT_sumPartial + lchromfixleng_iter.getObjetiveFunc(); 
	   }
	   ); 
	T_REAL lrt_aveFitness =
	  lrt_sumFitness / T_REAL(lvectorchrom_matingPool.size());

	auto lit_chromMax = std::max_element
	  (lvectorchrom_matingPool.begin(), 
	   lvectorchrom_matingPool.end(), 
	   [](gaencode::ChromCodeBook
	      <T_FEATURE,
	      T_CLUSTERIDX, //-1, 0, 1, .., K
	      T_INSTANCE_FREQUENCY,
	      T_INSTANCES_CLUSTER_K,    
	      T_FEATURE_SUM,
	      T_REAL> & x, 
	      gaencode::ChromCodeBook
	      <T_FEATURE,
	      T_CLUSTERIDX, //-1, 0, 1, .., K
	      T_INSTANCE_FREQUENCY,
	      T_INSTANCES_CLUSTER_K,    
	      T_FEATURE_SUM,
	      T_REAL>& y
	      ) 
	   {
	     return x.getObjetiveFunc() < y.getObjetiveFunc(); 
	   }
	   );
       
	T_REAL lrt_maxFitness = (*lit_chromMax).getObjetiveFunc();
	
       	gaiterator::crossover
	  (lvectorchrom_matingPool.begin(),
	   lvectorchrom_matingPool.end(),
	   lvectorchrom_population.begin(),
	   lvectorchrom_population.end(),
	   [&](gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>&
	       aichrom_parent1,
	       gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>&
	       aichrom_parent2,
	       gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>&
	       aochrom_child1, 
	       gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>&
	       aochrom_child2
	       )
	   {

	    /* Probability of crossover is  adaptive
	     */

	     T_REAL lT_fprime =
	       (aichrom_parent1.getObjetiveFunc()
		> aichrom_parent2.getObjetiveFunc() )?
	       aichrom_parent1.getObjetiveFunc()
	       :aichrom_parent2.getObjetiveFunc();
	 
	     T_REAL lT_probCrossover = 
	       prob::adaptiveProb
	       ((T_REAL) 1.0,  /*k1*/
		(T_REAL) 1.0,  /*k3*/
		lrt_maxFitness,
		lrt_aveFitness, 
		lT_fprime
		);

	     /* 
		Another option evaluated was with a probability 
	       of crossing with a fixed
	       if ( uniformdis_real01(gmt19937_eng) <
	    	  aiinp_inParamPcPmVk.getProbCrossover()  ) {

	       It checks if the crossover operator 
		is carried out
	      */
	  if ( uniformdis_real01(gmt19937_eng) < lT_probCrossover ) {

	    /* The BLX-alpha operator is used for the value of k''_v 
	       of the daughter chromosomes, which is uniformly randomly 
	       selected in the interval:

	                   [Cmin -I * alpha, Cmax + I * alpha  
	    */

	    const float lf_alpha =  0.5f;
	    
	       T_CLUSTERIDX lcidx_minNumCluster = (T_CLUSTERIDX)
		 std::min
		 (aichrom_parent1.getCodeBook().getNumRows(),
		  aichrom_parent2.getCodeBook().getNumRows()
		  );
	       
	       T_CLUSTERIDX lcidx_maxNumCluster = (T_CLUSTERIDX)
		 std::max
		 (aichrom_parent1.getCodeBook().getNumRows(),
		  aichrom_parent2.getCodeBook().getNumRows()
		  );

	       T_CLUSTERIDX lcidx_intNumCluster = (T_CLUSTERIDX)
		 std::round( (float)  lf_alpha *
			     (lcidx_maxNumCluster - lcidx_minNumCluster));

	       T_CLUSTERIDX lcidx_minNumClusterInterval =
		 (lcidx_minNumCluster-lcidx_intNumCluster) <
		 aiinp_inParamPcPmVk.getNumClusterKMinimum()?
		 aiinp_inParamPcPmVk.getNumClusterKMinimum():
		 lcidx_minNumCluster-lcidx_intNumCluster;
	       
	       T_CLUSTERIDX lcidx_maxNumClusterInterval =
		 (lcidx_maxNumCluster+lcidx_intNumCluster) >
		 aiinp_inParamPcPmVk.getNumClusterKMaximum()?
		 aiinp_inParamPcPmVk.getNumClusterKMaximum():
		 lcidx_maxNumCluster+lcidx_intNumCluster;
	       
	       if ((lcidx_minNumCluster + lcidx_maxNumCluster) <
		   lcidx_maxNumClusterInterval )
		 lcidx_maxNumClusterInterval =
		   lcidx_minNumCluster + lcidx_maxNumCluster;
	       
	       std::uniform_int_distribution<T_CLUSTERIDX>
	       uniformdis_mmcidxChild
	       (lcidx_minNumClusterInterval, 
		lcidx_maxNumClusterInterval
		);
										
	       T_CLUSTERIDX lcidx_numClusterChild1 =
		 uniformdis_mmcidxChild(gmt19937_eng);
	       T_CLUSTERIDX lcidx_numClusterChild2 =
		 uniformdis_mmcidxChild(gmt19937_eng);

#ifdef __VERBOSE_YES
	       ++geiinparam_verbose;
	       if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		 std::cout
		   << " PARENT NUM CLUSTER: "
		   << aichrom_parent1.getCodeBook().getNumRows()
		   << " - "  << aichrom_parent2.getCodeBook().getNumRows()
		   << " MIN MAX " << lcidx_minNumCluster << " : "
		   << lcidx_maxNumCluster
		   << " INT  " << lcidx_minNumClusterInterval  << " : "
		   << lcidx_maxNumClusterInterval
		   << " CHILD1 " << lcidx_numClusterChild1
		   << " CHILD2 " << lcidx_numClusterChild2
		   << std::endl;
		 std::cout
		   << "CROSS PAREN1  ITER: " << std::endl;
		 aichrom_parent1.print();
		 std::cout
		   << std::endl;
		 std::cout
		   << "CROSS PAREN2  ITER: " << std::endl;
		 aichrom_parent2.print();
		 std::cout
		   << std::endl;
	       }
	       --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	       
	       
	       gaclusteringop::crossPNNnew
		 (aochrom_child1,
		  aichrom_parent1,
		  aichrom_parent2,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lcidx_numClusterChild1,
		  aifunc2p_dist
		  );

	       gaclusteringop::crossPNNnew
		 (aochrom_child2,
		  aichrom_parent1,
		  aichrom_parent2,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lcidx_numClusterChild2,
		  aifunc2p_dist
		  );


	       /*Check if the first child chromosome is valid.
		*/
	       std::pair<T_REAL,bool>  lpair_distortion =
		 um::distortion
		 (aochrom_child1.getCodeBook(),
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  aochrom_child1.getPartition().getMembersShip(),
		  aifunc2p_dist,
		  [](const data::Instance<T_FEATURE>* aiinst_iter ) -> T_REAL
		  {
		    data::InstanceFreq
		      <T_FEATURE,
		       T_INSTANCE_FREQUENCY
		       >
		      *lptinstfreq_iter =
		      (data::InstanceFreq
		       <T_FEATURE,
		       T_INSTANCE_FREQUENCY
		       >*)
		      aiinst_iter;

		    return  (T_REAL) lptinstfreq_iter->getFrequency();
		  }
		  );

#ifdef __VERBOSE_YES
	       ++geiinparam_verbose;
	       if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		 std::cout
		   << "CHILD1 DISTORTION: "
		   << lpair_distortion.first << " : "
		   <<  lpair_distortion.second << std::endl;
		 std::cout
		   << "CROSS CHILD 1 ITER: "
		   <<  aochrom_child1.getCodeBook().getNumRows()
		   << std::endl;
		 aochrom_child1.print();
		 std::cout << std::endl;
	       }
	       --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	       	     
	     if ( lpair_distortion.second == false ) {
	       
	       std::pair<bool,T_REAL> lpair_generate  =   
		  clusteringop::reassignCluster
		  (aochrom_child1.getPartition(),
		   aochrom_child1.getCodeBook(),
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   aifunc2p_dist
		   );
	       
	       lpair_distortion.first = lpair_generate.second;
	       lpair_distortion.second = lpair_generate.first;
	       
	      }
	     
#if defined(__FITNESS_DB_INDEX__)     
	     T_REAL lor_objetiveValue = 1.0 / measuare_undefDBindex(T_REAL); 
#endif //__FITNESS_DB_INDEX__
      
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)     
	     T_REAL lor_objetiveValue = measuare_undefSilhouette(T_REAL); 
#endif //__FITNESS_SIMPLIFIED_SILHOUETTE__
      
#if defined(__FITNESS_SILHOUETTE__)
	     T_REAL lor_objetiveValue = measuare_undefSilhouette(T_REAL); 
#endif //__FITNESS_SILHOUETTE__

#if defined(__FITNESS_INDEX_I__)      
	     T_REAL lor_objetiveValue = measuare_undefIndexI(T_REAL); 
#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
	     T_REAL lor_objetiveValue = 1.0/measuare_undefWBIndex(T_REAL); 
#endif //__FITNESS_WB_INDEX__
	     
#if defined(__FITNESS_VRC__)
	     T_REAL lor_objetiveValue = measuare_undefVRC(T_REAL); 
#endif //__FITNESS_VRC__)

	     
	     if ( lpair_distortion.second &&
		  aochrom_child1.getCodeBook().getNumRows() >= 2 ) {

	       lpair_distortion.first *=
		 ((T_REAL) lconstui_numInstances *  
		 (T_REAL) data::Instance<T_FEATURE>::getNumDimensions());
		
	       mat::MatrixRow<T_FEATURE> aimatrixt_centroids
		 (aochrom_child1.getCodeBook().getNumRows(),
		  data::Instance<T_FEATURE>::getNumDimensions(),
		  aochrom_child1.getCodeBook().toArray()
		  );


#if defined(__FITNESS_DB_INDEX__)

	       partition::PartitionLabel
		 <T_CLUSTERIDX>
		 lpartitionlabel_clusters
		 (aochrom_child1.getPartition().getMembersShip(),
		  lconstui_numInstances,
		  (T_CLUSTERIDX) aochrom_child1.getCodeBook().getNumRows() 
		  );
	
	       lor_objetiveValue =
		 1.0 / um::dbindex
		 (aimatrixt_centroids,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lpartitionlabel_clusters,
		  aifunc2p_dist
		  );

#endif //__FITNESS_DB_INDEX__
       
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)

	       partition::PartitionLabel
		 <T_CLUSTERIDX>
		 lpartitionlabel_clusters
		 (aochrom_child1.getPartition().getMembersShip(),
		  lconstui_numInstances,
		  (T_CLUSTERIDX) aochrom_child1.getCodeBook().getNumRows()
		  );
	
	       std::vector<T_REAL>&&  lvectort_partialSilhouette =
		 um::simplifiedSilhouette
		 (aimatrixt_centroids,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lpartitionlabel_clusters,
		  aochrom_child1.getPartition().getNumInstancesClusterK(),
		  aifunc2p_dist
		  );

	       T_REAL lmetrict_sumPartialSilhouette = 
		 interfacesse::sum
		 (lvectort_partialSilhouette.data(),
		  (uintidx) lvectort_partialSilhouette.size()
		  );

	       lor_objetiveValue = 
		 (lvectort_partialSilhouette.size() != 0)?
		 lmetrict_sumPartialSilhouette / (T_REAL) lvectort_partialSilhouette.size():
		 measuare_undefSilhouette(T_REAL);
       	
#endif // __FITNESS_SIMPLIFIED_SILHOUETTE__
	

		
#if defined(__FITNESS_SILHOUETTE__)

	       lor_objetiveValue = 
		 um::silhouette
		 (lmatrixtriagT_dissimilarity,
		  aochrom_child1.getPartition()
		  );
	
#endif //__FITNESS_SILHOUETTE__
	
		
#if defined(__FITNESS_INDEX_I__)
		
	       T_REAL lmetric_dk =
		 um::maxDistCjCjp
		 (aimatrixt_centroids, 	    
		  aifunc2p_dist
		  );
	       lor_objetiveValue = ( lpair_distortion.first > 0.0 )?
		 (( lmetric_e1 / lpair_distortion.first ) *  lmetric_dk )
		 / T_REAL(aimatrixt_centroids.getNumRows())
		 : measuare_undefIndexI(T_REAL);
  
	       lor_objetiveValue = std::pow(lor_objetiveValue,airt_p);
		
#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
	
	       T_REAL lmetrict_SSb =
		 um::ssb
		 (aimatrixt_centroids,
		  larray_centroid1,
		  aochrom_child1.getPartition().getNumInstancesClusterK(),
		  aifunc2p_dist
		  );

	       lor_objetiveValue =
		 ( aimatrixt_centroids.getNumRows() >
		   0.0 && lpair_distortion.first > 0.0 )?
		 lmetrict_SSb/(T_REAL(aimatrixt_centroids.getNumRows())
			       * lpair_distortion.first):
		 1.0/measuare_undefWBIndex(T_REAL);
       
#endif //__FITNESS_WB_INDEX__
	

#if defined(__FITNESS_VRC__)

	       T_REAL lmetrict_SSb =
		 um::ssb
		 (aimatrixt_centroids,
		  larray_centroid1,
		  aochrom_child1.getPartition().getNumInstancesClusterK(),
		  aifunc2p_dist
		  );

	       lor_objetiveValue =
		 ( lpair_distortion.first > 0.0 )?
		 (lmetrict_SSb  * (T_REAL(lconstui_numInstances) -
				   T_REAL(aimatrixt_centroids.getNumRows())))
		 / (lpair_distortion.first
		    * (T_REAL(aimatrixt_centroids.getNumRows()) -1) )
		 : measuare_undefIndexI(T_REAL);
	 
#endif //__FITNESS_VRC__
	
	     }
       
	     aochrom_child1.setValidString(lpair_distortion.second);
	     aochrom_child1.setObjetiveFunc(lor_objetiveValue);

	     if ( aochrom_child1.getValidString() == false ) 
	       ++ll_invalidOffspring;


#ifdef __VERBOSE_YES
	     ++geiinparam_verbose;
	     if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	       std::cout
		 << "CROSS CHILD 1 ITER: "
		 <<  aochrom_child1.getCodeBook().getNumRows()
		 << " : " << lor_objetiveValue << std::endl;
	       aochrom_child1.print();
	       std::cout << std::endl;
	     }
	     --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	      
	      lpair_distortion = 
		um::distortion
		(aochrom_child2.getCodeBook(),
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 aochrom_child2.getPartition().getMembersShip(),
		 aifunc2p_dist,
		 [](const data::Instance<T_FEATURE>* aiinst_iter ) -> T_REAL
		 {
		   data::InstanceFreq
		     <T_FEATURE,
		      T_INSTANCE_FREQUENCY
		      >
		     *lptinstfreq_iter =
		     (data::InstanceFreq
		      <T_FEATURE,
		      T_INSTANCE_FREQUENCY
		      >*)
		     aiinst_iter;

		   return  (T_REAL) lptinstfreq_iter->getFrequency();
		   
		 }
		 );

	      if ( lpair_distortion.second == false ) {
		std::pair<bool,T_REAL> lpair_generate  =   
		  clusteringop::reassignCluster
		  (aochrom_child2.getPartition(),
		   aochrom_child2.getCodeBook(),
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   aifunc2p_dist
		   );
		
		lpair_distortion.first = lpair_generate.second;
		lpair_distortion.second = lpair_generate.first;
		
	      }
	      
#if defined(__FITNESS_DB_INDEX__)     
	      lor_objetiveValue = 1.0 / measuare_undefDBindex(T_REAL); 
#endif //__FITNESS_DB_INDEX__
	      
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)     
	      lor_objetiveValue = measuare_undefSilhouette(T_REAL); 
#endif //__FITNESS_SIMPLIFIED_SILHOUETTE__
      
#if defined(__FITNESS_SILHOUETTE__)
	      lor_objetiveValue = measuare_undefSilhouette(T_REAL); 
#endif //__FITNESS_SILHOUETTE__

#if defined(__FITNESS_INDEX_I__)      
	      lor_objetiveValue = measuare_undefIndexI(T_REAL); 
#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
	      lor_objetiveValue = 1.0/measuare_undefWBIndex(T_REAL); 
#endif //__FITNESS_WB_INDEX__
      
#if defined(__FITNESS_VRC__)
	      lor_objetiveValue = measuare_undefVRC(T_REAL); 
#endif //__FITNESS_VRC__)

	      if ( lpair_distortion.second &&
		   aochrom_child2.getCodeBook().getNumRows() >= 2 ) {

		lpair_distortion.first *=
		  ((T_REAL) lconstui_numInstances *
		   (T_REAL) data::Instance<T_FEATURE>::getNumDimensions());
		 
		mat::MatrixRow<T_FEATURE> aimatrixt_centroids
		  (aochrom_child2.getCodeBook().getNumRows(),
		   data::Instance<T_FEATURE>::getNumDimensions(),
		   aochrom_child2.getCodeBook().toArray()
		   );


#if defined(__FITNESS_DB_INDEX__)

	       partition::PartitionLabel
		 <T_CLUSTERIDX>
		 lpartitionlabel_clusters
		 (aochrom_child2.getPartition().getMembersShip(),
		  lconstui_numInstances,
		  (T_CLUSTERIDX) aochrom_child2.getCodeBook().getNumRows() 
		  );
	
	       lor_objetiveValue =
		 1.0 / um::dbindex
		 (aimatrixt_centroids,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lpartitionlabel_clusters,
		  aifunc2p_dist
		  );

#endif //__FITNESS_DB_INDEX__
		
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)

		partition::PartitionLabel
		  <T_CLUSTERIDX>
		  lpartitionlabel_clusters
		  (aochrom_child2.getPartition().getMembersShip(),
		   lconstui_numInstances,
		   (T_CLUSTERIDX) aochrom_child2.getCodeBook().getNumRows() 
		   );
	
		std::vector<T_REAL>&&  lvectort_partialSilhouette =
		  um::simplifiedSilhouette
		  (aimatrixt_centroids,
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   lpartitionlabel_clusters,
		   aochrom_child2.getPartition().getNumInstancesClusterK(),
		   aifunc2p_dist
		   );

		T_REAL lmetrict_sumPartialSilhouette = 
		  interfacesse::sum
		  (lvectort_partialSilhouette.data(),
		   (uintidx) lvectort_partialSilhouette.size()
		   );

		lor_objetiveValue = 
		  (lvectort_partialSilhouette.size() != 0)?
		  lmetrict_sumPartialSilhouette / (T_REAL) lvectort_partialSilhouette.size():
		  measuare_undefSilhouette(T_REAL);
       	
#endif // __FITNESS_SIMPLIFIED_SILHOUETTE__
	
		

#if defined(__FITNESS_SILHOUETTE__)

		lor_objetiveValue = 
		  um::silhouette
		  (lmatrixtriagT_dissimilarity,
		   aochrom_child2.getPartition()
		   );
	
#endif //__FITNESS_SILHOUETTE__
	

#if defined(__FITNESS_INDEX_I__)
	
		T_REAL lmetric_dk =
		  um::maxDistCjCjp
		  (aimatrixt_centroids, 	    
		   aifunc2p_dist
		   );
		lor_objetiveValue = ( lpair_distortion.first > 0.0 )?
		  (( lmetric_e1 / lpair_distortion.first ) *  lmetric_dk )
		  / T_REAL(aimatrixt_centroids.getNumRows())
		  : measuare_undefIndexI(T_REAL);
  
		lor_objetiveValue = std::pow(lor_objetiveValue,airt_p);
	

#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
	
		T_REAL lmetrict_SSb =
		  um::ssb
		  (aimatrixt_centroids,
		   larray_centroid1,
		   aochrom_child2.getPartition().getNumInstancesClusterK(),
		   aifunc2p_dist
		   );

		lor_objetiveValue =
		  ( aimatrixt_centroids.getNumRows() >
		    0.0 && lpair_distortion.first > 0.0 )?
		  lmetrict_SSb /( T_REAL(aimatrixt_centroids.getNumRows())
				  * lpair_distortion.first):
		  1.0/measuare_undefWBIndex(T_REAL);
       
#endif //__FITNESS_WB_INDEX__
		
#if defined(__FITNESS_VRC__)

		T_REAL lmetrict_SSb =
		  um::ssb
		  (aimatrixt_centroids,
		   larray_centroid1,
		   aochrom_child2.getPartition().getNumInstancesClusterK(),
		   aifunc2p_dist
		   );

		lor_objetiveValue =
		  ( lpair_distortion.first > 0.0 )?
		  (lmetrict_SSb  * (T_REAL(lconstui_numInstances) -
				    T_REAL(aimatrixt_centroids.getNumRows())))
		  / (lpair_distortion.first
		     *  (T_REAL(aimatrixt_centroids.getNumRows()) -1) )
		  : measuare_undefIndexI(T_REAL);
	
		
#endif //__FITNESS_VRC__

	      }
	       
	      aochrom_child2.setValidString(lpair_distortion.second);
	      aochrom_child2.setObjetiveFunc(lor_objetiveValue);

#ifdef __VERBOSE_YES
	      ++geiinparam_verbose;
	      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		std::cout
		  << "CROSS CHILD 2 ITER: "
		  <<  aochrom_child2.getCodeBook().getNumRows()
		  << " : "  << lor_objetiveValue << std::endl;
		aochrom_child2.print();
		std::cout << std::endl;
	      }
	      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	       
	      if ( aochrom_child2.getValidString() == false ) 
		++ll_invalidOffspring; 	       
	       
	     } //if  Crossover
	     else {

	       aochrom_child1 =  aichrom_parent1;
	       aochrom_child2 =  aichrom_parent2;
	   
	   }
	 }
	 );

	aoop_outParamGAC.sumTotalInvalidOffspring(ll_invalidOffspring);
	
#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << " OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
      } /*END CROSSOVER*/


   { /*BEGIN PRESERVING THE BEST STRING
       */
#ifdef __VERBOSE_YES
     const char *geverbosepc_labelstep = "ELITISM PRESERVING THE BEST CHROMOSOME";
     ++geiinparam_verbose;
     if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       std::cout
	 << geverbosepc_labelstep
	 << ":  IN(" << geiinparam_verbose << ')'
	 << std::endl;
     }
#endif /*__VERBOSE_YES*/

      auto lit_chromMax = std::max_element
	(lvectorchrom_population.begin(), 
	 lvectorchrom_population.end(), 
	 [](gaencode::ChromCodeBook
	    <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_INSTANCE_FREQUENCY,
	    T_INSTANCES_CLUSTER_K,    
	    T_FEATURE_SUM,
	    T_REAL> & x, 
            gaencode::ChromCodeBook
	    <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_INSTANCE_FREQUENCY,
	    T_INSTANCES_CLUSTER_K,    
	    T_FEATURE_SUM,
	    T_REAL>& y
	    ) 
	 {
	   return x.getObjetiveFunc() < y.getObjetiveFunc(); 
	 }
	 );
      
      
      if ( lochrom_best.getObjetiveFunc() < (*lit_chromMax).getObjetiveFunc() ) {
	lochrom_best = *lit_chromMax;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION 
	 */
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  lochrom_best.print();
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      } /*END IF*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END PRESERVING THE BEST STRING*/      

      
      /*The mutation operator allows exploring different 
	numbers of groups splitting groups, merging groups 
	or keeping groups with a new distribution. 
	It has been called:

	         SMoK-kV mutation (Split, Merged or Kept groups for k variable).
		 
      */
      { /*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "MUTATION";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << geverbosepc_labelstep  
	    << ":  IN(" << geiinparam_verbose << ')'
	    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	long ll_invalidOffspring = 0;

	T_REAL lrt_sumFitness = 
	  std::accumulate
	  (lvectorchrom_population.begin(),
	   lvectorchrom_population.end(),
	   T_REAL(0.0), //VALUE INITIAL
	   [&](T_REAL lT_sumPartial,
	       gaencode::ChromCodeBook
	       <T_FEATURE,
	       T_CLUSTERIDX, //-1, 0, 1, .., K
	       T_INSTANCE_FREQUENCY,
	       T_INSTANCES_CLUSTER_K,    
	       T_FEATURE_SUM,
	       T_REAL>& lchromfixleng_iter)
	   {
	     return lT_sumPartial + lchromfixleng_iter.getObjetiveFunc(); 
	   }
	   ); 
	T_REAL lrt_aveFitness =
	  lrt_sumFitness / T_REAL(lvectorchrom_matingPool.size());

	auto lit_chromMax = std::max_element
	  (lvectorchrom_population.begin(), 
	   lvectorchrom_population.end(), 
	   [](gaencode::ChromCodeBook
	      <T_FEATURE,
	      T_CLUSTERIDX, //-1, 0, 1, .., K
	      T_INSTANCE_FREQUENCY,
	      T_INSTANCES_CLUSTER_K,    
	      T_FEATURE_SUM,
	      T_REAL> & x, 
	      gaencode::ChromCodeBook
	      <T_FEATURE,
	      T_CLUSTERIDX, //-1, 0, 1, .., K
	      T_INSTANCE_FREQUENCY,
	      T_INSTANCES_CLUSTER_K,    
	      T_FEATURE_SUM,
	      T_REAL>& y
	      ) 
	   {
	     return x.getObjetiveFunc() < y.getObjetiveFunc(); 
	   }
	   );
       
	T_REAL lrt_maxFitness = (*lit_chromMax).getObjetiveFunc();
	
	
	for ( auto &&liter_chromCodeBook: lvectorchrom_population ) {

	  T_REAL lT_probMutation = 
	  prob::adaptiveProb
	  (T_REAL(0.5),  /*k2*/
	   T_REAL(0.5),  /*k4*/
	   lrt_maxFitness,
	   lrt_aveFitness, 
	   liter_chromCodeBook.getObjetiveFunc()
	   );

	  /* 
	       Another option evaluated was with a probability 
	       of mutation with a fixed

	       if ( uniformdis_real01(gmt19937_eng) <
	           aiinp_inParamPcPmVk.getProbMutation() ) {
	    
	       It checks if the mutation operator 
		is carried out
	      */

	if ( uniformdis_real01(gmt19937_eng) < lT_probMutation )  
	  { //IF BEGIN  MUTATION

	    

	    //Number of clusters k of the current chromosome
	    T_REAL   lrt_geneNumCluster = (T_REAL) liter_chromCodeBook.getCodeBook().getNumRows(); 

	    /* The Mühlenbein procedure improved determines the number of clusters k'
	      for the new individual.
	    */
	    garealop::muhlenbeinMutation
	      (lrt_geneNumCluster,    //Gene to Mutate,
	       aifeact_ai,
	       aifeact_bi,
	       llfh_listFuntionHist.getDomainUpperBound(), 
	       aiinp_inParamPcPmVk.getNumMaxGenerations(),
	       2.0,
	       8.0
	       );
	    //Conversion to integer  
	    uintidx  lconstui_numClusterFk  = (uintidx) lrt_geneNumCluster;
	    std::pair<bool,T_REAL>  lpair_validPartCVIDistortion;
	    
	    //MERGE CLUSTERS
	    if ( lconstui_numClusterFk
		 < liter_chromCodeBook.getCodeBook().getNumRows() ) {
  
	      mat::MatrixRow<T_FEATURE> lmatrixrowt_chromcbgaM1
		(liter_chromCodeBook.getCodeBook().getNumRows(),
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );

	      interfacesse::copy
		(lmatrixrowt_chromcbgaM1.toArray(),
		 liter_chromCodeBook.getCodeBook().toArray(),
		 lmatrixrowt_chromcbgaM1.getNumElems()
		 );

	      std::vector<T_INSTANCES_CLUSTER_K>
		lvectort_numInstClusterK
		(liter_chromCodeBook.getPartition().getNumInstancesClusterK());
	      
	      do {

		//Get maximum number of Kmeans iterations
		uintidx lui_numKmeansMaxIter =
		  getMaximumNumKmeansIter( lmatrixrowt_chromcbgaM1.getNumRows() );
		
		lmatrixrowt_chromcbgaM1 =
		  gaclusteringop::decreaseM1
		  (lmatrixrowt_chromcbgaM1,
		   liter_chromCodeBook.getPartition().getMembersShip(),
		   lvectort_numInstClusterK,
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   aifunc2p_dist,
		   lui_numKmeansMaxIter
		   );
		
	      } while ( lmatrixrowt_chromcbgaM1.getNumRows() > lconstui_numClusterFk);

	      mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
		*lmatresiablerow_mutation =
		new mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
		(lmatrixrowt_chromcbgaM1.getNumRows(),
		 data::Instance<T_FEATURE>::getNumDimensions(),
		 lmatrixrowt_chromcbgaM1.getNumRows()
		 );
	  
	      interfacesse::copy
		(lmatresiablerow_mutation->toArray(),
		 lmatrixrowt_chromcbgaM1.toArray(),
		 lmatrixrowt_chromcbgaM1.getNumElems()
		 );
	      
	      ds::PartitionLinkedStats
		<T_FEATURE,
		 T_CLUSTERIDX,
		 T_INSTANCE_FREQUENCY,
		 T_INSTANCES_CLUSTER_K,
		 T_FEATURE_SUM
		 >* lpartlinkstats_mutationM3 = 		
		new ds::PartitionLinkedStats
		<T_FEATURE,
		 T_CLUSTERIDX,
		 T_INSTANCE_FREQUENCY,
		 T_INSTANCES_CLUSTER_K,
		 T_FEATURE_SUM
		 >
		(lmatresiablerow_mutation->getNumRows(), 
		 lconstui_numInstances,
		 data::Instance<T_FEATURE>::getNumDimensions(), 
		 lmatresiablerow_mutation->getNumRows() 
		 );

	      lpair_validPartCVIDistortion =
		clusteringop::reassignCluster
		(*lpartlinkstats_mutationM3,
		 liter_chromCodeBook.getPartition().getMembersShip(),
		 *lmatresiablerow_mutation,
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 aifunc2p_dist
		 );

	      liter_chromCodeBook.changeCodeBook(lmatresiablerow_mutation);
	      liter_chromCodeBook.changePartition(lpartlinkstats_mutationM3);
	      liter_chromCodeBook.setOptimalityCBGA( gaencode::OPT_CB );

	    }
	    //SPLIT CLUSTERS
	    else if ( liter_chromCodeBook.getCodeBook().getNumRows()
		      < lconstui_numClusterFk ) {
	      
	       mat::MatrixRow<T_FEATURE> lmatrixrowt_chromcbgaM1
		(liter_chromCodeBook.getCodeBook().getNumRows(),
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );
	       

	      interfacesse::copy
		(lmatrixrowt_chromcbgaM1.toArray(),
		 liter_chromCodeBook.getCodeBook().toArray(),
		 lmatrixrowt_chromcbgaM1.getNumElems()
		 );

	      std::vector<T_INSTANCES_CLUSTER_K>
		lvectort_numInstClusterK
		(liter_chromCodeBook.getPartition().getNumInstancesClusterK());
	      lvectort_numInstClusterK.reserve( lconstui_numClusterFk );

	      do {
	      uintidx lui_numKmeansMaxIter =
		getMaximumNumKmeansIter( lmatrixrowt_chromcbgaM1.getNumRows() );
	     
	      lmatrixrowt_chromcbgaM1 =
		  gaclusteringop::increaseM1
		  (lmatrixrowt_chromcbgaM1,
		   liter_chromCodeBook.getPartition().getMembersShip(),
		   lvectort_numInstClusterK,
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   aifunc2p_dist,
		   lui_numKmeansMaxIter
		   );

	      } while ( lmatrixrowt_chromcbgaM1.getNumRows() < lconstui_numClusterFk);

	      mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
		*lmatresiablerow_mutation =
		new mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
		(lmatrixrowt_chromcbgaM1.getNumRows(),
		 data::Instance<T_FEATURE>::getNumDimensions(),
		 lmatrixrowt_chromcbgaM1.getNumRows()
		 );
	  
	      interfacesse::copy
		(lmatresiablerow_mutation->toArray(),
		 lmatrixrowt_chromcbgaM1.toArray(),
		 lmatrixrowt_chromcbgaM1.getNumElems()
		 );

	      ds::PartitionLinkedStats
		<T_FEATURE,
		 T_CLUSTERIDX,
		 T_INSTANCE_FREQUENCY,
		 T_INSTANCES_CLUSTER_K,
		 T_FEATURE_SUM
		 >* lpartlinkstats_mutationM3 = 		
		new ds::PartitionLinkedStats
		<T_FEATURE,
		 T_CLUSTERIDX,
		 T_INSTANCE_FREQUENCY,
		 T_INSTANCES_CLUSTER_K,
		 T_FEATURE_SUM
		 >
		(lmatresiablerow_mutation->getNumRows(), 
		 lconstui_numInstances, 
		 data::Instance<T_FEATURE>::getNumDimensions(), 
		 lmatresiablerow_mutation->getNumRows() 
		 );

	      lpair_validPartCVIDistortion =
		clusteringop::reassignCluster
		(*lpartlinkstats_mutationM3,
		 liter_chromCodeBook.getPartition().getMembersShip(),
		 *lmatresiablerow_mutation,
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 aifunc2p_dist
		 );

	      liter_chromCodeBook.changeCodeBook(lmatresiablerow_mutation);
	      liter_chromCodeBook.changePartition(lpartlinkstats_mutationM3);
	      liter_chromCodeBook.setOptimalityCBGA( gaencode::OPT_CB );

	      
	      
	    }
	    /*
	      They are equal, direct the search with the same number
	      of groups with other participation
	     */
	    else { //IS EQUAL if (lconstui_numClusterFk == liter_chromCodeBook.getCodeBook().getNumRows() ) {

	      uintidx lui_numKmeansMaxIter =
		getMaximumNumKmeansIter( lconstui_numClusterFk );
	      
	      {  //BEGIN MUTATION K-1 ---------------------------------------------
		
		static std::uniform_int_distribution<uintidx>
		  luniformdis_uiidxInstances0n
		  (0, (uintidx) std::distance(aiiterator_instfirst,aiiterator_instlast) -1);
		uintidx  lui_idxInsRand =
		  luniformdis_uiidxInstances0n(gmt19937_eng);
		//NEW GENE 
		data::Instance<T_FEATURE>* liter_newGeneInstance =
		  *std::next(aiiterator_instfirst,lui_idxInsRand);
		   
		mat::MatrixRow<T_FEATURE> lmatrixrowt_centroidsChrom
		(liter_chromCodeBook.getCodeBook().getNumRows(),
		 data::Instance<T_FEATURE>::getNumDimensions(),
		 liter_chromCodeBook.getCodeBook().toArray()
		 );

		mat::MatrixRow<T_FEATURE> llmatrixrowt_sumInstancesCluster
		    (liter_chromCodeBook.getCodeBook().getNumRows(),
		     data::Instance<T_FEATURE>::getNumDimensions(),
		     liter_chromCodeBook.getPartition().getSumInstancesClusterK().toArray()
		     );
		llmatrixrowt_sumInstancesCluster.initialize(); 

		
		std::vector<T_INSTANCES_CLUSTER_K>&  
		  lvectort_numInstClusterK =
		  liter_chromCodeBook.getPartition().getNumInstancesClusterK(); 
	      
		T_REAL       tmetric_sseKmeans;  
		T_CLUSTERIDX lmcidx_numClusterNull = 0;
	
		const std::vector<T_REAL>&&
		  lvectorT_probDistRouletteWheel =
		  prob::makeDistRouletteWheel
		  (lvectort_numInstClusterK.begin(),
		   lvectort_numInstClusterK.end(),
		   [](const T_INSTANCES_CLUSTER_K&
		      lui_numnstClusterK) -> T_REAL
		   {
		     return  T_REAL(1.0/ ((T_REAL) lui_numnstClusterK));
		   }
		   );

		uintidx lui_geneToMutate = 
		  gaselect::getIdxRouletteWheel
		  (lvectorT_probDistRouletteWheel,
		   uintidx(0)
		   );

		if ( lconstui_numClusterFk > 2 ) {

		  mat::MatrixRow<T_FEATURE> 
		    lmatrixrowt_centroidsChromM1
		    (lconstui_numClusterFk-1, 
		     data::Instance<T_FEATURE>::getNumDimensions()     
		     );

		  if ( lui_geneToMutate == 0 ) {
		    interfacesse::copy
		      (lmatrixrowt_centroidsChromM1.getRow(0),
		       lmatrixrowt_centroidsChrom.getRow(1),
		       (lconstui_numClusterFk-1)
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		  } else if ( lui_geneToMutate == (lconstui_numClusterFk-1) ) {
		    interfacesse::copy
		      (lmatrixrowt_centroidsChromM1.getRow(0),
		       lmatrixrowt_centroidsChrom.getRow(0),
		       (lconstui_numClusterFk-1)
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		  } else { 
		    interfacesse::copy
		      (lmatrixrowt_centroidsChromM1.getRow(0),
		       lmatrixrowt_centroidsChrom.getRow(0),
		       lui_geneToMutate
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		    interfacesse::copy
		      (lmatrixrowt_centroidsChromM1.getRow(lui_geneToMutate),
		       lmatrixrowt_centroidsChrom.getRow(lui_geneToMutate+1),
		       (lconstui_numClusterFk - lui_geneToMutate - 1)
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		  }

		  std::vector<T_INSTANCES_CLUSTER_K> 
		    lvectort_numInstClusterKM1
		    (lconstui_numClusterFk-1,
		     T_INSTANCES_CLUSTER_K(0)
		     );

		  mat::MatrixRow<T_FEATURE_SUM>      
		    lmat_sumInsCentroidsM1
		    (lmatrixrowt_centroidsChromM1.getNumRows(),
		     lmatrixrowt_centroidsChromM1.getNumColumns()
		     );
	    
		  T_REAL tmetric_sseKmeansM1;  
		  clusteringop::kmeansAlreadyInitResample
		    (tmetric_sseKmeansM1,
		     lmatrixrowt_centroidsChromM1,
		     lmat_sumInsCentroidsM1,
		     liter_chromCodeBook.getPartition().getMembersShip(),
		     lvectort_numInstClusterKM1,
		     aiiterator_instfirst,
		     aiiterator_instlast,
		     aifunc2p_dist,
		     (lui_numKmeansMaxIter-1)==0?1:lui_numKmeansMaxIter-1
		     );

		  //------------------------------------------------------------

		  const std::vector<T_REAL>&& lvectort_probDistRouletteWheelM1 =
		    prob::makeDistRouletteWheel
		    (lvectort_numInstClusterKM1.begin(),
		     lvectort_numInstClusterKM1.end(),
		     [](const T_INSTANCES_CLUSTER_K&
			lui_numnstClusterK) -> T_REAL
		     {
		       return  (T_REAL) (lui_numnstClusterK);
		     }
		     );

		  uintidx lui_geneToMutateSplit1 = 
		    gaselect::getIdxRouletteWheel
		    (lvectort_probDistRouletteWheelM1,
		     uintidx(0)
		     );
	  
		  //------------------------------------------------------------

		  uintidx lui_geneToMutateSplit2 = 0;
		  //COPY BACK
		  if ( lui_geneToMutate == 0 ) {
		    interfacesse::copy
		      (lmatrixrowt_centroidsChrom.getRow(1),
		       lmatrixrowt_centroidsChromM1.getRow(0),
		       (lconstui_numClusterFk-1)
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		    lui_geneToMutateSplit2 = lui_geneToMutateSplit1 +1;
		  } else if ( lui_geneToMutate == (lconstui_numClusterFk-1) ) {
		    interfacesse::copy
		      (lmatrixrowt_centroidsChrom.getRow(0),
		       lmatrixrowt_centroidsChromM1.getRow(0),
		       (lconstui_numClusterFk-1)
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		    lui_geneToMutateSplit2 = lui_geneToMutateSplit1;
		  } else {
		    interfacesse::copy
		      (lmatrixrowt_centroidsChrom.getRow(0),
		       lmatrixrowt_centroidsChromM1.getRow(0),
		       lui_geneToMutate
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		    interfacesse::copy
		      (lmatrixrowt_centroidsChrom.getRow(lui_geneToMutate+1),
		       lmatrixrowt_centroidsChromM1.getRow(lui_geneToMutate),
		       (lconstui_numClusterFk - lui_geneToMutate - 1)
		       * data::Instance<T_FEATURE>::getNumDimensions()
		       );
		    if (  lui_geneToMutateSplit1 <  lui_geneToMutate )
		      lui_geneToMutateSplit2 = lui_geneToMutateSplit1;
		    else
		      lui_geneToMutateSplit2 = lui_geneToMutateSplit1 +1;
		  }

		  interfacesse::copy
		    (lmatrixrowt_centroidsChrom.getRow(lui_geneToMutate),
		     lmatrixrowt_centroidsChromM1.getRow(lui_geneToMutateSplit1),
		     data::Instance<T_FEATURE>::getNumDimensions()
		     );


		  std::unordered_set<uintidx>
		    lunordset_iuiCluster;
		  lunordset_iuiCluster.insert(lui_geneToMutate);
		  lunordset_iuiCluster.insert(lui_geneToMutateSplit2); 

		  mat::MatrixRow<T_FEATURE> llmatrixrowt_sumInstancesCluster
		    (liter_chromCodeBook.getCodeBook().getNumRows(),
		     data::Instance<T_FEATURE>::getNumDimensions(),
		     liter_chromCodeBook.getPartition().getSumInstancesClusterK().toArray()
		     );
		  llmatrixrowt_sumInstancesCluster.initialize(); 
	      
		  tmetric_sseKmeans = 
		   clusteringop::updateCentroidsSqrtN
		    (lmcidx_numClusterNull,
		     liter_chromCodeBook.getPartition().getMembersShip(),
		     lmatrixrowt_centroidsChrom,
		     llmatrixrowt_sumInstancesCluster,
		     lvectort_numInstClusterK,
		     lunordset_iuiCluster,
		     aiiterator_instfirst,
		     aiiterator_instlast,
		     aifunc2p_dist
		     );
		}  // END if ( lconstui_numClusterFk > 2 )
		else { //CASE FOR lconstui_numClusterFk == 2
		  //CHROMOSOME IS ONLY TWO GENES
		  uintidx lui_geneTwo = (lui_geneToMutate == 0)?1:0;
		  
		  interfacesse::copy
		    (lmatrixrowt_centroidsChrom.getRow(lui_geneTwo),
		     larray_centroid1, 
		     data::Instance<T_FEATURE>::getNumDimensions()
		     );
	    
		  interfacesse::copy
		    (lmatrixrowt_centroidsChrom.getRow(lui_geneToMutate),
		     liter_newGeneInstance->getFeatures(),
		     data::Instance<T_FEATURE>::getNumDimensions()
		     );

		  lmcidx_numClusterNull =
		  clusteringop::kmeansAlreadyInitResample
		    (tmetric_sseKmeans,
		     lmatrixrowt_centroidsChrom,
		     llmatrixrowt_sumInstancesCluster,
		     liter_chromCodeBook.getPartition().getMembersShip(),
		     lvectort_numInstClusterK,
		     aiiterator_instfirst,
		     aiiterator_instlast,
		     aifunc2p_dist,
		     lui_numKmeansMaxIter
		     );
	    
		} //END MUTATION ANOTHER CASE----------------------------------
	  
		
	      }  //END  MUTATION K-1 ------------------------------------------  
		
		ds::PartitionLinkedStats
		  <T_FEATURE,
		   T_CLUSTERIDX,
		   T_INSTANCE_FREQUENCY,
		   T_INSTANCES_CLUSTER_K,
		   T_FEATURE_SUM
		   >* lpartlinkstats_mutationM3 = 		
		  new ds::PartitionLinkedStats
		  <T_FEATURE,
		   T_CLUSTERIDX,
		   T_INSTANCE_FREQUENCY,
		   T_INSTANCES_CLUSTER_K,
		   T_FEATURE_SUM
		   >
		  (lconstui_numClusterFk, 
		   lconstui_numInstances, 
		   data::Instance<T_FEATURE>::getNumDimensions(), 
		   lconstui_numClusterFk  
		  );

		lpair_validPartCVIDistortion =
		  clusteringop::reassignCluster
		  (*lpartlinkstats_mutationM3,
		   liter_chromCodeBook.getPartition().getMembersShip(),
		   liter_chromCodeBook.getCodeBook(),
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   aifunc2p_dist
		   );

		liter_chromCodeBook.changePartition(lpartlinkstats_mutationM3);
		liter_chromCodeBook.setOptimalityCBGA( gaencode::OPT_CB );
		
	    } // END if 
	  
#ifdef __VERBOSE_YES
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      std::cout
		<< "MUTATION DISTORTION: "
		<< lpair_validPartCVIDistortion.first  << " : "
		<< lpair_validPartCVIDistortion.second << std::endl;
	      std::cout
		<< "MUTATION ITER: "
		<<  liter_chromCodeBook.getCodeBook().getNumRows()
		<< std::endl;
	      liter_chromCodeBook.print();
	      std::cout << std::endl;
	    }
	    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	       
	    if ( lpair_validPartCVIDistortion.first == false ) {
	       
	      lpair_validPartCVIDistortion =
		  clusteringop::reassignCluster
		  (liter_chromCodeBook.getPartition(),
		   liter_chromCodeBook.getCodeBook(),
		   aiiterator_instfirst,
		   aiiterator_instlast,
		   aifunc2p_dist
		   );
	       
	      }

#if defined(__FITNESS_DB_INDEX__)     
	     T_REAL lor_objetiveValue = 1.0 / measuare_undefDBindex(T_REAL); 
#endif //__FITNESS_DB_INDEX__
	     
#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)     
	     T_REAL lor_objetiveValue = measuare_undefSilhouette(T_REAL); 
#endif //__FITNESS_SIMPLIFIED_SILHOUETTE__
      
#if defined(__FITNESS_SILHOUETTE__)
	     T_REAL lor_objetiveValue = measuare_undefSilhouette(T_REAL); 
#endif //__FITNESS_SILHOUETTE__

#if defined(__FITNESS_INDEX_I__)      
	     T_REAL lor_objetiveValue = measuare_undefIndexI(T_REAL); 
#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
      T_REAL lor_objetiveValue = 1.0/measuare_undefWBIndex(T_REAL); 
#endif //__FITNESS_WB_INDEX__	     
	     
#if defined(__FITNESS_VRC__)
	     T_REAL lor_objetiveValue = measuare_undefVRC(T_REAL); 
#endif //__FITNESS_VRC__)
	     
	     if ( lpair_validPartCVIDistortion.first
		  && liter_chromCodeBook.getCodeBook().getNumRows() >= 2 ) {

	       lpair_validPartCVIDistortion.second *=
		 ((T_REAL) lconstui_numInstances *  
		(T_REAL) data::Instance<T_FEATURE>::getNumDimensions());
	       mat::MatrixRow<T_FEATURE> aimatrixt_centroids
		 (liter_chromCodeBook.getCodeBook().getNumRows(),
		  data::Instance<T_FEATURE>::getNumDimensions(),
		  liter_chromCodeBook.getCodeBook().toArray()
		  );

#if defined(__FITNESS_DB_INDEX__)

	       partition::PartitionLabel
		 <T_CLUSTERIDX>
		 lpartitionlabel_clusters
		 (liter_chromCodeBook.getPartition().getMembersShip(),
		  lconstui_numInstances,
		  (T_CLUSTERIDX) liter_chromCodeBook.getCodeBook().getNumRows() 
		  );
	
	       lor_objetiveValue =
		 1.0 / um::dbindex
		 (aimatrixt_centroids,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lpartitionlabel_clusters,
		  aifunc2p_dist
		  );

#endif //__FITNESS_DB_INDEX__
	       

#if defined(__FITNESS_SIMPLIFIED_SILHOUETTE__)

	       partition::PartitionLabel
		 <T_CLUSTERIDX>
		 lpartitionlabel_clusters
		 (liter_chromCodeBook.getPartition().getMembersShip(),
		  lconstui_numInstances,
		  (T_CLUSTERIDX) liter_chromCodeBook.getCodeBook().getNumRows() 
		  );
	
	       std::vector<T_REAL>&&  lvectort_partialSilhouette =
		 um::simplifiedSilhouette
		 (aimatrixt_centroids,
		  aiiterator_instfirst,
		  aiiterator_instlast,
		  lpartitionlabel_clusters,
		  liter_chromCodeBook.getPartition().getNumInstancesClusterK(),
		  aifunc2p_dist
		  );

	       T_REAL lmetrict_sumPartialSilhouette = 
		 interfacesse::sum
		 (lvectort_partialSilhouette.data(),
		  (uintidx) lvectort_partialSilhouette.size()
		  );

	       lor_objetiveValue = 
		 (lvectort_partialSilhouette.size() != 0)?
		 lmetrict_sumPartialSilhouette /
		 (T_REAL) lvectort_partialSilhouette.size()
		 :measuare_undefSilhouette(T_REAL);
       	
#endif // __FITNESS_SIMPLIFIED_SILHOUETTE__

#if defined(__FITNESS_SILHOUETTE__)

	       lor_objetiveValue = 
		 um::silhouette
		 (lmatrixtriagT_dissimilarity,
		  liter_chromCodeBook.getPartition()
		  );
	
#endif //__FITNESS_SILHOUETTE__
	
#if defined(__FITNESS_INDEX_I__)
		
	       T_REAL lmetric_dk =
		 um::maxDistCjCjp
		 (aimatrixt_centroids, 	    
		  aifunc2p_dist
		  );
	       lor_objetiveValue = ( lpair_validPartCVIDistortion.second > 0.0 )?
		 (( lmetric_e1 / lpair_validPartCVIDistortion.second )
		  *  lmetric_dk )
		 / T_REAL(aimatrixt_centroids.getNumRows())
		 : measuare_undefIndexI(T_REAL);
  
	       lor_objetiveValue = std::pow(lor_objetiveValue,airt_p);

#endif //__FITNESS_INDEX_I__

#if defined(__FITNESS_WB_INDEX__)
	
	       T_REAL lmetrict_SSb =
		 um::ssb
		 (aimatrixt_centroids,
		  larray_centroid1,
		  liter_chromCodeBook.getPartition().getNumInstancesClusterK(),
		  aifunc2p_dist
		  );

	       lor_objetiveValue =
		 ( aimatrixt_centroids.getNumRows() > 0.0
		   && lpair_validPartCVIDistortion.second > 0.0 )?
		 lmetrict_SSb /( T_REAL(aimatrixt_centroids.getNumRows())
				 * lpair_validPartCVIDistortion.second)
		 :1.0/measuare_undefWBIndex(T_REAL);
       
#endif //__FITNESS_WB_INDEX__
	
	       
#if defined(__FITNESS_VRC__)

	       T_REAL lmetrict_SSb =
		 um::ssb
		 (aimatrixt_centroids,
		  larray_centroid1,
		  liter_chromCodeBook.getPartition().getNumInstancesClusterK(),
		  aifunc2p_dist
		  );

	       lor_objetiveValue =
		 ( lpair_validPartCVIDistortion.second > 0.0 )?
		 (lmetrict_SSb  * (T_REAL(lconstui_numInstances)
				   - T_REAL(aimatrixt_centroids.getNumRows())))
		 / (lpair_validPartCVIDistortion.second *
		    (T_REAL(aimatrixt_centroids.getNumRows()) -1) )
		 : measuare_undefIndexI(T_REAL);

#endif //__FITNESS_VRC__
	
	     }

	     liter_chromCodeBook.setValidString(lpair_validPartCVIDistortion.first);
	     liter_chromCodeBook.setObjetiveFunc(lor_objetiveValue);

	     if ( liter_chromCodeBook.getValidString() == false ) 
	       ++ll_invalidOffspring;
	   
	  }
	}

	aoop_outParamGAC.sumTotalInvalidOffspring
	  (ll_invalidOffspring);
	
#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << " OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    

      } /*END MUTATION */
	
    } /*END GENETIC OPERATIONS*/

   { /*BEGIN PRESERVING THE BEST STRING
       */
#ifdef __VERBOSE_YES
      const char *geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/


      auto lit_chromMax = std::max_element
	(lvectorchrom_population.begin(), 
	 lvectorchrom_population.end(), 
	 [](gaencode::ChromCodeBook
	    <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_INSTANCE_FREQUENCY,
	    T_INSTANCES_CLUSTER_K,    
	    T_FEATURE_SUM,
	    T_REAL> & x, 
            gaencode::ChromCodeBook
	    <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_INSTANCE_FREQUENCY,
	    T_INSTANCES_CLUSTER_K,    
	    T_FEATURE_SUM,
	    T_REAL>& y
	    ) 
	 {
	   return x.getObjetiveFunc() < y.getObjetiveFunc();
	 }
	 );
      
      if ( lochrom_best.getObjetiveFunc() < (*lit_chromMax).getObjetiveFunc() ) {
	lochrom_best = *lit_chromMax;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION 
	 */
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  lochrom_best.print();
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      } /*END IF*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END PRESERVING THE BEST STRING*/
    
    
    /*MEASUREMENT BEST: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT
    if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
      lofh_DBindex->setValue(lochrom_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectort_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectort_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

         
  } /*END EVOLUTION While*/ 

  uintidx lui_numClusterKBest =
    lochrom_best.getCodeBook().getNumRows();

  delete [] larray_centroid1;
  
  runtime::stop(let_executionTime);
  aoop_outParamGAC.setNumClusterK
    ((T_CLUSTERIDX)lui_numClusterKBest);
  aoop_outParamGAC.setMetricFuncRun
    (lochrom_best.getObjetiveFunc());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamPcPmVk,
       aoop_outParamGAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout
      << lpc_labelAlgGA 
      << ": OUT(" << geiinparam_verbose << ")\n";
    lochrom_best.print();
    std::cout << std::endl;

  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lochrom_best; 
 
} /* END gasgo_vkcentroid */

} /*END gasgo */

#endif /*__GASGO_VKCENTROID_HPP__*/

