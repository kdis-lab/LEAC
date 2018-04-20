/*! \file tgca_vkcentroid.hpp
 *
 * \brief TGCA \cite He:Tan:GAclusteringVarK:TGCA:2012
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the TGCA algorithm based on the paper:\n
 * Hong He and Yonghong Tan. A two-stage genetic algorithm for\n 
 * automatic clustering. Neurocomput., 81:49–59, April 2012.\n
 * <a href="http://dx.doi.org/10.1016/j.neucom.2011.11.001">doi:http://dx.doi.org/10.1016/j.neucom.2011.11.001</a>.\n
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


#ifndef __TGCA_VKCENTROID_HPP__
#define __TGCA_VKCENTROID_HPP__

#include <stdio.h>

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <map>

#include <leac.hpp>
#include "clustering_operator_hierarchical.hpp"
#include "vector_utils.hpp"
#include "container_out.hpp"
#include "inparam_tgca.hpp"
#include "outparam_eaclustering.hpp"

#include "plot_runtime_function.hpp"


/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of genetic and evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace eac {
 

/*! \fn auto tgca_getKSegments(const CONTAINER_ITERINST aiiterator_instfirst, const CONTAINER_ITERINST aiiterator_instlast, partition::PartitionDisjSets<T_CLUSTERIDX> &aipartitionDisjSets_clusters, FUNCTION_GETATTRIBUTE function_getAttribute)
  \brief tgca_getKSegments \cite He:Tan:GAclusteringVarK:TGCA:2012
  \details The attribute with maximum range is partitioned into ki segments by hierarchical agglomerative clustering algorithm. Then initial values of cluster centers are produced through uniformly random selection in every segment. 
  \returns the coordinate of a centroid.
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartitionDisjSets_clusters a partition::PartitionDisjSets of the attribute with maximum range maximum ki segments is partitioned into  by hierarchical agglomerative clustering algorithm.
  \param function_getAttribute a function to obtain the attributes of an instance
*/
template < typename CONTAINER_ITERINST,
	   typename T_CLUSTERIDX,
	   typename FUNCTION_GETATTRIBUTE
	   >
auto
tgca_getKSegments
(const CONTAINER_ITERINST       aiiterator_instfirst,
 const CONTAINER_ITERINST       aiiterator_instlast,
 partition::PartitionDisjSets   
 <T_CLUSTERIDX>                 &aipartitionDisjSets_clusters,
 FUNCTION_GETATTRIBUTE          function_getAttribute
 ) -> mat::MatrixRow<decltype(function_getAttribute(*aiiterator_instfirst))>
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "tgca_getKSegments";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  aiiterator_instfirst[" <<  *aiiterator_instfirst << "]\n"
	      << " input  aiiterator_instlast[" <<  *aiiterator_instlast << "]\n"
	      << " input aipartitionDisjSets_clusters[" << &aipartitionDisjSets_clusters << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES

  typedef decltype(function_getAttribute(*aiiterator_instfirst)) ResultType;
  
  uintidx lui_slinkK =
    uintidx(aipartitionDisjSets_clusters.getNumCluster());
	
  mat::MatrixRow<ResultType> 
    lomatrixrowt_newMinMaxKSegments
    (lui_slinkK,
     uintidx(2)
     );

  /*CALCULATE  KSEGMENTS
   */ 
  for (uintidx lui_iRow = 0;
       lui_iRow < lomatrixrowt_newMinMaxKSegments.getNumRows();
       lui_iRow++)
    {
      lomatrixrowt_newMinMaxKSegments(lui_iRow,0) = std::numeric_limits<ResultType>::max();
      lomatrixrowt_newMinMaxKSegments(lui_iRow,1) = -std::numeric_limits<ResultType>::max();
    }

  uintidx lui_idxInst = 0;
  for  (auto linst_iter = aiiterator_instfirst;
	linst_iter != aiiterator_instlast;
	linst_iter++)
    {
      const ResultType lt_valueDim = function_getAttribute(*linst_iter); 
      if (lomatrixrowt_newMinMaxKSegments
	  (aipartitionDisjSets_clusters.getClusterIdx(lui_idxInst),0) > lt_valueDim) {
	lomatrixrowt_newMinMaxKSegments(aipartitionDisjSets_clusters.getClusterIdx(lui_idxInst),0)
	  = lt_valueDim; 
      }
      if (lomatrixrowt_newMinMaxKSegments
	  (aipartitionDisjSets_clusters.getClusterIdx(lui_idxInst),1) < lt_valueDim) {
	lomatrixrowt_newMinMaxKSegments(aipartitionDisjSets_clusters.getClusterIdx(lui_idxInst),1)
	  = lt_valueDim; 
      }
      ++lui_idxInst;
    }
	
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelKiSegments;
    lostrstream_labelKiSegments << "<KISEGMENTS:" << lpc_labelFunc;
    lomatrixrowt_newMinMaxKSegments.print
      (std::cout,lostrstream_labelKiSegments.str().c_str(),',',';');
    
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lomatrixrowt_newMinMaxKSegments;
} /* End */


/*! \fn gaencode::ChromVariableLength<T_FEATURE,T_METRIC> tgca_vkcentroid(inout::OutParamEAClustering<T_METRIC,
 T_CLUSTERIDX> &aoop_outParamEAC,inout::InParamTGCA<T_CLUSTERIDX,T_METRIC,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamTGCA, const INPUT_ITERATOR  aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist)
  \brief TGCA \cite He:Tan:GAclusteringVarK:TGCA:2012
  \details Implementation of the TGCA algorithm based on \cite He:Tan:GAclusteringVarK:TGCA:2012. Which automatically finds K cluster using the Variance Ratio Criterion (VRC).
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid. Base to following equation:
  \f[
  x_i \in C_j \leftrightarrow  \| x_i - \mu_j \|  \begin{array}{c}min\\  k \end{array}
  \| x_i - \mu_k \|,\; j=1,2,..k,
  \f]
  where \f$m_j\f$, represents the medoid of cluster \f$C_j\f$
  \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
  \param aiinp_inParamTGCA a inout::InParamTGCA parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_METRIC,
           typename T_FEATURE,      //T_STRING,	   
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromVariableLength<T_FEATURE,T_METRIC> 
tgca_vkcentroid
(inout::OutParamEAClustering
 <T_METRIC,
 T_CLUSTERIDX>                         &aoop_outParamEAC,
 inout::InParamTGCA
 <T_CLUSTERIDX,
 T_METRIC,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                &aiinp_inParamTGCA,
 const INPUT_ITERATOR                  aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{

  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  if ( aiinp_inParamTGCA.getNumClusterKMaximum() == 
       INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED )
    aiinp_inParamTGCA.setNumClusterKMaximum
      (std::round(std::sqrt((double)lui_numInstances)));
  
  gaencode::ChromVariableLength<T_FEATURE,T_METRIC> 
    lochrom_best;
    
#ifdef __VERBOSE_YES
  
  /*ID PROC
   */
  geverboseui_idproc = 0;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "tgca_vkcentroid:";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << "  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output inout::OutParamEAClustering&: aoop_outParamEAC[" 
	      << &aoop_outParamEAC << "]\n"
	      << "\t input  InParamClusteringGaProbFk&: aiinp_inParamTGCA[" 
	      << &aiinp_inParamTGCA << "]\n"
              << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t\tPopulation size = " 
	      << aiinp_inParamTGCA.getSizePopulation()
	      << "\n\t\tProbCrossover = " 
	      << aiinp_inParamTGCA.getProbCrossover()
	      << "\n\t\tk-minimum  = " 
	      << aiinp_inParamTGCA.getNumClusterKMinimum()
	      << "\n\t\tk-maximum  = " 
	      << aiinp_inParamTGCA.getNumClusterKMaximum()
	      << "\n\t\trandom-seed = "
	      << aiinp_inParamTGCA.getRandomSeed()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 


  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamTGCA.getNumMaxGenerations(), "Iterations", "Clustering metrics");

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream                           lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_METRIC> *lofh_objetiveFunc = NULL;
  runtime::RuntimeFunctionStat<T_METRIC>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_METRIC>                   lvectort_statfuncObjetiveFunc;
  
  if ( aiinp_inParamTGCA.getWithPlotStatObjetiveFunc() ) {  
    
    lvectort_statfuncObjetiveFunc.reserve
      ( aiinp_inParamTGCA.getSizePopulation());
    //DEFINE FUNCTION
    lofh_objetiveFunc  = new runtime::RuntimeFunctionValue<T_METRIC>
      ("DB", 
       aiinp_inParamTGCA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_objetiveFunc);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_METRIC>
	( (char) li_i,
	  aiinp_inParamTGCA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamEAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamTGCA.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamTGCA.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open
      (aoop_outParamEAC.getFileNameOutPlotStatObjetiveFunc().c_str(),
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
  aoop_outParamEAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/

  runtime::ExecutionTime let_executionTime = runtime::start();

  std::uniform_int_distribution<T_CLUSTERIDX>
    uniformdis_mmcidxMinMaxK
    (aiinp_inParamTGCA.getNumClusterKMinimum(), 
     aiinp_inParamTGCA.getNumClusterKMaximum()
     );
  std::uniform_real_distribution<T_METRIC> uniformdis_real01(0, 1);
  std::uniform_int_distribution<uintidx> uniformdis_idxInstances
    (0,lui_numInstances-1);
    
  /*VARIABLE NEED FOR POPULATION
   */
  std::vector<gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* >  
    lvectorchrom_population;
  
  T_FEATURE *larray_maxFeactures =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  T_FEATURE *larray_minFeactures =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  stats::maxFeatures
    (larray_maxFeactures,
     aiiterator_instfirst,
     aiiterator_instlast
     );

  stats::minFeatures
    (larray_minFeactures,
     aiiterator_instfirst,
     aiiterator_instlast
     ); 
  
  /*1. INITIALIZATION-----------------------------------------------------------
    1. Initialization
    a) Input object number M, iteration number Gm , crossover
    probability Pc , evolutionary population N, and array of the
    numbers of clusters, i.e. K = {k1,k2,...,kN}.
    b) Adopt the maximum attribute range partition method to
    choose initial cluster centers for N times to form N initial
    individuals.
  */

#ifndef __INITIALIZATION_RANDOM_SAMPLING__
  
  uintidx  lui_idxAttMaxRange;
  uintidx* larrayiu_pi = NULL;
  T_METRIC*  larrayrt_lambda = NULL;

  std::vector<std::pair<uintidx,T_METRIC> > lvector_pairIdxLambda;
  std::map<T_CLUSTERIDX,mat::MatrixRow<T_FEATURE> > lmap_kSegments;
  
  try {
    larrayiu_pi     = new uintidx[lui_numInstances];
    larrayrt_lambda = new T_METRIC[lui_numInstances];
  }
  catch (std::bad_alloc& ba) {
    std::cerr << "bad_alloc caught in alg_slink_sibson1973.hpp: " << ba.what() << '\n';
  }

#endif /*__INITIALIZATION RANDOM SAMPLING__*/

  
  
  {/*BEGIN INITIALIZATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "1. POPULATION INITIALIZATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< ',' << "population_size = " << lvectorchrom_population.size()
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    /*CREATE SPACE FOR STORE POPULATION
     */
    lvectorchrom_population.reserve
      (aiinp_inParamTGCA.getSizePopulation());
     
    /*CHROMOSOME INITIALIZATION
     */
#ifdef __INITIALIZATION_RANDOM_SAMPLING__

    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinp_inParamTGCA.getSizePopulation(); 
	 luintidx_i++) 
      {
	/*Generate a number Ki in range Kmin to Kmax
	 */
	uintidx lui_krand =
	  uintidx(uniformdis_mmcidxMinMaxK(gmt19937_eng));
	
	gaencode::ChromVariableLength<T_FEATURE,T_METRIC> *lchrom_new=
	  new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>
	  ( lui_krand * data::Instance<T_FEATURE>::getNumDimensions() );
	
	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroidsChrom
	  (lui_krand,
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   lchrom_new->getString()
	   );

	std::unordered_set<uintidx>&& lunorderedset_idxRandInstances =
	  prob::getWithoutRepeatsSet
	  ( lmatrixrowt_centroidsChrom.getNumRows(),
	    [&]() -> uintidx
	    {
	      return uniformdis_idxInstances(gmt19937_eng);
	    }
	    );
	 
	centroids::centroidsInitialized
	  (lmatrixrowt_centroidsChrom, 
	   aiiterator_instfirst,
	   lunorderedset_idxRandInstances.begin()
	   );

	lchrom_new->setFitness(measuare_undefVRC(T_METRIC));	
        lchrom_new->setObjetiveFunc(measuare_undefVRC(T_METRIC));

	lvectorchrom_population.push_back( lchrom_new );
	
      }
    
#else /*PROPOSED BY THE AUTHOR*/

      /*attribute with maximum range is vm , and its maximum range Rm
	is [vm1,vm2].
      */
    lui_idxAttMaxRange = 0;
    T_FEATURE lrt_maxRange =
      std::abs(larray_maxFeactures[0] - larray_minFeactures[0]);
    for (uintidx lui_l = 1;
	 lui_l < data::Instance<T_FEATURE>::getNumDimensions(); ++lui_l)
      {
	T_FEATURE lrt_lRange =
	  std::abs(larray_maxFeactures[lui_l] - larray_minFeactures[lui_l]);
	if ( lrt_lRange >  lrt_maxRange) {
	  lui_idxAttMaxRange = lui_l;
	  lrt_maxRange = lrt_lRange;
	}
      }

    clusteringop::slink
      (larrayiu_pi,
       larrayrt_lambda,
       aiiterator_instfirst,
       aiiterator_instlast,
       [&](const data::Instance<T_FEATURE>* aiinst_a, const data::Instance<T_FEATURE>* aiinst_b)
       {
	 return std::abs
	   (aiinst_a->getAttribute(lui_idxAttMaxRange) -
	    aiinst_b->getAttribute(lui_idxAttMaxRange)
	    );
       }
       );

    lvector_pairIdxLambda = 
      clusteringop::sortLambda
      (larrayiu_pi,
       larrayrt_lambda,
       lui_numInstances
       );

    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinp_inParamTGCA.getSizePopulation(); 
	 luintidx_i++) 
      {
	/*Generate a number Ki in range Kmin to Kmax
	 */
	T_CLUSTERIDX lmmidx_krand = 
	  uniformdis_mmcidxMinMaxK(gmt19937_eng);
       
	auto lmapitem_kSegments = lmap_kSegments.find(lmmidx_krand);
	
	if ( !(lmapitem_kSegments != lmap_kSegments.end()) ) {

	  partition::PartitionDisjSets<T_CLUSTERIDX>
	    lpartitionDisjSets_clusters
	    (clusteringop::pointerToDisjSets
	     (lvector_pairIdxLambda,
	      larrayiu_pi,
	      (uintidx) lmmidx_krand 
	      )
	     );
	   
	  mat::MatrixRow<T_FEATURE>  
	    lmatrixrowt_newMinMaxKSegments = 
	    tgca_getKSegments
	    (aiiterator_instfirst,
	     aiiterator_instlast,
	     lpartitionDisjSets_clusters,
	     [&](const data::Instance<T_FEATURE>* aiinst_iter)
	     {
	       return aiinst_iter->getAttribute(lui_idxAttMaxRange);
	     }
	     );
	  
	  lmap_kSegments.emplace(lmmidx_krand,lmatrixrowt_newMinMaxKSegments);

	  lmapitem_kSegments = lmap_kSegments.find(lmmidx_krand);	  
	}
	gaencode::ChromVariableLength<T_FEATURE,T_METRIC> *lchrom_new=
	  gaclusteringop::newChromosome
	  ((*lmapitem_kSegments).second,
	   measuare_undefVRC(T_METRIC),
	   measuare_undefVRC(T_METRIC)
	   );
	lvectorchrom_population.push_back( lchrom_new );
      }

#endif /*__INITIALIZATION RANDOM SAMPLING__*/
    

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< ',' << "population_size = " << lvectorchrom_population.size();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END INITIALIZATION*/


  while ( 1 ) {

    /*2. Evaluation of individuals
      a) Obtain new cluster centers by k-means;
      b) Cluster the objects according to new cluster centers and
      calculate the fitness of the individuals on the basis of Eq. (3);
    */

    {/*BEGIN EVALUATION OF INDIVIDUALS*/
      
#ifdef __VERBOSE_YES
      /*ID PROC
       */
      ++geverboseui_idproc;
      geverbosepc_labelstep = "2. EVALUATION OF INDIVIDUALS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << ',' << "population_size = " << lvectorchrom_population.size()
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/


#ifdef __DELETE_EMPTY_CLUSTER__

      std::vector<gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* >  
	lvectorchrom_tmpPopulation;

      lvectorchrom_tmpPopulation.reserve
	( lvectorchrom_population.size() );

#endif /*__DELETE_EMPTY_CLUSTER__*/
                   
      for (auto lchrom_iter: lvectorchrom_population) {
	
        uintidx lui_numClusterK = 
	  lchrom_iter->getStringSize() / data::Instance<T_FEATURE>::getNumDimensions();
	
	if (lui_numClusterK > 1 ) {
	     
	  /*DECODE CHROMOSOME*/
	  mat::MatrixRow<T_FEATURE> 
	    lmatrixrowt_centroidsChrom
	    (lui_numClusterK, 
	     data::Instance<T_FEATURE>::getNumDimensions(),
	     lchrom_iter->getString()
	     );

	  /*a) Obtain new cluster centers by k-means;
	   */
	       
	  std::vector<T_CLUSTERIDX> lvectormcidx_memberShip
	    (lui_numInstances,
	     T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN)
	     );

	  partition::PartitionLabelVector<T_CLUSTERIDX>
	    lpartitionLabelVector_clusters
	    (lvectormcidx_memberShip,
	     (T_CLUSTERIDX) lui_numClusterK 
	     );
	       
	  T_CLUSTERIDX  lmcidx_numClusterNull;
	  uintidx       luintidx_numThreshold;

	  mat::MatrixRow<T_FEATURE_SUM>       
	    lmatrixrowt_sumInstCluster
	    (lui_numClusterK,
	     data::Instance<T_FEATURE>::getNumDimensions(),
	     T_FEATURE_SUM(0)
	     );
	
	  std::vector<T_INSTANCES_CLUSTER_K> 
	    lvectort_numInstClusterK
	    (lui_numClusterK,
	     T_INSTANCES_CLUSTER_K(0)
	     );

	  COMMON_IDOMAIN  lit_numMaxIter = 0;

	  do {
		 
	    /*STEP ASSIGNMENT
	     */
	    luintidx_numThreshold = 
	      clusteringop::reassignCluster 
	      (lvectormcidx_memberShip.data(),
	       lmatrixrowt_centroidsChrom,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );

	    
	    /*STEP UPDATE
	     */

	    lmcidx_numClusterNull =
	      clusteringop::getCentroids
	      (lmatrixrowt_centroidsChrom,
	       lmatrixrowt_sumInstCluster,
	       lvectort_numInstClusterK,
	       lpartitionLabelVector_clusters,
	       aiiterator_instfirst,
	       aiiterator_instlast
	       );
       
	    ++lit_numMaxIter;
		           
	  } while
	      ((luintidx_numThreshold > aiinp_inParamTGCA.getKmeansMinThreshold()) &&
	       ( lit_numMaxIter  <  aiinp_inParamTGCA.getKmeansNumMaxIter() )
	       );
	  
	  lchrom_iter->setValidString(lmcidx_numClusterNull > 0?false:true);


#ifdef __DELETE_EMPTY_CLUSTER__

	  if ( lmcidx_numClusterNull > 0 ) {
	    
	    uintidx lui_numClusterKWithoutNull =
	      lui_numClusterK - (uintidx) lmcidx_numClusterNull;
	    
	    gaencode::ChromVariableLength<T_FEATURE,T_METRIC> *lchrom_withoutNull =
	      new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>
	      ( lui_numClusterKWithoutNull * data::Instance<T_FEATURE>::getNumDimensions());
	    
	    /*DECODE CHROMOSOME*/
	    mat::MatrixRow<T_FEATURE> 
	      lmatrixrowt_centroidsChromWithoutNull
	      ( lui_numClusterKWithoutNull, 
		data::Instance<T_FEATURE>::getNumDimensions(),
		lchrom_withoutNull->getString()
		);
  
	    uintidx _lui_iRowWithoutNull = 0;
	    for (uintidx _lui_iRow = 0;
		 _lui_iRow  < lmatrixrowt_centroidsChrom.getNumRows();
		 ++_lui_iRow)
	      {
		if ( lvectort_numInstClusterK.at(_lui_iRow) != 0 ) {
		  lmatrixrowt_centroidsChromWithoutNull.copyRow
		    (_lui_iRowWithoutNull,
		     lmatrixrowt_centroidsChrom.getRow(_lui_iRow)
		     );
		  ++_lui_iRowWithoutNull;
		}
	      }

	    lchrom_withoutNull->setValidString
	      (lmatrixrowt_centroidsChromWithoutNull.getNumRows() > 1?true:false);
	    lvectorchrom_tmpPopulation.push_back(lchrom_withoutNull);
	    delete lchrom_iter;
       
	  } 
	  else {
	    lvectorchrom_tmpPopulation.push_back(lchrom_iter);
	  } /*if ( lmcidx_numClusterNull > 0 )*/

	    
#endif /*__DELETE_EMPTY_CLUSTER__*/
	  
	} /*if (lui_numClusterK > 1 )*/
	
      } /* for */

#ifdef __DELETE_EMPTY_CLUSTER__

      lvectorchrom_population.clear();
      lvectorchrom_population = std::move(lvectorchrom_tmpPopulation);

#endif /*__DELETE_EMPTY_CLUSTER__*/
      
      
      /* b) Cluster the objects according to new cluster centers and
	 calculate the fitness of the individuals on the basis of Eq. (3);
      */

      for (auto lchrom_iter: lvectorchrom_population) {
	//BEGIN VRC----------------------------------------------

	uintidx lui_numClusterK = 
	  lchrom_iter->getStringSize() / data::Instance<T_FEATURE>::getNumDimensions();
	
	if (  lchrom_iter->getValidString() &&  lui_numClusterK > 1 ) {
	     
	  /*DECODE CHROMOSOME*/
	  mat::MatrixRow<T_FEATURE> 
	    lmatrixrowt_centroidsChrom
	    (lui_numClusterK, 
	     data::Instance<T_FEATURE>::getNumDimensions(),
	     lchrom_iter->getString()
	     );

	  auto lpartition_clusters = 
	    partition::makePartition
	    (lmatrixrowt_centroidsChrom,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     T_CLUSTERIDX(lmatrixrowt_centroidsChrom.getNumRows()),
	     aifunc2p_dist
	     );
	  
	  T_METRIC lT_VRC;
	   
	  lT_VRC =
	    um::VRC
	    (lmatrixrowt_centroidsChrom,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     lpartition_clusters,
	     aifunc2p_dist
	     );

	  if (lT_VRC != measuare_undefVRC(T_METRIC) ) {
	    lchrom_iter->setValidString(true);
	  }
	  else {
	    lchrom_iter->setValidString(false);
	    aoop_outParamEAC.incTotalInvalidOffspring();
	  }
	   
	  lchrom_iter->setObjetiveFunc(lT_VRC); 
	  lchrom_iter->setFitness(lT_VRC);
	    
	}
	     
#ifndef __WITHOUT_PLOT_STAT
	lvectort_statfuncObjetiveFunc.push_back(lchrom_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
      }

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ")\n";
	inout::containerprint
	  (lvectorchrom_population.begin(),
	   lvectorchrom_population.end(),
	   [](const gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* lchrom_iter) -> std::string
	   {
	     std::ostringstream lostrstream_label;
	     uintidx lui_chromiterK =   lchrom_iter->getStringSize()
	       / data::Instance<T_FEATURE>::getNumDimensions();

	     lostrstream_label << lchrom_iter->getObjetiveFunc() << ' ' << lui_chromiterK;
	     return lostrstream_label.str();
	   },
	   std::cout,
	   "<POPULATIONOBJETIVEFUNC",
	   ','
	   );
	std::cout << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END EVALUATION OF INDIVIDUALS
       */

    /*c) Keep the best individuals to the next generation for con-
      tinuous evolution.
    */
    { /*BEGIN PRESERVING THE BEST STRING*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      const auto lit_chromMax =
	*(std::max_element
	  (lvectorchrom_population.begin(), 
	   lvectorchrom_population.end(), 
	   [](const gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* x, 
	      const gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* y
	      ) 
	{  return x->getFitness() < y->getFitness(); }
	   )
	  );
     
      if ( lochrom_best.getFitness() < lit_chromMax->getFitness() ) {
	lochrom_best = *lit_chromMax;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION
	 */
	aoop_outParamEAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamEAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {

	  /*DECODE CHROMOSOME*/
	  mat::MatrixRow<T_FEATURE> 
	    lmatrixrow_centroidsChromBest 
	    ( lochrom_best.getStringSize() / data::Instance<T_FEATURE>::getNumDimensions(),
	      data::Instance<T_FEATURE>::getNumDimensions(),
	      lochrom_best.getString()
	      );

	  std::ostringstream lostrstream_labelCentroids;
	  lostrstream_labelCentroids
	    << "<CENTROIDSCLUSTER:" << geverbosepc_labelstep
	    << ":generation " <<  llfh_listFuntionHist.getDomainUpperBound()
	    << ":lmatrixrow_centroidsChromBest["
	    << &lmatrixrow_centroidsChromBest << ']';
	  lmatrixrow_centroidsChromBest.print
	    (std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
	  std::cout << std::endl;

	  auto lpartitionCentroids_clustersChromBest = 
	    partition::makePartition
	    (lmatrixrow_centroidsChromBest,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     T_CLUSTERIDX(lmatrixrow_centroidsChromBest.getNumRows()),
	     aifunc2p_dist
	     );
	
	  std::ostringstream lostrstream_labelShipBest;
	  lostrstream_labelShipBest
	    << "<MEMBERCLUSTER:" << geverbosepc_labelstep //lpc_labelStep
	    << ":generation "    <<  llfh_listFuntionHist.getDomainUpperBound()
	    << ":lpartitionCentroids_clustersChromBest<>["
	    << &lpartitionCentroids_clustersChromBest << ']';
	  lpartitionCentroids_clustersChromBest.print
	    (std::cout,lostrstream_labelShipBest.str().c_str(),',');
	  
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
    if ( aiinp_inParamTGCA.getWithPlotStatObjetiveFunc() ) {  
      lofh_objetiveFunc->setValue(lochrom_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectort_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectort_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

    
    /*3. Generation of new individuals by genetic operations
      Each individual is subjected to the following conditions:
    */
   
    { /*BEGIN GENETIC OPERATIONS
       */
      
      T_METRIC  lrt_kcon = 0.0;
      
      /*CALCULATE n_sk maximum number of individuals with the 
	same number of clusters in the current population.
      */
      
      { /*BEGIN CALCULATE n_sk*/
	
	uintidx        lui_nsk  = 0;
	
#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "CALCULATE FOR SELECTION n_sk(k,NChrom)";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": IN(" << geiinparam_verbose << ')'
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/
	
	std::map<T_CLUSTERIDX,uintidx> lmap_chromSameNumclusterK;
      	
	for (auto lchrom_iter: lvectorchrom_population) {
	  
	  T_CLUSTERIDX lmmidx_numClusterK =
	    T_CLUSTERIDX
	    (lchrom_iter->getStringSize() / data::Instance<T_FEATURE>::getNumDimensions()
	     );      
	     
	  auto  litermap_numSamesClusterK = 
	    lmap_chromSameNumclusterK.find(lmmidx_numClusterK);
	  if  ( litermap_numSamesClusterK != lmap_chromSameNumclusterK.end() ) {
	    litermap_numSamesClusterK->second++;
	  }
	  else { //NEW NUM CLUSTERK
	    lmap_chromSameNumclusterK.insert
	      (std::pair<T_CLUSTERIDX,uintidx>(lmmidx_numClusterK,uintidx(1)));
	  }
	}

	auto litermap_max = std::max_element
	  (
	   std::begin(lmap_chromSameNumclusterK),
	   std::end(lmap_chromSameNumclusterK),
	   [] (const std::pair<T_CLUSTERIDX,uintidx> &litermap_it1,
	       const std::pair<T_CLUSTERIDX,uintidx> &litermap_it2)
	   {
	     return litermap_it1.second < litermap_it2.second;
	   }
	   );
	
	lui_nsk = litermap_max->second;
	lrt_kcon =
	  (T_METRIC) lui_nsk / (T_METRIC) lvectorchrom_population.size();

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": OUT(" << geiinparam_verbose << ')'
		    << ',' << "lui_nsk = "  << lui_nsk
		    << ',' << "lrt_kcon = " << lrt_kcon
		    << ',' << "lrt_alpha = " << std::exp(-lrt_kcon)
		    << ',' << "Kmode = " << litermap_max->first
		    << '\n';

	  inout::containerprint
	    (lmap_chromSameNumclusterK.begin(),
	     lmap_chromSameNumclusterK.end(),
	     [](const std::pair<T_CLUSTERIDX,uintidx>& lpairmap_sameNumclusterK)
	     -> std::string
	     {
	       std::ostringstream lostrstream_label;
	       lostrstream_label
		 <<  lpairmap_sameNumclusterK.first
		 << ","
		 << lpairmap_sameNumclusterK.second;
	       
	       return lostrstream_label.str();
	     },
	     std::cout,
	     "<CHROMOSEMESAMENUMCLUSTERK",
	     ';'
	     );
	  std::cout << std::endl;
	  
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	      
      } /*END CALCULATE n_sk
	 */
      
      /*2.4.1. Two-stage selection
	The selection process plays an important role of focusing the
	search effort on promising regions in the data space and control-
	ling the speed of convergence. In the clustering process, if the best
	number of clusters is unknwon, Z b can be obtained only when k b
	is found. Therefore, the search of k b should be the first essential
	task to be accomplished. In view of this point, the TGCA adopts
	two-stage selection operation.
      */

      std::vector<gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* > lvectorchrom_matingPool;
      lvectorchrom_matingPool.reserve
	(aiinp_inParamTGCA.getSizePopulation());
      
      {/*BEGIN 2.4.1. TWO-STAGE SELECTION
	*/

	std::vector<T_METRIC> lvectorT_probDistRouletteWheel;
	
	if ( lrt_kcon < 1.0 ) { /*FIRST-STAGE SELECTION*/
	
#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "FIRST-STAGE SELECTION";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": IN(" << geiinparam_verbose << ')'
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << ',' << "stringpool_size = " << lvectorchrom_matingPool.size()
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/

#ifdef __VERBOSE_YES
	  const char *lpc_labelSubStep = "CALCULATE SELECTION PROBABILITY Ps(i) OF Ki";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << lpc_labelSubStep
		      << ": IN(" << geiinparam_verbose << ')'
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/

	  std::map<T_METRIC,uintidx> lmap_chromSameFitness;

	  T_METRIC lrt_alpha = std::exp(-lrt_kcon);
	  T_METRIC lrt_sumfx = 0.0;

	  for (auto lchrom_iter: lvectorchrom_population)   {
	  
	    lrt_sumfx += std::pow(lchrom_iter->getFitness(),lrt_alpha);
	    
	    auto  litermap_numSamesFitness = 
	      lmap_chromSameFitness.find(lchrom_iter->getFitness());
	    
	    if  ( litermap_numSamesFitness != lmap_chromSameFitness.end() ) {
	      litermap_numSamesFitness->second++;
	    }
	    else { /*NEW INDIVIDUAL HAVE THE SAME FITNESS*/
	      lmap_chromSameFitness.insert
		(std::pair<T_METRIC,uintidx>(lchrom_iter->getFitness(),uintidx(1)));
	    }
	  }
  
#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {

	    T_METRIC lrt_testsumfx = 
	      std::accumulate
	      (lmap_chromSameFitness.begin(),
	       lmap_chromSameFitness.end(),
	       T_METRIC(0.0), //VALUE INITIAL
	       [&](T_METRIC airt_sumPartial,
		   const std::pair<T_METRIC,uintidx> &liter_sumDifFitness
		   )
	       {  
		 return airt_sumPartial + std::pow(liter_sumDifFitness.first,lrt_alpha)
		 * (T_METRIC)  liter_sumDifFitness.second;
	       }
	       );
	    
	    bool lbtest_sumFx = (lrt_sumfx == lrt_testsumfx);
	    
	    std::cout << lpc_labelSubStep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << ',' << "lrt_alpha = "  << lrt_alpha
		      <<  "\nCHECK:" << lpc_labelSubStep 
		      << " lrt_sumfx = " <<  lrt_sumfx
		      << " lrt_testsumfx = " <<  lrt_testsumfx << ' ' << lbtest_sumFx
		      << '\n';

	    std::ostringstream lostrstream_labelSameFitness;
	    lostrstream_labelSameFitness
	      << "<POPULATIONSAMEFITNESS:"
	      << lpc_labelSubStep
	      << ":sum," << lrt_sumfx;
	      
	    inout::containerprint
	      (lmap_chromSameFitness.begin(),
	       lmap_chromSameFitness.end(),
	       [](const std::pair<T_METRIC,uintidx>& lpairmap_sameFitness)
	       -> std::string
	       {
		 std::ostringstream lostrstream_label;
		 lostrstream_label
		   << lpairmap_sameFitness.first << ','
		   << lpairmap_sameFitness.second;
		 
		 return lostrstream_label.str();
	       },
	       std::cout,
	       lostrstream_labelSameFitness.str().c_str(),
	       ';'
	       );
	    std::cout << std::endl;
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES
	 END CALCULATE SELECTION PROBABILITY Ps(i) OF Ki";
       */
	  
	  lvectorT_probDistRouletteWheel =
	    prob::makeDistRouletteWheel
	    (lvectorchrom_population.begin(),
	     lvectorchrom_population.end(),
	     [&](const gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* lchrom_iter) -> T_METRIC
	     {
	       auto litermap_numSamesFitness = 
		 lmap_chromSameFitness.find(lchrom_iter->getFitness());
	       if  ( litermap_numSamesFitness == lmap_chromSameFitness.end() )
		 throw std::runtime_error
		   ("tgca_vkcentroid:"
		    " the algorithm failed in two-stage selection");
	       
	       T_METRIC lrt_sumfXjnotinSs =
		 lrt_sumfx -
		 ( std::pow(litermap_numSamesFitness->first,lrt_alpha)
		   * (T_METRIC) litermap_numSamesFitness->second );

	       return (lchrom_iter->getFitness()!= 0.0)? 1.0 /
		 ( T_METRIC(litermap_numSamesFitness->second) +
		   lrt_sumfXjnotinSs/std::pow(lchrom_iter->getFitness(),lrt_alpha)):0.0;
	     }
	     );

	  /*ELITISMO
	   */
	  lvectorchrom_matingPool.push_back
	    (new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>(lochrom_best));

	  for (uintidx luintidx_i = 1; 
	       luintidx_i < aiinp_inParamTGCA.getSizePopulation(); 
	       luintidx_i++) 
	    {      
	      uintidx luintidx_chrom = 
		gaselect::getIdxRouletteWheel
		(lvectorT_probDistRouletteWheel,
		 uintidx(0)
		 );
	      lvectorchrom_matingPool.push_back
		(new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>
		 (*lvectorchrom_population.at(luintidx_chrom))
		 );

	    }

	  for (uintidx lui_i = 0; lui_i < lvectorchrom_population.size(); ++lui_i) {
	    delete  lvectorchrom_population[lui_i];
	  }
	 
	  lvectorchrom_population.clear();
       
#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << ',' << "stringpool_size = " << lvectorchrom_matingPool.size()
		      << '\n';
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES
	 END FIRST-STAGE SELECTION
       */
		
	} else {
	
#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "SECOND-STAGE SELECTION";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": IN(" << geiinparam_verbose << ')'
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << ',' << "stringpool_size = " << lvectorchrom_matingPool.size()
	      //<< ',' << "lui_nsk = " << lui_nsk
		      << ',' << "lrt_kcon = " << lrt_kcon
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/

	  lvectorT_probDistRouletteWheel =
	    prob::makeDistRouletteWheel
	    (lvectorchrom_population.begin(),lvectorchrom_population.end(),
	     [](const gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* lchromfixleng_iter) -> T_METRIC
	     {
	       return lchromfixleng_iter->getFitness();
	     }
	     );


	  /*ELITISMO
	   */
	  lvectorchrom_matingPool.push_back
	    (new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>(lochrom_best));

	  for (uintidx luintidx_i = 1; 
	       luintidx_i < aiinp_inParamTGCA.getSizePopulation(); 
	       luintidx_i++) 
	    {
	    
	      uintidx luintidx_chrom = 
		gaselect::getIdxRouletteWheel
		(lvectorT_probDistRouletteWheel,
		 uintidx(0)
		 );
	    
	      lvectorchrom_matingPool.push_back
		(new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>
		 (*lvectorchrom_population.at(luintidx_chrom))
		 );

	    }

	  for (uintidx lui_i = 0; lui_i < lvectorchrom_population.size(); ++lui_i) {
	    delete  lvectorchrom_population[lui_i];
	  }
	 
	  lvectorchrom_population.clear();


#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << ',' << "stringpool_size = " << lvectorchrom_matingPool.size()
		      << '\n';
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES
	 END FIRST-STAGE SELECTION
       */

	} /*ELSE SECOND-STAGE SELECTION*/ 

      
      } /*END 2.4.1. TWO-STAGE SELECTION
	 */

      
      /*2.4.2. Parallel crossover
       */
      { /*BEGIN CROSSOVER*/

     
#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "PARALLEL CROSSOVER";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": IN(" << geiinparam_verbose << ')'
		    << ',' << "stringpool_size = " << lvectorchrom_matingPool.size()
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	auto &&lvectoriter_subpopulations = 
	  vectorutils::partition
	  (lvectorchrom_matingPool.begin(),
	   lvectorchrom_matingPool.end(),
	   aiinp_inParamTGCA.getNumSubpopulationsCross()
	   );

	auto liter_beginSubpopulation = lvectoriter_subpopulations.begin();
	auto liter_endSubpopulation   = lvectoriter_subpopulations.begin();
      
	for ( ++liter_endSubpopulation;
	      liter_endSubpopulation != lvectoriter_subpopulations.end();
	      ++liter_beginSubpopulation, ++liter_endSubpopulation)
	  {
	    /*BEGIN: Crossover subpopulation
	     */
	    auto lt_numChromSubPopulation =
	      std::distance(*liter_beginSubpopulation,*liter_endSubpopulation);
	  	  
	    auto jchrom_matingPool = *liter_beginSubpopulation;

	    if ( ( lt_numChromSubPopulation % 2 ) != 0 ) {
	      ++jchrom_matingPool;
	    }
	    while ( jchrom_matingPool !=  *liter_endSubpopulation ) { 

	      auto lchrom_individualXi  = jchrom_matingPool;
	      ++jchrom_matingPool;
	      auto lchrom_individualXj  = jchrom_matingPool;
	      ++jchrom_matingPool;

	      if ( uniformdis_real01(gmt19937_eng) //if  Crossover
		   < aiinp_inParamTGCA.getProbCrossover() ) {
	     
		gagenericop::onePointCrossover
		  (*(*lchrom_individualXi),
		   *(*lchrom_individualXj)
		   );
		(*lchrom_individualXi)->setObjetiveFunc(measuare_undefObjetiveFunc(T_METRIC));
		(*lchrom_individualXj)->setObjetiveFunc(measuare_undefObjetiveFunc(T_METRIC));

	      } //if  Crossover
	 
	    } /*while Crossover*/
	  
	    /*END: Crossover subpopulation
	     */

	  } /*for subpopulation*/

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": OUT(" << geiinparam_verbose << ')'
		    << ',' << "stringpool_size = " << lvectorchrom_matingPool.size()
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
      } /*END CROSSOVER*/

      /*2.4.3. Two-stage mutation
	Find the optimal partition in S, the mutation probabilty pm 
	varies with ki and the variation tendency of p m can be also 
	divided into two stages.
      */
    
      { /*(C). TWO-STAGE MUTATION*/

      
	/*CALCULATE n_sk maximum number of individuals with the 
	  same number of clusters in the current population.
	*/
	{ /*BEGIN CALCULATE n_sk*/
	  uintidx        lui_nsk  = 0;
#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "CALCULATE FOR MUTATION n_sk(k,NChrom)";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": IN(" << geiinparam_verbose << ')'
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/
	
	  std::map<T_CLUSTERIDX,uintidx> lmap_chromSameNumclusterK;
      
	  for (auto lchrom_iter: lvectorchrom_matingPool ) {
	    T_CLUSTERIDX lmmidx_numClusterK =
	      T_CLUSTERIDX
	      (lchrom_iter->getStringSize() / data::Instance<T_FEATURE>::getNumDimensions()
	       );      
	    auto litermap_numSamesClusterK = 
	      lmap_chromSameNumclusterK.find(lmmidx_numClusterK);
	    if  ( litermap_numSamesClusterK != lmap_chromSameNumclusterK.end() ) {
	      litermap_numSamesClusterK->second++;
	    }
	    else { //NEW NUM CLUSTERK
	      lmap_chromSameNumclusterK.insert
		(std::pair<T_CLUSTERIDX,uintidx>(lmmidx_numClusterK,uintidx(1)));
	    }
	  }

	  auto litermap_max = std::max_element
	    (
	     std::begin(lmap_chromSameNumclusterK),
	     std::end(lmap_chromSameNumclusterK),
	     [] (const std::pair<T_CLUSTERIDX,uintidx> &litermap_it1,
		 const std::pair<T_CLUSTERIDX,uintidx> &litermap_it2)
	     {return litermap_it1.second < litermap_it2.second; }
	     );
	  lui_nsk = litermap_max->second;
	  lrt_kcon =
	    (T_METRIC) lui_nsk / (T_METRIC) lvectorchrom_matingPool.size();

#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << ',' << "lui_nsk = "  << lui_nsk
		      << ',' << "lrt_kcon = " << lrt_kcon
		      << ',' << "Kmode = " << litermap_max->first
		      << '\n';
	 
	    auto  litermap_it = lmap_chromSameNumclusterK.begin();
	    std::cout << "<CHROMOSEMESAMENUMCLUSTERK:"
		      << geverbosepc_labelstep  << ",length," << lmap_chromSameNumclusterK.size() << '>';
	    if ( litermap_it != lmap_chromSameNumclusterK.end()) {
	      std::cout << litermap_it->first << "," << litermap_it->second;
	      ++litermap_it;
	    }
	    for (;litermap_it != lmap_chromSameNumclusterK.end(); ++litermap_it) {
	      std::cout << ';' << litermap_it->first << "," << litermap_it->second;
	    }
	    std::cout << std::endl;
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	      
	} /*END CALCULATE n_sk*/

      
	if (lrt_kcon <= 0.9) { /*FIRST-STAGE MUTATION*/

	
	  T_METRIC lrt_pm = 0.1 * ( 1.0 - lrt_kcon );
	
#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "(C). FIRST-STAGE MUTATION";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ":  IN(" << geiinparam_verbose << ')'
		      << ',' << "stringpoolcross_size = " << lvectorchrom_matingPool.size()
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << ',' << "lrt_pm = " << lrt_pm
		      << ',' << "lrt_kcon = " << lrt_kcon
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/

	
	  for (auto lchrom_iter: lvectorchrom_matingPool) {

	    if ( uniformdis_real01(gmt19937_eng) < lrt_pm ) {
	    
	      T_CLUSTERIDX lmmidx_krandNew = 
		uniformdis_mmcidxMinMaxK(gmt19937_eng);
	    
#ifdef __INITIALIZATION_RANDOM_SAMPLING__
	    
	      gaencode::ChromVariableLength<T_FEATURE,T_METRIC> *lchrom_new =
		new gaencode::ChromVariableLength<T_FEATURE,T_METRIC>
		( uintidx(lmmidx_krandNew) * data::Instance<T_FEATURE>::getNumDimensions());

	      /*DECODE CHROMOSOME*/
	      mat::MatrixRow<T_FEATURE> 
		lmatrixrowt_centroidsChromNew
		( lchrom_new->getStringSize() / data::Instance<T_FEATURE>::getNumDimensions(),
		  data::Instance<T_FEATURE>::getNumDimensions(),
		  lchrom_new->getString()
		  );

	      std::unordered_set<uintidx>&& lunorderedset_idxRandInstances =
		prob::getWithoutRepeatsSet
		( lmatrixrowt_centroidsChromNew.getNumRows(),
		  [&]() -> uintidx
		  {
		    return uniformdis_idxInstances(gmt19937_eng);
		  }
		  );
	    
	      clusteringop::centroidsInitialized
		(lmatrixrowt_centroidsChromNew, 
		 aivectorptinst_instances,
		 lunorderedset_idxRandInstances.begin()
		 );

	    
#else /*PROPOSED BY THE AUTHOR*/

	      auto lmapitem_kSegments = lmap_kSegments.find(lmmidx_krandNew);
	
	      if ( !(lmapitem_kSegments != lmap_kSegments.end()) ) {

		partition::PartitionDisjSets<T_CLUSTERIDX>
		  lpartitionDisjSets_clusters
		  (clusteringop::pointerToDisjSets
		   (lvector_pairIdxLambda,
		    larrayiu_pi,
		    (uintidx) lmmidx_krandNew
		    )
		   );
	   
		mat::MatrixRow<T_FEATURE>  
		  lmatrixrowt_newMinMaxKSegments = 
		  tgca_getKSegments
		  (aiiterator_instfirst,
		   aiiterator_instlast,
		   lpartitionDisjSets_clusters,
		   [&](const data::Instance<T_FEATURE>* aiinst_iter)
		   {
		     return aiinst_iter->getAttribute(lui_idxAttMaxRange);
		   }
		   );
	  
		lmap_kSegments.emplace(lmmidx_krandNew,lmatrixrowt_newMinMaxKSegments);

		lmapitem_kSegments = lmap_kSegments.find(lmmidx_krandNew);
	    
	      } //THERE IS NO SEGMENT Ki

	      gaencode::ChromVariableLength<T_FEATURE,T_METRIC> *lchrom_new=
		gaclusteringop::newChromosome
		((*lmapitem_kSegments).second,
		 measuare_undefVRC(T_METRIC),
		 measuare_undefVRC(T_METRIC)
		 );
	     
#endif /*__INITIALIZATION RANDOM SAMPLING__*/

	      lvectorchrom_population.push_back(lchrom_new);  
	      delete lchrom_iter;
	  
	    } else {
	      lvectorchrom_population.push_back(lchrom_iter);
	    } //NOT MUTATION
	 
	  }
	
	  lvectorchrom_matingPool.clear();

#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << ',' << "stringpoolcross_size = " << lvectorchrom_matingPool.size()
		      << ',' << "population_size = " << lvectorchrom_population.size();
	    std::cout << std::endl;
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
			
	} else { /*SECOND STAGE MUTATION*/

	  for (auto lchrom_iter: lvectorchrom_population) {
	  
	    if ( measuare_undefObjetiveFunc(T_METRIC) == lchrom_iter->getObjetiveFunc() ) {

	      /*DECODE CHROMOSOME*/
	      mat::MatrixRow<T_FEATURE> 
		lmatrixrowt_centroidsChrom
		(lchrom_iter->getStringSize() / data::Instance<T_FEATURE>::getNumDimensions(),
		 data::Instance<T_FEATURE>::getNumDimensions(),
		 lchrom_iter->getString()
		 );
	    
	      auto lpartition_clusters = 
		partition::makePartition
		(lmatrixrowt_centroidsChrom,
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 T_CLUSTERIDX(lmatrixrowt_centroidsChrom.getNumRows()),
		 aifunc2p_dist
		 );
	    
	      T_METRIC lT_VRC;
	   
   
	      lT_VRC =
		um::VRC
		(lmatrixrowt_centroidsChrom,
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 lpartition_clusters,
		 aifunc2p_dist
		 );

	      if (lT_VRC != measuare_undefVRC(T_METRIC) ) {
		lchrom_iter->setValidString(true);
	      }
	      else {
		lchrom_iter->setValidString(false);
		aoop_outParamEAC.incTotalInvalidOffspring();
	      }
	   
	      lchrom_iter->setObjetiveFunc(lT_VRC); 
	      lchrom_iter->setFitness(lT_VRC);  
	    }
	  }
	  /*The article does not establish if the fitness function should be 
	    recalculated when the chromosomes crossed. crossing implies that 
	    others are obtained at this time is not being recalculated
	  */
	  T_METRIC lrt_sumfx = 
	    std::accumulate
	    (lvectorchrom_matingPool.begin(),
	     lvectorchrom_matingPool.end(),
	     T_METRIC(0.0),
	     [&](T_METRIC airt_sumPartial, const gaencode::ChromVariableLength<T_FEATURE,T_METRIC>* lchrom_iter)
	     {  
	       return airt_sumPartial + lchrom_iter->getFitness();
	     }
	     );
	
#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "(C). SECOND-STAGE MUTATION";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ":  IN(" << geiinparam_verbose << ')'
		      << ',' << "stringpoolcross_size = " << lvectorchrom_matingPool.size()
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << ',' << "lrt_kcon = " << lrt_kcon
		      << ',' << "lrt_sumfx = " << lrt_sumfx
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/

	  for (auto lchrom_iter: lvectorchrom_matingPool) {
	    T_METRIC lrt_pm = std::exp(-lchrom_iter->getFitness() / lrt_sumfx);
	    if ( uniformdis_real01(gmt19937_eng) < lrt_pm ) {
	      gaclusteringop::biDirectionMutation
		(*lchrom_iter,
		 larray_minFeactures,
		 larray_maxFeactures
		 ); 
	    } //END IF MUTATION
	    lvectorchrom_population.push_back(lchrom_iter);
	  }
	  lvectorchrom_matingPool.clear();
	
#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << ',' << "stringpoolcross_size = " << lvectorchrom_matingPool.size()
		      << ',' << "population_size = " << lvectorchrom_population.size()
		      << std::endl;
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
	} /*END ELSE*/


      } /*END MUTATION */

    	
    } /*END GENETIC OPERATIONS*/
   
   
    { /*BEGIN TERMINATION CRITERION*/
#ifdef __VERBOSE_YES
      
      /*ID PROC
       */
      ++geverboseui_idproc;
      
      const char *lpc_labelStep = "TERMINATION CRITERION ATTAINED:";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << lpc_labelStep
		  << "generation " << llfh_listFuntionHist.getDomainUpperBound()
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
      if ( !(llfh_listFuntionHist.getDomainUpperBound() 
	     < aiinp_inParamTGCA.getNumMaxGenerations() ) 
	   )
	break;
      llfh_listFuntionHist.increaseDomainUpperBound();
    } /*END 3.1.5 TERMINATION CRITERION*/
    
  } /*END EVOLUTION While*/ 

  for (auto  liter_chrom: lvectorchrom_population)
    delete liter_chrom;
  
  runtime::stop(let_executionTime);
  aoop_outParamEAC.setNumClusterK
    (T_CLUSTERIDX(lochrom_best.getStringSize() / data::Instance<T_FEATURE>::getNumDimensions())
     );
  aoop_outParamEAC.setMetricFuncRun
    (lochrom_best.getObjetiveFunc());
  aoop_outParamEAC.setFitness
    (lochrom_best.getFitness());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamTGCA.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamTGCA,
       aoop_outParamEAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES

  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA 
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochrom_best.print();
    std::cout << std::endl;
 
    /*DECODE CHROMOSOME*/
    mat::MatrixRow<T_FEATURE> 
      lmatrixrow_centroidsChromBest 
      ( lochrom_best.getStringSize() / data::Instance<T_FEATURE>::getNumDimensions(),
	data::Instance<T_FEATURE>::getNumDimensions(),
	lochrom_best.getString()
	);

    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDSCLUSTER: " << lpc_labelAlgGA
			       << ": generation " <<  llfh_listFuntionHist.getDomainUpperBound()
			       << ": lmatrixrow_centroidsChromBest["
			       << &lmatrixrow_centroidsChromBest << ']';
    lmatrixrow_centroidsChromBest.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    std::cout << std::endl;

    auto lpartitionCentroids_clustersChromBest = 
      partition::makePartition
      (lmatrixrow_centroidsChromBest,
       aiiterator_instfirst,
       aiiterator_instlast,
       T_CLUSTERIDX(lmatrixrow_centroidsChromBest.getNumRows()),
       aifunc2p_dist
       );
    
    std::ostringstream lostrstream_labelShipBest;
    lostrstream_labelShipBest
      << "<MEMBERCLUSTER:" << lpc_labelAlgGA
      << "generation "    <<  llfh_listFuntionHist.getDomainUpperBound()
      << ":lpartitionCentroids_clustersChromBest<>["
      << &lpartitionCentroids_clustersChromBest << ']';
    lpartitionCentroids_clustersChromBest.print
      (std::cout,lostrstream_labelShipBest.str().c_str(),',');
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

#ifndef __INITIALIZATION_RANDOM_SAMPLING__
  if ( larrayiu_pi != NULL )
    delete [] larrayiu_pi;
  if ( larrayiu_pi != NULL )
    delete [] larrayrt_lambda;
#endif /*__INITIALIZATION RANDOM SAMPLING__*/
  
  /*DELETE VARIBLES EXTRAS OF THE ALGORITHM
   */
  delete[] larray_minFeactures;
  delete[] larray_maxFeactures;

  return lochrom_best; 
 
} /* END tgca_vkcentroid
   */

} /*END eac */

#endif /*__TGCA_VKCENTROID_HPP__*/
