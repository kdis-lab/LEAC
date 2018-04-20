/*! \file gaclustering_fklabel.hpp
 *
 * \brief GA Clustering \cite Murthy:Chowdhury:GAclustering:GA:1996
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GA algorithm based on the paper:\n
 * C. A. Murthy and Nirmalya Chowdhury. In search of optimal\n 
 * clusters using genetic algorithms. Pattern Recogn. Lett.,\n
 * 17(8):825–832, 1996.\n
 * <a href="http://dx.doi.org/10.1016/0167-8655(96)00043-8">doi:http://dx.doi.org/10.1016/0167-8655(96)00043-8</a>\n.
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

#ifndef __GA_MURTHY_CHOWDHURY_1996_HPP__
#define __GA_MURTHY_CHOWDHURY_1996_HPP__

#include <vector>
#include <algorithm>
#include <random>
#include <unordered_set>

#include <leac.hpp>

#include "inparam_pcpmfk.hpp"
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
  
/*! \fn gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> gaclustering_fklabel(inout::OutParamEAClustering<T_REAL,T_CLUSTERIDX> &aoop_outParamEAC, inout::InParamGAClusteringPcPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinpcgaprobfixedk_inParamGA, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR  aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist)
 \brief GA Clustering \cite Murthy:Chowdhury:GAclustering:GA:1996
 \details Implementation of GA algorithm based on  \cite Murthy:Chowdhury:GAclustering:GA:1996. 
 \returns A partition of a data set, encoded on a chromosome where each gene is the index of a cluster to which the instance belongs.
 \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
 \param aiinpcgaprobfixedk_inParamGA a inparam::InParamGAClusteringPcPmFk parameters required by the algorithm
 \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
 \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
 \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,   
	   typename T_REAL, //OBJETIVE FUNCTION,  FITNESS AND CLUSTERING_METRIC       
	   typename T_FEATURE,        
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename INPUT_ITERATOR
	   >
gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> 
gaclustering_fklabel
(inout::OutParamEAClustering
 <T_REAL,
 T_CLUSTERIDX>                   &aoop_outParamEAC,
 inout::InParamPcPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,         
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>        &aiinpcgaprobfixedk_inParamGA,
 const INPUT_ITERATOR          aiiterator_instfirst,
 const INPUT_ITERATOR          aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist
 )
{
  
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));

  /*ASIGNED CHROMOSOME SIZE
   */
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::setStringSize
    (lui_numInstances);
  
  /*DISTRIBUTION FOR RANDOM NUMBER
   */
  std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_mmcidx0K
    (0,aiinpcgaprobfixedk_inParamGA.getNumClusterK()-1);
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
  /*CHROMOSOME BEST 
   */
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> lochromfixleng_best;

  /*SPACE FOR STORE POPULATION
   */
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* >  
    lvectorchromfixleng_population;

  /*SPACE FOR STORE MATINGPOOL
   */
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* >  
    lvectorchromfixleng_matingPool;


#ifdef __VERBOSE_YES

   /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gaclustering_fklabel"; 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\t\n(output gaencode::ChromFixedLength: lochromfixleng_best[" 
	      << &lochromfixleng_best << "]\n"
	      << "\t output inout::OutParamEAClustering&: aoop_outParamEAC[" 
	      << &aoop_outParamEAC << "]\n"
	      << "\t input  InParamClusteringGaProbFk&: aiinpcgaprobfixedk_inParamGA[" 
	      << &aiinpcgaprobfixedk_inParamGA << "]\n"
	      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\nGA parameters: "
	      << "\tPopulation size = " << aiinpcgaprobfixedk_inParamGA.getSizePopulation()
	      << "\tProbCrossover = "   << aiinpcgaprobfixedk_inParamGA.getProbCrossover() 
	      << "\tProbMutation  = "   << aiinpcgaprobfixedk_inParamGA.getProbMutation()
	      << "\tKFind = "           << aiinpcgaprobfixedk_inParamGA.getNumClusterK()
	      << "\n\t)"
	      << std::endl;
  }
  #endif /*__VERBOSE_YES*/ 


  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinpcgaprobfixedk_inParamGA.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

 
  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  
  std::ofstream lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;
  
  if ( aiinpcgaprobfixedk_inParamGA.getWithPlotStatObjetiveFunc() ) {  

    lvectorT_statfuncObjetiveFunc.reserve
      (aiinpcgaprobfixedk_inParamGA.getSizePopulation());  
    //DEFINE FUNCTION
    lofh_SSE  = new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinpcgaprobfixedk_inParamGA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat<T_REAL>
	( (char) li_i,
	  aiinpcgaprobfixedk_inParamGA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION 
    aoop_outParamEAC.setFileNameOutPlotStatObjetiveFunc
      (aiinpcgaprobfixedk_inParamGA.getFileNamePlotStatObjetiveFunc(),
       aiinpcgaprobfixedk_inParamGA.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open   
      (aoop_outParamEAC.getFileNameOutPlotStatObjetiveFunc().c_str(), 
       std::ios::out | std::ios::app
       );
    lfileout_plotStatObjetiveFunc.precision(COMMON_COUT_PRECISION);

    //FUNCTION HEADER
    lfileout_plotStatObjetiveFunc 
      << llfh_listFuntionHist.getHeaderFuntions() 
      << "\n";
  }

#endif /*__WITHOUT_PLOT_STAT*/

 
  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING*/
  aoop_outParamEAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/
  runtime::ExecutionTime let_executionTime = runtime::start();


  /*POPULATION CREATE-----------------------------------------------------------
   */
  lvectorchromfixleng_population.reserve
    (aiinpcgaprobfixedk_inParamGA.getSizePopulation());
  for (uintidx luintidx_i = 0; 
       luintidx_i < aiinpcgaprobfixedk_inParamGA.getSizePopulation(); 
       luintidx_i++) 
    {
      lvectorchromfixleng_population.push_back
	( new gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>() );
    }
 
  /*CREATE SPACE FOR STORE MATINGPOOL-------------------------------------------
   */
  lvectorchromfixleng_matingPool.reserve
    (aiinpcgaprobfixedk_inParamGA.getSizePopulation());
  for (uintidx luintidx_i = 0; 
       luintidx_i < aiinpcgaprobfixedk_inParamGA.getSizePopulation(); 
       luintidx_i++) 
    {
      lvectorchromfixleng_matingPool.push_back
	( new gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>() );
    }
 
  /*POPULATION INITIAL----------------------------------------------------------
   */
  { /*BEGIN INITIALIZATION*/
   
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "(0) POPULATION INITIAL";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/  

    for (auto lchromfixleng_iter :lvectorchromfixleng_population) {
      
      gagenericop::initializeGenes
	(lchromfixleng_iter->begin(),
	 lchromfixleng_iter->end(),
	 [&]() 
	 {
	   return uniformdis_mmcidx0K(gmt19937_eng);
	 }
	 );
      
      lchromfixleng_iter->setObjetiveFunc(std::numeric_limits<T_REAL>::max());
      lchromfixleng_iter->setFitness(-std::numeric_limits<T_REAL>::max());
      
    }
    

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END INITIALIZATION*/

 
  /*COMPUTING CLUSTERING METRIC AND FITNESS COMPUTATION-------------------------
   */

  { /*BEGIN FITNESS COMPUTATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "COMPUTING CLUSTERING METRIC AND FITNESS COMPUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

    long ll_invalidOffspring = 0;

    mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroids
	((uintidx) aiinpcgaprobfixedk_inParamGA.getNumClusterK(),
	 data::Instance<T_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<T_FEATURE_SUM>       
	lmatrixrowt_sumInstCluster
	((uintidx) aiinpcgaprobfixedk_inParamGA.getNumClusterK(), 
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
	
      std::vector<T_INSTANCES_CLUSTER_K> 
	lvectort_numInstClusterK
	((uintidx) aiinpcgaprobfixedk_inParamGA.getNumClusterK());
      
      for (auto lchromfixleng_iter :lvectorchromfixleng_population) {

	partition::PartitionLabel<T_CLUSTERIDX>
	  lpartition_clusters
	  (lchromfixleng_iter->getString(),
	   lui_numInstances,
	   aiinpcgaprobfixedk_inParamGA.getNumClusterK()
	   );
    
	/*TRANSFORM OF LABEL TO CENTROID
	 */
	T_CLUSTERIDX lmcidx_numClusterNull =
	  clusteringop::getCentroids
	  (lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   lpartition_clusters,
	   aiiterator_instfirst,
	   aiiterator_instlast
	   );

	/*EVALUATE OBJETIVE FUNCTION
	 */ 
	T_REAL lT_objetiveFunc;
	if ( lmcidx_numClusterNull == 0 ) { 
	  lT_objetiveFunc = 
	    um::SSE
	    (lmatrixrowt_centroids,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     lchromfixleng_iter->getString(),
	     aifunc2p_dist
	     );
	}
	else {
	  lT_objetiveFunc = std::numeric_limits<T_REAL>::max();
	}

	lchromfixleng_iter->setObjetiveFunc(lT_objetiveFunc);
	lchromfixleng_iter->setValidString(lmcidx_numClusterNull == 0);
	
	lchromfixleng_iter->setFitness(1.0 / lchromfixleng_iter->getObjetiveFunc());
	
	if ( lchromfixleng_iter->getValidString() == false )
	  ++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

      }
      
      aoop_outParamEAC.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END FITNESS COMPUTATION*/
     
  /*(a) COPY THE BEST STRING TO S0----------------------------------------------
   */

  { /*BEGIN (a) COPY THE BEST STRING TO S0*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "(a) COPY THE BEST STRING TO S0";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* lchromfixleng_max = 
      *(std::max_element
	(lvectorchromfixleng_population.begin(), 
	 lvectorchromfixleng_population.end(), 
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* x, 
	    const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* y
	    ) 
	 {  return x->getFitness() < y->getFitness(); }
	 ));
      
     if ( lochromfixleng_best.getFitness() <  lchromfixleng_max->getFitness() ) {
       
       lochromfixleng_best  = *lchromfixleng_max;
    
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	aoop_outParamEAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamEAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

      } //END IF lochromfixleng_best < lchrom_maxFitness
    

#ifdef __VERBOSE_YES
     if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
 
  } /*END (a) COPY THE BEST STRING TO S0*/


  /*INITIAL MEASUREMENT: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT  
  if ( aiinpcgaprobfixedk_inParamGA.getWithPlotStatObjetiveFunc() ) {  
    lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
    functionhiststat_evaluateAll
      (lofhs_statObjectiveFunc,
       lvectorT_statfuncObjetiveFunc
       );
    lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
    lvectorT_statfuncObjetiveFunc.clear();
  }
#endif  /*__WITHOUT_PLOT_STAT */


 
  /*ITERATION: STEPS (B), (C), AND (D)
   */
  while ( (llfh_listFuntionHist.getDomainUpperBound() < 
	   aiinpcgaprobfixedk_inParamGA.getNumMaxGenerations()) &&
	  (runtime::elapsedTime(let_executionTime) < aiinpcgaprobfixedk_inParamGA.getMaxExecutiontime() ) ) {
   
    llfh_listFuntionHist.increaseDomainUpperBound();

    /*(B) PERFORM SELECTION, CROSSOVER AND MUTATION OBTAIN NEW Q1---------------
     */

    /*SELECTION--
     */

    { //BEGIN SELECTION

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(B) SELECTION: COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* lchromfixleng_iter) -> T_REAL
	 {
	   return lchromfixleng_iter->getFitness();
	 }
	 );
	      
      /*COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL--
       */
      for ( auto lchromfixleng_iter: lvectorchromfixleng_matingPool) {

	uintidx luiidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );
	
	*lchromfixleng_iter = *lvectorchromfixleng_population.at(luiidx_chrom);
      }
         
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


    } /*END SELECTION*/


    /*CROSSOVER-----------------------------------------------------------------
     */
    { /*BEGIN CROSSOVER*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(B) CROSSOVER";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      long ll_invalidOffspring = 0;

     gaiterator::crossover
      (lvectorchromfixleng_matingPool.begin(),
       lvectorchromfixleng_matingPool.end(),
       lvectorchromfixleng_population.begin(),
       lvectorchromfixleng_population.end(),
       [&](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* aichrom_parent1,
	   const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* aichrom_parent2,
	   gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>*  aochrom_child1, 
	   gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>*  aochrom_child2
	   )
       {

	 if ( uniformdis_real01(gmt19937_eng) <
	      aiinpcgaprobfixedk_inParamGA.getProbCrossover()  ) {
	   
	    gagenericop::onePointCrossover
	      (*aochrom_child1,
	       *aochrom_child2,
	       *aichrom_parent1,
	       *aichrom_parent2
	       );
	    aochrom_child1->setObjetiveFunc
	      (std::numeric_limits<T_REAL>::max());
	    aochrom_child2->setObjetiveFunc
	      (std::numeric_limits<T_REAL>::max());
	   
	    /*CHECK IF STRINGS ARE VALID
	     */
	    bool  lb_isNotValidChild1 = 
	      gaintegerop::isNotValid
	      (*aochrom_child1,
	       (T_CLUSTERIDX) 0,
	       aiinpcgaprobfixedk_inParamGA.getNumClusterK()
	       );
	    bool  lb_isNotValidChild2 =
	      gaintegerop::isNotValid
	      (*aochrom_child2,
	       (T_CLUSTERIDX) 0,
	       aiinpcgaprobfixedk_inParamGA.getNumClusterK()
	       );
	 
	    if (lb_isNotValidChild1) ++ll_invalidOffspring;
	    if (lb_isNotValidChild2) ++ll_invalidOffspring;

	    uintidx  lst_countlimitReachedCrossover = 0;

	    /* 100 IS THE CONSTANT PROPOSED BY THE AUTHOR
	      */
	    while( (lst_countlimitReachedCrossover < (uintidx) 100) && 
		   (lb_isNotValidChild1 || lb_isNotValidChild2) ) {
	   
	      gagenericop::onePointCrossover
		(*aochrom_child1,
		 *aochrom_child2,
		 *aichrom_parent1,
		 *aichrom_parent2
		 );
	      aochrom_child1->setObjetiveFunc
		(std::numeric_limits<T_REAL>::max());
	      aochrom_child2->setObjetiveFunc
		(std::numeric_limits<T_REAL>::max());
	      lb_isNotValidChild1 = 
		gaintegerop::isNotValid
		(*aochrom_child1,
		 (T_CLUSTERIDX) 0,
		 aiinpcgaprobfixedk_inParamGA.getNumClusterK()
		 );
	      lb_isNotValidChild2 =
		gaintegerop::isNotValid
		(*aochrom_child2,
		 (T_CLUSTERIDX) 0,
		 aiinpcgaprobfixedk_inParamGA.getNumClusterK()
		 );
	      
	      if (lb_isNotValidChild1) ++ll_invalidOffspring;
	      if (lb_isNotValidChild2) ++ll_invalidOffspring;

	      ++lst_countlimitReachedCrossover;
	   
	    } /*END WHILE FINDING A VALID CHROMOSOME
	       */

	    if ( lb_isNotValidChild1 ) {
	      if ( uniformdis_real01(gmt19937_eng) < (T_REAL) 0.5 ) {
		
		*aochrom_child1 = *aichrom_parent1;
	      }
	      else {
		*aochrom_child1 = *aichrom_parent2;
	      }

#ifdef __VERBOSE_YES
	      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		std::cout << "limit is reached without finding a child1 [" 
			  << *aochrom_child1 << "]\n";
	      }
#endif //__VERBOSE_YES
      
	    }

	    if ( lb_isNotValidChild2 ) {
	      if ( uniformdis_real01(gmt19937_eng) < (T_REAL) 0.5 ) {
		*aochrom_child2 = *aichrom_parent1;
	      }
	      else {
		*aochrom_child2 = *aichrom_parent2;
	      }

#ifdef __VERBOSE_YES
	      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		std::cout << "limit is reached without finding a child1 [" 
			  << *aochrom_child1 << "]\n";
	      }
#endif //__VERBOSE_YES

	    }
		  
	 } /*END IF CROSSOVER
	    */
	 else {
	   *aochrom_child1 = *aichrom_parent1;
	   *aochrom_child2 = *aichrom_parent2;
	 }
       }
       );
	
      aoop_outParamEAC.sumTotalInvalidOffspring(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
      
  } //END CROSSOVER

    /*MUTATION------------------------------------------------------------------
     */
    { /*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(B) MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      long ll_invalidOffspring = 0;
      
      T_REAL lrt_pcj =
	  aiinpcgaprobfixedk_inParamGA.getProbMutation() +
	  ((T_REAL) llfh_listFuntionHist.getDomainUpperBound() /
	   (T_REAL) aiinpcgaprobfixedk_inParamGA.getNumMaxGenerations())
	* ( (1.0 / (T_REAL) gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize() )
	    - aiinpcgaprobfixedk_inParamGA.getProbMutation() );
 
      for (auto lchromfixleng_iter :lvectorchromfixleng_population) {
	
	if ( uniformdis_real01(gmt19937_eng) <  lrt_pcj ) 
	  {                    
	    std::pair<uintidx,T_CLUSTERIDX> lpair_previousGene =  
	    gaintegerop::mutation
	      (*lchromfixleng_iter,
	       [&](const T_CLUSTERIDX aicidx_previousGene) -> T_CLUSTERIDX
	       {	 
		 T_CLUSTERIDX  lit_rand;
		
		 do {
		   lit_rand = uniformdis_mmcidx0K(gmt19937_eng);
		 } while (lit_rand == aicidx_previousGene); 

		 return lit_rand;
	       }
	       );
	       
	    /*IS VERIFIED IF CORRECT THE MUTATION
	     */
	    T_CLUSTERIDX 
	      *lptmcidx_find = 
	      std::find
	      (lchromfixleng_iter->getString(),
	       lchromfixleng_iter->getString()+lchromfixleng_iter->getStringSize(),
	       lpair_previousGene.second
	       );
	    while ( lptmcidx_find == 
		    (lchromfixleng_iter->getString()+lchromfixleng_iter->getStringSize()) ) 
	      {
		//RECOVER THE ORIGINAL CHROMOSOME
		lchromfixleng_iter->setGene(lpair_previousGene.first,lpair_previousGene.second);
		//APPLY FOR NEW ACCOUNT MUTATION
	        lpair_previousGene =  
		gaintegerop::mutation
		  (*lchromfixleng_iter,
		   [&](const T_CLUSTERIDX aicidx_previousGene) -> T_CLUSTERIDX
		   {
		     T_CLUSTERIDX  lit_rand;
		     do {
		       lit_rand = uniformdis_mmcidx0K(gmt19937_eng);
		     } while (lit_rand == aicidx_previousGene); 

		     return lit_rand;
		   }
		   );
		
		++ll_invalidOffspring;
		lptmcidx_find = 
		  std::find
		  (lchromfixleng_iter->getString(),
		   lchromfixleng_iter->getString()+lchromfixleng_iter->getStringSize(),
		   lpair_previousGene.second
		   );
	      } /*while*/
	
	    lchromfixleng_iter->setObjetiveFunc
	      (std::numeric_limits<T_REAL>::max());
	       
	     } /*END IF MUTATION*/
      }
      
      aoop_outParamEAC.sumTotalInvalidOffspring(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END MUTATION*/

    
    /*COMPUTING CLUSTERING METRIC AND FITNESS COMPUTATION--------------------------
     */
    { /*BEGIN FITNESS COMPUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "COMPUTING CLUSTERING METRIC AND FITNESS COMPUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      long ll_invalidOffspring = 0;

      mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroids
	((uintidx) aiinpcgaprobfixedk_inParamGA.getNumClusterK(),
	 data::Instance<T_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<T_FEATURE_SUM>       
	lmatrixrowt_sumInstCluster
	((uintidx) aiinpcgaprobfixedk_inParamGA.getNumClusterK(), 
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
	
      std::vector<T_INSTANCES_CLUSTER_K> 
	lvectort_numInstClusterK
	((uintidx) aiinpcgaprobfixedk_inParamGA.getNumClusterK()
	 );
      
      for ( auto lchromfixleng_iter :lvectorchromfixleng_population ) {
	
	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (lchromfixleng_iter->getString(),
	   lui_numInstances,
	   aiinpcgaprobfixedk_inParamGA.getNumClusterK()
	   );

	T_CLUSTERIDX lmcidx_numClusterNull =
	  clusteringop::getCentroids
	  (lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   lpartition_clusters,
	   aiiterator_instfirst,
	   aiiterator_instlast
	   );
	   
	T_REAL lT_objetiveFunc;
	if ( lmcidx_numClusterNull == 0 ) { 
	  lT_objetiveFunc = 
	    um::SSE
	    (lmatrixrowt_centroids,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     lchromfixleng_iter->getString(),
	     aifunc2p_dist
	     );
	}
	else {
	  lT_objetiveFunc = std::numeric_limits<T_REAL>::max();
	}
	lchromfixleng_iter->setObjetiveFunc(lT_objetiveFunc);
	  
	lchromfixleng_iter->setFitness(1.0 / lchromfixleng_iter->getObjetiveFunc());
	lchromfixleng_iter->setValidString(lmcidx_numClusterNull == 0);
	
	if ( lchromfixleng_iter->getValidString() == false )
	  ++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

      }
      aoop_outParamEAC.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END FITNESS COMPUTATION*/


    /*(c) COMPARE THE WORST STRING IN Q1 IS FOUND TO WORSE THEN REPLACE---------
     */

    { /*BEGIN (c) COMPARE THE WORST STRING IN Q1 IS FOUND TO WORSE THEN REPLACE*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(C) COMPARE THE WORST STRING IN Q1 IS FOUND TO WORSE THEN REPLACE";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* lchromfixleng_iterMinFitness = 
       *(std::min_element
	(lvectorchromfixleng_population.begin(), 
	 lvectorchromfixleng_population.end(), 
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* x, 
	    const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* y
	    ) 
	 {  return x->getFitness() < y->getFitness(); }
	 )
	 );
      
      if ( lchromfixleng_iterMinFitness->getFitness() < lochromfixleng_best.getFitness() ) {
	
	*lchromfixleng_iterMinFitness = lochromfixleng_best;  
      
	} /*IF*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END (c) COMPARE THE WORST STRING IN Q1 IS FOUND TO WORSE THEN REPLACE*/


    /*(d) FIND THE BEST STRING IN Q1 (SAY S2) AND REPLACE S0 BY S2--------------
     */
    
    { /*BEGIN (d) FIND THE BEST STRING IN Q1 (SAY S2) AND REPLACE S0 BY S2*/
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(D) FIND THE BEST STRING IN Q1 (SAY S2) AND REPLACE S0 BY S2";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* lchromfixleng_max = 
	*(std::max_element
	  (lvectorchromfixleng_population.begin(), 
	   lvectorchromfixleng_population.end(), 
	   [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* x, 
	      const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* y
	      ) 
	{  return x->getFitness() < y->getFitness(); }
	   ));

      if ( lochromfixleng_best.getFitness() < lchromfixleng_max->getFitness() ) {
	
	lochromfixleng_best  = *lchromfixleng_max;
	
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	aoop_outParamEAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamEAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));
      }

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END (d) FIND THE BEST STRING IN Q1 (SAY S2) AND REPLACE S0 BY S2*/

   
    /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT  
    if ( aiinpcgaprobfixedk_inParamGA.getWithPlotStatObjetiveFunc() ) {  
      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    } 

#endif /*__WITHOUT_PLOT_STAT */

#ifdef __VERBOSE_YES
    
    /*ID PROC
     */
    ++geverboseui_idproc;
      
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "END ITERATION STEPS (B), (C), AND (D): "   << llfh_listFuntionHist.getDomainUpperBound()
		<< "\tobjetivoFunc = " << lochromfixleng_best.getObjetiveFunc() 
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END While*/

  
  {/*BEGIN FREE MEMORY OF POPULATION*/ 
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "DELETEPOPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    for (uintidx lui_i = 0; lui_i < lvectorchromfixleng_population.size(); ++lui_i) {
      delete lvectorchromfixleng_population[lui_i];
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }/*END FREE MEMORY OF POPULATION*/
    
  
  {/*BEGIN FREE MEMORY OF STRINGPOOL*/ 
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "DELETESTRINGPOOL";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    for (uintidx lui_i = 0; lui_i < lvectorchromfixleng_matingPool.size(); ++lui_i) {
      delete lvectorchromfixleng_matingPool[lui_i];
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END FREE MEMORY OF STRINGPOOL*/ 

  std::unordered_set<T_CLUSTERIDX> lounorderedset_numClusterK;
  lounorderedset_numClusterK.reserve(aiinpcgaprobfixedk_inParamGA.getNumClusterK());
  
  for ( auto liter_gene =  lochromfixleng_best.begin();
	liter_gene !=lochromfixleng_best.end();
	++liter_gene)
    {
      lounorderedset_numClusterK.insert(*liter_gene);
    }
  
  runtime::stop(let_executionTime);
  
  aoop_outParamEAC.setNumClusterK
    ((T_CLUSTERIDX)lounorderedset_numClusterK.size() );
  aoop_outParamEAC.setEndingCondition
    ( aiinpcgaprobfixedk_inParamGA.getNumClusterK() == aoop_outParamEAC.getNumClusterK() );
  aoop_outParamEAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamEAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
        
  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinpcgaprobfixedk_inParamGA.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinpcgaprobfixedk_inParamGA,
       aoop_outParamEAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA 
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochromfixleng_best.print();
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  return lochromfixleng_best;

} /*END  gaclustering_fklabel */

} /*END eac */

#endif /*__GA_MURTHY_CHOWDHURY_1996_HPP__*/
 
