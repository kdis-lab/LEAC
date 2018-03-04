/*! \file cga_vklabel.hpp
 *
 * \brief CGA \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003 \cite Hruschka:etal:GAclusteringLabelKVar:CGAII:2004 \cite Hruschka:etal:GAClusteringLabelKVar:EAC:2006
 * \details This file is part of the LEAC.\n\n
 * Implementation of the CGA algorithm based on the paper:\n  
 * E. R. Hruschka and N. F. F. Ebecken. A genetic algorithm for\n 
 * cluster analysis. Intell. Data Anal., 7(1):15–25, January 2003.\n 
 * <a href="http://dl.acm.org/citation.cfm?id=1293920.1293922">http://dl.acm.org/citation.cfm?id=1293920.1293922</a>.\n
 * \n
 * Library Evolutionary Algorithms for Clustering (LEAC) is a library\n
 * for the implementation of evolutionary and genetic algorithms\n
 * focused on the partition type clustering problem. Based on the\n
 * current standards of the <a href="http://en.cppreference.com">C++</a> language, as well as on Standard\n
 * Template Library <a href="http://en.cppreference.com/w/cpp/container">STL</a> 
 * and also  <a href="http://www.openblas.net/">OpenBLAS</a> to have a better performance.\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> \endlink license
 */

#ifndef __CGA_VKLABEL_HPP__
#define __CGA_VKLABEL_HPP__

#include <vector>
#include <algorithm>
#include <cmath>

#include <leac.hpp>

#include "inparam_probcprobm_rangek.hpp"
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
  
/*! \fn gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> cga_vklabel(inout::OutParamEAClustering <T_REAL, T_CLUSTERIDX> &aoop_outParamEAC, inout::InParamPcPmRk <T_CLUSTERIDX, T_REAL, T_FEATURE, T_FEATURE_SUM, T_INSTANCES_CLUSTER_K> &aiinpcgaprobfixedk_inParamGA, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
 \brief CGA \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003 
 \details Implementation of CGA algorithm based on \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003. 
 \returns A partition of a data set, encoded on a chromosome where each gene is the index of a cluster to which the instance belongs.
 \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
 \param aiinpcgaprobfixedk_inParamGA a inout::InParamPcPmRk parameters required by the algorithm
 \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
 \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
 \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,  //DATATYPE OF CHROMOSOME*/
	   typename T_REAL,        //OBJETIVE FUNCTION, FITNESS AND CLUSTERING_METRIC 
	   typename T_FEATURE,     //T_CENTROIDS
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename INPUT_ITERATOR
	   >
gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> 
cga_vklabel
(inout::OutParamEAClustering
 <T_REAL,
 T_CLUSTERIDX>                      &aoop_outParamEAC,
 inout::InParamPcPmRk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>             &aiinpcgaprobfixedk_inParamGA,
 const INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  if ( aiinpcgaprobfixedk_inParamGA.getNumClusterKMaximum() == 
       INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED )
    aiinpcgaprobfixedk_inParamGA.setNumClusterKMaximum
      ((T_CLUSTERIDX) lui_numInstances/2 );
     
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::setStringSize
    ((uintidx) lui_numInstances + 1);

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> lochromfixleng_best;
  lochromfixleng_best.setObjetiveFunc(-std::numeric_limits<T_REAL>::max());
  
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* >  
    lvectorchromfixleng_population;

  /*SPACE FOR STORE MATINGPOOL*/
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* >  
    lvectorchromfixleng_matingPool;

 
#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "cga_vklabel"; 
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
	      << "\t input aiiterator_instlast[" <<  *aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\nGA parameters: "
	      << "\tPopulation size = " << aiinpcgaprobfixedk_inParamGA.getSizePopulation()
	      << "\tProbCrossover = "   << aiinpcgaprobfixedk_inParamGA.getProbCrossover() 
	      << "\tProbMutation  = "   << aiinpcgaprobfixedk_inParamGA.getProbMutation()
	      << "\tKinitial = " << aiinpcgaprobfixedk_inParamGA.getNumClusterKMaximum()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 


  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinpcgaprobfixedk_inParamGA.getNumMaxGenerations(), "Iterations", "Clustering metrics");

 
  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
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
  

  /*POPULATION CREATE------------------------------------------------------------
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
 
  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
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
 
  /*POPULATION INITIAL-----------------------------------------------------------
   */
  { /*BEGIN INITIALIZATION*/
   
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "1. INITIALIZE A POPULATION OF RANDOM GENOTYPES;";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/  

    std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_kMinMax
      ( aiinpcgaprobfixedk_inParamGA.getNumClusterKMinimum(),
	aiinpcgaprobfixedk_inParamGA.getNumClusterKMaximum()
      );

    
    for (auto lchromfixleng_iter :lvectorchromfixleng_population) {

      T_CLUSTERIDX lcidx_Kini = 
	uniformdis_kMinMax(gmt19937_eng);

      std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidxKini
	(0,lcidx_Kini-1);
      
      gagenericop::initializeGenes
	(lchromfixleng_iter->begin(),
	 lchromfixleng_iter->end()-1,
	 [&]() 
	 {
	   return uniformdis_cidxKini(gmt19937_eng);
	 }
	 );
      
      lchromfixleng_iter->setGene
	(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
	 lcidx_Kini
	 );
				   
      lchromfixleng_iter->setObjetiveFunc(-std::numeric_limits<T_REAL>::max());
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


  while ( ( llfh_listFuntionHist.getDomainUpperBound() <= 
	    aiinpcgaprobfixedk_inParamGA.getNumMaxGenerations() )
          && ( runtime::elapsedTime(let_executionTime)
	       < aiinpcgaprobfixedk_inParamGA.getMaxExecutiontime() )
	  ) {
   
    llfh_listFuntionHist.increaseDomainUpperBound();


 
  /*COMPUTING CLUSTERING METRIC AND FITNESS COMPUTATION------------------------
   */

  { /*BEGIN FITNESS COMPUTATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "2. EVALUATE EACH GENOTYPE ACCORDING TO ITS SILHOUETTE;";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

  
    for (auto lchromfixleng_iter :lvectorchromfixleng_population) {

      if ( lchromfixleng_iter->getObjetiveFunc() == -std::numeric_limits<T_REAL>::max() ) {

	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (lchromfixleng_iter->getString(),
	   (uintidx) lui_numInstances,
	   lchromfixleng_iter->getGene
	   (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
	   );

	ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>
	  &&lpartlink_memberShip =
	  ds::getPartitionLinkedNumInst
	  (aiiterator_instfirst,
	   aiiterator_instlast,
	   lpartition_clusters,
	   [&](data::Instance<T_FEATURE>* liter_inst)
	   {
	     return T_INSTANCES_CLUSTER_K(1);
	   }
	   );

	auto li_clusterNull =
	  std::count_if
	  (lpartlink_memberShip.getVectorNumInstClusterK().begin(),
	   lpartlink_memberShip.getVectorNumInstClusterK().end(),
	   [](T_INSTANCES_CLUSTER_K liter_numClusterK)
	   {return liter_numClusterK == 0;}
	   );

	if ( li_clusterNull != 0 ) {
	  lchromfixleng_iter->setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
	  lchromfixleng_iter->setValidString(false);
	  aoop_outParamEAC.incTotalInvalidOffspring();
	}
	else {
	
	  lchromfixleng_iter->setObjetiveFunc
	    (um::silhouette
	     (aiiterator_instfirst,
	      lpartlink_memberShip,
	      lpartlink_memberShip.getVectorNumInstClusterK(),
	      aifunc2p_dist
	      )
	     );
	  lchromfixleng_iter->setValidString(true);
	}   
      }
      
#ifndef __WITHOUT_PLOT_STAT
      lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

    }
      
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END FITNESS COMPUTATION*/
     
  /*(a) COPY THE BEST STRING TO S0-----------------------------------------------
   */

  { /*BEGIN (a) COPY THE BEST STRING TO S0*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "COPY THE BEST STRING";
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
	       {  return x->getObjetiveFunc() < y->getObjetiveFunc(); }
	   )
	  );
    
      if ( lochromfixleng_best.getObjetiveFunc() <  lchromfixleng_max->getObjetiveFunc() ) {
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


  /*INITIAL MEASUREMENT: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM-------
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

  
    /*SELECTION------------------------------------------------------------------
     */

    { //BEGIN SELECTION

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "3. APPLY A LINEAR NORMALIZATION (RANKING);";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

	prob::linearNormalization
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* liter_iChrom) -> T_REAL
	 {
	   return liter_iChrom->getObjetiveFunc();
	 },
	 [](gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* liter_iChrom,
	    T_REAL airt_funcFitnessLineal
	    )
	 {
	   liter_iChrom->setFitness(airt_funcFitnessLineal);
	 },
	 T_REAL(1.0)
	 );

      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* liter_iChrom) -> T_REAL
	 {
	   return liter_iChrom->getFitness();

	 }
	 );

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "4. SELECT GENOTYPES BY PROPORTIONAL SELECTION;";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
	      
      /*COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL--------------------------
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


    /*CROSSOVER------------------------------------------------------------------
     */
    { /*BEGIN CROSSOVER*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "5.APPLY CROSSOVER";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      //long ll_invalidOffspring = 0;
      auto ichrom_newoffspring = lvectorchromfixleng_population.begin();
      auto jchrom_matingPool = lvectorchromfixleng_matingPool.begin();

      if ( ( lvectorchromfixleng_population.size() % 2 ) != 0 ) {
	*(*ichrom_newoffspring) = *(*jchrom_matingPool);
	++ichrom_newoffspring;
	++jchrom_matingPool;
      }
      while ( (ichrom_newoffspring != lvectorchromfixleng_population.end() )
	      && (jchrom_matingPool != lvectorchromfixleng_matingPool.end()) )
	{
	  auto lchrom_child1  = ichrom_newoffspring;
	  ++ichrom_newoffspring;
	  auto lchrom_child2  = ichrom_newoffspring;
	  ++ichrom_newoffspring;

	  auto lchrom_parent1  = jchrom_matingPool;
	  ++jchrom_matingPool;
	  auto lchrom_parent2  = jchrom_matingPool;
	  ++jchrom_matingPool;

	  if ( uniformdis_real01(gmt19937_eng) 
	       < aiinpcgaprobfixedk_inParamGA.getProbCrossover() ) {

	    mat::MatrixRow<T_FEATURE>          lmatrixrowt_centroids;
	    mat::MatrixRow<T_FEATURE_SUM>    lmatrixrowt_sumInstCluster;
	    std::vector<T_INSTANCES_CLUSTER_K> lvectort_numInstClusterK;
	      
	    gaclusteringop::crossoverCGA
	      (*(*lchrom_child1),
	       *(*lchrom_parent1),
	       *(*lchrom_parent2),
	       lmatrixrowt_centroids,
	       lmatrixrowt_sumInstCluster,
	       lvectort_numInstClusterK,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );

	    gaclusteringop::crossoverCGA
	      (*(*lchrom_child2),
	       *(*lchrom_parent2),
	       *(*lchrom_parent1),
	       lmatrixrowt_centroids,
	       lmatrixrowt_sumInstCluster,
	       lvectort_numInstClusterK,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );
	    
	    (*lchrom_child1)->setObjetiveFunc
	      (-std::numeric_limits<T_REAL>::max());
	    (*lchrom_child2)->setObjetiveFunc
	      (-std::numeric_limits<T_REAL>::max());
	   
	 	  
	  } //if  Crossover
	  else {
	    *(*lchrom_child1) = *(*lchrom_parent1);
	    *(*lchrom_child2) = *(*lchrom_parent2);
	   
	  }
	} /*END FOR*/
	
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
    } /*END CROSSOVER*/

    /*MUTATION-------------------------------------------------------------------
     */
    { /*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "5.APPLY MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

       
      for (auto lchromfixleng_iter :lvectorchromfixleng_population) {
	
	if ( uniformdis_real01(gmt19937_eng) <
	     aiinpcgaprobfixedk_inParamGA.getProbMutation()  ) 
	  {

	    mat::MatrixRow<T_FEATURE>          lmatrixrowt_centroids;
	    mat::MatrixRow<T_FEATURE_SUM>      lmatrixrowt_sumInstCluster;
	    std::vector<T_INSTANCES_CLUSTER_K> lvectort_numInstClusterK;
	    
	    gaclusteringop::MO1
	      (*lchromfixleng_iter,
	       lmatrixrowt_centroids,
	       lmatrixrowt_sumInstCluster,
	       lvectort_numInstClusterK,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );
	   
	    lchromfixleng_iter->setObjetiveFunc
	      (-std::numeric_limits<T_REAL>::max());
	    
	  } /*IF FOR MUTATION*/

	if ( uniformdis_real01(gmt19937_eng) <
	     aiinpcgaprobfixedk_inParamGA.getProbMutation()  ) 
	  {
	    gaclusteringop::MO2
	      (*lchromfixleng_iter,
	       aiiterator_instfirst,
	       aifunc2p_dist
	       );
	   
	    lchromfixleng_iter->setObjetiveFunc
	      (-std::numeric_limits<T_REAL>::max());
	    
	  } /*IF FOR MUTATION*/
	   	       
      } /*END FOR MUTATION*/
      
     
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END MUTATION*/

#ifdef __VERBOSE_YES
    
    /*ID PROC
     */
    ++geverboseui_idproc;
      
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "END ITERATION: "   << llfh_listFuntionHist.getDomainUpperBound()
		<< "\tobjetivoFunc = " << lochromfixleng_best.getObjetiveFunc() 
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END While*/

  /*FREE MEMORY*/
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
    
  runtime::stop(let_executionTime);
  aoop_outParamEAC.setNumClusterK
    (lochromfixleng_best.getGene
     (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1));
  aoop_outParamEAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamEAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
        
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

} /*END cga_vklabel */

} /*END eac */  

#endif /*__CGA_VKLABEL_HPP__*/
 
