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

#ifndef __CGA_VKLABEL_HPP__
#define __CGA_VKLABEL_HPP__

#include <vector>
#include <algorithm>
#include <cmath>

#include <leac.hpp>

#include "inparam_withoutpcpmvk.hpp"
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
  
/*! \fn gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> cga_vklabel(inout::OutParamGAC <T_REAL, T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamWithoutPcPmVk <T_CLUSTERIDX, T_REAL, T_FEATURE, T_FEATURE_SUM, T_INSTANCES_CLUSTER_K> &aiinp_inParamWithoutPcPmVk, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief CGA \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003 
  \details Implementation of CGA algorithm based on \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003. 
  \returns A partition of a data set, encoded on a chromosome where each gene is the index of a cluster to which the instance belongs.
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinp_inParamWithoutPcPmVk a inout::InParamWithoutPcPmVk parameters required by the algorithm
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
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                      &aoop_outParamGAC,
 inout::InParamWithoutPcPmVk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>             &aiinp_inParamWithoutPcPmVk,
 const INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  if ( aiinp_inParamWithoutPcPmVk.getNumClusterKMaximum() == 
       INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED )
    aiinp_inParamWithoutPcPmVk.setNumClusterKMaximum
      ((T_CLUSTERIDX) lui_numInstances/2 );
     
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::setStringSize
    ((uintidx) lui_numInstances + 1);
  
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> lochromfixleng_best;
  lochromfixleng_best.setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
  
  /*POPULATION CREATE------------------------------------------------------------
   */
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> >  
    lvectorchromfixleng_population(aiinp_inParamWithoutPcPmVk.getSizePopulation());

  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
   */
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> >  
    lvectorchromfixleng_matingPool(aiinp_inParamWithoutPcPmVk.getSizePopulation());
 
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
	      << "\t output inout::OutParamGAC&: aoop_outParamGAC[" 
	      << &aoop_outParamGAC << "]\n"
	      << "\t input  InParamClusteringGaProbFk&: aiinp_inParamWithoutPcPmVk[" 
	      << &aiinp_inParamWithoutPcPmVk << "]\n"
              << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  *aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\nGA parameters: "
	      << "\tPopulation size = " << aiinp_inParamWithoutPcPmVk.getSizePopulation()
	      << "\tKinitial = " << aiinp_inParamWithoutPcPmVk.getNumClusterKMaximum()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 


  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamWithoutPcPmVk.getNumMaxGenerations(), "Iterations", "Clustering metrics");

 
  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT
  
  std::ofstream lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL> lvectorT_statfuncObjetiveFunc;
  
  if ( aiinp_inParamWithoutPcPmVk.getWithPlotStatObjetiveFunc() ) {  

    lvectorT_statfuncObjetiveFunc.reserve
      (aiinp_inParamWithoutPcPmVk.getSizePopulation());  
    //DEFINE FUNCTION
    lofh_SSE  = new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinp_inParamWithoutPcPmVk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat<T_REAL>
	( (char) li_i,
	  aiinp_inParamWithoutPcPmVk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION 
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamWithoutPcPmVk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamWithoutPcPmVk.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open   
      (aoop_outParamGAC.getFileNameOutPlotStatObjetiveFunc().c_str(), 
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
  aoop_outParamGAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/
  
  runtime::ExecutionTime let_executionTime = runtime::start();

  /*calculate matrix dissimilarity
   */
  mat::MatrixTriang<T_REAL>&&
    lmatrixtriagT_dissimilarity = 
    dist::getMatrixDissimilarity
    (aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  { /*BEGIN INITIALIZE A POPULATION*/
   
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
      ( aiinp_inParamWithoutPcPmVk.getNumClusterKMinimum(),
	aiinp_inParamWithoutPcPmVk.getNumClusterKMaximum()
	);

    for (auto& lchromfixleng_iter :lvectorchromfixleng_population) {

      T_CLUSTERIDX lcidx_Kini = 
	uniformdis_kMinMax(gmt19937_eng);

      std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidxKini
	(0,lcidx_Kini-1);
      
      gagenericop::initializeGenes
	(lchromfixleng_iter.begin(),
	 lchromfixleng_iter.end()-1,
	 [&]() 
	 {
	   return uniformdis_cidxKini(gmt19937_eng);
	 }
	 );
      
      lchromfixleng_iter.setGene
	(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
	 lcidx_Kini
	 );
				   
      lchromfixleng_iter.setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL)); // -std::numeric_limits<T_REAL>::max());
      lchromfixleng_iter.setFitness(0.0);
      
    }
       
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END INITIALIZE A POPULATION**/


  while ( ( llfh_listFuntionHist.getDomainUpperBound() <= 
	    aiinp_inParamWithoutPcPmVk.getNumMaxGenerations() )
          && ( runtime::elapsedTime(let_executionTime)
	       < aiinp_inParamWithoutPcPmVk.getMaxExecutiontime() )
	  ) {
   
    llfh_listFuntionHist.increaseDomainUpperBound();


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

  
      for (auto& lchromfixleng_iter :lvectorchromfixleng_population) {


	if ( (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1) >= 2 ) {

	  partition::PartitionLabel
	    <T_CLUSTERIDX>
	    lpartition_clusters
	    (lchromfixleng_iter.getString(),
	     (uintidx) lui_numInstances,
	     lchromfixleng_iter.getGene
	     (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
	     );

	  ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>
	    &&lpartlinknuminst_memberShip =
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
	    (lpartlinknuminst_memberShip.getVectorNumInstClusterK().begin(),
	     lpartlinknuminst_memberShip.getVectorNumInstClusterK().end(),
	     [](T_INSTANCES_CLUSTER_K liter_numClusterK)
	     {return liter_numClusterK == 0;}
	     );

	  if ( li_clusterNull != 0 ) {
	    lchromfixleng_iter.setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
	    lchromfixleng_iter.setFitness(0.0);
	    lchromfixleng_iter.setValidString(false);
	    aoop_outParamGAC.incTotalInvalidOffspring();
	  }
	  else {
	
	    lchromfixleng_iter.setObjetiveFunc
	      (um::silhouette
	       (lmatrixtriagT_dissimilarity,
		lpartlinknuminst_memberShip
		)
	       );

	    lchromfixleng_iter.setFitness(lchromfixleng_iter.getObjetiveFunc()+1.0);
	    lchromfixleng_iter.setValidString(true);

	  }   
	} else {

	  lchromfixleng_iter.setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
	  lchromfixleng_iter.setFitness(0.0);
      
	}
      
#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter.getObjetiveFunc());
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
     
  
    { /*BEGIN GET BEST CHROMOSOME*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "GET BEST CHROMOSOME";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      auto lchromfixleng_max = 
	std::max_element
	(lvectorchromfixleng_population.begin(), 
	 lvectorchromfixleng_population.end(), 
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& x, 
	    const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& y
	    ) 
	 {  return x.getObjetiveFunc() < y.getObjetiveFunc(); }
	 );
    
      if ( lochromfixleng_best.getObjetiveFunc() <  lchromfixleng_max->getObjetiveFunc() ) {
	lochromfixleng_best  = *lchromfixleng_max;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
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
 
    } /*END GET BEST CHROMOSOME*/


    /*INITIAL MEASUREMENT: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM-------
     */
#ifndef __WITHOUT_PLOT_STAT  
    if ( aiinp_inParamWithoutPcPmVk.getWithPlotStatObjetiveFunc() ) {  
      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif  /*__WITHOUT_PLOT_STAT */

  
    { //BEGIN SELECTION    
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "3. APPLY A LINEAR NORMALIZATION (RANKING);";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif //__VERBOSE_YES
      
      prob::linearNormalization
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& liter_iChrom) -> T_REAL
	 {
	   return liter_iChrom.getObjetiveFunc();
	 },
	 [](gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& liter_iChrom,
	    T_REAL airt_funcFitnessLineal
	    )
	 {
	   liter_iChrom.setFitness(airt_funcFitnessLineal);
	 },
	 T_REAL(1.0)
	 );

      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& liter_iChrom) -> T_REAL
	 {
	   return liter_iChrom.getFitness();
	 }
	 );

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
      

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "4. SELECT GENOTYPES BY PROPORTIONAL SELECTION;";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      /*ELITIST STRATEGY THE BEST INDIVIDUAL IS COPIED INTO THE SUCCEEDING 
	GENERATION
       */
      auto lchromfixleng_iterMatilPool = lvectorchromfixleng_matingPool.begin();

      *lchromfixleng_iterMatilPool = lochromfixleng_best;
      lchromfixleng_iterMatilPool++;

      /*COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL--------------------------
       */
      for ( ;lchromfixleng_iterMatilPool != lvectorchromfixleng_matingPool.end();
	    lchromfixleng_iterMatilPool++) {
	uintidx luiidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );
	
	*lchromfixleng_iterMatilPool = lvectorchromfixleng_population.at(luiidx_chrom);
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

    const auto it_matilPoolG2 = std::next(lvectorchromfixleng_matingPool.begin(),lvectorchromfixleng_matingPool.size() / 2);
    const auto it_populationG2 = std::next(lvectorchromfixleng_population.begin(),lvectorchromfixleng_population.size() / 2);
    const auto it_matilPool3G4 = std::next(lvectorchromfixleng_matingPool.begin(),3 * lvectorchromfixleng_matingPool.size() / 4);

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
            
      gaiterator::crossoverFirstLast
	(lvectorchromfixleng_matingPool.begin(),
	 it_matilPoolG2,
	 lvectorchromfixleng_population.begin(),
	 it_populationG2,
	 [&](gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>&
	     aichrom_parent1,
	     gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>&
	     aichrom_parent2,
	     gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>&
	     aochrom_child1, 
	     gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>&
	     aochrom_child2
	     )
	 {
	   mat::MatrixRow<T_FEATURE>          lmatrixrowt_centroids;
	   mat::MatrixRow<T_FEATURE_SUM>    lmatrixrowt_sumInstCluster;
	   std::vector<T_INSTANCES_CLUSTER_K> lvectort_numInstClusterK;
	      
	   gaclusteringop::crossoverCGA
	     (aochrom_child1,
	      aichrom_parent1,
	      aichrom_parent2,
	      lmatrixrowt_centroids,
	      lmatrixrowt_sumInstCluster,
	      lvectort_numInstClusterK,
	      aiiterator_instfirst,
	      aiiterator_instlast,
	      aifunc2p_dist
	      );

	   gaclusteringop::crossoverCGA
	     (aochrom_child2,
	      aichrom_parent2,
	      aichrom_parent1,
	      lmatrixrowt_centroids,
	      lmatrixrowt_sumInstCluster,
	      lvectort_numInstClusterK,
	      aiiterator_instfirst,
	      aiiterator_instlast,
	      aifunc2p_dist
	      );
	  
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
      
    } /*END CROSSOVER*/

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

      auto lchromfixleng_iterPopulation = it_populationG2;
      auto lchromfixleng_iterMatilPool  = it_matilPoolG2;
 
      for ( ;lchromfixleng_iterMatilPool != it_matilPool3G4; 
	    lchromfixleng_iterMatilPool++, lchromfixleng_iterPopulation++) {       
      
	*lchromfixleng_iterPopulation = *lchromfixleng_iterMatilPool;
	
	mat::MatrixRow<T_FEATURE>          lmatrixrowt_centroids;
	mat::MatrixRow<T_FEATURE_SUM>      lmatrixrowt_sumInstCluster;
	std::vector<T_INSTANCES_CLUSTER_K> lvectort_numInstClusterK;
	    
	gaclusteringop::MO1
	  (*lchromfixleng_iterPopulation,
	   lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
      }

      for ( ;lchromfixleng_iterMatilPool != lvectorchromfixleng_matingPool.end();
	    lchromfixleng_iterMatilPool++, lchromfixleng_iterPopulation++) {

	*lchromfixleng_iterPopulation = *lchromfixleng_iterMatilPool;

	gaclusteringop::MO2
	  (*lchromfixleng_iterPopulation,
	   aiiterator_instfirst,
	   aifunc2p_dist
	   );
    
      } /*FOR*/
	   	        
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
    
  runtime::stop(let_executionTime);
  aoop_outParamGAC.setNumClusterK
    (lochromfixleng_best.getGene
     (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1));
  aoop_outParamGAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamGAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
        
#ifndef __WITHOUT_PLOT_STAT
  if ( aiinp_inParamWithoutPcPmVk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamWithoutPcPmVk,
       aoop_outParamGAC
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
 
