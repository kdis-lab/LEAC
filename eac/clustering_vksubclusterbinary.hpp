/*! \file clustering_vksubclusterbinary.hpp
 *
 * \brief CLUSTERING \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the CLUSTERING algorithm based on the paper:\n 
 * Lin Yu Tseng and Shiueng Bien Yang. A genetic approach to the\n
 * automatic clustering problem. Pattern Recognition, 34(2):415–424,2001.\n
 * <a href="http://dx.doi.org/10.1016/S0031-3203(00)00005-4">doi:http://dx.doi.org/10.1016/S0031-3203(00)00005-4</a>\n
 * \n
 * Library Evolutionary Algorithms for Clustering (LEAC) is a library\n
 * for the implementation of evolutionary algorithms\n
 * focused on the partition type clustering problem. Based on the\n
 * current standards of the <a href="http://en.cppreference.com">C++</a> language, as well as on Standard\n
 * Template Library <a href="http://en.cppreference.com/w/cpp/container">STL</a> 
 * and also  <a href="http://www.openblas.net/">OpenBLAS</a> to have a better performance.\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CLUSTERING_VKSUBCLUSTERBINARY_HPP__
#define __CLUSTERING_VKSUBCLUSTERBINARY_HPP__

#include <sstream>
#include <vector>
#include <algorithm>
#include <tuple>

#include <leac.hpp>

#include "plot_runtime_function.hpp"

#include "inparam_subclusterbinaryvk.hpp"
#include "outparam_gac.hpp"


/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace eac {

/*! \fn gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> clustering_genetic(inout::OutParamGAC<T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, const mat::MatrixRow<T_FEATURE> &aimatrixrowt_Vi, const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstBi, const T_REAL aitr_w, const partition::PartitionDisjSets <T_CLUSTERIDX> &aimembclassdisjsets_Bi, inout::InParamSubClusterBinaryVk<T_REAL,T_BITSIZE,T_FEATURE,T_FEATURE_SUM, T_INSTANCES_CLUSTER_K> &aiinp_inParamSubClusterBinVk, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist, const COMMON_IDOMAIN aii_numrunalg, const runtime::ExecutionTime aiet_executionTime) 
 \brief CLUSTERING GENETIC \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
 \details
 \param aoop_outParamGAC a OutParamClusteringGA<T_REAL,T_CLUSTERIDX>
 \param aimatrixrowt_Vi
 \param aivectort_numInstBi
 \param aitr_w
 \param aimembclassdisjsets_Bi
 \param aiinp_inParamSubClusterBinVk a inparam::InParamGAClusteringPcPmFk
 \param aifunc2p_dist a dist::Dist<T_REAL,T_FEATURE>
 \param aii_numrunalg
 \param aiet_executionTime
*/  
template < typename T_REAL,
           typename T_BITSIZE,
           typename T_FEATURE,  
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX
	   >
gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
clustering_genetic
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                                 &aoop_outParamGAC,
 const mat::MatrixRow<T_FEATURE>               &aimatrixrowt_Vi,
 const std::vector<T_INSTANCES_CLUSTER_K>      &aivectort_numInstBi,
 const T_REAL                                  aitr_w,
#ifdef __VERBOSE_YES
 const partition::PartitionDisjSets
 <T_CLUSTERIDX>                                &aimembclassdisjsets_Bi,
#endif /*__VERBOSE_YES*/
 inout::InParamSubClusterBinaryVk
 <T_REAL,
 T_BITSIZE,
 T_CLUSTERIDX,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                        &aiinp_inParamSubClusterBinVk,
 const dist::Dist<T_REAL,T_FEATURE>            &aifunc2p_dist,
 const COMMON_IDOMAIN                          aii_numrunalg,
 const runtime::ExecutionTime                  aiet_executionTime
 )
{
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   */
 
  gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    lochrom_best( aimatrixrowt_Vi.getNumRows() );

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0,1);
  
#ifdef __VERBOSE_YES
  const char* lpc_labelAlgGA = "clustering_genetic";
  geverbosepc_labelstep = lpc_labelAlgGA;
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA 
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output inout::OutParamGAC&: aoop_outParamGAC[" 
	      << &aoop_outParamGAC << "]\n"
	      << "\t input  InParamClusteringGAProb&: aiinp_inParamSubClusterBinVk[" 
	      << &aiinp_inParamSubClusterBinVk << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t\tPopulation size = " 
	      << aiinp_inParamSubClusterBinVk.getSizePopulation()
	      << "\n\t\tProbCrossover = " 
	      << aiinp_inParamSubClusterBinVk.getProbCrossover() 
	      << "\n\t\tProbMutation  = " 
	      << aiinp_inParamSubClusterBinVk.getProbMutation()
	      << "\n\t\tw  = " 
	      << aitr_w
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 


  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/
  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamSubClusterBinVk.getNumMaxGenerations(), "Iterations", "Clustering metrics");

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream              lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL>  *lofh_VRC = NULL;
  runtime::RuntimeFunctionStat<T_REAL>   *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>        lvectorT_statfuncObjetiveFunc;
  
  if ( aiinp_inParamSubClusterBinVk.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      ( aiinp_inParamSubClusterBinVk.getSizePopulation());
    //DEFINE FUNCTION
    lofh_VRC  = new runtime::RuntimeFunctionValue<T_REAL>
      ("VRC", 
       aiinp_inParamSubClusterBinVk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_VRC);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamSubClusterBinVk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamSubClusterBinVk.getFileNamePlotStatObjetiveFunc(),
       aitr_w
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
  /*aoop_outParamGAC.setTotalInvalidOffspring(0);
    executiontime_initialize(let_executionTime);
    executiontime_start(let_executionTime);
  */

#ifdef _INITIALIZED_ONLY_ONCE
  
  /*The genetic algorithm consists of an initializa-
    tion step and the iterative generations with three phases
    in each generation. They are described in the following.
   
    Initialization step: 
    A population of N strings is randomly generated. 
    The length of each string is m, which is
    the number of the sets obtained in the first stage.
    N strings are generated in such a way that the number of
    1's in the strings uniformly distributes within [1, m]. Each
    string represents a subset of {B_1, B_2,...,B_m}. If B_i is in
    this subset, the ith position of the string will be 1; other-
    wise, it will be 0. Each B_i in the subset is used as a seed to
    generate a cluster.
  */  
  
  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> >  lvectorchrom_population;
  
  {/*BEGIN INITIALIZE POPULATION*/

     
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "POPULATION INITIALIZATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    std::uniform_int_distribution<uintidx> uniformdis_ui2Vi
    (2,(uintidx) (aimatrixrowt_Vi.getNumRows()-1) );
 
    /*CREATE SPACE FOR STORE POPULATION-----------------------------------------
     */
    lvectorchrom_population.reserve
      ( aiinp_inParamSubClusterBinVk.getSizePopulation() );

    for (uintidx lui_i = 0; 
	 lui_i < aiinp_inParamSubClusterBinVk.getSizePopulation(); 
	 lui_i++) 
      {
	lvectorchrom_population.push_back
	  (gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>(aimatrixrowt_Vi.getNumRows()));
  
      }
    
    for (auto&& lchrom_iter: lvectorchrom_population ) {

      uintidx luintidx_krand = uniformdis_ui2Vi(gmt19937_eng);

      std::unordered_set<uintidx>&& lunorderedset_idxRandVi =
	prob::getWithoutRepeatsSet
	(luintidx_krand
	 ,[&]() -> uintidx
	 {
	   return uniformdis_ui0N(gmt19937_eng);
	 }
	 );

      lchrom_iter.initialize();

      std::for_each
	(lunorderedset_idxRandVi.begin(),
	 lunorderedset_idxRandVi.end(),
	 [&] (const uintidx &luiidx_vi)
	 {
	   lchrom_iter.setBit(luiidx_vi);
	 }
	 );	  
    }
   
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
        
  } /*END INITIALIZE POPULATION*/

#endif /*_INITIALIZED_ONLY_ONCE*/


#ifndef _INITIALIZED_ONLY_ONCE
  
  /*The genetic algorithm consists of an initializa-
    tion step and the iterative generations with three phases
    in each generation. They are described in the following.
   
    Initialization step: 
    A population of N strings is randomly generated. 
    The length of each string is m, which is
    the number of the sets obtained in the first stage.
    N strings are generated in such a way that the number of
    1's in the strings uniformly distributes within [1, m]. Each
    string represents a subset of {B_1, B_2,...,B_m}. If B_i is in
    this subset, the ith position of the string will be 1; other-
    wise, it will be 0. Each B_i in the subset is used as a seed to
    generate a cluster.
  */  
  
  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> >  lvectorchrom_population;
  
  {/*BEGIN INITIALIZE POPULATION*/

     
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "POPULATION INITIALIZATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
 
    /*CREATE SPACE FOR STORE POPULATION-----------------------------------------
     */
    lvectorchrom_population.reserve
      ( aiinp_inParamSubClusterBinVk.getSizePopulation() );

    for (uintidx lui_i = 0; 
	 lui_i < aiinp_inParamSubClusterBinVk.getSizePopulation(); 
	 lui_i++) 
      {
	lvectorchrom_population.push_back
	  (gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>(aimatrixrowt_Vi.getNumRows()));  
      }

    std::uniform_int_distribution<int>     uniformdis_01(0,1);

    for (auto&& lchrom_iter: lvectorchrom_population ) {
      gabinaryop::initializeGenes
	(lchrom_iter,
	 [&]() 
	 {
	   return uniformdis_01(gmt19937_eng);
	 }
	 );
    }
   
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
        
  } /*END INITIALIZE POPULATION*/

#endif /*_INITIALIZED_ONLY_ONCE*/
  
  while ( 1 ) {

    /*Reproduction phase: Let C be one of the clusters gener-
      ated by string R. We define Dintra to represent the intra-
    */
    {/*BEGIN REPRODUCTION PHASE*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "REPRODUCTION PHASE";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": IN" << '(' << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      //uintidx luintidx_invalidOffspring = 0;

      for (auto&& lchrom_iter: lvectorchrom_population ) {

	if ( lchrom_iter.getNumBitOn() > 1 ) {
	   
	  mat::MatrixRow<T_FEATURE>              lmatrixrowt_S;
	  partition::PartitionDisjSets
	    <T_CLUSTERIDX>                       lmembclassdisjsets_BkinCi;
	  std::vector<T_INSTANCES_CLUSTER_K>     lvectort_numInstCi;
	  std::tie(lmatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
	    clusteringop::getClusters
	    (lchrom_iter,
	     aimatrixrowt_Vi,
	     aivectort_numInstBi,
	     aifunc2p_dist,
	     (T_CLUSTERIDX) 0
	     );

#ifdef __VERBOSE_YES
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    
	    std::vector<T_CLUSTERIDX>
	      lvectormmidx_memberShip =
	      clusteringop::bkinCiToMemberShip
	      <T_CLUSTERIDX>
	      (aimembclassdisjsets_Bi,
	       lmembclassdisjsets_BkinCi
	       );
	    std::ostringstream lostrstream_labelMember;
	    lostrstream_labelMember << "<MEMBERCLUSTERCLUSTERCj";
	    inout::containerprint
	      (lvectormmidx_memberShip.begin(),
	       lvectormmidx_memberShip.end(),
	       std::cout,
	       lostrstream_labelMember.str().c_str()
	       ,','
	       );
	    std::cout << std::endl;
	     
	    std::ostringstream lostrstream_labelNumInstCi;
	    lostrstream_labelNumInstCi
	      << "<NUMINSTANCES CLUSTER Cj"; 
	    inout::containerprint
	      (lvectort_numInstCi.begin(),
	       lvectort_numInstCi.end(),
	       std::cout,
	       lostrstream_labelNumInstCi.str().c_str()
	       ,','
	       );
	    std::cout << std::endl;
	     
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/   
	   
	  /*Dintra(Ci)*/
	     
	  T_REAL lt_sumDintra =  
	    um::Dintra
	    (lmatrixrowt_S,
	     aimatrixrowt_Vi,
	     aivectort_numInstBi,
	     lmembclassdisjsets_BkinCi,
	     aifunc2p_dist
	     );
	     
	  T_REAL lt_sumDinter =
	    um::Dinter
	    (lmatrixrowt_S,
	     aimatrixrowt_Vi,
	     aivectort_numInstBi,
	     lmembclassdisjsets_BkinCi,
	     aifunc2p_dist
	     );    
	 
	  lchrom_iter.setObjetiveFunc(lt_sumDinter * aitr_w - lt_sumDintra); 
	  lchrom_iter.setFitness(lchrom_iter.getObjetiveFunc());
	  lchrom_iter.setValidString(true); 
	}
	else {
	  lchrom_iter.setObjetiveFunc(T_REAL(0)); 
	  lchrom_iter.setFitness(T_REAL(0)); 
	  lchrom_iter.setValidString(false); 
	  aoop_outParamGAC.incTotalInvalidOffspring();
	} 

#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(lchrom_iter.getObjetiveFunc());
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

      
    } /*END REPRODUCTION PHASE*/

    /*The best solution was retained
     */
    { /*BEGIN PRESERVING THE BEST STRING*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
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
	 [](const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>& x, 
	    const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>& y
	    ) 
	 {  return x.getFitness() < y.getFitness(); }
	 );
     
      if ( lochrom_best.getFitness() < (*lit_chromMax).getFitness() ) {
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION
	 */
	lochrom_best = *lit_chromMax;

	aoop_outParamGAC.setIterationGetsBest
	  (aii_numrunalg  * aiinp_inParamSubClusterBinVk.getNumMaxGenerations()
	   +  llfh_listFuntionHist.getDomainUpperBound()
	   );
	aoop_outParamGAC.setRunTimeGetsBest
	(runtime::elapsedTime(aiet_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  const char* lpc_labelStep = "CHROMOSOME ONE WAS FOUND IN THIS ITERATION: ";

	  mat::MatrixRow<T_FEATURE>                lmatrixrowt_S;
	  partition::PartitionDisjSets
	    <T_CLUSTERIDX>                         lmembclassdisjsets_BkinCi;

	  std::vector<T_INSTANCES_CLUSTER_K>       lvectort_numInstCi;
	  std::tie(lmatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
	    clusteringop::getClusters
	    (lochrom_best,
	     aimatrixrowt_Vi,
	     aivectort_numInstBi,
	     aifunc2p_dist,
	     (T_CLUSTERIDX) 0
	     );

	  std::ostringstream lostrstream_labelCentroids;
	  lostrstream_labelCentroids << "<CENTROIDS Cj:" << lpc_labelStep;
	  lmatrixrowt_S.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
	  std::cout << std::endl;
	  
	  std::ostringstream lostrstream_labelBkinCi;
	  lostrstream_labelBkinCi << "<MEMBER BkinCi: " << lpc_labelStep;
	  lmembclassdisjsets_BkinCi.print(std::cout,lostrstream_labelBkinCi.str().c_str(),',');
	  std::cout << std::endl;
	     
	  std::vector<T_CLUSTERIDX>
	    lvectormmidx_memberShip =
	    clusteringop::bkinCiToMemberShip
	    <T_CLUSTERIDX>
	    (aimembclassdisjsets_Bi,
	     lmembclassdisjsets_BkinCi
	     );
	  
	  std::ostringstream lostrstream_labelMember;
	  lostrstream_labelMember << "<MEMBERCLUSTER CLUSTER Cj: " << lpc_labelStep;
	  
	 inout::containerprint
	    (lvectormmidx_memberShip.begin(),
	     lvectormmidx_memberShip.end(),
	     std::cout,
	     lostrstream_labelMember.str().c_str(),
	     ','
	     );
	  
	  std::cout << std::endl;
	     
	  std::ostringstream lostrstream_labelNumInstCi;
	  lostrstream_labelNumInstCi << "<NUM INSTANCES CLUSTER Cj: " <<  lpc_labelStep;
	  
	  inout::containerprint
	    (lvectort_numInstCi.begin(),
	     lvectort_numInstCi.end(),
	     std::cout,
	     lostrstream_labelNumInstCi.str().c_str(),
	     ','
	     );
	  std::cout << std::endl; 	      
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
    if ( aiinp_inParamSubClusterBinVk.getWithPlotStatObjetiveFunc() ) {  
      lofh_VRC->setValue(lochrom_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

    /*100 generations were run
     */
      
#ifdef __VERBOSE_YES
    /*ID PROC
     */
    ++geverboseui_idproc;
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "TERMINATION CRITERION ATTAINED?: NumMaxGenerations =  " 
		<< llfh_listFuntionHist.getDomainUpperBound() 
		<< std::endl; 
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    if ( (llfh_listFuntionHist.getDomainUpperBound() >= aiinp_inParamSubClusterBinVk.getNumMaxGenerations() ) )
      break;
   
    /*Selection
     */

    /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
     */
    std::vector
      <gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> > 
      lvectorchrom_stringPool;

    lvectorchrom_stringPool.reserve
      (aiinp_inParamSubClusterBinVk.getSizePopulation());

    {/*BEGIN SELECTION*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "STEP SELECTION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      

      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchrom_population.begin(),lvectorchrom_population.end(),
	 [](const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>& liter_iChrom) -> T_REAL
	 {
	   return T_REAL(liter_iChrom.getFitness());
	 }
	 );
     
      
      for (uintidx lui_i = 0; 
	   lui_i < aiinp_inParamSubClusterBinVk.getSizePopulation(); 
	   lui_i++) 
	{      
	  uintidx lstidx_chrom = 
	    gaselect::getIdxRouletteWheel
	    (lvectorT_probDistRouletteWheel,
	     uintidx(0)
	     );
	  lvectorchrom_stringPool.push_back
	    (gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
	     (lvectorchrom_population.at(lstidx_chrom))
	     );
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


      /*Crossover phase: 
	If a pair of strings R and Q are
	chosen for applying the crossover operator, two random
	numbers p and q in [1, m] are generated to decide
	which pieces of the strings are to be interchanged.
	Suppose p(q, the bits from position p to position q of
	string R will be interchanged with those bits that are in
	the same position of string Q. For each chosen pair
	of strings, the crossover operator is done with prob-
	ability p.
      */

    { /*BEGIN CROSSOVER*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "STEP CROSSOVER";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      auto ichrom_population = lvectorchrom_population.begin();
      auto jchrom_matingPool = lvectorchrom_stringPool.begin();

      if ( ( lvectorchrom_population.size() % 2 ) != 0 ) {
	*ichrom_population = *jchrom_matingPool;
	++ichrom_population;
	++jchrom_matingPool;
      }
      while ( (ichrom_population != lvectorchrom_population.end() )
	      && (jchrom_matingPool != lvectorchrom_stringPool.end()) )
	{
	  auto lchrom_child1  = ichrom_population;
	  ++ichrom_population;
	  auto lchrom_child2  = ichrom_population;
	  ++ichrom_population;

	  auto lchrom_parent1  = jchrom_matingPool;
	  ++jchrom_matingPool;
	  auto lchrom_parent2  = jchrom_matingPool;
	  ++jchrom_matingPool;

	  if ( uniformdis_real01(gmt19937_eng) 
	       < aiinp_inParamSubClusterBinVk.getProbCrossover() ) {
	  
	    gabinaryop::onePointDistCrossover
	      (*lchrom_child1,
	       *lchrom_child2,
	       *lchrom_parent1,
	       *lchrom_parent2   
	       );

	    (*lchrom_child1).setFitness(-std::numeric_limits<T_REAL>::max());  
	    (*lchrom_child1).setObjetiveFunc(-std::numeric_limits<T_REAL>::max());
	    (*lchrom_child2).setFitness(-std::numeric_limits<T_REAL>::max());  
	    (*lchrom_child2).setObjetiveFunc(-std::numeric_limits<T_REAL>::max());
	  
	  } //if  Crossover
	  else {
	    *lchrom_child1 = *lchrom_parent1;
	    *lchrom_child2 = *lchrom_parent2;	
	  }
	  
	} /*while*/

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
      
    } /*END CROSSOVER*/


      /*Mutation phase: 
	In the mutation phase, bits of the
	strings in the population will be chosen with probability
	p . Each chosen bit will be changed from 0 to 1 or from
	1 to 0.
      */

    { /*BEGIN MUTATION*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "STEP MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for (auto &&liter_iChrom: lvectorchrom_population) {
	
	   if ( uniformdis_real01(gmt19937_eng) 
		< aiinp_inParamSubClusterBinVk.getProbMutation() ) 
	     { //IF BEGIN  MUTATION
	       gabinaryop::bitMutation(liter_iChrom);
	       liter_iChrom.setFitness
		 (-std::numeric_limits<T_REAL>::max());  
	       liter_iChrom.setObjetiveFunc
		 (-std::numeric_limits<T_REAL>::max());
		 	     
	     } //END BEGIN  MUTATION
	 }


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


    } /*END MUTATION */
	
      /*END GENETIC OPERATIONS*/
   
    llfh_listFuntionHist.increaseDomainUpperBound();

  } /*END EVOLUTION While*/ 

  aoop_outParamGAC.incNumGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamSubClusterBinVk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamSubClusterBinVk,
       aoop_outParamGAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/
   
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA
	      << ": OUT(" << geiinparam_verbose << ")\n";

    const char* lpc_labelAlgGAOut = "OUT";
    
    for (auto&& lchrom_iter: lvectorchrom_population ) {
      lchrom_iter.print(std::cout,lpc_labelAlgGAOut);
      std::cout << '\n';
    }
    
    std::ostringstream lostrstream_labelChrom;
    lostrstream_labelChrom << lpc_labelAlgGAOut
			   << ":w," << aitr_w ;
    lochrom_best.print(std::cout,lostrstream_labelChrom.str().c_str());
    std::cout << std::endl;
      
      mat::MatrixRow<T_FEATURE>          lmatrixrowt_S;
      partition::PartitionDisjSets
      	<T_CLUSTERIDX>                   lmembclassdisjsets_BkinCi;
      std::vector<T_INSTANCES_CLUSTER_K> lvectort_numInstCi;
      std::tie(lmatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
	clusteringop::getClusters
	(lochrom_best,
	 aimatrixrowt_Vi,
	 aivectort_numInstBi,
	 aifunc2p_dist,
	 (T_CLUSTERIDX) 0
	 );
      
      std::ostringstream lostrstream_labelCentroids;
      lostrstream_labelCentroids << "<CENTROIDS Sj:"
				 << lpc_labelAlgGAOut
				 << ":w," << aitr_w; 
      lmatrixrowt_S.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
      std::cout << std::endl;

      std::vector<T_CLUSTERIDX>
	lvectormmidx_memberShip =
	clusteringop::bkinCiToMemberShip
	<T_CLUSTERIDX>
	(aimembclassdisjsets_Bi,
	 lmembclassdisjsets_BkinCi
	 );

      std::ostringstream lostrstream_labelMember;
      lostrstream_labelMember << "<MEMBERCLUSTER CLUSTER Cj:"
			      << lpc_labelAlgGAOut
	                      << ":w," << aitr_w;      
      inout::containerprint
	(lvectormmidx_memberShip.begin(),
	 lvectormmidx_memberShip.end(),
	 std::cout,
	 lostrstream_labelMember.str().c_str(),
	 ','
	 );
      std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lochrom_best;  
 
} /* END clustering_genetic */


/*! \fn std::tuple<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>, mat::MatrixRow<T_FEATURE>, std::vector<T_INSTANCES_CLUSTER_K>, partition::PartitionDisjSets<T_CLUSTERIDX> > clustering_vksubclusterbinary(inout::OutParamGAC <T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamSubClusterBinaryVk<T_REAL,T_BITSIZE,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamSubClusterBinVk, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
 \brief CLUSTERING \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
 \details
 \details Implementation of the CBGA algorithm based on \cite Franti:etal:GAclustering:gafranti:1997.  
 \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid or a codebook.
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinParam_CBGA a inout::InParamGAClusteringPcPmFk parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param aoop_outParamGAC a OutParamClusteringGA<T_REAL,T_CLUSTERIDX>
  \param aiinpcgaprobfixedk_inParamGA a inparam::InParamGAClusteringPcPmFk<T_CLUSTERIDX,T_REAL>
  \param aivectorptinst_instances a std::vector<data::Instance<T_FEATURE>* >
  \param aifunc2p_dist a dist::Dist<T_REAL,T_FEATURE>
*/  
template < typename T_BITSIZE,
           typename T_FEATURE, 
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_REAL,
	   typename INPUT_ITERATOR
	   >
std::tuple
<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>,
 mat::MatrixRow<T_FEATURE>,
 std::vector<T_INSTANCES_CLUSTER_K>,
 partition::PartitionDisjSets<T_CLUSTERIDX>
 >
clustering_vksubclusterbinary
(inout::OutParamGAC
 <T_REAL,T_CLUSTERIDX>                          &aoop_outParamGAC,
 inout::InParamSubClusterBinaryVk
 <T_REAL,
 T_BITSIZE,
 T_CLUSTERIDX,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                         &aiinp_inParamSubClusterBinVk,
 const INPUT_ITERATOR                           aiiterator_instfirst,
 const INPUT_ITERATOR                           aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>                   &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  const char* lpc_labelAlgHeuristic = "clustering_vksubclusterbinary"; 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgHeuristic 
	      << ": IN(" << geiinparam_verbose << ")\n"
	      << "\t(output inout::OutParamGAC&: aoop_outParamGAC[" 
	      << &aoop_outParamGAC << "]\n"
	      << "\t input  InParamClusteringGAProb&: aiinp_inParamSubClusterBinVk[" 
	      << &aiinp_inParamSubClusterBinVk << "]\n"
              << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t\tPopulation size = " 
	      << aiinp_inParamSubClusterBinVk.getSizePopulation()
	      << "\n\t\tProbCrossover = " 
	      << aiinp_inParamSubClusterBinVk.getProbCrossover() 
	      << "\n\t\tProbMutation  = " 
	      << aiinp_inParamSubClusterBinVk.getProbMutation()
	      << "\n\t\tu  = " 
	      << aiinp_inParamSubClusterBinVk.getU()
	      << "\n\t\tw1  = " 
	      << aiinp_inParamSubClusterBinVk.getW1()
	      << "\n\t\tw2  = " 
	      << aiinp_inParamSubClusterBinVk.getW2()
	      << "\n\t\tlambda  = " 
	      << aiinp_inParamSubClusterBinVk.getLambda()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 
  
  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING
   */
  aoop_outParamGAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/
  runtime::ExecutionTime let_executionTime = runtime::start();
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "NEAREST-NEIGHBOR ALGORITHM";
  geverbosepc_labelstep = lpc_labelFunc;
#endif /*__VERBOSE_YES*/
      
  std::vector<T_CLUSTERIDX>  lovectormmidx_memberShip;
  mat::MatrixRow<T_FEATURE>  lomatrixrowt_S;

  gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> lochrom_best;
 
  partition::PartitionDisjSets
    <T_CLUSTERIDX>            lmembclassdisjsets_BkinCi;
  std::vector<T_INSTANCES_CLUSTER_K> lvectort_numInstCi;
  
  partition::PartitionDisjSets<T_CLUSTERIDX>
    lmembclassdisjsets_Bi =
    graph::nearestNeighbor
    (aiinp_inParamSubClusterBinVk.getU(),
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist,
     T_CLUSTERIDX(0)
     );

  mat::MatrixRow<T_FEATURE> //centroids
    lmatrixrowt_Vi
    ( uintidx(lmembclassdisjsets_Bi.getNumCluster()),
      data::Instance<T_FEATURE>::getNumDimensions()
      );
  std::vector<T_INSTANCES_CLUSTER_K>
    lvectort_numInstBi(lmatrixrowt_Vi.getNumRows(),T_INSTANCES_CLUSTER_K(0));
   
  if ( lmembclassdisjsets_Bi.getNumCluster() > 1) {
 
    std::vector<T_REAL> lvectorr_meanRadiusBk;
  
    {/*BEGIN NEAREST-NEIGHBOR ALGORITHM
       The first stage
     */
#ifdef __VERBOSE_YES
      // "NEAREST-NEIGHBOR ALGORITHM";
      geverbosepc_labelstep = lpc_labelFunc;
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      mat::MatrixRow<T_FEATURE_SUM>       
	lmatrixrowt_sumInstancesClusterK
	(lmatrixrowt_Vi.getNumRows(),
	 data::Instance<T_FEATURE>::getNumDimensions(),
	 T_FEATURE_SUM(0)
	 );

      clusteringop::getCentroids
	(lmatrixrowt_Vi,
	 lmatrixrowt_sumInstancesClusterK,
	 lvectort_numInstBi,
	 lmembclassdisjsets_Bi,
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );
       
      lvectorr_meanRadiusBk =
	um::avgRadiusClusterK
	(lmatrixrowt_Vi,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 lmembclassdisjsets_Bi,
	 aifunc2p_dist
	 );
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	lpc_labelFunc = geverbosepc_labelstep;
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ")\n";
	std::ostringstream lostrstream_labelVi;
	lostrstream_labelVi << "<CENTROIDS Vi:";
	lmatrixrowt_Vi.print(std::cout,lostrstream_labelVi.str().c_str(),',',';');
	std::cout << '\n';

	std::ostringstream lostrstream_labelBi;
	lostrstream_labelBi << "<MEMBERCLUSTER Bi:";
	lmembclassdisjsets_Bi.print(std::cout,lostrstream_labelBi.str().c_str(),',');
	std::cout << '\n';
	
	std::ostringstream lostrstream_labelNumInstBi;
	lostrstream_labelNumInstBi << "<NUM INSTANCES Bi:";
	
	inout::containerprint
	  (lvectort_numInstBi.begin(),
	   lvectort_numInstBi.end(),
	   std::cout,
	   lostrstream_labelNumInstBi.str().c_str(),
	   ','
	   );
	std::cout << std::endl;
	  
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END NEAREST-NEIGHBOR ALGORITHM*/
  
    /*4. The heuristic strategy to find a good clustering
     */ 
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "4. THE HEURISTIC STRATEGY TO FIND A GOOD CLUSTERING:  IN"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    /*Basically, the second stage is a genetic algorithm, which will
      merge some of these B_i's if they are close enough to one
      another. 
    */
  
    /*Step 1:
      Initially, let variables wS and wL indicate, respec-
      tively, the smallest value and the largest value
      within the given range, that is, ws = w1 and wL= w2. 
      Use CLUSTERING with the parameter wL 
      to cluster the data set. Use CLUSTERING
      with the parameter w to cluster the data set.
    */
      
#ifdef __VERBOSE_YES
    lpc_labelFunc =  "STEP 1: INITIALLY, LET VARIABLES Ws AND WL";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      geverbosepc_labelstep = lpc_labelFunc;
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< "\n\twS = " << aiinp_inParamSubClusterBinVk.getW1()
		<< "\n\twL = " << aiinp_inParamSubClusterBinVk.getW2()
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
      

    aoop_outParamGAC.setNumTotalGenerations(0); 

    COMMON_IDOMAIN lii_numrunalg = 0;
      
    T_REAL ltr_wS = aiinp_inParamSubClusterBinVk.getW1();
    T_REAL ltr_wL = aiinp_inParamSubClusterBinVk.getW2();

    /*Ws
     */
    gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>&& lchrom_wS = 
      clustering_genetic
      (aoop_outParamGAC,
       lmatrixrowt_Vi,
       lvectort_numInstBi,
       ltr_wS,     
#ifdef __VERBOSE_YES
       lmembclassdisjsets_Bi,
#endif /*__VERBOSE_YES*/
       aiinp_inParamSubClusterBinVk,
       aifunc2p_dist,
       lii_numrunalg++,
       let_executionTime
       );

   
    std::tie(lomatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
      clusteringop::getClusters
      (lchrom_wS,
       lmatrixrowt_Vi,
       lvectort_numInstBi,
       aifunc2p_dist,
       (T_CLUSTERIDX) 0
       );

    T_REAL ltr_D1Ws = 
      um::minDistCjCjp //D1(w)
      (lomatrixrowt_S,
       aifunc2p_dist
       );
  
    T_REAL ltr_D2Ws =
      um::D2
      (lomatrixrowt_S,
       lmembclassdisjsets_BkinCi,
       lmembclassdisjsets_Bi,
       lmatrixrowt_Vi,
       lvectorr_meanRadiusBk,
       aiiterator_instfirst,
       aiiterator_instlast,
       aifunc2p_dist
       );
   
    /*Wl
     */
    gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>&& lchrom_wL = 
      clustering_genetic
      (aoop_outParamGAC,
       lmatrixrowt_Vi,
       lvectort_numInstBi,
       ltr_wL,
#ifdef __VERBOSE_YES
       lmembclassdisjsets_Bi,
#endif /*__VERBOSE_YES*/
       aiinp_inParamSubClusterBinVk,
       aifunc2p_dist,
        lii_numrunalg++,
       let_executionTime
       );
    
    std::tie(lomatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
      clusteringop::getClusters
      (lchrom_wL,
       lmatrixrowt_Vi,
       lvectort_numInstBi,
       aifunc2p_dist,
       (T_CLUSTERIDX) 0
       );
    
    T_REAL ltr_D1Wl = 
      um::minDistCjCjp //D1(w)
      (lomatrixrowt_S,
       aifunc2p_dist
       );

    T_REAL ltr_D2Wl =
      um::D2
      (lomatrixrowt_S,
       lmembclassdisjsets_BkinCi,
       lmembclassdisjsets_Bi,
       lmatrixrowt_Vi,
       lvectorr_meanRadiusBk,
       aiiterator_instfirst,
       aiiterator_instlast,
       aifunc2p_dist
       ); 

#ifdef __VERBOSE_YES
    
    gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> lchrom_wLtmp(lchrom_wL);
    gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> lchrom_WStmp(lchrom_wS);
    
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      geverbosepc_labelstep = lpc_labelFunc;
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
      
    /*Step 2:
      Do while $wL - wS > lambda$  ( lambda is 0.125 in our experiments)
    */

#ifdef __VERBOSE_YES
    lpc_labelFunc =  "STEP 2: DO WHILE w_L - w_S > lambda";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      geverbosepc_labelstep =  lpc_labelFunc;
      std::cout << geverbosepc_labelstep 
		<< ": IN(" << geiinparam_verbose << ") "
		<< ltr_wL -ltr_wS << " > " << aiinp_inParamSubClusterBinVk.getLambda()  
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    T_REAL ltr_D2Wm;
    bool lb_saveFisrt = true;

    while ( (ltr_wL -ltr_wS) > aiinp_inParamSubClusterBinVk.getLambda() ) { 

      /*Wm
       */
      T_REAL ltr_wm = (ltr_wS + ltr_wL) / 2.0;

      gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> lchrom_Wm = 
	clustering_genetic
	(aoop_outParamGAC,
	 lmatrixrowt_Vi,
	 lvectort_numInstBi,
	 ltr_wm,
#ifdef __VERBOSE_YES
	 lmembclassdisjsets_Bi,
#endif /*__VERBOSE_YES*/
	 aiinp_inParamSubClusterBinVk,
	 aifunc2p_dist,
	 lii_numrunalg++,
	 let_executionTime
	 );

      std::tie(lomatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
	clusteringop::getClusters
	(lchrom_Wm,
	 lmatrixrowt_Vi,
	 lvectort_numInstBi,
	 aifunc2p_dist,
	 (T_CLUSTERIDX) 0
	 );

      T_REAL ltr_D1Wm = 
	um::minDistCjCjp //D1(w)
	(lomatrixrowt_S,
	 aifunc2p_dist
	 );

      if (lb_saveFisrt) {
	ltr_D2Wm =
	  um::D2
	  (lomatrixrowt_S,
	   lmembclassdisjsets_BkinCi,
	   lmembclassdisjsets_Bi,
	   lmatrixrowt_Vi,
	   lvectorr_meanRadiusBk,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
	lb_saveFisrt = false;
      }

      T_REAL ltr_ratiosD1WmWs = ltr_D1Wm /ltr_D1Ws;
      T_REAL ltr_ratiosD1WlWm = ltr_D1Wl /ltr_D1Wm;

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "<HEURISTIC STRATEGY D1:>"
		  << "\nWs\t"         << ltr_wS
		  << "\tk\t"          << lchrom_wS.getNumBitOn()
		  << "\tD1Ws\t"       << ltr_D1Ws
		  << "\tfitness wS\t" << lchrom_wS.getFitness()
		  << "\nWl\t"         << ltr_wL
		  << "\tk\t"          << lchrom_wL.getNumBitOn()
		  << "\tD1Wl\t"       << ltr_D1Wl
		  << "\tfitness wL\t" << lchrom_wL.getFitness()
		  << "\nWm\t"         << ltr_wm
		  << "\tk\t"          << lchrom_Wm.getNumBitOn()
		  << "\tD1Wm\t"       << ltr_D1Wm
		  << "\tfitness wm\t" << lchrom_Wm.getFitness()
		  << "\n\tD1(Wm)/D1(wS)\t" << ltr_ratiosD1WmWs
		  << "\n\tD1(Wl)/D1(Wm)\t" << ltr_ratiosD1WlWm
		  << std::endl;
      }
      //--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

      if ( ltr_ratiosD1WmWs > ltr_ratiosD1WlWm ) {
	ltr_wL   = ltr_wm;
	ltr_D1Wl = ltr_D1Wm;
	lchrom_wL = lchrom_Wm;
      }
      else {
	ltr_wS   = ltr_wm;
	ltr_D1Ws = ltr_D1Wm;
	lchrom_wS = lchrom_Wm;
      }

#ifdef __VERBOSE_YES
      // ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "\t[" << ltr_wS << "," << ltr_wL << ']'
		  << std::endl;
      }
      //--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
    } /* end while */


    T_REAL ltr_wp = ltr_wL;
    gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> lchrom_wp(lchrom_wL); 
    
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	geverbosepc_labelstep = lpc_labelFunc;
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << "\n\twS = " << ltr_wS
		  << "\twL = " << ltr_wL
		  << "\tw' = " << ltr_wp
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES "STEP 2: DO WHILE w_L - w_S > lambda"*/


      /* Step 3: Among all subranges within the whole range
	 [w_1 , w_2], find the subrange [w_a, w_b] that has the
	 largest ratio of D2(w_b)/D2(w_a). Let w_s = w_a and
	 w_L = w_b.
      */

#ifdef __VERBOSE_YES
      lpc_labelFunc =
	"STEP 3: FIND SUBRANGE [w_a, w_b] LARGEST RATIO D2(w_b)/D2(w_a)";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	geverbosepc_labelstep = lpc_labelFunc;
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      ltr_wS = aiinp_inParamSubClusterBinVk.getW1();
      ltr_wL = aiinp_inParamSubClusterBinVk.getW2();
      T_REAL ltr_wm = (ltr_wS + ltr_wL) / 2.0;
      
      T_REAL ltr_ratiosD2WmWs = ltr_D2Wm / ltr_D2Ws;
      T_REAL ltr_ratiosD2WlWm = ltr_D2Wl / ltr_D2Wm;
        
      if ( ltr_ratiosD2WmWs > ltr_ratiosD2WlWm ) {
	ltr_wL   = ltr_wm;
	ltr_D2Wl = ltr_D2Wm;
      }
      else {
	ltr_wS   = ltr_wm;
	ltr_D2Ws = ltr_D2Wm;
      }

#ifdef __VERBOSE_YES
      lchrom_wL = lchrom_wLtmp;
      lchrom_wS = lchrom_WStmp;
     
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	geverbosepc_labelstep = lpc_labelFunc;
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << "\n\tD2(Wm)/D2(wS)\t" << ltr_ratiosD2WmWs
		  << "\n\tD2(Wl)/D2(Wm)\t" << ltr_ratiosD2WlWm
		  << "\n\twS = " << ltr_wS
		  << "\n\twL = " << ltr_wL
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    /*Step 4: Do while w_L - w_s >  lambda '
      Begin
      Let w_m = (w_s + w_L )/2. Use CLUSTERING
      with the parameter w to cluster the data set.
      Calculate the ratios D_2(w_m)/D_2(w_s) and
      D_2(w_L)/D_2(w_m). Among all subranges within
      the whole range [w_1, w_2], find the subrange
      [w_a, w_b] that has the largest ratio of
      D_2(w_b)/D_2(w_a). Let w_s = w_a and w_L = w_b.
      End
      w" = wL.
      w''' = w_S.
    */

#ifdef __VERBOSE_YES
      lpc_labelFunc = "STEP 4: DO WHILE w_L - w_S > lambda";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	geverbosepc_labelstep = lpc_labelFunc;
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      while ( (ltr_wL -ltr_wS) > aiinp_inParamSubClusterBinVk.getLambda() ) { 
    
	ltr_wm = (ltr_wS + ltr_wL) / 2.0;

	gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> lchrom_Wm = 
	  clustering_genetic
	  (aoop_outParamGAC,
	   lmatrixrowt_Vi,
	   lvectort_numInstBi,
	   ltr_wm,
#ifdef __VERBOSE_YES
	   lmembclassdisjsets_Bi,
#endif /*__VERBOSE_YES*/
	   aiinp_inParamSubClusterBinVk,
	   aifunc2p_dist,
	   lii_numrunalg++,
	   let_executionTime
	   );

	std::tie(lomatrixrowt_S,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
	  clusteringop::getClusters
	  (lchrom_Wm,
	   lmatrixrowt_Vi,
	   lvectort_numInstBi,
	   aifunc2p_dist,
	   (T_CLUSTERIDX) 0
	   );

	T_REAL ltr_D2Wm = 
	  um::D2
	  (lomatrixrowt_S,
	   lmembclassdisjsets_BkinCi,
	   lmembclassdisjsets_Bi,
	   lmatrixrowt_Vi,
	   lvectorr_meanRadiusBk,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );

	T_REAL ltr_ratiosD2WmWs = ltr_D2Wm / ltr_D2Ws;
	T_REAL ltr_ratiosD2WlWm = ltr_D2Wl / ltr_D2Wm;

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "<HEURISTIC STRATEGY D2:>"
		    << "\nWs\t"         << ltr_wS
		    << "\tk\t"          << lchrom_wS.getNumBitOn()
		    << "\tD2Ws\t"       << ltr_D2Ws
		    << "\tfitness wS\t" << lchrom_wS.getFitness()
		    << "\nWl\t"         << ltr_wL
		    << "\tk\t"          << lchrom_wL.getNumBitOn()
		    << "\tD2Wl\t"       << ltr_D2Wl
		    << "\tfitness wL\t" << lchrom_wL.getFitness()
		    << "\nWm\t"         << ltr_wm
		    << "\tk\t"          << lchrom_Wm.getNumBitOn()
		    << "\tD2Wm\t"       << ltr_D2Wm
		    << "\tfitness wm\t" << lchrom_Wm.getFitness()
		    << "\n\tD2(Wm)/D2(wS)\t" << ltr_ratiosD2WmWs
		    << "\n\tD2(Wl)/D2(Wm)\t" << ltr_ratiosD2WlWm
		    << std::endl;
	}
	//--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
        
	if ( ltr_ratiosD2WmWs > ltr_ratiosD2WlWm ) {
	  ltr_wL   = ltr_wm;
	  ltr_D2Wl = ltr_D2Wm;

	  lchrom_wL = lchrom_Wm;
	}
	else {
	  ltr_wS = ltr_wm;
	  ltr_D2Ws = ltr_D2Wm;
	  lchrom_wS = lchrom_Wm;
	}

#ifdef __VERBOSE_YES
	// ++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "\t[" << ltr_wS << "," << ltr_wL << ']'
		    << std::endl;
	}
	//--geiinparam_verbose;
#endif /*__VERBOSE_YES*/


      } /* end while */
  
      T_REAL ltr_wpp  = ltr_wL;
      T_REAL ltr_wppp = ltr_wS;

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	geverbosepc_labelstep = lpc_labelFunc;
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << "\n\tw'' = "  << ltr_wpp
		  << "\n\tw''' = " << ltr_wppp  
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


      /* Step 5: If w'' >= w' Then Output the clustering obtained
	 with the parameter w'.
	 If w'' < w' Then Output the clustering obtained
	 with the parameter w'''.
      */


#ifdef __VERBOSE_YES
      lpc_labelFunc = "STEP 5: OUTPUT THE CLUSTERING";
      geverbosepc_labelstep = lpc_labelFunc;
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      T_REAL lotr_outW = 0;
      if ( ltr_wpp  >= ltr_wp  ) {
	lotr_outW = ltr_wp;

      }
      else { //if ( ltr_wpp  < ltr_wp )
	lotr_outW = ltr_wppp;
      }
#endif /*__VERBOSE_YES*/
      
      if ( ltr_wpp  >= ltr_wp  ) {
	lochrom_best = lchrom_wp;      
      }
      else {
	lochrom_best = lchrom_wL;
      }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       geverbosepc_labelstep = lpc_labelFunc;
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< " outW = " << lotr_outW
		<< std::endl;
      lochrom_best.print();
      std::cout	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "4. THE HEURISTIC STRATEGY TO FIND A GOOD CLUSTERING: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/ 
  }
  else {
    aoop_outParamGAC.setEndingCondition(false);
    aoop_outParamGAC.setRuntimeMessage
      ("componet set Bi is one in the dataset impossible to group");
  }

  runtime::stop(let_executionTime);
  aoop_outParamGAC.setNumClusterK
    ((T_CLUSTERIDX) lochrom_best.getNumBitOn() ); 
  aoop_outParamGAC.setMetricFuncRun
    ( lochrom_best.getObjetiveFunc() );
  aoop_outParamGAC.setFitness
    (lochrom_best.getFitness());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax )  {
    geverbosepc_labelstep = lpc_labelAlgHeuristic;
    std::cout << lpc_labelAlgHeuristic
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochrom_best.print();
    std::cout << std::endl;
	   	     
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

   return std::make_tuple
    (lochrom_best,
     lmatrixrowt_Vi,
     lvectort_numInstBi,
     lmembclassdisjsets_Bi
     );
    
} /* END clustering_vksubclusterbinary */

} /*END eac */
  
#endif /*__CLUSTERING_VKSUBCLUSTERBINARY_HPP__*/
