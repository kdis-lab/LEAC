/*! \file hka_fkmedoid.hpp
 *
 * \brief HKA \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the HKA algorithm based on the paper:\n
 * Weiguo Sheng and Xiaohui Liu. A hybrid algorithm for\n 
 * k-medoid clustering of large data sets. In Evolutionary\n 
 * Computation, 2004. CEC2004. Congress on, volume 1,pages 77–82\n 
 * Vol.1, June 2004.\n 
 * <a href="http://dx.doi.org/10.1109/CEC.2004.1330840">doi:http://dx.doi.org/10.1109/CEC.2004.1330840</a>.\n
 * \n
 * Library Evolutionary Algorithms for Clustering (LEAC) is a library\n
 * for the implementation of evolutionary and genetic algorithms\n
 * focused on the partition type clustering problem. Based on the\n
 * current standards of the <a href="http://en.cppreference.com">C++</a> language, as well as on Standard\n
 * Template Library <a href="http://en.cppreference.com/w/cpp/container">STL</a> 
 * and also  <a href="http://www.openblas.net/">OpenBLAS</a> to have a better performance.\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> \endlink license
 */

#ifndef __HKA_FKMEDOID_SHENG_LUI_2004_HPP__
#define __HKA_FKMEDOID_SHENG_LUI_2004_HPP__

#include <vector>
#include <algorithm>    // std::sort

#include <leac.hpp>
#include "inparam_gamedoidclustering_hka.hpp"
#include "outparam_gaclustering.hpp"

#include "plot_runtime_function.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of genetic and evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace eac {
  
/*! \fn gaencode::ChromFixedLength<uintidx,T_REAL> hka_fkmedoid (inout::OutParamGAClustering<T_REAL,T_CLUSTERIDX> &aoopcga_outParamClusteringGA, inout::InParamGAMedoidClusteringHKA<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiipcgapc_inParamHKA, mat::MatrixTriang<T_REAL> &aimatrixtriagrt_dissimilarity)
  \brief HKA \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004
  \details Implementation of the HKA algorithm based on \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004, Find a partition based on the most representative instances, also called prototypes 
  \returns A chromosome with k genes, where each gene is the index of the most representative instance of each cluster
  \param aoopcga_outParamClusteringGA a inout::OutParamGAClustering with the output parameters of the algorithm
  \param aiipcgapc_inParamHKA a inout::InParamGAMedoidClusteringHKA parameters required by the algorithm
  \param aimatrixtriagrt_dissimilarity a triangular matrix with the distances between the instances
*/
template <typename T_CLUSTERIDX,
	  typename T_REAL,
	  typename T_FEATURE,         
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K
	  >
gaencode::ChromFixedLength<uintidx,T_REAL> 
hka_fkmedoid
(inout::OutParamGAClustering
 <T_REAL,
 T_CLUSTERIDX>                          &aoopcga_outParamClusteringGA,
 inout::InParamGAMedoidClusteringHKA
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                 &aiipcgapc_inParamHKA,
 mat::MatrixTriang<T_REAL>              &aimatrixtriagrt_dissimilarity
 )
{  
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC*/
  gaencode::ChromFixedLength<uintidx,T_REAL>::setStringSize
    ( (uintidx) aiipcgapc_inParamHKA.getNumClusterK() );

  gaencode::ChromFixedLength<uintidx,T_REAL> lochromfixleng_best;

  std::vector<gaencode::ChromFixedLength<uintidx,T_REAL>* >  
    lvectorchromfixleng_population;
  
  /*SPACE FOR STORE MATINGPOOL*/  
  std::vector<gaencode::ChromFixedLength<uintidx,T_REAL>* >
    lvectorchromfixleng_matingPool;
  lvectorchromfixleng_matingPool.reserve(aiipcgapc_inParamHKA.getSizePopulation()/2);

  std::vector<gaencode::ChromFixedLength<uintidx,T_REAL>* >
    lvectorchromfixleng_newoffspring; 

#ifdef __VERBOSE_YES
  
  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "hka_fkmedoid"; 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ")\n"
       	      << "\t(output Chromosome: lochromfixleng_best[" 
	      << &lochromfixleng_best << "]\n"
	      << "\t output inout::OutParamGAClustering&: aoopcga_outParamClusteringGA[" 
	      << &aoopcga_outParamClusteringGA << "]\n"
	      << "\t input  InParamGAMedoidClusteringHKA&: aiipcgapc_inParamHKA[" 
	      << &aiipcgapc_inParamHKA << ']'
	      << "\n\t\tNumClusterK = " 
	      << aiipcgapc_inParamHKA.getNumClusterK()
	      << "\n\t\tPopulation size = " 
	      << aiipcgapc_inParamHKA.getSizePopulation()
	      << "\n\t\tProbRecombination p_r = " 
	      << aiipcgapc_inParamHKA.getProbCrossover() 
	      << "\n\t\tProb_m.point  = " 
	      << aiipcgapc_inParamHKA.getProbMutation()
	      << "\n\t\tProbMixMutation  = " 
	      << aiipcgapc_inParamHKA.getProbMixMutation()
	      << "\n\t input  mat::MatrixTriang<T_REAL>: &aimatrixtriagrt_dissimilarity[" 
	      <<  &aimatrixtriagrt_dissimilarity << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 
 
  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiipcgapc_inParamHKA.getNumMaxGenerations(), "Iterations", "Clustering metrics"); 

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT
  
  std::ofstream lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;
  
  if ( aiipcgapc_inParamHKA.getWithPlotStatObjetiveFunc() ) {  

    lvectorT_statfuncObjetiveFunc.reserve
      (aiipcgapc_inParamHKA.getSizePopulation());  
    //DEFINE FUNCTION
    lofh_SSE  = new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiipcgapc_inParamHKA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);
    
    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat<T_REAL>
	( (char) li_i,
	  aiipcgapc_inParamHKA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION 
    aoopcga_outParamClusteringGA.setFileNameOutPlotStatObjetiveFunc
      (aiipcgapc_inParamHKA.getFileNamePlotStatObjetiveFunc(),
       aiipcgapc_inParamHKA.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open   
      (aoopcga_outParamClusteringGA.getFileNameOutPlotStatObjetiveFunc().c_str(), 
       std::ios::out | std::ios::app
       );
    lfileout_plotStatObjetiveFunc.precision(COMMON_COUT_PRECISION);

    //FUNCTION HEADER
    lfileout_plotStatObjetiveFunc 
      << llfh_listFuntionHist.getHeaderFuntions() 
      << "\n";
  }

#endif /*__WITHOUT_PLOT_STAT*/

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  std::uniform_int_distribution<uintidx> uniformdis_idxInstances
    (0,aimatrixtriagrt_dissimilarity.getNumRows()-1);
  std::uniform_int_distribution<uintidx> uniformdis_uiMutation0N
    (0,gaencode::ChromFixedLength<uintidx,T_REAL>::stcgetStringSize()-1);

  runtime::ExecutionTime let_executionTime = runtime::start();
  
  /*POPULATION CREATE------------------------------------------------------------
   */
  lvectorchromfixleng_population.reserve
    (aiipcgapc_inParamHKA.getSizePopulation()+1);
  for (uintidx luintidx_i = 0; 
       luintidx_i < aiipcgapc_inParamHKA.getSizePopulation(); 
       luintidx_i++) 
    {
      lvectorchromfixleng_population.push_back
	(new gaencode::ChromFixedLength<uintidx,T_REAL>());
    }
  
  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
   */
  lvectorchromfixleng_newoffspring.reserve
    ( aiipcgapc_inParamHKA.getSizePopulation()/2 + 1 );
  
  for (uintidx luintidx_i = 0; 
       luintidx_i < aiipcgapc_inParamHKA.getSizePopulation()/2; 
       luintidx_i++) 
    {
      lvectorchromfixleng_newoffspring.push_back
	(new gaencode::ChromFixedLength<uintidx,T_REAL>());
    }
 
  /*Step 1. Randomly initialize a population with P
    individuals. Each individual uniquely encodes a
    candidate solution of the k-medoid clustering
    problem.
  */

  { /*STEP 1. INITIAL A POPULATION*/
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "STEP 1. INITIAL A POPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    const uintidx  luintidx_numClusterK =
      (uintidx) aiipcgapc_inParamHKA.getNumClusterK();
    
    for (auto ichrom_population :lvectorchromfixleng_population) {
       
      std::unordered_set<uintidx> &&lunorderedset_medoids =
	prob::getWithoutRepeatsSet
	( luintidx_numClusterK,
	  [&]() -> uintidx
	  {
	    return uniformdis_idxInstances(gmt19937_eng);
	  }
	  );

      std::copy_n
	(lunorderedset_medoids.begin(),
	 lunorderedset_medoids.size(),
	 ichrom_population->getString()   
	 );

      ichrom_population->setObjetiveFunc(std::numeric_limits<T_REAL>::max());
      ichrom_population->setFitness(-std::numeric_limits<T_REAL>::max());

    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END STEP 1. INITIAL A POPULATION*/

  /*Step 2. Calculate SED value according to equation (2)
    for each individual in the initial population and
    set the fitness value as $f = 1/SED$.
    \[ SED = \sum_{i=1}^{n}\sum_{j=1,x_i\in C_j}^{k} d(x_i,m_j) (2) \]
  */
  
  { /*BEGIN STEP 2. CALCULATE SED VALUE ACCORDING TO EQUATION (2) */

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "STEP 2. CALCULATE SED VALUE";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

      
    for ( auto ichrom_population:  lvectorchromfixleng_population ) {
	 	 
      T_REAL lT_objetiveFunc = 
	um::SSEMedoid
	(ichrom_population->getString(),
	 aiipcgapc_inParamHKA.getNumClusterK(),
	 aimatrixtriagrt_dissimilarity
	 );
      ichrom_population->setObjetiveFunc(lT_objetiveFunc);
      ichrom_population->setFitness(1.0 / ichrom_population->getObjetiveFunc());
      ichrom_population->setValidString(true);

#ifndef __WITHOUT_PLOT_STAT
      lvectorT_statfuncObjetiveFunc.push_back(ichrom_population->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
      
    }
      
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep 
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END STEP 2. CALCULATE SED VALUE ACCORDING TO EQUATION (2) */

 
  { /*BEGIN FIND THE BEST CANDIDATE SOLUTION*/
      
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "FIND THE BEST CANDIDATE SOLUTION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    gaencode::ChromFixedLength<uintidx,T_REAL> *lchrom_maxFitness =
      *(std::max_element
	( lvectorchromfixleng_population.begin(),
	  lvectorchromfixleng_population.end(),
	  [&](const gaencode::ChromFixedLength<uintidx,T_REAL>* x, 
	      const gaencode::ChromFixedLength<uintidx,T_REAL>* y
	      ) 
      {  return x->getFitness() < y->getFitness(); } 
	  )
	);

    if ( lochromfixleng_best.getFitness() <  lchrom_maxFitness->getFitness() ) {
      /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
      lochromfixleng_best = *lchrom_maxFitness;	
      aoopcga_outParamClusteringGA.setIterationGetsBest
	(llfh_listFuntionHist.getDomainUpperBound());
      aoopcga_outParamClusteringGA.setRunTimeGetsBest
	(runtime::elapsedTime(let_executionTime));
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep 
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END FIND THE BEST CANDIDATE SOLUTION*/ 

    /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT  
  if ( aiipcgapc_inParamHKA.getWithPlotStatObjetiveFunc() ) {  
    lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
    functionhiststat_evaluateAll
      (lofhs_statObjectiveFunc,
       lvectorT_statfuncObjetiveFunc
       );
    lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
    lvectorT_statfuncObjetiveFunc.clear();
  }
#endif /*__WITHOUT_PLOT_STAT */


  /*Step 3. Repeat (a) to (e) until the stopping criterion is
    met.
  */

  do {

    /*a) Use 2-fold tournament selection to select P/2
      parent pairs for reproduction.
    */

    { /*BEGIN TOURNAMENT SELECTION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "A). USE 2-FOLD TOURNAMENT SELECTION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      lvectorchromfixleng_matingPool.clear();
  
      for (uintidx luintidx_i = 0; 
	   luintidx_i < lvectorchromfixleng_population.size()/2; 
	   luintidx_i++) {
      
	gaencode::ChromFixedLength<uintidx,T_REAL>* lchrom_selectTour = 
	  *(gaselect::tournament
	    (lvectorchromfixleng_population.begin(),
	     lvectorchromfixleng_population.end(),
	     aiipcgapc_inParamHKA.getOrderTournament(),
	     [&](const gaencode::ChromFixedLength<uintidx,T_REAL>* iter_chrom)
	     {
	        return iter_chrom->getFitness();
	     }
	     )
	    );
	lvectorchromfixleng_matingPool.push_back(lchrom_selectTour);
      }

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END TOURNAMENT SELECTION*/
    
      /*b) Do mix subset recombination crossover on
	parent pairs and then perform flip mutation
	on each feature of the offspring with some
	probability.

	The crossover is applied on each paired parent with
	probability of $P_c$.
      */

    { /*BEGIN B) DO MIX SUBSET RECOMBINATION CROSSOVER*/
     
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "B) DO MIX SUBSET RECOMBINATION CROSSOVER";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      gaiterator::crossover
	(lvectorchromfixleng_matingPool.begin(),
	 lvectorchromfixleng_matingPool.end(),
	 lvectorchromfixleng_newoffspring.begin(),
	 lvectorchromfixleng_newoffspring.end(),
	 [&](const gaencode::ChromFixedLength<uintidx,T_REAL>*
	     aichrom_parent1,
	     const gaencode::ChromFixedLength<uintidx,T_REAL>*
	     aichrom_parent2,
	     gaencode::ChromFixedLength<uintidx,T_REAL>*
	     aochrom_child1, 
	     gaencode::ChromFixedLength<uintidx,T_REAL>*
	     aochrom_child2
	     )
	 {
	 
	   if ( uniformdis_real01(gmt19937_eng) <
		aiipcgapc_inParamHKA.getProbCrossover()  ) {
	      
	     gaintegerop::recombinationD_MX
	       (*aochrom_child1,
		*aochrom_child2,
		*aichrom_parent1,
		*aichrom_parent2,
		aiipcgapc_inParamHKA.getNumInstances(),
		aiipcgapc_inParamHKA.getProbMixMutation()
		);

	   } //if  Crossover
	   else {
	     *aochrom_child1 = *aichrom_parent1;
	     *aochrom_child2 = *aichrom_parent2;
	   }
	 }
	 );
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END   B) DO MIX SUBSET RECOMBINATION CROSSOVER*/
    
      /*After crossover, the probability $P_m$, of flip mutation will
	be applied to the offspring. Flip mutation replaces the
	chosen feature by another randomly generated feature,
	subject to the restriction that the new feature is not
	presented in the current genotype of the individual.
      */

    { /*BEGIN B) FLIP MUTATE*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "B) FLIP MUTATE";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto ichrom_newoffspring: lvectorchromfixleng_newoffspring ) {
	
	if ( uniformdis_real01(gmt19937_eng) 
	     < aiipcgapc_inParamHKA.getProbMutation() ) 
	  {

	    std::unordered_set<uintidx>
	      lunorderedset_idxInstances
	      (ichrom_newoffspring->begin(),
	       ichrom_newoffspring->end()
	       );

	    uintidx luintidx_newMedoid =
	      prob::getRandSetUnlike
	      (lunorderedset_idxInstances,
	       [&]()
	       {
		 return uniformdis_idxInstances(gmt19937_eng);
	       }
	       );

	    uintidx lui_positionGene =
	      uniformdis_uiMutation0N(gmt19937_eng);
	    
	    ichrom_newoffspring->setGene(lui_positionGene,luintidx_newMedoid);
	    
	    
	  } 
      }
       
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "B) FLIP MUTATE: OUT"
		  << '(' << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END B) FLIP MUTATE */

      /*c) Run one step of the local search heuristic on
	the new offspring with some probability and
	update the offspring.
      */
    
    { /*BEGIN C) LOCAL SEARCH HEURISTIC*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "C) LOCAL SEARCH HEURISTIC";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      
      for ( auto ichrom_newoffspring: lvectorchromfixleng_newoffspring ) {
	
	if ( uniformdis_real01(gmt19937_eng)  
	     < aiipcgapc_inParamHKA.getProbSearchHeuristic()  ) {

	  clusteringop::updateMedoids
	    (ichrom_newoffspring->getString(),
	     aiipcgapc_inParamHKA.getNumClusterK(),
	     aiipcgapc_inParamHKA.getNearestNeighbors(),
	     aimatrixtriagrt_dissimilarity
	     );
	} 
      }

      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END C) LOCAL SEARCH HEURISTIC*/

      /*d) Calculate SED value according to equation
	(2) for each offspring member and set fitness
	of the member equal to 1/SED,
      */

    { /*BEGIN D) CALCULATE SED VALUE FOR EACH OFFSPRING*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "D) CALCULATE SED VALUE FOR EACH OFFSPRING";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto&& ichrom_newoffspring: lvectorchromfixleng_newoffspring ) {
     
	T_REAL lT_objetiveFunc = 
	  um::SSEMedoid
	  (ichrom_newoffspring->getString(),
	   aiipcgapc_inParamHKA.getNumClusterK(),
	   aimatrixtriagrt_dissimilarity
	   );
	ichrom_newoffspring->setObjetiveFunc(lT_objetiveFunc);
	ichrom_newoffspring->setFitness(1.0 / ichrom_newoffspring->getObjetiveFunc());
	ichrom_newoffspring->setValidString(true);

#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(ichrom_newoffspring->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
	 
      }
       

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END STEP 2. CALCULATE SED VALUE ACCORDING TO EQUATION (2) */

      /*e) Create the new generation of size P from the
	one best member of the previous generation
	and the best offspring.
      */
    { /*BEGIN E) CREATE THE NEW GENERATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "E) CREATE THE NEW GENERATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif //__VERBOSE_YES
      
      std::vector<gaencode::ChromFixedLength<uintidx,T_REAL>* >
	lvectorchromfixleng_select;
      lvectorchromfixleng_select.reserve( aiipcgapc_inParamHKA.getSizePopulation() );

      std::sort
	( lvectorchromfixleng_population.begin(),
	  lvectorchromfixleng_population.end(),
	  [](const gaencode::ChromFixedLength<uintidx,T_REAL>* x, 
	     const gaencode::ChromFixedLength<uintidx,T_REAL>* y
	     ) 
	  {  return x->getFitness() > y->getFitness(); } 
	  );
      
      std::sort
	( lvectorchromfixleng_newoffspring.begin(),
	  lvectorchromfixleng_newoffspring.end(),
	  [](const gaencode::ChromFixedLength<uintidx,T_REAL>* x, 
	     const gaencode::ChromFixedLength<uintidx,T_REAL>* y
	     ) 
	  {  return x->getFitness() > y->getFitness(); } 
	  );
       
      //CENTINEL CROMOSOMES
      gaencode::ChromFixedLength<uintidx,T_REAL> *lchrombase_centinelaLeft =
	new gaencode::ChromFixedLength<uintidx,T_REAL>();
      gaencode::ChromFixedLength<uintidx,T_REAL> *lchrombase_centinelaRight =
	new gaencode::ChromFixedLength<uintidx,T_REAL>();

      lvectorchromfixleng_population.push_back(lchrombase_centinelaLeft); 
      lvectorchromfixleng_newoffspring.push_back(lchrombase_centinelaRight);
       
      uintidx luintidx_l = 0;
      uintidx luintidx_r = 0;

      for (uintidx lui_i = 0;
	   lui_i < aiipcgapc_inParamHKA.getSizePopulation();
	   lui_i++)
	{

	  if ( lvectorchromfixleng_population[luintidx_l]->getFitness() <
	       lvectorchromfixleng_newoffspring[luintidx_r]->getFitness() )
	    {
	 
#ifdef __VERBOSE_YES
	      ++geiinparam_verbose;
	      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		std::cout
		  << " lvectorchromfixleng_select[" << lui_i <<  ']' 
		  << " <--  lvectorchromfixleng_newoffspring[" << luintidx_r << ']' 
		  << '[' << & lvectorchromfixleng_newoffspring[luintidx_r] << ']' 
		  << " Obj: "
		  <<  lvectorchromfixleng_newoffspring[luintidx_r]->getObjetiveFunc() 
		  << '\n';
	      }
	      --geiinparam_verbose;
#endif //__VERBOSE_YES

	      lvectorchromfixleng_select.push_back
		( new gaencode::ChromFixedLength<uintidx,T_REAL>
		  (*lvectorchromfixleng_newoffspring[luintidx_r])
		  );
	      ++luintidx_r;
	    }
	  else {

#ifdef __VERBOSE_YES
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      std::cout
		<< " lvectorchromfixleng_select[" << lui_i <<  ']' 
		<< " <--  lvectorchromfixleng_population[" << luintidx_l << ']' 
		<< "[" << & lvectorchromfixleng_population[luintidx_l] << ']' 
		<< " Obj: "
		<<  lvectorchromfixleng_population[luintidx_l]->getObjetiveFunc() 
		<< '\n';
	    }
	    --geiinparam_verbose;	     
#endif //__VERBOSE_YES

	    lvectorchromfixleng_select.push_back
	      ( new gaencode::ChromFixedLength<uintidx,T_REAL>
		(*lvectorchromfixleng_population[luintidx_l])
		);
	    ++luintidx_l;
	  }
	}

      lvectorchromfixleng_population.pop_back();
      lvectorchromfixleng_newoffspring.pop_back();
       
      delete lchrombase_centinelaLeft;
      delete lchrombase_centinelaRight;

       
      for ( auto  ichrom_population = lvectorchromfixleng_population.begin()
	      , jchrom_select = lvectorchromfixleng_select.begin()
	      ; ichrom_population != lvectorchromfixleng_population.end()
	      ; ++ichrom_population, ++jchrom_select
	    )
	{
	   
	  *(*ichrom_population) = *(*jchrom_select);
	     
#ifndef __WITHOUT_PLOT_STAT
	  lvectorT_statfuncObjetiveFunc.push_back
	    ( (*ichrom_population)->getObjetiveFunc());
#endif //__WITHOUT_PLOT_STAT

	}
      
      //DELETE SPACE SELECT
      for (uintidx lui_i = 0; lui_i < lvectorchromfixleng_select.size(); ++lui_i) {
	delete lvectorchromfixleng_select[lui_i];
      }
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
      
      
    } /*END E) CREATE THE NEW GENERATION*/


    { /*BEGIN FIND THE BEST CANDIDATE SOLUTION*/
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "FIND THE BEST CANDIDATE SOLUTION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
	lvectorchromfixleng_population[0]->print();
	std::cout << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      if ( lochromfixleng_best.getFitness() <
	   lvectorchromfixleng_population[0]->getFitness() )
	{
	  /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	  lochromfixleng_best = *lvectorchromfixleng_population[0];
	 
	  aoopcga_outParamClusteringGA.setIterationGetsBest
	    (llfh_listFuntionHist.getDomainUpperBound());
	  aoopcga_outParamClusteringGA.setRunTimeGetsBest
	    (runtime::elapsedTime(let_executionTime));
	}

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END FIND THE BEST CANDIDATE SOLUTION*/

    llfh_listFuntionHist.increaseDomainUpperBound();

    /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT  
    if ( aiipcgapc_inParamHKA.getWithPlotStatObjetiveFunc() ) {  
      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/


#ifdef __VERBOSE_YES
    /*ID PROC
     */
    ++geverboseui_idproc;
 
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< "END ITERATION: "
	<< llfh_listFuntionHist.getDomainUpperBound()
	<< "\tobjetivoFunc = "
	<< lochromfixleng_best.getObjetiveFunc() 
	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  } while ( llfh_listFuntionHist.getDomainUpperBound() <
	    aiipcgapc_inParamHKA.getNumMaxGenerations() );  
  /* END  REPEAT */

  
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

    for (uintidx lui_i = 0; lui_i <  lvectorchromfixleng_population.size(); ++lui_i) {
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
    geverbosepc_labelstep = "DELETESNEWOFFSPRING";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    for (uintidx lui_i = 0; lui_i <  lvectorchromfixleng_newoffspring.size(); ++lui_i) {
      delete lvectorchromfixleng_newoffspring[lui_i];
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }/*END FREE MEMORY OF STRINGPOOL*/
  
  runtime::stop(let_executionTime);
  aoopcga_outParamClusteringGA.setNumClusterK
    (aiipcgapc_inParamHKA.getNumClusterK());
  aoopcga_outParamClusteringGA.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoopcga_outParamClusteringGA.setFitness
    (lochromfixleng_best.getFitness());
  aoopcga_outParamClusteringGA.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoopcga_outParamClusteringGA.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
 
  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiipcgapc_inParamHKA.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiipcgapc_inParamHKA,
       aoopcga_outParamClusteringGA
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout
      << lpc_labelAlgGA 
      << ": OUT("
      << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelChrom;
    lostrstream_labelChrom
      << "<CHROMOSMEBEST:"
      << lpc_labelAlgGA
      << ":generation: "
      <<  llfh_listFuntionHist.getDomainUpperBound();
    lochromfixleng_best.print(std::cout,lostrstream_labelChrom.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lochromfixleng_best;

} /*END  hka_fkmedoid  */

} /*END eac */

#endif /*__HKA_FKMEDOID_SHENG_LUI_2004_HPP__*/
