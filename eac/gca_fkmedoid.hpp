/*! \file gca_fkmedoid.hpp
 *
 * \brief GCA \cite Lucasius:etal:GAclusteringMedoid:GCA:1993
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GCA algorithm based on the paper:\n
 * C.B. Lucasius, A.D. Dane, and G. Kateman. On k-medoid\n
 * clustering of large data sets with the aid of a genetic\n 
 * algorithm: background, feasibility and comparison. Analytica\n
 * Chimica Acta, 282:647–669, 1993.\n
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


#ifndef __GCA_FKMEDOID_LUCASIUS_ETAL1993_HPP__
#define __GCA_FKMEDOID_LUCASIUS_ETAL1993_HPP__

#include <vector>
#include <algorithm>
#include <iterator>

#include <leac.hpp>
#include "outparam_eaclustering.hpp"
#include "inparam_gca.hpp"

#include "plot_runtime_function.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of genetic and evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace eac {
  
/*! \fn gaencode::ChromFixedLength<uintidx,T_REAL> gca_fkmedoid(inout::OutParamEAClustering<T_REAL,T_CLUSTERIDX> &aoop_outParamEAC, inout::InParamGCA<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiipcgagca_inParamGCA, const mat::MatrixTriang<T_REAL> &aimatrixtriagrt_dissimilarity)
  \brief  GCA \cite Lucasius:etal:GAclusteringMedoid:GCA:1993  
  \details Implementation of the GCA algorithm based on \cite Lucasius:etal:GAclusteringMedoid:GCA:1993
  \returns A chromosome with k genes, where each gene is the index of the most representative instance of each cluster
  \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
  \param aiipcgagca_inParamGCA a inout::InParamGCA parameters required by the algorithm
  \param aimatrixtriagrt_dissimilarity a triangular matrix with the distances between the instances
*/
template <typename T_CLUSTERIDX,
	  typename T_REAL,
	  typename T_FEATURE,         
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K
	  >
gaencode::ChromFixedLength<uintidx,T_REAL> 
gca_fkmedoid
(inout::OutParamEAClustering
 <T_REAL,
 T_CLUSTERIDX>                         &aoop_outParamEAC,
 inout::InParamGCA
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                 &aiinp_inParamGCA,
 const mat::MatrixTriang<T_REAL>        &aimatrixtriagrt_dissimilarity
 )
{
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC*/
  
  gaencode::ChromFixedLength<uintidx,T_REAL>::setStringSize
    ((uintidx) aiinp_inParamGCA.getNumClusterK() );

  gaencode::ChromFixedLength<uintidx,T_REAL> lochromfixleng_best;

  std::vector<gaencode::ChromFixedLength<uintidx,T_REAL>* >  
    lvectorchromfixleng_population;
  
  /*SPACE FOR STORE MATINGPOOL*/  
  std::vector<gaencode::ChromFixedLength<uintidx,T_REAL>* >
    lvectorchromfixleng_stringPool;
    
#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gca_fkmedoid";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA
      << "  IN(" << geiinparam_verbose << ")\n"
      << "\t(output Chromosome: lochromfixleng_best[" 
      << &lochromfixleng_best << "]\n"
      << "\t output inout::OutParamEAClustering&: aoop_outParamEAC[" 
      << &aoop_outParamEAC << "]\n"
      << "\t input  InParamClusteringGALucasius&: aiinp_inParamGCA[" 
      << &aiinp_inParamGCA << ']'
      << "\n\t\tPopulation size = " 
      << aiinp_inParamGCA.getSizePopulation()
      << "\n\t\tProbRecombination p_r = " 
      << aiinp_inParamGCA.getProbCrossover() 
      << "\n\t\tProb_m.point  = " 
      << aiinp_inParamGCA.getProbMutation()
      << "\n\t\tProbMixMutation  = " 
      << aiinp_inParamGCA.getProbMixMutation()
      << "\n\t input  mat::MatrixTriang<T_REAL>: &aimatrixtriagrt_dissimilarity[" 
      <<  &aimatrixtriagrt_dissimilarity << ']'
      << "\n\t)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 
  
  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamGCA.getNumMaxGenerations(), "Iterations", "Clustering metrics");
 
  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT
  
  std::ofstream lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;
  
  if ( aiinp_inParamGCA.getWithPlotStatObjetiveFunc() ) {  

    lvectorT_statfuncObjetiveFunc.reserve
      (aiinp_inParamGCA.getSizePopulation());  
    //DEFINE FUNCTION
    lofh_SSE  = new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinp_inParamGCA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);
    
    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat<T_REAL>
	( (char) li_i,
	  aiinp_inParamGCA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION 
    aoop_outParamEAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamGCA.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamGCA.getTimesRunAlgorithm()
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

  runtime::ExecutionTime let_executionTime  = runtime::start();

  std::uniform_int_distribution<uintidx> uniformdis_idxInstances
    (0,aimatrixtriagrt_dissimilarity.getNumRows()-1);
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  std::uniform_int_distribution<uintidx> uniformdis_uiMutation0N
    (0,gaencode::ChromFixedLength<uintidx,T_REAL>::stcgetStringSize()-1);
  
  /*POPULATION CREATE------------------------------------------------------------
   */
  lvectorchromfixleng_population.reserve
    (aiinp_inParamGCA.getSizePopulation());
  for (uintidx luintidx_i = 0; 
       luintidx_i < aiinp_inParamGCA.getSizePopulation(); 
       luintidx_i++) 
    {
      lvectorchromfixleng_population.push_back
	(new gaencode::ChromFixedLength<uintidx,T_REAL>());
    }

  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
   */
  lvectorchromfixleng_stringPool.reserve
    (aiinp_inParamGCA.getSizePopulation());
  for (uintidx luintidx_i = 0; 
       luintidx_i < aiinp_inParamGCA.getSizePopulation(); 
       luintidx_i++) 
    {
      lvectorchromfixleng_stringPool.push_back
	(new gaencode::ChromFixedLength<uintidx,T_REAL>());
    }

  /*-1. Initiate strings. The initial population is
    created. Each string in the population uniquely
    encodes a candidate solution.
  */
  
  { /*BEGIN INITIAL STRINGS*/
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep =  "-1. INITIAL STRINGS";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    const uintidx  luintidx_numClusterK =
      (uintidx) aiinp_inParamGCA.getNumClusterK();
    
    for (auto lchromfixleng_iter :lvectorchromfixleng_population) {
       
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
	 lchromfixleng_iter->getString()   
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

  } /*END INITIAL STRINGS*/

  
  while ( true ) {

    /*0. Evaluate strings. Each string in the popula-
      tion is decoded to obtain its actual meaning. This
      is then passed on to the objective function 
      $(f_0 = 1/\phi)$.
    */
    { /*BEGIN EVALUATE STRINGS */

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "0. EVALUATE STRINGS:  IN"
		  << '(' << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      std::for_each
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [&](gaencode::ChromFixedLength<uintidx,T_REAL>* lchromfixleng_iter)
	 {
	   T_REAL lT_objetiveFunc = 
	     um::SSEMedoid
	     (lchromfixleng_iter->getString(),
	      aiinp_inParamGCA.getNumClusterK(),
	      aimatrixtriagrt_dissimilarity
	      );
	   lchromfixleng_iter->setObjetiveFunc(lT_objetiveFunc);
	   lchromfixleng_iter->setFitness(1.0 / lchromfixleng_iter->getObjetiveFunc());
	   lchromfixleng_iter->setValidString(true);

#ifndef __WITHOUT_PLOT_STAT
	   lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
	
	 }
	 );

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "0. EVALUATE STRINGS: OUT"
		  << '(' << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END EVALUATE STRINGS */
    
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

      gaencode::ChromFixedLength<uintidx,T_REAL> *lchromfixleng_maxFitness  = 
	*(std::max_element
	  (lvectorchromfixleng_population.begin(),
	   lvectorchromfixleng_population.end(),
	   [&](const gaencode::ChromFixedLength<uintidx,T_REAL>* x, 
	       const gaencode::ChromFixedLength<uintidx,T_REAL>* y
	       ) 
	{  return x->getFitness() < y->getFitness(); } 
	   )
	  );

      if ( lochromfixleng_best.getFitness() <  lchromfixleng_maxFitness->getFitness() ) {

	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	lochromfixleng_best = *lchromfixleng_maxFitness;
    
	aoop_outParamEAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamEAC.setRunTimeGetsBest
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
    if ( aiinp_inParamGCA.getWithPlotStatObjetiveFunc() ) {  
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
      std::cout
	<< "END ITERATION: " << llfh_listFuntionHist.getDomainUpperBound()
	<< "\tobjetivoFunc = " << lochromfixleng_best.getObjetiveFunc() 
	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    /*GCA the termination criterion is based on a predetermined 
      maximum genetation
    */
    if ( llfh_listFuntionHist.getDomainUpperBound() >=
	 aiinp_inParamGCA.getNumMaxGenerations()) 
      break;
   
    llfh_listFuntionHist.increaseDomainUpperBound();
 
    /*1. Scale fitnesses.
     */

    /*2. Selectively reproduce strings.
     */

    { /*2. BEGIN SELECTIVELY REPRODUCE STRINGS*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "2. SELECTIVELY REPRODUCE STRINGS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      
      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchromfixleng_population.begin(),lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<uintidx,T_REAL>* lchromfixleng_iter) -> T_REAL
	 {
	   return lchromfixleng_iter->getFitness();
	 }
	 );
      
      
      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL--------------------------
       */
      for ( auto lchromfixleng_iter: lvectorchromfixleng_stringPool) {

	uintidx luiidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );
	
	*lchromfixleng_iter = *lvectorchromfixleng_population.at(luiidx_chrom);
      }
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


    } /*END 2. SELECTIVELY REPRODUCE STRINGS*/
    
    /*3. Pair strings.
     */

    {/*BEGIN A. RECOMBINE STRINGS*/

      /*a. Recombine strings.
       */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "A. RECOMBINE STRINGS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      std::uniform_int_distribution<uintidx> uniformdis_uintidx0PoolSize
	(0,lvectorchromfixleng_stringPool.size()-1);
      
      for (auto lchromfixleng_iter = lvectorchromfixleng_population.begin(); 
	   lchromfixleng_iter != lvectorchromfixleng_population.end();
	   ++lchromfixleng_iter)
	{
	  //GET TWO RANDOM
	  std::pair<uintidx,uintidx>
	    lpair_idxChrom =
	    prob::getRandPairUnlike
	    ([&]() -> uintidx
	     {
	       return uniformdis_uintidx0PoolSize(gmt19937_eng);
	     }
	     );
	  
	  gaencode::ChromFixedLength<uintidx,T_REAL>* lchromfixleng_child1 = *lchromfixleng_iter;
	  
	  if ( ++lchromfixleng_iter == lvectorchromfixleng_population.end()) {
	    *lchromfixleng_child1 = *lvectorchromfixleng_stringPool.at(lpair_idxChrom.first); 
	    break;
	  }

	  gaencode::ChromFixedLength<uintidx,T_REAL>* lchromfixleng_child2 = *lchromfixleng_iter;
	  
	  if ( uniformdis_real01(gmt19937_eng) //if  Crossover
	       < aiinp_inParamGCA.getProbCrossover() ) {

	    gaintegerop::recombinationD_MX
	      (*lchromfixleng_child1,
	       *lchromfixleng_child2,
	       *lvectorchromfixleng_stringPool.at(lpair_idxChrom.first),
	       *lvectorchromfixleng_stringPool.at(lpair_idxChrom.second),
	       aiinp_inParamGCA.getNumInstances(), 
	       aiinp_inParamGCA.getProbMixMutation()
	       );
	  	   
	  } //END if  Crossover
	  else {
	    *lchromfixleng_child1 = *lvectorchromfixleng_stringPool.at(lpair_idxChrom.first);
	    *lchromfixleng_child2 = *lvectorchromfixleng_stringPool.at(lpair_idxChrom.second);	
	  }
	  
	} /*for*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    }/*END A. RECOMBINE STRINGS*/

    /*b. Mutate strings.
     */

    { /*B. MUTATE STRINGS*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "B. MUTATE STRINGS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      /*Mutation. 
	D_PM
	C ----> C'
      */
      for (auto lchromfixleng_iter :lvectorchromfixleng_population) {

	if ( uniformdis_real01(gmt19937_eng) 
	     < aiinp_inParamGCA.getProbMutation() ) 
	  {	 
	    std::unordered_set<uintidx> lunorderedset_medoids
	      (lchromfixleng_iter->begin(),
	       lchromfixleng_iter->end()
	       );
	    uintidx luintidx_medoidRand =
	      prob::getRandSetUnlike
	      (lunorderedset_medoids,
	       [&]()
	       {
		 return uniformdis_idxInstances(gmt19937_eng);
	       }
	       );

	    uintidx lui_positionGene = uniformdis_uiMutation0N(gmt19937_eng); 
     
	    lchromfixleng_iter->setGene(lui_positionGene,luintidx_medoidRand);
	    
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

    } /*END MUTATION */

    /*4. Selectively replace strings.
     */

    { /*BEGIN 4. SELECTIVELY REPLACE STRINGS*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "4. SELECTIVELY REPLACE STRINGS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END 4. SELECTIVELY REPLACE STRINGS*/


  } /*END EVOLUTION While*/ 

  
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
    
    for (uintidx lui_i = 0; lui_i < lvectorchromfixleng_stringPool.size(); ++lui_i) {
      delete lvectorchromfixleng_stringPool[lui_i];
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
  aoop_outParamEAC.setNumClusterK
    (aiinp_inParamGCA.getNumClusterK());
  aoop_outParamEAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());


#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamGCA.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamGCA,
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

  
} /*END  gca_fkmedoid  */

} /*END eac */

#endif /*__GCA_FKMEDOID_LUCASIUS_ETAL1993_HPP__*/
