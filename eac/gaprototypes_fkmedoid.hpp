/*! \file gaprototypes_fkmedoid.hpp
 *
 * \brief GA-Prototypes \cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GA-Prototypes algorithm based on the paper:\n
 * L. I. Kuncheva and J. C. Bezdek. Selection of cluster\n
 * prototypes from data by a genetic algorithm. In inProc.\n 
 * 5th Eur. Congr. Intell. Tech. Soft Comput., pages\n
 * 1683–1688, 1997.\n
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
 
#ifndef __GAPROTOTYPES_FKMEDOID_HPP__
#define __GAPROTOTYPES_FKMEDOID_HPP__

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <leac.hpp>
#include "inparam_gaprototypesfk.hpp"
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
  

/*! \fn gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> gaprototypes_fkmedoid(inout::OutParamGAC<T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamGAPrototypes<T_BITSIZE,T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamGAPrototypes, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE>   &aifunc2p_dist)
  \brief GA-Prototypes \cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997
  \details Implementation of the GA-Prototypes algorithm based on \cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997.
  \returns A partition of a data set, encoded in a binary chromosome, if the value of the gene is "1" in the i-th position the i-th instance is a medoid
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinp_inParamGAPrototypes a inout::InParamGAPrototypes parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distance
*/
template < typename T_BITSIZE,
           typename T_REAL,// J1,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,  //0, 1, .., N
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
gaprototypes_fkmedoid
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                  &aoop_outParamGAC,
 inout::InParamGAPrototypesFk
 <T_BITSIZE,
 T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>         &aiinp_inParamGAPrototypes,
 const INPUT_ITERATOR           aiiterator_instfirst,
 const INPUT_ITERATOR           aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>   &aifunc2p_dist
 )
{
  const uintidx  luintidx_numClusterK =
    (uintidx) aiinp_inParamGAPrototypes.getNumClusterK();
  const T_REAL   lrt_numClusterK = (T_REAL) aiinp_inParamGAPrototypes.getNumClusterK();
  const uintidx  luintidx_numIntances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));

  /*CONVERT INSTANCES TO FORMAT MATRIX
   */
  mat::MatrixRow<T_FEATURE>&& lmatrixt_y =
    data::toMatrixRow
    (aiiterator_instfirst,
     aiiterator_instlast);
    
  gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    lochrombitarray_best( luintidx_numIntances );

  lochrombitarray_best.setFitness(std::numeric_limits<T_REAL>::max());

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0,1);
  
  if ( aiinp_inParamGAPrototypes.getPini() < 0) {
    aiinp_inParamGAPrototypes.setPini
      ( (T_REAL)luintidx_numClusterK / (T_REAL)luintidx_numIntances); 
  }
  
#ifdef __VERBOSE_YES
  /*ID PROC
   */
  geverboseui_idproc = 1;

  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gaprototypes_fkmedoid";  
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA
      << "  IN(" << geiinparam_verbose << ")\n"
      << "\t(output lochrombitarray_best[" 
      << &lochrombitarray_best << "]\n"
      << "\t output inout::OutParamGAC&: aoop_outParamGAC[" 
      << &aoop_outParamGAC << "]\n"
      << "\t input  InParamClusteringGABezdek1994&: aiinp_inParamGAPrototypes[" 
      << &aiinp_inParamGAPrototypes << "]\n"
      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
      << "\t input  dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist[" 
      << &aifunc2p_dist << ']'
      << "\n\t\tPopulation size = " 
      << aiinp_inParamGAPrototypes.getSizePopulation()
      << "\n\t\t Generations  = "
      << aiinp_inParamGAPrototypes.getNumMaxGenerations()
      << "\n\t\t Pini  = "
      << aiinp_inParamGAPrototypes.getPini()
      << "\n\t\t alpha  = "
      << aiinp_inParamGAPrototypes.getAlpha()
      << "\n\t\t ProbCrossover = "
      << aiinp_inParamGAPrototypes.getProbCrossover() 
      << "\n\t\t ProbMutation  = "
      << aiinp_inParamGAPrototypes.getProbMutation()
      << "\n\t\trandom-seed = "
      << aiinp_inParamGAPrototypes.getRandomSeed()
      << "\n\t)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  runtime::ListRuntimeFunction<COMMON_IDOMAIN>
    llfh_listFuntionHist
    (aiinp_inParamGAPrototypes.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream               lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_J1 = NULL;
  runtime::RuntimeFunctionValue
    <T_INSTANCES_CLUSTER_K>   *lofh_errorBezdek = NULL; /*function extra*/
  runtime::RuntimeFunctionStat<T_REAL>
    *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;
  
  if ( aiinp_inParamGAPrototypes.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      ( aiinp_inParamGAPrototypes.getSizePopulation());
    //DEFINE FUNCTION
    lofh_J1  = new runtime::RuntimeFunctionValue<T_REAL>
      ("J1", 
       aiinp_inParamGAPrototypes.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_J1);

    if ( aiinp_inParamGAPrototypes.getClassInstanceColumn()  ) {
      lofh_errorBezdek = new runtime::RuntimeFunctionValue<T_INSTANCES_CLUSTER_K>
	("ErrorBezdek", 
	 aiinp_inParamGAPrototypes.getAlgorithmoName(),
	 RUNTIMEFUNCTION_NOT_STORAGE
	 );
      llfh_listFuntionHist.addFuntion(lofh_errorBezdek);
    }

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamGAPrototypes.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamGAPrototypes.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamGAPrototypes.getTimesRunAlgorithm()
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

  runtime::ExecutionTime let_executionTime = runtime::start();
  
  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* >  lvectorchrom_population;
  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* >  lvectorchrom_matingPool;
  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* >  lvectorchrom_setO;

  /*SPACE FOR STORE POPULATION
   */
  lvectorchrom_population.reserve
    ( aiinp_inParamGAPrototypes.getSizePopulation() );

  /*SPACE FOR STORE MATINGPOOL
   */
  lvectorchrom_matingPool.reserve
    ( aiinp_inParamGAPrototypes.getSizePopulation() );

  lvectorchrom_setO.reserve
    ( aiinp_inParamGAPrototypes.getSizePopulation() );

  
  {/*BEGIN 1. INITIALIZE POPULATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "POPULATION INITIALIZATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    long ll_invalidOffspring = 0;

    /*CREATE SPACE FOR STORE POPULATION-----------------------------------------
     */
    for (uintidx lui_i = 0; 
	 lui_i < aiinp_inParamGAPrototypes.getSizePopulation(); 
	 lui_i++) 
      {
	lvectorchrom_population.push_back
	  (new gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>(luintidx_numIntances));

	lvectorchrom_matingPool.push_back
	  (new gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>(luintidx_numIntances));

	lvectorchrom_setO.push_back
	  (new gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>(luintidx_numIntances));
	
      }

    for ( auto lchrombitarray_iter: lvectorchrom_population) {
      
      gabinaryop::initializeGenes
	(*lchrombitarray_iter,
	 [&]() 
	 {
	   return (uniformdis_real01(gmt19937_eng)
		   <= aiinp_inParamGAPrototypes.getPini());
	 }
	 );

      std::list<uintidx>&& listui_idxInst =
	lchrombitarray_iter->getIdxWithBitOn();

      T_REAL lT_j1 =  std::numeric_limits<T_REAL>::max();

      if ( luintidx_numClusterK != (uintidx) listui_idxInst.size() ) 
	++ll_invalidOffspring;
	     
      if ( listui_idxInst.size() > 1 ) {
	mat::MatrixRow<T_FEATURE> 
	  lmatrixt_v                         
	  ( (uintidx) listui_idxInst.size(),
	    data::Instance<T_FEATURE>::getNumDimensions()
	    );

	clusteringop::initialize
	  (lmatrixt_v,
	   aiiterator_instfirst,
	   listui_idxInst.begin()
	   );

	mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>
	  lbcrispmatrix_w(lmatrixt_v.getNumRows(),luintidx_numIntances);

	clusteringop::getPartition
	  (lbcrispmatrix_w,
	   lmatrixt_y,
	   lmatrixt_v,
	   aifunc2p_dist
	   );

	lT_j1 = 
	  um::j1
	  (lbcrispmatrix_w,
	   lmatrixt_v, 
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );

	T_REAL lrt_pmod = (T_REAL) lmatrixt_v.getNumRows();
   
	lT_j1 +=  aiinp_inParamGAPrototypes.getAlpha()
	  *(lrt_pmod-lrt_numClusterK)*(lrt_pmod-lrt_numClusterK);
      }
      
      lchrombitarray_iter->setFitness(lT_j1);
      
    }

    aoop_outParamGAC.sumTotalInvalidOffspring
      (ll_invalidOffspring);

#ifdef __VERBOSE_YES 
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END 1. INITIALIZE POPULATION*/


  while ( (llfh_listFuntionHist.getDomainUpperBound() < 
	   aiinp_inParamGAPrototypes.getNumMaxGenerations()) &&
	  (runtime::elapsedTime(let_executionTime) <
	   aiinp_inParamGAPrototypes.getMaxExecutiontime() ) ) {
   
    llfh_listFuntionHist.increaseDomainUpperBound();

    
    {/*BEGIN 2. FORMING THE MATING SET M. IN THE CURRENT 
       IMPLEMENTATION M CONICIDES WITH PI
     */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep =
	"2. FORMING THE MATING SET M. IN THE CURRENT IMPLEMENTATION M CONICIDES WITH PI";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
    
      for (uintidx lui_i = 0; lui_i < lvectorchrom_population.size(); ++lui_i) {
	*lvectorchrom_matingPool[lui_i] = *lvectorchrom_population[lui_i]; 
      }
        
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    
    } /*END 2. FORMING THE MATING SET M. IN THE CURRENT 
	IMPLEMENTATION M CONICIDES WITH PI
      */

    {/*BEGIN 3. CROSSOVER*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "CROSSOVER";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      gaiterator::crossoverRandSelect
	(lvectorchrom_matingPool.begin(),
	 lvectorchrom_matingPool.end(),
	 lvectorchrom_setO.begin(),
	 lvectorchrom_setO.end(),
	 [&](gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aichrom_parent1,
	     gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aichrom_parent2,
	     gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aochrom_child1, 
	     gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aochrom_child2
	     )
	 {

	   gabinaryop::uniformCrossover
	     (*aochrom_child1,
	      *aochrom_child2,
	      *aichrom_parent1,
	      *aichrom_parent2,
	      aiinp_inParamGAPrototypes.getProbCrossover()
	      );
	  
	   aochrom_child1->setObjetiveFunc(std::numeric_limits<double>::max());
	   aochrom_child2->setObjetiveFunc(std::numeric_limits<double>::max());
	   
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

    } /*END 3. CROSSOVER*/

    { /*BEGIN MUTATION EACH BIT OF EACH OFFSPRING CHROMOSOME ALTERNATES (MUTATES)*/
    
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto lchrombitarray_iter: lvectorchrom_setO) {
      
	gabinaryop::eachBitArrayMutation
	  (*lchrombitarray_iter,
	   aiinp_inParamGAPrototypes.getProbMutation()
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
    
    } /*END MUTATION EACH BIT OF EACH OFFSPRING CHROMOSOME 
	ALTERNATES (MUTATES)
      */

    { /*BEGIN ALL ELEMENTS OF O ARE THEN EVALUATED BY THE 
	FITNESS FUNCTION
      */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "EVALUATED ALL ELEMENTS OF O";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      long ll_invalidOffspring = 0;
    
      for ( auto lchrombitarray_iter: lvectorchrom_setO) {
      
	std::list<uintidx>&& listui_idxInst =
	  lchrombitarray_iter->getIdxWithBitOn();

	T_REAL lT_j1 =  std::numeric_limits<T_REAL>::max();

	if ( luintidx_numClusterK != (uintidx) listui_idxInst.size() ) 
	  ++ll_invalidOffspring;
      
	if ( listui_idxInst.size() > 1 ) {
	
	  mat::MatrixRow<T_FEATURE> 
	    lmatrixt_v                         
	    ( (uintidx) listui_idxInst.size(),
	      data::Instance<T_FEATURE>::getNumDimensions()
	      );

	  clusteringop::initialize
	    (lmatrixt_v,
	     aiiterator_instfirst,
	     listui_idxInst.begin()
	     );

	  mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>
	    lbcrispmatrix_w(lmatrixt_v.getNumRows(),luintidx_numIntances);

	  clusteringop::getPartition
	    (lbcrispmatrix_w,
	     lmatrixt_y,
	     lmatrixt_v,
	     aifunc2p_dist
	     );

	  lT_j1 = 
	    um::j1
	    (lbcrispmatrix_w,
	     lmatrixt_v, 
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     aifunc2p_dist
	     );

	  T_REAL lrt_pmod = (T_REAL) lmatrixt_v.getNumRows();
   
	  lT_j1 +=  aiinp_inParamGAPrototypes.getAlpha()
	    *(lrt_pmod-lrt_numClusterK)*(lrt_pmod-lrt_numClusterK);
	}
      
	lchrombitarray_iter->setFitness(lT_j1);
      
      }

      aoop_outParamGAC.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    
    } /*END ALL ELEMENTS OF O ARE THEN EVALUATED BY 
	THE FITNESS FUNCTION
      */
  
    {/*BEGIN 5. RECOMBINATION PI AND O ARE POOLED AND 
       THE BEST Npop INDIVIDUALS SURVIVE SORT
     */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "RECOMBINATION PI AND O";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* >
	lvectorchrom_tmpPooled;

      lvectorchrom_tmpPooled.reserve
	( lvectorchrom_matingPool.size() + lvectorchrom_setO.size() );

      for ( auto lchrombitarray_iter: lvectorchrom_matingPool) {
      	lvectorchrom_tmpPooled.push_back
	  (lchrombitarray_iter);
      }
      for ( auto lchrombitarray_iter: lvectorchrom_setO) {
      	lvectorchrom_tmpPooled.push_back
	  (lchrombitarray_iter);
      }

      std::sort
	(lvectorchrom_tmpPooled.begin(),
	 lvectorchrom_tmpPooled.end(),
	 [](const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* x, 
	    const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* y
	    ) 
	 {  return x->getFitness() < y->getFitness(); } 
	 );

      for (uintidx lui_i = 0; lui_i < lvectorchrom_population.size(); ++lui_i) {
	*lvectorchrom_population[lui_i] = *lvectorchrom_tmpPooled[lui_i];


	if ( lochrombitarray_best.getFitness() > lvectorchrom_population[0]->getFitness() ) {

#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "COPY BEST";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep  
		      << ":  IN(" << geiinparam_verbose << ')'
		      << std::endl;
	  }
#endif /*__VERBOSE_YES*/
      
	  lochrombitarray_best = *lvectorchrom_population[0];
	
	  /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	  aoop_outParamGAC.setIterationGetsBest
	    (llfh_listFuntionHist.getDomainUpperBound());
	  aoop_outParamGAC.setRunTimeGetsBest
	    (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << geverbosepc_labelstep
		      << ": OUT(" << geiinparam_verbose << ')'
		      << std::endl;
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
	}

	
#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back
	  (lvectorchrom_population[lui_i]->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

      }

      lvectorchrom_tmpPooled.clear();
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END RECOMBINATION PI AND O ARE POOLED AND THE BEST 
	Npop INDIVIDUALS SURVIVE SORT
      */

      /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL 
	AND METRIC OF THE ALGORITHM
      */
#ifndef __WITHOUT_PLOT_STAT  
    if ( aiinp_inParamGAPrototypes.getWithPlotStatObjetiveFunc() ) {  
      lofh_J1->setValue(lvectorchrom_population[0]->getFitness());
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
      std::cout << "END ITERATION: "   << llfh_listFuntionHist.getDomainUpperBound()
		<< "\tobjetivoFunc = " << lvectorchrom_population[0]->getFitness()
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END While*/

  
  lochrombitarray_best = *lvectorchrom_population[0];

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

    for (uintidx lui_i = 0; lui_i < lvectorchrom_population.size(); ++lui_i) {
      delete lvectorchrom_population[lui_i];
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
    
    for (uintidx lui_i = 0; lui_i < lvectorchrom_matingPool.size(); ++lui_i) {
      delete lvectorchrom_matingPool[lui_i];
      delete lvectorchrom_setO[lui_i];
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
  aoop_outParamGAC.setNumClusterK
    (T_CLUSTERIDX(lochrombitarray_best.getNumBitOn()));
  aoop_outParamGAC.setMetricFuncRun
    (lochrombitarray_best.getFitness());
  aoop_outParamGAC.setFitness
    (lochrombitarray_best.getFitness());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
  

#ifdef __VERBOSE_YES
  geverbosepc_labelstep = lpc_labelAlgGA;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA 
	      << " OUT(" << geiinparam_verbose << ")\n";
    lochrombitarray_best.print();
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lochrombitarray_best;

} /*END algGA_BezdekEtAl1994 */

} /*END eac */

#endif /*__GAPROTOTYPES_FKMEDOID_HPP__*/

