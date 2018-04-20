/*! \file gaclustering_fkcrispmatrix.hpp
 * 
 * \brief GA CLUSTERING \cite Bezdek:etal:GAclustering:GA:1994
 * 
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GA algorithm based on the paper:\n
 * J.C. Bezdek, S. Boggavarapu, L.O. Hall, and A. Bensaid.\n 
 * Genetic algorithm guided clustering. In Evolutionary Computation,\n
 * 1994. IEEE World Congress on Computational Intelligence., Proceed-\n
 * ings of the First IEEE Conference on, pages 34–39 vol.1, Jun 1994.\n
 * <a href="http://dx.doi.org/10.1109/ICEC.1994.350046">doi:10.1109/ICEC.1994.350046</a>\n.
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

#ifndef __GACLUSTERING_FKCRISPMATRIX_HPP__
#define __GACLUSTERING_FKCRISPMATRIX_HPP__

#include <iostream>
#include <iomanip>
#include <vector>

#include <leac.hpp>

#include "plot_runtime_function.hpp"
#include "inparam_withoutpcpmfk.hpp"
#include "outparam_gac.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace eac {

/*! \fn gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL> gaclustering_fkcrispmatrix(inout::OutParamGAC<T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamWithoutPcPm<T_CLUSTERIDX,T_BITSIZE,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinpkbezdekga_inParam, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief  gaclustering_fkcrispmatrix 
  \details GA clustering based on \cite Bezdek:etal:GAclustering:GA:1994. 
  \returns A crisp matrix, which encodes a partition of a data set, for a defined k. 
  \param aoop_outParamGAC a inout::OutParamGAC that contains information relevant to program execution
  \param aiinpkbezdekga_inParam a inout::InParamWithoutPcPm with the input parameters for the program configuration  
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_BITSIZE,
           typename T_REAL,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, //0, 1, .., N
	   typename T_CLUSTERIDX,          //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL> 
gaclustering_fkcrispmatrix
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                &aoop_outParamGAC,
 inout::InParamWithoutPcPmFk
 <T_CLUSTERIDX,
 T_BITSIZE,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>       &aiinp_inParamWithoutPcPmFk,
 const INPUT_ITERATOR         aiiterator_instfirst,
 const INPUT_ITERATOR         aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  /*ID PROC
   */
  geverboseui_idproc = 1;

  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gaclustering_fkcrispmatrix";  
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA
      << "  IN(" << geiinparam_verbose << ")\n"
      << "\t(output inout::OutParamGAC&: aoop_outParamGAC[" 
      << &aoop_outParamGAC << "]\n"
      << "\t input  InParamWithoutPcPm&: aiinp_inParamWithoutPcPmFk[" 
      << &aiinp_inParamWithoutPcPmFk << "]\n"
      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
      << "\t input  dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist[" 
      << &aifunc2p_dist << ']'
      << "\n\t\tPopulation size = " 
      << aiinp_inParamWithoutPcPmFk.getSizePopulation()
      << "\n\t\tMatingPool size = " 
      << aiinp_inParamWithoutPcPmFk.getSizeMatingPool()
      << "\n\t\t Generations  = "
      << aiinp_inParamWithoutPcPmFk.getNumMaxGenerations()
      << "\n\t\trandom-seed = "
      << aiinp_inParamWithoutPcPmFk.getRandomSeed()
      << "\n\t)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  const uintidx  luintidx_numClusterK =
    (uintidx) aiinp_inParamWithoutPcPmFk.getNumClusterK();
  const uintidx  luintidx_numIntances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  /*CONVERT INSTANCES TO FORMAT MATRIX
   */
  mat::MatrixRow<T_FEATURE>&& lmatrixt_y =
    data::toMatrixRow
    (aiiterator_instfirst,
     aiiterator_instlast
     );

  std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_mmcidx0K
    (0,aiinp_inParamWithoutPcPmFk.getNumClusterK()-1);

  gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>
    lochrombitcrispmatrix_best(luintidx_numClusterK,luintidx_numIntances);
  
  /*STL container for storing the chromosome population
   */
  std::vector<gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* >  
    lvectorchrombitcrispmatrix_population;

  /*Vector for matingpool
   */
  std::vector<gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* >
    lvectorchrombitcrispmatrix_matingPool;

  /*Vector for temporary storage when applying generic operators
   */
  std::vector<gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* >
    lvectorchromfixleng_childR;
  
  if ( aiinp_inParamWithoutPcPmFk.getSizePopulation()
       <= aiinp_inParamWithoutPcPmFk.getSizeMatingPool() )
    throw std::invalid_argument
      ("gaclustering_fkcrispmatrix: "
       "size population should be greater than size matingpool"
       );
  
  runtime::ListRuntimeFunction<COMMON_IDOMAIN>
    llfh_listFuntionHist
    (aiinp_inParamWithoutPcPmFk.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

  /*Declaration of variables: computing statistical 
    and metric of the algorithm
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream               lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_J1 = NULL;
  runtime::RuntimeFunctionValue<T_INSTANCES_CLUSTER_K>
    *lofh_misclassified = NULL; /*function extra*/
  runtime::RuntimeFunctionStat<T_REAL>
    *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;
  
  if ( aiinp_inParamWithoutPcPmFk.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      ( aiinp_inParamWithoutPcPmFk.getSizePopulation());
    //Variable to monitor in the execution of the program
    lofh_J1  = new runtime::RuntimeFunctionValue<T_REAL>
      ("J1", 
       aiinp_inParamWithoutPcPmFk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_J1);

    if ( aiinp_inParamWithoutPcPmFk.getClassInstanceColumn()  ) {
      lofh_misclassified =
	new runtime::RuntimeFunctionValue<T_INSTANCES_CLUSTER_K>
	("Misclassified", 
	 aiinp_inParamWithoutPcPmFk.getAlgorithmoName(),
	 RUNTIMEFUNCTION_NOT_STORAGE
	 );
      llfh_listFuntionHist.addFuntion(lofh_misclassified);
    }

    //Statistics of variable J1 in runtime
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamWithoutPcPmFk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamWithoutPcPmFk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamWithoutPcPmFk.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open
      (aoop_outParamGAC.getFileNameOutPlotStatObjetiveFunc().c_str(),
       std::ios::out | std::ios::app
       );

    lfileout_plotStatObjetiveFunc.precision(COMMON_COUT_PRECISION);

    //Header function
    lfileout_plotStatObjetiveFunc 
      <<  llfh_listFuntionHist.getHeaderFuntions() 
      << "\n";
  }
#endif /*__WITHOUT_PLOT_STAT*/

  runtime::ExecutionTime let_executionTime  = runtime::start();

  /*Create space for store population
   */
  lvectorchrombitcrispmatrix_population.reserve
    (aiinp_inParamWithoutPcPmFk.getSizePopulation() + 1);
  for (uintidx lui_i = 0; 
       lui_i < aiinp_inParamWithoutPcPmFk.getSizePopulation(); 
       lui_i++) 
    {
      lvectorchrombitcrispmatrix_population.push_back
	(new gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>
	 (luintidx_numClusterK,luintidx_numIntances)
	 );
    }

  /*Space for store matingpool
   */
  lvectorchrombitcrispmatrix_matingPool.reserve
    (aiinp_inParamWithoutPcPmFk.getSizeMatingPool());

  /*Space for chromosomes R  
   */
  lvectorchromfixleng_childR.reserve
    (aiinp_inParamWithoutPcPmFk.getSizeMatingPool() + 1 );

  /*Initialization of population

    Initial population of size P, consisting of U matrices is pseudo randomly 
    generated such that each has one at least one 1 in every row 
    (\sum_{j=1}^{n} U_{ij} >= 1 \all_i) and each column sums to 1, i.e. 
    \sum_{i=1}^{c} U_{ij} = 1,\forall j. 

    The partly random initialization is obtained as follows. For each cluster 
    center v_i, we choose the kth element of the cluster center to be the kth 
    feature of a randomly chosen pattern to be clustered. This is done for each 
    of the s elements of a cluster center. The process is repeated for each 
    cluster center. An initial U matrix is then generated from the cluster 
    centers. For a GA, population P (the population size) U matrices are 
    generated in this manner.
  */
  {/*BEGIN INITIALIZE POPULATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "POPULATION INITIALIZATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep  
	<< ":  IN(" << geiinparam_verbose << ')'
	<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    mat::MatrixRow<T_FEATURE> 
      lmatrixt_v                         
      ( luintidx_numClusterK, 
	data::Instance<T_FEATURE>::getNumDimensions()
	);

    for ( auto lchrombitcrispmatrix_iter: lvectorchrombitcrispmatrix_population) {
     	 
      clusteringop::randomInitialize
	(lmatrixt_v,
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );

      clusteringop::getPartition
	(*lchrombitcrispmatrix_iter,
	 lmatrixt_y,
	 lmatrixt_v,
	 aifunc2p_dist
	 );

      T_REAL lT_j1 = 
	um::j1
	(*lchrombitcrispmatrix_iter, 
	 lmatrixt_v, 
	 aiiterator_instfirst,
	 aiiterator_instlast, 
	 aifunc2p_dist
	 );

      lchrombitcrispmatrix_iter->setObjetiveFunc(lT_j1);
	 
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

  } /*END INITIALIZE POPULATION*/

  /*Population sort by J_1
 
    The U matrices are sorted by J_1 value. and the R with the lowest J_1, 
    values are choses to reproduce. 
  */
  {/*BEGIN POPULATION SORT BY J_1*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "SORT POPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep  
	<< ":  IN(" << geiinparam_verbose << ')'
	<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    
    std::sort
      (lvectorchrombitcrispmatrix_population.begin(),
       lvectorchrombitcrispmatrix_population.end(),
       [](const gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* x, 
	  const gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* y
	  ) 
       {  return x->getObjetiveFunc() < y->getObjetiveFunc(); } 
       );

	   
#ifdef __VERBOSE_YES
      
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {

      for ( auto lchrombitcrispmatrix_iter: lvectorchrombitcrispmatrix_population) {
				      
	lchrombitcrispmatrix_iter->print
	  (std::cout,
	   geverbosepc_labelstep,
	   ',',
	   ';'
	   );
	std::cout << '\n';
      }
    }   
    --geiinparam_verbose;

    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep
	<< ": OUT(" << geiinparam_verbose << ')'
	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END POPULATION SORT BY J_1*/


  while( true ) {

    {/*BEGIN PRESERVING THE CHROMOSOME BEST
      */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ":  IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      if ( lvectorchrombitcrispmatrix_population[0]->getObjetiveFunc()
	   < lochrombitcrispmatrix_best.getObjetiveFunc() ) {
	lochrombitcrispmatrix_best =
	  *lvectorchrombitcrispmatrix_population[0];
	 /*A better chromosome is found in this iteration
	  */
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));
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

    } /*END PRESERVING THE CHROMOSOME BEST*/

    /*COMPUTING STATISTICAL OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT

    if ( aiinp_inParamWithoutPcPmFk.getWithPlotStatObjetiveFunc() ) {  

      for ( auto lchrombitcrispmatrix_iter:
	      lvectorchrombitcrispmatrix_population) {
	lvectorT_statfuncObjetiveFunc.push_back
	  (lchrombitcrispmatrix_iter->getObjetiveFunc());
      }

      lofh_J1->setValue
	(lvectorchrombitcrispmatrix_population[0]->getObjetiveFunc());

      if ( lofh_misclassified != NULL ) {

	partition::PartitionCrispMatrix
	  <T_BITSIZE,T_CLUSTERIDX>
	  lpartitionCrispMatrix_classifierU
	  (*lvectorchrombitcrispmatrix_population[0]);

	sm::ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>&&
	  lmatchmatrix_confusion =
	  sm::getConfusionMatrix
	  (aiiterator_instfirst,
	   aiiterator_instlast,
	   lpartitionCrispMatrix_classifierU,
	   [](const data::Instance<T_FEATURE>* aiinst_iter )
	   -> T_INSTANCES_CLUSTER_K
	   {
	     return  T_INSTANCES_CLUSTER_K(1);
	   },
	   [](const data::Instance<T_FEATURE>* aiinst_iter )
	   -> T_CLUSTERIDX
	   {
	     data::InstanceClass
	       <T_FEATURE,
		T_INSTANCES_CLUSTER_K,
		T_CLUSTERIDX>
	       *linstclass_iter = 
	       (data::InstanceClass
		<T_FEATURE,
		T_INSTANCES_CLUSTER_K,
		T_CLUSTERIDX>*)
	       aiinst_iter;
	       
	     return linstclass_iter->getClassIdx();
	       
	   }
	   );
	lofh_misclassified->setValue
	  (lmatchmatrix_confusion.getMisclassified());
      
      }
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
	<< lochrombitcrispmatrix_best.getObjetiveFunc() 
	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    /*Termination criterion attained?
     */
    if ( (llfh_listFuntionHist.getDomainUpperBound()
	  >= aiinp_inParamWithoutPcPmFk.getNumMaxGenerations()) ||
	 (runtime::elapsedTime(let_executionTime) >
	  aiinp_inParamWithoutPcPmFk.getMaxExecutiontime())
	 )
      break;
  
    /*Selection
      R matrices with the lowest J_1, values are choses to reproduce. 
     */
    
    {/*BEGIN SELECTION
      */
      auto ichrom_population = lvectorchrombitcrispmatrix_population.begin();
   
      for (uintidx lui_i = 0; 
	   lui_i < aiinp_inParamWithoutPcPmFk.getSizeMatingPool(); 
	   lui_i++)  {
	lvectorchrombitcrispmatrix_matingPool.push_back
	  (*ichrom_population);
	++ichrom_population;
      }
      
    }/*END SELECTION*/

    /*Crossover operator
      
      The crossover point and number of columns in the two U matrices
      chosen for reproduction are randomly chosen. The columns of the 
      matrices are combined t o create the children matrices.
    */
    {/*BEGIN CROSSOVER OPERATORS*/
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "CROSSOVER OPERATORS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ":  IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for (uintidx lui_i = 0; 
	   lui_i < aiinp_inParamWithoutPcPmFk.getSizeMatingPool(); 
	   lui_i++) {
	lvectorchromfixleng_childR.push_back
	  (new gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>
	   (luintidx_numClusterK,luintidx_numIntances)
	   );
      }

      gaiterator::crossoverRandSelect
      (lvectorchrombitcrispmatrix_matingPool.begin(),
       lvectorchrombitcrispmatrix_matingPool.end(),
       lvectorchromfixleng_childR.begin(),
       lvectorchromfixleng_childR.end(),
       [&](gaencode::ChromosomeCrispMatrix
	   <T_BITSIZE,T_CLUSTERIDX,T_REAL>* aichrom_parent1,
	   gaencode::ChromosomeCrispMatrix
	   <T_BITSIZE,T_CLUSTERIDX,T_REAL>* aichrom_parent2,
	   gaencode::ChromosomeCrispMatrix
	   <T_BITSIZE,T_CLUSTERIDX,T_REAL>* aochrom_child1, 
	   gaencode::ChromosomeCrispMatrix
	   <T_BITSIZE,T_CLUSTERIDX,T_REAL>* aochrom_child2
	   )
       {
	 
	 gabinaryop::onePointDistCrossover
	   (*aochrom_child1,
	    *aochrom_child2,
	    *aichrom_parent1, 
	    *aichrom_parent2
	    );
	 
	 aochrom_child1->setObjetiveFunc(std::numeric_limits<T_REAL>::max());
	 aochrom_child2->setObjetiveFunc(std::numeric_limits<T_REAL>::max());
	   
       }
       );

      lvectorchrombitcrispmatrix_matingPool.clear();
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    }/*END CROSSOVER OPERATORS*/


    /*Mutation consists of randomly choosing an element of a
      column to have the value 1, such that it is a different 
      element than the one currently having a value of 1.
    */
    
    {/*BEGIN MUTATION OPERATOR*/ 
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "MUTATION OPERATOR";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ":  IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto ichrom_childR: lvectorchromfixleng_childR ) {
	gabinaryop::bitMutation(*ichrom_childR);
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

    } /*END MUTATION OPERATOR*/   


    {/*BEGIN EVALUATE J1 FOR CHILDR*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "EVALUATE J1 FOR CHILDR"; 
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ":  IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      mat::MatrixRow<T_FEATURE> 
	lmatrixt_v 
	(luintidx_numClusterK, 
	 data::Instance<T_FEATURE>::getNumDimensions() 
	 );
      
      mat::MatrixRow<T_FEATURE_SUM>
	lmatrixT_sumWX
	(lmatrixt_v.getNumRows(),
	 lmatrixt_v.getNumColumns()
	 );
      std::vector<T_INSTANCES_CLUSTER_K>
	lvectorT_sumWik(lmatrixt_v.getNumRows());
      
      for ( auto ichrom_childR: lvectorchromfixleng_childR ) {

	 /*Calculate the centroid associated with U_i
	  */
	clusteringop::getCentroids
	  (lmatrixt_v,
	   lmatrixT_sumWX,
	   lvectorT_sumWik,
	   *ichrom_childR, 
	   lmatrixt_y
	   );
	   
	T_REAL lT_j1 = 
	  um::j1
	  (*ichrom_childR,
	   lmatrixt_v, 
	   aiiterator_instfirst,
	   aiiterator_instlast, 
	   aifunc2p_dist
	   );

	ichrom_childR->setObjetiveFunc(lT_j1);
	 
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

    } /*EVALUATE J1 FOR CHILDR*/

    /*The R chuild U matrices are added to the population 
      with the P-R U matrices with the greatest J1 values 
      dropped from the population. 
    */
    { /*BEGIN ADD P-R U MATRICES TO POPULATION*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ADD P-R U MATRICES TO POPULATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ":  IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      std::sort
	(lvectorchromfixleng_childR.begin(),
	 lvectorchromfixleng_childR.end(),
	 [](const gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* x, 
	    const gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* y
	    ) 
	 {  return x->getObjetiveFunc() < y->getObjetiveFunc(); } 
	 );
      
      std::vector<gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>* >  
	lvectorchrombitcrispmatrix_tmpL;

      lvectorchrombitcrispmatrix_tmpL.swap(lvectorchrombitcrispmatrix_population);

      /*Insert a sentinel to merge the two vectors
       */
      lvectorchrombitcrispmatrix_tmpL.push_back
	(new gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>());
      lvectorchromfixleng_childR.push_back
	(new gaencode::ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_REAL>());
      
      lvectorchrombitcrispmatrix_population.reserve
	(aiinp_inParamWithoutPcPmFk.getSizePopulation() + 1);

      uintidx luintidx_l = 0;
      uintidx luintidx_r = 0;
      
      for (uintidx lui_i = 0;
	   lui_i < aiinp_inParamWithoutPcPmFk.getSizePopulation();
	   lui_i++)
	{

	  if (  lvectorchrombitcrispmatrix_tmpL[luintidx_l]->getObjetiveFunc() <
		lvectorchromfixleng_childR[luintidx_r]->getObjetiveFunc() )
	    {
	 
#ifdef __VERBOSE_YES
	      ++geiinparam_verbose;
	      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
		std::cout
		  << " lvectorchrombitcrispmatrix_population[" << lui_i <<  ']' 
		  << " <--  lvectorchrombitcrispmatrix_tmpL[" << luintidx_l << ']' 
		  << '[' << & lvectorchrombitcrispmatrix_population[luintidx_l] << ']' 
		  << " Fitness: "
		  <<  lvectorchrombitcrispmatrix_tmpL[luintidx_l]->getObjetiveFunc() 
		  << '\n';
	      }
	      --geiinparam_verbose;
#endif //__VERBOSE_YES

	      lvectorchrombitcrispmatrix_population.push_back
		(lvectorchrombitcrispmatrix_tmpL[luintidx_l]);
	      lvectorchrombitcrispmatrix_tmpL[luintidx_l] = NULL;
	      ++luintidx_l;
	      
	    }
	  else {

#ifdef __VERBOSE_YES
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      std::cout
		<< " lvectorchrombitcrispmatrix_population[" << lui_i <<  ']' 
		<< " <--  lvectorchromfixleng_childR[" << luintidx_r << ']' 
		<< "[" << & lvectorchrombitcrispmatrix_population[luintidx_r] << ']' 
		<< " Fitness: "
		<<  lvectorchromfixleng_childR[luintidx_r]->getObjetiveFunc()
		<< '\n';
	    }
	    --geiinparam_verbose;
#endif //__VERBOSE_YES
	    
	    lvectorchrombitcrispmatrix_population.push_back
	      (lvectorchromfixleng_childR[luintidx_r]);
	    lvectorchromfixleng_childR[luintidx_r] = NULL;
	    ++luintidx_r;
	     
	  }
	
	}

      for (uintidx lui_i = 0;
	   lui_i < lvectorchromfixleng_childR.size();
	   ++lui_i) {
	if ( lvectorchromfixleng_childR[lui_i] != NULL )
	  delete lvectorchromfixleng_childR[lui_i];
      }
      lvectorchromfixleng_childR.clear();

      for (uintidx lui_i = 0;
	   lui_i < lvectorchrombitcrispmatrix_tmpL.size();
	   ++lui_i) {
	if ( lvectorchrombitcrispmatrix_tmpL[lui_i] != NULL )
	  delete lvectorchrombitcrispmatrix_tmpL[lui_i];
      }
      lvectorchrombitcrispmatrix_tmpL.clear();
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    } /*END ADD P-R U MATRICES TO POPULATION*/

    /*The reproduction and survival of fittest process 
      continues for some set number of generations
    */

    llfh_listFuntionHist.increaseDomainUpperBound();

  } /*while*/

  /*FREE MEMORY*/
  {/*BEGIN FREE MEMORY OF POPULATION*/ 
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "DELETEPOPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<<  geverbosepc_labelstep
	<< ":  IN(" << geiinparam_verbose << ')'
	<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    for (uintidx lui_i = 0;
	 lui_i < lvectorchrombitcrispmatrix_population.size();
	 ++lui_i) {
      delete lvectorchrombitcrispmatrix_population[lui_i];
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

  }/*END FREE MEMORY OF POPULATION*/

  runtime::stop(let_executionTime);
  aoop_outParamGAC.setNumClusterK
    (aiinp_inParamWithoutPcPmFk.getNumClusterK());
  aoop_outParamGAC.setMetricFuncRun
    (lochrombitcrispmatrix_best.getObjetiveFunc());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
 
  aoop_outParamGAC.setFitness
    (lochrombitcrispmatrix_best.getObjetiveFunc());
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamWithoutPcPmFk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamWithoutPcPmFk,
       aoop_outParamGAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

  

#ifdef __VERBOSE_YES
  geverbosepc_labelstep = lpc_labelAlgGA;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA 
      << " OUT(" << geiinparam_verbose << ")\n";
    std::setprecision(COMMON_COUT_PRECISION);
    
    mat::MatrixRow<T_FEATURE> 
      lmatrixt_vBestChrom
      ( luintidx_numClusterK, 
	data::Instance<T_FEATURE>::getNumDimensions() 
	);

    mat::MatrixRow<T_FEATURE_SUM>
      lmatrixT_sumWX
      (lmatrixt_vBestChrom.getNumRows(),
       lmatrixt_vBestChrom.getNumColumns()
       );
      
    std::vector<T_INSTANCES_CLUSTER_K>
      lvectorT_sumWik(lmatrixt_vBestChrom.getNumRows());
	  
    clusteringop::getCentroids
      (lmatrixt_vBestChrom,
       lmatrixT_sumWX,
       lvectorT_sumWik,
       lochrombitcrispmatrix_best, 
       lmatrixt_y
       );

    lochrombitcrispmatrix_best.print
      (std::cout,
       geverbosepc_labelstep,
       ',',
       ';'
       );
    
    std::cout << '\n';
    um::j1
      (lochrombitcrispmatrix_best, 
       lmatrixt_vBestChrom, 
       aiiterator_instfirst,
       aiiterator_instlast, 
       aifunc2p_dist
       );
   
    std::cout << std::endl;
    
    std::setprecision(COMMON_VERBOSE_COUT_PRECISION);
    
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lochrombitcrispmatrix_best;

} /*END gaclustering_fkcrispmatrix */

} /*END namespace alg*/

#endif /*__GACLUSTERING_FKCRISPMATRIX_HPP__*/

