/*! \file gaclustering_vktreebinary.hpp
 *
 * \brief GA Clustering \cite Casillas:etal:GAclusteringVarK:GA:2003
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GACLUSTERING algorithm based on the paper:\n
 * A. Casillas, M.T. Gonzalez de Lena, and R. Martinez. Document\n
 * clustering into an unknown number of clusters using a genetic\n
 * algorithm. In V ́aclav Matouˇsek and Pavel Mautner, editors, Text,\n 
 * Speech and Dialogue, volume 2807 of Lecture Notes in Computer\n
 * Science, pages 43–49. Springer Berlin Heidelberg, 2003.\n
 * <a href="http://dx.doi.org/10.1007/978-3-540-39398-6_7">doi:http://dx.doi.org/10.1007/978-3-540-39398-6_7</a>\n.
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
 

#ifndef __GACLUSTERING_VKTREEBINARY_HPP__
#define __GACLUSTERING_VKTREEBINARY_HPP__

#include <vector>
#include <algorithm>   // std::iter_swap
#include <utility>      // std::pair
#include <leac.hpp>
#include "inparam_gaclustering_vktreebinary.hpp"
#include "outparam_gaclustering.hpp"
#include "container_out.hpp"
#include "plot_runtime_function.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of genetic and evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace eac {

/*! \fn std::pair<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>,std::vector<uintidx> > gaclustering_vktreebinary(inout::OutParamGAClustering <T_REAL, T_CLUSTERIDX> &aoopcga_outParamClusteringGA, inout::InParamGAClusteringVKTreeBinary<T_BITSIZE, T_CLUSTERIDX, T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinpcga_inParamGAClustering, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief GA Clustering \cite Casillas:etal:GAclusteringVarK:GA:2003
  \details Implementation of GA algorithm based on \cite Casillas:etal:GAclusteringVarK:GA:2003. 
 \returns A partition of a data set, encoded on a chromosome with n − 1 binary genes and a minimum spanning Tree (MST) where each gene with value “0” means that this  edge remains and “1” means that this edge is eliminated. The number of elements with value “1” represents the value of k − 1.
  \param aoopcga_outParamClusteringGA a OutParamClusteringGA<T_REAL,T_CLUSTERIDX>
  \param aiinpcgaprobfixedk_inParamGA a inparam::InParamGAClusteringProbCProbMFixedK<T_CLUSTERIDX,T_REAL>
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
 \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist a dist::Dist<T_REAL,T_FEATURE>
*/
template < typename T_BITSIZE,
	   typename T_REAL,
           typename T_FEATURE,     //T_STRING,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
std::pair<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>,
	  std::vector<uintidx> >
gaclustering_vktreebinary
(inout::OutParamGAClustering
 <T_REAL,
 T_CLUSTERIDX>                          &aoopcga_outParamClusteringGA,
 inout::InParamGAClusteringVKTreeBinary
 <T_BITSIZE,
 T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                 &aiinpcga_inParamGAClustering,
 const INPUT_ITERATOR                   aiiterator_instfirst,
 const INPUT_ITERATOR                   aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>           &aifunc2p_dist
 )
{
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  if ( aiinpcga_inParamGAClustering.getNumMaxGenerations() == 0 )  {
    aiinpcga_inParamGAClustering.setNumMaxGenerations
      ( (COMMON_IDOMAIN) lui_numInstances );
  }

  if ( aiinpcga_inParamGAClustering.getSizePopulation() == 0 ) {
    aiinpcga_inParamGAClustering.setSizePopulation
      ( uintidx( 10 * lui_numInstances ) );
  }
    
#ifdef _INITIATES_KMIN_KMAX_POPULATION_
  if ( aiinpcga_inParamGAClustering.getNumClusterKMaximum() == 
       INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED )
    aiinpcga_inParamGAClustering.setNumClusterKMaximum
      (T_CLUSTERIDX(((lui_numInstances -1) /2) + 1));
#endif // _INITIATES_KMIN_KMAX_POPULATION_

 
 
#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gaclustering_vktreebinary";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output OutParamClusteringGA&: aoopcga_outParamClusteringGA[" 
	      << &aoopcga_outParamClusteringGA << "]\n"
	      << "\t input  InParamGAClusteringVKTreeBinary&: aiinpcga_inParamGAClustering[" 
	      << &aiinpcga_inParamGAClustering << "]\n"
	      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t\tPopulation size = " 
	      << aiinpcga_inParamGAClustering.getSizePopulation()
	      << "\n\t\tProbCrossover = " 
	      << aiinpcga_inParamGAClustering.getProbCrossover() 
	      << "\n\t\tProbMutation  = " 
	      << aiinpcga_inParamGAClustering.getProbMutation()
	      << "\n\t\tNumMaxGenerations = "
	      << aiinpcga_inParamGAClustering.getNumMaxGenerations()
#ifdef _INITIATES_KMIN_KMAX_POPULATION_  
	      << "\n\t\tk-minimum  = " 
	      << aiinpcga_inParamGAClustering.getNumClusterKMinimum()
	      << "\n\t\tk-maximum  = " 
	      << aiinpcga_inParamGAClustering.getNumClusterKMaximum()
#endif //_INITIATES_RANDOM_POPULATION_      
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 

  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   */
  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* >
    lvectorchrom_population;

  std::vector<gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* >
    lvectorchrom_stringPool;
   
  gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    lochrom_best(lui_numInstances-1 );

  std::uniform_int_distribution<int>     uniformdis_01(0,1);
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
   
  runtime::ListRuntimeFunction<COMMON_IDOMAIN>
    llfh_listFuntionHist
    (aiinpcga_inParamGAClustering.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );
  
  COMMON_IDOMAIN  lit_iterNotChange =
    aiinpcga_inParamGAClustering.getNumNotChangeStop();
  //GA_CLUSTERINGVARK_NOT_CHANGE_STOP;

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream               lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_VRC = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;
  
  if ( aiinpcga_inParamGAClustering.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      ( aiinpcga_inParamGAClustering.getSizePopulation());
    //DEFINE FUNCTION
    lofh_VRC  = new runtime::RuntimeFunctionValue<T_REAL>
      ("VRC", 
       aiinpcga_inParamGAClustering.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_VRC);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinpcga_inParamGAClustering.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoopcga_outParamClusteringGA.setFileNameOutPlotStatObjetiveFunc
      (aiinpcga_inParamGAClustering.getFileNamePlotStatObjetiveFunc(),
       aiinpcga_inParamGAClustering.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open
      (aoopcga_outParamClusteringGA.getFileNameOutPlotStatObjetiveFunc().c_str(),
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
  aoopcga_outParamClusteringGA.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/

  runtime::ExecutionTime let_executionTime = runtime::start();
    
  /*An n × n distance tringular  matrix is then calculated. 
   */

  mat::MatrixTriang<T_REAL>&& 
    lmatrixtriagT_dissimilarity =
    dist::getMatrixDissimilarity
    (aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  std::vector<uintidx> lvectorinstidx_piMST;
  
  /*Next, the method needs to calculate the 
    Minimum Spanning Tree (MST), so that the 
    enormous number of possible partitions of 
    a set of points is reduced to those which 
    are obtainable by splitting the MST
  */
  {/*BEGIN CALCULATE MINIMUM SPANNING TREE (MST)*/

#ifdef __VERBOSE_YES
    const char *lpc_labelStep = "CALCULATE MINIMUM SPANNING TREE";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelStep;
      std::cout << ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    lvectorinstidx_piMST =
      graph::prim(lmatrixtriagT_dissimilarity);

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelStep
		<< ": OUT(" << geiinparam_verbose << ")\n";
	  
      std::vector<std::list<uintidx> >&& lvectorlist_tmpgraphpi =
	graph::fromPiToGraph(lvectorinstidx_piMST);
      std::ostringstream lostrstream_labelGraphPi;
      lostrstream_labelGraphPi << "<GRAPH: " << lpc_labelStep
			       << ":lvectorlist_tmpgraphpi["
			       << &lvectorlist_tmpgraphpi << ']';

      inout::vectorlistprint
	(lvectorlist_tmpgraphpi,
	 std::cout,
	 lostrstream_labelGraphPi.str().c_str(),
	 ';'
	 );
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
  } /*END CALCULATE MINIMUM SPANNING TREE (MST)*/
  
  /*3.1 Population Representation
    A vector with n − 1 binary elements, value “0” means that this 
    edge remains and “1” means that this edge is eliminated. The 
    number of elements with value “1” represents the value of k − 1.
  */
  
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
      ( aiinpcga_inParamGAClustering.getSizePopulation() );

    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinpcga_inParamGAClustering.getSizePopulation(); 
	 luintidx_i++) 
      {	
	lvectorchrom_population.push_back
	  (new gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
	   (lui_numInstances-1)
	   );
      }

#ifdef _INITIATES_KMIN_KMAX_POPULATION_
    
    std::uniform_int_distribution<T_CLUSTERIDX>
      uniformdis_iMinMaxK
      (aiinpcga_inParamGAClustering.getNumClusterKMinimum(), 
       aiinpcga_inParamGAClustering.getNumClusterKMaximum()
       );

    std::uniform_int_distribution<uintidx> uniformdis_uiInitialize0N_1
      (0,(lui_numInstances-2));
  
    for (auto&& lchrom_iter: lvectorchrom_population ) {
      T_CLUSTERIDX li_krand = uniformdis_iMinMaxK(gmt19937_eng);
      --li_krand;
      lchrom_iter->initialize();
      for ( T_CLUSTERIDX _li_k = 1; _li_k <= li_krand; ++_li_k) {
	lchrom_iter->setBit(uniformdis_uiInitialize0N_1(gmt19937_eng));
      }
    }
    
#endif //_INITIATES_RANDOM_POPULATION_
    
#ifdef _INITIATES_RANDOM_POPULATION_
    for (auto&& lchrom_iter: lvectorchrom_population ) {
    
      gabinaryop::initializeGenes
	(*lchrom_iter,
	 [&]() 
	 {
	   return uniformdis_01(gmt19937_eng);
	 }
	 );
      
    }
#endif //_INITIATES_RANDOM_POPULATION_

    /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
     */
    lvectorchrom_stringPool.reserve
      ( aiinpcga_inParamGAClustering.getSizePopulation() );

    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinpcga_inParamGAClustering.getSizePopulation(); 
	 luintidx_i++) 
      {	
	lvectorchrom_stringPool.push_back
	  (new gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
	   (lui_numInstances-1));
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


  llfh_listFuntionHist.increaseDomainUpperBound();
  
  while ( 1 ) {

    /*2. Fitness computation
      Based on the Calinski and Harabasz Stopping Rule.
    */

    {/*BEGIN FITNESS COMPUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "FITNESS COMPUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      long ll_invalidOffspring = 0;
      for (auto&& lchrom_iter: lvectorchrom_population ) {
	
	partition::PartitionDisjSets
	  <T_CLUSTERIDX>
	  lpartitionDisjSets_clusters
	  (graph::component
	   (lvectorinstidx_piMST,
	    *lchrom_iter
	    )
	   );
		
	uintidx lui_numClusterK =
	  uintidx(lpartitionDisjSets_clusters.getNumCluster());

	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroids
	  (lui_numClusterK,
	   data::Instance<T_FEATURE>::getNumDimensions() 
	   );

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
    
	T_CLUSTERIDX lmcidx_numClusterNull =
	  clusteringop::getCentroids
	  (lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   lpartitionDisjSets_clusters,
	   aiiterator_instfirst,
	   aiiterator_instlast
	   );
	
        T_REAL lT_VRC;
	
	if ( lmcidx_numClusterNull == 0 ) {

	  lT_VRC =
	    um::VRC
	    (lmatrixrowt_centroids,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     lpartitionDisjSets_clusters,
	     aifunc2p_dist
	     );

	  if (lT_VRC != measuare_undefVRC(T_REAL) ) {
	    lchrom_iter->setValidString(true);
	  }
	  else {
	    lchrom_iter->setValidString(false);
	    aoopcga_outParamClusteringGA.incTotalInvalidOffspring();
	  }

	}
	else {
	  lT_VRC = 
	    measuare_undefVRC(DATATYPE_REAL);
	  lchrom_iter->setValidString(false);
	}
	  
	lchrom_iter->setObjetiveFunc(lT_VRC); 
	lchrom_iter->setFitness(lT_VRC);
	    
#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(lchrom_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
	   
      }
     
      aoopcga_outParamClusteringGA.sumTotalInvalidOffspring
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
	 [](const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* x, 
	    const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* y
	    ) 
	 {  return x->getFitness() < y->getFitness(); }
	 );

      if ( lochrom_best.getFitness() < (*lit_chromMax)->getFitness() ) {
	
	lochrom_best = *(*lit_chromMax);
	
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION
	 */
	lit_iterNotChange =
	  aiinpcga_inParamGAClustering.getNumNotChangeStop();
	aoopcga_outParamClusteringGA.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoopcga_outParamClusteringGA.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES   
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {

	  partition::PartitionDisjSets
	    <T_CLUSTERIDX>
	    lpartitionDisjSets_clustersChromBest
	    (graph::component
	     (lvectorinstidx_piMST,
	      lochrom_best
	      )
	     );

	  uintidx lui_chromBestNumClusterK =
	    uintidx(lpartitionDisjSets_clustersChromBest.getNumCluster());

	  mat::MatrixRow<T_FEATURE> 
	    lmatrixrowt_centroids
	    (lui_chromBestNumClusterK,
	     data::Instance<T_FEATURE>::getNumDimensions() 
	     );

	  mat::MatrixRow<T_FEATURE_SUM>       
	    lmatrixrowt_sumInstCluster
	    (lui_chromBestNumClusterK, 
	     data::Instance<T_FEATURE>::getNumDimensions(),
	     T_FEATURE_SUM(0)
	     );
	
	  std::vector<T_INSTANCES_CLUSTER_K> 
	    lvectort_numInstClusterK
	    (lui_chromBestNumClusterK,
	     T_INSTANCES_CLUSTER_K(0)
	     );
    
	  //T_CLUSTERIDX lmcidx_numClusterNull =
	  clusteringop::getCentroids
	    (lmatrixrowt_centroids,
	     lmatrixrowt_sumInstCluster,
	     lvectort_numInstClusterK,
	     lpartitionDisjSets_clustersChromBest,
	     aiiterator_instfirst,
	     aiiterator_instlast
	     );

	  std::ostringstream lostrstream_labelCentroids;
	  lostrstream_labelCentroids << "<CENTROIDSCLUSTER:";
	  lmatrixrowt_centroids.print
	    (std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
	  std::cout << std::endl;

	  std::ostringstream lostrstream_labelShipBest;
	  lostrstream_labelShipBest << "<MEMBERCLUSTERDISJSETS:";
	  lpartitionDisjSets_clustersChromBest.print
	    (std::cout,lostrstream_labelShipBest.str().c_str(),',');
	  std::cout << std::endl;

	  std::vector<std::list<uintidx> >&& lvectorlist_tmpgraphpi = 
	    graph::fromPiToGraph<uintidx,T_BITSIZE>(lvectorinstidx_piMST,lochrom_best);
	  std::ostringstream lostrstream_labelGraphPi;
	  lostrstream_labelGraphPi << "<GRAPH:";

	  inout::vectorlistprint
	    (lvectorlist_tmpgraphpi,
	     std::cout,
	     lostrstream_labelGraphPi.str().c_str(),
	     ';'
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
    if ( aiinpcga_inParamGAClustering.getWithPlotStatObjetiveFunc() ) {  
      lofh_VRC->setValue(lochrom_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

    /*3.6 Stopping criterion
      There is no stopping criterion in the relevant literature 
      which ensures the convergence of a GA to an optimal solution. 
      We have used the two most usual criteria.
      Our GA stops when:
      – After a number x of iterations, the best chromosome does 
      not change. We have fixed x = 3.
      – The maximum number of generations is reached. We chose 
      this number to be n, the number of documents.
    */

#ifdef __VERBOSE_YES

    /*ID PROC
     */
    ++geverboseui_idproc;
    
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "TERMINATION CRITERION ATTAINED?: NumMaxGenerations =  " 
		<< llfh_listFuntionHist.getDomainUpperBound()
		<< " < "
		<< aiinpcga_inParamGAClustering.getNumMaxGenerations()
		<< "  OR Not Change " << lit_iterNotChange
		<< std::endl; 
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


    if ( (llfh_listFuntionHist.getDomainUpperBound() >= aiinpcga_inParamGAClustering.getNumMaxGenerations() ) || 
	 (lit_iterNotChange == 0 ) ||
	 (runtime::elapsedTime(let_executionTime) > aiinpcga_inParamGAClustering.getMaxExecutiontime())
	 )
      break;
   
    /*3.3 Selection
      The selection operator mimics the selection concept of natural genetic systems:
      the best chromosome survives. The probability of selection of chromosome is di-
      rectly proportional to the fitness value (VRC formula for us). The chromosomes
      with the highest VRC values have more chances of reproducing and generating
      new chromosomes.
    */

    /*BEGIN GENETIC OPERATIONS*/

    { /*BEGIN SELECTION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "SELECTION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchrom_population.begin(),lvectorchrom_population.end(),
	 [](const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* liter_iChrom) -> T_REAL
	 {
	   return liter_iChrom->getFitness();
	 }
	 );

      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL--------------------------
       */

      for ( auto&& lchrom_iter: lvectorchrom_stringPool) {
	uintidx luiidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );
	*lchrom_iter = *lvectorchrom_population.at(luiidx_chrom);
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


    /*ELITIST STRATEGY IS FOUND TO WORSE THEN REPLACE---------------------------
     */

    { /*BEGIN ELITIST STRATEGY IS FOUND TO WORSE THEN REPLACE*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITIST STRATEGY IS FOUND TO WORSE THEN REPLACE";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      auto lit_chromMin = std::min_element
	(lvectorchrom_stringPool.begin(), 
	 lvectorchrom_stringPool.end(), 
	 [](const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* x, 
	    const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* y
	    ) 
	 {  return x->getFitness() < y->getFitness(); }
	 );
      
      if ( (*lit_chromMin)->getFitness() < lochrom_best.getFitness() ) {
	
        *(*lit_chromMin) = lochrom_best;  
      
      } /*IF*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END ELITIST STRATEGY IS FOUND TO WORSE THEN REPLACE*/


      /*3.4 Crossover
	Once two parents are selected, two offspring are generated. 
	These offspring will receive information from both parents. 
	The classical crossover method uses only a crossing point 
	chosen at random. This crossing point marks the position 
	where the vector will be cut in order to exchange information 
	between the parents.
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

      gaiterator::crossoverRandSelect
      (lvectorchrom_stringPool.begin(),
       lvectorchrom_stringPool.end(),
       lvectorchrom_population.begin(),
       lvectorchrom_population.end(),
       [&](gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aichrom_parent1,
	   gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aichrom_parent2,
	   gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aochrom_child1, 
	   gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aochrom_child2
	   )
       {

	 if ( uniformdis_real01(gmt19937_eng) //if  Crossover
	     < aiinpcga_inParamGAClustering.getProbCrossover() ) {
	  
	  gabinaryop::onePointCrossover
	    (*aochrom_child1,
	     *aochrom_child2,
	     *aichrom_parent1,
	     *aichrom_parent2
	     );

	    aochrom_child1->setFitness(-std::numeric_limits<T_REAL>::max());
	    aochrom_child1->setObjetiveFunc(-std::numeric_limits<T_REAL>::max());
	    aochrom_child2->setFitness(-std::numeric_limits<T_REAL>::max());
	    aochrom_child2->setObjetiveFunc(-std::numeric_limits<T_REAL>::max());
      
	} //if  Crossover
	else {
	  *aochrom_child1 = *aichrom_parent1;
	  *aochrom_child2 = *aichrom_parent2;
	}
	  
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

      /*3.5 Mutation
	Each chromosome is subjected to a low probability of change or mutation. 
	In order to guarantee that all the search space can be explored, our 
	GA uses a mutation probability of 0.008%.
      */
    { /*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for (auto&& liter_iChrom :lvectorchrom_population) {
	if ( uniformdis_real01(gmt19937_eng)
	     <  aiinpcga_inParamGAClustering.getProbMutation() ) 
	  {
	    gabinaryop::bitMutation(*liter_iChrom);
	  }
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
	
      /*END GENETIC OPERATIONS
       */
   
    --lit_iterNotChange;
    llfh_listFuntionHist.increaseDomainUpperBound();

  } /*END EVOLUTION While*/ 

 
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
    
    for (uintidx lui_i = 0; lui_i < lvectorchrom_stringPool.size(); ++lui_i) {
      delete lvectorchrom_stringPool[lui_i];
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
  aoopcga_outParamClusteringGA.setNumClusterK
    ( (T_CLUSTERIDX) lochrom_best.getNumBitOn()+1);
  aoopcga_outParamClusteringGA.setMetricFuncRun
    (lochrom_best.getObjetiveFunc());
  aoopcga_outParamClusteringGA.setFitness
    (lochrom_best.getFitness());
  aoopcga_outParamClusteringGA.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoopcga_outParamClusteringGA.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinpcga_inParamGAClustering.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinpcga_inParamGAClustering ,
       aoopcga_outParamClusteringGA
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA 
	      << ": OUT(" << geiinparam_verbose << ")\n";

    partition::PartitionDisjSets
      <T_CLUSTERIDX>
      lpartitionDisjSets_tmpClusters
      (graph::component
       (lvectorinstidx_piMST,
	lochrom_best
	)
       );
     
    uintidx lui_numClusterK
      = lochrom_best.getNumBitOn()+1;

    mat::MatrixRow<T_FEATURE> 
      lmatrixrowt_centroids
      (lui_numClusterK,
       data::Instance<T_FEATURE>::getNumDimensions() 
       );

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

    clusteringop::getCentroids
      (lmatrixrowt_centroids,
       lmatrixrowt_sumInstCluster,
       lvectort_numInstClusterK,
       lpartitionDisjSets_tmpClusters,
       aiiterator_instfirst,
       aiiterator_instlast
       );

    lochrom_best.print();
    std::cout << std::endl;
	  
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids
      << "<CENTROIDSCLUSTER: "
      << lpc_labelAlgGA
      << ": generation "
      <<  llfh_listFuntionHist.getDomainUpperBound()
      << ": lmatrixrowt_centroids["
      << &lmatrixrowt_centroids << ']';
    lmatrixrowt_centroids.print
      (std::cout,lostrstream_labelCentroids.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;

    std::ostringstream lostrstream_labelShipBest;
    lostrstream_labelShipBest
      << "<MEMBERCLUSTER: "
      << lpc_labelAlgGA
      << ": generation "
      <<  llfh_listFuntionHist.getDomainUpperBound()
      << ": lpartitionDisjSets_tmpClusters<>["
      << &lpartitionDisjSets_tmpClusters << ']';
    lpartitionDisjSets_tmpClusters.print
      (std::cout,lostrstream_labelShipBest.str().c_str(),',');
    std::cout << std::endl;

    std::vector<std::list<uintidx> >&& lvectorlist_tmpgraphpi = 
      graph::fromPiToGraph<uintidx,T_BITSIZE>
      (lvectorinstidx_piMST,lochrom_best);
    std::ostringstream lostrstream_labelGraphPi;
    lostrstream_labelGraphPi
      << "<GRAPH: "
      << lpc_labelAlgGA
      << ": generation "
      <<  llfh_listFuntionHist.getDomainUpperBound()
      << ": lvectorlist_tmpgraphpi<>["
      << &lvectorlist_tmpgraphpi << ']';
    inout::vectorlistprint
      (lvectorlist_tmpgraphpi,
       std::cout,
       lostrstream_labelGraphPi.str().c_str(),
       ';'
       );
    std::cout << std::endl;
        
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return std::make_pair(lochrom_best,lvectorinstidx_piMST);
  
} /* END gaclustering_vktreebinary */

} /*END eac */

#endif /*__GACLUSTERING_VKTREEBINARY_HPP__*/
