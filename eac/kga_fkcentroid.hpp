/*! \file kga_fkcentroid.hpp
 *
 * \brief KGA \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the KGA algorithm based on the paper:\n
 * S. Bandyopadhyay and U. Maulik. An evolutionary technique based\n
 * on k-means algorithm for optimal clustering in rn. Inf. Sci. Appl.,\n 
 * 146(1-4):221–237, 2002.\n
 * URL: http://www.sciencedirect.com/science/article/pii/S0020025502002086,\n 
 * <a href="http://dx.doi.org/10.1016/S0020-0255(02)00208-6">doi:http://dx.doi.org/10.1016/S0020-0255(02)00208-6</a>\n.
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

#ifndef __KGA_FKCENTROID_HPP__
#define __KGA_FKCENTROID_HPP__

#include <vector>
#include <algorithm>

#include <leac.hpp>
#include "inparam_pcpmfk.hpp"
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
  
/*! \fn gaencode::ChromFixedLength<T_FEATURE,T_REAL> kga_fkcentroid (inout::OutParamGAC<T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamPcPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamPcPmFk, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief  KGA \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002
  \details Implementation of the KGA algorithm based on \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002. 
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid. Base to following equation:
  \f[
  x_i \in C_j \leftrightarrow  \| x_i - \mu_j \|  \begin{array}{c}min\\  k \end{array}
  \| x_i - \mu_k \|,\; j=1,2,..k,
  \f]
  where \f$mu_j\f$, represents the centroid of cluster \f$C_j\f$
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinp_inParamPcPmFk a inout::InParamPcPmFk parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_REAL,  
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromFixedLength<T_FEATURE,T_REAL> 
kga_fkcentroid
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                       &aoop_outParamGAC,
 inout::InParamPcPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>              &aiinp_inParamPcPmFk,
 const INPUT_ITERATOR                aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist
 )
{  
   uintidx lconstui_numClusterFk =
    (uintidx) aiinp_inParamPcPmFk.getNumClusterK();

  /*ASSIGN SIZE FOR ALL CHROMOSOMES
   */
  gaencode::ChromFixedLength<T_FEATURE,T_REAL>::setStringSize
    ( lconstui_numClusterFk * data::Instance<T_FEATURE>::getNumDimensions() );

  gaencode::ChromFixedLength<T_FEATURE,T_REAL> lochromfixleng_best;
  
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   */
   
  /*POPULATION CREATE
   */
  std::vector<gaencode::ChromFixedLength<T_FEATURE,T_REAL> >
    lvectorchromfixleng_population
    (aiinp_inParamPcPmFk.getSizePopulation());

  /*CREATE SPACE FOR STORE MATINGPOOL
   */
  std::vector<gaencode::ChromFixedLength<T_FEATURE,T_REAL> >
    lvectorchromfixleng_matingPool
    (aiinp_inParamPcPmFk.getSizePopulation());

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  

#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "kga_fkcentroids"; 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t(output Chromosome: lochromfixleng_best[" 
      << &lochromfixleng_best << "]\n"
      << "\t output inout::OutParamGAC&: "
      << "aoop_outParamGAC[" 
      << &aoop_outParamGAC << "]\n"
      << "\t input  InParamPcPmFk&: "
      <<"aiinp_inParamPcPmFk[" 
      << &aiinp_inParamPcPmFk << "]\n"
      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      << &aifunc2p_dist << ']'
      << "\n\t\tPopulation size = " 
      << aiinp_inParamPcPmFk.getSizePopulation()
      << "\n\t\tProbCrossover = " 
      << aiinp_inParamPcPmFk.getProbCrossover() 
      << "\n\t\tProbMutation  = " 
      << aiinp_inParamPcPmFk.getProbMutation() 
      << "\n\t)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 

  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamPcPmFk.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream               lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>
    *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>         lvectorT_statfuncObjetiveFunc;

  if ( aiinp_inParamPcPmFk.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      ( aiinp_inParamPcPmFk.getSizePopulation());
    //DEFINE FUNCTION
    lofh_SSE  = new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinp_inParamPcPmFk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamPcPmFk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamPcPmFk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamPcPmFk.getTimesRunAlgorithm()
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
  aoop_outParamGAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/
  runtime::ExecutionTime let_executionTime = runtime::start();
    
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
  
  /*3.1.1 Population initialization 
    Chosen distict points from the data set are used to initialize 
    the K cluster centers encoded in each choromosome.
    This is similar to the initialization od the centers 
    in K-Means algorithm. This process es repeat for each chromosome 
    in the population 
    \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
  */
  {/*BEGIN INITIALIZE POPULATION P(t)*/
      
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "(0) POPULATION INITIAL";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< geverbosepc_labelstep  
	<< ": IN(" << geiinparam_verbose << ')'
	<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    for ( auto& lchromfixleng_iter: lvectorchromfixleng_population ) {
      
      /*DECODE CHROMOSOME
       */
      mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroidsChrom
	(lconstui_numClusterFk,
	 data::Instance<T_FEATURE>::getNumDimensions(),
	 lchromfixleng_iter.getString()
	 );
 	 
      clusteringop::randomInitialize
	(lmatrixrowt_centroidsChrom,
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );

      lchromfixleng_iter.setFitness
	(-std::numeric_limits<T_REAL>::max());
      lchromfixleng_iter.setObjetiveFunc
	(std::numeric_limits<T_REAL>::max());
	 
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

  } /*END INITIALIZE POPULATION P(t)*/

  while ( 1 ) {
    
    /*BEGIN ITERATION
     */
    llfh_listFuntionHist.increaseDomainUpperBound();

    /*CLUSTERING
      In this step, the cluster are formed according to the center 
      encoded in the chromosome. This is done by assigning each 
      point x_i,i=1,2,...,n to one of the clusters Cj with center z_i* 
      such that ||x_i-z_j|| <= ||x_i-z_p||, p=1,2,...,K, and p != j.
      All ties are resolved arbitrarily. As like the K-Means algorithm, 
      for each cluster
      C_i , its new center z^* is computed as

      z_i^* = 1/n_i \sum_{x_j\inC_j} x_j, i=1,2,...,K,
      where ni is the number of points in cluster C_i . These z^* 
      now replace the preious zi's in the chromosome.
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
    */
    {/*BEGIN CLUSTERING*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "A. THE CLUSTERS ARE FORMED";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto& liter_iChrom: lvectorchromfixleng_population ) {
	
	/*DECODE CHROMOSOME*/
	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroidsChrom
	  (lconstui_numClusterFk,
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   liter_iChrom.getString()
	   );

	mat::MatrixRow<T_FEATURE_SUM>       
	  llmatrixrowt_sumInstancesCluster
	  (lconstui_numClusterFk,
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   T_FEATURE_SUM(0)
	   );
	    
	std::vector<T_INSTANCES_CLUSTER_K> 
	  lvectort_numInstancesInClusterK
	  (lconstui_numClusterFk,
	   T_INSTANCES_CLUSTER_K(0)
	   );
	    
	T_CLUSTERIDX lmcidx_numClusterNull;
	    
	clusteringop::updateCentroids
	  (lmcidx_numClusterNull,
	   lmatrixrowt_centroidsChrom,
	   llmatrixrowt_sumInstancesCluster,
	   lvectort_numInstancesInClusterK,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
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

    } /*END CLUSTERING*/

    
    /*FITNESS FUNCTION
      3.1.3 FITNESS COMPUTATION 
      For each chromosome, the clusters formed in the previous 
      step are utilized computing the clustering metric, M, as 
      follows: M =  \sum_{i=1}{x}\sum_{x_j\inC_j} ||x_j-z_i}||
      For finding the appropriate clusters M has to be minimized. 
      The fitness function of a chromosome is defined as 1/M. 
      Therefore, maximization of the fitness function will lead 
      to minimization of the clustering metric M.
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002} 
    */
    { /*BEGIN COMPUTED METRIC M AND FITNESS*/
     
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "B. COMPUTED METRIC M AND FITNESS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ": IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      long ll_invalidOffspring = 0;

      for ( auto& lchromfixleng_iter: lvectorchromfixleng_population ) {
	 
	/*DECODE CHROMOSOME*/
	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroidsChrom
	  (lconstui_numClusterFk,
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   lchromfixleng_iter.getString()
	   );

	std::pair<T_REAL,bool> lpair_SSE =
	  um::SSE 
	  (lmatrixrowt_centroidsChrom,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
	   
	lchromfixleng_iter.setObjetiveFunc(lpair_SSE.first);
	lchromfixleng_iter.setFitness(1.0 / lpair_SSE.first);
	lchromfixleng_iter.setValidString(lpair_SSE.second);
	    
	if ( lchromfixleng_iter.getValidString() == false )
	  ++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back
	  (lchromfixleng_iter.getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

      } //End for
        
      aoop_outParamGAC.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    } /*END COMPUTED METRIC M AND FITNESS*/


    /*ELITISM
      Elitism has been implemented in each generation by 
      replacing the worst chromosome of the population with 
      the best one seen up to the previous generation.
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
    */
    { /*BEGIN ELITISM REPLACING THE WORST CHROMOSOME*/ 

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITISM REPLACING THE WORST CHROMOSOME";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      auto lit_chromMin =
	std::min_element   
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& x,
	    const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& y
	    )
	 { return x.getFitness() < y.getFitness(); }
	 );
     
      if ( lit_chromMin->getFitness() < lochromfixleng_best.getFitness() ) {

	*lit_chromMin = lochromfixleng_best;

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

    } /*END ELITISM REPLACING THE WORST CHROMOSOME*/

    
    /*The best string or chromosome seen up to the last generation 
      provides the solution to the clustering problem. 
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
    */
    { /*BEGIN PRESERVING THE BEST STRING*/

      auto lchromfixleng_iterMax  =
	std::max_element
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& x, 
	    const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& y
	    ) 
	 {  return x.getFitness() < y.getFitness(); }
	 );

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ": IN(" << geiinparam_verbose << ")\tmax fitness = "
	  << lchromfixleng_iterMax->getFitness()  
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      if ( lochromfixleng_best.getFitness() <
	   lchromfixleng_iterMax->getFitness() ) {

	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	lochromfixleng_best  = *lchromfixleng_iterMax;
       
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

    } /*END PRESERVING THE BEST STRING*/

  
    /*MEASUREMENT BEST: COMPUTING STATISTICAL AND METRIC OF THE 
      ALGORITHM
    */
#ifndef __WITHOUT_PLOT_STAT
    if ( aiinp_inParamPcPmFk.getWithPlotStatObjetiveFunc() ) {  

      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());

      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

    /*TERMINATION CRITERION
      3.1.5 TERMINATION CRITERION 
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
    */
#ifdef __VERBOSE_YES
    /*ID PROC
     */
    ++geverboseui_idproc;
      
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< "TERMINATION CRITERION ATTAINED?: " 
	<< llfh_listFuntionHist.getDomainUpperBound()
	<< std::endl; 
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
   
    if ( !(llfh_listFuntionHist.getDomainUpperBound() 
	   < aiinp_inParamPcPmFk.getNumMaxGenerations() ) 
	 )
      break;

    /*3.1.4 GENETIC OPERATIONS
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
    */

    /*SELECTION
      Selection. The selection process selects chromosomes from the 
      mating pool directed by the survival of the fittest concept 
      of natural genetic systems. In the proportional selection 
      strategy adopted in this paper, a chromosome is assigned a 
      number of copies, which is proportional to its fitness 
      in the population.
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
    */
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
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>&
	    lchromfixleng_iter) -> T_REAL
	 {
	   return lchromfixleng_iter.getFitness();
	 }
	 );
      
      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL
       */ 
      for ( auto& lchromfixleng_iter: lvectorchromfixleng_matingPool ) {
	
	uintidx lstidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );

	lchromfixleng_iter = lvectorchromfixleng_population.at(lstidx_chrom);
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

    } /*END SELECTION*/

    /*CROSSOVER
      Crossover. Crossover is a probabilistic process that 
      exchanges information between two parent chromosomes 
      for generating two offspring. Here, single-point 
      crossover with a fixed crossover probability of $\mu_c$ 
      is used. For chromosomes of length $l$ $(l= NK)$, a random 
      integer, called the crossover point, is generated in the 
      range $[1,l-1]$. The portions of the chromosomes lying to
      the right of the crossover point are exchanged to produce 
      two offspring. \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
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

      long ll_invalidOffspring = 0;
 
      gaiterator::crossover
	(lvectorchromfixleng_matingPool.begin(),
	 lvectorchromfixleng_matingPool.end(),
	 lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>&
	     aichrom_parent1,
	     const gaencode::ChromFixedLength<T_FEATURE,T_REAL>&
	     aichrom_parent2,
	     gaencode::ChromFixedLength<T_FEATURE,T_REAL>&
	     aochrom_child1, 
	     gaencode::ChromFixedLength<T_FEATURE,T_REAL>&
	     aochrom_child2
	     )
	 {
	 
	   if ( uniformdis_real01(gmt19937_eng) <
		aiinp_inParamPcPmFk.getProbCrossover()  ) {
	      
	     gagenericop::onePointCrossover
	       (aochrom_child1,
		aochrom_child2,
		aichrom_parent1,
		aichrom_parent2
		);

	     /*DECODE CHROMOSOME CHILD1*/
	     mat::MatrixRow<T_FEATURE> 
	       lmatrixrowt_centroidsChromChild1
	       (lconstui_numClusterFk,
		data::Instance<T_FEATURE>::getNumDimensions(),
		aochrom_child1.getString()
		);

	     std::pair<T_REAL,bool> 
	       lpair_SSE1 =
	       um::SSE
	       (lmatrixrowt_centroidsChromChild1,
		aiiterator_instfirst,
		aiiterator_instlast, 
		aifunc2p_dist
		);
	     aochrom_child1.setObjetiveFunc(lpair_SSE1.first);
	     aochrom_child1.setFitness(1.0 / lpair_SSE1.first);
	     aochrom_child1.setValidString(lpair_SSE1.second);
	    
	     if ( aochrom_child1.getValidString() == false )
	       ++ll_invalidOffspring;

	     /*DECODE CHROMOSOME CHILD1*/
	     mat::MatrixRow<T_FEATURE> 
	       lmatrixrowt_centroidsChromChild2
	       (lconstui_numClusterFk,
		data::Instance<T_FEATURE>::getNumDimensions(),
		aochrom_child2.getString()
		);

	     std::pair<T_REAL,bool> 
	       lpair_SSE2 =
	       um::SSE 
	       (lmatrixrowt_centroidsChromChild2,
		aiiterator_instfirst,
		aiiterator_instlast,
		aifunc2p_dist
		);
	    
	     aochrom_child2.setObjetiveFunc(lpair_SSE2.first);
	     aochrom_child2.setFitness(1.0 / lpair_SSE2.first);
	     aochrom_child2.setValidString(lpair_SSE2.second);
	    
	     if ( aochrom_child2.getValidString() == false )
	       ++ll_invalidOffspring;

	   } //if  Crossover
	   else {

	     aochrom_child1 =  aichrom_parent1;
	     aochrom_child2 =  aichrom_parent2;
	   
	   }
	 }
	 );
	 
      aoop_outParamGAC.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
    } /*END CROSSOVER*/
    
      /*MUTATION
	Mutation. Each chromosome undergoes mutation with 
        a fixed probability $\mu_c$. Let $M_min$ and $M_max$ 
        be the minimum and maximum values of the clustering 
        metric, respectively, in the current population. For 
        mutating a chromosome whose clustering metric is $M$, 
        a number... \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
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

      auto lchrom_minObjFunc  =
	std::min_element
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& x, 
	    const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& y
	    )
	 {  return x.getObjetiveFunc() < y.getObjetiveFunc(); }
	 );
      T_REAL  lrt_minClusteringMetric =
	lchrom_minObjFunc->getObjetiveFunc();
      
      auto lchrom_maxObjFunc  =
	std::max_element
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& x,
	    const gaencode::ChromFixedLength<T_FEATURE,T_REAL>& y
	    )
	 {  return x.getObjetiveFunc() < y.getObjetiveFunc(); }
	    
	 );
      
      T_REAL lrt_maxClusteringMetric =
	lchrom_maxObjFunc->getObjetiveFunc();
      
      for ( auto& lchromfixleng_iter: lvectorchromfixleng_population ) {
	if ( uniformdis_real01(gmt19937_eng)
	     < aiinp_inParamPcPmFk.getProbMutation() ) 
	  { //IF MUTATION	
	    gaclusteringop::biDirectionHMutation
	      (lchromfixleng_iter,
	       lrt_minClusteringMetric,
	       lrt_maxClusteringMetric,
	       larray_minFeactures,
	       larray_maxFeactures
	       );
	    lchromfixleng_iter.setFitness
	      (-std::numeric_limits<T_REAL>::max());  
	    lchromfixleng_iter.setObjetiveFunc
	      (std::numeric_limits<T_REAL>::max());
	  } //END IF MUTATION
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
      
    } /*END MUTATION*/

      
  } /*END EVOLUTION While*/ 

  /*FREE MEMORY
   */
  delete [] larray_maxFeactures;
  delete [] larray_minFeactures;
    
  runtime::stop(let_executionTime);
  aoop_outParamGAC.setNumClusterK
    (aiinp_inParamPcPmFk.getNumClusterK());
  aoop_outParamGAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamGAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamPcPmFk.getWithPlotStatObjetiveFunc() ) {  
    runtime::plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamPcPmFk,
       aoop_outParamGAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout
      << lpc_labelAlgGA 
      << ": OUT(" << geiinparam_verbose << ")\n";
    lochromfixleng_best.print();
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  return lochromfixleng_best;
 
} /* END kga_fkcentroid */

} /*END eac */
  
#endif /*__KGA_FKCENTROID_HPP__*/
