/*! \file gas_fkcentroid.hpp
 *
 * \brief GAS \cite Maulik:Bandyopadhyay:GAclustering:GAS:2000
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GAS algorithm based on the paper:\n
 * U. Maulik and S. Bandyopadhyay. Genetic algorithm-based clustering\n
 * technique. Pattern Recognition, 33(9):1455–1465, 2000.\n
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
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
 

#ifndef __GAS_FKCENTROID_HPP__
#define __GAS_FKCENTROID_HPP__

#include <vector>
#include <algorithm>

#include <leac.hpp>
#include "inparam_probcprobm_fixedk.hpp"
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
  
/*! \fn gaencode::ChromFixedLength<T_FEATURE,T_REAL> gas_fkcentroid (inout::OutParamEAClustering<T_REAL,T_CLUSTERIDX> &aoop_outParamEAC, inout::InParamPcPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinpcgaprobfixedk_inParamKGA, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief GAS \cite Maulik:Bandyopadhyay:GAclustering:GAS:2000
  \details Implementation of GAS algorithm based on \cite Maulik:Bandyopadhyay:GAclustering:GAS:2000. 
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid. Base to following equation:
  \f[
  x_i \in C_j \leftrightarrow  \| x_i - \mu_j \|  \begin{array}{c}min\\  k \end{array}
  \| x_i - \mu_k \|,\; j=1,2,..k,
  \f]
  where \f$m_j\f$, represents the medoid of cluster \f$C_j\f$
  \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
  \param aiinpcgaprobfixedk_inParamKGA a inout::InParamPcPmFk parameters required by the algorithm
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
gas_fkcentroid
(inout::OutParamEAClustering
 <T_REAL,
 T_CLUSTERIDX>                      &aoop_outParamEAC,
 inout::InParamPcPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>              &aiinpcgaprobfixedk_inParamKGA,
 const INPUT_ITERATOR                aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist
 )
{  
  const uintidx lconstui_numClusterFk =
    (uintidx) aiinpcgaprobfixedk_inParamKGA.getNumClusterK();

  /*ASSIGN SIZE FOR ALL CHROMOSOMES
   */
  gaencode::ChromFixedLength<T_FEATURE,T_REAL>::setStringSize
    ( lconstui_numClusterFk * data::Instance<T_FEATURE>::getNumDimensions() );

  gaencode::ChromFixedLength<T_FEATURE,T_REAL> lochromfixleng_best;
  
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   
   *POPULATION CREATE------------------------------------------------------------
   */
  std::vector<gaencode::ChromFixedLength<T_FEATURE,T_REAL> >
    lvectorchromfixleng_population
    (aiinpcgaprobfixedk_inParamKGA.getSizePopulation());

  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
   */
  std::vector<gaencode::ChromFixedLength<T_FEATURE,T_REAL> >
    lvectorchromfixleng_matingPool
    (aiinpcgaprobfixedk_inParamKGA.getSizePopulation());

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  

#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gas_fkcentroids"; 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelAlgGA
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t(output Chromosome: lochromfixleng_best[" 
      << &lochromfixleng_best << "]\n"
      << "\t output inout::OutParamEAClustering&: "
      << "aoop_outParamEAC[" 
      << &aoop_outParamEAC << "]\n"
      << "\t input  InParamPcPmFk&: "
      <<"aiinpcgaprobfixedk_inParamKGA[" 
      << &aiinpcgaprobfixedk_inParamKGA << "]\n"
      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      << &aifunc2p_dist << ']'
      << "\n\t\tPopulation size = " 
      << aiinpcgaprobfixedk_inParamKGA.getSizePopulation()
      << "\n\t\tProbCrossover = " 
      << aiinpcgaprobfixedk_inParamKGA.getProbCrossover() 
      << "\n\t\tProbMutation  = " 
      << aiinpcgaprobfixedk_inParamKGA.getProbMutation() 
      << "\n\t)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 

  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinpcgaprobfixedk_inParamKGA.getNumMaxGenerations(),
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

  if ( aiinpcgaprobfixedk_inParamKGA.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      ( aiinpcgaprobfixedk_inParamKGA.getSizePopulation());
    //DEFINE FUNCTION
    lofh_SSE  = new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinpcgaprobfixedk_inParamKGA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinpcgaprobfixedk_inParamKGA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamEAC.setFileNameOutPlotStatObjetiveFunc
      (aiinpcgaprobfixedk_inParamKGA.getFileNamePlotStatObjetiveFunc(),
       aiinpcgaprobfixedk_inParamKGA.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open
      (aoop_outParamEAC.getFileNameOutPlotStatObjetiveFunc().c_str(),
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
  aoop_outParamEAC.setTotalInvalidOffspring(0);

  /*OUT: GENETIC ALGORITHM CHARACTERIZATION*/
  runtime::ExecutionTime let_executionTime = runtime::start();
    

  /*INITIALIZE POPULATION--------------------------------------------------------
    3.2.2 POPULATION INITIALIZATION    
    The K cluster centres encoded in each chromosome 
    are initialized to K randomly chosen points from the
    data set. This process es repeated for each of the P 
    chromosomes in the population, where P is the size of the
    population. \cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
  */
  {/*BEGIN INITIALIZE POPULATION P(t)*/     
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "(0) POPULATION INITIAL";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
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

      lchromfixleng_iter.setFitness(-std::numeric_limits<T_REAL>::max());
      lchromfixleng_iter.setObjetiveFunc(std::numeric_limits<T_REAL>::max());
	 
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

    /*BEGIN ITERATION*/
    llfh_listFuntionHist.increaseDomainUpperBound();

    /*CLUSTERING----------------------------------------------------------
      The fitness computation process consists of two
      phases. In the "rst phase, the clusters are formed accord-
      ing to the centres encoded in the chromosome under
      consideration. This is done by assigning each point
      x_i, i =1,2,...,n, to one of the clusters C_i with centre
      z_i
      \cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
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
   

    /*FITNESS FUNCTION----------------------------------------------------------
      Subsequently, the clustering metric M is computed as
      follows:
      M = \sum_{i=1}^{k} M_i
      M_i = \sum_{x_j\inC_i} ||x_j-zi||
      The fitness fuction is defined as f= 1/M, so that maxi-
      mization of the "tness function leads to minimization
      of M.
      \cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
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

      } 
        
      aoop_outParamEAC.sumTotalInvalidOffspring
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

    /*ELITISM-------------------------------------------------------------------
      The best string seen
      upto the last generation provides the solution to the
      lustering problem. We have implemented elitism at each
      generation by preserving the best string seen upto that
      generation in a location outside the population. Thus on
      ermination, this location contains the centres of the fnal
      lusters.
      \cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
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
       
	aoop_outParamEAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamEAC.setRunTimeGetsBest
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
    if ( aiinpcgaprobfixedk_inParamKGA.getWithPlotStatObjetiveFunc() ) {  

      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());

      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

    /*TERMINATION CRITERION-----------------------------------------------------
      3.2.7 TERMINATION CRITERION 
      \cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
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
	   < aiinpcgaprobfixedk_inParamKGA.getNumMaxGenerations() ) 
	 )
      break;

    /*SELECTION------------------------------------------------------------------
      The selection process selects chromosomes from the
      mating pool directed by the survival of the fittest concept
      of natural genetic systems. In the proportional selection
      strategy adopted in this article, a chromosome is assigned
      a number of copies, which is proportional to its fitness in
      the population, that go into the mating pool for further
      genetic operations. Roulette wheel selection is one com-
      mon technique that implements the proportional selec-
      tion strategy.
      \cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
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
      
      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL--------------------------
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


    /*CROSSOVER------------------------------------------------------------------
      Crossover. Crossover is a probabilistic process that 
      exchanges information between two parent chromosomes 
      for generating two offspring. Here, single-point 
      crossover with a fixed crossover probability of $\mu_c$ 
      is used. For chromosomes of length $l$ $(l= NK)$, a random 
      integer, called the crossover point, is generated in the 
      range $[1,l-1]$. The portions of the chromosomes lying to
      the right of the crossover point are exchanged to produce 
      two offspring.
      \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
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
		aiinpcgaprobfixedk_inParamKGA.getProbCrossover()  ) {
	      
	     gagenericop::onePointCrossover
	       (aochrom_child1,
		aochrom_child2,
		aichrom_parent1,
		aichrom_parent2
		);

	   } //if  Crossover
	   else {
	     aochrom_child1 =  aichrom_parent1;
	     aochrom_child2 =  aichrom_parent2;
	   }
	 }
	 );
	 
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

    
      /*MUTATION----------------------------------------------------------------
       */
    {/*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ": IN(" << geiinparam_verbose << ')'
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto& lchromfixleng_iter: lvectorchromfixleng_population ) {
	if ( uniformdis_real01(gmt19937_eng) 
	     < aiinpcgaprobfixedk_inParamKGA.getProbMutation() ) 
	  { //IF BEGIN  MUTATION	
	    garealop::randomMutation(lchromfixleng_iter);
	    //THE FITNESS VALUE AND OBJECTIVE IS RESET
	    lchromfixleng_iter.setFitness
	      (-std::numeric_limits<T_REAL>::max());  
	    lchromfixleng_iter.setObjetiveFunc
	      (std::numeric_limits<T_REAL>::max());  
	  } //END BEGIN  MUTATION
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


    } /*END MUTATION */
      
  } /*END EVOLUTION While*/ 
    
  runtime::stop(let_executionTime);
  aoop_outParamEAC.setNumClusterK
    (aiinpcgaprobfixedk_inParamKGA.getNumClusterK());
  aoop_outParamEAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinpcgaprobfixedk_inParamKGA.getWithPlotStatObjetiveFunc() ) {  
    runtime::plot_funtionHist
      (llfh_listFuntionHist,
       aiinpcgaprobfixedk_inParamKGA,
       aoop_outParamEAC
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
 
} /* END gas_fkcentroid */

} /*END eac */
  
#endif /*__GAS_FKCENTROID_HPP__*/
