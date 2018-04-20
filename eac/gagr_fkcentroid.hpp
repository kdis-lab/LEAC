/*! \file gagr_fkcentroid.hpp
 *
 * \brief GAGR \cite Chang:etal:GAclustering:GAGR:2009

 * \details This file is part of the LEAC.\n\n
 * Implementation of the GAGR algorithm based on the paper:\n
 * Dong-Xia Chang, Xian-Da Zhang, and Chang-Wen Zheng. A genetic\n
 * algorithm with gene rearrangement for k-means clustering.\n
 * Pattern Recogn., 42(7):1210–1222, 2009.\n
 * <a href="http://dx.doi.org/10.1016/j.patcog.2008.11.006">doi:http://dx.doi.org/10.1016/j.patcog.2008.11.006</a>\n.
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

#ifndef __GAGR_FKCENTROID_HPP__
#define __GAGR_FKCENTROID_HPP__

#include <vector>

#include <leac.hpp>
#include "inparam_adaptivepcpmfk.hpp"
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
  
/*! \fn gaencode::ChromFixedLength<T_FEATURE,T_REAL> gagr_fkcentroid(inout::OutParamGAC <T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, const inout::InParamAdaptivePcPm<T_CLUSTERIDX,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamProbAdaptive, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief GAGR \cite Chang:etal:GAclustering:GAGR:2009  
  \details Implementation of GAGR algorithm based on \cite Chang:etal:GAclustering:GAGR:2009. 
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid. Base to following equation:
  \f[
  x_i \in C_j \leftrightarrow  \| x_i - \mu_j \|  \begin{array}{c}min\\  k \end{array}
  \| x_i - \mu_k \|,\; j=1,2,..k,
  \f]
  where \f$m_j\f$, represents the medoid of cluster \f$C_j\f$
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinp_inParamProbAdaptive a inout::InParamAdaptivePcPm parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE, //DATA TYPE OF A GENE,
	   typename T_REAL,  //OBJETIVE FUNCTION, FITNESS AND CLUSTERING METRIC
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromFixedLength<T_FEATURE,T_REAL>  
gagr_fkcentroid
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                     &aoop_outParamGAC,
 const inout::InParamAdaptivePcPmFk
 <T_CLUSTERIDX,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>             &aiinp_inParamProbAdaptive,
 const INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  
  /*ID PROC
   */
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  geverboseui_idproc = 1;
  T_CLUSTERIDX *larraymcidx_memberShipTmp =
    new T_CLUSTERIDX[lui_numInstances]; 
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gagr_fkcentroid";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output inout::OutParamGAC&: aoop_outParamGAC[" 
	      << &aoop_outParamGAC << "]\n"
	      << "\t input  InParamClusteringGAProbAdap&: aiinp_inParamProbAdaptive[" 
	      << &aiinp_inParamProbAdaptive << "]\n"
	      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t\tPopulation size = " 
	      << aiinp_inParamProbAdaptive.getSizePopulation()
	      << "\n\t\trandom-seed = "
	      << aiinp_inParamProbAdaptive.getRandomSeed()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 
  
  runtime::ListRuntimeFunction<COMMON_IDOMAIN>  
    llfh_listFuntionHist
    (aiinp_inParamProbAdaptive.getNumMaxGenerations(), "Iterations", "Clustering metrics");

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT 
  std::ofstream lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL>  *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>   *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>        lvectorT_statfuncObjetiveFunc;

  if ( aiinp_inParamProbAdaptive.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      (aiinp_inParamProbAdaptive.getSizePopulation());
    //DEFINE FUNCTION
    lofh_SSE  = 
      new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinp_inParamProbAdaptive.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);
    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamProbAdaptive.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamProbAdaptive.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamProbAdaptive.getTimesRunAlgorithm()
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

  auto lfuncobjgagr_functionObjetive =
    gafuncobj::makeGAFuncObjSSE
    (aiinp_inParamProbAdaptive.getNumClusterK(),
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   */
  const uintidx lconstui_numClusterFk =
    (uintidx) aiinp_inParamProbAdaptive.getNumClusterK();
  
  gaencode::ChromFixedLength<T_FEATURE,T_REAL>::setStringSize
    ( lconstui_numClusterFk * data::Instance<T_FEATURE>::getNumDimensions() );

  gaencode::ChromFixedLength<T_FEATURE,T_REAL> lochromfixleng_best;
  lochromfixleng_best.setObjetiveFunc(std::numeric_limits<T_REAL>::max());
  lochromfixleng_best.setFitness(-std::numeric_limits<T_REAL>::max());

  std::vector<gaencode::ChromFixedLength<T_FEATURE,T_REAL>* >
    lvectorchromfixleng_population;

  std::vector<gaencode::ChromFixedLength<T_FEATURE,T_REAL>* >
    lvectorchromfixleng_matingPool;

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING*/
  aoop_outParamGAC.setTotalInvalidOffspring(0);

  runtime::ExecutionTime let_executionTime = runtime::start();

  llfh_listFuntionHist.increaseDomainUpperBound();
         
  /*VARIBLES EXTRAS OF THE ALGORITHM*/
  
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


  /*POPULATION CREATE------------------------------------------------------------
   */
  lvectorchromfixleng_population.reserve
    ( aiinp_inParamProbAdaptive.getSizePopulation() );

  for (uintidx lui_i = 0; 
       lui_i < aiinp_inParamProbAdaptive.getSizePopulation(); 
       lui_i++) 
    {
      lvectorchromfixleng_population.push_back
	(new gaencode::ChromFixedLength<T_FEATURE,T_REAL>());
    }

  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
   */
  lvectorchromfixleng_matingPool.reserve
    (aiinp_inParamProbAdaptive.getSizePopulation());
  
  for (uintidx lui_i = 0; 
       lui_i < aiinp_inParamProbAdaptive.getSizePopulation(); 
       lui_i++) 
    {
      lvectorchromfixleng_matingPool.push_back
	(new gaencode::ChromFixedLength<T_FEATURE,T_REAL>());
    }
  
  /*1. Initialize a group of cluster centers with size of P, only valid
    chromosomes (that have at least one data point in each cluster)
    are taken into consideration. Each data point of the set is assigned
    to the cluster with closest cluster center using the Euclidean
    distance. \cite{Chang:etal:GAclustering:GAGR:2009}
  */
  { /*BEGIN 1. INITIALIZE A GROUP OF CLUSTER CENTERS
      INITIALIZE POPULATION
    */

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "1. INITIALIZE A GROUP OF CLUSTER CENTERS";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    for ( auto lchromfixleng_iter: lvectorchromfixleng_population ) {
      /*DECODE CHROMOSOME*/
      mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroidsChrom
	(lconstui_numClusterFk,
	 data::Instance<T_FEATURE>::getNumDimensions(),
	 lchromfixleng_iter->getString()
	 );
	
      clusteringop::randomInitialize
	(lmatrixrowt_centroidsChrom,
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );

      lchromfixleng_iter->setObjetiveFunc(std::numeric_limits<T_REAL>::max());
      lchromfixleng_iter->setFitness(-std::numeric_limits<T_REAL>::min());
    }
  
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END 1. INITIALIZE A GROUP OF CLUSTER CENTERS
      INITIALIZE POPULATION
    */

     
  {/*BEGIN EACH DATA POINT OF THE SET IS ASSIGNED TO CLUSTER
     USING EUCLIDEAN DISTANCE
   */
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "EACH DATA POINT OF THE SET IS ASSIGNED TO CLUSTER";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep 
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    for ( auto lchromfixleng_iter: lvectorchromfixleng_population ) {
      
      //DECODE CHROMOSOME
      mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroidsChrom
	(lconstui_numClusterFk,
	 data::Instance<T_FEATURE>::getNumDimensions(),
	 lchromfixleng_iter->getString()
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
	    
	 
#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	   
	std::ostringstream lostrstream_labelCentroids;
	lostrstream_labelCentroids
	  << "<CENTROIDSCLUSTER:"
	  << geverbosepc_labelstep  
	  << ":gaencode::ChromFixedLength:lchromfixleng_iter["
	  << geverboseui_idproc  << ':'  << lchromfixleng_iter
	  << ']';			      
	lmatrixrowt_centroidsChrom.print
	  (std::cout,
	   lostrstream_labelCentroids.str().c_str(),
	   ',',
	   ';'
	   );
	std::cout << '\n';
	
	clusteringop::reassignCluster
	  (larraymcidx_memberShipTmp,
	   lmatrixrowt_centroidsChrom,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
 		  
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    }

#ifdef __VERBOSE_YES 
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END EACH DATA POINT OF THE SET IS ASSIGNED TO CLUSTER
     */
  
    /*2. Evaluate each chromosome and copy the best chromosome say
      pbest of the initial population in a separate location.
      \cite{Chang:etal:GAclustering:GAGR:2009}
    */
  
  { /*BEGIN 2. EVALUATE EACH CROMOSOME
     */
#ifdef __VERBOSE_YES
    const char* lpc_labelStep = "2. EVALUATE EACH CROMOSOME";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelStep  
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    long ll_invalidOffspring = 0;
    
    for ( auto lchromfixleng_iter: lvectorchromfixleng_population ) {
	 
      std::pair<T_REAL,bool> lpair_SSE = 
	lfuncobjgagr_functionObjetive.getObjetiveFunc
	(lchromfixleng_iter->getString());
      lchromfixleng_iter->setObjetiveFunc(lpair_SSE.first);
      lchromfixleng_iter->setFitness
	(lfuncobjgagr_functionObjetive.getFitness(lpair_SSE.first));
      lchromfixleng_iter->setValidString(lpair_SSE.second);
	 
      if ( lchromfixleng_iter->getValidString() == false )
	++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
      lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

    } //END std::for_each ObjetiveFunc
       
    aoop_outParamGAC.sumTotalInvalidOffspring
      (ll_invalidOffspring);

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelStep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END 2. EVALUATE EACH CROMOSOME
     */
    
  { /* BEGIN COPY THE BEST CHROMOSOME SAY Pbes
     */
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "COPY THE BEST CHROMOSOME SAY Pbes";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep 
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchromfixleng_iterMax  = 
      *(std::max_element
	(lvectorchromfixleng_population.begin(), 
	 lvectorchromfixleng_population.end(), 
	 [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	     const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	     ) 
      {  return x->getFitness() < y->getFitness(); } 
	 )
	);

    if ( lochromfixleng_best.getFitness() <  lchromfixleng_iterMax->getFitness() ) {
    
      /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
      lochromfixleng_best  = *lchromfixleng_iterMax;
      
      aoop_outParamGAC.setIterationGetsBest
	(llfh_listFuntionHist.getDomainUpperBound());
      aoop_outParamGAC.setRunTimeGetsBest
	(runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	
	/*DECODE CHROMOSOME*/
	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroidsChromBest
	  (lconstui_numClusterFk,
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   lochromfixleng_best.getString()
	   );
	 
	std::ostringstream lostrstream_labelCentroids;
	lostrstream_labelCentroids
	  << "<CENTROIDSCLUSTER:"
	  <<  geverbosepc_labelstep //lpc_labelStep 
	  << ":Chromosome:lochromfixleng_best["
	  << geverboseui_idproc << ':' << &lochromfixleng_best
	  << ']';
	lmatrixrowt_centroidsChromBest.print
	  (std::cout,
	   lostrstream_labelCentroids.str().c_str(),
	   ',',
	   ';'
	   );
	std::cout << '\n';
	/*  <T_FEATURE,
	    T_CLUSTERIDX, //-1, 0, 1, .., K
	    T_REAL
	    >*/
	clusteringop::reassignCluster
	  (larraymcidx_memberShipTmp,
	   lmatrixrowt_centroidsChromBest,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
	
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
    } /* END if Best*/
      
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep //lpc_labelStep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    /*MEASUREMENT BEST: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
     */
#ifndef __WITHOUT_PLOT_STAT
    if ( aiinp_inParamProbAdaptive.getWithPlotStatObjetiveFunc() ) {  
      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

  } /*END 2. EVALUATE EACH CROMOSOME*/
  
    /*3. If the termination condition is not reached, go to Step 4. Other-
      wise, select the best individual from the population as the best
      cluster result.
    */
  
  while ( llfh_listFuntionHist.getDomainUpperBound() 
	  < aiinp_inParamProbAdaptive.getNumMaxGenerations() ) {

    /*4. Select individuals from the population for crossover and muta-
      tion.
    */
    
    {/*BEGIN SELECTION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "4. SELECTION INDIVIDUAL FROM THE POPULATION FOR CROSSOVER AND MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep //lpc_labelStep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchromfixleng_population.begin(),lvectorchromfixleng_population.end(),
	 [](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* lchromfixleng_iter)
	 -> T_REAL
	 {
	   return lchromfixleng_iter->getFitness();
	 }
	 );
           
      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL--------------------------
       */ 
      for ( auto lchromfixleng_iter: lvectorchromfixleng_matingPool ) {
	
	uintidx lstidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );

	*lchromfixleng_iter = *lvectorchromfixleng_population.at(lstidx_chrom);
	   
      }
      

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END SELECTION
       */

      /*5. Apply crossover operator to the selected individuals based on the
	crossover probability.
	3.4.1 CROSSOVER
      */
    
    {/*BEGIN APPLY CROSSOVER OPERATOR*/

      T_REAL lrt_sumFitness = 
	std::accumulate
	(lvectorchromfixleng_matingPool.begin(),
	 lvectorchromfixleng_matingPool.end(),
	 T_REAL(0.0), //VALUE INITIAL
	 [&](T_REAL lT_sumPartial,
	     const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* lchromfixleng_iter)
	 {
	   return lT_sumPartial + lchromfixleng_iter->getFitness();
	 }
	 ); 
      T_REAL lrt_aveFitness =
	lrt_sumFitness / T_REAL(lvectorchromfixleng_matingPool.size());

      gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchromfixleng_iterMax  =
	*(std::max_element
	  (lvectorchromfixleng_matingPool.begin(), 
	   lvectorchromfixleng_matingPool.end(), 
	   [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	       const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	       ) 
	{  return x->getFitness() < y->getFitness(); } 
	   ));
      T_REAL lrt_maxFitness = lchromfixleng_iterMax->getFitness();

      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "5. APPLY CROSSOVER OPERATOR";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep  
	  << ":  IN(" << geiinparam_verbose << ')'
	  << "\tlrt_maxFitness = " << lrt_maxFitness
	  << "\tlrt_aveFitness = " << lrt_aveFitness 
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/
 
      long ll_invalidOffspring = 0;

      auto ichrom_population = lvectorchromfixleng_population.begin();
      auto jchrom_matingPool = lvectorchromfixleng_matingPool.begin();
     
      if ( ( lvectorchromfixleng_population.size() % 2 ) != 0 ) {
	*(*ichrom_population) = *(*jchrom_matingPool);
	++ichrom_population;
	++jchrom_matingPool;
      }
      while ( (ichrom_population != lvectorchromfixleng_population.end() )
	      && (jchrom_matingPool != lvectorchromfixleng_matingPool.end()) )
	{
	  auto lchrom_child1  = ichrom_population;
	  ++ichrom_population;
	  auto lchrom_child2  = ichrom_population;
	  ++ichrom_population;

	  auto lchrom_parent1  = jchrom_matingPool;
	  ++jchrom_matingPool;
	  auto lchrom_parent2  = jchrom_matingPool;
	  ++jchrom_matingPool;
	  
	  T_REAL lT_fprime =
	    ((*lchrom_parent1)->getFitness() > (*lchrom_parent2)->getFitness())?
	    (*lchrom_parent1)->getFitness():(*lchrom_parent2)->getFitness();
	 
	  T_REAL lT_probCrossover = 
	    prob::adaptiveProb
	    ((T_REAL) 1.0,  /*k1*/
	     (T_REAL) 1.0,  /*k3*/
	     lrt_maxFitness,
	     lrt_aveFitness, 
	     lT_fprime
	     );

	  //if  Crossover
	  if ( uniformdis_real01(gmt19937_eng) < lT_probCrossover ) {
	 
	    uintidx luintidx_GAGRdistM1M2 = 
	      gagenericop::GAGRdist
	      (*(*lchrom_parent1),
	       *(*lchrom_parent2)
	       );
	    
	    //Distance smaller minGAGRdist =  2 pag 1213
	    if ( luintidx_GAGRdistM1M2  < (uintidx) 2 ) { 

	      garealop::heuristicCrossover
		(*(*lchrom_child1),
		 *(*lchrom_child2),
		 *(*lchrom_parent1),
		 *(*lchrom_parent2)
		 );
	       
	      std::pair<T_REAL,bool> lpair_SSE1 = 
		lfuncobjgagr_functionObjetive.getObjetiveFunc
		((*lchrom_child1)->getString());
	      (*lchrom_child1)->setObjetiveFunc(lpair_SSE1.first);
	      (*lchrom_child1)->setFitness
		(lfuncobjgagr_functionObjetive.getFitness(lpair_SSE1.first));
	      (*lchrom_child1)->setValidString(lpair_SSE1.second);
	      
	    }
	    else {
	       
	      gagenericop::pathCrossover
		(*(*lchrom_child1),
		 *(*lchrom_child2),
		 *(*lchrom_parent1),
		 *(*lchrom_parent2),
		 luintidx_GAGRdistM1M2,
		 //uniformdis_uintidx0GAGRdistM1M2(gmt19937_eng),
		 lfuncobjgagr_functionObjetive
		 );
	    }  
	    	    
	    if ( (*lchrom_child1)->getValidString() == false )
	      ++ll_invalidOffspring;
	    if ( (*lchrom_child2)->getValidString() == false )
	      ++ll_invalidOffspring;
	  	   
	  } //if  Crossover
	  else {
	    *(*lchrom_child1) = *(*lchrom_parent1);
	    *(*lchrom_child2) = *(*lchrom_parent2);
	  }
	} /*END  While
	   */
      
      aoop_outParamGAC.sumTotalInvalidOffspring(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END APPLY CROSSOVER OPERATOR*/


      /*6. Apply mutation operator to the selected individuals based on the
	mutation probability.
      */
    {/*BEGIN MUTATION*/

      /*Here,the mutation probability is also selected adaptively as [50]. The
	expression for mutation probability $p_m$ is given below
      */
      T_REAL lrt_sumFitness = 
	std::accumulate
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 T_REAL(0.0), //VALUE INITIAL
	 [&](T_REAL lT_sumPartial,
	     const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* lchromfixleng_iter
	     )
	 {
	   return lT_sumPartial + lchromfixleng_iter->getFitness();
	 }
	 );
      T_REAL lrt_aveFitness =
	lrt_sumFitness / T_REAL(lvectorchromfixleng_population.size());

      gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchrom_maxFitness  =
	*(std::max_element
	  (lvectorchromfixleng_population.begin(), 
	   lvectorchromfixleng_population.end(), 
	   [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	       const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	       ) 
	{  return x->getFitness() < y->getFitness(); } 
	   )
	  );
      T_REAL lrt_maxFitness = lchrom_maxFitness->getFitness();
     

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "6. APPLY MUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ":  IN(" << geiinparam_verbose << ')'
		  << "\tlrt_maxFitness = " << lrt_maxFitness
		  << "\tlrt_aveFitness = " << lrt_aveFitness 
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
           
      /* \cite{Chang:etal:GAclustering:GAGR:2009}
	 The mutation process adopted in this paper is the same as that
	 used in [36] which will be described as below. Let $f_min$ and $f_max$ be
	 the minimum and maximum fitness values in the current population,
	 respectively. For an individual with fitness value $f$ , a number $\delta$ in
	 the range [-R, +R] is generated with uniform distribution, where

	 En el articulo esta mal planteado por que a $R$ lo define en función  de la
	 fitness y debe ser un función de la M como la hace 
	 \cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
      */

      gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchrom_minObjFunc =
	*(std::min_element
	  (lvectorchromfixleng_population.begin(), 
	   lvectorchromfixleng_population.end(), 
	   [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	       const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	       ) 
	{  return x->getObjetiveFunc() < y->getObjetiveFunc(); } 
	   )
	  );
            
      T_REAL lT_minClusteringMetric =
	lchrom_minObjFunc->getObjetiveFunc();

      gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchrom_maxObjFunc =
	*(std::max_element
	  (lvectorchromfixleng_population.begin(), 
	   lvectorchromfixleng_population.end(), 
	   [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	       const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	       ) 
	{  return x->getObjetiveFunc() < y->getObjetiveFunc(); } 
	   )
	  );
          
      T_REAL lT_maxClusteringMetric =
	lchrom_maxObjFunc->getObjetiveFunc();
      
      for ( auto lchromfixleng_iter: lvectorchromfixleng_population ) {
	
	T_REAL lT_probMutation = 
	  prob::adaptiveProb
	  (T_REAL(0.5),  /*k2*/
	   T_REAL(0.5),  /*k4*/
	   lrt_maxFitness,
	   lrt_aveFitness, 
	   lchromfixleng_iter->getFitness()
	   );

	if ( uniformdis_real01(gmt19937_eng) < lT_probMutation )  
	  { //IF BEGIN  MUTATION
	    
	    gaclusteringop::biDirectionHMutation
	      (*lchromfixleng_iter,
	       lT_minClusteringMetric,
	       lT_maxClusteringMetric,
	       larray_minFeactures,
	       larray_maxFeactures
	       );
	    std::pair<T_REAL,bool> lpair_SSE = 
	      lfuncobjgagr_functionObjetive.getObjetiveFunc
	      (lchromfixleng_iter->getString());
	    lchromfixleng_iter->setObjetiveFunc(lpair_SSE.first);
	    lchromfixleng_iter->setFitness
	      (lfuncobjgagr_functionObjetive.getFitness(lpair_SSE.first));
	    lchromfixleng_iter->setValidString(lpair_SSE.second);
	    
	  } //IF END MUTATION
      }
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END MUTATION*/

      /*7. Evaluate the newly generated candidates.
       */
    {/*BEGIN EVALUATE*/

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "7. EVALUATE THE NEWLY GENERATED CANDIDATES:  IN"
		  << '(' << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      /*After applying the operators, the fitness function is evaluated, 
	so it is not necessary to apply this step
       */

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "7. EVALUATE THE NEWLY GENERATED CANDIDATES: OUT"
		  << '(' << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END EVALUATE*/

      /*8. Compare the worst chromosome in the new population with $P_best$
	in term of their fitness values. If the former is worse than the
	later, then replace it by pbest.
      */

    { /*BEGIN REPLACE THE WORST CHROMOSOME BY  Pbest*/ 

      gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchrom_minFitness  =
	*(std::min_element
	  (lvectorchromfixleng_population.begin(), 
	   lvectorchromfixleng_population.end(), 
	   [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	       const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	       ) 
	{  return x->getFitness() < y->getFitness(); } 
	   )
	  );
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "8. REPLACE THE WORST CHROMOSOME BY Pbest";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ":  IN(" << geiinparam_verbose << ')'
		  << "\tlchrom_minFitnes = " << lchrom_minFitness->getFitness()
		  << "\tbest fitness = " << lochromfixleng_best.getFitness()
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
    
      if (  lchrom_minFitness->getFitness() < lochromfixleng_best.getFitness() )  {	
	*lchrom_minFitness = lochromfixleng_best;
      }


#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END REPLACE THE WORST CHROMOSOME BY  Pbest*/ 

    
      /*9. Find the best chromosome in the new population and replace
	Pbest.
      */
    { /*BEGIN FIND THE BEST CHROMOSOME*/

      gaencode::ChromFixedLength<T_FEATURE,T_REAL> *lchrom_maxFitness =
	*(std::max_element
	  (lvectorchromfixleng_population.begin(), 
	   lvectorchromfixleng_population.end(), 
	   [&](const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* x, 
	       const gaencode::ChromFixedLength<T_FEATURE,T_REAL>* y
	       ) 
	{  return x->getFitness() < y->getFitness(); } 
	   )
	  );

#ifdef __VERBOSE_YES
      geverbosepc_labelstep =  "9. FIND THE BEST CHROMOSOME"; 
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep
		  << ": IN(" << geiinparam_verbose << ") max fitness = "
		  << lchrom_maxFitness->getFitness()  
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      if (lochromfixleng_best.getFitness() < lchrom_maxFitness->getFitness() ) {
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	lochromfixleng_best = *lchrom_maxFitness;
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  /*DECODE CHROMOSOME*/
	  mat::MatrixRow<T_FEATURE> 
	    lmatrixrowt_centroidsChromBest
	    (lconstui_numClusterFk,
	     data::Instance<T_FEATURE>::getNumDimensions(),
	     lochromfixleng_best.getString()
	     );
	 
	  std::ostringstream lostrstream_labelCentroids;
	  lostrstream_labelCentroids
	    << "<CENTROIDSCLUSTER:"
	    << geverbosepc_labelstep
	    << ":Chromosome:lochromfixleng_best["
	    <<  geverboseui_idproc << ':' << &lochromfixleng_best
	    << ']';
	  lmatrixrowt_centroidsChromBest.print
	    (std::cout,
	     lostrstream_labelCentroids.str().c_str(),
	     ',',
	     ';'
	     );
	  std::cout << '\n';
	  
	  clusteringop::reassignCluster
	    (larraymcidx_memberShipTmp,
	     lmatrixrowt_centroidsChromBest,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     aifunc2p_dist
	     );
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      }
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
 
    } /*BEGIN FIND THE BEST CHROMOSOME*/

      /*10. For the new population, select the best chromosome as a
	reference, which other chromosomes might fall into the gene
	rearrangement if needed.
      */
    { /*BEGIN REARRANGEMENT*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "10. REARRANGEMENT CHROMOSOME";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep 
		  << ": IN(" << geiinparam_verbose << ")\n";
	std::ostringstream lostrstream_labelChromBest;
	lostrstream_labelChromBest << "BEST";
	lochromfixleng_best.print(std::cout,lostrstream_labelChromBest.str().c_str());
	std::cout	<< std::endl;
      }
#endif /*__VERBOSE_YES*/

      /*DECODE CHROMOSOME*/
      mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroidsBestChrom
	(lconstui_numClusterFk,
	 data::Instance<T_FEATURE>::getNumDimensions(),
	 lochromfixleng_best.getString()
	 );

      for ( auto lchromfixleng_iter: lvectorchromfixleng_population ) {
	    
	gaclusteringop::rearrangedCluster
	  (*lchromfixleng_iter,
	   lmatrixrowt_centroidsBestChrom,
	   aifunc2p_dist
	   );
	   
      }
	
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout <<  geverbosepc_labelstep 
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    } /*END REARRANGEMENT*/


      /*MEASUREMENT BEST: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT  
    if ( aiinp_inParamProbAdaptive.getWithPlotStatObjetiveFunc() ) {

      for ( auto lchromfixleng_iter: lvectorchromfixleng_population ) {
	lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter->getObjetiveFunc());
      }
      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }

#endif /*__WITHOUT_PLOT_STAT*/


    llfh_listFuntionHist.increaseDomainUpperBound();
    
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


  } /*END While*/
  
    /*FREE MEMORY*/

  delete [] larray_maxFeactures;
  delete [] larray_minFeactures;
  
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
  aoop_outParamGAC.setNumClusterK
    (aiinp_inParamProbAdaptive.getNumClusterK());
  aoop_outParamGAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));

  aoop_outParamGAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
 
  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */

#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamProbAdaptive.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamProbAdaptive,
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
    std::cout << '\n';
    
    /*DECODE CHROMOSOME*/
    mat::MatrixRow<T_FEATURE> 
      lmatrixrowt_centroidsChromBest
      (lconstui_numClusterFk,
       data::Instance<T_FEATURE>::getNumDimensions(),
       lochromfixleng_best.getString()
       );
    
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids
      << "<CENTROIDSCLUSTER:"
      << lpc_labelAlgGA
      << ":generation: " <<  llfh_listFuntionHist.getDomainUpperBound()
      << ":Chromosome:lochromfixleng_best["
      << geverboseui_idproc << ':' << &lochromfixleng_best
      << ']';
    lmatrixrowt_centroidsChromBest.print
      (std::cout,
       lostrstream_labelCentroids.str().c_str(),
       ',',
       ';'
       );
    std::cout << '\n';
    /* <T_FEATURE,
       T_CLUSTERIDX, //-1, 0, 1, .., K
       T_REAL
       >*/
    clusteringop::reassignCluster
      (larraymcidx_memberShipTmp,
       lmatrixrowt_centroidsChromBest,
       aiiterator_instfirst,
       aiiterator_instlast,
       aifunc2p_dist
       );
  }
  --geiinparam_verbose;

  delete [] larraymcidx_memberShipTmp;
  
#endif /*__VERBOSE_YES*/

  return lochromfixleng_best;
 
} /* END  gagr_fkcentroid */

} /*END eac */

#endif /*__GAGR_FKCENTROID_HPP__*/
