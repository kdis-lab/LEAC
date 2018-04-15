/*! \file gka_fklabel.hpp
 *
 * \brief GKA Clustering \cite Krishna:Murty:GAClustering:GKA:1999. 
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GKA algorithm based on the paper:\n
 * K. Krishna and M. Narasimha Murty. Genetic k-means algorithm.\n 
 * IEEE Transactions on Systems, Man, and Cybernetics, Part B\n
 * (Cybernetics), 29(3):433–439, Jun 1999.\n
 * doi:http://dx.doi.org/10.1109/3477.764879.\n
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
 

#ifndef __GKA_FKLABEL_KRISHNA_AND_MURTY1999_HPP__
#define __GKA_FKLABEL_KRISHNA_AND_MURTY1999_HPP__

#include <vector>
#include <algorithm>

#include <leac.hpp>
#include "inparam_pmfk.hpp"
#include "outparam_eaclustering.hpp"
#include "ga_integer_operator.hpp"

#include "plot_runtime_function.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of genetic and evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace eac {

/*! \fn  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>  gka_fklabel (inout::OutParamEAClustering<T_REAL, T_CLUSTERIDX> &aoop_outParamEAC, inout::InParamPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE, T_FEATURE_SUM, T_INSTANCES_CLUSTER_K> &aiinp_inParamPmFk, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
 \brief GKA 
 \details Implementation of GKA algorithm based on \cite Krishna:Murty:GAClustering:GKA:1999. 
 \returns a partition of a data set, encoded on a chromosome where each gene is the index of a cluster to which the instance belongs.
 \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
 \param aiinp_inParamPmFk a inout::InParamPmFk parameters required by the algorithm
 \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
 \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
 \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX, //T_CLUSTERIDX,-1,0,1,..,K DATATYPE OF CHROMOSOME
           typename T_REAL,       //T_METRIC, T_FITNESS, T_CLUSTERING_METRIC
	   typename T_FEATURE,    //T_CENTROIDS
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename INPUT_ITERATOR
	   >
gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> 
gka_fklabel
(inout::OutParamEAClustering
 <T_REAL,
 T_CLUSTERIDX>                   &aoop_outParamEAC,  
 inout::InParamPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>          &aiinp_inParamPmFk,
 const INPUT_ITERATOR            aiiterator_instfirst,
 const INPUT_ITERATOR            aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>    &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  /*ID PROC
   */
  geverboseui_idproc = 1;
 
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gka_fklabel";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(\t output inout::OutParamEAClustering&: aoop_outParamEAC[" 
	      << &aoop_outParamEAC << "]\n"
	      << "\t input  InParamPmFk&: aiinp_inParamPmFk[" 
	      << &aiinp_inParamPmFk << "]\n"
	      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist[" << &aifunc2p_dist << "]\n"
	      << "GA parameters: "
	      << "\tPopulation size = " << aiinp_inParamPmFk.getSizePopulation()
	      << "\tProbMutation  = "   << aiinp_inParamPmFk.getProbMutation()
	      << "\tK = " << aiinp_inParamPmFk.getNumClusterK()
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamPmFk.getNumMaxGenerations(), "Iterations", "Clustering metrics");
  
  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;

  if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {  

    lvectorT_statfuncObjetiveFunc.reserve
      (aiinp_inParamPmFk.getSizePopulation());
    //DEFINE FUNCTION
    lofh_SSE  = 
      new runtime::RuntimeFunctionValue<T_REAL>
      ("SSE", 
       aiinp_inParamPmFk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_SSE);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat<T_REAL>
	( (char) li_i,
	  aiinp_inParamPmFk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamEAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamPmFk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamPmFk.getTimesRunAlgorithm()
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
#endif  /*__WITHOUT_PLOT_STAT*/
 
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::setStringSize
    (uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)));
  
  gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> lochromfixleng_best;

  /*POPULATION CREATE-----------------------------------------------------------
   */
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> >  
    lvectorchromfixleng_population(aiinp_inParamPmFk.getSizePopulation());

  /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
   */ 
  std::vector<gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> >  
    lvectorchromfixleng_stringPool(aiinp_inParamPmFk.getSizePopulation());

  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING*/
  aoop_outParamEAC.setTotalInvalidOffspring(0);

  runtime::ExecutionTime let_executionTime = runtime::start();
 
  /*1) CODING:------------------------------------------------------------------ 
   */
  
  /*2) INITIALIZATION: THE INITIAL POPULATION P(0) SELCTED RANDOMY--------------
   */

  { /*BEGIN INITIALIZATION*/
   
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "2) INITIALIZATION: THE INITIAL POPULATION P(0) SELCTED RANDOMY";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_mmcidx0K
      (0,aiinp_inParamPmFk.getNumClusterK()-1);
   
    for (auto &lchromfixleng_iter :lvectorchromfixleng_population) {
      
      gagenericop::initializeGenes
	(lchromfixleng_iter.begin(),
	 lchromfixleng_iter.end(),
	 [&]() 
	 {
	   return uniformdis_mmcidx0K(gmt19937_eng);
	 }
	 );
      
      lchromfixleng_iter.setObjetiveFunc(std::numeric_limits<T_REAL>::max());
      lchromfixleng_iter.setFitness(-std::numeric_limits<T_REAL>::max());
      
    }
     
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END INITIALIZATION*/

 
  /*EVALUTE OBJECTIVE FUNCTION---------------------------------------------------
   */
 
  { /*BEGIN EVALUTE OBJECTIVE FUNCTION*/
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "EVALUTE OBJECTIVE FUNCTION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    long ll_invalidOffspring = 0;

    mat::MatrixRow<T_FEATURE> 
      lmatrixrowt_centroids
      ((uintidx) aiinp_inParamPmFk.getNumClusterK(),
       data::Instance<T_FEATURE>::getNumDimensions() 
       );

    mat::MatrixRow<T_FEATURE_SUM>       
      lmatrixrowt_sumInstCluster
      ((uintidx) aiinp_inParamPmFk.getNumClusterK(), 
       data::Instance<T_FEATURE>::getNumDimensions()
       );
	
    std::vector<T_INSTANCES_CLUSTER_K> 
      lvectort_numInstClusterK
      ((uintidx) aiinp_inParamPmFk.getNumClusterK()
       );
     
    for (auto &lchromfixleng_iter :lvectorchromfixleng_population) {
		
      partition::PartitionLabel
	<T_CLUSTERIDX>
	lpartition_clusters
	(lchromfixleng_iter.getString(),
	 lchromfixleng_iter.getStringSize(),
	 aiinp_inParamPmFk.getNumClusterK()
	 );
     
      T_CLUSTERIDX lmcidx_numClusterNull =
	clusteringop::getCentroids
	(lmatrixrowt_centroids,
	 lmatrixrowt_sumInstCluster,
	 lvectort_numInstClusterK,
	 lpartition_clusters,
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );
    
      T_REAL   lrt_objetiveFunc;
      if ( lmcidx_numClusterNull == 0 ) { 
	lrt_objetiveFunc = 
	  um::SSE
	  (lmatrixrowt_centroids,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   lchromfixleng_iter.getString(),
	   aifunc2p_dist
	   );
      }
      else {
	lrt_objetiveFunc = std::numeric_limits<T_REAL>::max();
      }

      lchromfixleng_iter.setObjetiveFunc(lrt_objetiveFunc);
      lchromfixleng_iter.setValidString(lmcidx_numClusterNull == 0);      

      if ( lchromfixleng_iter.getValidString() == false ) ++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
      lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter.getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::ostringstream lostrstream_labelChrom;
	lostrstream_labelChrom  << "<GRAPH: " << geverbosepc_labelstep;
	lchromfixleng_iter.print(std::cout,lostrstream_labelChrom.str().c_str(),',',';');
	std::cout << std::endl;
      }
#endif /*__VERBOSE_YES*/

    }

    aoop_outParamEAC.sumTotalInvalidOffspring
      (ll_invalidOffspring);

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END EVALUTE OBJECTIVE FUNCTION*/

  /*How are illegal chromosomes handled?
   */

  /*COPY CROMOSOME 0
    s* = P1; (Pi IS THE ith STRING IN P) 
  */
#ifdef __VERBOSE_YES
  geverbosepc_labelstep = "s* = P1; (Pi IS THE ith STRING IN P)";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << geverbosepc_labelstep  
	      << ": IN(" << geiinparam_verbose << ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  lochromfixleng_best  = lvectorchromfixleng_population.at(0); 
 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << geverbosepc_labelstep
	      << ": OUT(" << geiinparam_verbose << ')'
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

   
  /*INITIAL MEASUREMENT: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {
  
    lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
    functionhiststat_evaluateAll
      (lofhs_statObjectiveFunc,
       lvectorT_statfuncObjetiveFunc
       );
    lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
    lvectorT_statfuncObjetiveFunc.clear();
  }
  
#endif  /*__WITHOUT_PLOT_STAT */

  while ( llfh_listFuntionHist.getDomainUpperBound() < 
	  aiinp_inParamPmFk.getNumMaxGenerations()) {
    
    llfh_listFuntionHist.increaseDomainUpperBound();


    /*CALCULATE FITNESS VALUES OF STRINGS IN P:
      fitness value of a solution string s_W depends on the 
      total within-cluster variation S(W). 
    */
   
    { /*BEGIN CALCULATE FITNESS*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "CALCULATE FITNESS VALUES OF STRINGS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
 
      T_REAL lrt_objFuncAverage =
	std::accumulate
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 0, 
	 [&](T_REAL lT_sumPartial,
	     const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& lchromfixleng_iter
	     )
	 {
	   return lT_sumPartial + lchromfixleng_iter.getObjetiveFunc();
	 }
	 );
   
      /*AVERAGE*/
      lrt_objFuncAverage /= (T_REAL ) lvectorchromfixleng_population.size();
      /*STDDESV*/
      T_REAL lrt_objFuncStdDesv =
	std::accumulate
	(lvectorchromfixleng_population.begin(),
	 lvectorchromfixleng_population.end(),
	 0, 
	 [&](T_REAL lT_sumPartial,
	     const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& lchromfixleng_iter
	     )
	 {
	   T_REAL ltr_diff  = lchromfixleng_iter.getObjetiveFunc() - lrt_objFuncAverage;
	 
	   return lT_sumPartial + ltr_diff*ltr_diff;
	 }
	 );
      lrt_objFuncStdDesv /= std::sqrt( (T_REAL ) lvectorchromfixleng_population.size() -1.0);
     
      const T_REAL  lT_varianceTimes = 2;
      T_REAL lT_aveMcstdDesv = lrt_objFuncAverage - lT_varianceTimes * lrt_objFuncStdDesv;
   
      for (auto &lchromfixleng_iter :lvectorchromfixleng_population) {
       
	T_REAL lf_gsw = lchromfixleng_iter.getObjetiveFunc() - lT_aveMcstdDesv;
	lchromfixleng_iter.setFitness((lf_gsw >= 0.0)?lf_gsw:0.0);
      }
    

#ifdef __VERBOSE_YES    
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << "\n<" << geverbosepc_labelstep  << ':' << geverboseui_idproc << '>'
		  << "lrt_objFuncAverage" << ',' << lrt_objFuncAverage
		  << ",lrt_objFuncStdDesv" << ',' << lrt_objFuncStdDesv
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    }/*END CALCULATE FITNESS*/

    /*3) SELECTION--------------------------------------------------------------
      ^P = Selection (P);
    */ 
    {/*BEGIN SELECTION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(3) SELECTION: ^P = Selection(P)";
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
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& lchromfixleng_iter) -> T_REAL
	 {
	   return T_REAL(lchromfixleng_iter.getFitness());
	 }
	 );
      
      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL--------------------------
       */ 
      for (auto &lchromfixleng_iter: lvectorchromfixleng_stringPool) {
       
	uintidx luiidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );
	  
	lchromfixleng_iter = lvectorchromfixleng_population.at(luiidx_chrom);
	
      }


#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    
    }/*END SELECTION*/

    
    /*4) MUTATION
      for i = 1 to N , Pi = Mutation(^Pi);
    */
    { /*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(4) MUTATION for i = 1 to N , Pi = Mutation(^Pi)";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      /* std::vector<T_REAL> lvector_distINSTiCLUSTER1k
	( (uintidx) aiinp_inParamPmFk.getNumClusterK() );
      */
      
      auto lchromfixleng_iterPop_j = lvectorchromfixleng_population.begin();
      auto lchromfixleng_iterPool_i = lvectorchromfixleng_stringPool.begin();


      mat::MatrixRow<T_FEATURE> 
	lmatrixrowt_centroids
	((uintidx) aiinp_inParamPmFk.getNumClusterK(),
	 data::Instance<T_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<T_FEATURE_SUM>       
	lmatrixrowt_sumInstCluster
	((uintidx) aiinp_inParamPmFk.getNumClusterK(), 
	 data::Instance<T_FEATURE>::getNumDimensions(),
	 T_FEATURE_SUM(0)
	 );
	
      std::vector<T_INSTANCES_CLUSTER_K> 
	lvectort_numInstClusterK
	((uintidx) aiinp_inParamPmFk.getNumClusterK(),
	 T_INSTANCES_CLUSTER_K(0)
	 );
		  
      for ( ;lchromfixleng_iterPool_i != lvectorchromfixleng_stringPool.end(); 
	    ++lchromfixleng_iterPool_i, ++lchromfixleng_iterPop_j) 
	{
	  gaintegerop::mutationgka
	    (*lchromfixleng_iterPool_i,
	     aiinp_inParamPmFk.getProbMutation(),
	     lmatrixrowt_centroids,
	     lmatrixrowt_sumInstCluster,
	     lvectort_numInstClusterK,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     aifunc2p_dist
	     );
	  	  
	  *lchromfixleng_iterPop_j = *lchromfixleng_iterPool_i;
	}
     
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


    } /*END MUTATION*/

    /*5) K-Means Operator and evaluate S(Pi)
     */
    {//BEGIN K-MEANS OPERATOR
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "(5) K-Means Operator and evaluate S(Pi)";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      const uintidx  lui_numclusterK =  uintidx(aiinp_inParamPmFk.getNumClusterK());
      
      for (auto &lchromfixleng_iter :lvectorchromfixleng_population) {
	  
	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroids
	  (lui_numclusterK,
	   data::Instance<T_FEATURE>::getNumDimensions() 
	   );

	mat::MatrixRow<T_FEATURE_SUM>       
	  lmatrixrowt_sumInstCluster
	  (lui_numclusterK,
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   T_FEATURE_SUM(0)
	   );
	
	std::vector<T_INSTANCES_CLUSTER_K> 
	  lvectort_numInstClusterK
	  (lui_numclusterK,
	   T_INSTANCES_CLUSTER_K(0)
	   );
 	  
	T_CLUSTERIDX 
	  lmcidx_numClusterNull =
	  clusteringop::kmeansoperator
	  (lchromfixleng_iter.getString(),
	   lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   aiiterator_instfirst,
	   aiiterator_instlast,	     
	   aifunc2p_dist
	   );
	  
	lchromfixleng_iter.setValidString(lmcidx_numClusterNull == 0);

      }

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } //END K-MEANS OPERATOR

   
    /*EVALUTE OBJECTIVE FUNCTION---------------------------------------------------
     */
    { /*BEGIN EVALUTE OBJECTIVE FUNCTION*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "EVALUTE OBJECTIVE FUNCTION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
      
      long ll_invalidOffspring = 0;
      for (auto &lchromfixleng_iter :lvectorchromfixleng_population) {
	  
	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (lchromfixleng_iter.getString(),
	   lchromfixleng_iter.getStringSize(),
	   aiinp_inParamPmFk.getNumClusterK()
	   );

	mat::MatrixRow<T_FEATURE> 
	  lmatrixrowt_centroids
	  ((uintidx) aiinp_inParamPmFk.getNumClusterK(),
	   data::Instance<T_FEATURE>::getNumDimensions() 
	   );

	mat::MatrixRow<T_FEATURE_SUM>       
	  lmatrixrowt_sumInstCluster
	  ((uintidx) aiinp_inParamPmFk.getNumClusterK(), 
	   data::Instance<T_FEATURE>::getNumDimensions(),
	   T_FEATURE_SUM(0)
	   );
	
	std::vector<T_INSTANCES_CLUSTER_K> 
	  lvectort_numInstClusterK
	  ((uintidx) aiinp_inParamPmFk.getNumClusterK(),
	   T_INSTANCES_CLUSTER_K(0)
	   );

	T_CLUSTERIDX lmcidx_numClusterNull =
	  clusteringop::getCentroids
	  (lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   lpartition_clusters,
	   aiiterator_instfirst,
	   aiiterator_instlast
	   );
	
	T_REAL lrt_objetiveFunc;
	if ( lmcidx_numClusterNull == 0 ) { 
	  lrt_objetiveFunc = 
	    um::SSE
	    (lmatrixrowt_centroids,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     lchromfixleng_iter.getString(),
	     aifunc2p_dist
	     );
	}
	else {
	  lrt_objetiveFunc = std::numeric_limits<T_REAL>::max();
	}

	lchromfixleng_iter.setObjetiveFunc(lrt_objetiveFunc);
	lchromfixleng_iter.setValidString(lmcidx_numClusterNull == 0);      

	if ( lchromfixleng_iter.getValidString() == false ) ++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
	lvectorT_statfuncObjetiveFunc.push_back(lchromfixleng_iter.getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

      }
       
      aoop_outParamEAC.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END EVALUTE OBJECTIVE FUNCTION*/

    /*IF (S(Ws*) > S(W)), s* = s
     */
    { /*BEGIN (S(Ws*) > S(W)), s* = s */
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "IF (S(Ws*) > S(W)), s* = s";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
     
      auto lchromfixleng_iterMinObjFun =
	std::min_element
	(lvectorchromfixleng_population.begin(), 
	 lvectorchromfixleng_population.end(), 
	 [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& x, 
	    const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>& y
	    ) 
	 { return x.getObjetiveFunc() < y.getObjetiveFunc(); }
	 );
    
      if ( lchromfixleng_iterMinObjFun->getObjetiveFunc() < lochromfixleng_best.getObjetiveFunc() ) {
	lochromfixleng_best = *lchromfixleng_iterMinObjFun;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	aoop_outParamEAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamEAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));
      }
   
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
    } /*END (S(Ws*) > S(W)), s* = s */

    /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT
    if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {
  
      lofh_SSE->setValue(lochromfixleng_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectorT_statfuncObjetiveFunc.clear();
    }
#endif  /*__WITHOUT_PLOT_STAT */

#ifdef __VERBOSE_YES
   
    /*ID PROC
     */
    ++geverboseui_idproc;
      
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ITERATION: OUT" 
		<< '(' << geiinparam_verbose << ')'
		<< "\titeration = " << llfh_listFuntionHist.getDomainUpperBound()
		<< "\tobjetivoFunc = " << lochromfixleng_best.getObjetiveFunc() 
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }/*END While*/

  runtime::stop(let_executionTime);
  aoop_outParamEAC.setNumClusterK
    (aiinp_inParamPmFk.getNumClusterK());
  aoop_outParamEAC.setMetricFuncRun
    (lochromfixleng_best.getObjetiveFunc());
  aoop_outParamEAC.setFitness
    (lochromfixleng_best.getFitness());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
  
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamPmFk,
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

} /*END gka_fklabel*/

} /*END eac */

#endif /*__GKA_FKLABEL_KRISHNA_AND_MURTY1999_HPP__*/
