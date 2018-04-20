/*! \file igka_fklabel.hpp
 *
 * \brief IGKA \cite Lu:etal:GAclusteringLabel:IGKA:2004 and FGKA \cite Lu:etal:GAclusteringLabel:FGKA:2004 
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the IGKA and FGKA algorithm based on the papers:\n
 * Yi Lu, Shiyong Lu, Farshad Fotouhi, Youping Deng, and Susan J.\n
 * Brown.Incremental genetic k-means algorithm and its application\n
 * in gene expression data analysis. BMC Bioinformatics, 5:172, 2004.\n
 * \n
 * Yi Lu, Shiyong Lu, Farshad Fotouhi, Youping Deng, and Susan J. Brown.\n 
 * Fgka: a fast genetic k-means clustering algorithm. In Proceedings\n
 * of the 2004 ACM symposium on Applied computing, SAC ’04,\n
 * pages 622–623, New York, NY, USA, 2004. ACM.\n
 * <a href="http://doi.acm.org/10.1145/967900.968029">doi:http://doi.acm.org/10.1145/967900.968029</a>.\n
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

#ifndef __IGKA_FKLABEL_LU_ETA2004__
#define __IGKA_FKLABEL_LU_ETA2004__

#include <vector>
#include <algorithm>
#include <utility>      // std::pair
#include <limits>       // std::numeric_limits

#include <leac.hpp>
#include "outparam_gac.hpp"
#include "inparam_pmfk.hpp"

#include "plot_runtime_function.hpp"


/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace eac {

/*! \fn  gaencode::ChromosomeIGKA <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> igka_fklabel(inout::OutParamGAC<T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamPmFk, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief IGKA \cite Lu:etal:GAclusteringLabel:IGKA:2004 and FGKA \cite Lu:etal:GAclusteringLabel:FGKA:2004 
  \details 
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinp_inParamPmFk a inout::InParamPmFk parameters required by the algorithm
  \param aivectorptinst_instances a std::vector<data::Instance<T_FEATURE>* >
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,
           typename T_REAL,
	   typename T_FEATURE, //CENTROIDS
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
gaencode::ChromosomeIGKA 
<T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
igka_fklabel
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                   &aoop_outParamGAC, 
 inout::InParamPmFk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>          &aiinp_inParamPmFk,
 std::vector
 <data::Instance<T_FEATURE>* >   &aivectorptinst_instances,
 dist::Dist<T_REAL,T_FEATURE>    &aifunc2p_dist
 )
{

#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
  
  gaencode::ChromosomeIGKA
    <T_CLUSTERIDX,
     T_REAL,
     T_FEATURE,
     T_FEATURE_SUM,
     T_INSTANCES_CLUSTER_K>
    ::setStringSize
    ((uintidx)aivectorptinst_instances.size());
 
  gaencode::ChromosomeIGKA
    <T_CLUSTERIDX,
     T_REAL,
     T_FEATURE,
     T_FEATURE_SUM,
     T_INSTANCES_CLUSTER_K
     >
    lochromigka_best
    ((uintidx) aiinp_inParamPmFk.getNumClusterK()
     );

  std::vector
    <gaencode::ChromosomeIGKA
     <T_CLUSTERIDX,
      T_REAL,
      T_FEATURE,
      T_FEATURE_SUM,
      T_INSTANCES_CLUSTER_K
      >* 
     >
    lvectorchromigka_populationS;
 
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
 
#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004

  gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>::setStringSize
    ((uintidx) aivectorptinst_instances.size());

  gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>::setNumClusterK
    (aiinp_inParamPmFk.getNumClusterK() );

  gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL> lochromigka_best;

  std::vector<gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>* >  
    lvectorchromigka_populationS;

  
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/

  std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_mmcidx0K
    (0,aiinp_inParamPmFk.getNumClusterK()-1);
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
    
  T_REAL lrt_TWCVMax = 0.0;

#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "igka_fklabel";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\t\n(output gaencode::ChromosomeIGKA: lochromigka_best[" 
	      << &lochromigka_best << "]\n"
	      << "\t output inout::OutParamGAC&: aoop_outParamGAC[" 
	      << &aoop_outParamGAC << "]\n"
	      << "\t input  InParamPmFk&: aiinp_inParamPmFk[" 
	      << &aiinp_inParamPmFk << "]\n"
	      << "\t input  std::vector<Instance>: aivectorptinst_instances[" 
	      <<  &aivectorptinst_instances << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\nGA parameters: "
	      << "\tPopulation size = " << aiinp_inParamPmFk.getSizePopulation()
	      << "\tProbMutation  = " << aiinp_inParamPmFk.getProbMutation()
	      << "\tKFind = " << aiinp_inParamPmFk.getNumClusterK()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamPmFk.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT

  std::ofstream lfileout_plotStatObjetiveFunc;
  //char          lc_plotStatSeparator = '\t';
  runtime::RuntimeFunctionValue<T_REAL> *lofh_SSE = NULL;
  runtime::RuntimeFunctionValue<runtime::ExecutionTime> *lofh_time = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>                   lvectorT_statfuncObjetiveFunc;
  //runtime::ExecutionTime               lexetime_iteration;

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
    lofh_time  = 
      new runtime::RuntimeFunctionValue<runtime::ExecutionTime>
      ("time", 
       aiinp_inParamPmFk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );
    llfh_listFuntionHist.addFuntion(lofh_SSE);
    llfh_listFuntionHist.addFuntion(lofh_time);

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
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamPmFk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamPmFk.getTimesRunAlgorithm()
       );

    lfileout_plotStatObjetiveFunc.open
      (aoop_outParamGAC.getFileNameOutPlotStatObjetiveFunc().c_str(),  
       std::ios::out | std::ios::app
       );
    lfileout_plotStatObjetiveFunc.precision(COMMON_COUT_PRECISION);

    //FUNCTION HEADER
    lfileout_plotStatObjetiveFunc 
      <<  llfh_listFuntionHist.getHeaderFuntions() 
      << '\n';
  }

#endif  /*__WITHOUT_PLOT_STAT*/

  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING*/
  aoop_outParamGAC.setTotalInvalidOffspring(0);

  runtime::ExecutionTime lexetime_time = runtime::start();
 
  /*POPULATION CREATE AND INITIALIZE------------------------------------------------------------
   */

  { /*BEGIN INITIALIZATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "INITIALIZE POPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    long ll_invalidOffspring = 0;
 
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
    const uintidx luintidx_numClusterK =
      (uintidx) aiinp_inParamPmFk.getNumClusterK();
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
   
    lvectorchromigka_populationS.reserve
      (aiinp_inParamPmFk.getSizePopulation());
   
    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinp_inParamPmFk.getSizePopulation(); 
	 luintidx_i++) 
      {

#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
       
	lvectorchromigka_populationS.push_back
	  (new
	   gaencode::ChromosomeIGKA
	   <T_CLUSTERIDX,
	   T_REAL,
	   T_FEATURE,
	   T_FEATURE_SUM,
	   T_INSTANCES_CLUSTER_K
	   >
	   (luintidx_numClusterK)
	   );

#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004

	lvectorchromigka_populationS.push_back
	  ( new gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>() );
       
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
       
      }
    for (auto lchromigka_iter: lvectorchromigka_populationS) {
      
      gagenericop::initializeGenes
	(lchromigka_iter->begin(),
	 lchromigka_iter->end(),
	 [&]() 
	 {
	   return uniformdis_mmcidx0K(gmt19937_eng);
	 }
	 );

      //EVALUATE ObjetiveFunc
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
      lchromigka_iter->initialize
	(aivectorptinst_instances,
	 aifunc2p_dist
	 );
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/


#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
      
      partition::PartitionLabel
	<T_CLUSTERIDX>
	lpartition_clusters
	(lchromigka_iter->getString(),
	 lchromigka_iter->getStringSize(),
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
    
      /*TRANSFORM OF LABEL TO CENTROID
       */
      T_CLUSTERIDX lmcidx_numClusterNull =
	clusteringop::getCentroids
	(lmatrixrowt_centroids,
	 lmatrixrowt_sumInstCluster,
	 lvectort_numInstClusterK,
	 lpartition_clusters,
	 aivectorptinst_instances.begin(),
	 aivectorptinst_instances.end()
	 );
	
      T_REAL   lrt_objetiveFunc =
	um::SSE
	(lmatrixrowt_centroids,
	 aivectorptinst_instances.begin(),
	 aivectorptinst_instances.end(),
	 lchromigka_iter->getString(),
	 aifunc2p_dist
	 );
      lchromigka_iter->setObjetiveFunc(lrt_objetiveFunc);
      lchromigka_iter->setNumClusterNotNull
	(aiinp_inParamPmFk.getNumClusterK()- lmcidx_numClusterNull);
      lchromigka_iter->setValidString( lmcidx_numClusterNull == 0);

#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/      
     
      if ( lrt_TWCVMax <  lchromigka_iter->getObjetiveFunc() )
	lrt_TWCVMax =  lchromigka_iter->getObjetiveFunc();

      if ( lchromigka_iter->getValidString() == false )
	++ll_invalidOffspring;
	
#ifndef __WITHOUT_PLOT_STAT
      lvectorT_statfuncObjetiveFunc.push_back(lchromigka_iter->getObjetiveFunc());
#endif //__WITHOUT_PLOT_STAT
      
      lchromigka_iter->setFitness(-std::numeric_limits<T_REAL>::max());
      
    }
   
    aoop_outParamGAC.sumTotalInvalidOffspring
      (ll_invalidOffspring);
  
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')' << " lrt_TWCVMax = " << lrt_TWCVMax
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  } /*END INITIALIZATION*/


  /* s* = P1; (Pi IS THE ith STRING IN P)---------------------------------------- 
     COPY CROMOSOME 0
  */
 
  { //BEGIN s* = P1; (Pi IS THE ith STRING IN P)

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "s* = P1; (Pi IS THE ith STRING IN P)";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif //__VERBOSE_YES
    

    lochromigka_best = *lvectorchromigka_populationS.at(0);


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< " lochromigka_best.getFitness() = " << lochromigka_best.getFitness()
		<< " lochromigka_best.getObjetiveFunc()  = " << lochromigka_best.getObjetiveFunc() 
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  } //END s* = P1; (Pi IS THE ith STRING IN P)
 
  /*INITIAL MEASUREMENT: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {
   
    lofh_SSE->setValue(lochromigka_best.getObjetiveFunc());
    lofh_time->setValue(runtime::elapsedTime(lexetime_time));
   
    functionhiststat_evaluateAll
      (lofhs_statObjectiveFunc,
       lvectorT_statfuncObjetiveFunc
       );

    lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      
    lvectorT_statfuncObjetiveFunc.clear();
    //lexetime_iteration = runtime::start();
  }
  
#endif  /*__WITHOUT_PLOT_STAT */

  /*BEGIN EVOLUTIONARY LOOP------------------------------------------------------ 
   */
  while ( llfh_listFuntionHist.getDomainUpperBound() <  aiinp_inParamPmFk.getNumMaxGenerations()) {
    
    llfh_listFuntionHist.increaseDomainUpperBound();
   
     
    /*SELECTION---------------------------------------------------------------- 
     */ 
    { /*BEGIN CALCULATE FITNESS*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "CALCULATE FITNESS VALUES OF STRINGS IN P";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
    
      T_REAL  lT_Fmin = std::numeric_limits<T_REAL>::max();
      bool    lb_flagExistFmin = false;

      for ( auto lchromigka_iter: lvectorchromigka_populationS )  {
	if ( lchromigka_iter->getLegalityRation() == 1.0 ) {
	  lb_flagExistFmin = true;
	  T_REAL lT_Fitness = 1.5 * lrt_TWCVMax - lchromigka_iter->getObjetiveFunc();
	  lchromigka_iter->setFitness(lT_Fitness);
	  if ( lT_Fmin > lT_Fitness ) {
	    lT_Fmin = lT_Fitness;
	  } 
	}
      }

      if (lb_flagExistFmin == false )
	lT_Fmin = 1.0;

      for ( auto lchromigka_iter: lvectorchromigka_populationS )  {
	if ( lchromigka_iter->getLegalityRation() < 1.0 ) {
	  T_REAL lT_Fitness = lchromigka_iter->getLegalityRation() * lT_Fmin;
	  lchromigka_iter->setFitness(lT_Fitness);
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

    }/*END CALCULATE FITNESS*/
 
    {/*BEGIN SELECTION  P = Selection(P)*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "P = SELECTION(P)";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
       */
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
    
      std::vector
	<gaencode::ChromosomeIGKA
	 <T_CLUSTERIDX,
	  T_REAL,
	  T_FEATURE,
	  T_FEATURE_SUM,
	  T_INSTANCES_CLUSTER_K
	  >*
	 >   
	lvectorchrombase_populationSprime;
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/


#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004

      std::vector<gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>* >  
	lvectorchrombase_populationSprime;
     
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
     
      lvectorchrombase_populationSprime.reserve(lvectorchromigka_populationS.size());
     
      for ( auto lchromigka_iter: lvectorchromigka_populationS )  {
	/*COPY REFERENCES S TO Sprime*/
	lvectorchrombase_populationSprime.push_back(lchromigka_iter);
      }
     
      const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorchromigka_populationS.begin(),
	 lvectorchromigka_populationS.end(),
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
	 [](const gaencode::ChromosomeIGKA
	    <T_CLUSTERIDX,
	    T_REAL,
	    T_FEATURE,
	    T_FEATURE_SUM,
	    T_INSTANCES_CLUSTER_K
	    >* iter_chrom) -> T_REAL
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
	 [](const gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>*  iter_chrom) -> T_REAL
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
	 {
	   return T_REAL(iter_chrom->getFitness());
	 }
	 );

     
      /*COPY POPULATION TO STRING POOL FOR ROULETTE WHEEL AND BIG CHROMOSOME---
       */ 
      std::vector<std::pair<uintidx,uintidx> > lvectorpairst_match;

      for (uintidx luintidx_i = 0; luintidx_i < lvectorchromigka_populationS.size(); luintidx_i++) {

	uintidx luintidx_chrom = 
	  gaselect::getIdxRouletteWheel
	  (lvectorT_probDistRouletteWheel,
	   uintidx(0)
	   );
 
	auto  lochromfgka_sPrime =
	  lvectorchrombase_populationSprime.at(luintidx_chrom);
	 
	if ( lochromfgka_sPrime->getSelected() == false ) {
	  lvectorchromigka_populationS.at(luintidx_i) = lochromfgka_sPrime;
	  lochromfgka_sPrime->setSelected(true);
	}
	else {
	  lvectorpairst_match.push_back(std::make_pair(luintidx_i,luintidx_chrom));
	}
      } /*END for luintidx_i */
     
      uintidx luintidx_j = 0;
      for (uintidx luintidx_i = 0; luintidx_i < lvectorchromigka_populationS.size(); luintidx_i++) {

	auto lochromfgka_sPrime_i =
	  lvectorchrombase_populationSprime.at(luintidx_i);
	 
	if (lochromfgka_sPrime_i->getSelected() == false && 
	    luintidx_j < lvectorpairst_match.size() ) {
	 
	  auto 
	    lochromfgka_sPrime_j =
	    lvectorchrombase_populationSprime.at(lvectorpairst_match.at(luintidx_j).second);;
	         
	  *lochromfgka_sPrime_i  = *lochromfgka_sPrime_j;
	  lvectorchromigka_populationS.at(lvectorpairst_match.at(luintidx_j).first) 
	    =  lochromfgka_sPrime_i; 
	  lvectorchromigka_populationS.at(lvectorpairst_match.at(luintidx_j).first)->setSelected(false);
	  ++luintidx_j;
	}
	lochromfgka_sPrime_i->setSelected(false);
      
      }

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')';
	std::cout << "\nlvectorchromigka_populationS [" << &lvectorchromigka_populationS << "]:\n" ;
	for ( auto  liter_iChrom: lvectorchromigka_populationS ) {
	  liter_iChrom->print();
	  std::cout << std::endl;
	}
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
                  
    } /*END SELECTION*/


    /*MUTATION ----------------------------------------------------------------
      for i = 1 to N , Pi = Mutation(^Pi);
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
      T_INSTANCES_CLUSTER_K linstck_numAllelsMutated;
#endif /*__VERBOSE_YES*/
          
      std::vector<T_REAL> lvectorrt_probDistDistK(aiinp_inParamPmFk.getNumClusterK());
     
      for ( auto lchromigka_iter: lvectorchromigka_populationS )  {

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "chromosome mutation:  IN"
		    << '(' << geiinparam_verbose << ')'
		    << std::endl;
	}
	linstck_numAllelsMutated = 0;
#endif /*__VERBOSE_YES*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004

	/*Calculate cluster centers,cj's corresponding to sw:*/
	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (lchromigka_iter->getString(),
	   gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize(),
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
    
	/*TRANSFORM OF LABEL TO CENTROID
	 */
	clusteringop::getCentroids
	  (lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   lpartition_clusters,
	   aivectorptinst_instances.begin(),
	   aivectorptinst_instances.end()
	   );
	 
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
	  
	for ( uintidx luintidx_i = 0; luintidx_i < aivectorptinst_instances.size(); luintidx_i++) {
	    
	  if ( uniformdis_real01(gmt19937_eng) < aiinp_inParamPmFk.getProbMutation() ) {
	    /*IF BEGIN PROBABILITY*/
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
	      
	    for ( uintidx lui_k = 1; 
		  lui_k < lvectorrt_probDistDistK.size();
		  ++lui_k) 
	      {
		lvectorrt_probDistDistK[lui_k] = lvectorrt_probDistDistK[lui_k-1]
		  +  
		  ( (1.5 * lchromigka_iter->getDistMaxINSTiCLUSTEk(luintidx_i) - 
		     lchromigka_iter->getDistINSTiCLUSTE1k(luintidx_i,lui_k) + 0.5) 
		    / lchromigka_iter->getDistSumINSTiCLUSTE1k(luintidx_i));   
	      }
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
		
#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
		
	    std::vector<T_REAL>&& lvectorrt_distK = 
	      nearest::dist
	      (lmatrixrowt_centroids,
	       lvectort_numInstClusterK,
	       aivectorptinst_instances.at(luintidx_i)->getFeatures(),
	       aifunc2p_dist
	       );
	       
	    auto lit_dmax = 
	      std::max_element
	      (lvectorrt_distK.begin(),
	       lvectorrt_distK.end()
	       );

	    T_REAL lrt_dmax = std::distance(lvectorrt_distK.begin(), lit_dmax);
		
	    std::vector<T_REAL>&& lvectorrt_probDistDistK =
	      prob::makeDistRouletteWheel
	      (lvectorrt_distK.begin(),
	       lvectorrt_distK.end(),
	       [&](T_REAL airt_dXnck )
	       {
		 return 1.5 * lrt_dmax - airt_dXnck +0.5;
	       }
	       );
		
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
		
	    T_CLUSTERIDX lmcidx_oldAllen = 
	      lchromigka_iter->getGene(luintidx_i);
		
	    T_CLUSTERIDX lmcidx_newAllen =
	      gaselect::getIdxRouletteWheel
	      (lvectorrt_probDistDistK,
	       T_CLUSTERIDX(0)
	       );
		
	    if ( lmcidx_oldAllen != lmcidx_newAllen ) {
	      lchromigka_iter->setGene(luintidx_i,lmcidx_newAllen);
		  
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
		  
	      lchromigka_iter->accumulateUpdate
		(lmcidx_newAllen, 
		 lmcidx_oldAllen,
		 aivectorptinst_instances.at(luintidx_i)
		 );
	     
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
		  
	      --lvectort_numInstClusterK[lmcidx_oldAllen];
	      ++lvectort_numInstClusterK[lmcidx_newAllen];

	      interfacesse::axpy
		(lmatrixrowt_sumInstCluster.getRow(lmcidx_oldAllen),
		 (T_FEATURE) -1,
		 aivectorptinst_instances.at(luintidx_i)->getFeatures(),
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );

	      interfacesse::axpy
		(lmatrixrowt_sumInstCluster.getRow(lmcidx_newAllen),
		 (T_FEATURE) 1,
		 aivectorptinst_instances.at(luintidx_i)->getFeatures(), 
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );
	      T_CLUSTERIDX  lmcidx_numClusterNull;

	      clusteringop::meanCentroids
		(lmcidx_numClusterNull,
		 lmatrixrowt_centroids, 
		 lmatrixrowt_sumInstCluster, 
		 lvectort_numInstClusterK
		 );
    
	      lmcidx_numClusterNull =
		(T_CLUSTERIDX)
		std::count_if
		(lvectort_numInstClusterK.begin(),
		 lvectort_numInstClusterK.end(),
		 [] (const uintidx aiit_num) {return aiit_num == 0;}
		 );	  
	      lchromigka_iter->setNumClusterNotNull
		(aiinp_inParamPmFk.getNumClusterK()- lmcidx_numClusterNull);
	      lchromigka_iter->setValidString( lmcidx_numClusterNull == 0);
	     
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
		    
#ifdef __VERBOSE_YES
	      ++linstck_numAllelsMutated;
#endif /*__VERBOSE_YES*/
	      /*AccumulateUpdate(m, n, olda, S[m].a[n]);
	       */
	    }
	  } /*IF BEGIN PROBABILITY*/
	    
	} /*FOR*/
	  
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
	  
	lchromigka_iter->incrementalUpdate
	  (aivectorptinst_instances,
	   aifunc2p_dist
	   );

	if ( lrt_TWCVMax <  lchromigka_iter->getObjetiveFunc() )
	  lrt_TWCVMax =  lchromigka_iter->getObjetiveFunc();
	  
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
	    
#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "chromosome mutation: OUT"
		    << '(' << geiinparam_verbose << ")\n";
	  lchromigka_iter->print();
	  float lf_percentAllelsMutated 
	    = linstck_numAllelsMutated / (float) lchromigka_iter->getStringSize();  
	  std::cout << "allels changed: " 
		    << linstck_numAllelsMutated  
		    << ", percent: " << lf_percentAllelsMutated
		    << '\n';
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
    
       
    } /*END MUTATION*/

    { /*BEGIN K-Means Operator
       */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "K-Means Operator";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      long ll_invalidOffspring = 0;
      
      for ( auto lchromigka_iter: lvectorchromigka_populationS )  {
      
	if ( lchromigka_iter->getNumClusterNotNull() ==  aiinp_inParamPmFk.getNumClusterK() ) {

#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
	  
	  for ( uintidx luintidx_i = 0; luintidx_i < aivectorptinst_instances.size(); luintidx_i++) {
	    T_CLUSTERIDX lmcidx_oldAllen = 
	      lchromigka_iter->getGene(luintidx_i);
	    T_CLUSTERIDX lmcidx_minAllen = 
	      lchromigka_iter->getIdxMinINSTiCLUSTEk(luintidx_i);
	      
	    if (lmcidx_minAllen != lmcidx_oldAllen ) {
	      lchromigka_iter->setGene(luintidx_i,lmcidx_minAllen);
	      lchromigka_iter->accumulateUpdate
		(lmcidx_minAllen, 
		 lmcidx_oldAllen, 
		 aivectorptinst_instances.at(luintidx_i)
		 );
	    }
	  }
	  lchromigka_iter->incrementalUpdate
	    (aivectorptinst_instances,
	     aifunc2p_dist
	     ); 
	
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004

	  const uintidx  lui_numclusterK =  uintidx(aiinp_inParamPmFk.getNumClusterK());
	  
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
 	  
	  clusteringop::kmeansoperator
	    (lchromigka_iter->getString(),
	     lmatrixrowt_centroids,
	     lmatrixrowt_sumInstCluster,
	     lvectort_numInstClusterK,
	     aivectorptinst_instances.begin(),
	     aivectorptinst_instances.end(),	     
	     aifunc2p_dist
	     );
	
	
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/

#ifdef __VERBOSE_YES
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    lchromigka_iter->print();
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

	}//IF
	else {
	  ++ll_invalidOffspring;
	}
	
      }//FOR 

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
    
    }  /*END K-Means Operator
	*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004

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

      const uintidx  lui_numclusterK =  uintidx(aiinp_inParamPmFk.getNumClusterK());
     
      for (auto lchromigka_iter: lvectorchromigka_populationS) {

	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (lchromigka_iter->getString(),
	   lchromigka_iter->getStringSize(),
	   aiinp_inParamPmFk.getNumClusterK() 
	   );

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
    
	/*TRANSFORM OF LABEL TO CENTROID
	 */
	T_CLUSTERIDX lmcidx_numClusterNull =
	  clusteringop::getCentroids
	  (lmatrixrowt_centroids,
	   lmatrixrowt_sumInstCluster,
	   lvectort_numInstClusterK,
	   lpartition_clusters,
	   aivectorptinst_instances.begin(),
	   aivectorptinst_instances.end()
	   );
	 
	T_REAL   lrt_objetiveFunc =
	  um::SSE
	  (lmatrixrowt_centroids,
	   aivectorptinst_instances.begin(),
	   aivectorptinst_instances.end(),
	   lchromigka_iter->getString(),
	   aifunc2p_dist
	   );
	lchromigka_iter->setObjetiveFunc(lrt_objetiveFunc);
	lchromigka_iter->setNumClusterNotNull
	  (aiinp_inParamPmFk.getNumClusterK()- lmcidx_numClusterNull);
	lchromigka_iter->setValidString( lmcidx_numClusterNull == 0);
     
	if ( lrt_TWCVMax <  lchromigka_iter->getObjetiveFunc() )
	  lrt_TWCVMax =  lchromigka_iter->getObjetiveFunc();

	if ( lchromigka_iter->getValidString() == false )
	  ++ll_invalidOffspring;
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

    } /*END EVALUTE OBJECTIVE FUNCTION*/

   
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/

     
    /*IF (S(Ws*) > S(W)), s* = s  --------------------------------------------
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
      
      auto   *loch_minChromosomeIGKA  =
	*(std::min_element
	  (lvectorchromigka_populationS.begin(),
	   lvectorchromigka_populationS.end(),
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
	   [](const gaencode::ChromosomeIGKA
	      <T_CLUSTERIDX,
	      T_REAL,
	      T_FEATURE,
	      T_FEATURE_SUM,
	      T_INSTANCES_CLUSTER_K>* x,
	      const gaencode::ChromosomeIGKA
	      <T_CLUSTERIDX,
	      T_REAL,
	      T_FEATURE,
	      T_FEATURE_SUM,
	      T_INSTANCES_CLUSTER_K>* y
	      )
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/
#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
	   [](const gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>* x,
	      const gaencode::ChromosomeFGKA<T_CLUSTERIDX,T_REAL>* y
	      )
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/
	{  return x->getObjetiveFunc() < y->getObjetiveFunc(); }
	   ));
      if ( loch_minChromosomeIGKA->getObjetiveFunc() < lochromigka_best.getObjetiveFunc() )
	{/*BEGIN IF FIN std::min_element*/
	  lochromigka_best = *loch_minChromosomeIGKA;
	  /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	  aoop_outParamGAC.setIterationGetsBest
	    (llfh_listFuntionHist.getDomainUpperBound());
	  aoop_outParamGAC.setRunTimeGetsBest
	    (runtime::elapsedTime(lexetime_time));

	} /*END IF FIN std::min_element*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << " loch_minChromosomeIGKA->getFitness() = " << loch_minChromosomeIGKA->getFitness()
		  << " loch_minChromosomeIGKA->getObjetiveFunc() = " << loch_minChromosomeIGKA->getObjetiveFunc()
		  << " lochromigka_best.getFitness() = " << lochromigka_best.getFitness()
		  << " lochromigka_best.getObjetiveFunc()  = " << lochromigka_best.getObjetiveFunc() 
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    

    } /*END (S(Ws*) > S(W)), s* = s */

#ifdef __VERBOSE_YES

    /*ID PROC
     */
    ++geverboseui_idproc;
    
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "END ITERATION: " << llfh_listFuntionHist.getDomainUpperBound()
		<< "\tobjetivoFunc = " << lochromigka_best.getObjetiveFunc() 
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT  

    if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {
      //runtime::stop(lexetime_iteration);

   
      for ( auto lchromigka_iter: lvectorchromigka_populationS )  {
	lvectorT_statfuncObjetiveFunc.push_back(lchromigka_iter->getObjetiveFunc());
      }
   
      lofh_SSE->setValue(lochromigka_best.getObjetiveFunc());
      lofh_time->setValue(runtime::elapsedTime(lexetime_time)); 
      //lofh_time->setValue(lexetime_iteration);

      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectorT_statfuncObjetiveFunc
	 );

      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
  
      lvectorT_statfuncObjetiveFunc.clear();
      //lexetime_iteration = runtime::start();
    }

#endif /*__WITHOUT_PLOT_STAT */


  } /*END While*/
  /*END EVOLUTIONARY LOOP------------------------------------------------------ 
   */

  /*FREE MEMORY*/

#ifdef __VERBOSE_YES
  geverbosepc_labelstep = "DELETEPOPULATION";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout <<  geverbosepc_labelstep
	      << ":  IN(" << geiinparam_verbose << ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  for (uintidx lui_i = 0; lui_i < lvectorchromigka_populationS.size(); ++lui_i) {
    delete lvectorchromigka_populationS[lui_i];
  }
 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << geverbosepc_labelstep
	      << ": OUT(" << geiinparam_verbose << ')'
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
 
  runtime::stop(lexetime_time);
  aoop_outParamGAC.setNumClusterK
    (aiinp_inParamPmFk.getNumClusterK());
  aoop_outParamGAC.setMetricFuncRun
    ( lochromigka_best.getObjetiveFunc() );
  aoop_outParamGAC.setFitness
    (lochromigka_best.getFitness());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(lexetime_time));
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());
   
#ifndef __WITHOUT_PLOT_STAT
  if ( aiinp_inParamPmFk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamPmFk,
       aoop_outParamGAC
       );  
  }
#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA 
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochromigka_best.print();
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

   
  return lochromigka_best;
}

} /*END eac */

#endif /*__IGKA_FKLABEL_LU_ETA2004__*/
