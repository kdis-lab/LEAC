/*! \file gga_vklabel.hpp
 *
 * \brief GGA \cite Agustin:etal:GAclusteringVarK:GGA:2012
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GGA algorithm based on the paper:\n
 * L. E. Agustın-Blas, S. Salcedo-Sanz, S. Jimenez-Fernandez,\n 
 * L. Carro-Calvo, J. Del Ser, and J. A. Portilla-Figueras.\n
 * A new grouping genetic algorithm for clustering problems.\n
 * Expert Syst. Appl., 39(10):9695–9703, August 2012.\n
 * <a href="http://dx.doi.org/10.1016/j.eswa.2012.02.149">doi:http://dx.doi.org/10.1016/j.eswa.2012.02.149</a>.\n
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

#ifndef __GGA_VKLABEL_HPP__
#define __GGA_VKLABEL_HPP__

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <assert.h> 
#include <utility>      // std::move

#include <leac.hpp>
#include "inparam_gga.hpp"
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

  
/*! \fn gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> gga_vklabel(inout::OutParamEAClustering<T_REAL,T_CLUSTERIDX> &aoop_outParamEAC,inout::InParamGGA<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamGGA, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief GGA \cite Agustin:etal:GAclusteringVarK:GGA:2012
  \details Implementation of the GGA algorithm based on \cite Agustin:etal:GAclusteringVarK:GGA:2012. Which automatically finds K cluster using the Silhouette and  Davies–Bouldin index.
  \returns A partition of a data set, encoded on a chromosome where each gene is the index of a cluster to which the instance belongs.
  \param aoop_outParamEAC a inout::OutParamEAClustering with the output parameters of the algorithm
  \param aiinp_inParamGGA a inout::InParamGGA parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX, //DATATYPE OF CHROMOSOME
           typename T_REAL,       //T_METRIC, T_FITNESS,
	   typename T_FEATURE,    //T_CENTROIDS
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename INPUT_ITERATOR
	   >
gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> 
gga_vklabel
(inout::OutParamEAClustering
 <T_REAL,
 T_CLUSTERIDX>                        &aoop_outParamEAC,
 inout::InParamGGA
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>               &aiinp_inParamGGA,
 const INPUT_ITERATOR                 aiiterator_instfirst,
 const INPUT_ITERATOR                 aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>   &aifunc2p_dist
 )
{
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  if ( aiinp_inParamGGA.getNumClusterKMaximum() == 
       INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED )
    aiinp_inParamGGA.setNumClusterKMaximum
      (std::round(std::sqrt((double)lui_numInstances)));

  gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>
    ::setElementSize((uintidx) lui_numInstances);

  gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> lochrom_best
    ( uintidx
      ((aiinp_inParamGGA.getNumClusterKMaximum() +
	aiinp_inParamGGA.getNumClusterKMinimum() )/ 2
       )
      );

  std::uniform_int_distribution<uintidx>
    uniformdis_uiMinMaxK
    ((uintidx) aiinp_inParamGGA.getNumClusterKMinimum(), 
     (uintidx) aiinp_inParamGGA.getNumClusterKMaximum()
     );
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0,1);

  std::uniform_int_distribution<uintidx> uniformdis_ui0SubPopSize
    (0, aiinp_inParamGGA.getSubPopulationSize()-1);

  std::uniform_int_distribution<uintidx> uniformdis_ui0NumIsland
    (0, aiinp_inParamGGA.getNumIsland()-1 );

  
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012	      
  lochrom_best.setFitness(-measuare_undefDBindex(T_REAL)); 
  lochrom_best.setObjetiveFunc(measuare_undefDBindex(T_REAL));
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012

#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
  lochrom_best.setFitness(measuare_lowerValueSilhouette(T_REAL));
  lochrom_best.setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	 
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gga_vklabel";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output inout::OutParamEAClustering&: aoop_outParamEAC[" 
	      << &aoop_outParamEAC << "]\n"
              << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  *aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << "\t input  InParamClusteringGaProbFk&: aiinp_inParamGGA[" 
	      << &aiinp_inParamGGA << "]\n";
    aiinp_inParamGGA.print();
    std::cout << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 


  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamGGA.getNumMaxGenerations(), "Iterations", "Clustering metrics");

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream             lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_objetiveFunc = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectort_statfuncObjetiveFunc;
  
  if ( aiinp_inParamGGA.getWithPlotStatObjetiveFunc() ) {  
    
    lvectort_statfuncObjetiveFunc.reserve
      ( aiinp_inParamGGA.getSubPopulationSize());
    //DEFINE FUNCTION
    lofh_objetiveFunc  = new runtime::RuntimeFunctionValue<T_REAL>
      ("DB", 
       aiinp_inParamGGA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_objetiveFunc);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamGGA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamEAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamGGA.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamGGA.getTimesRunAlgorithm()
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
  
  int  li_opMutation = 1;
  
  /*VARIABLE NEED FOR POPULATION
   */
  std::vector<std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >* >
    lvectorvector_subPopulation;
  lvectorvector_subPopulation.reserve(aiinp_inParamGGA.getNumIsland());
  
  std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >  lvectorchrom_bestIsland;


#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
 
  /*calculate matrix dissimilarity
   */
  mat::MatrixTriang<T_REAL>&&
    lmatrixtriagT_dissimilarity = 
    dist::getMatrixDissimilarity
    (aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	    
  
  {/*BEGIN INITIALIZATION BEST ISLAND*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "0. INITIALIZATION BEST ISLAND";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< ',' << "number_island = "
		<< aiinp_inParamGGA.getNumIsland() 
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    uintidx  lui_numKMean =
      uintidx((aiinp_inParamGGA.getNumClusterKMaximum() +
	       aiinp_inParamGGA.getNumClusterKMinimum() )/ 2);

    for (uintidx lui_l = 0; lui_l < aiinp_inParamGGA.getNumIsland(); lui_l++) {
      lvectorchrom_bestIsland.push_back
	(new gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>(lui_numKMean));
      gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> *lchom_bestIslandIter =
	lvectorchrom_bestIsland.back();

#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012	      
      lchom_bestIslandIter->setFitness(-measuare_undefDBindex(T_REAL)); 
      lchom_bestIslandIter->setObjetiveFunc(measuare_undefDBindex(T_REAL));
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012

#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
      lchom_bestIslandIter->setFitness(measuare_lowerValueSilhouette(T_REAL));
      lchom_bestIslandIter->setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
       
#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	lchom_bestIslandIter->print();
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
        
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')';
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    
  }/*END INITIALIZATION BEST ISLAND*/
  
  /*1. INITIALIZATION-----------------------------------------------------------
    1. Initialization
    a) Input object number M, iteration number Gm , crossover
    probability Pc , evolutionary population N, and array of the
    numbers of clusters, i.e. K = {k 1 ,k 2 ,y,k N }.
    b) Adopt the maximum attribute range partition method to
    choose initial cluster centers for N times to form N initial
    individuals.
  */
  
  {/*BEGIN INITIALIZATION POPULATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "1. POPULATION INITIALIZATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< ',' << "number_island = "
		<< aiinp_inParamGGA.getNumIsland()
		<< ',' << "sub_population_size = "
		<< aiinp_inParamGGA.getSubPopulationSize() 
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
 
    /*CREATE SPACE FOR STORE POPULATION
     */    
    for (uintidx lui_l = 0;
	 lui_l < aiinp_inParamGGA.getNumIsland();
	 lui_l++)
      {
	std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	  lvectorchrom_subpopulation =
	  new std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >();
	
	lvectorchrom_subpopulation->reserve
	  ( aiinp_inParamGGA.getSubPopulationSize() );

	lvectorvector_subPopulation.push_back(lvectorchrom_subpopulation);
	
#ifdef __VERBOSE_YES
	const char* lpc_labelFunc  = "CREATE SUBPOPULATION";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << lpc_labelFunc  
	    << ": IN(" << geiinparam_verbose << ')'
	    << ','
	    << "lvectorchrom_subpopulation["
	    << &lvectorchrom_subpopulation << ']'
	    << std::endl;
	}
#endif /*__VERBOSE_YES*/
   
	for (uintidx lui_i = 0; 
	     lui_i < aiinp_inParamGGA.getSubPopulationSize(); 
	     lui_i++) 
	  {
	    /*Generate a number Ki in range Kmin to Kmax
	     */
	    uintidx luintidx_krand = uniformdis_uiMinMaxK(gmt19937_eng);
	    gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* liter_iChrom =
	      new gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>( luintidx_krand );
	    /*CHROMOSOME INITIALIZATION
	     */
	    T_CLUSTERIDX lmmidx_chromK = T_CLUSTERIDX( luintidx_krand );
	    std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_mmcidx0K
	      (0,lmmidx_chromK-1);
	
	    gagenericop::initializeGenes
	      (liter_iChrom->begin(),
	       liter_iChrom->end(),
	       [&]() 
	       {
		 return uniformdis_mmcidx0K(gmt19937_eng);
	       }
	       );
	    liter_iChrom->initializeGroupSec();
	    
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012	      
	    liter_iChrom->setFitness(-measuare_undefDBindex(T_REAL)); 
	    liter_iChrom->setObjetiveFunc(measuare_undefDBindex(T_REAL));
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012

#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	    liter_iChrom->setFitness(measuare_lowerValueSilhouette(T_REAL));
	    liter_iChrom->setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012

	    lvectorchrom_subpopulation->push_back(liter_iChrom);
	    
#ifdef __VERBOSE_YES
	    
	    assert(liter_iChrom->isValid() ||
		   assert_msg( geiinparam_verbose << ": Chromosoma is invalid"));
		   
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      liter_iChrom->print();
	      std::cout << std::endl;
	    }
	    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	    
	  } /*for SubPopulation*/

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << lpc_labelFunc
		    << ": OUT(" << geiinparam_verbose << ')'
		    << ',' << "lvectorchrom_subpopulation[" << &lvectorchrom_subpopulation << ']'
		    << ", size = " << lvectorchrom_subpopulation->size()
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

      } /* for NumIsland */
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< ',' << "subPopulation_size = " << lvectorvector_subPopulation.size();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*END INITIALIZATION*/


  while ( 1 ) {

    /*2. Evaluation of individuals
      a) Obtain new cluster centers by k-means;
      b) Cluster the objects according to new cluster centers and
      calculate the fitness of the individuals on the basis of Eq. (3);
    */

    {/*BEGIN EVALUATION OF INDIVIDUALS OBJETIVE FUNCTION*/
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "2.2 CLUSTERING EVALUATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << ',' << "population_size = " << lvectorvector_subPopulation.size()
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for (auto lvector_iSubPopulation: lvectorvector_subPopulation) {
	for (auto liter_iChrom: *lvector_iSubPopulation) {
	  //DECODE CHROMOSOME
	  partition::PartitionLabel
	    <T_CLUSTERIDX>
	    lpartition_clusters
	    (liter_iChrom->getString(),
	     gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize(),
	     liter_iChrom->getNumClusterK()
	     );

#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012

	  uintidx lui_numClusterK = uintidx(liter_iChrom->getNumClusterK());
	    
	  mat::MatrixRow<T_FEATURE> 
	    lmatrixrowt_centroids
	    (lui_numClusterK,
	     data::Instance<T_FEATURE>::getNumDimensions() 
	     );

	  mat::MatrixRow<T_FEATURE_SUM>       
	    lmatrixrowt_sumInstCluster
	    (lui_numClusterK,
	     data::Instance<T_FEATURE>::getNumDimensions()
	     );
	
	  std::vector<T_INSTANCES_CLUSTER_K> 
	    lvectort_numInstClusterK
	    (lui_numClusterK);

	  T_CLUSTERIDX  lmcidx_numClusterNull =
	    clusteringop::getCentroids
	    (lmatrixrowt_centroids,
	     lmatrixrowt_sumInstCluster,
	     lvectort_numInstClusterK,
	     lpartition_clusters,
	     aiiterator_instfirst,
	     aiiterator_instlast
	     );
	  
	  if ( lmcidx_numClusterNull == 0) {

	    std::vector<T_REAL>
	      lvectorrt_avgRadiusClusterK = 
	      um::avgRadiusClusterK
	      (lmatrixrowt_centroids,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       lpartition_clusters,
	       aifunc2p_dist
	       );
	      
	    liter_iChrom->setObjetiveFunc
	      (um::dbindex
	       (lmatrixrowt_centroids,
		lvectorrt_avgRadiusClusterK,
		aifunc2p_dist
		)
	       );
	    liter_iChrom->setValidString(true);
	  }
	  else {
	    liter_iChrom->setObjetiveFunc(measuare_undefDBindex(T_REAL));
	    liter_iChrom->setFitness(-measuare_undefDBindex(T_REAL)); 
	    liter_iChrom->setValidString(false);
	    aoop_outParamEAC.incTotalInvalidOffspring();
	  }
	  
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012

#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012

	  ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>
	    &&lpartlinknuminst_memberShip =
	    ds::getPartitionLinkedNumInst
	    (aiiterator_instfirst,
	     aiiterator_instlast,
	     lpartition_clusters,
	     [&](data::Instance<T_FEATURE>* liter_inst)
	     {
	       return T_INSTANCES_CLUSTER_K(1);
	     }
	     );

	  auto li_clusterNull =
	    std::count_if
	    (lpartlinknuminst_memberShip.getVectorNumInstClusterK().begin(),
	     lpartlinknuminst_memberShip.getVectorNumInstClusterK().end(),
	     [](T_INSTANCES_CLUSTER_K liter_numClusterK)
	     {return liter_numClusterK == 0;}
	     );

	  if ( li_clusterNull != 0 ) {
	    liter_iChrom->setObjetiveFunc(measuare_lowerValueSilhouette(T_REAL));
	    liter_iChrom->setFitness(measuare_lowerValueSilhouette(T_REAL));
	    liter_iChrom->setValidString(false);
	    aoop_outParamEAC.incTotalInvalidOffspring();
	  }
	  else {

	    liter_iChrom->setObjetiveFunc
	      (um::silhouette
	       (lmatrixtriagT_dissimilarity,
		lpartlinknuminst_memberShip
		)
	       );
	    
	    liter_iChrom->setFitness(liter_iChrom->getObjetiveFunc());
	    liter_iChrom->setValidString(true);
	  }   
	  
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012

	}
	
      } /*END All individual of the island determinate 
	  the quality clustering (lvectorvector_subPopulation)  
	*/
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')';
	std::cout << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END EVALUATION OF INDIVIDUALS
       */

    
    { /*BEGIN PRESERVING THE BEST STRING*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
	std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	  lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l); //.get();
	       
	const auto lit_chromMaxIsland =
	  *(std::max_element
	    (lvectorchrom_subpopulation->begin(), 
	     lvectorchrom_subpopulation->end(), 
	     [](const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* x, 
		const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* y
		) 
	  {
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	    return x->getObjetiveFunc() > y->getObjetiveFunc();
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	     
#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	    return x->getObjetiveFunc() < y->getObjetiveFunc();
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	  }
	     )
	    );
	
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	if ( lvectorchrom_bestIsland[lui_l]->getObjetiveFunc() >
	     lit_chromMaxIsland->getObjetiveFunc() )
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	     
#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	  if ( lvectorchrom_bestIsland[lui_l]->getObjetiveFunc() <
	       lit_chromMaxIsland->getObjetiveFunc() )
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	
	    {
	      /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION*/
	      *lvectorchrom_bestIsland[lui_l] = *lit_chromMaxIsland;
	    }
      }
      
      const auto lit_chromMax =
	*(std::max_element
	  (lvectorchrom_bestIsland.begin(), 
	   lvectorchrom_bestIsland.end(), 
	   [](const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* x, 
	      const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* y
	      ) 
	{
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	  return x->getObjetiveFunc() > y->getObjetiveFunc();
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	     
#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	  return x->getObjetiveFunc() < y->getObjetiveFunc();
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012	   
	}
	   )
	  );
    
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
      if ( lochrom_best.getObjetiveFunc() > lit_chromMax->getObjetiveFunc() )
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	     
#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	if ( lochrom_best.getObjetiveFunc() < lit_chromMax->getObjetiveFunc() )
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
      
	  {
	    /*CHROMOSOME ONE WAS FOUND IN THIS ITERATION
	     */
	    lochrom_best = *lit_chromMax;
	 
	    aoop_outParamEAC.setIterationGetsBest
	      (llfh_listFuntionHist.getDomainUpperBound());
	    aoop_outParamEAC.setRunTimeGetsBest
	      (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {

	      partition::PartitionLabel
		<T_CLUSTERIDX>
		lpartition_clusters
		(lochrom_best.getString(),
		 gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize(),
		 lochrom_best.getNumClusterK()
		 );

	      uintidx lui_numClusterK = uintidx(lochrom_best.getNumClusterK());

	      mat::MatrixRow<T_FEATURE> 
		lmatrixrow_centroidsChromBest
		(lui_numClusterK,
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );

	      mat::MatrixRow<T_FEATURE_SUM>       
		lmatrixrowt_sumInstClusterChromBest
		(lui_numClusterK,
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );
	
	      std::vector<T_INSTANCES_CLUSTER_K> 
		lvectort_numInstClusterKChromBest
		(lui_numClusterK
		 );
	      
	      clusteringop::getCentroids
		(lmatrixrow_centroidsChromBest,
		 lmatrixrowt_sumInstClusterChromBest,
		 lvectort_numInstClusterKChromBest,
		 lpartition_clusters,
		 aiiterator_instfirst,
		 aiiterator_instlast
		 );
	      
	      std::ostringstream lostrstream_labelCentroids;
	      lostrstream_labelCentroids
		<< "<CENTROIDSCLUSTER: " << geverbosepc_labelstep 
		<< ": generation " <<  llfh_listFuntionHist.getDomainUpperBound()
		<< ": lmatrixrow_centroidsChromBest["
		<< &lmatrixrow_centroidsChromBest << ']';
	      lmatrixrow_centroidsChromBest.print
		(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
	      std::cout << std::endl;
       
	    }
	    --geiinparam_verbose;
#endif // __VERBOSE_YES
	
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
    if ( aiinp_inParamGGA.getWithPlotStatObjetiveFunc() ) {
      for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
	std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	  lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);//.get();
	lofh_objetiveFunc->setValue(lvectorchrom_bestIsland.at(lui_l)->getObjetiveFunc());
	std::for_each
	  (lvectorchrom_subpopulation->begin(),
	   lvectorchrom_subpopulation->end(),
	   [&](gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* liter_iChrom)
	   {
	     lvectort_statfuncObjetiveFunc.push_back(liter_iChrom->getObjetiveFunc());
	   }
	   );
      }
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectort_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectort_statfuncObjetiveFunc.clear();
    }
#endif //__WITHOUT_PLOT_STAT


    /*TERMINATION CRITERION 
     */
    { /*TERMINATION CRITERION*/
      
#ifdef __VERBOSE_YES
      
      /*ID PROC
       */
      ++geverboseui_idproc;
      
      const char *lpc_labelStep = "TERMINATION CRITERION ATTAINED:";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << lpc_labelStep
	  << "generation "
	  << llfh_listFuntionHist.getDomainUpperBound()
	  << std::endl;
      }
      --geiinparam_verbose;
    
#endif /*__VERBOSE_YES*/

    
      if ( (llfh_listFuntionHist.getDomainUpperBound() >=
	    aiinp_inParamGGA.getNumMaxGenerations() ) ||
	   (runtime::elapsedTime(let_executionTime) >
	    aiinp_inParamGGA.getMaxExecutiontime())
	   )
	break;
      llfh_listFuntionHist.increaseDomainUpperBound();
    } /*END 3.1.5 TERMINATION CRITERION*/
 

    /*3.6. Local search
      We use a local search procedure to try to find local optimums in a
      close neighborhood of an individual. The local search proposed is
      based on slight modifications of the current individual, as far as they
      produce an increase of the associated objective function.
    */
    { /*BEGIN LOCAL SEARCH*/
      
      T_REAL lrt_pbj =
	aiinp_inParamGGA.getPbi() +
	((T_REAL) llfh_listFuntionHist.getDomainUpperBound() /
	 (T_REAL) aiinp_inParamGGA.getNumMaxGenerations())
	* ( aiinp_inParamGGA.getPbf() - aiinp_inParamGGA.getPbi() );
     	
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "LOCAL SEARCH";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ":  IN(" << geiinparam_verbose << ')'
		  << ',' << "lrt_pbj = " << lrt_pbj
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
	  
	std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	  lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);
	  
	for (uintidx lui_i = 0; lui_i < lvectorchrom_subpopulation->size(); lui_i++) {
	  
	  if ( uniformdis_real01(gmt19937_eng) < lrt_pbj ) {

	    //DECODE CHROMOSOME
	    gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* liter_iChrom =
	      lvectorchrom_subpopulation->at(lui_i);
	   
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	      	      
	    partition::PartitionLabel
	      <T_CLUSTERIDX>
	      lpartition_clusters
	      (liter_iChrom->getString(),
	       gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize(),
	       liter_iChrom->getNumClusterK()
	       );

	    uintidx lui_numClusterK = uintidx(liter_iChrom->getNumClusterK());

	    mat::MatrixRow<T_FEATURE> 
	      lmatrixrowt_centroids
	      (lui_numClusterK,
	       data::Instance<T_FEATURE>::getNumDimensions() 
	       );

	    mat::MatrixRow<T_FEATURE_SUM>       
	      lmatrixrowt_sumInstCluster
	      (lui_numClusterK,
	       data::Instance<T_FEATURE>::getNumDimensions()
	       );
	
	    std::vector<T_INSTANCES_CLUSTER_K> 
	      lvectort_numInstClusterK(lui_numClusterK);

	    clusteringop::getCentroids
	      (lmatrixrowt_centroids,
	       lmatrixrowt_sumInstCluster,
	       lvectort_numInstClusterK,
	       lpartition_clusters,
	       aiiterator_instfirst,
	       aiiterator_instlast
	       );
	 	
	    clusteringop::reassignCluster 
	      (liter_iChrom->getString(),
	       lmatrixrowt_centroids,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );

	    lmatrixrowt_sumInstCluster.initialize();
    
	    interfacesse::copya
	      (lvectort_numInstClusterK.data(),
	       T_INSTANCES_CLUSTER_K(0),
	       (uintidx) lvectort_numInstClusterK.size()
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
	   
	    if ( lmcidx_numClusterNull != 0 ) {
	      std::vector<uintidx> lvectormcidx_clustersKeep;
	      lvectormcidx_clustersKeep.reserve
		(lvectort_numInstClusterK.size() - (uintidx) lmcidx_numClusterNull);
	      for ( uintidx lui_i = 0;
		    lui_i < lvectort_numInstClusterK.size();
		    lui_i++) {
		if ( lvectort_numInstClusterK[lui_i] != 0 ) {
		  lvectormcidx_clustersKeep.push_back(lui_i);
		}
	      }
	      liter_iChrom->decrementGroupSecSize((uintidx) lmcidx_numClusterNull);
	      gaintegerop::labelKeep
		<T_CLUSTERIDX,T_REAL>
		(*liter_iChrom,
		 lvectormcidx_clustersKeep
		 );
		  
	      lpartition_clusters.setNumCluster
		( liter_iChrom->getNumClusterK() );
		  
	      mat::MatrixRow<T_FEATURE>  lmatrixrowt_centroidswithoutNull
		((uintidx)liter_iChrom->getNumClusterK(),
		 data::Instance<T_FEATURE>::getNumDimensions()
		 );

	      uintidx lui_iRow = 0;
	      for ( uintidx lui_i = 0;
		    lui_i < lvectort_numInstClusterK.size();
		    lui_i++) {
		if ( lvectort_numInstClusterK[lui_i] != 0 ) {
		  lmatrixrowt_centroidswithoutNull.copyRow
		    (lui_iRow,
		     lmatrixrowt_centroids.getRow(lui_i)
		     );
		  ++lui_iRow;
		}
	      }
	      lmatrixrowt_centroids = std::move(lmatrixrowt_centroidswithoutNull);
	    }
	    std::vector<T_REAL>
	      lvectorrt_avgRadiusClusterK = 
	      um::avgRadiusClusterK
	      (lmatrixrowt_centroids,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       lpartition_clusters,
	       aifunc2p_dist
	       );
	    
	    T_REAL lrt_dbindex =
	      um::dbindex
	      (lmatrixrowt_centroids,
	       lvectorrt_avgRadiusClusterK,
	       aifunc2p_dist
	       );
	     
#ifdef __VERBOSE_YES
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      std::cout
		<< geverbosepc_labelstep
		<< "liter_iChrom->objetiveFunc = "
		<< liter_iChrom->getObjetiveFunc() 
		<< ", lrt_dbindex = " << lrt_dbindex
		<< std::endl;
	    }
	    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	    liter_iChrom->setObjetiveFunc(lrt_dbindex);
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	      
#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
 
	    partition::PartitionLabel
	      <T_CLUSTERIDX>
	      lpartition_clusters
	      (liter_iChrom->getString(),
	       gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize(),
	       liter_iChrom->getNumClusterK()
	       );

	    ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>
	      &&lpartlinknuminst_memberShip =
	      ds::getPartitionLinkedNumInst
	      (aiiterator_instfirst,
	       aiiterator_instlast,
	       lpartition_clusters,
	       [&](data::Instance<T_FEATURE>* liter_inst)
	       {
		 return T_INSTANCES_CLUSTER_K(1);
	       }
	       );

#ifdef __VERBOSE_YES
	    T_REAL lrt_origSilhouette;
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      lrt_origSilhouette =
		(um::silhouette
		 (lmatrixtriagT_dissimilarity,
		  lpartlinknuminst_memberShip
		  )
		 );	      
	    }
	    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	      
	    T_CLUSTERIDX *larraymmidx_iChrom = liter_iChrom->getString();

	    for (uintidx lui_i = 0;
		 lui_i <  gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize();
		 lui_i++)
	      {
		T_CLUSTERIDX lmmidx_geneMaxObj  = -1;
		T_REAL lrt_geneSilhouetteMax =
		  measuare_lowerValueSilhouette(T_REAL);
		lpartlinknuminst_memberShip.subInstanceFromCluster
		  (*larraymmidx_iChrom,lui_i);
		
		for ( T_CLUSTERIDX lmmidx_iterK = 0;
		      lmmidx_iterK < liter_iChrom->getNumClusterK();
		      ++lmmidx_iterK ) 
		  {

		    T_REAL lrt_geneSilhouette  = 
		      um::silhouette
		      (lui_i,
		       lmmidx_iterK,
		       lmatrixtriagT_dissimilarity,
		       lpartlinknuminst_memberShip,
		       lpartlinknuminst_memberShip.getVectorNumInstClusterK()
		       );

		    if ( lrt_geneSilhouette > lrt_geneSilhouetteMax ) {
		      lmmidx_geneMaxObj = lmmidx_iterK;
		      lrt_geneSilhouetteMax = lrt_geneSilhouette;
		    }
		  } /*FOR K*/
		
		*larraymmidx_iChrom = lmmidx_geneMaxObj;
		lpartlinknuminst_memberShip.addInstanceToCluster
		  (lmmidx_geneMaxObj,lui_i);
		++larraymmidx_iChrom;
	      }

	    T_REAL lrt_searchLocalSilhouette =
	      (um::silhouette
	       (lmatrixtriagT_dissimilarity,
		lpartlinknuminst_memberShip
		)
	       );
	    for ( T_CLUSTERIDX lmmc_k = (liter_iChrom->getNumClusterK()-1) ;
		  lmmc_k >= 0; --lmmc_k) {
	      if ( lpartlinknuminst_memberShip.getVectorNumInstClusterK().at(lmmc_k)
		   ==  0 )		
		liter_iChrom->deleteGroupNull(lmmc_k);
	    }
	    
#ifdef __VERBOSE_YES
	    ++geiinparam_verbose;
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
  
	      partition::PartitionLabel
		<T_CLUSTERIDX>
		lpartitionLabel_reclassifierC
		(liter_iChrom->getString(),
		 gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize(),
		 liter_iChrom->getNumClusterK()
		 );

	      ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>
		&&lpartlinkNumInst_rememberShip =
		ds::getPartitionLinkedNumInst
		(aiiterator_instfirst,
		 aiiterator_instlast,
		 lpartitionLabel_reclassifierC,
		 [&](data::Instance<T_FEATURE>* liter_inst)
		 {
		   return T_INSTANCES_CLUSTER_K(1);
		 }
		 );
	    
	      T_REAL lrt_clacSilhouette =
		(um::silhouette
		 (lmatrixtriagT_dissimilarity,
		  lpartlinkNumInst_rememberShip
		  )
		 );
		 
	      std::cout
		<< geverbosepc_labelstep << ":"
		<< geverboseui_idproc
		<< ": liter_iChrom->objetiveFunc = "
		<< liter_iChrom->getObjetiveFunc()
		<< ", lrt_origSilhouette = " << lrt_origSilhouette
		<< ", lrt_searchLocalSilhouette = "
		<< lrt_searchLocalSilhouette
		<< ", lrt_clacSilhouette = "
		<<  lrt_clacSilhouette
		<<  '\n';
	      liter_iChrom->print();
	      std::cout << std::endl;
	    }
	    if ( liter_iChrom->isValid() == false ) {
	      std::cout << "lrt_searchLocalSilhouette = "
			<< lrt_searchLocalSilhouette
			<< std::endl;
	      liter_iChrom->print();
	      std::cout << std::endl;
	      inout::containerprint
		(lpartlinknuminst_memberShip.getVectorNumInstClusterK().begin(),
		 lpartlinknuminst_memberShip.getVectorNumInstClusterK().end(),
		 std::cout,
		 geverbosepc_labelstep,
		 ','
		 );
	      std::cout << std::endl;
	      assert( false ||
		      assert_msg("LOCAL SEARCH: Chromosoma is invalid"));
	    }
	    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

	    liter_iChrom->setObjetiveFunc(lrt_searchLocalSilhouette); 
	      
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	      
	  } /* if */
	}
      }
      	
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << ',' << "subPopulation_size = "
	  << lvectorvector_subPopulation.size()
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END LOCAL SEARCH
       */ 

    /*3.7 AN ISLAND MODEL TO IMPROVE  
     */
    { /*BEGIN ISLAND MODEL TO IMPROVE  
       */
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "3.7 ISLAND MIGRATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": IN(" << geiinparam_verbose << ')'
	  << ',' << "subpopulation_size = "
	  << lvectorvector_subPopulation.size()
	  << std::endl;
      }
#endif /*__VERBOSE_YES*/
	
      for (uintidx lui_l = 0; lui_l < lvectorchrom_bestIsland.size(); lui_l++) {
	
	if ( uniformdis_real01(gmt19937_eng) < aiinp_inParamGGA.getPe() ) {
	  /*2. Randomly choose the island toward each individual 
	    will migrate.
	  */
	  uintidx lui_islandToward = 
	    prob::getRandUnlike
	    (lui_l,
	     [&]() 
	     {
	       return uniformdis_ui0NumIsland(gmt19937_eng);
	     }
	     );
	  
	  std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	    lvectorchrom_subpopulation =
	    lvectorvector_subPopulation.at(lui_islandToward);
	  
	  /*3. Randomly choose an individual in the destiny island 
	    and change it by the migrating individual.
	  */
	  uintidx lui_idxIndividual =
	    uniformdis_ui0SubPopSize(gmt19937_eng);
	  
	  *lvectorchrom_subpopulation->at(lui_idxIndividual) =
	    *lvectorchrom_bestIsland[lui_l];
	}
      }
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << ": OUT(" << geiinparam_verbose << ')'
	  << ',' << "subPopulation_size = "
	  << lvectorvector_subPopulation.size()
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
    } /* END ISLAND MODEL TO IMPROVE  
       */
 
    /*3. GENETIC OPERATIONS
     */
    { /*BEGIN GENETIC OPERATIONS*/
      
      /*Evaluate fitness value by rank of the individual
	First, the individuals are sorted in a list based 
	on their quality. The position of the individuals 
	in the list is called rank of the individual,
	A fitness value associated to each individual is 
	then defined, as follows:
	f_i = \frac{ 2 \cdot R_i}{\xi \cdot ( \xi + 1 )}
      */
      { /*BEGIN EVALUATE FITNESS VALUE BY RANK
	 */
	
#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "EVALUATE FITNESS VALUE BY RANK";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": IN(" << geiinparam_verbose << ')'
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
	  std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	    lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);
      	    
	  T_REAL lrt_scaleProbDist = 2.0 /
	    ( (T_REAL) lvectorchrom_subpopulation->size() *
	      ((T_REAL) lvectorchrom_subpopulation->size() + 1.0));
	  
#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	    
	  prob::linearNormalization
	    (lvectorchrom_subpopulation->begin(),
	     lvectorchrom_subpopulation->end(),
	     [](const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* lchrom_iter)
	     -> T_REAL
	     {
	       return 1.0 / lchrom_iter->getObjetiveFunc();
	     },
	     [](gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* lchrom_iter,
		T_REAL airt_funcFitnessLineal)
	     {
	       lchrom_iter->setFitness(airt_funcFitnessLineal);
	     },
	     lrt_scaleProbDist
	     );
	  
#endif //ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
	     
#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	  
	  prob::linearNormalization
	    (lvectorchrom_subpopulation->begin(),
	     lvectorchrom_subpopulation->end(),
	     [](const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* lchromfeac_iter)
	     -> T_REAL
	     {
	       return lchromfeac_iter->getObjetiveFunc();
	     },
	     [](gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* lchromfeac_iter,
		T_REAL airt_funcFitnessLineal)
	     {
	       lchromfeac_iter->setFitness(airt_funcFitnessLineal);
	     },
	     lrt_scaleProbDist
	     );
	  
#endif //ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
	  
	}

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

	
      } /*END EVALUATE FITNESS VALUE BY RANK
	 */

      std::vector
	<std::unique_ptr
	 <std::vector<std::unique_ptr<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> > > > >
	lvectorvector_subPoolString;
      lvectorvector_subPoolString.reserve( lvectorvector_subPopulation.size() );
      
      
      /*3.2 Selection operator
       */
      {/*BEGIN SELECTION*/

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "3.2 SELECTION";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << geverbosepc_labelstep
	    << ": IN(" << geiinparam_verbose << ')'
	    << ',' << "subpopulation_size = "
	    << lvectorvector_subPopulation.size()
	    << ',' << "subpoolstring_size = "
	    << lvectorvector_subPoolString.size()
	    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {

	  std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	    lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);

	  const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	    prob::makeDistRouletteWheel
	    (lvectorchrom_subpopulation->begin(),lvectorchrom_subpopulation->end(),
	     [](const gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* liter_iChrom)
	     -> T_REAL
	     {
	       return liter_iChrom->getFitness();
	     }
	     );
	
	  lvectorvector_subPoolString.push_back
	    (std::unique_ptr
	     <std::vector< std::unique_ptr<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> > > >
	     (new std::vector<std::unique_ptr< gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> > >())
	     ); 
	  std::vector<std::unique_ptr<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> > >
	    *lvectorchrom_subPoolString =
	    lvectorvector_subPoolString.back().get();	  
	  lvectorchrom_subPoolString->reserve
	    ( lvectorchrom_subpopulation->size() );
	
	  /*ELITISMO
	   */
	  lvectorchrom_subPoolString->push_back
	    ( std::unique_ptr<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> >
	      (new gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>
	       ( *lvectorchrom_bestIsland[lui_l] )
	       )
	      );
	  for (uintidx lui_i = 1; 
	       lui_i < lvectorchrom_subpopulation->size(); 
	       lui_i++)
	    {
	      uintidx luiidx_chrom =
		gaselect::getIdxRouletteWheel
		(lvectorT_probDistRouletteWheel,
		 uintidx(0)
		 );
	     
	      lvectorchrom_subPoolString->push_back
		( std::unique_ptr<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> >
		  (new gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>
		   ( *lvectorchrom_subpopulation->at(luiidx_chrom) )
		   )
		  );
	    }
	  
	}
	
#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << geverbosepc_labelstep
	    << ": OUT(" << geiinparam_verbose << ')'
	    << ','
	    << "subPopulation_size = "
	    << lvectorvector_subPopulation.size()
	    << ','
	    << "subPoolString_size = "
	    << lvectorvector_subPoolString.size()
	    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      } /*END SELECTION*/

      /*3.3 Crossover operator
       */    
      { /*BEGIN CROSSOVER*/

	T_REAL lrt_pcj =
	  aiinp_inParamGGA.getPci() +
	  ((T_REAL) llfh_listFuntionHist.getDomainUpperBound() /
	   (T_REAL) aiinp_inParamGGA.getNumMaxGenerations())
	  * ( aiinp_inParamGGA.getPcf() - aiinp_inParamGGA.getPci() );
	

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "3.3 CROSSOVER OPERATOR";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << geverbosepc_labelstep
	    << ": IN(" << geiinparam_verbose << ')'
	    << ','
	    << "subPopulation_size = "
	    << lvectorvector_subPopulation.size()
	    << ','
	    << "subPoolString_size = "
	    << lvectorvector_subPoolString.size()
	    << ',' << "lrt_pcj = " << lrt_pcj
	    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
	  
	  std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	    lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);
	  
	  std::vector<std::unique_ptr<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> > >
	    *lvectorchrom_subPoolString = lvectorvector_subPoolString.at(lui_l).get();
	  
	  uintidx luintidx_idxCrossChromi = 0;
	  while ( luintidx_idxCrossChromi < lvectorchrom_subpopulation->size() ) {
	    //GET TWO RANDOM
	    std::pair<uintidx,uintidx>
	      lpair_idxChrom =
	      prob::getRandPairUnlike
	      ([&]() -> uintidx
	       {
		 return uniformdis_ui0SubPopSize(gmt19937_eng);
	       }
	       );

	    if ( uniformdis_real01(gmt19937_eng) < lrt_pcj ) {

	      *lvectorchrom_subpopulation->at(luintidx_idxCrossChromi) =
		gaclusteringop::mergeCrossover
		<T_CLUSTERIDX,T_REAL>
		(*lvectorchrom_subPoolString->at(lpair_idxChrom.first).get(),
		 *lvectorchrom_subPoolString->at(lpair_idxChrom.second).get()
		 );
	      ++luintidx_idxCrossChromi;
	    }
	    else {
	      *lvectorchrom_subpopulation->at(luintidx_idxCrossChromi) =
		*lvectorchrom_subPoolString->at(lpair_idxChrom.first).get();
		
	      if ( ++luintidx_idxCrossChromi >= lvectorchrom_subpopulation->size() ) 
		break;

	      *lvectorchrom_subpopulation->at(luintidx_idxCrossChromi) =
		*lvectorchrom_subPoolString->at(lpair_idxChrom.second).get();
	      ++luintidx_idxCrossChromi;
	
	    } //END ELSE  Crossover
	 
	  } /*while Crossover*/
	}

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout
	    << geverbosepc_labelstep
	    << ": OUT(" << geiinparam_verbose << ')'
	    << ','
	    << "subPopulation_size = "
	    << lvectorvector_subPopulation.size()
	    << ','
	    << "subPoolString_size = "
	    << lvectorvector_subPoolString.size()
	    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
      } /*END CROSSOVER*/

      
      {/*BEGIN MUTATION*/

	T_REAL lrt_pmj =
	  aiinp_inParamGGA.getPmi() +
	  ((T_REAL) llfh_listFuntionHist.getDomainUpperBound() /
	   (T_REAL) aiinp_inParamGGA.getNumMaxGenerations())
	  * ( aiinp_inParamGGA.getPmf() - aiinp_inParamGGA.getPmi() );

		
#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "3.4 MUTATION OPERATOR";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ":  IN(" << geiinparam_verbose << ')'
		    << ',' << "lrt_pmj = " << lrt_pmj
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/

	for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
	  
	  std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	    lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);
	  
	  for (uintidx lui_i = 0; lui_i < lvectorchrom_subpopulation->size(); lui_i++) {

	    if ( uniformdis_real01(gmt19937_eng)  < lrt_pmj ) {
	      //DECODE CHROMOSOME
	      if ( li_opMutation == 1 ) {

		gaclusteringop::splittingMutation
		  (*lvectorchrom_subpopulation->at(lui_i));
		
	      }
	      else {
		
		gaclusteringop::mergeMutation
		  (*lvectorchrom_subpopulation->at(lui_i));
		
	      }
	      li_opMutation *= -1;
	    }


	  } /* FOR END All individual of the island
	     */
	
	} /* FOR END All individual of the island
	   */
      
#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << ": OUT(" << geiinparam_verbose << ')';
	  std::cout << std::endl;
	}
	--geiinparam_verbose;
	geverbosepc_labelstep = "END GENETIC OPERATIONS";	
#endif /*__VERBOSE_YES*/

      } /*BEGIN MUTATION*/
       
    } /*END GENETIC OPERATIONS*/
    
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
#endif //__VERBOSE_YES

    for (uintidx lui_l = 0; lui_l < lvectorvector_subPopulation.size(); lui_l++) {
      std::vector<gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* >*
	lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_l);//.get();
      for (uintidx lui_i = 0; lui_i < lvectorchrom_subpopulation->size(); lui_i++) {
	gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>* liter_iChrom =
	  lvectorchrom_subpopulation->at(lui_i);
	
	delete liter_iChrom;
      }
      delete lvectorchrom_subpopulation;
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
    
  
  {/*BEGIN FREE MEMORY OF BESTISLAND*/ 
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "DELETEBESTISLAND";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  geverbosepc_labelstep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    for (uintidx lui_i = 0; lui_i < lvectorchrom_bestIsland.size(); ++lui_i) {
      delete lvectorchrom_bestIsland[lui_i];
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
  aoop_outParamEAC.setNumClusterK
    (lochrom_best.getNumClusterK());
  aoop_outParamEAC.setMetricFuncRun
    (lochrom_best.getObjetiveFunc());
  aoop_outParamEAC.setFitness
    (lochrom_best.getFitness());
  aoop_outParamEAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamEAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  
  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamGGA.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamGGA,
       aoop_outParamEAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
  
  geverbosepc_labelstep = lpc_labelAlgGA;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA 
	      << " OUT(" << geiinparam_verbose << ")\n";
    
    lochrom_best.print();
    std::cout << std::endl;

    partition::PartitionLabel
      <T_CLUSTERIDX>
      lpartition_clusters
      (lochrom_best.getString(),
       gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>::getElementSize(),
       lochrom_best.getNumClusterK()
       );
    
    uintidx lui_numClusterK = uintidx(lochrom_best.getNumClusterK());

    mat::MatrixRow<T_FEATURE> 
      lmatrixrow_centroidsChromBest
      (lui_numClusterK,
       data::Instance<T_FEATURE>::getNumDimensions()
       );

    mat::MatrixRow<T_FEATURE_SUM>       
      lmatrixrowt_sumInstClusterChromBest
      (lui_numClusterK,
       data::Instance<T_FEATURE>::getNumDimensions()
       );
	
    std::vector<T_INSTANCES_CLUSTER_K> 
      lvectort_numInstClusterKChromBest
      (lui_numClusterK
       );

    clusteringop::getCentroids
      (lmatrixrow_centroidsChromBest,
       lmatrixrowt_sumInstClusterChromBest,
       lvectort_numInstClusterKChromBest,
       lpartition_clusters,
       aiiterator_instfirst,
       aiiterator_instlast
       );
			      
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids
      << "<CENTROIDSCLUSTER: " << lpc_labelAlgGA
      << ": generation " <<  llfh_listFuntionHist.getDomainUpperBound()
      << ": lmatrixrow_centroidsChromBest["
      << &lmatrixrow_centroidsChromBest << ']';
    lmatrixrow_centroidsChromBest.print
      (std::cout,
       lostrstream_labelCentroids.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
   
  return lochrom_best; 
 
} /* END gga_vklabel
   */

} /*END eac */

#endif /*__GGA_VKLABEL_HPP__*/
