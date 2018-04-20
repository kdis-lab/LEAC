/*! \file gcuk_vkcentroid.hpp
 *
 * \brief GCUK \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the GCUK algorithm based on the paper:\n
 * S. Bandyopadhyay and U. Maulik. Genetic clustering for automatic\n 
 * evolution of clusters and application to image classification.\n
 * Pattern Recognition, 35(6):1197 – 1208, 2002.\n
 * <a href="http://www.sciencedirect.com/science/article/pii/S003132030100108X">http://www.sciencedirect.com/science/article/pii/S003132030100108X</a>,\n
 * <a href="http://dx.doi.org/10.1016/S0031-3203(01)00108-X">doi:http://dx.doi.org/10.1016/S0031-3203(01)00108-X</a>.\n
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

#ifndef __GCUK_VKCENTROID_HPP__
#define __GCUK_VKCENTROID_HPP__

#include <vector>
#include <algorithm>

#include <leac.hpp>
#include "inparam_pcpmvk.hpp"
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
  
/*! \fn mat::MatrixRow<T_FEATURE> gcuk_vkcentroid(inout::OutParamGAC<T_REAL,T_CLUSTERIDX> &aoop_outParamGAC, inout::InParamPcPmVk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K> &aiinp_inParamPcPmVk, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist) 
  \brief GCUK \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002
  \details Implementation of the GCUK algorithm based on \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002. Which automatically finds K cluster using the Davies–Bouldin index.
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid. Base to following equation:
  \f[
  x_i \in C_j \leftrightarrow  \| x_i - \mu_j \|  \begin{array}{c}min\\  k \end{array}
  \| x_i - \mu_k \|,\; j=1,2,..k,
  \f]
  where \f$mu_j\f$, represents the centroid of cluster \f$C_j\f$
  \param aoop_outParamGAC a inout::OutParamGAC with the output parameters of the algorithm
  \param aiinp_inParamPcPmVk a inout::InParamPcPmVk parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,      //T_STRING,
	   typename T_REAL,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, 
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
mat::MatrixRow<T_FEATURE>
gcuk_vkcentroid
(inout::OutParamGAC
 <T_REAL,
 T_CLUSTERIDX>                                &aoop_outParamGAC,
 inout::InParamPcPmVk
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>                       &aiinp_inParamPcPmVk,
 const INPUT_ITERATOR                         aiiterator_instfirst,
 const INPUT_ITERATOR                         aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>                 &aifunc2p_dist
 )
{

  gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>  lochrom_best;

  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC
   */
  std::vector<gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL> >
    lvectorchrom_population;
  
  std::vector<gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL> >
    lvectorchrom_matingPool;
   
  std::uniform_int_distribution<uintidx>
    uniformdis_uiMinMaxK
    ((uintidx) aiinp_inParamPcPmVk.getNumClusterKMinimum(), 
     (uintidx) aiinp_inParamPcPmVk.getNumClusterKMaximum()
     );
    
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "gcuk_vkcentroid:";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << "  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output inout::OutParamGAC&: aoop_outParamGAC[" 
	      << &aoop_outParamGAC << "]\n"
	      << "\t input  InParamClusteringGaProbFk&: aiinp_inParamPcPmVk[" 
	      << &aiinp_inParamPcPmVk << "]\n"
	      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t\tPopulation size = " 
	      << aiinp_inParamPcPmVk.getSizePopulation()
	      << "\n\t\tProbCrossover = " 
	      << aiinp_inParamPcPmVk.getProbCrossover() 
	      << "\n\t\tProbMutation  = " 
	      << aiinp_inParamPcPmVk.getProbMutation()
	      << "\n\t\tk-minimum  = " 
	      << aiinp_inParamPcPmVk.getNumClusterKMinimum()
	      << "\n\t\tk-maximum  = " 
	      << aiinp_inParamPcPmVk.getNumClusterKMaximum()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 


  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinp_inParamPcPmVk.getNumMaxGenerations(),
     "Iterations",
     "Clustering metrics"
     );

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  std::ofstream  lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_DBindex = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL> lvectort_statfuncObjetiveFunc;
  
  if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
    
    lvectort_statfuncObjetiveFunc.reserve
      ( aiinp_inParamPcPmVk.getSizePopulation());
    //DEFINE FUNCTION
    lofh_DBindex  = new runtime::RuntimeFunctionValue<T_REAL>
      ("DB", 
       aiinp_inParamPcPmVk.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_DBindex);

    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = 
	new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinp_inParamPcPmVk.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoop_outParamGAC.setFileNameOutPlotStatObjetiveFunc
      (aiinp_inParamPcPmVk.getFileNamePlotStatObjetiveFunc(),
       aiinp_inParamPcPmVk.getTimesRunAlgorithm()
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
  
  {/*BEGIN CREATE SPACE FOR STORE POPULATION*/

#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "CREATE SPACE FOR STORE POPULATION";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    /*CREATE SPACE FOR STORE POPULATION-----------------------------------------
     */

    lvectorchrom_population.reserve
      (aiinp_inParamPcPmVk.getSizePopulation());
  
    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinp_inParamPcPmVk.getSizePopulation(); 
	 luintidx_i++) 
      {
	/*Generate a number Ki in range Kmin to Kmax
	 */
	lvectorchrom_population.push_back
	  (gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>
	   (uniformdis_uiMinMaxK(gmt19937_eng),
	    (uintidx) aiinp_inParamPcPmVk.getNumClusterKMaximum(),
	    data::Instance<T_FEATURE>::getNumDimensions())
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
    
  } /*END CREATE SPACE FOR STORE POPULATION*/


  /*2.1.2. Population initialization
    For each string i in the population (i = 1;...;P, where
    P is the size of the population), a random number Ki
    in the range [Kmin;Kmax] is generated.
  */    
  { /*BEGIN POPULATION INITIALIZATION*/
      
#ifdef __VERBOSE_YES
    const char *geverbosepc_labelstep = "POPULATION INITIALIZATION:";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< " IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
   
    
    for ( auto &&liter_iChrom: lvectorchrom_population ) {
      
      /*Chose Ki point randomly from the data
       */
      clusteringop::randomInitialize
	(liter_iChrom,
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );
      
      std::shuffle
	(liter_iChrom.toArray(),
	 liter_iChrom.toArray()+liter_iChrom.getNumRowsMax(),
	 gmt19937_eng
	 );
      liter_iChrom.setObjetiveFunc(std::numeric_limits<T_REAL>::max());
      liter_iChrom.setFitness(-std::numeric_limits<T_REAL>::max());

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	liter_iChrom.print();
	std::cout << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    }


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< " OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  } /*END POPULATION INITIALIZATION*/


  while ( 1 ) {

    /*2.1.3. Fitness computation
      The fitness of a chromosome is computed using the
      Davies–Bouldin index. 
    */

    {/*BEGIN FITNESS COMPUTATION*/
      
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "FITNESS COMPUTATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto&& liter_iChrom : lvectorchrom_population ) {
      
	if ( (liter_iChrom.getNumRows() > 1) && 
	     (liter_iChrom.getFitness() == -std::numeric_limits<T_REAL>::max()) ) 
	  {
	    /*DECODE CHROMOSOME*/
	    mat::MatrixRow<T_FEATURE>&& 
	      lmatrixrow_centroidsChrom = 
	      liter_iChrom.getMatrix();

	    mat::MatrixRow<T_FEATURE_SUM>       
	      llmatrixrowt_sumInstancesCluster
	      (lmatrixrow_centroidsChrom.getNumRows(),
	       data::Instance<T_FEATURE>::getNumDimensions(),
	       T_FEATURE_SUM(0)
	       );
	    
	    std::vector<T_INSTANCES_CLUSTER_K> 
	      lvectort_numInstancesInClusterK
	      (lmatrixrow_centroidsChrom.getNumRows(),
	       T_INSTANCES_CLUSTER_K(0)
	       );
	  
	    T_CLUSTERIDX lmcidx_numClusterNull; 
	    clusteringop::updateCentroids
	      (lmcidx_numClusterNull,
	       lmatrixrow_centroidsChrom,
	       llmatrixrowt_sumInstancesCluster,
	       lvectort_numInstancesInClusterK,
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );
	    liter_iChrom.setMatrix(lmatrixrow_centroidsChrom);
	     
	    if (lmcidx_numClusterNull == 0 ) {
		 
	      auto  lpartitionCentroids_clusters = 
		partition::makePartition
		(lmatrixrow_centroidsChrom,
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 T_CLUSTERIDX(lmatrixrow_centroidsChrom.getNumRows()),
		 aifunc2p_dist
		 );
	    
	      T_REAL lrt_dbindex = 
		um::dbindex
		(lmatrixrow_centroidsChrom,
		 aiiterator_instfirst,
		 aiiterator_instlast,
		 lpartitionCentroids_clusters,
		 aifunc2p_dist
		 );

	      liter_iChrom.setObjetiveFunc(lrt_dbindex); 
	      liter_iChrom.setFitness(1.0 / lrt_dbindex);

	      if ( lrt_dbindex < measuare_undefDBindex(T_REAL) ) {
		liter_iChrom.setValidString(true);
	      }
	      else {
		aoop_outParamGAC.incTotalInvalidOffspring();
		liter_iChrom.setValidString(false);
	      }
		 
	     
#ifndef __WITHOUT_PLOT_STAT
	      lvectort_statfuncObjetiveFunc.push_back(liter_iChrom.getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
	   
	    }
	    else {
	      T_REAL lrt_dbindex =  measuare_undefDBindex(T_REAL);
	      liter_iChrom.setObjetiveFunc(lrt_dbindex); 
	      liter_iChrom.setFitness(1.0 / lrt_dbindex);
	      aoop_outParamGAC.incTotalInvalidOffspring();
	      liter_iChrom.setValidString(false); 
	    }
	  }
	else {
	  if ( liter_iChrom.getNumRows() < 2 ) {
	    T_REAL lrt_dbindex = std::numeric_limits<T_REAL>::max();
	    liter_iChrom.setObjetiveFunc(lrt_dbindex); 
	    liter_iChrom.setFitness(1.0 / lrt_dbindex);
	    aoop_outParamGAC.incTotalInvalidOffspring();
	    liter_iChrom.setValidString(false); 
	  }
	}
      }
   
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << " OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    

    } /*END FITNESS COMPUTATION*/
   
    { /*BEGIN PRESERVING THE BEST STRING
       */
#ifdef __VERBOSE_YES
      const char *geverbosepc_labelstep = "ELITISM PRESERVING THE BEST";
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
	 [](const gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>& x, 
	    const gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>& y
	    ) 
	 {  return x.getFitness() < y.getFitness(); }
	 );
      
      if ( lochrom_best.getFitness() < (*lit_chromMax).getFitness() ) {
	lochrom_best = *lit_chromMax;
	/*CHROMOSOME ONE WAS FOUND IN THIS ITERATION
	 */
	aoop_outParamGAC.setIterationGetsBest
	  (llfh_listFuntionHist.getDomainUpperBound());
	aoop_outParamGAC.setRunTimeGetsBest
	  (runtime::elapsedTime(let_executionTime));

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	
	  /*DECODE CHROMOSOME*/
	  mat::MatrixRow<T_FEATURE>&& 
	    lmatrixrow_centroidsChromBest = 
	    lochrom_best.getMatrix();

	  auto lpartitionCentroids_clustersChromBest = 
	    partition::makePartition
	    (lmatrixrow_centroidsChromBest,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     T_CLUSTERIDX(lmatrixrow_centroidsChromBest.getNumRows()),
	     aifunc2p_dist
	     );
	  
	  std::ostringstream lostrstream_labelShipBest;
	  lostrstream_labelShipBest
	    << "<MEMBERCLUSTER:"
	    << geverbosepc_labelstep
	    << "generation "
	    <<  llfh_listFuntionHist.getDomainUpperBound()
	    << ":lpartitionCentroids_clustersChromBest<>["
	    << &lpartitionCentroids_clustersChromBest << ']';
	  lpartitionCentroids_clustersChromBest.print
	    (std::cout,lostrstream_labelShipBest.str().c_str(),',');
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
    if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
      lofh_DBindex->setValue(lochrom_best.getObjetiveFunc());
      functionhiststat_evaluateAll
	(lofhs_statObjectiveFunc,
	 lvectort_statfuncObjetiveFunc
	 );
      lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
      lvectort_statfuncObjetiveFunc.clear();
    }
#endif /*__WITHOUT_PLOT_STAT*/

    /*2.1.4 GENETIC OPERATIONS
     */

    { /*BEGIN GENETIC OPERATIONS*/

      /*Selection: Conventional proportional selection is 
	applied on the population of strings.
      */
      { /*BEGIN SELECTION*/

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "SELECTION";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep  
		    << ":  IN(" << geiinparam_verbose << ')'
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/
    
	const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
	  prob::makeDistRouletteWheel
	  (lvectorchrom_population.begin(),lvectorchrom_population.end(),
	   [](const gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>& lchrommatrixwrn_iter)
	   -> T_REAL
	   {
	     return lchrommatrixwrn_iter.getFitness();
	   }
	   );
	
	/*COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL--------------------------
	 */ 
	lvectorchrom_matingPool.reserve
	  (aiinp_inParamPcPmVk.getSizePopulation());

	/*ELITISMO
	 */
	lvectorchrom_matingPool.push_back
	  (gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>(lochrom_best));

	for (uintidx luintidx_i = 1; 
	     luintidx_i < aiinp_inParamPcPmVk.getSizePopulation(); 
	     luintidx_i++) 
	  {      
	    uintidx luintidx_chrom = 
	      gaselect::getIdxRouletteWheel
	      (lvectorT_probDistRouletteWheel,
	       uintidx(0)
	       );
	    
	    lvectorchrom_matingPool.push_back
	      (gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>
	       (lvectorchrom_population.at(luintidx_chrom))
	       );
	  }

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << " OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

      } /*END SELECTION*/


      /*Crossover: During crossover each cluster centre is con-
	sidered to be an indivisible gene. Single point crossover,
	applied stochastically with probability c,
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

	lvectorchrom_population.clear();

	auto&& jchrom_matingPool = lvectorchrom_matingPool.begin();
	
	if ( ( lvectorchrom_matingPool.size() % 2 ) != 0 ) {
	  ++jchrom_matingPool;
	}
	while ( jchrom_matingPool != lvectorchrom_matingPool.end() ) {
	  
	  auto lchrom_parent1  = jchrom_matingPool;
	  ++jchrom_matingPool;
	  auto lchrom_parent2  = jchrom_matingPool;
	  ++jchrom_matingPool;
	  
	  if ( uniformdis_real01(gmt19937_eng)
	       < aiinp_inParamPcPmVk.getProbCrossover() ) {

	    gaclusteringop::onePointCrossover
	      ( *lchrom_parent1,
		*lchrom_parent2
		);
	  }
	  
	} /*END While
	   */

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << " OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
      } /*END CROSSOVER*/
 
      /*Mutation: Each valid position (i.e., which is not ‘#’)
	in a chromosome is mutated with probability $\mu_m$.
      */
      { /*BEGIN MUTATION*/

#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "MUTATION";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep  
		    << ":  IN(" << geiinparam_verbose << ')'
		    << std::endl;
	}
#endif /*__VERBOSE_YES*/
    
    
	for ( auto&& liter_iChrom : lvectorchrom_matingPool ) {
	  if ( uniformdis_real01(gmt19937_eng) 
	       < aiinp_inParamPcPmVk.getProbMutation() ) {
	    //IF BEGIN  MUTATION
	    gaclusteringop::randomMutation
	      (liter_iChrom);	   
	    liter_iChrom.setFitness(-std::numeric_limits<T_REAL>::max());  
	    liter_iChrom.setObjetiveFunc(std::numeric_limits<T_REAL>::max());
	  } //END BEGIN  MUTATION
	}

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << geverbosepc_labelstep
		    << " OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    

      } /*END MUTATION */
	
    } /*END GENETIC OPERATIONS*/

    lvectorchrom_population.swap(lvectorchrom_matingPool);
     
    /*3.1.5 TERMINATION CRITERION 
     */
    { /*BEGIN 3.1.5 TERMINATION CRITERION*/
#ifdef __VERBOSE_YES

      /*ID PROC
       */
      ++geverboseui_idproc;
    
      const char *geverbosepc_labelstep = "TERMINATION CRITERION ATTAINED:";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout
	  << geverbosepc_labelstep
	  << "generation "
	  << llfh_listFuntionHist.getDomainUpperBound()
	  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
      if ( !(llfh_listFuntionHist.getDomainUpperBound() 
	     < aiinp_inParamPcPmVk.getNumMaxGenerations() ) 
	   )
	break;
      llfh_listFuntionHist.increaseDomainUpperBound();
    } /*END 3.1.5 TERMINATION CRITERION*/
    
  } /*END EVOLUTION While*/ 

  /*DECODE CHROMOSOME*/
  mat::MatrixRow<T_FEATURE>&& 
    lomatrixrow_centroidsChromBest = 
    lochrom_best.getMatrix();
  
  runtime::stop(let_executionTime);
  aoop_outParamGAC.setNumClusterK
    ((T_CLUSTERIDX)lochrom_best.getNumRows());
  aoop_outParamGAC.setMetricFuncRun
    (lochrom_best.getObjetiveFunc());
  aoop_outParamGAC.setFitness
    (lochrom_best.getFitness());
  aoop_outParamGAC.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoop_outParamGAC.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */ 
#ifndef __WITHOUT_PLOT_STAT

  if ( aiinp_inParamPcPmVk.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinp_inParamPcPmVk,
       aoop_outParamGAC
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA 
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochrom_best.print();
    std::cout << std::endl;

    auto lpartitionCentroids_clustersChromBest = 
      partition::makePartition
      (lomatrixrow_centroidsChromBest,
       aiiterator_instfirst,
       aiiterator_instlast,
       T_CLUSTERIDX(lomatrixrow_centroidsChromBest.getNumRows()),
       aifunc2p_dist
       );
    	     
    std::ostringstream lostrstream_labelShipBest;
    lostrstream_labelShipBest
      << "<MEMBERCLUSTER:" << lpc_labelAlgGA
      << "generation "
      <<  llfh_listFuntionHist.getDomainUpperBound()
      << ":lpartitionCentroids_clustersChromBest<>["
      << &lpartitionCentroids_clustersChromBest << ']';
    lpartitionCentroids_clustersChromBest.print
      (std::cout,lostrstream_labelShipBest.str().c_str(),',');
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lomatrixrow_centroidsChromBest; 
 
} /* END gcuk_vkcentroid */

} /*END eac */

#endif /*__GCUK_VKCENTROID_HPP__*/
