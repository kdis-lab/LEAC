/*! \file cbga_fkcentroid.hpp
 *
 * \brief  CBGA \cite Franti:etal:GAclustering:gafranti:1997
 *
 * \details  This file is part of the LEAC.\n\n
 * Implementation of the CBGA algorithm based on the paper:\n
 * Pasi Franti, Juha Kivijarvi, Timo Kaukoranta, and Olli Nevalainen.\n 
 * Genetic algorithms for large scale clustering problems.\n
 * Comput. J, 40:547–554, 1997
 *
 * Library Evolutionary Algorithms for Clustering (LEAC) is a library\n
 * for the implementation of evolutionary and genetic algorithms\n
 * focused on the partition type clustering problem. Based on the\n
 * current standards of the <a href="http://en.cppreference.com">C++</a> language, as well as on Standard\n
 * Template Library <a href="http://en.cppreference.com/w/cpp/container">STL</a> 
 * and also  <a href="http://www.openblas.net/">OpenBLAS</a> to have a better performance.\n
 *
 * \brief  
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CBGA_FKCENTROID_HPP__
#define __CBGA_FKCENTROID_HPP__

#include <iostream>
#include <functional>   // std::less
#include <algorithm>    // std::sort

#include <leac.hpp>
#include "inparam_cbga.hpp"
#include "outparam_gaclustering.hpp"

#include "plot_runtime_function.hpp"

/*! \namespace eac
  \brief Evolutionary Algorithms for Clustering
  \details Implementation of genetic and evolutionary algorithms used to solve the clustering problem 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace eac {

#define GACLUSTERING_CBGA_SELECMETH_ROULETTE 0
#define GACLUSTERING_CBGA_SELECMETH_ELITIST1 1
#define GACLUSTERING_CBGA_SELECMETH_ELITIST2 2
#define GACLUSTERING_CBGA_SELECMETH_ZIGZAG   3

#define GACLUSTERING_CBGA_SELECMETH {"roulette", "elitist1", "elitist2", "zigzag", (char *) NULL }

template <typename T_FEATURE,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename T_REAL,
	  typename INPUT_ITERATOR
	  >
void 
iterateGLAAux
(gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM, 
 T_REAL>                            &aochrom_cbga,
 INPUT_ITERATOR                     aiiterator_instfirst,


 const INPUT_ITERATOR               aiiterator_instlast,
 const int                          aii_numIteration, 
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 );


/*! \fn uintidx getCrossSetSizeAux (const int aii_opSelectMethod, const uintidx aist_populationSize)
   \brief  Cross size 
   \details Calculate the amount of chromosomes to crossover
   \return an integer number with the number of individuals to cross
   \param aii_opSelectMethod ￼type of selected method: GACLUSTERING_CBGA_SELECMETH_ROULETTE, GACLUSTERING_CBGA_SELECMETH_ELITIST1, GACLUSTERING_CBGA_SELECMETH_ZIGZAG or GACLUSTERING_CBGA_SELECMETH_ELITIST2
   \param aist_populationSize an integer with the size of the population
*/
uintidx
getCrossSetSizeAux
(const int     aii_opSelectMethod, 
 const uintidx aist_populationSize
 )
{
  uintidx lost_crossSetSize = 0;

  switch( aii_opSelectMethod ) {
  case GACLUSTERING_CBGA_SELECMETH_ROULETTE:
    lost_crossSetSize = aist_populationSize;
    break;
  case GACLUSTERING_CBGA_SELECMETH_ELITIST1:
  case GACLUSTERING_CBGA_SELECMETH_ZIGZAG:
    lost_crossSetSize=1; 
    while ( (lost_crossSetSize*(lost_crossSetSize-1)/2) < aist_populationSize ) 
      lost_crossSetSize++;
    break;
  case GACLUSTERING_CBGA_SELECMETH_ELITIST2:
    lost_crossSetSize=1; 
    while ( (lost_crossSetSize*(lost_crossSetSize+1)/2) < aist_populationSize ) 
      lost_crossSetSize++;
    break;
  default:
    break;
  }

  return lost_crossSetSize;

}

/*! \fn uintidx getSurvivorsAux (int aii_opSelectMethod, uintidx aist_populationSize)  
   \brief   Get survivors 
   \details Calculate the number of selected individuals in an elitist manner
   \param aii_opSelectMethod type of selected method
   \param aist_populationSize an integer with the size of the population
*/
uintidx
getSurvivorsAux
(int     aii_opSelectMethod, 
 uintidx aist_populationSize
 )
{
  uintidx lost_survivors = 0;

  switch( aii_opSelectMethod ) {
  case GACLUSTERING_CBGA_SELECMETH_ROULETTE:
    lost_survivors = 1;
    break;
  case GACLUSTERING_CBGA_SELECMETH_ELITIST1:
  case GACLUSTERING_CBGA_SELECMETH_ZIGZAG:
    lost_survivors = 1;
    break;
  case GACLUSTERING_CBGA_SELECMETH_ELITIST2:
    lost_survivors = 1; 
    while ( (lost_survivors*(lost_survivors+1)/2) < aist_populationSize ) 
      lost_survivors++;
    break;
  default:
    break;
  }

  return lost_survivors;
}


/*! \fn void changeOptimalityAux (gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX, T_INSTANCE_FREQUENCY, T_INSTANCES_CLUSTER_K, T_FEATURE_SUM, T_REAL> &aochrom_cbga, gaencode::OptimalityCBGA aiooptimalitycbga_a, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist)
  \brief Change optimality 
  \details Does a local search depending on the status of the chromosome
  \param aochrom_cbga a chromosome of type gaencode::ChromosomeCBGA
  \param aiooptimalitycbga_a gaencode::OptimalityCBGA
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distance
*/
template <typename T_FEATURE,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM, 
	  typename T_REAL,
	  typename INPUT_ITERATOR
	  >
void 
changeOptimalityAux
(gaencode::ChromosomeCBGA
 <T_FEATURE, 
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL>                             &aochrom_cbga,
 gaencode::OptimalityCBGA            aiooptimalitycbga_a,
 INPUT_ITERATOR                      aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char   *las_optimalityCBGA[] = OPTIMALITY_CBGA_STRING;

  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "eac::changeOptimalityAux:  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "(output gaencode::ChromosomeCBGA: aochrom_cbga["
	      << &aochrom_cbga << "]\n"
	      << " input gaencode::ChromosomeCBGA.optimalityCBGA = " 
	      <<  las_optimalityCBGA[ aochrom_cbga.getOptimalityCBGA() ] <<  '\n'
	      << " input  OptimalityCBGA:  aiooptimalitycbga_a = " 
	      << las_optimalityCBGA[ aiooptimalitycbga_a ] << '\n'
              << " input  aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      	      << &aifunc2p_dist << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  if ((aochrom_cbga.getOptimalityCBGA() != aiooptimalitycbga_a) ||
      (aochrom_cbga.getOptimalityCBGA() != gaencode::OPT_BOTH) 
      )  {

    switch ( aiooptimalitycbga_a )
      {
      case gaencode::OPT_NONE:
	break;
      case gaencode::OPT_CB:
	{
	  aochrom_cbga.optimalCodebook();
	  aochrom_cbga.setOptimalityCBGA( gaencode::OPT_CB );
	  break;
	}
      case gaencode::OPT_PA:
	{
	  clusteringop::reassignCluster
	    (aochrom_cbga.getPartition(),
	     aochrom_cbga.getCodeBook(),
	     aiiterator_instfirst,
             aiiterator_instlast,
	     aifunc2p_dist 
	     );
	  
	  aochrom_cbga.setOptimalityCBGA( gaencode::OPT_PA );
	  break;
	}
      case gaencode::OPT_BOTH:
	{ 
	  iterateGLAAux
	    (aochrom_cbga,
	     aiiterator_instfirst,
             aiiterator_instlast,
	     GENETICOPCLUSTER_MAX_GLA,
	     aifunc2p_dist
	     );	
	  aochrom_cbga.setOptimalityCBGA( gaencode::OPT_BOTH );
	  break;
	}
      }
  }
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "eac::changeOptimalityAux: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t input  OptimalityCBGA:  aiooptimalitycbga_a = " 
	      << las_optimalityCBGA[ aochrom_cbga.getOptimalityCBGA() ] << '\n';
	
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
}

/*! \fn void iterateGLAAux (gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_FREQUENCY_SUM,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM,T_REAL> &aochromfr_chromosome, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const int aii_numIteration, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief iterateGLAAux 
  \details
  \param aochromfr_chromosome a gaencode::ChromosomeCBGA
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aii_numIteration
  \param aifunc2p_dist an object of type dist::Dist to calculate distance
*/
template <typename T_FEATURE, 
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename T_REAL,
          typename INPUT_ITERATOR	  
	  >
void 
iterateGLAAux
(gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM, 
 T_REAL>                            &aochrom_cbga,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const int                          aii_numIteration, 
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{
  int           li_round = 0;
  T_CLUSTERIDX  *lcidx_newindices;   
  uintidx       lst_countNewIndex = 1; 


#ifdef __VERBOSE_YES
  const char   *las_optimalityCBGA[] = OPTIMALITY_CBGA_STRING;

  const char* lpc_labelFunc = "eac::iterateGLA";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
  	      << "(output gaencode::ChromosomeCBGA: aochrom_cbga["
	      << &aochrom_cbga << "]\n"
	      << " input ChromosomeCBGA.optimalityCBGA = " 
	      <<  las_optimalityCBGA[ aochrom_cbga.getOptimalityCBGA() ] <<  '\n'
	      << " input  aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input  int: aii_numIteration: " << aii_numIteration << '\n'
	      << " input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      	      << &aifunc2p_dist << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  if (aii_numIteration < 1) {

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "eac::iterateGLA: OUT (aii_numIteration < 1)\n"
		<< "li_round: " << li_round 
		<< "\taii_numIteration: " << aii_numIteration  
		<< "\tlst_countNewIndex: " << lst_countNewIndex
		<< '\n'
		<< "\toutput gaencode::ChromosomeCBGA: aochrom_cbga["
		<< &aochrom_cbga << "]\n";
    }
    return;
#endif /*__VERBOSE_YES*/

  }
  
  switch ( aochrom_cbga.getOptimalityCBGA() )
    {
    case gaencode::OPT_NONE: /* start by generating an optimal partitioning */
    case gaencode::OPT_CB:
      {
	changeOptimalityAux
	  (aochrom_cbga, 
	   gaencode::OPT_PA,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
	
	iterateGLAAux
	  (aochrom_cbga,
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aii_numIteration-1,
	   aifunc2p_dist
	   );
	
	changeOptimalityAux
	  (aochrom_cbga, 
	   gaencode::OPT_CB, 
	   aiiterator_instfirst,
	   aiiterator_instlast, 
	   aifunc2p_dist
	   );
	
	clusteringop::fillEmptyPartitions
	  (aochrom_cbga.getPartition(),
	   aochrom_cbga.getCodeBook(),
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   aifunc2p_dist
	   );
	
	break;
      }
    case gaencode::OPT_PA: //K-MEAN
      {  
	lcidx_newindices = 
	  new T_CLUSTERIDX[aochrom_cbga.getCodeBook().getNumRows()];
	li_round = 0;
	while ( (li_round < aii_numIteration)  && (lst_countNewIndex > 0) ) {
	  
	  mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
	    lmT_codeBook(aochrom_cbga.getCodeBook());

	  T_CLUSTERIDX lcidx_numClusterNull; 
	  
	  clusteringop::meanCentroids
	    (lcidx_numClusterNull,
	     aochrom_cbga.getCodeBook(),
	     aochrom_cbga.getPartition().getSumInstancesClusterK(),
	     aochrom_cbga.getPartition().getNumInstancesClusterK()
	     );
	  
	  clusteringop::fillEmptyPartitions
	    (aochrom_cbga.getPartition(),
	     aochrom_cbga.getCodeBook(),
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     aifunc2p_dist
	     );
	  
	  aochrom_cbga.checkChangedCodeVectors
	    (lcidx_newindices,
	     &lst_countNewIndex,
	     lmT_codeBook 
	     );
	  
	  clusteringop::reassignCluster
	    (aochrom_cbga.getPartition(),
	     aochrom_cbga.getCodeBook(),
	     lcidx_newindices,
	     lst_countNewIndex,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     aifunc2p_dist
	     );
	  
	  li_round++;
	  
	} /*END WHILE*/
	delete[] lcidx_newindices;
	
	break;
      }
    case gaencode::OPT_BOTH:
      break;
    }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "\tli_round: " << li_round 
	      << "\taii_numIteration: " << aii_numIteration  
	      << "\tlst_countNewIndex: " << lst_countNewIndex
	      << " gaencode::ChromosomeCBGA.optimalityCBGA = " 
	      <<  las_optimalityCBGA[ aochrom_cbga.getOptimalityCBGA() ] <<  '\n';
    aochrom_cbga.print(std::cout,lpc_labelFunc,',',';');
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*! \fn gaencode::ChromosomeCBGA <T_FEATURE, T_CLUSTERIDX, T_INSTANCE_FREQUENCY, _INSTANCES_CLUSTER_K, T_FEATURE_SUM, T_REAL>  cbga_fkcentroid (inout::OutParamGAClustering<T_REAL,T_CLUSTERIDX> &aoopcga_outParamClusteringGA, const inout::InParamCBGA<T_CLUSTERIDX,T_REAL,T_FEATURE, T_FEATURE_SUM, T_INSTANCES_CLUSTER_K, T_INSTANCE_FREQUENCY> &aiinParam_CBGA, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief CBGA \cite Franti:etal:GAclustering:gafranti:1997
  \details Implementation of the CBGA algorithm based on \cite Franti:etal:GAclustering:gafranti:1997.  
  \returns A partition of a data set, encoded on a chromosome where each gene is the coordinate of a centroid or a codebook.
  \param aoopcga_outParamClusteringGA a inout::OutParamGAClustering with the output parameters of the algorithm
  \param aiinParam_CBGA a inout::InParamGAClusteringProbCProbMFixedK parameters required by the algorithm
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances  
*/
template < typename T_FEATURE,
	   typename T_REAL,
	   typename T_INSTANCE_FREQUENCY,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K, //INSTANCES FREQUENCY SUM
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
gaencode::ChromosomeCBGA
<T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL
 >
cbga_fkcentroid
(inout::OutParamGAClustering
 <T_REAL,
 T_CLUSTERIDX>                      &aoopcga_outParamClusteringGA,
 const inout::InParamCBGA
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K,
 T_INSTANCE_FREQUENCY>              &aiinParam_CBGA,
 const INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{ /* BEGIN cbga_fkcentroid */

  const uintidx lconstui_numClusterFixedK =
    (uintidx) aiinParam_CBGA.getNumClusterK();

  const uintidx lconstui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));

  std::uniform_int_distribution<uintidx> uniformdis_idxInstances
    (0,lconstui_numInstances-1);
  
  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
  /*VARIABLE NEED FOR POPULATION AND MATINGPOOL GENETIC*/
  std::vector   
    <gaencode::ChromosomeCBGA
     <T_FEATURE,
      T_CLUSTERIDX, //-1, 0, 1, .., K
      T_INSTANCE_FREQUENCY,
      T_INSTANCES_CLUSTER_K,    
      T_FEATURE_SUM,
      T_REAL>* >                           
    lvectorchromcbga_populationNew;
  
  std::vector   
    <gaencode::ChromosomeCBGA
     <T_FEATURE,
      T_CLUSTERIDX, //-1, 0, 1, .., K
      T_INSTANCE_FREQUENCY,
      T_INSTANCES_CLUSTER_K,    
      T_FEATURE_SUM,
      T_REAL>* >                           
    lvectorchromcbga_populationOld;
  
  gaencode::ChromosomeCBGA 
    <T_FEATURE,
     T_CLUSTERIDX,
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,    
     T_FEATURE_SUM,
     T_REAL> 
    lochromcbga_best
    (aiinParam_CBGA.getNumClusterK(),
     aiinParam_CBGA.getNumClusterK(),
     aiinParam_CBGA.getNumDimensionsInstances(),
     aiinParam_CBGA.getNumInstances()
     );

  uintidx luintidx_crossSetSize = 
    getCrossSetSizeAux
    (aiinParam_CBGA.getOpSelectMethod(),
     aiinParam_CBGA.getSizePopulation()
     );
  uintidx luintidx_survivors = 
    getSurvivorsAux
    (aiinParam_CBGA.getOpSelectMethod(),
     aiinParam_CBGA.getSizePopulation()
     );


#ifdef __VERBOSE_YES

  /*ID PROC
   */
  geverboseui_idproc = 1;
  
  ++geiinparam_verbose;
  const char* lpc_labelAlgGA = "cbga_fkcentroid" ; 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelAlgGA
	      << ":  IN(" << geiinparam_verbose << ")\n"
              << "\t(output gaencode::ChromosomeCBGA: lochromcbga_best[" 
	      << &lochromcbga_best << "]\n"
	      << "\t output inout::OutParamGAClustering&: aoopcga_outParamClusteringGA[" 
	      << &aoopcga_outParamClusteringGA << "]\n"
	      << "\t input  InParamCBGA&: aiinParam_CBGA[" 
	      << &aiinParam_CBGA << "]\n"
              << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  &aiiterator_instlast << "]\n"
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t input GA parameters: "
	      << "\n\t\tPopulation size = " 
	      << aiinParam_CBGA.getSizePopulation()
	      << "\n\t\tCross set size = " << luintidx_crossSetSize 
	      << "\n\t\tSurvivors      = " << luintidx_survivors
	      << "\n\t\tMutation       = " << aiinParam_CBGA.getProbMutation()
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  
  runtime::ListRuntimeFunction<COMMON_IDOMAIN> 
    llfh_listFuntionHist
    (aiinParam_CBGA.getNumMaxGenerations(), "Iterations", "Clustering metrics");

  /*DECLARATION OF VARIABLES: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */
#ifndef __WITHOUT_PLOT_STAT
  
  std::ofstream             lfileout_plotStatObjetiveFunc;
  runtime::RuntimeFunctionValue<T_REAL> *lofh_distortion = NULL;
  runtime::RuntimeFunctionStat<T_REAL>  *lofhs_statObjectiveFunc[STATISTICAL_ALL_MEASURES];
  std::vector<T_REAL>       lvectorT_statfuncObjetiveFunc;

  if ( aiinParam_CBGA.getWithPlotStatObjetiveFunc() ) {  
    
    lvectorT_statfuncObjetiveFunc.reserve
      (aiinParam_CBGA.getSizePopulation());
    //DEFINE FUNCTION
    lofh_distortion  = new runtime::RuntimeFunctionValue<T_REAL>
      ("Distortion", 
       aiinParam_CBGA.getAlgorithmoName(),
       RUNTIMEFUNCTION_NOT_STORAGE
       );

    llfh_listFuntionHist.addFuntion(lofh_distortion);
    //DEFINE FUNCTION STATISTICAL
    for  (int li_i = 0; li_i < STATISTICAL_ALL_MEASURES; li_i++) {
      lofhs_statObjectiveFunc[li_i] = new runtime::RuntimeFunctionStat
	<T_REAL>
	( (char) li_i,
	  aiinParam_CBGA.getAlgorithmoName(),
	  RUNTIMEFUNCTION_NOT_STORAGE
	  );
      llfh_listFuntionHist.addFuntion(lofhs_statObjectiveFunc[li_i]);
    }
  
    //OPEN FILE STRORE FUNCTION
    aoopcga_outParamClusteringGA.setFileNameOutPlotStatObjetiveFunc
      (aiinParam_CBGA.getFileNamePlotStatObjetiveFunc(),
       aiinParam_CBGA.getTimesRunAlgorithm()
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
   
  /*WHEN CAN MEASURE STARTS AT ZERO INVALID OFFSPRING*/
  aoopcga_outParamClusteringGA.setTotalInvalidOffspring(0);

  runtime::ExecutionTime let_executionTime = runtime::start();

  { /*BEGIN 1. GENERATE S RANDOM SOLUTIONS FOR THE INITIAL GENERATION
      GENERATE CODEBOOK  AND PARTITIONING
    */
  
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "1. GENERATE S RANDOM SOLUTIONS";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    /*BEGIN POPULATION CREATE
     */
    lvectorchromcbga_populationNew.reserve
      (aiinParam_CBGA.getSizePopulation());

    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinParam_CBGA.getSizePopulation(); 
	 luintidx_i++) 
      {
      
	lvectorchromcbga_populationNew.push_back
	  (new gaencode::ChromosomeCBGA
	   <T_FEATURE,
	   T_CLUSTERIDX, //-1, 0, 1, .., K
	   T_INSTANCE_FREQUENCY,
	   T_INSTANCES_CLUSTER_K,    
	   T_FEATURE_SUM,
	   T_REAL> 
	   (lconstui_numClusterFixedK,
	    2 * lconstui_numClusterFixedK,
	    data::Instance<T_FEATURE>::getNumDimensions(),
	    lconstui_numInstances
	    )
	   );
      
      }
       
    /*CREATE SPACE FOR STORE MATINGPOOL--------------------------------------------
     */
    lvectorchromcbga_populationOld.reserve
      (aiinParam_CBGA.getSizePopulation());
    for (uintidx luintidx_i = 0; 
	 luintidx_i < aiinParam_CBGA.getSizePopulation(); 
	 luintidx_i++) 
      {
	lvectorchromcbga_populationOld.push_back
	  (new gaencode::ChromosomeCBGA
	   <T_FEATURE,
	   T_CLUSTERIDX, //-1, 0, 1, .., K
	   T_INSTANCE_FREQUENCY,
	   T_INSTANCES_CLUSTER_K,    
	   T_FEATURE_SUM,
	   T_REAL>
	   (lconstui_numClusterFixedK,
	    lconstui_numClusterFixedK,
	    data::Instance<T_FEATURE>::getNumDimensions(),
	    lconstui_numInstances
	    )
	   );
      }  
    /*END POPULATION CREATE
     */

    /*BEGIN 1. GENERATE S RANDOM SOLUTIONS
     */
    long ll_invalidOffspring = 0;
    for ( auto literchrom_cbga: lvectorchromcbga_populationNew ) {
       
      clusteringop::randomInitialize
	(literchrom_cbga->getCodeBook(), 
	 aiiterator_instfirst,
	 aiiterator_instlast
	 );

      /*Generate optimal Partitioning
       */
      std::pair<bool,T_REAL> lpair_generate  =  
	clusteringop::reassignCluster
	(literchrom_cbga->getPartition(),
	 literchrom_cbga->getCodeBook(),
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 aifunc2p_dist
	 );
      literchrom_cbga->setValidString(lpair_generate.first);
      literchrom_cbga->setObjetiveFunc(lpair_generate.second);

      if ( literchrom_cbga->getValidString() == false ) 
	++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
      lvectorT_statfuncObjetiveFunc.push_back(literchrom_cbga->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/
	
    }    
    /*METRIC INVALID SOLUTION
     */
    aoopcga_outParamClusteringGA.sumTotalInvalidOffspring
      (ll_invalidOffspring);

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
   
    /*END 1. GENERATE S RANDOM SOLUTIONS
     */
    
  } /*END 1. GENERATE S RANDOM SOLUTIONS FOR THE INITIAL GENERATION
      GENERATE CODEBOOK  AND PARTITIONING
    */
  

  /*SORT POPULATION
   */
  {/*BEGIN SORT POPULATION*/
    
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "GENERATE S RANDOM AND SORT";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    std::sort
      (lvectorchromcbga_populationNew.begin(),
       lvectorchromcbga_populationNew.end(),
       [](const typename decltype(lvectorchromcbga_populationNew)::value_type literchrom_cbgaA,
	  const typename decltype(lvectorchromcbga_populationNew)::value_type literchrom_cbgaB
	  )
       {
	 return ( literchrom_cbgaA->getObjetiveFunc() < literchrom_cbgaB->getObjetiveFunc() );
       }
       );
     
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;    
      for ( auto literchrom_cbga: lvectorchromcbga_populationNew ) {

	//TEST CALCULATE METRIC----------------------------------------------------------
	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (literchrom_cbga->getPartition().getMembersShip(),
	   lconstui_numInstances,
	   aiinParam_CBGA.getNumClusterK()
	   );
	
	std::vector<T_REAL>      lvectorrt_sumDistEuSqInstInClusterK =
	  um::sumDistInstCentInK
	  (literchrom_cbga->getCodeBook(),
	   aiiterator_instfirst,
	   aiiterator_instlast,
	   lpartition_clusters,
	   aifunc2p_dist,
	   [](const data::Instance<T_FEATURE>* aiinst_iter ) -> T_REAL
	   {
	     data::InstanceFreq
	       <T_FEATURE,
		T_INSTANCE_FREQUENCY
		>
	       *lptinstfreq_iter =
	       (data::InstanceFreq
		<T_FEATURE,
		T_INSTANCE_FREQUENCY
		>*)
	       aiinst_iter;

	     return  (T_REAL) lptinstfreq_iter->getFrequency(); 
	   }
	   );
	T_REAL lmetrict_SSENormEuSq = 
	  interfacesse::sum
	  (lvectorrt_sumDistEuSqInstInClusterK.data(),
	   (uintidx) lvectorrt_sumDistEuSqInstInClusterK.size()
	   );

	T_REAL lrt_distortion =
	  lmetrict_SSENormEuSq / 
	  ((T_REAL) lconstui_numInstances * 
	   (T_REAL) data::Instance<T_FEATURE>::getNumDimensions());
        //TEST CALCULATE METRIC----------------------------------------------------------
	
	std::cout << geverbosepc_labelstep
		  << ":gaencode::ChromosomeCBGA"  << ":id[" << geverboseui_idproc
		  << ':' << literchrom_cbga << ']'
		  << ":objetive function: " << literchrom_cbga->getObjetiveFunc()
		  << " :lrt_distortion: " << lrt_distortion
		  << '\n';
      } //for
    } //if 
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }/*END SORT POPULATION*/

  /*ONLY FOR MEASURE NOT IS PART ALGORITHM*/
  auto lT_bestObjeticeFunc = lvectorchromcbga_populationNew.at(0)->getObjetiveFunc();
  aoopcga_outParamClusteringGA.setIterationGetsBest
    (llfh_listFuntionHist.getDomainUpperBound());
  aoopcga_outParamClusteringGA.setRunTimeGetsBest
    (runtime::elapsedTime(let_executionTime));

  /*MEASUREMENT NEW GENERATION: COMPUTING STATISTICAL AND METRIC OF THE--------- 
    ALGORITHM*/
#ifndef __WITHOUT_PLOT_STAT  
  if ( aiinParam_CBGA.getWithPlotStatObjetiveFunc() ) {  
    lofh_distortion->setValue
      (lvectorchromcbga_populationNew.at(0)->getObjetiveFunc());
    functionhiststat_evaluateAll
      (lofhs_statObjectiveFunc,
       lvectorT_statfuncObjetiveFunc
       );
    lfileout_plotStatObjetiveFunc << llfh_listFuntionHist;
    lvectorT_statfuncObjetiveFunc.clear();
  }
#endif  /*__WITHOUT_PLOT_STAT */


  while ( llfh_listFuntionHist.getDomainUpperBound() < 
	  aiinParam_CBGA.getNumMaxGenerations()) {
   
    llfh_listFuntionHist.increaseDomainUpperBound();
    
    
    {/* COPY NEW GENERATION */

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "COPY NEW GENERATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      auto iIter = lvectorchromcbga_populationNew.begin();
      auto jIter = lvectorchromcbga_populationOld.begin();
  
      for (;iIter != lvectorchromcbga_populationNew.end();++iIter, ++jIter) {
	*(*jIter) = *(*iIter);
      }

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
      
    } /* END COPY */

    /*BEGIN: GENERATE NEXT GENERATION
     */

    /*CROSSOVER-----------------------------------------------------------------
     */
    { /*BEGIN: CROOSSOVER
       */
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "2.2 GENERATE NEW SOLUTIONS BY CROSSOVER";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << " luintidx_survivors = "
		  << luintidx_survivors
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
    
      std::vector<uintidx> lvectorst_idxChromPairsCross;
      lvectorst_idxChromPairsCross.reserve(2);

      for ( auto literchrom_cbgaNew =
	      std::next(lvectorchromcbga_populationNew.begin(),luintidx_survivors);
	    literchrom_cbgaNew != lvectorchromcbga_populationNew.end();
	    ++literchrom_cbgaNew)
	{
	  std::pair<uintidx,uintidx> 
	    lpairst_pairCross =
	    prob::getElitistPairCross
	    (&lvectorst_idxChromPairsCross,
	     (uintidx) 0,
	     luintidx_crossSetSize
	     );
	  auto lchf_chromOld1 = 
	    lvectorchromcbga_populationOld.at(lpairst_pairCross.first);
	  auto lchf_chromOld2 = 
	    lvectorchromcbga_populationOld.at(lpairst_pairCross.second);
	   
	  gaclusteringop::crossPNNnew
	    (*(*literchrom_cbgaNew),
	     *lchf_chromOld1,
	     *lchf_chromOld2,
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     aiinParam_CBGA.getNumClusterK(),
	     aifunc2p_dist
	     );
	}
      
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
        
    } /*END: CROOSSOVER*/
 
    /*MUTATE--------------------------------------------------------------------
     */
    { /*BEGIN: MUTATE*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "2.3 GENETATE MUTATIONS TO THE SOLUTIONS";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << " luintidx_survivors = "
		  << luintidx_survivors
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/

      for ( auto literchrom_cbgaNew =
	      std::next(lvectorchromcbga_populationNew.begin(),luintidx_survivors);
	    literchrom_cbgaNew != lvectorchromcbga_populationNew.end();
	    ++literchrom_cbgaNew)
	{
	  if ( uniformdis_real01(gmt19937_eng) < aiinParam_CBGA.getProbMutation() ) {

	    clusteringop::randomInitialize
	      ((*literchrom_cbgaNew)->getCodeBook(), 
	       aiiterator_instfirst,
	       aiiterator_instlast
	       );
	    (*literchrom_cbgaNew)->setOptimalityCBGA(gaencode::OPT_NONE);
	  }
	 
	  if ( aiinParam_CBGA.getNumGLAIterations() > 0 )  { 
	    iterateGLAAux
	      (*(*literchrom_cbgaNew),
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aiinParam_CBGA.getNumGLAIterations(),
	       aifunc2p_dist
	       );    
	  } /*END IF*/
	}
    
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
        
    } /*END: MUTATE*/


    { /*BEGIN  SURVIVORS*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "EVALUATE SSE";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": IN(" << geiinparam_verbose << ')'
		  << " luintidx_survivors = "
		  << luintidx_survivors
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
    
      long ll_invalidOffspring = 0;

      for ( auto literchrom_cbga =
	      std::next(lvectorchromcbga_populationNew.begin(),luintidx_survivors);
	    literchrom_cbga != lvectorchromcbga_populationNew.end();
	    ++literchrom_cbga)
	{	
	  std::pair<T_REAL,bool> 
	    lpair_distortion =
	    um::distortion
	    ((*literchrom_cbga)->getCodeBook(),
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     (*literchrom_cbga)->getPartition().getMembersShip(),
	     aifunc2p_dist,
	     [](const data::Instance<T_FEATURE>* aiinst_iter ) -> T_REAL
	     {
	       data::InstanceFreq
		 <T_FEATURE,
		  T_INSTANCE_FREQUENCY
		  >
		 *lptinstfreq_iter =
		 (data::InstanceFreq
		  <T_FEATURE,
		  T_INSTANCE_FREQUENCY
		  >*)
		 aiinst_iter;

	       return  (T_REAL) lptinstfreq_iter->getFrequency();
	   
	     }
	     );

	  if ( measuare_undefObjetiveFunc(T_REAL) == lpair_distortion.first ) {
	    std::pair<bool,T_REAL> lpair_generate  =  
	      clusteringop::reassignCluster
	      ((*literchrom_cbga)->getPartition(),
	       (*literchrom_cbga)->getCodeBook(),
	       aiiterator_instfirst,
	       aiiterator_instlast,
	       aifunc2p_dist
	       );
	    (*literchrom_cbga)->setValidString(lpair_generate.first);
	    (*literchrom_cbga)->setObjetiveFunc(lpair_generate.second);
	  }
	  else {
       
	    (*literchrom_cbga)->setObjetiveFunc(lpair_distortion.first);
	    (*literchrom_cbga)->setValidString(lpair_distortion.second);

	  }
		  
	  if ( (*literchrom_cbga)->getValidString() == false ) 
	    ++ll_invalidOffspring;

#ifndef __WITHOUT_PLOT_STAT
	  lvectorT_statfuncObjetiveFunc.push_back((*literchrom_cbga)->getObjetiveFunc());
#endif /*__WITHOUT_PLOT_STAT*/

	}

      /*METRIC INVALID SOLUTION*/
      aoopcga_outParamClusteringGA.sumTotalInvalidOffspring
	(ll_invalidOffspring);

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    } /*END  SURVIVOTS*/


    /*SORT POPULATION
     */
    {/*BEGIN SORT POPULATION*/

#ifdef __VERBOSE_YES
      geverbosepc_labelstep = "SORT POPULATION";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep  
		  << ": IN(" << geiinparam_verbose << ')'
		  << std::endl;
      }
#endif /*__VERBOSE_YES*/
    
      std::sort
	(lvectorchromcbga_populationNew.begin(),
	 lvectorchromcbga_populationNew.end(),
	 [](const typename decltype(lvectorchromcbga_populationNew)::value_type literchrom_cbgaA,
	    const typename decltype(lvectorchromcbga_populationNew)::value_type literchrom_cbgaB
	    )
	 {
	   return ( literchrom_cbgaA->getObjetiveFunc() < literchrom_cbgaB->getObjetiveFunc() );
	 }
	 );
   
#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << geverbosepc_labelstep
		  << ": OUT(" << geiinparam_verbose << ')'
		  << std::endl;    
	for ( auto literchrom_cbga: lvectorchromcbga_populationNew ) {
	  //TEST CALCULATE METRIC----------------------------------------------------------
	  partition::PartitionLabel
	    <T_CLUSTERIDX>
	    lpartition_clusters
	    (literchrom_cbga->getPartition().getMembersShip(),
	     lconstui_numInstances,
	     aiinParam_CBGA.getNumClusterK()
	     );
	  std::vector<T_REAL>      lvectorrt_sumDistEuSqInstInClusterK =
	    um::sumDistInstCentInK
	    (literchrom_cbga->getCodeBook(),
	     aiiterator_instfirst,
	     aiiterator_instlast,
	     lpartition_clusters,
	     aifunc2p_dist,
	     [](const data::Instance<T_FEATURE>* aiinst_iter ) -> T_REAL
	     {
	       data::InstanceFreq
		 <T_FEATURE,
		  T_INSTANCE_FREQUENCY
		  >
		 *lptinstfreq_iter =
		 (data::InstanceFreq
		  <T_FEATURE,
		  T_INSTANCE_FREQUENCY
		  >*)
		 aiinst_iter;

	       return  (T_REAL) lptinstfreq_iter->getFrequency(); 
	     }
	     );
	  T_REAL lmetrict_SSENormEuSq = 
	    interfacesse::sum
	    (lvectorrt_sumDistEuSqInstInClusterK.data(),
	     (uintidx) lvectorrt_sumDistEuSqInstInClusterK.size()
	     );

	  T_REAL lrt_distortion =
	    lmetrict_SSENormEuSq / 
	    ((T_REAL) lconstui_numInstances *
	     (T_REAL) data::Instance<T_FEATURE>::getNumDimensions());
	  //TEST CALCULATE METRIC-----------------------------------------------
	  std::cout << geverbosepc_labelstep
		    << ":gaencode::ChromosomeCBGA"
		    << ":id[" << geverboseui_idproc << ':'
		    << literchrom_cbga << ']'
		    << ":objetive function: " << literchrom_cbga->getObjetiveFunc()
		    << " :lrt_distortion: " << lrt_distortion
		    << '\n';
	}	
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    } /*END SORT POPULATION*/
    
    /*ONLY FOR MEASURE NOT IS PART ALGORITHM*/
    if ( lvectorchromcbga_populationNew.at(0)->getObjetiveFunc() < lT_bestObjeticeFunc) {
      lT_bestObjeticeFunc =  lvectorchromcbga_populationNew.at(0)->getObjetiveFunc();
      aoopcga_outParamClusteringGA.setIterationGetsBest
	(llfh_listFuntionHist.getDomainUpperBound());
      aoopcga_outParamClusteringGA.setRunTimeGetsBest
	(runtime::elapsedTime(let_executionTime));
    }

    /*MEASUREMENT BEST: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM*/

#ifndef __WITHOUT_PLOT_STAT  
    if ( aiinParam_CBGA.getWithPlotStatObjetiveFunc() ) {  
      lofh_distortion->setValue
	(lvectorchromcbga_populationNew.at(0)->getObjetiveFunc());
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
      std::cout << "END ITERATION: " << llfh_listFuntionHist.getDomainUpperBound()
		<< ":gaencode::ChromosomeCBGA"  << ":id[" << geverboseui_idproc
		<< ':' << lvectorchromcbga_populationNew.at(0) << ']'
		<< ":objetive function: " << lvectorchromcbga_populationNew.at(0)->getObjetiveFunc()
		<< std::endl;
    }
    
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /* END  While */

  /*FREE MEMORY*/

#ifdef __VERBOSE_YES
  geverbosepc_labelstep = "3. OUTPUT THE BEST SOLUTION OF THE FINAL GENERATION";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << geverbosepc_labelstep  
	      << ": IN(" << geiinparam_verbose << ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  lochromcbga_best = *(lvectorchromcbga_populationNew.at(0));

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout  << geverbosepc_labelstep
	       << ": OUT(" << geiinparam_verbose << ")\n";
    lochromcbga_best.print();
    std::cout << std::endl;

    //TEST CALCULATE METRIC-----------------------------------------------------------------------------------------
    partition::PartitionLabel
      <T_CLUSTERIDX>
      lpartition_clusters
      (lochromcbga_best.getPartition().getMembersShip(),
       lconstui_numInstances, 
       aiinParam_CBGA.getNumClusterK()
       );

    std::vector<T_REAL>      lvectorrt_sumDistEuSqInstInClusterK =
      um::sumDistInstCentInK
      (lochromcbga_best.getCodeBook(),
       aiiterator_instfirst,
       aiiterator_instlast,
       lpartition_clusters,
       aifunc2p_dist,
       [](const data::Instance<T_FEATURE>* aiinst_iter ) -> T_REAL
       {
	 data::InstanceFreq
	   <T_FEATURE,
	    T_INSTANCE_FREQUENCY
	    >
	   *lptinstfreq_iter =
	   (data::InstanceFreq
	    <T_FEATURE,
	    T_INSTANCE_FREQUENCY
	    >*)
	   aiinst_iter;

	 return  (T_REAL) lptinstfreq_iter->getFrequency(); 
       }
       );
    T_REAL lmetrict_SSENormEuSq = 
      interfacesse::sum
      (lvectorrt_sumDistEuSqInstInClusterK.data(),
       (uintidx) lvectorrt_sumDistEuSqInstInClusterK.size()
       );
    T_REAL lrt_distortion =
      lmetrict_SSENormEuSq / 
      ((T_REAL) lconstui_numInstances *  
       (T_REAL) data::Instance<T_FEATURE>::getNumDimensions());

    std::cout << "cbga_fkcentroid.hpp "
	      << "objetive function best = "
	      << lochromcbga_best.getObjetiveFunc()
	      << " lrt_distortion  ="
	      << lrt_distortion
	      << std::endl;
      
    //TEST CALCULATE METRIC-----------------------------------------------------
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
          
  for ( auto  literchrom_cbga: lvectorchromcbga_populationNew ) 
    delete literchrom_cbga;

  for ( auto  literchrom_cbga: lvectorchromcbga_populationOld ) 
    delete literchrom_cbga;
         		
  runtime::stop(let_executionTime);
  aoopcga_outParamClusteringGA.setNumClusterK
    (aiinParam_CBGA.getNumClusterK());
  aoopcga_outParamClusteringGA.setAlgorithmRunTime
    (runtime::getTime(let_executionTime));
  aoopcga_outParamClusteringGA.setMetricFuncRun
    (lochromcbga_best.getObjetiveFunc());

  aoopcga_outParamClusteringGA.setFitness
    (OUTPARAMCLUSTERING_FITNESS_NaN);
  aoopcga_outParamClusteringGA.setNumTotalGenerations
    (llfh_listFuntionHist.getDomainUpperBound());

  /*FREE: COMPUTING STATISTICAL AND METRIC OF THE ALGORITHM
   */

#ifndef __WITHOUT_PLOT_STAT

  if ( aiinParam_CBGA.getWithPlotStatObjetiveFunc() ) {  
    plot_funtionHist
      (llfh_listFuntionHist,
       aiinParam_CBGA,
       aoopcga_outParamClusteringGA
       );  
  }

#endif /*__WITHOUT_PLOT_STAT*/


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {    
    geverbosepc_labelstep = lpc_labelAlgGA;
    std::cout << lpc_labelAlgGA
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochromcbga_best.print(); 
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


  return  lochromcbga_best; 
 
} /* END cbga_fkcentroid */

} /*END eac */

#endif /*__CBGA_FKCENTROID_HPP__*/
