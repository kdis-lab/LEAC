/*! \file supervised_measures.hpp
 *
 * \brief Supervised measures
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef SUPERVISED_MEASURE_HPP
#define SUPERVISED_MEASURE_HPP


#include "matching_matrix.hpp"
#include "verbose_global.hpp"

/*! \namespace sm
  \brief  Supervised measures 
  \details For a given algorithm's are measures that evaluate the quality of the cluster for previously classified instances
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace sm {


#define measuare_undefRandIndex(T_METRIC) (T_METRIC) 0.0
  
/*! \fn inline T_METRIC randIndex(T_INSTANCES_CLUSTER_K aiT_agreementObjects, T_INSTANCES_CLUSTER_K aiT_numObjetos)
    \brief Get Rand Index \cite Rand:ClusterAnalysis:RandIndex:1971 \cite humber:arabie:clusteranalysis:randindex:1988 
    \details The Rand Index is defined for two partitions of the same data set \f$X\f$, \f$C\f$ in \f$k\f$ cluster and \f$R\f$ in \f$k'\f$ known clases.

\f[ 
\Omega(R,C)= {a+d \over a+b+c+d }
\f]

Where:
-\f$a\f$: Number of pairs of data objects belonging to the same class in \f$R\f$ and to the same cluster in \f$C\f$.
-\f$b\f$: Number of pairs of data objects belonging to the same class in \f$R\f$ yet to different clusters in \f$C\f$.
-\f$c\f$: Number of pairs of data objects belonging to different \f$R\f$ yet to the same cluster in \f$C\f$.
-\f$d\f$: Number of pairs of data objects belonging to different classes in \f$R\f$ and to different clusters in \f$C\f$.

    \param  aimatchmatrix_confusion a sm::ConfusionMatchingMatrix 
 */
template < typename T_METRIC, 
	   typename T_INSTANCES_CLUSTER_K 
	   >
T_METRIC
randIndex
(const ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "sm::randIndex";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "aimatchmatrix_confusion [ = " << &aimatchmatrix_confusion << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  T_INSTANCES_CLUSTER_K lit_pairsSomeClassUSomeClusterV_a =
    aimatchmatrix_confusion.getSomeClassUSomeClusterV();
  T_INSTANCES_CLUSTER_K lit_pairsDiffClassUDiffClusterV_d =
    aimatchmatrix_confusion.getDiffClassUDiffClusterV();
  T_INSTANCES_CLUSTER_K lit_agreementObjects =  lit_pairsSomeClassUSomeClusterV_a + lit_pairsDiffClassUDiffClusterV_d;

  T_INSTANCES_CLUSTER_K lit_numObjetos =
    aimatchmatrix_confusion.getNumObjetos();
  
  T_METRIC lrt_randIndex =
    ((lit_numObjetos != 0) && (lit_agreementObjects>=0) )?          //BINOMIAL COEFFICIENT 
    T_METRIC(lit_agreementObjects) / (T_METRIC( (lit_numObjetos-1) * lit_numObjetos / 2.0)):
    measuare_undefRandIndex(T_METRIC);

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << "\tlit_agreementObjects = " << lit_agreementObjects  << "\tlrt_randIndex = " << lrt_randIndex 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lrt_randIndex;
   
}

/*Jaccard index is quite similar to Rand index, but
	  this measure does not consider the number of correct assign-
	  ments when two elements are assigned to different clusters.
 */
template < typename T_METRIC, 
	   typename T_INSTANCES_CLUSTER_K 
	   >
T_METRIC
jaccardIndex
(const ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
{
  
  T_INSTANCES_CLUSTER_K lit_pairsSomeClassUSomeClusterV_a =
    aimatchmatrix_confusion.getSomeClassUSomeClusterV();
  T_INSTANCES_CLUSTER_K lit_pairsDiffClassUDiffClusterV_d =
    aimatchmatrix_confusion.getDiffClassUDiffClusterV();

  T_INSTANCES_CLUSTER_K lit_numObjetos =
    aimatchmatrix_confusion.getNumObjetos();
  
  return ( lit_numObjetos != 0 && lit_pairsSomeClassUSomeClusterV_a>=0)?
    T_METRIC(lit_pairsSomeClassUSomeClusterV_a) / (T_METRIC( (lit_numObjetos-1) * lit_numObjetos / 2) - lit_pairsDiffClassUDiffClusterV_d ): measuare_undefRandIndex(T_METRIC);
  
}


template < typename T_METRIC, 
	   typename T_INSTANCES_CLUSTER_K 
	   >
T_METRIC
precision
(const ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
{
  T_INSTANCES_CLUSTER_K lit_pairsSomeClassUSomeClusterV_a =
    aimatchmatrix_confusion.getSomeClassUSomeClusterV();
	
  T_INSTANCES_CLUSTER_K lit_pairsDiffClassUSomeClusterV_c =
    aimatchmatrix_confusion.getDiffClassUSomeClusterV();
	
  return ((lit_pairsSomeClassUSomeClusterV_a + lit_pairsDiffClassUSomeClusterV_c) != 0)? 
    T_METRIC(lit_pairsSomeClassUSomeClusterV_a) /
    T_METRIC(lit_pairsSomeClassUSomeClusterV_a + lit_pairsDiffClassUSomeClusterV_c)
    :measuare_undefRandIndex(T_METRIC);
}



template < typename T_METRIC, 
	   typename T_INSTANCES_CLUSTER_K 
	   >
T_METRIC
recall
(const ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
{
  T_INSTANCES_CLUSTER_K lit_pairsSomeClassUSomeClusterV_a =
    aimatchmatrix_confusion.getSomeClassUSomeClusterV();
  T_INSTANCES_CLUSTER_K lit_pairsSameClassUDiffClusterV_b =
    aimatchmatrix_confusion.getSomeClassUDiffClusterV();
  
  return ((lit_pairsSomeClassUSomeClusterV_a + lit_pairsSameClassUDiffClusterV_b) != 0)? 
    T_METRIC(lit_pairsSomeClassUSomeClusterV_a) /
    T_METRIC(lit_pairsSomeClassUSomeClusterV_a + lit_pairsSameClassUDiffClusterV_b)
    :measuare_undefRandIndex(T_METRIC);
}


template < typename T_METRIC, 
	   typename T_INSTANCES_CLUSTER_K 
	   >
T_METRIC
fmeasure
(const ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
{
  T_INSTANCES_CLUSTER_K lit_pairsSomeClassUSomeClusterV_a =
    aimatchmatrix_confusion.getSomeClassUSomeClusterV();
  T_INSTANCES_CLUSTER_K lit_pairsDiffClassUSomeClusterV_c =
    aimatchmatrix_confusion.getDiffClassUSomeClusterV();
  T_INSTANCES_CLUSTER_K lit_pairsSameClassUDiffClusterV_b =
    aimatchmatrix_confusion.getSomeClassUDiffClusterV();
	
  return  ( (lit_pairsSomeClassUSomeClusterV_a + lit_pairsSameClassUDiffClusterV_b +
	     lit_pairsDiffClassUSomeClusterV_c) != 0 )? 
    (2. * T_METRIC(lit_pairsSomeClassUSomeClusterV_a)) /
    (T_METRIC
     (2*lit_pairsSomeClassUSomeClusterV_a+
      lit_pairsSameClassUDiffClusterV_b+ lit_pairsDiffClassUSomeClusterV_c))
    :measuare_undefRandIndex(T_METRIC);
}


    
/*! \fn T_METRIC purity(ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
    \brief purity
    \details
    \param
 */
template < typename T_METRIC,
	   typename T_INSTANCES_CLUSTER_K
	   >
T_METRIC
purity
(const ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>  &aimatchmatrix_confusion)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::purity";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\naimatchmatrix_confusion[" << &aimatchmatrix_confusion << ']' 
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  T_METRIC                lort_purity;
  T_INSTANCES_CLUSTER_K   lit_numObjetosMax;
  T_INSTANCES_CLUSTER_K   lit_numObjetos;
  uintidx                 list_numClassV;
  uintidx                 list_numClassU;

  list_numClassV  = aimatchmatrix_confusion.getNumRows()-1;
  list_numClassU  = aimatchmatrix_confusion.getNumColumns()-1;
  lit_numObjetos  = aimatchmatrix_confusion(list_numClassV,list_numClassU);

  //lort_purity = T_METRIC(0);
  lit_numObjetosMax = 0;
  for (uintidx li_j = 0; li_j < list_numClassU; li_j++) {
    T_INSTANCES_CLUSTER_K linstclusterK_max = aimatchmatrix_confusion(0,li_j);
    for (uintidx lui_i = 1; lui_i < list_numClassV; lui_i++) {
      if ( linstclusterK_max < aimatchmatrix_confusion(lui_i,li_j) )
	linstclusterK_max = aimatchmatrix_confusion(lui_i,li_j);
    }
    lit_numObjetosMax += linstclusterK_max;
  }

  lort_purity = (lit_numObjetos!=0)?T_METRIC(lit_numObjetosMax) / T_METRIC( lit_numObjetos )
    :measuare_undefPurity(T_METRIC);

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << "\tlit_numObjetosMax = " << lit_numObjetosMax 
	      << "\tlort_purity = " << lort_purity 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lort_purity;
}

} /*END namespace sm clustering measure*/

#endif /*UNSUPERVISED_MEASURE_HPP*/
