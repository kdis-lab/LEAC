/*! \file ga_clustering_operator.hpp
 *
 * \brief genetic operators oriented to clustering
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __GA_CLUSTERING_OPERATOR_HPP 
#define __GA_CLUSTERING_OPERATOR_HPP  

#include <climits>
#include <algorithm>  
#include <cmath> //std::round
#include <iterator>
#include "common.hpp"
#include "chromosome_cbga.hpp"
#include "chromosome_withrownull.hpp"
#include "chromosome_feac.hpp"
#include "chromosome_gga.hpp"
#include "chromosome_variablelength.hpp"
#include "linear_algebra_level1.hpp"
#include "linear_algebra_level2.hpp"
#include "bit_array.hpp"
#include "probability_selection.hpp"
#include "vector_utils.hpp"
#include "probability_distribution.hpp"
#include "ga_selection.hpp"
#include "ga_integer_operator.hpp"
#include "partition_label.hpp"
#include "probability_distribution.hpp"
#include "probability_selection.hpp"
#include "ga_function_objective.hpp"

#define  GENETICOPCLUSTER_MAX_GLA  INT_MAX

#include "verbose_global.hpp"

extern StdMT19937       gmt19937_eng;


/*! \namespace gaclusteringop
  \brief Genetic cluster-oriented operators
  \details Cluster-oriented operators mean the operators that are task-dependent, such as operators that copy, split, merge, and eliminate clusters of data objects, in contrast to conventional evolutionary operators that just exchange or switch bits without any regard to their task-dependent meaning \cite Hruschka:etal:GAclustering:survey:2009
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaclusteringop {
  
/*! \fn gaencode::ChromVariableLength<T_FEATURE,T_REAL>* newChromosome(mat::MatrixRow<T_FEATURE>& aimatrixrowt_minMaxKSegments, T_REAL airt_fitnessInitial, T_REAL airt_objetiveInitial)
  \brief Create new chromosome \cite He:Tan:GAclusteringVarK:TGCA:2012 
  \details Create new chromosome initial values of cluster centers are produced through uniformly random selection in ki segments by hierarchical agglomerative clustering algorithm.
  \param aimatrixrowt_minMaxKSegments a matrix with ki segments
  \param airt_fitnessInitial real number with value of the fitness initial
  \param airt_objetiveInitial real number with value of the function objective initial
*/
template <typename T_FEATURE, 
	   typename T_REAL  //TYPE OF OBJETIVE FUNCTION,
	   >
gaencode::ChromVariableLength<T_FEATURE,T_REAL>*
newChromosome
(mat::MatrixRow<T_FEATURE>& aimatrixrowt_minMaxKSegments,
 T_REAL                     airt_fitnessInitial,
 T_REAL                     airt_objetiveInitial
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::newChromosome";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( input aimatrixrowt_minMaxKSegments["  
	      << &aimatrixrowt_minMaxKSegments << "]\n"
	      << "  input airt_fitnessInitial = "  << airt_fitnessInitial
	      << "\n input airt_objetiveInitial = " << airt_objetiveInitial
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  gaencode::ChromVariableLength<T_FEATURE,T_REAL> *lochromvarlength_new =
    new gaencode::ChromVariableLength<T_FEATURE,T_REAL>
    ( aimatrixrowt_minMaxKSegments.getNumRows() *
      data::Instance<T_FEATURE>::getNumDimensions() );
	 
  mat::MatrixRow<T_FEATURE> 
    lmatrixrowt_centroidsChrom
    ( aimatrixrowt_minMaxKSegments.getNumRows(), 
      data::Instance<T_FEATURE>::getNumDimensions(),
      lochromvarlength_new->getString()
      );

  for (uintidx _lui_i =  0;
       _lui_i < lmatrixrowt_centroidsChrom.getNumRows();
       _lui_i ++) {

#if  DATATYPE_CENTROIDS_ROUND == 0
    std::uniform_real_distribution<T_FEATURE>
      uniformdis_kSegments
      (aimatrixrowt_minMaxKSegments(_lui_i,0),
       aimatrixrowt_minMaxKSegments(_lui_i,1)
       );
#else
    std::uniform_int_distribution<T_FEATURE>
      uniformdis_kSegments
      (aimatrixrowt_minMaxKSegments(_lui_i,0),
       aimatrixrowt_minMaxKSegments(_lui_i,1)
       );
#endif //DATATYPE_CENTROIDS_ROUND
	   
    for (uintidx _lui_j =  0; _lui_j < lmatrixrowt_centroidsChrom.getNumColumns(); _lui_j ++) {
      lmatrixrowt_centroidsChrom(_lui_i,_lui_j) = uniformdis_kSegments(gmt19937_eng);
    }
	  
  }

  lochromvarlength_new->setFitness(airt_fitnessInitial); 
  lochromvarlength_new->setObjetiveFunc(airt_objetiveInitial);

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochromvarlength_new->print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lochromvarlength_new;
  
}
    
  
/*CROSSOVER---------------------------------------------------------------------
 */

/* BEGIN crossPNNnew----------------------------------------------
 */

/*! \fn void merge(gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM,T_REAL> &aiochromcbga_offspring, gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM,T_REAL>  &aiochromcbga_part1, gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM,T_REAL> &aiochromcbga_part2, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist)
    \brief Merge two chomosomas 
    \details
    \param aiochromcbga_offspring
    \param aiochromcbga_part1
    \param aiochromcbga_part2 
    \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template <typename T_FEATURE, 
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename T_REAL,
	  typename INPUT_ITERATOR
	  >
void merge
(gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL>                             &aiochromcbga_offspring,
 gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL>                             &aiochromcbga_part1,
 gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL>                             &aiochromcbga_part2,
 INPUT_ITERATOR                      aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES

  const char* lpc_labelFunc = "gaclusteringop::merge";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labeoffspring;
    lostrstream_labeoffspring << lpc_labelFunc << ":aiochromcbga_offspring";
    aiochromcbga_offspring.print(std::cout,lostrstream_labeoffspring.str().c_str(),',',';');
    std::cout << '\n';
    std::ostringstream lostrstream_labelPart1;
    lostrstream_labelPart1 << lpc_labelFunc << ":aiochromcbga_part1";
    aiochromcbga_part1.print(std::cout,lostrstream_labelPart1.str().c_str(),',',';');
    std::cout << '\n';
    std::ostringstream lostrstream_labelPart2;
    lostrstream_labelPart2 << lpc_labelFunc << ":aiochromcbga_part2";
    aiochromcbga_part2.print(std::cout,lostrstream_labelPart2.str().c_str(),',',';');
    
    std::cout << "\ninput aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  ds::PartitionLinkedStats
    <T_FEATURE,
     T_CLUSTERIDX, //-1, 0, 1, .., K
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,
     T_FEATURE_SUM
     >  &lpartlinkstats_offspring  =  aiochromcbga_offspring.getPartition();

  ds::PartitionLinkedStats
    <T_FEATURE,
     T_CLUSTERIDX, //-1, 0, 1, .., K
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,
     T_FEATURE_SUM
     >  &aipartition_chromPart1  =  aiochromcbga_part1.getPartition();

  ds::PartitionLinkedStats
    <T_FEATURE,
     T_CLUSTERIDX, //-1, 0, 1, .., K
     T_INSTANCE_FREQUENCY,
     T_INSTANCES_CLUSTER_K,
     T_FEATURE_SUM
     >  &aipartition_chromPart2  =  aiochromcbga_part2.getPartition();
   
  mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
    &aiomatrixresizerow_offspring = aiochromcbga_offspring.getCodeBook();

  mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
    &aiomatrixresizerow_chromPart1 = aiochromcbga_part1.getCodeBook();

  mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
    &aiomatrixresizerow_chromPart2 = aiochromcbga_part2.getCodeBook();

  aiomatrixresizerow_offspring.merge
    (aiomatrixresizerow_chromPart1,
     aiomatrixresizerow_chromPart2
     );

  lpartlinkstats_offspring = aiochromcbga_part1.getPartition();
  
  lpartlinkstats_offspring.resize(aiomatrixresizerow_offspring.getNumRows());

  for (uintidx lidxinst_i = 0;
       aiiterator_instfirst != aiiterator_instlast;
       ++aiiterator_instfirst, lidxinst_i++)
    {
     
      T_CLUSTERIDX lIdxK_memberClusterPart1 =
	aipartition_chromPart1.getMemberShip(lidxinst_i);
	 
      T_CLUSTERIDX lIdxK_memberClusterPart2 =
	aipartition_chromPart2.getMemberShip(lidxinst_i);
      
      data::InstanceFreq
	<T_FEATURE,
	 T_INSTANCE_FREQUENCY> *literinstfo_iInstance  = 
	(data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
	*aiiterator_instfirst;
	 
      T_REAL lTnCentInst_distChrom1 = aifunc2p_dist
	(aiomatrixresizerow_chromPart1.getRow(lIdxK_memberClusterPart1),
	 literinstfo_iInstance->getFeatures(),
	 literinstfo_iInstance->getNumDimensions()
	 );
	    
      T_REAL lTnCentInst_distChrom2 = 
	aifunc2p_dist
	(aiomatrixresizerow_chromPart2.getRow(lIdxK_memberClusterPart2),
	 literinstfo_iInstance->getFeatures(),
	 literinstfo_iInstance->getNumDimensions()
	 );
	 
      if( lTnCentInst_distChrom2 < lTnCentInst_distChrom1 )  {
	
	lpartlinkstats_offspring.changeMemberShip
	  (aiomatrixresizerow_chromPart1.getNumRows() + lIdxK_memberClusterPart2, 
	   lidxinst_i,
	   literinstfo_iInstance->getFeatures(),
	   literinstfo_iInstance->getFrequency()
	   );
      }
    } 
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aiochromcbga_offspring.print(std::cout,lpc_labelFunc,',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

} /*END merge*/
  

/*! \fn void crossPNNnew (gaencode::ChromosomeCBGA <T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_FREQUENCY_SUM,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM, T_REAL> &aochromcbga_new, gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_FREQUENCY_SUM,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM,T_REAL> &aichromcbga_old1, gaencode::ChromosomeCBGA<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_FREQUENCY_SUM,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM, T_REAL> &aichromcbga_old2, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinstfo_instances, const T_CLUSTERIDX aicidx_numclusterKToReduce, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)

  \brief crossPNNnew \cite Franti:etal:GAclustering:gafranti:1997
  \details
  \param aochromcbga_new  OUT a Codebook encode by gaencode::ChromosomeCBGA
  \param aichromcbga_old1 IN  a Codebook encode by gaencode::ChromosomeCBGA
  \param aichromcbga_old2
  \param aivectorptinstfo_instances
  \param aicidx_numclusterKToReduce a numberof cluster to reduce
  \param aifunc2p_dist
 */
template <typename T_FEATURE,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,  //0, 1, .., N
	  typename T_FEATURE_SUM,
	  typename T_CLUSTERIDX, //-1, 0, 1, .., K
	  typename T_REAL,
	  typename INPUT_ITERATOR
	  >
void
crossPNNnew
(gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM, 
 T_REAL>                            &aochromcbga_new,
 gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM,
 T_REAL>                            &aichromcbga_old1,
 gaencode::ChromosomeCBGA
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM, 
 T_REAL>                            &aichromcbga_old2,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const T_CLUSTERIDX                 aicidx_numclusterKToReduce,
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "gaclusteringop::crossPNNnew:  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t(output gaencode::ChromosomeCBGA: aochromcbga_new["
	      << &aochromcbga_new << "]\n"
	      << "\t input  gaencode::ChromosomeCBGA: aichromcbga_old1["
	      << &aichromcbga_old1 << "]\n"
	      << "\t input  gaencode::ChromosomeCBGA: aichromcbga_old2["
	      << &aichromcbga_old2 << "]\n"
	      << " input  aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "\t input  T_CLUSTERIDX:     aicidx_numclusterKToReduce " 
	      << aicidx_numclusterKToReduce << '\n'
	      << "\t input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      	      << &aifunc2p_dist << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  /* Combine chromosomes: 
         aochromcbga_new <-- aichromcbga_old1 + aichromcbga_old2, 
   */
  merge
    (aochromcbga_new,
     aichromcbga_old1,
     aichromcbga_old2,
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  /* Calculate optimal centroids (2K) for combined partitioning (mean) 
   */
  aochromcbga_new.optimalCodebook(); 
  
  /* Reduce clustering size from 2K -> K using PNN algorithm 
   */
  clusteringop::pnnFast
    (aochromcbga_new.getPartition(),
     aochromcbga_new.getCodeBook(),
     aiiterator_instfirst,
     aiiterator_instlast,
     aicidx_numclusterKToReduce,
     aifunc2p_dist
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "gaclusteringop::crossPNNnew: OUT"
	      << '(' << geiinparam_verbose << ")\n";
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*END crossPNNnew*/


/*! \fn void onePointIndivisibleCrossover(gaencode::ChromFixedLength<T_GENE,T_METRIC> &aochrom_child1, gaencode::ChromFixedLength<T_GENE,T_METRIC> &aochrom_child2, const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent1, const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent2) 
  \brief single-point crossover (Michalewicz, 1992).
  \details
  \param aochrom_child1 a gaencode::ChromFixedLength<T_GENE,T_METRIC>
  \param aochrom_child2 a gaencode::ChromFixedLength<T_GENE,T_METRIC>
  \param aichrom_parent1 a gaencode::ChromFixedLength<T_GENE,T_METRIC>
  \param aichrom_parent2 a ChromFixedLength<T_GENE,T_METRIC> 
 */
template < typename T_GENE,
	   typename T_METRIC
	   >
void
onePointIndivisibleCrossover
(gaencode::ChromVariableLength<T_GENE,T_METRIC>       &aochrom_child1,
 gaencode::ChromVariableLength<T_GENE,T_METRIC>       &aochrom_child2, 
 const gaencode::ChromVariableLength<T_GENE,T_METRIC> &aichrom_parent1, 
 const gaencode::ChromVariableLength<T_GENE,T_METRIC> &aichrom_parent2
 )
{
  const long liu_k1 =
    long((aichrom_parent1.getStringSize() / data::Instance<T_GENE>::getNumDimensions()));
  const long liu_k2 =
    long((aichrom_parent2.getStringSize() / data::Instance<T_GENE>::getNumDimensions()));
  long lui_c1 = 0;
  if ( liu_k1 >= 2 ) {
    std::uniform_int_distribution<long> uniformdis_uiRandC1(0,liu_k1-1);
    // C1 = rand() mod K1
    lui_c1 = uniformdis_uiRandC1(gmt19937_eng); 
  }
  const long lui_lbc2 = std::min<long>(2,std::max<long>(0,2-(liu_k1-lui_c1)));
  const long lui_ubc2 = liu_k2 - std::max<long>(0,2-lui_c1);
  long lui_c2 = 0;
  if ( liu_k2 >= 2 ) {
    if ( lui_lbc2 < (lui_ubc2-lui_lbc2-1) ) {
      //C2 = LB(C2) + rand()mod (UB(C2) - LB(C2))
      std::uniform_int_distribution<long> uniformdis_uiRandC2(lui_lbc2,lui_ubc2-lui_lbc2-1);
      lui_c2 =  uniformdis_uiRandC2(gmt19937_eng);
    }
    else {
      lui_c2 = lui_lbc2;
    }
  }
 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::onePointIndivisibleCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(\noutput gaencode::ChromVariableLength<T_GENE>: aochrom_child1[" 
	      << &aochrom_child1 << "]\n"
	      << "output gaencode::ChromVariableLength<T_GENE>: aochrom_child2[" 
	      << &aochrom_child2 << "]\n";
    std::ostringstream lostrstream_labelparent1;
    lostrstream_labelparent1 << lpc_labelFunc << ":input parent1";
    aichrom_parent1.print(std::cout,lostrstream_labelparent1.str().c_str());
    std::cout << '\n';
    
    std::ostringstream lostrstream_labelparent2;
    lostrstream_labelparent2 << lpc_labelFunc << ":input parent2";
    aichrom_parent2.print(std::cout,lostrstream_labelparent2.str().c_str());
    std::cout << '\n';
    std::cout
      << "liu_k1 = " <<  liu_k1 << '\n'
      << "liu_k2 = " <<  liu_k2 << '\n'
      << "lui_lbc2 = " << lui_lbc2 << '\n'
      << "lui_ubc2 = " << lui_ubc2 << '\n'
      << "lui_c1 = " << lui_c1 << '\n'
      << "lui_c2 = " << lui_c2 << std::endl;
  }
#endif /*__VERBOSE_YES*/
   
  aochrom_child1.resize( uintidx((lui_c1 + liu_k2 - lui_c2) * data::Instance<T_GENE>::getNumDimensions()) );
  aochrom_child2.resize( uintidx((lui_c2 + liu_k1 - lui_c1) * data::Instance<T_GENE>::getNumDimensions()) );
 
  interfacesse::copy
    (aochrom_child1.getString(), 
     aichrom_parent1.getString(), 
     uintidx(lui_c1 * data::Instance<T_GENE>::getNumDimensions())
     );

  interfacesse::copy
    (aochrom_child1.getString() + lui_c1 * data::Instance<T_GENE>::getNumDimensions(), 
     aichrom_parent2.getString() + lui_c2 * data::Instance<T_GENE>::getNumDimensions(),
     uintidx((liu_k2 - lui_c2) * data::Instance<T_GENE>::getNumDimensions())
     );

  interfacesse::copy
    (aochrom_child2.getString(), 
     aichrom_parent2.getString(), 
     uintidx(lui_c2 * data::Instance<T_GENE>::getNumDimensions())
     );

  interfacesse::copy
    (aochrom_child2.getString() + lui_c2 * data::Instance<T_GENE>::getNumDimensions(), 
     aichrom_parent1.getString() + lui_c1 * data::Instance<T_GENE>::getNumDimensions(),
     uintidx((liu_k1 - lui_c1) * data::Instance<T_GENE>::getNumDimensions())
     );
 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aochrom_child1.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aochrom_child2.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*onePointIndivisibleCrossover*/


/*! \fn void onePointCrossover(gaencode::Chromosomemat::MatrixWithRowNull<T_GENE,T_METRIC> &aiochrommwrn_child1, gaencode::ChromosomeMatrixWithRowNull <T_GENE,T_METRIC> &aiochrommwrn_child2, const uintidx aiui_positionGene)
  \brief One point crossover \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002
  \details This operator uses two chromosomes that are combined as described in the example:

   \code
  Eeach cluster centre is considered to be an indivisible gene, is explained below with an example.

  Let two strings

  #   (20.4, 13.2) # # (15.8, 2.9)| # (10.0, 5.0) (22.7, 17.7) #    #
  (13.2, 15.6)   # # # (5.3, 13.7)| # (10.5, 16.2) (7.9, 15.3) # (18.3, 14.5)

  If the crossover position be 5 as shown above. Then the offspring are

  # (20.4, 13.2) # # (15.8, 2.9)  | # (10.5, 16.2) (7.9, 15.3) # (18.3, 14.5) #
  (13.2, 15.6)   # # # (5.3, 13.7)| # (10.0, 5.0) (22.7, 17.7) #  #
  \endcode

  \param aiochrommwrn_child1 a gaencode::ChromosomeMatrixWithRowNull offspring chromosome 
  \param aiochrommwrn_child2 a gaencode::ChromosomeMatrixWithRowNull offspring chromosome 
*/
template <typename T_GENE,
	   typename T_METRIC
	   >
inline
void
onePointCrossover
(gaencode::ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &aiochrommwrn_child1,
 gaencode::ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &aiochrommwrn_child2
 )
{
 uintidx lui_minStringLength =
    std::min(aiochrommwrn_child1.getNumRowsMax(),aiochrommwrn_child1.getNumRowsMax());
 std::uniform_int_distribution<uintidx> luniformdis_uiCrossover1KMax
   (1,lui_minStringLength-1);
 uintidx lui_positionGene = luniformdis_uiCrossover1KMax(gmt19937_eng);
 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::onePointCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
  	      << "\n(output gaencode::Chromosomemat::MatrixWithRowNull<>&: aiochrommwrn_child1[" 
	      << &aiochrommwrn_child1 << "]\n"
	      << " output gaencode::Chromosomemat::MatrixWithRowNull<>&: aiochrommwrn_child2[" 
	      << &aiochrommwrn_child2 << "]\n"
	      << " lui_positionGene =  " << lui_positionGene
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  if ( lui_minStringLength  == aiochrommwrn_child1.getNumRowsMax() ) {
    aiochrommwrn_child1.combination(lui_positionGene,aiochrommwrn_child2);
  }
  else {
    aiochrommwrn_child2.combination(lui_positionGene,aiochrommwrn_child1);
  }
  
  aiochrommwrn_child1.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aiochrommwrn_child1.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());
  aiochrommwrn_child2.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aiochrommwrn_child2.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aiochrommwrn_child1.print(std::cout,lpc_labelFunc,',',';');
    std::cout << '\n';
    aiochrommwrn_child2.print(std::cout,lpc_labelFunc,',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*! \fn gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC> mergeCrossover(gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>& aichrom_parent1, gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>& aichrom_parent2)
  \brief mergeCrossover \cite Agustin:etal:GAclusteringVarK:GGA:2012
  \details
  \param aichrom_parent1 a gaencode::ChromosomeGGA
  \param aichrom_parent2 a gaencode::ChromosomeGGA
 */
template <typename T_CLUSTERIDX,
	   typename T_METRIC
	   >
gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>
mergeCrossover
(gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>& aichrom_parent1, 
 gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>& aichrom_parent2
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::mergeCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "PARENT1:"
	      << aichrom_parent1
	      << "\nPARENT2:"
	      << aichrom_parent2
	      << std::endl;
  }
#endif //__VERBOSE_YES
 
  T_CLUSTERIDX lcidx_kParent1 = aichrom_parent1.getNumClusterK();
  T_CLUSTERIDX lcidx_kParent2 = aichrom_parent2.getNumClusterK();
  
  std::pair<T_CLUSTERIDX,T_CLUSTERIDX> lpaircidx_parent1(0,0);
  std::pair<T_CLUSTERIDX,T_CLUSTERIDX> lpaircidx_parent2(0,0);

  if ( lcidx_kParent1 > 1) { 
    std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_KParent1
      (0,lcidx_kParent1-1);
    lpaircidx_parent1 = 
      prob::getRandPairUnlikeInOrd
      ([&]() -> T_CLUSTERIDX
       {
	 return uniformdis_cidx_KParent1(gmt19937_eng);
       }
       );
  }
  

  if ( lcidx_kParent2 > 1) { 
    std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_KParent2
      (0,lcidx_kParent2-1);
    lpaircidx_parent2 = 
      prob::getRandPairUnlikeInOrd
      ([&]() -> T_CLUSTERIDX
       {
	 return uniformdis_cidx_KParent2(gmt19937_eng);
       }
       );
  }
    
   
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "lcidx_kParent1 = "  << lcidx_kParent1
	      << ",LimRandom[" << lpaircidx_parent1.first << ',' << lpaircidx_parent1.second << ']'
	      << ",lcidx_kParent2 = " << lcidx_kParent2
	      << ",LimRandom[" << lpaircidx_parent2.first << ',' << lpaircidx_parent2.second << ']'
	      << std::endl;
  }
   --geiinparam_verbose;
#endif //__VERBOSE_YES

 const T_CLUSTERIDX lcidx_kOffspring = 
    lpaircidx_parent1.second - lpaircidx_parent1.first + lpaircidx_parent2.second - lpaircidx_parent2.first + 2;

 std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_KOffspring
    (0,lcidx_kOffspring-1);
 
  gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC> lochrom_offspring((uintidx) lcidx_kOffspring );
  std::vector<uintidx> lvectort_numInstClusterK((uintidx) lcidx_kOffspring, 0);

  lochrom_offspring.initializeGroupSec();

  //INITIALIZE
  interfacesse::copya
    (lochrom_offspring.getString(),
     T_CLUSTERIDX(UNKNOWN_CLUSTER_IDX),
     lochrom_offspring.getStringSize()
     );

  T_CLUSTERIDX *larraycidx_parent    = aichrom_parent1.getString();
  T_CLUSTERIDX *larraycidx_offspring = lochrom_offspring.getString();
  for (uintidx lui_i = 0;
       lui_i < gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>::getElementSize();
       lui_i++)
    {
      if ( lpaircidx_parent1.first <= *larraycidx_parent && *larraycidx_parent <= lpaircidx_parent1.second) {
	*larraycidx_offspring = *larraycidx_parent - lpaircidx_parent1.first;
	++lvectort_numInstClusterK.at(*larraycidx_offspring);
      }
      ++larraycidx_parent;
      ++larraycidx_offspring;
    }

  T_CLUSTERIDX lcidx_beginKparent2 = lpaircidx_parent1.second - lpaircidx_parent1.first + 1;
  larraycidx_parent    = aichrom_parent2.getString();
  larraycidx_offspring = lochrom_offspring.getString();
  for (uintidx lui_i = 0;
       lui_i < gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>::getElementSize();
       lui_i++)
    {
    if ( (*larraycidx_offspring == UNKNOWN_CLUSTER_IDX) &&
	 (lpaircidx_parent2.first <= *larraycidx_parent && *larraycidx_parent <= lpaircidx_parent2.second)
	 )
      {
	*larraycidx_offspring = *larraycidx_parent - lpaircidx_parent2.first + lcidx_beginKparent2;
	++lvectort_numInstClusterK.at(*larraycidx_offspring);
	
      }
    ++larraycidx_parent;
    ++larraycidx_offspring;
  }
  
  larraycidx_offspring = lochrom_offspring.getString();
  for (uintidx lui_i = 0;
       lui_i < gaencode::ChromosomeGGA<T_CLUSTERIDX,T_METRIC>::getElementSize();
       lui_i++)
    {
    if ( (*larraycidx_offspring == UNKNOWN_CLUSTER_IDX) ) {
      T_CLUSTERIDX lcidx_randK =
	uniformdis_cidx_KOffspring(gmt19937_eng);
      *larraycidx_offspring = lcidx_randK;
      ++lvectort_numInstClusterK.at(*larraycidx_offspring);
    }
    ++larraycidx_offspring;
  }

  std::vector<uintidx> lvectorcidx_clustersKeep;

  lvectorcidx_clustersKeep.reserve(lvectort_numInstClusterK.size());
  for ( uintidx lui_i = 0; lui_i < lvectort_numInstClusterK.size(); lui_i++) {
    if ( lvectort_numInstClusterK[lui_i] != 0 ) {
      lvectorcidx_clustersKeep.push_back(lui_i);
    }
  }
  if (lvectort_numInstClusterK.size() !=  lvectorcidx_clustersKeep.size() ) {
    lochrom_offspring.decrementGroupSecSize
      (lvectort_numInstClusterK.size() - lvectorcidx_clustersKeep.size() );
    gaintegerop::labelKeep
      <T_CLUSTERIDX,T_METRIC>
      (lochrom_offspring,
       lvectorcidx_clustersKeep
       );

  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lochrom_offspring.print();
    std::cout << std::endl;
	
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lochrom_offspring;
} 


/*! \fn void crossoverCGA (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> &aochrom_childC, gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> &aichrom_parentA, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist )
  \brief crossoverCGAAux \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003 
  \details
  \param aochrom_childC
  \param aichrom_parentA
  \param aivectorptinst_instances
  \param aifunc2p_dist
 */
template <typename T_CLUSTERIDX,
	  typename T_REAL,
	  typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename INPUT_ITERATOR
	  >
void
crossoverCGA
(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>   &aochrom_childC,
 gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>   &aichrom_parentA,
 gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>   &aichrom_parentB,
 mat::MatrixRow<T_FEATURE>                         &aomatrixrowt_centroids,
 mat::MatrixRow<T_FEATURE_SUM>                     &aomatrixrowt_sumInstCluster,
 std::vector<T_INSTANCES_CLUSTER_K>                &aovectort_numInstClusterK,
 INPUT_ITERATOR                                    aiiterator_instfirst,
 const INPUT_ITERATOR                              aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>                &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::crossoverCGA";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output gaencode::ChromFixedLength<>: aochrom_childC[" 
	      << &aochrom_childC << "]\n";
    std::ostringstream lostrstream_labelchildA;
    lostrstream_labelchildA << "childA:" << lpc_labelFunc;
    aichrom_parentA.print(std::cout,lostrstream_labelchildA.str().c_str());
    std::cout  << "\ninput  gaencode::ChromFixedLength<>: aichrom_parentB[" 
	      << &aichrom_parentB << "]\n"
      	      << "input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      <<  "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  aochrom_childC = aichrom_parentB;
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "CROSSOVERCGA 1. RANDOMLY CHOOSES c in {1,...,k_1}:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << std::endl; 
  }
#endif /*__VERBOSE_YES*/
  
  std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_1K1
    (1,aichrom_parentA.getGene(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1));

  std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_0K1
    (0,aichrom_parentA.getGene(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)-1);

  T_CLUSTERIDX lcidx_c = uniformdis_cidx_1K1(gmt19937_eng);
  
  std::unordered_set<T_CLUSTERIDX>&& lunorderedset_cidx_0K1 =
    prob::getWithoutRepeatsSet
    (lcidx_c,
     [&]() -> T_CLUSTERIDX
     {
       return uniformdis_cidx_0K1(gmt19937_eng);
     }
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << "CROSSOVERCGA 1. RANDOMLY CHOOSES c in {1,...,k_1}: OUT"
      << '(' << geiinparam_verbose << ')'
      << "\nc = " << lcidx_c << " in {1,...,"
      << aichrom_parentA.getGene(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
      << '}'
      << std::endl; 
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  partition::PartitionLabel
    <T_CLUSTERIDX>
    lpartition_clustersParent1
    (aichrom_parentA.getString(),
     gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
     aichrom_parentA.getGene
     (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
     );
  
  ds::PartitionLinked<T_CLUSTERIDX>
    lpartlink_memberShipParent1 =
    ds::getPartitionlinked
    (lpartition_clustersParent1);

  ds::IteratorPartitionLinked<T_CLUSTERIDX>
    literpartlink_parent1(&lpartlink_memberShipParent1);

  std::unordered_set<T_CLUSTERIDX> lounorderedset_affected;
  lounorderedset_affected.reserve(lunorderedset_cidx_0K1.size());
  
  for (const auto lcidx_gene: lunorderedset_cidx_0K1) {
    for ( literpartlink_parent1.begin(lcidx_gene);
	  literpartlink_parent1.end();
	  literpartlink_parent1.next()
	  )
      {
	lounorderedset_affected.insert
	  ( aochrom_childC.getGene(literpartlink_parent1.getValue()));
      }  
  }


#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    inout::containerprint
    (lounorderedset_affected.begin(),
     lounorderedset_affected.end(),
     std::cout,
     "CROSSOVERCGA 2. GENES INDIRECTLY AFFECTED",
     ','
     );
  std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  for ( auto liter_gene = aochrom_childC.begin(); liter_gene != aochrom_childC.end()-1; liter_gene++)
    {
      auto  liter_setaffected = lounorderedset_affected.find(*liter_gene);
      if ( liter_setaffected != lounorderedset_affected.end() )
	*liter_gene = T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN);
    }

  for (const auto lcidx_gene: lunorderedset_cidx_0K1) {
    for ( literpartlink_parent1.begin(lcidx_gene);
	  literpartlink_parent1.end();
	  literpartlink_parent1.next()
	  )
      {
	aochrom_childC.setGene(literpartlink_parent1.getValue(),lcidx_gene);
      }  
  }


#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
   aochrom_childC.print(std::cout,"CROSSOVERCGA 3. AFFECTED NOW CHANGED TO 0");
   std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  //renumber
  T_CLUSTERIDX lcidx_consecK = 0;
    
  std::map<T_CLUSTERIDX,T_CLUSTERIDX> mymap;
  for ( auto liter_gene = aochrom_childC.begin(); liter_gene != aochrom_childC.end()-1; liter_gene++)
    {
      if ( *liter_gene != T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN) ) {
	auto liter_map = mymap.find(*liter_gene);
	if ( liter_map == mymap.end() )
	  mymap.insert
	    (std::pair<T_CLUSTERIDX,T_CLUSTERIDX>(*liter_gene,lcidx_consecK++) );
      }

    }
  
  for ( auto liter_gene = aochrom_childC.begin(); liter_gene != aochrom_childC.end()-1; liter_gene++)
    {
      if ( *liter_gene != T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN) ) {
	*liter_gene = mymap[*liter_gene];
      }
    }

  aochrom_childC.setGene
    (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
     (T_CLUSTERIDX) mymap.size()
     );
   
  partition::PartitionLabel
    <T_CLUSTERIDX>
    lpartition_clusters
    (aochrom_childC.getString(),
     gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
     aochrom_childC.getGene
     (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
     );

  mat::MatrixRow<T_FEATURE> 
    lmatrixrowt_centroids
    ((uintidx) mymap.size(), 
     data::Instance<T_FEATURE>::getNumDimensions() 
     );
  
  mat::MatrixRow<T_FEATURE_SUM>       
    lmatrixrowt_sumInstCluster
    ((uintidx) mymap.size(),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
	
  std::vector<T_INSTANCES_CLUSTER_K> 
    lvectort_numInstClusterK((uintidx) mymap.size());

  clusteringop::getCentroids
    (lmatrixrowt_centroids,
     lmatrixrowt_sumInstCluster,
     lvectort_numInstClusterK,
     lpartition_clusters,
     aiiterator_instfirst,
     aiiterator_instlast
     );
  
  clusteringop::setUpCuster
    (aochrom_childC.getString(),
     lmatrixrowt_centroids,
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  std::swap(aomatrixrowt_centroids,lmatrixrowt_centroids);
  std::swap(aomatrixrowt_sumInstCluster,lmatrixrowt_sumInstCluster);
  std::swap(aovectort_numInstClusterK,lvectort_numInstClusterK);
    	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')';
    std::ostringstream lostrstream_labelchildC;
    lostrstream_labelchildC << "childC:" << lpc_labelFunc;
    aochrom_childC.print(std::cout,lostrstream_labelchildC.str().c_str());
    std::cout << std::endl;
    std::cout << std::endl; 
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
} /*crossoverCGA*/

  
/*MUTATION---------------------------------------------------------------------
*/

  
/*! \fn void splittingMutation(gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>  &aiochromab_toMutate)
  \brief  Mutation by cluster splitting \cite Agustin:etal:GAclusteringVarK:GGA:2012
  \details Mutation by clusters merging it consists of merging two existing clusters, randomly selected, into just one.
  \param aiochromab_toMutate a gaencode::ChromosomeGGA a mutate
 */
template <typename T_CLUSTERIDX,
	   typename T_REAL
	   >
void
splittingMutation
(gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL>  &aiochromab_toMutate)
{ 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::splittingMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    aiochromab_toMutate.print(std::cout,lpc_labelFunc);
    std::cout << "\n)"
	      << std::endl;
   }
#endif //__VERBOSE_YES

  T_CLUSTERIDX lcidx_ktoMutate = aiochromab_toMutate.getNumClusterK();
  std::vector<uintidx> lvectorrt_numInstClusterK(lcidx_ktoMutate,0);
  T_CLUSTERIDX *larraycidx_toMutate = aiochromab_toMutate.getString();
  for (uintidx lui_i = 0; lui_i < aiochromab_toMutate.getElementSize(); lui_i++) {
    ++lvectorrt_numInstClusterK.at(*larraycidx_toMutate);
    ++larraycidx_toMutate;
  }

  const std::vector<T_REAL>&& lvectorprob_selectBySizeClusterK =
    prob::makeDistRouletteWheel
    (lvectorrt_numInstClusterK.begin(),lvectorrt_numInstClusterK.end(),
     [](const uintidx& lui_numInstClusterK) -> T_REAL
     {
       return T_REAL(lui_numInstClusterK);
     }
     );
  T_CLUSTERIDX lcidx_selectClusterK =
    gaselect::getIdxRouletteWheel<T_REAL>
     (lvectorprob_selectBySizeClusterK,
      T_CLUSTERIDX(0)
     );

  if ( lvectorrt_numInstClusterK[lcidx_selectClusterK] >= 2 ) {

    larraycidx_toMutate = aiochromab_toMutate.getString();
    bool lb_mutOpApp = false;
    for (uintidx lui_i = 0; lui_i < aiochromab_toMutate.getElementSize(); lui_i++) {
      if ( lb_mutOpApp &&  *larraycidx_toMutate == lcidx_selectClusterK ) {
	*larraycidx_toMutate = lcidx_ktoMutate;
	lb_mutOpApp = false;
      }
      else if ( *larraycidx_toMutate == lcidx_selectClusterK ) {
	lb_mutOpApp = true;
      }
      ++larraycidx_toMutate;
    }

    aiochromab_toMutate.increaseGroupSecSize(); //Add cluster
    
  }
  
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "lcidx_selectClusterK = " << lcidx_selectClusterK
	      << "\tlui_numInstClusterK = " << lvectorrt_numInstClusterK[lcidx_selectClusterK] << '\n';
    aiochromab_toMutate.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
   }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

}  /*gabinaryop::splittingMutation*/

  
/*! \fn void mergeMutation(gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> &aiochromab_toMutate, FUNCTION_RAND function_rand01) 
  \brief Mutation by clusters merging \cite Agustin:etal:GAclusteringVarK:GGA:2012
  \details Mutation by clusters merging it consists of merging two existing clusters, randomly selected, into just one.
  \param aiochromab_toMutate a gaencode::ChromosomeGGA to mutate
 */
template <typename T_CLUSTERIDX,
	   typename T_REAL
	   >
void
mergeMutation
(gaencode::ChromosomeGGA<T_CLUSTERIDX,T_REAL> &aiochromab_toMutate)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::mergeMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    aiochromab_toMutate.print(std::cout,lpc_labelFunc);
    std::cout << "\n)"
	      << std::endl;
   }
#endif //__VERBOSE_YES


  T_CLUSTERIDX lcidx_ktoMutate = aiochromab_toMutate.getNumClusterK();
  std::pair<T_CLUSTERIDX,T_CLUSTERIDX> lpaircidx_mergeK
    (UNKNOWN_CLUSTER_IDX,UNKNOWN_CLUSTER_IDX);
  
  if ( lcidx_ktoMutate > 2 ) {
    std::vector<uintidx> lvectorrt_numInstClusterK(lcidx_ktoMutate,0);
    T_CLUSTERIDX *larraycidx_toMutate = aiochromab_toMutate.getString();
    for (uintidx lui_i = 0; lui_i < aiochromab_toMutate.getElementSize(); lui_i++) {
      ++lvectorrt_numInstClusterK.at(*larraycidx_toMutate);
      ++larraycidx_toMutate;
    }

    const std::vector<T_REAL>&& lvectorprob_selectBySizeClusterK =
      prob::makeDistRouletteWheel
      (lvectorrt_numInstClusterK.begin(),lvectorrt_numInstClusterK.end(),
       [](const uintidx& lui_numInstClusterK) -> T_REAL
       {
	 return T_REAL(1.0) / T_REAL(lui_numInstClusterK);
       }
       );

    std::pair<uintidx,uintidx> lpair_mergeIdxK =   
      prob::getIdxPairUnlikeRoulette
      (lvectorprob_selectBySizeClusterK, 
       uintidx(0)
       );
     
    if (lpair_mergeIdxK.second < lpair_mergeIdxK.first) {
	std::swap(lpair_mergeIdxK.second,lpair_mergeIdxK.first);
    }

    lpaircidx_mergeK.first  = (T_CLUSTERIDX) lpair_mergeIdxK.first;
    lpaircidx_mergeK.second = (T_CLUSTERIDX) lpair_mergeIdxK.second;
    
    T_CLUSTERIDX *larraycidx_chromMerge  = aiochromab_toMutate.getString();
    for (uintidx lui_i = 0; lui_i < aiochromab_toMutate.getElementSize(); lui_i++) {
      if ( *larraycidx_chromMerge == lpaircidx_mergeK.second ) {
	*larraycidx_chromMerge = lpaircidx_mergeK.first;
      }
      else if ( *larraycidx_chromMerge > lpaircidx_mergeK.second ) {
	--(*larraycidx_chromMerge);
      }
      ++larraycidx_chromMerge;
    }

    aiochromab_toMutate.decrementGroupSecSize(1);
    
  }

 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "pair_mergeIdxK: " << lpaircidx_mergeK.first
	      << " and " << lpaircidx_mergeK.second
	      << '\n';
    aiochromab_toMutate.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
   }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
} //mutationMergingAgustin2012


/*! \fn void void randomMutation(gaencode::ChromosomeMatrixWithRowNull<T_FEATURE,T_REAL>& aiochromatrixwrn_toMutate)
    \brief Random mutation \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002
    \details Each valid position (i.e., which is not ‘#’) in a chromosome is mutated. A number \f$\delta\f$ in the range [0,1] is generated with uniform distribution. If the value at that position is \f$g_i\f$, then after mutation it becomes

    \f[ 
    \cases{
    g_i \cdot (1 \pm 2 \sigma ), &\quad $g_i \neq 0$ \cr 
    \pm 2 \sigma,                &\quad $g_i = 0$. \cr
    }
    \f]

    The ‘+’ or ‘−’ sign occurs with equal probability.

    \param aiochromatrixwrn_toMutate a gaencode::ChromosomeMatrixWithRowNull chromosome to mutate 
*/  
template <typename T_FEATURE,
	  typename T_REAL
	  > 
void 
randomMutation
(gaencode::ChromosomeMatrixWithRowNull
 <T_FEATURE,T_REAL>&                    aiochromatrixwrn_toMutate)
{

  static std::uniform_real_distribution<T_REAL> lsuniformdis_real01(0.0,1.0);
  const T_REAL lrt_sign2delta = 2.0 * ((lsuniformdis_real01(gmt19937_eng) < 0.5)?1.0:-1.0) * lsuniformdis_real01(gmt19937_eng);

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::randomMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
  	      << "(output gaencode::Chromosomemat::MatrixWithRowNull<>&: aiochromatrixwrn_toMutate["
	      << &aiochromatrixwrn_toMutate << "]\n"
	      << "lrt_sign2delta =  " <<  lrt_sign2delta
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  for ( uintidx lst_i = 0; lst_i < aiochromatrixwrn_toMutate.getNumRowsMax(); ++lst_i) {

    T_FEATURE *larrayt_string  = aiochromatrixwrn_toMutate.getRow(lst_i);
 
    if ( larrayt_string != NULL) {

      for ( uintidx lst_j = 0; lst_j < aiochromatrixwrn_toMutate.getNumColumns(); ++lst_j) {
	    
#if  DATATYPE_CENTROIDS_ROUND == 0
 
	larrayt_string[lst_j] = ( larrayt_string[lst_j] != 0.0)?
	  larrayt_string[lst_j]  * ( 1.0  + lrt_sign2delta )
	  :lrt_sign2delta;

#else

	larrayt_string[lst_j] = ( larrayt_string[lst_j] != 0)?
	  T_FEATURE( std::round((DATATYPE_REAL) larrayt_string[lst_j]  * ( 1.0  +  lrt_sign2delta)))
	  :T_FEATURE(std::round((DATATYPE_REAL) lrt_sign2delta));

#endif //DATATYPE_CENTROIDS_ROUND
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
  	      << aiochromatrixwrn_toMutate
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
}

  
/*! \fn void biDirectionHMutation(gaencode::ChromosomeString<T_GENE,T_METRIC> &aiochromstr_toMutate, const T_METRIC airt_minObjetiveFunc, const T_METRIC airt_maxObjetiveFunc, const T_GENE* aiarrayt_minFeatures, const T_GENE* aiarrayt_maxFeatures)
  \brief biDirectionHMutation \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002
  \details The mutation does the following
\f[
\hbox{mutate}(g_{jl}) = \cases{
    g_{jl} + \delta \times ( \hbox{max}(x_l) - g_{jl})  &\quad  if $ \delta \geq 0 $, for $j= 1,2,..,k$ and $l= 1,2,...d$, \cr
    g_{jl} + \delta \times (g_{jl} - \hbox{min}(x_l))   &\quad  if $ \delta < 0 $. \cr
    }
\f]       
Where \f$\delta\f$ is a random number in the interval \f$[-R,+R]\f$:
\f[
R = \cases{
    M - M_{min} \over M_{max} - M_{min} &\quad  if $ M_{max} > M $, \cr
    1                                   &\quad  if $ M_{min} = M_{max} $. \cr
    }
\f]
\f$M_{min}\f$ and \f$M_{max}\f$ be the minimum and maximum values of the clustering metric, respectively, in the current population. \f$M\f$ is the clustering metric value of the current chromosome that must be mutated.
    \param aiochromstr_toMutate a gaencode::ChromosomeString chromosome to mutate
    \param airt_minObjetiveFunc a real number with the minimum value of the objective function of the current population
    \param airt_maxObjetiveFunc a real number with the maximum value of the objective function of the current population
    \param aiarrayt_minFeatures an array with the minimum values of the dataset instances
    \param aiarrayt_maxFeatures an array with the maximum values of the dataset instances

  \note The authors of this operator do not assign a specific name for this, here it is called gaclusteringop::biDirectionHMutation
 */
template <typename T_GENE, 
	  typename T_METRIC
	 >
void 
biDirectionHMutation
(gaencode::ChromosomeString<T_GENE,T_METRIC>  &aiochromstr_toMutate,
 const T_METRIC                               airt_minObjetiveFunc,
 const T_METRIC                               airt_maxObjetiveFunc,
 const T_GENE*                                aiarrayt_minFeatures,
 const T_GENE*                                aiarrayt_maxFeatures
 )
{ //BEGIN biDirectionHMutation

  const uintidx lconstui_numClusterK =
    aiochromstr_toMutate.getStringSize() / data::Instance<T_GENE>::getNumDimensions();
    
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::biDirectionHMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output gaencode::Chromosome: &aiochromstr_toMutate[" 
	      << &aiochromstr_toMutate << "]\n";
    std::ostringstream lostrstream_labelChromCentroids;
    lostrstream_labelChromCentroids << lpc_labelFunc << ":inout aiochromstr_toMutate";
    aiochromstr_toMutate.print(std::cout,lostrstream_labelChromCentroids.str().c_str());
    std::cout << '\n';
    std::cout << "input  M = " << aiochromstr_toMutate.getObjetiveFunc()
	      << "\tairt_minObjetiveFunc = " << airt_minObjetiveFunc
              << "\tairt_maxObjetiveFunc = " << airt_maxObjetiveFunc
	      << "\tlconstui_numClusterK = " << lconstui_numClusterK
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
     
#if  DATATYPE_CENTROIDS_ROUND == 0
    
  T_GENE  lrt_R;
  T_GENE  lrt_delta;

  lrt_R     =  (airt_maxObjetiveFunc != airt_minObjetiveFunc)? 
    (aiochromstr_toMutate.getObjetiveFunc() - airt_minObjetiveFunc) 
    / (airt_maxObjetiveFunc - airt_minObjetiveFunc)
    :1.0;
  std::uniform_real_distribution<T_GENE> uniformdis_deltaRR(-lrt_R, lrt_R);
  lrt_delta = uniformdis_deltaRR(gmt19937_eng);

  if (lrt_delta >= 0.0) {

    interfacesse::aasxpa
      (-lrt_delta,
       aiochromstr_toMutate.getString(),
       lconstui_numClusterK,
       data::Instance<T_GENE>::getNumDimensions(),
       aiarrayt_maxFeatures
       );

  }
  else {

    interfacesse::aasxpa
      (lrt_delta,
       aiochromstr_toMutate.getString(),
       lconstui_numClusterK,
       data::Instance<T_GENE>::getNumDimensions(),
       aiarrayt_minFeatures
       );
  }

#else //DATATYPE_CENTROIDS_ROUND

  DATATYPE_REAL  lrt_R;
  DATATYPE_REAL  lrt_delta;

  lrt_R     = (airt_maxObjetiveFunc != airt_minObjetiveFunc)? 
    (aiochromstr_toMutate.getObjetiveFunc() - airt_minObjetiveFunc) 
    / (airt_maxObjetiveFunc - airt_minObjetiveFunc)
    :1.0;
  std::uniform_real_distribution<DATATYPE_REAL> uniformdis_deltaRR(-lrt_R, lrt_R);
  lrt_delta = uniformdis_deltaRR(gmt19937_eng);
  

  T_GENE* lchromstr_toMutate = aiochromstr_toMutate.getString();
  
  if (lrt_delta >= 0.0) {

      for ( uintidx liu_i = 0; liu_i < lconstui_numClusterK; ++liu_i) {
	
	interfacesse::aysxpy
	  (lchromstr_toMutate,
	   -lrt_delta,
	   aiarrayt_maxFeatures,
	   data::Instance<T_GENE>::getNumDimensions()
	   );
	lchromstr_toMutate +=  data::Instance<T_GENE>::getNumDimensions();
	
      }

  }
  else {

    for ( uintidx liu_i = 0; liu_i < lconstui_numClusterK; ++liu_i) {
	
	interfacesse::aysxpy
	  (lchromstr_toMutate,
	   lrt_delta,
	   aiarrayt_minFeatures,
	   data::Instance<T_GENE>::getNumDimensions()
	   );
	
	lchromstr_toMutate +=  data::Instance<T_GENE>::getNumDimensions();
	
    }
    
  }

#endif //DATATYPE_CENTROIDS_ROUND  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      <<  "\tR = " <<  lrt_R << "\tdelta = " << lrt_delta <<'\n';
    aiochromstr_toMutate.print(std::cout,lpc_labelFunc,',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

} //END mutationChromosomeBandyMaulik2002

  

/*! \fn void biDirectionMutation(gaencode::ChromosomeString<T_FEATURE,T_REAL>& aiochrom_centroids, const T_FEATURE *aiarrayt_minFeatures, const T_FEATURE *aiarrayt_maxFeatures, const T_REAL airt_sign, const T_REAL airt_delta) 
    \brief biDirectionMutation \cite He:Tan:GAclusteringVarK:TGCA:2012
    \details The mutation does the following
    \f[
    \hbox{mutate}(g_{jl}) = \cases{
    g_{jl} + \delta \times ( \hbox{max}(x_l) - g_{jl})  &\quad  if $ \delta_{jl} \geq 0 $, for $j= 1,2,..,k$ and $l= 1,2,...d$, \cr
    g_{jl} + \delta \times (g_{jl} - \hbox{min}(x_l))   &\quad  if $ \delta_{jl} < 0 $. \cr
    }
    \f]       
    Where \f$\delta_{jl}\f$ is a random number in the interval \f$[-1,+1]\f$:
    \param aiochromstr_toMutate a gaencode::ChromosomeString chromosome to mutate
    \param aiarrayt_minFeatures an array with the minimum values of the dataset instances
    \param aiarrayt_maxFeatures an array with the maximum values of the dataset instances
  */
template <typename T_FEATURE,
	  typename T_REAL
	  >
void

biDirectionMutation
(gaencode::ChromosomeString<T_FEATURE,T_REAL>& aiochromstr_toMutate,
 const T_FEATURE                               *aiarrayt_minFeatures,
 const T_FEATURE                               *aiarrayt_maxFeatures
 )
{ //BEGIN biDirectionMutation
     
  uintidx  lcidx_numRows = 
    (aiochromstr_toMutate.getStringSize() / data::Instance<T_FEATURE>::getNumDimensions());
    
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "gaclusteringop::biDirectionMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output gaencode::Chromosome: &aiochromstr_toMutate[" 
	      << &aiochromstr_toMutate << "]\n"
	      << " input lcidx_numRows = " 
	      << lcidx_numRows
	      << "\n)"
	      << std::endl; 
  }
#endif //__VERBOSE_YES

  static std::uniform_real_distribution<T_REAL> lsuniformdis_real11(-1.0,1.0);

  T_FEATURE *lrt_gene = aiochromstr_toMutate.getString();
  
#if  DATATYPE_CENTROIDS_ROUND == 0
  
  for ( uintidx liu_i = 0; liu_i < lcidx_numRows; ++liu_i) {

    for ( uintidx liu_j = 0; liu_j < data::Instance<T_FEATURE>::getNumDimensions(); ++liu_j) {

      T_REAL lrt_delta = lsuniformdis_real11(gmt19937_eng);

      if ( lrt_delta >= 0.0 ) {
	*lrt_gene += lrt_delta * (aiarrayt_maxFeatures[liu_j] - *lrt_gene);
      }
      else {
	*lrt_gene += lrt_delta * (*lrt_gene -  aiarrayt_minFeatures[liu_j]);
      }
      
      ++lrt_gene;
    }
    
  }

#else //DATATYPE_CENTROIDS_ROUND

   for ( uintidx liu_i = 0; liu_i < lcidx_numRows; ++liu_i) {

    for ( uintidx liu_j = 0; liu_j < data::Instance<T_FEATURE>::getNumDimensions(); ++liu_j) {

      T_REAL lrt_delta = lsuniformdis_real11(gmt19937_eng);

      if ( lrt_delta >= 0.0 ) {
	*lrt_gene += T_FEATURE(std::round(lrt_delta * (T_REAL)(aiarrayt_maxFeatures[liu_j] - *lrt_gene)));
      }
      else {
	*lrt_gene += T_FEATURE(std::round(lrt_delta * (T_REAL)(*lrt_gene -  aiarrayt_minFeatures[liu_j])));
      }
      
      ++lrt_gene;
    }
    
  }

#endif //DATATYPE_CENTROIDS_ROUND  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelFunc;
    lostrstream_labelFunc  << "<CENTROIDSCLUSTER:" << geverbosepc_labelstep << ':' << lpc_labelFunc; 
    aiochromstr_toMutate.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

} /*END biDirectionMutation
   */

 /*OTHERE OPERATOR ------------------------------------------------------------------------------
   */

/*! \fn bool rearrangedCluster(gaencode::ChromosomeString<T_CENTROIDS,T_REAL> &aochrom_YRearranged, const mat::MatrixRow<T_CENTROIDS> &aimatrixrowt_XReferenced, const dist::Dist<T_REAL,T_CENTROIDS> &aifunc2p_dist)
  \brief Rearranged cluster \cite Chang:etal:GAclustering:GAGR:2009
  \details For the new population, select the best chromosome as as a refrence, wich other chromosomes might fall into the gene rearrangement if needed 
  \param aochrom_YRearranged
  \param aimatrixrowt_XReferenced
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CENTROIDS,
	  typename T_REAL
	  >   
bool
rearrangedCluster
(gaencode::ChromosomeString<T_CENTROIDS,T_REAL> &aochrom_YRearranged,       
 const mat::MatrixRow<T_CENTROIDS>              &aimatrixrowt_XReferenced,
 const dist::Dist<T_REAL,T_CENTROIDS>           &aifunc2p_dist
 )
{
  bool      lob_isNeededRearrange;
  bool      lb_changeRow;
  uintidx   luintidx_i, luintidx_j;
  uintidx   luintid_minYiXj;
  T_REAL    lrt_distMinYiXj;  
  T_REAL    lrt_distYjXi;

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::rearrangedClusterGAGR";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelChrom;
    lostrstream_labelChrom << "input:YREARRANGED:" << lpc_labelFunc;
    aochrom_YRearranged.print(std::cout,lostrstream_labelChrom.str().c_str());
    std::cout << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  mat::BitArray<unsigned long> lbitarray_p2(aimatrixrowt_XReferenced.getNumRows());


  /*DECODE CHROMOSOME*/
  mat::MatrixRow<T_CENTROIDS> 
    lmatrixrowt_YRearranged
    (aimatrixrowt_XReferenced.getNumRows(),
     data::Instance<T_CENTROIDS>::getNumDimensions(),
     aochrom_YRearranged.getString()
     );

  lrt_distMinYiXj = 0.0;
  
  lbitarray_p2.initialize();

  lob_isNeededRearrange = false;
  luintid_minYiXj = 0;
  luintidx_i = 0;
  while (luintidx_i < lmatrixrowt_YRearranged.getNumRows()) {
    if (!lbitarray_p2.getBit(luintidx_i)) {
      luintidx_j = 0;
      lb_changeRow = false;
      while ( lbitarray_p2.getBit(luintidx_j) && (luintidx_j < lmatrixrowt_YRearranged.getNumRows()) )
	++luintidx_j;
      if  ( luintidx_j < lmatrixrowt_YRearranged.getNumRows() ) {
	
	lrt_distMinYiXj = 
	  aifunc2p_dist
	  (lmatrixrowt_YRearranged.getRow(luintidx_i), 
	   aimatrixrowt_XReferenced.getRow(luintidx_j), 
	   aimatrixrowt_XReferenced.getNumColumns()
	   ); 
	luintid_minYiXj = luintidx_j;
	lb_changeRow  = true;
      }
      for (luintidx_j = luintidx_j +1;
	   luintidx_j < lmatrixrowt_YRearranged.getNumRows();
	   luintidx_j++)
	{  
	if (!lbitarray_p2.getBit(luintidx_j)) { 
	  

	  lrt_distYjXi = 
	    aifunc2p_dist
	    (lmatrixrowt_YRearranged.getRow(luintidx_i), 
	     aimatrixrowt_XReferenced.getRow(luintidx_j), 
	     aimatrixrowt_XReferenced.getNumColumns()
	     );
	  if ( lrt_distYjXi <  lrt_distMinYiXj) {
	    lrt_distMinYiXj = lrt_distYjXi;
	    luintid_minYiXj  = luintidx_j;
	    lb_changeRow   = true;
	  }
	} 
      }
      
      if ( lb_changeRow && luintidx_i != luintid_minYiXj) {
	lmatrixrowt_YRearranged.swapRows(luintidx_i,luintid_minYiXj);
	lbitarray_p2.setBit(luintid_minYiXj); // = true;
	lob_isNeededRearrange = true;        
	luintidx_i = 0;
      }
    }
    ++luintidx_i;
  } //WHILE i    


#ifdef __VERBOSE_YES
   if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
	      
    std::ostringstream lostrstream_labelChrom;
    lostrstream_labelChrom << "output:YREARRANGED:NeededRearrange:" << lob_isNeededRearrange 
			   <<  ':' << lpc_labelFunc;
    aochrom_YRearranged.print(std::cout,lostrstream_labelChrom.str().c_str());
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lob_isNeededRearrange;
}

  
/*! \fn void kmeansfeac(gaencode::ChromosomeFEAC<T_CLUSTERIDX,T_REAL,T_FEATURE,T_INSTANCES_CLUSTER_K> &aochrom_feac, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, COMMON_IDOMAIN aiit_kmeansNumMaxIter, T_FEATURE airt_kmeansMaxDiffCent, dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist ) 
  \brief kmeansfeac \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
  \brief 
  \param aochrom_feac a gaencode::ChromosomeFEAC
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aiit_kmeansNumMaxIter an integer number
  \param airt_kmeansMaxDiffCent a real number
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename  T_CLUSTERIDX,
	   typename T_REAL,
	   typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_FEATURE_SUM,
	  typename T_COMMON_IDOMAIN,
	   typename INPUT_ITERATOR
	   > 
void
kmeansfeac
(gaencode::ChromosomeFEAC
 <T_CLUSTERIDX,
  T_REAL,
  T_FEATURE,
  T_FEATURE_SUM,
  T_INSTANCES_CLUSTER_K>             &aochrom_feac,
 INPUT_ITERATOR                      aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 T_COMMON_IDOMAIN                    aiit_kmeansNumMaxIter,
 T_FEATURE                           airt_kmeansMaxDiffCent,
 dist::Dist<T_REAL,T_FEATURE>        &aifunc2p_dist
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::kmeansfeac";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelFunc;
    lostrstream_labelFunc  << geverbosepc_labelstep << ':' << lpc_labelFunc; 
    aochrom_feac.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
    std::cout << "\ninput  aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << "input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
      	      << "input dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      <<  "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  mat::MatrixWithRowNull<T_FEATURE>         &lmatrixwrownull_centroids =
    aochrom_feac.getCentroids();
  T_CLUSTERIDX                              *larraycidx_idxMemberShip =
    aochrom_feac.getString();
  std::vector<T_INSTANCES_CLUSTER_K>        &lvectorT_numInstancesClusterK = 
    aochrom_feac.getNumInstancesClusterK();

  T_COMMON_IDOMAIN liT_newIteration = 1;
  mat::BitArray<unsigned long> lbitarray_recalculeteK(lmatrixwrownull_centroids.getNumRows());
  T_INSTANCES_CLUSTER_K *larrayt_instacesDiffNumClusterK =
    new T_INSTANCES_CLUSTER_K[lmatrixwrownull_centroids.getNumRows()];
  
  mat::MatrixRow<T_FEATURE_SUM>       
    lmatrixrowt_centroidsDiffSumInsta
    (lmatrixwrownull_centroids.getNumRows(), 
     lmatrixwrownull_centroids.getNumColumns()
     );
	      
  uintidx lui_numClusterNull = 0;
  
  while ( liT_newIteration == 1 ) {

    liT_newIteration = 0;
    
    lbitarray_recalculeteK.initialize();

    lmatrixrowt_centroidsDiffSumInsta.initialize();
    
    interfacesse::copya
      (larrayt_instacesDiffNumClusterK,
       T_INSTANCES_CLUSTER_K(0),
       lmatrixwrownull_centroids.getNumRows()
       );
    
    /*3. Assign each object to the cluster with the nearest centroid;
     */
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "KMEANS-FEAC 3. ASSIGN EACH OBJECT TO THE CLUSTER WITH THE NEAREST CENTROID:  IN"
              << '(' << geiinparam_verbose << ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

    mat::MatrixWithRowNull<T_FEATURE> lmatrixwrownullT_oldCentroids(lmatrixwrownull_centroids);
    for (uintidx lui_idxInsti = 0 ; aiiterator_instfirst != aiiterator_instlast;
	 aiiterator_instfirst++, lui_idxInsti++) {
      const T_FEATURE* liter_iInstance = (*aiiterator_instfirst)->getFeatures();

      T_REAL lrt_distMinCentInst;
      T_CLUSTERIDX lmgidx_j =
	nearest::checkNullCentroidsNN
	<T_CLUSTERIDX, 
	 T_FEATURE, 
	 T_REAL>
	(lrt_distMinCentInst,
	 lmatrixwrownull_centroids,
	 liter_iInstance,
	 aifunc2p_dist
	 );
      if ( lmgidx_j == NEARESTCENTROID_UNKNOWN )
	throw  std::range_error
	  ("gaclusteringop::kmeansfeac: instance without group, or centroids at infinity");

      if ( larraycidx_idxMemberShip[lui_idxInsti] == NEARESTCENTROID_UNKNOWN) {
	
        /*update cromosome*/
       	larraycidx_idxMemberShip[lui_idxInsti] = lmgidx_j;
	
	/*update diff centroide*/
	interfacesse::axpy
	  (lmatrixrowt_centroidsDiffSumInsta.getRow(lmgidx_j),
	   T_FEATURE(1),
	   liter_iInstance,
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	++larrayt_instacesDiffNumClusterK[lmgidx_j];
	
	lbitarray_recalculeteK.setBit(lmgidx_j);

      } /*end if NEARESTCENTROID_UNKNOWN*/
      else { /*instance in cluster k*/

	if ( lmgidx_j != larraycidx_idxMemberShip[lui_idxInsti] ) {

	  /*update diff centroide*/
	  interfacesse::axpy
	    (lmatrixrowt_centroidsDiffSumInsta.getRow(larraycidx_idxMemberShip[lui_idxInsti]),
	     T_FEATURE(-1),
	     liter_iInstance,
	     data::Instance<T_FEATURE>::getNumDimensions()
	     );
	  --larrayt_instacesDiffNumClusterK[larraycidx_idxMemberShip[lui_idxInsti]];
	  lbitarray_recalculeteK.setBit(larraycidx_idxMemberShip[lui_idxInsti]);

	  /*update cromosome*/
	  larraycidx_idxMemberShip[lui_idxInsti] = lmgidx_j;

	  /*update diff centroide*/
	  interfacesse::axpy
	    (lmatrixrowt_centroidsDiffSumInsta.getRow(lmgidx_j),
	     T_FEATURE(1),
	     liter_iInstance,
	     data::Instance<T_FEATURE>::getNumDimensions()
	     );
	  ++larrayt_instacesDiffNumClusterK[lmgidx_j];

	  lbitarray_recalculeteK.setBit(lmgidx_j);
	}
	
      } /*end else*/
    } /*end for clustering instances*/
   
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "KMEANS-FEAC 3. ASSIGN EACH OBJECT TO THE CLUSTER WITH THE NEAREST CENTROID:  OUT"
              << '(' << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDSDIFFSUM:" << lpc_labelFunc;
    lmatrixrowt_centroidsDiffSumInsta.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    std::cout <<'\n';
    std::ostringstream lostrstream_labelinstacesDiffNumClusterK;
    lostrstream_labelinstacesDiffNumClusterK
      << "<INSTACESDIFFNUMCLUSTERK:"
      << lpc_labelFunc;
    inout::containerprint
      (larrayt_instacesDiffNumClusterK,
       larrayt_instacesDiffNumClusterK+ lmatrixwrownull_centroids.getNumRows(),
       std::cout,
       lostrstream_labelinstacesDiffNumClusterK.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
    /*for moves between cluster instances
     */
    for (uintidx lui_j = 0; lui_j < lmatrixwrownull_centroids.getNumRows(); lui_j++) {

      if ( lbitarray_recalculeteK.getBit(lui_j) ) {
	
	interfacesse::scal
	  (lmatrixwrownull_centroids.getRow(lui_j),
	   lvectorT_numInstancesClusterK.at(lui_j),
	   lmatrixwrownull_centroids.getNumColumns()
	   );
	
	interfacesse::axpy
	  (lmatrixwrownull_centroids.getRow(lui_j), 
	   T_FEATURE(1),                            
	   lmatrixrowt_centroidsDiffSumInsta.getRow(lui_j),
	   lmatrixwrownull_centroids.getNumColumns()
	   );

	lvectorT_numInstancesClusterK.at(lui_j) += 
	  larrayt_instacesDiffNumClusterK[lui_j];

	/*What to do when a cluster has no instances assigned?
	 */
	if ( lvectorT_numInstancesClusterK[lui_j] != 0 ) {

	  interfacesse::scalInv
	    (lmatrixwrownull_centroids.getRow(lui_j),
	     lvectorT_numInstancesClusterK.at(lui_j),
	     lmatrixwrownull_centroids.getNumColumns()
	     );
	  
	}
	else {
	  
	  ++lui_numClusterNull;
	  interfacesse::copya
	    (lmatrixwrownull_centroids.getRow(lui_j),
	     data::Instance<T_FEATURE>::getInfiniteFeature(), 
	     lmatrixwrownull_centroids.getNumColumns()
	     );
	}
      }
    } /*for moves between cluster instances
       */
 
    /*calculated if another iteration is required
     */
    for ( uintidx lui_j = 0; lui_j < lmatrixwrownull_centroids.getNumRows(); lui_j++) {
      T_FEATURE *larrayrowt_centroidj = lmatrixwrownull_centroids.getRow(lui_j); 
      if (!data::Instance<T_FEATURE>::isInfiniteFeature(larrayrowt_centroidj[0]) &&  
	  aifunc2p_dist
	  (larrayrowt_centroidj, 
	   lmatrixwrownullT_oldCentroids.getRow(lui_j), 
	   lmatrixwrownull_centroids.getNumColumns()
	   ) > airt_kmeansMaxDiffCent
	  )
	liT_newIteration = 1;
     
    }
     
    --aiit_kmeansNumMaxIter;
    if (aiit_kmeansNumMaxIter == 0) 
      liT_newIteration = 0;

  } /*while*/

  if ( lui_numClusterNull > 0 ) {

    std::vector<uintidx>  
      lvectorstidx_clustersKeep;
    
    lvectorstidx_clustersKeep.reserve
      (lmatrixwrownull_centroids.getNumRows()); 

    for (uintidx lst_k = 0; lst_k < lvectorT_numInstancesClusterK.size(); lst_k++) 
      if ( lvectorT_numInstancesClusterK[lst_k] != 0) 
	lvectorstidx_clustersKeep.push_back( lst_k );
   
    gaintegerop::labelKeep
      <T_CLUSTERIDX,T_REAL>
      (aochrom_feac,
       lvectorstidx_clustersKeep
       );
    
    lmatrixwrownull_centroids.keepRows(lvectorstidx_clustersKeep);

    lvectorT_numInstancesClusterK =
      vectorutils::keepItems
      <T_INSTANCES_CLUSTER_K>
      (lvectorT_numInstancesClusterK,
       lvectorstidx_clustersKeep
       );

    std::vector<T_REAL>& lvectorT_genotypePartialFcC =
	aochrom_feac.getPartialFcC();

    lvectorT_genotypePartialFcC =
      vectorutils::keepItems
      (lvectorT_genotypePartialFcC,
       lvectorstidx_clustersKeep
       );
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
     std::cout << lpc_labelFunc
	       << ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelFunc;
      lostrstream_labelFunc  << geverbosepc_labelstep << ':' << lpc_labelFunc; 
      aochrom_feac.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  

  delete[] larrayt_instacesDiffNumClusterK;
}


/*! \fn void MO1(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> &aiochromfixlen_toMutate, mat::MatrixRow<T_FEATURE> &aomatrixrowt_centroids, mat::MatrixRow<T_FEATURE_SUM> &aomatrixrowt_sumInstCluster, std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstClusterK, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief MO1 mutation operator \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003
  \details It eliminates one selected clusters, placing its objects into the nearest remaining clusters (according to their centroids).To remove a group it is necessary to calculate the centroids associated with the chromosome. The parameters of aomatrixt_meanCentroids, aimatrixtrowt_sumInstancesClusterK, and aovectort_numInstClusterK are returned with the centroid, the sum of instances, and the number of instances per group, respectively.
  \param aiochromfixlen_toMutate a gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>
  \param aomatrixrowt_centroids a MatrixBase<T_CENTROIDS>
  \param aomatrixrowt_sumInstCluster a MatrixBase<T_FEATURE_SUM>
  \param aovectort_numInstClusterK a std::vector<T_INSTANCES_CLUSTER_K>
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX,
	  typename T_REAL,//DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename INPUT_ITERATOR
	  >
void
MO1
(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>   &aiochromfixlen_toMutate,
 mat::MatrixRow<T_FEATURE>                         &aomatrixrowt_centroids,
 mat::MatrixRow<T_FEATURE_SUM>                     &aomatrixrowt_sumInstCluster,
 std::vector<T_INSTANCES_CLUSTER_K>                &aovectort_numInstClusterK,
 INPUT_ITERATOR                                    aiiterator_instfirst,
 const INPUT_ITERATOR                              aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE>                &aifunc2p_dist
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::MO1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ":  IN(" << geiinparam_verbose << ')'
      << " K = "
      << aiochromfixlen_toMutate.getGene
      (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
      << '\n';
    std::ostringstream lostrstream_labeltoMutate;
    lostrstream_labeltoMutate << "IN:" << lpc_labelFunc << ":aiochromfixlen_toMutate:";
    aiochromfixlen_toMutate.print(std::cout,lostrstream_labeltoMutate.str().c_str());
    std::cout
     
      << "\n input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" << &aifunc2p_dist << ']'
      <<  "\n)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  T_CLUSTERIDX lcidx_Cs = T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN);
  if ( aiochromfixlen_toMutate.getGene
       (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1) > 2
       )
    {

    /*1.2.1 Randomly choose a cluster Cs encoded into g //
     */
      
    std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_0K1
      (0,aiochromfixlen_toMutate.getGene(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)-1);

    lcidx_Cs = uniformdis_cidx_0K1(gmt19937_eng);

    //Relabel genotype
    for ( auto liter_gene = aiochromfixlen_toMutate.begin();
	  liter_gene != aiochromfixlen_toMutate.end()-1;
	  liter_gene++)
      {
      if ( *liter_gene > lcidx_Cs ) {
	--(*liter_gene);
      }
      else if ( *liter_gene == lcidx_Cs ) { 
	*liter_gene = T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN);
      }
    }

    partition::PartitionLabel
      <T_CLUSTERIDX>
      lpartition_clusters
      (aiochromfixlen_toMutate.getString(),
       gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
       (aiochromfixlen_toMutate.getGene
	(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)-1)
       );

    mat::MatrixRow<T_FEATURE> 
      lmatrixrowt_centroids
      (uintidx
       (aiochromfixlen_toMutate.getGene
	(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)-1
	),
       data::Instance<T_FEATURE>::getNumDimensions() 
       );

    mat::MatrixRow<T_FEATURE_SUM>       
      lmatrixrowt_sumInstCluster
      (lmatrixrowt_centroids.getNumRows(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
	
    std::vector<T_INSTANCES_CLUSTER_K> 
      lvectort_numInstClusterK(lmatrixrowt_centroids.getNumRows()); 

    clusteringop::getCentroids
      (lmatrixrowt_centroids,
       lmatrixrowt_sumInstCluster,
       lvectort_numInstClusterK,
       lpartition_clusters,
       aiiterator_instfirst,
       aiiterator_instlast 
       );
  
    clusteringop::setUpCuster
      (aiochromfixlen_toMutate.getString(),
       lmatrixrowt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aifunc2p_dist
       );
    
    aiochromfixlen_toMutate.setGene
    (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
     aiochromfixlen_toMutate.getGene
	(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)-1
     );

    std::swap(aomatrixrowt_centroids,lmatrixrowt_centroids);
    std::swap(aomatrixrowt_sumInstCluster,lmatrixrowt_sumInstCluster);
    std::swap(aovectort_numInstClusterK,lvectort_numInstClusterK);
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " K = "
      << aiochromfixlen_toMutate.getGene
      (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
      << " lcidx_Cs = " << lcidx_Cs
      << '\n';
    std::ostringstream lostrstream_labeltoMutate;
    lostrstream_labeltoMutate << "OUT:" << lpc_labelFunc << ":aiochromfixlen_toMutate:";
    aiochromfixlen_toMutate.print(std::cout,lostrstream_labeltoMutate.str().c_str());
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*END MO1
   */

/*! \fn void MO2(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> &aiochromfixlen_toMutate, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief MO2 mutation operator \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003
  \details It splits one selected cluster, each of which into two new clusters.
  \param aiochromfixlen_toMutate  a gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>
  \param aivectorptinst_instances a std::vector<data::Instance<T_FEATURE>* >
  \param aifunc2p_dist a dist::Dist<T_REAL,T_FEATURE>
 */
template <typename T_CLUSTERIDX,
	  typename T_REAL,
	  typename T_FEATURE,
	  typename INPUT_ITERATOR
	  >
void
MO2
(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>   &aiochromfixlen_toMutate,
 const INPUT_ITERATOR                              aiiterator_instfirst,
 const dist::Dist<T_REAL,T_FEATURE>                &aifunc2p_dist
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::MO2";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ":  IN(" << geiinparam_verbose << ')'
      << " K = "
      << aiochromfixlen_toMutate.getGene
      (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
        << '\n';
    std::ostringstream lostrstream_labeltoMutate;
    lostrstream_labeltoMutate << "IN:" << lpc_labelFunc << ":aiochromfixlen_toMutate:";
    aiochromfixlen_toMutate.print(std::cout,lostrstream_labeltoMutate.str().c_str());
    std::cout
      << "\n input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" << &aifunc2p_dist << ']'
      <<  "\n)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  /*1.2.1 Randomly choose a cluster Cs encoded into g //
   */
  
  std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_cidx_0K1
    (0,aiochromfixlen_toMutate.getGene(gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)-1);

  T_CLUSTERIDX lcidx_Cs = uniformdis_cidx_0K1(gmt19937_eng);

  std::list<uintidx> llist_idxInstanceInCs;
  uintidx luintidx_iInst = 0;
  for ( auto liter_gene = aiochromfixlen_toMutate.begin();
	liter_gene != aiochromfixlen_toMutate.end()-1;
	liter_gene++
	)
    {
    if ( *liter_gene == lcidx_Cs ) {
      llist_idxInstanceInCs.push_back(luintidx_iInst);
    }
    ++luintidx_iInst;
  }
  
  if (llist_idxInstanceInCs.size() > 2 ) {
      
    std::uniform_int_distribution<uintidx> uniformdis_consecutiveS1
      (0, uintidx(llist_idxInstanceInCs.size()-1) );
      
    uintidx luintidx_S1 =
      uniformdis_consecutiveS1(gmt19937_eng);
    auto literlist_idxS1 = llist_idxInstanceInCs.begin();
    std::advance(literlist_idxS1,luintidx_S1);
    luintidx_S1  = *literlist_idxS1;

    uintidx luintidx_S2 = 
      nearest::farthestInstanceFromS1
      (luintidx_S1,
       llist_idxInstanceInCs.begin(),
       llist_idxInstanceInCs.end(),
       aiiterator_instfirst,
       aifunc2p_dist
       );

    T_CLUSTERIDX lcidx_newCluster = 
      aiochromfixlen_toMutate.getGene
      (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1);

    const T_FEATURE* linst_featuresS1 =
      (*std::next(aiiterator_instfirst,luintidx_S1))->getFeatures();
      
    const T_FEATURE* linst_featuresS2 =
      (*std::next(aiiterator_instfirst,luintidx_S2))->getFeatures();
      
    for (auto lui_idxInstanceInCs: llist_idxInstanceInCs) {

      T_REAL lrt_distS1toSi = 
	aifunc2p_dist 
	(linst_featuresS1,
	 (*std::next(aiiterator_instfirst,lui_idxInstanceInCs))->getFeatures(),
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );

      T_REAL lrt_distS2toSi = 
	aifunc2p_dist 
	(linst_featuresS2,
	 (*std::next(aiiterator_instfirst,lui_idxInstanceInCs))->getFeatures(),
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
    
      if ( lrt_distS2toSi < lrt_distS1toSi ) {
	aiochromfixlen_toMutate.setGene(lui_idxInstanceInCs,lcidx_newCluster);
      }
      
    }
    //Code in string number of cluster
    aiochromfixlen_toMutate.setGene
      (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1,
       lcidx_newCluster+1
       );
  }
    

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " K = "
      << aiochromfixlen_toMutate.getGene
      (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>::stcgetStringSize()-1)
      << " lcidx_Cs = " << lcidx_Cs
      << '\n';
    std::ostringstream lostrstream_labeltoMutate;
    lostrstream_labeltoMutate << "OUT:" << lpc_labelFunc << ":aiochromfixlen_toMutate:";
    aiochromfixlen_toMutate.print(std::cout,lostrstream_labeltoMutate.str().c_str());
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*END MO2
   */


/*! \fn void MO1(gaencode::ChromosomeFEAC<T_CLUSTERIDX,T_REAL,T_FEATURE,T_INSTANCES_CLUSTER_K> &aochrom_feac, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
  \brief MO1 mutation operator
  \details It eliminates one or more randomly selected clusters, placing its objects into the nearest remaining clusters (according to their centroids).
  \param aiochromfixlen_toMutate a gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>
  \param aivectorptinst_instances a std::vector<data::Instance<T_FEATURE>* >
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
 */
template <typename T_CLUSTERIDX,
	  typename T_REAL, //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  typename T_FEATURE,
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K,
	  typename INPUT_ITERATOR
	  >
void
MO1
(gaencode::ChromosomeFEAC
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>          &aochromfeac_toMutate,
 INPUT_ITERATOR                  aiiterator_instfirst,
 const INPUT_ITERATOR            aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>    &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::MO1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << " K = " << aochromfeac_toMutate.getNumClusterK()
	      << " fitness = " <<  aochromfeac_toMutate.getFitness() << '\n';
    aochromfeac_toMutate.print(std::cout,lpc_labelFunc,',',';');
    std::cout << "input const aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << "input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
      	      << "\n input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      <<  "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  
  if ( aochromfeac_toMutate.getNumClusterK() > 2 )  {

    std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
    mat::MatrixWithRowNull<T_FEATURE> &lmatrixwrownull_genotypeCentroids = 
      aochromfeac_toMutate.getCentroids();
  
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "MO1 1.1 RANDOMLY GENERATE A NUMBER n in {1,..., k_g-2}:  IN"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl; 
    }
#endif /*__VERBOSE_YES*/

   /*1.1 Randomly generate a number n \in {1,..., k_g-2}; 
      at least two clusters must remain in g 
    */
    std::uniform_int_distribution<uintidx> uniformdis_uiN
      (1,lmatrixwrownull_genotypeCentroids.getNumRows()- (uintidx) 2);
    uintidx luintidx_n = uniformdis_uiN(gmt19937_eng);
      
       
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "MO1 1.1 RANDOMLY GENERATE A NUMBER n in {1,..., k_g-2}: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< "\nn = " << luintidx_n << " in {1,...,"
		<< lmatrixwrownull_genotypeCentroids.getNumRows()-2  << '}'
		<< std::endl; 
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    
    /*1.2 For i = (1,...,n) do:
     */
    for (uintidx  lui_i = 1; lui_i <= luintidx_n; lui_i++) {
      
      /*1.2.1 Randomly choose a cluster Cs encoded into g 
	is lvectoruintidx_Cs[lui_i] for strategy without
	replacement
      */

      std::vector<T_REAL>& lvectorrt_genotypePartialFcC =
	aochromfeac_toMutate.getPartialFcC();

      // { BEGIN 1.2.1
	 	
#ifdef __VERBOSE_YES
      const char* lpc_labelFuncStep = "MO1 1.2.1 RANDOMLY CHOOSE A CLUSTER CS ENCODED INTO G";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<  lpc_labelFuncStep
		<< ":  IN(" << geiinparam_verbose << ')'
		<< std::endl; 
    }
#endif /*__VERBOSE_YES*/
      
     
#if defined(ALG_EAC_I_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006) || \
  defined(ALG_EAC_III_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006)

    /*\cite{Alves:etal:GAclusteringLabelKVar:FEAC:2006}
      EAC-I
      EAC-III
    */
    
    std::vector<T_REAL> lvectorT_partialFcCNorma01(lvectorrt_genotypePartialFcC.size(),1);
    /*1 - fc(Ci)
     */
    interfacesse::axpy
      (lvectorT_partialFcCNorma01.data(),
       -1.0,
       lvectorrt_genotypePartialFcC.data(),
       (uintidx) lvectorT_partialFcCNorma01.size()
       );
    /* \sum(1-fc(Ci))
     */
    T_REAL lrt_sum1_fc =
      interfacesse::sum
      (lvectorT_partialFcCNorma01.data(),
       (uintidx) lvectorT_partialFcCNorma01.size()
       );
    lrt_sum1_fc = 1.0 / lrt_sum1_fc;
  
    interfacesse::scal
      (lvectorT_partialFcCNorma01.data(),
       lrt_sum1_fc,
       (uintidx) lvectorT_partialFcCNorma01.size()
       );

    std::vector<T_REAL>&& lvectorrt_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorT_partialFcCNorma01.begin(),
	 lvectorT_partialFcCNorma01.end(),
	 [](const T_REAL& lt_metrici) -> T_REAL
	 {
	   return lt_metrici;
	 }
	 );
      
    T_CLUSTERIDX lcidx_Cs =
	gaselect::getIdxRouletteWheel
	(lvectorrt_probDistRouletteWheel,
	 T_CLUSTERIDX(0)
	 );
        
    /*END:
      ALG_EAC_I_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
      ALG_EAC_III_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
    */
    
#elif defined(ALG_EAC_II_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006) || \
  defined(ALG_FEAC_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006)  
      /*\cite{Alves:etal:GAclusteringLabelKVar:FEAC:2006}
	EAC-II
	F-EAC
      */
    
      std::vector<T_REAL>&& lvectorT_partialFcCLinearNorm =
	prob::linearNormalization
	(lvectorrt_genotypePartialFcC.begin(),
	 lvectorrt_genotypePartialFcC.end(),
	 [](const T_REAL& lt_metrici) -> T_REAL
	 {
	   return lt_metrici;
	 }
	 );
	
      std::vector<T_REAL>&& lvectorrt_probDistRouletteWheel =
	prob::makeDistRouletteWheel
	(lvectorT_partialFcCLinearNorm.begin(),
	 lvectorT_partialFcCLinearNorm.end(),
	 [](const T_REAL& lt_metrici) -> T_REAL
	 {
	   return lt_metrici;
	 }
	 );

      T_CLUSTERIDX lcidx_Cs =
	gaselect::getIdxRouletteWheel
	(lvectorrt_probDistRouletteWheel,
	T_CLUSTERIDX(0)
	);
      
      /*END:
	ALG_EAC_II_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
	ALG_FEAC_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006  
      */
      
#else /*BEGIN: ALG_EAC_CLUSTERINGLABELVARK_HRUSCHKA_CAMPELLO_CASTRO2006 || OTHER
       */

      /*\cite{Hruschka:Campello:Castro:GAClusteringLabelKVar:EAC:2006}
	EAC
      */
      std::uniform_int_distribution<T_CLUSTERIDX> uniformdis_0K
	(0,(T_CLUSTERIDX ) aochromfeac_toMutate.getNumClusterK()-1);

      T_CLUSTERIDX lcidx_Cs =
	uniformdis_0K(gmt19937_eng);
      
#endif /*END: 
	 ALG_EAC_CLUSTERINGLABELVARK_HRUSCHKA_CAMPELLO_CASTRO2006
	*/

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFuncStep
		<< ": OUT(" << geiinparam_verbose << ')'
		<< "\nCs = " << lcidx_Cs
		<< std::endl; 
    }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    
    // }  END 1.2.1
       

      
      /*1.2.2 Place the objects that belong to C_s into the
	nearest remaining clusters C_j \in g, j \neq s,
	according to the distances between objects
	and cluster centroids (mean vectors of C_j);
	C_s --> C_j
      */
      {/*BEGIN 1.2.2
	*/
#ifdef __VERBOSE_YES
	const char* lpc_labelFuncStep = "MO1 1.2.2 PLACE THE OBJECTS A CLUSTER CS";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout <<  lpc_labelFuncStep
		    << ":  IN(" << geiinparam_verbose << ')'
		    << std::endl; 
	}
#endif /*__VERBOSE_YES*/

	std::vector<T_INSTANCES_CLUSTER_K>& lvectorT_genotypeNumInstancesClusterK =
	  aochromfeac_toMutate.getNumInstancesClusterK();

	T_CLUSTERIDX lcidx_Cj =
	  nearest::findTwoNearestCentroids
	  (lcidx_Cs,
	   lmatrixwrownull_genotypeCentroids,
	   aifunc2p_dist
	   );
      
	T_INSTANCES_CLUSTER_K lT_denominator = 
	  lvectorT_genotypeNumInstancesClusterK.at(lcidx_Cj) + 
	  lvectorT_genotypeNumInstancesClusterK.at(lcidx_Cs);
	
	if ( lT_denominator != 0 ) {

	  interfacesse::scal
	    (lmatrixwrownull_genotypeCentroids.getRow(lcidx_Cj),
	     lvectorT_genotypeNumInstancesClusterK.at(lcidx_Cj),//lT_beta,
	     lmatrixwrownull_genotypeCentroids.getNumColumns()
	     );
	  
	  interfacesse::axpy
	    (lmatrixwrownull_genotypeCentroids.getRow(lcidx_Cj),
	     lvectorT_genotypeNumInstancesClusterK.at(lcidx_Cs),
	     lmatrixwrownull_genotypeCentroids.getRow(lcidx_Cs),
	     lmatrixwrownull_genotypeCentroids.getNumColumns()
	     );

	  interfacesse::scalInv
	    (lmatrixwrownull_genotypeCentroids.getRow(lcidx_Cj),
	     lT_denominator,
	     lmatrixwrownull_genotypeCentroids.getNumColumns()
	     );
	  
	  lvectorT_genotypeNumInstancesClusterK.at(lcidx_Cj) += 
	    lvectorT_genotypeNumInstancesClusterK.at(lcidx_Cs);
	  
	  std::vector<uintidx>  
	    lvectorstidx_clustersKeep;

	  lvectorstidx_clustersKeep.reserve
	    (lmatrixwrownull_genotypeCentroids.getNumRows() - (uintidx) 1);

	  for ( uintidx lui_k = 0; lui_k < lmatrixwrownull_genotypeCentroids.getNumRows(); lui_k++) 
	    
	    if ( lui_k != (uintidx) lcidx_Cs ) 
	      lvectorstidx_clustersKeep.push_back( lui_k );

	  /*UPDATE LABELS
	   */
	  gaintegerop::oneChangelabel
	    <T_CLUSTERIDX,T_REAL>
	    (aochromfeac_toMutate,
	     lcidx_Cs,
	     lcidx_Cj
	     );

	  gaintegerop::labelKeep
	    <T_CLUSTERIDX,T_REAL>
	    (aochromfeac_toMutate,
	     lvectorstidx_clustersKeep
	     );

	  /*DELETE CENTROIDE
	   */
	  lmatrixwrownull_genotypeCentroids.keepRows
	    (lvectorstidx_clustersKeep);

	  lvectorT_genotypeNumInstancesClusterK.erase
	    (lvectorT_genotypeNumInstancesClusterK.begin() + lcidx_Cs);
	  lvectorrt_genotypePartialFcC.erase
	    (lvectorrt_genotypePartialFcC.begin() + lcidx_Cs);

	  

#ifdef __FITNESS_SIMPLIFIED_SILHOUETTE__
      
	  partition::PartitionLabel
	    <T_CLUSTERIDX>
	    lpartition_clustersLabel
	    (aochromfeac_toMutate.getString(),
	     aochromfeac_toMutate.getStringSize(),
	     aochromfeac_toMutate.getNumClusterK()
	     );

	     std::vector<T_REAL>&&  lvectort_partialSilhouette =
	       um::simplifiedSilhouette
	       (aochromfeac_toMutate.getCentroids(),
		aiiterator_instfirst,
		aiiterator_instlast,
		lpartition_clustersLabel,
		aochromfeac_toMutate.getNumInstancesClusterK(),
		aifunc2p_dist
		);
	   

		T_REAL lmetrict_partialSilhouette = 
		  interfacesse::sum
		  (lvectort_partialSilhouette.data(),
		   (uintidx) lvectort_partialSilhouette.size()
		   );

		lmetrict_partialSilhouette /= (T_REAL) lvectort_partialSilhouette.size();

		aochromfeac_toMutate.setPartialFcC(lvectort_partialSilhouette);
		
		aochromfeac_toMutate.setObjetiveFunc(lmetrict_partialSilhouette);
		

#endif /*__FITNESS_SIMPLIFIED_SILHOUETTE__*/
      
	
	} /*End if*/
	
#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << lpc_labelFuncStep
		      << ": OUT(" << geiinparam_verbose << ")\n";
	    std::ostringstream lostrstream_labelSubStep;
	    lostrstream_labelSubStep << geverbosepc_labelstep
				     << ':' << lpc_labelFunc       
				     << ':' << lpc_labelFuncStep;  //":1.2.2 PLACE THE OBJECTS A CLUSTER CS";
	    aochromfeac_toMutate.print(std::cout,lostrstream_labelSubStep.str().c_str(),',',';');
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      } /*END 1.2.2
	 */
    } /*End for*/
    
  } /*Else g is not mutated*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " K = " << aochromfeac_toMutate.getNumClusterK() << '\n';        
    aochromfeac_toMutate.print(std::cout,lpc_labelFunc,',',';');
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*END MO1
   */

  
/*! \fn void MO2(gaencode::ChromosomeFEAC<T_CLUSTERIDX,T_REAL,T_FEATURE,T_INSTANCES_CLUSTER_K> &aochromfeac_toMutate, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist )
  \brief MO2 mutation operator
  \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
  \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
  \details It splits one or more randomly selected clusters, each of which into two new clusters.
  \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
  \param aochromfeac_toMutate a gaencode::ChromosomeFEAC
  \param aivectorptinst_instances a std::vector of pointers to data::Instance
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX,
	  typename T_REAL, //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  typename T_FEATURE,
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K,
	  typename INPUT_ITERATOR
	  >
void
MO2
(gaencode::ChromosomeFEAC
 <T_CLUSTERIDX,
 T_REAL,
 T_FEATURE,
 T_FEATURE_SUM,
 T_INSTANCES_CLUSTER_K>         &aochromfeac_toMutate,
 INPUT_ITERATOR                 aiiterator_instfirst,
 const INPUT_ITERATOR           aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>   &aifunc2p_dist
 )
{
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::MO2";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << " K = " << aochromfeac_toMutate.getNumClusterK()
	      << " fitness = " <<  aochromfeac_toMutate.getFitness()
	      << "\n(output gaencode::ChromosomeFEAC<>&: aochromfeac_toMutate["
	      << geverboseui_idproc << ':'<<  &aochromfeac_toMutate << ']'
      	      << "\n input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      <<  "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  
  mat::MatrixWithRowNull<T_FEATURE> &lmatrixwrownull_genotypeCentroids =  
    aochromfeac_toMutate.getCentroids();
  std::vector<T_INSTANCES_CLUSTER_K>& lvectorit_genotypeNumInstClusterK =
    aochromfeac_toMutate.getNumInstancesClusterK();
  std::vector<T_REAL>& lvectorrt_genotypePartialFcC =
    aochromfeac_toMutate.getPartialFcC();

  /*1. Randomly generate a number n \in {1,..., kg};
   */
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "MO2 1. RANDOMLY GENERATE A NUMBER n in {1,...,k_g}:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << std::endl; 
  }
#endif /*__VERBOSE_YES*/
  
  std::uniform_int_distribution<uintidx> uniformdis_ui1N
    (1, lmatrixwrownull_genotypeCentroids.getNumRows());

  uintidx luintidx_n = uniformdis_ui1N(gmt19937_eng); 
    

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "MO2 1. RANDOMLY GENERATE A NUMBER n in {1,...,k_g}: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\nn = " << luintidx_n << " in {1,...,"
	      << lmatrixwrownull_genotypeCentroids.getNumRows() << '}'
	      << std::endl; 
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
 

  /*2. For i = (1,...,n) do:
   */
  for (uintidx lui_i = 1; lui_i <= luintidx_n; lui_i++) {

    
    /*2.1 Randomly choose a cluster Cs encoded into g;
     */
    

#if defined(ALG_EAC_I_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006) || \
  defined(ALG_EAC_III_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006)

    /*\cite{Alves:etal:GAclusteringLabelKVar:FEAC:2006}
      EAC-I
      EAC-III
    */

    /*Alternatively, it is possible to use information concerning
      the quality of the clusters to select them for mutation. Our
      first proposal to enhance the EAC goes into this direction,
    
      basically changing the way in which the EAC mutation
      operators are applied (step 6 in Fig. 5). Particularly, we
      hypothesize that the better the cluster the smaller should be
      its probability of being mutated.Thus, good clusters tend to
      be maintained during the evolutionary process, whereas bad
      clusters are more likely to be mutated, hopefully towards
      improving the partition encoded into the genotype.

      In addition, let us consider that f_c(C_i) can be normalized
      within the interval [0,1], which is precisely the case when
      using the fitness functions to be described in Section IV.
    */
    
    std::vector<T_REAL> lvectorT_partialFcCNorma01(lvectorrt_genotypePartialFcC.size(),1);
    /*1 - fc(Ci)
     */
    interfacesse::axpy
      (lvectorT_partialFcCNorma01.data(),
       -1.0,
       lvectorrt_genotypePartialFcC.data(),
       (uintidx) lvectorT_partialFcCNorma01.size()
       );
    /* / \sum(1-fc(Ci))
     */
    T_REAL lrt_sum1_fc =
      interfacesse::sum
      (lvectorT_partialFcCNorma01.data(),
       (uintidx) lvectorT_partialFcCNorma01.size()
       );
    lrt_sum1_fc = 1.0 / lrt_sum1_fc;
  
    interfacesse::scal
      (lvectorT_partialFcCNorma01.data(),
       lrt_sum1_fc,
       (uintidx) lvectorT_partialFcCNorma01.size()
       );

    std::vector<T_REAL>&& lvectorrt_probDistRouletteWheel =
      prob::makeDistRouletteWheel
      (lvectorT_partialFcCNorma01.begin(),
       lvectorT_partialFcCNorma01.end(),
       [](const T_REAL& lt_metrici) -> T_REAL
       {
	 return lt_metrici;
       }
       );
   
    T_CLUSTERIDX lcidx_Cs = 
      gaselect::getIdxRouletteWheel
      (lvectorrt_probDistRouletteWheel,
       T_CLUSTERIDX(0)
       );
      
    /*ALG_EAC_I_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
      ALG_EAC_III_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
    */
 
#elif defined(ALG_EAC_II_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006) || \
  defined(ALG_FEAC_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006)

    /*\cite{Alves:etal:GAclusteringLabelKVar:FEAC:2006}
      EAC-II
      F-EAC
    */

    std::vector<T_REAL>&& lvectorT_partialFcCLinearNorm =
      prob::linearNormalization
      (lvectorrt_genotypePartialFcC.begin(),
       lvectorrt_genotypePartialFcC.end(),
       [](const T_REAL& lt_metrici) -> T_REAL
       {
	 return lt_metrici;
       }
       );

    std::vector<T_REAL>&& lvectorrt_probDistRouletteWheel =
      prob::makeDistRouletteWheel
      (lvectorT_partialFcCLinearNorm.begin(),
       lvectorT_partialFcCLinearNorm.end(),
       [](const T_REAL& lt_metrici) -> T_REAL
       {
	 return lt_metrici;
       }
       );
  
    T_CLUSTERIDX lcidx_Cs = 
      gaselect::getIdxRouletteWheel
      (lvectorrt_probDistRouletteWheel,
       T_CLUSTERIDX(0)
       );
    
    /*END ALG_EAC_II_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
      ALG_FEAC_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006
    */
#else
    /* ELSE
      ALG_EAC_CLUSTERINGLABELVARK_HRUSCHKA_CAMPELLO_CASTRO2006 || OTHER
      \cite{Hruschka:Campello:Castro:GAClusteringLabelKVar:EAC:2006}
      EAC
    */
    
    std::uniform_int_distribution<uintidx> uniformdis_0K
      (0,(uintidx) aochromfeac_toMutate.getNumClusterK()-1);

    T_CLUSTERIDX lcidx_Cs = uniformdis_0K(gmt19937_eng);


#endif  /* END  #if defined(ALG_EAC_I_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006) || \
	   defined(ALG_EAC_III_CLUSTERINGLABELVARK_ALVES_CAMPELLO_HRUSCHKA2006)
	*/

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "MO2 2.1 RANDOMLY CHOOSE A CLUSTER C_s in g,  IN"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl; 
    }
#endif /*__VERBOSE_YES*/


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "MO2 2.1 RANDOMLY CHOOSE A CLUSTER C_s in g, OUT"
		<< '(' << geiinparam_verbose << ')'
		<< "\nCs = " << lcidx_Cs 
		<< std::endl; 
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    /*Cs must have more than two objects to be
      eligible for this operator 
    */
    if  ( lvectorit_genotypeNumInstClusterK.at(lcidx_Cs) > 2 ) {
      
      T_CLUSTERIDX *larraycidx_memberShip =
	aochromfeac_toMutate.getString();
      
      uintidx luintidx_S1 = 0;
      
      {/* BEGIN *2.2.1 Randomly choose an object s1 \in Cs;
	*/
	
#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "MO2 2.2.1 RANDOMLY CHOOSE A OBJECT s_1 in C_s,  IN"
		    << '(' << geiinparam_verbose << ')' 
		    << std::endl; 
	}
#endif /*__VERBOSE_YES*/
	
	
	std::uniform_int_distribution<uintidx> uniformdis_consecutiveS1
	  (0, uintidx(lvectorit_genotypeNumInstClusterK.at(lcidx_Cs)-1));
	
	uintidx luintidx_ConsecutiveS1 =  /*Index object S1*/
	  uniformdis_consecutiveS1(gmt19937_eng);
	  
	uintidx luintidx_l = 0;
	do { 
	  if ( larraycidx_memberShip[luintidx_l] == lcidx_Cs )
	    ++luintidx_S1;
	} while ( (luintidx_S1  < luintidx_ConsecutiveS1) &&
		  (++luintidx_l <  lui_numInstances ) // aivectorptinst_instances.size())
		  );
	luintidx_S1 = luintidx_l;
	
	
#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "MO2 2.2.1 RANDOMLY CHOOSE A OBJECT s_1 in C_s, OUT"
		    << '(' << geiinparam_verbose << ')'
		    << "\nluintidx_ConsecutiveS1 = " << luintidx_ConsecutiveS1 
		    << "\ts_1 = " << luintidx_S1;
	  std::cout << std::endl; 
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
	
      } /* END *2.2.1 Randomly choose an object s1 \in Cs;
	 */
          
      /*2.2.2 Determine the farthest object s2 \in Cs from 
	s1;
      */
      uintidx luintidx_S2 = 
	nearest::farthestInstanceFromS1
	(luintidx_S1,      
	 aochromfeac_toMutate.getString(),
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 aifunc2p_dist
	 );
      
      /*2.2.3 Generate two new clusters Cs’ and Cs’’,
	placing the object(s) of Cs closer to s1 into
	Cs’ and the object(s) closer to s2 into Cs’’;
      */
      mat::MatrixRow<T_FEATURE>  lmatrixrowt_centroidsNewS
	(2,data::Instance<T_FEATURE>::getNumDimensions());  /*New centroids*/

      mat::MatrixRow<T_FEATURE_SUM>       
	llmatrixrowt_sumInstancesCluster
	(2,data::Instance<T_FEATURE>::getNumDimensions());
	
      
      lmatrixrowt_centroidsNewS.copyRow
	((uintidx) 0,
	 (*std::next(aiiterator_instfirst,luintidx_S1))->getFeatures()
	 );
      lmatrixrowt_centroidsNewS.copyRow
	((uintidx) 1,
	 (*std::next(aiiterator_instfirst,luintidx_S2))->getFeatures()
	 );

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {	
	std::ostringstream lostrstream_labelSubStep;
	lostrstream_labelSubStep
	  << "<CENTROIDS:"
	  <<  geverbosepc_labelstep << ':' << lpc_labelFunc 
	  << ":S1_S2[" << geverboseui_idproc << ':' << &aochromfeac_toMutate << "]:";
	lmatrixrowt_centroidsNewS.print
	  (std::cout,
	   lostrstream_labelSubStep.str().c_str(),
	   ',',
	   ';'
	   );
	std::cout << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  				 
      std::vector<T_INSTANCES_CLUSTER_K> 
	lvectort_numInstInClusterNew(2,0);

      std::vector<T_CLUSTERIDX>
	lvectorcidx_idxClustersNew;
      lvectorcidx_idxClustersNew.reserve(2);
      
      lvectorcidx_idxClustersNew.push_back( lcidx_Cs );
      lvectorcidx_idxClustersNew.push_back( aochromfeac_toMutate.getNumClusterK() );  

      mat::MatrixRow<T_FEATURE_SUM>       
	llmatrixt_sumInstancesCluster
	(lmatrixrowt_centroidsNewS.getNumRows(), 
	 lmatrixrowt_centroidsNewS.getNumColumns()
	 );
      
      clusteringop::updateClusterCj
	(lmatrixrowt_centroidsNewS, //Centroids average of each cluster
	 llmatrixt_sumInstancesCluster,
	 lvectort_numInstInClusterNew,
	 aochromfeac_toMutate.getString(),
	 lcidx_Cs,
	 lvectorcidx_idxClustersNew,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 aifunc2p_dist
	 );

      lmatrixwrownull_genotypeCentroids.copyRow
	(lcidx_Cs,lmatrixrowt_centroidsNewS.getRow(0));
      lmatrixwrownull_genotypeCentroids.addRow
	(lmatrixrowt_centroidsNewS.getRow(1));


#ifdef __VERBOSE_YES
     if ( lvectorit_genotypeNumInstClusterK[lcidx_Cs] !=
	  (lvectort_numInstInClusterNew[0] + lvectort_numInstInClusterNew[1]) )
       {
	    std::cout
	      << "\nERROR:"
	      << geverbosepc_labelstep
	      <<":number of cluster not valid, genotypeNumInstClusterK "
	      << lvectorit_genotypeNumInstClusterK[lcidx_Cs]
	      << " = numInstInClusterNew[0] + numInstInClusterNew[1] "
	      << lvectort_numInstInClusterNew[0] + lvectort_numInstInClusterNew[1];
      }
#endif /*__VERBOSE_YES*/
           
      lvectorit_genotypeNumInstClusterK[lcidx_Cs] = lvectort_numInstInClusterNew[0];
      lvectorit_genotypeNumInstClusterK.push_back(lvectort_numInstInClusterNew[1]);
      lvectorrt_genotypePartialFcC.push_back(-1.0);


#ifdef __FITNESS_RAND_INDEX__

      partition::PartitionLabel
	<T_CLUSTERIDX>
	lpartition_clusters
	(aochromfeac_toMutate.getString(),
	 aochromfeac_toMutate.getStringSize(),
	 aochromfeac_toMutate.getNumClusterK() 
	 );
      
      sm::ConfusionMatchingMatrix<T_INSTANCES_CLUSTER_K>&&
	lmatchmatrix_confusion = 
	sm::getConfusionMatrix
	(aiiterator_instfirst,
	 aiiterator_instlast,
	 lpartition_clusters,
	 [](const data::Instance<T_FEATURE>* aiinst_iter )
	 -> T_INSTANCES_CLUSTER_K
	 {
	   return  T_INSTANCES_CLUSTER_K(1);
	 },
	 [](const data::Instance<T_FEATURE>* aiinst_iter )
	 -> T_CLUSTERIDX
	 {
	   data::InstanceClass
	     <T_FEATURE,
	      T_INSTANCES_CLUSTER_K,
	      T_CLUSTERIDX>
	     *linstclass_iter = 
	     (data::InstanceClass
	      <T_FEATURE,
	      T_INSTANCES_CLUSTER_K,
	      T_CLUSTERIDX>*)
	     aiinst_iter;
	       
	   return linstclass_iter->getClassIdx();
	       
	 }
	 );

      std::vector<T_REAL>&&  lvectort_partialFitness =
	sm::partialRandIndex<T_REAL,T_INSTANCES_CLUSTER_K>(lmatchmatrix_confusion);

/*#endif __FITNESS_RAND_INDEX__*/      

#else /* __FITNESS_SIMPLIFIED_SILHOUETTE__ */
      
      partition::PartitionLabel<T_CLUSTERIDX>
	lpartition_clustersLabel
	(aochromfeac_toMutate.getString(),
	 aochromfeac_toMutate.getStringSize(),
	 aochromfeac_toMutate.getNumClusterK()
	 );

      std::vector<T_REAL>&&  lvectort_partialFitness =
	um::simplifiedSilhouette
	(aochromfeac_toMutate.getCentroids(),
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 lpartition_clustersLabel,
	 aochromfeac_toMutate.getNumInstancesClusterK(),
	 aifunc2p_dist
	 );

#endif /*__FITNESS_SIMPLIFIED_SILHOUETTE__*/
	   
      T_REAL lmetrict_partialSilhouette = 
	interfacesse::sum
	(lvectort_partialFitness.data(),
	 (uintidx) lvectort_partialFitness.size()
	 );

      lmetrict_partialSilhouette /= (T_REAL) lvectort_partialFitness.size();

      aochromfeac_toMutate.setPartialFcC(lvectort_partialFitness);
      aochromfeac_toMutate.setObjetiveFunc(lmetrict_partialSilhouette);
      
#ifdef __VERBOSE_YES
      const char* lpc_labelFuncStep = "2. FOR END";
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {	
	std::ostringstream lostrstream_labelFunStep;
	lostrstream_labelFunStep <<  geverbosepc_labelstep
				 << ':' << lpc_labelFunc
				 << ':' <<  lpc_labelFuncStep; 
	aochromfeac_toMutate.print
	  (std::cout,
	   lostrstream_labelFunStep.str().c_str(),
	   ',',
	   ';'
	   );
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/ 
      

    } /*Else do not split Cs.*/

 
  } /*End For*/
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " K = " << aochromfeac_toMutate.getNumClusterK() << '\n'; 
    std::ostringstream lostrstream_labelFunc;
    lostrstream_labelFunc
      <<  geverbosepc_labelstep
      << ':' << lpc_labelFunc;
    aochromfeac_toMutate.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

} /*END namespace gaclusteringop*/

#endif /*__GA_CLUSTERING_OPERATOR_HPP*/

  
