/*! \file clustering_operator_centroids.hpp
 *
 * \brief clustering operator centroids  
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __CLUSTERING_OPERATOR_CENTROIDS_HPP
#define __CLUSTERING_OPERATOR_CENTROIDS_HPP

#include <iostream>
#include <stdio.h>
#include <vector>
#include <utility>      // std::move
#include <stdexcept>
#include <unordered_set>

#include "random_ext.hpp"
#include "instance_frequency.hpp"
#include "matrix.hpp"
#include "nearestinstance_operator.hpp"
#include "linear_algebra_level1.hpp"
#include "partition_linked_stats.hpp"
#include "partition_disjsets.hpp"
#include "bit_array.hpp"
#include "stats_instances.hpp"
#include "probability_selection.hpp"
#include "nearestcentroids_operator.hpp"
#include "partition_label.hpp"

#include "verbose_global.hpp"

extern StdMT19937  gmt19937_eng;

/*! \namespace clusteringop
  \brief Clustering operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace clusteringop {
  

/*! \fn void meanCentroids(T_CLUSTERIDX &aocidx_numClusterNull, MatrixBase<T_CENTROIDS> &aomatrixt_meanCentroids, MatrixBase<T_FEATURE_SUM>  &aimatrixtrowt_sumInstancesClusterK, std::vector<T_INSTANCES_CLUSTER_K>  &aovectort_numInstClusterK)
    \brief Mean of the centroids
    \details Calculate the centroids with the sum of the instances and their number per cluster
    \param aocidx_numClusterNull a integer with number of null cluster 
    \param aomatrixt_meanCentroids a MatrixBase<T_CENTROIDS>
    \param aimatrixtrowt_sumInstancesClusterK a MatrixBase<T_FEATURE_SUM>
    \param aovectort_numInstClusterK a std::vector<T_INSTANCES_CLUSTER_K>
 */
template < typename T_CENTROIDS, 
	   typename T_FEATURE_SUM, 
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
void
meanCentroids
(T_CLUSTERIDX                        &aocidx_numClusterNull,
 mat::MatrixBase<T_CENTROIDS>        &aomatrixt_meanCentroids, 
 mat::MatrixBase<T_FEATURE_SUM>      &aimatrixtrowt_sumInstancesClusterK, 
 std::vector<T_INSTANCES_CLUSTER_K>  &aovectort_numInstClusterK
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::meanCentroids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output mat::MatrixRow<T_CENTROIDS>: aomatrixt_meanCentroids["  
	      << &aomatrixt_meanCentroids << "]\n"
	      << " input  mat::MatrixRow<T_FEATURE_SUM>: aimatrixtrowt_sumInstancesClusterK[" 
	      << &aimatrixtrowt_sumInstancesClusterK << "]\n"
	      << " input  std::vector<T_INSTANCES_CLUSTER_K>: aovectort_numInstClusterK["  
	      << &aovectort_numInstClusterK << "]\n"
	      << ")\n";
  }
#endif //__VERBOSE_YES

  aocidx_numClusterNull = 0;

  for ( uintidx lui_i = 0; lui_i < aomatrixt_meanCentroids.getNumRows(); lui_i++) { 
    if ( aovectort_numInstClusterK[lui_i] == 0) {
      ++aocidx_numClusterNull;
    }
    else {
      stats::meanVector
	(aomatrixt_meanCentroids.getRow(lui_i),
	 aovectort_numInstClusterK[lui_i],
	 aimatrixtrowt_sumInstancesClusterK.getRow(lui_i),
	 lui_i
	 );
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "aocidx_numClusterNull = " << aocidx_numClusterNull << '\n';
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS MEAN:" << lpc_labelFunc;
    aomatrixt_meanCentroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    if ( aocidx_numClusterNull > 0 ) {
	std::cout << "\nERROR: solution is not valid, cluster without instances";
    }
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
}


/*! \fn void initialize(mat::MatrixBase<T_FEATURE> &aomatrixt_centroids, const INPUT_ITERATOR aiiterator_instfirst, INPUT_ITERATOR iterator_idxInstanceRandFirst)
    \brief Centroids initialized
    \details Initialize a centroids randomly or based on some prior knowledge
    \param aomatrixt_centroids a mat::MatrixBase with the centroids of each cluster
    \param aiiterator_instfirst a iterator with the instances
    \param iterator_idxInstanceRandFirst a iterator with index of the instances to copy
 */
template < typename T_FEATURE,
	   typename INPUT_ITERATOR_IDX,
	   typename INPUT_ITERATOR
	   > 
void
initialize
(mat::MatrixBase<T_FEATURE> &aomatrixt_centroids,
 const INPUT_ITERATOR       aiiterator_instfirst,
 INPUT_ITERATOR_IDX         aiiterator_idxInstanceRandFirst
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::initialize";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output mat::MatrixRow<T_FEATURE>&: aomatrixt_centroids["
	      << &aomatrixt_centroids << "]\n"
	      << "input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << "input aiiterator_idxInstanceRandFirst[" << *aiiterator_idxInstanceRandFirst << "]\n"
	      <<  ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  /*Copy rows*/
  for (uintidx li_j = 0; li_j < aomatrixt_centroids.getNumRows(); li_j++) { 
    T_FEATURE* lmrT_centroids = aomatrixt_centroids.getRow(li_j);

    data::Instance<T_FEATURE>* liter_iInstance =
      *std::next(aiiterator_instfirst,*aiiterator_idxInstanceRandFirst);
    
    interfacesse::copy
      (lmrT_centroids,
       liter_iInstance->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    
    ++aiiterator_idxInstanceRandFirst;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDSCLUSTER:"<< lpc_labelFunc;
    
    aomatrixt_centroids.print
      (std::cout,
       lostrstream_labelCentroids.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  

} /*initialized*/
  
  
/*! \fn void randomInitialize(mat::MatrixBase<T_FEATURE> &aomatrixt_centroids, const INPUT_ITERATOR aiiterator_instfirst, INPUT_ITERATOR iterator_idxInstanceRandFirst, const INPUT_ITERATOR aiiterator_instlast)
    \brief Initialize the centroids by selecting random instances.
    \details 
    \param aomatrixt_centroids a mat::MatrixBase with the centroids of each cluster
    \param aiiterator_instfirst a iterator with the instances
    \param aiiterator_instlast a const input iterator of the instances
 */
template < typename T_FEATURE,
	   typename INPUT_ITERATOR
	   > 
void
randomInitialize
(mat::MatrixBase<T_FEATURE> &aomatrixt_centroids,
 const INPUT_ITERATOR       aiiterator_instfirst,
 const INPUT_ITERATOR       aiiterator_instlast
 )
{

  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::randomInitialize";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output mat::MatrixRow<T_FEATURE>&: aomatrixt_centroids["
	      << &aomatrixt_centroids << "]\n"
              <<  aomatrixt_centroids << '\n'
	      << "input const aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << "input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "lui_numInstances =  "  << lui_numInstances << '\n'
	      <<  ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  /*Copy rows*/
  std::uniform_int_distribution<uintidx> luniformdis_uiidxInstances0n
    (0,lui_numInstances-1);
  
  std::unordered_set<uintidx>&& lunorderedset_idxRandInstances =
    prob::getWithoutRepeatsSet
    ( aomatrixt_centroids.getNumRows()
      ,[&]() -> uintidx
      { 
	return luniformdis_uiidxInstances0n(gmt19937_eng);
      }
      );

  auto literator_idxInstanceRandFirst = lunorderedset_idxRandInstances.begin();
  
  for (uintidx li_j = 0; li_j < aomatrixt_centroids.getNumRows(); li_j++) { 
    T_FEATURE* lmrT_centroids = aomatrixt_centroids.getRow(li_j);

    data::Instance<T_FEATURE>* liter_iInstance =
      *std::next(aiiterator_instfirst,*literator_idxInstanceRandFirst);
    
    interfacesse::copy
      (lmrT_centroids,
       liter_iInstance->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    
    ++literator_idxInstanceRandFirst;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDSCLUSTER:"<< lpc_labelFunc;
    
    aomatrixt_centroids.print
      (std::cout,
       lostrstream_labelCentroids.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  

} /*initialize*/



/*! \fn void randomInitialize(mat::MatrixBase<T_FEATURE> &aomatrixt_centroids, T_CLUSTERIDX *aoarraycidx_memberShip, std::vector<T_INSTANCES_CLUSTER_K> &aovectorit_numInstClusterK, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast)
    \brief Initialize the centroids and labels by selecting random instances.
    \details 
    \param aomatrixt_centroids a mat::MatrixBase with the centroids of each cluster
    \param aoarraycidx_memberShip a T_CLUSTERIDX
    \param aovectorit_numInstClusterK a std::vector
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a  const input iterator of the instances
 */
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX,
	   typename INPUT_ITERATOR
	   > 
void
randomInitialize
(mat::MatrixBase<T_FEATURE>         &aomatrixt_centroids,
 T_CLUSTERIDX                       *aoarraycidx_memberShip,
 std::vector<T_INSTANCES_CLUSTER_K> &aovectorit_numInstClusterK,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast
 )
{
const uintidx  lui_numInstances =
  uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::randomInitialize";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output mat::MatrixBase<T_FEATURE>&: aomatrixt_centroids["
	      << &aomatrixt_centroids << ']'
              << "\n output std::vector<T_INSTANCES_CLUSTER_K>&: aovectorit_numInstClusterK[" 
	      << &aovectorit_numInstClusterK << ']'
              << "\n output T_CLUSTERIDX*: aoarraycidx_memberShip[" 
	      << aoarraycidx_memberShip << "]\n"
	      << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
              << "lui_numInstances = " << lui_numInstances 
	      <<  "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::uniform_int_distribution<uintidx> luniformdis_uiidxInstances0n
    (0,lui_numInstances-1);

  std::unordered_set<uintidx>&& lunorderedset_idxRandInstances =
    prob::getWithoutRepeatsSet
    ( aomatrixt_centroids.getNumRows()
      ,[&]() -> uintidx
      {
	return luniformdis_uiidxInstances0n(gmt19937_eng);
      }
      );

  auto literator_idxInstanceRandFirst = lunorderedset_idxRandInstances.begin();
  
  interfacesse::copya
    (aoarraycidx_memberShip,
     T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN),
     lui_numInstances
     );

  /*Copy rows*/
  for (uintidx li_j = 0; li_j < aomatrixt_centroids.getNumRows(); li_j++) { 
    T_FEATURE* larrayrowt_centroids = aomatrixt_centroids.getRow(li_j);
    data::Instance<T_FEATURE>* liter_iInstance =
      *std::next(aiiterator_instfirst,*literator_idxInstanceRandFirst);
    interfacesse::copy 
      (larrayrowt_centroids,
       liter_iInstance->getFeatures(),
       liter_iInstance->getNumDimensions()
       );
    aoarraycidx_memberShip[*literator_idxInstanceRandFirst] =
      (T_CLUSTERIDX) li_j;
    aovectorit_numInstClusterK[li_j] = 1;
    ++literator_idxInstanceRandFirst;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aomatrixt_centroids.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:"
				<< lpc_labelFunc;
    inout::containerprint
      (aoarraycidx_memberShip,
       aoarraycidx_memberShip+lui_numInstances,
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelNumInst;
    lostrstream_labelNumInst << "<NUMINSTANCES:"
			     << lpc_labelFunc;
 
    inout::containerprint
      (aovectorit_numInstClusterK.begin(),
       aovectorit_numInstClusterK.end(),
       std::cout,
       lostrstream_labelNumInst.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  

}

  
/*! \fn  uintidx reassignCluster(T_CLUSTERIDX *aioarraycidx_memberShip, const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist) 
    \brief Change each instance to the nearest cluster
    \details Each data instance of the container is assigned to the cluster with the nearest neighbor rule
    \param aioarraycidx_memberShip an array of indexes belonging to previously assigned or UNKNOWN_CLUSTER_IDX. Also returns how many cluster changed
    \param aimatrixt_centroids a mat::MatrixRow with the centroids of each cluster
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a const input iterator of the instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template < typename T_FEATURE,
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_DIST, 
	   typename INPUT_ITERATOR
	   >
uintidx 
reassignCluster
(T_CLUSTERIDX                       *aioarraycidx_memberShip,
 const mat::MatrixRow<T_FEATURE>    &aimatrixt_centroids,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist
 ) 
{
  uintidx      louintidx_threshold = 0; 
  T_DIST       lT_distMinCentInst;
  
#ifdef __VERBOSE_YES
  const uintidx  lui_numInstances =
      uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  const char* lpc_labelFunc = "clusteringop::reassignCluster";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output T_CLUSTERIDX*: aioarraycidx_memberShip[" 
	      << aioarraycidx_memberShip << "]\n"
	      << "\t input  mat::MatrixRow<T_CENTROIDS>&: aimatrixt_centroids["
	      << &aimatrixt_centroids << "]\n"
	      << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "\t input  aifunc2p_dista\n"
	      << "\t)\n";
  }
#endif //__VERBOSE_YES 

  for (; aiiterator_instfirst != aiiterator_instlast;
       aiiterator_instfirst++, aioarraycidx_memberShip++)
    {
    data::Instance<T_FEATURE>* linst_inter =
      (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    T_CLUSTERIDX lmgidx_j = 
      nearest::NN
      <T_CLUSTERIDX,T_FEATURE,T_DIST>
      (lT_distMinCentInst,
       aimatrixt_centroids,
       linst_inter->getFeatures(),
       aifunc2p_dist
       );
    if ( *aioarraycidx_memberShip != lmgidx_j ) {
      *aioarraycidx_memberShip = lmgidx_j;
      ++louintidx_threshold;
    }
  } //END FOR

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::cout << "louintidx_threshold = " << louintidx_threshold 
	      << '\n';
    
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:" << lpc_labelFunc;
    inout::containerprint
      (aioarraycidx_memberShip,
       aioarraycidx_memberShip + lui_numInstances,
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout << std::endl;
 
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return louintidx_threshold; 
} /*reassignCluster*/


/*! \fn void setUpCuster(T_CLUSTERIDX *aioarraycidx_memberShip, const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)  
    \brief  Set up cluster index
    \details Assigns the instances of the container with an unknown cluster using the nearest neighbor rule
    \param aioarraycidx_memberShip an array of indexes only with the seed instances as the clustering center. The other instances with UNKNOWN_CLUSTER_IDX
    \param aimatrixt_centroids a mat::MatrixRow with the centroids of each cluster
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a const input iterator of the instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template < typename T_FEATURE,
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
void
setUpCuster
(T_CLUSTERIDX                       *aioarraycidx_memberShip,
 const mat::MatrixRow<T_FEATURE>    &aimatrixt_centroids,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist
 ) 
{ 
  T_DIST       lT_distMinCentInst;
  
#ifdef __VERBOSE_YES
  const uintidx  lui_numInstances =
      uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  const char* lpc_labelFunc = "clusteringop::setUpCuster";
  const T_CLUSTERIDX  *aioarraycidx_memberShipVerbose = aioarraycidx_memberShip;
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output T_CLUSTERIDX*: aioarraycidx_memberShip[" 
	      << aioarraycidx_memberShipVerbose << "]\n";
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:" << lpc_labelFunc;
    inout::containerprint
      (aioarraycidx_memberShipVerbose,
       aioarraycidx_memberShipVerbose + lui_numInstances,
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout  << "\n input  mat::MatrixRow<T_CENTROIDS>&: aimatrixt_centroids["
	      << &aimatrixt_centroids << "]\n"
              << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "input  aifunc2p_dista\n"
	      << ")\n";
  }
#endif //__VERBOSE_YES 

  for ( ; aiiterator_instfirst != aiiterator_instlast;
       aiiterator_instfirst++, aioarraycidx_memberShip++)
    {
    if (*aioarraycidx_memberShip == NEARESTCENTROID_UNKNOWN ) {
      data::Instance<T_FEATURE>* linst_inter =
	(data::Instance<T_FEATURE>*) *aiiterator_instfirst;
      T_CLUSTERIDX lmgidx_j = 
	nearest::NN
	<T_CLUSTERIDX,T_FEATURE,T_DIST>
	(lT_distMinCentInst,
	 aimatrixt_centroids,
	 linst_inter->getFeatures(),
	 aifunc2p_dist
	 );
      *aioarraycidx_memberShip = lmgidx_j;
    }
  } //End for

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:" << lpc_labelFunc;
    inout::containerprint
      (aioarraycidx_memberShipVerbose,
       aioarraycidx_memberShipVerbose + lui_numInstances,
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout << std::endl;
  
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
 
} /*setUpCuster*/

  
/*! \fn T_CLUSTERIDX  getCentroids(const mat::MatrixRow<T_FEATURE>  &aomatrixt_centroids, mat::MatrixRow<T_FEATURE_SUM> &aomatrixt_sumInstancesCluster, std::vector<T_INSTANCES_CLUSTER_K>  &aovectort_numInstancesInClusterK, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast)
    \brief Calculate centroids
    \details Calculate centroids, instances sum and instances number for a for a given instance partition
    \param aomatrixt_centroids a mat::MatrixRow with centroids of each cluster
    \param aomatrixt_sumInstancesCluster a mat::MatrixRow with the sum of instances per cluster
    \param aovectort_numInstancesInClusterK a std::vector with the number of instances per cluster
    \param aipartition_clusters a partition::Partition of instances in clusters
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a  const input iterator of the instances
 */
template < typename T_FEATURE, 
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,  //-1, 0, 1, .., N
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename INPUT_ITERATOR
	   >
T_CLUSTERIDX  
getCentroids
(mat::MatrixRow<T_FEATURE>          &aomatrixt_centroids,
 mat::MatrixRow<T_FEATURE_SUM>      &aomatrixt_sumInstancesCluster,
 std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstancesInClusterK,
 partition::Partition<T_CLUSTERIDX> &aipartition_clusters,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast
 ) 
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::getCentroids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input mat::MatrixRow<T_FEATURE>  aomatrixt_centroids"
	      <<  aomatrixt_centroids.getNumRows() << "x"
	      <<  aomatrixt_centroids.getNumColumns()
	      << "[" <<  &aomatrixt_centroids << "]\n"
	      << "MatrixRow<T_FEATURE_SUM> aomatrixt_sumInstancesCluster"
	      <<  aomatrixt_sumInstancesCluster.getNumRows() << "x"
	      <<  aomatrixt_sumInstancesCluster.getNumColumns()
	      << "[" <<  &aomatrixt_sumInstancesCluster << "]\n"
	      << "(input partition::Partition<>: aipartition_clusters[" 
	      << &aipartition_clusters << "]\tk = "
	      << aipartition_clusters.getNumCluster() << "\n"
	      << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES 

  const T_CLUSTERIDX  lcidx_numClusterK = aipartition_clusters.getNumCluster();  
  const T_FEATURE lT_alpha = T_FEATURE(1);

  aomatrixt_sumInstancesCluster.initialize(T_FEATURE(0.0));
  
  interfacesse::copya
    (aovectort_numInstancesInClusterK.data(),
     T_INSTANCES_CLUSTER_K(0),
     aomatrixt_sumInstancesCluster.getNumRows()
     );
  
  for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    
    T_CLUSTERIDX lcidx_instInClusterJ = aipartition_clusters.next();
    
    if ( 0 <= lcidx_instInClusterJ  && lcidx_instInClusterJ <  lcidx_numClusterK  ) {
      const T_FEATURE* linst_inter =
      ((data::Instance<T_FEATURE>*) *aiiterator_instfirst)->getFeatures();
      T_FEATURE_SUM* larrarrowt_sumInstancesCluster = 
	aomatrixt_sumInstancesCluster.getRow(lcidx_instInClusterJ);
      interfacesse::axpy
	(larrarrowt_sumInstancesCluster,
	 lT_alpha,
	 linst_inter,
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
      
      aovectort_numInstancesInClusterK[lcidx_instInClusterJ]++;
      
    }
  } //End for
  
  T_CLUSTERIDX locidx_numClusterNull = 0;
  
  for ( uintidx lui_i = 0; lui_i < aomatrixt_centroids.getNumRows(); lui_i++) { 
    if ( aovectort_numInstancesInClusterK[lui_i] == 0) {
      ++locidx_numClusterNull;
    }
    else {
      stats::meanVector
	(aomatrixt_centroids.getRow(lui_i),
	 aovectort_numInstancesInClusterK[lui_i],
	 aomatrixt_sumInstancesCluster.getRow(lui_i),
	 lui_i
	 );
    }
  }


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << "\n\tlocidx_numClusterNull = " << locidx_numClusterNull << '\n';

    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS:" << lpc_labelFunc;
      aomatrixt_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    std::cout << std::endl;

    std::ostringstream lostrstream_labelSumInst;
    lostrstream_labelSumInst << "<INSTANCES SUM:" << lpc_labelFunc;
   
    aomatrixt_sumInstancesCluster.print(std::cout,lostrstream_labelSumInst.str().c_str(),',',';');
    std::cout << std::endl;

    std::ostringstream lostrstream_labelInstClusterK;
    lostrstream_labelInstClusterK << "<INSTANCES CLUSTER K: " << lpc_labelFunc
				  << ":aovectort_numInstancesInClusterK["
				  << &aovectort_numInstancesInClusterK << ']';
    inout::containerprint
      (aovectort_numInstancesInClusterK.begin(),
       aovectort_numInstancesInClusterK.end(),
       std::cout,
       lostrstream_labelInstClusterK.str().c_str(),
       ','
       );
    std::cout << std::endl;
   
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return locidx_numClusterNull; 
  
}


/*! \fn std::tuple<mat::MatrixRow<T_FEATURE>,partition::PartitionDisjSets<T_CLUSTERIDX>,std::vector<T_INSTANCES_CLUSTER_K> > getClusters(const mat::BitArray<T_BITSIZE> &aibirarray_clusterSeedIdx, const mat::MatrixRow<T_FEATURE> &aimatrixrowt_Vi, const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstBi, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
    \brief Calculate centroids \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
    \details generates a set of cluster centroids from other centroids used as seeds
    \param aibirarray_clusterSeedIdx a mat::BitArray 
    \param aimatrixrowt_Vi a mat::MatrixRow with centroid seeds
    \param aivectort_numInstBi a std::vector with the number of instances in the seed clusters
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
    \param aiit_datatypeClusterIdx to specify the datatype of the cluster labels
 */
template < typename T_BITSIZE,
	   typename T_FEATURE,
	   // typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_REAL  
	   >
std::tuple
<mat::MatrixRow<T_FEATURE>,
 partition::PartitionDisjSets<T_CLUSTERIDX>,
 std::vector<T_INSTANCES_CLUSTER_K> 
 >
getClusters
(const mat::BitArray<T_BITSIZE>            &aibirarray_clusterSeedIdx,
 const mat::MatrixRow<T_FEATURE>           &aimatrixrowt_Vi,
 const std::vector<T_INSTANCES_CLUSTER_K>  &aivectort_numInstBi,
 const dist::Dist<T_REAL,T_FEATURE>        &aifunc2p_dist,
 const T_CLUSTERIDX                        aiit_datatypeClusterIdx
 )
{
  mat::MatrixRow<T_FEATURE>       
    lomatrixrowt_S
    (aibirarray_clusterSeedIdx.getNumBitOn(),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
 
  ds::DisjSets lodisjsets_BkinCi( aibirarray_clusterSeedIdx.size() );
  std::vector<uintidx> lvectoruintidx_idxFirstInstVi(lomatrixrowt_S.getNumRows());  
  std::vector<T_INSTANCES_CLUSTER_K> lovectort_numInstCi( lomatrixrowt_S.getNumRows() );
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::getCentroids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n";
    
    std::ostringstream lostrstream_labelChrom;
    lostrstream_labelChrom << "INPUT:" << lpc_labelFunc;
    aibirarray_clusterSeedIdx.print(std::cout,lostrstream_labelChrom.str().c_str());
    std::cout << std::endl;

    std::cout << "mat::MatrixRow<T_FEATURE>: aimatrixrowt_Vi["
	      << &aimatrixrowt_Vi << ']' << " NumRows: " << aimatrixrowt_Vi.getNumRows() 
	      << "\ninput  std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstBi[" 
	      << &aivectort_numInstBi << ']' << " size: " << aivectort_numInstBi.size()
	      << "\ninput  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  if ( lomatrixrowt_S.getNumRows() > 0 ) {
    
    /*INITIALIZE: First Vi
     */
    uintidx luintidx_cj = 0;
    
    for ( uintidx luintidx_instii = 0;
	  luintidx_instii < aibirarray_clusterSeedIdx.size();
	  luintidx_instii++ )
      {
     
      if ( aibirarray_clusterSeedIdx.getBit(luintidx_instii) ) {

	lomatrixrowt_S.copyRow
	  (luintidx_cj,
	   aimatrixrowt_Vi.getRow(luintidx_instii)
	   );
	lovectort_numInstCi[luintidx_cj] = aivectort_numInstBi[luintidx_instii];
	lvectoruintidx_idxFirstInstVi.at(luintidx_cj) = luintidx_instii;
	++luintidx_cj;
	
      }
    }
	    
    /*Eqs(3) and Eqs(4)
     */
    for (uintidx luintidx_instii = 0;
	 luintidx_instii < aibirarray_clusterSeedIdx.size();
	 luintidx_instii++ )
      {
     
      if ( !aibirarray_clusterSeedIdx.getBit(luintidx_instii) ) {

	T_REAL  lnormt_distMinBiCj;
	uintidx luintidx_insticj =
	  nearest::NN
	  <uintidx,
	   T_FEATURE,
	   T_REAL
	   >
	  (lnormt_distMinBiCj,
	   lomatrixrowt_S,
	   aimatrixrowt_Vi.getRow(luintidx_instii),
	   aifunc2p_dist
	   );

	lodisjsets_BkinCi.merge(lvectoruintidx_idxFirstInstVi[luintidx_insticj],luintidx_instii);

	//#if  DATATYPE_CENTROIDS_ROUND == 0
  
	/*FOR INSTANCES NUMBER REAL DATATYPE 
	  T_FEATURE IS EQUAL TO T_FEATURE_SUM
	*/		
	/*	T_FEATURE lt_sumInstCiBi =
	  (T_FEATURE) (lovectort_numInstCi[luintidx_insticj] + aivectort_numInstBi[luintidx_instii]);

	T_FEATURE lt_alpha = 
	(T_FEATURE) (lovectort_numInstCi[luintidx_insticj]) / lt_sumInstCiBi; 		  */

	interfacesse::scal
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lovectort_numInstCi[luintidx_insticj],
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	interfacesse::axpy
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   aivectort_numInstBi[luintidx_instii],
	   aimatrixrowt_Vi.getRow(luintidx_instii),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	interfacesse::scalInv
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lovectort_numInstCi[luintidx_insticj]+aivectort_numInstBi[luintidx_instii],
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	/*OK
	interfacesse::scal
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lt_alpha,
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	lt_alpha = 
	  (T_FEATURE) (aivectort_numInstBi[luintidx_instii]) / lt_sumInstCiBi;
		
	interfacesse::axpy
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lt_alpha,
	   aimatrixrowt_Vi.getRow(luintidx_instii),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	ok*/

	//#else

	/*T_FEATURE lt_sumInstCiBi =
	  (lovectort_numInstCi[luintidx_insticj] + aivectort_numInstBi[luintidx_instii]);

	double lt_alpha = 
	(double) (lovectort_numInstCi[luintidx_insticj]) / lt_sumInstCiBi; 		  */

	/*	interfacesse::scal
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lovectort_numInstCi[luintidx_insticj],
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	interfacesse::axpy
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   aivectort_numInstBi[luintidx_instii],
	   aimatrixrowt_Vi.getRow(luintidx_instii),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	interfacesse::scalInv
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lovectort_numInstCi[luintidx_insticj]+aivectort_numInstBi[luintidx_instii],
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	*/

	/*	lt_alpha = 
	  (T_FEATURE) (aivectort_numInstBi[luintidx_instii]) / lt_sumInstCiBi;
	*/	
	  

	/*T_FEATURE_SUM *larrayt_Stmp = 
	  new T_FEATURE_SUM[data::Instance<T_FEATURE>::getNumDimensions()];
	*/
	/*Stmp = 0
	 */
	/*interfacesse::copya
	  (larrayt_Stmp, 
	   T_FEATURE_SUM(0),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	*/

	 /*Stmp = Stmp +  Sj*|Cj|
	 */
	/*interfacesse::axpy
	  (larrayt_Stmp,
	   lovectort_numInstCi[luintidx_insticj],
	   lomatrixrowt_S.getRow(luintidx_insticj),  
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	*/
	/*Stmp = Stmp + Vi*|Bi|
	 */
	/*interfacesse::axpy
	  (larrayt_Stmp,
	   lovectort_numInstCi[luintidx_insticj],
	   lomatrixrowt_S.getRow(luintidx_insticj),  
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	*/	
	/*	interfacesse::copya
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   T_FEATURE(0),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	DATATYPE_REAL lt_alpha = 1.0 /
	  ((DATATYPE_REAL)
	   (lovectort_numInstCi[luintidx_insticj] + aivectort_numInstBi[luintidx_instii]));

	interfacesse::axpy
	  (lomatrixrowt_S.getRow(luintidx_insticj),  
	   lt_alpha,
	   larrayt_Stmp,
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
		
	delete[] larrayt_Stmp;
	*/
		 
	//#endif /*DATATYPE_CENTROIDS_ROUND*/
	
	lovectort_numInstCi[luintidx_insticj] += aivectort_numInstBi[luintidx_instii];
      } /*if*/

    } /*for*/


  } /*if ( lomatrixrowt_S.getNumRows() > 0 ) */

  std::map<uintidx,T_CLUSTERIDX> lmapset_clusterk;
  std::vector<T_CLUSTERIDX>      lovectormmidx_memberBkinCi(lodisjsets_BkinCi.size());
      
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

    std::ostringstream lostrstream_labelSj;
    lostrstream_labelSj << "<CENTROIDS Sj:" << lpc_labelFunc;
    lomatrixrowt_S.print(std::cout,lostrstream_labelSj.str().c_str(),',',';');
    std::cout << std::endl;
	  
    std::ostringstream lostrstream_labelBkinCi;
    lostrstream_labelBkinCi
      <<':'  <<lpc_labelFunc << ":lodisjsets_BkinCi";
    lodisjsets_BkinCi.print
      (std::cout,
       lostrstream_labelBkinCi.str().c_str(),
       ',');
    std::cout << std::endl;
  
    if ( lodisjsets_BkinCi.getNumSet() != lomatrixrowt_S.getNumRows()  ) {
      std::cout << lpc_labelFunc << " warning number of sets and different groups "
		<< lodisjsets_BkinCi.getNumSet() << " <> " << lomatrixrowt_S.getNumRows()
		<< std::endl;
    }
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
 
  return std::make_tuple
    (lomatrixrowt_S,
     partition::PartitionDisjSets<T_CLUSTERIDX>
     (std::move(lodisjsets_BkinCi),lvectoruintidx_idxFirstInstVi),
     lovectort_numInstCi
     );

} /*tsengyang2001_getClusters*/

  

/*! \fn T_CLUSTERIDX updateCentroids(mat::MatrixRow<T_FEATURE> &aomatrixt_centroids, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Update centroid 
    \details Update centroid with themselves  the nearest neighbor rule. Also obtains the sum the sum of instances and number of instances per cluster. 
    \param aomatrixt_centroids a mat::MatrixRow<T_FEATURE> with the centroids of each cluster
    \param aomatrixt_sumInstancesCluster a mat::MatrixRow with the sum of instances per cluster
    \param aovectort_numInstancesInClusterK a std::vector with the number of instances per cluster
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast a InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
  template < typename INPUT_ITERATOR, 
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,  //0, 1, .., N
	   typename T_CLUSTERIDX,    //-1, 0, 1, .., K
	   typename T_DIST
	   >
void
updateCentroids
(T_CLUSTERIDX                        &aocidx_numClusterNull, 
 mat::MatrixRow<T_FEATURE>           &aiomatrixt_centroids, /*Mean for each cluster*/
 mat::MatrixRow<T_FEATURE_SUM>       &aomatrixt_sumInstancesCluster,
 std::vector<T_INSTANCES_CLUSTER_K>  &aovectort_numInstancesInClusterK,
 INPUT_ITERATOR                      aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
)       
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::updateCentroids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output mat::MatrixRow<T_FEATURE>: aimatrixt_centroids[" 
	      << &aiomatrixt_centroids << ']'
	      << "\n\t input aiiterator_instfirst[" << *aiiterator_instfirst << ']'
	      << "\n\t input aiiterator_instlast[" <<  *aiiterator_instlast << ']'
	      << "\n\t input  dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
 
  T_DIST  lT_distMinCentInst;
 
  const T_FEATURE lT_alpha = T_FEATURE(1);
  
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    
    const T_FEATURE* linst_inter =
      ((data::Instance<T_FEATURE>*) *aiiterator_instfirst)->getFeatures();

    T_CLUSTERIDX lmgidx_j = 
      nearest::NN
      <T_CLUSTERIDX,T_FEATURE,T_DIST>
      (lT_distMinCentInst,
       aiomatrixt_centroids,
       linst_inter,
       aifunc2p_dist
       );
    aovectort_numInstancesInClusterK[lmgidx_j]++;
    T_FEATURE_SUM   *larrarrowt_sumInstancesCluster = 
      aomatrixt_sumInstancesCluster.getRow(lmgidx_j);
    interfacesse::axpy
      (larrarrowt_sumInstancesCluster,
       lT_alpha,
       linst_inter,
       data::Instance<T_FEATURE>::getNumDimensions()
       );
  }
   
  meanCentroids
    (aocidx_numClusterNull,
     aiomatrixt_centroids,
     aomatrixt_sumInstancesCluster,
     aovectort_numInstancesInClusterK 
     ); 

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*! \fn void updateClusterCj(mat::MatrixRow<T_FEATURE> &aiomatrixt_centroids, std::vector<T_INSTANCES_CLUSTER_K> &aovectorit_numInstClusterK, T_CLUSTERIDX *aioarraycidx_memberShip, T_CLUSTERIDX aicidx_Cs, std::vector<T_CLUSTERIDX> &aivectorcidx_clustersNew, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Update cluster members Cj
    \details Update cluster members Cj, using the nearby centroid rule. This function is used to split a cluster
    \param aiomatrixt_centroids a mat::MatrixRow
    \param aovectorit_numInstClusterK a std::vector
    \param aioarraycidx_memberShip a an array with the labels belong to the cluster
    \param aicidx_Cs a cluster index to update
    \param aivectorcidx_clustersNew a std::vector<T_CLUSTERIDX>
    \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template < typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,  //0, 1, .., N
	   typename T_CLUSTERIDX,    //-1, 0, 1, .., K
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
void
updateClusterCj
(mat::MatrixRow<T_FEATURE>       &aiomatrixt_centroids, 
 mat::MatrixRow<T_FEATURE_SUM>   &aiomatrixt_sumInstCluster,
 std::vector
 <T_INSTANCES_CLUSTER_K>         &aovectorit_numInstClusterK,
 T_CLUSTERIDX                    *aioarraycidx_memberShip,
 T_CLUSTERIDX                    aicidx_Cj,  
 std::vector<T_CLUSTERIDX>       &aivectorcidx_clustersNew,
 INPUT_ITERATOR                  aiiterator_instfirst,
 const INPUT_ITERATOR            aiiterator_instlast,
 dist::Dist<T_DIST,T_FEATURE>    &aifunc2p_dist
 )       
{
#ifdef __VERBOSE_YES
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  const char* lpc_labelFunc = "clusteringop::generateNewClusters";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(inout mat::MatrixRow<>: aimatrixt_centroids[" 
      	      << &aiomatrixt_centroids << ']'   
	      << "\ninout std::vector<>: &aovectorit_numInstClusterK["
	      << &aovectorit_numInstClusterK << ']'
	      << "\ninout T_CLUSTERIDX aioarraycidx_memberShip[" 
	      << aioarraycidx_memberShip << ']'
	      << "\ninput T_CLUSTERIDX  aicidx_Cj = " << aicidx_Cj
	      << "\ninput std::vector<> aivectorcidx_clustersNew[" 
	      << &aivectorcidx_clustersNew << ']'
              << "input const aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "\ninput dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
 
  aiomatrixt_sumInstCluster.initialize();

  const T_FEATURE lT_alpha = T_FEATURE(1);
  T_DIST  lT_distMinCentInst;
 
  for (uintidx lui_idxInsti = 0 ; aiiterator_instfirst != aiiterator_instlast;
       aiiterator_instfirst++, lui_idxInsti++) {
    if ( aicidx_Cj == aioarraycidx_memberShip[lui_idxInsti] ) { 
      const T_FEATURE* linst_inter = (*aiiterator_instfirst)->getFeatures();
 
      T_CLUSTERIDX lmgidx_j = 
	nearest::NN
	<T_CLUSTERIDX,T_FEATURE,T_DIST>
	(lT_distMinCentInst,
	 aiomatrixt_centroids,
	 linst_inter,
	 aifunc2p_dist
	 );
      aovectorit_numInstClusterK[lmgidx_j]++;
      aioarraycidx_memberShip[lui_idxInsti] = aivectorcidx_clustersNew[lmgidx_j]; 
    
      T_FEATURE_SUM *larrarrowt_sumInstancesCluster = 
	aiomatrixt_sumInstCluster.getRow(lmgidx_j);
      interfacesse::axpy
	(larrarrowt_sumInstancesCluster,
	 lT_alpha,
	 linst_inter,
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
    }
  }
   
  T_CLUSTERIDX  lcidx_numClusterNull;

  meanCentroids
    (lcidx_numClusterNull,
     aiomatrixt_centroids,
     aiomatrixt_sumInstCluster,
     aovectorit_numInstClusterK 
     ); 

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS:" << lpc_labelFunc;
    aiomatrixt_centroids.print
      (std::cout,
       lostrstream_labelCentroids.str().c_str(),
       ',',
       ';'
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip
      << "<MEMBERCLUSTER:"
      << geverbosepc_labelstep
      << ':' << lpc_labelFunc
      << ':' << geverboseui_idproc
      << ":aioarraycidx_memberShip[" 
      << aioarraycidx_memberShip << ']';
    inout::containerprint
      (aioarraycidx_memberShip,
       aioarraycidx_memberShip + lui_numInstances,
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelNumInst;
    lostrstream_labelNumInst
      << "<NUMINSTANCES:"
      << lpc_labelFunc
      << ":aovectorit_numInstClusterK[" 
      << &aovectorit_numInstClusterK << ']';
    inout::containerprint
      (aovectorit_numInstClusterK.begin(),
       aovectorit_numInstClusterK.end(),
       std::cout,
       lostrstream_labelNumInst.str().c_str(),
       ','
       );
    std::cout << std::endl;
      
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

  
/*! \fn T_CLUSTERIDX kmeansoperator(T_CLUSTERIDX *aioarraycidx_memberShip, mat::MatrixRow<T_FEATURE> &aomatrixt_centroids, mat::MatrixRow<T_FEATURE_SUM>  &aomatrixt_sumInstancesCluster, std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstancesInClusterK, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief K-Means Operator \cite Krishna:Murty:GAClustering:GKA:1999
    \details For the given partition in an array of labels, it performs an update using the k-means algorithm. Returns the membership tags, centroids, sum of instances and number of instances in each cluster.
    \param aioarraycidx_memberShip  an array with the labels belong to the cluster
    \param aomatrixt_centroids a mat::MatrixRow with the centroids of each cluster
    \param aomatrixt_sumInstancesCluster a  mat::MatrixRow<T_FEATURE_SUM>
    \param aovectort_numInstancesInClusterK a std::vector<T_INSTANCES_CLUSTER_K>    
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast a InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX, //-1, 0, 1, .., K
	  typename T_FEATURE,
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_DIST,
	  typename INPUT_ITERATOR
	  >
T_CLUSTERIDX 
kmeansoperator
(T_CLUSTERIDX                       *aioarraycidx_memberShip,
 mat::MatrixRow<T_FEATURE>          &aomatrixt_centroids,
 mat::MatrixRow<T_FEATURE_SUM>      &aomatrixt_sumInstancesCluster,
 std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstancesInClusterK,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 dist::Dist<T_DIST,T_FEATURE>       &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::kmeansoperator";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output T_CLUSTERIDX*: aioarraycidx_memberShip["
	      << aioarraycidx_memberShip << "]\n"
	      << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "input numclusterK = " << aomatrixt_centroids.getNumRows() << '\n'
	      << "input dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist[" << &aifunc2p_dist << "]\n"
	      << ")\n";
  }
#endif /*__VERBOSE_YES*/

  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  partition::PartitionLabel<T_CLUSTERIDX>
    lpartition_clustersLabel
    (aioarraycidx_memberShip,
     lui_numInstances,
     aomatrixt_centroids.getNumRows()
     );

  T_CLUSTERIDX locidx_numClusterNull =
    getCentroids
    (aomatrixt_centroids,
     aomatrixt_sumInstancesCluster,
     aovectort_numInstancesInClusterK,
     lpartition_clustersLabel,
     aiiterator_instfirst,
     aiiterator_instlast
     );

  reassignCluster 
    (aioarraycidx_memberShip,
     aomatrixt_centroids, 
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:" << lpc_labelFunc;
      inout::containerprint
      (aioarraycidx_memberShip,
       aioarraycidx_memberShip + lui_numInstances, 
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
      
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return locidx_numClusterNull;
}


/*! \fn std::pair<bool,T_DIST> reassignCluster(ds::PartitionLinkedStats<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K, T_FEATURE_SUM> &aoipartlinkstats_partition, const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist)
    \brief Change each instance to the nearest cluster
    \details
    \param aoipartlinkstats_partition a ds::PartitionLinkedStats for specifying cluster membership
    \param aimatrixt_centroids a mat::MatrixBase with the centroids of each cluster
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a const input iterator of the instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_INSTANCE_FREQUENCY,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_FEATURE_SUM,
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
std::pair<bool,T_DIST>
reassignCluster
(ds::PartitionLinkedStats
 <T_FEATURE,
  T_CLUSTERIDX,
  T_INSTANCE_FREQUENCY,
  T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM>                      &aoipartlinkstats_partition,
 const mat::MatrixBase<T_FEATURE>    &aimatrixt_centroids,
 INPUT_ITERATOR                      aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
 )
{
  bool    lob_validPartitioning = false;
  T_DIST  lort_distortionDist   = 0.0;
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::reassignCluster";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output PartitionLinkedStats: aoipartlinkstats_partition["
	      << &aoipartlinkstats_partition << "]\n";

    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    
    std::cout << " input  mat::MatrixBase<T_FEATURE>: aimatrixt_centroids[" 
	      <<  &aimatrixt_centroids << "]\n"
              << " input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input  dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES
    
  uintidx lidxinstT_i = 0;

  T_DIST                lrt_distMinCentInst;
  T_INSTANCES_CLUSTER_K lit_totalInstFreq = 0;
  
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    
       data::InstanceFreq
	 <T_FEATURE,
	  T_INSTANCE_FREQUENCY> *liter_instFreq =
	 (data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
	 *aiiterator_instfirst;
	 
       T_CLUSTERIDX lmgidx_j =
	 nearest::NN
	 (lrt_distMinCentInst,
	  aimatrixt_centroids,
	  aoipartlinkstats_partition.getMemberShip(lidxinstT_i),
	  liter_instFreq->getFeatures(),
	  aifunc2p_dist
	  );
       
       if ( lmgidx_j != aoipartlinkstats_partition.getMemberShip(lidxinstT_i) ) {
	 
	 aoipartlinkstats_partition.changeMemberShip
	   (lmgidx_j,
	    lidxinstT_i,
	    liter_instFreq->getFeatures(),
	    liter_instFreq->getFrequency()
	    );
	 
       } 
       ++lidxinstT_i;
       lit_totalInstFreq += liter_instFreq->getFrequency();
       lort_distortionDist += lrt_distMinCentInst * liter_instFreq->getFrequency();
	 
     }
     
  /*Check is valid parttition*/
  lob_validPartitioning =
    ( aoipartlinkstats_partition.numClusterEmpty() == 0 );

  lort_distortionDist /=
    ( (T_DIST)  lit_totalInstFreq * data::Instance<T_FEATURE>::getNumDimensions() );
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " validPartitioning  =  "  << lob_validPartitioning
              << ", distortionDistance =  "  << lort_distortionDist << '\n';
    
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    
    std::cout << '\n'
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_pair(lob_validPartitioning,lort_distortionDist);
    
} /*END reassignCluster
   */
  

/*! \fn void reassignCluster(ds::PartitionLinkedStats<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM> &aoipartlinkstats_partition, const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, T_CLUSTERIDX *aiarraycidx_newIndex, uintidx aiui_countNewIndex, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Reassign cluster
    \details Reassign cluster based on a local search. Returns if the partition is valid and the distance distortion metric
    \param aoipartlinkstats_partition a ds::PartitionLinkedStats for specifying cluster membership
    \param aimatrixt_centroids a mat::MatrixBase with the centroids of each cluster
    \param aiarraycidx_newIndex an array with the new indexes
    \param aiui_countNewIndex number of new indices
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a const input iterator of the instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template <typename T_FEATURE,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM, 
	  typename T_DIST,
	  typename INPUT_ITERATOR
	  >
void reassignCluster
(ds::PartitionLinkedStats
 <T_FEATURE,
  T_CLUSTERIDX,
  T_INSTANCE_FREQUENCY,
  T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM>                     &aoipartlinkstats_partition,
 const mat::MatrixBase<T_FEATURE>   &aimatrixt_centroids,
 T_CLUSTERIDX                       *aiarraycidx_newIndex,
 uintidx                            aiui_countNewIndex,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::reassignCluster";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output PartitionLinkedStats: aoipartlinkstats_partition["
	      << &aoipartlinkstats_partition << "]\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    std::cout << " input  mat::MatrixBase<T_FEATURE>: aimatrixt_centroids[" 
	      <<  &aimatrixt_centroids << "]\n"
	      << " input  T_CLUSTERIDX*: aiarraycidx_newIndex[" << aiarraycidx_newIndex << "]\n"
	      << " input  uintidx: aiui_countNewIndex = " << aiui_countNewIndex << '\n'
	      << " input  aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input  dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  T_CLUSTERIDX lcidx_codeBookNumClusterK = (T_CLUSTERIDX) (aimatrixt_centroids.getNumRows());

  for (uintidx lui_idxInsti = 0 ; aiiterator_instfirst != aiiterator_instlast;
       aiiterator_instfirst++, lui_idxInsti++) {
    
    T_CLUSTERIDX lcidx_best = -1;
    data::InstanceFreq
      <T_FEATURE,
       T_INSTANCE_FREQUENCY> *literinstfo_iInstance  =
      (data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
      *aiiterator_instfirst;
      
    T_DIST  lrt_olddist = 
      aifunc2p_dist
      (aimatrixt_centroids.getRow(aoipartlinkstats_partition.getMemberShip(lui_idxInsti)),
       literinstfo_iInstance->getFeatures(),
       literinstfo_iInstance->getNumDimensions()
       );
      
    bool  lb_changeInstMemberShip = false;

    //Check if the instance is in the changing cluster
    for ( uintidx lui_j = 0; lui_j < aiui_countNewIndex; lui_j++)
      if ( aoipartlinkstats_partition.getMemberShip(lui_idxInsti) == aiarraycidx_newIndex[lui_j] ) {
	 
	lb_changeInstMemberShip = true;
	break;
      }
    if ( lb_changeInstMemberShip ) {  

      //If you belong look for the cluster with minimum distance or -1 if it remains the same
      for ( T_CLUSTERIDX lcidx_j = 0; 
	    lcidx_j < lcidx_codeBookNumClusterK; 
	    lcidx_j++) {
	T_DIST lrt_newdist = 
	  aifunc2p_dist
	  (aimatrixt_centroids.getRow( lcidx_j ),
	   literinstfo_iInstance->getFeatures(),
	   literinstfo_iInstance->getNumDimensions()
	   );
	  
	if ( lrt_newdist < lrt_olddist ) {
	  lcidx_best    = lcidx_j;
	  lrt_olddist = lrt_newdist;
	}
      }
    }
    else {  //Belong to one that does not change
      for (uintidx lui_j = 0; lui_j < aiui_countNewIndex; lui_j++ ) {
	T_DIST lrt_newdist = 
	  aifunc2p_dist
	  (aimatrixt_centroids.getRow(aiarraycidx_newIndex[lui_j]),
	   literinstfo_iInstance->getFeatures(),
	   literinstfo_iInstance->getNumDimensions()
	   );
	 
	if ( lrt_newdist < lrt_olddist ) {
	  lcidx_best    = aiarraycidx_newIndex[lui_j];
	  lrt_olddist = lrt_newdist;
	}
      }
    }
    
    if ( lcidx_best > -1 ) 
      aoipartlinkstats_partition.changeMemberShip
	(lcidx_best,
	 lui_idxInsti,
	 literinstfo_iInstance->getFeatures(),
	 literinstfo_iInstance->getFrequency()
	 );
      
  } //For lui_idxInsti 

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*END reassignCluster */


/*! \fn void  joinPartition (ds::PartitionLinkedStats<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM> &aoipartlinkstats_partition, const T_CLUSTERIDX aicidx_partitionTo, const T_CLUSTERIDX aicidx_partitionFrom, const T_CLUSTERIDX aicidx_numClusterK, const INPUT_ITERATOR aiiterator_instfirst)
    \brief Join partition
    \details
    \param aoipartlinkstats_partition a ds::PartitionLinkedStats to join, this also specifies membership in the cluster
    \param aicidx_partitionFrom index of the partition to move to another
    \param aicidx_partitionTo index of the partition to which they join
    \param aiiterator_instfirst a iterator with the instances
*/
template < typename T_FEATURE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_INSTANCE_FREQUENCY,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_FEATURE_SUM,
	   typename INPUT_ITERATOR
>
void  joinPartition
(ds::PartitionLinkedStats
 <T_FEATURE,
 T_CLUSTERIDX,
 T_INSTANCE_FREQUENCY,
 T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM>             &aoipartlinkstats_partition,
 const T_CLUSTERIDX         aicidx_partitionFrom,
 const T_CLUSTERIDX         aicidx_partitionTo,
 const INPUT_ITERATOR       aiiterator_instfirst
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::joinPartition";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    std::cout << "aicidx_partitionFrom: " << aicidx_partitionFrom
	      << "\taicidx_partitionTo: " <<  aicidx_partitionTo
	      << "\ninput const aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_j(&aoipartlinkstats_partition);
   
  for ( literpart_j.begin(aicidx_partitionFrom); literpart_j.end(); literpart_j.next() ) {
    
    data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>* lptinstfo_changeInstance = 
      (data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
      *std::next(aiiterator_instfirst,literpart_j.getValue());

    aoipartlinkstats_partition.setMemberShip
      (literpart_j.getValue(),
       aicidx_partitionTo
       );
    aoipartlinkstats_partition.changeSumInstances
      (aicidx_partitionFrom,
       aicidx_partitionTo,
       lptinstfo_changeInstance->getFeatures(),
       lptinstfo_changeInstance->getFrequency()
       );    
  }
 
  aoipartlinkstats_partition.joinCluster(aicidx_partitionFrom,aicidx_partitionTo);

  const T_CLUSTERIDX  lcidx_numClusterK =  aoipartlinkstats_partition.getNumPartitions() -1;
  
  if ( aicidx_partitionFrom != (lcidx_numClusterK) ) {
 
    while ( aoipartlinkstats_partition.getFirstInstClusterK(lcidx_numClusterK)
	    != UINTIDX_NIL )
      {

	data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>* lptinstfo_changeInstance =
	  (data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
	  *std::next
	  (aiiterator_instfirst,
	   aoipartlinkstats_partition.getFirstInstClusterK(lcidx_numClusterK)
	   );
      
	aoipartlinkstats_partition.changeMemberShip 
	  (aicidx_partitionFrom,
	   aoipartlinkstats_partition.getFirstInstClusterK(lcidx_numClusterK),
	   lptinstfo_changeInstance->getFeatures(),
	   lptinstfo_changeInstance->getFrequency()
	   );	
      }
  }

  aoipartlinkstats_partition.decreasePartition();
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << std::endl;  
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
}


/*! \fn void fillEmptyPartitions (ds::PartitionLinkedStats<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY, T_INSTANCES_CLUSTER_K, T_FEATURE_SUM> &aiopartlinkstats_partition, mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K> &aiomatrixresizerow_centroids, const INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist)
    \brief Fill empty partitions
    \details For each empty partition adds a new centroid
    \param aiopartlinkstats_partition a ds::PartitionLinkedStats to fill empty clusters
    \param aiomatrixresizerow_centroids a mat::MatrixResizableRow with the centroids of each cluster
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a const input iterator of the instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template <typename T_FEATURE,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM, 
	  typename T_DIST,
	  typename INPUT_ITERATOR
	  >
void fillEmptyPartitions
(ds::PartitionLinkedStats
 <T_FEATURE,
  T_CLUSTERIDX,
  T_INSTANCE_FREQUENCY,
  T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM>                      &aiopartlinkstats_partition,
 mat::MatrixResizableRow
 <T_FEATURE,T_INSTANCES_CLUSTER_K>   &aiomatrixresizerow_centroids,
 const INPUT_ITERATOR                aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::fillEmptyPartitions";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aiopartlinkstats_partition";
    aiopartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS:aiomatrixresizerow_centroids" << lpc_labelFunc;
    aiomatrixresizerow_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    std::cout << '\n';
    
    std::cout << " input  aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
              << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input  dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << ")"
	      << std::endl;		
  }
#endif //__VERBOSE_YES

  T_CLUSTERIDX       lcidx_nextNew;
  uintidx            *larrayst_newIndex = NULL; 
  T_DIST             *larrayNorm_error = NULL; 
  T_DIST             lnorm_newerror;     
  T_CLUSTERIDX       lcidx_numClusterEmpty = 0;
  T_CLUSTERIDX       lcidx_oldClusterK;
  T_CLUSTERIDX       lcidx_codeBookNumClusterK;
  
  lcidx_codeBookNumClusterK = (T_CLUSTERIDX) aiomatrixresizerow_centroids.getNumRows();

  lcidx_numClusterEmpty =
    aiopartlinkstats_partition.numClusterEmpty();
   
  if (lcidx_numClusterEmpty == 0) {

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    return;
  }

  //Allocate temporary memory 
  larrayst_newIndex  =  new uintidx[(uintidx) (lcidx_numClusterEmpty + 1)];
  larrayNorm_error   =  new T_DIST[(uintidx) (lcidx_numClusterEmpty + 1) ];
   
  //Find new vectors
  lcidx_nextNew = 0;

  INPUT_ITERATOR literator_instfirst = aiiterator_instfirst;
  for (uintidx  luintidx_i = 0;
       literator_instfirst != aiiterator_instlast;
       ++literator_instfirst, luintidx_i++)
    {
      data::InstanceFreq
	<T_FEATURE,
	 T_INSTANCE_FREQUENCY> *liter_instFreq  = 
	(data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
	*aiiterator_instfirst;
     
      lnorm_newerror =
	larrayNorm_error[lcidx_nextNew] = 
	aifunc2p_dist
	(aiomatrixresizerow_centroids.getRow(aiopartlinkstats_partition.getMemberShip(luintidx_i)),
	 liter_instFreq->getFeatures(),
	 liter_instFreq->getNumDimensions()
	 );
     
      T_CLUSTERIDX lcidx_n = lcidx_nextNew;
      while( lcidx_n > 0 && larrayNorm_error[lcidx_n] > larrayNorm_error[lcidx_n - 1] ) {
	larrayst_newIndex[lcidx_n] = larrayst_newIndex[lcidx_n - 1];
	larrayNorm_error[lcidx_n] = larrayNorm_error[lcidx_n - 1];
	lcidx_n--;
      }
      
      larrayst_newIndex[lcidx_n] = luintidx_i;
      larrayNorm_error[lcidx_n] = lnorm_newerror;
      if( lcidx_nextNew < lcidx_numClusterEmpty  ) {
	lcidx_nextNew++;
      }
      
    }

  for( T_CLUSTERIDX lcidx_j = 0, lcidx_n = 0; 
       lcidx_n < lcidx_numClusterEmpty && lcidx_j <  lcidx_codeBookNumClusterK; 
       lcidx_j++ ) 
    { //Begin for
	
      if ( (aiopartlinkstats_partition.getNumInstancesClusterKi(lcidx_j) 
	    == T_INSTANCES_CLUSTER_K(0)) )  
	{  
	  lcidx_oldClusterK = aiopartlinkstats_partition.getMemberShip( larrayst_newIndex[lcidx_n] );
	    
	  T_FEATURE* larrayrow_centroids = 
	    aiomatrixresizerow_centroids.getRow((uintidx)lcidx_j);
	  data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>* liter_instFreq = 
	    (data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>*)
	    *std::next(aiiterator_instfirst,larrayst_newIndex[lcidx_n]);
	  
	  interfacesse::copy
	    (larrayrow_centroids,
	     liter_instFreq->getFeatures(),
	     liter_instFreq->getNumDimensions()
	     );
	   
	  aiopartlinkstats_partition.changeMemberShip
		 (lcidx_j,
		  larrayst_newIndex[lcidx_n],
		  liter_instFreq->getFeatures(),
		  liter_instFreq->getFrequency()
		  );
	    
	  if( aiopartlinkstats_partition.getNumInstancesClusterKi(lcidx_oldClusterK) > 0 ) {
	     
	    stats::meanVector  
	      (aiomatrixresizerow_centroids.getRow(lcidx_oldClusterK),
	       aiopartlinkstats_partition.getNumInstancesClusterKi(lcidx_oldClusterK),
	       aiopartlinkstats_partition.getSumInstancesClusterKi(lcidx_oldClusterK),
	       (uintidx) lcidx_j
	       );
	  }
	   
	  lcidx_n++; 
	} 
    } //End for
  delete[] larrayst_newIndex;
  delete[] larrayNorm_error;

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aiopartlinkstats_partition";
    aiopartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS:aiomatrixresizerow_centroids" << lpc_labelFunc;
    aiomatrixresizerow_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

} /*fillEmptyPartitions*/


/*! \fn void removeEmptyPartitions(ds::PartitionLinkedStats<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM> &aoipartlinkstats_partition, mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K> &aomatrixresizerow_centroids, const INPUT_ITERATOR aiiterator_instfirst, const T_CLUSTERIDX aicidx_numclusterKTopMin
 )
    \brief Remove empty partitions
    \details
    \param aoipartlinkstats_partition a ds::PartitionLinkedStats to join, this also specifies membership in the cluster
    \param aomatrixresizerow_centroids a mat::MatrixResizableRow with the centroids of each cluster
    \param aiiterator_instfirst a input iterator of the instances
    \param aicidx_numclusterKTopMin a minimum number of clusters to maintain
*/
template < typename T_FEATURE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_INSTANCE_FREQUENCY,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_FEATURE_SUM,
	   typename INPUT_ITERATOR
>
void removeEmptyPartitions
(ds::PartitionLinkedStats
 <T_FEATURE,
  T_CLUSTERIDX, //-1, 0, 1, .., K
  T_INSTANCE_FREQUENCY,
  T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM>                     &aoipartlinkstats_partition,
 mat::MatrixResizableRow
 <T_FEATURE,T_INSTANCES_CLUSTER_K>  &aomatrixresizerow_centroids,
 const INPUT_ITERATOR               aiiterator_instfirst,
 const T_CLUSTERIDX                 aicidx_numclusterKTopMin
 )
{
#ifdef __VERBOSE_YES
   /*Parameters*/
  bool  lb_verbose_removeEmptyPartitions = false;
  const char* lpc_labelFunc = "clusteringop::removeEmptyPartitions";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output ds::PartitionLinkedStats: aoipartlinkstats_partition["
	      << &aoipartlinkstats_partition << "]\n"
	      << " input const aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input T_CLUSTERIDX: aicidx_numclusterKTopMin "
	      << aicidx_numclusterKTopMin << '\n'
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
 
  T_CLUSTERIDX lcidx_codeBookNumClusterK =
    (T_CLUSTERIDX) aomatrixresizerow_centroids.getNumRows();
  
  for ( T_CLUSTERIDX lcidx_j = 0; lcidx_j < lcidx_codeBookNumClusterK; lcidx_j++) {
    
    if ( (aoipartlinkstats_partition.getNumInstancesClusterKi(lcidx_j) 
	  == T_INSTANCES_CLUSTER_K(0))  &&
	 lcidx_codeBookNumClusterK > aicidx_numclusterKTopMin ){

#ifdef __VERBOSE_YES
      lb_verbose_removeEmptyPartitions = true;
#endif /*__VERBOSE_YES*/

      T_CLUSTERIDX lcidx_clusterLast =  lcidx_codeBookNumClusterK -1;
      
      aomatrixresizerow_centroids.swapRows
	((uintidx) lcidx_j, (uintidx) lcidx_clusterLast);

      joinPartition
	(aoipartlinkstats_partition,
	 lcidx_clusterLast,
	 lcidx_j,
	 aiiterator_instfirst
	 );
      aomatrixresizerow_centroids.decreaseRow();
      lcidx_codeBookNumClusterK = (T_CLUSTERIDX) aomatrixresizerow_centroids.getNumRows();
      --lcidx_j;
    } 	 
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    if (  lb_verbose_removeEmptyPartitions ) {
      std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    }
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

} /* END removeEmptyPartitions */

  
/*! \fn void pnnFast (ds::PartitionLinkedStats<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K, T_FEATURE_SUM> &aoipartlinkstats_partition, mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K> &aiomatrixresizerow_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATORaiiterator_instlast, const T_CLUSTERIDX aicidx_numclusterKToReduce, const dist::Dist<T_DIST,T_FEATURE>&aifunc2p_dist)
  \brief PNN fast algorithm \cite Franti:etal:GAclustering:gafranti:1997
  \details Reduce the number of partitions to a number specified by the PNN algorithm fast
  \param aoipartlinkstats_partition a ds::PartitionLinkedStats for specifying cluster membership
  \param aiomatrixresizerow_centroids a mat::MatrixResizableRow with the centroids of each cluster
  \param aiiterator_instfirst a input iterator of the instances
  \param aiiterator_instlast a const input iterator of the instances
  \param aicidx_numclusterKToReduce a number of cluster to reduce
*/
template <typename T_FEATURE, 
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename T_DIST,
	  typename INPUT_ITERATOR
	  >
void pnnFast
(ds::PartitionLinkedStats
 <T_FEATURE,
  T_CLUSTERIDX,
  T_INSTANCE_FREQUENCY,
  T_INSTANCES_CLUSTER_K,
 T_FEATURE_SUM>                      &aoipartlinkstats_partition,
 mat::MatrixResizableRow
 <T_FEATURE,T_INSTANCES_CLUSTER_K>   &aiomatrixresizerow_centroids,
 INPUT_ITERATOR                      aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const T_CLUSTERIDX                  aicidx_numclusterKToReduce,
 const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES


  const char* lpc_labelFunc = "gaclusteringop::pnnFast";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS:aiomatrixresizerow_centroids" << lpc_labelFunc;
    aiomatrixresizerow_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    
    std::cout << " input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input T_CLUSTERIDX: aicidx_numclusterKToReduce "
	      << aicidx_numclusterKToReduce << '\n'
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  clusteringop::removeEmptyPartitions
    (aoipartlinkstats_partition,
     aiomatrixresizerow_centroids,
     aiiterator_instfirst,
     aicidx_numclusterKToReduce
     );
  
  T_CLUSTERIDX lcidx_numClusters =
    (T_CLUSTERIDX) aiomatrixresizerow_centroids.getNumRows();
  
  std::vector<nearest::NearestCentroids<T_CLUSTERIDX,T_DIST> >
    lvector_nearestCentroidsDist(lcidx_numClusters);
  
  nearest::findAllPairCentroidsClosest 
    (lvector_nearestCentroidsDist,
     aiomatrixresizerow_centroids,
     aoipartlinkstats_partition.getNumInstancesClusterK(),
     aifunc2p_dist
     );

  while ( lcidx_numClusters >  aicidx_numclusterKToReduce ) {

    std::pair<T_CLUSTERIDX,T_CLUSTERIDX> lpair_centroidIdxToMerged =
      nearest::findNearestCentroids(lvector_nearestCentroidsDist);

    /*Begin merge
     */    
    nearest::markRecalculationNearest
      (lvector_nearestCentroidsDist,
       lpair_centroidIdxToMerged
       );

    aiomatrixresizerow_centroids.mergeTwoRow
      (lpair_centroidIdxToMerged,
       aoipartlinkstats_partition.getNumInstancesClusterKi( lpair_centroidIdxToMerged.first  ),
       aoipartlinkstats_partition.getNumInstancesClusterKi( lpair_centroidIdxToMerged.second  )
       );
    
    joinPartition
      (aoipartlinkstats_partition,
       lpair_centroidIdxToMerged.second,
       lpair_centroidIdxToMerged.first,
       aiiterator_instfirst
       );
    aiomatrixresizerow_centroids.removeRow((uintidx) lpair_centroidIdxToMerged.second);
   
    nearest::deleteNearestCentroid
      (lvector_nearestCentroidsDist,
       lpair_centroidIdxToMerged.second
       );
    /*End merge
     */
      
    nearest::recalculateDistNearestCentroids
      (lvector_nearestCentroidsDist,
       aiomatrixresizerow_centroids,
       aoipartlinkstats_partition.getNumInstancesClusterK(),
       aifunc2p_dist
       );

    --lcidx_numClusters;

  } /* While */

  aoipartlinkstats_partition.resize( (uintidx) aicidx_numclusterKToReduce);
  aiomatrixresizerow_centroids.resize( (uintidx) aicidx_numclusterKToReduce);

#ifdef __VERBOSE_YES

  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition << lpc_labelFunc << ":aoipartlinkstats_partition";
    aoipartlinkstats_partition.print(std::cout,lostrstream_labelPartition.str().c_str(),',');
    std::cout << '\n';
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS:aiomatrixresizerow_centroids" << lpc_labelFunc;
    aiomatrixresizerow_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
}

} /*END namespace clusteringop*/
  
#endif /*__CLUSTERING_OPERATOR_CENTROIDS_HPP*/
