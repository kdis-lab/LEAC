/*! \file nearestinstance_operator.hpp
 *
 * \brief nearest instance operator
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __NEAREST_INSTANCES_OPERATOR_HPP
#define __NEAREST_INSTANCES_OPERATOR_HPP

#include "instance.hpp"
#include "matrix.hpp"
#include "matrix_triangular.hpp"
#include "dist.hpp"
#include "vector_utils.hpp"

#include "verbose_global.hpp"

/*! \namespace nearest
  \brief Function to find nearest instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace nearest {

  
#define NEARESTCENTROID_UNKNOWN -1
#define NEARESTCENTROID_NOISE   -2


/*! \fn uintidx farthestInstanceFromS1(const uintidx aiuintidx_S1, ITERATOR_IDXINST aiiterator_instfirst, ITERATOR_IDXINST aiiterator_instlast, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Farthest instance from S1 \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003
    \details Find the instance furthest from s1. The sequence is the set of indexes of intents that belong to a cluster 
    \param aiuintidx_S1 a index object S1
    \param aiiterator_instfirst a iterator to the initial positions of the sequence to search
    \param aiiterator_instlast  a iterator to the final positions of the sequence to search
    \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template < typename ITERATOR_IDXINST,//-1, 0, 1, .., K
           typename T_FEATURE, 
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
uintidx
farthestInstanceFromS1
(const uintidx                          aiuintidx_S1,  //Index object S1
 ITERATOR_IDXINST                       aiiterator_idxInstfirst,
 ITERATOR_IDXINST                       aiiterator_idxInstlast,
 const INPUT_ITERATOR                   aiiterator_instfirst,
 const dist::Dist<T_DIST,T_FEATURE>     &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "nearest::farthestInstanceFromS1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')' << " S1 = " <<  aiuintidx_S1
	      << std::endl;
  }
#endif //__VERBOSE_YES

  uintidx louintidx_S2 = aiuintidx_S1;
  T_DIST  lrt_maxDistS1toS2 = T_DIST(0);
  
  const T_FEATURE* linst_featuresS1 =
    (*std::next(aiiterator_instfirst,aiuintidx_S1))->getFeatures();
  
  for (; aiiterator_idxInstfirst != aiiterator_idxInstlast; ++aiiterator_idxInstfirst) {

    T_DIST lrt_distS1toSi = 
      aifunc2p_dist 
      (linst_featuresS1, 
       (*std::next(aiiterator_instfirst,*aiiterator_idxInstfirst))->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    if ( lrt_maxDistS1toS2 < lrt_distS1toSi ) {
      louintidx_S2 = *aiiterator_idxInstfirst;
      lrt_maxDistS1toS2 = lrt_distS1toSi;
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')' << " S2 = " << louintidx_S2
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return louintidx_S2;

}


/*! \fn uintidx farthestInstanceFromS1(const uintidx aiuintidx_S1, const T_CLUSTERIDX *aiarraycidx_idxMemberShip, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)                      
    \brief Farthest instance from S1 \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006 
    \details Find the instance furthest from s1. The sequence is the set of indexes of intents that belong to a cluster. The belonging to a group is defined by a string of labels
    \param aiuintidx_S1 a Index of instance s1
    \param aiarraycidx_idxMemberShip an array with the labels belong to the cluster
    \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_FEATURE, 
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
uintidx
farthestInstanceFromS1
(const uintidx                 aiuintidx_S1,
 const T_CLUSTERIDX            *aiarraycidx_idxMemberShip,
 INPUT_ITERATOR                aiiterator_instfirst,
 const INPUT_ITERATOR          aiiterator_instlast,
 dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "nearest::farthestInstanceFromS1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')' << " S1 = " <<  aiuintidx_S1
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  T_CLUSTERIDX lcidx_S1inKj =  aiarraycidx_idxMemberShip[aiuintidx_S1];
  uintidx louintidx_S2     = aiuintidx_S1;
  T_DIST lrt_maxDistS1toS2 = T_DIST(0);

  const T_FEATURE* linst_featuresS1 =
    (*std::next(aiiterator_instfirst,aiuintidx_S1))->getFeatures();
  
  for (uintidx lui_i = 0 ; aiiterator_instfirst != aiiterator_instlast;
	 aiiterator_instfirst++, lui_i++) {
    
    if ( lcidx_S1inKj == aiarraycidx_idxMemberShip[lui_i] ) {
      T_DIST lrt_distS1toSi = 
	aifunc2p_dist 
	(linst_featuresS1,
	 (*aiiterator_instfirst)->getFeatures(),
	
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
      if ( lrt_maxDistS1toS2 < lrt_distS1toSi ) {
	louintidx_S2 = lui_i;
	lrt_maxDistS1toS2 = lrt_distS1toSi;
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')' << " S2 = " << louintidx_S2
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return louintidx_S2;

}


/*! \fn T_CLUSTERIDX NN(T_DIST &aort_distMinCentInst, const mat::MatrixRow<T_FEATURE> &aimatrixrowt_centroids, const T_FEATURE *aiat_instance, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief NN Find the centroid closest to an instance 
    \details 
    \param aort_distMinCentInst output distance between centroid and instance 
    \param aimatrixrowt_centroids matrix with the centroids
    \param aiat_instance a array with the instance attributes
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX, 
	  typename T_FEATURE, 
	  typename T_DIST
	  >
T_CLUSTERIDX
NN
(T_DIST                             &aort_distMinCentInst,
 const mat::MatrixRow<T_FEATURE>    &aimatrixrowt_centroids,
 const T_FEATURE                    *aiat_instance,
 const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist
 )
{
  T_CLUSTERIDX oidxT_cj = NEARESTCENTROID_UNKNOWN;
  T_CLUSTERIDX lidxT_K;
  T_DIST              lrt_distTiCi;
  
  oidxT_cj = 0;
  aort_distMinCentInst =  
    aifunc2p_dist 
    (aimatrixrowt_centroids.getRow(0), 
     aiat_instance, 
     aimatrixrowt_centroids.getNumColumns()
     );
  lidxT_K = T_CLUSTERIDX(aimatrixrowt_centroids.getNumRows());
  for (T_CLUSTERIDX li_j = 1; li_j < lidxT_K; li_j++) {
    lrt_distTiCi = 
      aifunc2p_dist
      (aimatrixrowt_centroids.getRow(li_j),aiat_instance, aimatrixrowt_centroids.getNumColumns());
    if (aort_distMinCentInst > lrt_distTiCi) {
      oidxT_cj = li_j;
      aort_distMinCentInst = lrt_distTiCi;
    }
  }

  return oidxT_cj;
}


/*! \fn T_CLUSTERIDX checkNullCentroidsNN(T_DIST &aort_distMinCentInst, const mat::MatrixBase<T_FEATURE> &aimatrixrowt_centroids, const T_FEATURE *aiat_instance, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Find the centroid closest to an instance and check null centroids
    \details 
    \param aort_distMinCentInst output distance between centroid and instance
    \param aiat_instance characteristics of the instance
    \param aimatrixrowt_centroids a matrix with centroids
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX, 
	  typename T_FEATURE, 
	  typename T_DIST
	  >
T_CLUSTERIDX
checkNullCentroidsNN
(T_DIST                           &aort_distMinCentInst,
 const mat::MatrixBase<T_FEATURE> &aimatrixrowt_centroids,
 const T_FEATURE                  *aiat_instance,
 dist::Dist<T_DIST,T_FEATURE>     &aifunc2p_dist
 )
{
  T_CLUSTERIDX oidxT_cj = NEARESTCENTROID_UNKNOWN;
  T_CLUSTERIDX lidxT_K;
  T_DIST              lrt_distTiCi;
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "nearest::checkNullCentroidsNN:  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "(output T_CLUSTERIDX: oidxT_cj = " << oidxT_cj << '\n'
	      << " output T_DIST: aort_distMinCentInst = " << aort_distMinCentInst << "]\n"
	      << " input  mat::MatrixRow<T_CENTROIDS>: aimatrixrowt_v[" <<  &aimatrixrowt_centroids << "]\n"
	      << " input  T_FEATURE *aiat_instance[" << aiat_instance << "]\n"
	      << " input  aifunc2p_dist[" << &aifunc2p_dist << "]\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/


  const T_FEATURE* larrayT_centroidj = aimatrixrowt_centroids.getRow(0);

  if ( data::Instance<T_FEATURE>::isInfiniteFeature(larrayT_centroidj[0]) ) {
    aort_distMinCentInst = std::numeric_limits<T_DIST>::max();
    oidxT_cj = NEARESTCENTROID_UNKNOWN;
  }
  else { 
    aort_distMinCentInst =
      aifunc2p_dist 
      (larrayT_centroidj,
       aiat_instance,
       aimatrixrowt_centroids.getNumColumns()
       );
    oidxT_cj = 0;
  }
  
  lidxT_K = T_CLUSTERIDX(aimatrixrowt_centroids.getNumRows());
  for (T_CLUSTERIDX li_j = 1; li_j < lidxT_K; li_j++) {
    const T_FEATURE* larrayT_centroidj = aimatrixrowt_centroids.getRow(li_j);
    lrt_distTiCi =
      (data::Instance<T_FEATURE>::isInfiniteFeature(larrayT_centroidj[0]))
      ?  std::numeric_limits<T_DIST>::max()
      : aifunc2p_dist
      (aimatrixrowt_centroids.getRow(li_j),
       aiat_instance, 
       aimatrixrowt_centroids.getNumColumns()
       );
    if ( aort_distMinCentInst > lrt_distTiCi ) {
      oidxT_cj = li_j;
      aort_distMinCentInst = lrt_distTiCi;
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "nearest::checkNullCentroidsNN: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "T_CLUSTERIDX: oidxT_cj = " << oidxT_cj << '\n'
	      << "T_DIST: aort_distMinCentInst = " << aort_distMinCentInst
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return oidxT_cj;
}


/*! \fn T_CLUSTERIDX medoidsNN(T_DIST_INSTANCESINSTANCES   &aort_distMinMedoidInst, const uintidx aiuidx_instance, const uintidx *aiarrayuidx_medoids, const T_CLUSTERIDX aicidx_numMedoisK, const mat::MatrixTriang<T_DIST_INSTANCESINSTANCES> &aimatrixtriagt_dissimilarity //DIST INSTANCES)
    \brief medoidsNN Find the medoid closest to an instance
    \details 
    \param aort_distMinMedoidInst output distance between medoid and instance 
    \param aiuidx_instance a index of the instance
    \param aiarrayuidx_medoids  a array with index of the medoids
    \param aimatrixtriagt_dissimilarity  a triangular matrix with distances between instances 
 */
template <typename T_CLUSTERIDX, 
	  typename T_DIST_INSTANCESINSTANCES
	  >
T_CLUSTERIDX
medoidsNN
(T_DIST_INSTANCESINSTANCES   &aort_distMinMedoidInst,
 const uintidx               aiuidx_instance,      //FOR THIS INSTANCE FIND NEAR MEDOIDS
 const uintidx               *aiarrayuidx_medoids,
 const T_CLUSTERIDX   aicidx_numMedoisK,
 const mat::MatrixTriang
 <T_DIST_INSTANCESINSTANCES> &aimatrixtriagt_dissimilarity //DIST INSTANCES
 )
{
  T_CLUSTERIDX        oidxt_cj;
  T_DIST_INSTANCESINSTANCES  lrt_distTiCi;
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "nearest::medoidsNN";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
  	      << "\t(output T_DIST: aort_distMinMedoidInst = " << aort_distMinMedoidInst << '\n'
	      << "\t input  uintidx: aiuidx_instance = " 
	      << aiuidx_instance << "]\n"
	      << "\t input  uintidx: aiarrayuidx_medoids[" 
	      <<  aiarrayuidx_medoids << "]\n"
	      << "\t input  T_CLUSTERIDX: aicidx_numMedoisK = " 
	      << aicidx_numMedoisK << "]\n"
              << "\t input  mat::MatrixTriang<T_DIST> aimatrixtriagt_dissimilarity[" 
	      << &aimatrixtriagt_dissimilarity << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  oidxt_cj = 0;
  aort_distMinMedoidInst = ( aiarrayuidx_medoids[0] > aiuidx_instance )?
    aimatrixtriagt_dissimilarity(aiarrayuidx_medoids[0],aiuidx_instance):
    aimatrixtriagt_dissimilarity(aiuidx_instance,aiarrayuidx_medoids[0]);
  for ( T_CLUSTERIDX li_j = 1; li_j < aicidx_numMedoisK; li_j++) {
    lrt_distTiCi = (aiarrayuidx_medoids[li_j] > aiuidx_instance)?
      aimatrixtriagt_dissimilarity(aiarrayuidx_medoids[li_j],aiuidx_instance):
      aimatrixtriagt_dissimilarity(aiuidx_instance,aiarrayuidx_medoids[li_j]);
    
    if ( lrt_distTiCi < aort_distMinMedoidInst ) {
      oidxt_cj = li_j;
      aort_distMinMedoidInst = lrt_distTiCi;
    }
    
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
  	      << " aort_distMinMedoidInst = " << aort_distMinMedoidInst
	      << " oidxt_cj = " << oidxt_cj
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return oidxt_cj;
  
}


/*! \fn T_CLUSTERIDX NN(T_DIST &aort_distMinCentInst, const mat::MatrixBase<T_FEATURE> &aimatrixrowt_centroids, T_CLUSTERIDX aimcIdx_currentCentroidK, const T_FEATURE *aiat_instance, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist) 
    \brief NN Find the centroid closest to an instance starting from a defined centroid
    \details 
    \param aort_distMinCentInst output distance between centroid and instance 
    \param aimatrixrowt_centroids matrix with the centroids
    \param aimcIdx_currentCentroidK a current centroid that the instance belongs to
    \param aiat_instance a array with the instance attributes
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX,
	  typename T_FEATURE,
	  typename T_DIST
	  >
T_CLUSTERIDX
NN
(T_DIST                             &aort_distMinCentInst,
 const mat::MatrixBase<T_FEATURE>   &aimatrixrowt_centroids,
 T_CLUSTERIDX                       aicidx_currentCentroidK,
 const T_FEATURE                    *aiat_instance,
 const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist
 )
{
  T_CLUSTERIDX oidxt_cj;
  T_CLUSTERIDX lcidx_numClustersK;
  T_DIST       lrt_distTiCi;

  lcidx_numClustersK = T_CLUSTERIDX(aimatrixrowt_centroids.getNumRows());
  if ( ( 0 > aicidx_currentCentroidK) || (aicidx_currentCentroidK >= lcidx_numClustersK) )
    aicidx_currentCentroidK = 0;

  oidxt_cj = aicidx_currentCentroidK;
  aort_distMinCentInst =  
    aifunc2p_dist 
    (aimatrixrowt_centroids.getRow(aicidx_currentCentroidK), 
     aiat_instance, 
     aimatrixrowt_centroids.getNumColumns()
     );
  
  for (T_CLUSTERIDX li_j = 0; li_j < lcidx_numClustersK; li_j++) {
    if ( li_j != aicidx_currentCentroidK ) {
      lrt_distTiCi = 
	aifunc2p_dist
	(aimatrixrowt_centroids.getRow(li_j),
	 aiat_instance,
	 aimatrixrowt_centroids.getNumColumns()
	 );
      if (aort_distMinCentInst > lrt_distTiCi) {
	oidxt_cj = li_j;
	aort_distMinCentInst = lrt_distTiCi;
      }
    }
  }

  return oidxt_cj;

}


/*! \fn  std::vector<T_DIST> dist(const mat::MatrixBase<T_FEATURE> &aimatrixrowt_centroids, const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstClusterK, const T_FEATURE *aiat_instance, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Calculate the distances of the instance to the centroids \cite Lu:etal:GAclusteringLabel:FGKA:2004
    \details
    \param aimatrixrowt_centroids a mat::MatrixBase with the centroids
    \param aivectort_numInstClusterK a std::vector with the number of instances
    \param aiat_instance a array with the instance attributes
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template<typename T_FEATURE,
	 typename T_INSTANCES_CLUSTER_K,
	 typename T_DIST
	 >
std::vector<T_DIST>
dist
(const mat::MatrixBase<T_FEATURE>         &aimatrixrowt_centroids,
 const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstClusterK,
 const T_FEATURE                          *aiat_instance,
 const dist::Dist<T_DIST,T_FEATURE>       &aifunc2p_dist
 )
{
  std::vector<T_DIST>  aovectorrt_distK(aimatrixrowt_centroids.getNumRows());
  
  for ( uintidx lui_j = 0; lui_j < aimatrixrowt_centroids.getNumRows(); lui_j++) {
    
    aovectorrt_distK.at(lui_j) =
      (aivectort_numInstClusterK.at(lui_j) != 0)?  
      aifunc2p_dist
      (aimatrixrowt_centroids.getRow(lui_j),
       aiat_instance,
       aimatrixrowt_centroids.getNumColumns()
       )
      :0.0;
      
  }

  return aovectorrt_distK;

}

} /*END namespace nearest*/

#endif /*__NEAREST_INSTANCES_OPERATOR_HPP*/
