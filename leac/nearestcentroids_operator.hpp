/*! \file nearestcentroids_operator.hpp
 *
 *
 * \brief  nearest centroids operator
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __NEAREST_CENTROIDS_OPERATOR_HPP
#define __NEAREST_CENTROIDS_OPERATOR_HPP

#include <iostream>
#include <vector>
#include "nearest_centroids.hpp"
#include "unsupervised_measures.hpp"
#include "dist.hpp"

/*! \namespace nearest
  \brief Function to find nearest instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace nearest {


/*! \fn T_CLUSTERIDX findTwoNearestCentroids(T_CLUSTERIDX aicidx_centroidCj, const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Find the centroid Cjp closest to the centroid Cj 
    \details
    \param aicidx_centroidCj an index of the centroid. From this look for the nearest
    \param aimatrixt_centroids a mat::MatrixRow with centroids
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template <typename T_CLUSTERIDX, 
	  typename T_FEATURE, 
	  typename T_DIST
	  >
T_CLUSTERIDX
findTwoNearestCentroids 
(T_CLUSTERIDX                     aicidx_centroidCj,
 const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids,
 dist::Dist<T_DIST,T_FEATURE>     &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "nearest::findTwoNearestCentroids:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(input T_CLUSTERIDX: aicidx_centroidCj = " << aicidx_centroidCj
	      << "\n\t input  mat::MatrixBase<>: aimatrixt_centroids["
	      <<  &aimatrixt_centroids << ']'
              << "\n\t input  dist::Dist<> aifunc2p_dist[" << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  T_CLUSTERIDX locidx_nearestCentroid; /*INDEX NEAREST CENTROID*/
  T_CLUSTERIDX lidxT_numCluster;

  locidx_nearestCentroid = NEARESTCENTROID_UNKNOWN;
  lidxT_numCluster = T_CLUSTERIDX(aimatrixt_centroids.getNumRows());

  if ( lidxT_numCluster > 1 ) {
    
    T_CLUSTERIDX lcidx_l = ( aicidx_centroidCj == 0 )?1:0; /*ITERATOR OF CENTROIDS*/
    T_DIST lrt_minDistCjCnearest =  
      aifunc2p_dist 
      (aimatrixt_centroids.getRow(aicidx_centroidCj), 
       aimatrixt_centroids.getRow(lcidx_l), 
       aimatrixt_centroids.getNumColumns()
       );
    locidx_nearestCentroid = lcidx_l;

    for ( ; lcidx_l < lidxT_numCluster; lcidx_l++) {
      if ( lcidx_l != aicidx_centroidCj) {
	T_DIST lrt_distCjCl = 
	  aifunc2p_dist 
	  (aimatrixt_centroids.getRow(aicidx_centroidCj), 
	   aimatrixt_centroids.getRow(lcidx_l), 
	   aimatrixt_centroids.getNumColumns()
	   );
	if ( lrt_distCjCl < lrt_minDistCjCnearest ) {
	  lrt_minDistCjCnearest = lrt_distCjCl;
	  locidx_nearestCentroid  = lcidx_l;
	}
	
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "nearest::findTwoNearestCentroids: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "T_CLUSTERIDX: locidx_nearestCentroid = " << locidx_nearestCentroid
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return locidx_nearestCentroid;
}

/*! \fn void findTwoNearestCentroids (NearestCentroids<T_CLUSTERIDX,T_DIST> &aonearcent_centroidCjp, mat::MatrixBase<T_CENTROIDS> &aimatrixt_centroids, std::vector<T_INSTANCES_CLUSTER_K> &aivectorit_numInstClusterK, T_CLUSTERIDX aicidx_centroidCj, const dist::Dist<T_DIST,T_CENTROIDS>  &aifunc2p_dist)
  \brief Find two nearest centroids
  \details
  \param aonearcent_centroidCjp a nearest::NearestCentroids stores the index and distance to the nearest centroid
  \param aimatrixt_centroids a mat::MatrixRow with centroids
  \param aivectorit_numInstClusterK a vector with the number of instances in each cluster
  \param aicidx_centroidCj an index of the centroid. From this look for the nearest
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CENTROIDS,
	   typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_DIST 
	   > 
void findTwoNearestCentroids
(NearestCentroids
 <T_CLUSTERIDX,
  T_DIST>                              &aonearcent_centroidCjp,
 mat::MatrixBase<T_CENTROIDS>          &aimatrixt_centroids,
 std::vector<T_INSTANCES_CLUSTER_K>    &aivectorit_numInstClusterK,
 T_CLUSTERIDX                          aicidx_centroidCj,
 const dist::Dist<T_DIST,T_CENTROIDS>  &aifunc2p_dist
 )
{  
  T_CLUSTERIDX     lcidx_iteratorClusterK;
  T_CLUSTERIDX     lcidx_numClusterK;
  T_DIST           lrt_distCjCi;
  
  lcidx_iteratorClusterK = (aicidx_centroidCj == 0)?1:0;
  lcidx_numClusterK  = (T_CLUSTERIDX) aimatrixt_centroids.getNumRows();
  
  lrt_distCjCi = 
    um::distortion
    (aimatrixt_centroids.getRow( (uintidx) aicidx_centroidCj),
     aivectorit_numInstClusterK[aicidx_centroidCj],
     aimatrixt_centroids.getRow( (uintidx) lcidx_iteratorClusterK),
     aivectorit_numInstClusterK[lcidx_iteratorClusterK],
     aifunc2p_dist
     );
  aonearcent_centroidCjp.setDistCentroidCentroid(lrt_distCjCi);
  aonearcent_centroidCjp.setNearestClusterK(lcidx_iteratorClusterK);
  
  for(lcidx_iteratorClusterK = lcidx_iteratorClusterK+1 ; 
      lcidx_iteratorClusterK < lcidx_numClusterK;  
      lcidx_iteratorClusterK++ ) {
    if (  aicidx_centroidCj != lcidx_iteratorClusterK ) {
      lrt_distCjCi = 
	um::distortion
	(aimatrixt_centroids.getRow( (uintidx) aicidx_centroidCj),
	 aivectorit_numInstClusterK[aicidx_centroidCj],
	 aimatrixt_centroids.getRow( (uintidx) lcidx_iteratorClusterK),
	 aivectorit_numInstClusterK[lcidx_iteratorClusterK],
	 aifunc2p_dist
	 );
     
      if ( lrt_distCjCi < aonearcent_centroidCjp.getDistCentroidCentroid() ) {
	aonearcent_centroidCjp.setDistCentroidCentroid(lrt_distCjCi);
	aonearcent_centroidCjp.setNearestClusterK(lcidx_iteratorClusterK);
	//aonearcent_centroidCjp.setDistRecalculate(false);
      }
    } 
  }
}

/*! \fn void findAllPairCentroidsClosest (std::vector<NearestCentroids<T_CLUSTERIDX, T_DIST> > &aovector_nearestCentroids, mat::MatrixBase<T_CENTROIDS> &aimatrixt_centroids, std::vector<T_INSTANCES_CLUSTER_K> &aivectorit_numInstClusterK, const dist::Dist<T_DIST,T_CENTROIDS> &aifunc2p_dist)
   \brief Find all pair of centroids closest
   \details The pairs of centroids are determined by the index of the vector and the index stored in nearest::NearestCentroids
   \param aovector_nearestCentroids a std::vector of nearest::NearestCentroids where the nearest pairs of centroids are recorded
   \param aimatrixt_centroids a mat::MatrixRow with centroids
   \param aivectorit_numInstClusterK a vector with numero number of instances in each cluster
   \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CENTROIDS,
	   typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_DIST 
	   > 
void 
findAllPairCentroidsClosest
(std::vector
 <NearestCentroids
  <T_CLUSTERIDX,
   T_DIST> >                          &aovector_nearestCentroids,
 mat::MatrixBase<T_CENTROIDS>         &aimatrixt_centroids,
 std::vector<T_INSTANCES_CLUSTER_K>   &aivectorit_numInstClusterK,
 const dist::Dist<T_DIST,T_CENTROIDS> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "clusteringop::findAllPairCentroidsClosest:  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\t(output std::vector<NearestCentroids>: aovector_nearestCentroids[" 
		<< &aovector_nearestCentroids << "]\n"
		<< "\t input mat::MatrixRow<T_CENTROIDS>: aimatrixt_centroids["
		<< &aimatrixt_centroids << "]\n"
		<< "\t input  aifunc2p_dist\n"
		<< "\t)\n";
    }
#endif /*__VERBOSE_YES*/

    T_CLUSTERIDX        lcidx_j = T_CLUSTERIDX(0); 
    for ( auto & lnearestCentroidsDist_iter :  aovector_nearestCentroids ) {
    findTwoNearestCentroids
      (lnearestCentroidsDist_iter,
	aimatrixt_centroids,
	aivectorit_numInstClusterK,
	lcidx_j,
	aifunc2p_dist
	);
     lnearestCentroidsDist_iter.setDistRecalculate(false);
    ++lcidx_j;
  }

 #ifdef __VERBOSE_YES
     if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       std::cout << "clusteringop::findAllPairCentroidsClosest: OUT"
		 << '(' << geiinparam_verbose << ")\n"
		 << "\t(output std::vector<NearestCentroids>: aovector_nearestCentroids[" 
		 << &aovector_nearestCentroids << "]:\n"
		 << "[";

       for ( uintidx _lui_k = 0; _lui_k < aovector_nearestCentroids.size(); _lui_k++) {
	 std::cout << ' ' << _lui_k << ':'; 
	 aovector_nearestCentroids.at(_lui_k).print();
       }
       std::cout << ']' <<  std::endl;
     }
     --geiinparam_verbose;
 #endif /*__VERBOSE_YES*/

}

/*! \fn  void recalculateDistNearestCentroids (std::vector<NearestCentroids<T_CLUSTERIDX,T_DIST> > &aovector_nearestCentroids, mat::MatrixBase<T_CENTROIDS> &aimatrixt_centroids, std::vector<T_INSTANCES_CLUSTER_K> &aivectorit_numInstClusterK, const dist::Dist<T_DIST,T_CENTROIDS> &aifunc2p_dist)
  \brief Recalculates distance between nearest::NearestCentroids
  \details
  \param aovector_nearestCentroids a std::vector of nearest::NearestCentroids where the nearest pairs of centroids are recorded
  \param aimatrixt_centroids a mat::MatrixRow with centroids
  \param aivectorit_numInstClusterK a vector with numero number of instances in each cluster
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CENTROIDS,
	   typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_DIST 
	   > 
void 
recalculateDistNearestCentroids
(std::vector
 <NearestCentroids
  <T_CLUSTERIDX,
   T_DIST> >                          &aovector_nearestCentroids,
 mat::MatrixBase<T_CENTROIDS>         &aimatrixt_centroids,
 std::vector<T_INSTANCES_CLUSTER_K>   &aivectorit_numInstClusterK,
 const dist::Dist<T_DIST,T_CENTROIDS> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax  ) {
      std::cout << "clusteringop::recalculateDistNearestCentroids:  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "(<std::vector<NearestTypeFrant>: aovector_nearestCentroids[" 
		<< &aovector_nearestCentroids << "]\n"
		<< ">";
      for ( uintidx _lui_k = 0; _lui_k < aovector_nearestCentroids.size(); _lui_k++) {
	 std::cout << ' ' << _lui_k << ':'; 
	 aovector_nearestCentroids.at(_lui_k).print();
       } 
      std::cout << "\n input  mat::MatrixRow<T_CENTROIDS>: aimatrixt_centroids["
		<<  &aimatrixt_centroids << "]\n"
		<< " input  T_DIST_CENTROIDCENTROID: aifunc2p_dist\n"
		<< ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
  
  T_CLUSTERIDX lcidx_j = T_CLUSTERIDX(0);
 
  for ( auto & lnearestCentroidsDist_iter :  aovector_nearestCentroids ) {
    if ( lnearestCentroidsDist_iter.getDistRecalculate() ) {
      findTwoNearestCentroids
	(lnearestCentroidsDist_iter,
	 aimatrixt_centroids,
	 aivectorit_numInstClusterK,
	 lcidx_j,
	 aifunc2p_dist
	 );
      lnearestCentroidsDist_iter.setDistRecalculate(false);
    }
    ++lcidx_j;
  }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "clusteringop::recalculateDistNearestCentroids: OUT"
		<< '(' << geiinparam_verbose << ")\n"
		<< "<std::vector<NearestTypeFrant>: aovector_nearestCentroids[" 
		<< &aovector_nearestCentroids << "]\n"
		<< ">";
      for ( uintidx _lui_k = 0; _lui_k < aovector_nearestCentroids.size(); _lui_k++) {
	 std::cout << ' ' << _lui_k << ':'; 
	 aovector_nearestCentroids.at(_lui_k).print();
       } 
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

/*! \fn  std::pair<T_CLUSTERIDX,T_CLUSTERIDX> findNearestCentroids (const std::vector<NearestCentroids<T_CLUSTERIDX,T_DIST> >   &aivector_nearestCentroidsDist)
 \brief Find the nearest pair of centroids
 \details
 \param aivector_nearestCentroidsDist a std::vector of nearest::NearestCentroids with the nearest pairs
*/
template < typename T_CLUSTERIDX,
	   typename T_DIST
	   > 
std::pair<T_CLUSTERIDX,T_CLUSTERIDX>
findNearestCentroids
(const std::vector
 <NearestCentroids
  <T_CLUSTERIDX,
   T_DIST> >        &aivector_nearestCentroidsDist
 )
{
  std::pair<T_CLUSTERIDX,T_CLUSTERIDX> lopair_merged;

  T_DIST lT_minDistCjCi =
    aivector_nearestCentroidsDist.at(0).getDistCentroidCentroid();
  lopair_merged.first  = T_CLUSTERIDX(0);
  lopair_merged.second = aivector_nearestCentroidsDist.at(0).getNearestClusterK();
  
  for (uintidx luintidx_i = 1; luintidx_i < aivector_nearestCentroidsDist.size(); luintidx_i++) {
    if (aivector_nearestCentroidsDist.at(luintidx_i).getDistCentroidCentroid() < lT_minDistCjCi ) {
      lT_minDistCjCi = 
	aivector_nearestCentroidsDist.at(luintidx_i).getDistCentroidCentroid();
      lopair_merged.first  = T_CLUSTERIDX(luintidx_i);
      lopair_merged.second = 
	aivector_nearestCentroidsDist.at(luintidx_i).getNearestClusterK();
    }
  } 
  if ( lopair_merged.first > lopair_merged.second  ) { 
    T_CLUSTERIDX  lcidx_tmp = lopair_merged.first;
    lopair_merged.first = lopair_merged.second;
    lopair_merged.second = lcidx_tmp;
  }

  return lopair_merged;
}


/*! \fn  void markRecalculationNearest (std::vector<NearestCentroids<T_CLUSTERIDX,T_DIST> > &aovector_nearestCentroids, const std::pair<T_CLUSTERIDX,T_CLUSTERIDX> &aipair_merged)
 \brief Mark pairs of centroids that depend on the pair to be joined 
 \details
 \param aovector_nearestCentroids a std::vector of nearest::NearestCentroids with the nearest pairs
 \param aipair_merged a std::pair with cluster indexes to Join
*/
template < typename T_CLUSTERIDX,
	   typename T_DIST
	   >
void
markRecalculationNearest
(std::vector
 <NearestCentroids
 <T_CLUSTERIDX,
 T_DIST> >           &aovector_nearestCentroids,
 const std::pair
 <T_CLUSTERIDX,
 T_CLUSTERIDX>       &aipair_merged
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::markRecalculationNearest";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output std::vector<NearestCentroids>: *aovector_nearestCentroids["  
	      << &aovector_nearestCentroids << "]\n"
	      << " input  std::pair<>: aipair_merged("
	      << aipair_merged.first << ',' << aipair_merged.second << ")\n"
	      << ")\n";
  }
#endif //__VERBOSE_YES
  
  for ( auto & lnearestCentroidsDist_iter :  aovector_nearestCentroids ) {
    lnearestCentroidsDist_iter.setDistRecalculate
      ( (lnearestCentroidsDist_iter.getNearestClusterK() == aipair_merged.first ||
       lnearestCentroidsDist_iter.getNearestClusterK() == aipair_merged.second )?
      true:
      false
	);
  }
  aovector_nearestCentroids.at(aipair_merged.first).setDistRecalculate(true);

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "<aovector_nearestCentroids:" <<  lpc_labelFunc
	      << ":size" << aovector_nearestCentroids.size() << '>';
    
    for ( uintidx _lui_k = 0; _lui_k < aovector_nearestCentroids.size(); _lui_k++) {
	 std::cout << ' ' << _lui_k << ':'; 
	 aovector_nearestCentroids.at(_lui_k).print();
    }	    
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
}

/*! \fn void deleteNearestCentroid(std::vector<NearestCentroids<T_CLUSTERIDX,T_DIST> > &aovector_nearestCentroids, const T_CLUSTERIDX aicidx_clusterToDelete)
 \brief Delete a nearest::NearestCentroids 
 \details
 \param aovector_nearestCentroids a std::vector of nearest::NearestCentroids with the nearest pairs
 \param aicidx_clusterToDelete cluster index to delete 
*/  
template < typename T_CLUSTERIDX,
	   typename T_DIST
	   >
void deleteNearestCentroid  
(std::vector
 <NearestCentroids
 <T_CLUSTERIDX,
 T_DIST> >               &aovector_nearestCentroids,
 const T_CLUSTERIDX      aicidx_clusterToDelete
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "nearest::deleteNearestCentroid";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n";
  }
#endif /*__VERBOSE_YES*/

  T_CLUSTERIDX lidxK_idxClusterKLast =
    (T_CLUSTERIDX) aovector_nearestCentroids.size() - 1;

  if ( aicidx_clusterToDelete != lidxK_idxClusterKLast ) {

    std::swap
      (aovector_nearestCentroids[lidxK_idxClusterKLast],
       aovector_nearestCentroids[aicidx_clusterToDelete]
       );
  
    for ( auto & lnearestCentroidsDist_iter :  aovector_nearestCentroids ) {
      if (lnearestCentroidsDist_iter.getNearestClusterK() == lidxK_idxClusterKLast )
	lnearestCentroidsDist_iter.setNearestClusterK( aicidx_clusterToDelete );
    }
   
  }
  aovector_nearestCentroids.pop_back();

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*End deleteNearestCentroid*/

} /*END namespace nearest*/

#endif /* __NEAREST_CENTROIDS_OPERATOR_HPP */
