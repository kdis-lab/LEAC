/*! \file unsupervised_measures.hpp
 *
 * \brief Unsupervised measures
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __UNSUPERVISED_MEASURE_HPP
#define __UNSUPERVISED_MEASURE_HPP

#include <algorithm>    // std::find
#include <iterator>     // std::distance
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <limits>
#include <cmath>        //std::isnan 
#include <math.h>
#include "matrix.hpp"
#include "bit_matrix.hpp"
#include "partition.hpp"
#include "partition_linked_numinst.hpp"
#include "vector_utils.hpp"
#include "container_out.hpp"
#include "leac_utils.hpp"
#include "stats_instances.hpp"

#include "verbose_global.hpp"


/*! \namespace um
  \brief  Unsupervised measures for clustering analysis
  \details Evaluate tries to determine the quality of a given obtained partition of the data without any external information available. 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace um {


/*When the algorithms do not end with a correct value, for example 
  cluster number k = 0 or k = 1, metrics are assigned to special 
  values to each according to his definition
*/
#define measuare_undefObjetiveFunc(T_METRIC) (T_METRIC) -1.0

#define measuare_undefSSE(T_METRIC)    std::numeric_limits<T_METRIC>::max()
#define measuare_undefIndexI(T_METRIC) (T_METRIC) 0.0
#define measuare_undefVRC(T_METRIC)    (T_METRIC) 0.0
#define measuare_undefCS(T_METRIC)      std::numeric_limits<T_METRIC>::max()
#define measuare_undefDunnIndex(T_METRIC) (T_METRIC) 0.0
#define measuare_undefDBindex(T_METRIC) std::numeric_limits<T_METRIC>::max()
#define measuare_undefSilhouette(T_METRIC) (T_METRIC) 0.0 // -1.0
#define measuare_lowerValueSilhouette(T_METRIC) (T_METRIC) -1.0

#define measuare_undefPurity(T_METRIC) (T_METRIC) 0.0
#define measuare_undefVRC(T_METRIC) (T_METRIC) 0.0
#define measuare_undefCS(T_METRIC) std::numeric_limits<T_METRIC>::max()

#define measuare_undefScoreFunction(T_METRIC) (T_METRIC) 0.0
#define measuare_undefXieBeniIndex(T_METRIC) std::numeric_limits<T_METRIC>::max()
#define measuare_undefEntropy(T_METRIC) (T_METRIC) std::numeric_limits<T_METRIC>::max()
#define measuare_undefPartitionCoefficient(T_METRIC) (T_METRIC) 0.0

#define measuare_undefWBIndex(T_METRIC) (T_METRIC) std::numeric_limits<T_METRIC>::max()

#define measuare_undefOverlap(T_METRIC) (T_METRIC) -std::numeric_limits<T_METRIC>::max()
  
/*! \fn T_METRIC maxDistCjCjp(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Find the maximum distance beteen the centroids
  \details 
  \param aimatrixt_centroids a mat::MatrixRow with centroids
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_METRIC  
	   >
T_METRIC 
maxDistCjCjp
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::maxDistCjCjp";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(\n\t input  mat::MatrixRow<T_FEATURE>&: aimatrixt_centroids[" 
	      << &aimatrixt_centroids << ']'
	      << "\n\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lot_DkMax = T_METRIC(0);
 
  for(uintidx lui_i = 0; lui_i < aimatrixt_centroids.getNumRows(); lui_i++) {
    for(uintidx lui_j = lui_i+1; lui_j < aimatrixt_centroids.getNumRows(); lui_j++) {
     
      T_METRIC  lot_distSiSj = 
	aifunc2p_dist
	(aimatrixt_centroids.getRow(lui_i),
	 aimatrixt_centroids.getRow(lui_j),
	 aimatrixt_centroids.getNumColumns()
	 );
      if ( lot_DkMax < lot_distSiSj ) 
	lot_DkMax = lot_distSiSj;
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " T_METRIC lot_DkMax = " << lot_DkMax
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lot_DkMax;

}


/*! \fn T_METRIC minDistCjCjp(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Find the minimum distance between centroids
  \details
  \f[
  \begin{array}{c} min\\j \in k, j \neq j' \end{array}  D(\mu_j,\mu_{j'}),
  \f]    
  where \f$D\f$ is a distance function
  \para aimatrixt_centroids a mat::MatrixRow with centroids
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_METRIC  
	   >
T_METRIC
minDistCjCjp
(const mat::MatrixRow<T_FEATURE>         &aimatrixt_centroids,
 const dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::minDistCjCjp";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t(input  mat::MatrixRow<T_FEATURE>&: aimatrixt_centroids[" << &aimatrixt_centroids << ']'
      << "\n\tinput  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist["  << &aifunc2p_dist << ']'
      << "\n\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lot_minDistCjCjp   = std::numeric_limits<T_METRIC>::max();
  
  for(uintidx lui_i = 0; lui_i < aimatrixt_centroids.getNumRows(); lui_i++) {
    for(uintidx lui_j = lui_i+1; lui_j < aimatrixt_centroids.getNumRows(); lui_j++) {
   
      T_METRIC  lot_distSiSj = 
	aifunc2p_dist
	(aimatrixt_centroids.getRow(lui_i),
	 aimatrixt_centroids.getRow(lui_j),
	 aimatrixt_centroids.getNumColumns()
	 );
      if ( lot_minDistCjCjp > lot_distSiSj ) 
	lot_minDistCjCjp = lot_distSiSj;
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " T_METRIC lot_minDistCjCjp = " << lot_minDistCjCjp
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lot_minDistCjCjp;

}


/*! \fn T_METRIC minDistCjCjp(const T_CLUSTERIDX aicidx_Clusterj, const mat::MatrixRow<T_FEATURE>   &aimatrixt_centroids, dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist) 
  \brief Find the minimum distance between aicidx_Clusterj and all centroids,
  \details
  \f[
  \begin{array}{c} min\\j \in k, j \neq j' \end{array}  D(\mu_j,\mu_{j'}),
  \f]    
  where  \f$D\f$ is a distance function
  \param aicidx_Clusterj a const integer centroid index 
  \param aimatrixt_centroids a mat::MatrixBase with centroids
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,
           typename T_FEATURE,
	   typename T_METRIC  
	   >
T_METRIC
minDistCjCjp
(const  T_CLUSTERIDX                   aicidx_Clusterj,
 const mat::MatrixBase<T_FEATURE>      &aimatrixt_centroids,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::minDistCjCjp";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t( input  mat::MatrixRow<T_FEATURE>&: aimatrixt_centroids[" 
	      << &aimatrixt_centroids << ']'
	      << "\n\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lot_minDistCjCjp   = std::numeric_limits<T_METRIC>::max();

  uintidx lui_j =  uintidx(aicidx_Clusterj);
  
  for(uintidx lui_jp = 0; lui_jp < aimatrixt_centroids.getNumRows(); lui_jp++) {
    
    if ( lui_j  != lui_jp ) {
      T_METRIC  lot_distSiSj = 
	aifunc2p_dist
	(aimatrixt_centroids.getRow(lui_jp),
	 aimatrixt_centroids.getRow(lui_j),
	 aimatrixt_centroids.getNumColumns()
	 );
      if ( lot_minDistCjCjp > lot_distSiSj ) 
	lot_minDistCjCjp = lot_distSiSj;
    }
   
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " T_METRIC lot_minDistCjCjp = " << lot_minDistCjCjp
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lot_minDistCjCjp;

}
  
/*! \fn T_METRIC sumMinCjCjp(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist) 
  \brief Calculate minimum sum of distances between all centroids
  \details
  \f[
  \sum_{j=1}^{k} \left\{ \begin{array}{c} min\\j \in k, j \neq j' \end{array}  \left\{ D(\mu_j,\mu_{j'}) \right\}\
  \right\}
  \f]    
  where \f$D\f$ is a distance function
  \param aimatrixt_centroids a mat::MatrixRow with centroids
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_METRIC  
	   >
T_METRIC
sumMinCjCjp
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{ 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::sumMinCjCjp";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\n\t input  mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lot_sumMinCjCjp = T_METRIC(0.0);
  
  for ( uintidx luintidx_Ci = 0; luintidx_Ci < aimatrixt_centroids.getNumRows(); luintidx_Ci++)
    {
      T_METRIC lt_minDistCjCjp = std::numeric_limits<T_METRIC>::max();  
      for ( uintidx luintidx_Cj = 0; luintidx_Cj < aimatrixt_centroids.getNumRows(); luintidx_Cj++)
	{
	  if ( luintidx_Ci != luintidx_Cj ) {
	 
	    T_METRIC lt_distViVj =
	      aifunc2p_dist
	      (aimatrixt_centroids.getRow(luintidx_Ci), 
	       aimatrixt_centroids.getRow(luintidx_Cj),
	       aimatrixt_centroids.getNumColumns()
	       );     
	
	    if ( lt_distViVj < lt_minDistCjCjp )
	      lt_minDistCjCjp = lt_distViVj;
	  }
	}
      if ( lt_minDistCjCjp < std::numeric_limits<T_METRIC>::max() ) {
	lot_sumMinCjCjp += (T_METRIC) lt_minDistCjCjp;
      }
    }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lot_sumMinCjCjp = " << lot_sumMinCjCjp
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lot_sumMinCjCjp;
  
} //END sumMinCjCjp

/*! \fn std::vector<T_METRIC> avgRadiusClusterK(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Calculate the average radius of clusters
  \details
  \param aimatrixt_centroids a matrix with centroids
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
std::vector<T_METRIC>
avgRadiusClusterK
(const mat::MatrixBase<T_FEATURE>     &aimatrixt_centroids,
 INPUT_ITERATOR                       aiiterator_instfirst,
 const INPUT_ITERATOR                 aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>   &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::avgRadiusClusterK";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<T_METRIC> lovectort_sumDistXZi
    (aimatrixt_centroids.getNumRows(),T_METRIC(0));

  std::vector<uintidx> lovectoui_numInstClusterK
    (aimatrixt_centroids.getNumRows(),T_METRIC(0));

  if ( aimatrixt_centroids.getNumRows() >= 1 ) {
    const T_CLUSTERIDX  lcidx_numClusterK = 
      T_CLUSTERIDX(aimatrixt_centroids.getNumRows());

    for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
      T_CLUSTERIDX lcidx_xinK = aipartition_clusters.next();
      
      if ( 0 <= lcidx_xinK  && lcidx_xinK <  lcidx_numClusterK  ) {

	data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
	
	lovectort_sumDistXZi[lcidx_xinK] += 
	  aifunc2p_dist
	  (aimatrixt_centroids.getRow(lcidx_xinK),
	   linst_inter->getFeatures(),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	++lovectoui_numInstClusterK[lcidx_xinK];
      }
    }

    for ( uintidx _lui_k = 0; _lui_k < lovectort_sumDistXZi.size(); _lui_k++)
      {
	if ( lovectoui_numInstClusterK[_lui_k] != 0 ){
	  lovectort_sumDistXZi[_lui_k] /= (T_METRIC) lovectoui_numInstClusterK[_lui_k];
	}
	else {
	  lovectort_sumDistXZi[_lui_k] = std::numeric_limits<T_METRIC>::quiet_NaN();
	}
      }
    
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelRadium;
    lostrstream_labelRadium << "<AVGRADADIUSCLUSTERK:" << lpc_labelFunc;
    inout::containerprint
      (lovectort_sumDistXZi.begin(),
       lovectort_sumDistXZi.end(),
       std::cout,
       lostrstream_labelRadium.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectort_sumDistXZi;
    
}

  
/*! \fn std::vector<T_METRIC> maxRadiusClusterK(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist)
  \brief Calculates the maximum radius of a cluster for a partition 
  \details
  \param aimatrixt_centroids a mat::MatrixRow with centroids
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
std::vector<T_METRIC>
maxRadiusClusterK
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::maxRadiusClusterK";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<T_METRIC> lovectort_maxRadiusClusterK
    (aimatrixt_centroids.getNumRows(),T_METRIC(0));

  if ( aimatrixt_centroids.getNumRows() >= 1 ) {
    const T_CLUSTERIDX  lcidx_numClusterK = 
      T_CLUSTERIDX(aimatrixt_centroids.getNumRows());

    for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
      T_CLUSTERIDX lcidx_xinK =
	aipartition_clusters.next();
      
      if ( 0 <= lcidx_xinK  && lcidx_xinK <  lcidx_numClusterK  ) {

	data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
	
	T_METRIC lt_distInstCentroid = 
	  aifunc2p_dist
	  (aimatrixt_centroids.getRow(lcidx_xinK),
	   linst_inter->getFeatures(),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	if ( lovectort_maxRadiusClusterK[lcidx_xinK] <  lt_distInstCentroid )
	  lovectort_maxRadiusClusterK[lcidx_xinK] = lt_distInstCentroid;
	  
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelRadium;
    lostrstream_labelRadium << "<MAXRADIUSCLUSTERK: " << lpc_labelFunc
			    << "lovectort_maxRadiusClusterK[" << &lovectort_maxRadiusClusterK << ']';
    inout::containerprint
      (lovectort_maxRadiusClusterK.begin(),
       lovectort_maxRadiusClusterK.end(),
       std::cout,
       lostrstream_labelRadium.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectort_maxRadiusClusterK;
    
}

/*! \fn  T_METRIC diameterClusterKj(const T_CLUSTERIDX aicidx_Clusterj, ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, INPUT_ITERATOR aiiterator_instfirst, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist) 
  \brief  Diameter of Clusterj
  \details Find the diameter of cluster j, the maximum distance between two instances of the cluster 
  \param aicidx_Clusterj an integer index of the cluster
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,
           typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_METRIC
	   >
T_METRIC
diameterClusterKj
(const T_CLUSTERIDX                      aicidx_Clusterj,
 const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 INPUT_ITERATOR                          aiiterator_instfirst,
 const dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::diameterClusterKj";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t( input const T_CLUSTERIDX  aicidx_Clusterj[" << aicidx_Clusterj << ']'
      << "\t  input PartitionLinked&: aipartlink_memberShip[" << &aipartlink_memberShip << "]\n"
      << "\t  input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t  input dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" << &aifunc2p_dist << "]\n"
      << "\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_j(&aipartlink_memberShip);
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_k(&aipartlink_memberShip);
  T_METRIC  lot_delta = T_METRIC(0);

  for ( literpart_j.begin(aicidx_Clusterj); literpart_j.end(); literpart_j.next() ) {
    for ( literpart_k.begin(aicidx_Clusterj); literpart_k.end(); literpart_k.next() ) {
      data::Instance<T_FEATURE>* linst_j = *std::next(aiiterator_instfirst,literpart_j.getValue());
      data::Instance<T_FEATURE>* linst_k = *std::next(aiiterator_instfirst,literpart_k.getValue());
      T_METRIC lt_distk =
	aifunc2p_dist
	(linst_j->getFeatures(),
	 linst_k->getFeatures(),
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
      if ( lot_delta < lt_distk )
	lot_delta = lt_distk;
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " T_METRIC lot_delta = " << lot_delta
      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lot_delta;
}


/*! \fn  T_METRIC diameterClusterKj(const T_CLUSTERIDX aicidx_Clusterj, const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity) 
  \brief  Diameter of Clusterj
  \details Find the diameter of cluster j, the maximum distance between two instances of the cluster 
  \param aicidx_Clusterj an integer index of the cluster
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aimatrixtriagt_dissimilarity a matrix of distances
*/
template < typename T_CLUSTERIDX,
	   typename T_METRIC
	   >
T_METRIC
diameterClusterKj
(const T_CLUSTERIDX                      aicidx_Clusterj,
 const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 const mat::MatrixTriang<T_METRIC>       &aimatrixtriagt_dissimilarity
) 
{

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::diameterClusterKj";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t(input const T_CLUSTERIDX  aicidx_Clusterj[" << aicidx_Clusterj << "]\n"
      << "\t input PartitionLinked&: aipartlink_memberShip[" << &aipartlink_memberShip << "]\n"
      << "\t input  mat::MatrixTriang<T_METRIC>: aimatrixtriagt_dissimilarity[" 
      << &aimatrixtriagt_dissimilarity << "]\n"
      << "\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES
    
 ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_j(&aipartlink_memberShip);
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_k(&aipartlink_memberShip);
  
  T_METRIC  lot_delta = T_METRIC(0);

  for ( literpart_j.begin(aicidx_Clusterj); literpart_j.end(); literpart_j.next() ) {
    for ( literpart_k.begin(aicidx_Clusterj); literpart_k.end(); literpart_k.next() ) {
      T_METRIC lt_distk =
	aimatrixtriagt_dissimilarity(literpart_j.getValue(),literpart_k.getValue());
      if ( lot_delta < lt_distk )
	lot_delta = lt_distk;
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " T_METRIC lot_delta = " << lot_delta 
      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lot_delta;
}


/*! \fn  T_METRIC radiusClusterKj (const T_CLUSTERIDX aicidx_Clusterj, const T_FEATURE* ait_centroid, const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,INPUT_ITERATOR aiiterator_instfirst, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief  Radius of Clusterj
  \details Find the radius of cluster j, the maximum distance between an instance of the cluster and the centroid 
  \param aicidx_Clusterj an integer index of the cluster
  \param ait_centroid the centroid of the Clusterj
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,
           typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_METRIC
	   >
T_METRIC
radiusClusterKj
(const T_CLUSTERIDX                      aicidx_Clusterj,
 const T_FEATURE*                        ait_centroid, 
 const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 INPUT_ITERATOR                          aiiterator_instfirst,
 const dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::radiusClusterKj";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t(input  const T_CLUSTERIDX aicidx_Clusterj[" << aicidx_Clusterj << ']'
      << "\n\t const T_FEATURE* ait_centroid,[" << ait_centroid << ']'    
      << "\n\t input ds::PartitionLinked<T_CLUSTERIDX>&: aipartlink_memberShip[" 
      << &aipartlink_memberShip << ']'
      << "\n\t input aiiterator_instfirst[" << *aiiterator_instfirst << ']'
      << "\n\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist["  << &aifunc2p_dist << ']'
      << "\n\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_j(&aipartlink_memberShip);

  T_METRIC  lot_delta = T_METRIC(0);

  for ( literpart_j.begin(aicidx_Clusterj); literpart_j.end(); literpart_j.next() ) {
    
    data::Instance<T_FEATURE>* linst_j = *std::next(aiiterator_instfirst,literpart_j.getValue());

    T_METRIC lt_distk =
      aifunc2p_dist
      (linst_j->getFeatures(),
       ait_centroid,
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    if ( lot_delta < lt_distk )
      lot_delta = lt_distk;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " T_METRIC lot_delta = " << lot_delta
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lot_delta;
}



/*! \fn T_METRIC D2(const mat::MatrixRow<T_FEATURE> &aimatrixrowt_S, const partition::Partition<T_CLUSTERIDX> &aipartition_clustersBkinCi, const partition::Partition<T_CLUSTERIDX> &aipartition_clustersBi, const mat::MatrixRow<T_FEATURE> &aimatrixrowt_Vi, const std::vector<T_METRIC> &aivectorr_meanRadiusBk, const std::vector<data::Instance<T_FEATURE>* >  &aivectorptinst_instances, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief D2 \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
  \details
  \f[
  D_{2}(w) = \begin{array}{c}max\\ 1 \leq i \leq k \end{array}
  \begin{array}{c}max\\ B_l \subset C_j \end{array}
  \frac{\sum_{x_i \in C_j} \| x_i - S_j \| / | C_j |}
  {\sum_{x_i \in B_l} \| x_i - V_i \| / | B_l |}
  \f]    
  where   \f$D\f$ is a distance function
  \param aimatrixrowt_S 
  \param aipartition_clustersBkinCi 
  \param aipartition_clustersBi
  \param aimatrixrowt_Vi
  \param aivectorr_meanRadiusBk
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_METRIC,
	   typename T_FEATURE,
           typename T_CLUSTERIDX,
	   typename INPUT_ITERATOR
	   >
T_METRIC
D2
(const mat::MatrixRow<T_FEATURE>          &aimatrixrowt_S,     //Ci
 const partition::Partition<T_CLUSTERIDX> &aipartition_clustersBkinCi, //Centroids
 const partition::Partition<T_CLUSTERIDX> &aipartition_clustersBi,      //OjinVk
 const mat::MatrixRow<T_FEATURE>          &aimatrixrowt_Vi,
 const std::vector<T_METRIC>              &aivectorr_meanRadiusBk,
 INPUT_ITERATOR                           aiiterator_instfirst,
 const INPUT_ITERATOR                     aiiterator_instlast,
 const dist::Dist<T_METRIC,T_FEATURE>     &aifunc2p_dist
 )
{
  std::vector<T_METRIC>   lvectort_sumDistOjSi(aimatrixrowt_S.getNumRows(),0);
  std::vector<uintidx>    lvectort_numOjinCi(aimatrixrowt_S.getNumRows(),uintidx(0));
 

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::D2";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n"   
 	      << "(input  mat::MatrixRow<T_FEATURE>&: aimatrixrowt_S[" 
	      << &aimatrixrowt_S << ']' << " NumRows: " << aimatrixrowt_S.getNumRows()
	      << "\ninput  partition::Partition<>&: aipartition_clustersBkinCi" 
	      << &aipartition_clustersBkinCi << "]\n"
	      << "input  partition::Partition<>&: aipartition_clustersBi["
	      << &aipartition_clustersBi  << ']'
	      << "\nmat::MatrixRow<T_FEATURE>: aimatrixrowt_Vi["
	      << &aimatrixrowt_Vi << ']' << " NumRows: " << aimatrixrowt_Vi.getNumRows() << '\n'
	      << "input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  uintidx  luintidx_Oj = 0;
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst, luintidx_Oj++) {
    T_CLUSTERIDX luintidx_OjinBk =
      aipartition_clustersBi.getClusterIdx(luintidx_Oj);
    T_CLUSTERIDX luintidx_OjinCi =
      aipartition_clustersBkinCi.getClusterIdx(luintidx_OjinBk);
    
    data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;

    lvectort_sumDistOjSi.at(luintidx_OjinCi) +=
      aifunc2p_dist
      (aimatrixrowt_S.getRow(luintidx_OjinCi), 
       linst_inter->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    ++lvectort_numOjinCi.at(luintidx_OjinCi);
    
  }

  /*The mean radius of Ci
   */  
  std::vector<T_METRIC>  lvectorr_meanRadiusCi(aimatrixrowt_S.getNumRows());
  for ( uintidx  luintidx_Ci = 0; luintidx_Ci < aimatrixrowt_S.getNumRows(); luintidx_Ci++) {
    lvectorr_meanRadiusCi[luintidx_Ci] =
      ((T_METRIC) lvectort_sumDistOjSi[luintidx_Ci] / (T_METRIC) lvectort_numOjinCi[luintidx_Ci]);
  }
  
    
  T_METRIC lot_maxD2 = 0; 
  
  for ( uintidx  luintidx_Bk = 0; luintidx_Bk < aimatrixrowt_Vi.getNumRows(); luintidx_Bk++) {
    
    T_CLUSTERIDX luintidx_BkinCi =
      aipartition_clustersBkinCi.getClusterIdx(luintidx_Bk);

    T_METRIC  lt_quotientRadius = (aivectorr_meanRadiusBk[luintidx_Bk]>0)?
      lvectorr_meanRadiusCi[luintidx_BkinCi] / aivectorr_meanRadiusBk[luintidx_Bk]:0.0;
    
    if ( lt_quotientRadius > lot_maxD2 )
      lot_maxD2 = lt_quotientRadius;
     
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')' <<  " lot_maxD2 = " << lot_maxD2 << '\n';

    std::ostringstream lostrstream_labelmeanRadiusCi;
    lostrstream_labelmeanRadiusCi << "<MEANRADIUSCi:" << lpc_labelFunc;

    inout::containerprint
      (lvectorr_meanRadiusCi.begin(),
       lvectorr_meanRadiusCi.end(),
       std::cout,
       lostrstream_labelmeanRadiusCi.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lot_maxD2;

}


/*! \fn T_METRIC Dintra(const mat::MatrixRow<T_FEATURE> &aimatrixrowt_S, const mat::MatrixRow<T_FEATURE> &aimatrixrowt_Vi, const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstBi, const partition::Partition<T_CLUSTERIDX> &aipartition_clustersBkinCi, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Dintra \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
  \param aimatrixrowt_S a matrix with centroid S mat::MatrixRow<T_FEATURE>
  \param aimatrixrowt_Vi a matrix with centroid Vi 
  \param aivectort_numInstBi a vector with numero number of instances in Bi
  \param aipartition_clustersBkinCi a partition of B in the groups Ci
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
Dintra
(const mat::MatrixRow<T_FEATURE>           &aimatrixrowt_S,
 const mat::MatrixRow<T_FEATURE>           &aimatrixrowt_Vi,
 const std::vector<T_INSTANCES_CLUSTER_K>  &aivectort_numInstBi,
 const partition::Partition<T_CLUSTERIDX>  &aipartition_clustersBkinCi,
 const dist::Dist<T_METRIC,T_FEATURE>      &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::Dintra<>  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n(input  mat::MatrixRow<T_FEATURE>&: aimatrixrowt_S[" 
	      << &aimatrixrowt_S << ']' << " size: " << aimatrixrowt_S.getNumRows()  << '\n' 
	      << "input  mat::MatrixRow<T_FEATURE>&: aimatrixrowt_Vi[" 
	      << &aimatrixrowt_Vi << ']' << " size: " << aimatrixrowt_Vi.getNumRows()  << "\n"
              << "input  std::vector<T_INSTANCES_CLUSTER_K>&: aivectort_numInstBi[" 
	      << &aivectort_numInstBi << ']' << " size: " << aivectort_numInstBi.size()  << "\n"
              << "input  partition::Partition<T_CLUSTERIDX>&: aipartition_clustersBkinCi[" 
	      << &aipartition_clustersBkinCi << "]\n"
	      << "input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lort_Dintra = 0.0;
  T_CLUSTERIDX   lcidx_numClusterK = (T_CLUSTERIDX) aimatrixrowt_S.getNumRows();

  for ( uintidx luintidx_Vi = 0; luintidx_Vi < aimatrixrowt_Vi.getNumRows(); luintidx_Vi++) {
    
    T_CLUSTERIDX lcidx_ViinCj =
      aipartition_clustersBkinCi.getClusterIdx(luintidx_Vi);

    if ( 0 <= lcidx_ViinCj  && lcidx_ViinCj < lcidx_numClusterK ) {
      
      lort_Dintra += 
	(T_METRIC)
	aifunc2p_dist
	(aimatrixrowt_S.getRow(lcidx_ViinCj),
	 aimatrixrowt_Vi.getRow(luintidx_Vi),
	 data::Instance<T_FEATURE>::getNumDimensions()
	 )
	* (T_METRIC) aivectort_numInstBi.at(luintidx_Vi);
      
    }
    else {
      std::cout
	<< "um::Dintra: warning out of range number of clusters [0,"
	<< lcidx_numClusterK -1
	<< "] the value is " << lcidx_ViinCj
	<< std::endl;  
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::Dintra<>  OUT"
	      << '(' << geiinparam_verbose << ')' << " lort_Dintra = " << lort_Dintra 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lort_Dintra;

}


/*! \fn T_METRIC Dinter(const mat::MatrixRow<T_FEATURE> &aimatrixrowt_S, const mat::MatrixRow<T_FEATURE> &aimatrixrowt_Vi, const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstBi, const partition::Partition<T_CLUSTERIDX> &aipartition_clustersBkinCi, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist) 
  \brief Dinter \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
  \param aimatrixrowt_S a matrix with centroid S mat::MatrixRow<T_FEATURE>
  \param aimatrixrowt_Vi a matrix with centroid Vi 
  \param aivectort_numInstBi a vector with numero number of instances in Bi
  \param aipartition_clustersBkinCi a partition of B in the groups Ci
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
Dinter
(const mat::MatrixRow<T_FEATURE>           &aimatrixrowt_S,
 const mat::MatrixRow<T_FEATURE>           &aimatrixrowt_Vi,
 const std::vector<T_INSTANCES_CLUSTER_K>  &aivectort_numInstBi,
 const partition::Partition<T_CLUSTERIDX>  &aipartition_clustersBkinCi,
 const dist::Dist<T_METRIC,T_FEATURE>      &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::Dinter<>  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n(input  mat::MatrixRow<T_FEATURE>&: aimatrixrowt_S[" 
	      << &aimatrixrowt_S << "]\n"
	      << "input  mat::MatrixRow<T_FEATURE>&: aimatrixrowt_Vi[" 
	      << &aimatrixrowt_Vi << "]\n"
              << "input  std::vector<T_INSTANCES_CLUSTER_K>&: aivectort_numInstBi[" 
	      << &aivectort_numInstBi << "]\n"
              << "input  partition::Partition<T_CLUSTERIDX>&: aipartition_clustersBkinCi[" 
	      << &aipartition_clustersBkinCi << "]\n"
	      << "input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lort_Dinter = 0.0;
  T_CLUSTERIDX lcidx_numClusterK = (T_CLUSTERIDX) aimatrixrowt_S.getNumRows();

  if ( aimatrixrowt_S.getNumRows() > 1 ) {

    for ( uintidx luintidx_Vi = 0; luintidx_Vi < aimatrixrowt_Vi.getNumRows(); luintidx_Vi++) {
 
      T_METRIC  lt_distMinVkSj = std::numeric_limits<T_METRIC>::max();
 
      T_CLUSTERIDX lcidx_ViinCj =
	aipartition_clustersBkinCi.getClusterIdx(luintidx_Vi);
      if ( 0 <= lcidx_ViinCj  && lcidx_ViinCj < lcidx_numClusterK  ) {
	for ( T_CLUSTERIDX  lcidx_Cj = 0; lcidx_Cj < lcidx_numClusterK; lcidx_Cj++) { 
	  if ( lcidx_Cj != lcidx_ViinCj ) {
	    T_METRIC lt_distVkSj = 
	      aifunc2p_dist
	      (aimatrixrowt_S.getRow(lcidx_Cj),
	       aimatrixrowt_Vi.getRow(luintidx_Vi),
	      
	       data::Instance<T_FEATURE>::getNumDimensions()
	       ); 
	    if ( lt_distVkSj < lt_distMinVkSj )
	      lt_distMinVkSj = lt_distVkSj;
	  }
	} //for cluster
	lort_Dinter += (T_METRIC) lt_distMinVkSj * (T_METRIC) aivectort_numInstBi.at(luintidx_Vi); 
      }
    }    
  } //END IF 

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::Dinter<>  OUT"
	      << '(' << geiinparam_verbose << ')' << " lort_Dinter = " << lort_Dinter 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lort_Dinter;

}


/*! \fn T_METRIC e1(const T_FEATURE *aiarrayt_meanInstances, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist) 
  \brief e1 \cite Maulik:Bandyopadhyay:GAclustering:IndexI:2002 
  \details
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aiarrayt_meanInstances an array with the average of the attributes of the instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_METRIC  
	   >
T_METRIC
e1
(const T_FEATURE                      *aiarrayt_meanInstances,
 INPUT_ITERATOR                       aiiterator_instfirst,
 const INPUT_ITERATOR                 aiiterator_instlast,
 const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::e1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "input dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lrt_e1 = T_METRIC(0);

  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    data::Instance<T_FEATURE>* linst_inter =
      (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    lrt_e1 += 
      aifunc2p_dist
      (aiarrayt_meanInstances,
       linst_inter->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       ); 
  }
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')' << " lrt_e1 = " << lrt_e1
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lrt_e1;
    
}

  
/*! \fn T_METRIC e1 (const T_FEATURE *aiarrayt_meanInstances, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist, const FUNCINSTFREQUENCY func_instfrequency)
  \brief e1 \cite Maulik:Bandyopadhyay:GAclustering:IndexI:2002  for frequent instances 
  \details
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aiarrayt_meanInstances an array with the average of the attributes of the instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param func_instfrequency a function that returns the frequency of the instance
    
  \code{.cpp}

  DATATYPE_REAL lmetric_e1;

  lmetric_e1 =
  um::e1
  (larray_meanFeactures,
  lvectorptinst_instances.begin(),
  lvectorptinst_instances.end(),
  *pfunc2p_dist, 
  [](const data::Instance<DATATYPE_INST_CENT>* aiinst_iter ) -> DATATYPE_REAL
  {
  data::InstanceFreq
  <DATATYPE_INST_CENT,
  DATATYPE_INSTANCE_FREQUENCY
  >
  *lptinstfreq_iter =
  (data::InstanceFreq
  <DATATYPE_INST_CENT,
  DATATYPE_INSTANCE_FREQUENCY
  >*)
  aiinst_iter;

  return  (DATATYPE_REAL) lptinstfreq_iter->getFrequency(); 
  }
  );

  \endcode
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename FUNCINSTFREQUENCY, //RETURN INSTANCE_FREQUENCY AS TYPEDATA T_METRIC
	   typename T_METRIC  
	   >
T_METRIC
e1
(const T_FEATURE                       *aiarrayt_meanInstances,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist,
 const FUNCINSTFREQUENCY               func_instfrequency
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::e1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(input const const *aiarrayt_meanInstances["
	      << aiarrayt_meanInstances << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lrt_e1 = T_METRIC(0);

  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    lrt_e1 += 
      aifunc2p_dist
      (aiarrayt_meanInstances,
       linst_inter->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       )
      * func_instfrequency(linst_inter);
  }
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')' << " lrt_e1 = " << lrt_e1
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lrt_e1;
    
}

/*! \fn std::tuple<std::vector<T_METRIC>,std::vector<uintidx> > sumDistInstCentInK(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist)
  \brief Computes the overall between-cluster variance \f$SS_B\f$ and the number of instances in each cluster 
  \details 
  \f[ 
  SS_W = \sum_{j=1}^k \sum_{x_i \in C_j} \| x_i -\mu_j\|^2 
  \f]
  \param aimatrixt_centroids a mat::MatrixBase with centroids
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
std::tuple
<std::vector<T_METRIC>,
 std::vector<uintidx>
 >
sumDistInstCentInK
(const mat::MatrixBase<T_FEATURE>      &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::sumDistInstCentInK";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<T_METRIC> lovectort_sumDistXZi
    (aimatrixt_centroids.getNumRows(),T_METRIC(0));

  std::vector<uintidx> 
    lovectort_numInstClusterK
    (aimatrixt_centroids.getNumRows(),
     uintidx(0)
     );

  if ( aimatrixt_centroids.getNumRows() >= 1 ) {
    const T_CLUSTERIDX  lcidx_numClusterK = 
      T_CLUSTERIDX(aimatrixt_centroids.getNumRows());

    for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
      T_CLUSTERIDX lcidx_xinK =
	aipartition_clusters.next();
      
      if ( 0 <= lcidx_xinK  && lcidx_xinK <  lcidx_numClusterK  ) {

	data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
	
	lovectort_sumDistXZi.at(lcidx_xinK) += 
	  aifunc2p_dist
	  (aimatrixt_centroids.getRow(lcidx_xinK),
	   linst_inter->getFeatures(),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	++lovectort_numInstClusterK.at(lcidx_xinK);

      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelRadium;
    lostrstream_labelRadium << "<SUMDISTXZ: " << lpc_labelFunc
			    << "lovectort_sumDistXZi[" << &lovectort_sumDistXZi << ']';
    inout::containerprint
      (lovectort_sumDistXZi.begin(),
       lovectort_sumDistXZi.end(),
       std::cout,
       lostrstream_labelRadium.str().c_str(),
       ','
       );

    std::ostringstream lostrstream_labelNumInst;
    lostrstream_labelNumInst << "\n<NUMINSTCLUSTERK: " << lpc_labelFunc
			     << "lovectort_numInstClusterK[" << &lovectort_numInstClusterK << ']';
    inout::containerprint
      (lovectort_numInstClusterK.begin(),
       lovectort_numInstClusterK.end(),
       std::cout,
       lostrstream_labelNumInst.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_tuple
    (lovectort_sumDistXZi,
     lovectort_numInstClusterK
     );
    
}

/*! \fn std::vector<T_METRIC> sumDistInstCentInK(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast,partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist, const FUNCINSTFREQUENCY func_instfrequency)
  \brief Computes the overall between-cluster variance \f$SS_B\f$ for instances for frequent instances
  \details 
  \param aimatrixt_centroids a mat::MatrixBase with centroids
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param func_instfrequency a function that returns the frequency of an instance
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename FUNCINSTFREQUENCY, //RETURN INSTANCE_FREQUENCY AS TYPEDATA T_METRIC
	   typename T_METRIC  
	   >
std::vector<T_METRIC>
sumDistInstCentInK
(const mat::MatrixBase<T_FEATURE>           &aimatrixt_centroids,
 INPUT_ITERATOR                             aiiterator_instfirst,
 const INPUT_ITERATOR                       aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>  &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>       &aifunc2p_dist,
 const FUNCINSTFREQUENCY                    func_instfrequency
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::sumDistInstCentInK";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<T_METRIC> lovectort_sumDistXZi
    (aimatrixt_centroids.getNumRows(),T_METRIC(0));

  if ( aimatrixt_centroids.getNumRows() >= 1 ) {
    const T_CLUSTERIDX  lcidx_numClusterK = 
      T_CLUSTERIDX(aimatrixt_centroids.getNumRows());

    for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
      T_CLUSTERIDX lcidx_xinK = aipartition_clusters.next();
      
      if ( 0 <= lcidx_xinK  && lcidx_xinK <  lcidx_numClusterK  ) {

	data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
		
	lovectort_sumDistXZi[lcidx_xinK] += 
	  aifunc2p_dist
	  (aimatrixt_centroids.getRow(lcidx_xinK),
	   linst_inter->getFeatures(),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   )
	  * func_instfrequency(linst_inter);
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelRadium;
    lostrstream_labelRadium << "<SUMDISTXZ: " << lpc_labelFunc
			    << "lovectort_sumDistXZi[" << &lovectort_sumDistXZi << ']';
    inout::containerprint
      (lovectort_sumDistXZi.begin(),
       lovectort_sumDistXZi.end(),
       std::cout,
       lostrstream_labelRadium.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectort_sumDistXZi;
    
}

  
/*! \fn T_METRIC ssb(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, const T_FEATURE *aiarrayt_meanInstances, const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstancesClusterK, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief The overall between-cluster variance ssb
  \details
  \f[ 
  SS_B = \sum_{j=i}^k |C_j| \| \mu_j - M \|^2 
  \f]
  Where \f$\mu_j\f$ is the centroid of cluster \f$j\f$, \f$M\f$ is the overall mean of the instances

  \param aimatrixt_centroids a mat::MatrixRow with centroids
  \param aiarrayt_meanInstances an array with the mean of the attributes
  \param aivectort_numInstancesClusterK a vector with the number of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
           typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC  
	   >
T_METRIC
ssb
(const mat::MatrixRow<T_FEATURE>          &aimatrixt_centroids,
 const T_FEATURE                          *aiarrayt_meanInstances,
 const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstancesClusterK,
 const dist::Dist<T_METRIC,T_FEATURE>     &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::ssb";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ')'
      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
      << &aimatrixt_centroids << "]\n"
      << "\t input  T_FEATURE  *aiarrayt_meanInstances[" 
      << aiarrayt_meanInstances << "]\n"
      << "\t std::vector aivectort_numInstancesCluster[" 
      << aivectort_numInstancesClusterK << "]\n"
      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
      << &aifunc2p_dist << ']'
      << "\n\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES
  
  T_METRIC lt_SSb = T_METRIC(0.0);
  
  for ( uintidx luintidx_Ci = 0; luintidx_Ci < aimatrixt_centroids.getNumRows(); luintidx_Ci++) {
    lt_SSb += (T_METRIC) aivectort_numInstancesClusterK[luintidx_Ci] *
      aifunc2p_dist
      (aimatrixt_centroids.getRow(luintidx_Ci),
       aiarrayt_meanInstances,
       aimatrixt_centroids.getNumColumns()
       );    
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " lt_SSb = " << lt_SSb
      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lt_SSb;
}


template < typename INPUT_ITERATOR,
           typename T_FEATURE,
           typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC  
	   >
T_METRIC
ssbreeval
(const mat::MatrixRow<T_FEATURE>          &aimatrixt_centroids,
 INPUT_ITERATOR                           aiiterator_instfirst,
 const INPUT_ITERATOR                     aiiterator_instlast,
 const std::vector<T_INSTANCES_CLUSTER_K> &aivectort_numInstancesClusterK,
 const dist::Dist<T_METRIC,T_FEATURE>     &aifunc2p_dist
 )
{
  T_METRIC lt_SSb = T_METRIC(0.0);
  
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    *larray_sumFeatureTmp =
    new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    [data::Instance<T_FEATURE>::getNumDimensions()];

  stats::sumFeactures
    (larray_sumFeatureTmp,
     aiiterator_instfirst,
     aiiterator_instlast,
     T_FEATURE(0)
     );
  
  stats::meanVector
    (larray_centroid1,
     lui_numInstances,
     larray_sumFeatureTmp
     );

  for ( uintidx luintidx_Ci = 0; luintidx_Ci < aimatrixt_centroids.getNumRows(); luintidx_Ci++) {
    lt_SSb += (T_METRIC) aivectort_numInstancesClusterK[luintidx_Ci] *
      aifunc2p_dist
      (aimatrixt_centroids.getRow(luintidx_Ci),
       larray_centroid1,
       aimatrixt_centroids.getNumColumns()
       );    
  }

  delete [] larray_centroid1;
  delete [] larray_sumFeatureTmp;

  return lt_SSb;
}


/*! \fn T_METRIC T_METRIC VRC(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist)
  \brief  Variance ratio criterion (VRC) 

  \details VRC is sometimes called the Calinski-Harabasz criterion
  \f[ 
  VRC_k=\frac{SS_B}{SS_W} \times \frac{(n−k)}{(k−1)},   
  \f]
  Where \f$SS_B\f$ is the overall between-cluster variance, \f$SS_W\f$ is the overall  within-cluster variance, \f$k\f$ is the number of clusters, and $n$ is the number of instances. The larger the \f$VRC_k\f$ ratio, the better the data partition.
  \cite Calinski:Harabasz:Metricclustering:VRC:1974

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_squaredDist an object of type dist::Dist to calculate distances
  \param aib_withNullK a bool it is possible that for testing you have null clusters and VRC should be calculated omitting cluster nulls

  \code{.cpp}
    
  lT_VRC =
  um::VRC
  (lmatrixrowt_centroids,
  larray_centroid1,
  aivectorptinst_instances.begin(),
  aivectorptinst_instances.end(),
  lipartitionDisjSets_clusters,
  aifunc2p_squaredDist
  );

  \endcode
*/   
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
VRC
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_squaredDist,
 const bool                            aib_withNullK = false
 )
{  
  const T_CLUSTERIDX lcidx_numClusterK = T_CLUSTERIDX(aimatrixt_centroids.getNumRows());
  const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));

  static T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  static utils::RunOnce runOnce ([&]() {
      decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	*larray_sumFeatureTmp =
	new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	[data::Instance<T_FEATURE>::getNumDimensions()];

      stats::sumFeactures
	(larray_sumFeatureTmp,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 T_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_centroid1,
	 lui_numInstances,
	 larray_sumFeatureTmp
	 );

      delete [] larray_sumFeatureTmp;
    }
    );

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::VRC";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ')'
      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
      << &aimatrixt_centroids << "]\n"
      << "\t input  partition::Partition<>&: aipartition_clusters[" 
      << &aipartition_clusters << "]\n"
      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist[" 
      << &aifunc2p_squaredDist << ']'
      << "\n\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lometric_VRC = measuare_undefVRC(T_METRIC);
  T_METRIC  lmetrict_SSb = 0.0;
  T_METRIC  lmetrict_SSw = 0.0;
  
  if ( lcidx_numClusterK >= 2 ) {

    std::vector<T_METRIC> lvectorrt_sumDistInstCentInK;
    std::vector<uintidx>  lvectorui_numInstClusterK;
  
    std::tie
      (lvectorrt_sumDistInstCentInK,
       lvectorui_numInstClusterK) =
      sumDistInstCentInK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_squaredDist
       );

    const T_CLUSTERIDX lcidx_numClusterKNull =
      T_CLUSTERIDX
      (count(lvectorui_numInstClusterK.begin(),
       lvectorui_numInstClusterK.end(),
       0));

    if  ( (lcidx_numClusterK - lcidx_numClusterKNull) >= 2 ) {

      if ( lcidx_numClusterKNull == 0 || aib_withNullK  ) {
	
	lmetrict_SSw = 
	  interfacesse::sum
	  (lvectorrt_sumDistInstCentInK.data(),
	   (uintidx) lvectorrt_sumDistInstCentInK.size()
	   );

	if  ( lmetrict_SSw > 0.0 ) {
      
	  lmetrict_SSb =
	    ssb
	    (aimatrixt_centroids,
	     larray_centroid1,
	     lvectorui_numInstClusterK,
	     aifunc2p_squaredDist
	     );

	  lometric_VRC =
	    (lmetrict_SSb  *
	     ( (T_METRIC) lui_numInstances -
	       (T_METRIC) (lcidx_numClusterK - lcidx_numClusterKNull) ))
	    / (lmetrict_SSw * (T_METRIC) ( lcidx_numClusterK - lcidx_numClusterKNull -1 ));

	}
      }
    }
  }
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " lmetrict_SSb = " << lmetrict_SSb
      << " lmetrict_SSw = " << lmetrict_SSw
      << " lometric_VRC = " << lometric_VRC
      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_VRC;
    
}


/*! \fn T_METRIC T_METRIC VRCreeval(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist)
  \brief  Variance ratio criterion (VRC), recalculate the average centroid on each call.

  \details VRC is sometimes called the Calinski-Harabasz criterion, 

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_squaredDist an object of type dist::Dist to calculate distances
  \param aib_withNullK a bool it is possible that for testing you have null clusters and VRC should be calculated omitting cluster nulls
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
VRCreeval
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_squaredDist,
 const bool                            aib_withNullK = false
 )
{

  const T_CLUSTERIDX lcidx_numClusterK = T_CLUSTERIDX(aimatrixt_centroids.getNumRows());
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::VRCreeval";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist[" 
	      << &aifunc2p_squaredDist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lometric_VRC = measuare_undefVRC(T_METRIC);
  T_METRIC  lmetrict_SSb = 0.0;
  T_METRIC  lmetrict_SSw = 0.0;

  if ( lcidx_numClusterK >= 2 ) {

    const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
    T_FEATURE *larray_centroid1 =
      new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

    decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
      *larray_sumFeatureTmp =
      new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
      [data::Instance<T_FEATURE>::getNumDimensions()];

    stats::sumFeactures
      (larray_sumFeatureTmp,
       aiiterator_instfirst,
       aiiterator_instlast,
       T_FEATURE(0)
       );
  
    stats::meanVector
      (larray_centroid1,
       lui_numInstances,
       larray_sumFeatureTmp
       );
 
    std::vector<T_METRIC> lvectorrt_sumDistInstCentInK;
    std::vector<uintidx>  lvectorui_numInstClusterK;
  
    std::tie
      (lvectorrt_sumDistInstCentInK,
       lvectorui_numInstClusterK) =
      sumDistInstCentInK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_squaredDist
       );

    const T_CLUSTERIDX lcidx_numClusterKNull =
      T_CLUSTERIDX
      (count(lvectorui_numInstClusterK.begin(),
       lvectorui_numInstClusterK.end(),
       0));
    
    if  ( (lcidx_numClusterK - lcidx_numClusterKNull) >= 2 ) {

      if ( lcidx_numClusterKNull == 0 || aib_withNullK  ) {
	
	lmetrict_SSw = 
	  interfacesse::sum
	  (lvectorrt_sumDistInstCentInK.data(),
	   (uintidx) lvectorrt_sumDistInstCentInK.size()
	   );

	if  ( lmetrict_SSw > 0.0 ) {
    
	  lmetrict_SSb =
	    ssb
	    (aimatrixt_centroids,
	     larray_centroid1,
	     lvectorui_numInstClusterK,
	     aifunc2p_squaredDist
	     );

	   lometric_VRC =
	    (lmetrict_SSb  *
	     ( (T_METRIC) lui_numInstances -
	       (T_METRIC) (lcidx_numClusterK - lcidx_numClusterKNull) ))
	    / (lmetrict_SSw * (T_METRIC) ( lcidx_numClusterK - lcidx_numClusterKNull -1 ));
	   
	}
      }
    }
  
    delete [] larray_centroid1;
    delete [] larray_sumFeatureTmp;
    
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " lmetrict_SSb = " << lmetrict_SSb
      << " lmetrict_SSw = " << lmetrict_SSw
      << " lometric_VRC = " << lometric_VRC
      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lometric_VRC;
    
}


/*! \fn T_METRIC T_METRIC WBIndex(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist)
  \brief  WB-Index \cite Zhao:Xu:Franti:ClusterMeasure:WBIndex:2009 \cite Zhao:Franti:ClusterMeasure:WBIndex:2014
  \details WB-Index 
  \f[ 
  WBIndex_k= k \times \frac{SS_B}{SS_W}
  \f]
  Where \f$SS_B\f$ is the overall between-cluster variance, \f$SS_W\f$ is the overall  within-cluster variance, \f$k\f$ is the number of clusters, and $n$ is the number of instances. The minimum value of \f$WBIndex_k\f$ determine the best cluster number.
  \cite Zhao:Xu:Franti:ClusterMeasure:WBIndex:2009 \cite Zhao:Franti:ClusterMeasure:WBIndex:2014

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_squaredDist an object of type dist::Dist to calculate distances

  \code{.cpp}
    
  lT_WBIndex =
  um::WBIndex
  (lmatrixrowt_centroids,
  larray_centroid1,
  aivectorptinst_instances.begin(),
  aivectorptinst_instances.end(),
  lipartitionDisjSets_clusters,
  aifunc2p_squaredDist
  );

  \endcode
*/   
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
WBIndex
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_squaredDist
 )
{
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  static T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  static utils::RunOnce runOnce ([&]() {
      decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	*larray_sumFeatureTmp =
	new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	[data::Instance<T_FEATURE>::getNumDimensions()];

      stats::sumFeactures
	(larray_sumFeatureTmp,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 T_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_centroid1,
	 lui_numInstances,
	 larray_sumFeatureTmp
	 );

      delete [] larray_sumFeatureTmp;
    }
    );
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::WBIndex";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist[" 
	      << &aifunc2p_squaredDist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lometric_WBIndex = measuare_undefWBIndex(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 1 ) { 

    std::vector<T_METRIC> lvectorrt_sumDistInstCentInK;
    std::vector<uintidx>  lvectorui_numInstClusterK;
  
    std::tie
      (lvectorrt_sumDistInstCentInK,
       lvectorui_numInstClusterK) =
      sumDistInstCentInK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_squaredDist
       );

    T_METRIC lmetrict_SSb =
      ssb
      (aimatrixt_centroids,
       larray_centroid1,
       lvectorui_numInstClusterK,
       aifunc2p_squaredDist
       );

    if  ( lmetrict_SSb > 0.0 ) {

      T_METRIC lmetrict_SSw = 
	interfacesse::sum
	(lvectorrt_sumDistInstCentInK.data(),
	 (uintidx) lvectorrt_sumDistInstCentInK.size()
	 );

      lometric_WBIndex = (lmetrict_SSw / lmetrict_SSb) *
	(T_METRIC) ( aimatrixt_centroids.getNumRows() -1 );    
    }
    
  }
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_WBIndex = " << lometric_WBIndex
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_WBIndex;
    
}



template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
WBIndexreeval
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_squaredDist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::WBIndexreeval";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_squaredDist[" 
	      << &aifunc2p_squaredDist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lometric_WBIndex = measuare_undefWBIndex(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 1 ) {

    const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
    T_FEATURE *larray_centroid1 =
      new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

    decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
      *larray_sumFeatureTmp =
      new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
      [data::Instance<T_FEATURE>::getNumDimensions()];

    stats::sumFeactures
      (larray_sumFeatureTmp,
       aiiterator_instfirst,
       aiiterator_instlast,
       T_FEATURE(0)
       );
  
    stats::meanVector
      (larray_centroid1,
       lui_numInstances,
       larray_sumFeatureTmp
       );
                
    std::vector<T_METRIC> lvectorrt_sumDistInstCentInK;
    std::vector<uintidx>  lvectorui_numInstClusterK;
  
    std::tie
      (lvectorrt_sumDistInstCentInK,
       lvectorui_numInstClusterK) =
      sumDistInstCentInK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_squaredDist
       );

    T_METRIC lmetrict_SSb =
      ssb
      (aimatrixt_centroids,
       larray_centroid1,
       lvectorui_numInstClusterK,
       aifunc2p_squaredDist
       );

    if  ( lmetrict_SSb > 0.0 ) {

      T_METRIC lmetrict_SSw = 
	interfacesse::sum
	(lvectorrt_sumDistInstCentInK.data(),
	 (uintidx) lvectorrt_sumDistInstCentInK.size()
	 );

      lometric_WBIndex = (lmetrict_SSw / lmetrict_SSb) *
	(T_METRIC) ( aimatrixt_centroids.getNumRows() -1 );    
    }

    delete [] larray_centroid1;
    delete [] larray_sumFeatureTmp;

  }
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_WBIndex = " << lometric_WBIndex
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_WBIndex;
    
}


/*! \fn T_METRIC scoreFunction(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Score function (\f$SF\f$)
  \cite Sandro:Benny:IanF:ClusteringMeasure:2007
  \details  Score function is based on inter-cluster and intra-cluster distances. The score function is used for two purposes: i) to estimate the number of clusters and ii) to evaluate the quality of the clustering results \cite Sandro:Benny:IanF:ClusteringMeasure:2007. The higher the value of the Score function, the more suitable the number of clusters.
  \f[ 
  SF = 1 - \frac{1}{e^{e^{bcd-wcd}}} 
  \f]
  Where
  \f[ 
  bcd = \frac{\sum_{j=1}^k \| \mu_k -M \| \cdot |C_j|}{n \cdot k}, 
  \f]

  \f[ 
  wcd = \sum_{j=1}^k \left( \frac{1}{ |C_j|} \sum_{x_i \in C_j} \| x_i - \mu_j \| \right)
  \f]

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition::Partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
scoreFunction
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  static T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  static utils::RunOnce runOnce ([&]() {
      decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	*larray_sumFeatureTmp =
	new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	[data::Instance<T_FEATURE>::getNumDimensions()];

      stats::sumFeactures
	(larray_sumFeatureTmp,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 T_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_centroid1,
	 lui_numInstances,
	 larray_sumFeatureTmp
	 );

      delete [] larray_sumFeatureTmp;
    }
    );
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::scoreFunction";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lometric_scoreFunction = measuare_undefScoreFunction(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 0 ) {
    
    std::vector<T_METRIC> lvectorrt_sumDistInstCentInK;
    std::vector<uintidx>  lvectorui_numInstClusterK;
  
    std::tie
      (lvectorrt_sumDistInstCentInK,
       lvectorui_numInstClusterK) =
      sumDistInstCentInK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_dist
       );
   
    T_METRIC lmetrict_wcd =
      std::inner_product
      (lvectorrt_sumDistInstCentInK.begin(),
       lvectorrt_sumDistInstCentInK.end(),
       lvectorui_numInstClusterK.begin(),
       T_METRIC(0.0),
       [](T_METRIC airt_partialsum, T_METRIC airt_prod) 
       {return airt_partialsum + airt_prod; },
       [](T_METRIC airt_sumDistInstCentInK,
	  uintidx  aiit_numInstClusterK )
       {return (aiit_numInstClusterK!=0)?
	airt_sumDistInstCentInK / T_METRIC(aiit_numInstClusterK)
	:T_METRIC(0.0);
       }
       );

    /* SS_B  == bcd */
    T_METRIC lmetrict_bcd =
      ssb
      (aimatrixt_centroids,
       larray_centroid1,
       lvectorui_numInstClusterK,
       aifunc2p_dist
       )
      / ( T_METRIC(lui_numInstances)  * T_METRIC(aimatrixt_centroids.getNumRows()) );
 

    lometric_scoreFunction =
      1.0 - 1.0 /(std::exp(std::exp(lmetrict_bcd-lmetrict_wcd)));
    
  }
    
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_scoreFunction = " << lometric_scoreFunction
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_scoreFunction;
    
} /*END scoreFunction
   */
  

template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
scoreFunctionreeval
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::scoreFunctionreeval";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];

  decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    *larray_sumFeatureTmp =
    new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    [data::Instance<T_FEATURE>::getNumDimensions()];

  stats::sumFeactures
    (larray_sumFeatureTmp,
     aiiterator_instfirst,
     aiiterator_instlast,
     T_FEATURE(0)
     );
  
  stats::meanVector
    (larray_centroid1,
     lui_numInstances,
     larray_sumFeatureTmp
     );

  T_METRIC  lometric_scoreFunction = measuare_undefScoreFunction(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 0 ) {
    
    std::vector<T_METRIC> lvectorrt_sumDistInstCentInK;
    std::vector<uintidx>  lvectorui_numInstClusterK;
  
    std::tie
      (lvectorrt_sumDistInstCentInK,
       lvectorui_numInstClusterK) =
      sumDistInstCentInK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_dist
       );
   
    T_METRIC lmetrict_wcd =
      std::inner_product
      (lvectorrt_sumDistInstCentInK.begin(),
       lvectorrt_sumDistInstCentInK.end(),
       lvectorui_numInstClusterK.begin(),
       T_METRIC(0.0),
       [](T_METRIC airt_partialsum, T_METRIC airt_prod) 
       {return airt_partialsum + airt_prod; },
       [](T_METRIC airt_sumDistInstCentInK,
	  uintidx  aiit_numInstClusterK )
       {return (aiit_numInstClusterK!=0)?
	airt_sumDistInstCentInK / T_METRIC(aiit_numInstClusterK)
	:T_METRIC(0.0);
       }
       );

    /* SS_B  == bcd */
    T_METRIC lmetrict_bcd =
      ssb
      (aimatrixt_centroids,
       larray_centroid1,
       lvectorui_numInstClusterK,
       aifunc2p_dist
       )
      / ( T_METRIC(lui_numInstances)  * T_METRIC(aimatrixt_centroids.getNumRows()) );
 

    lometric_scoreFunction =
      1.0 - 1.0 /(std::exp(std::exp(lmetrict_bcd-lmetrict_wcd)));
    
  }
    
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_scoreFunction = " << lometric_scoreFunction
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  delete [] larray_centroid1;
  delete [] larray_sumFeatureTmp;

  return lometric_scoreFunction;
    
} /*END scoreFunction
   */

/*! \fn T_METRIC distMinIntraClusterKj(const T_CLUSTERIDX aicidx_Clusterj, const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, INPUT_ITERATOR  aiiterator_instfirst, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief  Distance min intra cluster
  \details 
  \param aicidx_Clusterj an integer index of the cluster
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_CLUSTERIDX,
           typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_METRIC
	   >
T_METRIC
distMinIntraClusterKj
(const T_CLUSTERIDX                       aicidx_Clusterj,
 const ds::PartitionLinked<T_CLUSTERIDX>  &aipartlink_memberShip,
 INPUT_ITERATOR                           aiiterator_instfirst,
 const dist::Dist<T_METRIC,T_FEATURE>     &aifunc2p_dist
 )
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::distMinIntraClusterKj";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t( input const T_CLUSTERIDX  aicidx_Clusterj[" << aicidx_Clusterj << "]\n"
      << "\t  input PartitionLinked&: aipartlink_memberShip[" << &aipartlink_memberShip << "]\n"
      << "\t  input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
      << "\t  input dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" << &aifunc2p_dist << "]\n"
      << "\t  lcidx_numClusterK = " << lcidx_numClusterK << '\n'
      << "\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_i(&aipartlink_memberShip);
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_ip(&aipartlink_memberShip);
  
  T_METRIC  lot_distMinIntraClusterKj = std::numeric_limits<T_METRIC>::max();

  for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {
    
    if ( aicidx_Clusterj != lcidx_Ck  ) {

      for ( literpart_i.begin(aicidx_Clusterj); literpart_i.end(); literpart_i.next() ) {
	
	for ( literpart_ip.begin(lcidx_Ck); literpart_ip.end(); literpart_ip.next() ) {
	  
	  data::Instance<T_FEATURE>* linst_i = *std::next(aiiterator_instfirst,literpart_i.getValue());
	  data::Instance<T_FEATURE>* linst_ip = *std::next(aiiterator_instfirst,literpart_ip.getValue());
	  T_METRIC lt_distk =
	    aifunc2p_dist
	    (linst_i->getFeatures(),
	     linst_ip->getFeatures(),
	     data::Instance<T_FEATURE>::getNumDimensions()
	     );
	  if ( lt_distk < lot_distMinIntraClusterKj )
	    lot_distMinIntraClusterKj = lt_distk;
	}
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " T_METRIC lot_distMinIntraClusterKj = " << lot_distMinIntraClusterKj
      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lot_distMinIntraClusterKj;
}


/*! \fn T_METRIC distMinIntraClusterKj(const T_CLUSTERIDX aicidx_Clusterj, const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity) 
  \brief  Distance min intra cluster
  \details 
  \param aicidx_Clusterj an integer index of the cluster
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aimatrixtriagt_dissimilarity a matrix of distances
*/
template < typename T_CLUSTERIDX,
	   typename T_METRIC
	   >
T_METRIC
distMinIntraClusterKj
(const T_CLUSTERIDX                      aicidx_Clusterj,
 const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 const mat::MatrixTriang<T_METRIC>       &aimatrixtriagt_dissimilarity
) 
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::distMinIntraClusterKj";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc 
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t( input const T_CLUSTERIDX  aicidx_Clusterj[" << aicidx_Clusterj << "]\n"
      << "\t  input PartitionLinked&: aipartlink_memberShip[" << &aipartlink_memberShip << "]\n"
      << "\t input  mat::MatrixTriang<T_METRIC>: aimatrixtriagt_dissimilarity[" 
      << &aimatrixtriagt_dissimilarity << "]\n"
      << "\t  lcidx_numClusterK = " << lcidx_numClusterK << '\n'
      << "\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES  
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_i(&aipartlink_memberShip);
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_ip(&aipartlink_memberShip);
  
  T_METRIC  lot_distMinIntraClusterKj = std::numeric_limits<T_METRIC>::max();

  for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {
    
    if ( aicidx_Clusterj != lcidx_Ck  ) {

      for ( literpart_i.begin(aicidx_Clusterj); literpart_i.end(); literpart_i.next() ) {
	
	for ( literpart_ip.begin(lcidx_Ck); literpart_ip.end(); literpart_ip.next() ) {
	  
	  T_METRIC lt_distk =
	    aimatrixtriagt_dissimilarity(literpart_i.getValue(),literpart_ip.getValue());

	  if ( lt_distk < lot_distMinIntraClusterKj )
	    lot_distMinIntraClusterKj = lt_distk;
	}
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ": OUT(" << geiinparam_verbose << ')'
      << " T_METRIC lot_distMinIntraClusterKj = " << lot_distMinIntraClusterKj
      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lot_distMinIntraClusterKj;
}

/*! \fn std::vector<T_METRIC> diameterClusterK(INPUT_ITERATOR aiiterator_instfirst, ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Calculate the cluster diameters for a partition
  \details
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
	   typename T_CLUSTERIDX,
	   typename T_FEATURE,
	   typename T_METRIC 
	   >
std::vector<T_METRIC>
diameterClusterK
(INPUT_ITERATOR                          aiiterator_instfirst,
 const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::diameterClusterK";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "(\n input  PartitionLinked&: aipartlink_memberShip[" 
	      << &aipartlink_memberShip << "]\n"
	      << " input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "] K = " << aipartlink_memberShip.getNumPartitions()
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  const  T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();
    
  std::vector<T_METRIC>        lvectort_diameterClusterK(lcidx_numClusterK,T_METRIC(0));
    
  for ( T_CLUSTERIDX lcidx_Ci = 0; lcidx_Ci < lcidx_numClusterK; lcidx_Ci++) {
      
    lvectort_diameterClusterK[lcidx_Ci] = 
      diameterClusterKj
      (lcidx_Ci,
       aipartlink_memberShip,
       aiiterator_instfirst,
       aifunc2p_dist
       );
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelDiameter;
    lostrstream_labelDiameter << "<DIAMETER: " << lpc_labelFunc
			      << "lvectort_diameterClusterK[" << &lvectort_diameterClusterK << ']';
    inout::containerprint
      (lvectort_diameterClusterK.begin(),
       lvectort_diameterClusterK.end(),
       std::cout,
       lostrstream_labelDiameter.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
     
  return lvectort_diameterClusterK;
  
} //END diameterClusterK



/*! \fn T_METRIC overlap (const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist)

  \brief Overlap
  \details To measure the overlap, we calculated the distance from every point to its centroid (d1) and to its nearest point in another cluster (d2). If this nearest point is closer than its own centroid (d1 > d2), the point is evidence of overlap \cite Franti:Sieranoja:ClusterMeasure:Benchmark:2018. 

  \f[
  Overlap = \frac{1}{n}\sum ov(d_1,d_2)
  \f]

  where
  \f[
  ov(d_1,d_2) = 
  \left\{
        \begin{array}{ll}
               1,  & \mbox{if } d_1 > d_2  \\
               0,  & \mbox{otherwise}
        \end{array}
  \right
  \f]

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/ 
template < typename INPUT_ITERATOR,
	   typename T_METRIC, 
	   typename T_FEATURE,
	   typename T_CLUSTERIDX /*-1,0,..,K*/
	   >
T_METRIC
overlap
(const mat::MatrixRow<T_FEATURE>         &aimatrixt_centroids,
 const INPUT_ITERATOR                    aiiterator_instfirst,
 const INPUT_ITERATOR                    aiiterator_instlast,
 const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist
 )
{  

  T_CLUSTERIDX lcidx_numClusterK = 
    (T_CLUSTERIDX) aimatrixt_centroids.getNumRows();

  /*CALCULATE RADIUS CLUSTERS
   */
  std::vector<T_METRIC>
    lvector_radiusClusterK(aimatrixt_centroids.getNumRows());
  for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {
    lvector_radiusClusterK.at(lcidx_Ck) = 
      radiusClusterKj
      (lcidx_Ck,
       aimatrixt_centroids.getRow(lcidx_Ck),
       aipartlink_memberShip,
       aiiterator_instfirst,
       aifunc2p_dist
       );
  }

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  const  char* lpc_labelFunc = "um::overlap";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout 
      << lpc_labelFunc
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\n\t input  mat::MatrixRow<T_FEATURE>&: aimatrixt_centroids[" 
      << &aimatrixt_centroids << ']'
      << "\n\t input aiiterator_instfirst[" << *aiiterator_instfirst << ']'
      << "\n\t input aiiterator_instlast[" << *aiiterator_instlast << ']'
      << "\n\t input ds::PartitionLinked<T_CLUSTERIDX>&: aipartlink_memberShip[" 
      << &aipartlink_memberShip << ']'
      << "\n\t input  aifunc2p_dist[" << &aifunc2p_dist << ']'
      << "\n\t lvector_radiusClusterK[" << lvector_radiusClusterK << ']'
      << "\n\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES
  
  uintidx  lui_overlap = 0;
  T_METRIC lot_overlap  = measuare_undefOverlap(T_METRIC);

  if  ( lcidx_numClusterK >= 1 ) {
    
    ds::IteratorPartitionLinked <T_CLUSTERIDX>
      literpart_j(&aipartlink_memberShip);
      
    for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {
      
      for ( literpart_j.begin(lcidx_Ck); literpart_j.end(); literpart_j.next() ) {
    
	data::Instance<T_FEATURE>* linst_j = *std::next(aiiterator_instfirst,literpart_j.getValue());

	T_METRIC lt_distCinK =
	  aifunc2p_dist
	  (linst_j->getFeatures(),
	   aimatrixt_centroids.getRow(lcidx_Ck),
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );

	bool lbool_isOverlap = false;
	for ( T_CLUSTERIDX lcidx_Ckp = 0; lcidx_Ckp < lcidx_numClusterK; lcidx_Ckp++) {

	  if ( lcidx_Ck != lcidx_Ckp ) {

	    T_METRIC lt_distCOtherK =
	      aifunc2p_dist
	      (linst_j->getFeatures(),
	       aimatrixt_centroids.getRow(lcidx_Ckp),
	       data::Instance<T_FEATURE>::getNumDimensions()
	       );
	    if  ( std::abs(lt_distCOtherK - lvector_radiusClusterK.at(lcidx_Ckp)) <= lt_distCinK ) {

	      ds::IteratorPartitionLinked <T_CLUSTERIDX>
		literpart_jp(&aipartlink_memberShip);

	      for ( literpart_jp.begin(lcidx_Ckp); literpart_jp.end(); literpart_jp.next() ) {

		data::Instance<T_FEATURE>* linst_jp = 
		  *std::next(aiiterator_instfirst,literpart_jp.getValue());

		T_METRIC lt_distInstInst =
		  aifunc2p_dist
		  (linst_j->getFeatures(),
		   linst_jp->getFeatures(),
		   data::Instance<T_FEATURE>::getNumDimensions()
		   );
		if ( lt_distInstInst < lt_distCinK ) {
		  lbool_isOverlap = true;
		  break;
		}
	      }
	    } //if distance
	  } //if Other cluster
	  if ( lbool_isOverlap ) {
	    break;
	  }
	}
	if ( lbool_isOverlap ) {
	  ++lui_overlap;
	}
      } //Next instance of cluster
    }
    const uintidx  lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
    
    lot_overlap =(lui_numInstances != 0 )? 
      T_METRIC(lui_overlap) / T_METRIC(lui_numInstances)
      :measuare_undefOverlap(T_METRIC);;
  }  
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lot_overlap = "   << lot_overlap
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lot_overlap;

} /*END overlap*/



/*! \fn T_METRIC DunnIndex(mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity, ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const bool aib_withNullK = false) 
  \brief  The Dunn Index (DI) \cite Dunn:ClusterMeasure:CS:1974 \cite Zhang:Cao:KernelclusteringLabelKVar:2011 
  \details The Dunn Index determines the minimum ratio between inter-cluster distance and cluster diameter for a given partitioning. It captures the notion that, in a good clustering solution, data elements within one cluster should be much closer to each other than to elements within different clusters. It is defined as

  \f[
  DI(C) = {  \min  \atop j \in C} \left\{ { \min \it \atop j' \in C, j' \neq j } \left\{ { \delta(C_{j},C_{j'}) \over  { \max  \atop j'' \in C}  \{\Delta(C_{j''})\} }  \right\} \right\}  
  \f]

  Where

  \f[
  \delta(C_{j},C_{j'}) =   \min   \left\{  d(x_{i},x_{i''}) | x_i \in C_j, x_{i''} \in C_{j''}  \right\}
  \f]

  \f[
  \Delta(C_{j}) =   \max   \left\{  d(x_{i},x_{i''}) | x_i , x_{i''} \in C_{j}  \right\}
  \f]

  \param aimatrixtriagt_dissimilarity a matrix of distances
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aib_withNullK a bool it is possible that for testing you have null clusters and Dunn Index should be calculated omitting cluster nulls
*/
template < typename T_CLUSTERIDX,
	   typename T_METRIC
	   >
T_METRIC
DunnIndex
(mat::MatrixTriang<T_METRIC>       &aimatrixtriagt_dissimilarity,
 ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip,
 const bool                        aib_withNullK = false
 ) 
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();
  const T_CLUSTERIDX lcidx_numClusterKNull = aipartlink_memberShip.getNumPartitionsNull(); 
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::DunnIndex";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout
      << lpc_labelFunc
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "\t(input  mat::MatrixTriang<T_METRIC>: aimatrixtriagt_dissimilarity[" 
      << &aimatrixtriagt_dissimilarity << "]\n"
      << "\t input  PartitionLinked&: aipartlink_memberShip[" << &aipartlink_memberShip << "]\n"
      << "\t lcidx_numClusterK =  " << lcidx_numClusterK << '\n'
      << "\t lcidx_numClusterKNull = " << lcidx_numClusterKNull << '\n'
      << "\t)"
      << std::endl;
  }
#endif //__VERBOSE_YES


   T_METRIC lot_DunnIndex = measuare_undefDunnIndex(T_METRIC);
 
  if  ( (lcidx_numClusterK - lcidx_numClusterKNull) >= 2 ) {

    if ( lcidx_numClusterKNull == 0 || aib_withNullK  ) {
      
      lot_DunnIndex = std::numeric_limits<T_METRIC>::max();
      
      for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {

	if ( !aipartlink_memberShip.isPartitionsNull(lcidx_Ck) ) {  
	  T_METRIC lrt_delta =
	    distMinIntraClusterKj
	    (lcidx_Ck,
	     aipartlink_memberShip,
	     aimatrixtriagt_dissimilarity
	     );
      
	  T_METRIC lrt_Delta =
	    diameterClusterKj
	    (lcidx_Ck,
	     aipartlink_memberShip,
	     aimatrixtriagt_dissimilarity
	     );

	  T_METRIC lrt_partilDI  = (lrt_Delta != 0)?
	    lrt_delta / lrt_Delta:measuare_undefDunnIndex(T_METRIC);

	  if ( lrt_partilDI < lot_DunnIndex) 
	    lot_DunnIndex = lrt_partilDI;
	}
      }
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lot_DunnIndex = "   << lot_DunnIndex
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lot_DunnIndex;
  
}


/*! \fn T_METRIC DunnIndex(INPUT_ITERATOR aiiterator_instfirst, const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist,  const bool aib_withNullK = false)
  \brief  The Dunn Index (DI) \cite Dunn:ClusterMeasure:CS:1974 \cite Zhang:Cao:KernelclusteringLabelKVar:2011
  \details The Dunn Index determines the minimum ratio between inter-cluster distance and cluster diameter for a given partitioning. It captures the notion that, in a good clustering solution, data elements within one cluster should be much closer to each other than to elements within different clusters.  It is defined as

  \f[
  DI(C) = {  \min  \atop j \in C} \left\{ { \min \it \atop j' \in C, j' \neq j } \left\{ { \delta(C_{j},C_{j'}) \over  { \max  \atop j'' \in C}  \{\Delta(C_{j''})\} }  \right\} \right\}  
  \f]

  Where

  \f[
  \delta(C_{j},C_{j'}) =   \min   \left\{  d(x_{i},x_{i''}) | x_i \in C_j, x_{i''} \in C_{j''}  \right\}
  \f]

  \f[
  \Delta(C_{j}) =   \max   \left\{  d(x_{i},x_{i''}) | x_i , x_{i''} \in C_{j}  \right\}
  \f]
      
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param aib_withNullK a bool it is possible that for testing you have null clusters and Dunn Index should be calculated omitting cluster nulls
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC
	   >
T_METRIC
DunnIndex
(INPUT_ITERATOR                        aiiterator_instfirst,
 ds::PartitionLinked<T_CLUSTERIDX>     &aipartlink_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist,
 const bool                            aib_withNullK = false
 )
{ //BEGIN DunnIndex
  
  const T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();
  const T_CLUSTERIDX lcidx_numClusterKNull = aipartlink_memberShip.getNumPartitionsNull(); 
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::DunnIndex";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << " input &aipartlink_memberShip[" << &aipartlink_memberShip << "]\n";
    aipartlink_memberShip.print();
    std::cout << "\nlcidx_numClusterK =  " << lcidx_numClusterK << '\n'
	      << "lcidx_numClusterKNull = " << lcidx_numClusterKNull << '\n'
              << " input const &aifunc2p_dist[" << &aifunc2p_dist << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lot_DunnIndex = measuare_undefDunnIndex(T_METRIC);
 
  if  ( (lcidx_numClusterK - lcidx_numClusterKNull) >= 2 ) {

    if ( lcidx_numClusterKNull == 0 || aib_withNullK  ) {
      
      lot_DunnIndex = std::numeric_limits<T_METRIC>::max();
      
      for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {

	if ( !aipartlink_memberShip.isPartitionsNull(lcidx_Ck) ) {  
	  T_METRIC lrt_delta =
	    distMinIntraClusterKj
	    (lcidx_Ck,
	     aipartlink_memberShip,
	     aiiterator_instfirst,
	     aifunc2p_dist
	     );
      
	  T_METRIC lrt_Delta =
	    diameterClusterKj
	    (lcidx_Ck,
	     aipartlink_memberShip,
	     aiiterator_instfirst,
	     aifunc2p_dist
	     );

	  T_METRIC lrt_partilDI  = (lrt_Delta != 0)?
	    lrt_delta / lrt_Delta:measuare_undefDunnIndex(T_METRIC);

	  if ( lrt_partilDI < lot_DunnIndex) 
	    lot_DunnIndex = lrt_partilDI;
	}
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lot_DunnIndex = "   << lot_DunnIndex
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lot_DunnIndex;
  
} //END DunnIndex


/*! \fn T_METRIC simplifiedDunnIndex(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist, const bool aib_withNullK = false)
  \brief  The  Simplified Dunn Index (SDI)
  \details To simplify the calculations of the the Dunn Index, simplifiedDunnIndex uses the centroids to calculate the radius and decrease the computational complexity
      
  \param aimatrixt_centroids a matrix with the centroids of the clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param aib_withNullK a bool it is possible that for testing you have null clusters and Dunn Index should be calculated omitting cluster nulls
*/  
template < typename T_FEATURE,
	   typename INPUT_ITERATOR,
           typename T_CLUSTERIDX,
	   typename T_METRIC
	   >
T_METRIC 
simplifiedDunnIndex
(const mat::MatrixBase<T_FEATURE>         &aimatrixt_centroids, 
 INPUT_ITERATOR                           aiiterator_instfirst,
 ds::PartitionLinked<T_CLUSTERIDX>        &aipartlink_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>     &aifunc2p_dist,
 const bool                               aib_withNullK = false
 )
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();
  const T_CLUSTERIDX lcidx_numClusterKNull = aipartlink_memberShip.getNumPartitionsNull();
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::simplifiedDunnIndex";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << " input  aipartlink_memberShip[" 
	      << &aipartlink_memberShip << "]\n";
    aipartlink_memberShip.print();
    std::cout << "\nlcidx_numClusterK =  " << lcidx_numClusterK << '\n'
	      << "lcidx_numClusterKNull = " << lcidx_numClusterKNull << '\n'
	      << " input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lot_DunnIndex = measuare_undefDunnIndex(T_METRIC);
  
  if  ( (lcidx_numClusterK - lcidx_numClusterKNull) >= 2 ) {

    if ( lcidx_numClusterKNull == 0 || aib_withNullK  ) {
      
      lot_DunnIndex = std::numeric_limits<T_METRIC>::max();
      
      for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) {

	if ( !aipartlink_memberShip.isPartitionsNull(lcidx_Ck) ) {
	  
	  T_METRIC lrt_delta =
	    minDistCjCjp
	    (lcidx_Ck,
	     aimatrixt_centroids,
	     aifunc2p_dist
	     );
      
	  T_METRIC lrt_Delta =
	    radiusClusterKj
	    (lcidx_Ck,
	     aimatrixt_centroids.getRow(lcidx_Ck),
	     aipartlink_memberShip,
	     aiiterator_instfirst,
	     aifunc2p_dist
	     );

	  T_METRIC lrt_partilDI  =
	    (lrt_Delta != 0)? lrt_delta / (2.0 * lrt_Delta):measuare_undefDunnIndex(T_METRIC);

	  if ( lrt_partilDI < lot_DunnIndex) 
	    lot_DunnIndex = lrt_partilDI;
	}
      }
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lot_DunnIndex = "   << lot_DunnIndex
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lot_DunnIndex;

}

  
/*! \fn T_METRIC T_METRIC CSmeasure(INPUT_ITERATOR aiiterator_instfirst, const mat::MatrixRow<T_FEATURE>  &aimatrixt_centroids, const ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &aipartlinknuminst_memberShip, const dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist)
  \brief CS measure \cite Chou:Su:Lai:ClusterMeasure:CS:2004 \cite Das:Abraham:Konar:GAclusteringLabelKVar:ACDE:2008
  \details CS measure is a function of the ratio of the sum of within-cluster scatter to between-cluster separation. The smallest \f$CS(C)\f$ indicates a valid optimal partition.

  \f[
  CS(C) =
  {
  {1 \over k} \sum_{j=1}^{k}
  \left\{
  {1 \over |C_j| } \sum_{x_i \in C_j} { \max \atop x_{i'} \in C_j }
  \left\{ D(x_j,x_{j'}) \right\} \right\} 
  \over
  {1 \over k} \sum_{j=1}^{k} \left\{ { \min \atop j \in k, j \neq j'}
  \left\{ D(\mu_j,\mu_{j'}) \right\} \right\}
  }
  \f]

  Where \f$D\f$ is a distance function.

  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aipartlinknuminst_memberShip a clusters partition in a ds::PartitionLinkedNumInst data structure   
  \param aifunc2p_dist an object of type dist::Dist to calculate distances 
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE,
	   typename T_CLUSTERIDX,
           typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC
	   >
T_METRIC
CSmeasure
(INPUT_ITERATOR                         aiiterator_instfirst,
 const mat::MatrixRow<T_FEATURE>        &aimatrixt_centroids,
 const ds::PartitionLinkedNumInst
 <T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>   &aipartlinknuminst_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist
 )
{ //BEGIN CSmeasure
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::CSmeasure";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << " input const aimatrixt_centroids[" << &aimatrixt_centroids << "]\n"
	      << " input const &aipartlinknuminst_memberShip[" << &aipartlinknuminst_memberShip << "]\n"
              << " input const &aifunc2p_dist[" << &aifunc2p_dist << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lort_CSmeasure = measuare_undefCS(T_METRIC);

  T_METRIC lrt_sumMinCjCjp =
    sumMinCjCjp
    (aimatrixt_centroids,	  
     aifunc2p_dist
     );

  if ( lrt_sumMinCjCjp != 0  ) {

    const std::vector<T_METRIC>&&
      lvectort_diameterClusterK =
      diameterClusterK
      (aiiterator_instfirst,
       aipartlinknuminst_memberShip,
       aifunc2p_dist
       );

    const std::vector<T_INSTANCES_CLUSTER_K>  &lvectorit_numInstClusterK =
      aipartlinknuminst_memberShip.getVectorNumInstClusterK();
    
    lort_CSmeasure = 0.0;
    
    for ( uintidx lui_i = 0; lui_i < lvectort_diameterClusterK.size(); lui_i++) {

      if ( lvectorit_numInstClusterK.at(lui_i) != 0 ) {

	lort_CSmeasure += lvectort_diameterClusterK[lui_i]
	  / (T_METRIC) lvectorit_numInstClusterK[lui_i];
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " T_METRIC lort_CSmeasure = " << lort_CSmeasure
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lort_CSmeasure;

} //END CSmeasure

   
/*! \fn T_METRIC silhouette(mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity, ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &aipartlinknuminst_memberShip)
  \brief Silhouette \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006 
  \details Calculate the average of Silhouette for the \f$x_i\f$ instances
  \f[
  s(x_i)=\frac{b(x_i)-a(x_i)}{\max\left\{a(x_i),b(x_i)\right\} }
  \f]
  \param aimatrixtriagt_dissimilarity a triangular matrix with distances between instances
  \param aipartlinknuminst_memberShip a clusters partition in a ds::PartitionLinkedNumInst data structure

  \note Reduces the complexity of \f$O(d \cdot n^2)\f$ to \f$O(n^2)\f$, for data set d large better performance
*/
template < typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC
	   >
T_METRIC
silhouette
(mat::MatrixTriang<T_METRIC>                                    &aimatrixtriagt_dissimilarity,
 ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &aipartlinknuminst_memberShip
 ) 
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlinknuminst_memberShip.getNumPartitions();
  const std::vector<T_INSTANCES_CLUSTER_K>  &lvectorit_numInstClusterK =
    aipartlinknuminst_memberShip.getVectorNumInstClusterK();
  const T_CLUSTERIDX lcidx_numNullCluster =
    (T_CLUSTERIDX)
    std::count_if
    (lvectorit_numInstClusterK.begin(),
     lvectorit_numInstClusterK.end(),
     [] (const T_INSTANCES_CLUSTER_K aiit_num) {return aiit_num == T_INSTANCES_CLUSTER_K(0); }
     );
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::silhouette";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  PartitionLinked&: aipartlinknuminst_memberShip[" 
	      << &aipartlinknuminst_memberShip << "]\n"
	      << "lvectorit_numInstClusterK[" << &lvectorit_numInstClusterK << "]\n"; 
    inout::containerprint
      (lvectorit_numInstClusterK.begin(),
       lvectorit_numInstClusterK.end(),
       std::cout,
       lpc_labelFunc,
       ','
       );
   
    std::cout << "lcidx_numClusterK =  " << lcidx_numClusterK
	      << "\nlcidx_numNullCluster =  " << lcidx_numNullCluster 
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
   
  T_METRIC  lort_silhouette = measuare_undefSilhouette(T_METRIC);
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_i(&aipartlinknuminst_memberShip);
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_ip(&aipartlinknuminst_memberShip);

  if ( (lcidx_numClusterK - lcidx_numNullCluster) > 1)  {
    
    T_METRIC lrt_sumSilhouette = T_METRIC(0.0);

    for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) { 
      
      if ( lvectorit_numInstClusterK.at(lcidx_Ck) > 1 ) {
	  
	for ( literpart_i.begin(lcidx_Ck); literpart_i.end(); literpart_i.next() ) {
      
	  T_METRIC lt_sumai = T_METRIC(0);
	  
	  for ( literpart_ip.begin(lcidx_Ck); literpart_ip.end(); literpart_ip.next() ) {
	   
	    lt_sumai +=
	      aimatrixtriagt_dissimilarity
	      (literpart_i.getValue(),
	       literpart_ip.getValue()
	       );
	 
	  }
	  T_METRIC lt_ai = lt_sumai / T_METRIC(lvectorit_numInstClusterK.at(lcidx_Ck) - 1);


	  T_METRIC lt_bi = std::numeric_limits<T_METRIC>::max();  
	  for ( T_CLUSTERIDX lcidx_Ckp = 0; lcidx_Ckp < lcidx_numClusterK; lcidx_Ckp++) {
	    if ( lcidx_Ck != lcidx_Ckp ) {
	      T_METRIC lt_sumbi = T_METRIC(0);
	      for ( literpart_ip.begin(lcidx_Ckp); literpart_ip.end(); literpart_ip.next() ) {
		
		lt_sumbi +=
		  aimatrixtriagt_dissimilarity
		  (literpart_i.getValue(),
		   literpart_ip.getValue()
		   );
	      }
	      T_METRIC lt_diC = lt_sumbi / T_METRIC(lvectorit_numInstClusterK.at(lcidx_Ckp));
	      if ( lt_diC < lt_bi ) lt_bi = lt_diC;
	    }
	  }

	  T_METRIC lt_max =  std::max(lt_ai,lt_bi);

#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose < 0 &&  geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << lpc_labelFunc << ": i = " << literpart_i.getValue() << " a " <<  lt_ai << " b " << lt_bi << std::endl;
	  }
#endif //__VERBOSE_YES
 
	  lrt_sumSilhouette += (lt_max == 0.0)?0.0:(lt_bi - lt_ai)/ lt_max;
      
    
	} //For Object_i

	lort_silhouette = lrt_sumSilhouette / (T_METRIC) aipartlinknuminst_memberShip.getNumInstances();
	
      }  

    } // For all cluster

  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lort_silhouette = "   << lort_silhouette
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lort_silhouette;
  
}


/*! \fn T_METRIC silhouette (INPUT_ITERATOR aiiterator_instfirst, ds::PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &aipartlinknuminst_memberShip, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Silhouette \cite Kaufman:Rousseeuw:Book:ClusterAnalysis:1990
  \details Calculate the average of Silhouette for the instances
  \f[
  s(x_i)=\frac{b(x_i)-a(x_i)}{\max\left\{a(x_i),b(x_i)\right\} }
  \f]
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aipartlinknuminst_memberShip a clusters partition in a ds::PartitionLinkedNumInst data structure
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
	   typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC,
	   typename T_FEATURE
	   >
T_METRIC
silhouette
(INPUT_ITERATOR                          aiiterator_instfirst,
 ds::PartitionLinkedNumInst
 <T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>    &aipartlinknuminst_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>    &aifunc2p_dist
 ) 
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlinknuminst_memberShip.getNumPartitions();

  const T_CLUSTERIDX lcidx_numNullCluster =
    aipartlinknuminst_memberShip.getNumNullCluster();
  /*   (T_CLUSTERIDX)
       std::count_if
       (aivectorit_numInstClusterK.begin(),
       aivectorit_numInstClusterK.end(),
       [] (const T_INSTANCES_CLUSTER_K aiit_num) {return aiit_num == T_INSTANCES_CLUSTER_K(0);}
       );
  */
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::silhouette";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout 
      << lpc_labelFunc
      << ":  IN(" << geiinparam_verbose << ")\n"
      << "(input  PartitionLinked&: aipartlinknuminst_memberShip[" 
      << &aipartlinknuminst_memberShip << "]\n";
    aipartlinknuminst_memberShip.print();
    std::cout 
      << "\ninput  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
      << &aifunc2p_dist << "]\n"
      << "lcidx_numClusterK =  " << lcidx_numClusterK
      << "\nlcidx_numNullCluster =  " << lcidx_numNullCluster 
      << "\n)"
      << std::endl;
  }
#endif /*__VERBOSE_YES*/
   
  T_METRIC  lort_silhouette = measuare_undefSilhouette(T_METRIC); 
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_i(&aipartlinknuminst_memberShip);
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_ip(&aipartlinknuminst_memberShip);

  if ( (lcidx_numClusterK - lcidx_numNullCluster) > 1)  {
   
    T_METRIC lrt_sumSilhouette = T_METRIC(0.0);

    //FOR ALL CLUSTER
    for ( T_CLUSTERIDX lcidx_Ck = 0; lcidx_Ck < lcidx_numClusterK; lcidx_Ck++) { 
      
      /*INSTANCE NOT IS OUTLIER*/
      //if ( aivectorit_numInstClusterK.at(lcidx_Ck) > 1 ) {
      if ( aipartlinknuminst_memberShip.getNumInstClusterK(lcidx_Ck) > 1 ) {
	  
	for ( literpart_i.begin(lcidx_Ck); literpart_i.end(); literpart_i.next() ) {
      
	  data::Instance<T_FEATURE>* linst_i
	    = *std::next(aiiterator_instfirst,literpart_i.getValue());
	  T_METRIC lt_sumai = T_METRIC(0);
	  for ( literpart_ip.begin(lcidx_Ck); literpart_ip.end(); literpart_ip.next() ) {
	
	    data::Instance<T_FEATURE>* linst_l
	      = *std::next(aiiterator_instfirst,literpart_ip.getValue());
	    
	    T_METRIC lt_distai =
	      aifunc2p_dist
	      (linst_i->getFeatures(),
	       linst_l->getFeatures(),
	       data::Instance<T_FEATURE>::getNumDimensions()
	       );

	    lt_sumai += lt_distai;

#ifdef __VERBOSE_YES
	    if ( geiinparam_verbose < 0 &&  geiinparam_verbose <= geiinparam_verboseMax ) {
	      std::cout << lpc_labelFunc << ": D( " <<  literpart_i.getValue() << ", " <<  literpart_ip.getValue() << ") = "
			<< lt_distai << std::endl;
	    }
#endif //__VERBOSE_YES
	      
	  }
	  T_METRIC lt_ai = lt_sumai / T_METRIC(aipartlinknuminst_memberShip.getNumInstClusterK(lcidx_Ck) -1);


	  T_METRIC lt_bi = std::numeric_limits<T_METRIC>::max();  
	  for ( T_CLUSTERIDX lcidx_Ckp = 0; lcidx_Ckp < lcidx_numClusterK; lcidx_Ckp++) {
	    if ( lcidx_Ck != lcidx_Ckp ) {
	      T_METRIC lt_sumbi = T_METRIC(0);
	      for ( literpart_ip.begin(lcidx_Ckp); literpart_ip.end(); literpart_ip.next() ) {
		data::Instance<T_FEATURE>* linst_l
		  = *std::next(aiiterator_instfirst,literpart_ip.getValue());
		lt_sumbi +=
		  aifunc2p_dist
		  (linst_i->getFeatures(),
		   linst_l->getFeatures(),
		   data::Instance<T_FEATURE>::getNumDimensions()
		   );
	      }
	      T_METRIC lt_diC = lt_sumbi / aipartlinknuminst_memberShip.getNumInstClusterK(lcidx_Ckp);
	      if ( lt_diC < lt_bi ) lt_bi = lt_diC;
	    }	
	  }

#ifdef __VERBOSE_YES
	  if ( geiinparam_verbose < 0 &&  geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << lpc_labelFunc << ": i = " << literpart_i.getValue() << " a " <<  lt_ai << " b " << lt_bi << std::endl; 
	  }
#endif //__VERBOSE_YES
	  
	  T_METRIC lt_max =  std::max(lt_ai,lt_bi);
      
	  lrt_sumSilhouette += (lt_max == 0.0)?0.0:(lt_bi - lt_ai)/ lt_max;
      
    
	} //FOR OBJETO I

	lort_silhouette = lrt_sumSilhouette / (T_METRIC) aipartlinknuminst_memberShip.getNumInstances();
	
      }  

    } //ALL CLUSTER  

  }
   
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lort_silhouette = "   << lort_silhouette
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lort_silhouette;
  
}


/*! \fn T_METRIC silhouette(const uintidx aiui_idxInstance, const T_CLUSTERIDX aicidx_instInClusterJ, mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity, ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_memberShip, const std::vector<T_INSTANCES_CLUSTER_K> &aivectorit_numInstClusterK)
  \brief Silhouette \cite Kaufman:Rousseeuw:Book:ClusterAnalysis:1990
  \details Calculate Silhouette for an instance
  \f[
  s(x_i)=\frac{b(x_i)-a(x_i)}{\max\left\{a(x_i),b(x_i)\right\} }
  \f]
    
  \param aiui_idxInstance Instance index to calculate \f$ s(x_i) \f$
  \param aicidx_instInClusterJ Index of the cluster to which the aiui_idxInstance belongs
  \param aimatrixtriagt_dissimilarity a matrix of distances
  \param aipartlink_memberShip a clusters partition in a ds::PartitionLinked data structure
  \param aivectorit_numInstClusterK a vector with the number of instances per cluster
*/
template < typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC
	   >
T_METRIC
silhouette
(const uintidx                             aiui_idxInstance,
 const T_CLUSTERIDX                        aicidx_instInClusterJ,
 mat::MatrixTriang<T_METRIC>               &aimatrixtriagt_dissimilarity,
 ds::PartitionLinked<T_CLUSTERIDX>         &aipartlink_memberShip,
 const std::vector<T_INSTANCES_CLUSTER_K>  &aivectorit_numInstClusterK
 ) 
{
  const T_CLUSTERIDX lcidx_numClusterK = aipartlink_memberShip.getNumPartitions();

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::silhouette";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  PartitionLinked&: aipartlink_memberShip[" 
	      << &aipartlink_memberShip << "]\n"
	      << "aivectorit_numInstClusterK[" << &aivectorit_numInstClusterK << "]\n"; 
    inout::containerprint
      (aivectorit_numInstClusterK.begin(),
       aivectorit_numInstClusterK.end(),
       std::cout,
       lpc_labelFunc,
       ','
       );

    std::cout << "aiui_idxInstance =  " << aiui_idxInstance
	      << "\taicidx_instInClusterJ =  " << aicidx_instInClusterJ
	      << "\tlcidx_numClusterK =  " << lcidx_numClusterK
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
   
  T_METRIC  lort_silhouette = measuare_undefSilhouette(T_METRIC);
  
  ds::IteratorPartitionLinked <T_CLUSTERIDX>
    literpart_ip(&aipartlink_memberShip);
      
  if ( aivectorit_numInstClusterK.at(aicidx_instInClusterJ) > 1 ) {

    T_METRIC lt_sumai = T_METRIC(0);
	  
    for ( literpart_ip.begin(aicidx_instInClusterJ); literpart_ip.end(); literpart_ip.next() ) {
		    
      lt_sumai +=
	aimatrixtriagt_dissimilarity
	(aiui_idxInstance,
	 literpart_ip.getValue()
	 );
	 
    }
    T_METRIC lt_ai = lt_sumai / ( aivectorit_numInstClusterK.at(aicidx_instInClusterJ) + 1);
    T_METRIC lt_bi = std::numeric_limits<T_METRIC>::max();  
    for ( T_CLUSTERIDX lcidx_Ckp = 0; lcidx_Ckp < lcidx_numClusterK; lcidx_Ckp++) {
      if ( aicidx_instInClusterJ != lcidx_Ckp ) {
	T_METRIC lt_sumbi = T_METRIC(0);
	for ( literpart_ip.begin(lcidx_Ckp); literpart_ip.end(); literpart_ip.next() ) {
	  lt_sumbi +=
	    aimatrixtriagt_dissimilarity
	    (aiui_idxInstance,
	     literpart_ip.getValue()
	     );
	}
	T_METRIC lt_diC = lt_sumbi / aivectorit_numInstClusterK.at(lcidx_Ckp);
	if ( lt_diC < lt_bi ) lt_bi = lt_diC;
      }		
    } //FOR OTHER CLUSTER

    T_METRIC lt_max =  std::max(lt_ai,lt_bi);
    lort_silhouette = (lt_max == 0.0)?0.0:(lt_bi - lt_ai)/ lt_max;  

  } //IF 

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "T_METRIC lort_silhouette = "   << lort_silhouette
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lort_silhouette;
  
}


/*! \fn std::vector<T_METRIC> simplifiedSilhouette(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const std::vector<T_INSTANCES_CLUSTER_K> &aivectorit_numInstClusterK, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief Simplified Silhouette \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006 
  \details The silhouette proposed in \cite Kaufman:Rousseeuw:Book:ClusterAnalysis:1990 depends on the computation of all distances between objects, leading to a computational cost of \f$O(n^2)\f$, which is often not sufficiently efficient for real-world clustering applications (e.g. data mining, text mining, gene-expression data analysis). To circumvent this limitation, a simplified silhouette can be employed. The simplified silhouette is based on the computation of distances between objects and cluster centroids, which are the mean vectors of the clusters \cite Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003.
  \param aimatrixt_centroids a mat::MatrixBase with centroids
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition of instances in clusters
  \param aivectorit_numInstClusterK a vector with the number of instances per cluster
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename INPUT_ITERATOR,
           typename T_CLUSTERIDX,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC
	   >
std::vector<T_METRIC> 
simplifiedSilhouette
(const mat::MatrixBase<T_FEATURE>          &aimatrixt_centroids, 
 INPUT_ITERATOR                            aiiterator_instfirst,
 const INPUT_ITERATOR                      aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>        &aipartition_clusters,
 const std::vector<T_INSTANCES_CLUSTER_K>  &aivectorit_numInstClusterK,
 const dist::Dist<T_METRIC,T_FEATURE>      &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::simplifiedSilhouette";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
 
  std::vector<T_METRIC>  lovectort_partialSilhouette
    (aimatrixt_centroids.getNumRows(),T_METRIC(0));

  interfacesse::copya
    (lovectort_partialSilhouette.data(), 
     (T_METRIC) 0.0, 
     lovectort_partialSilhouette.size()
     );

  const T_CLUSTERIDX lcidx_numClusterK =
    (T_CLUSTERIDX) aimatrixt_centroids.getNumRows();
  const T_CLUSTERIDX lcidx_numNullCluster =
    (T_CLUSTERIDX)
    std::count_if
    (aivectorit_numInstClusterK.begin(),
     aivectorit_numInstClusterK.end(),
     [] (const T_INSTANCES_CLUSTER_K aiit_num) {return aiit_num == T_INSTANCES_CLUSTER_K(0);}
     );

  if ( (lcidx_numClusterK - lcidx_numNullCluster) > 1)  {
      
    for (aipartition_clusters.begin();
	 aiiterator_instfirst != aiiterator_instlast;
	 ++aiiterator_instfirst)
      {
      
	T_CLUSTERIDX lcidx_xinK = aipartition_clusters.next();
	
	/*INSTANCE NOT IS OUTLIER
	 */
	if ( (0 <= lcidx_xinK)  && ( lcidx_xinK < lcidx_numClusterK)  ) {
      
	  /* IF  NOT Ci IS A SINGLETON
	   */ 
	  if ( aivectorit_numInstClusterK[lcidx_xinK] > 1 ) {

	    data::Instance<T_FEATURE>* linst_inter =
	      (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
	  
	    T_METRIC lt_a = 
	      aifunc2p_dist
	      (aimatrixt_centroids.getRow(lcidx_xinK),
	       linst_inter->getFeatures(),
	       data::Instance<T_FEATURE>::getNumDimensions()
	       );
	  
	    T_METRIC lt_b = 
	      std::numeric_limits<T_METRIC>::max();
	    for (T_CLUSTERIDX lcidx_k= 0; lcidx_k < lcidx_numClusterK; lcidx_k++)  {
	      if (lcidx_k !=  lcidx_xinK  ) {
		T_METRIC lt_distik =
		  aifunc2p_dist
		  (aimatrixt_centroids.getRow(lcidx_k),
		   linst_inter->getFeatures(), 
		   data::Instance<T_FEATURE>::getNumDimensions()
		   );
		if (lt_b > lt_distik)
		  lt_b = lt_distik;	
	      }
	    }
	  
	    lovectort_partialSilhouette[lcidx_xinK] += 
	      (T_METRIC) (lt_b - lt_a) / (T_METRIC) std::max(lt_a,lt_b);
	  }  else if  ( aivectorit_numInstClusterK[lcidx_xinK] == 1 ) {
	    lovectort_partialSilhouette[lcidx_xinK] = (T_METRIC) 0.0;
	  }
	}
      } /*End for lui_i*/
  
    for (uintidx lui_j = 0; lui_j <aimatrixt_centroids.getNumRows(); lui_j++)  {
      if ( aivectorit_numInstClusterK[lui_j] > 0 ) {
	lovectort_partialSilhouette[lui_j] /= 
	  (T_METRIC) aivectorit_numInstClusterK[lui_j];
      }     
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

    std::ostringstream lostrstream_labelPartialSilhouette;
    lostrstream_labelPartialSilhouette  << "<PARTIALSILHUETTE: " << lpc_labelFunc
					<< "lovectort_partialSilhouette["
					<< &lovectort_partialSilhouette << ']';
    inout::containerprint
      (lovectort_partialSilhouette.begin(),
       lovectort_partialSilhouette.end(),
       std::cout,
       lostrstream_labelPartialSilhouette.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lovectort_partialSilhouette;
}
  

/*! \fn T_METRIC dbindex(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, const std::vector<T_METRIC>  &aivectort_scatter, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief DB Index \cite Davies:Bouldin:Metricclustering:DB:1979 \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002
  \details This index is a function of the ratio of the sum of within-cluster scatter to between-cluster separation
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aivectort_scatter a vector whith aivectort_scatter of the clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_METRIC
	   >
T_METRIC
dbindex
(const mat::MatrixRow<T_FEATURE>      &aimatrixt_centroids,
 const std::vector<T_METRIC>          &aivectort_scatter, //Radius of Clusters
 const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::dbindex<>  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << ']' << ": numrows "
	      <<  aimatrixt_centroids.getNumRows() << '\n'
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  T_METRIC loT_dbindex = measuare_undefDBindex(T_METRIC);
  uintidx lui_numClusterK = 0;
  
  
  if ( aimatrixt_centroids.getNumRows() > 1 ) {
           
    T_METRIC lrt_sumdbindex = 0.0;
      
    for (uintidx lui_i = 0; lui_i < aimatrixt_centroids.getNumRows(); lui_i++) {
	
      if ( std::isnan(aivectort_scatter[lui_i]) )
	continue;
      ++lui_numClusterK;
      uintidx lui_j = 0;
      while ( (lui_j < aivectort_scatter.size())
	      && ((lui_i == lui_j) || std::isnan(aivectort_scatter[lui_j]))
	      ) ++lui_j;
      if ( lui_j >= aivectort_scatter.size() )
	continue;
      T_METRIC lt_dijt = (T_METRIC) aifunc2p_dist
	(aimatrixt_centroids.getRow(lui_i),
	 aimatrixt_centroids.getRow(lui_j),
	 aimatrixt_centroids.getNumColumns()
	 );
      if ( lt_dijt != 0.0) {
	T_METRIC lt_rMax = ( aivectort_scatter[lui_i] + aivectort_scatter[lui_j] ) / lt_dijt;
	  
	while ( ++lui_j < aimatrixt_centroids.getNumRows() ) {
	  if ( lui_i != lui_j && !std::isnan(aivectort_scatter[lui_j]) ) {
	      
	    lt_dijt =
	      aifunc2p_dist
	      (aimatrixt_centroids.getRow(lui_i),
	       aimatrixt_centroids.getRow(lui_j),
	       aimatrixt_centroids.getNumColumns()
	       );
	    if ( lt_dijt != 0.0) {
	      T_METRIC lt_r = ( aivectort_scatter[lui_i] + aivectort_scatter[lui_j] ) / lt_dijt;
		
	      if ( lt_r > lt_rMax ) {
		lt_rMax =  lt_r;
	      }
		
	    } 
		
	  } 
	    
	} /* while j */ 
	lrt_sumdbindex += lt_rMax;
      }
    } /* while i*/
    loT_dbindex  = (lrt_sumdbindex > 0.0)?  
      lrt_sumdbindex / (T_METRIC) lui_numClusterK: 
      measuare_undefDBindex(T_METRIC); 
  }
     
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::dbindex<>  OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\tlui_numClusterK = " << lui_numClusterK 
	      << "\tloT_dbindex = " << loT_dbindex 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return loT_dbindex;
}



/*! \fn T_METRIC dbindex(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief DB Index \cite Davies:Bouldin:Metricclustering:DB:1979 \cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002
  \details This index is a function of the ratio of the sum of within-cluster scatter to between-cluster separation
  \param aimatrixt_centroids a mat::MatrixBase with centroids
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/  
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
dbindex
(const mat::MatrixBase<T_FEATURE>     &aimatrixt_centroids,
 INPUT_ITERATOR                       aiiterator_instfirst,
 const INPUT_ITERATOR                 aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>   &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::dbindex";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC loT_dbindex = measuare_undefDBindex(T_METRIC);
  uintidx lui_numClusterK = 0;
    
  if ( aimatrixt_centroids.getNumRows() > 1 ) {

    std::vector<T_METRIC>&& lvectort_scatter =
      avgRadiusClusterK
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_dist
       );
                 
    T_METRIC lrt_sumdbindex = 0.0;
      
    for (uintidx lui_i = 0; lui_i < aimatrixt_centroids.getNumRows(); lui_i++) {
	
      if ( std::isnan(lvectort_scatter[lui_i]) )
	continue;
      ++lui_numClusterK;
      uintidx lui_j = 0;
      while ( (lui_j < lvectort_scatter.size())
	      && ((lui_i == lui_j) || std::isnan(lvectort_scatter[lui_j]))
	      ) ++lui_j;
      if ( lui_j >= lvectort_scatter.size() )
	continue;
      T_METRIC lt_dijt = (T_METRIC) aifunc2p_dist
	(aimatrixt_centroids.getRow(lui_i),
	 aimatrixt_centroids.getRow(lui_j),
	 aimatrixt_centroids.getNumColumns()
	 );
      if ( lt_dijt != 0.0) {
	T_METRIC lt_rMax = ( lvectort_scatter[lui_i] + lvectort_scatter[lui_j] ) / lt_dijt;
	  
	while ( ++lui_j < aimatrixt_centroids.getNumRows() ) {
	  if ( lui_i != lui_j && !std::isnan(lvectort_scatter[lui_j]) ) {
	      
	    lt_dijt =
	      aifunc2p_dist
	      (aimatrixt_centroids.getRow(lui_i),
	       aimatrixt_centroids.getRow(lui_j),
	       aimatrixt_centroids.getNumColumns()
	       );
	    if ( lt_dijt != 0.0) {
	      T_METRIC lt_r = ( lvectort_scatter[lui_i] + lvectort_scatter[lui_j] ) / lt_dijt;
		
	      if ( lt_r > lt_rMax ) {
		lt_rMax =  lt_r;
	      }
		
	    } 
		
	  } 
	    
	} /* while j */ 
	lrt_sumdbindex += lt_rMax;
      }
    } /* while i*/
    loT_dbindex  = (lrt_sumdbindex > 0.0)?  
      lrt_sumdbindex / (T_METRIC) lui_numClusterK: 
      measuare_undefDBindex(T_METRIC); 
  }
     
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::dbindex<>  OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\tlui_numClusterK = " << lui_numClusterK 
	      << "\tloT_dbindex = " << loT_dbindex 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return loT_dbindex;
}


/*! \fn std::pair<T_METRIC,bool> SSE(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief SSE \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002 \cite Chang:etal:GAclustering:GAGR:2009. A boolean is also returned to indicate if the partition is valid
  \details  SSE A common clustering criterion or quality indicator is the sum of squared error SSE  measure

  \f[
  SSE=\sum_{C_j}\sum_{x_i\in C_j}(x_i-\mu_j)^{T}(x_i-\mu_j)=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i-\mu_j\Vert^{2}.
  \f]

  Or with some slight variation Sum of Euclidean Distance (SED):

  \f[
  SED=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i-\mu_j\Vert
  \f]

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances

  \note The membership of an instance to a cluster for this function is determined by the nearest object rule
*/
template < typename INPUT_ITERATOR,
           typename T_METRIC, 
	   typename T_FEATURE
       	   >
std::pair<T_METRIC,bool>
SSE
(const mat::MatrixRow<T_FEATURE>          &aimatrixt_centroids,
 INPUT_ITERATOR                           aiiterator_instfirst,
 const INPUT_ITERATOR                     aiiterator_instlast,
 const dist::Dist<T_METRIC,T_FEATURE>     &aifunc2p_dist
 )
{  
  T_METRIC              loT_SSE;
  T_METRIC              lT_distMinCentInst;

  loT_SSE = T_METRIC(0);
  std::vector<uintidx> 
    lvectorT_numInstancesClusterK(aimatrixt_centroids.getNumRows(),0);

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::SSEWithOutMember:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t input  mat::MatrixRow<T_FEATURE>: aimatrixt_centroids[" 
	      << &aimatrixt_centroids << ']'
	      << "\n\t input aiiterator_instfirst[" << *aiiterator_instfirst << ']'
	      << "\n\t input aiiterator_instlast[" <<  *aiiterator_instlast << ']'
	      << "\t input  aifunc2p_dist\n"
	      << "\t)\n";
  }
#endif //__VERBOSE_YES
  
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {

    data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    
    intidx lcidx_instInClusterJ = 
      nearest::NN
      <intidx,
       T_FEATURE,
       T_METRIC
       >
      (lT_distMinCentInst,
       aimatrixt_centroids,
       linst_inter->getFeatures(),
       aifunc2p_dist
       );
    loT_SSE += (T_METRIC) lT_distMinCentInst;
    lvectorT_numInstancesClusterK[lcidx_instInClusterJ]++;
       

  }
  
  intidx  locidx_numClusterNull = (intidx)
    std::count
    (lvectorT_numInstancesClusterK.begin(),
     lvectorT_numInstancesClusterK.end(),
     uintidx(0)
     );
   
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::SSEWithOutMember: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t output T_METRIC: loT_SSE = " << loT_SSE << '\n'
	      << "\t output T_CLUSTERIDX  locidx_numClusterNull = " 
	      << locidx_numClusterNull << '\n';
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_pair(loT_SSE,locidx_numClusterNull==intidx(0));
} /*END SSEWithOutMember*/



/*! \fn std::pair<T_METRIC,bool> SSE (const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist)

  \brief SSE \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002 \cite Chang:etal:GAclustering:GAGR:2009. A boolean is also returned to indicate if the partition is valid
  \details  SSE A common clustering criterion or quality indicator is the sum of squared error SSE  measure

  \f[
  SSE=\sum_{C_j}\sum_{x_i\in C_j}(x_i-\mu_j)^{T}(x_i-\mu_j)=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i-\mu_j\Vert^{2}.
  \f]

  Or with some slight variation Sum of Euclidean Distance (SED):

  \f[
  SED=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i-\mu_j\Vert
  \f]

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/ 
template < typename INPUT_ITERATOR,
	   typename T_METRIC, 
	   typename T_FEATURE,
	   typename T_CLUSTERIDX /*-1,0,..,K*/
	   >
std::pair<T_METRIC,bool>
SSE
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::SSEWithMember:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t input  mat::MatrixRow<T_FEATURE>&: aimatrixt_centroids[" 
	      << &aimatrixt_centroids << ']'
	      << "\n\t input  partition::Partition<T_CLUSTERIDX>&: aipartition_clusters[" 
	      << &aipartition_clusters << ']'
	      << "\n\t input aiiterator_instfirst[" << *aiiterator_instfirst << ']'
	      << "\n\t input aiiterator_instlast[" << *aiiterator_instlast << ']'
	      << "\n\t input  aifunc2p_dist[" << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<uintidx>
    lvectorT_numInstancesClusterK(aimatrixt_centroids.getNumRows(),0);
  
  T_CLUSTERIDX lcidx_numClusterK = 
    (T_CLUSTERIDX) aimatrixt_centroids.getNumRows();

  T_METRIC  loT_SSE =  T_METRIC(0);
  
  for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    T_CLUSTERIDX lcidx_instInClusterJ = aipartition_clusters.next();
  
    if ( 0 <= lcidx_instInClusterJ  && lcidx_instInClusterJ <  lcidx_numClusterK  ) {

      data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
      
      T_METRIC  lT_distMinCentInst =  
	aifunc2p_dist
	(aimatrixt_centroids.getRow(lcidx_instInClusterJ),
	 linst_inter->getFeatures(),
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );

      loT_SSE += lT_distMinCentInst;

      lvectorT_numInstancesClusterK[lcidx_instInClusterJ]++;

    }
    else {
      loT_SSE = measuare_undefSSE(T_METRIC);
      break;
    }
    //++luintidx_insti;
  }

  T_CLUSTERIDX  locidx_numClusterNull =
    (T_CLUSTERIDX) std::count
    (lvectorT_numInstancesClusterK.begin(),
     lvectorT_numInstancesClusterK.end(),
     uintidx(0)
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::SSEWithMember: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t output T_METRIC: loT_SSE = " << loT_SSE << '\n'
	      << "\t output T_CLUSTERIDX  locidx_numClusterNull = " 
	      << locidx_numClusterNull
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_pair(loT_SSE,locidx_numClusterNull==T_CLUSTERIDX(0));

} /*END SSEWithOutMember*/



/*! \fn T_METRIC SSE (const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, T_CLUSTERIDX *aiarraycidx_memberShip, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist)
  \brief SSE \cite Bandyopadhyay:Maulik:GAclustering:KGA:2002 \cite Chang:etal:GAclustering:GAGR:2009. Also known as Total within cluster variation (TWCV) \cite Murthy:Chowdhury:GAclustering:GA:1996 \cite Lu:etal:GAclusteringLabel:FGKA:2004
  \details  SSE A common clustering criterion or quality indicator is the sum of squared error SSE  measure

  \f[
  SSE=\sum_{C_j}\sum_{x_i\in C_j}(x_i-\mu_j)^{T}(x_i-\mu_j)=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i-\mu_j\Vert^{2}.
  \f]

  Or with some slight variation Sum of Euclidean Distance (SED):

  \f[
  SED=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i-\mu_j\Vert
  \f]

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aiarraycidx_memberShip an array with the labels belong to the cluster
  \param aifunc2p_dist an object of type dist::Dist to calculate distances

  \note The membership of an instance to a cluster for this function is determined by aiarraycidx_memberShip
*/
template < typename INPUT_ITERATOR,
           typename T_CLUSTERIDX,
	   typename T_FEATURE,
	   typename T_METRIC
	   >
T_METRIC
SSE
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                       aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 T_CLUSTERIDX                          *aiarraycidx_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::SSE";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(input  mat::MatrixRow<T_CENTROIDS>: aimatrixt_centroids["
	      <<  &aimatrixt_centroids << "]\n"
	      << "\t input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "\t input  T_CLUSTERIDX*: aiarraycidx_memberShip[" 
	      << aiarraycidx_memberShip << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  T_METRIC           lot_SSE;
  
  lot_SSE = T_METRIC(0.0);
       
  for (; aiiterator_instfirst != aiiterator_instlast;
       ++aiiterator_instfirst, ++aiarraycidx_memberShip)
    {
      data::Instance<T_FEATURE>* linst_inter =
	(data::Instance<T_FEATURE>*) *aiiterator_instfirst;
  
      const T_FEATURE* lmatrowt_centroid =  aimatrixt_centroids.getRow(*aiarraycidx_memberShip);
    
      lot_SSE += 
	aifunc2p_dist
	(lmatrowt_centroid,
	 linst_inter->getFeatures(),
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
    }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
  	      << " T_METRIC: lot_SSE = " 
	      << lot_SSE
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lot_SSE;
  
}

/*! \fn std::pair<T_METRIC,T_CLUSTERIDX> distortion(const mat::MatrixBase<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, T_CLUSTERIDX  aiarraycidx_memberShip, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist, const FUNCINSTFREQUENCY func_instfrequency )
  \brief Distortion for instances with frequency. A boolean is also returned to indicate if the partition is valid. if there is inconsistency of instances when they belong to a cluster, returns an undefined value for the metri. 
  \details Similar to SSE but considers the number of instances and attributes

  \f[
  distortion(C) = \frac{1}{n \cdot d} \sum_{i =j}^n D(x_i,f_C(x_i))^2 
  \f]

  Where \f$D\f$ is the Euclidean distance, \f$f_C(x_i)\f$ be a mapping which gives the 
  closest centroid in solution \f$C\f$ for a instance  \f$x_i\f$, \f$n\f$ is the number 
  of instances, and \f$d\f$ number of attributes of the instances.

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aiarraycidx_memberShip an array with the labels belong to the cluster
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param func_instfrequency a function to obtain the frequency of the instance
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename FUNCINSTFREQUENCY, 
	   typename T_METRIC  
	   >
std::pair<T_METRIC,bool>
distortion
(const mat::MatrixBase<T_FEATURE>      &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 T_CLUSTERIDX                          *aiarraycidx_memberShip,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist,
 const FUNCINSTFREQUENCY               func_instfrequency
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::distortion";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(input  mat::MatrixRow<T_CENTROIDS>: aimatrixt_centroids["
	      <<  &aimatrixt_centroids << "]\n"
	      << "\t input  INPUT_ITERATOR    aiiterator_instfirst[" 
	      <<  &aiiterator_instfirst << "]\n"
	      << "\t input  T_CLUSTERIDX*: aiarraycidx_memberShip[" 
	      << aiarraycidx_memberShip << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::vector<bool> lvectorbool_haveInstClusterK(aimatrixt_centroids.getNumRows(),false);
  T_METRIC          lort_SSE = T_METRIC(0);
  T_METRIC          lT_numInstances  = T_METRIC(0);
  
  uintidx lui_i = 0;
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {

    data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    T_METRIC          lT_distMinCentInst;
      
    T_CLUSTERIDX lcidx_instInClusterJ = 
      nearest::NN
      (lT_distMinCentInst,
       aimatrixt_centroids,
       aiarraycidx_memberShip[lui_i],
       linst_inter->getFeatures(),
       aifunc2p_dist
       );

    if ( lcidx_instInClusterJ != aiarraycidx_memberShip[lui_i]) {

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << lpc_labelFunc
		  << ": OUT(" << geiinparam_verbose << ')'
		  << " Warning: Inconsistency of instances when they belong to a cluster" 
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
      return std::make_pair(measuare_undefObjetiveFunc(T_METRIC),false);

    }
          
    T_METRIC lrt_instfrequency = func_instfrequency(linst_inter);
    lort_SSE += lT_distMinCentInst * lrt_instfrequency; 
    lT_numInstances += lrt_instfrequency;     
    lvectorbool_haveInstClusterK[lcidx_instInClusterJ] =  true;
    ++lui_i;
    
  } /*END FOR*/

  T_CLUSTERIDX  locidx_numClusterNull =
    (T_CLUSTERIDX)
    std::count
    (lvectorbool_haveInstClusterK.begin(),
     lvectorbool_haveInstClusterK.end(),
     false);
  
  lort_SSE = ( lT_numInstances > 0 )?
    lort_SSE /  (  T_METRIC(lT_numInstances) * T_METRIC(data::Instance<T_FEATURE>::getNumDimensions() ) )
    :measuare_undefSSE(T_METRIC);
   
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
  	      << " T_METRIC: lort_SSE = " 
	      << lort_SSE
	      << "\tT_CLUSTERIDX  locidx_numClusterNull = " 
	      << locidx_numClusterNull 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  return std::make_pair(lort_SSE,locidx_numClusterNull== T_CLUSTERIDX(0));
  
}


/*! \fn T_METRIC distortion (T_FEATURE *aiarrayt_centroid1, T_INSTANCES_CLUSTER_K aiit_numInstClusterK1, T_FEATURE *aiarrayt_centroid2, T_INSTANCES_CLUSTER_K aiit_numInstClusterK2, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief  Distortion metric distance between two centroids \cite Franti:etal:GAclustering:gafranti:1997.
  \details 
  \param aiarrayt_centroid1 an array with the coordinates of the centroid of the first cluster
  \param aiit_numInstClusterK1 number of instances of the first group
  \param aiarrayt_centroid2 an array with the coordinates of the centroid of the second cluster
  \param aiit_numInstClusterK2 number of instances of the second group
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_METRIC
	   >
T_METRIC       
distortion
(T_FEATURE                            *aiarrayt_centroid1,
 T_INSTANCES_CLUSTER_K                aiit_numInstClusterK1,
 T_FEATURE                            *aiarrayt_centroid2,
 T_INSTANCES_CLUSTER_K                aiit_numInstClusterK2,
 const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 )
{
  T_METRIC lrt_distCjCi;
  T_METRIC lrt_factor;
  
  if ( aiit_numInstClusterK1 == 0 || aiit_numInstClusterK2 == 0 ) {
    lrt_distCjCi = 0.0;
  }
  else {
    lrt_distCjCi = 
      aifunc2p_dist
      (aiarrayt_centroid1,
       aiarrayt_centroid2,
       data::Instance<T_FEATURE>::getNumDimensions() 
       );
    lrt_factor  = T_METRIC(aiit_numInstClusterK1 * aiit_numInstClusterK2);
    lrt_factor /=  T_METRIC(aiit_numInstClusterK1 + aiit_numInstClusterK2);
    lrt_distCjCi *= lrt_factor; 
  }
  
  return lrt_distCjCi;
  
}
  
/*! \fn T_METRIC j1(BitMatrix<T_BITSIZE> &aicrispmatrix_partition, mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast,, dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief  \f$J_1\f$ \cite Ball:Hall:ClusterAnalysis:ISODATA:1967
  \details Least-squared errors functional for crisp partition or hard partition

  \f[
  J_1(U,\mu)=\sum_{i=1}^{n}\sum_{j=1}^{k}u_{ji}D_{Ind}(x_{i}-\mu_{i})
  \f]

  Where

  \f$U\in M_{c}\f$ crisp partition of \f$X\f$;

  \f[ 
  M_c = \{U_{k\times n} | u_{ji} \in \{0,1\}; \sum_{i=1}^{n}u_{ji} > 0, \hbox{for all } j \hbox{, }
  \sum_{i=1}^{n}u_{ji} = 1 \hbox{, for all } i \}
  \f]

  \f$\mu=[\mu_{1},\mu_{2},...,\mu_{k}]\f$ centroids,

  \f$D_{Ind}(x_{k},\mu_{i})\f$  is one of the induced distances from \f$x_{i}\f$ 
  to  \f$\mu_{j}\f$. You can also use another distance that inherits from dist::Dist

  \param aicrispmatrix_partition a  mat::CrispMatrix partition of instances in clusters
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances

*/
template < typename INPUT_ITERATOR,
           typename T_METRIC, 
	   typename T_FEATURE,
	   typename T_BITSIZE
	   >
T_METRIC
j1
(mat::BitMatrix<T_BITSIZE>      &aicrispmatrix_partition, 
 mat::MatrixRow<T_FEATURE>      &aimatrixt_centroids,
 INPUT_ITERATOR                 aiiterator_instfirst,
 const INPUT_ITERATOR           aiiterator_instlast,
 dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 )
{
  T_METRIC loT_J1;

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::j1";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
  	      << "\t(input  mat::BitMatrix<T_BITSIZE>: aicrispmatrix_partition[" << &aicrispmatrix_partition << "]\n"
	      << "\t(input  mat::MatrixRow<T_FEATURE>: aimatrixt_centroids["
	      << &aimatrixt_centroids << "]\n"
	      << "\t input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES*/

  loT_J1 = T_METRIC(0);
  for (uintidx lui_i = 0;
       aiiterator_instfirst != aiiterator_instlast;
       ++aiiterator_instfirst, lui_i++)
    {
      data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
      for (uintidx lui_j = 0; lui_j < aimatrixt_centroids.getNumRows(); lui_j++)  {
	if ( aicrispmatrix_partition(lui_j,lui_i) ) {
	  T_FEATURE* lmatrowt_centroid = aimatrixt_centroids.getRow(lui_j);
	  loT_J1 += 
	    aifunc2p_dist
	    (lmatrowt_centroid,
	     linst_inter->getFeatures(),
	     data::Instance<T_FEATURE>::getNumDimensions()
	     );
	}
      }
    }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") J1 =  " << loT_J1
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
 
  return loT_J1;
}


/*! \fn T_METRIC jm(mat::MatrixRow<T_METRIC> &aimatrixt_u, mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, T_METRIC aif_m, dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist )
  \brief Least-squared errors functional \f$(J_m)\f$ \cite Bezdek:ClusterAnalysis:FCM:1974
  \details For fuzzy c-partitions in \f$X\f$, get Least-squared errors functional
  \f[
  J_m(U,\mu)=\sum_{i=1}^{n}\sum_{j=1}^{k}u_{ji}^mD_{Ind}(x_{i}-\mu_{j}),
  \f]
  Where
  \$fU\in M_{fc}\f$ fuzzy c-partition of \$fX\f$;
  \$f\mu=[\mu_{1},\mu_{2},...,\mu_{k}]\f$ centroids,
  \$fm\f$ weighting exponent; \$f1 \leq m < \infty\f$
  \$fD_{Ind}(x_{k},\mu_{i})\f$  is one of the induced distances from \$fx_{i}\f$ to  \$f\mu_{j}\f$.
  \param &aimatrixt_u a fuzzy c-partition mat::MatrixRow
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aif_m a real number with the weighting exponent
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_METRIC,
	   typename T_FEATURE
	   >
T_METRIC 
jm
(mat::MatrixRow<T_METRIC>         &aimatrixt_u,
 mat::MatrixRow<T_FEATURE>        &aimatrixt_centroids,
 INPUT_ITERATOR                   aiiterator_instfirst,
 const INPUT_ITERATOR             aiiterator_instlast,
 T_METRIC                         aif_m,
 dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist
 )
{
  T_METRIC  lort_Jm;
  T_METRIC  lrt_dist;
  
  lort_Jm = T_METRIC(0);
  for (uintidx lui_i = 0;
       aiiterator_instfirst != aiiterator_instlast;
       ++aiiterator_instfirst, lui_i++)
    {
      data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
   
      for (uintidx lui_j = 0; lui_j < aimatrixt_centroids.getNumRows(); lui_j++)  {
	T_FEATURE* lmatrowt_centroid = aimatrixt_centroids.getRow(lui_j);
	lrt_dist = 
	  aifunc2p_dist
	  (lmatrowt_centroid,
	   linst_inter->getFeatures(),
	   data::Instance<T_FEATURE>::getNumDimensions() 
	   );
	lort_Jm  +=  std::pow(aimatrixt_u(lui_j, lui_i), aif_m) * lrt_dist;
      }
    }
  return lort_Jm;
}


/*! \fn T_METRIC SSEMedoid (const uintidx *aiarrayidxinst_medoids, const T_CLUSTERIDX aicidx_numClustersK, const mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity) 
  \brief SSEMedoid  \cite Lucasius:etal:GAclusteringMedoid:GCA:1993
  \details The cost function sum of Euclidean distances to the most representative instance
  \f[
  SED_{medoid}=\sum_{C_j}\sum_{x_i\in C_j}\Vert x_i - m_j \Vert
  \f]

  where \f$m_j\f$ represents the \a medoid or \a prototype of cluster \f$C_j\f$

  \param aiarrayidxinst_medoids an array of indexes of more representative instances
  \param aicidx_numClustersK Number of medoids or clusters
  \param aimatrixtriagt_dissimilarity a matrix of distances 
*/
template <typename T_METRIC,
          typename T_CLUSTERIDX
	  >
T_METRIC   
SSEMedoid
(const uintidx                     *aiarrayidxinst_medoids,
 const T_CLUSTERIDX                aicidx_numClustersK,
 const mat::MatrixTriang<T_METRIC> &aimatrixtriagt_dissimilarity
 ) 
{
  T_METRIC lort_SSE;
  T_METRIC lrt_distMinMedoidsInst;
   
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::SSEMedoid:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(input  uintidx: *aiarrayidxinst_medoids[" 
	      <<  aiarrayidxinst_medoids << "]\n"
	      << "\t input  T_CLUSTERIDX: aicidx_numClustersK = " 
	      <<  aicidx_numClustersK << "\n"
	      << "\t input  mat::MatrixTriang<T_METRIC>: aimatrixtriagt_dissimilarity[" 
	      << &aimatrixtriagt_dissimilarity << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  lort_SSE = T_METRIC(0);
  for (uintidx lui_i = 0; lui_i < aimatrixtriagt_dissimilarity.getNumRows(); lui_i++) {
    nearest::medoidsNN
      (lrt_distMinMedoidsInst,
       lui_i,
       aiarrayidxinst_medoids,
       aicidx_numClustersK,
       aimatrixtriagt_dissimilarity
       );
  
    lort_SSE += (T_METRIC) lrt_distMinMedoidsInst;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "um::SSEMedoid: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\noutput  T_METRIC:  lort_SSE = " << lort_SSE 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lort_SSE;  
}


/*! \fn T_METRIC indexI (mat::MatrixRow<T_METRIC> &aimatrixt_u, mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist, const T_METRIC airt_p)

  \brief  Validity index I or index I for a fuzzy c-partition  \cite Maulik:Bandyopadhyay:GAclustering:IndexI:2002 \cite Bandyopadhuay:Maulik:GAclustering:MO:2007.

  \details The \a validity \a index \a I or simply \a Index \a I. It is used as a metric to measure clustering performance. It was proposed as a measure to indicate the (goodness) validity of the solution in the cluster. It is defined as follows:

  \f[
  I(k)=\left(  {1 \over  k} \cdot {E_{1} \over E_{k} } \cdot D_{k}\right)^{p},
  \f]

  Where \f$k\f$ is the number of clusters 

  \f[
  E_{k}=\sum_{j=1}^{k}\sum_{i=1}^{n}u_{ji}\left\Vert x_{i}-\mu_{j}\right\Vert,
  \f]

  and

  \f[
  D_{k}=\max_{j,j'=1}^{k}\left\Vert \mu_{j'}-\mu_{j}\right\Vert
  \f]

  \f$n\f$ is the total number of objects. \f$U(X)=\left[u_{kj}\right]_{k \times n}\f$
  is a partition matrix of the objects and \f$\mu_{j}\f$  is the centroid es el centro
  \f$j^{th}\f$. The value of \f$k\f$ that maximizes \f$I(k)\f$ is considered
  the correct number of clusters.

  \param &aimatrixt_u a fuzzy c-partitions mat::MatrixRow
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param airt_p the power used to calculate the index, by default 2.0 
*/
template < typename INPUT_ITERATOR,
           typename T_METRIC,
	   typename T_FEATURE
	   >
T_METRIC 
indexI
(mat::MatrixRow<T_METRIC>         &aimatrixt_u,
 mat::MatrixRow<T_FEATURE>        &aimatrixt_centroids,
 INPUT_ITERATOR                   aiiterator_instfirst,
 const INPUT_ITERATOR             aiiterator_instlast,
 dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist,
 const T_METRIC                   airt_p = 2.0
 )
{
static T_METRIC lmetric_e1;
 
  static utils::RunOnce runOnce ([&]() {

      const uintidx  lui_numInstances = 
	uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
      
      T_FEATURE *larray_centroid1 =
	new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];
      
      decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	*larray_sumFeatureTmp =
	new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	[data::Instance<T_FEATURE>::getNumDimensions()];

      stats::sumFeactures
	(larray_sumFeatureTmp,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 T_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_centroid1,
	 lui_numInstances,
	 larray_sumFeatureTmp
	 );

      lmetric_e1 =
	e1
	(larray_centroid1,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 aifunc2p_dist
	 );

      delete [] larray_sumFeatureTmp;
      delete [] larray_centroid1;
    }
    );
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::IndexI";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_METRIC>"
	      << &aimatrixt_u << "]\n"
	      << aimatrixt_u << "\n"
	      << "\t input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << aimatrixt_centroids << "\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lometric_indexI = measuare_undefIndexI(T_METRIC);
  T_METRIC lmetric_ek = T_METRIC(0);
  T_METRIC lmetric_dk = T_METRIC(-1);

  if ( aimatrixt_centroids.getNumRows() > 0 ) {

       
    for (uintidx lui_i = 0;
	 aiiterator_instfirst != aiiterator_instlast;
	 ++aiiterator_instfirst, lui_i++)
      {
	data::Instance<T_FEATURE>* linst_inter =  
	  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
   
	for (uintidx lui_j = 0; lui_j < aimatrixt_centroids.getNumRows(); lui_j++)  {
	  T_FEATURE* lmatrowt_centroid = aimatrixt_centroids.getRow(lui_j);
	  T_METRIC  lrt_dist = 
	    aifunc2p_dist
	    (lmatrowt_centroid,
	     linst_inter->getFeatures(),
	     data::Instance<T_FEATURE>::getNumDimensions() 
	     );
	  	 
	  lmetric_ek  += aimatrixt_u(lui_j, lui_i) * lrt_dist;
	  
	}
      } //End for

      lmetric_dk =
      maxDistCjCjp
      (aimatrixt_centroids,	    
       aifunc2p_dist
       );
      
    lometric_indexI = (( lmetric_e1 / lmetric_ek ) *  lmetric_dk )
      / T_METRIC(aimatrixt_centroids.getNumRows());
  
    lometric_indexI = std::pow(lometric_indexI,airt_p); 

  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lmetric_e1 " << lmetric_e1
              << " lmetric_ek " << lmetric_ek
              << " lmetric_dk " << lmetric_dk
	      << " lometric_indexI = " << lometric_indexI 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_indexI;
  
}


/*! \fn T_METRIC indexI (const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX>    &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist, const T_METRIC airt_p = 2.0)   
  \brief Validity index I for a crisp partition \cite Maulik:Bandyopadhyay:GAclustering:IndexI:2002 \cite Bandyopadhuay:Maulik:GAclustering:MO:2007.
  \details The \a validity \a index \a I or simply \a Index \a I. It is used as a metric to measure clustering performance. It was proposed as a measure to indicate the (goodness) validity of the solution in the cluster. It is defined as follows:

  \f[
  I(k)=\left(  {1 \over  k} \cdot {E_{1} \over E_{k} } \cdot D_{k}\right)^{p},
  \f]

  Where \f$k\f$ is the number of clusters 

  \f[
  E_{k}=\sum_{j=1}^{k}\sum_{i=1}^{n}u_{ji}\left\Vert x_{i}-\mu_{j}\right\Vert,
  \f]

  and

  \f[
  D_{k}=\max_{j,j'=1}^{k}\left\Vert \mu_{j'}-\mu_{j}\right\Vert
  \f]

  \f$n\f$ is the total number of objects. \f$U(X)=\left[u_{kj}\right]_{k \times n}\f$
  is a partition matrix of the objects and \f$\mu_{j}\f$  is the centroid es el centro
  \f$j^{th}\f$. The value of \f$k\f$ that maximizes \f$I(k)\f$ is considered
  the correct number of clusters.

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters 
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a crisp partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances 
  \param airt_p the power used to calculate the index, by default 2.0 
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
indexI
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist,
 const T_METRIC                        airt_p = 2.0
 )
{

  static T_METRIC lmetric_e1;
 
  static utils::RunOnce runOnce ([&]() {

      const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
      
      T_FEATURE *larray_centroid1 =
	new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];
      
      decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	*larray_sumFeatureTmp =
	new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
	[data::Instance<T_FEATURE>::getNumDimensions()];

      stats::sumFeactures
	(larray_sumFeatureTmp,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 T_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_centroid1,
	 lui_numInstances,
	 larray_sumFeatureTmp
	 );

      lmetric_e1 =
	e1
	(larray_centroid1,
	 aiiterator_instfirst,
	 aiiterator_instlast,
	 aifunc2p_dist
	 );

      delete [] larray_sumFeatureTmp;
      delete [] larray_centroid1;
    }
    );
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::IndexI";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lometric_indexI = measuare_undefIndexI(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 0 ) {
    std::pair<T_METRIC,bool> lpair_ek = 
      SSE
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_dist
       );

    T_METRIC
      lmetric_dk =
      maxDistCjCjp
      (aimatrixt_centroids,	    
       aifunc2p_dist
       );

    lometric_indexI = (( lmetric_e1 / lpair_ek.first ) *  lmetric_dk )
      / T_METRIC(aimatrixt_centroids.getNumRows());
  
    lometric_indexI = std::pow(lometric_indexI,airt_p);
  }
  
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_indexI = " << lometric_indexI 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_indexI;
  
}


/*! \fn T_METRIC indexIreeval (mat::MatrixRow<T_METRIC> &aimatrixt_u, mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist, const T_METRIC airt_p)

  \brief  Validity index I or index I for a fuzzy c-partition  \cite Maulik:Bandyopadhyay:GAclustering:IndexI:2002 \cite Bandyopadhuay:Maulik:GAclustering:MO:2007. Use this function when the data set changes in the same execution, parameter \f$E_1\f$ is recalculated.

  \details The \a validity \a index \a I or simply \a Index \a I. It is used as a metric to measure clustering performance. It was proposed as a measure to indicate the (goodness) validity of the solution in the cluster. It is defined as follows:

  \f[
  I(k)=\left(  {1 \over  k} \cdot {E_{1} \over E_{k} } \cdot D_{k}\right)^{p},
  \f]

  Where \f$k\f$ is the number of clusters 

  \f[
  E_{k}=\sum_{j=1}^{k}\sum_{i=1}^{n}u_{ji}\left\Vert x_{i}-\mu_{j}\right\Vert,
  \f]

  and

  \f[
  D_{k}=\max_{j,j'=1}^{k}\left\Vert \mu_{j'}-\mu_{j}\right\Vert
  \f]

  \f$n\f$ is the total number of objects. \f$U(X)=\left[u_{kj}\right]_{k \times n}\f$
  is a partition matrix of the objects and \f$\mu_{j}\f$  is the centroid es el centro
  \f$j^{th}\f$. The value of \f$k\f$ that maximizes \f$I(k)\f$ is considered
  the correct number of clusters.

  \param &aimatrixt_u a fuzzy c-partitions mat::MatrixRow
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
  \param airt_p the power used to calculate the index, by default 2.0 
*/
template < typename INPUT_ITERATOR,
           typename T_METRIC,
	   typename T_FEATURE
	   >
T_METRIC 
indexIreeval
(mat::MatrixRow<T_METRIC>         &aimatrixt_u,
 mat::MatrixRow<T_FEATURE>        &aimatrixt_centroids,
 INPUT_ITERATOR                   aiiterator_instfirst,
 const INPUT_ITERATOR             aiiterator_instlast,
 dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist,
 const T_METRIC                   airt_p = 2.0
 )
{
  T_METRIC lmetric_e1;
 
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
      
  T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];
      
  decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    *larray_sumFeatureTmp =
    new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    [data::Instance<T_FEATURE>::getNumDimensions()];

  stats::sumFeactures
    (larray_sumFeatureTmp,
     aiiterator_instfirst,
     aiiterator_instlast,
     T_FEATURE(0)
     );
  
  stats::meanVector
    (larray_centroid1,
     lui_numInstances,
     larray_sumFeatureTmp
     );

  lmetric_e1 =
    e1
    (larray_centroid1,
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  delete [] larray_sumFeatureTmp;
  delete [] larray_centroid1;
 
  
  T_METRIC lometric_indexI = measuare_undefIndexI(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 0 ) {

    T_METRIC lmetric_ek = T_METRIC(0);
    
    for (uintidx lui_i = 0;
	 aiiterator_instfirst != aiiterator_instlast;
	 ++aiiterator_instfirst, lui_i++)
      {
	data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
   
	for (uintidx lui_j = 0; lui_j < aimatrixt_centroids.getNumRows(); lui_j++)  {
	  T_FEATURE* lmatrowt_centroid = aimatrixt_centroids.getRow(lui_j);
	  T_METRIC  lrt_dist = 
	    aifunc2p_dist
	    (lmatrowt_centroid,
	     linst_inter->getFeatures(),
	     data::Instance<T_FEATURE>::getNumDimensions() 
	     );
	  T_METRIC lrt_uji = aimatrixt_u(lui_j, lui_i);
	  
	  lmetric_ek  += lrt_uji * lrt_dist;
	  
	}
      } //End for

    T_METRIC
      lmetric_dk =
      maxDistCjCjp
      (aimatrixt_centroids,	    
       aifunc2p_dist
       );

    lometric_indexI = (( lmetric_e1 / lmetric_ek ) *  lmetric_dk )
      / T_METRIC(aimatrixt_centroids.getNumRows());
  
    lometric_indexI = std::pow(lometric_indexI,airt_p); 

  }

  return lometric_indexI;
  
}

  
/*! \fn T_METRIC indexIreeval (const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX>    &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist, const T_METRIC airt_p = 2.0)   
  \brief Validity index I for a crisp partition \cite Maulik:Bandyopadhyay:GAclustering:IndexI:2002 \cite Bandyopadhuay:Maulik:GAclustering:MO:2007. Use this function when the data set changes in the same execution, parameter \f$E_1\f$ is recalculated.

  \details The \a validity \a index \a I or simply \a Index \a I. It is used as a metric to measure clustering performance. It was proposed as a measure to indicate the (goodness) validity of the solution in the cluster. It is defined as follows:

  \f[
  I(k)=\left(  {1 \over  k} \cdot {E_{1} \over E_{k} } \cdot D_{k}\right)^{p},
  \f]

  Where \f$k\f$ is the number of clusters 

  \f[
  E_{k}=\sum_{j=1}^{k}\sum_{i=1}^{n}u_{ji}\left\Vert x_{i}-\mu_{j}\right\Vert,
  \f]

  and

  \f[
  D_{k}=\max_{j,j'=1}^{k}\left\Vert \mu_{j'}-\mu_{j}\right\Vert
  \f]

  \f$n\f$ is the total number of objects. \f$U(X)=\left[u_{kj}\right]_{k \times n}\f$
  is a partition matrix of the objects and \f$\mu_{j}\f$  is the centroid es el centro
  \f$j^{th}\f$. The value of \f$k\f$ that maximizes \f$I(k)\f$ is considered
  the correct number of clusters.

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters 
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a crisp partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances 
  \param airt_p the power used to calculate the index, by default 2.0 
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
indexIreeval
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>    &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist,
 const T_METRIC                        airt_p = 2.0
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::IndexIreeval";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC lmetric_e1;
 
  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
      
  T_FEATURE *larray_centroid1 =
    new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];
      
  decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    *larray_sumFeatureTmp =
    new decltype(utils::InstanceDataType().sum(data::Instance<T_FEATURE>::type()))
    [data::Instance<T_FEATURE>::getNumDimensions()];

  stats::sumFeactures
    (larray_sumFeatureTmp,
     aiiterator_instfirst,
     aiiterator_instlast,
     T_FEATURE(0)
     );
  
  stats::meanVector
    (larray_centroid1,
     lui_numInstances,
     larray_sumFeatureTmp
     );

  lmetric_e1 =
    e1
    (larray_centroid1,
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist
     );

  delete [] larray_sumFeatureTmp;
  delete [] larray_centroid1;
   


  T_METRIC lometric_indexI = measuare_undefIndexI(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 0 ) {
    std::pair<T_METRIC,bool> lpair_ek = 
      SSE
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_dist
       );

    T_METRIC
      lmetric_dk =
      maxDistCjCjp
      (aimatrixt_centroids,	    
       aifunc2p_dist
       );

    lometric_indexI = (( lmetric_e1 / lpair_ek.first ) *  lmetric_dk )
      / T_METRIC(aimatrixt_centroids.getNumRows());
  
    lometric_indexI = std::pow(lometric_indexI,airt_p);
  }
  
	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_indexI = " << lometric_indexI 
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_indexI;
  
}
  
  
/*! \fn T_METRIC xb(mat::MatrixRow<T_METRIC> &aimatrixt_u, mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, T_METRIC aif_m, dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist )

  \brief  The index of Xie-Beni (XB) 
  \cite Xie:Beni:Metricclustering:XieBeniIndex:1991
  \details Get the index of Xie-Beni for a fuzzy c-partition

  \f[
  XB =  {\sigma  \over n \cdot (d_{min})^2} 
  \f]

  In detail
 
  \f[
  \sigma = \sum_{j=1}^k\sum_{x_i}^n u_{ji}^2 \Vert x_i - \mu_j \Vert^2  \hbox{,}
  \f]

  \f[
  d_{min}  =\min_{j,j'=1,j\neq j'}^{k}\left\Vert \mu_{j}-\mu_{j'}\right \Vert 
  \f]

  \f[
  XB = {\sum_{j=1}^k\sum_{x_i}^n u_{ji}^2 \Vert x_i - \mu_j \Vert^2 \over n \cdot (d_{min})^2}
  \f]

  \param &aimatrixt_u a fuzzy c-partitions mat::MatrixRow
  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_METRIC,
	   typename T_FEATURE
	   >
T_METRIC 
xb
(mat::MatrixRow<T_METRIC>         &aimatrixt_u,
 mat::MatrixRow<T_FEATURE>        &aimatrixt_centroids,
 INPUT_ITERATOR                   aiiterator_instfirst,
 const INPUT_ITERATOR             aiiterator_instlast,
 dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dist
 )
{

  T_METRIC  lometric_xb = measuare_undefXieBeniIndex(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 1 ) {

    const T_METRIC  lrt_numInstances =
      T_METRIC(std::distance(aiiterator_instfirst,aiiterator_instlast));
     
    T_METRIC  lrt_sigma =  T_METRIC(0);
     
    for (uintidx lui_i = 0;
	 aiiterator_instfirst != aiiterator_instlast;
	 ++aiiterator_instfirst, lui_i++)
      {
	data::Instance<T_FEATURE>* linst_inter =  (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
   
	for (uintidx lui_j = 0; lui_j < aimatrixt_centroids.getNumRows(); lui_j++)  {
	  T_FEATURE* lmatrowt_centroid = aimatrixt_centroids.getRow(lui_j);
	  T_METRIC  lrt_dist = 
	    aifunc2p_dist
	    (lmatrowt_centroid,
	     linst_inter->getFeatures(),
	     data::Instance<T_FEATURE>::getNumDimensions() 
	     );
	  T_METRIC lrt_uji = aimatrixt_u(lui_j, lui_i);
	  
	  lrt_sigma  +=  lrt_uji * lrt_uji * lrt_dist;
	  
	}
      } //End for

    T_METRIC
      lmetrict_dmin =
      minDistCjCjp
      (aimatrixt_centroids,	    
       aifunc2p_dist
       );

    if ( lmetrict_dmin > 0 ) 
      lometric_xb =  lrt_sigma / ( lrt_numInstances *  lmetrict_dmin );
      
  }
  
  return lometric_xb;
  
}
  
/*! \fn T_METRIC xb(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const partition::Partition<T_CLUSTERIDX> &aipartition_clusters, const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist)
  \brief  The index of Xie-Beni (XB) 
  \cite Xie:Beni:Metricclustering:XieBeniIndex:1991
  \details Get the index of Xie-Beni extended to a crisp partition

  \param aimatrixt_centroids a mat::MatrixRow with centroids clusters
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
  \param aipartition_clusters a partition of instances in clusters
  \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE,
	   typename T_CLUSTERIDX,
	   typename T_METRIC  
	   >
T_METRIC
xb
(const mat::MatrixRow<T_FEATURE>      &aimatrixt_centroids,
 INPUT_ITERATOR                       aiiterator_instfirst,
 const INPUT_ITERATOR                 aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>   &aipartition_clusters,
 const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "um::XB";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n\t(input mat::MatrixRow<T_FEATURE>& aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << "\t input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_METRIC  lometric_xb = measuare_undefXieBeniIndex(T_METRIC);

  if ( aimatrixt_centroids.getNumRows() > 1 ) {

    const T_METRIC  lrt_numInstances =
      T_METRIC(std::distance(aiiterator_instfirst,aiiterator_instlast));
    
    std::pair<T_METRIC,bool> lpair_sigma = 
      SSE
      (aimatrixt_centroids,
       aiiterator_instfirst,
       aiiterator_instlast,
       aipartition_clusters,
       aifunc2p_dist
       );

    T_METRIC
      lmetrict_dmin =
      minDistCjCjp
      (aimatrixt_centroids,	    
       aifunc2p_dist
       );

    if ( lmetrict_dmin > 0 ) 
      lometric_xb =  lpair_sigma.first / ( lrt_numInstances *  lmetrict_dmin );
  
  }
  	  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lometric_xb = " << lometric_xb
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lometric_xb;
    
}


/*! \fn T_METRIC  entropy(const mat::MatrixRow<T_METRIC>  &aimatrixt_u)

  \brief Get entropy \f$ H_c(U) \f$ \cite Bezdek:etal:ClusterAnalysis:FCM:1984
  \details 

  \f[
  H_c(U) =  -\sum_{i=1}^n \sum_{j=1}^k (u_{ji} \log_a( u_{ji} )) / n.
  \f]

  \param aimatrixt_u a fuzzy c-partitions mat::MatrixRow

*/
template < typename T_METRIC>
T_METRIC  entropy(const mat::MatrixRow<T_METRIC>  &aimatrixt_u)
{

  T_METRIC  lometric_entropy = measuare_undefEntropy(T_METRIC);

  uintidx lui_numElems  = aimatrixt_u.getNumElems();  
  
  if ( lui_numElems > 0 ) {

    lometric_entropy = 
      std::accumulate
      (aimatrixt_u.toArray(),
       aimatrixt_u.toArray()+lui_numElems,
       T_METRIC(0.0),
       [](T_METRIC airt_sumPart, T_METRIC airt_uij)
       {
	 return (airt_uij) * std::log(airt_uij) + airt_sumPart;
       }
       );
    lometric_entropy /= -T_METRIC(aimatrixt_u.getNumColumns());
  }
    
  return lometric_entropy;
  
}

/*! \fn T_METRIC  partitioncoefficient(const mat::MatrixRow<T_METRIC>  &aimatrixt_u)

  \brief Get partition coefficient \f$ F_c(U) \f$ \cite Bezdek:etal:ClusterAnalysis:FCM:1984
  \details 

  \f[
  F_c(U) =  \sum_{i=1}^n \sum_{j=1}^k (u_{ji})^2 / n.
  \f]

  \param aimatrixt_u a fuzzy c-partitions mat::MatrixRow

*/
template < typename T_METRIC>
T_METRIC  partitioncoefficient(const mat::MatrixRow<T_METRIC>  &aimatrixt_u)
{

  T_METRIC  lometric_partitioncoefficient = measuare_undefPartitionCoefficient(T_METRIC);

  uintidx lui_numElems  = aimatrixt_u.getNumElems();  
  
  if ( lui_numElems > 0 ) {

    lometric_partitioncoefficient  = 
      std::accumulate
      (aimatrixt_u.toArray(),
       aimatrixt_u.toArray()+lui_numElems,
       T_METRIC(0.0),
       [](T_METRIC airt_sumPart, T_METRIC airt_uij)
       {
	 return (airt_uij) * airt_uij + airt_sumPart;
       }
       );
    lometric_partitioncoefficient /= T_METRIC(aimatrixt_u.getNumColumns());
  }
    
  return lometric_partitioncoefficient;
  
}
  
  
} /*END namespace um Unsupervised measures*/

#endif /*__UNSUPERVISED_MEASURE_HPP*/
