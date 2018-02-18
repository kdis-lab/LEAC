/*! \file partition_centroids.hpp
 *
 * \brief Partition based on centroids
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MEMBERSHIP_PARTITION_CENTROIDS_HPP
#define MEMBERSHIP_PARTITION_CENTROIDS_HPP

#include <cassert>
#include <iterator>     // std::distance
#include <cassert>

#include "partition.hpp"
#include "matrix.hpp"
#include "nearestinstance_operator.hpp"

#include "verbose_global.hpp"

/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace partition {

/*! \class PartitionCentroids
  \brief Partition of instances in base to following equation:
  \f[
  x_i \in C_j \leftrightarrow  \| x_i - \mu_j \|  \begin{array}{c}min\\  k \end{array}
  \| x_i - \mu_k \|,\; j=1,2,..k,
  \f]
  where \f$m_j\f$, represents the medoid of cluster \f$C_j\f$
*/  
template < typename T_FEATURE,
           typename T_CLUSTERIDX,
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
class  PartitionCentroids:
  public Partition<T_CLUSTERIDX>  
{
public:

  PartitionCentroids()
    : Partition<T_CLUSTERIDX>()
    , _ptmatrixt_centroids(NULL)
    , _ptfunc2p_dist(NULL)
  {}
  
  PartitionCentroids
  (mat::MatrixRow<T_FEATURE>           &aimatrixt_centroids,
   const INPUT_ITERATOR                aiiterator_instfirst,
   const INPUT_ITERATOR                aiiterator_instlast,
   const T_CLUSTERIDX                  aicidx_numCluster,
   const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
   )
    : Partition<T_CLUSTERIDX>()
    , _ptmatrixt_centroids(&aimatrixt_centroids)
    , _iterator_instfirst(aiiterator_instfirst)
    , _iterator_instlast(aiiterator_instlast)
    , _iterator_next(aiiterator_instfirst)
    , _ptfunc2p_dist(&aifunc2p_dist)
  {
#ifdef __VERBOSE_YES
    assert(aicidx_numCluster ==  T_CLUSTERIDX(_ptmatrixt_centroids->getNumRows()));
#endif //__VERBOSE_YES
  }

  
  ~PartitionCentroids() {}

    //copy constructor
  PartitionCentroids(const PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>& aipartcentroids_clustering)
    : _ptmatrixt_centroids(aipartcentroids_clustering._ptmatrixt_centroids)
    , _iterator_instfirst(aipartcentroids_clustering._iterator_instfirst)
    , _iterator_instlast(aipartcentroids_clustering._iterator_instlast)
    , _iterator_next(aipartcentroids_clustering._iterator_next)
    , _ptfunc2p_dist(aipartcentroids_clustering._ptfunc2p_dist)
  {
  }
    

  //move constructor
  PartitionCentroids(PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>&& aipartcentroids_clustering)
    : _ptmatrixt_centroids(aipartcentroids_clustering._ptmatrixt_centroids)
    , _iterator_instfirst(aipartcentroids_clustering._iterator_instfirst)
    , _iterator_instlast(aipartcentroids_clustering._iterator_instlast)
    , _iterator_next(aipartcentroids_clustering._iterator_next)
    , _ptfunc2p_dist(aipartcentroids_clustering._ptfunc2p_dist)
  {
    
    aipartcentroids_clustering._ptmatrixt_centroids = NULL;
    aipartcentroids_clustering._ptfunc2p_dist = NULL;
  }

  PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>&
  operator=(const PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>& aipartcentroids_clustering)
  {
    if ( this != &aipartcentroids_clustering ) {
      _ptmatrixt_centroids = aipartcentroids_clustering._ptmatrixt_centroids;
      _iterator_instfirst = aipartcentroids_clustering._iterator_instfirst;
      _iterator_instlast = aipartcentroids_clustering._iterator_instlast;
      _iterator_next = aipartcentroids_clustering._iterator_next;
     _ptfunc2p_dist = aipartcentroids_clustering._ptfunc2p_dist;
    }

    return *this;
  }

  PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>&
  operator=(PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>&& aipartcentroids_clustering)
  {
    if ( this != &aipartcentroids_clustering ) {
      _ptmatrixt_centroids = aipartcentroids_clustering._ptmatrixt_centroids;
      _iterator_instfirst = aipartcentroids_clustering._iterator_instfirst;
      _iterator_instlast = aipartcentroids_clustering._iterator_instlast;
      _iterator_next = aipartcentroids_clustering._iterator_next;
     _ptfunc2p_dist = aipartcentroids_clustering._ptfunc2p_dist;

     aipartcentroids_clustering._ptmatrixt_centroids = NULL;
     aipartcentroids_clustering._iterator_instfirst = NULL;
     aipartcentroids_clustering._iterator_instlast = NULL;
     aipartcentroids_clustering._iterator_next = NULL;
     aipartcentroids_clustering._ptfunc2p_dist = NULL;
    }

    return (*this);
  }
  
  void begin()
  {
    _iterator_next = _iterator_instfirst;  
  }

  const T_CLUSTERIDX next()
  {     
    T_CLUSTERIDX  locidx_cluster;
    T_DIST               lT_distMinCentInst;

    data::Instance<T_FEATURE>* linst_inter =
      (data::Instance<T_FEATURE>*) *_iterator_next;
    locidx_cluster = 
      nearest::NN
      <T_CLUSTERIDX,
       T_FEATURE,
       T_DIST>
      (lT_distMinCentInst,
       *_ptmatrixt_centroids,
       linst_inter->getFeatures(),
       *_ptfunc2p_dist
       );
    //this->_ptvectorptinst_instances->at(aiuintidx_instanceIdx)->getFeatures(),
    ++_iterator_next;
    
    return locidx_cluster;
    
  }
 
  
  bool end() const
  {
    return ( _iterator_next != _iterator_instlast);
  }
  
  const T_CLUSTERIDX 
  getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {     
    T_CLUSTERIDX  locidx_cluster;
    T_DIST               lT_distMinCentInst;

    assert
      (0 <= aiuintidx_instanceIdx &&
       aiuintidx_instanceIdx < this->getNumInstances()
       );
    
    data::Instance<T_FEATURE>* linst_next = *std::next(_iterator_instfirst,aiuintidx_instanceIdx);
    
    locidx_cluster = 
      nearest::NN
      <T_CLUSTERIDX,
       T_FEATURE,
       T_DIST>
      (lT_distMinCentInst,
       *_ptmatrixt_centroids,
       linst_next->getFeatures(),
       *_ptfunc2p_dist
       );

    //this->_ptvectorptinst_instances->at(aiuintidx_instanceIdx)->getFeatures(),
    
    return locidx_cluster;
    
  }

  inline const uintidx getNumInstances() const
  {
    return uintidx(std::distance(_iterator_instfirst,_iterator_instlast));
  }
  
  const T_CLUSTERIDX getNumCluster() const
  {
    return T_CLUSTERIDX(_ptmatrixt_centroids->getNumRows());
  }

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label     = "",
   const char  aic_delimCoef  = ','
   )
  {
    os << aipc_label << aic_delimCoef << "_outK" << aic_delimCoef << this->getNumCluster();
    
#if defined(__VERBOSE_YES)
    os << aic_delimCoef
       << "PartitionCentroids[" << geverboseui_idproc << ':' << this << ']';
#else
    os << aic_delimCoef
       << "PartitionCentroids[" << this << ']';
#endif
    
    os << aic_delimCoef << "length" << aic_delimCoef << getNumInstances() << '>';
    this->begin();
    os << this->next();
    while ( this->end() ) { 
      os << aic_delimCoef <<  this->next();
    }

    os << std::endl;
  }

protected:
  
  mat::MatrixRow<T_FEATURE>           *_ptmatrixt_centroids;
  const INPUT_ITERATOR                _iterator_instfirst;
  const INPUT_ITERATOR                _iterator_instlast;
  INPUT_ITERATOR                      _iterator_next;
   
  const dist::Dist<T_DIST,T_FEATURE>  *_ptfunc2p_dist;
  
};
  

template < typename T_FEATURE,
           typename T_CLUSTERIDX,
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>
makePartition
(mat::MatrixRow<T_FEATURE>           &aimatrixt_centroids,
 const INPUT_ITERATOR                aiiterator_instfirst,
 const INPUT_ITERATOR                aiiterator_instlast,
 const T_CLUSTERIDX                  aicidx_numClusters,
 const dist::Dist<T_DIST,T_FEATURE>  &aifunc2p_dist
)
{

  return PartitionCentroids<T_FEATURE,T_CLUSTERIDX,T_DIST,INPUT_ITERATOR>(aimatrixt_centroids,aiiterator_instfirst,aiiterator_instlast,aicidx_numClusters,aifunc2p_dist);
  
}
    
} /* END namespace partition*/


#endif  /* MEMBERSHIP_PARTITION_CENTROIDS_HPP */
