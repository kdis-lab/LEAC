/*! \file partition_disjsets.hpp
 *
 * \brief Partition based on structure of disjoint joint data
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef PARTITION_DISJSETS_HPP
#define PARTITION_DISJSETS_HPP

#include <map>
#include "partition.hpp"
#include "disjsets.hpp"
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


/*! \class PartitionDisjSets
  \brief Partition of instances in base DisjSets data structure
*/
template < class T_CLUSTERIDX >
class  PartitionDisjSets
  : public Partition<T_CLUSTERIDX>  
  , public ds::DisjSets 
{
public:

  PartitionDisjSets()
    : Partition<T_CLUSTERIDX>()
    , ds::DisjSets()
  {
  }
 
  PartitionDisjSets(ds::DisjSets   &&aids_disjset)
    : Partition<T_CLUSTERIDX>()
    , ds::DisjSets(aids_disjset)
  {
    T_CLUSTERIDX lcidx_labelClustersK = 0;
    
    for (uintidx luintidx_i = 0; luintidx_i <  this->_vectorst_parent.size(); ++luintidx_i) {
      if ( luintidx_i == this->_vectorst_parent[luintidx_i] ) {
	_mapset_clusterk.insert
	  ( std::pair<uintidx,T_CLUSTERIDX>(luintidx_i,lcidx_labelClustersK++) );
      }
    }
    
  }

  PartitionDisjSets
  (ds::DisjSets               &&aids_disjset,
   const std::vector<uintidx> &aivectoruintidx_idxToMap
   )
    : Partition<T_CLUSTERIDX>()
    , DisjSets(aids_disjset)
  {
    
    T_CLUSTERIDX lcidx_labelClustersK = 0;

    for (uintidx luintidx_i = 0; luintidx_i <  aivectoruintidx_idxToMap.size(); luintidx_i++) {
      uintidx luintidx_mapIdx = constfind(aivectoruintidx_idxToMap[luintidx_i]);

      _mapset_clusterk.insert                  
	( std::pair<uintidx,T_CLUSTERIDX>(luintidx_mapIdx,lcidx_labelClustersK) );
      ++lcidx_labelClustersK;
    }
    
  }

  PartitionDisjSets
  (const PartitionDisjSets<T_CLUSTERIDX> &aimembclassdisjset_b)
    : Partition<T_CLUSTERIDX>(aimembclassdisjset_b)
    , DisjSets(aimembclassdisjset_b)
    , _mapset_clusterk(aimembclassdisjset_b._mapset_clusterk)
  {
  }
  
  PartitionDisjSets
  (PartitionDisjSets<T_CLUSTERIDX> &&aimembclassdisjset_b)
    : Partition<T_CLUSTERIDX>(aimembclassdisjset_b)
    , DisjSets(aimembclassdisjset_b)
    , _mapset_clusterk(std::move(aimembclassdisjset_b._mapset_clusterk))
  {
  }

  ~PartitionDisjSets() {}

  void begin()
  {
    _ui_iteratornext = 0;
  }

  inline const T_CLUSTERIDX next()
  {
    const uintidx luintidx_clusterK = this->constfind(_ui_iteratornext);
   
    
    typename std::map<uintidx,T_CLUSTERIDX>::const_iterator litmapset_clusterk =
      _mapset_clusterk.find(luintidx_clusterK);

    T_CLUSTERIDX locidx_clusterK =
      ( litmapset_clusterk != _mapset_clusterk.end() )?
      litmapset_clusterk->second:NEARESTCENTROID_UNKNOWN;

    _ui_iteratornext++;
      
    return locidx_clusterK;
    
  }

  inline bool end() const
  {
    return ( _ui_iteratornext < this->_vectorst_parent.size() );
  }
  
  PartitionDisjSets<T_CLUSTERIDX>&
  operator=(const PartitionDisjSets<T_CLUSTERIDX> &aimembclassdisjset_b)
  {
    if( this !=  &aimembclassdisjset_b ) {
      DisjSets::operator=(aimembclassdisjset_b);
      this->_mapset_clusterk = std::move(aimembclassdisjset_b._mapset_clusterk);
    }

    return *this;
  }

  PartitionDisjSets<T_CLUSTERIDX>&
  operator=(PartitionDisjSets<T_CLUSTERIDX> &&aimembclassdisjset_b)
  {
    if( this !=  &aimembclassdisjset_b ) {
      DisjSets::operator=(aimembclassdisjset_b);
      this->_mapset_clusterk = std::move(aimembclassdisjset_b._mapset_clusterk);
    }

    return *this;
  }

  const T_CLUSTERIDX getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {
    const uintidx luintidx_clusterK = this->constfind(aiuintidx_instanceIdx);
   
    typename std::map<uintidx,T_CLUSTERIDX>::const_iterator litmapset_clusterk =
      _mapset_clusterk.find(luintidx_clusterK);

    T_CLUSTERIDX locidx_clusterK =
      ( litmapset_clusterk != _mapset_clusterk.end() )?
      litmapset_clusterk->second:NEARESTCENTROID_UNKNOWN;
   
    return locidx_clusterK;
   
  }

  inline const uintidx getNumInstances() const
  {
    return uintidx(this->_vectorst_parent.size());
  }
  
  const T_CLUSTERIDX getNumCluster() const
  {
    return (T_CLUSTERIDX) this->_m;
  }

  void  print
  (std::ostream &os=std::cout,
   const char   *aipc_label     = "",
   const char   aic_delimCoef   =','
   )
  {

#if defined(__VERBOSE_YES)
    os << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << this << ']'
#else
    os << aipc_label
#endif
       << ":number of set," << this->_m << '>';
    
    if ( this->_vectorst_parent.size() > 0 ) {
      uintidx  luintidx_last = this->_vectorst_parent.size() - 1;   
      for (uintidx luintidx_i = 0; luintidx_i < luintidx_last ; ++luintidx_i ) {
	os << getClusterIdx(luintidx_i) << aic_delimCoef;
      }
      os << getClusterIdx(luintidx_last);
   
    }
    os << std::endl;
  }
   
protected:

  std::map<uintidx,T_CLUSTERIDX> _mapset_clusterk;
  uintidx                        _ui_iteratornext;
 
};

} /* END namespace partition*/

#endif  /* PARTITION_DISJSETS_HPP */
