/*! \file partition_label.hpp
 *
 * \brief Partition based on the label of the instances
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MEMBERSHIP_PARTITION_LABEL_HPP
#define MEMBERSHIP_PARTITION_LABEL_HPP

#include <cassert>
#include "partition.hpp"
#include "container_out.hpp"
#include "common.hpp"


/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace partition {


/*! \class PartitionLabel
  \brief Partition of instances with label string assigned previous
*/
template < class T_CLUSTERIDX >
class  PartitionLabel:
  public 
Partition<T_CLUSTERIDX>  {
public:

  PartitionLabel
  (T_CLUSTERIDX       *aiarrayT_idxMemberShip,
   const uintidx      aiui_numInstances,
   const T_CLUSTERIDX aicidx_numClustersK
   )
    : Partition<T_CLUSTERIDX>()
    , _arrayit_idxMemberShip(aiarrayT_idxMemberShip)
    , _ui_numInstances(aiui_numInstances)
    , _cidx_numClustersK(aicidx_numClustersK)
    //, _ptcidx_iteratornext(aiarrayT_idxMemberShip)
    , _ptcidx_iteratorlast(aiarrayT_idxMemberShip + aiui_numInstances -1)    
  {}

  ~PartitionLabel() {} 

  void begin() 
  {

    _ptcidx_iteratornext = _arrayit_idxMemberShip;
    --_ptcidx_iteratornext;

  }

  inline const T_CLUSTERIDX next()
  {
    //T_CLUSTERIDX locidx_cluster = *_ptcidx_iteratornext;
   
    ++_ptcidx_iteratornext;
    
    return *_ptcidx_iteratornext; //locidx_cluster; 
  }

  inline bool end() const
  {
    return ( _ptcidx_iteratornext != _ptcidx_iteratorlast);
  }
  
  virtual const T_CLUSTERIDX 
  getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {
    assert(0 <= aiuintidx_instanceIdx && aiuintidx_instanceIdx < _ui_numInstances);
    return _arrayit_idxMemberShip[aiuintidx_instanceIdx];
  }

  inline const uintidx getNumInstances() const
  {
    return _ui_numInstances;
  }
  
  const T_CLUSTERIDX getNumCluster() const
  {
    return  _cidx_numClustersK; 
  }

  void setNumCluster(const T_CLUSTERIDX aicidx_numClustersK)
  {
    _cidx_numClustersK = aicidx_numClustersK; 
  }

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label     = "",
   const char  aic_delimCoef  = ','
   ) //const
  {
    std::ostringstream lostrstream_labelMember;

    lostrstream_labelMember << aipc_label << aic_delimCoef << "_outK" << aic_delimCoef << this->getNumCluster();

#if defined(__VERBOSE_YES)
    lostrstream_labelMember << aic_delimCoef
			    << "PartitionLabel[" << geverboseui_idproc << ':' << this << ']';
#else
    lostrstream_labelMember << aic_delimCoef
			    << "PartitionLabel[" << this << ']';
#endif

    inout::containerprint
      (_arrayit_idxMemberShip,
       _arrayit_idxMemberShip + _ui_numInstances,
       os,
       lostrstream_labelMember.str().c_str(),
       aic_delimCoef
       );
    os << std::endl;
  }

protected:
  
  T_CLUSTERIDX*       _arrayit_idxMemberShip;
  const uintidx       _ui_numInstances;
  T_CLUSTERIDX        _cidx_numClustersK;
  T_CLUSTERIDX*       _ptcidx_iteratornext;
  const T_CLUSTERIDX* _ptcidx_iteratorlast;
  
}; /*PartitionLabel*/

  
} /* END namespace partition*/

#endif  /* MEMBERSHIP_PARTITION_LABEL_HPP */
