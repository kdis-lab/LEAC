/*! \file partition_labelvector.hpp
 *
 * \brief Partition based on the label vector
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __PARTITION_LABEL_VECTOR_HPP
#define __PARTITION_LABEL_VECTOR_HPP

#include <vector>
#include "partition.hpp"
#include "common.hpp"
#include "container_out.hpp"

/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace partition {


/*! \class PartitionLabelVector
  \brief Partition of instances with label vector assigned previous
*/
template <class T_CLUSTERIDX>
class  PartitionLabelVector:
  public Partition<T_CLUSTERIDX>   
{
public:

  PartitionLabelVector
  (std::vector<T_CLUSTERIDX> &avectorcidx_idxMemberShip,
   T_CLUSTERIDX              aicidx_numClustersK
   )
    : Partition<T_CLUSTERIDX>()
    , _vectorcidx_idxMemberShip(avectorcidx_idxMemberShip)
    , _cidx_numClustersK(aicidx_numClustersK)
    
  {}

  ~PartitionLabelVector() {} 

  void begin() 
  {
    _itvec_iteratornext = _vectorcidx_idxMemberShip.begin();
    --_itvec_iteratornext;
  }

   inline const T_CLUSTERIDX next()
  {
    ++_itvec_iteratornext;
    return *_itvec_iteratornext;
    
  }

  inline bool end() const
  {
    return ( _itvec_iteratornext != _vectorcidx_idxMemberShip.end() );
  }
  
  virtual const T_CLUSTERIDX 
  getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {
    return _vectorcidx_idxMemberShip.at(aiuintidx_instanceIdx);
  }

  const T_CLUSTERIDX getNumCluster() const
  {
    return  _cidx_numClustersK; 
  }

  const uintidx getNumInstances() const
  {
    return uintidx(_vectorcidx_idxMemberShip.size());
  }
    
 virtual void print 
  (std::ostream &os=std::cout,
   const char* aipc_label     = "",
   const char  aic_delimCoef  = ','
   )
  {
    std::ostringstream lostrstream_labelMember;

    lostrstream_labelMember << aipc_label << aic_delimCoef << "_outK" << aic_delimCoef << this->getNumCluster();
    
#if defined(__VERBOSE_YES)
    lostrstream_labelMember << aic_delimCoef
			    << "PartitionLabelVector[" << geverboseui_idproc << ':' << this << ']';
#else
    lostrstream_labelMember << aic_delimCoef << "PartitionLabelVector[" << this << ']';
#endif

    inout::containerprint
      (_vectorcidx_idxMemberShip.begin(),
       _vectorcidx_idxMemberShip.end(),
       os,
       lostrstream_labelMember.str().c_str(),
       aic_delimCoef
       );
    os << std::endl;
    
  }
  
protected:
  
  std::vector<T_CLUSTERIDX> &_vectorcidx_idxMemberShip;
  T_CLUSTERIDX                    _cidx_numClustersK;

  typename std::vector<T_CLUSTERIDX>::iterator _itvec_iteratornext;
  
}; /*PartitionLabelVector*/

  
} /* END namespace partition*/


#endif  /* __PARTITION_LABEL_VECTOR_HPP */
