/*! \file partition_linked_numinst.hpp
 *
 * \brief Linked partition data structure with number of instances
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef PARTITION_LINKED_NUMINST_HPP
#define PARTITION_LINKED_NUMINST_HPP

#include <algorithm>    // std::find
#include <vector>
#include <iterator>     // std::distance
#include "partition.hpp"
#include "nearestinstance_operator.hpp"  // NEARESTCENTROID_UNKNOWN
#include "common.hpp"
#include "partition_linked.hpp"
#include "container_out.hpp"

#include "verbose_global.hpp"

/*! \namespace ds
  \brief Data structure
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace ds {

/*! \class PartitionLinkedNumInst
  \brief A data structure for managing an instance partition with the number of instances per cluster
  \details 
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_INSTANCES_CLUSTER_K
	   > 
class PartitionLinkedNumInst : public PartitionLinked<T_CLUSTERIDX>
{
public:
  PartitionLinkedNumInst()
    : PartitionLinked<T_CLUSTERIDX>()
    , _vectorit_numInstClusterK(0)
  {}
  
  PartitionLinkedNumInst
  (const uintidx aiuintidx_numInstances, 
   const uintidx aiuintidx_numClusterK 
   ): PartitionLinked<T_CLUSTERIDX>
    (aiuintidx_numInstances, aiuintidx_numClusterK)
    , _vectorit_numInstClusterK(aiuintidx_numClusterK,0)
 {
   
  }

  //move constructor
  PartitionLinkedNumInst
  (PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &&aipartlink_b)
    : PartitionLinked<T_CLUSTERIDX>(aipartlink_b)
    , _vectorit_numInstClusterK(aipartlink_b._vectorit_numInstClusterK)
  {
  }
  
  //copy constructor
  PartitionLinkedNumInst
  (const PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &aipartlink_b)
    : PartitionLinked<T_CLUSTERIDX>(aipartlink_b)
    , _vectorit_numInstClusterK(aipartlink_b._vectorit_numInstClusterK)
  {
  }
  
  
  ~PartitionLinkedNumInst() 
  {
  }

  
  inline const T_INSTANCES_CLUSTER_K getNumInstClusterK(T_CLUSTERIDX aicidx_clusterK) const
  {
    return _vectorit_numInstClusterK[aicidx_clusterK];
  }

  inline const std::vector<T_INSTANCES_CLUSTER_K>& getVectorNumInstClusterK()  const
  {
    return _vectorit_numInstClusterK;
  }
  
  inline const bool haveNullCluster() const
  {
    return std::find
	  (_vectorit_numInstClusterK.begin(),
	   _vectorit_numInstClusterK.end(),
	   0) != _vectorit_numInstClusterK.end(); 
  }

  inline const T_CLUSTERIDX getNumNullCluster() const
  {
    return (T_CLUSTERIDX)
      std::count_if
       (_vectorit_numInstClusterK.begin(),
	_vectorit_numInstClusterK.end(),
	[] (const T_INSTANCES_CLUSTER_K aiit_num) {return aiit_num == 0;}
	);
  }
  
  PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>& 
  operator=(const PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>  &aipartlink_b)
  { 
    if ( this != &aipartlink_b ) {
      PartitionLinked<T_CLUSTERIDX>::operator=(aipartlink_b);
      _vectorit_numInstClusterK = aipartlink_b._vectorit_numInstClusterK;       
    }
    
    return *this;
  }

  PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>&
  operator=(PartitionLinkedNumInst<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> &&aipartlink_b)
  { 
    if ( this != &aipartlink_b ) {
      PartitionLinked<T_CLUSTERIDX>::operator=(aipartlink_b);
      _vectorit_numInstClusterK = aipartlink_b._vectorit_numInstClusterK;       
    }
    
    return *this;
  }
 
  inline void initialize()
  {
    PartitionLinked<T_CLUSTERIDX>::initialize();
    interfacesse::copya
      (_vectorit_numInstClusterK.data(),
       T_INSTANCES_CLUSTER_K(0),
       (uintidx) _vectorit_numInstClusterK.size()
       );
  }

 
  inline void addInstanceToCluster
  (const T_CLUSTERIDX   aiT_clusterIdx, 
   const uintidx               aiidxiT_instanceI,
   const T_INSTANCES_CLUSTER_K aiit_numInstFrec = 1
   )
  {
    PartitionLinked<T_CLUSTERIDX>::addInstanceToCluster(aiT_clusterIdx,aiidxiT_instanceI);
    _vectorit_numInstClusterK[aiT_clusterIdx] += aiit_numInstFrec;
  }
  

  void subInstanceFromCluster
  (T_CLUSTERIDX         aiidxmcT_clusterIdx,
   uintidx                     aiidxinstT_instanceI, //size__t aiidxinstT_instanceI)
   const T_INSTANCES_CLUSTER_K aiit_numInstFrec = 1
   ) 
  {
    PartitionLinked<T_CLUSTERIDX>::subInstanceFromCluster
      (aiidxmcT_clusterIdx,aiidxinstT_instanceI);
    if ( aiidxmcT_clusterIdx != NEARESTCENTROID_UNKNOWN ) {
      _vectorit_numInstClusterK[aiidxmcT_clusterIdx] -= aiit_numInstFrec;
    }
    
  }

  void changeMemberShip
  (T_CLUSTERIDX  aicidx_nearestClusterK,
   uintidx              aiidxinstT_instanceI,
   const T_INSTANCES_CLUSTER_K aiit_numInstFrec = 1
   )
  {
    this->subInstanceFromCluster(aiidxinstT_instanceI,aiit_numInstFrec);
    if ( aicidx_nearestClusterK != NEARESTCENTROID_UNKNOWN ) {
      this->addInstanceToCluster(aicidx_nearestClusterK,aiidxinstT_instanceI,aiit_numInstFrec);
    }
  }
  
  
  void  print
  (std::ostream &os=std::cout,
   const char*  aipc_label   = "",
   const char   aic_delimCoef =',',
   const char   aic_delimClusterK =';'
   )
  {
    PartitionLinked<T_CLUSTERIDX>
      ::print(os,aipc_label,aic_delimCoef,aic_delimClusterK);
    std::ostringstream lostrstream_labelNumInstClusterK;
    lostrstream_labelNumInstClusterK << "<INSTANCESCLUSTERK:" << aipc_label
				    << "_vectorit_numInstClusterK["
				    << &_vectorit_numInstClusterK << ']';
    inout::containerprint
      (_vectorit_numInstClusterK.begin(),
       _vectorit_numInstClusterK.end(),
       os,
       lostrstream_labelNumInstClusterK.str().c_str(),
       aic_delimCoef
       );
  }
 
protected:
 
  std::vector<T_INSTANCES_CLUSTER_K> _vectorit_numInstClusterK;
 
}; /*PartitionLinkedNumInst*/ 


template < typename T_CLUSTERIDX,
	   typename INPUT_ITERATOR,
	   typename FUNCINSTFREQUENCY
	   >
auto 
getPartitionLinkedNumInst
(INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX> &aipartition_clusters,
 const FUNCINSTFREQUENCY            func_instfrequency
 ) -> PartitionLinkedNumInst<T_CLUSTERIDX,decltype(func_instfrequency(*aiiterator_instfirst))>
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "ds::getPartitionLinkedNumInst";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "]\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  typedef decltype(func_instfrequency(*aiiterator_instfirst)) ResultType;
  
  const T_CLUSTERIDX lcidx_numClusterK = aipartition_clusters.getNumCluster();
  PartitionLinkedNumInst
    <T_CLUSTERIDX,
     ResultType
     > lopartlink_memberShip
    ((uintidx) std::distance(aiiterator_instfirst,aiiterator_instlast), 
     (uintidx) lcidx_numClusterK 
     );

  uintidx lui_i = 0;
  for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    T_CLUSTERIDX lcidx_xinK =
      aipartition_clusters.next();
    if ( 0 <= lcidx_xinK  && lcidx_xinK <  lcidx_numClusterK  ) {
      lopartlink_memberShip.addInstanceToCluster
	(lcidx_xinK,
	 lui_i,
	 func_instfrequency(*aiiterator_instfirst)
	 );
    }
    ++lui_i;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelPartLink;
    lostrstream_labelPartLink 
      << geverbosepc_labelstep << ':' << lpc_labelFunc; 
    lopartlink_memberShip.print
      (std::cout,
       lostrstream_labelPartLink.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lopartlink_memberShip;
  
} /* END getPartitionLinkedNumInst
   */

} /*END namespace ds*/

#endif /*PARTITION_LINKED_NUMINST_HPP*/
