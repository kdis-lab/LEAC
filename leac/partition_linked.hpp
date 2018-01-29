/*! \file partition_linked.hpp
 *
 * \brief Linked partition data structure
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef PARTITION_LINKED_HPP
#define PARTITION_LINKED_HPP

#include <vector>
#include <utility>
#include "nearestinstance_operator.hpp" //NEARESTCENTROID_UNKNOWN
#include "partition.hpp"
#include "linear_algebra_level1.hpp"
#include "common.hpp"

#include "verbose_global.hpp"

/*! \namespace ds
  \brief Data structure
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace ds {

  
template < typename T_CLUSTERIDX >  class IteratorPartitionLinked; 

/*! \class PartitionLinked
  \brief Data structure to define the instances that belong to a cluster consecutively
  \details 
*/
template < typename T_CLUSTERIDX >   //-1, 0, 1, .., K 
class PartitionLinked  {
public:
  PartitionLinked()
  : _vectorInstIdx_firstClusterK()
  , _vectorInstIdx_nextClusterK()
  {}
  
  PartitionLinked
  (const uintidx aiui_numInstances, 
   const uintidx aiui_numClusterK 
   ) : _vectorInstIdx_firstClusterK(aiui_numClusterK,UINTIDX_NIL)
     , _vectorInstIdx_nextClusterK(aiui_numInstances,UINTIDX_NIL)
  {
  }

  //copy constructor
  PartitionLinked(const PartitionLinked<T_CLUSTERIDX> &aipartlink_b)
    : _vectorInstIdx_firstClusterK(aipartlink_b._vectorInstIdx_firstClusterK)
    , _vectorInstIdx_nextClusterK(aipartlink_b._vectorInstIdx_nextClusterK)
  {
  }
  
  //move constructor
  PartitionLinked(PartitionLinked<T_CLUSTERIDX>&& aipartlink_b)
    : _vectorInstIdx_firstClusterK(aipartlink_b._vectorInstIdx_firstClusterK)
    , _vectorInstIdx_nextClusterK(aipartlink_b._vectorInstIdx_nextClusterK)
  {
  }
    
  ~PartitionLinked() 
  {}

  void resize(const uintidx  aiui_numClusterK)
  {
    if ( _vectorInstIdx_firstClusterK.size() < aiui_numClusterK) {
      _vectorInstIdx_firstClusterK.resize(aiui_numClusterK,UINTIDX_NIL);
    }
    else if ( _vectorInstIdx_firstClusterK.size() > aiui_numClusterK) {
      _vectorInstIdx_firstClusterK.resize(aiui_numClusterK);
    }
  }
  
  inline const T_CLUSTERIDX getNumPartitions() const
  {
    return T_CLUSTERIDX(_vectorInstIdx_firstClusterK.size());
  }

  inline const uintidx getFirstInstClusterK(T_CLUSTERIDX aicidx_clusterK)
  {  
    return  _vectorInstIdx_firstClusterK.at(aicidx_clusterK);
  }
  
  inline const uintidx getNumInstances() const
  {
    return (uintidx) _vectorInstIdx_nextClusterK.size();
  }

  
  PartitionLinked<T_CLUSTERIDX>& 
  operator=(const PartitionLinked<T_CLUSTERIDX>& aipartlink_b)
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "PartitionLinked::operator(copy)=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
      		<< "(output PartitionLinked: this[" << this << "]\n"
	        << " input  PartitionLinked: aipartlink_b["
		<< &aipartlink_b << "]\n" 
		<< ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    if ( this != &aipartlink_b ) {

      if (  _vectorInstIdx_nextClusterK.size() == aipartlink_b._vectorInstIdx_nextClusterK.size())  {
	interfacesse::copy
	  (_vectorInstIdx_nextClusterK.data(),
	   aipartlink_b._vectorInstIdx_nextClusterK.data(),
	   (uintidx) _vectorInstIdx_nextClusterK.size()
	   );
      }
      else {
	_vectorInstIdx_nextClusterK = aipartlink_b._vectorInstIdx_nextClusterK;
      }

      if ( _vectorInstIdx_firstClusterK.size() == aipartlink_b._vectorInstIdx_firstClusterK.size() ) {
	interfacesse::copy
	  (_vectorInstIdx_firstClusterK.data(),
	   aipartlink_b._vectorInstIdx_firstClusterK.data(),
	   (uintidx) _vectorInstIdx_firstClusterK.size()
	   );
      }
      else {
	_vectorInstIdx_firstClusterK = aipartlink_b._vectorInstIdx_firstClusterK;
      }

    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      PartitionLinked<T_CLUSTERIDX>::print
	(std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return *this;
    
  }

  PartitionLinked<T_CLUSTERIDX>&
  operator=(PartitionLinked<T_CLUSTERIDX>&& aipartlink_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "PartitionLinked::operator(move)=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
      		<< "(output PartitionLinked: this[" << this << "]\n"
	        << " input  PartitionLinked: aipartlink_b["
		<< &aipartlink_b << "]\n" 
		<< ')'
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
    
    if ( this != &aipartlink_b ) {
      _vectorInstIdx_nextClusterK = aipartlink_b._vectorInstIdx_nextClusterK;
      _vectorInstIdx_firstClusterK = aipartlink_b._vectorInstIdx_firstClusterK;       
    }


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      PartitionLinked<T_CLUSTERIDX>::print
	(std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    return *this;
   
  }

  inline void initialize()
  {
    interfacesse::copya
      (_vectorInstIdx_nextClusterK.data(), 
       UINTIDX_NIL,
       (uintidx) _vectorInstIdx_nextClusterK.size()
       );
    
    interfacesse::copya
      (_vectorInstIdx_firstClusterK.data(),
       UINTIDX_NIL,
       (uintidx) _vectorInstIdx_firstClusterK.size()
       );
  }

  inline void addInstanceToCluster
  (T_CLUSTERIDX  aiT_clusterIdx, 
   uintidx              aiidxiT_instanceI
   )
  {
    this->_vectorInstIdx_nextClusterK[aiidxiT_instanceI] 
      = this->_vectorInstIdx_firstClusterK[aiT_clusterIdx];
    this->_vectorInstIdx_firstClusterK[aiT_clusterIdx] = aiidxiT_instanceI;
  }

  void subInstanceFromCluster
  (T_CLUSTERIDX aiidxmcT_clusterIdx,
   uintidx             aiidxinstT_instanceI 
   ) 
  {
    uintidx   lidxinstT_iterator;

    if ( aiidxmcT_clusterIdx != NEARESTCENTROID_UNKNOWN ) { //IF
      if ( this->_vectorInstIdx_firstClusterK[aiidxmcT_clusterIdx] == aiidxinstT_instanceI )
	{
	  this->_vectorInstIdx_firstClusterK[aiidxmcT_clusterIdx] 
	    = this->_vectorInstIdx_nextClusterK[aiidxinstT_instanceI];
	  this->_vectorInstIdx_nextClusterK[aiidxinstT_instanceI] 
	    = UINTIDX_NIL; 
	}
      else {
        lidxinstT_iterator = 
	  this->_vectorInstIdx_firstClusterK[aiidxmcT_clusterIdx];
	while ( this->_vectorInstIdx_nextClusterK[lidxinstT_iterator] 
		!= aiidxinstT_instanceI) {
	  lidxinstT_iterator = 
	    this->_vectorInstIdx_nextClusterK[lidxinstT_iterator];
	}

	this->_vectorInstIdx_nextClusterK[lidxinstT_iterator] 
	  = this->_vectorInstIdx_nextClusterK[aiidxinstT_instanceI];
      
	this->_vectorInstIdx_nextClusterK[aiidxinstT_instanceI] = UINTIDX_NIL; 
      }
      
    } //IF 
   
  }
  
  void changeMemberShip
  (T_CLUSTERIDX  aicidx_nearestClusterK,
   uintidx              aiidxinstT_instanceI)
  {
    this->subInstanceFromCluster(aiidxinstT_instanceI);
    if ( aicidx_nearestClusterK != NEARESTCENTROID_UNKNOWN ) {
      this->addInstanceToCluster(aicidx_nearestClusterK,aiidxinstT_instanceI);
    }
  }
 
  /* getInstanceIdxOfClusterK:
     Searching for instance m of a cluster and returns the index of 
     the instance. If the instance is out of range returns UINTIDX_NIL
  */
  uintidx getInstanceIdxOfClusterK
  (T_CLUSTERIDX  aicidx_clusterK,
   const uintidx aiuintidx_instanceNi
   ) 
  {
    uintidx luintidx_contInstance = 0;
    uintidx lidxinstT_iterator = this->_vectorInstIdx_firstClusterK[aicidx_clusterK];

    while ( lidxinstT_iterator != UINTIDX_NIL && 
	    luintidx_contInstance < aiuintidx_instanceNi ) {
      ++luintidx_contInstance;
      lidxinstT_iterator  = this->_vectorInstIdx_nextClusterK[lidxinstT_iterator];
    }
    
    return lidxinstT_iterator;
  }
  
  void joinCluster
  (const T_CLUSTERIDX aiIdxK_clusterKFrom,
   const T_CLUSTERIDX aiIdxK_clusterKTo
   )
  {
   uintidx  luintidx_iterInst = this->_vectorInstIdx_firstClusterK[aiIdxK_clusterKTo];

    if ( luintidx_iterInst == UINTIDX_NIL ) {
      this->_vectorInstIdx_firstClusterK[aiIdxK_clusterKTo] 
	= this->_vectorInstIdx_firstClusterK[aiIdxK_clusterKFrom];
    }
    else {
      while ( this->_vectorInstIdx_nextClusterK[luintidx_iterInst] 
	      != UINTIDX_NIL )  {
	luintidx_iterInst = this->_vectorInstIdx_nextClusterK[luintidx_iterInst];
      }
      this->_vectorInstIdx_nextClusterK[luintidx_iterInst] 
	= this->_vectorInstIdx_firstClusterK[aiIdxK_clusterKFrom];
    }
    this->_vectorInstIdx_firstClusterK[aiIdxK_clusterKFrom] =  UINTIDX_NIL;
  }

  void  print
  (std::ostream &os=std::cout,
   const char*  aipc_label        = "",
   const char   aic_delimCoef     = ',',
   const char   aic_delimClusterK = ';'
   )
  {
    T_CLUSTERIDX lcidx_numClusterK = this->getNumPartitions();
      
    os  << "<PARTITIONLINKED:" << aipc_label
	<< ":length" << aic_delimCoef  << lcidx_numClusterK << '>';

    IteratorPartitionLinked <T_CLUSTERIDX> literpart_ip(this);
    
    for ( T_CLUSTERIDX lcidx_Ckp = 0; lcidx_Ckp < lcidx_numClusterK; lcidx_Ckp++) {
      literpart_ip.begin(lcidx_Ckp);
      if ( literpart_ip.end() ) {
	os << literpart_ip.getValue();
	
	literpart_ip.next();
      
	for ( ; literpart_ip.end(); literpart_ip.next() ) {
	  os << aic_delimCoef << literpart_ip.getValue();
	}
      }
      if ( (lcidx_Ckp + 1) != lcidx_numClusterK )
	os << aic_delimClusterK;
    }
    
  }
  
protected:
                         
  std::vector<uintidx> _vectorInstIdx_firstClusterK;  //UINTIDX_NIL, 0,...N
  std::vector<uintidx> _vectorInstIdx_nextClusterK;   //UINTIDX_NIL, 0,...N

  friend class IteratorPartitionLinked<T_CLUSTERIDX>;  
 
}; /*PartitionLinked*/ 


/*! \class IteratorPartitionLinked
  \brief  An iterator for a PartitionLinked
  \details 
*/
template < typename T_CLUSTERIDX > //-1, 0, 1, .., K  
class IteratorPartitionLinked {
public:
  IteratorPartitionLinked
  (const  PartitionLinked<T_CLUSTERIDX>   *aippartl_partitionLinked ):
    _ppartl_partitionLinked(aippartl_partitionLinked),
    _uintidx_instI(UINTIDX_NIL)
  {}

  ~IteratorPartitionLinked() 
  {}

  inline void begin(T_CLUSTERIDX aicidx_clusterK)
  {
    _uintidx_instI = 
      _ppartl_partitionLinked->_vectorInstIdx_firstClusterK[aicidx_clusterK];
  }

  inline void next() 
  {
    _uintidx_instI = 
      _ppartl_partitionLinked->_vectorInstIdx_nextClusterK[_uintidx_instI];
  }

  inline bool end() 
  {
    return _uintidx_instI != UINTIDX_NIL;
  }

  inline uintidx getValue() 
  {
    return this->_uintidx_instI;
  }

protected:
  const PartitionLinked<T_CLUSTERIDX>*  _ppartl_partitionLinked;  
  uintidx                               _uintidx_instI;
};

/* \fn template < typename T_CLUSTERIDX > PartitionLinked<T_CLUSTERIDX> getPartitionlinked(const partition::Partition<T_CLUSTERIDX> &aipartition_clusters)
   \brief Gets a PartitionLinked from encoding on a partition of instances in clusters
   \details 
   \param aipartition_clusters a partition of instances in clusters
*/
template < typename T_CLUSTERIDX >
PartitionLinked<T_CLUSTERIDX>
getPartitionlinked
(partition::Partition<T_CLUSTERIDX> &aipartition_clusters)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "getPartitionlinked";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  partition::Partition<>&: aipartition_clusters[" 
	      << &aipartition_clusters << "] K = "   << aipartition_clusters.getNumCluster()
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  const T_CLUSTERIDX lcidx_numClusterK = aipartition_clusters.getNumCluster();
  
  PartitionLinked<T_CLUSTERIDX> lopartlink_memberShip
    (aipartition_clusters.getNumInstances(),
     (uintidx) lcidx_numClusterK 
     );
  
  aipartition_clusters.begin();
  for (uintidx lui_i = 0; aipartition_clusters.end(); lui_i++) {
    T_CLUSTERIDX lcidx_xinK = aipartition_clusters.next();
  
    if ( 0 <= lcidx_xinK  && lcidx_xinK <  lcidx_numClusterK  ) {
      lopartlink_memberShip.addInstanceToCluster(lcidx_xinK,lui_i);
    }
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
  
} /*getPartitionlinked
   */

} /*END namespace ds*/

#endif /*PARTITION_LINKED_HPP*/
