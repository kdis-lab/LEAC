/*! \file partition_linked_stats.hpp
 *
 * \brief Linked partition data structure with statistics
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __PARTITION_LINKED_STATS_HPP
#define __PARTITION_LINKED_STATS_HPP

#include <algorithm>    // std::count 
#include "partition_linked.hpp"
#include "matrix_resizablerow.hpp"
#include "linear_algebra_level1.hpp"

#include "verbose_global.hpp"
#include "container_out.hpp"

/*! \namespace ds
  \brief Data structure
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace ds {

/*! \class PartitionLinkedStats
  \brief Partition with an array of membership 
  \details 
*/
template < typename T_FEATURE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_INSTANCE_FREQUENCY,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_FEATURE_SUM
	   > 
class PartitionLinkedStats: 
  public  PartitionLinked<T_CLUSTERIDX> 
{
public:
  PartitionLinkedStats
  (const uintidx  aiui_numClusterK,
   const uintidx  aiui_numInstances,
   const uintidx  aiui_numDimensionsInstances, 
   const uintidx  aiui_numClusterKSpace // 0, 1, .., K, .., NumClusterKSpace
   ): PartitionLinked<T_CLUSTERIDX>(aiui_numInstances,aiui_numClusterKSpace)
    , _arrayidxk_memberShip(new T_CLUSTERIDX[aiui_numInstances])
    , _matrixresizerow_m
      (aiui_numClusterK,
       aiui_numDimensionsInstances,
       aiui_numClusterKSpace
       )
    , _vectorit_n(aiui_numClusterKSpace,0)
  {
  
    interfacesse::copya
      (this->_arrayidxk_memberShip, 
       T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN),
       aiui_numInstances
       );

    this->_matrixresizerow_m.initialize();
    
  }

  //copy constructor
  PartitionLinkedStats
  (const PartitionLinkedStats
   <T_FEATURE, 
   T_CLUSTERIDX, 
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM
   >                         &aipartlinkstats_b)
    : PartitionLinked<T_CLUSTERIDX>(aipartlinkstats_b)
    , _arrayidxk_memberShip(new T_CLUSTERIDX[aipartlinkstats_b.getNumInstances()])
    , _matrixresizerow_m(aipartlinkstats_b._matrixresizerow_m)
    , _vectorit_n(aipartlinkstats_b._vectorit_n)
  {
    
    interfacesse::copy
      (_arrayidxk_memberShip,
       aipartlinkstats_b._arrayidxk_memberShip,
       aipartlinkstats_b.getNumInstances()
       );
  }

  //move constructor
  PartitionLinkedStats
  (PartitionLinkedStats
   <T_FEATURE, 
   T_CLUSTERIDX, 
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM
   >                         &&aipartlinkstats_b)
    : PartitionLinked<T_CLUSTERIDX>(aipartlinkstats_b)
    , _arrayidxk_memberShip(aipartlinkstats_b._arrayidxk_memberShip)
    , _matrixresizerow_m(aipartlinkstats_b._matrixresizerow_m)
    , _vectorit_n(aipartlinkstats_b._vectorit_n)
  {
    aipartlinkstats_b._arrayidxk_memberShip = NULL;
  }
  
  ~PartitionLinkedStats() 
  {
    if ( _arrayidxk_memberShip != NULL ) 
      delete[] _arrayidxk_memberShip;
  }
 
  inline void initialize()   
  {
    PartitionLinked<T_CLUSTERIDX>
      ::initialize(); 
     interfaceclapackext_copya
      (this->_arrayidxk_memberShip, 
       T_CLUSTERIDX(NEARESTCENTROID_UNKNOWN),
       this->getNumInstances()
       );
  }

  inline const T_CLUSTERIDX getNumPartitions() const
  {
    return T_CLUSTERIDX(_matrixresizerow_m.getNumRows());
  }
  
  inline void decreasePartition()
  {
    _matrixresizerow_m.decreaseRow();
  }

  inline void increasesPartition()
  {
    _matrixresizerow_m.increasesRow();
  }

  void resize(const uintidx  aiui_numClusterK)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "PartitionLinkedStats::resize";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n";
      PartitionLinkedStats
	<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM>
	::print
	(std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout << "\n input  aiui_numClusterK: " << aiui_numClusterK
		<< '\n' 
		<< ')'
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    PartitionLinked<T_CLUSTERIDX>::resize(aiui_numClusterK);
    if ( this->_matrixresizerow_m.getNumRowsMaximum() < aiui_numClusterK )  {
      _vectorit_n.resize(aiui_numClusterK,T_INSTANCES_CLUSTER_K(0));  
    }
    else if ( this->_matrixresizerow_m.getNumRowsMaximum() > aiui_numClusterK )  {
      _vectorit_n.resize(aiui_numClusterK);  
    }
    this->_matrixresizerow_m.resize(aiui_numClusterK);
     
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      PartitionLinkedStats
	<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM>
	::print
	(std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }
  

  inline
  T_CLUSTERIDX 
  numClusterEmpty()
  {
    return (T_CLUSTERIDX) 
      std::count
      (&this->_vectorit_n[0],
       &this->_vectorit_n[_matrixresizerow_m.getNumRows()], 
       T_INSTANCES_CLUSTER_K(0)
       );
  }

  inline 
  mat::MatrixResizableRow<T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>&
  getSumInstancesClusterK()
  {
    return this->_matrixresizerow_m;
  }

  inline
  T_FEATURE_SUM* 
  getSumInstancesClusterKi(T_CLUSTERIDX aicidx_clusterK)
  {
    return this->_matrixresizerow_m.getRow((uintidx) aicidx_clusterK);
  }

  inline
  std::vector<T_INSTANCES_CLUSTER_K>& 
  getNumInstancesClusterK()
  {
    return this->_vectorit_n;
  }

  inline
  T_INSTANCES_CLUSTER_K
  getNumInstancesClusterKi(T_CLUSTERIDX aicidx_clusterK)
  {
    return this->_vectorit_n.at(aicidx_clusterK);
  }
  
  inline 
  T_CLUSTERIDX* 
  getMembersShip()
  {
    return this->_arrayidxk_memberShip;
  }

  inline 
  T_CLUSTERIDX 
  getMemberShip(const uintidx aiui_idxinstanceI ) const
  {
    return this->_arrayidxk_memberShip[aiui_idxinstanceI];
  }

  inline 
  void 
  setMemberShip(const uintidx aiui_idxinstanceI, const T_CLUSTERIDX aiidxk_clusterK)
  {
    _arrayidxk_memberShip[aiui_idxinstanceI] = aiidxk_clusterK;
  }
  
  PartitionLinkedStats
   <T_FEATURE, 
   T_CLUSTERIDX, 
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM
   >&
  operator=
  (const PartitionLinkedStats
   <T_FEATURE, 
   T_CLUSTERIDX, 
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM
   >                         &aipartlinkstats_b
   )
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "PartitionLinkedStats::(copy)=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
      		<< "(output PartitionLinkedStats: this[" << this << "]\n"
	        << " input  PartitionLinkedStats: aipartlinkstats_b["
		<< &aipartlinkstats_b << "]\n" 
		<< ')'
		<< std::endl;
    }
#endif //__VERBOSE_YES

    if ( this != &aipartlinkstats_b ) {

      PartitionLinked<T_CLUSTERIDX>::
	operator=(aipartlinkstats_b);
    
      interfacesse::copy
	(this->_arrayidxk_memberShip,
	 aipartlinkstats_b._arrayidxk_memberShip,
	 this->getNumInstances()
	 );

         if ( _vectorit_n.size() == aipartlinkstats_b._vectorit_n.size() ) {
	interfacesse::copy
	  (_vectorit_n.data(),
	   aipartlinkstats_b._vectorit_n.data(),
	   (uintidx) _vectorit_n.size()
	   );
      }
      else {
	_vectorit_n = aipartlinkstats_b._vectorit_n;
      }  
      _matrixresizerow_m = aipartlinkstats_b._matrixresizerow_m;
    
    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      PartitionLinkedStats
	<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM>
	::print
	(std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    return *this;
  }

  PartitionLinkedStats
  <T_FEATURE, 
   T_CLUSTERIDX, 
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM
   >&
  operator=
  (PartitionLinkedStats
   <T_FEATURE, 
   T_CLUSTERIDX, 
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM
   >                         &&aipartlinkstats_b
   )
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "PartitionLinkedStats::(move)=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
      		<< "(output PartitionLinkedStats: this[" << this << "]\n"
	        << " input  PartitionLinkedStats: aipartlinkstats_b["
		<< &aipartlinkstats_b << "]\n" 
		<< ')'
		<< std::endl;
    }
#endif //__VERBOSE_YES

    if ( this !=  &aipartlinkstats_b ) {
      
      PartitionLinked<T_CLUSTERIDX>::
	operator=(aipartlinkstats_b);
      this->_arrayidxk_memberShip = aipartlinkstats_b._arrayidxk_memberShip;
      
      aipartlinkstats_b._arrayidxk_memberShip = NULL;

      _vectorit_n =  aipartlinkstats_b._vectorit_n;
      _matrixresizerow_m = aipartlinkstats_b._matrixresizerow_m;
       
    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      PartitionLinkedStats
	<T_FEATURE,T_CLUSTERIDX,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_FEATURE_SUM>
	::print
	(std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    return *this;
 
  }
  
  
  void addInstanceToCluster
  (T_CLUSTERIDX               aicidx_toClusterIdx, 
   uintidx                    aiui_idxInstance,
   const T_FEATURE            *aiarrayt_feature,
   const T_INSTANCE_FREQUENCY aiit_instfrequency = 1
   )  
  {
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionLinkedStats::addInstanceToCluster:  IN"
		<< '(' << geiinparam_verbose << ")\n"
      		<< "(output PartitionLinkedStats: this[" << this << "]\n"
		<< " input  T_INSTANCES_IDX aiui_idxInstance = " 
		<<  aiui_idxInstance << '\n'  
		<< " input  aiarrayt_feature: [" 
		<<  aiarrayt_feature << "]\n"
		<< ")"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    PartitionLinked<T_CLUSTERIDX>
      ::addInstanceToCluster(aicidx_toClusterIdx,aiui_idxInstance);
  
    this->_arrayidxk_memberShip[aiui_idxInstance] = aicidx_toClusterIdx;

    T_FEATURE_SUM *larrayrowt_sumInstancesClusterK  = 
      this->_matrixresizerow_m.getRow( (uintidx)aicidx_toClusterIdx );
    T_FEATURE lT_instanceFrec = (T_FEATURE) aiit_instfrequency;
    interfacesse::axpy
      (larrayrowt_sumInstancesClusterK,
       lT_instanceFrec,
       aiarrayt_feature,
       data::Instance<T_FEATURE>::getNumDimensions()
       ); 
    ++_vectorit_n.at(aicidx_toClusterIdx);
   
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionLinkedStats::addInstanceToCluster: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*addInstanceToCluster*/

  void subInstanceFromCluster
  (uintidx                    aiui_idxInstance,
   const T_FEATURE            *aiarrayt_feature,
   const T_INSTANCE_FREQUENCY aiit_instfrequency = 1
   )
  {

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionLinkedStats::subInstanceFromCluster:  IN"
		<< '(' << geiinparam_verbose << ")\n"
      		<< "(output PartitionLinkedStats: this[" << this << "]\n"
		<< " input uintidx aiui_idxInstance = " 
		<<  aiui_idxInstance << '\n'  
		<< " input  aiarrayt_feature [" 
		<<  aiarrayt_feature << "]\n"
		<< ")"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    T_CLUSTERIDX lT_clusterIdx  
      = this->_arrayidxk_memberShip[aiui_idxInstance];
   
    if ( lT_clusterIdx != NEARESTCENTROID_UNKNOWN ) { //IF 
      PartitionLinked<T_CLUSTERIDX>::subInstanceFromCluster
      (lT_clusterIdx,aiui_idxInstance);
      
      this->_arrayidxk_memberShip[aiui_idxInstance] = NEARESTCENTROID_UNKNOWN;

      T_FEATURE_SUM *larrayrowt_sumInstancesClusterK =
      this->_matrixresizerow_m.getRow((uintidx) lT_clusterIdx);
      T_FEATURE lT_instanceFrec = 
      (T_FEATURE) aiit_instfrequency;
    interfacesse::axpy
      (larrayrowt_sumInstancesClusterK,
       -lT_instanceFrec,             
       aiarrayt_feature,
       data::Instance<T_FEATURE>::getNumDimensions()
       ); 
    --_vectorit_n.at(lT_clusterIdx);
    } //IF 
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionLinkedStats::subInstanceFromCluster: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  void changeSumInstances
  (const T_CLUSTERIDX         aicidx_clusterKFrom,
   const T_CLUSTERIDX         aicidx_clusterKTo,
   const T_FEATURE            *aiarrayt_feature,
   const T_INSTANCE_FREQUENCY aiit_instfrequency = 1
   )
  {
      T_FEATURE_SUM *larrayrowt_sumInstancesClusterK =
      this->_matrixresizerow_m.getRow((uintidx) aicidx_clusterKFrom);

      T_FEATURE lT_instanceFrec = (T_FEATURE) aiit_instfrequency;
      
    interfacesse::axpy
      (larrayrowt_sumInstancesClusterK,
       -lT_instanceFrec,             
       aiarrayt_feature,
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    
    --_vectorit_n.at(aicidx_clusterKFrom);

    larrayrowt_sumInstancesClusterK  = 
      this->_matrixresizerow_m.getRow( (uintidx) aicidx_clusterKTo );

    interfacesse::axpy
      (larrayrowt_sumInstancesClusterK,
       lT_instanceFrec,
       aiarrayt_feature,
       data::Instance<T_FEATURE>::getNumDimensions()
       );
    
    ++_vectorit_n.at(aicidx_clusterKTo);
  }
      
  void changeMemberShip 
  (T_CLUSTERIDX               aicidx_toClusterIdx,
   uintidx                    aiui_idxInstance,
   const T_FEATURE            *aiarrayt_feature,
   const T_INSTANCE_FREQUENCY aiit_instfrequency = 1
  )   
  {

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionLinkedStats::changeMemberShip:  IN"
		<< '(' << geiinparam_verbose << ")\n"
      		<< "\t(output PartitionLinkedStats: this[" << this << "]\n"
		<< "\t input  T_CLUSTERIDX aicidx_toClusterIdx = " 
		<<  aicidx_toClusterIdx << '\n'
		<< "\t input  uintidx aiui_idxInstance = " 
		<<  aiui_idxInstance << '\n'  
		<< "\t input  aiarrayt_feature[" 
		<<  aiarrayt_feature << "]\n"
		<< "\t)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    this->subInstanceFromCluster
      (aiui_idxInstance,
       aiarrayt_feature,
       aiit_instfrequency
       );
    if ( aicidx_toClusterIdx != NEARESTCENTROID_UNKNOWN ) {
      this->addInstanceToCluster
	(aicidx_toClusterIdx,
	 aiui_idxInstance,
	 aiarrayt_feature,
	 aiit_instfrequency
	 ); 
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionLinkedStats::changeMemberShip: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  void  print
  (std::ostream &os=std::cout,
   const char*  aipc_label    = "",
   const char   aic_delimCoef = ','
   )
  {

#ifdef __VERBOSE_YES
    IteratorPartitionLinked <T_CLUSTERIDX> literpart_ip(this);
    
    T_CLUSTERIDX lcidx_numClusterK = this->getNumPartitions();
    
    for ( T_CLUSTERIDX lcidx_Ckp = 0; lcidx_Ckp < lcidx_numClusterK; lcidx_Ckp++) {
      T_INSTANCES_CLUSTER_K  lti_numInstOk = 0;	
      for (literpart_ip.begin(lcidx_Ckp); literpart_ip.end(); literpart_ip.next() ) {
	if ( this->_arrayidxk_memberShip[literpart_ip.getValue()] == lcidx_Ckp ) {
	  ++lti_numInstOk;
	}
      }
      if ( lti_numInstOk != _vectorit_n[lcidx_Ckp]) {
	std::cout << "\nERROR: PARTITION_MEMBERCLUSTER lui_numInstOk: K= "
		  << lcidx_Ckp << ": " <<lti_numInstOk  << '/'  << this->getNumInstances()
		  << std::endl;
      }	
    }    
#endif /*__VERBOSE_YES*/

    std::ostringstream lostrstream_labelPartition;
    lostrstream_labelPartition
      << aipc_label
      << ':' << "NumPartition" << aic_delimCoef << this->getNumPartitions();
      
    PartitionLinked<T_CLUSTERIDX>::print(os,lostrstream_labelPartition.str().c_str(),aic_delimCoef);

    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "\n<PARTITION_MEMBERCLUSTER:" << aipc_label;
    inout::containerprint
      (this->_arrayidxk_memberShip,
       this->_arrayidxk_memberShip + this->getNumInstances(),
       os,
       lostrstream_labelMemberShip.str().c_str(),
       aic_delimCoef
       );
    os << '\n';
     std::ostringstream lostrstream_labelNumInstClusterK;
    lostrstream_labelNumInstClusterK << "<NUMBERINSTANCESCLUSTERK:" << aipc_label;
    inout::containerprint
      (_vectorit_n.begin(),
       _vectorit_n.end(),
       os,
       lostrstream_labelNumInstClusterK.str().c_str(),
       aic_delimCoef
       );
    os << '\n';
    std::ostringstream lostrstream_labelSumInstClusterK;
    lostrstream_labelSumInstClusterK << "<SUMINSTANCESCLUSTER:" << aipc_label;
    this->_matrixresizerow_m.print
      (os,
       lostrstream_labelSumInstClusterK.str().c_str(),
       ',',
       ';'
       );
   
  }
  
protected:
  
  T_CLUSTERIDX     *_arrayidxk_memberShip;

  /*Contain sum of instances or point per cluster d x k
   */
  mat::MatrixResizableRow<T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>  _matrixresizerow_m;
  /*Number of instances or point per cluster  k x 1  VALUES 
    -1,0,1,..,K,.., NumClusterKSpace 
  */
  std::vector<T_INSTANCES_CLUSTER_K> _vectorit_n;
  
}; /*PartitionLinkedStats*/

} /*END namespace ds*/

#endif /*__PARTITION_LINKED_STATS_HPP*/
