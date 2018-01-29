/*! \file instance_interfazclass.hpp
 *
 * \brief Instance with class
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INSTANCE_INTERFAZCLASS_HPP
#define INSTANCE_INTERFAZCLASS_HPP


#include "common_clustering.hpp"
#include "count_label.hpp"


/*! \namespace data
  \brief Module for the model object, point or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class InstanceIterfazClass
  \brief Instance with class
*/
template < class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
class InstanceIterfazClass
{
public:

  typedef std::map
    <const std::string, 
     CountLabel
     <T_CLUSTERIDX,
     T_INSTANCES_CLUSTER_K> >    MapClassLabel;

  InstanceIterfazClass()
    :  _cidx_classIdx(UNKNOWN_CLUSTER_IDX)
  {} 

  //copy constructor
  InstanceIterfazClass
  (const InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiinstiterclass_b)
    : _cidx_classIdx(aiinstiterclass_b._cidx_classIdx)
  { }

  //move constructor
  InstanceIterfazClass
  (InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &&aiinstiterclass_b)
    : _cidx_classIdx(aiinstiterclass_b._cidx_classIdx)
  {
    aiinstiterclass_b._cidx_classIdx =  UNKNOWN_CLUSTER_IDX;
  }

  virtual ~InstanceIterfazClass() 
  {
  }


  InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>&
  operator=(const InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiinstiterclass_b)
  {

    if( this != &aiinstiterclass_b )
      _cidx_classIdx = aiinstiterclass_b._cidx_classIdx;
    
    return *this;
  }

  InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>&
  operator=(InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &&aiinstiterclass_b)
  {

    if( this != &aiinstiterclass_b )
      _cidx_classIdx = aiinstiterclass_b._cidx_classIdx;
    
    return *this;
  }

  
  static void initialize()
  {
    _stcvector_classLabel.clear();
    _stmap_instanceClass.clear();
    _idxmcT_consecutiveClass = 0;
  }
           
  void setClassIdx(std::string &aistr_keyMapClass)
  {  
    typename MapClassLabel::iterator limap_iterClassLabel = 
      _stmap_instanceClass.find( aistr_keyMapClass );
    if  ( limap_iterClassLabel != _stmap_instanceClass.end() ) {
      limap_iterClassLabel->second.increasesCountLabel();
      _cidx_classIdx = limap_iterClassLabel->second.getIdxLabel();
    }
    else { /*NEW CLASS LABEL*/
      _cidx_classIdx = _idxmcT_consecutiveClass;
      _stmap_instanceClass.insert
	(std::pair
	 <std::string,
	  CountLabel<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> 
	  >
	(aistr_keyMapClass,
	 CountLabel
	 <T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>
	 (aistr_keyMapClass,_idxmcT_consecutiveClass)
	 ) 
	 );
      ++_idxmcT_consecutiveClass;
      
    } /* IF END CLASS LABEL*/

  }

  static void setVectorClassLabel()
  {
     _stcvector_classLabel.resize(_stmap_instanceClass.size());
    
    for (typename std::map<const std::string, 
	   CountLabel<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K> >::iterator it = 
	   _stmap_instanceClass.begin(); 
	 it != _stmap_instanceClass.end(); ++it) 
      {
	_stcvector_classLabel.at(it->second.getIdxLabel()) = &it->second;
      }
   }

  static const std::vector
  <CountLabel
   <T_CLUSTERIDX,
    T_INSTANCES_CLUSTER_K>* >& getVectorClassLabel()
  {
            
    return  _stcvector_classLabel;
  }

 
  T_CLUSTERIDX 
  getClassIdx()
  {
    return this->_cidx_classIdx;
  }

  static std::vector<std::string>  
  getLabelClass
  (uintidx  aiui_numClass)
  {
    std::vector<std::string>   lovectorstr_labelClass;
    
    if (_stcvector_classLabel.size() > 0 ) { 
      lovectorstr_labelClass.reserve(_stcvector_classLabel.size());
      for ( auto  liter_classLabel: _stcvector_classLabel ) {
	lovectorstr_labelClass.push_back
	  (liter_classLabel->getLabel());
      }
    }
    else {
      std::stringstream lss_labelClass;

      lovectorstr_labelClass.reserve(aiui_numClass);
      for (uintidx  luintidx_j = 0; luintidx_j < aiui_numClass; luintidx_j++) {  
	lss_labelClass.str("");
	lss_labelClass << "cluster_c" << luintidx_j;
	
	lovectorstr_labelClass.push_back( lss_labelClass.str() ); 
      }
    } 
   
    return lovectorstr_labelClass;

  }

  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const
  {
    os  << aic_delim
	<< InstanceIterfazClass
	<T_INSTANCES_CLUSTER_K,
	 T_CLUSTERIDX>
	::_stcvector_classLabel.at
	(InstanceIterfazClass
	 <T_INSTANCES_CLUSTER_K,
	  T_CLUSTERIDX>
	 ::_cidx_classIdx
	 )->getLabel();
      }

protected:
  
  T_CLUSTERIDX          _cidx_classIdx;
  static MapClassLabel  _stmap_instanceClass;
  static T_CLUSTERIDX   _idxmcT_consecutiveClass;

  static std::vector
  <CountLabel<T_CLUSTERIDX,T_INSTANCES_CLUSTER_K>* > _stcvector_classLabel;

}; /*END CLASS InstanceIterfazClass
    */

template < class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
typename InstanceIterfazClass
<T_INSTANCES_CLUSTER_K, 
 T_CLUSTERIDX>::MapClassLabel
InstanceIterfazClass
<T_INSTANCES_CLUSTER_K,
T_CLUSTERIDX>::_stmap_instanceClass;

template < class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
T_CLUSTERIDX InstanceIterfazClass
<T_INSTANCES_CLUSTER_K,
T_CLUSTERIDX>::_idxmcT_consecutiveClass = 0;

template < class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
std::vector
<CountLabel
 <T_CLUSTERIDX,
  T_INSTANCES_CLUSTER_K>* > 
InstanceIterfazClass
<T_INSTANCES_CLUSTER_K,
 T_CLUSTERIDX>::_stcvector_classLabel;

} /*END namespace data 
   */

#endif /*INSTANCE_INTERFAZCLASS_HPP*/
