/*! \file instance_class.hpp
 *
 * \brief instance with class 
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef INSTANCE_CLASS_HPP
#define INSTANCE_CLASS_HPP

#include "instance.hpp"
#include "instance_interfazclass.hpp"

/*! \namespace data
  \brief Module for the model object, point or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class InstanceClass
  \brief Patterns, usually vectors in a multidimensional space, previamente clasificado
*/
template < class T_FEATURE,
	   class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
class InstanceClass
  : public Instance<T_FEATURE>
  , public InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InstanceClass()
  : Instance<T_FEATURE>()
  , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {
  }

  InstanceClass
  (const InstanceClass<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiinstclass_b)
    : Instance<T_FEATURE>(aiinstclass_b)
    , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>(aiinstclass_b)
  {
  }

  InstanceClass
  (InstanceClass<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &&aiinstclass_b)
    : Instance<T_FEATURE>(aiinstclass_b)
    , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>(aiinstclass_b)
  {
  }
  
  virtual ~InstanceClass() {}

  InstanceClass<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>&
  operator=(const InstanceClass<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiinstclass_b)
  {
    if( this != &aiinstclass_b ) {
      Instance<T_FEATURE>::operator=(aiinstclass_b);
      InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::operator=(aiinstclass_b);
    }
     
    return *this;
  }

  InstanceClass<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>&
  operator=(InstanceClass<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &&aiinstclass_b)
  {

    if( this != &aiinstclass_b ) {
      Instance<T_FEATURE>::operator=(aiinstclass_b);
      InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::operator=(aiinstclass_b);
    }
     
    return *this;
    
  }
  
  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const 
  {
    Instance<T_FEATURE>::print(os,aic_delim);
    InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(os,aic_delim);
  }

  virtual std::string
  getToString(const char aic_delim='\t')
  {
    std::stringstream lss_instance;
    
    lss_instance <<  Instance<T_FEATURE>::getToString(aic_delim);
    
      lss_instance 
	<< aic_delim 
	<< InstanceIterfazClass
	<T_INSTANCES_CLUSTER_K,
	 T_CLUSTERIDX>
	::_stcvector_classLabel.at
	(InstanceIterfazClass
	 <T_INSTANCES_CLUSTER_K,
	  T_CLUSTERIDX>
	 ::_cidx_classIdx
	 )->getLabel(); 

    return lss_instance.str();
  }
  
}; /*END CLASS InstanceClass-------------------------------------------------------
    */

} /*END namespace data 
   */


#endif /*INSTANCE_CLASS_HPP*/
