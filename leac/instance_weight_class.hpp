/*! \file instance_weight_class.hpp
 *
 * \brief instance with weight and class
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef INSTANCE_WEIGHT_CLASS_HPP
#define INSTANCE_WEIGHT_CLASS_HPP

#include "instance_weight.hpp"
#include "instance_interfazclass.hpp"

/*! \namespace data
  \brief Module for the model object, point or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class InstanceWeightClass
  \brief Patterns, usually vectors in a multidimensional space
*/
template < class T_FEATURE,
	   class T_WEIGHT,
	   class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
class InstanceWeightClass
  : public InstanceWeight<T_FEATURE,T_WEIGHT>
  , public InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InstanceWeightClass()
  : InstanceWeight<T_FEATURE,T_WEIGHT>()
  , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}

  virtual ~InstanceWeightClass() {}

  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const 
  {
    
    InstanceWeight<T_FEATURE,T_WEIGHT>::print(os,aic_delim);
    os << aic_delim 
      // << '\"'
       << InstanceIterfazClass
      <	T_INSTANCES_CLUSTER_K,
	T_CLUSTERIDX>
      ::_stcvector_classLabel.at
      (InstanceIterfazClass
       <T_INSTANCES_CLUSTER_K,
       T_CLUSTERIDX>
       ::_membidx_memberCluster
       )->getLabel()
      // << '\"'
      ;  
  }

  virtual std::string
  getToString(const char aic_delim='\t')
  {
    std::stringstream lss_instance;
    
    lss_instance << InstanceWeight<T_FEATURE,T_WEIGHT>::getToString(aic_delim);
    lss_instance 
      << aic_delim 
      << InstanceIterfazClass
      <T_INSTANCES_CLUSTER_K,
	T_CLUSTERIDX>
      ::_stcvector_classLabel.at
      (InstanceIterfazClass
       <T_INSTANCES_CLUSTER_K,
       T_CLUSTERIDX>
       ::_membidx_memberCluster
       )->getLabel(); 

    return lss_instance.str();
  }

protected:

}; /*class InstanceWeightClass*/ 

 
} /*END namespace data 
   */

#endif /*INSTANCE_WEIGHT_CLASS_HPP*/
