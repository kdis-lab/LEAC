/*! \file instance_class_frequency.hpp
 *
 * \brief instance class with frequency
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef INSTANCE_CLASS_FREQUENCY_HPP
#define INSTANCE_CLASS_FREQUENCY_HPP

#include "instance_frequency.hpp"
#include "instance_interfazclass.hpp"

/*! \namespace data
  \brief Module for the model object, point or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class InstanceClassFreq
  \brief Patterns, usually vectors in a multidimensional space, which is classified and with a number of occurrences
*/
template < class T_FEATURE,
	   class T_INSTANCE_FREQUENCY,
	   class T_INSTANCES_CLUSTER_K,
	   class T_CLUSTERIDX
	   >
class InstanceClassFreq
  : public InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>
  , public InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InstanceClassFreq()
    : InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>()
    , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}

  //copy constructor
  InstanceClassFreq
  (const InstanceClassFreq<T_FEATURE,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiinstclassfreq_b)
    : InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>(aiinstclassfreq_b)
    , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>(aiinstclassfreq_b)
  {
  }

  //move constructor
  InstanceClassFreq
  (InstanceClassFreq<T_FEATURE,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &&aiinstclassfreq_b)
    : InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>(aiinstclassfreq_b)
    , InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>(aiinstclassfreq_b)
  {
  }

  virtual ~InstanceClassFreq() {}

  InstanceClassFreq<T_FEATURE,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>&
  operator=(const InstanceClassFreq<T_FEATURE,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiinstclassfreq_b)
  {
    if( this != &aiinstclassfreq_b ) {
      InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>::operator=(aiinstclassfreq_b);
      InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::operator=(aiinstclassfreq_b);
    }
    
    return *this;
  }

  InstanceClassFreq<T_FEATURE,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>&
  operator=(InstanceClassFreq<T_FEATURE,T_INSTANCE_FREQUENCY,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &&aiinstclassfreq_b)
  {
     if( this != &aiinstclassfreq_b ) {
       InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>::operator=(aiinstclassfreq_b);
       InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::operator=(aiinstclassfreq_b);
    }
    
    return *this;
  }

  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const 
  {
    InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>::print(os,aic_delim);
    InstanceIterfazClass<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(os,aic_delim);
  }

  virtual std::string
  getToString(const char aic_delim='\t')
  {
    std::stringstream lss_instance;

    lss_instance
      <<  
      InstanceFreq
      <T_FEATURE,
       T_INSTANCE_FREQUENCY>
      ::getToString(aic_delim)
      << aic_delim
      << InstanceIterfazClass
      <T_INSTANCES_CLUSTER_K,
       T_CLUSTERIDX>
      ::_stcvector_classLabel.at
      (InstanceIterfazClass
       <T_INSTANCES_CLUSTER_K,
       T_CLUSTERIDX>
       ::_cidx_classIdx
       )->getLabel()
      ;
	
    return lss_instance.str();
  }

protected:

}; /*InstanceClassFreq*/

} /*END namespace data 
   */

#endif /*INSTANCE_CLASS_FREQUENCY_HPP*/
