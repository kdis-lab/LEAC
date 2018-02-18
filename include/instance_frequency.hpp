/*! \file instance_frequency.hpp
 *
 * \brief instance with frequency 
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INSTANCE_FREQUENCY_HPP
#define INSTANCE_FREQUENCY_HPP

#include "instance.hpp"

/*! \namespace data
  \brief Module for the model object, point or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class InstanceFreq
  \brief Patterns, usually vectors in a multidimensional space, with an attribute to store occurrences of instances
*/
template < class T_FEATURE,
	   class T_INSTANCE_FREQUENCY
	   >
class InstanceFreq:
  public Instance<T_FEATURE>
{
public:

  InstanceFreq()
    : Instance<T_FEATURE>()
    , _t_frequency(1)
  {
  }
  
  InstanceFreq
  (T_INSTANCE_FREQUENCY ait_frequency):
    Instance<T_FEATURE>(),
    _t_frequency(ait_frequency)
  {
  }

  //copy constructor
  InstanceFreq
  (const InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY> &aiinstfreq_b)
    : Instance<T_FEATURE>(aiinstfreq_b)
    , _t_frequency(aiinstfreq_b._t_frequency)
  {
  }
  
  //move constructor
  InstanceFreq
  (InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY> &aiinstfreq_b)
    : Instance<T_FEATURE>(aiinstfreq_b)
    , _t_frequency(aiinstfreq_b._t_frequency)
  {
    _t_frequency = 0;
  }

  virtual ~InstanceFreq() { }

  InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>&
  operator=(const InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY> &aiinstfreq_b)
  {
    if( this != &aiinstfreq_b ) {
      Instance<T_FEATURE>::operator=(aiinstfreq_b);
      _t_frequency = aiinstfreq_b._t_frequency;
    }
    
    return *this;
  }

  InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY>&
  operator=(InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY> &&aiinstfreq_b)
  {
    if( this != &aiinstfreq_b ) {
      Instance<T_FEATURE>::operator=(aiinstfreq_b);
      _t_frequency = aiinstfreq_b._t_frequency;
      
      aiinstfreq_b._t_frequency = 0;
    }
    
    return *this;
  }

  void setFrequency(T_INSTANCE_FREQUENCY ait_frequency)
  {
    _t_frequency = ait_frequency;
  }

  const T_INSTANCE_FREQUENCY getFrequency() const
  {
    return _t_frequency;
  } 

  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const 
  {
    Instance<T_FEATURE>::print(os,aic_delim);
    os << aic_delim << _t_frequency; 
  }

  virtual std::string
  getToString(const char aic_delim='\t')
  {
    std::stringstream lss_instance;

    lss_instance << Instance<T_FEATURE>::getToString(aic_delim);
    lss_instance << aic_delim << _t_frequency; 

    return lss_instance.str();
  }

protected:

  T_INSTANCE_FREQUENCY _t_frequency;
  
}; /*END CLASS InstanceFreq*/

} /*END namespace data*/

#endif /*INSTANCE_FREQUENCY_HPP*/
