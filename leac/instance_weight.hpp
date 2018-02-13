/*! \file instance_weight.hpp
 *
 * \brief instance with weight
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef INSTANCE_WEIGHT_HPP
#define INSTANCE_WEIGHT_HPP

#include "instance.hpp"

/*! \namespace data
  \brief Module for the model object, point or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class InstanceWeight
  \brief Patterns, usually vectors in a multidimensional space
*/
template < class T_FEATURE,
	   class T_WEIGHT
	   >
class InstanceWeight:
  public Instance<T_FEATURE>
{
public:
  InstanceWeight()
  : Instance<T_FEATURE>()
  {}

  virtual ~InstanceWeight() { }


  void setWeight(const T_WEIGHT aiT_weight)
  {
    _t_weight = aiT_weight;
  }

  T_WEIGHT getWeight() const
  {
    return _t_weight;
  } 

  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const 
  {
    Instance<T_FEATURE>::print(os,aic_delim);
    os << aic_delim << _t_weight; 
    //return os;
  }

  virtual std::string
  getToString(const char aic_delim='\t')
  {
    std::stringstream lss_instance;

    lss_instance << Instance<T_FEATURE>::getToString(aic_delim);
    lss_instance << aic_delim << _t_weight; 

    return lss_instance.str();
  }

protected:
  T_WEIGHT _t_weight;
}; 

 
} /*END namespace data 
   */

#endif /*INSTANCE_WEIGHT_HPP*/
