/*! \file chromosome_base.hpp
 *
 * \brief  chromosome base
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CHROMOSOMEBASE_HPP
#define CHROMOSOMEBASE_HPP

#include <iostream>
#include "common.hpp" //uintidx
#include "verbose_global.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {


/*! \class ChromosomeBase
  \brief Chromosome Base define basic attributes for a chromosome 
*/  
template <class T_METRIC //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  >
class ChromosomeBase
{   
public:
  ChromosomeBase():
    _b_stringInvalid(false),
    _t_objetiveFunc(std::numeric_limits<T_METRIC>::max()),
    _t_fitness(-std::numeric_limits<T_METRIC>::max())
  { }

  ChromosomeBase(const T_METRIC airt_objetiveFunc, const T_METRIC airt_fitness)
    : _b_stringInvalid(false)
    , _t_objetiveFunc(airt_objetiveFunc)
    , _t_fitness(airt_fitness)
  { }

 
  //move constructor 
  ChromosomeBase(ChromosomeBase<T_METRIC> &&aich_chromosome):
    _b_stringInvalid(aich_chromosome._b_stringInvalid),
    _t_objetiveFunc(aich_chromosome._t_objetiveFunc),
    _t_fitness(aich_chromosome._t_fitness)
  {    
  }

  //copy constructor
  ChromosomeBase
  (const ChromosomeBase<T_METRIC> &aich_chromosome):
    _b_stringInvalid(aich_chromosome._b_stringInvalid),
    _t_objetiveFunc(aich_chromosome._t_objetiveFunc),
    _t_fitness(aich_chromosome._t_fitness)
  {
  }

  virtual ~ChromosomeBase() 
  { 
  }
  
 ChromosomeBase<T_METRIC>& 
  operator=(const ChromosomeBase<T_METRIC> &aich_chromosome)
  {
    if ( this != &aich_chromosome ) { 
      this->_b_stringInvalid = aich_chromosome._b_stringInvalid;
      this->_t_objetiveFunc  = aich_chromosome._t_objetiveFunc;
      this->_t_fitness  = aich_chromosome._t_fitness; 
    }

    return *this;
  }

  ChromosomeBase<T_METRIC>& 
  operator=(ChromosomeBase<T_METRIC> &&aich_chromosome)
  {
    if ( this != &aich_chromosome ) {
     
      this->_b_stringInvalid = aich_chromosome._b_stringInvalid;
      this->_t_objetiveFunc = aich_chromosome._t_objetiveFunc;
      this->_t_fitness  = aich_chromosome._t_fitness; 
    }

    return *this;
  }

  inline void setValidString(bool aib_stringInvalid) 
  {
    this->_b_stringInvalid = aib_stringInvalid;
  }
  
  inline const bool getValidString() const 
  {
    return this->_b_stringInvalid;
  }

  inline void setObjetiveFunc(T_METRIC ait_objetiveFunc) 
  {
    this->_t_objetiveFunc = ait_objetiveFunc;
  }

  inline const T_METRIC getObjetiveFunc() const 
  {
    return this->_t_objetiveFunc;
  }
  

  inline const T_METRIC getFitness() const 
  {
    return this->_t_fitness;
  }

  inline void setFitness(T_METRIC ait_fitness) 
  {
    this->_t_fitness = ait_fitness;
  }


  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label = "",
   const char aic_delimCoef=',',
   const char aic_delimRow=';'
   ) const
  {

#if defined(__VERBOSE_YES)
    os << "<CROMOSOME:"
       << geverbosepc_labelstep
       << ':' <<  aipc_label
       << ":objetive," << this->_t_objetiveFunc << aic_delimCoef
       << "fitness," << _t_fitness
       << ",id[" << geverboseui_idproc << '-' << this << ']';
#else
   os  << aipc_label
       << ":objetive," << this->_t_objetiveFunc << aic_delimCoef
       << "fitness," << _t_fitness;

#endif 
       
  }

  friend std::ostream& operator<<
  (std::ostream& os, const ChromosomeBase<T_METRIC> &aich_chromosome)
  {
    aich_chromosome.print(os); 
    
    return os;
  }

protected:

   bool       _b_stringInvalid;  /*!< an bool for check is valid chromosome */
   T_METRIC   _t_objetiveFunc;  /*!< an real number for objetive function of solution  */
   T_METRIC   _t_fitness;       /*!< an real number for chromosome fitness*/

}; /*End ChromosomeBase*/


} /*END namespace gaencode*/

#endif  /*CHROMOSOMEBASE_HPP*/
