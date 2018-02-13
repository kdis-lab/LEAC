/*! \file chromosome_string.hpp
 *
 * \brief chromosome string
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef CHROMOSOME_STRING_HPP
#define CHROMOSOME_STRING_HPP

#include <iostream>
#include <stdexcept>
#include <limits>
#include <assert.h>        /* assert */
#include "outfilename.hpp" //outparam::OutFileName::getDelim()
#include "chromosome_base.hpp"
#include "common.hpp" //uintidx


/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {
  

/*! \class ChromosomeString
  \brief Chhromosome String define Chromosome encode for an string 
*/  
template <class T_GENE,
	  class T_METRIC //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  >
class ChromosomeString: 
    public ChromosomeBase<T_METRIC>
{   
public:
   /*Constructed a null chromosome
   */
  ChromosomeString()
    :  ChromosomeBase<T_METRIC>()
    
  { }

  ChromosomeString(const T_METRIC airt_objetiveFunc, const T_METRIC airt_fitness)
    :  ChromosomeBase<T_METRIC>(airt_objetiveFunc, airt_fitness)
  { }

  //move constructor 
  ChromosomeString(ChromosomeString<T_GENE,T_METRIC> &&aichrom_b)
    :  ChromosomeBase<T_METRIC>(aichrom_b)
  {}

  //copy constructor
  ChromosomeString
  (const ChromosomeString<T_GENE,T_METRIC> &aichrom_b)
    :  ChromosomeBase<T_METRIC>(aichrom_b)
  {}

  virtual ~ChromosomeString() 
  {}

  ChromosomeString<T_GENE,T_METRIC>& 
  operator=(const ChromosomeString<T_GENE,T_METRIC> &aichrom_b)
  {
     if ( this != &aichrom_b ) {
       ChromosomeBase<T_METRIC>::operator=(aichrom_b);
     }

     return *this;
  }

  ChromosomeString<T_GENE,T_METRIC>& 
   operator=(ChromosomeString<T_GENE,T_METRIC> &&aichrom_b)
  {
    if ( this != &aichrom_b ) {
      ChromosomeBase<T_METRIC>::operator=(aichrom_b);
    }
    return *this;
  }
  
  virtual const uintidx getStringSize() const = 0;

  virtual void setString(T_GENE *aips_string) = 0;
    
  virtual T_GENE* getString() = 0;
 
  virtual const T_GENE* getString() const = 0;

  virtual const T_GENE getGene(const uintidx aiuintidx_idxGene) const = 0;
 
  virtual void setGene(uintidx aiuintidx_idxGene, T_GENE aiT_newGene) = 0; 
   
  
protected:


}; //End ChromosomeString


} /*END namespace gaencode*/

#endif  /*CHROMOSOME_STRING_HPP*/
