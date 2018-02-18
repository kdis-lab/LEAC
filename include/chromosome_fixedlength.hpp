/*! \file chromosome_fixedlength.hpp
 *
 * \brief chromosome fixed length
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef CHROMOSOME_FIXEDLENGTH_HPP
#define CHROMOSOME_FIXEDLENGTH_HPP

#include <iostream>
#include <stdexcept>
#include <limits>
#include <assert.h>        
#include "outfilename.hpp" //outparam::OutFileName::getDelim()
#include "chromosome_string.hpp"
#include "linear_algebra_level1.hpp"

#include "verbose_global.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {


/*! \class ChromFixedLength
  \brief Chromosome with fixed length string 
*/  
template <class T_GENE,
	  class T_METRIC 
	  >
class ChromFixedLength: 
  public ChromosomeString<T_GENE,T_METRIC>
{   
public:
  /*Constructed a null chromosome
   */
  ChromFixedLength()
    :  ChromosomeString<T_GENE,T_METRIC>()
    , _pts_string(new T_GENE[_stcui_stringSize])
  { }

  ChromFixedLength(const T_METRIC airt_objetiveFunc, const T_METRIC airt_fitness)
    :  ChromosomeString<T_GENE,T_METRIC>(airt_objetiveFunc, airt_fitness)
    ,  _pts_string(new T_GENE[_stcui_stringSize])
  { }
  
  //move constructor 
  ChromFixedLength(ChromFixedLength<T_GENE,T_METRIC> &&aichrom_b)
    :  ChromosomeString<T_GENE,T_METRIC>(aichrom_b)
    ,  _pts_string(aichrom_b._pts_string)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "gaencode::ChromFixedLength::ChromFixedLength:move";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  gaencode::ChromFixedLength<>: &&aichrom_b[" <<  &aichrom_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromFixedLength::print(std::cout,lpc_labelFunc);
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    aichrom_b._pts_string = NULL;
    
  }

  //copy constructor
  ChromFixedLength
  (const ChromFixedLength<T_GENE,T_METRIC> &aichrom_b)
    :  ChromosomeString<T_GENE,T_METRIC>(aichrom_b)
    , _pts_string(new T_GENE[_stcui_stringSize])
  {
    if (_pts_string != NULL) {
      interfacesse::copy
	(this->_pts_string, aichrom_b._pts_string, _stcui_stringSize);
    }
  }

  virtual ~ChromFixedLength() 
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "gaencode::ChromFixedLength::ChromFixedLength~";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ':' << geverbosepc_labelstep
		<< '[' << geverboseui_idproc << ':' << this << ']'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    if ( this->_pts_string != NULL ) {
      delete[] _pts_string; 
    }
  }

  
  ChromFixedLength<T_GENE,T_METRIC>& 
  operator=(const ChromFixedLength<T_GENE,T_METRIC> &aichrom_b)
  { 
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "gaencode::ChromFixedLength::operator:copy=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  gaencode::ChromFixedLength<>: &&aichrom_b[" <<  &aichrom_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichrom_b ) {
      ChromosomeString<T_GENE,T_METRIC>::operator=(aichrom_b);	
      interfacesse::copy
	(this->_pts_string, 
	 aichrom_b._pts_string, 
	 _stcui_stringSize
	 );
    }
    

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromFixedLength::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    
    return *this;
  }

  ChromFixedLength<T_GENE,T_METRIC>& 
  operator=(ChromFixedLength<T_GENE,T_METRIC> &&aichrom_b)
  {
    if ( this != &aichrom_b ) {
      ChromosomeString<T_GENE,T_METRIC>::operator=(aichrom_b);
      if ( this->_pts_string != NULL ) {
	delete[] _pts_string; 
      }
      this->_pts_string = aichrom_b._pts_string;
      aichrom_b._pts_string = NULL;
    }

    return *this;
  }

  inline static void setStringSize(uintidx aiuintidx_stringSize) 
  {
    _stcui_stringSize = aiuintidx_stringSize;
  }
   
  inline static const uintidx stcgetStringSize()
  {
    return _stcui_stringSize;
  }

  virtual const uintidx getStringSize() const 
  {
    return _stcui_stringSize;
  }

  virtual void setString(T_GENE *aips_string) 
  {
    interfacesse::copy
      (this->_pts_string, 
       aips_string, 
       _stcui_stringSize
       );
  }
  
  virtual T_GENE* getString() 
  {
    return this->_pts_string;
  }

  virtual const T_GENE* getString() const
  {
    return this->_pts_string;
  }

  inline T_GENE* begin()
  {
    return this->_pts_string;
  }

  inline T_GENE* end()
  {
    return this->_pts_string + _stcui_stringSize;
  }
  
  virtual const T_GENE getGene(const uintidx aiuintidx_idxGene) const 
  {   
    assert(0 <= aiuintidx_idxGene && aiuintidx_idxGene < _stcui_stringSize );
    return this->_pts_string[aiuintidx_idxGene];
  }

  virtual void setGene(uintidx aiuintidx_idxGene, T_GENE ait_newGene) 
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "gaencode::ChromFixedLength::setGene";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input gaencode::ChromFixedLength<>: this[" <<  this << ']'
		<< "\n input aiuintidx_idxGene = " << aiuintidx_idxGene
		<< "\tait_newGene = " << ait_newGene 
		<< "\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    assert(0 <= aiuintidx_idxGene && aiuintidx_idxGene < _stcui_stringSize );
    this->_pts_string[aiuintidx_idxGene] = ait_newGene;

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ')';
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }
   
  virtual void  print
  (std::ostream &os=std::cout,
   const char* aipc_label   = "",
   const char aic_delimCoef = ',',
   const char aic_delimRow  = ';'
   ) const 
  {

#if defined(__VERBOSE_YES)
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    os   << '>';
#else
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    os  << ':';
#endif
    
    if ( this->_pts_string != NULL ) {
      for(uintidx li_j = 0; li_j < ( this->_stcui_stringSize - 1 ); li_j++) {
	os << this->_pts_string[li_j] << aic_delimCoef;
      }
      os << this->_pts_string[(this->_stcui_stringSize-1)];
    }
  }

protected:

  T_GENE         *_pts_string;
  static uintidx _stcui_stringSize;

}; //End ChromFixedLength

template <class T_GENE, class T_METRIC>
uintidx gaencode::ChromFixedLength<T_GENE,T_METRIC>::_stcui_stringSize = 0;


} /*END namespace gaencode*/

#endif  /*CHROMOSOME_FIXEDLENGTH_HPP*/
