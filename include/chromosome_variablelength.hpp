/*! \file chromosome_variablelength.hpp
 *
 * \brief chromosome variable length
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef CHROMOSOME_VARIABLELENGTH_HPP
#define CHROMOSOME_VARIABLELENGTH_HPP

#include <iostream>
#include <stdexcept>
#include <limits>
#include <assert.h>        
#include "outfilename.hpp" //outparam::OutFileName::getDelim()
#include "chromosome_string.hpp"


/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {
  

/*! \class ChromVariableLength
  \brief Chromosome with variable length string 
*/
template <class T_GENE,
	  class T_METRIC //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  >
class ChromVariableLength: 
  public ChromosomeString<T_GENE,T_METRIC>
{   
public:
  ChromVariableLength()
    :  ChromosomeString<T_GENE,T_METRIC>()
    ,  _pts_string(NULL)
    ,  _uintidx_stringSize((uintidx) 0)
  { }

  ChromVariableLength(const T_METRIC airt_objetiveFunc, const T_METRIC airt_fitness)
    :  ChromosomeString<T_GENE,T_METRIC>(airt_objetiveFunc, airt_fitness)
    ,  _pts_string(NULL)
    ,  _uintidx_stringSize((uintidx) 0)
  { }
  
  ChromVariableLength(const uintidx aiuintidx_stringSize)
  : ChromosomeString<T_GENE,T_METRIC>()
  , _pts_string(new T_GENE[aiuintidx_stringSize])
  , _uintidx_stringSize(aiuintidx_stringSize)
  { }
  
  //move constructor 
  ChromVariableLength(ChromVariableLength<T_GENE,T_METRIC> &&aichromvarlength_b)
    :  ChromosomeString<T_GENE,T_METRIC>(aichromvarlength_b)
    ,  _pts_string(aichromvarlength_b._pts_string)
    ,  _uintidx_stringSize(aichromvarlength_b._uintidx_stringSize)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromVariableLength::Move";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromVariableLength<>: &&aichromvarlength_b[" <<  &aichromvarlength_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromVariableLength<T_GENE,T_METRIC>::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    aichromvarlength_b._pts_string = NULL;
    aichromvarlength_b._uintidx_stringSize = 0;
    
  }

  //copy constructor
  ChromVariableLength
  (const ChromVariableLength<T_GENE,T_METRIC> &aichromvarlength_b)
    :  ChromosomeString<T_GENE,T_METRIC>(aichromvarlength_b)
    , _pts_string(new T_GENE[aichromvarlength_b._uintidx_stringSize])
    , _uintidx_stringSize(aichromvarlength_b._uintidx_stringSize)
  {
    
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "ChromVariableLength::copy";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( aichromvarlength_b ["  
	      << &aichromvarlength_b << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
    if (_pts_string != NULL) {
      interfacesse::copy
	(this->_pts_string, aichromvarlength_b._pts_string, this->_uintidx_stringSize);
    }
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    ChromVariableLength<T_GENE,T_METRIC>::print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

    
  }

  virtual ~ChromVariableLength() 
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromVariableLength::ChromVariableLength~";
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

  void resize(const uintidx aiui_newStringSize)
  {
    if ( _uintidx_stringSize != aiui_newStringSize )  {
	if ( this->_pts_string != NULL )  delete[] _pts_string;
	_uintidx_stringSize = aiui_newStringSize;
	_pts_string = new T_GENE[_uintidx_stringSize];
    }
  }
 
  ChromVariableLength<T_GENE,T_METRIC>& 
  operator=(const ChromVariableLength<T_GENE,T_METRIC> &aichromvarlength_b)
  { 
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromVariableLength::operator=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromVariableLength<>: &&aichromvarlength_b[" <<  &aichromvarlength_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichromvarlength_b ) { 
      ChromosomeString<T_GENE,T_METRIC>::operator=(aichromvarlength_b);
      if ( _uintidx_stringSize != aichromvarlength_b._uintidx_stringSize )  {
	if ( this->_pts_string != NULL )  delete[] _pts_string;
	_uintidx_stringSize = aichromvarlength_b._uintidx_stringSize;
	_pts_string = new T_GENE[_uintidx_stringSize];
      }
      interfacesse::copy
	(this->_pts_string, 
	 aichromvarlength_b._pts_string, 
	 this->_uintidx_stringSize
	 );
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromVariableLength<T_GENE,T_METRIC>::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    
    return *this;
  }

   ChromVariableLength<T_GENE,T_METRIC>& 
   operator=(ChromVariableLength<T_GENE,T_METRIC> &&aichromvarlength_b)
  {
    if ( this != &aichromvarlength_b ) {
      ChromosomeString<T_GENE,T_METRIC>::operator=(aichromvarlength_b);
      if ( this->_pts_string != NULL ) {
	delete[] _pts_string; 
      }
      this->_pts_string = aichromvarlength_b._pts_string;
      this->_uintidx_stringSize = aichromvarlength_b._uintidx_stringSize;
 
      aichromvarlength_b._uintidx_stringSize = 0;
      aichromvarlength_b._pts_string = NULL;
    }

    return *this;
  }

  virtual const uintidx getStringSize() const 
  {
    return this->_uintidx_stringSize;
  }

  virtual void setString(T_GENE *aips_string) 
  {
    interfacesse::copy
      (this->_pts_string, 
       aips_string, 
       this->_uintidx_stringSize
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

  virtual const T_GENE getGene(const uintidx aiuintidx_idxGene) const 
  {
    assert(0 <= aiuintidx_idxGene && aiuintidx_idxGene < _uintidx_stringSize );
    return this->_pts_string[aiuintidx_idxGene];
  }

  virtual void setGene(uintidx aiuintidx_idxGene, T_GENE ait_newGene) 
  {
    assert(0 <= aiuintidx_idxGene && aiuintidx_idxGene < _uintidx_stringSize );
    this->_pts_string[aiuintidx_idxGene] = ait_newGene;
  }

  inline T_GENE* begin()
  {
    return this->_pts_string;
  }

  inline T_GENE* end()
  {
    return this->_pts_string + _uintidx_stringSize;
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
    os << ",length," << _uintidx_stringSize
       << '>';
#else
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    os << ",length," << _uintidx_stringSize
       << '>';
#endif
    
    if ( this->_pts_string != NULL ) {
      for(uintidx li_j = 0; li_j < ( this->_uintidx_stringSize - 1 ); li_j++) {
	os << this->_pts_string[li_j] << aic_delimCoef;
      }
      os << this->_pts_string[(this->_uintidx_stringSize-1)];
    }
  }

protected:

  T_GENE   *_pts_string;
  uintidx     _uintidx_stringSize;

}; //End ChromVariableLength

} /*END namespace gaencode*/

#endif  /*CHROMOSOME_VARIABLELENGTH_HPP*/
