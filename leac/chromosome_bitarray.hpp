/*! \file chromosome_bitarray.hpp
 *
 * \brief chromosome bitarray
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CHROMOSOME_BITARRAY_HPP
#define CHROMOSOME_BITARRAY_HPP

#include "chromosome_base.hpp"
#include "bit_array.hpp"

#include "verbose_global.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {


/*! \class ChromosomeBitArray
  \brief Chromosome encoding string of bits
  \details 
*/ 
template < class T_BITSIZE,
	   class T_METRIC
	   >
class ChromosomeBitArray
  : public ChromosomeBase<T_METRIC> 
  , public mat::BitArray<T_BITSIZE> 
{
public:
  ChromosomeBitArray() 
    : ChromosomeBase<T_METRIC>()
    , mat::BitArray<T_BITSIZE>()  
    {}
  
  ChromosomeBitArray
  (const uintidx aiintidx_numBits)
    : ChromosomeBase<T_METRIC>()
    , mat::BitArray<T_BITSIZE>(aiintidx_numBits)
  { }

  //copy constructor
  ChromosomeBitArray
  (const ChromosomeBitArray<T_BITSIZE,T_METRIC> &aichrombitarray_b)
    : ChromosomeBase<T_METRIC>(aichrombitarray_b)
    , mat::BitArray<T_BITSIZE>(aichrombitarray_b)
  {

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ChromosomeBitArray::ChromosomeBitArray:  IN"
		<< '(' << geiinparam_verbose << ')'
	        << "\n\t input  ChromosomeBitArray<>: &aichrombitarray_b[" 
		<<  &aichrombitarray_b << "]\n"
		<< "\t)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ChromosomeBitArray::ChromosomeBitArray: OUT"
		<< '(' << geiinparam_verbose << ")\n"
	        << *this
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  //move constructor
  ChromosomeBitArray
  (ChromosomeBitArray<T_BITSIZE,T_METRIC> &&aichrombitarray_b)
    : ChromosomeBase<T_METRIC>(aichrombitarray_b)
    , mat::BitArray<T_BITSIZE>(aichrombitarray_b)
  {
    
  }

  virtual ~ChromosomeBitArray()
  {
    
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeBitArray::ChromosomeBitArray~";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ':' << geverbosepc_labelstep
		<< '[' << geverboseui_idproc << ':' << this << ']'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }
  
 
  ChromosomeBitArray<T_BITSIZE,T_METRIC>& 
   operator=(const ChromosomeBitArray<T_BITSIZE,T_METRIC> &aichrombitarray_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeBitArray::operator=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeBitArray<>: &&aichrombitarray_b[" <<  &aichrombitarray_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichrombitarray_b  ) {
      ChromosomeBase<T_METRIC>::operator=(aichrombitarray_b);
      mat::BitArray<T_BITSIZE>::operator=(aichrombitarray_b);
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeBitArray<T_BITSIZE,T_METRIC>::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    return *this;
  }

  ChromosomeBitArray<T_BITSIZE,T_METRIC>&
   operator=(ChromosomeBitArray<T_BITSIZE,T_METRIC> &&aichrombitarray_b)
  {
    if ( this != &aichrombitarray_b  ) {
      ChromosomeBase<T_METRIC>::operator=(aichrombitarray_b);
      mat::BitArray<T_BITSIZE>::operator=(aichrombitarray_b);
    }

    return *this;
  }

  
  virtual void  print
  (std::ostream &os=std::cout,
   const char*  aipc_label    = "",
   const char   aic_delimCoef = ',',
   const char   aic_delimRow  = ';'
   ) const 
  {
    
#if defined(__VERBOSE_YES)
    std::string  lpc_labelChrom("BITARRAY");
    ChromosomeBase<T_METRIC>::print(os,aipc_label);
    os << ',';
    mat::BitArray<T_BITSIZE>::print(os,lpc_labelChrom.c_str());
#else
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    os << ',';
    mat::BitArray<T_BITSIZE>::print(os);
#endif
    
  }

  friend std::ostream& operator<<
  (std::ostream& os, 
   const ChromosomeBitArray<T_BITSIZE,T_METRIC> &aichrombitarray_b
   )
  {
    aichrombitarray_b.print(os,""); 
    
    return os;
  }

protected:
 
}; /*ChromosomeBitArray*/

} /*END namespace gaencode*/

#endif  /*CHROMOSOME_BITARRAY_HPP*/
