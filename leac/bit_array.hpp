/*! \file bit_array.hpp
 *
 * \brief bit array  
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __BIT_ARRAY_HPP
#define __BIT_ARRAY_HPP

#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include "bit_container.hpp"
#include "common.hpp"

#include <stdio.h>

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {

/*! \class BitArray
  \brief Bit array
  \details 
*/ 
template < class T_BITSIZE >
class BitArray: public ds::BitContainer<T_BITSIZE> {
public:
  BitArray():
    ds::BitContainer<T_BITSIZE>(),
    _uintidx_numBits(0)
  {}

  BitArray(const uintidx aiintidx_numBits):
    ds::BitContainer<T_BITSIZE>(ds::BitContainer<T_BITSIZE>::_getNumWords(aiintidx_numBits)),
  _uintidx_numBits(aiintidx_numBits)
  {}
  
  BitArray(const uintidx aiintidx_numBits, T_BITSIZE* aiarrayT_data):
    ds::BitContainer<T_BITSIZE>
    (ds::BitContainer<T_BITSIZE>::_getNumWords(aiintidx_numBits),
     aiarrayT_data),
    _uintidx_numBits(aiintidx_numBits)
  {}

  BitArray(const BitArray<T_BITSIZE>& aibitarray_B):
  ds::BitContainer<T_BITSIZE>(aibitarray_B),
    _uintidx_numBits(aibitarray_B._uintidx_numBits)
  {}

  //move constructor
  BitArray(const BitArray<T_BITSIZE>&& aibitarray_B):
  ds::BitContainer<T_BITSIZE>(aibitarray_B),
    _uintidx_numBits(aibitarray_B._uintidx_numBits)
  {
    aibitarray_B._uintidx_numBits = 0;
  }

  BitArray<T_BITSIZE>& operator=(BitArray<T_BITSIZE> &&aibitarray_B)
  {
    if ( this !=  &aibitarray_B ){
      ds::BitContainer<T_BITSIZE>::operator=(aibitarray_B);
      _uintidx_numBits = aibitarray_B._uintidx_numBits;
      aibitarray_B._uintidx_numBits  = 0;
    }

    return *this;
  }

  BitArray<T_BITSIZE>& operator=(const BitArray<T_BITSIZE> &aibitarray_B)
  {
    if( this !=  &aibitarray_B ){
      ds::BitContainer<T_BITSIZE>::operator=(aibitarray_B);
      _uintidx_numBits = aibitarray_B._uintidx_numBits;
    }

    return *this;
  }

  ~BitArray() 
  {
  }

  /*initializeOn:
     set all elements of data to one
  */
  inline void initializeOn()
  {
    ds::BitContainer<T_BITSIZE>::setAll();
    _maskTopWord();
  }

  void toggleAll()
  {
    ds::BitContainer<T_BITSIZE>::toggleAll();
    _maskTopWord();
  }

  
  inline 
  uintidx size() const
  {
    return _uintidx_numBits; 
  }

  //const char delim='') const
  void print(std::ostream &os=std::cout, const char *aipc_label = "") const 
  {

#if defined(__VERBOSE_YES)
    int li_tmpVerboseMax =  geiinparam_verboseMax; // desactive verbose: getNumBitOn
    geiinparam_verboseMax = -1;
    os << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << this << ']'
       << ":length," << _uintidx_numBits << ",biton," << this->getNumBitOn() << '>';
    geiinparam_verboseMax = li_tmpVerboseMax;
#else
      os << aipc_label;
#endif
   
      //os << aipc_label << ":length," << _uintidx_numBits << '>';
    for(uintidx i = 0; i < this->_uintidx_numBits; i++) {
      os << std::noboolalpha << this->getBit(i);
    }
    
  }
  
  friend std::ostream& operator<<(std::ostream& os, const BitArray<T_BITSIZE>& aibitarray_B) 
  {
    aibitarray_B.print(os);
    
    return os;
  }

private:

  void _maskTopWord()
  {
    uintidx luintidx_bitsActive = this->_bits_in_top_word(this->_uintidx_numBits);
    if ( luintidx_bitsActive > 0 )  {
      this->_arrayT_data[this->_uintidx_numWords-1] &= this->_bitMask(luintidx_bitsActive);
    }
  }


protected:
  
  /*_bit_array_get_word:
   Get and set words (internal use only -- no bounds checking)
  */
  uintidx        _uintidx_numBits;
  
} /* BitArray */;


} /*END namespace mat*/


#endif /*__BIT_ARRAY_HPP*/
