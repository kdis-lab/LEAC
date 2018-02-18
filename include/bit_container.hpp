/*! \file bit_container.hpp
 *
 * \brief  bit container
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef BIT_CONTAINER_HPP
#define BIT_CONTAINER_HPP

#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <list>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include "common.hpp"

#include "verbose_global.hpp"

namespace mat {
  template<class T_BITSIZE> class BitMatrix;
}


/*! \namespace ds
  \brief Data structure
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace ds {
  
/*! \class BitContainer
  \brief Bit container
  \details 
*/ 
template < class T_BITSIZE >
class BitContainer {
public:

  static const unsigned short int _stusi_numBitsWords = 
    (unsigned short int) (sizeof(T_BITSIZE) * 8);
  static const T_BITSIZE  t_wordMax = (~(T_BITSIZE)0);

  BitContainer()
  : _uintidx_numWords(0)
  , _arrayuit_data(NULL)
  ,  b_externalData(false)
  { }

  BitContainer(const uintidx  aiuintidx_numWords)
  : _uintidx_numWords(aiuintidx_numWords)
  , b_externalData(false)
  {
    _arrayuit_data = (T_BITSIZE*) calloc(this->_uintidx_numWords, sizeof(T_BITSIZE));
    if(_arrayuit_data == NULL) {
      throw  std::range_error
	("BitContainer::BitContainer(uintidx):" 
	 "error - could not allocate enough memory");
    }
    this->initialize();
  }

  BitContainer(const uintidx  aiuintidx_numWords, T_BITSIZE* aiarrayit_data):
    _uintidx_numWords(aiuintidx_numWords),
    _arrayuit_data(aiarrayit_data),
    b_externalData(true)
  {
    
  }
 
  BitContainer(const ds::BitContainer<T_BITSIZE>& aibitcontainer_B):
    _uintidx_numWords(aibitcontainer_B._uintidx_numWords),
    b_externalData(false)
  {
    _arrayuit_data = (T_BITSIZE*) calloc(this->_uintidx_numWords, sizeof(T_BITSIZE));
    if(_arrayuit_data == NULL) {
      throw  std::range_error
	("BitContainer::BitContainer(const BitContainer&):" 
	 "error - could not allocate enough memory");
    }
    memcpy
      (_arrayuit_data,
       aibitcontainer_B._arrayuit_data, 
       aibitcontainer_B._uintidx_numWords * sizeof(T_BITSIZE)
       );    
  }

  //move constructor
  BitContainer(ds::BitContainer<T_BITSIZE>&& aibitcontainer_B)
    :  _uintidx_numWords(aibitcontainer_B._uintidx_numWords)
    ,  _arrayuit_data(aibitcontainer_B._arrayuit_data)
    ,  b_externalData(aibitcontainer_B.b_externalData)
  {
    aibitcontainer_B._uintidx_numWords    = 0;
    aibitcontainer_B._arrayuit_data   = NULL;
    aibitcontainer_B.b_externalData = true;
  }

  virtual ~BitContainer() 
  {
    if( _arrayuit_data != NULL && this->b_externalData == false )
      free(_arrayuit_data);
  }


  ds::BitContainer<T_BITSIZE>& operator=(const ds::BitContainer<T_BITSIZE>& aibitcontainer_B)
  {
    if ( this != &aibitcontainer_B ) {
    
      if ( _uintidx_numWords != aibitcontainer_B._uintidx_numWords ) {
	if ( _arrayuit_data != NULL ) free(_arrayuit_data);
	_uintidx_numWords = aibitcontainer_B._uintidx_numWords;
        _arrayuit_data = (T_BITSIZE*) calloc(this->_uintidx_numWords, sizeof(T_BITSIZE));
      }
    
      if(_arrayuit_data == NULL) {
	throw  std::range_error
	  ("BitContainer::BitContainer(uintidx):" 
	   "error - could not allocate enough memory");
      }
      memcpy
	(_arrayuit_data,
	 aibitcontainer_B._arrayuit_data, 
	 aibitcontainer_B._uintidx_numWords * sizeof(T_BITSIZE)
	 );
    }
    return *this;
  }

  ds::BitContainer<T_BITSIZE>& operator=(ds::BitContainer<T_BITSIZE>&& aibitcontainer_B)
  {
    if ( this != &aibitcontainer_B ) {
      if ( _arrayuit_data != NULL )
	free(_arrayuit_data);
      _uintidx_numWords    = aibitcontainer_B._uintidx_numWords;
      _arrayuit_data   = aibitcontainer_B._arrayuit_data;
      b_externalData = aibitcontainer_B.b_externalData;
      
      aibitcontainer_B._uintidx_numWords    = 0;
      aibitcontainer_B._arrayuit_data   = NULL;
      aibitcontainer_B.b_externalData = true;
    }
    return *this;
  }

  inline void setArray(T_BITSIZE* aiarrayit_data)
  {
    _arrayuit_data = aiarrayit_data;
  }

  inline T_BITSIZE* toArray()
  { 
    return _arrayuit_data; 
  }

  //! initialize
  /*!
    set all elements of data to zero
  */
  inline void initialize()
  {
    memset(_arrayuit_data, 0, this->_uintidx_numWords * sizeof(T_BITSIZE));
  }

  //!initializeOn:
  /*! 
    set all elements of data to one
  */
  void initializeOn()
  {
    uintidx num_of_bytes = this->_uintidx_numWords * sizeof(T_BITSIZE);
    memset(_arrayuit_data, 0xFF, num_of_bytes);
  }

  //!toggleAll
  /*!
    Set all 1 bits to 0, and all 0 bits to 1. AKA flip
  */
  void toggleAll()
  { 
    for(uintidx i = 0; i < this->_uintidx_numWords; i++) {
      _arrayuit_data[i] ^= ds::BitContainer<T_BITSIZE>::t_wordMax; /*WORD_MAX;*/
    }    
  }

   /*clearBit:
    Set all 1 bits to 0, and all 0 bits to 1. AKA flip
  */
  inline void clearBit(uintidx i) 
  {
    _arrayuit_data[BitContainer<T_BITSIZE>::_getWordIndexBit(i)] &= 
      ~((T_BITSIZE)1 << BitContainer<T_BITSIZE>::_getBitOffset(i));
  }

  inline void setBit(uintidx i) 
  {
    _arrayuit_data[BitContainer<T_BITSIZE>::_getWordIndexBit(i)] |=  
      ((T_BITSIZE)1 << BitContainer<T_BITSIZE>::_getBitOffset(i));
  }

  inline void toggleBit(uintidx i) 
  {
    _arrayuit_data[BitContainer<T_BITSIZE>::_getWordIndexBit(i)] ^=  
      (T_BITSIZE)1 << BitContainer<T_BITSIZE>::_getBitOffset(i);
  }

  inline void assignBit(uintidx i, bool aib_bit)
  {
    ((aib_bit)?setBit(i):clearBit(i)); 
  }

  inline bool getBit(uintidx i) const 
  {
    return (this->_arrayuit_data[BitContainer<T_BITSIZE>::_getWordIndexBit(i)] >>  
	    (BitContainer<T_BITSIZE>::_getBitOffset(i)) & 0x1);
  }
  
  uintidx getNumBitOn() const 
  {
    uintidx    louintidx_numBitOn;
    T_BITSIZE lt_word;

    louintidx_numBitOn = 0;

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ds::BitContainer<T_BITSIZE>::getNumBitOn  IN" 
		<< '(' << geiinparam_verbose << ')'
		<< "\n( input  ds::BitContainer<T_BITSIZE>: this[" << this << "]\n"
		<< ")\n";
    }
#endif /*__VERBOSE_YES*/

    for(uintidx i = 0; i < this->_uintidx_numWords; i++) {
      lt_word = _arrayuit_data[i];
      for(unsigned short int j = 0; j < ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords; j++) {
	if ( lt_word & 0x1) ++louintidx_numBitOn;  
	lt_word >>= 1;
      }
    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ds::BitContainer<T_BITSIZE>::getNumBitOn OUT" 
		<< '(' << geiinparam_verbose << ')'
		<< "\tuintidx: louintidx_numBitOn = " << louintidx_numBitOn
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return louintidx_numBitOn;
  }


  std::list<uintidx> getIdxWithBitOn() const 
  {
    std::list<uintidx> lostui_idxInstance;
    T_BITSIZE lt_word;
    
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ds::BitContainer<T_BITSIZE>::getIdxWithBitOn()";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
   		<< "( input  ds::BitContainer<T_BITSIZE>: this[" << this << "]\n"
		<< ")\n";
    }
#endif /*__VERBOSE_YES*/

    for(uintidx i = 0; i < this->_uintidx_numWords; i++) {
      lt_word = _arrayuit_data[i];
      for(unsigned short int j = 0; j < ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords; j++) {
	if ( lt_word & 0x1 )
	  lostui_idxInstance.push_back((uintidx)(ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords * i + j));
	lt_word >>= 1;
      }
    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      inout::containerprint
	(lostui_idxInstance.begin(),
	 lostui_idxInstance.end(),
	 std::cout,
	 lpc_labelFunc,
	 ','
	 );
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return lostui_idxInstance;
  }
 

  //! opXor
  /*! 
    xor
  */
  void opXor(const ds::BitContainer<T_BITSIZE>   &aibitc_bits)
  { 
    for(uintidx i = 0; i < this->_uintidx_numWords; i++) {
      _arrayuit_data[i] ^= aibitc_bits._arrayuit_data[i]; 
    }    
  }
  
  void copyAligned
  (const ds::BitContainer<T_BITSIZE> &aibitc_source, 
   const uintidx                     aiuintidx_idxTarget, 
   uintidx                           aiuintidx_length
   )
  {
    const ds::BitContainer<T_BITSIZE> *lbitc_tmpSource;
 
#ifdef __VERBOSE_YES
    uintidx  luintidx_length = aiuintidx_length;
    const char* lpc_labelFunc = "ds::BitContainer<T_BITSIZE>::copyAligned"; 
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n";
       std::ostringstream lostrstream_aibitc_this;
       lostrstream_aibitc_this << lpc_labelFunc << ":aibitc_this";
       ds::BitContainer<T_BITSIZE>::print(std::cout,lostrstream_aibitc_this.str().c_str());
       std::cout <<  '\n';
       std::ostringstream lostrstream_aibitc_source;
       lostrstream_aibitc_source << lpc_labelFunc << ":aibitc_source";
       aibitc_source.print(std::cout,lostrstream_aibitc_source.str().c_str());
       std::cout <<  '\n';

       std::cout << " input  uintidx: aiuintidx_idxTarget = " << aiuintidx_idxTarget << '\n'
		 << " input  uintidx: aiuintidx_length = " << aiuintidx_length << '\n'
		 << ")\n";
    }
#endif /*__VERBOSE_YES*/

    if ( this == &aibitc_source )  {
      lbitc_tmpSource = new ds::BitContainer<T_BITSIZE>(aibitc_source);
    }
    else {
      lbitc_tmpSource = &aibitc_source;
    }
    /*THE PART FIRST*/
    
    unsigned short int lusi_boffset = this->_getBitOffset(aiuintidx_idxTarget);
    uintidx luintidx_wordIndexBegin = this->_getWordIndexBit(aiuintidx_idxTarget);
    
    if ( lusi_boffset != 0 ) {
      unsigned short int lusi_maxBitsCopyBoffset =
	ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords - lusi_boffset;
      unsigned short int lusi_maxBitsCopy = ( aiuintidx_length >  lusi_maxBitsCopyBoffset)?
	lusi_maxBitsCopyBoffset: (unsigned short int) aiuintidx_length;
      _arrayuit_data[luintidx_wordIndexBegin] = 
	this->_setWordAligned
	(_arrayuit_data[luintidx_wordIndexBegin],
	 lbitc_tmpSource->_arrayuit_data[luintidx_wordIndexBegin], 
	 lusi_boffset, 
	 lusi_maxBitsCopy);
      aiuintidx_length -=  lusi_maxBitsCopy;
      ++luintidx_wordIndexBegin;
    }
    uintidx luintidx_numFullWords = this->_getWordIndexBit(aiuintidx_length);
    
    if ( luintidx_numFullWords > 0  ) {

      memcpy
	(&_arrayuit_data[luintidx_wordIndexBegin], 
	 &lbitc_tmpSource->_arrayuit_data[luintidx_wordIndexBegin], 
	 luintidx_numFullWords * sizeof(T_BITSIZE)
	 ); 
      aiuintidx_length -= (luintidx_numFullWords *  ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords);
      luintidx_wordIndexBegin += luintidx_numFullWords;
    } 	   
    
    if ( aiuintidx_length > 0) {

      _arrayuit_data[luintidx_wordIndexBegin] = 
	this->_setWordAligned
	(_arrayuit_data[luintidx_wordIndexBegin],
	 lbitc_tmpSource->_arrayuit_data[luintidx_wordIndexBegin], 
	 0, 
	 aiuintidx_length);
      
    }

    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::string lstr_termination = "OK";
      for ( uintidx  lui_j = aiuintidx_idxTarget; lui_j < (aiuintidx_idxTarget+luintidx_length);lui_j++) {
	if ( this->getBit(lui_j) != aibitc_source.getBit(lui_j) ) { 
	  lstr_termination = "FAILURE";
	  break;
	}
      }
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ')' << " Termination:" << lstr_termination << '\n';
      ds::BitContainer<T_BITSIZE>::print(std::cout,lpc_labelFunc);
      std::cout <<  std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    if ( this == &aibitc_source )  {
      delete lbitc_tmpSource;
    }
    
  }

  void copyUnaligned
  (uintidx         aiuintidx_idxTarget, 
   BitContainer
   <T_BITSIZE>    &aibitc_source, 
   uintidx         aiuintidx_idxSource, 
   uintidx         aiuintidx_length
   )
  {
    uintidx                 luintidx_idxWordSrc;
    uintidx                 luintidx_idxWordTarget;
    unsigned short int      lusi_idxBitTarget;
    unsigned short int      lusi_idxBitWordSrc;
    unsigned short int      lusi_maxBitsCopy;
    ds::BitContainer<T_BITSIZE> *lbitc_tmpSource;

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ds::BitContainer<T_BITSIZE>::copyUnaligned  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\t(output ds::BitContainer<T_BITSIZE>: this[" << this << "]\n"
		<< "\t input  uintidx: aiuintidx_idxTarget = " << aiuintidx_idxTarget << '\n'
		<< "\t input  ds::BitContainer<T_BITSIZE>: &aibitc_source[" << &aibitc_source << "]\n"
		<< "\t input  uintidx: aiuintidx_idxSource = " << aiuintidx_idxSource << '\n'
		<< "\t input  uintidx: aiuintidx_length = " << aiuintidx_length << '\n'
		<< "\t)\n";
    }
#endif /*__VERBOSE_YES*/

    if ( this == &aibitc_source )  {
      lbitc_tmpSource = new ds::BitContainer<T_BITSIZE>(aibitc_source);
    }
    else {
      lbitc_tmpSource = &aibitc_source;
    }

    while ( aiuintidx_length != 0) {
      luintidx_idxWordSrc  = this->_getWordIndexBit(aiuintidx_idxSource);
      lusi_idxBitWordSrc = _getBitOffset(aiuintidx_idxSource);
      luintidx_idxWordTarget  = this->_getWordIndexBit(aiuintidx_idxTarget);
      lusi_idxBitTarget   = _getBitOffset(aiuintidx_idxTarget);

      lusi_maxBitsCopy = 
	std::min
	(ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords - lusi_idxBitTarget,
	 ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords - lusi_idxBitWordSrc);
      
      if ( aiuintidx_length < lusi_maxBitsCopy )
	lusi_maxBitsCopy = aiuintidx_length % ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords; 
    
      _arrayuit_data[luintidx_idxWordTarget] = 
	_setWordUnaligned
	(_arrayuit_data[luintidx_idxWordTarget],
	 lusi_idxBitTarget,
	 lbitc_tmpSource->_arrayuit_data[luintidx_idxWordSrc],
	 lusi_idxBitWordSrc,
	 lusi_maxBitsCopy
	 );
    
      aiuintidx_length    -= lusi_maxBitsCopy;
      aiuintidx_idxSource += lusi_maxBitsCopy;
      aiuintidx_idxTarget += lusi_maxBitsCopy;
    
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ds::BitContainer<T_BITSIZE>::copyUnaligned: OUT"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\toutput ds::BitContainer<T_BITSIZE>: this[" << this << "]\n";
      ds::BitContainer<T_BITSIZE>::print();
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    if ( this == &aibitc_source )  {
      delete lbitc_tmpSource;
    }
     
  }

  //!_getNumWords
  /*!
    Number of words required to store a given number of bits
    0 -> 0
    1..WORD_SIZE -> 1
    WORD_SIZE+1..2*WORD_SIZE -> 2 etc.
  */
  static const uintidx _getNumWords(uintidx aiuintidx_numBits)                
  { 
    return (aiuintidx_numBits + ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords - 1) 
      / ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords;
  }

  void  print(std::ostream &os=std::cout, const char   *aipc_label = "") const 
  {

#if defined(__VERBOSE_YES)
    os << "<BitContainer:"
       << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << this << ']'
       << ":_uintidx_numWords," << _uintidx_numWords << '>';
#else
    os << aipc_label;
#endif
    for(uintidx i = 0; i < this->_uintidx_numWords; i++) {
      this->print(_arrayuit_data[i],os);
    }
  }

private:

  //!_setWordAligned
  /*!
    common internal functions
  */
  T_BITSIZE _setWordAligned
  (T_BITSIZE aiT_wordTarget, 
   T_BITSIZE aiT_wordSource, 
   unsigned short int at, 
   unsigned short int numbits
   )
  {
  
    T_BITSIZE mask = 
      (ds::BitContainer<T_BITSIZE>::t_wordMax>>
       (ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords-numbits))<<at;
  
    return (aiT_wordTarget&~mask)|((aiT_wordSource)&mask);
  }

  T_BITSIZE _setWordUnaligned
  (T_BITSIZE             aiT_wordTarget,
   unsigned short int    dstindx,
   T_BITSIZE             aiT_wordSource,
   unsigned short int    aiuintidx_idxSource,
   unsigned short int    numbits
   )
  {
  
    if ( dstindx < aiuintidx_idxSource)
      aiT_wordSource >>= (aiuintidx_idxSource- dstindx);
    else if ( aiuintidx_idxSource < dstindx)
      aiT_wordSource <<= (dstindx - aiuintidx_idxSource);
  
    T_BITSIZE mask = (ds::BitContainer<T_BITSIZE>::t_wordMax>>(ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords-numbits))<<dstindx;
  
    return (aiT_wordTarget&~mask)|((aiT_wordSource)&mask);
  }

  
  void  print(T_BITSIZE ai_word, std::ostream &os=std::cout) const //, char delim='\t')
  {  
    for(unsigned short int i = 0; i < ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords; i++) {
      std::cout << (ai_word & 0x1);
      ai_word >>= 1;
    }
  }

protected:
 
  inline unsigned short int _getBitOffset(uintidx b) const
  { 
    return (unsigned short int) ((b) % ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords);
  }

  uintidx _bits_in_top_word(uintidx b)
  {
    return (b == 0 ? 0 : _getBitOffset(b - 1) + 1);
  }

  inline uintidx _getWordIndexBit(uintidx b) const 
  {
    return ((b) / ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords);
  }

  inline T_BITSIZE _bitMask(uintidx length)
  {                                    //WORD_MAX
    T_BITSIZE lt_mask = ds::BitContainer<T_BITSIZE>::t_wordMax;
    lt_mask  >>= (ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords-(length));

    return lt_mask;
  }

  uintidx    _uintidx_numWords;
  T_BITSIZE  *_arrayuit_data;
  bool       b_externalData;

  friend class  mat::BitMatrix<T_BITSIZE>; 

}; /*BitContainer*/

} /*END namespace ds*/

#endif /*BIT_CONTAINER_HPP*/
