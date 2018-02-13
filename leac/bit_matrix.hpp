/*! \file bit_matrix.hpp
 *
 * \brief bit matrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> licenset
 */


#ifndef BIT_MATRIX_HPP
#define BIT_MATRIX_HPP

#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include "bit_container.hpp"
#include "bit_array.hpp"
#include "vector_utils.hpp"
#include <stdio.h>

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {
  
/*! \class BitMatrix
  \brief Bit matrix
  \details 
*/ 
template < class T_BITSIZE >
class BitMatrix { 
public:
  BitMatrix() 
    : _bitContainer_data()
    , _ui_numRows(0) 
    , _ui_numColumns(0)
    , _st_numWordsRows( ds::BitContainer<T_BITSIZE>::_getNumWords(0) )
  {}

  BitMatrix
  (const uintidx aiui_numRows, 
   const uintidx aiui_numColumns
   ): 
    _bitContainer_data
    (aiui_numRows * ds::BitContainer<T_BITSIZE>::_getNumWords(aiui_numColumns)),
    _ui_numRows(aiui_numRows), 
    _ui_numColumns(aiui_numColumns),
    _st_numWordsRows( ds::BitContainer<T_BITSIZE>::_getNumWords(aiui_numColumns) )
  { }
  
  BitMatrix
  (const uintidx aiui_numRows, 
   const uintidx aiui_numColumns,
   T_BITSIZE*   aiarrayT_data
   ): 
    _bitContainer_data
    (aiui_numRows * ds::BitContainer<T_BITSIZE>::_getNumWords(aiui_numColumns),aiarrayT_data),
    _ui_numRows(aiui_numRows), 
    _ui_numColumns(aiui_numColumns),
    _st_numWordsRows( ds::BitContainer<T_BITSIZE>::_getNumWords(aiui_numColumns) )
  { }

  BitMatrix(const BitMatrix<T_BITSIZE>& aibitmatrixb_B) :
    _bitContainer_data(aibitmatrixb_B._bitContainer_data),
    _ui_numRows(aibitmatrixb_B._ui_numRows),
    _ui_numColumns(aibitmatrixb_B._ui_numColumns), 
    _st_numWordsRows( aibitmatrixb_B._st_numWordsRows )
  { }

  //move constructor
  BitMatrix(BitMatrix<T_BITSIZE> &&aibitmatrixb_B) :
    _bitContainer_data(aibitmatrixb_B._bitContainer_data),
    _ui_numRows(aibitmatrixb_B._ui_numRows),
    _ui_numColumns(aibitmatrixb_B._ui_numColumns), 
    _st_numWordsRows( aibitmatrixb_B._st_numWordsRows )
  {
    aibitmatrixb_B._ui_numRows = 0;
    aibitmatrixb_B._ui_numColumns = 0; 
    aibitmatrixb_B._st_numWordsRows = 0;
  }

  BitMatrix<T_BITSIZE>& operator=(BitMatrix<T_BITSIZE> &&aibitmatrixb_B)
  {
    if( this !=  &aibitmatrixb_B ){
      _bitContainer_data = aibitmatrixb_B._bitContainer_data;
      _ui_numRows = aibitmatrixb_B._ui_numRows;
      _ui_numColumns = aibitmatrixb_B._ui_numColumns;
      _st_numWordsRows= aibitmatrixb_B._st_numWordsRows;
      
      aibitmatrixb_B._ui_numRows  = 0;
      aibitmatrixb_B._ui_numColumns = 0; 
      aibitmatrixb_B._st_numWordsRows = 0;
    }

    return *this;
  }


  BitMatrix<T_BITSIZE>& operator=(const BitMatrix<T_BITSIZE> &aibitmatrixb_B)
  {
    if( this !=  &aibitmatrixb_B ){
       _bitContainer_data = aibitmatrixb_B._bitContainer_data;
      _ui_numRows = aibitmatrixb_B._ui_numRows;
      _ui_numColumns = aibitmatrixb_B._ui_numColumns;
      _st_numWordsRows= aibitmatrixb_B._st_numWordsRows;
      /*aibitmatrixb_B._ui_numRows  = 0;
      aibitmatrixb_B._ui_numColumns = 0; 
      aibitmatrixb_B._st_numWordsRows = 0;*/
    }

    return *this;
  }

  ~BitMatrix() { }

  inline void  initialize()
  {
    _bitContainer_data.initialize();
  }
  
  inline const uintidx getNumRows()    const 
  { 
    return this->_ui_numRows; 
  }

  inline const uintidx getNumColumns() const 
  { 
    return this->_ui_numColumns; 
  }
  
  inline const uintidx getNumElems()   const 
  { 
    return (this->_ui_numRows)*(this->_ui_numColumns); 
  }

  /*setAll:
     set all elements of data to one
  */
  inline void setAll()
  {
    _bitContainer_data.initializeOn();
    _maskTopWord();
  }

  void toggleAll()
  {
    _bitContainer_data.toggleAll();
    _maskTopWord();
  }

  inline void clearBit(uintidx i, uintidx j) 
  {
    _bitContainer_data._arrayuit_data[(i)*_st_numWordsRows+ _bitContainer_data._getWordIndexBit(j) ] &= 
      ~((T_BITSIZE)1 << _bitContainer_data._getBitOffset(j));
  }

  inline void setBit(uintidx i, uintidx j) 
  {
    _bitContainer_data._arrayuit_data[(i)*_st_numWordsRows+ _bitContainer_data._getWordIndexBit(j) ] |=  
      ((T_BITSIZE)1 << _bitContainer_data._getBitOffset(j) );
  }

  inline const T_BITSIZE* getRow (const uintidx i) const  
  { 
    return &_bitContainer_data._arrayuit_data[i*_st_numWordsRows]; 
  }

  inline T_BITSIZE* getRow ( const  uintidx i)   
  { 
    return &_bitContainer_data._arrayuit_data[i*_st_numWordsRows]; 
    }

  inline bool  operator() (uintidx i, uintidx j) const
  { 
    return ((_bitContainer_data._arrayuit_data
	     [(i) * _st_numWordsRows+(j)/ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords ] 
	     >> ((j) % ds::BitContainer<T_BITSIZE>::_stusi_numBitsWords)) & 0x1);
  }

  void copyAligned
  (BitMatrix<T_BITSIZE>  &aibitmatrix_source, 
   const uintidx               aiui_idxTarget, 
   const uintidx               aiui_length
   )
  {

    ds::BitContainer<T_BITSIZE> lbitcontainer_rowSrc(_st_numWordsRows,NULL);
    ds::BitContainer<T_BITSIZE> lbitcontainer_rowThis(_st_numWordsRows,NULL);

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "BitMatrix<T_BITSIZE>::copyAligned";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n";
   
      std::ostringstream lostrstream_this;
      lostrstream_this << lpc_labelFunc << ":this";
      this->print
	(std::cout,
	 lostrstream_this.str().c_str(),
	 '\0',
	 '\n'
	 );
      std::cout << '\n';
      std::ostringstream lostrstream_source;
      lostrstream_source << lpc_labelFunc << ":aibitmatrix_source";
      aibitmatrix_source.print
	(std::cout,
	 lostrstream_source.str().c_str(),
	 '\0',
	 '\n'
	 );
      std::cout	<< "\ninput  uintidx: aiui_idxTarget = " << aiui_idxTarget << '\n'
		<< " input  uintidx: aiui_length = " << aiui_length << '\n'
		<< ")\n";
    }
#endif /*__VERBOSE_YES*/


    for ( uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
      lbitcontainer_rowThis.setArray(this->getRow(lui_i));
      lbitcontainer_rowSrc.setArray(aibitmatrix_source.getRow(lui_i));
      lbitcontainer_rowThis.copyAligned(lbitcontainer_rowSrc,aiui_idxTarget,aiui_length);
    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      
       std::string lstr_termination = "OK";
       for ( uintidx  lui_i = 0; lui_i < this->getNumRows();  lui_i++) { 
	 for ( uintidx  lui_j = aiui_idxTarget; lui_j <
		 (aiui_idxTarget+aiui_length);
	       lui_j++)
	   {
	   if ( (*this)(lui_i,lui_j) != aibitmatrix_source(lui_i,lui_j) ) {
	     lstr_termination = "FAILURE";
	     break;
	   }
	     
	 }
       }
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ')' << " Termination:" << lstr_termination << "\n";
      this->print
      (std::cout,
       lpc_labelFunc,
        '\0',//',',
       '\n' //';'
       );
     std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  inline uintidx getRowNumBitOn(uintidx aiui_rowi) 
  {
    uintidx    lost_numBitOn;

    ds::BitContainer<T_BITSIZE> lbContainer_row(_st_numWordsRows,NULL);
     
    #ifdef __VERBOSE_YES 
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "BitMatrix<T_BITSIZE>::getRowNumBitOn  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\t(output BitMatrix<T_BITSIZE>: this[" << this << "]\n"
		<< "\t input  uintidx: aiui_rowi = " << aiui_rowi << '\n'
		<< "\t)\n";
    }
#endif /*__VERBOSE_YES*/

    lbContainer_row.setArray( this->getRow(aiui_rowi) );

    lost_numBitOn = lbContainer_row.getNumBitOn();

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "BitMatrix<T_BITSIZE>::getRowNumBitOn: OUT"
		<< '(' << geiinparam_verbose << ")\n"
		<<  "uintidx: lost_numBitOn = " << lost_numBitOn << '\n';
      
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
   
    return lost_numBitOn;
  }
  
  uintidx maxNorm()
  {
    ds::BitContainer<T_BITSIZE> lbitcontainer_row(_st_numWordsRows,NULL);
    uintidx lost_maxNorm = 0;

#ifdef __VERBOSE_YES 
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "BitMatrix<T_BITSIZE>::maxNorm  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\t(output BitMatrix<T_BITSIZE>: this[" << this << "]\n"
		<< "\t)\n";
    }
#endif /*__VERBOSE_YES*/
   
    for ( uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
      lbitcontainer_row.setArray(this->getRow(lui_i));
      uintidx lost_numBitOn = lbitcontainer_row.getNumBitOn();
      if ( lost_numBitOn > lost_maxNorm )
	lost_maxNorm = lost_numBitOn;
      
    } 

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "BitMatrix<T_BITSIZE>::maxNorm: OUT"
		<< '(' << geiinparam_verbose << ")\n"
		<< "uintidx: lost_maxNorm = " << lost_maxNorm << '\n';
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return lost_maxNorm;

  }

  void copyUnaligned
  (uintidx            aiui_idxTarget,  
   BitMatrix
   <T_BITSIZE>       &aibitmatrix_source, 
   uintidx            aiui_idxSource, 
   uintidx            aiui_length
   )
  {
    ds::BitContainer<T_BITSIZE> lbitcontainer_rowSrc(_st_numWordsRows,NULL);
    ds::BitContainer<T_BITSIZE> lbitcontainer_rowThis(_st_numWordsRows,NULL);
        
#ifdef __VERBOSE_YES 
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "BitMatrix<T_BITSIZE>::copyUnaligned  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\t(output BitMatrix<T_BITSIZE>: this[" << this << "]\n"
		<< "\t input  uintidx: aiui_idxTarget = " << aiui_idxTarget << '\n'
		<< "\t input  BitMatrix<T_BITSIZE>: &aibitmatrix_source[" 
		<< &aibitmatrix_source << "]\n"
		<< "\t input  uintidx: aiui_idxSource = " << aiui_idxSource << '\n'
		<< "\t input  uintidx: aiui_length = " << aiui_length << '\n'
		<< "\t)\n";
    }
#endif /*__VERBOSE_YES*/

    for ( uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
      lbitcontainer_rowThis.setArray(this->getRow(lui_i));
      lbitcontainer_rowSrc.setArray(aibitmatrix_source.getRow(lui_i));
      lbitcontainer_rowThis.copyUnaligned
	(aiui_idxTarget,lbitcontainer_rowSrc,aiui_idxSource,aiui_length);
    } 

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ds::BitContainer<T_BITSIZE>::copyUnaligned: OUT"
		<< '(' << geiinparam_verbose << ")\n";
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  void print_listAdj
  (std::ostream &os=std::cout,
   const char   *aipc_label = "",
   const char   aic_delimVert   = ',',
   const char   aic_delimAdj    =  ';'
   ) const
  {

#if defined(__VERBOSE_YES)
    os << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << this << ']'
       << aic_delimVert << "vertex"
       << aic_delimVert << this->_ui_numRows
       << '>';
#else
    os << aipc_label;
#endif 

    for(uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
      const T_BITSIZE* lbitarray_row = this->getRow(lui_i);
      os << lui_i;
      for(uintidx lui_j=0; lui_j < this->_ui_numColumns; lui_j++){
	if (lbitarray_row[_bitContainer_data._getWordIndexBit(lui_j)] >> 
	    (_bitContainer_data._getBitOffset(lui_j)) & 0x1 )
	  os << aic_delimVert << lui_j;
      }
      if ( lui_i <  this->_ui_numRows-1) 
	os << aic_delimAdj;
    }
    os <<'\n';
  }


  void  printByColumns
  (std::ostream    &os=std::cout,
   const  uintidx  aiui_numRows = 80,
   const char      aic_delimRow='\n'
   ) const
  {
    //const uintidx lui_numColumnsMinusOne  = this->getNumColumns() -1;
    //const uintidx lui_numRowMinusOne = this->getNumRows() -1;
    uintidx lui_limMin = 0;
    uintidx lui_limMax = aiui_numRows;
    lui_limMax = (lui_limMax > this->getNumColumns())?this->getNumColumns():lui_limMax;
    if ( this->_ui_numRows * this->_ui_numColumns > 0 ) {
      while ( lui_limMin < lui_limMax  ) {
	if ( this->getNumColumns() > aiui_numRows ) {
	  os << "Columns "  << lui_limMin << " through " << (lui_limMax-1) << ":"  << std::endl;
	}
	for(uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
	  const T_BITSIZE* lbitarray_row = this->getRow(lui_i);
	  for(uintidx lui_j=lui_limMin; lui_j < lui_limMax; lui_j++){
	    os << (lbitarray_row[_bitContainer_data._getWordIndexBit(lui_j)] >> 
		   (_bitContainer_data._getBitOffset(lui_j)) & 0x1);
	  }
	  os << aic_delimRow;
	}
	lui_limMin = lui_limMax;
	lui_limMax += aiui_numRows;
	lui_limMax = (lui_limMax > this->getNumColumns())?this->getNumColumns():lui_limMax;
	
      } /* While end  */
    }
    else {
      os << "NULL";
    }    
  }

  //! opXor
  /*! 
    xor
  */
  inline void opXor(const BitMatrix<T_BITSIZE> &aibitmatrixb_B)
  {
    _bitContainer_data.opXor(aibitmatrixb_B._bitContainer_data);    
  }
  
  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label = "",
   const char aic_delimCoef=',',
   const char aic_delimRow='\n'
   ) const
  {
    const uintidx lui_numColumnsMinusOne  = this->getNumColumns() -1;
    const uintidx lui_numRowMinusOne = this->getNumRows() -1;

    std::string lstr_delimCoef;

    if ( aic_delimCoef== '\0' ) {
      lstr_delimCoef = "";
    }
    else {
      lstr_delimCoef = aic_delimCoef;
    }
    
    os << aipc_label
       << ',' << "rows"
       << ',' << this->_ui_numRows
       << ',' << "columns"
       << ',' << this->_ui_numColumns
       << '>';
    if ( this->_ui_numRows * this->_ui_numColumns > 0 ){
      for(uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
	const T_BITSIZE* lbitarray_row = this->getRow(lui_i);
	for(uintidx lui_j=0; lui_j < lui_numColumnsMinusOne; lui_j++){
	  os << (lbitarray_row[_bitContainer_data._getWordIndexBit(lui_j)] >> 
		 (_bitContainer_data._getBitOffset(lui_j)) & 0x1)
	     << lstr_delimCoef; //aic_delimCoef;
	}
	os << (lbitarray_row[_bitContainer_data._getWordIndexBit(lui_numColumnsMinusOne)] >>
	       (_bitContainer_data._getBitOffset(lui_numColumnsMinusOne)) & 0x1);
	if ( lui_i < lui_numRowMinusOne )
	  os << aic_delimRow;
      }
    }
    else {
      os << "NULL";
    }    
  }

  friend std::ostream& operator<<(std::ostream& os, const BitMatrix<T_BITSIZE> &aibitmatrix_A)
  {
    aibitmatrix_A.print(os);
    
    return os;
  }
private:

  void _maskTopWord()
  {
    uintidx lui_bitsActive = this->_bits_in_top_word(this->_ui_numColumns);
    if ( lui_bitsActive > 0 )  {
      for ( uintidx lui_i = 0; lui_i < this->_ui_numRows; lui_i++) {
	T_BITSIZE* lT_matrixRow = this->getRow(lui_i);
	lT_matrixRow[_st_numWordsRows-1] &= this->_bitMask(lui_bitsActive);
      }
    }
  }

protected:

  inline void  _setNumRows(uintidx aiui_numRows) 
  {
    this->_ui_numRows = aiui_numRows;
  }

  inline void  _setNumColumns(uintidx aiui_numColumns) 
  {
    this->_ui_numColumns = aiui_numColumns;
  }
  
  ds::BitContainer<T_BITSIZE> _bitContainer_data;
  uintidx  _ui_numRows;        //m;
  uintidx  _ui_numColumns;     //n;
  uintidx  _st_numWordsRows;

}; /*BitMatrix*/


/*getRowsNumBitOn
 */
template < typename T_NUMBITON, 
	   typename T_BITSIZE 
	   >
void getRowsNumBitOn
(std::vector<T_NUMBITON>  &aovectort_numNumBitOn,
 BitMatrix<T_BITSIZE>     &aibmatrix_w
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "mat::getRowsNumBitOn";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output std::vector<T_NUMBITON>&: aovectort_numNumBitOn[" 
	      << &aovectort_numNumBitOn << "]\n"
	      << "\t  input BitMatrix<T_BITSIZE>: aibmatrix_w[" << &aibmatrix_w << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES*/
   
  for ( uintidx li_i = 0; li_i < aibmatrix_w.getNumRows(); li_i++) {   
    aovectort_numNumBitOn[li_i] = 
      (T_NUMBITON) aibmatrix_w.getRowNumBitOn(li_i);    
  }


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelVectorBitOn;
    lostrstream_labelVectorBitOn
      << "<VECTORBITON:"
      << geverbosepc_labelstep << ':' << lpc_labelFunc
      << ":aovectort_numNumBitOn["
      << geverboseui_idproc << ':'
      << &aovectort_numNumBitOn
      << ']';
    inout::containerprint
      (aovectort_numNumBitOn.begin(),
       aovectort_numNumBitOn.end(),
       std::cout,
       lostrstream_labelVectorBitOn.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

} /*END namespace ds*/

#endif /*BIT_MATRIX_HPP*/
