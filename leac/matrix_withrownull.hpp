/*! \file matrix_withrownull.hpp
 *
 * \brief Matrix with null of rows
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MATRIX_WITHROWNULL_HPP
#define MATRIX_WITHROWNULL_HPP

#include <iostream>
#include <typeinfo>
#include <stdexcept>
#include <algorithm> // std::copy
#include "outfilename.hpp"
#include "matrix_base.hpp"


/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {

#define MATRIX_WITHROWNULL_CHAR '#'


/*! \class MatrixWithRowNull
  \brief MatrixWithRowNull
*/
template <class T>
class MatrixWithRowNull: public MatrixBase<T>
{
public:

  MatrixWithRowNull()
    : MatrixBase<T>()
    , _ui_numRowsMax(0)
    , _arrayt_row(NULL) 
  { }

  MatrixWithRowNull
  (const uintidx aiui_numRows,
   const uintidx aiui_numRowsMax,
   const uintidx aiui_numColumns
   )
    : MatrixBase<T>(aiui_numRows,aiui_numColumns)
    , _ui_numRowsMax(aiui_numRowsMax)  
    , _arrayt_row(new T*[aiui_numRowsMax]) 
  { 
    for (uintidx li_i = 0; li_i < aiui_numRows; ++li_i) { 
      _arrayt_row[li_i] = new T[aiui_numColumns];
    }
    for (uintidx li_i =  aiui_numRows; li_i < aiui_numRowsMax; ++li_i) 
      _arrayt_row[li_i] = NULL;
  }

  //copy constructor
  MatrixWithRowNull(const mat::MatrixWithRowNull<T> &aimatrixwrownull_b)
    :  MatrixBase<T>(aimatrixwrownull_b)
    , _ui_numRowsMax(aimatrixwrownull_b._ui_numRowsMax)
    , _arrayt_row(new T*[aimatrixwrownull_b._ui_numRowsMax])
  {  
    for (uintidx li_i = 0; li_i < _ui_numRowsMax; ++li_i) { 
      if ( aimatrixwrownull_b._arrayt_row[li_i] != NULL ) { 
	_arrayt_row[li_i] = new T[this->getNumColumns()];
	interfacesse::copy
	  (_arrayt_row[li_i],
	   aimatrixwrownull_b._arrayt_row[li_i], 
	   this->getNumColumns()
	   );
      }
      else {
	_arrayt_row[li_i] = NULL;
      }
    }
  }

  //move constructor
  MatrixWithRowNull(mat::MatrixWithRowNull<T> &&aimatrixwrownull_b)
    : MatrixBase<T>(aimatrixwrownull_b)
    , _ui_numRowsMax(aimatrixwrownull_b.aiui_numRowsMax)
    , _arrayt_row(aimatrixwrownull_b._arrayt_row)
  {
    for (uintidx li_i = 0; li_i < _ui_numRowsMax; ++li_i) { 
      _arrayt_row[li_i] = aimatrixwrownull_b._arrayt_row[li_i]; 
    }
    aimatrixwrownull_b._ui_numRowsMax = 0;
    aimatrixwrownull_b._arrayt_row = NULL;
  }

  ~MatrixWithRowNull() 
  {
    if (  _arrayt_row !=  NULL ) { 
      for (uintidx li_i =  0; li_i < _ui_numRowsMax; ++li_i) { 
	if (_arrayt_row[li_i] != NULL )
	  delete [] _arrayt_row[li_i];
      }
      delete[] _arrayt_row;
    }
  }

  mat::MatrixRow<T>  getMatrix() 
  {
    mat::MatrixRow<T>       
      lomatrixT_extract
      (this->getNumRows(),
       this->getNumColumns()
       );

    uintidx lui_j = 0;
    for ( uintidx lui_i = 0; lui_i < _ui_numRowsMax; ++lui_i) {
      const T* lmr_dataRow = _arrayt_row[lui_i];
      if ( lmr_dataRow != NULL) {
	lomatrixT_extract.copyRow(lui_j++,lmr_dataRow);
      }
    }
    return lomatrixT_extract;
  }


  void setMatrix(const mat::MatrixRow<T> &aimatrix_B) 
  {
    uintidx lui_j = 0;
    
    if ( this->getNumRows() != aimatrix_B.getNumRows() ) 
      throw  std::range_error
	("mat::MatrixWithRowNull<>::setMatrix: order of the matrices is different");
    for ( uintidx lui_i = 0; lui_i < _ui_numRowsMax; ++lui_i) {
      T* lmr_dataRow = _arrayt_row[lui_i];
      if ( lmr_dataRow != NULL) {
	interfacesse::copy
	  (lmr_dataRow,
	   aimatrix_B.getRow(lui_j++),
	   this->getNumColumns()
	   );
      }
    } //for
    if ( this->getNumRows() != lui_j ) 
      throw  std::range_error
	("mat::MatrixWithRowNull<>::setMatrix: order of the matrix this is different");
  }

  void combination(uintidx aiui_numrow, mat::MatrixWithRowNull<T> &aimatrixwrownull_b)
  {
    uintidx lui_thisNumRows = this->getNumRows();
    uintidx lui_bNumRows    = aimatrixwrownull_b.getNumRows();

    for (uintidx lui_i = aiui_numrow; lui_i < _ui_numRowsMax; ++lui_i) {
      if (_arrayt_row[lui_i] == NULL && aimatrixwrownull_b._arrayt_row[lui_i] != NULL) {
	++lui_thisNumRows; --lui_bNumRows;
      }
      if (_arrayt_row[lui_i] != NULL && aimatrixwrownull_b._arrayt_row[lui_i] == NULL) {
	--lui_thisNumRows; ++lui_bNumRows;
      }
      std::swap(_arrayt_row[lui_i],aimatrixwrownull_b._arrayt_row[lui_i]);
    }
    this->_setNumRows(lui_thisNumRows);
    aimatrixwrownull_b._setNumRows(lui_bNumRows);
  }

  T* getRow ( uintidx i) 
  {
    return _arrayt_row[i];
  }

  const T* getRow ( uintidx i) const
  {
    return _arrayt_row[i];
  }
  
  void copyRow(const uintidx aiui_idxRow, const T* aiT_row)
  {
    interfacesse::copy
      (_arrayt_row[aiui_idxRow],aiT_row,this->getNumColumns());
  }

  inline const uintidx getNumRowsMax()    const 
  { 
    return this->_ui_numRowsMax; 
  }

  inline T** toArray() 
  { 
    return this->_arrayt_row; 
  }

  mat::MatrixWithRowNull<T>& operator=(const mat::MatrixWithRowNull<T>& aimatrixwrownull_b)
  {
    if( this != &aimatrixwrownull_b){
      MatrixBase<T>::operator=(aimatrixwrownull_b);
      if (  _arrayt_row !=  NULL ) { 
	for (uintidx li_i =  0; li_i < _ui_numRowsMax; ++li_i) { 
	  if (_arrayt_row[li_i] != NULL )
	    delete [] _arrayt_row[li_i];
	}
	delete[] _arrayt_row;
      }
      _ui_numRowsMax = aimatrixwrownull_b._ui_numRowsMax;
      _arrayt_row = new T*[_ui_numRowsMax];
      for (uintidx li_i = 0; li_i < _ui_numRowsMax; ++li_i) { 
	if ( aimatrixwrownull_b._arrayt_row[li_i] != NULL ) { 
	  _arrayt_row[li_i] = new T[this->getNumColumns()];
	  interfacesse::copy
	    (_arrayt_row[li_i],
	     aimatrixwrownull_b._arrayt_row[li_i], 
	     this->getNumColumns()
	     );
	}
	else {
	  _arrayt_row[li_i] = NULL;
     
	}
      }
    }
    return *this;
  }

  mat::MatrixWithRowNull<T>& operator=(mat::MatrixWithRowNull<T> &&aimatrixwrownull_b)
  {
    if( this !=  &aimatrixwrownull_b ){
      MatrixBase<T>::operator=(aimatrixwrownull_b);
      if (  _arrayt_row !=  NULL ) { 
	for (uintidx li_i =  0; li_i < _ui_numRowsMax; ++li_i) { 
	  if (_arrayt_row[li_i] != NULL )
	    delete [] _arrayt_row[li_i];
	}
	delete[] _arrayt_row;
      }
      _ui_numRowsMax = aimatrixwrownull_b._ui_numRowsMax;
      _arrayt_row    = aimatrixwrownull_b._arrayt_row;
      for (uintidx li_i = 0; li_i < _ui_numRowsMax; ++li_i) { 
	_arrayt_row[li_i] = aimatrixwrownull_b._arrayt_row[li_i]; 
      }
      aimatrixwrownull_b._ui_numRowsMax = 0;
      aimatrixwrownull_b._arrayt_row = NULL;
    }

    return *this;
  }

  /*keepRows: passes all the rows of the matrix B to a new array 
    and it voids at B, without repeated indices
  */
  void 
  keepRows(std::vector<uintidx> &aivectorT_idxRow)
  {
    uintidx lui_increasesRows = (uintidx) aivectorT_idxRow.size() * (uintidx) 2; 
    uintidx lui_numRowMax = (lui_increasesRows < _ui_numRowsMax)
      ? lui_increasesRows
      : _ui_numRowsMax;

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "MatrixWithRowNull::keepRows";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ')'
		<< "\n(output mat::MatrixWithRowNull<>&: this[" << this << ']'
		<< "\ninput   vector<uintidx> &aivectorT_idxRow[" << &aivectorT_idxRow << ']'
		<< "\n)\n";
    }
#endif /*__VERBOSE_YES*/

    T** larrayT_row = new T*[lui_numRowMax];
    for (uintidx lui_i = 0; lui_i < lui_numRowMax; lui_i++)
      larrayT_row[lui_i] = NULL;
    for ( uintidx lui_i = 0; lui_i < aivectorT_idxRow.size(); lui_i++) {
      std::swap(larrayT_row[lui_i],_arrayt_row[aivectorT_idxRow[lui_i]]);
    }
    for (uintidx li_i =  0; li_i < _ui_numRowsMax; ++li_i) { 
      if (_arrayt_row[li_i] != NULL )
	delete [] _arrayt_row[li_i];
    }
    delete[] _arrayt_row;

    _arrayt_row = larrayT_row;
    this->_setNumRows(aivectorT_idxRow.size());
    _ui_numRowsMax = lui_numRowMax;

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelMatrix;
      lostrstream_labelMatrix << lpc_labelFunc
			      << ':' << "MatrixWithRowNull:[" << this  << ']';
      this->print
      (std::cout,
       lostrstream_labelMatrix.str().c_str(),
       ',',
       ';'
       );
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  }

  void addRow(const T* aiT_row)
  {
    const uintidx lui_increasesRows = this->getNumRows() + (uintidx) 1;
    uintidx lui_numRowMax = (lui_increasesRows <= _ui_numRowsMax)
      ? _ui_numRowsMax
      : _ui_numRowsMax * (uintidx) 2;
 
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "MatrixWithRowNull::addRow";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ')'
		<< "\n(output mat::MatrixWithRowNull<>&: this[" << this << ']'
		<< "\n)\n";
    }
#endif /*__VERBOSE_YES*/

    if ( lui_numRowMax != _ui_numRowsMax ) {
      T** larrayT_row = new T*[lui_numRowMax];
      for (uintidx lui_i = 0; lui_i < lui_numRowMax; lui_i++)
	larrayT_row[lui_i] = NULL;
      for ( uintidx lui_i = 0; lui_i < this->getNumRows(); lui_i++) {
	std::swap(larrayT_row[lui_i],_arrayt_row[lui_i]);
      }
      for (uintidx li_i =  0; li_i < _ui_numRowsMax; ++li_i) { 
	if (_arrayt_row[li_i] != NULL )
	  delete [] _arrayt_row[li_i];
      }
      delete[] _arrayt_row;

      _arrayt_row = larrayT_row;
      _ui_numRowsMax = lui_numRowMax;
    }

    _arrayt_row[this->getNumRows()] = new T[this->getNumColumns()];
    interfacesse::copy
      (_arrayt_row[this->getNumRows()],aiT_row,this->getNumColumns());

    this->_setNumRows( lui_increasesRows );

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelMatrix;
      lostrstream_labelMatrix
	<< lpc_labelFunc
	<< ':' << "mat::MatrixWithRowNull<>&: this[" << this  << ']';
      this->print
	(std::cout,
	 lostrstream_labelMatrix.str().c_str(),
	 ',',
	 ';'
	 );
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label    = "",
   const char  aic_delimCoef = '\t',
   const char  aic_delimRow  = '\n'
   ) const
  {
    if ( 0 <  this->_ui_numRowsMax )  {

      uintidx  lui_numColumns = this->getNumColumns() - 1;
      uintidx  lui_numRowsMax = this->_ui_numRowsMax - 1;
      const T* lmr_dataRow;
      
      os << aipc_label
	 << ":rows,"    << this->getNumRows()
	 << ",columns,"  << this->getNumColumns()
	<< ",rowsmax,"  << _ui_numRowsMax
	 << '>';
      for(uintidx lui_i=0; lui_i< lui_numRowsMax; lui_i++) {
	lmr_dataRow = _arrayt_row[lui_i];
  
	if ( lmr_dataRow != NULL) {
	  for(uintidx lui_j=0; lui_j < lui_numColumns; lui_j++){
	    os << lmr_dataRow[lui_j] <<  aic_delimCoef; 
	  } 
	  os << lmr_dataRow[lui_numColumns] << aic_delimRow;
	}
	else {
	  os << MATRIX_WITHROWNULL_CHAR << aic_delimRow;
	}
      } /*end for*/
      lmr_dataRow = _arrayt_row[lui_numRowsMax];
      if ( lmr_dataRow != NULL) {
	for(uintidx lui_j=0; lui_j < lui_numColumns; lui_j++){
	  os << lmr_dataRow[lui_j] <<  aic_delimCoef;
	} 
	os << lmr_dataRow[lui_numColumns];
      }
      else {
	os << MATRIX_WITHROWNULL_CHAR;
      }
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const mat::MatrixWithRowNull<T> &aimatrixwrownull_a)
  {
    aimatrixwrownull_a.print(os,"",inout::OutFileName::getDelim(),'\n');
    
    return os;
  }

protected:

  uintidx  _ui_numRowsMax;
  T        **_arrayt_row;

}; /*MatrixWithRowNull*/

} /*END namespace mat*/

#endif /*MATRIX_WITHROWNULL_HPP*/
