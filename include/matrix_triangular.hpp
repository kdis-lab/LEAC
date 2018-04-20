/*! \file matrix_triangular.hpp
 *
 * \brief matrix triangular
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MATRIX_TRIANG_HPP
#define MATRIX_TRIANG_HPP

#include <iostream>
#include "common.hpp"

#define MATRIX_TRIANG_OUT_SEPARATOR_DEFAULT ','

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {


/*! \class MatrixTriang
  \brief MatrixTriang
*/

template <class T_FEATURE>
class MatrixTriang 
{
public:
  MatrixTriang(): 
    _ui_numRows(0), 
    _t_data(NULL)
  { }

  MatrixTriang(const uintidx aiuintidx_numRows)
    : _ui_numRows(aiuintidx_numRows)
    , _t_data(new T_FEATURE*[aiuintidx_numRows]) 
  {
    for(uintidx li_i = 0; li_i < aiuintidx_numRows; ++li_i) {
      _t_data[li_i] = new T_FEATURE[li_i+1];
    }
  }

  //copy constructor
  MatrixTriang(const MatrixTriang<T_FEATURE>& B)
    : _ui_numRows(B._ui_numRows)
    , _t_data(new T_FEATURE*[B._ui_numRows]) 
  {
    for(uintidx lui_i = 0; lui_i < B._ui_numRows; ++lui_i) {
      _t_data[lui_i] = new T_FEATURE[lui_i+1];
      interfacesse::copy(&_t_data[lui_i][0],&B._t_data[lui_i][0],lui_i+1);
    }
  }

  //move constructor
  MatrixTriang(MatrixTriang<T_FEATURE> &&B)
    : _ui_numRows(B._ui_numRows)
    , _t_data(B._t_data)
  {  
    B._ui_numRows = true;
    B._t_data = NULL;
  }
  
  ~MatrixTriang()
  {
     if ( _t_data != NULL ) {
       for(uintidx li_i = 0; li_i < _ui_numRows; ++li_i) {
	 delete[] _t_data[li_i];
       }
       delete[] _t_data;
     }
  }

  MatrixRow<T_FEATURE>& operator=(const MatrixRow<T_FEATURE>& B)
  {
    if( this != &B ){
      if  ( _ui_numRows != B._ui_numRows ) {
	if  ( _t_data != NULL ) {
	  for(uintidx li_i = 0; li_i < _ui_numRows; ++li_i) {
	    delete[] _t_data[li_i];
	  }
	  delete[] _t_data;
	}
	_t_data = new T_FEATURE*[B._ui_numRows];
	for(uintidx lui_i = 0; lui_i < B._ui_numRows; ++lui_i) {
	  _t_data[lui_i] = new T_FEATURE[lui_i+1];
	}	
      }
      for(uintidx lui_i = 0; lui_i < B._ui_numRows; ++lui_i) {
	  interfacesse::copy(&_t_data[lui_i][0],&B._t_data[lui_i][0],lui_i+1);
      } 
    }
    
    return *this;
  }
  

  MatrixRow<T_FEATURE>& operator=(MatrixRow<T_FEATURE> &&B)
  {
    if( this !=  &B ){
      if  ( _t_data != NULL ) {
	for(uintidx li_i = 0; li_i < _ui_numRows; ++li_i) {
	  delete[] _t_data[li_i];
	}
	delete[] _t_data;
      }
      _ui_numRows   = B._ui_numRows;
      this->_t_data = B._t_data;
      
      B._ui_numRows = 0;
      B._t_data     = NULL;
    }

    return *this;
  }

  
  inline const uintidx getNumRows()    const 
  { 
    return this->_ui_numRows; 
  }
 
  inline T_FEATURE& operator() (uintidx i, uintidx j)
  {
    assert(0 <= i && i < this->getNumRows() && 0 <= j && j < this->getNumRows());
    return (i>j)?this->_t_data[i][j]:this->_t_data[j][i];
  }

  inline T_FEATURE  operator() (uintidx i, uintidx j) const
  {
    assert(0 <= i && i < this->getNumRows() && 0 <= j && j < this->getNumRows());
    return (i>j)?this->_t_data[i][j]:this->_t_data[j][i];
  }
  
  const T_FEATURE* getRow (const uintidx i) const  
  { 
    assert(0 <= i && i < this->getNumRows());
    return &this->_t_data[i]; 
  } 

  T_FEATURE* getRow (const  uintidx i)  
  {
    assert(0 <= i && i < this->getNumRows());
    return this->_t_data[i]; 
  }

  
  void print
  (std::ostream &os=std::cout,
   const char   aic_delimCoef = ',',
   const char   aic_delimRow  = '\n'
   ) const
  {
    for(uintidx luintidx_i = 0; luintidx_i < this->_ui_numRows; luintidx_i++) {
      std::cout << luintidx_i  << ":\t";
      for(uintidx luintidx_j = 0; luintidx_j < luintidx_i; luintidx_j++) {
	std::cout <<  this->_t_data[luintidx_i][luintidx_j] << aic_delimCoef;
      }
      std::cout << this->_t_data[luintidx_i][luintidx_i] << aic_delimRow;
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const mat::MatrixTriang<T_FEATURE> &aimatrix_a)
  {
    aimatrix_a.print(os,MATRIX_TRIANG_OUT_SEPARATOR_DEFAULT,'\n');
    
    return os;
  }

protected:
  uintidx    _ui_numRows;
  T_FEATURE  **_t_data;
} /*MatrixTriang*/;

} /*END namespace mat*/
  
#endif /* MATRIX_TRIANG_HPP */
