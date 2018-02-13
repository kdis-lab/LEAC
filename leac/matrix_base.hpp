/*! \file matrix_base.hpp
 *
 * \brief matrix base
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MATRIX_BASE_HPP
#define MATRIX_BASE_HPP

#include <iostream>
#include <sstream>      // std::ostringstream
#include <typeinfo>
#include <stdexcept>
#include <utility>
#include <cassert>
#include <string.h> 
#include <vector>
#include "common.hpp"

#include "verbose_global.hpp"

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {

/*! \class MatrixBase
  \brief MatrixBase
*/
template <class T_FEATURE>
class MatrixBase
{public:
  MatrixBase()
    : _ui_numRows(0)
    , _ui_numColumns(0)
  { }

  MatrixBase(const uintidx aiui_numRows)
    : _ui_numRows(aiui_numRows)
    , _ui_numColumns(1)
  { }

  MatrixBase(const uintidx aiui_numRows, const uintidx aiui_numColumns)
    : _ui_numRows(aiui_numColumns > 0?aiui_numRows:0) 
    , _ui_numColumns(aiui_numRows > 0?aiui_numColumns:0) 
  { }

  //move constructor
  MatrixBase(MatrixBase<T_FEATURE> &&B) 
    : _ui_numRows(B._ui_numRows)
    , _ui_numColumns(B._ui_numColumns)
  {
    B._ui_numRows    = 0;
    B._ui_numColumns = 0;
  }
  
  //copy constructor
  MatrixBase(const MatrixBase<T_FEATURE>& B) 
    : _ui_numRows(B._ui_numRows)
    , _ui_numColumns(B._ui_numColumns)
  { }

  virtual ~MatrixBase() {}

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
    return ( this->getNumRows() * this->getNumColumns() ); 
  }

  virtual const T_FEATURE* getRow ( const uintidx i) const  = 0;

  virtual T_FEATURE* getRow ( const uintidx i) = 0;

  
  MatrixBase<T_FEATURE>& operator=(const MatrixBase<T_FEATURE>& B)
  {
    if( this != &B ){
      this->_ui_numRows = B._ui_numRows;
      this->_ui_numColumns = B._ui_numColumns;
    }
    return *this;
  }

  MatrixBase<T_FEATURE>& operator=(MatrixBase<T_FEATURE> &&B)
  {
    if( this !=  &B ){
      this->_ui_numRows = B._ui_numRows;
      this->_ui_numColumns = B._ui_numColumns;
      B._ui_numRows = 0;
      B._ui_numColumns = 0;
    }

    return *this;
  }

  
  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label = "",
   const char aic_delimCoef=',',
   const char aic_delimRow='\n'
   ) const = 0;
 
  //virtual void print(std::ostream &os=std::cout,const char aic_delimCoef='\t',const char aic_delimRow='\n') const = 0;

  friend std::ostream& operator<<(std::ostream& os, const MatrixBase<T_FEATURE> &aimatrix_a)
  {
    aimatrix_a.print(os,"",',',';');
    //aimatrix_a.print(os,"",outparam::OutFileName::getDelim(), '\n');
    
    return os;
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
    
  uintidx  _ui_numRows;        //m;
  uintidx  _ui_numColumns;     //n;

}; /*MatrixBase*/


} /*END namespace mat*/


#endif /* MATRIX_BASE_HPP */

