/*! \file matrix_resizablerow.hpp
 *
 * \brief Matrix with resizable number of rows
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MATRIX_RESIZABLEROW_HPP
#define MATRIX_RESIZABLEROW_HPP

#include <limits>
#include <cmath> /*round*/
#include "matrix_base.hpp"


/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {


/*! \class MatrixResizableRow
  \brief MatrixResizableRow
  \cite Franti:etal:GAclustering:gafranti:1997
*/

template <class T_FEATURE,
	  class T_WEIGHT
	  >
class MatrixResizableRow: public MatrixBase<T_FEATURE> {
public:
  MatrixResizableRow()
    : MatrixBase<T_FEATURE>()
    , _ui_numRowsMaximum(0)
    , _arrayt_data(NULL)
  { }
  
  MatrixResizableRow
  (const uintidx aiui_numRows,         
   const uintidx aiui_numColumns,
   const uintidx aiui_numRowsMaximum 
   )
    : MatrixBase<T_FEATURE>(aiui_numRows,aiui_numColumns)
    , _ui_numRowsMaximum(aiui_numRowsMaximum)
    , _arrayt_data
    ((aiui_numRows > 0 && aiui_numColumns > 0)?
     new T_FEATURE[aiui_numRowsMaximum * aiui_numColumns]:NULL)
  { } 

  //copy constructor
  MatrixResizableRow
  (const MatrixResizableRow<T_FEATURE,T_WEIGHT> &aimatrixresizerow_B)
    : MatrixBase<T_FEATURE>(aimatrixresizerow_B)
    , _ui_numRowsMaximum(aimatrixresizerow_B._ui_numRowsMaximum)
    , _arrayt_data
    ((aimatrixresizerow_B._ui_numRowsMaximum > 0
      && aimatrixresizerow_B.getNumColumns() > 0)?
     new T_FEATURE[aimatrixresizerow_B._ui_numRowsMaximum * aimatrixresizerow_B.getNumColumns()]:NULL)
  {
    interfacesse::copy
      (this->_arrayt_data, 
       aimatrixresizerow_B._arrayt_data, 
       aimatrixresizerow_B._ui_numRowsMaximum *  aimatrixresizerow_B.getNumColumns()
       );		
  }

  //move constructor
  MatrixResizableRow
  (MatrixResizableRow<T_FEATURE,T_WEIGHT> &&aimatrixresizerow_B)
    : mat::MatrixRow<T_FEATURE>(aimatrixresizerow_B)
    , _ui_numRowsMaximum(aimatrixresizerow_B._ui_numRowsMaximum)
    , _arrayt_data(aimatrixresizerow_B._arrayt_data)
  {  
    aimatrixresizerow_B._ui_numRowsMaximum = 0;
    aimatrixresizerow_B._arrayt_data = NULL;
  }

  virtual ~MatrixResizableRow()
  {
    if( this->_arrayt_data != NULL ) 
      delete[] _arrayt_data;
  }
  
  MatrixResizableRow<T_FEATURE,T_WEIGHT>& operator=
  (const MatrixResizableRow<T_FEATURE,T_WEIGHT> &B)
  {
    if( this != &B ){
      uintidx lui_thisNumElems = this->getNumElemsMaximum();
      uintidx lui_BNumElems = B.getNumElemsMaximum();
      MatrixBase<T_FEATURE>::operator=(B);
      _ui_numRowsMaximum = B._ui_numRowsMaximum;
      if  (lui_thisNumElems != lui_BNumElems ) {
	delete[] _arrayt_data;
	_arrayt_data = new T_FEATURE[lui_BNumElems];
      }
      interfacesse::copy(this->_arrayt_data, B._arrayt_data,B.getNumElems());
    }
    
    return *this;
  }
  
  MatrixResizableRow<T_FEATURE,T_WEIGHT>& operator=
  (MatrixResizableRow<T_FEATURE,T_WEIGHT> &&aimatrixresizerow_B)
  {
    if( this !=  &aimatrixresizerow_B ){
      MatrixBase<T_FEATURE>::operator=(aimatrixresizerow_B);
      delete[] _arrayt_data;
      _ui_numRowsMaximum = aimatrixresizerow_B._ui_numRowsMaximum;    
      this->_arrayt_data = aimatrixresizerow_B._arrayt_data;
      aimatrixresizerow_B._ui_numRowsMaximum = 0;
      aimatrixresizerow_B._arrayt_data     = NULL;
    }

    return *this;
  }

  void resize(const uintidx aiui_numRowsMaximum)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "MatrixResizableRow::resize";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
      		<< "(output MatrixResizableRow: this[" << this << "]\n"
	        << " input  aiui_numRowsMaximum: " <<  aiui_numRowsMaximum
	        << " _ui_numRowsMaximum: " << _ui_numRowsMaximum
	        << " this->getNumRows(): " << this->getNumRows()
	        << " this->getNumColumns(): " << this->getNumColumns()
		<< '\n' 
		<< ')'
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( aiui_numRowsMaximum != this->_ui_numRowsMaximum ) {
      T_FEATURE *larrayt_datatmp =
	(aiui_numRowsMaximum > 0 && this->getNumColumns() > 0)?
	new T_FEATURE[ aiui_numRowsMaximum * this->getNumColumns()]:NULL;
    
      if ( this->getNumRows() <  aiui_numRowsMaximum ) {
	interfacesse::copy
	  (larrayt_datatmp,
	   _arrayt_data, 
	   this->getNumRows() * this->getNumColumns() 
	   );
	 interfacesse::copya
	  (larrayt_datatmp + this->getNumRows() * this->getNumColumns(), 
	   T_FEATURE(0), 
	   (aiui_numRowsMaximum - this->getNumRows()) * this->getNumColumns()
	   );
      }
      else if (aiui_numRowsMaximum < this->_ui_numRowsMaximum ) {
	interfacesse::copy
	  (larrayt_datatmp, 
	   _arrayt_data, 
	   aiui_numRowsMaximum * this->getNumColumns() 
	   );
      }
      if ( larrayt_datatmp != NULL) { 
	delete [] _arrayt_data;
	_arrayt_data = larrayt_datatmp;
      }
      this->_ui_numRowsMaximum = aiui_numRowsMaximum;
    }
    this->_setNumRows(aiui_numRowsMaximum);
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }
  
  void merge
  (const MatrixResizableRow<T_FEATURE,T_WEIGHT> &aimatrixresizerow_B,
   const MatrixResizableRow<T_FEATURE,T_WEIGHT> &aimatrixresizerow_C
   )
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "MatrixResizableRow::merge";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
      		<< "(output MatrixResizableRow: this[" << this << "]\n";
      std::ostringstream lostrstream_labelmatrixresizerow_B;
      lostrstream_labelmatrixresizerow_B
	<< "<MatrixResizableRow:"
	<< lpc_labelFunc  
	<< ":aimatrixresizerow_B[" << &aimatrixresizerow_B << ']';
      aimatrixresizerow_B.print(std::cout,lostrstream_labelmatrixresizerow_B.str().c_str(),',',';');
      std::cout << '\n';
      std::ostringstream lostrstream_labelmatrixresizerow_C;
      lostrstream_labelmatrixresizerow_C
	<< "<MatrixResizableRow:"
	<< lpc_labelFunc  
	<< ":aimatrixresizerow_C[" << &aimatrixresizerow_C << ']';
      aimatrixresizerow_C.print(std::cout,lostrstream_labelmatrixresizerow_C.str().c_str(),',',';');
      std::cout << std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this->getNumColumns() != aimatrixresizerow_B.getNumColumns() ||
	 this->getNumColumns() != aimatrixresizerow_C.getNumColumns() )
      throw  std::range_error
	("mat::MatrixResizableRow::merge: error - the number of columns in the matrices is different"); 
    if (  _ui_numRowsMaximum <
	  aimatrixresizerow_B.getNumRows() + aimatrixresizerow_C.getNumRows() ) {
      _ui_numRowsMaximum =
	aimatrixresizerow_B.getNumRows() + aimatrixresizerow_C.getNumRows();
      delete[] _arrayt_data;
      _arrayt_data =
	new T_FEATURE[aimatrixresizerow_B.getNumElems() + aimatrixresizerow_C.getNumElems()];
    }
     
    T_FEATURE  *larray_tmp = _arrayt_data;
   
    interfacesse::copy
      (larray_tmp, 
       aimatrixresizerow_B.toArray(), 
       aimatrixresizerow_B.getNumElems()
       );
    
    larray_tmp += aimatrixresizerow_B.getNumElems();
   
    interfacesse::copy
      (larray_tmp,
       aimatrixresizerow_C.toArray(), 
       aimatrixresizerow_C.getNumElems()
       );

    this->_setNumRows
      (aimatrixresizerow_B.getNumRows()  +
       aimatrixresizerow_C.getNumRows()
       );
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelthis;
      lostrstream_labelthis
	<< "<MatrixResizableRow:"
	<< lpc_labelFunc  
	<< ":this[" << this << ']';
      this->print(std::cout,lostrstream_labelthis.str().c_str(),',',';');
      std::cout << std::endl;
      std::cout<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }
  
  inline const uintidx getNumElemsMaximum() const 
  { 
    return (this->_ui_numRowsMaximum)*(this->getNumColumns()); 
  }

  inline const uintidx getNumRowsMaximum() const
  {
    return this->_ui_numRowsMaximum;
  }

  inline const T_FEATURE *toArray() const
  { 
    return this->_arrayt_data; 
  }

  inline T_FEATURE* toArray() 
  { 
    return this->_arrayt_data; 
  }
  
  inline void  initialize() 
  {
    interfacesse::copya
      (this->_arrayt_data, T_FEATURE(0), this->getNumElemsMaximum());	
  }

  const T_FEATURE* getRow (const uintidx i) const  
  {
    assert(0 <= i && i < this->getNumRowsMaximum());

    return &this->_arrayt_data[i*(this->getNumColumns())]; 
  } 

  T_FEATURE* getRow (const  uintidx i)  
  {
    assert(0 <= i && i < getNumRowsMaximum());
    
    return &this->_arrayt_data[i*(this->getNumColumns())]; 
  } 

  inline void decreaseRow() 
  {
    uintidx lst_numRows = this->getNumRows();
    --lst_numRows;
    assert(0 <= lst_numRows &&  lst_numRows < this->getNumRowsMaximum());
    this->_setNumRows(lst_numRows);
  }

  inline void increasesRow() 
  {
    uintidx lst_numRows = this->getNumRows();
    ++lst_numRows;
    assert(0 <= lst_numRows &&  lst_numRows < this->getNumRowsMaximum());
    this->_setNumRows(lst_numRows);
  }
  
  void mergeTwoRow
  (std::pair<uintidx,uintidx> aipairst_idxRow,
   T_WEIGHT                 aiT_weightRowFirst,
   T_WEIGHT                 aiT_weightRowSecond
   )
  {
    T_FEATURE* larrayrowt_first  = this->getRow( aipairst_idxRow.first );
    T_FEATURE* larrayrowt_second = this->getRow( aipairst_idxRow.second );

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "MatrixResizableRow::mergeTwoRow";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
		<< "\t(output MatrixResizableRow: this[" << this << "]\n"
		<< "\t input  pair: aipairst_idxRow: (" 
		<< aipairst_idxRow.first << ", " << aipairst_idxRow.second << ")\n"
		<< "\t)\n";
    
      std::ostringstream lostrstream_labelRowFirst;
      lostrstream_labelRowFirst << "<ARRAYROW_FIRST:"
				<< lpc_labelFunc
				<< ':' << aipairst_idxRow.first << ':';
      inout::containerprint
	(larrayrowt_first,
	 larrayrowt_first + this->getNumColumns(),
	 std::cout,
	 lostrstream_labelRowFirst.str().c_str(),
	 ','
	 );
      std::cout << '\n';
      std::ostringstream lostrstream_labelRowSecond;
      lostrstream_labelRowSecond << "<ARRAYROW_SECOND:"
				<< lpc_labelFunc
				<< ':' << aipairst_idxRow.second << ':';
      inout::containerprint
	(larrayrowt_second,
	 larrayrowt_second + this->getNumColumns(),
	 std::cout,
	 lostrstream_labelRowSecond.str().c_str(),
	 ','
	 );
      std::cout << std::endl;
    }
#endif /*__VERBOSE_YES*/

    /*
                 (WeightRowFirst * RowFirst +  WeightRowSecond *  RowSecond ) 
      RowFirst = ------------------------------------------------------------
                            (WeightRowFirst + WeightRowSecond)

       
                  (WeightRowSecond  *  RowSecond  + WeightRowFirst * RowFirst 
      RowFirst = ------------------------------------------------------------
                            (WeightRowFirst + WeightRowSecond)

      y = ax + by

                 WeightRowSecond
      a   = ---------------------------------
             (WeightRowFirst + WeightRowSecond)

                   WeightRowFirst
      b   = ---------------------------------
                (WeightRowFirst + WeightRowSecond)

      x  =  RowSecond

      y  =  RowFirst
       
    */
 
    T_WEIGHT lT_denominator  = aiT_weightRowFirst + aiT_weightRowSecond;
    
    if ( lT_denominator != 0 ) {
      interfacesse::scal
	(larrayrowt_first,
	 aiT_weightRowFirst,
	 this->getNumColumns()
	 );

      interfacesse::axpy
	(larrayrowt_first,
	 aiT_weightRowSecond,
	 larrayrowt_second,
	 this->getNumColumns()
	 );

      interfacesse::scalInv
	(larrayrowt_first,
	 lT_denominator,
	 this->getNumColumns()
	 );
	 	
    } //IF ( lT_denominator != 0 ) 


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelRowFirst;
      lostrstream_labelRowFirst << "<ARRAYROW_FIRST:"
				<< lpc_labelFunc
				<< ':' << aipairst_idxRow.first << ':';
      inout::containerprint
	(larrayrowt_first,
	 larrayrowt_first + this->getNumColumns(),
	 std::cout,
	 lostrstream_labelRowFirst.str().c_str(),
	 ','
	 );
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  }

  void copyRow(const uintidx aiui_idxFrom, const uintidx aiui_idxTo)
  {
    assert(0 <= aiui_idxFrom && aiui_idxFrom < this->getNumRows());
    assert(0 <= aiui_idxTo && aiui_idxTo < this->getNumRows());
    interfacesse::copy
      (this->getRow(aiui_idxTo),this->getRow(aiui_idxFrom),this->getNumColumns());
  }
  
  inline 
  void 
  copyRows(const MatrixResizableRow<T_FEATURE,T_WEIGHT>& B, const uintidx aiui_idxRows) 
  {
    assert(1 <= aiui_idxRows && aiui_idxRows <= this->getNumRows());
    interfacesse::copy
      ( this->_arrayt_data, B._arrayt_data, aiui_idxRows * this->getNumColumns());		
  }

  void  swapRows(const uintidx i, const uintidx j)
  {
    assert(0 <= i && i < this->getNumRows());
    assert(0 <= j && j < this->getNumRows());
    
    interfacesse::swap
      (this->getRow(i),this->getRow(j),this->getNumColumns());
  }

  void removeRow(const uintidx aiui_idxRow) 
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "MatrixResizableRow::removeRow";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n";
    }
#endif /*__VERBOSE_YES*/

    assert(0 <= aiui_idxRow &&  aiui_idxRow < this->getNumRows());
    uintidx lui_idxLastRow = this->getNumRows() -1;

    if ( aiui_idxRow != lui_idxLastRow ) {
      this->copyRow(lui_idxLastRow, aiui_idxRow);
    }
    this->decreaseRow();

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  } /*removeLastCentroid END*/

  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label = "",
   const char aic_delimCoef=',',
   const char aic_delimRow='\n'
   ) const
  {
    const uintidx lui_numColumns     = this->getNumColumns() -1;
    const uintidx lui_numRowMinusOne = this->getNumRows() -1;

#if defined(__VERBOSE_YES)
    os << aipc_label  << ':'
       <<  geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << this << ']'
#else
    os << aipc_label
#endif 
       << aic_delimCoef << "rows"
       << aic_delimCoef << this->getNumRows()
       << aic_delimCoef << "columns"
       << aic_delimCoef << this->getNumColumns()
       << aic_delimCoef << "rowsmaximum"
       << aic_delimCoef << _ui_numRowsMaximum
       << '>';
    if ( this->getNumElems() > 0  ) {
      for(uintidx lui_i=0; lui_i < this->getNumRows(); lui_i++) {
	const T_FEATURE* larrayrow_data = this->getRow(lui_i);
	for(uintidx lui_j=0; lui_j < lui_numColumns; lui_j++){
	  os << larrayrow_data[lui_j] << aic_delimCoef;
	} 
	os << larrayrow_data[lui_numColumns];
	if ( lui_i < lui_numRowMinusOne )
	  os << aic_delimRow;
      } 
    }
    else {
      os << "NULL";
    }
  }

protected:
  uintidx    _ui_numRowsMaximum; 
  T_FEATURE  *_arrayt_data;
 
}; /*MatrixResizableRow*/ 

} /*END namespace mat*/
  
#endif /* MATRIX_RESIZABLEROW_HPP */
