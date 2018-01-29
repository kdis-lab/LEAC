/*! \file matrix.hpp
 *
 * \brief matrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <sstream>      // std::ostringstream
#include <iomanip>      //std::setw
#include <typeinfo>
#include <stdexcept>
#include <utility>
#include <cassert>
#include <string>
#include <string.h> 
#include <vector>
#include "matrix_base.hpp"
#include "linear_algebra_level1.hpp"
#include "line_split.hpp"
#include "outfilename.hpp"

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/


namespace mat {


/*! \class MatrixRow
  \brief Matrix store items for row
*/
template <class T_FEATURE>
class MatrixRow 
  : public MatrixBase<T_FEATURE>
{
public:
  MatrixRow()
    : MatrixBase<T_FEATURE>()
    , _b_externalData(false)
    , _arrayt_data(NULL)
  { }

  MatrixRow(const uintidx aiuintidx_numRows)
    : MatrixBase<T_FEATURE>(aiuintidx_numRows)
    ,  _b_externalData(false)
    ,  _arrayt_data((aiuintidx_numRows > 0)?new T_FEATURE[aiuintidx_numRows]:NULL)
  { }


  MatrixRow
  (const uintidx aiuintidx_numRows,
   const uintidx aiuintidx_numColumns
   ) :  MatrixBase<T_FEATURE>(aiuintidx_numRows,aiuintidx_numColumns)
    ,  _b_externalData(false)
    , _arrayt_data((aiuintidx_numRows > 0 && aiuintidx_numColumns > 0)?new T_FEATURE[aiuintidx_numRows * aiuintidx_numColumns]:NULL) 
  { }

  MatrixRow
  (const uintidx aiuintidx_numRows,
   const uintidx aiuintidx_numColumns,
   const T_FEATURE ait_initialValue 
   ) :  MatrixBase<T_FEATURE>(aiuintidx_numRows,aiuintidx_numColumns)
    ,  _b_externalData(false)
    , _arrayt_data((aiuintidx_numRows > 0 && aiuintidx_numColumns > 0)?new T_FEATURE[aiuintidx_numRows * aiuintidx_numColumns]:NULL) 
  {
    this->initialize(ait_initialValue);
  }
  

  MatrixRow
  (const uintidx aiuintidx_numRows,
   const uintidx aiuintidx_numColumns,
   T_FEATURE* aiT_data
   )
    :  MatrixBase<T_FEATURE>(aiuintidx_numRows,aiuintidx_numColumns)
    , _b_externalData(true)
    , _arrayt_data(aiT_data)
  { }
  
  //copy constructor
  MatrixRow(const MatrixRow<T_FEATURE>& B) 
    : MatrixBase<T_FEATURE>(B)
    , _b_externalData(false)
    , _arrayt_data((B.getNumElems()>0)?new T_FEATURE[B.getNumElems()]:NULL)
  {
    interfacesse::copy(this->_arrayt_data, B._arrayt_data, B.getNumElems());
  }

  //move constructor
  MatrixRow(MatrixRow<T_FEATURE> &&B) 
    : MatrixBase<T_FEATURE>(B)
    , _b_externalData(B._b_externalData)
    , _arrayt_data(B._arrayt_data)
  {
    B._b_externalData = true;
    B._arrayt_data    = NULL;
  }

  virtual ~MatrixRow()
  {
    if( this->_arrayt_data != NULL && this->_b_externalData == false ) 
      delete[] _arrayt_data; 
  }


  inline const T_FEATURE *toArray() const
  { 
    return this->_arrayt_data; 
  }

  inline T_FEATURE* toArray() 
  { 
    return this->_arrayt_data; 
  }

  inline void  initialize(const T_FEATURE ait_value=T_FEATURE(0)) 
  {                    
    interfacesse::copya( this->_arrayt_data, T_FEATURE(ait_value), this->getNumElems() );  
  }

  inline T_FEATURE& operator[] (uintidx i)       { return _arrayt_data[i]; }
  inline T_FEATURE  operator[] (uintidx i) const { return _arrayt_data[i]; }

  MatrixRow<T_FEATURE>& operator=(const MatrixRow<T_FEATURE>& B)
  {
     
    if( this != &B ){
      uintidx lui_thisNumElems = this->getNumElems();
      uintidx lui_BNumElems    = B.getNumElems();
      MatrixBase<T_FEATURE>::operator=(B);
      //if  ( (_arrayt_data != NULL) && (lui_thisNumElems != lui_BNumElems) ) {
      if  ( lui_thisNumElems != lui_BNumElems ) {
	if  ( _arrayt_data != NULL ) delete[] _arrayt_data;
	//this->_b_externalData = false;
	_arrayt_data = new T_FEATURE[lui_BNumElems];
      }
      interfacesse::copy(this->_arrayt_data, B._arrayt_data,lui_BNumElems);
    }
    
    return *this;
  }
  

  MatrixRow<T_FEATURE>& operator=(MatrixRow<T_FEATURE> &&B)
  {
    if( this !=  &B ){
      MatrixBase<T_FEATURE>::operator=(B);
      if( this->_arrayt_data != NULL) 
	delete[] this->_arrayt_data;
      this->_b_externalData = B._b_externalData;
      this->_arrayt_data = B._arrayt_data;
      B._b_externalData  = true;
      B._arrayt_data     = NULL;
    }

    return *this;
  }


  //MATRIX row
  inline T_FEATURE& operator() (const uintidx i, const uintidx j)
  { return this->_arrayt_data[i*(this->getNumColumns())+j]; }
  
  inline T_FEATURE  operator() (const uintidx i, const uintidx j) const
  { return this->_arrayt_data[i*(this->getNumColumns())+j]; }
 
  
  inline const T_FEATURE* getRow (const uintidx i) const  
  {
    assert(0 <= i && i < this->getNumRows());

    return &this->_arrayt_data[i*(this->getNumColumns())]; 
  } 

  inline T_FEATURE* getRow (const  uintidx i)  
  {
    assert(0 <= i && i < this->getNumRows());

    return &this->_arrayt_data[i*(this->getNumColumns())]; 
  } 

  
  void  swapRows(const uintidx i, const uintidx j)
  {
    assert(0 <= i && i < this->getNumRows());
    assert(0 <= j && j < this->getNumRows());

    interfacesse::swap
      (this->getRow(i),this->getRow(j),this->getNumColumns());
  }

  void copyRow(const uintidx aiuintidx_idxFrom, const uintidx aiuintidx_idxTo)
  {
    assert(0 <= aiuintidx_idxFrom && aiuintidx_idxFrom < this->getNumRows());
    assert(0 <= aiuintidx_idxTo && aiuintidx_idxTo < this->getNumRows());

    interfacesse::copy
      (this->getRow(aiuintidx_idxTo),this->getRow(aiuintidx_idxFrom),this->getNumColumns());
  }

  void copyRow(const uintidx aiuintidx_idxRow, const T_FEATURE* aiT_row)
  {
    assert(0 <= aiuintidx_idxRow && aiuintidx_idxRow < this->getNumRows());

    interfacesse::copy
      (this->getRow(aiuintidx_idxRow),aiT_row,this->getNumColumns());
  }

  void copyRow(const uintidx aiuintidx_idxRow, const T_FEATURE* aiT_row, const uintidx aiuintidx_numColumns)
  {
    assert((0 <= aiuintidx_idxRow && aiuintidx_idxRow < this->getNumRows()) && (aiuintidx_numColumns <= this->getNumColumns()) );

    interfacesse::copy
      (this->getRow(aiuintidx_idxRow),aiT_row,aiuintidx_numColumns);
  }
    
  void cpyLoweToUpper()
  { 
   
    if ( this->getNumRows() != this->getNumColumns() )
      throw  std::range_error
	("MatrixRow<>::cpyLoweToUpper(): number of rows is different from the number of columns"); 
    for (uintidx li_i = 0; li_i < this->getNumRows(); li_i++) {
      T_FEATURE* larrayrow_data = this->getRow(li_i);
      for (uintidx li_j = li_i+1; li_j < this->getNumColumns(); li_j++) {
	larrayrow_data[li_j] = (*this)(li_j,li_i); 
      }
    }
  }

  void cpyUpper(const MatrixRow<T_FEATURE>& B)
  {
    T_FEATURE* larrayrow_data;

    if ( ( this->getNumRows() != B.getNumRows() ) &&  ( this->getNumColumns() != this->getNumColumns() ) )
      throw  std::range_error("MatrixRow<T_FEATURE>::cpyUpper(): order of the matrices is different");
    for (uintidx li_i = 0; li_i < this->getNumColumns(); li_i++) {
      larrayrow_data = this->getRow(li_i);
      for (uintidx  li_j = 0; li_j <= li_i; li_j++) {
	larrayrow_data[li_j] = larrayrow_data[li_j];
      }
    }
  }
 
  void cpySubMatrix
  (MatrixRow<T_FEATURE>& B, 
   const uintidx     aiuintidx_beginRow, 
   const uintidx     aiuintidx_beginCol, 
   const uintidx     aiuintidx_numRow, 
   const uintidx     aiuintidx_numCol
   )
  {

    for (uintidx luintidx_iDest = 0, luintidx_iSource = aiuintidx_beginRow; 
	 luintidx_iDest < aiuintidx_numRow;  
	 luintidx_iDest++, luintidx_iSource++) {
      T_FEATURE* lmatrixrow_row   = this->getRow(luintidx_iDest);
      T_FEATURE* lmr_sourceRow = B.getRow(luintidx_iSource);
      for (uintidx luintidx_jDest = 0, luintidx_jSource = aiuintidx_beginCol; 
	   luintidx_jDest < aiuintidx_numCol; luintidx_jDest++, luintidx_jSource++) {
	lmatrixrow_row[luintidx_jDest] = lmr_sourceRow[luintidx_jSource];
      }
    }
  }

   //these (+-*/) won't work for complex
  MatrixRow<T_FEATURE>& operator+=(const MatrixRow<T_FEATURE> &B)
  {
    if ( this->getNumRows() != B.getNumRows() && 
	 this->getNumColumns() != B.getNumColumns() )
      throw  std::range_error
	("MatrixRow<>::operator+=: order of the matrices is different");
    /*this = 1 * B + this*/
    interfacesse::axpy
      (this->_arrayt_data,
       T_FEATURE(1.0),
       B._arrayt_data,
       this->getNumElems()
       ); 

    return *this;
  }

  //these (+-*/) won't work for complex
  MatrixRow<T_FEATURE>& operator-=(const MatrixRow<T_FEATURE> &B)
  {
    if ( this->getNumRows() != B.getNumRows() && 
	 this->getNumColumns() != B.getNumColumns() )
      throw  std::range_error
	("MatrixRow<>::operator+=: order of the matrices is different");
    /*this = -1 * B + this*/
    interfacesse::axpy
      (this->_arrayt_data,
       T_FEATURE(-1.0),
       B._arrayt_data,
       this->getNumElems()
       ); 

    return *this;
    
  }
  
  std::pair<uintidx,uintidx> getMaximumElemIJ()
  {
    std::pair<uintidx,uintidx> lopair_maxElemIJ(0,0);
    T_FEATURE  lT_maxElem = this->_arrayt_data[0];
    
    for (uintidx li_i = 0; li_i < this->getNumRows(); li_i++) {
      T_FEATURE* larrayrow_data = this->getRow(li_i);
      for (uintidx li_j = 0; li_j < this->getNumColumns(); li_j++) {
	if ( larrayrow_data[li_j] > lT_maxElem ) {
	  lT_maxElem = larrayrow_data[li_j];
	  lopair_maxElemIJ.first  = li_i;
	  lopair_maxElemIJ.second = li_j;
	} 
      }
    }

    return lopair_maxElemIJ;
  }

  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label = "",
   const char aic_delimCoef=',',
   const char aic_delimRow='\n'
   ) const
  {
    const uintidx luintidx_numColumnsMinusOne  = this->getNumColumns() -1;
    const uintidx luintidx_numRowMinusOne = this->getNumRows() -1;

#if defined(__VERBOSE_YES)
    os << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << this << ']'
#else
    os << aipc_label
#endif 
       << ":rows,"
       << this->getNumRows()
       << ",columns,"
       << this->getNumColumns()
       << '>';
    if ( this->getNumElems() > 0  ) {
      for(uintidx luintidx_i=0; luintidx_i < this->getNumRows(); luintidx_i++) {
	const T_FEATURE* larrayrow_data = this->getRow(luintidx_i);
	for(uintidx luintidx_j=0; luintidx_j < luintidx_numColumnsMinusOne; luintidx_j++){
	  os << larrayrow_data[luintidx_j] << aic_delimCoef;
	} 
	os << larrayrow_data[luintidx_numColumnsMinusOne];
	if ( luintidx_i < luintidx_numRowMinusOne )
	  os << aic_delimRow;
      } 
    }
    else {
      os << "NULL";
    }
  }

  void r8mat_print
  (const char   *aipc_label       = "",
   const char   *aipc_labelColums = '\0',
   const int    aii_numColumns    = 5,
   std::ostream &os=std::cout
  ) const
//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in row-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, T A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
  {
    //# define INCX 5

 
    int m = (int) this->getNumRows();
    int n = (int) this->getNumColumns();
    
    //T_FEATURE a[];
    int ilo = 1;
    int jlo = 1;
    int ihi = (int) this->getNumRows();
    int jhi = (int) this->getNumColumns();
  
  
    int i;
    int i2hi;
    int i2lo;
    int j;
    int j2hi;
    int j2lo;

    os << "\n";
    os << aipc_label << "\n";

    std::string      lstr_colNames(aipc_labelColums);
    std::string      lstr_separator(1,inout::OutFileName::getDelim());
    inout::LineSplit lsplit_colNames;
    
    lsplit_colNames.setSeparator(lstr_separator);
   
    uintidx     lui_numColName = lsplit_colNames.split(lstr_colNames);
      
    if ( m <= 0 || n <= 0 )
      {
	os << "\n";
	os << "  (None)\n";
	return;
      }
    //
    //  Print the columns of the matrix, in strips of 5.
    //
    for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + aii_numColumns )
      {
	j2hi = j2lo + aii_numColumns - 1;
	if ( n < j2hi )
	  {
	    j2hi = n;
	  }
	if ( jhi < j2hi )
	  {
	    j2hi = jhi;
	  }
	os << "\n";
	//
	//  For each column J in the current range...
	//
	//  Write the header.
	//
	os << "  Col:    ";
	for ( j = j2lo; j <= j2hi; j++ )
	  {
	    if ( lui_numColName == 0 ) {
	      os << std::left <<std::setw(7) << j - 1 << "       ";
	    }
	    else {
	      os << std::left << std::setw(7) <<  lsplit_colNames.getItem(j) << "       ";
	    }
	  }
	os << "\n";
	os << "  Row\n";
	os << "\n";
	//
	//  Determine the range of the rows in this strip.
	//
	if ( 1 < ilo )
	  {
	    i2lo = ilo;
	  }
	else
	  {
	    i2lo = 1;
	  }
	if ( ihi < m )
	  {
	    i2hi = ihi;
	  }
	else
	  {
	    i2hi = m;
	  }

	for ( i = i2lo; i <= i2hi; i++ )
	  {
	    //
	    //  Print out (up to) 5 entries in row I, that lie in the current strip.
	    //
	    os << std::right << std::setw(5) << i - 1 << ": ";
	    for ( j = j2lo; j <= j2hi; j++ )
	      {
		os <<  std::left  << std::setw(12) << this->_arrayt_data[j-1+(i-1)*n] << "  ";
	      }
	    os << "\n";
	  }
      }

    return;
    //# undef INCX
    
  }
    
protected:

  bool          _b_externalData;
  T_FEATURE     *_arrayt_data;
  
};  /*MatrixRowT*/

  
template <class T_FEATURE>
MatrixRow<T_FEATURE> 
getIdentity(const uintidx aiui_dimension)
{
  MatrixRow<T_FEATURE> lomatrix_identity(aiui_dimension,aiui_dimension);

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "mat::getIdentity";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  uintidx: aiui_dimension = " << aiui_dimension << '\n'
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  lomatrix_identity.initialize();
  for (uintidx li_i = 0; li_i < lomatrix_identity.getNumRows();  li_i++) {
    lomatrix_identity(li_i, li_i) = (T_FEATURE) 1; 
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelIdentity;
    lostrstream_labelIdentity << "<IDENTITY:" << lpc_labelFunc  
			       << ", lomatrix_identity[" << &lomatrix_identity << ']';
    lomatrix_identity.print(std::cout,lostrstream_labelIdentity.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lomatrix_identity;
  
}

template <class T_FEATURE>
MatrixRow<T_FEATURE> 
getTranspose(MatrixRow<T_FEATURE>& aimatrix_B)
{
  MatrixRow<T_FEATURE> lomatrix_transpose(aimatrix_B.getNumColumns(), aimatrix_B.getNumRows());
   
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "matrix::getTranspose";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( input const MatrixRow<T_FEATURE>: &aimatrix_B[" << &aimatrix_B << ']'
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES


  T_FEATURE* larrayt_columnB = new T_FEATURE[aimatrix_B.getNumColumns()];
  T_FEATURE* lpt_columnB  = aimatrix_B.toArray();
  T_FEATURE* lpt_rowTranspose     = lomatrix_transpose.toArray();
  for (uintidx li_j = 0; li_j < aimatrix_B.getNumColumns();  li_j++) {
    T_FEATURE* lpt_itemColumn = lpt_columnB;
    for (uintidx li_i = 0; li_i < aimatrix_B.getNumRows();  li_i++) {
      larrayt_columnB[li_i] = *lpt_itemColumn;
      lpt_itemColumn  += aimatrix_B.getNumColumns();
    }
    interfacesse::copy
      (lpt_rowTranspose,larrayt_columnB,lomatrix_transpose.getNumColumns());
    ++lpt_columnB;
    lpt_rowTranspose += lomatrix_transpose.getNumColumns();
  }
  delete [] larrayt_columnB;
      
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelTranspose;
    lostrstream_labelTranspose << "<IDENTITY:" << lpc_labelFunc  
			      << ", lomatrix_transpose[" << &lomatrix_transpose << ']';
    lomatrix_transpose.print(std::cout,lostrstream_labelTranspose.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lomatrix_transpose;
  
}


template <class T_FEATURE,
	  class T_IDX
	  >
MatrixRow<T_FEATURE>
keepRows
(const MatrixRow<T_FEATURE>     &aimatrixT_b,
 std::vector<T_IDX>             &aovectorT_idxRow
 )
{
  MatrixRow<T_FEATURE> lomatrix_row
    (aovectorT_idxRow.size(),
     aimatrixT_b.getNumColumns()
     );

  for ( uintidx luintidx_i = 0; luintidx_i < aovectorT_idxRow.size(); luintidx_i++) {
    lomatrix_row.copyRow(luintidx_i,aimatrixT_b.getI(aovectorT_idxRow.at(luintidx_i)));
  }
    
  return lomatrix_row;
}

  
} /*END namespace mat*/


#endif /* MATRIX_HPP */

