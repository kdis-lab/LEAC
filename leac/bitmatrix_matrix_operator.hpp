/*! \file bitmatrix_matrix_operator.hpp
 *
 * \brief  bitmatrix matrix operator
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __BITMATRIX_MATRIX_OPERATOR_HPP
#define __BITMATRIX_MATRIX_OPERATOR_HPP

#include "bit_array.hpp"
#include "bit_matrix.hpp"
#include "matrix.hpp"
#include "linear_algebra_level1.hpp"

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/


namespace mat {

  
/*! \fn void mulRowsIColumns0N(T_C *aoarrayT_rowC, uintidx aiuintidx_idxBegin, uintidx aiuintidx_idxEnd, mat::BitArray<T_BITSIZE>  &aibarray_rowA, MatrixRow<T_B> &aimatrixrowt_B
 )
    \brief mulRowsIColumns0N:
    \details
    \param aoarrayT_rowC a array of type T_C
    \param aiuintidx_idxBegin a uintidx 
    \param aiuintidx_idxEnd a uintidx 
    \param aibarray_rowA a mat::BitArray<T_BITSIZE>
    \param aimatrixrowt_B a MatrixRow<T_B>
 */
template < typename T_C, 
	   typename T_BITSIZE,
           typename T_B
	   >
void
mulRowsIColumns0N
(T_C                      *aoarrayT_rowC,
 uintidx                  aiuintidx_idxBegin,
 uintidx                  aiuintidx_idxEnd,
 mat::BitArray<T_BITSIZE>  &aibarray_rowA,
 MatrixRow<T_B>           &aimatrixrowt_B
 )
{
  interfacesse::copya
    ( aoarrayT_rowC, 
      T_C(0), 
      aimatrixrowt_B.getNumColumns()
      );  
  for ( uintidx li_i = aiuintidx_idxBegin; li_i < aiuintidx_idxEnd; li_i++) {
    if ( aibarray_rowA.getBit(li_i) ) {
      
      interfacesse::axpy
	(aoarrayT_rowC,
	 1.0,
	 aimatrixrowt_B.getRow(li_i),
	 aimatrixrowt_B.getNumColumns()
	 );
    }
  }
}


/*! \fn void bitgemm(MatrixRow<T_C>  &aomatrixrowt_C, mat::BitMatrix<T_BITSIZE> &aibitmatrixT_A, MatrixRow<T_B> &aimatrixrowt_B)
    \brief bitgemm multiply a matrix of bits by another of real or integers
    \details 
    \f[
     C = A  \times B,
    \f]
    \param aomatrixrowt_C a The resulting matrix  mat::MatrixRow<T_C>
    \param aibitmatrixT_A a matrix of bits mat::BitMatrix
    \param aimatrixrowt_B a matrix of integers or real MatrixRow<T_B>
 */
template < typename T_C,
	   typename T_BITSIZE,
	   typename T_B
	   >
void 
bitgemm
(MatrixRow<T_C>               &aomatrixrowt_C,  
 mat::BitMatrix<T_BITSIZE>     &aibitmatrixT_A, 
 MatrixRow<T_B>               &aimatrixrowt_B
 )
{
  mat::BitArray<T_BITSIZE> lbarray_row(aibitmatrixT_A.getNumColumns(),NULL);

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "mat::bitgemm";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output MatrixRow<T_C>: aomatrixrowt_C[" << &aomatrixrowt_C << "]\n"
	      << " input  MatrixByte: aibitmatrixT_A[" <<  &aibitmatrixT_A << "]\n"
	      << " input  MatrixRow<T_B>: aimatrixrowt_B[" << &aimatrixrowt_B << "]\n"
	      << ")\n";
  }
#endif /*__VERBOSE_YES*/

    if ( aibitmatrixT_A.getNumColumns() != aimatrixrowt_B.getNumRows()  )
      throw  std::range_error
	("mat::bitgemm:"
	 "the number of columns of the first matrix is different than the second"
	 );

  for ( uintidx luintidx_i = 0; luintidx_i < aibitmatrixT_A.getNumRows(); luintidx_i++) {
    lbarray_row.setArray(aibitmatrixT_A.getRow(luintidx_i));
    mulRowsIColumns0N
      (aomatrixrowt_C.getRow(luintidx_i),
       (uintidx) 0,
       (uintidx) aimatrixrowt_B.getNumRows(),
       lbarray_row,
       aimatrixrowt_B
       );
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMatrix;
    lostrstream_labelMatrix << "<MATRIX:" << lpc_labelFunc
			       << "aomatrixrowt_C[" << &aomatrixrowt_C << ']';
    aomatrixrowt_C.print(std::cout,lostrstream_labelMatrix.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

} /*END namespace mat*/
  
#endif  /* __BITMATRIX_MATRIX_OPERATOR_HPP */
