/*! \file interface_level2.hpp
 *
 * \brief Interface_level 2 
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INTERFACE_LEVEL2_HPP
#define INTERFACE_LEVEL2_HPP

#include <stdexcept>
#include "common_interfacelapack.hpp"
#include "matrix.hpp"


/*! \namespace interfacesse
  \brief Interface to Streaming SIMD Extensions (SSE) high Performance Level Low Functions
  \details Functions base in SSE instructions operate on packed double-precision floating-point values contained in XMM registers and on packed integers contained in MMX and XMM registers \cite progguide:intel10

  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace interfacesse {

  
/*! \fn void gemv(double *aoarrayt_y, const mat::MatrixRow<double>  &aimatrixrowt_a, const double *aiarrayt_x)
  \brief  gemv of blas operation 
  \details
  \param aoarrayt_y a array of  number double
  \param aimatrixrowt_a a const mat::MatrixRow<double>
  \param aiarrayt_x a array of  number double
*/

template < typename T>
void
gemv
(T                        *aoarrayt_y,
 const mat::MatrixRow<T>  &aimatrixrowt_a,
 const T                  *aiarrayt_x
 )
{
  const uintidx lii_m =  aimatrixrowt_a.getNumRows();
  const uintidx lii_n =  aimatrixrowt_a.getNumColumns();
  for ( uintidx lui_i = 0; lui_i < lii_m;  lui_i++) {
    const T *lt_rowarray_a = aimatrixrowt_a.getRow(lui_i); 
    *aoarrayt_y = std::inner_product(lt_rowarray_a,lt_rowarray_a+lii_n, aiarrayt_x, T(0));
    ++aoarrayt_y;
  }
}


/*! \fn void gemm(mat::MatrixRow<double> &aomatrixrowt_C, const mat::MatrixRow<double> &aimatrixrowt_A, const mat::MatrixRow<double> &aimatrixrowt_B, const enum CBLAS_TRANSPOSE aienum_transA = CblasNoTrans, const enum CBLAS_TRANSPOSE aienum_transB = CblasNoTrans, const double airt_alpha = 1.0, const double airt_beta  = 0.0)
  \brief  gemm of blas operation 
  \details
\f[
C \leftarrow \alpha op(A)op(B) +  \beta C
\f]
  \param aomatrixrowt_C a mat::MatrixRow<double>
  \param aimatrixrowt_A a const mat::MatrixRow<double>
  \param aienum_transA a enum CBLAS_TRANSPOSE default No transpose
  \param aienum_transB a  enum CBLAS_TRANSPOSE default No transpose
  \param airt_alpha a real number double 
  \param airt_beta a real number double 
*/  
template < typename T>
void
gemm
(mat::MatrixRow<T>         &aomatrixrowt_C,
 const mat::MatrixRow<T>   &aimatrixrowt_A,
 const mat::MatrixRow<T>   &aimatrixrowt_B,
 const enum mat::TRANSPOSE aienum_transA = mat::CblasNoTrans,
 const enum mat::TRANSPOSE aienum_transB = mat::CblasNoTrans,
 const T                   airt_alpha = 1.0,
 const T                   airt_beta  = 0.0
 )
{
  if ( (aienum_transA == mat::CblasNoTrans) && (aienum_transB == mat::CblasNoTrans) ) {
    if (aimatrixrowt_A.getNumColumns() != aimatrixrowt_B.getNumRows()) { // A B 
      std::string lstr_error("gemm: Number columns of matrix A are different from rows of B");
      //lstr_error += aistr_fileInstance;
      throw  std::invalid_argument(lstr_error);
    }
    for (uintidx li_i = 0; li_i < aomatrixrowt_C.getNumRows(); li_i++) {
      T *ltpt_c =  aomatrixrowt_C.getRow(li_i);
      //T *ltpt_a =  aimatrixrowt_A.getRow(li_i);
      
      for (uintidx li_j = 0; li_j < aomatrixrowt_C.getNumColumns(); li_j++) {
	//*ltpt_c *= airt_beta;
	T sum = 0.;
	const T *ltpt_a =  aimatrixrowt_A.getRow(li_i);
	for (uintidx li_k = 0; li_k < aimatrixrowt_A.getNumColumns(); li_k++) {
	  sum += *ltpt_a * aimatrixrowt_B(li_k,li_j);
	  ++ltpt_a;
	  //sum += aimatrixrowt_A(li_i,li_k) * aimatrixrowt_B(li_k,li_j);
	}
	*ltpt_c *= airt_beta;
	*ltpt_c += airt_alpha* sum;
	++ltpt_c;
	//aomatrixrowt_C(li_i,li_j) = sum;
	/*const T *ltpt_cEnd = ltpt_c + aomatrixrowt_C.getNumRows();
	  while (ltpt_c != ltpt_cEnd ) {
	
	++ltpt_c;
	}*/
      }
    }
  }
  else if ( (aienum_transA == mat::CblasNoTrans) && (aienum_transB == mat::CblasTrans) ) {
    if (aimatrixrowt_A.getNumColumns() != aimatrixrowt_B.getNumColumns()) { //A B'
      std::string lstr_error("gemm: Number columns of matrix A are different from Columns of B'");
      //lstr_error += aistr_fileInstance;
      throw  std::invalid_argument(lstr_error);
    }
    for (uintidx li_i = 0; li_i < aomatrixrowt_C.getNumRows(); li_i++) {
      T *ltpt_c =  aomatrixrowt_C.getRow(li_i);  
      for (uintidx li_j = 0; li_j < aomatrixrowt_C.getNumColumns(); li_j++) {
	//*ltpt_c *= airt_beta;
	T sum = 0.;
	const T *ltpt_a =  aimatrixrowt_A.getRow(li_i);
	const T *ltpt_b =  aimatrixrowt_B.getRow(li_j);
	for (uintidx li_k = 0; li_k < aimatrixrowt_A.getNumColumns(); li_k++) {
	  sum += *ltpt_a * (*ltpt_b);
	  ++ltpt_a;
	  ++ltpt_b;
	  //sum += aimatrixrowt_A(li_i,li_k) * aimatrixrowt_B(li_k,li_j);
	}
	*ltpt_c *= airt_beta;
	*ltpt_c += airt_alpha* sum;
	++ltpt_c;
	//aomatrixrowt_C(li_i,li_j) = sum;
	/*const T *ltpt_cEnd = ltpt_c + aomatrixrowt_C.getNumRows();
	  while (ltpt_c != ltpt_cEnd ) {
	
	++ltpt_c;
	}*/
      }
    }
  } // (aienum_transA == mat::CblasNoTrans) && (aienum_transB == mat::CblasTrans)
  else if ( (aienum_transA == mat::CblasTrans) && (aienum_transB == mat::CblasTrans) ) {
    if (aimatrixrowt_A.getNumRows() != aimatrixrowt_B.getNumColumns()) { //A' B'
      std::string lstr_error("gemm: Number columns of matrix A are different from Columns of B'");
      throw  std::invalid_argument(lstr_error);
    }
    for (uintidx li_i = 0; li_i < aomatrixrowt_C.getNumRows(); li_i++) {
      T *ltpt_c =  aomatrixrowt_C.getRow(li_i);  
      for (uintidx li_j = 0; li_j < aomatrixrowt_C.getNumColumns(); li_j++) {
	//*ltpt_c *= airt_beta;
	T sum = 0.;
	//const T *ltpt_a =  aimatrixrowt_A.getRow(li_i);
	const T *ltpt_b =  aimatrixrowt_B.getRow(li_j);
	for (uintidx li_k = 0; li_k < aimatrixrowt_A.getNumRows(); li_k++) {
	  //sum += *ltpt_a * (*ltpt_b);
	  //++ltpt_a;
	  //
	  sum += aimatrixrowt_A(li_k,li_i) *(*ltpt_b); // aimatrixrowt_B(li_j,li_k);
	  ++ltpt_b;
	  //sum += aimatrixrowt_A(li_i,li_k) * aimatrixrowt_B(li_k,li_j);
	}
	*ltpt_c *= airt_beta;
	*ltpt_c += airt_alpha* sum;
	++ltpt_c;
	//aomatrixrowt_C(li_i,li_j) = sum;
	/*const T *ltpt_cEnd = ltpt_c + aomatrixrowt_C.getNumRows();
	  while (ltpt_c != ltpt_cEnd ) {
	
	++ltpt_c;
	}*/
      }
    }
  } 
  else if ( (aienum_transA == mat::CblasTrans) && (aienum_transB == mat::CblasNoTrans) ) {
    if (aimatrixrowt_A.getNumRows() != aimatrixrowt_B.getNumRows()) { //A' B
      std::string lstr_error("gemm: Number columns of matrix A are different from Columns of B'");
      throw  std::invalid_argument(lstr_error);
    }
    for (uintidx li_i = 0; li_i < aomatrixrowt_C.getNumRows(); li_i++) {
      T *ltpt_c =  aomatrixrowt_C.getRow(li_i);  
      for (uintidx li_j = 0; li_j < aomatrixrowt_C.getNumColumns(); li_j++) {
	//*ltpt_c *= airt_beta;
	T sum = 0.;
	//const T *ltpt_a =  aimatrixrowt_A.getRow(li_i);
	//const T *ltpt_b =  aimatrixrowt_B.getRow(li_j);
	for (uintidx li_k = 0; li_k < aimatrixrowt_A.getNumRows(); li_k++) {
	  //sum += *ltpt_a * (*ltpt_b);
	  //++ltpt_a;
	  //
	  sum += aimatrixrowt_A(li_k,li_i) * aimatrixrowt_B(li_k,li_j);
	  //++ltpt_b;
	  //sum += aimatrixrowt_A(li_i,li_k) * aimatrixrowt_B(li_k,li_j);
	}
	*ltpt_c *= airt_beta;
	*ltpt_c += airt_alpha* sum;
	++ltpt_c;
	//aomatrixrowt_C(li_i,li_j) = sum;
	/*const T *ltpt_cEnd = ltpt_c + aomatrixrowt_C.getNumRows();
	  while (ltpt_c != ltpt_cEnd ) {
	
	++ltpt_c;
	}*/
      }
    }
  }

  
}
  
} /*END namespace interfacesse*/


#endif  /* INTERFACE_LEVEL2_HPP */
