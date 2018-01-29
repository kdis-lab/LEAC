/*! \file  interface_sse_level2_64bits.hpp
 *
 * \brief interface sse level2 64bits
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INTERFACE_CLAPACK_LEVEL2_HPP
#define INTERFACE_CLAPACK_LEVEL2_HPP

#include "common_interfacelapack.hpp"
#include "matrix.hpp"

/*! \namespace interfacesse
  \brief Interface to Streaming SIMD Extensions (SSE) kernel function
  \details Functions base in SSE instructions operate on packed double-precision floating-point values contained in XMM registers and on packed integers contained in MMX and XMM registers \cite progguide:intel10.

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
inline
void
gemv
(double                        *aoarrayt_y ,
 const mat::MatrixRow<double>  &aimatrixrowt_a,
 const double                  *aiarrayt_x
 )
{
  const int64_t lii_inc = (int64_t) 1;
  const int64_t lii_m = (int64_t) aimatrixrowt_a.getNumRows();
  const int64_t lii_n = (int64_t) aimatrixrowt_a.getNumColumns();
  const double  lrt_alpha = 1.0;
  const double  lrt_beta = 0.0;

  cblas_dgemv
    (CblasRowMajor,
     CblasNoTrans,
     lii_m,
     lii_n,
     lrt_alpha,
     aimatrixrowt_a.toArray(),
     lii_n,
     aiarrayt_x,
     lii_inc,
     lrt_beta,
     aoarrayt_y,
     lii_inc
     );
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
inline
void
gemm
(mat::MatrixRow<double>       &aomatrixrowt_C,
 const mat::MatrixRow<double> &aimatrixrowt_A,
 const mat::MatrixRow<double> &aimatrixrowt_B,
 const enum mat::TRANSPOSE    aienum_transA = mat::CblasNoTrans,
 const enum mat::TRANSPOSE    aienum_transB = mat::CblasNoTrans,
 const double                 airt_alpha = 1.0,
 const double                 airt_beta  = 0.0
 )
{
  const int lii_ka = (aienum_transA == mat::CblasNoTrans)
    ?(int) aimatrixrowt_A.getNumColumns()
    :(int) aimatrixrowt_A.getNumRows();
  const int lii_kb = (aienum_transB == mat::CblasNoTrans)
    ?(int) aimatrixrowt_B.getNumRows()
    :(int) aimatrixrowt_B.getNumColumns();
  enum CBLAS_TRANSPOSE lienum_transA = (CBLAS_TRANSPOSE) aienum_transA;
  enum CBLAS_TRANSPOSE lienum_transB = (CBLAS_TRANSPOSE) aienum_transB;
    
  if ( lii_ka != lii_kb )  
    throw  std::range_error
    ("Matrix<double>::operator*: order of the matrices is different");
  
  const int64_t lii_m = (int64_t)  aomatrixrowt_C.getNumRows(); 
  const int64_t lii_n = (int64_t)  aomatrixrowt_C.getNumColumns();
  
  const int64_t lii_lda = (int64_t) aimatrixrowt_A.getNumColumns();
  const int64_t lii_ldb = (int64_t) aimatrixrowt_B.getNumColumns(); 
  const int64_t lii_ldc = (int64_t) aomatrixrowt_C.getNumColumns(); 

  cblas_dgemm
    (CblasRowMajor,
     lienum_transA,
     lienum_transB,
     lii_m,
     lii_n,
     lii_ka,
     airt_alpha,
     aimatrixrowt_A.toArray(),
     lii_lda,
     aimatrixrowt_B.toArray(),
     lii_ldb,
     airt_beta,
     aomatrixrowt_C.toArray(),
     lii_ldc
     );  
}

} /*END namespace interfacesse*/

#endif  /* INTERFACE_CLAPACK_LEVEL2_HPP */
