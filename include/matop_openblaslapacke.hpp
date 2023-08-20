/*! \file matop_openblaslapacke.hpp
 *
 * \brief Matrix operations with OpenBlas and Lapacke
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef  MATOP_OPENBLAS_LAPACKE_HPP
#define  MATOP_OPENBLAS_LAPACKE_HPP

#include <tuple>
#include "matrix.hpp"
#include <iostream>
#include <stdlib.h>

#include "verbose_global.hpp"

#include "lapacke.h"

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {

/*! \fn float lange(char aic_norm, mat::MatrixRow<float> &aimatrixt_a)
    \brief Calculates the norm of a matrix for float
    \details Return the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absoute value of a float matrix A
    \param aimatrixt_a  a  mat::MatrixRow to calculate the norm
    \param aic_norm specifies the type of norm: 'M', '1', 'I' or 'F' 
    max(abs(A(i,j))), NORM = 'M' or 'm' 
    norm1(A),         NORM = '1', 'O' or 'o' 
    normI(A),         NORM = 'I' or 'i' 
    normF(A),         NORM = 'F', 'f', 'E' or 'e' 
    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
    normF  denotes the  Frobenius norm of a matrix (square root of sum of squares).  

    Note that  max(abs(A(i,j)))  is not a consistent matrix norm. 
*/
float lange
(char                  aic_norm, 
 mat::MatrixRow<float> &aimatrixt_a
 )
{
  const lapack_int lii_m = (lapack_int)   aimatrixt_a.getNumRows();
  const lapack_int lii_n = (lapack_int)   aimatrixt_a.getNumColumns();
  const lapack_int lii_lda = (lapack_int) aimatrixt_a.getNumColumns();

  float lolpdr_normMax =
    LAPACKE_slange
    (CblasRowMajor,
     aic_norm,
     lii_m,
     lii_n,
     aimatrixt_a.toArray(),
     lii_lda
     );
 
  return lolpdr_normMax;
}


/*! \fn float lange(char aic_norm, mat::MatrixRow<double> &aimatrixt_a)
    \brief Calculates the norm of a matrix for double
    \details Return the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absoute value of a float matrix A
    \param aimatrixt_a  a  mat::MatrixRow to calculate the norm
    \param aic_norm specifies the type of norm: 'M', '1', 'I' or 'F' 
    max(abs(A(i,j))), NORM = 'M' or 'm' 
    norm1(A),         NORM = '1', 'O' or 'o' 
    normI(A),         NORM = 'I' or 'i' 
    normF(A),         NORM = 'F', 'f', 'E' or 'e' 
    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
    normF  denotes the  Frobenius norm of a matrix (square root of sum of squares).  

    Note that  max(abs(A(i,j)))  is not a consistent matrix norm. 
*/
double lange
(char                   aic_norm, 
 mat::MatrixRow<double> &aimatrixt_a
 )
{
  const lapack_int lii_m = (lapack_int) aimatrixt_a.getNumRows();
  const lapack_int lii_n = (lapack_int) aimatrixt_a.getNumColumns();
  const lapack_int lii_lda = (lapack_int) aimatrixt_a.getNumColumns();

  float lolpdr_normMax =
    LAPACKE_dlange
    (CblasRowMajor,
     aic_norm,
     lii_m,
     lii_n,
     aimatrixt_a.toArray(),
     lii_lda
     );
  
  return lolpdr_normMax;
}


/*! \fn std::tuple<std::vector<float>,mat::MatrixRow<float> > syevd(const mat::MatrixRow<float> &aimatrixt_symatrix)
    \brief  syevd computes all eigenvalues and eigenvectors
    \details syevd computes all eigenvalues and eigenvectors of a float symmetric matrix A. If eigenvectors are desired, it uses a divide and conquer algorithm.
    \param aimatrixt_symatrix a mat::MatrixRow with the data
 */
std::tuple<std::vector<float>,mat::MatrixRow<float> >
syevd(const mat::MatrixRow<float> &aimatrixt_symatrix)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "mat::syevd";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( input const mat::MatrixRow<float>: &aimatrixt_symatrix[" << &aimatrixt_symatrix << ']'
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  
  std::vector<float>  lovector_eigenVal(aimatrixt_symatrix.getNumRows());
  mat::MatrixRow<float>    aomatrixcolumnt_eigenVec(aimatrixt_symatrix);
  char       lc_jobz = 'V';
  char       lc_uplo = 'L';
  lapack_int n = (lapack_int) aimatrixt_symatrix.getNumColumns(); 
  lapack_int lda = (lapack_int) aimatrixt_symatrix.getNumColumns();
  
  lapack_int info =
    LAPACKE_ssyev
    (CblasRowMajor,
     lc_jobz,
     lc_uplo,
     n,
     aomatrixcolumnt_eigenVec.toArray(),
     lda,
     lovector_eigenVal.data()
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelEgenVal;
    lostrstream_labelEgenVal << "<EIGENVAL: " << lpc_labelFunc << "info," << info  
			     << ",lovector_eigenVal[" << &lovector_eigenVal << ']';
    inout::containerprint
      (lovector_eigenVal.begin(),
       lovector_eigenVal.end(),
       std::cout,lostrstream_labelEgenVal.str().c_str(),
       ','
       );
    std::cout << std::endl;
    std::ostringstream lostrstream_labelEigenVec;
    lostrstream_labelEigenVec << "<EIGENVEC:" << lpc_labelFunc << "info," << info  
			       << ", aomatrixcolumnt_eigenVec[" << &aomatrixcolumnt_eigenVec << ']';
    aomatrixcolumnt_eigenVec.print(std::cout,lostrstream_labelEigenVec.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  if ( info > 0) 
    throw std::runtime_error("mat::syevd: the algorithm failed");
  
  return std::make_tuple(lovector_eigenVal,aomatrixcolumnt_eigenVec);

}


/*! \fn std::tuple<std::vector<double>,mat::MatrixRow<float> > syevd(const mat::MatrixRow<double> &aimatrixt_symatrix)
    \brief  syevd computes all eigenvalues and eigenvectors
    \details syevd computes all eigenvalues and eigenvectors of a double symmetric matrix A. If eigenvectors are desired, it uses a divide and conquer algorithm.
    \param aimatrixt_symatrix a mat::MatrixRow with the data
 */
std::tuple<std::vector<double>,mat::MatrixRow<double> >
syevd(const mat::MatrixRow<double> &aimatrixt_symatrix)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "matop_openblaslapacke.hpp::mat::syevd";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( input const mat::MatrixRow<double>: &aimatrixt_symatrix[" << &aimatrixt_symatrix << ']'
	      << ")";
    std::ostringstream lostrstream_labelsymatrix;
    lostrstream_labelsymatrix << "<SYMATRIX:" << lpc_labelFunc  
			      << ",aimatrixt_symatrix[" << &aimatrixt_symatrix << ']';
    aimatrixt_symatrix.print(std::cout,lostrstream_labelsymatrix.str().c_str(),',',';');
    std::cout << std::endl;
  }
#endif //__VERBOSE_YES
  
  
  std::vector<double>  lovector_eigenVal(aimatrixt_symatrix.getNumRows());
  mat::MatrixRow<double>    aomatrixcolumnt_eigenVec(aimatrixt_symatrix);
  char       lc_jobz = 'V';
  char       lc_uplo = 'L';
  lapack_int n = (lapack_int) aimatrixt_symatrix.getNumColumns(); 
  lapack_int lda = (lapack_int) aimatrixt_symatrix.getNumColumns();
  
  lapack_int info =
    LAPACKE_dsyev
    (CblasRowMajor,
     lc_jobz,
     lc_uplo,
     n,
     aomatrixcolumnt_eigenVec.toArray(),
     lda,
     lovector_eigenVal.data()
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelEgenVal;
    lostrstream_labelEgenVal << "<EIGENVAL: " << lpc_labelFunc << "info," << info  
			     << ",lovector_eigenVal[" << &lovector_eigenVal << ']';
    inout::containerprint
      (lovector_eigenVal.begin(),
       lovector_eigenVal.end(),
       std::cout,lostrstream_labelEgenVal.str().c_str(),
       ','
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelEigenVec;
    lostrstream_labelEigenVec
      << "<EIGENVEC:aomatrixcolumnt_eigenVec:"
      << lpc_labelFunc
      << "info," << info;
    aomatrixcolumnt_eigenVec.print(std::cout,lostrstream_labelEigenVec.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  if ( info > 0) 
    throw std::runtime_error("mat::syevd: the algorithm failed");
  
  return std::make_tuple(lovector_eigenVal,aomatrixcolumnt_eigenVec);

}

/* \fn int inverse(mat::MatrixRow<double>& aiomatrix_a)
   \brief Compute the inverse of a matrix
   \details Compute the inverse of a matrix using the LU factorization for double
   \param aiomatrix_a a mat::MatrixRow to invest
*/
int
inverse(mat::MatrixRow<double>& aiomatrix_a)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "matop_openblaslapacke.hpp::mat::inverse";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelInverse;
      lostrstream_labelInverse << "<INOUT:" << lpc_labelFunc  
			       << ":aiomatrix_a";
      aiomatrix_a.print(std::cout,lostrstream_labelInverse.str().c_str(),',',';');
      std::cout	<< std::endl;
  }
#endif //__VERBOSE_YES
  
  const lapack_int lii_m = (lapack_int) aiomatrix_a.getNumRows();
  const lapack_int lii_n = (lapack_int) aiomatrix_a.getNumColumns();
  const lapack_int lii_lda = (lapack_int) aiomatrix_a.getNumColumns();

  lapack_int *IPIV = new lapack_int[aiomatrix_a.getNumRows()+1];
  
  lapack_int INFO =
    LAPACKE_dgetrf
    (LAPACK_ROW_MAJOR,
     lii_m,
     lii_n,
     aiomatrix_a.toArray(),
     lii_lda,
     IPIV
     );

  INFO =
    LAPACKE_dgetri
    (LAPACK_ROW_MAJOR,
     lii_n,
     aiomatrix_a.toArray(),
     lii_lda,
     IPIV
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelInverse;
    lostrstream_labelInverse << "<INVERSE:" << lpc_labelFunc  
			       << ":aiomatrix_a";
    aiomatrix_a.print(std::cout,lostrstream_labelInverse.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  
  delete[] IPIV;

  return INFO;
}

} /*END namespace mat */
  
#endif /*MATOP_OPENBLAS_LAPACKE_HPP*/
