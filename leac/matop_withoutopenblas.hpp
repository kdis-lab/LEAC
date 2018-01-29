/*! \file matop_withoutopenblas.hpp
 *
 * \brief Matrix operations
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef  MATOP_HPP
#define  MATOP_HPP

#include <vector>
#include <tuple>
#include <ctype.h>
#include <cmath>


#include "matrix_inverse_lup.hpp"
#include "jacobi_eigenvalue.hpp"
#include "interface_level1.hpp"

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {


/*! \fn T lange(char aic_norm, mat::MatrixRow<T> &aimatrixt_a)
    \brief Calculates the norm of a matrix for T
    \details Return the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absoute value of a T matrix A
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
template <class T>
T lange
(const char              aic_norm,
 const mat::MatrixRow<T> &aimatrixt_a
)
{
  T lt_value(std::numeric_limits<T>::min());

  switch (aic_norm) {
  case 'm':
  case 'M':
    {
      const T* larray_matrix = aimatrixt_a.toArray();
      const T* larray_matrixEnd = larray_matrix + aimatrixt_a.getNumElems();
      while ( larray_matrix != larray_matrixEnd ) {
	T lt_temp = std::abs(*larray_matrix);
	if ( lt_temp >  lt_value)
	  lt_value = lt_temp;
	++larray_matrix;
      }  
    }
    break;
  case '1':
  case 'O':
  case 'o':
    {
      //norm1 denotes the  one norm of a matrix (maximum column sum)
      std::vector<T> lvector_tmp(aimatrixt_a.getNumColumns(),T(0));
      for (uintidx lui_i = 0; lui_i < aimatrixt_a.getNumRows(); lui_i++) {
	const T *lptt_row = aimatrixt_a.getRow(lui_i);
	for ( auto liter_col = lvector_tmp.begin();
	      liter_col != lvector_tmp.end();
	      ++liter_col, ++lptt_row)
	  {
	    *liter_col += std::abs(*lptt_row);
	  }
      }
      lt_value = *std::max_element(lvector_tmp.begin(),lvector_tmp.end());  
    }
    break;
  case 'I':
  case 'i':
    {
      //normI  denotes the  infinity norm  of a matrix  (maximum row sum)
      for (uintidx lui_i = 0; lui_i < aimatrixt_a.getNumRows(); lui_i++) {
	T lt_temp =
	  std::accumulate
	  (aimatrixt_a.getRow(lui_i),
	   aimatrixt_a.getRow(lui_i)+aimatrixt_a.getNumColumns(),
	   T(0),
	   [](T ait_sum, T ait_elem) {
	    return ait_sum + std::abs(ait_elem); 
	  }
	   );
	if ( lt_temp >  lt_value)
	  lt_value = lt_temp;
      }
    }
    break;
  case 'F':
  case 'f':
  case 'E':
  case 'e':
    {
      //Frobenius norm of a matrix (square root of sum of squares)
      lt_value = 0;
      const T *larray_matrix = aimatrixt_a.toArray();
      const T *larray_matrixEnd = larray_matrix + aimatrixt_a.getNumElems();
      while ( larray_matrix != larray_matrixEnd ) {
	lt_value +=   (*larray_matrix) * (*larray_matrix);
	++larray_matrix;
      }
      lt_value = std::sqrt(lt_value);
    }
    break;   
  default:
    break;
  }

  return lt_value;
}
  
/*! \fn int getrf(mat::MatrixRow<T>& aiomatrix_a)
  \brief Computes an L * U * P factorization of a general matrix aiomatrix_a
  \details This function is an interface to LUPdecompose. 
  If successful, the L and the U are stored in aiomatrix_a, and information 
  about the pivot in aovectorui_p. The diagonal elements of 'L' are all 1, 
  and therefore they are not stored.
  \param aiomatrix_a a mat::MatrixRow<T>
  \param aovectorui_p a vector with information about the pivot.
*/
template <class T>
int getrf(mat::MatrixRow<T>& aiomatrix_a, uintidx *aovectorui_p)
{
  int INFO =  LUPdecompose(aiomatrix_a,aovectorui_p);
  
  return INFO;
}


/*! \fn std::tuple<std::vector<T>,mat::MatrixRow<T> > syevd(const mat::MatrixRow<T> &aimatrixt_symatrix)
    \brief  syevd computes all eigenvalues and eigenvectors
    \details syevd computes all eigenvalues and eigenvectors of a T symmetric matrix A. If eigenvectors are desired, it uses a divide and conquer algorithm.
    \param aimatrixt_symatrix a mat::MatrixRow with the data
 */
template <class T>
std::tuple<std::vector<T>,mat::MatrixRow<T> >
syevd(const mat::MatrixRow<T> &aimatrixt_symatrix)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "matop_withoutopenblas.hpp::mat::syevd";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( input const mat::MatrixRow<T>: &aimatrixt_symatrix["
		<< &aimatrixt_symatrix << ']'
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  mat::MatrixRow<T> lmatrixt_symatrix(aimatrixt_symatrix);
  std::vector<T>    lovector_eigenVal(aimatrixt_symatrix.getNumRows());
  mat::MatrixRow<T> lomatrix_eigenVec
    (aimatrixt_symatrix.getNumRows(),aimatrixt_symatrix.getNumRows());
    
  int li_itmax = 100;
  int li_n     = (int) aimatrixt_symatrix.getNumRows();
  int li_rotnum;
  int li_itnum;
  
  jacobi_eigenvalue
    (li_n,
     lmatrixt_symatrix.toArray(),
     li_itmax,
     lomatrix_eigenVec.toArray(),
     lovector_eigenVal.data(),
     li_itnum,
     li_rotnum
     );

  lomatrix_eigenVec = getTranspose(lomatrix_eigenVec);
  
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelEgenVal;
    lostrstream_labelEgenVal << "<EIGENVAL: " << lpc_labelFunc 
			     << ",lovector_eigenVal[" << &lovector_eigenVal << ']';
    inout::containerprint
      (lovector_eigenVal.begin(),
       lovector_eigenVal.end(),
       std::cout,
       lostrstream_labelEgenVal.str().c_str(),
       ','
       );
    std::cout << std::endl;
    std::ostringstream lostrstream_labelEigenVec;
    lostrstream_labelEigenVec << "<EIGENVEC:" << lpc_labelFunc  
			       << ", lomatrix_eigenVec[" << &lomatrix_eigenVec << ']';
    lomatrix_eigenVec.print(std::cout,lostrstream_labelEigenVec.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return std::make_tuple(lovector_eigenVal,lomatrix_eigenVec);

}


/*! \fn int inverse(mat::MatrixRow<T>& aiomatrix_a)
    \brief Compute the inverse of a matrix
    \details Compute the inverse of a matrix using the LU factorization
    \param aiomatrix_a a mat::MatrixRow to invest
 */
template <class T>
int inverse(mat::MatrixRow<T>& aiomatrix_a)
{
  uintidx *lvectorui_p = new uintidx[aiomatrix_a.getNumRows()];

  int INFO =  LUPdecompose(aiomatrix_a,lvectorui_p);

  if ( INFO == 0 )
    LUPinverse(aiomatrix_a,lvectorui_p);
  
  delete[] lvectorui_p;

  return INFO;
  
}

} /*END namespace mat */

#endif /*MATOP_HPP*/
