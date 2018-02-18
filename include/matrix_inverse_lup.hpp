/*! \file matrix_inverse_lup.hpp
 *
 * \brief Reverse of a matrix by the LUP method
 *
 * \details Functions LUPdecompose and  LUPinverse are part of\n
 * Copyright 2015 Chandra Shekhar (chandraiitk AT yahoo DOT co DOT in).\n
 * Homepage: https://sites.google.com/site/chandraacads\n
 * \n
 * chandraacads  is free software: you can redistribute it and/or modify\n
 * it under the terms of the GNU General Public License as published by\n
 * the Free Software Foundation, either version 3 of the License, or\n
 * any later version.\n
 * \n
 * This program is distributed in the hope that it will be useful,\n
 * but WITHOUT ANY WARRANTY; without even the implied warranty of\n
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n
 * GNU General Public License for more details.\n
 * \n
 * You should have received a copy of the GNU General Public License\n
 * along with this program. If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef MATRIX_INVERSE_LUP_HPP
#define MATRIX_INVERSE_LUP_HPP

#include <stdio.h>
#include "matrix.hpp"
#include "vector_utils.hpp"


/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {

/* This function decomposes the matrix 'A' into L, U, and P. If successful,
 * the L and the U are stored in 'A', and information about the pivot in 'P'.
 * The diagonal elements of 'L' are all 1, and therefore they are not stored. */
template <class T_FEATURE>
int LUPdecompose(mat::MatrixRow<T_FEATURE> &aiomatrixrt_a, uintidx *aovectorui_p)
{
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "matrix_inverse_lup.hpp::mat::LUPdecompose";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
		<< "(inout mat::MatrixRow<T_FEATURE>: &aiomatrixrt_a[" << &aiomatrixrt_a << ']'
		<< "\n(inout uintidx *aovectorui_p[" << aovectorui_p << ']'
		<< "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  //int i, j, k, kd = 0, T;
  //sankhya p, t;

  int lio_info = 0;
  
  uintidx lui_kd = 0;
  /* Finding the pivot of the LUP decomposition. */
  for(uintidx lui_i=0; lui_i<aiomatrixrt_a.getNumRows(); lui_i++) aovectorui_p[lui_i] = lui_i; //Initializing.

  for(uintidx lui_k=0; lui_k<aiomatrixrt_a.getNumRows()-1; lui_k++)
    {
      T_FEATURE p = 0;
      for(uintidx lui_i=lui_k; lui_i<aiomatrixrt_a.getNumRows(); lui_i++)
	{
	  T_FEATURE lrt_t= aiomatrixrt_a(lui_i,lui_k);
	  if(lrt_t < 0) lrt_t*= -1; //aiomatrixrt_abosolute value of 't'.
	  if(lrt_t > p)
	    {
	      p = lrt_t;
	      lui_kd = lui_i;
	    }
	}

      if(p == 0)
	{
	  printf("\nLUPdecompose(): ERROR: aiomatrixrt_a singular matrix is supplied.\n"\
		 "\tRefusing to proceed any further.\n");
	  lio_info = -1;
	  break;
	  //return -1;
	}

      /* Exchanging the rows according to the pivot determined above. */
      uintidx lui_t = aovectorui_p[lui_kd];
      aovectorui_p[lui_kd] = aovectorui_p[lui_k];
      aovectorui_p[lui_k] = lui_t;
      for(uintidx lui_i=0; lui_i<aiomatrixrt_a.getNumRows(); lui_i++)
	{
	  T_FEATURE t = aiomatrixrt_a(lui_kd,lui_i);
	  aiomatrixrt_a(lui_kd,lui_i) = aiomatrixrt_a(lui_k,lui_i);
	  aiomatrixrt_a(lui_k,lui_i) = t;
	}

      //Performing substraction to decompose aiomatrixrt_a as LU.
      for(uintidx lui_i=lui_k+1; lui_i<aiomatrixrt_a.getNumRows(); lui_i++) 
	{
	  aiomatrixrt_a(lui_i,lui_k) = aiomatrixrt_a(lui_i,lui_k)/aiomatrixrt_a(lui_k,lui_k);
	  for(uintidx lui_j=lui_k+1; lui_j<aiomatrixrt_a.getNumRows(); lui_j++)
	    aiomatrixrt_a(lui_i,lui_j) -= aiomatrixrt_a(lui_i,lui_k)*aiomatrixrt_a(lui_k,lui_j);
	}
    } //Now, 'aiomatrixrt_a' contains the L (without the diagonal elements, which are all 1)
  //and the U.
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelP;
    lostrstream_labelP << "<P: " << lpc_labelFunc << "info," << lio_info  
		       << ",aovectorui_p[" << &aovectorui_p << ']';
    inout::containerprint
      (aovectorui_p,
       aovectorui_p + aiomatrixrt_a.getNumRows(),
       std::cout,
       lostrstream_labelP.str().c_str(),
       ','
       );
    std::ostringstream lostrstream_labelLU;
    lostrstream_labelLU << "<LU:" << lpc_labelFunc << "info," << lio_info 
			       << ", aiomatrixrt_a[" << &aiomatrixrt_a << ']';
    aiomatrixrt_a.print(std::cout,lostrstream_labelLU.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lio_info;
  
}

/* This function calculates the inverse of the LUP decomposed matrix 'LU' and pivoting
 * information stored in 'P'. The inverse is returned through the matrix 'LU' itselt.
*/
template <class T_FEATURE>
int LUPinverse(mat::MatrixRow<T_FEATURE> &aiomatrixrt_lu, uintidx *aivectorui_p)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "matrix_inverse_lup.hpp::mat::LUPinverse";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelInverse;
      lostrstream_labelInverse << "<INOUT:" << lpc_labelFunc  
			       << ":aiomatrixrt_lu";
      aiomatrixrt_lu.print(std::cout,lostrstream_labelInverse.str().c_str(),',',';');
      std::cout	<< "\ninout uintidx *aivectorui_p[" << aivectorui_p << ']'
		<< "\n)"
		<< std::endl;
  }
#endif //__VERBOSE_YES
  
  /* 'lmatrix_b', 'lvectorrt_x', and 'lvectorrt_y' are used as temporary spaces. 
   */
  mat::MatrixRow<T_FEATURE> lmatrix_b(aiomatrixrt_lu.getNumRows(),aiomatrixrt_lu.getNumRows());
  T_FEATURE *lvectorrt_x = new T_FEATURE[aiomatrixrt_lu.getNumRows()]; 
  T_FEATURE *lvectorrt_y = new T_FEATURE[aiomatrixrt_lu.getNumRows()];
  //int i, j, n, m;
  //sankhya t;

  //Initializing lvectorrt_x and lvectorrt_y.
  for(uintidx n=0; n<aiomatrixrt_lu.getNumRows(); n++) lvectorrt_x[n] = lvectorrt_y[n] = 0;

  /* Solving LUX = Pe, in order to calculate the inverse of 'A'. Here, 'e' is a column
   * vector of the identity matrix of size 'size-1'. Solving for all 'e'. */
  for(uintidx i=0; i<aiomatrixrt_lu.getNumRows(); i++)
    {
      //Storing elements of the i-th column of the identity matrix in i-th row of 'B'.
      for(uintidx j = 0; j<aiomatrixrt_lu.getNumRows(); j++)
  	lmatrix_b(i,j) = 0;
      lmatrix_b(i,i) = 1;

      //Solving Ly = Pb.
      for(uintidx n=0; n<aiomatrixrt_lu.getNumRows(); n++)
  	{
  	  T_FEATURE lrt_t = 0;
  	  for(uintidx m=0; (n != 0) && (m <= n-1); m++) lrt_t += aiomatrixrt_lu(n,m)*lvectorrt_y[m];
  	  lvectorrt_y[n] = lmatrix_b(i,aivectorui_p[n])-lrt_t;
  	}
      
      //Solving Ux = y.
      //for(uintidx n=aiomatrixrt_lu.getNumRows()-1; n >= 0; n--)
      for(uintidx n=aiomatrixrt_lu.getNumRows()-1; n < aiomatrixrt_lu.getNumRows(); n--)
      	{
      	    T_FEATURE lrt_t = 0;
	    for(uintidx m = n+1; m < aiomatrixrt_lu.getNumRows(); m++)
	      lrt_t += aiomatrixrt_lu(n,m)*lvectorrt_x[m];
	    lvectorrt_x[n] = (lvectorrt_y[n]-lrt_t)/aiomatrixrt_lu(n,n);
      	}//Now, lvectorrt_x contains the solution.
     
      //Copying 'lvectorrt_x' into the same row of 'B'.
      for(uintidx j = 0; j<aiomatrixrt_lu.getNumRows(); j++)
    	lmatrix_b(i,j) = lvectorrt_x[j]; 
    } //Now, 'B' the transpose of the inverse of 'A'.

  /* Copying transpose of 'B' into 'LU', which would the inverse of 'A'. */
  for(uintidx  i=0; i<aiomatrixrt_lu.getNumRows(); i++)
    for(uintidx j=0; j<aiomatrixrt_lu.getNumRows(); j++)
      aiomatrixrt_lu(i,j) = lmatrix_b(j,i);

  delete [] lvectorrt_x;
  delete [] lvectorrt_y;

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelP;
    lostrstream_labelP << "<P:" << lpc_labelFunc  
		       << ",aivectorui_p[" << &aivectorui_p << ']';
    inout::containerprint
      (aivectorui_p,
       aivectorui_p + aiomatrixrt_lu.getNumRows(),
       std::cout,lostrstream_labelP.str().c_str(),
       ','
       );
    std::cout << '\n';
    std::ostringstream lostrstream_labelInverse;
    lostrstream_labelInverse << "<INVERSE:" << lpc_labelFunc  
			       << ":aiomatrixrt_lu";
    aiomatrixrt_lu.print(std::cout,lostrstream_labelInverse.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return 0;
}
  
} /*END namespace mat*/
  
#endif /* MATRIX_INVERSE_LUP_HPP */
