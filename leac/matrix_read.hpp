/*! \file matrix_read.hpp
 *
 * \brief read a matrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MATRIX_READ_HPP
#define MATRIX_READ_HPP

#include <iostream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <utility>      // std::pair, std::make_pair
#include "line_split.hpp"
#include "matrix.hpp"
#include "matrix_withrownull.hpp"


/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {

template <class T_FEATURE>
mat::MatrixRow<T_FEATURE>
matrixread_get(const std::string &aistr_matrix, const bool aib_homogeneousCoord = false)
{
  using namespace std;
  istringstream liss_stringstream;

#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "matrixread_get  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\ninput  std::string &aistr_matrix[" << &aistr_matrix << ']' 
	      <<  '\n' << aistr_matrix
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  inout::LineSplit
    lls_row
    (std::string(";"),
     std::string("")
     );

  inout::LineSplit
    lls_columns
    (std::string(","),
     std::string("")
     );

  std::string lstr_lineData;
  std::size_t lst_found =  aistr_matrix.find_last_of('>');
  if (lst_found != std::string::npos) {
    lstr_lineData = aistr_matrix.substr(lst_found+1);
  }
  else {
    lstr_lineData = aistr_matrix;
  }
  const uintidx luintidx_numRowsNull = lls_row.split(lstr_lineData);
  uintidx luintidx_numRows    = 0;
  uintidx luintidx_numColumns = 0;

  if ( luintidx_numRowsNull > 0 ) {
    for (uintidx luintidx_i = 1; luintidx_i <= luintidx_numRowsNull; luintidx_i++) {
      std::string lstr_matrixRow = lls_row.getItem(luintidx_i);
      if ( lstr_matrixRow.at(0) != MATRIX_WITHROWNULL_CHAR ) {
	++luintidx_numRows;
	if ( luintidx_numColumns == 0 ) {
	  luintidx_numColumns = lls_columns.split(lstr_matrixRow);
	}
      }
    } 
  }
  /*cout << "luintidx_numRows " << luintidx_numRows << std::endl;
    cout << "luintidx_numColumns " << luintidx_numColumns << std::endl;
  */
					    
  mat::MatrixRow<T_FEATURE> lomatrixrow_read
    (luintidx_numRows,
     (aib_homogeneousCoord)?luintidx_numColumns +1 : luintidx_numColumns
     );
  // (lls_row.split(lstr_lineData),
  //  lls_columns.split(lls_row.getItem(1))
  //  );

  if ( luintidx_numRows > 0 ) {
    /*T *larrayrow_data = lomatrixrow_read.getRowCol(0);
      for (uintidx luintidx_j = 1; luintidx_j <= luintidx_numColumns; luintidx_j++) {
      liss_stringstream.clear();
      liss_stringstream.str(lls_columns.getItem(luintidx_j));
      liss_stringstream >>  larrayrow_data[luintidx_j-1]; 
      }
      if ( aib_homogeneousCoord ) {
      larrayrow_data[luintidx_numColumns] = T(1);
      }*/

    uintidx luintidx_k=0;
    for(uintidx luintidx_i=1; luintidx_i <= luintidx_numRowsNull; luintidx_i++) {
      std::string lstr_matrixRow = lls_row.getItem(luintidx_i);
      if ( lstr_matrixRow.at(0) != MATRIX_WITHROWNULL_CHAR ) {
	lls_columns.split(lstr_matrixRow);
	T_FEATURE* larrayrow_data = lomatrixrow_read.getRow(luintidx_k);
	for (uintidx luintidx_j = 1; luintidx_j <= luintidx_numColumns; luintidx_j++) {
	  liss_stringstream.clear();
	  liss_stringstream.str(lls_columns.getItem(luintidx_j));
	  liss_stringstream >>  larrayrow_data[luintidx_j-1]; 
	}
	if (aib_homogeneousCoord) {
	  larrayrow_data[luintidx_numColumns] = T_FEATURE(1);
	}
	++luintidx_k;
      }
    } /*for*/
  }
 
#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "matrixread_get OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\noutput mat::MatrixRow<T_FEATURE> lomatrixrow_read [" << &lomatrixrow_read << ']' 
	      << '\n' << lomatrixrow_read
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  
  return lomatrixrow_read;

}

} /*END namespace mat*/

#endif /*MATRIX_READ_HPP*/
