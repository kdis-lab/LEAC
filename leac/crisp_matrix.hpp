/*! \file crisp_matrix.hpp
 *
 * \brief  crisp matrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef BIT_CRISP_MATRIX_HPP
#define BIT_CRISP_MATRIX_HPP

#include "bit_matrix.hpp"

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {
  
/*! \class CrispMatrix
  \brief Bit crisp matrix
  \details 
*/ 
template < class T_BITSIZE,
	   class T_CLUSTERIDX
	   >
class CrispMatrix: public BitMatrix<T_BITSIZE> {
public:
  CrispMatrix()
    : BitMatrix<T_BITSIZE>()
  { }
  
  CrispMatrix
  (const uintidx aiuintidx_numRows, 
   const uintidx aiuintidx_numColumns
   ): 
    BitMatrix<T_BITSIZE>(aiuintidx_numRows,aiuintidx_numColumns)
  { }
  
  CrispMatrix
  (const uintidx aiuintidx_numRows, 
   const uintidx aiuintidx_numColumns,
   T_BITSIZE    *aiarrayT_data
   ):
    BitMatrix<T_BITSIZE>(aiuintidx_numRows,aiuintidx_numColumns, aiarrayT_data)
  { }

  //copy constructor
  CrispMatrix
  (const CrispMatrix
   <T_BITSIZE,
    T_CLUSTERIDX>& aibitcrispmatrix_B
   ):
    BitMatrix<T_BITSIZE>(aibitcrispmatrix_B)
  { }

  //move constructor
  CrispMatrix
  (CrispMatrix
   <T_BITSIZE,
    T_CLUSTERIDX> &&aibitcrispmatrix_B
   )
    : BitMatrix<T_BITSIZE>(aibitcrispmatrix_B)
  { }

  virtual ~CrispMatrix() { }

  //move copy
  mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>& operator=
  (const mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aibitcrispmatrix_B)
  {
    if( this !=  &aibitcrispmatrix_B ){
      BitMatrix<T_BITSIZE>::operator=(aibitcrispmatrix_B);
    }

    return *this;
  }
 
  //move asigned
  mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>& operator=
  (CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &&aibitcrispmatrix_B)
  {
    if( this !=  &aibitcrispmatrix_B ){
      BitMatrix<T_BITSIZE>::operator=(aibitcrispmatrix_B);
    }

    return *this;
  }
   
  const T_CLUSTERIDX getMember(uintidx aiuintidx_instanceIdx) const
  {      
    for (uintidx luintidx_i = 0; luintidx_i < this->getNumRows(); luintidx_i++) {
      if ( (*this)(luintidx_i, aiuintidx_instanceIdx) ) {
	return (T_CLUSTERIDX) luintidx_i;  
      }
    }
    return (T_CLUSTERIDX) -1;
  }

  bool setMember(uintidx aiuintidx_instanceIdx, T_CLUSTERIDX aimgidxT_memberCluster)
  {
    T_CLUSTERIDX  lmgidx_prevMemberCluster;
    bool loi_threshold;

    loi_threshold = !(*this)(aimgidxT_memberCluster, aiuintidx_instanceIdx);
    if ( loi_threshold ) {
      lmgidx_prevMemberCluster = this->getMember(aiuintidx_instanceIdx);
      if (lmgidx_prevMemberCluster != -1 ) {
	this->clearBit(lmgidx_prevMemberCluster, aiuintidx_instanceIdx);
      }
      this->setBit(aimgidxT_memberCluster, aiuintidx_instanceIdx);
    }
    return loi_threshold;
  }

private:

}; /*CrispMatrix*/

} /*END namespace ds*/

#endif /*BIT_CRISP_MATRIX_HPP*/
