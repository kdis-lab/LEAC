/*! \file dist_euclidean.hpp
 *
 * \brief Euclidean distance
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __DIST_EUCLIDEAN_HPP
#define __DIST_EUCLIDEAN_HPP

#include <iostream>
#include <cmath>
#include <stdexcept>
#include "dist.hpp"
#include "matrix.hpp"
#include "instance.hpp"
#include "distance_operation.hpp"
#include "matrix_operation.hpp"
#include "linear_algebra_level1.hpp"
#include "linear_algebra_level2.hpp"
#include  "common.hpp"


/*! \namespace dist
  \brief Module for definition of distance between objects or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  dist {

/*! \struct Euclidean
  \brief Euclidean distance
*/
template < class T_DIST,
	   class T_FEATURE
	   >
struct Euclidean: public Dist<T_DIST,T_FEATURE> {
  inline T_DIST operator() (const T_FEATURE *aiarrayT_p, const T_FEATURE* aiarrayT_q, const  uintidx uintidx_length) const
  {
    return kernelEuclidean(aiarrayT_p,aiarrayT_q,uintidx_length);
  }
}; /* Dist */


/*! \struct EuclideanSquared
  \brief Euclidean distance sqrt
*/
template < class T_DIST,
	   class T_FEATURE
	   >
struct EuclideanSquared: public Dist<T_DIST,T_FEATURE> {
  inline T_DIST operator() (const T_FEATURE *aiarrayT_p, const T_FEATURE* aiarrayT_q, const uintidx uintidx_length) const 
  {
    return kernelEuclideanSquared(aiarrayT_p,aiarrayT_q,uintidx_length);
  }
}; /* EuclideanSquared */


/*INDUCE  DIST --------------------------------------------------------------------
*/  



#if  DATATYPE_CENTROIDS_ROUND == 0

/*! \struct Induced
  \brief Dist induced for real number  \cite Bezdek:ClusterAnalysis:FCM:1974 \cite Bezdek:etal:ClusterAnalysis:FCM:1984 \cite Bezdek:etal:GAclustering:GA:1994 
*/
template < class T_DIST,
	   class T_FEATURE
	   >
struct Induced: public Dist<T_DIST,T_FEATURE> {
  Induced()
    : _matrix_weight(NULL)
    , _arrayt_xt(NULL)
    , _arrayt_x(NULL)
  {}
  Induced(const mat::MatrixRow<T_DIST>& aimatrix_weight)
    : _matrix_weight(aimatrix_weight)
    , _arrayt_xt(new T_DIST[aimatrix_weight.getNumColumns()])
    , _arrayt_x(new T_DIST[aimatrix_weight.getNumColumns()])
  {}

  
 virtual  ~Induced()
  {
    if ( _arrayt_xt != NULL ) {
      delete []_arrayt_xt;
      delete []_arrayt_x;
    }
  }
 
  T_DIST operator() (const T_FEATURE *aiarrayT_p, const T_FEATURE* aiarrayT_q, const uintidx uintidx_length) const 
  {

    T_DIST  loT_dist;

    
    interfacesse::copy
      (_arrayt_x,
       aiarrayT_p,
       uintidx_length
       );

    interfacesse::axpy
      (_arrayt_x,
       T_FEATURE(-1),
       aiarrayT_q,
       uintidx_length
       );
    
    interfacesse::gemv
      (_arrayt_xt,
       _matrix_weight,
       _arrayt_x
       );
      
    loT_dist = 
      interfacesse::dot
      (_arrayt_xt,
       _arrayt_x,
       uintidx_length
       );

    return loT_dist;
    
  }

  mat::MatrixRow<T_DIST> _matrix_weight;
  T_DIST*           _arrayt_xt;
  T_DIST*           _arrayt_x;

}; /* Induced */


/*! \struct InducedSquareRoot
  \brief Dist induced sqrt for real number  \cite Bezdek:ClusterAnalysis:FCM:1974 \cite Bezdek:etal:ClusterAnalysis:FCM:1984 \cite Bezdek:etal:GAclustering:GA:1994
*/
template < class T_DIST,
	   class T_FEATURE
	   >
struct InducedSquareRoot: public Induced<T_DIST,T_FEATURE> {
  InducedSquareRoot()
    : Induced<T_DIST,T_FEATURE>()
  {
  }
  
  InducedSquareRoot(const mat::MatrixRow<T_DIST>& aimatrix_weight)
    : Induced<T_DIST,T_FEATURE>(aimatrix_weight)
  {
  }
  
 virtual  ~InducedSquareRoot()
  {
  }
 
  T_DIST operator() (const T_FEATURE *aiarrayT_p, const T_FEATURE* aiarrayT_q, const uintidx uintidx_length) const 
  {
    return std::sqrt(Induced<T_DIST,T_FEATURE>::operator()(aiarrayT_p,aiarrayT_q,uintidx_length)); 
  }
  
};/* InducedSquareRoot */

#endif //DATATYPE_CENTROIDS_ROUND

  
} /*END namespace dist 
   */

#endif /*__DIST_EUCLIDEAN_HPP*/



