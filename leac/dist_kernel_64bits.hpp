/*! \file dist_kernel_64bits.hpp
 *
 * \brief distance kernel 64bits 
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef DIST_KERNEL_HPP
#define DIST_KERNEL_HPP

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "common_ssekernel_64bits.hpp"


/*! \namespace dist
  \brief Module for definition of distance between objects or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  dist {

/*! \fn double kernelEuclidean(const double* aiarrayd_p, const double* aiarrayd_q, const uintidx uintidx_length)
  \brief Gets the Euclidean distance of two double n-dimensional points with SSE
  \details
  \param aiarrayd_p a double array 
  \param aiarrayd_q a double array
  \param uintidx_length  a length of the array
 */
inline
double
kernelEuclidean
(const double* aiarrayd_p, const double* aiarrayd_q, const uintidx uintidx_length)
{
  double loT_dist = 0.0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
  
  loT_dist = ddnrm2_k(lit_n, aiarrayd_p, lit_inc, aiarrayd_q, lit_inc);

  return  std::sqrt(loT_dist);  
}


/*! \fn double kernelEuclideanSquared(const double* aiarrayd_p, const double* aiarrayd_q, const uintidx uintidx_length)
  \brief Gets the square Euclidean distance of two double n-dimensional points with SSE
  \details
  \param aiarrayd_p a double array 
  \param aiarrayd_q a double array
  \param uintidx_length length of the array 
 */
inline
double
kernelEuclideanSquared
(const double* aiarrayd_p, const double* aiarrayd_q, const uintidx uintidx_length)
{
  double loT_dist = 0.0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
  
  loT_dist = ddnrm2_k(lit_n, aiarrayd_p, lit_inc, aiarrayd_q, lit_inc);

  return  loT_dist;  
} 


/*! \fn float kernelEuclidean(const float* aiarrayd_p, const float* aiarrayd_q, const  uintidx uintidx_length)
  \brief Gets the Euclidean distance of two float n-dimensional points with SSE
  \details
  \param aiarrayd_p a float array 
  \param aiarrayd_q a float array
  \param uintidx_length length of the array
 
 */
inline
float
kernelEuclidean
(const float* aiarrayd_p, const float* aiarrayd_q, const  uintidx uintidx_length)
{
  float   loT_dist = 0.0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
 

  loT_dist = ssnrm2_k(lit_n, aiarrayd_p, lit_inc, aiarrayd_q, lit_inc);

  return  std::sqrt(loT_dist);
  
} 
 

/*! \fn float kernelEuclideanSquared(const float* aiarrayf_p, const float* aiarrayf_q, const uintidx uintidx_length)
  \brief Gets the square Euclidean distance of two float n-dimensional points with SSE
  \details
  \param aiarrayf_p a float array 
  \param aiarrayf_q a float array
  \param uintidx_length length of the array
 */
inline
float
kernelEuclideanSquared
(const float* aiarrayf_p, const float* aiarrayf_q, const uintidx uintidx_length)
{
  float   loT_dist = 0.0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
 

  loT_dist = ssnrm2_k(lit_n, aiarrayf_p, lit_inc, aiarrayf_q, lit_inc);

  return  loT_dist;  
} 

  
/*! \fn double kernelEuclidean(const int* aiarrayi_p, const int* aiarrayi_q, const uintidx uintidx_length)
  \brief Kernel euclidean
  \details Gets the distance of two n-dimensional points with SSE for int 
  \param aiarrayi_p a int array 
  \param aiarrayi_q a int array
  \param uintidx_length length of the array
 */
inline
double
kernelEuclidean
(const int* aiarrayi_p, const int* aiarrayi_q, const uintidx uintidx_length)
{
  unsigned long int loT_dist = 0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
    
  loT_dist = inrm2_k(lit_n, aiarrayi_p, lit_inc, aiarrayi_q, lit_inc);

  return  std::sqrt((double) loT_dist);
  
}


/*! \fn double kernelEuclideanSquared(const int* aiarrayi_p, const int* aiarrayi_q, uintidx uintidx_length)
  \brief Gets the square Euclidean distance of two int n-dimensional points with SSE
  \details
  \param aiarrayi_p a int array 
  \param aiarrayi_q a int array
  \param uintidx_length length of the array
 */
inline
double
kernelEuclideanSquared
(const int* aiarrayi_p, const int* aiarrayi_q, uintidx uintidx_length)
{
  unsigned long int loT_dist = 0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
    
  loT_dist = inrm2_k(lit_n, aiarrayi_p, lit_inc, aiarrayi_q, lit_inc);

  return  (double) loT_dist;
  
} 

/*! \fn double kernelEuclidean(const short* aiarrays_p, const short* aiarrays_q, const uintidx uintidx_length)
  \brief Gets the Euclidean distance of two short integer  n-dimensional points with SSE 
  \details
  \param aiarrays_p a short int array 
  \param aiarrays_q a short int array
  \param uintidx_length length of the array
 */
inline
double
kernelEuclidean
(const short* aiarrays_p, const short* aiarrays_q, const uintidx uintidx_length)
{
  unsigned int loT_dist = 0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
  
  
  loT_dist = shnrm2_k(lit_n, aiarrays_p, lit_inc, aiarrays_q, lit_inc);
  
  return  std::sqrt((double) loT_dist);
  
}

/*! \fn double  kernelEuclideanSquared(const short* aiarrays_p,const  short* aiarrays_q, const uintidx uintidx_length)
  \brief Gets the square Euclidean distance of two short integer n-dimensional points with SSE
  \details
  \param aiarrays_p a short int array 
  \param aiarrays_q a short int array
  \param uintidx_length length of the array
 */
inline
double 
kernelEuclideanSquared
(const short* aiarrays_p, const  short* aiarrays_q, const uintidx uintidx_length)
{
  unsigned long int loT_dist = 0;
  int64_t lit_n;
  int64_t lit_inc;
  
  lit_n    = (int64_t) uintidx_length;
  lit_inc  = (int64_t) 1;
    
  loT_dist = shnrm2_k(lit_n, aiarrays_p, lit_inc, aiarrays_q, lit_inc);
  
  return  (double) loT_dist;
  
} 

} /*END namespace dist
   */

#endif /*DIST_KERNEL_HPP*/
