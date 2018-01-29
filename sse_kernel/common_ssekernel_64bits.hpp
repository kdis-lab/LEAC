/*! \file common_ssekernel_64bits.hpp
 * This file is part of the LEAC.
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */

#ifndef __COMMON_SSEKERNEL_HPP
#define __COMMON_SSEKERNEL_HPP

#include <stdint.h> //int64_t


#ifdef __cplusplus
extern "C" {
#endif

  /* x = a
   */
  void   uicopya_k(int64_t n, unsigned int alpha, unsigned int *x, int64_t incx);

  void   icopya_k(int64_t n, int alpha, int *x, int64_t incx);

  void   lcopya_k(int64_t n, long int alpha, long int *x, int64_t incx);

  void   scopya_k(int64_t n, float alpha, float  *x, int64_t incx);

  void   dcopya_k(int64_t n, double alpha, double  *x, int64_t incx);
  
  /* x <--> y
   */
  int    iswap_kh(const int64_t n, int64_t, int64_t, int, int *x, const int64_t incx, int *y, const int64_t incy, int*, int64_t);

  /*x <-- a x
   */
  int   iriscal_kh(const int64_t n, int64_t, int64_t, const double alpha, int  *x, const int64_t incx, int*, int64_t, int*, int64_t);

  /*y <-- x
   */
  void  icopy_k(const int64_t n, const int *x, const int64_t incx, int *y, const int64_t incy);

  void  lcopy_k(const int64_t n, const long int *x, const int64_t incx, long int *y, const int64_t incy);

  /* y = y + a
   */
  int   dtrans_kh(const int64_t n, int64_t, int64_t, const double alpha, double *y, const int64_t incy, double *, const int64_t , double*, int64_t);

  int  itrans_kh(const int64_t n, int64_t, int64_t, const int alpha, int *y, const int64_t incy, int *, const int64_t , int*, int64_t);
  
  /*y <-- ax + y
   */
  int   ilaxpy_kh(const int64_t n, int64_t, int64_t, const double alpha, const long int *x, const int64_t incx, int *y, const int64_t incy, long int*, int64_t);

  int   laxpy_kh(const int64_t n, int64_t, int64_t, const int alpha, const int *x, const int64_t incx, long int *y, const int64_t incy, long int*, int64_t);

  int   iiaxpy_kh(const int64_t n, int64_t, int64_t, const int alpha, const int *x, const int64_t incx, int *y, const int64_t incy,  int*, int64_t);
  
  /*y = y + a(y-x)
   */
  int   daysxpy_kh(const int64_t n, int64_t, int64_t, const double alpha, const double *x, const int64_t incx, double *y, const int64_t incy, double*, int64_t);

  int   iaysxpy_kh(const int64_t n, int64_t, int64_t, const double alpha, const int *x, const int64_t incx, int *y, const int64_t incy, int*, int64_t);

  /*dot <-- x^Ty
   */
  int         shdot_kh(const int64_t n, const short *x, const int64_t incx, const short *y, const int64_t incy);

  long int    idot_kh(const int64_t n, const int *x, const int64_t incx, const int *y, const int64_t incy);

  /*nom2 <-- ||x||2 
   */
  float ssnrm2_k(const int64_t n, const float *x, const  int64_t incx, const float *y, const int64_t incy);

  double ddnrm2_k(const int64_t n, const double  *x, const  int64_t incx, const double  *y, const int64_t incy);

  unsigned int shnrm2_k(const int64_t n, const short *x, const int64_t incx, const short *y, const int64_t incy);

  unsigned long int inrm2_k(const int64_t n, const int *x, const int64_t incx, const int *y, const int64_t incy);

  /* sum <-- x
   */
  float ssum_kh(const float   *aiarrayfloat_y,const int64_t  aiint64_t_n);
  
  double dsum_kh(const int64_t n, const double  *x, const int64_t incx, const int*, const int64_t);

  long int isum_kh(const int64_t n, const int *x, const int64_t incx, const int*, const int64_t);

  /*a_ij = a_ij + alpha(a_ij - xj)
   */
  int daasxpa_kh(const int64_t m, const int64_t n, const int64_t, const double alpha, const double *a, const int64_t lda, const double *x, int64_t, double*, int64_t);

  void saasxpa_kh(const float aif_alpha, float *aimatrixrowfloat_a, const int64_t  aiint64_numRows, const int64_t  aiint64_numColumns, const float    *aiarrayfloat_x);
  
#ifdef __cplusplus
}
#endif /*__cplusplus*/


#endif /*__COMMON_SSEKERNEL_HPP*/
