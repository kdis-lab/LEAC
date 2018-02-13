/*! \file interface_sse_level1_64bits.hpp
 *
 * \brief interface sse level1 64bits
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INTERFACE_SSE_LEVEL1_HPP
#define INTERFACE_SSE_LEVEL1_HPP

#include <string.h>   /*memcpy*/
#include <stdint.h> //int64_t

#include "common.hpp" /*uintidx*/
#include "common_ssekernel_64bits.hpp"
#include "cblas.h"


/*------------------------------------------------------------------------------
 INTERFACE SSE FUNCTION
------------------------------------------------------------------------------ */

/*! \namespace interfacesse
  \brief Interface to Streaming SIMD Extensions (SSE) high Performance Level Low Functions
  \details Functions base in SSE instructions operate on packed double-precision floating-point values contained in XMM registers and on packed integers contained in MMX and XMM registers \cite progguide:intel10

  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace interfacesse {

/*! \fn void copya(unsigned int *aoarrayuint_x, const unsigned int aiui_alpha, const uintidx aiui_lengthArray)
  \brief copy array aid_alpha in aoarrayuint_x of unsigned int
  \details
\f[
   x \leftarrow \alpha
\f]
  \param aoarrayuint_x a array of unsigned int
  \param aiui_alpha a constant unsigned
  \param aiui_lengthArray a unsigned integer with aoarrayuint_x length 
*/
inline
void
copya
(unsigned int        *aoarrayuint_x,
 const unsigned int  aiui_alpha,
 const uintidx       aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incx = (int64_t) 1;
  
  uicopya_k(lit_n, aiui_alpha, aoarrayuint_x, lit_incx);
}



/*! \fn void copya(int *aoarrayint_x, const int aii_alpha, const uintidx aiui_lengthArray)
  \brief  copy array aii_alpha in aoarrayint_x of int
  \details
\f[
   x \leftarrow \alpha
\f]
  \param aoarrayint_x a array of int
  \param aii_alpha a constant int
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void
copya
(int           *aoarrayint_x,
 const int     aii_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incx = (int64_t) 1;
  
  icopya_k(lit_n, aii_alpha ,aoarrayint_x, lit_incx);
}


/*! \fn void copya(long int *aoarraylong_x, const long int ail_alpha, const uintidx  aiui_lengthArray)
  \brief  copy array aid_alpha in aoarraylong_x of long int
  \details
\f[
   x \leftarrow \alpha
\f]
  \param aoarraylong_x a array of long int
  \param ail_alpha a long int 
  \param aiui_lengthArray a unsigned integer with aoarraylong_x length 
*/
void
copya
(long int       *aoarraylong_x,
 const long int ail_alpha,
 const uintidx  aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incx = (int64_t) 1;
  
  lcopya_k(lit_n, ail_alpha ,aoarraylong_x, lit_incx);
}


/*! \fn void copya(float *aoarrayfloat_x, float aif_alpha, uintidx aiui_lengthArray) 
  \brief  copy array aid_alpha in aoarrayfloat_x of float
  \details
\f[
   x \leftarrow \alpha
\f]
  \param aoarrayfloat_x a array of float
  \param aif_alpha a float
  \param aiui_lengthArray a unsigned integer with aoarrayfloat_x length 
*/
inline
void
copya
(float   *aoarrayfloat_x,
 float   aif_alpha,
 uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incx = (int64_t) 1;
  
  scopya_k(lit_n, aif_alpha ,aoarrayfloat_x, lit_incx);
}


/*! \fn void void copya(double *aoarraydouble_x, double  aid_alpha, uintidx aiui_lengthArray) 
  \brief  copy array aid_alpha in aoarraydouble_x of double
  \details
\f[
   x \leftarrow \alpha
\f]
  \param aoarraydouble_x a array of double
  \param aid_alpha a double
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/
inline
void
copya
(double  *aoarraydouble_x,
 double  aid_alpha,
 uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incx = (int64_t) 1;
  
  dcopya_k(lit_n, aid_alpha ,aoarraydouble_x, lit_incx);
}


/*! \fn void swap(int *aioarrayint_y, int *aioarrayint_x, const uintidx aiui_lengthArray) 
  \brief  swap items of aioarrayint_y to aioarrayint_x
  \details
\f[
   y \leftrightarrow x
\f]
  \param aioarrayint_y a array of int
  \param aioarrayint_x a array of int
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void 
swap
(int           *aioarrayint_y,
 int           *aioarrayint_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_incx;
  int64_t lbi_incy;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_incx = (int64_t) 1;
  lbi_incy = (int64_t) 1;

  iswap_kh(lbi_n, 0, 0, 0, aioarrayint_x, lbi_incx, aioarrayint_y, lbi_incy, NULL,0);
}


/*! \fn void scalInv(int *aoarrayint_x, const long ail_alpha, const uintidx aiui_lengthArray)
  \brief  scal items of the array int
  \details
\f[
   x \leftarrow \frac{1}{\alpha} x
\f]
  \param aoarrayint_x a array of int 
  \param ail_alpha a long scale factor
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void
scalInv
(int           *aoarrayint_x,  
 const long    ail_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  double ld_alpha = (ail_alpha ==0)?0.0d:(1.0d/(double)ail_alpha);
   
  iriscal_kh(lbi_n, 0, 0, ld_alpha, aoarrayint_x, lbi_inc, NULL, 0, NULL, 0);
}

/*! \fn void scal(int *aoarrayint_x, const long ail_alpha, const uintidx aiui_lengthArray)
  \brief  scal items of the array int
  \details
\f[
   x \leftarrow \alpha x
\f]
  \param aoarrayint_x a array of int 
  \param ail_alpha a long scale factor
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void
scal
(int           *aoarrayint_x,  
 const long    ail_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  double ld_alpha((double)ail_alpha);
  
  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;
   
  iriscal_kh(lbi_n, 0, 0, ld_alpha, aoarrayint_x, lbi_inc, NULL, 0, NULL, 0);
}
 

/*! \fn void copy(unsigned int *aoarrayuint_y, const unsigned int *aiarrayuint_x, const uintidx aiui_lengthArray)
  \brief  copy array aiarrayuint_x in aoarrayuint_y of unsigned int
  \details
\f[
   y \leftarrow x
\f]
  \param aoarrayuint_y a array unsigned int
  \param aiarrayuint_x a array unsigned int
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
inline
void
copy
(unsigned int       *aoarrayuint_y,
 const unsigned int *aiarrayuint_x,
 const uintidx      aiui_lengthArray
)
{
  memcpy(aoarrayuint_y,aiarrayuint_x,aiui_lengthArray*sizeof(unsigned int)); //FALTA
}


/*! \fn void copy(int *aoarrayint_y, const int *aiarrayint_x, const uintidx aiui_lengthArray) 
  \brief  copy array aiarrayint_x in aoarrayint_y of int
  \details
\f[
   y \leftarrow x
\f]
  \param aoarrayint_y a array int
  \param aiarrayint_x a array int
  \param aiui_lengthArray a unsigned integer with aiarrayint_x length 
*/  
inline
void
copy
(int          *aoarrayint_y,
 const int    *aiarrayint_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  icopy_k(lbi_n,aiarrayint_x,lbi_inc,aoarrayint_y,lbi_inc);
}



/*! \fn void copy(long int *aoarraylong_y, const long int *aiarraylong_x, const  uintidx aiui_lengthArray)
  \brief  copy array aoarrayt_x in aoarraylong_y of long int
  \details
\f[
   y \leftarrow x
\f]
  \param aoarraylong_y a array long int
  \param aiarraylong_x a array long int
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
inline
void
copy
(long int       *aoarraylong_y,
 const long int *aiarraylong_x,
 const  uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  lcopy_k(lbi_n,aiarraylong_x,lbi_inc,aoarraylong_y,lbi_inc);
}


/*! \fn  long int long int dot (const int *aiarrayint_x, const int *aiarrayint_y, const uintidx aiui_lengthArray)
  \brief  dot product 
  \details
\f[
   dot \leftarrow  x^T y 
\f]
  \param aiarrayint_x a array of int
  \param aiarrayint_y a array of int
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
inline
long int   
dot
(const int     *aiarrayint_x,
 const int     *aiarrayint_y,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t) 1;
  lit_incy  = (int64_t) 1;
   
  return idot_kh(lit_n, aiarrayint_x, lit_incx, aiarrayint_y, lit_incy);
}


/*! \fn void axpyInv (int *aoarraint_y, const long int  ail_alpha, const long int  *aiarraylong_x, const uintidx aiui_lengthArray)
  \brief  axpyInv
  \details
\f[
   y \leftarrow \frac{1}{\alpha} x + y
\f]
  \param aoarraint_y a array of int
  \param ail_alpha a long integer scale factor
  \param aiarraylong_x a array of long int
  \param aiui_lengthArray a unsigned integer with aiarraylong_x length 
*/
inline
void
axpyInv
(int             *aoarraint_y,
 const long int  ail_alpha,
 const long int  *aiarraylong_x,
 const uintidx   aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;

  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;

  double ld_alpha = (ail_alpha ==0)?0.0d:(1.0d/(double)ail_alpha);
   
  ilaxpy_kh(lit_n, 0, 0, ld_alpha, aiarraylong_x, lit_incx, aoarraint_y, lit_incy, NULL, 0);
}
  

/*! \fn void axpy(int *aoarrayint_y, const int aii_alpha, const long int  *aiarraylong_x, const uintidx aiui_lengthArray)
  \brief xAXPY  operation of blas for int with long 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarrayint_y a array of long int
  \param aii_alpha a const double
  \param aiarraylong_x a array of long int
  \param aiui_lengthArray a unsigned integer with aiarraylong_x length 
*/
inline
void
axpy
(int             *aoarrayint_y,
 const int       aii_alpha,
 const long int  *aiarraylong_x,
 const uintidx   aiui_lengthArray
) 
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  double  ld_alpha((double)aii_alpha);
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  ilaxpy_kh(lit_n, 0, 0, ld_alpha, aiarraylong_x, lit_incx, aoarrayint_y, lit_incy, NULL, 0);
}

  
/*! \fn void axpy (int *aoarrayint_y, const long int  ail_alpha, const int *aiarrayint_x, const uintidx aiui_lengthArray)
  \brief xAXPY operation of blas for int with long and  int 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarrayint_y a array of long int
  \param ail_alpha a const long int 
  \param aiarrayint_x a array of long int
  \param aiui_lengthArray a unsigned integer with aiarrayint_x length 
*/
inline
void
axpy
(int             *aoarrayint_y,
 const long int  ail_alpha,
 const int       *aiarrayint_x,
 const uintidx   aiui_lengthArray
) 
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  int     li_alpha((int)ail_alpha);
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  iiaxpy_kh(lit_n, 0, 0, li_alpha, aiarrayint_x, lit_incx, aoarrayint_y, lit_incy, NULL, 0);
}
  
  
/*! \fn void axpy (long int *aoarraylong_y, const int aii_alpha, const int *aiarrayint_x, const uintidx aiui_lengthArray)
  \brief  xAXPY operation of blas for integer
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarraylong_y a array of long int
  \param aii_alpha a const double
  \param aiarrayint_x a array of long int
  \param aiui_lengthArray a unsigned integer with length of array 
*/
inline
void
axpy
(long int      *aoarraylong_y,
 const int     aii_alpha,
 const int     *aiarrayint_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  laxpy_kh(lit_n, 0, 0, aii_alpha, aiarrayint_x, lit_incx, aoarraylong_y, lit_incy, NULL, 0);
}


/*! \fn void aysxpy(double *aioarraydouble_y, const double aid_alpha, const double *aiarraydouble_x, const uintidx aiui_lengthArray)
  \brief  extension of xAXPY operation for double
  \details
\f[
   y \leftarrow y + \alpha (y - x)
\f]
  \param aioarraydouble_y a array of double
  \param aid_alpha  a const double
  \param aiarraydouble_x a array of double
  \param aiui_lengthArray a unsigned integer with aiarraydouble_x length 
*/
inline
void
aysxpy
(double        *aioarraydouble_y,
 const double  aid_alpha,
 const double  *aiarraydouble_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  daysxpy_kh(lit_n, 0, 0, aid_alpha, aiarraydouble_x, lit_incx, aioarraydouble_y, lit_incy, NULL, 0);
  
}


/*! \fn void aysxpy(int *aioarrayint_y, const double  aid_alpha, const int *aiarrayint_x, const uintidx aiui_lengthArray)
  \brief  extension of xAXPY operation for int result
  \details
\f[
   y \leftarrow y + \alpha (y - x)
\f]
  \param aioarrayint_y a array of int
  \param aid_alpha  a const double
  \param aiarrayint_x a array of double
  \param aiui_lengthArray a unsigned integer with aiarrayint_x length 
*/
inline
void
aysxpy
(int           *aioarrayint_y,
 const double  aid_alpha,
 const int     *aiarrayint_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  iaysxpy_kh(lit_n, 0, 0, aid_alpha, aiarrayint_x, lit_incx, aioarrayint_y, lit_incy, NULL, 0);

}


/*! \fn void transy (int *aioarrayint_y, const int ait_alpha, const uintidx aiui_lengthArray)
  \brief  translate aioarrayint_y array for int
  \details
\f[
   y \leftarrow y +  \alpha
\f]
  \param aioarrayint_y a array of int
  \param ait_alpha a int 
  \param aiui_lengthArray a unsigned integer with aioarrayint_y length 
*/
inline
void
transy
(int          *aioarrayint_y,
 const int     ait_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incy  = (int64_t)  1;
   
  itrans_kh(lit_n, 0, 0, ait_alpha, aioarrayint_y, lit_incy, NULL, 0, NULL, 0);
}

/*! \fn void transy(double *aoarraydouble_y, const double  aid_alpha, const uintidx aiui_lengthArray)
  \brief  translate aoarraydouble_y array for double
  \details
\f[
   y \leftarrow y +  \alpha
\f]
  \param aoarraydouble_y a array of int
  \param aid_alpha a double 
  \param aiui_lengthArray a unsigned integer with aoarraydouble_y length 
*/
inline
void
transy
(double        *aoarraydouble_y,
 const double  aid_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incy  = (int64_t)  1;
   
  dtrans_kh(lit_n, 0, 0, aid_alpha, aoarraydouble_y, lit_incy, NULL, 0, NULL, 0);
}


/*! \fn long int sum(const int *aiarrayint_y, const uintidx aiui_lengthArray)
  \brief  return the sum array items int
  \details
\f[
    \sum x_i 
\f]
  \param aiarrayint_y a array of int
  \param aiui_lengthArray a unsigned integer with aiarrayint_y length 
*/
inline
long int
sum
(const int    *aiarrayint_y,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incy;
  long int  lt_out;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incy = (int64_t) 1;
  
  lt_out = isum_kh(lit_n, aiarrayint_y, lit_incy, NULL, 0);

  return lt_out;
}


/*! \fn  float sum(const float *aiarrayfloat_y, const uintidx aiui_lengthArray)
  \brief  return the sum array items float
  \details
\f[
    \sum x_i 
\f]
  \param aiarrayfloat_y a array of float
  \param aiui_lengthArray a unsigned integer with aiarrayfloat_y length 
*/
inline
float
sum
(const float   *aiarrayfloat_y,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incy;
  float lt_out;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incy = (int64_t) 1;
  
  lt_out = ssum_kh(lit_n, aiarrayfloat_y, lit_incy, NULL, 0);

  return lt_out;
}

/*! \fn  double sum(const double *aiarraydouble_y, const uintidx aiui_lengthArray)
  \brief  return the sum array items double
  \details
\f[
    \sum x_i 
\f]
  \param aiarraydouble_y a array of double
  \param aiui_lengthArray a unsigned integer with aiarraydouble_y length 
*/
inline
double
sum
(const double *aiarraydouble_y,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incy;
  double  lt_out;

  lit_n    = (int64_t) aiui_lengthArray;
  lit_incy = (int64_t) 1;
  
  lt_out = dsum_kh(lit_n, aiarraydouble_y, lit_incy, NULL, 0);

  return lt_out;
}


/*! \fn void aasxpa(const float aif_alpha, float *aimatrixrowfloat_a, const uintidx aiui_numRows, const uintidx aiui_numColumns, const float *aiarrayfloat_x)
  \brief  For a matrix each items change for constant and  a vector
  \details
\f[
   a_{ij} = a_{ij} + \alpha(a_{ij} - x{j})

\f]
  \param aif_alpha a float constant
  \param aimatrixrowfloat_a a matrix how a vector row
  \param aiui_numRows a unsigned integer numer rows of matrix 
  \param aiui_numColumns a unsigned integer numer columns  of matrix 
  \param aiarrayfloat_x an array of float
*/
inline 
void
aasxpa
(const float    aif_alpha,
 float          *aimatrixrowfloat_a,
 const uintidx  aiui_numRows,
 const uintidx  aiui_numColumns,
 const float    *aiarrayfloat_x
 )
{
  const int64_t m = (int64_t)  aiui_numColumns;
  const int64_t n = (int64_t) aiui_numRows;
  const int64_t lda = m;

  saasxpa_kh(m, n, 0, aif_alpha, aimatrixrowfloat_a, lda, aiarrayfloat_x, 0, NULL, 0);
}

/*! \fn void aasxpa(const double aid_alpha, double *aimatrixrowdouble_a, const uintidx aiui_numRows, const uintidx aiui_numColumns, const double *aoarrayt_x)
  \brief  For a matrix each items change for constant and  a vector
  \details
\f[
   a_{ij} = a_{ij} + \alpha(a_{ij} - x{j})

\f]
  \param aid_alpha a double constant
  \param aimatrixrowdouble_a a matrix how a vector row
  \param aiui_numRows a unsigned integer numer rows of matrix 
  \param aiui_numColumns a unsigned integer numer columns  of matrix 
  \param aoarrayt_x a vector 
*/
inline
void
aasxpa
(const double    aid_alpha,
 double          *aimatrixrowdouble_a,
 const uintidx   aiui_numRows,
 const uintidx   aiui_numColumns,
 const double    *aoarrayt_x
 )
{
  const int64_t m = (int64_t)  aiui_numColumns;
  const int64_t n = (int64_t) aiui_numRows;
  const int64_t lda = m;

  daasxpa_kh(m, n, 0, aid_alpha, aimatrixrowdouble_a, lda, aoarrayt_x, 0, NULL, 0);

}


/*------------------------------------------------------------------------------
  Openblas
------------------------------------------------------------------------------ */

/*! \fn void  swap(float *aioarrayfloat_y, float *aioarrayfloat_x, const uintidx aiui_lengthArray) 
  \brief  swap items of aioarrayfloat_y to aioarrayfloat_x for float
  \details
\f[

 y \leftrightarrow x
\f]
  \param aioarrayfloat_y a array of float
  \param aioarrayfloat_x a array of float
  \param aiui_lengthArray a unsigned integer with aioarrayfloat_x length 
*/
inline
void 
swap
(float *aioarrayfloat_y,
 float *aioarrayfloat_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_incx;
  int64_t lbi_incy;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_incx = (int64_t) 1;
  lbi_incy = (int64_t) 1;

  cblas_sswap(lbi_n, aioarrayfloat_x, lbi_incx, aioarrayfloat_y, lbi_incy);
  
}


/*! \fn void void swap(double *aioarraydouble_y, double *aioarraydouble_x, const uintidx aiui_lengthArray)
  \brief  swap items of aioarraydouble_y to aioarraydouble_x for double 
  \details
\f[
 y \leftrightarrow x
\f]
  \param aioarraydouble_y a array of double
  \param aioarraydouble_x a array of double
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/
inline
void 
swap
(double *aioarraydouble_y,
 double *aioarraydouble_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_incx;
  int64_t lbi_incy;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_incx = (int64_t) 1;
  lbi_incy = (int64_t) 1;

  cblas_dswap(lbi_n, aioarraydouble_x, lbi_incx, aioarraydouble_y, lbi_incy);
}


/*! \fn void scal(float *aoarrayfloat_x, const float aif_alpha, const uintidx aiui_lengthArray) 
  \brief scal items of the array float
  \details
\f[
   x \leftarrow \alpha x
\f]
  \param aoarrayfloat_x a array of float 
  \param aif_alpha a float scale factor constant
  \param aiui_lengthArray a unsigned integer with aoarrayfloat_x length 
*/  
inline
void
scal
(float         *aoarrayfloat_x,  
 const float   aif_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;
   
  cblas_sscal(lbi_n, aif_alpha, aoarrayfloat_x, lbi_inc);
}

/*! \fn void scalInv(int *aoarrayint_x, const long ail_alpha, const uintidx aiui_lengthArray)
  \brief  scal items of the array int
  \details
\f[
   x \leftarrow \frac{1}{\alpha} x
\f]
  \param aoarrayint_x a array of int 
  \param ail_alpha a long scale factor
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void
scalInv
(float            *aoarrayfloat_x,  
 const long int   ail_alpha,
 const uintidx    aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  float lf_alpha = (ail_alpha ==0)?0.0f:(1.0f/(float)ail_alpha);
   
  cblas_sscal(lbi_n, lf_alpha, aoarrayfloat_x, lbi_inc);
}


/*! \fn void scal(double *aoarraydouble_x, const double aid_alpha, const uintidx aiui_lengthArray) 
  \brief scal items of the array float
  \details
\f[
   x \leftarrow \alpha x
\f]
  \param aoarraydouble_x a array of float 
  \param aid_alpha a float scale factor constant
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/  
inline
void
scal
(double        *aoarraydouble_x,  
 const double  aid_alpha,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;
   
  cblas_dscal(lbi_n, aid_alpha, aoarraydouble_x, lbi_inc);
}

/*! \fn void scal(double *aoarraydouble_x, const long int ail_alpha,const uintidx  aiui_lengthArray)
  \brief  scal items of the array double
  \details
\f[
   x \leftarrow \frac{1}{\alpha} x
\f]
  \param aoarraydouble_x a array of double 
  \param ail_alpha a long scale factor
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void
scal
(double         *aoarraydouble_x,  
 const long int ail_alpha,
 const uintidx  aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  double ld_alpha = (double) ail_alpha;
  
  cblas_dscal(lbi_n, ld_alpha, aoarraydouble_x, lbi_inc);
}

/*! \fn void scalInv(int *aoarraydouble_x, const long ail_alpha, const uintidx aiui_lengthArray)
  \brief  scal items of the array double
  \details
\f[
   x \leftarrow \frac{1}{\alpha} x
\f]
  \param aoarraydouble_x a array of double 
  \param ail_alpha a long scale factor
  \param aiui_lengthArray a unsigned integer with aoarrayint_x length 
*/
inline
void
scalInv
(double         *aoarraydouble_x,  
 const long int ail_alpha,
 const uintidx  aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  double ld_alpha = (ail_alpha ==0)?0.0d:(1.0d/(double)ail_alpha);
  
  cblas_dscal(lbi_n, ld_alpha, aoarraydouble_x, lbi_inc);
}


/*! \fn void copy(float *aoarrayfloat_y, const float *aiarrayfloat_x, const uintidx aiui_lengthArray) 
  \brief  copy array aoarrayfloat_x in aoarrayfloat_y of float
  \details
\f[
   y \leftarrow x
\f]
  \param aoarrayfloat_y a array of float
  \param aiarrayfloat_x a array of float
  \param aiui_lengthArray a unsigned integer with aiarrayfloat_x length 
*/
inline
void
copy
(float         *aoarrayfloat_y,
 const float   *aiarrayfloat_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  cblas_scopy(lbi_n, aiarrayfloat_x, lbi_inc, aoarrayfloat_y, lbi_inc);
}

/*! \fn  void copy (double *aoarraydouble_y, const double  *aiarraydouble_x, const uintidx aiui_lengthArray)
  \brief  copy array aoarraydouble_x in aoarraydouble_y of double
  \details
\f[
   y \leftarrow x
\f]
  \param aoarraydouble_y a array of double
  \param aiarraydouble_x a array of double
  \param aiui_lengthArray a unsigned integer with aiarraydouble_x length 
*/  
inline
void
copy
(double        *aoarraydouble_y,
 const double  *aiarraydouble_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lbi_n;
  int64_t lbi_inc;

  lbi_n    = (int64_t) aiui_lengthArray;
  lbi_inc  = (int64_t) 1;

  cblas_dcopy(lbi_n,aiarraydouble_x,lbi_inc,aoarraydouble_y,lbi_inc);
}


/*! \fn float dot(const float *aiarrayfloat_x, const float *aiarrayfloat_y, const uintidx aiui_lengthArray)  
  \brief  dot product for float vectors
  \details
\f[
   dot \leftarrow  x^T y 
\f]
  \param aiarrayfloat_x a array of float
  \param aiarrayfloat_y a array of float
  \param aiui_lengthArray a unsigned integer with aoarrayfloat_x length 
*/  
inline
float   
dot
(const float   *aiarrayfloat_x,
 const float   *aiarrayfloat_y,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  return cblas_sdot(lit_n, aiarrayfloat_x, lit_incx, aiarrayfloat_y, lit_incy);
}


/*! \fn double dot(const double *aiarraydouble_x, const double  *aiarraydouble_y, const uintidx aiui_lengthArray)
  \brief  dot product for double vectors
  \details
\f[
   dot \leftarrow  x^T y 
\f]
  \param aiarraydouble_x a array of float
  \param aiarraydouble_y a array of float
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/  
inline
double   
dot
(const double  *aiarraydouble_x,
 const double  *aiarraydouble_y,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  return cblas_ddot(lit_n, aiarraydouble_x, lit_incx, aiarraydouble_y, lit_incy);
}



/*! \fn void axpy(float *aoarrayfloat_y, const float ait_alpha, const float *aiarrayfloat_x, const uintidx aiui_lengthArray)
  \brief  function xAXPY of blas for float vectors
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarrayfloat_y a array of float
  \param ait_alpha a scale factor constant
  \param aiarrayfloat_x a array of float
  \param aiui_lengthArray a unsigned integer with aoarrayfloat_x length 
*/
inline
void
axpy
(float         *aoarrayfloat_y,
 const float   ait_alpha,
 const float   *aiarrayfloat_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  cblas_saxpy(lit_n, ait_alpha, aiarrayfloat_x, lit_incx, aoarrayfloat_y, lit_incy);
  
}

/*! \fn void axpy(float *aoarrayfloat_y, const long ail_alpha, const float *aiarrayfloat_x, const uintidx aiui_lengthArray) 
  \brief  function xAXPY of blas for floats vectors 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarrayfloat_y a array of float
  \param ail_alpha a scale factor constant of long int
  \param aiarrayfloat_x a array of float
  \param aiui_lengthArray a unsigned integer with aoarrayfloat_x length 
*/
inline
void
axpy
(float         *aoarrayfloat_y,
 const long    ail_alpha,
 const float   *aiarrayfloat_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;

  float lf_alpha((float)ail_alpha);
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  cblas_saxpy(lit_n, lf_alpha, aiarrayfloat_x, lit_incx, aoarrayfloat_y, lit_incy);
  
}

/*! \fn void axpyInv(float *aoarrayfloat_y, const long ail_alpha, const float *aiarrayfloat_x, const uintidx aiui_lengthArray) 
  \brief  function xAXPY-inverse of blas for floats vectors 
  \details
\f[
   y \leftarrow \frac{1}{\alpha} x + y
\f]
  \param aoarrayfloat_y a array of float
  \param ail_alpha a scale factor constant of long int
  \param aiarrayfloat_x a array of float
  \param aiui_lengthArray a unsigned integer with aoarrayfloat_x length 
*/
inline
void
axpyInv
(float         *aoarrayfloat_y,
 const long    ail_alpha,
 const float   *aiarrayfloat_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  float  lf_alpha = (ail_alpha == 0)?0.0f: 1.0f/(float) ail_alpha;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
   
  cblas_saxpy(lit_n, lf_alpha, aiarrayfloat_x, lit_incx, aoarrayfloat_y, lit_incy);
  
}


/*! \fn void axpy(double *aoarraydouble_y, const double aid_alpha, const double *aiarraydouble_x, const uintidx aiui_lengthArray) 
  \brief  function xAXPY of blas for doubles vectors 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarraydouble_y a array of double
  \param aid_alpha a scale factor constant
  \param aiarraydouble_x a array of double
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/
inline
void
axpy
(double       *aoarraydouble_y,
 const double aid_alpha,
 const double *aiarraydouble_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
  
  cblas_daxpy(lit_n, aid_alpha, aiarraydouble_x, lit_incx, aoarraydouble_y, lit_incy);
  
}

/*! \fn void axpy(double *aoarraydouble_y, const long int ail_alpha, const double *aiarraydouble_x, const uintidx aiui_lengthArray) 
  \brief  function xAXPY of blas for doubles vectors 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarraydouble_y a array of double
  \param ail_alpha a scale factor constant long int
  \param aiarraydouble_x a array of double
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/
inline
void
axpy
(double         *aoarraydouble_y,
 const long int ail_alpha,
 const double   *aiarraydouble_x,
 const uintidx  aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;

  double ld_alpha((double)ail_alpha);
  
  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
  
  cblas_daxpy(lit_n, ld_alpha, aiarraydouble_x, lit_incx, aoarraydouble_y, lit_incy);
  
}


/*! \fn void axpyInv(double *aoarraydouble_y, const long ail_alpha, const double *aiarraydouble_x, const uintidx aiui_lengthArray) 
  \brief  function xAXPY-inverse of blas for doubles vectors 
  \details
\f[
   y \leftarrow \frac{1}{\alpha} x + y
\f]
  \param aoarraydouble_y a array of double
  \param ail_alpha a scale factor constant of long int
  \param aiarraydouble_x a array of double
  \param aiui_lengthArray a unsigned integer with aoarraydouble_x length 
*/
inline
void
axpyInv
(double       *aoarraydouble_y,
 const long   ail_alpha,
 const double *aiarraydouble_x,
 const uintidx aiui_lengthArray
)
{
  int64_t lit_n;
  int64_t lit_incx;
  int64_t lit_incy;
  double  ld_alpha = (ail_alpha == 0)?0.0d: 1.0d/(double) ail_alpha;

  lit_n     = (int64_t) aiui_lengthArray;
  lit_incx  = (int64_t)  1;
  lit_incy  = (int64_t)  1;
  

  cblas_daxpy(lit_n, ld_alpha, aiarraydouble_x, lit_incx, aoarraydouble_y, lit_incy);
  
}

} /*END namespace interfacesse*/


#endif  /* INTERFACE_SSE_LEVEL1_HPP */

