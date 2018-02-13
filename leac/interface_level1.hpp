/*! \file interface_level1.hpp
 *
 * \brief Interface_level 1
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INTERFACE_LEVEL1_HPP
#define INTERFACE_LEVEL1_HPP

#include <algorithm>    // std::copy
#include <numeric>
#include <limits>
#include <cmath>

/*! \namespace interfacesse
  \brief Interface to Streaming SIMD Extensions (SSE) high Performance Level Low Functions
  \details Functions base in SSE instructions operate on packed double-precision floating-point values contained in XMM registers and on packed integers contained in MMX and XMM registers \cite progguide:intel10

  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace interfacesse {



/*! \fn void copya(T *aoarrayt_x, const T ait_alpha, const uintidx aiui_lengthArray) 
  \brief Initializes the array with a constant
\f[
   x \leftarrow \alpha
\f]
  \param aoarrayt_x a array of type T
  \param ait_alpha a constant of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
inline
void
copya
(T             *aoarrayt_x,
 const T       ait_alpha,
 const uintidx aiui_lengthArray
)
{
  std::fill_n(aoarrayt_x, aiui_lengthArray, ait_alpha);
}
  
  
/*! \fn void void copy(T *aoarrayt_y, const T *aiarrayt_x, const uintidx aiui_lengthArray)
  \brief  copy array aoarrayt_x in aoarrayt_y of type T
  \details
\f[
   y \leftarrow x
\f]
  \param aoarrayt_y a array of type T
  \param aiarrayt_x a array of type T 
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
inline
void
copy
(T             *aoarrayt_y,
 const T       *aiarrayt_x,
 const uintidx aiui_lengthArray
)
{
  std::copy(aiarrayt_x, aiarrayt_x+aiui_lengthArray, aoarrayt_y);
}
  
/*! \fn  void swap(T *aoarrayt_y, T *aiarrayt_x, const uintidx aiui_lengthArray)
  \brief  swap items of aoarrayt_y to aiarrayt_x
  \details
\f[
   y \leftrightarrow x
\f]
  \param aoarrayt_y a array of type T
  \param aiarrayt_x a array of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
inline
void 
swap
(T   *aoarrayt_y,
 T   *aiarrayt_x,
 const uintidx aiui_lengthArray
)
{
  std::swap_ranges(aoarrayt_y, aoarrayt_y+aiui_lengthArray,aiarrayt_x);
}


/*! \fn  T dot(T *aoarrayt_y,T *aiarrayt_x, const uintidx aiui_lengthArray)
  \brief  swap items of aoarrayt_y to aiarrayt_x
  \details
\f[
   y \leftrightarrow x
\f]
  \param aoarrayt_y a array of type T
  \param aiarrayt_x a array of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
inline
T
dot
(T   *aoarrayt_y,
 T   *aiarrayt_x,
 const uintidx aiui_lengthArray
)
{
  return  std::inner_product(aoarrayt_y,aoarrayt_y+aiui_lengthArray, aiarrayt_x, T(0));
}


/*! \fn T sum(const T *aoarrayt_y, const uintidx aiui_lengthArray)
  \details
\f[
    \sum x_i 
\f]
  \param aoarrayt_y a array of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
inline
T
sum
(const T       *aoarrayt_y,
 const uintidx aiui_lengthArray
)
{
  return std::accumulate(aoarrayt_y,aoarrayt_y+aiui_lengthArray, T(0));
}

/*! \fn void axpy(T *aoarrayt_y, const T ait_alpha, const T *aiarrayt_x, const uintidx aiui_lengthArray)
  \brief  function xAXPY of blas for doubles vectors 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarrayt_y a array of type T
  \param ait_alpha a scale factor constant
  \param aiarrayt_x a array of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T_INT1, typename T_INT2, typename T_INT3 >
void
axpy
(T_INT1        *aoarrayt_y,  //INT
 const T_INT2  ait_alpha,    //LONG
 const T_INT3  *aiarrayt_x,  //INT
 const uintidx aiui_lengthArray
)
{
  const T_INT1 *last1  = aoarrayt_y + aiui_lengthArray;

  while (aoarrayt_y != last1) {
    *aoarrayt_y = (T_INT1) (ait_alpha * (*aiarrayt_x) + (*aoarrayt_y));
    ++aoarrayt_y;
    ++aiarrayt_x;
  }
}
  
/*! \fn void axpy(long int* aoarraylongint_y, const int ait_alpha, const int* aiarrayint_x, const uintidx aiui_lengthArray)
  \brief  function xAXPY of blas for doubles vectors 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarraylongint_y a array of long int
  \param ait_alpha a scale factor constant of type int
  \param aiarrayint_x a array of int 
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
/*void
axpy
(long int*      aoarraylongint_y,
 const int      ait_alpha,
 const int*     aiarrayint_x,
 const uintidx  aiui_lengthArray
)
{
  long int *last1  = aoarraylongint_y + aiui_lengthArray;

  while (aoarraylongint_y != last1) {
    *aoarraylongint_y += ait_alpha * (*aiarrayint_x); 
    ++aoarraylongint_y;
    ++aiarrayint_x;
  }
  
}
*/

/*! \fn void axpy (int* aoarrayint_y, const double aid_alpha, const long int* aiarraylongint_x, const uintidx aiui_lengthArray)
  \brief  function xAXPY of blas for doubles vectors 
  \details
\f[
   y \leftarrow \alpha x + y
\f]
  \param aoarrayint_y a array of int
  \param aid_alpha a constant with a double-scale factor
  \param aiarraylongint_x a array of long int
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
/*void
axpy
(int*            aoarrayint_y,
 const double    aid_alpha,
 const long int* aiarraylongint_x,
 const uintidx   aiui_lengthArray
)
{ 
  int *last1  = aoarrayint_y + aiui_lengthArray;

  while (aoarrayint_y != last1) {
    *aoarrayint_y +=  int(std::round(aid_alpha * (*aiarraylongint_x))); 
    ++aoarrayint_y;
    ++aiarraylongint_x;
  } 
}
*/

/*! \fn void axpyInv (T_INT1 *aoarrayint1_y, const T_INT2  aiint2_alpha, const T_INT3  *aiarrayint3_x, const uintidx aiui_lengthArray)
  \brief  function xAXPY of blas for doubles vectors 
  \details
\f[
   y \leftarrow \frac{1}{\alpha} x + y
\f]
  \param aoarrayint1_y a array of integer 
  \param aiint2_alpha a constant with an integer scaling factor
  \param aiarrayint3_x a array of integer
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T_INT1, typename T_INT3, typename T_INT2 >
void
axpyInv
(T_INT1        *aoarrayint1_y,
 const T_INT2  aiint2_alpha,
 const T_INT3  *aiarrayint3_x,
 const uintidx aiui_lengthArray
)
{
  T_INT1 *last1  = aoarrayint1_y + aiui_lengthArray;

  double ld_alpha = (aiint2_alpha ==0)?0.0:
    ((1.0+10.0*std::numeric_limits<double>::epsilon())/(double)aiint2_alpha);
  
  while (aoarrayint1_y != last1) {
    *aoarrayint1_y += (T_INT1) (ld_alpha * (double) (*aiarrayint3_x));
    ++aoarrayint1_y;
    ++aiarrayint3_x;
  }
}


/*! \fn  void scal(T *aoarrayt_x, const T ait_alpha, const uintidx aiui_lengthArray)
  \brief scal items of the array of type T
  \details
\f[
   x \leftarrow \alpha x
\f]
  \param aoarrayt_x a array of type  T
  \param ait_alpha a scale factor constant of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
void
scal
(T            *aoarrayt_x,  
 const T       ait_alpha,
 const uintidx aiui_lengthArray
)
{
  T *last1  = aoarrayt_x + aiui_lengthArray;
  while (aoarrayt_x != last1) {
    *aoarrayt_x *= ait_alpha;
    ++aoarrayt_x;
  } 
}

/*! \fn  void scal(T_INT1 *aoarrayit_x, const T_INT2 ait_alpha, const uintidx aiui_lengthArray)
  \brief scal items of the array of type integer
  \details
\f[
   x \leftarrow \alpha x
\f]
  \param aoarrayit_x a array of type integer
  \param ait_alpha a integer scale factor constant
  \param aiui_lengthArray a unsigned integer with aoarrayit_x length 
*/
template < typename T_INT1, typename T_INT2>
void
scal
(T_INT1        *aoarrayint1_x,  
 const T_INT2  ait_alpha,
 const uintidx aiui_lengthArray
)
{
  T_INT1 *last1  = aoarrayint1_x + aiui_lengthArray;
  while (aoarrayint1_x != last1) {
    *aoarrayint1_x = (T_INT1) ((T_INT2)(*aoarrayint1_x) * (ait_alpha));
    ++aoarrayint1_x;
  }
  
}


/*! \fn  void scalInv(T_INT *aoarrayint1_x, const T_INT2 ait_alpha, const uintidx  aiui_lengthArray)
  \brief scal items of the array of type integer
  \details
\f[
   x \leftarrow \frac{1}{\alpha} x
\f]
  \param aoarrayint1_x a array of type integer
  \param ait_alpha a constant with an integer scaling factor
  \param aiui_lengthArray a unsigned integer with aoarrayint1_x length 
*/
template < typename T_INT1, typename T_INT2 >
void
scalInv
(T_INT1         *aoarrayint1_x,  
 const T_INT2   ait_alpha,
 const uintidx  aiui_lengthArray
)
{
  double ld_alpha = (ait_alpha ==0)?0.0:
    ((1.0+10.0*std::numeric_limits<double>::epsilon())/(double)ait_alpha);

  T_INT1 *last1  = aoarrayint1_x + aiui_lengthArray;
  while (aoarrayint1_x != last1) {
    *aoarrayint1_x = (T_INT1) (ld_alpha * (double) (*aoarrayint1_x));
    ++aoarrayint1_x;
  }
}

/*! \fn void transy(T *aoarrayt_y, const T ait_alpha, const uintidx aiui_lengthArray)
  \brief  translate the  aoarrayt_y array
  \details
\f[
   y \leftarrow y +  \alpha
\f]
  \param aoarrayt_y a array of T
  \param ait_alpha a T
  \param aiui_lengthArray a unsigned integer with aoarrayt_y length 
*/
template < typename T>
void
transy
(T             *aoarrayt_y,
 const T       ait_alpha,
 const uintidx aiui_lengthArray
)
{

  T *last1  = aoarrayt_y + aiui_lengthArray;
  while (aoarrayt_y != last1) {
    *aoarrayt_y += ait_alpha;
    ++aoarrayt_y;
  }
  
}

/*! \fn void aasxpa(const T ait_alpha, T *aimatrixrowt_a, const uintidx aiui_numRows, const uintidx aiui_numColumns, const T *aoarrayt_x)
  \brief  For a matrix each items change for constant and  a vector
  \details
\f[
   a_{ij} = a_{ij} + \alpha(a_{ij} - x{j})

\f]
  \param ait_alpha a constant with a T scaling factor
  \param aimatrixrowt_a a matrix how a vector row
  \param aiui_numRows a unsigned integer numer rows of matrix 
  \param aiui_numColumns a unsigned integer numer columns  of matrix 
  \param aoarrayt_x a vector
*/
template < typename T>
void
aasxpa
(const T       ait_alpha,
 T             *aimatrixrowt_a,
 const uintidx aiui_numRows,
 const uintidx aiui_numColumns,
 const T       *aoarrayt_x
 )
{
  const T *last1 = aoarrayt_x + aiui_numColumns;

  for ( uintidx liu_i = 0; liu_i < aiui_numRows; ++liu_i) {
    const T *larrayt_xIter = aoarrayt_x;
    while (larrayt_xIter != last1) {
      *aimatrixrowt_a += ait_alpha * ((*aimatrixrowt_a) - (*larrayt_xIter));
      ++aimatrixrowt_a;
      ++larrayt_xIter;
    }
    
  }
}
  
/*! \fn void aysxpy(T *aoarrayt_y, const T ait_alpha, const T *aiarrayt_x, const uintidx aiui_lengthArray)
  \brief  extension of xAXPY operation for T result
  \details
\f[
   y \leftarrow y + \alpha (y - x)
\f]
  \param aoarrayt_y a array of type T
  \param ait_alpha  a const of type T
  \param aiarrayt_x a array of type T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
template < typename T>
void
aysxpy
(T             *aoarrayt_y,
 const T       ait_alpha,
 const T       *aiarrayt_x,
 const uintidx aiui_lengthArray
)
{
  T *last1  = aoarrayt_y + aiui_lengthArray;
  T lit_alpha = T(ait_alpha);
  while (aoarrayt_y != last1) {
    *aoarrayt_y += lit_alpha * ((*aoarrayt_y) - (*aiarrayt_x));
    ++aoarrayt_y;
    ++aiarrayt_x;
  }
}

/*! \fn void aysxpy (int *aoarrayi_y, const double aid_alpha, const int *aiarrayi_x, const uintidx  aiui_lengthArray)
  \brief  extension of xAXPY operation 
  \details
\f[
   y \leftarrow y + \alpha (y - x)
\f]
  \param aoarrayt_y a array of T
  \param T  a const T
  \param aoarrayt_x a array of T
  \param aiui_lengthArray a unsigned integer with aoarrayt_x length 
*/
void
aysxpy
(int            *aoarrayi_y,
 const double   aid_alpha,
 const int      *aiarrayi_x,
 const uintidx  aiui_lengthArray
)
{
  int  *last1  = aoarrayi_y + aiui_lengthArray;
  
  while (aoarrayi_y != last1) {
    *aoarrayi_y += int(std::round(aid_alpha * ((*aoarrayi_y) - (*aiarrayi_x))));
    ++aoarrayi_y;
    ++aiarrayi_x;
  }
}
  
} /*END namespace interfacesse*/


#endif  /* INTERFACE_LEVEL1_HPP */
