/*! \file dist_template.hpp
 *
 * \brief distance template
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef DIST_TEMPLATE_HPP
#define DIST_TEMPLATE_HPP

/*! \namespace dist
  \brief Module for definition of distance between objects or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  dist {


/*! \fn T_FEATURE kernelEuclidean(const T_FEATURE* aiad_p, const T_FEATURE* aiad_q, const uintidx aiuintidx_length)
  \brief Kernel euclidean distance for T_FEATURE
  \details
  \param aiad_p an array of T_FEATURE 
  \param aiad_q an array of T_FEATURE
  \param aiuintidx_length length of the array
 */
template < typename T_FEATURE>
T_FEATURE
kernelEuclidean
(const T_FEATURE* aiad_p, const T_FEATURE* aiad_q, const uintidx aiuintidx_length)
{
  T_FEATURE lt_sumDiff = T_FEATURE(0);
  T_FEATURE lt_diff = T_FEATURE(0);
  
  for ( uintidx lui_i = 0; lui_i < aiuintidx_length; lui_i++ ) {
    lt_diff  = aiad_p[lui_i] - aiad_q[lui_i];
    lt_sumDiff += lt_diff*lt_diff;
  }
  
  return   std::sqrt(lt_sumDiff);  
}


/*! \fn T_FEATURE kernelEuclideanSquared(const T_FEATURE* aiad_p, const T_FEATURE* aiad_q, const uintidx aiuintidx_length)
  \brief Kernel euclidean distance for T_FEATURE
  \details
  \param aiad_p an array of T_FEATURE 
  \param aiad_q an array of T_FEATURE
  \param aiuintidx_length length of the array
 */
template < typename T_FEATURE>
T_FEATURE
kernelEuclideanSquared
(const T_FEATURE* aiad_p, const T_FEATURE* aiad_q, const uintidx aiuintidx_length)
{
  T_FEATURE lt_sumDiff = T_FEATURE(0);
  T_FEATURE lt_diff = T_FEATURE(0);
  
  for ( uintidx lui_i = 0; lui_i < aiuintidx_length; lui_i++ ) {
    lt_diff  = aiad_p[lui_i] - aiad_q[lui_i];
    lt_sumDiff += lt_diff*lt_diff;
  }
  
  return  lt_sumDiff;

  //return  std::inner_product(aiad_p,aiad_p+aiuintidx_length,aiad_q, 0);
   
}

} /*END namespace dist 
   */

#endif /*DIST_TEMPLATE_HPP*/

