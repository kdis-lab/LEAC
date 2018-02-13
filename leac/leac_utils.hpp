/*! \file leac_utils.hpp
 *
 * \brief LEAC utils
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __LEAC_UTILS_HPP
#define __LEAC_UTILS_HPP

#include <utility>
#include <typeinfo>    // operator typeid
#include <type_traits> //std::common_type
#include <cmath>       //std::isless
#include <cassert>

/*! \namespace utils
  \brief General Purpose Utilities
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace utils {

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

/*! \struct InstanceDataType
  \brief Structure For a data type of an attribute of an instance, it defines the data type recommended for the sum of the instances
*/  
struct InstanceDataType
{
  long   sum(int)    const { return long(1); }
  double sum(double) const { return double(1); }
  float  sum(float) const { return float(1); }

  long   sumFrequency(int) const { return long(1); }
  
};

/*! \struct RunOnce
  \brief Structure to ensure that a code or function is executed only once 
*/  
struct RunOnce {
  template <typename T>
  RunOnce(T &&f) { f(); }
};

/*! \fn bool closeenough(A const & a, B const & b, typename std::common_type< A, B >::type const & epsilon)
    \brief Check if two real numbers are close enough or equal 
    \details 
    \param a  a real number
    \param b  a real number
    \param epsilon a real number
 */
template< typename A, typename B >
inline
bool closeenough(A const & a, B const & b,
		  typename std::common_type< A, B >::type const & epsilon)
{
  //using std::isless;
  assert(std::isless(0, epsilon)); // epsilon is a part of the whole quantity
  assert(std::isless(epsilon, 1));
  //using std::abs;
  auto const delta = std::abs(a - b);
  auto const x = std::abs(a);
  auto const y = std::abs(b);
  
  // comparable generally and |a - b| < eps * (|a| + |b|) / 2
  return std::isless(epsilon * y, x) &&
    std::isless(epsilon * x, y) &&
    std::isless((delta + delta) / (x + y), epsilon);
}
  

} /*END namespace utils 
   */


#endif /*__LEAC_UTILS_HPP*/
