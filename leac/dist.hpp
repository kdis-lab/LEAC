/*! \file dist.hpp
 *
 * \brief operations for the calculation of distances between objects
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef DIST_HPP
#define DIST_HPP

#include "common.hpp"

/*! \namespace dist
  \brief Module for definition of distance between objects or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  dist {

/*! \struct Dist
  \brief Generic dist
*/
template < class T_DIST,
	   class T_FEATURE
	   >
struct Dist {
  virtual ~Dist() {}
  virtual T_DIST operator() (const T_FEATURE*, const T_FEATURE*, const uintidx) const = 0;
}; /* Dist */


} /*END namespace dist 
   */

#endif /*DIST_HPP*/


