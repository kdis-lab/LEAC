/*! \file common.hpp
 *
 * \brief common
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef COMMON_HPP
#define COMMON_HPP

#include <stdint.h>
#include <limits>

#undef min
#undef max

/*! Data type for the addressing variables*/
typedef uint32_t uintidx;

typedef int32_t  intidx;


static const uintidx UINTIDX_NIL = -1;

#define COMMON_IDOMAIN          unsigned int 

#define COMMON_NUM_SEED  8

#define COMMON_COUT_PRECISION          std::numeric_limits<double>::digits10+1     
#define COMMON_VERBOSE_COUT_PRECISION  6

#endif /*COMMON_HPP*/
