/*! \file common_interfacelapack.hpp
 *
 * \brief common interfacelapack
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef COMMON_INTERFACELAPACK_HPP
#define COMMON_INTERFACELAPACK_HPP

/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat
{

typedef enum TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} TRANSPOSE;

}

#endif /* COMMON_INTERFACELAPACK_HPP */
