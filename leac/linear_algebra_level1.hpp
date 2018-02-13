/*! \file linear_algebra_level1.hpp
 *
 * \brief linear algebra level1
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __LINEAR_ALGEBRA_LEVEL1_HPP
#define __LINEAR_ALGEBRA_LEVEL1_HPP

#ifdef __WITH_OPEN_BLAS
#include "interface_sse_level1_64bits.hpp"
#else
#include "interface_level1.hpp"
#endif /*__WITH_OPEN_BLAS*/

#endif /*__LINEAR_ALGEBRA_LEVEL1_HPP*/
