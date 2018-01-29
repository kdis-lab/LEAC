/*! \file matrix_operation.hpp
 *
 * \brief matrix operation
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __MATRIX_OPERATION_HPP
#define __MATRIX_OPERATION_HPP

#ifdef __WITH_OPEN_BLAS
#include "matop_openblaslapacke.hpp"
#else
#include "matop_withoutopenblas.hpp"
#endif /*__WITH_OPEN_BLAS*/

#endif /*__MATRIX_OPERATION_HPP*/
