/*! \file verbose_global.hpp
 *
 * \brief verbose global
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __VERBOSE_GLOBAL_HPP
#define __VERBOSE_GLOBAL_HPP


#ifdef __VERBOSE_YES

#include <iostream>
#include <iomanip>
#include <sstream>

#include "container_out.hpp"

extern int        geiinparam_verbose;
extern int        geiinparam_verboseMax;
extern const char *geverbosepc_labelstep;
extern uintidx    geverboseui_idproc;

#endif /*__VERBOSE_YES*/


#endif /*__VERBOSE_GLOBAL_HPP*/

