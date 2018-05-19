/*! \file verbose.hpp
 *
 * \brief verbose
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __VERBOSE_HPP
#define __VERBOSE_HPP

#define  __VERBOSE_CHECK    -100
#define  __VERBOSE_NOTCHECK  0

/*Global variables to debug functions and programs to build reliable,
  correct and robust code
*/
#ifdef __VERBOSE_YES
/*Detail level to debug, -1 does not print debug messages.
 */
int        geiinparam_verboseMax  = -1;

/*Counter to print detail to debug
 */
int        geiinparam_verbose     =  0;

/*variable to keep identifier number of times
  that the task is being performed
*/
uintidx    geverboseui_idproc     =  0;
const char *geverbosepc_labelstep = "main";
#endif /*__VERBOSE_YES*/


#endif /*__VERBOSE_HPP*/

