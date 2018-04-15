/*! \file bar_progress.h
 *
 * \brief bar progress
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef BARPROGRESS_H
#define BARPROGRESS_H

#ifdef __cplusplus
extern "C" {
#endif

  void barprogress_update(int ai_numberOfRun, int ai_totalRuns);

#ifdef __cplusplus
}
#endif 

#endif /* BARPROGRESS_H */

