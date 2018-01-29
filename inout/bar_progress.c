/*! \file bar_progress.c
 *
 * \brief implementation bar progress
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#include <stdio.h>
#include "bar_progress.h"

#define BARPROGRESS_SIZE 50

void barprogress_update(int ai_numberOfRun, int ai_totalRuns)
{
  int li_i;
  int li_percentage;
  int li_currentSizeBar;

  li_percentage = 100 * ai_numberOfRun / ai_totalRuns;
  li_currentSizeBar = li_percentage * BARPROGRESS_SIZE  / 100;
  printf("Computing: %d%% [",li_percentage);
  for ( li_i = 0; li_i < li_currentSizeBar; li_i++)
    printf("#");
  for ( li_i = li_currentSizeBar; li_i < BARPROGRESS_SIZE; li_i++)
    printf(".");
  printf("] %d/%d",ai_numberOfRun,ai_totalRuns);
  if ( ai_numberOfRun < ai_totalRuns) {
    printf("\r");
    fflush(stdout);
      
  }
  else {
    printf("\n");
  }
}
