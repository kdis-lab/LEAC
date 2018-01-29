/*! \file execution_time.hpp
 *
 * \brief execution time
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef EXECUTION_TIME_HPP
#define EXECUTION_TIME_HPP

#include <omp.h>

/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace runtime {


/*! ExecutionTime variable for count execution time.
 */
typedef double ExecutionTime;

inline ExecutionTime initialize()
{
  return ExecutionTime(0.0);
}

inline ExecutionTime start()
{
  return omp_get_wtime();
}

inline void stop(ExecutionTime& aoexetime_time)
{
  aoexetime_time = omp_get_wtime() - aoexetime_time;
}

inline ExecutionTime elapsedTime(const ExecutionTime& aoexetime_time)
{
  return omp_get_wtime() - aoexetime_time;
}

inline ExecutionTime getTime(ExecutionTime& aoexetime_time)
{
  return aoexetime_time;
}

} /*END namespace runtime 
   */

#endif  /*EXECUTION_TIME_HPP*/

