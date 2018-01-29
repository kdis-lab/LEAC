/*! \file runtime_statistic.hpp
 *
 * \brief runtime statistic
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __RUNTIME_STATISTIC_HPP
#define __RUNTIME_STATISTIC_HPP

#include <vector>
#include <algorithm>    // std::min_element, std::max_element

#include "runtime_function_value.hpp"


/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  runtime {
  
  
#define STATISTICAL_ALL_MEASURES 6 
#define STATISTICAL_MINIMUM  0
#define STATISTICAL_MAXIMUM  1
#define STATISTICAL_AVERAGE  2
#define STATISTICAL_SUM      3
#define STATISTICAL_VARIANCE 4
#define STATISTICAL_STD_DEVIATION 5

const char* statistical_name[] = {"Min", "Max", "Avg", "Sum",  "Var", "StdD"};

/*! \class RuntimeFunctionStat
  \brief Function statistical for metrics in run time
*/
template <class T_VALUE>
class RuntimeFunctionStat:
  public RuntimeFunctionValue<T_VALUE> 
{
public:
  RuntimeFunctionStat
  (const char         aic_idStatistical, 
   const std::string& ais_algorithmoName, 
   char               aic_typeFunc
   ): RuntimeFunctionValue<T_VALUE>
      (statistical_name[(int) aic_idStatistical], ais_algorithmoName,aic_typeFunc)
    , c_idStatistical(aic_idStatistical)
  { }

  virtual ~RuntimeFunctionStat() {}

  void evaluate(std::vector<T_VALUE> &aivectorT_data) 
  {
    switch (this->c_idStatistical) {
    case STATISTICAL_MINIMUM:
      this->t_value = 
	*std::min_element(std::begin(aivectorT_data),std::end(aivectorT_data));
      break;
    case STATISTICAL_MAXIMUM:
      this->t_value = 
	*std::max_element(std::begin(aivectorT_data),std::end(aivectorT_data)); 
      break;
    case STATISTICAL_AVERAGE:
      this->t_value = _avg(aivectorT_data);
      break;
    case STATISTICAL_SUM:
      this->t_value = _sum(aivectorT_data);
      break;
    case STATISTICAL_VARIANCE:
      this->t_value = _var(aivectorT_data);
      break;
    case STATISTICAL_STD_DEVIATION:
      this->t_value = _devstd(aivectorT_data);
      break;
    default:
      break;
    }
  }
private:
  
  T_VALUE _sum(std::vector<T_VALUE> &aivectorT_data)
  {
    return std::accumulate(std::begin(aivectorT_data),std::end(aivectorT_data),0);
  } 

  T_VALUE _avg(std::vector<T_VALUE> &aivectorT_data)
  {
    return (_sum(aivectorT_data) / aivectorT_data.size()); 
  } 

  T_VALUE _var(std::vector<T_VALUE> &aivectorT_data)
  {
    T_VALUE lT_avg = _avg(aivectorT_data);
    T_VALUE lT_accum = 0.0;
    std::for_each (std::begin(aivectorT_data), std::end(aivectorT_data), [&](const T_VALUE liter_data) {
	lT_accum += (liter_data - lT_avg) * (liter_data - lT_avg);
      });
    return (lT_accum / (aivectorT_data.size()- 1.0)); 
  } 

  T_VALUE _devstd(std::vector<T_VALUE> &aivectorT_data)
  {
    return std::sqrt(_var(aivectorT_data));
  }

protected:
  int c_idStatistical;

};
  
/*! \fn void functionhiststat_evaluateAll(RuntimeFunctionStat<T_VALUE>  *aofhs_statObjectiveFunc[], std::vector<T_VALUE> &aivectorT_data) 
  \brief Calculate statistical for metrics in run time 
  \details 
  \param aofhs_statObjectiveFunc a array of RuntimeFunctionStat<T_VALUE>
  \param aivectorT_data a std::vector<T_VALUE>
*/
template <class T_VALUE>
void functionhiststat_evaluateAll
(RuntimeFunctionStat<T_VALUE>  *aofhs_statObjectiveFunc[],
 std::vector<T_VALUE>       &aivectorT_data
 )
{ //BEGIN functionhiststat_evaluateAll 
  aofhs_statObjectiveFunc[STATISTICAL_MINIMUM]->setValue
    (*std::min_element(std::begin(aivectorT_data),std::end(aivectorT_data)));
  aofhs_statObjectiveFunc[STATISTICAL_MAXIMUM]->setValue
    (*std::max_element(std::begin(aivectorT_data),std::end(aivectorT_data))); 

  T_VALUE lT_sum = 
    std::accumulate(std::begin(aivectorT_data),std::end(aivectorT_data),0);
  aofhs_statObjectiveFunc[STATISTICAL_SUM]->setValue(lT_sum);

  T_VALUE lT_avg = lT_sum / aivectorT_data.size();
  aofhs_statObjectiveFunc[STATISTICAL_AVERAGE]->setValue(lT_avg);
  
  T_VALUE lT_accum = 0.0;
  std::for_each
    (std::begin(aivectorT_data), 
     std::end(aivectorT_data), 
     [&](const T_VALUE liter_data) 
     {
      lT_accum += (liter_data - lT_avg) * (liter_data - lT_avg);
     }
     );

  T_VALUE lT_var = lT_accum / (aivectorT_data.size()- 1.0); 
  aofhs_statObjectiveFunc[STATISTICAL_VARIANCE]->setValue(lT_var);

  aofhs_statObjectiveFunc[STATISTICAL_STD_DEVIATION]->setValue(std::sqrt(lT_var));
 
} //ALL functionhiststat_evaluateAll 


} /*END namespace runtime 
   */

  
#endif /*__RUNTIME_STATISTIC_HPP*/
