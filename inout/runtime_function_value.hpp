/*! \file runtime_function_value.hpp
 *
 * \brief runtime function value
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __RUNTIME_FUNCTION_VALUE_HPP
#define __RUNTIME_FUNCTION_VALUE_HPP

#include "runtime_function.hpp"


/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace runtime {  
  

/*! \class RuntimeFunctionValue
  \brief Function store metrics in run time
*/
template <class T_VALUE>
class RuntimeFunctionValue:
  public RuntimeFunction
{ //BEGIN CLASS  RuntimeFunctionValue
public:
  RuntimeFunctionValue(const std::string& as_funcName, 
		   const std::string& ais_algorithmoName, 
		   char               aic_typeFunc
		   ) :
    RuntimeFunction(as_funcName,ais_algorithmoName,aic_typeFunc),
    t_value(0)
  { }

  virtual ~RuntimeFunctionValue() {}

  inline void setValue(T_VALUE aiT_value)
  {
    this->t_value = aiT_value;
  }
  inline T_VALUE getValue()
  {
    return this->t_value;
  }

  virtual std::ostream& print(std::ostream& out) const
  {
    out << t_value;

    return out;
  }

protected:
  T_VALUE   t_value;
}; //END CLASS  RuntimeFunctionValue


} /*END namespace executiontime 
   */

  
#endif /*__RUNTIME_FUNCTION_VALUE_HPP*/
