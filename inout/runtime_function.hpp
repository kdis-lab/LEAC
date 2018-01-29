/*! \file runtime_function.hpp
 *
 * \brief runtime function
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __RUNTIME_FUNCTION_HPP
#define __RUNTIME_FUNCTION_HPP

#include <iostream>
#include <string>

/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace runtime {
  
  
#define RUNTIMEFUNCTION_NOT_STORAGE        0
#define RUNTIMEFUNCTION_STORAGE_LINKEDLIST 1 
#define RUNTIMEFUNCTION_STORAGE_BUFFER     2

/*! \class RuntimeFunction
  \brief Function store metrics in run time
*/
class RuntimeFunction {
public:
  RuntimeFunction(const std::string& as_funcName, 
	       const std::string& as_algorithmoName, 
	       char               aic_typeFunc
	       )
    : _str_nameFunc(as_algorithmoName)
    , _c_typeFunc(aic_typeFunc)
  {
    _str_nameFunc.append(as_funcName);
  }

  virtual ~RuntimeFunction() {}

  inline std::string& getName() 
  {
    return this->_str_nameFunc;
  }

  inline bool isNotStorage()
  {
    return (this->_c_typeFunc == RUNTIMEFUNCTION_NOT_STORAGE ); 
  }

  inline bool isStorageBuffer()
  {
    return (this->_c_typeFunc ==  RUNTIMEFUNCTION_STORAGE_BUFFER); 
  }

  inline bool isStorageLinkedList() 
  {
    return (this->_c_typeFunc ==  RUNTIMEFUNCTION_STORAGE_LINKEDLIST); 
  }

  virtual std::ostream& print(std::ostream& os) const = 0;
  
protected:
  
  std::string _str_nameFunc;
  char        _c_typeFunc;
  
};

std::ostream& operator<<(std::ostream& os, const RuntimeFunction& airtf_runtimefunction)
{
  return airtf_runtimefunction.print(os);
}


} /*END namespace runtime 
   */
 
#endif  /*__RUNTIME_FUNCTION_HPP*/
