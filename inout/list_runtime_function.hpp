/*! \file list_runtime_function.hpp
 *
 * \brief list runtime function
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __LIST_RUNTIME_FUNCTION_HPP
#define __LIST_RUNTIME_FUNCTION_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "interval_positiveintegers.hpp"
#include "runtime_function.hpp"


/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace runtime {
  

/*! \class ListRuntimeFunction
  \brief ListRuntimeFunction  store metrics in run time and count genetations  
*/
template < typename T_INTEGERDOMAIN >
class ListRuntimeFunction {
public:
  ListRuntimeFunction
  (T_INTEGERDOMAIN aiiT_sizeInterval, 
   const std::string& ais_labelX, 
   const std::string& ais_labelY
   ): s_labelX(ais_labelX)

    , s_labelY(ais_labelY)
    , pIntervalZT_domain( new IntervalPositiveInteger<T_INTEGERDOMAIN>(aiiT_sizeInterval) )
    , _vector_ptRuntimeFunction()
    , _c_plotStatSeparator('\t')
  {}

  ~ListRuntimeFunction() 
  {
    for ( auto  liter_functionHist: _vector_ptRuntimeFunction ) 
      delete liter_functionHist;

    delete pIntervalZT_domain;
  }

  void  initialize()
  {
    pIntervalZT_domain->initialize();
    if ( _vector_ptRuntimeFunction.size() > 0 ) {
      for ( auto  liter_functionHist: _vector_ptRuntimeFunction ) 
	delete liter_functionHist;
       _vector_ptRuntimeFunction.clear();
    }
  }

  void addFuntion(RuntimeFunction* aifh_functionHistNew)
  {    
    this->_vector_ptRuntimeFunction.push_back(aifh_functionHistNew);
  }


  inline std::vector<RuntimeFunction* >&  getListRuntimeFunction()
  {
    return this->_vector_ptRuntimeFunction;
  }

  
  inline IntervalPositiveInteger<T_INTEGERDOMAIN>* getDomainInterval()	
  {
    return this->pIntervalZT_domain;
  }

  inline T_INTEGERDOMAIN getDomainSize()		
  {
    return this->pIntervalZT_domain->getSizeInterval();
  }

  inline T_INTEGERDOMAIN getDomainLowerBound()
  {
    return this->pIntervalZT_domain->getLowerBound();
  }

  inline const  T_INTEGERDOMAIN getDomainUpperBound() const  
  {
    return this->pIntervalZT_domain->getUpperBound(); 
  }

  inline void setDomainSize(T_INTEGERDOMAIN aiT_integerDomain)
  {
    this->pIntervalZT_domain->setSizeInterval(aiT_integerDomain);
  }

  inline void increaseDomainUpperBound()	
  {
    this->pIntervalZT_domain->increaseUpperBound();
  }

  inline std::string getLabelX()	
  {
    return this->s_labelX;
  }

  inline std::string getLabelY()	
  {
    return this->s_labelY;
  }
  
  inline uintidx getNumberFuntions()	
  {
    return (uintidx)  this->_vector_ptRuntimeFunction->size();
  }

  
  std::string getHeaderFuntions()
  {
    std::string los_headerFuntions(50,'\0');

    los_headerFuntions.assign(1,'#');
    los_headerFuntions.append(s_labelX);
    for (uintidx lui_i = 0; lui_i < _vector_ptRuntimeFunction.size(); lui_i++) {
      los_headerFuntions.append(1,_c_plotStatSeparator);
      los_headerFuntions.append(_vector_ptRuntimeFunction.at(lui_i)->getName());
    }
    
    return los_headerFuntions;
  }
  
  inline
  const  char getPlotStatSeparator()  const
  {
    return _c_plotStatSeparator;
  }
  
  template < typename U_INTEGERDOMAIN >
  friend std::ostream& operator<< 
  (std::ostream& out, const ListRuntimeFunction<U_INTEGERDOMAIN> &ailfhT_listFuntionHist);

private:
  
  std::string s_labelX;
  std::string s_labelY; 
  IntervalPositiveInteger<T_INTEGERDOMAIN> *pIntervalZT_domain; 
  std::vector<RuntimeFunction* > _vector_ptRuntimeFunction;
  char        _c_plotStatSeparator;
  
};
  

template<typename T_INTEGERDOMAIN>
std::ostream& operator<< (std::ostream& out, const ListRuntimeFunction<T_INTEGERDOMAIN> &ailfhT_listFuntionHist) 
{
  out << ailfhT_listFuntionHist.getDomainUpperBound();
  for ( auto  liter_functionHist: ailfhT_listFuntionHist._vector_ptRuntimeFunction ) {
    out << ailfhT_listFuntionHist.getPlotStatSeparator()
	<< *liter_functionHist;
  }
  out <<  '\n';
  return out;
}

} /*END namespace executiontime 
   */

#endif /*__LIST_RUNTIME_FUNCTION_HPP*/
