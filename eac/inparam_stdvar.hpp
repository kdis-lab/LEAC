/*! \file inparam_stdvar.hpp
 *
 * \brief Definition of standardization of variables program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_STDVAR_HPP
#define __IN_PARAM_STDVAR_HPP

#include "inparam.hpp"
#include "standardize_variable.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
/*! \class InParamStdVar
  \brief Input parameter for standardization of variables Algorithmo 
*/
class InParamStdVar: public InParamAlgorithmo
{
public: 
  InParamStdVar
  (std::string        ais_algorithmoName,
   std::string        ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut)
    : InParamAlgorithmo(ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut)
    , _i_standardizationVar(STD_VAR_Z0)
  {}

  ~InParamStdVar() {}

  inline void setStandardizationVar(int aii_standardizationVar)
  {
    _i_standardizationVar = aii_standardizationVar;
  }

  inline int getStandardizationVar() 
  { 
    return _i_standardizationVar;
  }
  
  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    const char  *las_optStandardizationVar[] = STD_VAR_NAME;
    
    InParamAlgorithmo::print(aipf_outFile,aic_separator);
      aipf_outFile
	<< aic_separator
	<< "_standardization of variables" 
     	<< aic_separator
	<< las_optStandardizationVar[this->_i_standardizationVar];
  }
  
protected:

  int          _i_standardizationVar;
  
}; /*InParamStdVar*/

} /* END namespace inout*/

#endif /*__IN_PARAM_STDVAR_HPP*/
