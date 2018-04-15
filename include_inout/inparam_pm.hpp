/*! \file inparam_probm.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_PROBM_HPP__
#define __IN_PARAM_PROBM_HPP__

#include "inparam_gaclustering.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
/*! \class InParamPm
  \brief Input parameter for GA with only probability mutation (Pm)
*/
  template <  typename T_REAL>
class InParamPm
    :  public InParamGAClustering
{
public:
  InParamPm
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   ):
    InParamGAClustering
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) {}

  ~InParamPm() {}

  inline void setProbMutation(T_REAL ait_probMutation) {
    this->_t_probMutation = ait_probMutation;
  }

  inline const T_REAL getProbMutation() const  {
    return this->_t_probMutation;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering::print(aipf_outFile,aic_separator);
   
    aipf_outFile << aic_separator << "_probability mutation" 
		 << aic_separator << this->_t_probMutation;
  }

protected:
  T_REAL  _t_probMutation;
}; 

} /* END namespace inout
   */

#endif /*__IN_PARAM_PROBM_HPP__*/
