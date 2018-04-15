/*! \file inparam_pcpm.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_GACLUSTERING_PROBC_PROBM_HPP__
#define __IN_PARAM_GACLUSTERING_PROBC_PROBM_HPP__

#include "inparam_gaclustering.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  

/*! \class InParamPcPm
  \brief Input parameter for GA with probability crossover (Pc) and  mutation (Pm)
*/
template < typename T_REAL>
class InParamPcPm
  : public InParamGAClustering
{
public:
  InParamPcPm
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamGAClustering
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm)
  {}

  ~InParamPcPm() {}

  inline void setProbCrossover(T_REAL aiT_probCrossover) 
  {
    this->t_probCrossover = aiT_probCrossover;
  }

  inline T_REAL getProbCrossover() {
    return this->t_probCrossover;
  }

  inline void setProbMutation(T_REAL aiT_probMutation) {
    this->t_probMutation = aiT_probMutation;
  }

  inline const T_REAL getProbMutation() const {
    return this->t_probMutation;
  }

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_probability crossover"   
		 << aic_separator << this->t_probCrossover;
    aipf_outFile << aic_separator << "_probability mutation" 
		 << aic_separator << this->t_probMutation;
  }

protected:
  
  T_REAL  t_probCrossover;
  T_REAL  t_probMutation;
  
}; /*InParamPcPm*/

} /* END namespace inout
   */

#endif /*__IN_PARAM_GACLUSTERING_PROBC_PROBM_HPP__*/
