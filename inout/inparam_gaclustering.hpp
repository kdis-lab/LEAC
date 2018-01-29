/*! \file inparam_gaclustering.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GACLUSTERING_HPP
#define IN_PARAM_GACLUSTERING_HPP

#include "inparam_clustering.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamGAClustering
  \brief Input parameter for genetic and evolutionary algorithm 
*/
class InParamGAClustering:
    public InParamClustering
{
public:
  InParamGAClustering
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm):
    InParamClustering
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm) {}
 
  ~InParamGAClustering() {}

  inline void setNumMaxGenerations(COMMON_IDOMAIN aiT_numMaxGenerations) 
  {
    this->t_numMaxGenerations = aiT_numMaxGenerations;
  }

  inline COMMON_IDOMAIN getNumMaxGenerations() const 
  {
    return this->t_numMaxGenerations;
  }

  inline void setSizePopulation(uintidx aist_sizePopulation) 
  {
    this->st_sizePopulation = aist_sizePopulation;
  }

  inline uintidx getSizePopulation() const 
  {
    return this->st_sizePopulation;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamClustering::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_size population"   
		 << aic_separator << this->st_sizePopulation;
    aipf_outFile << aic_separator << "_number maximum generations" 
		 << aic_separator << this->t_numMaxGenerations;
  }
protected:
  uintidx          st_sizePopulation;
  COMMON_IDOMAIN t_numMaxGenerations;
}; 


} /*END namespace inout 
   */

#endif /*IN_PARAM_GACLUSTERING_HPP*/
