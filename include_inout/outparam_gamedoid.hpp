/*! \file outparam_gamedoid.hpp
 *
 * \brief out parameters GA-Medoid
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __OUT_PARAM_GA_MEDOID_HPP__
#define __OUT_PARAM_GA_MEDOID_HPP__


#include "outparam_gac.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
/*! \class OutParamGAC
  \brief Output Output parameters of the evolutionary algorithm
*/
template < typename T_METRIC,
	   typename T_CLUSTERIDX
	  >
class OutParamGAMedoid:
  public OutParamGAC<T_METRIC,T_CLUSTERIDX> {
public:
  OutParamGAMedoid(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc):
    OutParamGAC<T_METRIC,T_CLUSTERIDX>(aienum_usedObjectiveFunc)   
  {
    this->initialize(-1);
  }

  virtual ~OutParamGAMedoid() {}


  inline void setStringMedoid(std::string aistr_medoid)  
  {
    _str_medoid = aistr_medoid;
  }

  inline  std::string getStringMedoid()  
  {
    return _str_medoid;
  }


  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    OutParamGAC<T_METRIC,T_CLUSTERIDX>::print(aipf_outFile);
   
    aipf_outFile 
      << aic_separator << "_medoid" 
      << aic_separator << _str_medoid;

  }

protected:  

  std::string            _str_medoid;


}; /*OutParamGAMedoid*/

} /*END namespace inout
   */

#endif /*__OUT_PARAM_GA_MEDOID_HPP__*/
