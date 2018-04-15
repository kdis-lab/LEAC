/*! \file inparam_gaprototypesfk.hpp
 *
 * \brief Definition of GA-Prototypes program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_GAPROTOTYPESFK_HPP__
#define __IN_PARAM_GAPROTOTYPESFK_HPP__

#include "inparam_pcpmfk.hpp"

#undef  __INPARAM_PCPMFK__
#define __INPARAM_GAPROTOTYPESFK__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

/*! \class InParamGAPrototypesFk
  \brief Input parameter for GA-Prototypes algorithm \cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997
*/
template < typename T_BITSIZE, //-1, 0, 1, .., K
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGAPrototypesFk
  : public InParamPcPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGAPrototypesFk
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   ) :
    InParamPcPmFk
    <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) {}

  ~InParamGAPrototypesFk() {}

  inline void setPini(T_REAL airt_pini) 
  {
    this->_rt_pini = airt_pini;
  }

  inline T_REAL getPini()
  {
    return this->_rt_pini;
  }
  
  inline void setAlpha(T_REAL airt_alpha) 
  {
    this->_rt_alpha = airt_alpha;
  }

  inline T_REAL getAlpha()
  {
    return this->_rt_alpha;
  }
  
  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamPcPmFk
      <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
      ::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_rt_pini" 
		 << aic_separator << _rt_pini;
    aipf_outFile << aic_separator << "_rt_alpha" 
		 << aic_separator << _rt_alpha;
  }
  
protected:
  
  T_REAL _rt_pini;
  T_REAL _rt_alpha;
  
};

  
} /* END namespace inout
   */

#endif /*__IN_PARAM_GAPROTOTYPESFK_HPP__*/
