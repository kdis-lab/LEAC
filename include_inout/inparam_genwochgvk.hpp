/*! \file inparam_genwochgvk.hpp
 *
 * \brief Definition of GA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_GENERATIONS_WITHOUT_CHANGE_VK_HPP__
#define __IN_PARAM_GENERATIONS_WITHOUT_CHANGE_VK_HPP__

#include "inparam_pcpmvk.hpp"

#undef  __INPARAM_PCPMVK__
#define __INPARAM_GENWOCHGVK__
 
/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamGenWOChgVk
  \brief Input parameter Document clustering into an unknown number of clusters using a genetic algorithm \cite Casillas:etal:GAclusteringVarK:GA:2003
*/
  template < typename T_BITSIZE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGenWOChgVk
  : public InParamPcPmVk
<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGenWOChgVk
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   ) : InParamPcPmVk
       <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) {}
   
  ~InParamGenWOChgVk() {}

  inline void setNumNotChangeStop(COMMON_IDOMAIN aiit_numNotChangeStop) 
  {
    this->_it_numNotChangeStop = aiit_numNotChangeStop;
  }

   inline COMMON_IDOMAIN getNumNotChangeStop() 
  {
    return this->_it_numNotChangeStop;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamPcPmVk
      <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
      ::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "NumNotChangeStop"   
		 << aic_separator << this->_it_numNotChangeStop; 
  }

protected:

 COMMON_IDOMAIN _it_numNotChangeStop;
  
}; /*InParamGenWOChgVk*/


} /*END namespace inout 
   */

#endif /*__IN_PARAM_GENERATIONS_WITHOUT_CHANGE_VK_HPP__*/
