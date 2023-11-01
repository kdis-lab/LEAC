/*! \file inparam_pcpmvk.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __INPARAM_PCPMFREQVK_HPP__
#define __INPARAM_PCPMFREQVK_HPP__

#include "inparam_pcpm.hpp"
#include "inparam_variablek.hpp"
#include "inparam_readinst.hpp"

#define __INPARAM_PCPMFREQVK__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

/*! \class InParamPcPmFreqVk
  \brief Input parameter for GA with probability Pc, Pm and Variable K 
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_INSTANCE_FREQUENCY
	   >
class InParamPcPmFreqVk
  : public InParamPcPm<T_REAL>
  , public InParamVk<T_CLUSTERIDX>
  , public InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY>
  //, public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamPcPmFreqVk
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   ) 
    : InParamPcPm<T_REAL>
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) 
    , InParamVk<T_CLUSTERIDX>()
    , InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY>()
    //, InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}

  ~InParamPcPmFreqVk() {}

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamPcPm<T_REAL>::print(aipf_outFile,aic_separator);
    InParamVk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY>::print(aipf_outFile,aic_separator);
    //InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }  
};

} /*END namespace inout 
   */

#endif /*__INPARAM_PCPMFREQVK_HPP__*/
