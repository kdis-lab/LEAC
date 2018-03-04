/*! \file inparam_withoutprobcprobm.hpp
 *
 * \brief Definition of GA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_WITHOUT_PCPM_HPP__
#define __IN_PARAM_WITHOUT_PCPM_HPP__

#include "inparam_fixedk.hpp"
#include "inparam_gaclustering.hpp"
#include "inparam_readinst.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamWithoutPcPm
  \brief Input parameter for GA without probability crossover (Pc) and mutation (Pm)\cite Bezdek:etal:GAclustering:GA:1994
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_BITSIZE,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	  > 
class InParamWithoutPcPm
  : public InParamGAClustering
  , public InParamFk<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamWithoutPcPm
  (const std::string&   ais_algorithmoName,
   const std::string&   ais_algorithmoAuthor,
   InParam_algTypeOut   aiato_algTypeOut,
   int                  aii_opNorm)
    :  InParamGAClustering
         (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) 
    ,  InParamFk<T_CLUSTERIDX>()
    ,  InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}
  
  inline void setSizeMatingPool(uintidx aist_sizeMatingPool) {
    this->st_sizeMatingPool = aist_sizeMatingPool;
  }

  inline uintidx getSizeMatingPool() 
  {
    return this->st_sizeMatingPool;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering
      ::print(aipf_outFile,aic_separator);
    InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_size mating pool" 
		 << aic_separator << this->st_sizeMatingPool;
  }
protected:
 
  uintidx st_sizeMatingPool;
  
};

} /* END namespace inout
   */

#endif /*__IN_PARAM_WITHOUT_PCPM_HPP__*/
