/*! \file inparam_adaptiveprobcprobm.hpp
 *
 * \brief Definition of GAGR program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_PADAPTIVE_HPP__
#define __IN_PARAM_PADAPTIVE_HPP__

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

/*! \class InParamAdaptivePcPm
  \brief Input parameter for EA with adaptive probabilities of crossover and mutation \cite Chang:etal:GAclustering:GAGR:2009
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamAdaptivePcPm
  : public InParamGAClustering
  , public InParamFk<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamAdaptivePcPm
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamGAClustering
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamFk<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}

  ~InParamAdaptivePcPm() {}

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering::print(aipf_outFile,aic_separator);
    InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }
protected:
};


} /* END namespace inout
   */

#endif /*__IN_PARAM_PADAPTIVE_HPP__*/
