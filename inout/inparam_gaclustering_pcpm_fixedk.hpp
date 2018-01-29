/*! \file inparam_gaclustering_pcpm_fixedk.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GACLUSTERING_PROB_FIXEDK_HPP
#define IN_PARAM_GACLUSTERING_PROB_FIXEDK_HPP

#include "inparam_gaclustering_pcpm.hpp"
#include "inparam_fixedk.hpp"
#include "inparam_definedatatypes.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {


/*! \class InParamGAClusteringProbCProbMFixedK
  \brief Input parameter for GA with probability Pc, Pm and Fixed K 
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGAClusteringProbCProbMFixedK
  : public InParamGAClusteringProbCProbM<T_REAL>
  , public InParamFixedK<T_CLUSTERIDX>
  , public InParamDefineFeatFeatSumInstK<T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGAClusteringProbCProbMFixedK
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   ) 
    : InParamGAClusteringProbCProbM<T_REAL>
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) 
    , InParamFixedK<T_CLUSTERIDX>()
    , InParamDefineFeatFeatSumInstK<T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>()
  {}

  ~InParamGAClusteringProbCProbMFixedK() {}

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClusteringProbCProbM<T_REAL>
      ::print(aipf_outFile,aic_separator);
    InParamFixedK<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }  
};

} /*END namespace inout 
   */

#endif /*IN_PARAM_GACLUSTERING_PROB_FIXEDK_HPP*/
