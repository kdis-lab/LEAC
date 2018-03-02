/*! \file inparam_gaclustering_gka.hpp
 *
 * \brief Definition of GKA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GACLUSTERING_GKA_HPP
#define IN_PARAM_GACLUSTERING_GKA_HPP

#include "inparam_fixedk.hpp"
#include "inparam_gaclustering_pm.hpp"
#include "inparam_readinst.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
    

/*! \class InParamGAClusteringGKA
  \brief Input parameter for GKA with only probability mutation (Pm) and Fixed K \cite Krishna:Murty:GAClustering:GKA:1999
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGAClusteringGKA
  : public InParamGAClusteringProbM<T_REAL>
  , public InParamFixedK<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamGAClusteringGKA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   )
    : InParamGAClusteringProbM<T_REAL>
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm)
    , InParamFixedK<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}

  ~InParamGAClusteringGKA() {}

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClusteringProbM<T_REAL>::print(aipf_outFile,aic_separator);
    InParamFixedK<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }

protected:
  
};

} /* END namespace inout
   */

#endif /*IN_PARAM_GACLUSTERING_GKA_HPP*/
