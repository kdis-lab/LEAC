/*! \file inparam_clustering_medoidspam.hpp
 *
 * \brief Definition of PAM program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef IN_PARAM_CLUSTERING_MEDOIDSPAM_HPP
#define IN_PARAM_CLUSTERING_MEDOIDSPAM_HPP

#include "inparam_clustering_max_iter.hpp"
#include "inparam_fixedk.hpp"
#include "inparam_readinst.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \author Hermes Robles Berumen
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

/*! \class InParamClusteringMedoidsPAM
  \brief Input parameter for Medoids PAM algorithm \cite Kaufman:Rousseeuw:Book:ClusterAnalysis:1990
*/
template <typename T_CLUSTERIDX,
	  typename T_FEATURE,
	  typename T_INSTANCES_CLUSTER_K
	 > 
class InParamClusteringMedoidsPAM
  : public InParamClusteringMaxIter
  , public InParamFk<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{ 
public:
  InParamClusteringMedoidsPAM
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm) 
    : InParamClusteringMaxIter
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamFk<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
   {}
  ~InParamClusteringMedoidsPAM() {}

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
      InParamClusteringMaxIter::print(aipf_outFile,aic_separator);
      InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator); 
      InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }
protected:
  
};

  
} /*END namespace inout 
   */

#endif /*IN_PARAM_CLUSTERING_MEDOIDSPAM_HPP*/
