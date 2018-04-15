/*! \file inparam_clustering_fcm.hpp
 *
 * \brief Definition of FCM program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_CLUSTERING_FCM_HPP
#define IN_PARAM_CLUSTERING_FCM_HPP

#include "inparam_clustering_max_iter.hpp"
#include "inparam_fixedk.hpp"
#include "inparam_readinst.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

/*! \class InParamClusteringFCM
  \brief Input parameter for FCM algorithm clustering \cite Bezdek:ClusterAnalysis:FCM:1974
*/
template < typename T_CLUSTERIDX,
	   typename T_REAL,
	   typename T_FEATURE,         
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamClusteringFCM
  : public InParamClusteringMaxIter
  , public InParamFk<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{ 
public:
  InParamClusteringFCM
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamClusteringMaxIter
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamFk<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , t_epsilon(INPARAMCLUSTERING_DEFAULT_EPSILON) 
    , t_weightingExponent(INPARAMCLUSTERING_DEFAULT_WEIGHTING_EXPONENT) 
  {}
  
  ~InParamClusteringFCM() {}

  inline void setWeightingExponent(T_REAL aiT_weightingExponent) 
  {
    this->t_weightingExponent = aiT_weightingExponent;
  }
  
  inline T_REAL getWeightingExponent() 
  {
    return this->t_weightingExponent;
  }

  inline void setEpsilon(T_REAL aiT_epsilon) 
  {
    this->t_epsilon = aiT_epsilon;
  }

  inline T_REAL getEpsilon() 
  {
    return this->t_epsilon;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',')  const
  {
    InParamClusteringMaxIter::print(aipf_outFile);
    InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_epsilon"   
		 << aic_separator << this->t_epsilon;
    aipf_outFile << aic_separator << "_weighting exponent"
		 << aic_separator << this->t_weightingExponent;
  }
protected:
  
  T_REAL   t_epsilon;
  T_REAL   t_weightingExponent;
  
}; 

  
} /*END namespace inout
   */

#endif /*IN_PARAM_CLUSTERING_FCM_HPP*/
