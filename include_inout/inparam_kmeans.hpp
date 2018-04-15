/*! \file inparam_kmeans.hpp
 *
 * \brief Definition of clustering clasic program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_KMEANS_HPP__
#define __IN_PARAM_KMEANS_HPP__

#include "inparam_clustering_max_iter.hpp"
#include "inparam_fixedk.hpp"
#include "inparam_readinst.hpp"

#define __INPARAM_KMEANS__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
/*! \class InParamKmeans
  \brief Input parameter for traditional algorithm clustering \cite MacQueen:ClusterAnalysis:KMeans:1967
*/
  template <typename T_CLUSTERIDX,
	    typename T_FEATURE,         
	    typename T_FEATURE_SUM,
	    typename T_INSTANCES_CLUSTER_K 
	    > 
class InParamKmeans
    : public InParamClusteringMaxIter
    , public InParamFk<T_CLUSTERIDX>
    , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamKmeans
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm) 
    : InParamClusteringMaxIter
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamFk<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , _ui_numMinThreshold(INPARAMCLUSTERING_DEFAULT_MIN_THRESHOLD)
  {}
  ~InParamKmeans() {}

  inline void setMinThreshold(uintidx aiui_numMinThreshold) 
  {
    this->_ui_numMinThreshold = aiui_numMinThreshold;
  }

  inline const uintidx getMinThreshold() const 
  {
    return this->_ui_numMinThreshold;
  }

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamClusteringMaxIter::print(aipf_outFile,aic_separator);
    InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_minimum threshold" 
		 << aic_separator <<  this->getMinThreshold();
  }
  
protected:
  uintidx _ui_numMinThreshold; 
}; 


} /*END namespace inparam 
   */

#endif /*__IN_PARAM_KMEANS_HPP__*/
