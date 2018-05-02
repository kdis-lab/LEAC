/*! \file outparam_clusteringalg.hpp
 *
 * \brief outparam clasic clustering
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __OUT_PARAM_CLUSTERING_ALG_HPP__
#define __OUT_PARAM_CLUSTERING_ALG_HPP__


#include "outparam_clustering.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace inout {

/*! \class OutParamDBSCAN
  \brief Output parameter for DBSCAN Algorithmo \cite Ester:Kriegel:Sander:Xu:Clustering:DBSCAN:1996
*/
template < typename T_METRIC,
	   typename T_MEMBERCLUSTER_IDX
	   >
class OutParamDBSCAN: 
  public OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX> {
public:
  OutParamDBSCAN(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc):
  OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX>
    ::OutParamClustering(aienum_usedObjectiveFunc)   
  {
    this->initialize(-1);
  }

  virtual ~OutParamDBSCAN() { }

  void initialize(int aii_numRunAlgorithm)
  {
    OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX>
    ::initialize(aii_numRunAlgorithm);
    this->_ui32t_numPointNoise = 0;
  }

  inline void setNumPointNoise(const  uint32_t aiui32t_numPointNoise)	
  {
    this->_ui32t_numPointNoise = aiui32t_numPointNoise;
  }

  inline const uintidx getNumPointNoise() const	
  {
    return this->_ui32t_numPointNoise;
  }

  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX>::print(aipf_outFile);
    aipf_outFile << aic_separator << "_number point noise" 
		 << aic_separator << this->_ui32t_numPointNoise;
  }
protected:

  uint32_t _ui32t_numPointNoise;

}; /*END class OutParamDBSCAN*/


/*! \class OutParamClusteringAlg
  \brief Output parameter for Kmeans Algorithmo \cite MacQueen:ClusterAnalysis:KMeans:1967
 */
template < typename T_METRIC,
	   typename T_MEMBERCLUSTER_IDX
	   >
class OutParamClusteringAlg: 
  public OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX> {
public:
  OutParamClusteringAlg(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc):
  OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX>
    ::OutParamClustering(aienum_usedObjectiveFunc)   
  {
    this->initialize(-1);
  }

  virtual ~OutParamClusteringAlg() { }

  void initialize(int aii_numRunAlgorithm)
  {
    OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX>
    ::initialize(aii_numRunAlgorithm);
    this->_uintidx_numThreshold = OUTPARAMCLUSTERING_INT_NaN;
  }

  inline void setNumThreshold(uintidx aist_numThreshold)	
  {
    this->_uintidx_numThreshold = aist_numThreshold;
  }

  inline uintidx& getNumThreshold()	
  {
    return this->_uintidx_numThreshold;
  }

  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    OutParamClustering<T_METRIC,T_MEMBERCLUSTER_IDX>::print(aipf_outFile);
    aipf_outFile << aic_separator << "_number threshold" 
		 << aic_separator << this->_uintidx_numThreshold;
  }
protected:

  uintidx  _uintidx_numThreshold;

}; /*END class OutParamClusteringAlg*/

} /*END namespace inout
   */

#endif /*__OUT_PARAM_CLUSTERING_ALG_HPP__*/
