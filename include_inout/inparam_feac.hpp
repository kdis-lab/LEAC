/*! \file inparam_feac.hpp
 *
 * \brief Definition of FEAC program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_FEAC_HPP__
#define __IN_PARAM_FEAC_HPP__

#include "inparam_gaclustering.hpp"
#include "inparam_rangek.hpp"
#include "inparam_readinst.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

#define INPARAMCLUSTERING_FEAC_KMEANSNUMMAXITER   5
#define INPARAMCLUSTERING_FEAC_KMEANSMAXDIFFCENT  0.001
#define INPARAMCLUSTERING_FEAC_DEFAULT_DESIABLEFITNESS 1
  
/*! \class InParamFEAC
  \brief Input parameter Fast Evolutionary Algorithm for Clustering F-EAC \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006 
*/
template < typename T_FEATURE,
	   typename T_REAL,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamFEAC
  : public InParamGAClustering
  , public InParamRk<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamFEAC
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamGAClustering
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamRk<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , _it_kmeansNumMaxIter(INPARAMCLUSTERING_FEAC_KMEANSNUMMAXITER)
    , _rT_kmeansMaxDiffCent(INPARAMCLUSTERING_FEAC_KMEANSMAXDIFFCENT)
    , _rT_desiableObjetiveFunc(INPARAMCLUSTERING_FEAC_DEFAULT_DESIABLEFITNESS)
  {}

  ~InParamFEAC() {}

  
  inline void setKmeansNumMaxIter(COMMON_IDOMAIN aiiT_numMaxIterKmeans) 
  {
    _it_kmeansNumMaxIter = aiiT_numMaxIterKmeans;
  }

  inline COMMON_IDOMAIN getKmeansNumMaxIter()
  {
    return _it_kmeansNumMaxIter;
  }

  inline T_FEATURE getKmeansMaxDiffCent()
  {
    return _rT_kmeansMaxDiffCent;
  }

  inline void setKmeansMaxDiffCent(T_FEATURE airT_kmeansMaxDiffCent) 
  {
    _rT_kmeansMaxDiffCent = airT_kmeansMaxDiffCent;
  }

  inline T_REAL getDesiableObjetiveFunc()
  {
    return _rT_desiableObjetiveFunc;
  }

  inline void setDesiableObjetiveFunc(T_REAL airT_desiableObjetiveFunc) 
  {
    _rT_desiableObjetiveFunc = airT_desiableObjetiveFunc;
  }
 

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering::print(aipf_outFile,aic_separator);
      InParamRk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
      InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
      aipf_outFile << aic_separator << "_k-means iterations"            
		   << aic_separator << _it_kmeansNumMaxIter; 
      aipf_outFile << aic_separator << "_maximum difference centroids iterations"
		   << aic_separator << _rT_kmeansMaxDiffCent;
      aipf_outFile << aic_separator << "_desiable fitness"
		   << aic_separator << _rT_desiableObjetiveFunc;
  }

protected:
  COMMON_IDOMAIN  _it_kmeansNumMaxIter;
  T_FEATURE       _rT_kmeansMaxDiffCent;
  //!_rT_desiableObjetiveFunc
  /*!
    Value desiable fitness for stop algorith
  */
  T_REAL          _rT_desiableObjetiveFunc; 
}; 

} /*END namespace inout 
   */

#endif /*__IN_PARAM_FEAC_HPP__*/
