/*! \file inparam_isodata.hpp
 *
 * \brief Definition of ISODATA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_ISODATA_HPP__
#define __IN_PARAM_ISODATA_HPP__

#include "inparam_clustering_max_iter.hpp"
#include "inparam_fixedk.hpp"
#include "inparam_readinst.hpp"

#define __INPARAM_ISODATA__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamIsoData
  \brief Input parameter for IsoData algorithm
*/
template < typename T_CLUSTERIDX,
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamIsoData
  : public InParamClusteringMaxIter
  , public InParamFk<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:

  InParamIsoData
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut, 
   int         aii_opNorm)
    : InParamClusteringMaxIter
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm)
    , InParamFk<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , b_optimalInitializacion(false)
    , t_epsilon(INPARAMCLUSTERING_DEFAULT_EPSILON) 
  {}
  
  ~InParamIsoData()
  {}

  inline void setIsOptimalInitializacion(bool aib_optimalInitializacion) 
  {
    this->b_optimalInitializacion = aib_optimalInitializacion;
  }

  inline bool getIsOptimalInitializacion() 
  {
    return this->b_optimalInitializacion;
  }

  inline void setEpsilon(T_REAL aifpT_epsilon) 
  {
    this->t_epsilon = aifpT_epsilon;
  }

  inline T_REAL getEpsilon() 
  {
    return this->t_epsilon;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamClusteringMaxIter::print(aipf_outFile);
    InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_epsilon" 
		 << aic_separator << this->t_epsilon;
  }
  
protected:
  bool       b_optimalInitializacion;
  T_REAL     t_epsilon;
};
  
} /*END namespace inparam 
   */

#endif /*__IN_PARAM_ISODATA_HPP__*/
