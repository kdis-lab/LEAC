/*! \file inparam_tgca.hpp
 *
 * \brief Definition of TGCA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_TGCA_HPP__
#define __IN_PARAM_TGCA_HPP__

#include "inparam_probcprobm_rangek.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
  
/*! \class InParamTGCA
  \brief Input parameter for TGCA a two-stage genetic algorithm \cite He:Tan:GAclusteringVarK:TGCA:2012
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamTGCA
  : public InParamPcPmRk
<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamTGCA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamPcPmRk
      <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , _it_kmeansNumMaxIter(INPARAMCLUSTERING_DEFAULT_MAX_ITER)
    , _ui_kmeansNumMinThreshold(INPARAMCLUSTERING_DEFAULT_MIN_THRESHOLD)
  {}

  ~InParamTGCA() {}
  
  inline void setKmeansNumMaxIter(COMMON_IDOMAIN aiiT_numMaxIterKmeans) 
  {
    _it_kmeansNumMaxIter = aiiT_numMaxIterKmeans;
  }

  inline COMMON_IDOMAIN getKmeansNumMaxIter()
  {
    return _it_kmeansNumMaxIter;
  }

  inline void setKmeansMinThreshold(uintidx aiui_numMinThreshold) 
  {
    this->_ui_kmeansNumMinThreshold = aiui_numMinThreshold;
  }

  inline const uintidx getKmeansMinThreshold() const 
  {
    return this->_ui_kmeansNumMinThreshold;
  }
  
  inline void setNumSubpopulationsCross(uintidx aiui_numSubpopulationsCross) 
  {
    _ui_numSubpopulationsCross = aiui_numSubpopulationsCross;
  }

  inline uintidx getNumSubpopulationsCross() const 
  {
    return _ui_numSubpopulationsCross;
  }
  
  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {

    InParamPcPmRk
      <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
      ::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_k-means iterations"            
		 << aic_separator << _it_kmeansNumMaxIter;
    aipf_outFile << aic_separator << "_minimum threshold" 
		 << aic_separator << _ui_kmeansNumMinThreshold;
    aipf_outFile << aic_separator << "_number subpopulations" 
		 << aic_separator << _ui_numSubpopulationsCross;
    
  }

protected:
  
  COMMON_IDOMAIN  _it_kmeansNumMaxIter;
  uintidx         _ui_kmeansNumMinThreshold;
  uintidx         _ui_numSubpopulationsCross;
  
};

}

#endif /*__IN_PARAM_TGCA_HPP__*/
