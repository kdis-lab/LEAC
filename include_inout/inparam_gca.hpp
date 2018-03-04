/*! \file inparam_gca.hpp
 *
 * \brief Definition of GCA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_GCA_HPP__
#define __IN_PARAM_GCA_HPP__

#include "inparam_probcprobm_fixedk.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
    
/*! \class InParamGCA 
  \brief Input parameter for GCA algorithm k-medoid clustering \cite Lucasius:etal:GAclusteringMedoid:GCA:1993
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGCA:
    public InParamPcPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGCA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   ) :
    InParamPcPmFk<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) {}
  ~InParamGCA() {}
  
  inline void setProbMixMutation(T_REAL aiT_probMixMutation) 
  {
    this->t_probMixMutation = aiT_probMixMutation;
  }

  inline T_REAL getProbMixMutation() 
  {
    return this->t_probMixMutation;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamPcPmFk
      <T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
      ::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_mix mutation probability" 
		 << aic_separator << this->t_probMixMutation;      
  }

protected:
  
  T_REAL  t_probMixMutation;
  
};

  
} /* END namespace inparam
   */

#endif /*__IN_PARAM_GCA_HPP__*/
