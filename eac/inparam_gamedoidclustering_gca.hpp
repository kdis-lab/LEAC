/*! \file inparam_gamedoidclustering_gca.hpp
 *
 * \brief Definition of GCA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GAMEDOIDCLUSTERING_GCA_HPP
#define IN_PARAM_GAMEDOIDCLUSTERING_GCA_HPP

#include "inparam_gaclustering_pcpm_fixedk.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
    
/*! \class InParamGAMedoidClusteringGCA 
  \brief Input parameter for GCA algorithm k-medoid clustering \cite Lucasius:etal:GAclusteringMedoid:GCA:1993
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGAMedoidClusteringGCA:
    public InParamGAClusteringProbCProbMFixedK<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGAMedoidClusteringGCA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   ) :
    InParamGAClusteringProbCProbMFixedK<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) {}
  ~InParamGAMedoidClusteringGCA() {}
  
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
    InParamGAClusteringProbCProbMFixedK
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

#endif /*IN_PARAM_GAMEDOIDCLUSTERING_GCA_HPP*/
