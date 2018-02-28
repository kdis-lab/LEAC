/*! \file inparam_gaclustering_withoutpcpm.hpp
 *
 * \brief Definition of GA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GACLUSTERING_WITHOUT_PCPM_HPP
#define IN_PARAM_GACLUSTERING_WITHOUT_PCPM_HPP

#include "inparam_fixedk.hpp"
#include "inparam_gaclustering.hpp"
#include "inparam_definedatatypes.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamGAClusteringWithoutProbCProbM
  \brief Input parameter for GA without probability crossover (Pc) and mutation (Pm)\cite Bezdek:etal:GAclustering:GA:1994
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_BITSIZE,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	  > 
class InParamGAClusteringWithoutProbCProbM
  : public InParamGAClustering
  , public InParamFixedK<T_CLUSTERIDX>
  , public InParamDefineBitSizeFeatFeatSumInstK<T_BITSIZE,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGAClusteringWithoutProbCProbM
  (const std::string&   ais_algorithmoName,
   const std::string&   ais_algorithmoAuthor,
   InParam_algTypeOut   aiato_algTypeOut,
   int                  aii_opNorm)
    :  InParamGAClustering
         (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm) 
    ,  InParamFixedK<T_CLUSTERIDX>()
    ,  InParamDefineBitSizeFeatFeatSumInstK<T_BITSIZE,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>()
  {}
  
  inline void setSizeMatingPool(uintidx aist_sizeMatingPool) {
    this->st_sizeMatingPool = aist_sizeMatingPool;
  }

  inline uintidx getSizeMatingPool() 
  {
    return this->st_sizeMatingPool;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering
      ::print(aipf_outFile,aic_separator);
    InParamFixedK<T_CLUSTERIDX>::print(aipf_outFile,aic_separator); 
    aipf_outFile << aic_separator << "_size mating pool" 
		 << aic_separator << this->st_sizeMatingPool;
  }
protected:
 
  uintidx st_sizeMatingPool;
  
};

} /* END namespace inout
   */

#endif /*IN_PARAM_GACLUSTERING_WITHOUT_PCPM_HPP*/
