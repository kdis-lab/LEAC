/*! \file inparam_gaclustering_padaptive.hpp
 *
 * \brief Definition of GAGR program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GACLUSTERING_PADAPTIVE_HPP
#define IN_PARAM_GACLUSTERING_PADAPTIVE_HPP

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

/*! \class InParamGAClusteringProbAdaptive
  \brief Input parameter for EA with adaptive probabilities of crossover and mutation \cite Chang:etal:GAclustering:GAGR:2009
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGAClusteringProbAdaptive
  : public InParamGAClustering
  , public InParamFixedK<T_CLUSTERIDX>
  , public InParamDefineFeatFeatSumInstK<T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamGAClusteringProbAdaptive
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamGAClustering
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamFixedK<T_CLUSTERIDX>()
    , InParamDefineFeatFeatSumInstK<T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>()
  {}

  ~InParamGAClusteringProbAdaptive() {}

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGAClustering::print(aipf_outFile,aic_separator);
    InParamFixedK<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }
protected:
};


} /* END namespace inout
   */

#endif /*IN_PARAM_GACLUSTERING_PADAPTIVE_HPP*/
