/*! \file inparam_hka.hpp
 *
 * \brief Definition of HKA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_HKA_HPP__
#define __IN_PARAM_HKA_HPP__

#include "inparam_gca.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamHKA 
  \brief Input parameter for HKA a hybrid algorithm for k-medoid clustering \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename T_REAL,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamHKA:
    public InParamGCA<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
{
public:
  InParamHKA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   ) :
    InParamGCA<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
    (ais_algorithmoName, ais_algorithmoAuthor, aiato_algTypeOut,aii_opNorm) {}
  ~InParamHKA() {}
  
  inline void setOrderTournament(uintidx aist_orderTournament) 
  {
    this->st_orderTournament = aist_orderTournament;
  }

  inline uintidx getOrderTournament() 
  {
    return this->st_orderTournament;
  }

  inline void setNearestNeighbors(uintidx aist_nearestNeighbors) 
  {
    this->st_nearestNeighbors = aist_nearestNeighbors;
  }

  inline uintidx getNearestNeighbors() 
  {
    return this->st_nearestNeighbors;
  }

  inline void setProbSearchHeuristic(T_REAL aiT_probSearchHeuristic) 
  {
    this->t_probSearchHeuristic = aiT_probSearchHeuristic;
  }

  inline T_REAL getProbSearchHeuristic() 
  {
    return this->t_probSearchHeuristic;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamGCA<T_CLUSTERIDX,T_REAL,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>
      ::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_order of tournament"
		 << aic_separator << this->st_orderTournament;
    aipf_outFile << aic_separator << "_nearest neighbors (p)" 
		 << aic_separator << st_nearestNeighbors;
    aipf_outFile << aic_separator << "_probability local search heuristic"
		 << aic_separator << this->t_probSearchHeuristic;
  }

protected:
  uintidx  st_orderTournament;
  uintidx  st_nearestNeighbors;
  T_REAL   t_probSearchHeuristic; 
};

} /* END namespace inout
   */

#endif /*__IN_PARAM_HKA_HPP__*/
