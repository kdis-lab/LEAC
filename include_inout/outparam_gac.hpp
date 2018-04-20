/*! \file outparam_gac.hpp
 *
 * \brief out parameters EA-clustering
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __OUT_PARAM_EA_HPP__
#define __OUT_PARAM_EA_HPP__


#include "outparam_clustering.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
/*! \class OutParamGAC
  \brief Output Output parameters of the evolutionary algorithm
*/
template < typename T_METRIC,
	   typename T_CLUSTERIDX
	  >
class OutParamGAC:
  public OutParamClustering<T_METRIC,T_CLUSTERIDX> {
public:
  OutParamGAC(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc):
  OutParamClustering<T_METRIC,T_CLUSTERIDX>
    ::OutParamClustering(aienum_usedObjectiveFunc)   
  {
    this->initialize(-1);
  }

  virtual ~OutParamGAC() {}

  void initialize
  (int             aii_numRunAlgorithm)
  {
    OutParamClustering<T_METRIC,T_CLUSTERIDX>
    ::initialize(aii_numRunAlgorithm);
   
    this->t_fitness = OUTPARAMCLUSTERING_FITNESS_NaN;
    this->t_numTotalGenerations = OUTPARAMCLUSTERING_INT_NaN; 
  }

  inline T_METRIC getFitness()
  {
    return this->t_fitness;
  }

  inline void setFitness(T_METRIC aiT_fitness) 
  {
    this->t_fitness = aiT_fitness;
  }

  inline void incNumGenerations(COMMON_IDOMAIN aiiT_numTotalGenerations) 
  {
    this->t_numTotalGenerations += aiiT_numTotalGenerations;
  }

  inline void setNumTotalGenerations(COMMON_IDOMAIN aiiT_numTotalGenerations) 
  {
    this->t_numTotalGenerations = aiiT_numTotalGenerations;
  }

  inline COMMON_IDOMAIN getNumTotalGenerations() 
  {
    return this->t_numTotalGenerations;
  }

  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    OutParamClustering<T_METRIC,T_CLUSTERIDX>::print(aipf_outFile);
   
    aipf_outFile << aic_separator << "_fitness" 
		 << aic_separator << this->t_fitness;
    aipf_outFile << aic_separator << "_number total generations" 
		 << aic_separator << this->t_numTotalGenerations;
  }

protected:  

  T_METRIC        t_fitness;             
  COMMON_IDOMAIN  t_numTotalGenerations; 

}; /*OutParamGAC*/

} /*END namespace inout
   */

#endif /*__OUT_PARAM_EA_HPP__*/
