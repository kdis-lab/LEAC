/*! \file outparam_dbscan.hpp
 *
 * \brief outparam DBSCAN
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __OUT_DBSCAN_HPP__
#define __OUT_DBSCAN_HPP__


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


} /*END namespace inout
   */

#endif /*__OUT_DBSCAN_HPP__*/
