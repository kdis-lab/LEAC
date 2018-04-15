/*! \file inparam_dbscan.hpp
 *
 * \brief Definition of DBSCAN program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __IN_PARAM_DBSCAN_HPP__
#define __IN_PARAM_DBSCAN_HPP__

#include "inparam_clustering.hpp"
#include "inparam_readinst.hpp"

#define __INPARAM_DBSCAN__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/


namespace  inout {
  
/*! \class InParamDBSCAN
  \brief Input parameter for DBSCAN algorithm \cite Ester:Kriegel:Sander:Xu:Clustering:DBSCAN:1996
*/
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	  > 
class InParamDBSCAN
  : public InParamClustering
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamDBSCAN
  (std::string        ais_algorithmoName,
   std::string        ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   ):
    InParamClustering
    (ais_algorithmoName,
     ais_algorithmoAuthor,
     aiato_algTypeOut,
     aii_opNorm)
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()    
  {}

  ~InParamDBSCAN() {}
  
  inline void setEps(double aid_eps) 
  {
    this->_d_eps = aid_eps;
  }

  inline const double getEps() const
  {
    return this->_d_eps;
  }

  inline void setMinPts(uint32_t aiui32t_minPts) 
  {
    this->_ui32t_minPts = aiui32t_minPts;
  }

  inline const uint32_t getMinPts() const
  {
    return this->_ui32t_minPts;
  }


  /*SpatialIndex
   */
  inline void setRTreeBufferCapacity(uint32_t airtreeui32t_bufferCapacity)
  {
    _rtreeui32t_bufferCapacity = airtreeui32t_bufferCapacity;
  }
  
  inline uint32_t getRTreeBufferCapacity()
  {
    return _rtreeui32t_bufferCapacity;
  }

  inline void setRTreeFillFactor(double airtreed_fillFactor)
  {
    _rtreed_fillFactor = airtreed_fillFactor;
  }

  inline double getRTreeFillFactor()
  {                
    return _rtreed_fillFactor;
  }

  inline void setRTreeCapacity(uint32_t airtreeui32t_indexCapacity)
  {  
    _rtreeui32t_indexCapacity = airtreeui32t_indexCapacity;
  }

  inline uint32_t getRTreeCapacity()
  {
    return _rtreeui32t_indexCapacity;
  }

  inline void setRTreeleafCapacity(uint32_t airtreeui32t_leafCapacity)
  {
    _rtreeui32t_leafCapacity = airtreeui32t_leafCapacity;
  }
  
  
  inline uint32_t getRTreeleafCapacity() 
  {
    return _rtreeui32t_leafCapacity;
  }


  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamClustering::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_eps" 
		 << aic_separator <<  this->_d_eps
		 << aic_separator << "_minpts" 
		 << aic_separator <<  this->_ui32t_minPts
		 << aic_separator << "_rtreeui32t_bufferCapacity"
		 << aic_separator << _rtreeui32t_bufferCapacity
		 << aic_separator << "_rtreed_fillFactor"
		 << aic_separator << _rtreed_fillFactor
		 << aic_separator << "_rtreeui32t_indexCapacity"
		 << aic_separator << _rtreeui32t_indexCapacity
		 << aic_separator << "_rtreeui32t_leafCapacity"
		 << aic_separator << _rtreeui32t_leafCapacity;
    //louintidx_numFields += (uintidx) 6; 
    //return louintidx_numFields;
  }
protected:
  
  double                _d_eps;
  uint32_t              _ui32t_minPts;

  /*SpatialIndex
   */
  uint32_t              _rtreeui32t_bufferCapacity;
  double                _rtreed_fillFactor;
  uint32_t              _rtreeui32t_indexCapacity;
  uint32_t              _rtreeui32t_leafCapacity;
  
};

  
} /*END namespace inout
   */

#endif /*__IN_PARAM_DBSCAN_HPP__*/
