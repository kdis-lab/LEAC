/*! \file partition_medoids.hpp
 *
 * \brief Partition based on centroids medoids
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MEMBERSHIP_PARTITION_MEDOIDS_HPP
#define MEMBERSHIP_PARTITION_MEDOIDS_HPP

#include "partition.hpp"
#include "matrix_triangular.hpp"
#include "nearestinstance_operator.hpp"

/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace partition {

/*! \class PartitionMedoids
  \brief Partition of instances with closest medoids
\f[
x_i \in C_j \leftrightarrow  \| x_i - m_j \|  \begin{array}{c}min\\  k \end{array}
\| x_i - m_k \|,\; j=1,2,..k,
\f]
where \f$m_j \f$, represents the medoid of cluster \f$C_j\f$
*/
template < class T_INSTANCES_IDX,
           class T_CLUSTERIDX,
	   class T_DIST 
	   >
class  PartitionMedoids:
  public Partition<T_CLUSTERIDX>
{
public:
  
  PartitionMedoids
  (T_INSTANCES_IDX             *aiarrayidxinstT_medoids,
   T_CLUSTERIDX         aicidx_numClusterK,
   mat::MatrixTriang<T_DIST>   &aipmatrixtriagT_dissimilarity
   )
    :  Partition<T_CLUSTERIDX>() 
    ,  _arrayidxinst_medoids(aiarrayidxinstT_medoids)
    ,  _cidx_numClusterK(aicidx_numClusterK)
    ,  _ptmatrixtriagt_dissimilarity(&aipmatrixtriagT_dissimilarity) 
  {}

  ~PartitionMedoids() {}

  inline const uintidx getNumInstances() const
  {
    return _ptmatrixtriagt_dissimilarity->getNumRows();
  }
  
  void begin() 
  {
    _ui_iteratornextInstIdx = 0;
  }

  inline const T_CLUSTERIDX next()
  {
    T_CLUSTERIDX  lodxT_clustering;
    T_DIST               lT_distMinCentInst;

    lodxT_clustering = 
      nearest::medoidsNN
      <T_CLUSTERIDX,T_DIST>
      (lT_distMinCentInst,
       _ui_iteratornextInstIdx,
       _arrayidxinst_medoids,
       _cidx_numClusterK,
       *_ptmatrixtriagt_dissimilarity
       );
    
    ++_ui_iteratornextInstIdx;

    return lodxT_clustering;
    
  }

  inline bool end() const
  {
    return ( _ui_iteratornextInstIdx != _ptmatrixtriagt_dissimilarity->getNumRows() );
  }
  
 const T_CLUSTERIDX getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {     
    T_CLUSTERIDX  lodxT_clustering;
    T_DIST               lT_distMinCentInst;
    
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionMedoids::getClusterIdx IN"
		<< '(' << geiinparam_verbose << ')'
		<< "\n(\t input  uintidx aiuintidx_instanceIdx = " <<  aiuintidx_instanceIdx << " )"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/
 
    lodxT_clustering = 
      nearest::medoidsNN
      <T_CLUSTERIDX,T_DIST>
      (lT_distMinCentInst,
       aiuintidx_instanceIdx,
       _arrayidxinst_medoids,
       _cidx_numClusterK,
       *_ptmatrixtriagt_dissimilarity
       );

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionMedoids::getClusterIdx OUT"
		<< '(' << geiinparam_verbose << ')'
		<< "\noutput T_CLUSTERIDX lodxT_clustering = " << lodxT_clustering 
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return lodxT_clustering;
  }

  const T_CLUSTERIDX getNumCluster() const
  {
    return _cidx_numClusterK;
  }

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label     = "",
   const char  aic_delimCoef  = ','
   ) 
  {
    /* uintidx lui_numInstances =
      (_ptmatrixtriagt_dissimilarity->getNumRows() > 0)?
      (uintidx) _ptmatrixtriagt_dissimilarity->getNumRows()-1:0;
    */
    os << aipc_label << aic_delimCoef << "_outK" << aic_delimCoef << this->getNumCluster();
    
#if defined(__VERBOSE_YES)
    os << aic_delimCoef
       << "PartitionMedoids[" << geverboseui_idproc << ':' << this << ']';
#else
    os << aic_delimCoef
       << "PartitionMedoids[" << this << ']';
#endif

    os  << aic_delimCoef << "length" << aic_delimCoef << _ptmatrixtriagt_dissimilarity->getNumRows()
	<< '>';
    this->begin();
    os << this->next();
    while ( this->end() ) {
      os << aic_delimCoef <<  this->next(); // << std::endl; //this->getClusterIdx(lui_i); 
    }
    /* for (uintidx  lui_i = 0; lui_i < lui_numInstances; lui_i++) {
      os << this->getClusterIdx(lui_i) << aic_delimCoef;
    }
    os << this->getClusterIdx(lui_numInstances);
    */
    os << std::endl;
  }

protected:
  T_INSTANCES_IDX            *_arrayidxinst_medoids;
  T_CLUSTERIDX        _cidx_numClusterK;
  mat::MatrixTriang<T_DIST>  *_ptmatrixtriagt_dissimilarity;
  uintidx                   _ui_iteratornextInstIdx;
};

  
} /* END namespace partition*/

#endif  /* MEMBERSHIP_PARTITION_MEDOIDS_HPP */
