/*! \file partition_fuzzy.hpp
 *
 * \brief Partition based on fuzzy belonged
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MEMBERSHIP_PARTITION_FUZZYPMATRIX_HPP
#define MEMBERSHIP_PARTITION_FUZZYPMATRIX_HPP

#include "partition.hpp"
#include "matrix.hpp"
#include "common.hpp"

/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace partition {


/*! \class PartitionFuzzy
  \brief Partition of instances with fuzzy matrix \cite Bezdek:ClusterAnalysis:FCM:1974 \cite Bezdek:etal:ClusterAnalysis:FCM:1984
*/
template< class T_U,
          class T_CLUSTERIDX
	  >
class PartitionFuzzy:
 public 
Partition<T_CLUSTERIDX>
{
public:

  PartitionFuzzy
  (const mat::MatrixRow<T_U>   &aimatrixt_crisp
   )
    : Partition<T_CLUSTERIDX>()
    , _ptmatrix_crisp(&aimatrixt_crisp)   
  { }
  
  ~PartitionFuzzy()
  { }

  void begin()
  {
    _ui_iteratornext = 0;
  }

  inline const T_CLUSTERIDX next()
  {
    return this->getClusterIdx(_ui_iteratornext++);	
  }

  inline bool end() const
  {
    return ( _ui_iteratornext < _ptmatrix_crisp->getNumColumns());
  }
  
  const T_CLUSTERIDX getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {
    uintidx luintidx_idxUjkMax;
    T_U    lT_ujkMax;

    luintidx_idxUjkMax = 0;
    lT_ujkMax = (*_ptmatrix_crisp)(0, aiuintidx_instanceIdx);
    for (uintidx li_i = 1; li_i < _ptmatrix_crisp->getNumRows(); li_i++) {
      if (  (*_ptmatrix_crisp)(li_i, aiuintidx_instanceIdx) > lT_ujkMax) {
	luintidx_idxUjkMax = li_i;
	lT_ujkMax    = (*_ptmatrix_crisp)(li_i, aiuintidx_instanceIdx);
      }
    }
    
    return T_CLUSTERIDX(luintidx_idxUjkMax);
  }

  inline const uintidx getNumInstances() const
  {
    return _ptmatrix_crisp->getNumColumns();
  }
  
  const T_CLUSTERIDX getNumCluster() const
  {
    return T_CLUSTERIDX( _ptmatrix_crisp->getNumRows() );  
  }

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label     = "",
   const char  aic_delimCoef  = ','
   )
  {
    uintidx lui_numInstances =
      (_ptmatrix_crisp->getNumColumns() > 0)?
      (uintidx) _ptmatrix_crisp->getNumColumns()-1:0;

    os << aipc_label << aic_delimCoef << "_outK" << aic_delimCoef << this->getNumCluster();
    
#if defined(__VERBOSE_YES)
    os << aic_delimCoef
       << "PartitionFuzzy[" << geverboseui_idproc << ':' << this << ']';
#else
    os << aic_delimCoef
       << "PartitionFuzzy[" << this << ']';
#endif

    os  << aic_delimCoef << "length" << aic_delimCoef << _ptmatrix_crisp->getNumColumns() << '>';
    for (uintidx  lui_i = 0; lui_i < lui_numInstances; lui_i++) {
      os << this->getClusterIdx(lui_i) << aic_delimCoef;
    }
    os << this->getClusterIdx(lui_numInstances);
    os << std::endl;
  }
  
protected:
  
  const mat::MatrixRow<T_U> *_ptmatrix_crisp;
   uintidx _ui_iteratornext;
  
};
  
  
} /* END namespace partition*/

#endif  /* MEMBERSHIP_PARTITION_FUZZYPMATRIX_HPP */
