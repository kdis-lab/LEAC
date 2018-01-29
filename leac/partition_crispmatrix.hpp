/*! \file partition_crispmatrix.hpp
 *
 * \brief Partition based on crispmatrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MEMBERSHIP_PARTITION_BITCRISPMATRIX_HPP
#define MEMBERSHIP_PARTITION_BITCRISPMATRIX_HPP

#include "partition.hpp"
#include "crisp_matrix.hpp"
#include "common.hpp"

/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace partition {


/*! \class PartitionCrispMatrix
  \brief Partition of instances with bit crisp matrix
*/
template< class T_BITSIZE,
          class T_CLUSTERIDX
	  >
class PartitionCrispMatrix: 
  public Partition<T_CLUSTERIDX>
{
public:

  PartitionCrispMatrix
  (const mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aibitcrisp_matrix)
    : Partition<T_CLUSTERIDX>() 
    , _ptcrisp_matrix(&aibitcrisp_matrix)
    , _ui_iteratornext(0)
  { }
  
 
  ~PartitionCrispMatrix()
  { }

  void begin()
  {
    _ui_iteratornext = 0;
  }

  inline const T_CLUSTERIDX next()
  {
    return this->_ptcrisp_matrix->getMember(_ui_iteratornext++);	
  }

  inline bool end() const
  {
    return ( _ui_iteratornext < _ptcrisp_matrix->getNumColumns());
  }
  
  virtual const T_CLUSTERIDX 
  getClusterIdx(uintidx aiuintidx_instanceIdx) const
  {
    return this->_ptcrisp_matrix->getMember(aiuintidx_instanceIdx);	
  }

  inline const uintidx getNumInstances() const
  {
    return _ptcrisp_matrix->getNumColumns();
  }

  const T_CLUSTERIDX getNumCluster() const
  {
    return T_CLUSTERIDX(this->_ptcrisp_matrix->getNumRows());  
  }

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label     = "",
   const char  aic_delimCoef  = ','
   )
  {
    os << aipc_label << aic_delimCoef << "_outK" << aic_delimCoef << this->getNumCluster();
    
#if defined(__VERBOSE_YES)
    os <<  aic_delimCoef
       << "PartitionCrispMatrix[" << geverboseui_idproc << ':' << this << ']';
#else
    os <<  aic_delimCoef
       << "PartitionCrispMatrix[" << this << ']';
#endif

    os  << aic_delimCoef << "length" << aic_delimCoef << _ptcrisp_matrix->getNumColumns() << '>';
    this->begin();
    os << this->next();
    while ( this->end() ) {
      os << aic_delimCoef << this->next();
    }
    os << std::endl;
  }

protected:

  const mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> *_ptcrisp_matrix;
  uintidx _ui_iteratornext;

}; /*PartitionCrispMatrix*/

  
} /* END namespace partition*/


#endif  /* MEMBERSHIP_PARTITION_BITCRISPMATRIX_HPP */
