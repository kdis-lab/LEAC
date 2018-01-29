/*! \file partition.hpp
 *
 * \brief partition
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef MEMBERSHIP_PARTITION_HPP
#define MEMBERSHIP_PARTITION_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


/*! \namespace partition
  \brief Gets the indexes of the group that each instance belongs to for different schemas representing a partition
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace partition {

/*! \class Partition
  \brief Definition of the abstract class Partition, to build a cluster of a dataset 
*/
template <class T_CLUSTERIDX>
class Partition {
public:

  Partition() { }
  
  Partition
  (const Partition<T_CLUSTERIDX>& aimembclassdisjset_b) { }

  Partition
  (Partition<T_CLUSTERIDX>&& aimembclassdisjset_b) { }
  
  ~Partition() {}
  
  virtual const T_CLUSTERIDX 
  getClusterIdx(uintidx aiuintidx_instanceIdx) const = 0;

  virtual void begin() = 0;
  virtual const T_CLUSTERIDX next() = 0;
  virtual bool end() const = 0;

  virtual const uintidx getNumInstances() const = 0;
  virtual const T_CLUSTERIDX getNumCluster() const = 0;

  virtual void print
  (std::ostream &os=std::cout,
   const char* aipc_label      = "",
   const char  aic_delimCoef   = ','
   ) = 0; 

protected:

}; /*Partition*/ 


} /* END namespace partition*/

#endif  /* MEMBERSHIP_PARTITION_HPP */



