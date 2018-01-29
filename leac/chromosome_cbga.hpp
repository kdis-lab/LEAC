/*! \file chromosome_cbga.hpp
 *
 * \brief chromosome CBGA
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CHROMOSOME_CBGA_HPP__
#define __CHROMOSOME_CBGA_HPP__

#include <vector>       // std::vector
#include <algorithm>    // std::equal
#include <stdexcept>
#include "partition_linked_stats.hpp"
#include "nearestcentroids_operator.hpp"
#include "clustering_operator_centroids.hpp"
#include "instance.hpp"

#include "verbose_global.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace gaencode {
  

typedef  enum { OPT_NONE=0, OPT_CB=1, OPT_PA=2, OPT_BOTH=3 } OptimalityCBGA;

#define OPTIMALITY_CBGA_STRING {"OPT_NONE", "OPT_CB", "OPT_PA", "OPT_BOTH", (char *) NULL }

/*! \class ChromosomeCBGA
  \brief A set of centroids stands for a codebook of the application and the partitions \cite Franti:etal:GAclustering:gafranti:1997 
*/
template <typename T_FEATURE,  //TYPE GENE
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY,
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_FEATURE_SUM,
	  typename T_REAL
	  >
class ChromosomeCBGA    
{   
public:
  ChromosomeCBGA
  (const uintidx aiui_numClusterK,
   const uintidx aiui_numClusterKSpace,
   const uintidx aiui_numDimensionsInstances,
   const uintidx aiui_numInstances 
   ) : _matrixresizerow_codeBook
       (new mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>
	(aiui_numClusterK,aiui_numDimensionsInstances,aiui_numClusterKSpace))
     , _partlinkstats_part
       (new ds::PartitionLinkedStats
	<T_FEATURE,
	T_CLUSTERIDX,
	T_INSTANCE_FREQUENCY,
	T_INSTANCES_CLUSTER_K,
	T_FEATURE_SUM
	>
	(aiui_numClusterK,
	 aiui_numInstances,
	 aiui_numDimensionsInstances,
	 aiui_numClusterKSpace)
	)  
     , _opf_optimalityCBGA(OPT_NONE)
     , _b_stringInvalid(false)
     , _t_objetiveFunc(std::numeric_limits<T_REAL>::max())
  { }
  
  //copy constructor 
  ChromosomeCBGA
  (const ChromosomeCBGA
   <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,    
   T_FEATURE_SUM,
   T_REAL>                &aichromcbga_src
   )
    : _matrixresizerow_codeBook
      (new mat::MatrixResizableRow
       <T_FEATURE,T_INSTANCES_CLUSTER_K>(*aichromcbga_src._matrixresizerow_codeBook)
       )
    , _partlinkstats_part
       (new ds::PartitionLinkedStats
	<T_FEATURE,
	T_CLUSTERIDX,
	T_INSTANCE_FREQUENCY,
	T_INSTANCES_CLUSTER_K,
	T_FEATURE_SUM
	>
	(*aichromcbga_src._partlinkstats_part)
	)
    , _opf_optimalityCBGA(aichromcbga_src._opf_optimalityCBGA)
    , _b_stringInvalid(aichromcbga_src._b_stringInvalid)
    , _t_objetiveFunc(aichromcbga_src._t_objetiveFunc)
  {
  }
  

  
  //move constructor 
  ChromosomeCBGA
  (ChromosomeCBGA
   <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,    
   T_FEATURE_SUM,
   T_REAL> &&aichromcbga_src)
    : _matrixresizerow_codeBook(aichromcbga_src._matrixresizerow_codeBook)
    , _partlinkstats_part(aichromcbga_src._partlinkstats_part)
    , _opf_optimalityCBGA(aichromcbga_src._opf_optimalityCBGA)
    , _b_stringInvalid(aichromcbga_src._b_stringInvalid)
    , _t_objetiveFunc(aichromcbga_src._t_objetiveFunc)
  {
    aichromcbga_src._matrixresizerow_codeBook = NULL;
    aichromcbga_src._partlinkstats_part = NULL;
  }
  
  
  virtual ~ChromosomeCBGA() 
  {
    if ( this->_matrixresizerow_codeBook != NULL ) 
      delete _matrixresizerow_codeBook;
    if ( this->_partlinkstats_part != NULL )
      delete _partlinkstats_part;
  }
         
  ChromosomeCBGA
  <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM,
   T_REAL>&
  operator=
  (const ChromosomeCBGA
   <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM,
   T_REAL>                       &aichromcbga_src
   )
  {
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ChromosomeCBGA::copy:  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\t(output ChromosomeCBGA: this[" << this << "]\n"
		<< "\t input  ChromosomeCBGA: aichromcbga_src[" 
		<<  &aichromcbga_src << "]\n"
		<< "\t)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

    if ( this != &aichromcbga_src) {
      
      *this->_matrixresizerow_codeBook = *aichromcbga_src._matrixresizerow_codeBook;
      *this->_partlinkstats_part =  *aichromcbga_src._partlinkstats_part;
     
      this->_opf_optimalityCBGA = aichromcbga_src._opf_optimalityCBGA;
      this->_b_stringInvalid = aichromcbga_src._b_stringInvalid;
      this->_t_objetiveFunc = aichromcbga_src._t_objetiveFunc;

    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "ChromosomeCBGA::copy: OUT"
		<< '(' << geiinparam_verbose << ")\n"
		<< "\toutput ChromosomeCBGA: this[" << this << "]\n";
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

  return *this;
  
  }

  
  ChromosomeCBGA
  <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM,
   T_REAL>& 
  operator=
  (ChromosomeCBGA
   <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,
   T_FEATURE_SUM,
   T_REAL>
   &&aichromcbga_src
   )
  {
    if ( this != &aichromcbga_src ) {
      
      delete this->_matrixresizerow_codeBook;
      delete this->_partlinkstats_part;

      this->_matrixresizerow_codeBook = aichromcbga_src._matrixresizerow_codeBook;
      this->_partlinkstats_part = aichromcbga_src._partlinkstats_part;
   
      this->_opf_optimalityCBGA = aichromcbga_src._opf_optimalityCBGA;
      this->_b_stringInvalid = aichromcbga_src._b_stringInvalid;
      this->_t_objetiveFunc = aichromcbga_src._t_objetiveFunc;

      aichromcbga_src._matrixresizerow_codeBook = NULL;
      aichromcbga_src._partlinkstats_part = NULL;
      aichromcbga_src._b_stringInvalid = true;
      aichromcbga_src._t_objetiveFunc = std::numeric_limits<T_REAL>::max();
      aichromcbga_src._opf_optimalityCBGA = false;
    }

    return *this;
  }
  
  inline 
  mat::MatrixResizableRow<T_FEATURE,T_INSTANCES_CLUSTER_K>&  getCodeBook()
  {
    return *(this->_matrixresizerow_codeBook);
  }
  
  inline 
  ds::PartitionLinkedStats
  <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,    
   T_FEATURE_SUM
   >&   getPartition()
  {
    return *(this->_partlinkstats_part);
  }
	
  inline OptimalityCBGA getOptimalityCBGA()
  { 
    return this->_opf_optimalityCBGA;
  }

  inline void setOptimalityCBGA(OptimalityCBGA aiooptimalitycbga_a)
  { 
    this->_opf_optimalityCBGA = aiooptimalitycbga_a;
  }

  inline 
  void optimalCodebook() /* Calculate Mean Centroids */ 
  {
    T_CLUSTERIDX  lcidx_numClusterNull;
   
    clusteringop::meanCentroids
      (lcidx_numClusterNull,
       *(this->_matrixresizerow_codeBook), 
       this->_partlinkstats_part->getSumInstancesClusterK(),
       this->_partlinkstats_part->getNumInstancesClusterK()
       );
    this->setOptimalityCBGA( OPT_CB );
    
  }
  
  void checkChangedCodeVectors
  (T_CLUSTERIDX             *aoarraycidx_newIndex,
   uintidx                  *aoui_countNewIndex,
   mat::MatrixResizableRow
   <T_FEATURE,
   T_INSTANCES_CLUSTER_K>   &aimatrixresize_oldCodeBook
   )
  {
    T_CLUSTERIDX lcidx_codeBookNumClusterK;

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeCBGA::checkChangedCodeVectors";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc 
		<< ":  IN(" << geiinparam_verbose << ")\n"
		<< "(output T_CLUSTERIDX*: aoarraycidx_newIndex[" 
		<< aoarraycidx_newIndex << "]\n"
		<< " output uintidx: aoui_countNewIndex = " 
		<<  *aoui_countNewIndex << '\n'
		<< " input  ChromosomeCBGA: this[" << this << "]\n";
      aimatrixresize_oldCodeBook.print
	(std::cout,
	 lpc_labelFunc,
	 ',',
	 ';'
	 );

      std::cout	<< ")" 
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    *aoui_countNewIndex = 0;  
    uintidx luintidx_arrayLength = 
      this->_matrixresizerow_codeBook->getNumColumns();
    lcidx_codeBookNumClusterK = (T_CLUSTERIDX) this->_matrixresizerow_codeBook->getNumRows();
    for ( T_CLUSTERIDX lcidx_j = 0; lcidx_j < lcidx_codeBookNumClusterK; lcidx_j++) {

      bool lb_isEqualArray=
	std::equal
	(this->_matrixresizerow_codeBook->getRow(lcidx_j),
	 this->_matrixresizerow_codeBook->getRow(lcidx_j) + luintidx_arrayLength, 
	 aimatrixresize_oldCodeBook.getRow(lcidx_j)
	 );	    
      
      if (!lb_isEqualArray) {
	
	aoarraycidx_newIndex[*aoui_countNewIndex] = lcidx_j;
	(*aoui_countNewIndex)++;
      }
    }
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelNewIndex;
      lostrstream_labelNewIndex << "<ARRAY_NEWINDEX:"
				<< lpc_labelFunc;
      
      inout::containerprint
	(aoarraycidx_newIndex,
	 aoarraycidx_newIndex+ *aoui_countNewIndex,
	 std::cout,
	 lostrstream_labelNewIndex.str().c_str(),
	 ','
	 );
      std::cout << std::endl;
   
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  } /*END checkChangedCodeVectors */

  inline void setValidString(bool aib_stringInvalid) 
  {
    this->_b_stringInvalid = aib_stringInvalid;
  }
  
  inline bool getValidString() 
  {
    return this->_b_stringInvalid;
  }

  inline void setObjetiveFunc(T_REAL ait_objetiveFunc) 
  {
    this->_t_objetiveFunc = ait_objetiveFunc;
  }

  inline T_REAL getObjetiveFunc() 
  {
    return this->_t_objetiveFunc;
  }

  inline void setNullObjetiveFunc() 
  {
    this->_t_objetiveFunc = -std::numeric_limits<T_REAL>::min();
  }

  inline  bool isNullObjetiveFunc() 
  {
    return (this->_t_objetiveFunc == -std::numeric_limits<T_REAL>::min());
  }

  virtual void  print
  (std::ostream &os=std::cout,
   const char* aipc_label   = "",
   const char aic_delimCoef =',',
   const char aic_delimRow  =';'
   ) const
  {
    
#if defined(__VERBOSE_YES)

    std::ostringstream lostrstream_labelCodeBook;
    lostrstream_labelCodeBook
      << "<CODEBOOK:CHROMOSOMECBGA:"
      << aipc_label
   
      << ":id[" << geverboseui_idproc << '-' << this << ']'
      << ":objetive function: " << this->_t_objetiveFunc;// << '\n';
#else
    std::ostringstream lostrstream_labelCodeBook;
    lostrstream_labelCodeBook << aipc_label;
#endif
    
    this->_matrixresizerow_codeBook->print
      (os,
       lostrstream_labelCodeBook.str().c_str(),
       aic_delimCoef,
       aic_delimRow
       );
    os << '\n';
 
    this->_partlinkstats_part->print
      (os,
       aipc_label,
       aic_delimCoef
       );
  }
  
protected:
	
  mat::MatrixResizableRow
  <T_FEATURE,
   T_INSTANCES_CLUSTER_K>           *_matrixresizerow_codeBook;
  ds::PartitionLinkedStats
  <T_FEATURE,
   T_CLUSTERIDX,
   T_INSTANCE_FREQUENCY,
   T_INSTANCES_CLUSTER_K,    
   T_FEATURE_SUM
   >                                *_partlinkstats_part;
  OptimalityCBGA                    _opf_optimalityCBGA;
  bool                              _b_stringInvalid;
  T_REAL                            _t_objetiveFunc;
	
}; /*END class ChromosomeCBGA*/

} /*END namespace gaencode*/


#endif /*__CHROMOSOME_CBGA_HPP__*/



