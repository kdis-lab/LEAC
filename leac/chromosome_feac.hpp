/*! \file chromosome_feac.hpp
 *
 * \brief chromosome FEAC \cite Hruschka:etal:GAClusteringLabelKVar:EAC:2006 and \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CHROMOSOME_FEAC_HPP
#define CHROMOSOME_FEAC_HPP

#include <vector>
#include "matrix_withrownull.hpp"
#include "vector_utils.hpp"
#include "chromosome_fixedlength.hpp"
#include "verbose_global.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {

typedef enum { FEAC_OPERATOR_UNKNOWN=0, 
	       FEAC_OPERATOR_MO1=1, 
	       FEAC_OPERATOR_MO2=2} 
EnumFeacOperatorApplied;


/*! \class ChromosomeFEAC
  \brief Chromosome encoded integer string and centroids 
  \details  A partition is encoded as an integer string (genotype) of N positions. Each string position corresponds to a dataset object, i.e., the i-th position represents the i-thobject of the dataset. Thus, each string component has avalue over the possible cluster labels {1,2,3,...,k} \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006. 
*/ 
template <class T_CLUSTERIDX,
	  class T_METRIC, //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  class T_FEATURE,
	  class T_FEATURE_SUM, //To apply kmeansfeac without template parameters
	  class T_INSTANCES_CLUSTER_K
	  >
class ChromosomeFEAC: public gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>
{   
public:
  ChromosomeFEAC()
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>() 
    , _matrixwrownull_centroids()
    , _vectorT_numInstancesClusterK()
    , _vectorT_partialFcC()
    , _enum_feacOperatorApplied(FEAC_OPERATOR_UNKNOWN)
    , _t_lastObjetiveFunc(-std::numeric_limits<T_METRIC>::max())
  {
    this->setObjetiveFunc(-std::numeric_limits<T_METRIC>::max());
  }

  ChromosomeFEAC
  (const uintidx aiuintidx_numClusterK,
   const uintidx aiuintidx_numDimensions
   )
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>() 
    , _matrixwrownull_centroids(aiuintidx_numClusterK, (uintidx) 2 * aiuintidx_numClusterK ,aiuintidx_numDimensions)
    , _vectorT_numInstancesClusterK(aiuintidx_numClusterK,0)
    , _vectorT_partialFcC(aiuintidx_numClusterK)
    , _enum_feacOperatorApplied(FEAC_OPERATOR_UNKNOWN)
    , _t_lastObjetiveFunc(-std::numeric_limits<T_METRIC>::max())
  {}

  //copy constructor
  ChromosomeFEAC
  (const ChromosomeFEAC
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
   T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K
    >                        &ichromfeac_b
   )
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>(ichromfeac_b)
    , _matrixwrownull_centroids(ichromfeac_b._matrixwrownull_centroids)
    , _vectorT_numInstancesClusterK(ichromfeac_b._vectorT_numInstancesClusterK)
    , _vectorT_partialFcC(ichromfeac_b._vectorT_partialFcC)
    , _enum_feacOperatorApplied(ichromfeac_b._enum_feacOperatorApplied)
    , _t_lastObjetiveFunc(ichromfeac_b._t_lastObjetiveFunc)
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "Copy:ChromosomeFEAC::ChromosomeFEAC";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ')'
	        << "\n\t input  ChromosomeFEAC<>: &ichromfeac_b[" 
		<<  &ichromfeac_b << "]\n"
		<< "\t)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< " OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelFunc;
      lostrstream_labelFunc <<  geverbosepc_labelstep << ':' << lpc_labelFunc;
      ChromosomeFEAC<T_CLUSTERIDX,T_METRIC,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>::print
	(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  }

  //move constructor
  ChromosomeFEAC
  (ChromosomeFEAC
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
   T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K
    >                       &&ichromfeac_b)
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>(ichromfeac_b)
    , _matrixwrownull_centroids(ichromfeac_b._matrixwrownull_centroids)
    , _vectorT_numInstancesClusterK(ichromfeac_b._vectorT_numInstancesClusterK)
    , _vectorT_partialFcC(ichromfeac_b._vectorT_partialFcC)
    , _enum_feacOperatorApplied(ichromfeac_b._enum_feacOperatorApplied)
    , _t_lastObjetiveFunc(ichromfeac_b._t_lastObjetiveFunc)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "Move:ChromosomeFEAC::ChromosomeFEAC";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ')'
	        << "\ninput  ChromosomeFEAC<>: &&ichromfeac_b[" 
		<<  &ichromfeac_b << "]\n)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelFunc;
      lostrstream_labelFunc <<  geverbosepc_labelstep << ':' << lpc_labelFunc;
      ChromosomeFEAC<T_CLUSTERIDX,T_METRIC,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>::print
      (std::cout,lostrstream_labelFunc.str().c_str(),',',';');
     
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  }

  virtual ~ChromosomeFEAC() {}

  ChromosomeFEAC
  <T_CLUSTERIDX,
   T_METRIC,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>& 
  operator=
  (const ChromosomeFEAC
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
    T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K
    >                        &ichromfeac_b)
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeFEAC::operator=Copy";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelFunc;
      lostrstream_labelFunc <<  geverbosepc_labelstep << ':' << lpc_labelFunc;
      ichromfeac_b.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
      std::cout	<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    if ( this != &ichromfeac_b ) {
      gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(ichromfeac_b);
      _matrixwrownull_centroids = ichromfeac_b._matrixwrownull_centroids;
      _vectorT_numInstancesClusterK = ichromfeac_b._vectorT_numInstancesClusterK;
      _vectorT_partialFcC = ichromfeac_b._vectorT_partialFcC;
      _enum_feacOperatorApplied = ichromfeac_b._enum_feacOperatorApplied;
      _t_lastObjetiveFunc = ichromfeac_b._t_lastObjetiveFunc;
    }
    
#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeFEAC<T_CLUSTERIDX,T_METRIC,T_FEATURE,T_FEATURE_SUM,T_INSTANCES_CLUSTER_K>::print
	(std::cout,lpc_labelFunc,',',';');  
     
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return *this;

  }

  ChromosomeFEAC
  <T_CLUSTERIDX,
   T_METRIC,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>& 
  operator=
  (const ChromosomeFEAC
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
    T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K
    >                    &&ichromfeac_b)
  {
    if ( this != &ichromfeac_b ) {
      gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(ichromfeac_b);
      _matrixwrownull_centroids = ichromfeac_b._matrixwrownull_centroids;
      _vectorT_numInstancesClusterK = ichromfeac_b._vectorT_numInstancesClusterK;
      _vectorT_partialFcC = ichromfeac_b._vectorT_partialFcC;
      _enum_feacOperatorApplied = ichromfeac_b._enum_feacOperatorApplied;
      _t_lastObjetiveFunc = ichromfeac_b._t_lastObjetiveFunc;
    }
    return *this;
  }

  inline
  mat::MatrixWithRowNull<T_FEATURE>& 
  getCentroids()
  {
    return _matrixwrownull_centroids;
  }

  inline
  std::vector<T_INSTANCES_CLUSTER_K>& 
  getNumInstancesClusterK()
  {
    return _vectorT_numInstancesClusterK;
  }

   inline
  T_CLUSTERIDX
  getNumClusterK() 
  {
    return (T_CLUSTERIDX) _matrixwrownull_centroids.getNumRows();
  }
  
  inline void setPartialFcC(const std::vector<T_METRIC>& aivectorrt_partialFcC)
  {
    interfacesse::copy
      (_vectorT_partialFcC.data(),
       aivectorrt_partialFcC.data(),
       (uintidx) aivectorrt_partialFcC.size() 
       );
  }

   inline
  std::vector<T_METRIC>& getPartialFcC()
  {
    return this->_vectorT_partialFcC;
  }

  inline 
  void setAppliedOperator
  (const EnumFeacOperatorApplied aienum_feacOperatorApplied) 
  {
    _enum_feacOperatorApplied = aienum_feacOperatorApplied;
  }

  inline 
  const EnumFeacOperatorApplied getAppliedOperator() const 
  {
    return _enum_feacOperatorApplied;
  }
  
  void saveLastObjetiveFunc() 
  {
    _t_lastObjetiveFunc = this->_t_objetiveFunc;
  }

  inline const T_METRIC getLastObjetiveFunc() const 
  {
    return _t_lastObjetiveFunc;
  }

  virtual void  print
  (std::ostream &os=std::cout,
   const char* aipc_label   = "",
   const char aic_delimCoef =',',
   const char aic_delimRow  =';'
   ) const 
  {
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << ":MEMBERCLUSTER:"
				<< aipc_label
				<< ":_t_lastObjetiveFunc:" <<  _t_lastObjetiveFunc
				<< ":_enum_feacOperatorApplied:" <<_enum_feacOperatorApplied
				<< ":ChromosomeFEAC:";
    gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::print
      (os,
       lostrstream_labelMemberShip.str().c_str(),
       aic_delimCoef,
       aic_delimRow
       );
    os << '\n';

    std::ostringstream lostrstream_labelCentroids;
#if defined(__VERBOSE_YES)    
    lostrstream_labelCentroids << "<CENTROIDSCLUSTER:"
			       << aipc_label
			       << ":ChromosomeFEAC[" << geverboseui_idproc << ':' <<  this << ']';
#else
    lostrstream_labelCentroids << "<CENTROIDSCLUSTER:"
			       << aipc_label
			       << ":ChromosomeFEAC[" <<  this << ']';
#endif 
     
    _matrixwrownull_centroids.print
      (os,
       lostrstream_labelCentroids.str().c_str(),
       aic_delimCoef,
       aic_delimRow
       );
    os << '\n';
    
    std::ostringstream lostrstream_labelNumInstancesK;
#if defined(__VERBOSE_YES)    
    lostrstream_labelNumInstancesK << "<NUMINSTANCESCLUSTERK:"
				   << aipc_label
				   << ":ChromosomeFEAC[" << geverboseui_idproc << ':' <<  this << ']';
#else
    lostrstream_labelNumInstancesK << "<NUMINSTANCESCLUSTERK:"
				   << aipc_label
				   << ":ChromosomeFEAC["  <<  this << ']';
#endif 
     
    inout::containerprint
      (_vectorT_numInstancesClusterK.begin(),
       _vectorT_numInstancesClusterK.end(),
       os,
       lostrstream_labelNumInstancesK.str().c_str(),
       aic_delimCoef
       );
    os << '\n';
    std::ostringstream lostrstream_labelPartialFcC;
#if defined(__VERBOSE_YES)    
    lostrstream_labelPartialFcC << "<PARTIALFcC:"
				<< aipc_label
				<< ":ChromosomeFEAC[" << geverboseui_idproc << ':' <<  this << ']';
#else
    lostrstream_labelPartialFcC << "<PARTIALFcC:"
				<< aipc_label
				<< ":ChromosomeFEAC[" <<  this << ']';
#endif
    
    inout::containerprint
      (_vectorT_partialFcC.begin(),
       _vectorT_partialFcC.end(),
       os,
       lostrstream_labelPartialFcC.str().c_str(),
       aic_delimCoef
       );
    os << std::endl;
  }

protected:

  mat::MatrixWithRowNull<T_FEATURE> _matrixwrownull_centroids;
  std::vector
   <T_INSTANCES_CLUSTER_K>       _vectorT_numInstancesClusterK;
  std::vector<T_METRIC>          _vectorT_partialFcC;
  EnumFeacOperatorApplied        _enum_feacOperatorApplied;
  T_METRIC                       _t_lastObjetiveFunc;

}; /*ChromosomeFEAC*/

} /*END namespace gaencode*/

#endif  /*CHROMOSOME_FEAC_HPP*/
