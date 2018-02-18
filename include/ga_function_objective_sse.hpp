/*! \file ga_function_objective_sse.hpp
 *
 * \brief function objective SSE
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __GA_FUNCTION_OBJECTIVE_SSE_HPP
#define __GA_FUNCTION_OBJECTIVE_SSE_HPP

#include <utility> /*std::pair*/
#include "ga_function_objective.hpp"
#include "matrix.hpp"
#include "unsupervised_measures.hpp"
#include "dist.hpp"

/*! \namespace gafuncobj
  \brief Definition of objective function
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace gafuncobj {
  
/*! \class GAFuncObjSSE
  \brief Definition of the class GAFuncObjSSE. \cite Chang:etal:GAclustering:GAGR:2009
*/
  
template < typename T_FEATURE,//DATATYPE T_GENE OF CHROMOSOME IS EQ T_CENTROIDS AND T_INSTANCES
           typename T_CLUSTERIDX,
	   typename T_METRIC,
	   typename INPUT_ITERATOR
 	   >
class GAFuncObjSSE: 
  public GAFunctionObjective<T_FEATURE,T_METRIC> 
{
public:
  GAFuncObjSSE
  (const T_CLUSTERIDX                   aicidx_numClustersK,
   const INPUT_ITERATOR                 aiiterator_instfirst,
   const INPUT_ITERATOR                 aiiterator_instlast,
   const dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist
   )
    :  GAFunctionObjective<T_FEATURE,T_METRIC>()
    , _cidx_numClustersK(aicidx_numClustersK)
    , _iterator_instfirst(aiiterator_instfirst)
    , _iterator_instlast(aiiterator_instlast)
    , _func2p_dist(&aifunc2p_dist)
  {}

  ~GAFuncObjSSE() { }

   //copy constructor
  GAFuncObjSSE(const GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>& aifuncobjgagr_b)
    : _cidx_numClustersK(aifuncobjgagr_b._cidx_numClustersK)
    , _iterator_instfirst(aifuncobjgagr_b._iterator_instfirst)
    , _iterator_instlast(aifuncobjgagr_b._iterator_instlast)
    , _func2p_dist(aifuncobjgagr_b._func2p_dist)
  {
  }
    

  //move constructor
  GAFuncObjSSE(GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>&& aifuncobjgagr_b)
    : _cidx_numClustersK(aifuncobjgagr_b._cidx_numClustersK)
    , _iterator_instfirst(aifuncobjgagr_b._iterator_instfirst)
    , _iterator_instlast(aifuncobjgagr_b._iterator_instlast)
    , _func2p_dist(aifuncobjgagr_b._func2p_dist)
  {
    aifuncobjgagr_b._func2p_dist = NULL;
  }

  GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>&
  operator=(const GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>& aifuncobjgagr_b)
  {
    if ( this != &aifuncobjgagr_b ) {
      _cidx_numClustersK = aifuncobjgagr_b._cidx_numClustersK;
      _iterator_instfirst = aifuncobjgagr_b._iterator_instfirst;
      _iterator_instlast = aifuncobjgagr_b._iterator_instlast;
     _func2p_dist = aifuncobjgagr_b._func2p_dist;
    }

    return *this;
  }

  GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>&
  operator=(GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>&& aifuncobjgagr_b)
  {
    if ( this != &aifuncobjgagr_b ) {
      _cidx_numClustersK = aifuncobjgagr_b._cidx_numClustersK;
      _iterator_instfirst = aifuncobjgagr_b._iterator_instfirst;
      _iterator_instlast = aifuncobjgagr_b._iterator_instlast;
     _func2p_dist = aifuncobjgagr_b._func2p_dist;

     aifuncobjgagr_b._cidx_numClustersK = 0;
     aifuncobjgagr_b._iterator_instfirst = NULL;
     aifuncobjgagr_b._iterator_instlast = NULL;
     aifuncobjgagr_b._func2p_dist = NULL;
    }

    return (*this);
  }
  
  /*getObjetiveFunc:
    based on SSE
  */ 
  std::pair<T_METRIC,bool> getObjetiveFunc(T_FEATURE* aistrt_string)
  {

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax  ) {
      std::cout << "GAFuncObjSSE::getObjetiveFunc:  IN"
		<< '(' << geiinparam_verbose << ")\n"
		<< "(input T_FEATURE*: aistrt_string[" << aistrt_string << ']'
		<< "\n)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    mat::MatrixRow<T_FEATURE>   
      lmatrixT_centroidsChrom
      ((uintidx) this->_cidx_numClustersK, 
       data::Instance<T_FEATURE>::getNumDimensions(),
       aistrt_string
       );
  
    std::pair<T_METRIC,bool> 
      lpair_SSE =
      um::SSE 
      (lmatrixT_centroidsChrom,
       _iterator_instfirst,
       _iterator_instlast, 
       *(this->_func2p_dist)
       );

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "GAFuncObjSSE::getObjetiveFunc: OUT"
		<< '(' << geiinparam_verbose << ")\n"
	        << "output lpair_SSE.first = " << lpair_SSE.first 
		<< "\nlpair_SSE.second  = "  << (lpair_SSE.second == 0)
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
    return std::make_pair(lpair_SSE.first,lpair_SSE.second == 0);

  }

  inline T_METRIC getFitness(T_METRIC airt_objetiveValue)
  {
    return 1.0 / airt_objetiveValue;
  }

protected:

  const T_CLUSTERIDX                      _cidx_numClustersK;
  const INPUT_ITERATOR                    _iterator_instfirst;
  const INPUT_ITERATOR                    _iterator_instlast;
  const dist::Dist<T_METRIC,T_FEATURE>    *_func2p_dist;
  
};

template < typename T_FEATURE,
           typename T_CLUSTERIDX,
	   typename T_METRIC,
	   typename INPUT_ITERATOR
	   >
GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>
makeGAFuncObjSSE
(const T_CLUSTERIDX                   aicidx_numClusters,
 const INPUT_ITERATOR                 aiiterator_instfirst,
 const INPUT_ITERATOR                 aiiterator_instlast,
 const dist::Dist<T_METRIC,T_FEATURE>  &aifunc2p_dist
)
{

  return GAFuncObjSSE<T_FEATURE,T_CLUSTERIDX,T_METRIC,INPUT_ITERATOR>(aicidx_numClusters,aiiterator_instfirst,aiiterator_instlast,aifunc2p_dist);
  
}

  
}  /*END namespace gafuncobj*/

#endif /*__GA_FUNCTION_OBJECTIVE_SSE_HPP*/


