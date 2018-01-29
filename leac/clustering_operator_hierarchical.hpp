/*! \file clustering_operator_hierarchical.hpp
 *
 * \brief  Clustering operator hierarchical
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CLUSTERING_OPERATOR_HIERARCHICAL_HPP
#define __CLUSTERING_OPERATOR_HIERARCHICAL_HPP

#include <iostream>
#include <new>
#include <iterator>
#include "disjsets.hpp"
#include "container_out.hpp"
#include "verbose_global.hpp"

/*! \namespace clusteringop
  \brief Clustering operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace clusteringop {

  
/*\cite{Sibson1973:HierarchicalCluster:SLINK:1973}
 */
template  <typename T_REAL>
std::vector<std::pair<uintidx,T_REAL> >
sortLambda
(const uintidx    *aiarrayiu_pi,
 const T_REAL     *aiarrayrt_lambda,
 const uintidx    aiui_numObj
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::sortLambda";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input aiarrayiu_pi[" << aiarrayiu_pi << "]\n"
	      << " input aiarrayrt_lambda[" << aiarrayrt_lambda << "]\n"
	      << " input aiui_numObj = " << aiui_numObj
      //<< "\naiui_numClusterK = " << aiui_numClusterK 
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<std::pair<uintidx,T_REAL> > lovector_pairIdxLambda;
  lovector_pairIdxLambda.reserve(aiui_numObj);

  uintidx luintidx_idxObj = 0;
  while ( luintidx_idxObj < aiui_numObj ) {
    lovector_pairIdxLambda.emplace_back(luintidx_idxObj++,*aiarrayrt_lambda);
    ++aiarrayrt_lambda;
  }

  std::sort
    (lovector_pairIdxLambda.begin(),
     lovector_pairIdxLambda.end(), 
     [](const std::pair<uintidx,T_REAL> &left, const std::pair<uintidx,T_REAL> &right) 
     {
       return left.second < right.second; //descend
     }
     );
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

    for (const auto& pairIdxFuncFitness: lovector_pairIdxLambda)
	{
	  std::cout << '(' <<pairIdxFuncFitness.first << "," << pairIdxFuncFitness.second  << ')';
	}
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lovector_pairIdxLambda;
  
}

/*\cite{Sibson1973:HierarchicalCluster:SLINK:1973}
 */
template  <typename T_REAL>
ds::DisjSets
pointerToDisjSets
(const std::vector<std::pair<uintidx,T_REAL> >& aivector_pairIdxLambda,
 const uintidx                                  *aiarrayiu_pi,
 /*const T_REAL     *aiarrayrt_lambda,
   const uintidx    aiui_numObj,*/
 const uintidx    aiui_numClusterK
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::pointerToDisjSets";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(aivector_pairIdxLambda[" << &aivector_pairIdxLambda << "]\n"
      //      << "  input aiarrayrt_lambda[" << aiarrayrt_lambda << ']'
      // << "\naiui_numObj = " << aiui_numObj
	      << " aiui_numClusterK = " << aiui_numClusterK 
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  /* std::vector<std::pair<uintidx,T_REAL> > lvectorT_tmpPairIdxLambda;
  lvectorT_tmpPairIdxLambda.reserve(aiui_numObj);

  uintidx luintidx_idxObj = 0;
  while ( luintidx_idxObj < aiui_numObj ) {
    lvectorT_tmpPairIdxLambda.emplace_back(luintidx_idxObj++,*aiarrayrt_lambda);
    ++aiarrayrt_lambda;
  }

  std::sort
    (lvectorT_tmpPairIdxLambda.begin(),
     lvectorT_tmpPairIdxLambda.end(), 
     [](const std::pair<uintidx,T_REAL> &left, const std::pair<uintidx,T_REAL> &right) 
     {
       return left.second < right.second; //descend
     }
     );

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "SORT FUNTION FITNESS:  IN"
	      << '(' << geiinparam_verbose << ')';
    std::cout << "\nSORT FUNCTION FITNESS: OUT"
	      << '(' << geiinparam_verbose << ")\n";
      for (const auto& pairIdxFuncFitness: lvectorT_tmpPairIdxLambda)
	{
	  std::cout << '(' <<pairIdxFuncFitness.first << "," << pairIdxFuncFitness.second  << ')';
	}
      std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  */

  const uintidx    liui_numObj =  (uintidx)aivector_pairIdxLambda.size();
  ds::DisjSets lodisjsets_memberCluster( liui_numObj );
  uintidx  lui_idxMinLambda = 0;

  //std::cout << "CLUSTER " << aiui_numClusterK << " lodisjsets_memberCluster.getNumSet(): " << lodisjsets_memberCluster.getNumSet() <<  std::endl;
  while ( (lui_idxMinLambda < liui_numObj) &&
	  (aiui_numClusterK < lodisjsets_memberCluster.getNumSet()) )
    {

    T_REAL  lrt_minLambda = aivector_pairIdxLambda.at(lui_idxMinLambda).second;

    while ( lrt_minLambda >=  aivector_pairIdxLambda.at(lui_idxMinLambda).second ) {
      
      lodisjsets_memberCluster.merge
	(aivector_pairIdxLambda.at(lui_idxMinLambda).first,
	 aiarrayiu_pi[aivector_pairIdxLambda.at(lui_idxMinLambda).first]
	 );
      ++lui_idxMinLambda;
    }
    
  }
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    lodisjsets_memberCluster.print(std::cout,lpc_labelFunc,',');  
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lodisjsets_memberCluster;
  
}

/*
  R. Sibson
  SLINK: An optimally efficient algorithm for the single-
  link cluster method
  \cite{Sibson1973:HierarchicalCluster:SLINK:1973}
*/
template  <typename T_REAL,
	   typename INPUT_ITERATOR,
	   typename FUNCTION_DISTANCE
	   //   typename T_FEATURE
	   >
void 
slink
(uintidx                        *aoarrayiu_pi,
 T_REAL                         *aoarrayrt_lambda,
 const INPUT_ITERATOR       aiiterator_instfirst,
 const INPUT_ITERATOR       aiiterator_instlast,
 FUNCTION_DISTANCE              function_distance
 //const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 //const long ail_i1
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::slink";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output aoarrayiu_pi[" << aoarrayiu_pi << "]\n"
	      << " output aoarrayrt_lambda[" << aoarrayrt_lambda << ']'
      //   << " input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
      //      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  const uintidx lui_nobj(uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)));
  T_REAL*       larrayrt_distance = NULL;

  try {
    larrayrt_distance = new T_REAL[lui_nobj];
  }
  catch (std::bad_alloc& ba) {
    std::cerr << "bad_alloc caught in clustering_operator_hierarchical.hpp: " << ba.what() << '\n';
  }

    /* initialize pi and lambda for a single point representation */

    aoarrayiu_pi[0] = 0;
    aoarrayrt_lambda[0] = std::numeric_limits<T_REAL>::max();

    /*
#ifdef __VERBOSE_YES
    geverbosepc_labelstep = "ITERATIVELY ADD THE REMAINING N-1 POINTD TO THE HIERARCHY";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << geverbosepc_labelstep  
		<< ": IN(" << geiinparam_verbose << ')'
		<< std::endl;
    }
#endif //__VERBOSE_YES
    */
    /*iteratively add the remaining N-1 points to the hierarchy */

    INPUT_ITERATOR liter_iInst = aiiterator_instfirst;
    
    for (uintidx lui_i = 1; lui_i < lui_nobj; lui_i++) {

      ++liter_iInst;
      
	aoarrayiu_pi[lui_i] = lui_i;
	aoarrayrt_lambda[lui_i] = std::numeric_limits<T_REAL>::max();

	/* calculate and store a row of the distance matrix for lui_i*/

	//INPUT_ITERATOR liter_jInst = aiiterator_instfirst;

	uintidx lui_j = 0;
	for  (auto liter_jInst = aiiterator_instfirst;
	      liter_jInst != liter_iInst;
	      liter_jInst++)
	  {
	  //for (uintidx lui_j = 0; lui_j < lui_i; lui_j++) {
	  larrayrt_distance[lui_j] = //calc_distance(i,lui_j);
	    function_distance(*liter_iInst,*liter_jInst);
	    /* aifunc2p_dist 
	    (aivectorptinst_instances.at(lui_i)->getFeatures(),
	     aivectorptinst_instances.at(lui_j)->getFeatures(),
	     data::Instance<T_FEATURE>::getNumDimensions()
	     );*/
	    
	  ++lui_j;
	}
	
	for (lui_j = 0; lui_j < lui_i; lui_j++) {
	    uintidx lui_next = aoarrayiu_pi[lui_j];
	    if (aoarrayrt_lambda[lui_j] < larrayrt_distance[lui_j]) {
	      larrayrt_distance[lui_next] = std::min(larrayrt_distance[lui_next],larrayrt_distance[lui_j]);
	    }
	    else {
	      larrayrt_distance[lui_next] = std::min(aoarrayrt_lambda[lui_j],larrayrt_distance[lui_next]);
		aoarrayiu_pi[lui_j] = lui_i;
		aoarrayrt_lambda[lui_j] = larrayrt_distance[lui_j];
	      }
	}

	/* relabel clusters if necessary */

	for (lui_j = 0; lui_j <lui_i; lui_j++) {
	    uintidx lui_next = aoarrayiu_pi[lui_j];
	    if (aoarrayrt_lambda[lui_next] < aoarrayrt_lambda [lui_j])
	      aoarrayiu_pi[lui_j] = lui_i;
	}
	
      }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

    std::ostringstream lostrstream_labelPi;
      lostrstream_labelPi << "<Pi:" << lpc_labelFunc;
       
      inout::containerprint
	(aoarrayiu_pi,
	 aoarrayiu_pi+lui_nobj,
	 std::cout,
	 lostrstream_labelPi.str().c_str(),
	 ','
	 );
      std::cout << '\n';

      std::ostringstream lostrstream_labelLambda;
      lostrstream_labelLambda << "<LAMBDA:" << lpc_labelFunc;
      inout::containerprint
	(aoarrayrt_lambda,
	 aoarrayrt_lambda+lui_nobj,
	 std::cout,
	 lostrstream_labelLambda.str().c_str(),
	 ','
	 );
    
      std::cout << std::endl;
      
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  

  if (larrayrt_distance != NULL ) 
    delete [] larrayrt_distance;
  
}

} /*END namespace clusteringop*/

#endif /*__CLUSTERING_OPERATOR_HIERARCHICAL_HPP*/
