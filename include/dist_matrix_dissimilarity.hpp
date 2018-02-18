/*! \file dist_matrix_dissimilarity.hpp
 *
 * \brief distance matrix between objects
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __DIST_MATRIX_DISSIMILARITY_HPP
#define __DIST_MATRIX_DISSIMILARITY_HPP

#include "matrix_triangular.hpp"
#include "dist.hpp"

/*! \namespace dist
  \brief Module for definition of distance between objects or instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  dist {
  
/*! \fn mat::MatrixTriang<T_DIST> getMatrixDissimilarity(INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Triangular distance matrix 
    \details Returns the triangular distance matrix using a specified distance measure.
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast a InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
*/
template < typename T_FEATURE,
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
mat::MatrixTriang<T_DIST>
getMatrixDissimilarity
(INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "dist::getMatrixDissimilarity";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "\t input  aifunc2p_dista\n"
	      << "\t)\n";
  }
#endif //__VERBOSE_YES

  const uintidx lui_numInstances =
    uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  mat::MatrixTriang<T_DIST>  lomatrixtriagt_dissimilarity(lui_numInstances);
  const INPUT_ITERATOR   lconstiterator_instfirst = aiiterator_instfirst;
    
  for(uintidx lui_i = 0; lui_i < lui_numInstances; lui_i++) {
    
    const T_FEATURE* linst_interi =
      ((data::Instance<T_FEATURE>*) *aiiterator_instfirst)->getFeatures();

    INPUT_ITERATOR  linstiterator_j = lconstiterator_instfirst;
    
    for(uintidx lui_j = 0; lui_j < lui_i; lui_j++) {
      
      const T_FEATURE* linst_interj =
	((data::Instance<T_FEATURE>*) *linstiterator_j)->getFeatures();
      
      lomatrixtriagt_dissimilarity(lui_i,lui_j) = 
	aifunc2p_dist
	(linst_interi,
	 linst_interj,
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
      ++linstiterator_j;
    }
    ++aiiterator_instfirst;
    lomatrixtriagt_dissimilarity(lui_i,lui_i) = T_DIST(0);
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "output mat::MatrixTriang: lomatrixtriagt_dissimilarity["
	      << &lomatrixtriagt_dissimilarity << "]\n";
    lomatrixtriagt_dissimilarity.print();
  }      
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lomatrixtriagt_dissimilarity;
}

  
} /*END namespace dist 
   */

#endif /*__DIST_MATRIX_DISSIMILARITY_HPP*/
