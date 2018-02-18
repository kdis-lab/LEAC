/*! \file container_instance.hpp
 *
 * \brief container utilities for instances
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CONTAINER_INSTANCE_HPP
#define CONTAINER_INSTANCE_HPP

#include <vector>
#include "instance.hpp"
#include "matrix.hpp"

#include "verbose_global.hpp"

/*! \namespace data
  \brief Module for the handling of instances or also called objects or points.
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \fn auto toMatrixRowTrans(INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, bool aiWithOutHomogeneousCoord ) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
    \brief Instaces matrix transpose.
    \details Gets a matrix from a container, where elements are the characteristics of instances stored by columns.
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast a InputIterator to the final positions of the sequence of instances
    \param aiWithOutHomogeneousCoord a  boolean to include another element and represent an instance in homogeneous coordinates
 */  
template < typename INPUT_ITERATOR>
auto 
toMatrixRowTrans
(INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR         aiiterator_instlast,
 bool                         aiWithOutHomogeneousCoord                                
 ) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "data::toMatrixRowTrans";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      <<  ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  const uintidx lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)); 
  const uintidx lui_numRows = (aiWithOutHomogeneousCoord)?
    Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions()-1:
    Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions();
  
  mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>  lomatrixrowt_instancesTans
    (lui_numRows,
     (uintidx) lui_numInstances
     );
  
  uintidx lui_i = 0;
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    const decltype((*aiiterator_instfirst)->getAttribute(0))*
      liter_iInstance =  (*aiiterator_instfirst)->getFeatures();
    for ( uintidx lui_j = 0; lui_j < lui_numRows; lui_j++) { 
      lomatrixrowt_instancesTans(lui_j,lui_i) = liter_iInstance[lui_j];
    }
    ++lui_i;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMatrixInstances;
    lostrstream_labelMatrixInstances
      << "<MATRIXINSTANCES:"
      << geverbosepc_labelstep << ':'
      << lpc_labelFunc
      << ":lomatrixrowt_instancesTans["
      << geverboseui_idproc << ':'
      << &lomatrixrowt_instancesTans << ']';
    lomatrixrowt_instancesTans.print
      (std::cout,
       lostrstream_labelMatrixInstances.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lomatrixrowt_instancesTans;
}

/*! \fn auto toMatrixRow(INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
    \brief Instaces matrix
    \details Gets a matrix from a container, where elements are the characteristics of instances stored by columns.
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast a InputIterator to the final positions of the sequence of instances
 */
template < typename INPUT_ITERATOR>
auto 
toMatrixRow
(INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR         aiiterator_instlast
 ) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
{

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "data::toMatrixRow";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
              << ":  IN(" << geiinparam_verbose << ")\n"
      	      << "(input INPUT_ITERATOR aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      <<  ')'
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  const uintidx lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)); 
  mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>  lomatrixrowt_instances
    (lui_numInstances,
     Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions()
     );
  
  uintidx lui_i = 0;
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    const Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>*
      liter_iInstance = *aiiterator_instfirst;
    decltype((*aiiterator_instfirst)->getAttribute(0))* lmatrixrowT_Ci =
      lomatrixrowt_instances.getRow(lui_i);
    interfacesse::copy 
      (lmatrixrowT_Ci,
       liter_iInstance->getFeatures(),
       Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions()
       );
    ++lui_i;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMatrixInstances;
    lostrstream_labelMatrixInstances
      << "<MATRIXINSTANCES:"
      << geverbosepc_labelstep << ':'
      << lpc_labelFunc
      << ":lomatrixrowt_instances["
      << geverboseui_idproc << ':'
      << &lomatrixrowt_instances << ']';
    lomatrixrowt_instances.print
      (std::cout,
       lostrstream_labelMatrixInstances.str().c_str(),
       ',',
       ';'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lomatrixrowt_instances;
}

} /*END data*/

#endif /*CONTAINER_INSTANCE_HPP*/

