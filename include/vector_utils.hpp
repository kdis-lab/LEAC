/*! \file vector_utils.hpp
 *
 * \brief Operations for vector manipulation
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __VECTOR_UTILS_HPP
#define __VECTOR_UTILS_HPP

#include <iostream>
#include <vector>
//#include <list>
#include <string>
#include <algorithm> 
//#include <iterator>
#include "common.hpp"

#include "verbose_global.hpp"


/*! \namespace vectorutils
  \brief Vector utils
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace vectorutils {


template<typename INPUT_ITERATOR,
	 typename T_SIZE>
std::vector<INPUT_ITERATOR>
partition
(INPUT_ITERATOR     aiiterator_first,
 INPUT_ITERATOR     aiiterator_last,
 const  T_SIZE      aist_numPartitions
 )
{
  auto lt_numItems = std::distance(aiiterator_first,aiiterator_last);
  auto lti_sizePartitions = lt_numItems / aist_numPartitions;
  auto lt_rest = lt_numItems %  aist_numPartitions;

  std::vector<INPUT_ITERATOR> lovector_iter;

  lovector_iter.reserve(aist_numPartitions+1);
  auto lit_endPartition   = aiiterator_first;

  lovector_iter.push_back(lit_endPartition);

  T_SIZE lt_i = 0;
  do {
    decltype(lt_numItems) lit_fraccion; 
    if ( lt_rest > 0) {
      lit_fraccion = 1;
      --lt_rest;
    }
    else {
      lit_fraccion = 0;
    }

    std::advance(lit_endPartition,lti_sizePartitions+lit_fraccion);
    
    lovector_iter.push_back(lit_endPartition);
    
    ++lt_i;
  }  while( lt_i < aist_numPartitions);
     
  return lovector_iter;
  
}

template <class T>
std::vector<T>
keepItems
(const std::vector<T>  &aivector_b,
 std::vector<uintidx>  &aivectorcidx_items
 )
{
  std::vector<T> lovector_keep
    (aivectorcidx_items.size());

  for ( uintidx luintidx_i = 0; luintidx_i < lovector_keep.size(); luintidx_i++) {
    lovector_keep[luintidx_i] = aivector_b.at(aivectorcidx_items.at(luintidx_i));
  }
   
  return lovector_keep;
}

/*template <class T>
std::vector<T>
vectorutil_normalized01
(const std::vector<T>  &aivectort_data)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "vectorutil_normalized01";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input aivectort_data[" << &aivectort_data << "]:"
	      << aivectort_data
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  T lT_min;
  T lT_max;
  T lT_delta;
  std::vector<T> lovector_norm;
  
  lovector_norm.reserve(aivectort_data.size());
  lT_min = 
    *std::min_element(std::begin(aivectort_data),std::end(aivectort_data));
  lT_max =
    *std::max_element(std::begin(aivectort_data),std::end(aivectort_data));
  lT_delta = lT_max - lT_min;

  std::for_each
    (aivectort_data.begin(),
     aivectort_data.end(),
     [&](const T& liter_elem)
     {
       lovector_norm.push_back((liter_elem - lT_min) / lT_delta);
     }
     ); 

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "lovector_norm[" << &lovector_norm << "]:"
	      << lovector_norm
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lovector_norm;
}

*/

std::vector<std::string> 
foundStrings
(const std::vector<std::string> &aivectorstr_stringsToSearch,
 const std::string              &aistr_stringToFound
 )
{
  std::vector<std::string> lovector_stringsFound(0);
  lovector_stringsFound.reserve(aivectorstr_stringsToSearch.size()/2);

  std::for_each 
    (std::begin(aivectorstr_stringsToSearch), 
     std::end(aivectorstr_stringsToSearch), 
     [&](const std::string literstr_nameFile) 
     {
       std::size_t found = literstr_nameFile.find(aistr_stringToFound);
       if ( found != std::string::npos) {
	 lovector_stringsFound.push_back(literstr_nameFile);
       }
     }
     );

  return lovector_stringsFound;

}

  
} /*END namespace vectorutils*/

#endif /*__VECTOR_UTILS_HPP*/
