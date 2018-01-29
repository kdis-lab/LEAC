/*! \file container_out.hpp
 *
 * \brief container out
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CONTAINER_OUT_HPP
#define __CONTAINER_OUT_HPP

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm> //std::copy
#include <iterator>
#include "common.hpp"

#include "insertion_operator.hpp"

#ifdef __VERBOSE_YES

#include <iostream>
#include <iomanip>
#include <sstream>

extern int        geiinparam_verbose;
extern int        geiinparam_verboseMax;
extern const char *geverbosepc_labelstep;
extern uintidx    geverboseui_idproc;

#endif /*__VERBOSE_YES*/

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace  inout {
  
template<typename INPUT_ITERATOR>
void 
containerprint
(INPUT_ITERATOR   iterator_first,
 INPUT_ITERATOR   iterator_last,
 std::ostream     &os=std::cout,
 const char       *aipc_label = "",
 const char       aic_delimCoef=','
 )
{  
  const uintidx lui_sizeContainer(uintidx(std::distance(iterator_first,iterator_last)));

#if defined(__VERBOSE_YES)
    os << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << &iterator_first << ']'
       << ":length" << aic_delimCoef << lui_sizeContainer << '>';
#else
    os << aipc_label
       << ":length" << aic_delimCoef << lui_sizeContainer << '>';
#endif
     
  if ( lui_sizeContainer > 0 ) {
    os << *iterator_first++;
  while( iterator_first != iterator_last ) {
      os << aic_delimCoef << *iterator_first++; 
    } 
  }
}


template<typename INPUT_ITERATOR, typename FUNCGETVALUE>
void 
containerprint
(INPUT_ITERATOR   iterator_first,
 INPUT_ITERATOR   iterator_last,
 FUNCGETVALUE     func_getValue,
 std::ostream     &os=std::cout,
 const char       *aipc_label = "",
 const char       aic_delimCoef=','
 )
{  
  const uintidx lui_sizeContainer(uintidx(std::distance(iterator_first,iterator_last)));

#if defined(__VERBOSE_YES)
    os << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << &iterator_first << ']'
       << ":length" << aic_delimCoef << lui_sizeContainer << '>';
#else
    os << aipc_label
       << ":length" << aic_delimCoef << lui_sizeContainer << '>';
#endif
     
  if ( lui_sizeContainer > 0 ) {
    os << func_getValue(*iterator_first);
    iterator_first++;			
  while( iterator_first != iterator_last ) {
    os << aic_delimCoef << func_getValue(*iterator_first);
      iterator_first++;
    } 
  }
}

/*! \fn void  maplistprint(const std::map<T,std::list<T> > &aimaplistadj_g, std::ostream &out=std::cout, const char *aipc_label = "", const char aic_delimCoef = ',', const char aic_delimRow = ';')
    \brief Prints graph of a map structure with list of adjacencies
    \details
    \param aimaplistadj_g a map de list
    \param out  a std::ostream
    \param aipc_label a label of graph
    \param aic_delimCoef a separtor item separator
    \param aic_delimRow a separtor item separator
 */
template<typename T>
void  maplistprint
(const std::map<T,std::list<T> > &aimaplistadj_g,
 std::ostream &out=std::cout,
 const char *aipc_label   = "",
 const char aic_delimCoef = ',',
 const char aic_delimRow  = ';'
 )
{
  //std::map<T,std::list<T> >::const_iterator
  auto lconstit_listAdj = aimaplistadj_g.begin();
  out << aipc_label;
  while ( lconstit_listAdj != aimaplistadj_g.end() ) {
    out << lconstit_listAdj->first << aic_delimCoef;
    out << lconstit_listAdj->second;
    if ( ++lconstit_listAdj != aimaplistadj_g.end() )
      out << aic_delimRow;
  }
  out << std::endl;
}


/*! \fn void  print(const std::vector<std::list<uintidx> > &aivectorlistadj_g, std::ostream &out=std::cout, const char *aipc_label   = "", const char aic_delimCoef = ',', const char aic_delimRow  = ';')
    \brief Prints graph of a vector structure with list of adjacencies
    \details
    \param aivectorlistadj_g a map de list
    \param out  a std::ostream
    \param aipc_label a label of graph
    \param aic_delimCoef a separtor item separator
    \param aic_delimRow a separtor item separator
 
void  print
(const std::vector<std::list<uintidx> > &aivectorlistadj_g,
 std::ostream &out=std::cout,
 const char *aipc_label   = "",
 const char aic_delimCoef = ',',
 const char aic_delimRow  = ';'
 )
{
  out << aipc_label
      << aivectorlistadj_g
      << std::endl;
  
}
*/  

 
  
template<typename T>
void  vectorlistprint
(const std::vector<std::list<T> > &aivectorlist_g,
 std::ostream &out=std::cout,
 const char *aipc_label   = "",
 const char aic_delimList  = ';'
 )
{
#if defined(__VERBOSE_YES)
   out << aipc_label << ':'
       << geverbosepc_labelstep
       << ":id[" << geverboseui_idproc << '-' << &aivectorlist_g << ']'
       << ":length," << aivectorlist_g.size() << '>';
#else
   out << aipc_label << '>';
#endif
      
  if (  aivectorlist_g.size() > 0 ) {
    uintidx lui_vectorSize =  (uintidx)  aivectorlist_g.size()-1;
    for (uintidx lui_i = 0; lui_i < lui_vectorSize; lui_i++) {
      out << lui_i;
      if ( aivectorlist_g[lui_i].size() > 0 ) {
	out << OutFileName::getDelim();
	::operator<< <T> (out,aivectorlist_g[lui_i]);
      }
      out << aic_delimList;
    }
    out << lui_vectorSize;
    if ( aivectorlist_g[lui_vectorSize].size() > 0 ) {
      out << OutFileName::getDelim();
      ::operator<< <T> (out,aivectorlist_g[lui_vectorSize]);
    }
  }
}
  
  
} /*END namespace inout*/

/*
template <class T>
std::vector<T>
vectorutil_keepItems
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

template <class T>
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

std::vector<std::string> 
vectorutil_foundStrings
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

*/
#endif /*__CONTAINER_OUT_HPP*/
