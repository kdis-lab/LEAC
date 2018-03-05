/*! \file container_out.hpp
 *
 * \brief container out
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
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


#endif /*__CONTAINER_OUT_HPP*/
