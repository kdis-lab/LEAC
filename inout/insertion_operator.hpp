/*! \file insertion_operator.hpp
 *
 * \brief insertion operator
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#include <iostream>
#include "outfilename.hpp"

#ifndef __INSERTION_OPERATOR_HPP
#define __INSERTION_OPERATOR_HPP


template< typename T
	  , template<typename ELEM, typename ALLOC=std::allocator<ELEM> > class Container
	  >
std::ostream& operator<< (std::ostream& out, const Container<T>& container)
{
  if ( container.size() > 0 ) {
    typename Container<T>::const_iterator beg = container.begin();

    if (beg != container.end())
      out << *beg;
    beg++;
    while(beg != container.end())
      {
	out << inout::OutFileName::getDelim() << *beg++;
      }
  }

  return out;
}

#endif /*__INSERTION_OPERATOR_HPP*/
