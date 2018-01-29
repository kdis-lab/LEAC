/*! \file disjsets.hpp
 *
 * \brief disjoint sets
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef DISJSETS_HPP
#define DISJSETS_HPP

#include <iostream>
#include <vector>
#include "common.hpp"

#include "verbose_global.hpp"

/*! \namespace ds
  \brief Data structure
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace ds {

/*! \class DisjSets
  \brief Disjoint-set
  \details A disjoint-set data structure, also called a union-find data structure or merge-find set, is a data structure that keeps track of a set of elements partitioned into a number of disjoint (nonoverlapping) subsets.
*/
class DisjSets
{  
public:

  DisjSets()
    : _m(0)
    , _vectorst_parent(0)
    , _vectorui_rank(0)
  {
  
  }
  
  DisjSets(const uintidx ai_numObj)
    : _m(ai_numObj)
    , _vectorst_parent(ai_numObj)
    , _vectorui_rank(ai_numObj,0)
  {
    
     for (uintidx i=0; i < ai_numObj; i++) {
      _vectorst_parent[i] = i;
     }
  }

  
  //move constructor
  DisjSets(DisjSets &&aids_disjset)
    : _m(aids_disjset._m)
    , _vectorst_parent(std::move( aids_disjset._vectorst_parent) )
    , _vectorui_rank(std::move (aids_disjset._vectorui_rank) ) 
  {
  }

  //copy constructor
  DisjSets(const DisjSets &aids_disjset)
    : _m(aids_disjset._m)
    , _vectorst_parent(aids_disjset._vectorst_parent)
    , _vectorui_rank(aids_disjset._vectorui_rank) 
  {
   
  }

  ~ DisjSets()
  {
  }

  DisjSets&
  operator=(DisjSets &&aids_disjset)
  {
    if ( this != &aids_disjset ) {
      _m = aids_disjset._m;
      _vectorst_parent = aids_disjset._vectorst_parent;
      _vectorui_rank = aids_disjset._vectorui_rank;

      aids_disjset._m = 0;
    }

    return *this;
   
  }


  DisjSets&
  operator=(const DisjSets &aids_disjset)
  {
    if ( this != &aids_disjset ) {
      _m = aids_disjset._m;
      _vectorst_parent = aids_disjset._vectorst_parent;
      _vectorui_rank = aids_disjset._vectorui_rank;
    }

    return *this;
   
  }
  
  
  inline void make_set(uintidx x)
  {
    _vectorst_parent[x] = x;
    _vectorui_rank[x] = 0;
  }
  
  // Find set that element i belongs to, represented as an index
  uintidx find(uintidx x)
  {
    if ( x != _vectorst_parent[x] )
      _vectorst_parent[x] = find(_vectorst_parent[x]);
    
    return _vectorst_parent[x];
  }

  const uintidx constfind(uintidx x) const 
  {

    while ( x != _vectorst_parent[x])
      x = _vectorst_parent[x];
    return x;
  }

  inline
  void merge(uintidx x, uintidx y)
  {
    link(find(x),find(y));
    --_m;
  }

  inline uintidx getNumSet()
  {
    return _m;
  }

  inline const  uintidx size() const 
  {
    return (uintidx) _vectorst_parent.size(); 
  }
  
  void  print
  (std::ostream &os=std::cout,
   const char *aipc_label = "",
   const char aic_delimCoef=','
   ) const
  {

#if defined(__VERBOSE_YES)
    os << "<DISJSETS:"
       << geverbosepc_labelstep
       << ':' <<  aipc_label
       << ":id[" << geverboseui_idproc << '-' << this << ']'
       << ":length," << _vectorst_parent.size()
#else
   os  << aipc_label
#endif
   
       << ":number of set," << _m << '>';
    
    if ( _vectorst_parent.size() > 0 ) {
      uintidx  luintidx_last = _vectorst_parent.size() - 1;   
      for (uintidx luintidx_i = 0; luintidx_i < luintidx_last ; ++luintidx_i ) {
	os << constfind(luintidx_i) << aic_delimCoef; 
      }
      os << constfind(luintidx_last);
    }

  }
    

protected:
  
  uintidx                _m;
  std::vector<uintidx>   _vectorst_parent;
  std::vector<unsigned> _vectorui_rank;
  
private:

  void link(uintidx x, uintidx y)
  {   
    if ( _vectorui_rank[x] > _vectorui_rank[y] )
      _vectorst_parent[y] = x;
    else {
      _vectorst_parent[x] = y;
      if (_vectorui_rank[x] == _vectorui_rank[y] )
	_vectorui_rank[y] += 1;
    }  
  }
  
};

} /*END namespace ds*/

#endif
