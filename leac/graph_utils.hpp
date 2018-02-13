/*! \file graph_utils.hpp
 *
 * \brief algorithms and utilities for graphs
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef GRAPH_UTILS_HPP
#define GRAPH_UTILS_HPP

#include <iostream>
#include <vector>
#include <utility>
#include <set> 
#include <list>
#include <map>
#include <limits> //std::numeric_limits
#include "matrix_triangular.hpp"
#include "instance.hpp"
#include "disjsets.hpp"
#include "common.hpp"

#include "insertion_operator.hpp"
#include "verbose_global.hpp"

/*! \namespace graph
  \brief Algorithms and utilities for graphs
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace graph {

  
typedef enum {VERTEXSTATE_WHITE=0, VERTEXSTATE_GRAY=1, VERTEXSTATE_BLACK=2} GraphStateVertex;


/*! \class AlgorithmsDFS
  \brief Depth-first search algorithms
  \details 
*/
class AlgorithmsDFS
{
public:
  AlgorithmsDFS(std::vector<std::list<uintidx> > &aivectorlistadj_g)
  : _vectorlistadj_g(aivectorlistadj_g)
  , _vectorenum_color(aivectorlistadj_g.size()) 
  , _vectorvertexidx_pi(aivectorlistadj_g.size())
  , f(aivectorlistadj_g.size())
  , d(aivectorlistadj_g.size())
  {
   
  }

  ~AlgorithmsDFS() {} 

  void resolve() 
  {
    uintidx luintidx_numVertex = (uintidx) _vectorlistadj_g.size();

    for ( uintidx luintidx_ui = 0; luintidx_ui < luintidx_numVertex; luintidx_ui++ ) { 
      _vectorenum_color[luintidx_ui]   = VERTEXSTATE_WHITE;
      _vectorvertexidx_pi[luintidx_ui] = UINTIDX_NIL;  
    }
    _uintidx_time = 0;
   
   for ( uintidx u = 0; u < luintidx_numVertex; u++ ) { 
      if ( _vectorenum_color[u] == VERTEXSTATE_WHITE )
	dfs_visit(u);
   }
  }
  
  void resolve(uintidx aiui_vertexU)
  {
    uintidx luintidx_numVertex = 
      (uintidx) _vectorlistadj_g.size();

    for ( uintidx luintidx_ui = 0; luintidx_ui < luintidx_numVertex; luintidx_ui++ ) { 
      _vectorenum_color[luintidx_ui]   = VERTEXSTATE_WHITE;
      _vectorvertexidx_pi[luintidx_ui] = UINTIDX_NIL; 
      f[luintidx_ui] = UINTIDX_NIL;  
      d[luintidx_ui] = UINTIDX_NIL;  
    }
    _uintidx_time = 0;
    dfs_visit(aiui_vertexU);
   
  }

  void dfs_visit(uintidx aiui_vertexU)
  {
    _vectorenum_color[aiui_vertexU] = VERTEXSTATE_GRAY;
    ++_uintidx_time;
    d[aiui_vertexU] = _uintidx_time;
   
    std::list<uintidx> &lvectorvertexidx_adj =  _vectorlistadj_g[aiui_vertexU]; 

    /*Obtain the list of vertices adjacent to aiui_vertexU
     lui_vertexV is a vertex adjacent to aiui_vertexU
    */
    for (auto lui_vertexV: lvectorvertexidx_adj) { 
      if ( _vectorenum_color[lui_vertexV] == VERTEXSTATE_WHITE ) {
	_vectorvertexidx_pi[lui_vertexV] = aiui_vertexU;
	dfs_visit(lui_vertexV);
      }
    }
    _vectorenum_color[aiui_vertexU] = VERTEXSTATE_BLACK;
    ++_uintidx_time;
    f[aiui_vertexU] = _uintidx_time;
  }
  
  /*getVertexReachable:
   */
  std::list<uintidx>
  getVertexReachable(uintidx aiui_vertexU)
  {

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "AlgorithmsDFS::getVertexReachable:  IN"
		<< '(' << geiinparam_verbose << ')'
		<< "\n\t( input uintidx aiui_vertexU: aiui_vertexU = "
		<< aiui_vertexU
		<< "\n\t)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    std::list<uintidx> lolistvertexidx_reachable;
    uintidx lui_numVertices = 
      (uintidx) _vectorlistadj_g.size();
    
    this->resolve(aiui_vertexU);
    for ( uintidx luintidx_ui = 0; luintidx_ui < lui_numVertices; luintidx_ui++ ) {
      if (d[luintidx_ui] != UINTIDX_NIL ) {
	lolistvertexidx_reachable.push_back(luintidx_ui);
      } 
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "AlgorithmsDFS::getVertexReachable: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< "\noutput std::list<uintidx> lolistvertexidx_reachable[" 
		<< &lolistvertexidx_reachable << "]\n"
		<< lolistvertexidx_reachable
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    return lolistvertexidx_reachable;

  }

  void print() 
  {
    for ( uintidx li_i = 0; li_i < _vectorvertexidx_pi.size(); li_i++) {
      std::cout << li_i << ":\t" << _vectorvertexidx_pi[li_i] 
		<< "\t" << d[li_i] << "\t" << f[li_i] << "\n";
    }
  }
private:

  std::vector<std::list<uintidx> >   &_vectorlistadj_g;
  //Establish if the vertex was visited
  std::vector<GraphStateVertex>      _vectorenum_color;
  //Stores parent vertices
  std::vector<uintidx>               _vectorvertexidx_pi;
  //Final time to know the vertex parent
  std::vector<uintidx>               f;
  //￼Initial time to know the vertex parent
  std::vector<uintidx>               d; 
  uintidx                            _uintidx_time;
    
}; /*class AlgorithmsDFS*/



/*! \class AlgorithmsDFSAdjMatrix
  \brief Depth-first search algorithms for a matrix of bits 
  \details 
*/
template < typename T_BITSIZE >
class AlgorithmsDFSAdjMatrix
{
public:
  AlgorithmsDFSAdjMatrix(mat::BitMatrix<T_BITSIZE> &aibitmatrixadj_g)
    : _bitmatrixadj_g(aibitmatrixadj_g)
    , _vectorenum_color(aibitmatrixadj_g.getNumRows()) 
    , _vectorvertexidx_pi(aibitmatrixadj_g.getNumRows())
    , f(aibitmatrixadj_g.getNumRows())
    , d(aibitmatrixadj_g.getNumRows())
  {
  }

  ~AlgorithmsDFSAdjMatrix() {} 
  
  void resolve() 
  {
    for ( uintidx luintidx_ui = 0;
	  luintidx_ui < _bitmatrixadj_g.getNumRows();
	  luintidx_ui++ )
      { 
	_vectorenum_color[luintidx_ui] = VERTEXSTATE_WHITE;
	_vectorvertexidx_pi[luintidx_ui] = UINTIDX_NIL; 
	f[luintidx_ui] = UINTIDX_NIL;  
	d[luintidx_ui] = UINTIDX_NIL;  
      }
    _uintidx_time = 0;
    for ( uintidx luintidx_ui = 0;
	  luintidx_ui < _bitmatrixadj_g.getNumRows();
	  luintidx_ui++ )
      { 
	if ( _vectorenum_color[luintidx_ui] == VERTEXSTATE_WHITE ) {
	  dfs_visit(luintidx_ui);
	}
      }
  }
  
  void dfs_visit(uintidx aiui_vertexU)
  {  
    _vectorenum_color[aiui_vertexU] = VERTEXSTATE_GRAY;
    ++_uintidx_time;
    d[aiui_vertexU] = _uintidx_time;

    //The list of adjacent to the u vertex is obtained 
   
    for ( uintidx luintidx_vi = 0;
	  luintidx_vi < _bitmatrixadj_g.getNumColumns();
	  luintidx_vi++ )
      {
      if ( _bitmatrixadj_g(aiui_vertexU,luintidx_vi) ) {
	if ( _vectorenum_color[luintidx_vi] == VERTEXSTATE_WHITE ) {
	  _vectorvertexidx_pi[luintidx_vi] = aiui_vertexU;
	  dfs_visit(luintidx_vi);
	}
      }
      // v is a vertex a u
    }
    _vectorenum_color[aiui_vertexU] = VERTEXSTATE_BLACK;
    ++_uintidx_time;
    f[aiui_vertexU] = _uintidx_time;
  }

  std::vector<uintidx>&  getPi()
  {
    return _vectorvertexidx_pi; 
  }
  
  
  void print() 
  {
    for ( uintidx li_i = 0; li_i < _vectorvertexidx_pi.size(); li_i++) {
      std::cout << li_i << ":\t" << _vectorvertexidx_pi[li_i] 
		<< "\t" << d[li_i] << "\t" << f[li_i] << '\n';
    }
  }
private:

  mat::BitMatrix<T_BITSIZE>           &_bitmatrixadj_g;
  //Establish if the vertex was visited
  std::vector<GraphStateVertex>      _vectorenum_color;
  //Stores parent vertices
  std::vector<uintidx>               _vectorvertexidx_pi;
  //Final time to know the vertex parent
  std::vector<uintidx>               f;
  //￼Initial time to know the vertex parent
  std::vector<uintidx>               d; 
  uintidx                            _uintidx_time;
}; /*class AlgorithmsDFSAdjMatrix*/


/*! \fn std::map< uintidx, std::list<uintidx> > fromPiToGraphMap(const std::vector<uintidx> &aivectorvertexidx_pi)
    \brief Constructs a list map graph from PI vertices for a Pi vector
    \details
    \param aivectorvertexidx_pi a vector from index Pi
 */
std::map< uintidx, std::list<uintidx> > 
fromPiToGraphMap
(const std::vector<uintidx> &aivectorvertexidx_pi 
)
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "graph::fromPiToGraphMap:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n( input std::vector<uintidx>: &aivectorvertexidx_pi["
	      << &aivectorvertexidx_pi << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::map<uintidx,std::list<uintidx> > lomaplistadj_g;
  
  uintidx lui_numVertices = (uintidx) aivectorvertexidx_pi.size();
  
  for (uintidx luintidx_u = 0; luintidx_u < lui_numVertices; luintidx_u++) {
    
    if ( aivectorvertexidx_pi[luintidx_u] != UINTIDX_NIL) {
      
      typename std::map<uintidx,std::list<uintidx> >::iterator lconstit_listAdj =
	lomaplistadj_g.find(aivectorvertexidx_pi[luintidx_u]);
      
      if ( lconstit_listAdj != lomaplistadj_g.end() ) {
	
	lconstit_listAdj->second.push_back(luintidx_u);
	
      }
      else {
	std::list<uintidx> llist_vertidx(1,luintidx_u); 
	lomaplistadj_g.insert
	  ( std::pair<uintidx,std::list<uintidx> >
	    (aivectorvertexidx_pi[luintidx_u],llist_vertidx)
	    );
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    const char* lpc_label = "graph::fromPiToGraphMap: OUT";    
    std::cout << lpc_label
	      << '(' << geiinparam_verbose << ')'
	      << "\nstd::map<uintidx,std::list<uintidx> > lomaplistadj_g[" 
	      << &lomaplistadj_g << "]"
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lomaplistadj_g;
  
}


/*! \fn std::vector<std::list<T_VERTEX_IDX> > fromPiToGraph(const std::vector<T_VERTEX_IDX> &aivectorvertexidx_pi)
    \brief Constructs a vector list graph from PI vertices for a Pi vector
    \details
    \param aivectorvertexidx_pi a vector from index Pi
 */
template < typename T_VERTEX_IDX >
std::vector<std::list<T_VERTEX_IDX> >
fromPiToGraph
(const std::vector<T_VERTEX_IDX> &aivectorvertexidx_pi)
{
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "graph::fromPiToGraph:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t( input std::vector<T_VERTEX_IDX>: &aivectorvertexidx_pi["
	      << &aivectorvertexidx_pi << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES*/
  
  std::vector<std::list<T_VERTEX_IDX> >
    lovectorlistadj_g
    (aivectorvertexidx_pi.size(),
     std::list<T_VERTEX_IDX>(0)
     );
  T_VERTEX_IDX lui_numVertices = (T_VERTEX_IDX) aivectorvertexidx_pi.size();

  for (T_VERTEX_IDX luintidx_u = 0; luintidx_u < lui_numVertices; luintidx_u++) {
    if ( aivectorvertexidx_pi[luintidx_u] != UINTIDX_NIL) { 
      lovectorlistadj_g[aivectorvertexidx_pi[luintidx_u]].push_back(luintidx_u);
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "graph::fromPiToGraph: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\noutput std::vector<std::list<T_VERTEX_IDX> > lovectorlistadj_g[" 
	      << &lovectorlistadj_g << ']'
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectorlistadj_g;
  
}


/*! \fn std::vector<std::list<T_VERTEX_IDX> > fromPiToGraph(const std::vector<T_VERTEX_IDX> &aivectorvertexidx_pi, const mat::BitArray<T_BITSIZE> &aibitarray_edgeRemains) 
    \brief Constructs a vector list graph from PI vertices for a Pi vector
    \details
    \param aivectorvertexidx_pi a vector from index Pi
    \param aibitarray_edgeRemains a bit array of edge remains
 */
template < typename T_VERTEX_IDX,
	   typename T_BITSIZE
	   >
std::vector<std::list<T_VERTEX_IDX> >
fromPiToGraph
(const std::vector<T_VERTEX_IDX>     &aivectorvertexidx_pi,
 const mat::BitArray<T_BITSIZE>      &aibitarray_edgeRemains
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "graph::fromPiToGraph:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t( input std::vector<T_VERTEX_IDX>: &aivectorvertexidx_pi["
	      << &aivectorvertexidx_pi << ']'
	      << "\n\t input  mat::BitArray<T_BITSIZE>: &aibitarray_edgeRemains [" 
	      << &aibitarray_edgeRemains << ']' 
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::vector<std::list<T_VERTEX_IDX> >
    lovectorlistadj_g
    (aivectorvertexidx_pi.size(),
     std::list<T_VERTEX_IDX>(0)
     );
  T_VERTEX_IDX lui_numVertices = (T_VERTEX_IDX) aivectorvertexidx_pi.size();
  
  uintidx  luintidx_numEdge = 0;
  for (T_VERTEX_IDX luintidx_u = 0; luintidx_u < lui_numVertices; luintidx_u++) {
    if ( aivectorvertexidx_pi[luintidx_u] != UINTIDX_NIL) { 
      if ( !aibitarray_edgeRemains.getBit(luintidx_numEdge) ) { 
	lovectorlistadj_g[aivectorvertexidx_pi[luintidx_u]].push_back(luintidx_u);
      }
      ++luintidx_numEdge;
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "graph::fromPiToGraph: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\noutput std::vector<std::list<T_VERTEX_IDX> > lovectorlistadj_g[" 
	      << &lovectorlistadj_g << "]\n"
	      << std::endl;
    }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectorlistadj_g;
  
}

/*! \fn ds::DisjSets component(const std::vector<uintidx>  &aivectoruintidx_vertexPi)
    \brief Find connected component of a graph
    \details A connected component (or just component) of an undirected graph is a subgraph in which any two vertices are connected to each other by paths \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001 
    \param aivectoruintidx_vertexPi a vector of indices Pi
*/
ds::DisjSets 
component
(const std::vector<uintidx>  &aivectoruintidx_vertexPi
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  const char* lpc_labelComponent = "graph::component";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelComponent
	      << ": IN(" << geiinparam_verbose << ')'
	      << "\n\t  input std::vector<uintidx> &aivectoruintidx_vertexPi[" 
	      << &aivectoruintidx_vertexPi << ']' 
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  uintidx  lui_numVertices = (uintidx) aivectoruintidx_vertexPi.size();
  uintidx  luintidx_numEdge;
  ds::DisjSets aiods_disjsets(lui_numVertices);

  luintidx_numEdge  = 0;
  for ( uintidx luintidx_i = 0; luintidx_i < lui_numVertices; ++luintidx_i) {
    if ( aivectoruintidx_vertexPi[luintidx_i] != UINTIDX_NIL )  { 
      if ( aiods_disjsets.find(aivectoruintidx_vertexPi[luintidx_i])
	   != aiods_disjsets.find(luintidx_i) )
	{
	  aiods_disjsets.merge(aivectoruintidx_vertexPi[luintidx_i], luintidx_i );
	}
      ++luintidx_numEdge;
    }
  }
  
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {    
    std::cout << lpc_labelComponent
	      << ": OUT(" << geiinparam_verbose << ")\n";   
    aiods_disjsets.print(std::cout,lpc_labelComponent,',');  
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return aiods_disjsets;
  
}


/*! \fn ds::DisjSets component(const std::vector<uintidx> &aivectoruintidx_vertexPi, const mat::BitArray<T_BITSIZE>  &aibitarray_edgeRemains)
    \brief Find connected component of a graph
    \details A connected component (or just component) of an undirected graph is a subgraph in which any two vertices are connected to each other by paths \cite Casillas:etal:GAclusteringVarK:GA:2003 
    \param aivectoruintidx_vertexPi a vector of indices Pi
    \param aibitarray_edgeRemains a bit array of edge remains
*/
template < typename T_BITSIZE >
ds::DisjSets 
component
(const std::vector<uintidx>      &aivectoruintidx_vertexPi,
 const mat::BitArray<T_BITSIZE>  &aibitarray_edgeRemains
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  const char* lpc_labelComponent = "graph::component";
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelComponent
	      << ": IN(" << geiinparam_verbose << ')'
	      << "\n\t  input std::vector<uintidx> &aivectoruintidx_vertexPi[" 
	      << &aivectoruintidx_vertexPi << ']' 
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  uintidx  lui_numVertices = (uintidx) aivectoruintidx_vertexPi.size();
  ds::DisjSets aiods_disjsets(lui_numVertices);

  uintidx  lui_numEdges  = 0;
  for ( uintidx luintidx_i = 0; luintidx_i < lui_numVertices; ++luintidx_i) {
    if ( aivectoruintidx_vertexPi[luintidx_i] != UINTIDX_NIL )  {
      if ( !aibitarray_edgeRemains.getBit(lui_numEdges) ) {
	if ( aiods_disjsets.find(aivectoruintidx_vertexPi[luintidx_i])
	     != aiods_disjsets.find(luintidx_i) )
	  aiods_disjsets.merge(aivectoruintidx_vertexPi[luintidx_i], luintidx_i );
      }
      ++lui_numEdges;
    }
  }
  
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
   
    std::cout << lpc_labelComponent
	      << ": OUT(" << geiinparam_verbose << ')'
	      << "\noutput std::vector<std::list<uintidx> > aiods_disjsets[" 
	      << &aiods_disjsets
	      << "]\n";    
    aiods_disjsets.print(std::cout,lpc_labelComponent,',');
        
    std::cout << std::endl;
    
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return aiods_disjsets;
  
}



/*! \fn std::vector<uintidx> prim(const mat::MatrixTriang<T_W> &aimatrixtriagt_w, const uintidx lui_startVertex = 0) 
    \brief Prim’s algorithm
    \details Prim’s algorithm. It is assumed that in the graph for this function all the vertices are connected to everyone, not the graph is required. Used to obtain component
    \param aimatrixtriagt_w a matrix of distances all the vertices are joined with an edge of weight equal to the distance
    \param lui_startVertex a starting vertex
 */
template < typename T_W >  
std::vector<uintidx> 
prim
(const mat::MatrixTriang<T_W> &aimatrixtriagt_w,
 const uintidx                lui_startVertex = 0
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "graph::prim:";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << "  IN(" << geiinparam_verbose << ')'
	      << "\n\t( input mat::MatrixTriang<T_W>: &aimatrixtriagt_w = ["
	      << &aimatrixtriagt_w << ']'
	      << "\n\t input  uintidx lui_startVertex = " << lui_startVertex
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::set<std::pair<T_W,uintidx> >    lset_Q;
  std::vector<uintidx>                 
    lovectoruintidx_pi(aimatrixtriagt_w.getNumRows(), UINTIDX_NIL);
  std::vector<std::pair<T_W,uintidx> > lvector_temp;

  lvector_temp.reserve(aimatrixtriagt_w.getNumRows());
  const uintidx lui_numVertices = (uintidx) aimatrixtriagt_w.getNumRows();
  
  for (uintidx luintidx_i = 0; luintidx_i < lui_numVertices; luintidx_i++) {
    lset_Q.insert(std::make_pair(std::numeric_limits<T_W>::max(),luintidx_i));
  }
  
  lset_Q.erase(lset_Q.find(std::make_pair(std::numeric_limits<T_W>::max(),lui_startVertex)));
  lset_Q.insert(std::make_pair(T_W(0),lui_startVertex));
  
  while ( !lset_Q.empty() ) {
    std::pair<uintidx,T_W> lpair_u = *lset_Q.begin();
    lset_Q.erase(lset_Q.begin()); 
      
    lvector_temp.clear();
    for ( auto  lpair_v: lset_Q ) {

      if ( aimatrixtriagt_w(lpair_u.second,lpair_v.second) < lpair_v.first ) {
	lvector_temp.push_back(lpair_v);
	lovectoruintidx_pi[lpair_v.second] = lpair_u.second;
      }
    }
    for ( auto  lpair_v: lvector_temp) {
      lset_Q.erase(lset_Q.find(lpair_v));
      lpair_v.first = aimatrixtriagt_w(lpair_u.second,lpair_v.second);
      lset_Q.insert(lpair_v);
    }

  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {

    std::cout << lpc_labelFunc
	      << " OUT(" << geiinparam_verbose << ')'
	      << "\noutput vector<uintidx> lovectoruintidx_pi["
	      << &lovectoruintidx_pi << "]\n";
  
    if ( lui_numVertices > 0 ) {
      for (uintidx luintidx_v = 0; luintidx_v < lui_numVertices-1; luintidx_v++) {
	std::cout << luintidx_v << ',' << lovectoruintidx_pi[luintidx_v] << ';';
      }
      std::cout << lui_numVertices-1 << ',' << lovectoruintidx_pi[lui_numVertices-1] << '\n';
    }
    std::ostringstream lostrstream_labelGraph;
   lostrstream_labelGraph << "<GRAPH:" << lpc_labelFunc;
  
    std::vector<std::list<uintidx> >&& lvectorlist_graph = 
      fromPiToGraph
      (lovectoruintidx_pi);

    inout::vectorlistprint
      (lvectorlist_graph,
       std::cout,
       lostrstream_labelGraph.str().c_str(),
       ';'
       );
    std::cout << std::endl;   
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lovectoruintidx_pi;

} /*prim*/


/*! \fn mat::BitMatrix<T_BITSIZE> getAdjacencyMatrix(const T_DIST ait_distAdj, const std::vector<data::Instance<T_FEATURE>* > &aivectorptinst_instances, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
    \brief Gets an adjacency matrix \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
    \details It obtains an adjacency matrix considering that if one vertex is adjacent to another if it is a less or equal distance specified by ait_distAdj
    \param ait_distAdj a real number 
    \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type Dist for example dist::Euclidean
*/ 
template < typename T_BITSIZE,
           typename T_FEATURE,
	   typename T_DIST,
	   typename INPUT_ITERATOR
	   >
mat::BitMatrix<T_BITSIZE> 
getAdjacencyMatrix
(const T_DIST                    ait_distAdj,
 INPUT_ITERATOR                  aiiterator_instfirst,
 const INPUT_ITERATOR            aiiterator_instlast,
 dist::Dist<T_DIST,T_FEATURE>    &aifunc2p_dist,
 const T_BITSIZE                 aii_datatypeBitSize = 0
)
{
 const uintidx  lui_numInstances =
   uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)); 
  mat::BitMatrix<T_BITSIZE> 
    lobitmatrix_adjacency
    ( lui_numInstances,
      lui_numInstances
      );

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "graph::getAdjacencyMatrix";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'      
	      << "\n(input  T_DIST ait_distAdj = " <<  ait_distAdj
	      << "input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "\n input  dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']' 
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  lobitmatrix_adjacency.initialize();

  uintidx lui_i = 0;
   for (INPUT_ITERATOR literator_i = aiiterator_instfirst;
	literator_i != aiiterator_instlast; literator_i++, lui_i++) {
     uintidx lui_j = 0;
    for (INPUT_ITERATOR literator_j = aiiterator_instfirst;
	 literator_j != aiiterator_instlast; literator_j++, lui_j++) {    
      T_DIST lt_dij =  
	aifunc2p_dist
	((*literator_i)->getFeatures(),
	 (*literator_j)->getFeatures(),
	 data::Instance<T_FEATURE>::getNumDimensions() 
	 );
      if ( lt_dij <= ait_distAdj )
	lobitmatrix_adjacency.setBit(lui_i,lui_j);
      
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
   	      << "\nmat::BitMatrix<T_BITSIZE> lobitmatrix_adjacency[" 
	      << &lobitmatrix_adjacency << "]\n"
	      << lobitmatrix_adjacency
	      << '\n';
    std::ostringstream lostrstream_labelGraph;
    lostrstream_labelGraph << "<GRAPH:" << lpc_labelFunc;
    
    lobitmatrix_adjacency.print_listAdj
      (std::cout,
       lostrstream_labelGraph.str().c_str(),
       ',',
       ';'
       );
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lobitmatrix_adjacency;

} /*END getAdjacencyMatrix
   */
 
/*! \fn partition::PartitionDisjSets<T_CLUSTERIDX> nearestNeighbor(const T_REAL aitr_parameterU, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist, const T_CLUSTERIDX aii_datatypeClusterIdx = 0)
    \brief  Nearest neighbor \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
    \details Make a grouping of the closest instances according to the \f$u\f$ parameter.
    \param aitr_parameterU a real number the parameter u
     \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param aifunc2p_dist an object of type Dist for example dist::Euclidean
    \param aii_datatypeClusterIdx 
 */
template < typename T_REAL,
           typename T_FEATURE,
           typename T_CLUSTERIDX, //-1, 0, 1, .., K
           typename INPUT_ITERATOR
	   >
partition::PartitionDisjSets<T_CLUSTERIDX>
nearestNeighbor
(const T_REAL                    aitr_parameterU,
 INPUT_ITERATOR                  aiiterator_instfirst,
 const INPUT_ITERATOR            aiiterator_instlast,
 dist::Dist<T_REAL,T_FEATURE>    &aifunc2p_dist,
 const T_CLUSTERIDX              aii_datatypeClusterIdx = 0
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "graph::nearestNeighbor";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(const input T_PARAMETER_REAL aitr_parameterU  = "
	      << aitr_parameterU
              << "input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "\n input  dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  /*Step 1: For each object O , find the distance between
    Oi and its nearest neighbor. That is,

    d_nn(O_i)  = min_{j\ne i} ||O_j - O_i||,
  */
  T_REAL lt_dav = T_REAL(0);
  
  for (INPUT_ITERATOR literator_i = aiiterator_instfirst;
       literator_i != aiiterator_instlast; literator_i++) {
    T_REAL lt_dnn = std::numeric_limits<T_REAL>::max();
    for (INPUT_ITERATOR literator_j = aiiterator_instfirst;
       literator_j != aiiterator_instlast; literator_j++) {
      if ( literator_i != literator_j ) {
	T_REAL lt_dij = 
	  aifunc2p_dist
	  ((*literator_i)->getFeatures(),
	   (*literator_j)->getFeatures(),
	   data::Instance<T_FEATURE>::getNumDimensions() 
	   );
	if (lt_dij < lt_dnn )
	  lt_dnn = lt_dij;
      }
    }
    lt_dav += lt_dnn;
  }
 
  /*Step 2: Compute d , the average of the nearest-neighbor 
    distances by using Eq. (1) as follows:
  */
  lt_dav /= (T_REAL) uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)); 

  T_REAL  lt_d = aitr_parameterU * lt_dav; 
  /*Step 3: View the n objects as nodes of a graph. Compute 
    the adjacency matrix A_{nxn} as follows:
  */
  auto &&lbitmatrix_adjacency =
    getAdjacencyMatrix
    (lt_d,
     aiiterator_instfirst,
     aiiterator_instlast,
     aifunc2p_dist,
     (unsigned int) 0
     );

  /*Step 4:
   */
  graph::AlgorithmsDFSAdjMatrix<unsigned int>
    lalgdfsadjmatrix(lbitmatrix_adjacency);

  lalgdfsadjmatrix.resolve();
 
  ds::DisjSets  lodisjsets_componentSetBi = 
    graph::component(lalgdfsadjmatrix.getPi());
    
  partition::PartitionDisjSets<T_CLUSTERIDX>
    lomembclassdisjsets_Bi(std::move(lodisjsets_componentSetBi));
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

    std::cout <<  "lt_dav = " << lt_dav << '\n';
    std::cout <<  "lt_d = "   << lt_d   << '\n';

    std::ostringstream lostrstream_labelBi;
    lostrstream_labelBi << "<MEMBER Bi:CLASSDISJSETS_Bi" << lpc_labelFunc;
    lomembclassdisjsets_Bi.print(std::cout,lostrstream_labelBi.str().c_str(),',');
    std::cout << std::endl;

    std::vector<std::list<uintidx> > lvectorlist_tmpgraphpi =
      graph::fromPiToGraph(lalgdfsadjmatrix.getPi());
    std::ostringstream lostrstream_labelGraphPi;
    lostrstream_labelGraphPi << "<GRAPH:PI:" << lpc_labelFunc;
    inout::vectorlistprint
      (lvectorlist_tmpgraphpi,
       std::cout,lostrstream_labelGraphPi.str().c_str(),
       ';'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lomembclassdisjsets_Bi;

} /*END nearestNeighbor */

} /*END namespace graph*/

#endif /*GRAPH_UTILS_HPP*/

