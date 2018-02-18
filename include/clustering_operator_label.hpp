/*! \file clustering_operator_label.hpp
 *
 * \brief  clustering operator label
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CLUSTERING_OPERATOR_LABEL_HPP
#define __CLUSTERING_OPERATOR_LABEL_HPP

#include <vector>
#include "partition_disjsets.hpp"

/*! \namespace clusteringop
  \brief Clustering operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace clusteringop {

/*The second stage is a genetic algorithm, which will
  merge some of these B_i's if they are close enough to one
  another. 

*/
template < typename T_CLUSTERIDX >
std::vector<T_CLUSTERIDX>
bkinCiToMemberShip
(const partition::PartitionDisjSets
 <T_CLUSTERIDX>              &aimembclassdisjsets_Bi,
 const partition::PartitionDisjSets
 <T_CLUSTERIDX>              &aimembclassdisjsets_BkinCi
 )
{
  /*#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "tsengyang2001_BkinCiToMemberShip:  IN"
    << '(' << geiinparam_verbose << ')'
    << "\n\t(input const partition::PartitionDisjSets<> &aimembclassdisjsets_Bi["
    << &aimembclassdisjsets_Bi << ']' << " size: " << aimembclassdisjsets_Bi.size() << '\n';
    aimembclassdisjsets_Bi.print();
    std::cout << "\n\t input const partition::PartitionDisjSets<> &aimembclassdisjsets_BkinCi[" 
    << &aimembclassdisjsets_BkinCi << ']' << " size: " << aimembclassdisjsets_BkinCi.size() << '\n';
    aimembclassdisjsets_BkinCi.print();
    std::cout << "\n\t)"
    << std::endl;
    }
    #endif //__VERBOSE_YES*/
  
  uintidx luintidx_numObj =
    aimembclassdisjsets_Bi.size();
  std::vector<T_CLUSTERIDX>
    lovectormmidx_memberShip
    (luintidx_numObj);
  

  // std::cout << "hola1" << std::endl;
  for (uintidx lui_i = 0; lui_i < luintidx_numObj; lui_i++) {
    lovectormmidx_memberShip.at(lui_i) =
      aimembclassdisjsets_BkinCi.getClusterIdx(aimembclassdisjsets_Bi.getClusterIdx(lui_i));
    //aivectormmidx_memberBkinCi[aimembclassdisjsets_Bi.getClusterIdx(lui_i)];

    //std::cout << aimembclassdisjsets_Bi.getClusterIdx(luintidx_numObj) << std::endl;
    //aimembclassdisjsets_Bi.constfind(lui_i)
  }

  //std::cout << "hola2"  <<std::endl;

  /*#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "tsengyang2001_BkinCiToMemberShip: OUT"
    << '(' << geiinparam_verbose << ')'
    << "\nstd::vector<T_CLUSTERIDX> lovectormmidx_memberShip["
    << &lovectormmidx_memberShip << "]\n";
    vectort_print(lovectormmidx_memberShip,std::cout,"<MEMBERCLUSTER",',');
    //<< "\noutput std::vector<std::list<uintidx> >&& lvectorlistcomponent_setBi[" 
    //	      << &lvectorlistcomponent_setBi << "]\n";
    }
    --geiinparam_verbose;
    #endif __VERBOSE_YES*/
  
  return lovectormmidx_memberShip;
}

} /*END namespace clusteringop*/
  
#endif /*__CLUSTERING_OPERATOR_LABEL_HPP*/
