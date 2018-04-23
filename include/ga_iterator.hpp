/*! \file ga_iterator.hpp
 *
 * \brief iterators for applying crossover operators
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef GA_ITERATOR_HPP
#define GA_ITERATOR_HPP

#include "probability_selection.hpp"

extern StdMT19937 gmt19937_eng;


/*! \namespace gaiterator
  \brief Iterators for applying crossover operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace gaiterator {
  
/*! \fn void crossover(INPUT_ITERATOR aiiterator_instfirstParent, const INPUT_ITERATOR aiiterator_instlastParent, INPUT_ITERATOR aiiterator_instfirstChild, const INPUT_ITERATOR aiiterator_instlastChild, const GENETIC_OPERATOR genetic_operator) 
    \brief Pairs crossover iterator 
    \details Function to iterate over two containers, selected pairs consecutively to apply a crossover operator
    \param aiiterator_instfirstParent an iterator
    \param aiiterator_instlastParent  a const iterator
    \param aiiterator_instfirstChild  an iterator
    \param aiiterator_instlastChild   a const iterator
    \param genetic_operator a crossover function  

    \code{.cpp}

    gaiterator::crossover
    (lvec_matingPool.begin(),
     lvec_matingPool.end(),
     lvec_population.begin(),
     lvec_population.end(),
     [&](const gaencode::ChromFixedLength<double,double>& aichrom_parent1,
	 const gaencode::ChromFixedLength<double,double>& aichrom_parent2,
	 gaencode::ChromFixedLength<double,double>&  aochrom_child1, 
	 gaencode::ChromFixedLength<double,double>&  aochrom_child2
	 ) -> bool
     {

       bool lb_appliesCrossover =
	 (uniformdis_real01(gmt19937_eng) < lr_crossoverProbability);
	 
       if ( lb_appliesCrossover  ) {
		     
	   gatplop::onePointCrossover
	     (aochrom_child1,
	      aochrom_child2,
	      aichrom_parent1,
	      aichrom_parent2,
	      uniformdis_uiCrossover1l_1(gmt19937_eng)
	      );
	}// End IF Crossover
       
       return lb_appliesCrossover; 
     }
     );
 */
template<typename INPUT_ITERATOR, typename GENETIC_OPERATOR>
void
crossover
(INPUT_ITERATOR         aiiterator_instfirstParent,
 const INPUT_ITERATOR   aiiterator_instlastParent,
 INPUT_ITERATOR         aiiterator_instfirstChild,
 const INPUT_ITERATOR   aiiterator_instlastChild,
 const GENETIC_OPERATOR genetic_operator
 )
{

  const uintidx lui_sizePopulation
    (uintidx(std::distance(aiiterator_instfirstParent,aiiterator_instlastParent)));

  if ( ( lui_sizePopulation % 2 ) != 0 ) {
    *aiiterator_instfirstChild = *aiiterator_instfirstParent;
    ++aiiterator_instfirstChild;
    ++aiiterator_instfirstParent;
  }
  while ( (aiiterator_instfirstChild != aiiterator_instlastChild )
	  && (aiiterator_instfirstParent != aiiterator_instlastParent) )
    {
	     
      INPUT_ITERATOR lchrom_parent1  = aiiterator_instfirstParent; 
      ++aiiterator_instfirstParent;
      INPUT_ITERATOR lchrom_parent2  = aiiterator_instfirstParent; 
      ++aiiterator_instfirstParent;

      INPUT_ITERATOR lchrom_child1  = aiiterator_instfirstChild; 
      ++aiiterator_instfirstChild;
      INPUT_ITERATOR lchrom_child2  = aiiterator_instfirstChild; 
      ++aiiterator_instfirstChild;

      genetic_operator
	(*lchrom_parent1,
	 *lchrom_parent2,
	 *lchrom_child1,
	 *lchrom_child2
	 );
	   
    }

}


/*! \fn void crossoverRandSelect(INPUT_ITERATOR aiiterator_instfirstParent, const INPUT_ITERATOR aiiterator_instlastParent, INPUT_ITERATOR aiiterator_instfirstChild, const INPUT_ITERATOR aiiterator_instlastChild, const GENETIC_OPERATOR genetic_operator) 
    \brief Random crossover iterator 
    \details Function to iterate over two containers, selected pairs randomly to apply a crossover operator
    \param aiiterator_instfirstParent an iterator
    \param aiiterator_instlastParent  a const iterator
    \param aiiterator_instfirstChild  an iterator
    \param aiiterator_instlastChild   a const iterator
    \param genetic_operator a crossover function 
*/
template<typename INPUT_ITERATOR, typename GENETIC_OPERATOR>
void
crossoverRandSelect
(const INPUT_ITERATOR   aiiterator_instfirstParent,
 const INPUT_ITERATOR   aiiterator_instlastParent,
 INPUT_ITERATOR         aiiterator_instfirstChild,
 const INPUT_ITERATOR   aiiterator_instlastChild,
 const GENETIC_OPERATOR genetic_operator
 )
{

  const uintidx lui_sizePopulation
    (uintidx(std::distance(aiiterator_instfirstParent,aiiterator_instlastParent)));
  
  std::uniform_int_distribution<uintidx> uniformdis_uintidx0PoolSize
    (0,lui_sizePopulation-1);

  if ( ( lui_sizePopulation % 2 ) != 0 ) {
    auto lchrom_parent1 =
      std::next(aiiterator_instfirstParent,uniformdis_uintidx0PoolSize(gmt19937_eng));
    *aiiterator_instfirstChild = *lchrom_parent1;
    ++aiiterator_instfirstChild;
  }
  while ( (aiiterator_instfirstChild != aiiterator_instlastChild )
	  && (aiiterator_instfirstParent != aiiterator_instlastParent) ) {

    std::pair<uintidx,uintidx>
      lpair_idxChrom =
      prob::getRandPairUnlike
      ([&]() -> uintidx
       {
	 return uniformdis_uintidx0PoolSize(gmt19937_eng);
       }
       );
	     
    auto  lchrom_parent1 = std::next(aiiterator_instfirstParent,lpair_idxChrom.first);
    auto  lchrom_parent2 = std::next(aiiterator_instfirstParent,lpair_idxChrom.second);

    INPUT_ITERATOR lchrom_child1  = aiiterator_instfirstChild;
    ++aiiterator_instfirstChild;
    INPUT_ITERATOR lchrom_child2  = aiiterator_instfirstChild;
    ++aiiterator_instfirstChild;

    genetic_operator
      (*lchrom_parent1,
       *lchrom_parent2,
       *lchrom_child1,
       *lchrom_child2
       );
    
  }

}

}  /*END namespace gaiterator*/

#endif /*GA_ITERATOR_HPP*/
