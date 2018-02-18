/*! \file ga_selection.hpp
 *
 * \brief chromosome selection method in genetic algorithms
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef GASELECTION_HPP
#define GASELECTION_HPP

#include <random>
#include <vector>
#include <iterator>
#include "common.hpp" //uintidx

#include "verbose_global.hpp"

extern std::mt19937       gmt19937_eng;

/*! \namespace gaselect
  \brief Chromosome selection method in genetic algorithms
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaselect {


/*! \fn  T_INTIDX getIdxRouletteWheel(const std::vector<T_PROBABILITY> &aivectorrt_probDist, const T_INTIDX aiuintidx_begin)
  \brief Selection a random index in [IdxBegin,n] 
  \details Selection a which is proportional to its probability dstribution by aivectorrt_probDist, utilizes for roulette wheel selection method
  \param aivectorrt_probDist a acumulate probability dstribution in [0,1]
  \param aiuintidx_begin a unsigned integer index where begin index to select 

  \code{.cpp}

  //COPY POPULATION TO MATING POOL FOR ROULETTE WHEEL--------------------------

  for ( auto&& lchrom_iter: lvectorchrom_matingPool) {

     uintidx luiidx_chrom = 
     gaselect::getIdxRouletteWheel
     (lvectorT_probDistRouletteWheel,
      uintidx(0)
     );

     lchrom_iter = lvectorchrom_population.at(luiidx_chrom);

  }
  \endcode

*/
template < typename T_PROBABILITY, typename T_INTIDX >
T_INTIDX
getIdxRouletteWheel
(const std::vector<T_PROBABILITY> &aivectorrt_probDist,
 const T_INTIDX                   aiuintidx_begin 
 )
{
  /* std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.1);
  */
  static std::uniform_real_distribution<T_PROBABILITY> lsuniformdis_real01(0.0,1.0);
  T_PROBABILITY lrt_numberRand01 = lsuniformdis_real01(gmt19937_eng);
  
  T_INTIDX  loiidx_rand = aiuintidx_begin; 
  T_INTIDX  liidx_sizeVectorProbDist = (T_INTIDX) aivectorrt_probDist.size();

  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaselect::getIdxRouletteWheel";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input  vector<T_PROBABILITY>&: aivectorrt_probDist[" << &aivectorrt_probDist << "]\n"
	      << " input  aiuintidx_begin = " << aiuintidx_begin << "\n"
	      << " T_PROBABILITY  lrt_numberRand01 = " << lrt_numberRand01
	      << "\n)\n";
  }
#endif /*__VERBOSE_YES*/

  while ( loiidx_rand < liidx_sizeVectorProbDist &&
	  lrt_numberRand01 > aivectorrt_probDist[loiidx_rand]
	  )
    loiidx_rand++;
  if ( loiidx_rand >= liidx_sizeVectorProbDist ) --loiidx_rand;
  loiidx_rand += aiuintidx_begin;
 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << ", T_INTIDX loiidx_rand = " << loiidx_rand
	      << std::endl;		
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return loiidx_rand;

}


/*! \fn INPUT_ITERATOR tournament(const INPUT_ITERATOR aiiterator_instfirst, const uintidx aiuintidx_orderTournament, const FITNESSOPERATION fitness_func)
  \brief Selection a item of container for method tournament
  \details
  \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of chromosome
  \param aiiterator_instlast a InputIterator to the final positions of the sequence of chromosome
  \param aiuintidx_orderTournament a unsigned integer the tournament order is set to he as low as possible, i.e. 2.
  \param fitness_func a function of fitness individuals
  \note use by
  \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004

  \code{.cpp}
     for (uintidx luintidx_i = 0; 
	   luintidx_i < lvectorchromfixleng_population.size()/2; 
	   luintidx_i++) {
      
	ChromFixedLength<uintidx,T_REAL>* lchrom_selectTour = 
	  *(gaselect::tournament
	    (lvectorchromfixleng_population.begin(),
	     aiipcgapc_inParamHKA.getOrderTournament(),
	     [&](const ChromFixedLength<uintidx,T_REAL>* iter_chrom)
	     {
	        return iter_chrom->getFitness();
	     }
	     )
	    );
	lvectorchromfixleng_matingPool.push_back(lchrom_selectTour);
      }
  \endcode
*/
template<typename INPUT_ITERATOR, typename FITNESSOPERATION>
INPUT_ITERATOR
tournament
(const INPUT_ITERATOR      aiiterator_instfirst,
 const INPUT_ITERATOR      aiiterator_instlast,
 const uintidx             aiuintidx_orderTournament,
 const FITNESSOPERATION    fitness_func
 )
{    
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::selectTournament";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
              << "(input aiiterator_instfirst[" << &aiiterator_instfirst << "]\n"
              << " input aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << " input  uintidx: aiuintidx_orderTournament = " << aiuintidx_orderTournament
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));
  
  std::uniform_int_distribution<uintidx> luniformdis_populationSize
	(0, (uintidx) (lui_numInstances-1));
  
  INPUT_ITERATOR loptchom_best = aiiterator_instfirst;
  std::advance(loptchom_best,luniformdis_populationSize(gmt19937_eng));
   
  for (uintidx li_i = 1; li_i < aiuintidx_orderTournament; li_i++) {

    INPUT_ITERATOR literator_select = aiiterator_instfirst;
    std::advance(literator_select,luniformdis_populationSize(gmt19937_eng));
  
    if ( fitness_func(*literator_select)  > fitness_func(*loptchom_best)  )
      loptchom_best = literator_select;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    (*loptchom_best)->print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return loptchom_best;
  
}

  
} /*END namespace stdgaselect*/

#endif /*GASELECTION_HPP*/
