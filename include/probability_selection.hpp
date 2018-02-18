/*! \file probability_selection.hpp
 *
 * \brief probability selection
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef PROBABILITY_SELECTION_HPP
#define PROBABILITY_SELECTION_HPP

#include <random>
#include <unordered_set>
#include <utility>        //std::pair
#include <algorithm>      // std::find

#include "ga_selection.hpp"
#include "verbose_global.hpp"
#include "insertion_operator.hpp"

//#include "container_out.hpp"

extern std::mt19937       gmt19937_eng;

/*! \namespace prob
  \brief functions for get randon number 
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace prob {
  
/*! \fn std::pair<uintidx,uintidx> getIdxPairUnlikeRoulette(const std::vector<T_PROBABILITY> &aivectorrt_probDist, uintidx aiuintidx_begins)
  \brief Selection a pair random index in [IdxBegin,n]
  \details Selection a pair random index proportional to its probability dstribution by aivectorrt_probDist
  \param aivectorrt_probDist a acumulate probability dstribution in [0,1]
  \param aiuintidx_begins a unsigned integer where begin number to selected
  \code{.cpp}
     const std::vector<T_REAL>&& lvectorprob_selectBySizeClusterK =
      prob::makeDistRouletteWheel
      (lvectorrt_numInstClusterK.begin(),lvectorrt_numInstClusterK.end(),
       [](const uintidx& lui_numInstClusterK) -> T_REAL
       {
	 return T_REAL(1.0) / T_REAL(lui_numInstClusterK);
       }
       );

    std::pair<uintidx,uintidx> lpair_mergeIdxK =   
      prob::getIdxPairUnlikeRoulette
      (lvectorprob_selectBySizeClusterK, 
       uintidx(0)
       );

  \endcode
 */
template < typename T_PROBABILITY>
std::pair<uintidx,uintidx>
getIdxPairUnlikeRoulette
(const std::vector<T_PROBABILITY> &aivectorrt_probDist,  
 uintidx                          aiuintidx_begins
 )
{
  uintidx loit_idxrand1;
  uintidx loit_idxrand2;

  loit_idxrand1 = 
    gaselect::getIdxRouletteWheel
    (aivectorrt_probDist,
     aiuintidx_begins
     );
  do {
    loit_idxrand2 =
      gaselect::getIdxRouletteWheel
      (aivectorrt_probDist,
       aiuintidx_begins
       );
  } while (loit_idxrand1 == loit_idxrand2); 

  return std::make_pair(loit_idxrand1,loit_idxrand2);
  
}


/*! \fn std::vector<uintidx> getIdxVectorRoulette(const std::vector<T_PROBABILITY> &aivectort_probDist, const bool aib_withRepeatedValues, const uintidx aiuintidx_numbersToGenerate, const uintidx aiuintidx_begins, FUNCTION_RAND function_rand01)
  \brief Selection a vector of random index in [IdxBegin,n]
  \details Selection a vector of random index  proportional to its probability dstribution by aivectorrt_probDist
  \param aivectort_probDist a const std::vector
  \param aib_withRepeatedValues a bool if is false get get with repeated values
  \param aiuintidx_numbersToGenerate a unsigned integer with numbers of items 
  \param aiuintidx_begins a unsigned integer where begin number to selected
  \param function_rand01 a random distribution function in [0,1]
*/
template < typename T_PROBABILITY, typename FUNCTION_RAND>
std::vector<uintidx>
getIdxVectorRoulette
(const std::vector<T_PROBABILITY>   &aivectort_probDist,  
 const bool                         aib_withRepeatedValues, 
 const uintidx                      aiuintidx_numbersToGenerate, 
 const uintidx                      aiuintidx_begins,
 FUNCTION_RAND                      function_rand01 
 )
{
  std::vector<uintidx> lovector_randomIdx;
  lovector_randomIdx.reserve(aiuintidx_numbersToGenerate);

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "prob::getIdxVectorUnlikeRoulette:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(input std::vector<T_PROBABILITY>: aivectort_probDist[" << &aivectort_probDist << ']'
	      << "\n\t input bool aib_withRepeatedValues = " << aib_withRepeatedValues
	      << "\n\t input uintidx: aiuintidx_numbersToGenerate = " << aiuintidx_numbersToGenerate
	      << "\n\t uintidx: aiuintidx_begins = " << aiuintidx_begins
	      << "\n\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  uintidx    luintidx_randDist;
  
  for (uintidx luintidx_i = 0; luintidx_i < aiuintidx_numbersToGenerate; luintidx_i++) {
    do {
      luintidx_randDist = 
	gaselect::getIdxRouletteWheel
	(aivectort_probDist,  
	 function_rand01(),
	 aiuintidx_begins
	 );
    } while 
	(!aib_withRepeatedValues 
	 && 
	 (std::find
	  (lovector_randomIdx.begin(),
	   lovector_randomIdx.end(),
	   luintidx_randDist) != lovector_randomIdx.end()) 
	 );
    lovector_randomIdx.push_back(luintidx_randDist);
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "prob::getIdxVectorUnlikeRoulette: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\nstd::vector<uintidx>: lovector_randomIdx[" << &lovector_randomIdx << "]:"
	      << lovector_randomIdx
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovector_randomIdx;

}

/*! \fn auto getRandPairUnlike(FUNCTION function_generate)
    \brief Get two random number diferent
    \param function_generate a random distribution function 
*/
template<typename FUNCTION>
auto
getRandPairUnlike
(FUNCTION function_generate) -> std::pair<decltype(function_generate()),decltype(function_generate())>
{
  typedef decltype(function_generate()) ResultType;
  ResultType lit_rand1;
  ResultType lit_rand2;

  lit_rand1 = function_generate();
  do {
    lit_rand2 = function_generate();
  } while (lit_rand1 == lit_rand2); 

  return std::make_pair(lit_rand1,lit_rand2);
}


/*! \fn auto getRandPairUnlikeInOrd(FUNCTION function_generate)
    \brief Get two random number diferent in order 
    \param function_generate a random distribution function 
*/
template<typename FUNCTION>
auto
getRandPairUnlikeInOrd
(FUNCTION function_generate) -> std::pair<decltype(function_generate()),decltype(function_generate())>
{
  typedef decltype(function_generate()) ResultType;
  ResultType lit_rand1;
  ResultType lit_rand2;

  lit_rand1 = function_generate();
  do {
    lit_rand2 = function_generate();
  } while (lit_rand1 == lit_rand2); 

  if (lit_rand2 < lit_rand1)
    std::swap(lit_rand1,lit_rand2);
  
  return std::make_pair(lit_rand1,lit_rand2);
}

/*! \fn auto getWithoutRepeatsSet(const uintidx aiuintidx_numbersToGenerate, FUNCTION function_generate)
    \brief Get a set of random number diferent
    \param aiuintidx_numbersToGenerate a unsigned integer with numbers to Generate
    \param function_generate a random distribution function 
*/
template<typename FUNCTION>
auto
getWithoutRepeatsSet
(const uintidx aiuintidx_numbersToGenerate,
 FUNCTION      function_generate
 ) ->  std::unordered_set<decltype(function_generate())>
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::getWithoutRepeatsSet";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "input aiuintidx_numbersToGenerate = " << aiuintidx_numbersToGenerate
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  typedef decltype(function_generate()) ResultType;

  std::unordered_set<ResultType> lounorderedset_rand;
  
  lounorderedset_rand.reserve(aiuintidx_numbersToGenerate);
 
  while ( lounorderedset_rand.size() <  aiuintidx_numbersToGenerate ) {
 
    lounorderedset_rand.insert( function_generate() );
    
  }
  
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_unorderedset;
      lostrstream_unorderedset << "<UNORDEREDSET:"  << lpc_labelFunc;

    inout::containerprint
      (lounorderedset_rand.begin(),
       lounorderedset_rand.end(),
       std::cout,
       lostrstream_unorderedset.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lounorderedset_rand;
}

/*! \fn auto getRandUnlike(const T_RAND  at_unlike, FUNCTION function_generate)

    \brief A random different than a given at_unlike
    \param at_unlike a number
    \param function_generate a random distribution function 

    \code{.cpp}

    // 3.7 Island migration

    std::uniform_int_distribution<uintidx> uniformdis_ui0NumIsland
    (0, aiinparam_GAProbAdapRangeK.getNumIsland()-1 );

    for (uintidx lui_l = 0; lui_l < lvectorchrom_bestIsland.size(); lui_l++) {
	
       if ( uniformdis_real01(gmt19937_eng) < aiinparam_GAProbAdapRangeK.getPe() ) {
	  
	  //2. Randomly choose the island toward each individual 
	  //  will migrate.
	  
	  uintidx lui_islandToward = 
	    prob::getRandUnlike
	    (lui_l,
	     [&]() 
	     {
	       return uniformdis_ui0NumIsland(gmt19937_eng);
	     }
	     );
	  
	  
	  std::vector<ChromAgustin2012<T_CLUSTERIDX,T_REAL>* >*
	    lvectorchrom_subpopulation = lvectorvector_subPopulation.at(lui_islandToward);
	  
	  //3. Randomly choose an individual in the destiny island 
	  //  and change it by the migrating individual.
	  
	  uintidx lui_idxIndividual =
	    uniformdis_ui0SubPopSize(gmt19937_eng);
	  
	  *lvectorchrom_subpopulation->at(lui_idxIndividual) =
	    *lvectorchrom_bestIsland[lui_l];
       }
    }

    \endcode
    \cite Agustin:etal:GAclusteringVarK:GGA:2012

*/
template < typename T_RAND, typename FUNCTION >
auto
getRandUnlike
(const T_RAND  at_unlike,
 FUNCTION      function_generate
 ) ->  decltype(function_generate())
{
  typedef decltype(function_generate()) ResultType;
  
  ResultType lit_rand;
  do {
    lit_rand = function_generate();
  } while (lit_rand == at_unlike); 

  return lit_rand;
}


/*! \fn T_RAND getRandSetUnlike(const std::unordered_set<T_RAND> &aiunorderedset_rand, FUNCTION function_generate)
  \brief Selection a random un like to set items
  \details this function is used for implementate a point mutation operator, D_PM \cite Lucasius:etal:GAclusteringMedoid:GCA:1993 
  \param aiunorderedset_rand IN a std::unordered_set 
  \param function_generate a random distribution function 
  
  \code{.cpp}
  
  for (auto lchromfixleng_iter :lvectorchromfixleng_population) {

	if ( uniformdis_real01(gmt19937_eng) 
	     < aiipcgapc_inParamGALucasius.getProbMutation() ) 
	  {	 

	    std::unordered_set<uintidx> lunorderedset_medoids
	      (lchromfixleng_iter->begin(),
	       lchromfixleng_iter->end()
	       );

	    uintidx luintidx_medoidRand =
	      prob::getRandSetUnlike
	      (lunorderedset_medoids,
	       [&]()
	       {
		 return uniformdis_idxInstance(gmt19937_eng);
	       }
	       );

	    uintidx lui_positionGene = uniformdis_uiMutation0N(gmt19937_eng); 
     
	    lchromfixleng_iter->setGene(lui_positionGene,luintidx_medoidRand);
	    
	  }
      }

  }
  \endcode

*/
template<typename T_RAND, typename FUNCTION>
T_RAND
getRandSetUnlike
(const std::unordered_set<T_RAND> &aiunorderedset_rand,
 FUNCTION                         function_generate
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::getRandSetUnlike";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output const std::unordered_set<>: &aiunorderedset_rand["  
	      << &aiunorderedset_rand << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES

  T_RAND lt_unlikeRand;
  lt_unlikeRand = function_generate();
  auto liter_searchMedoid = aiunorderedset_rand.find(lt_unlikeRand);
  while ( liter_searchMedoid  != aiunorderedset_rand.end() ) {
    lt_unlikeRand = function_generate();
    liter_searchMedoid = aiunorderedset_rand.find(lt_unlikeRand);
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lt_unlikeRand = " << lt_unlikeRand
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lt_unlikeRand;
}
  

/*! \fn std::pair<uintidx,uintidx>  getElitistPairCross(std::vector<uintidx> *aoptvectoruintidx_chromPairsCross, uintidx aiuintidx_idxChromTMinIdxSelect, uintidx aiuintidx_idxChromTMaxIdxSelect)
  \brief Selection a pair random index
  \details
  \param aoptvectoruintidx_chromPairsCross a std::vector<uintidx>
  \param aiuintidx_idxChromTMinIdxSelect a unsigned integer
  \param aiuintidx_idxChromTMaxIdxSelect a unsigned integer
  \cite Franti:etal:GAclustering:gafranti:1997
 */
std::pair<uintidx,uintidx>  
getElitistPairCross
(std::vector<uintidx> *aoptvectoruintidx_chromPairsCross,
 uintidx               aiuintidx_idxChromTMinIdxSelect,
 uintidx               aiuintidx_idxChromTMaxIdxSelect 
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "prob::getElitistPairCross: IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t(output std::vector<uintidx>: aoptvectoruintidx_chromPairsCross[" 
	      << aoptvectoruintidx_chromPairsCross << "]\n"
	      << "\t input  uintidx: aiuintidx_idxChromTMinIdxSelect = " 
	      << aiuintidx_idxChromTMinIdxSelect << "]\n"
	      << "\t input  uintidx: aiuintidx_idxChromTMaxIdxSelect = " 
	      << aiuintidx_idxChromTMaxIdxSelect << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  if ( aoptvectoruintidx_chromPairsCross->size() == 0)  {
    aoptvectoruintidx_chromPairsCross->push_back(aiuintidx_idxChromTMinIdxSelect);
    aoptvectoruintidx_chromPairsCross->push_back(aiuintidx_idxChromTMinIdxSelect+1);
  }
  else  {
    (*aoptvectoruintidx_chromPairsCross)[1] = (*aoptvectoruintidx_chromPairsCross)[1] + 1;
    if ( (*aoptvectoruintidx_chromPairsCross)[1] ==  aiuintidx_idxChromTMaxIdxSelect ) {
      (*aoptvectoruintidx_chromPairsCross)[0] = (*aoptvectoruintidx_chromPairsCross)[0] + 1;
      (*aoptvectoruintidx_chromPairsCross)[1] = (*aoptvectoruintidx_chromPairsCross)[0] + 1;
      
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "prob::getElitistPairCross: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output pairCross: " 
	      << '(' << (*aoptvectoruintidx_chromPairsCross)[0]  
	      << ',' << (*aoptvectoruintidx_chromPairsCross)[1] 
	      << ")\n";
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_pair((*aoptvectoruintidx_chromPairsCross)[0],(*aoptvectoruintidx_chromPairsCross)[1]);

}

}/*END namespace prob*/
  
#endif /*PROBABILITY_SELECTION_HPP*/
