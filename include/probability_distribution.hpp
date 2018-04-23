/*! \file probability_distribution.hpp
 *
 * \brief probability distribution
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef PROBABILITY_DISTRIBUTION_HPP
#define PROBABILITY_DISTRIBUTION_HPP

#include <iostream>
#include <list>
#include <utility>       //std::pair
#include <iterator>     // std::distance
#include <type_traits>
#include <typeinfo>
#include <random>
#include "probability_selection.hpp"
#include "vector_utils.hpp"
#include "container_out.hpp"
#include "linear_algebra_level1.hpp"
#include "leac_utils.hpp"

#include "verbose_global.hpp"

extern  StdMT19937 gmt19937_eng;

/*! \namespace prob
  \brief Method for calculate probability distributions
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace prob {

  
/*! \fn auto makeDistRouletteWheel(INPUT_ITERATOR aiiterator_instfirst, INPUT_ITERATOR aiiterator_instlast, const FITNESSOPERATION fitness_func)
  \brief
  \details 
  \param aiiterator_instfirst
  \param aiiterator_instlast
  \param fitness_func

  \code{.cpp}

  const std::vector<T_REAL>&& lvectorT_probDistRouletteWheel =
  prob::makeDistRouletteWheel
    (lvectorchromfixleng_population.begin(),
     lvectorchromfixleng_population.end(),
     [](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* lchromfixleng_iter) -> T_REAL
     {
         return lchromfixleng_iter->getFitness();
     }
    );

  \endcode

*/
template<typename INPUT_ITERATOR, typename FITNESSOPERATION>
auto
makeDistRouletteWheel
(INPUT_ITERATOR     aiiterator_instfirst,
 INPUT_ITERATOR     aiiterator_instlast,
 const FITNESSOPERATION fitness_func
 ) -> std::vector<decltype(fitness_func(*aiiterator_instfirst))>
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::makeDistRouletteWheel";
  ++geiinparam_verbose;
  auto literator_firstcopy =  aiiterator_instfirst;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labellinealNorm;
     lostrstream_labellinealNorm << "<FITNESS:"
				 << lpc_labelFunc;
     inout::containerprint
      (literator_firstcopy,
       aiiterator_instlast,
       fitness_func,
       std::cout,
       lostrstream_labellinealNorm.str().c_str(),
       ','
       );
     std::cout <<  "\n)"
	       << std::endl;
  }
#endif //__VERBOSE_YES

  typedef decltype(fitness_func(*aiiterator_instfirst)) ResultType;
  const uintidx lui_sizeProbDist(uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)));
  std::vector<ResultType> lovectorrt_probRouletteWheel;
  
  lovectorrt_probRouletteWheel.reserve(lui_sizeProbDist);

  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    lovectorrt_probRouletteWheel.push_back(fitness_func(*aiiterator_instfirst));
  }

  ResultType
    lt_sumRandVar =
    interfacesse::sum
    (lovectorrt_probRouletteWheel.data(),
     lui_sizeProbDist
     );

  if ( lt_sumRandVar > 0 ) {
     ResultType
     lt_scale = 1.0 / lt_sumRandVar;
    interfacesse::scal
      (lovectorrt_probRouletteWheel.data(),
       lt_scale,
       lui_sizeProbDist
       );
    for (uintidx li_i = 1; li_i < lui_sizeProbDist; li_i++) {
      lovectorrt_probRouletteWheel[li_i] += lovectorrt_probRouletteWheel[li_i-1];	
    }
  }
  else {
    ResultType
      lt_probProportional
      = 1.0 / (ResultType) lui_sizeProbDist; 
    interfacesse::copya
      (lovectorrt_probRouletteWheel.data(),
       lt_probProportional,
       lui_sizeProbDist
       );
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ") "
	      << "lt_sumRandVar =  " << lt_sumRandVar << '\n';
    inout::containerprint
      (lovectorrt_probRouletteWheel.begin(),
       lovectorrt_probRouletteWheel.end(),
       std::cout,
       lpc_labelFunc,
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lovectorrt_probRouletteWheel;
  
}
  

/*linearNormalization:
  linear normalization
  \cite{Alves:etal:GAclusteringLabelKVar:FEAC:2006}
  James et al 2007
  \cite{Agustin:etal:GAclusteringVarK:GGA:2012} 
*/

/*! \fn void linearNormalization(INPUT_ITERATOR aiiterator_first, const INPUT_ITERATOR aiiterator_last, const FITNESSFUNCTGET  fitness_funcGet, const FITNESSFUNCTSET  fitness_funcSet, const decltype(fitness_funcGet(*aiiterator_first)) lrt_scale = 1.0) 
  \brief Linear normalization
  \details
  \param aiiterator_first a Iterator vector
  \param aiiterator_last a  Iterator vector
  \param fitness_funcGet a  Fitness funtion get
  \param fitness_funcSet a  Fitness funtion a set
  \param lrt_scale a scaling factor the value by default  is 1.0
  \code{.cpp}
  
  prob::linearNormalization
    (lvectorchrom_population.begin(),
    lvectorchrom_population.end(),
    [](const ChromosomeFEAC<T_CLUSTERIDX,T_REAL,T_FEATURE,T_INSTANCES_CLUSTER_K>
    &lchromfeac_iter) -> T_REAL
    {
        return lchromfeac_iter.getObjetiveFunc();
    },
    [](ChromosomeFEAC<T_CLUSTERIDX,T_REAL,T_FEATURE,T_INSTANCES_CLUSTER_K>
    &lchromfeac_iter, T_REAL airt_funcFitnessLineal)
    {
        lchromfeac_iter.setFitness(airt_funcFitnessLineal);
    },
    T_REAL(1.0)
   );

  \endcode
*/
template<typename INPUT_ITERATOR, typename FITNESSFUNCTGET, typename FITNESSFUNCTSET>
void 
linearNormalization
(INPUT_ITERATOR          aiiterator_first,
 const INPUT_ITERATOR    aiiterator_last,
 const FITNESSFUNCTGET  fitness_funcGet, 
 const FITNESSFUNCTSET  fitness_funcSet, 
 const decltype(fitness_funcGet(*aiiterator_first)) lrt_scale = 1.0
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::linearNormalization";
  ++geiinparam_verbose;
  auto literator_firstcopy =  aiiterator_first; 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  typedef decltype(fitness_funcGet(*aiiterator_first)) ResultType;
  const uintidx lui_sizeProbDist(uintidx(std::distance(aiiterator_first,aiiterator_last)));
  const auto literator_firstbegin =  aiiterator_first; 
  
  std::vector<std::pair<uintidx,ResultType> > lvectorT_tmpPairIdxFuncVar;
  lvectorT_tmpPairIdxFuncVar.reserve(lui_sizeProbDist);

  uintidx luintidx_idxChrom = 0;
  for (; aiiterator_first != aiiterator_last; ++aiiterator_first) {
    lvectorT_tmpPairIdxFuncVar.emplace_back(luintidx_idxChrom++,fitness_funcGet(*aiiterator_first));
  }

  std::sort
    (lvectorT_tmpPairIdxFuncVar.begin(),
     lvectorT_tmpPairIdxFuncVar.end(), 
     [](const std::pair<uintidx,ResultType> &left, const std::pair<uintidx,ResultType> &right) 
     {
       return left.second > right.second; 
     }
     );

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "SORT FUNTION FITNESS:  IN"
	      << '(' << geiinparam_verbose << ')';
    std::cout << "\nSORT FUNCTION FITNESS: OUT"
	      << '(' << geiinparam_verbose << ")\n";
      for (const auto& pairIdxFuncFitness: lvectorT_tmpPairIdxFuncVar)
	{
	  std::cout << '(' <<pairIdxFuncFitness.first << "," << pairIdxFuncFitness.second  << ')';
	  
	}
      std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES   

  ResultType lt_probIni = ResultType(lui_sizeProbDist);
  for (const auto& pairIdxFuncFitness: lvectorT_tmpPairIdxFuncVar) {
    auto literator_inc =  literator_firstbegin;
    std::advance(literator_inc,pairIdxFuncFitness.first);
    fitness_funcSet(*literator_inc,lt_probIni * lrt_scale);
      lt_probIni -= 1.0;
  }
 

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

     std::ostringstream lostrstream_labellinealNorm;
     lostrstream_labellinealNorm << "<DISTRIBUTIONLINEALNORM:"
				 << lpc_labelFunc;
     inout::containerprint
      (literator_firstcopy,
       aiiterator_last,
       fitness_funcGet,
       std::cout,
       lostrstream_labellinealNorm.str().c_str(),
       ','
       );
     std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

}


/*! \fn auto linearNormalization(INPUT_ITERATOR aiiterator_first, const INPUT_ITERATOR aiiterator_last, const FITNESSOPERATION fitness_func, const decltype(fitness_func(*aiiterator_first)) lrt_scale = 1.0)
    \brief 
    \details  The values of the variable store in a container are normalized according to a linear normalization procedure \cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
  \param aiiterator_first a INPUT_ITERATOR
  \param aiiterator_last  a INPUT_ITERATOR
  \param fitness_func a fintness function 
  \param lrt_scale a number scale default 1.0

  \code{.cpp}

  std::vector<T_REAL>&& lvectorT_partialFcCLinearNorm =
  linearNormalization
  (lvectorrt_genotypePartialFcC.begin(),
   lvectorrt_genotypePartialFcC.end(),
   [](const T_REAL& lt_metrici) -> T_REAL
   {
     return lt_metrici;
   }
  );

  \endcode
*/
template<typename INPUT_ITERATOR, typename FITNESSOPERATION>
auto
linearNormalization
(INPUT_ITERATOR             aiiterator_first,
 const INPUT_ITERATOR       aiiterator_last,
 const FITNESSOPERATION    fitness_func,
 const decltype(fitness_func(*aiiterator_first)) lrt_scale = 1.0
 ) -> std::vector<decltype(fitness_func(*aiiterator_first))>
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::linearNormalization";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  typedef decltype(fitness_func(*aiiterator_first)) ResultType;
  const uintidx lui_sizeProbDist(uintidx(std::distance(aiiterator_first,aiiterator_last)));
  
  std::vector<ResultType> lovectorT_linealNorm(lui_sizeProbDist);
  
  std::vector<std::pair<uintidx,ResultType> > lvectorT_tmpPairIdxFuncVar;
  lvectorT_tmpPairIdxFuncVar.reserve(lui_sizeProbDist);

  uintidx luintidx_idxChrom = 0;
  for (; aiiterator_first != aiiterator_last; ++aiiterator_first) {
    lvectorT_tmpPairIdxFuncVar.emplace_back(luintidx_idxChrom++,fitness_func(*aiiterator_first));
  }

  std::sort
    (lvectorT_tmpPairIdxFuncVar.begin(),
     lvectorT_tmpPairIdxFuncVar.end(), 
     [](const std::pair<uintidx,ResultType> &left, const std::pair<uintidx,ResultType> &right) 
     {
       return left.second > right.second; //descend
     }
     );

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "SORT FUNTION FITNESS:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << '\n';
    for (const auto& pairIdxFuncFitness: lvectorT_tmpPairIdxFuncVar)
      {
	std::cout << '(' <<pairIdxFuncFitness.first << "," << pairIdxFuncFitness.second  << ')';
      }
    std::cout << "\nSORT FUNCTION FITNESS: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES   

  ResultType lt_probIni = ResultType(lui_sizeProbDist); 
  for (uintidx luintidx_i = 0; luintidx_i < lui_sizeProbDist; luintidx_i++) {
    lovectorT_linealNorm
      [ lvectorT_tmpPairIdxFuncVar.at(luintidx_i).first ] = lt_probIni * lrt_scale; 
    lt_probIni -= 1.0;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

     std::ostringstream lostrstream_labellinealNorm;
     lostrstream_labellinealNorm << "<DISTRIBUTIONLINEALNORM:"
				 << lpc_labelFunc;
     inout::containerprint
      (lovectorT_linealNorm.begin(),
       lovectorT_linealNorm.end(),
       std::cout,
       lostrstream_labellinealNorm.str().c_str(),
       ','
       );
     std::cout << std::endl;
     
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectorT_linealNorm;  
}


/*! \fn std::vector<T_REAL>  getDistGKA(const std::vector<T_REAL> &aivector_randVar, T_REAL airt_cm) 
  \brief GKA Probability Distribution \cite Krishna:Murty:GAClustering:GKA:1999
  \details Obtains a probability distribution based on the Euclidean distance for the gene mutation
  \param aivector_randVar a vector with the values of a random variable 
  \param airt_cm a constant usually \f$\geq\f$
*/
template < typename T_REAL >
std::vector<T_REAL>  
getDistGKA
(const std::vector<T_REAL> &aivector_randVar,
 const T_REAL              airt_cm
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::getDistGKA";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "( input const std::vector<T_REAL>: aivector_randVar["
	      << &aivector_randVar << "]\n"
	      << " T_REAL: airt_cm = " << airt_cm
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES */

 
  const uintidx lui_sizeProbDist = (uintidx) aivector_randVar.size();
 
  T_REAL  lrt_dmaxRandVar =
    *std::max_element(aivector_randVar.begin(),aivector_randVar.end());
  std::vector<T_REAL> lovectorrt_probDistGKA
    (lui_sizeProbDist, airt_cm * lrt_dmaxRandVar);

  interfacesse::axpy
    (lovectorrt_probDistGKA.data(),
     -1.0,
     aivector_randVar.data(),
     lui_sizeProbDist
     );
    
  T_REAL
    lrt_sumRandVar =
    interfacesse::sum
    (lovectorrt_probDistGKA.data(),
     lui_sizeProbDist
     );


  interfacesse::scal
    (lovectorrt_probDistGKA.data(),
     1.0 /  lrt_sumRandVar,
     lui_sizeProbDist
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")"
	      << " lt_sumRandVar =  " << lrt_sumRandVar
	      << " lrt_dmaxRandVar =  " << lrt_dmaxRandVar
	      << '\n';
    inout::containerprint
      (lovectorrt_probDistGKA.begin(),
       lovectorrt_probDistGKA.end(),
       std::cout,
       lpc_labelFunc,
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lovectorrt_probDistGKA;
  
} /*makeDistGKA*/



/*adaptiveProb:

  PROPOSED:
  \cite{srinivas:etal:94}

  \cite{Chang:etal:GAclustering:GAGRClustering:2009}
  GAGR clustering algorithm
  The crossover probability is selected adaptively as in ref \cite{srinivas:etal:94}
 
  POSTCONDITION:
  Return crossover probability is selected adaptively
*/

/*! \fn T_REAL adaptiveProb(T_REAL aiT_k1, T_REAL aiT_k3, T_REAL aiT_fmax, T_REAL aiT_faverage, T_REAL aiT_fprime) 
  \brief Adaptively probability
  \details 
  \param aiT_k1 a real number
  \param aiT_k3 a real number
  \param aiT_fmax a real number
  \param aiT_faverage a real number
  \param aiT_fprime a real number
 )
*/
template < typename T_REAL >
T_REAL
adaptiveProb
(T_REAL aiT_k1,
 T_REAL aiT_k3,
 T_REAL aiT_fmax,
 T_REAL aiT_faverage,
 T_REAL aiT_fprime
 )
{
  T_REAL  lort_pc;     /*Probability adaptively*/ 

  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "prob::adaptiveProb";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc << ':' << geverbosepc_labelstep 
              << ":  IN(" << geiinparam_verbose << ')';
    printf(" aiT_k1 = %g aiT_k3 = %g aiT_fmax = %2.30lf aiT_faverage = %2.30lf aiT_fprime = %2.30lf\n)",
	   aiT_k1,aiT_k3,aiT_fmax,aiT_faverage,aiT_fprime);
    std::cout << std::endl;
  }
#endif /*__VERBOSE_YES*/

  bool lb_eqFprimeFaverage =
    utils::closeenough(aiT_fprime,aiT_faverage,1e-8);

  if (!lb_eqFprimeFaverage && (aiT_fprime >  aiT_faverage)) {
    lort_pc = 
      aiT_k1 * (aiT_fmax - aiT_fprime) / (aiT_fmax - aiT_faverage); 
  }
  else  {
    lort_pc =  aiT_k3;
  }

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc << ':' << geverbosepc_labelstep
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lb_eqFprimeFaverage = " << lb_eqFprimeFaverage;
    printf(" lort_pc = %0.10lf",lort_pc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lort_pc;
}

  
}/*END namespace*/

#endif /*PROBABILITY_DISTRIBUTION_HPP*/
