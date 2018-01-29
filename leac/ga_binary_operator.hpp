/*! \file ga_binary_operator.hpp
 *
 * \brief traditional genetic binary operators
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef GA_BINARY_OPERATOR_HPP
#define GA_BINARY_OPERATOR_HPP

#include <vector>
#include <memory>
#include <limits>
#include <utility>
#include <iterator>     // std::distance
#include <typeinfo>
#include "chromosome_bitarray.hpp"
#include "crisp_matrix.hpp"
#include "probability_selection.hpp"

#include "verbose_global.hpp"

extern std::mt19937       gmt19937_eng;

/*! \namespace gabinaryop
  \brief Genetic binary operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gabinaryop {

/*! \fn void initializeGenes(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> &aobitarray_chrom, const FUNCTION function)
  \brief initialize genes where function evaluate to 1
  \details
  \param aobitarray_chrom OUT a gaencode::ChromosomeBitArray  
  \param function  IN  evalute in {0,1} where value is 1 set gene to 1 other 0

  \code{.cpp}

  for ( auto lchrombitarray_iter: lvectorchrom_population) {
      
      gabinaryop::initializeGenes
	(*lchrombitarray_iter,
	 [&]() 
	 {
	   return (uniformdis_real01(gmt19937_eng) <= aiinpkkunchevabezdek_inParamGA.getPini());
	 }
	 );
   }

   \endcode

 */
template < typename T_BITSIZE,typename T_REAL, typename FUNCTION >
void
initializeGenes
(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> &aobitarray_chrom,
 const FUNCTION                                 function
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::initializeGenes";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: this["
	      << &aobitarray_chrom << "]\n"
	      << ')'
	      << std::endl;		
  }
#endif /*__VERBOSE_YES*/
    

  aobitarray_chrom.initialize();
    
  for (uintidx li_j = 0; li_j < aobitarray_chrom.size(); li_j++) {
    if ( function() == 1 ) {
      aobitarray_chrom.setBit(li_j);
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aobitarray_chrom.print(std::cout,lpc_labelFunc);
    std::cout	<< std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
}

/*! \fn void onePointDistCrossover(BitMatrix<T_BITSIZE> &aobitmatrix_child1, mat::BitMatrix<T_BITSIZE> &aobitmatrix_child2, mat::BitMatrix<T_BITSIZE> &aibitmatrix_parent1, mat::BitMatrix<T_BITSIZE> &aibitmatrix_parent2) 
    \brief Crossover two matrices \cite Bezdek:etal:GAclustering:GA:1994 
    \details The crossover point and number of columns in the two two matrices. The columns of the matrices are combined to create the children matrices. 
    \param aobitmatrix_child1 OUT a mat::BitMatrix<T_BITSIZE> 
    \param aobitmatrix_child2 OUT a mat::BitMatrix<T_BITSIZE> 
    \parem aibitmatrix_parent1 IN a mat::BitMatrix<T_BITSIZE>
    \parem aibitmatrix_parent2 IN a mat::BitMatrix<T_BITSIZE>
    
    \code{.cpp}

    gaiterator::crossoverRandSelect
      (lvecchromcrispmatrix_matingPool.begin(),
       lvecchromcrispmatrix_matingPool.end(),
       lvecchromcrispmatrix_childR.begin(),
       lvecchromcrispmatrix_childR.end(),
       [&](gaencode::ChromosomeCrispMatrix<unsigned int,int,double>* aichrom_parent1,
	   gaencode::ChromosomeCrispMatrix<unsigned int,int,double>* aichrom_parent2,
	   gaencode::ChromosomeCrispMatrix<unsigned int,int,double>* aochrom_child1, 
	   gaencode::ChromosomeCrispMatrix<unsigned int,int,double>* aochrom_child2
	   )
       {
	 
	 gabinaryop::onePointDistBitCrossover
	   (*aochrom_child1,
	    *aochrom_child2,
	    *aichrom_parent1, 
	    *aichrom_parent2
	    );
	 
	 aochrom_child1->setFitness(std::numeric_limits<double>::max());
	 aochrom_child2->setFitness(std::numeric_limits<double>::max());
	   
       }
       );

    \endcode
*/
template <typename T_BITSIZE>
void
onePointDistCrossover
(mat::BitMatrix<T_BITSIZE>  &aobitmatrix_child1,
 mat::BitMatrix<T_BITSIZE>  &aobitmatrix_child2, 
 mat::BitMatrix<T_BITSIZE>  &aibitmatrix_parent1, 
 mat::BitMatrix<T_BITSIZE>  &aibitmatrix_parent2
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::onePointDistCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
  	      << "(output mat::BitMatrix<T_BITSIZE>: aobitmatrix_child1["
	      << &aobitmatrix_child1 << "]\n"
	      << " output mat::BitMatrix<T_BITSIZE>: aobitmatrix_child2["
	      << &aobitmatrix_child2 << "]\n";  
    aibitmatrix_parent1.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << '\n';
    aibitmatrix_parent2.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << std::endl;
    
  }
#endif /*__VERBOSE_YES*/

  static std::uniform_int_distribution<uintidx> lsuniformdis_uintidx0N
    (0, aobitmatrix_child1.getNumColumns()-1);

  const std::pair<uintidx,uintidx> lpair_pointCrossoverOrd = 
    prob::getRandPairUnlikeInOrd
    ([&]() -> uintidx
     {
       return lsuniformdis_uintidx0N(gmt19937_eng);
     }
     );
  
  const auto lui_distCrossover =
    lpair_pointCrossoverOrd.second - lpair_pointCrossoverOrd.first + 1;

  aobitmatrix_child1 = aibitmatrix_parent1;
  aobitmatrix_child2 = aibitmatrix_parent2;
  
  aobitmatrix_child1.copyAligned
    (aibitmatrix_parent2,
     lpair_pointCrossoverOrd.first,
     lui_distCrossover
     );
  aobitmatrix_child2.copyAligned
    (aibitmatrix_parent1,
     lpair_pointCrossoverOrd.first,
     lui_distCrossover
     );
 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {

    std::string lstr_termination = "OK";
    for ( uintidx  lui_i = 0; lui_i < aobitmatrix_child1.getNumRows();  lui_i++) { 
      for ( uintidx  lui_j = lpair_pointCrossoverOrd.first; lui_j <=lpair_pointCrossoverOrd.second;lui_j++) {
	if ( aobitmatrix_child1(lui_i,lui_j) != aibitmatrix_parent2(lui_i,lui_j) &&
	     aobitmatrix_child2(lui_i,lui_j) != aibitmatrix_parent1(lui_i,lui_j)) {
	  lstr_termination = "FAILURE";
	  break;
	}
	     
      }
      }
      
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')' << " Termination:" << lstr_termination
	      << ":lpair_pointCrossoverOrd["
	      << lpair_pointCrossoverOrd.first
	      << "," << lpair_pointCrossoverOrd.second
	      << "]\n";
    aobitmatrix_child1.print
      (std::cout,
       lpc_labelFunc,
       '\0',
       '\n'
       );
    std::cout << '\n';
    aobitmatrix_child2.print
      (std::cout,
       lpc_labelFunc,
       '\0',
       '\n'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*! \fn void onePointCrossover(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL> &aobitarray_chromChild1, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarray_chromChild2, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarray_chromParent1, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarray_chromParent2)
    \brief  One point crossover 
    \param aobitarray_chromChild1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aobitarray_chromChild1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aibitarray_chromParent1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aibitarray_chromParent2 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
*/
template < typename T_BITSIZE, typename T_REAL >
void
onePointCrossover
(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>        &aobitarray_chromChild1,
 gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>        &aobitarray_chromChild2,
 const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarray_chromParent1, 
 const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarray_chromParent2
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::onePointCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aobitarray_chromChild1[" 
	      << &aobitarray_chromChild1 << "]\n"
	      << " output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aobitarray_chromChild2[" 
	      << &aobitarray_chromChild2 << "]\n";
    
    aibitarray_chromParent1.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aibitarray_chromParent2.print(std::cout,lpc_labelFunc);
    
    std::cout << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

 
  static std::uniform_int_distribution<uintidx> lstuniformdis_uiCrossover1N
    (1,(uintidx)(aibitarray_chromParent1.size()-1));
   
  const uintidx lui_randPositionGene = lstuniformdis_uiCrossover1N(gmt19937_eng);
   
  uintidx uintidx_distCrossover = aibitarray_chromParent1.size() - lui_randPositionGene; 
  
  aobitarray_chromChild1.copyAligned
    (aibitarray_chromParent1,
     0,
     lui_randPositionGene
     );
  aobitarray_chromChild2.copyAligned
    (aibitarray_chromParent2,
     0,
     lui_randPositionGene
     );
  
  aobitarray_chromChild1.copyAligned
    (aibitarray_chromParent2,
     lui_randPositionGene,
     uintidx_distCrossover
     );
  
  aobitarray_chromChild2.copyAligned
    (aibitarray_chromParent1,
     lui_randPositionGene,
     uintidx_distCrossover
     );

  
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << "lui_randPositionGene = " << lui_randPositionGene
	      << '\n';
    aobitarray_chromChild1.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aobitarray_chromChild2.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

/*! \fn void onePointDistCrossover(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child1, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child2, const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent1, const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent2) 
    \brief One point distance crossover
    \param aobitarraychrom_child1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aobitarraychrom_child2 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aibitarraychrom_parent1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aibitarraychrom_parent2 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>

    \note Use by \cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001
*/
template < typename T_BITSIZE,  typename T_REAL >
void 
onePointDistCrossover
(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child1,
 gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child2,
 const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent1, 
 const gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent2
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::onePointDistCrossover";
    ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
    	      << "(output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aobitarraychrom_child1[" 
	      << &aobitarraychrom_child1 << "]\n"
	      << " output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aobitarraychrom_child2[" 
	      << &aobitarraychrom_child2 << "]\n";
	     
    aibitarraychrom_parent1.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << '\n';
   
   aibitarraychrom_parent2.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
       
    std::cout << std::endl;
  }
#endif /*__VERBOSE_YES*/

  static std::uniform_int_distribution<uintidx> lsuniformdis_uintidx0N
    (0, aibitarraychrom_parent1.size() -1 );

  const std::pair<uintidx,uintidx> lpair_pointCrossoverOrd = 
    prob::getRandPairUnlikeInOrd
    ([&]() -> uintidx
     {
       return lsuniformdis_uintidx0N(gmt19937_eng);
     }
     );
  
  const auto lui_distCrossover = lpair_pointCrossoverOrd.second - lpair_pointCrossoverOrd.first + 1;
  
  aobitarraychrom_child1 = aibitarraychrom_parent1;
  aobitarraychrom_child2 = aibitarraychrom_parent2;
  
  aobitarraychrom_child1.copyAligned
    (aibitarraychrom_parent2,
     lpair_pointCrossoverOrd.first,
     lui_distCrossover
     );
  aobitarraychrom_child2.copyAligned
    (aibitarraychrom_parent1,
     lpair_pointCrossoverOrd.first,
     lui_distCrossover
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "gabinaryop::onePointDistCrossover:  OUT"
	      << '(' << geiinparam_verbose << ')'
	      << " lpair_pointCrossoverOrd: ["
	      << lpair_pointCrossoverOrd.first
	      << "," << lpair_pointCrossoverOrd.second
	      << "]\n";
    aobitarraychrom_child1.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << '\n';
    aobitarraychrom_child2.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
   
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*! \fn void uniformCrossover(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child1, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child2, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent1, gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent2, const T_REAL airt_probCrossover) 
    \brief Uniform crossover
     \details The parent chromosomes swap their i-th bits with a certain probability, and i goes from 1 to n length chromosome. 
    \param aobitarraychrom_child1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aobitarraychrom_child2 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aibitarraychrom_parent1 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param aibitarraychrom_parent2 a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param airt_probCrossover a real number probability swap bits

    \code{.cpp}

    gaiterator::crossoverRandSelect
      (lvectorchrom_matingPool.begin(),
       lvectorchrom_matingPool.end(),
       lvectorchrom_setO.begin(),
       lvectorchrom_setO.end(),
       [&](gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aichrom_parent1,
	   gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aichrom_parent2,
	   gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aochrom_child1, 
	   gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>* aochrom_child2
	   )
       {

	 gabinaryop::uniformCrossover
	   (*lchrom_child1,
	    *lchrom_child2,
	    *aichrom_parent1,
	    *aichrom_parent2,
	    aiinpkkunchevabezdek_inParamGA.getProbCrossover()
	    );
	  
	 aochrom_child1->setObjetiveFunc(std::numeric_limits<double>::max());
	 aochrom_child2->setObjetiveFunc(std::numeric_limits<double>::max());
	   
       }
       );
    \endcode

    \note Use by \cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997 
*/
template < typename T_BITSIZE,  typename T_REAL>
void 
uniformCrossover
(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child1,
 gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aobitarraychrom_child2,
 gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent1, 
 gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aibitarraychrom_parent2,
 const T_REAL                                    airt_probCrossover
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::uniformCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aobitarraychrom_child1[" 
	      << &aobitarraychrom_child1 << "]\n"
	      << " output gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aobitarraychrom_child2[" 
	      << &aobitarraychrom_child2 << "]\n"
	      << " input  gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aibitarraychrom_parent1[" 
	      << &aibitarraychrom_parent1 << "]\n"
	      << " input  gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: aibitarraychrom_parent2[" 
	      << &aibitarraychrom_parent2 << "]\n";
    aibitarraychrom_parent1.print(std::cout,"PARENT1");
    std::cout << '\n';
    aibitarraychrom_parent2.print(std::cout,"PARENT2");
    std::cout << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  static std::uniform_real_distribution<T_REAL> lsuniformdis_real01(0, 1);

  aobitarraychrom_child1 = aibitarraychrom_parent1;
  aobitarraychrom_child2 = aibitarraychrom_parent2;

  for (uintidx li_j = 0; li_j < aobitarraychrom_child1.size(); li_j++) {
    if (aobitarraychrom_child1.getBit(li_j) !=  aobitarraychrom_child2.getBit(li_j) ) {
      if ( lsuniformdis_real01(gmt19937_eng) <= airt_probCrossover ) {
	aobitarraychrom_child1.toggleBit(li_j);
	aobitarraychrom_child2.toggleBit(li_j);
      }
    }
  }


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aobitarraychrom_child1.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aobitarraychrom_child2.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
}


/*! \fn void bitMutation(CrispMatrix<T_BITSIZE,T_CLUSTERIDX>  &aiobitcrispmatrix_chrom) 
    \brief  Crisp matrix mutation
    \details Consists of randomly choosing an element of a column to have the value 1, such thar it is a different element than the one 
    currently having a value of 1. 
    \param aiobitcrispmatrix_chrom a mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>
  
    \code{.cpp}

    for ( auto ichrom_childR: lvecchromcrispmatrix_childR ) {
      gabinaryop::bitMutation(*ichrom_childR)
    }

    \endcode

    \note Use by \cite Bezdek:etal:GAclustering:GA:1994
*/
template <typename T_BITSIZE,
	  typename T_CLUSTERIDX
	  >
void
bitMutation
(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>  &aiobitcrispmatrix_chrom)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::bitMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n";
    aiobitcrispmatrix_chrom.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  static std::uniform_int_distribution<uintidx> lsuniformdis_uintidx0N
    (0, aiobitcrispmatrix_chrom.getNumColumns()-1);

  static std::uniform_int_distribution<T_CLUSTERIDX> lsuniformdis_cidx_0K
    (0, T_CLUSTERIDX(aiobitcrispmatrix_chrom.getNumRows()-1));
  
  const uintidx lui_mutationRandPosGene(lsuniformdis_uintidx0N(gmt19937_eng));
  
  T_CLUSTERIDX lcidx_mutationMemberPrev = 
    aiobitcrispmatrix_chrom.getMember(lui_mutationRandPosGene);

  T_CLUSTERIDX lcidx_mutationNewMember =
    prob::getRandUnlike
    (lcidx_mutationMemberPrev,
     [&]() -> T_CLUSTERIDX {
      return lsuniformdis_cidx_0K(gmt19937_eng);
      }
     );
  
    //function(lcidx_mutationMemberPrev);
  
  aiobitcrispmatrix_chrom.setMember(lui_mutationRandPosGene,lcidx_mutationNewMember);
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
   	      << "uintidx: lui_mutationRandPosGene = " << lui_mutationRandPosGene
	      << " lcidx_mutationMemberPrev = " << lcidx_mutationMemberPrev
	      << " lcidx_mutationNewMember = " << lcidx_mutationNewMember
	      << '\n';
    
    aiobitcrispmatrix_chrom.print
      (std::cout,
       lpc_labelFunc,
       ',',
       ';'
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
} /*bitCrispMatrixMutation*/



/*! \fn void bitMutation(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aiobitarray_chrom)
    \brief  Change one bit
    \details Consists of randomly choosing an element of a  change the value 1 by 0 or 0 by 1.   
    \param aiobitarray_chrom a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    
    \code{.cpp}
      for (auto&& liter_iChrom :lvectorchrom_population) {

	if ( uniformdis_real01(gmt19937_eng) 
	     <  aiinpcga_inParamGACasillas.getProbMutation() ) 
	  {
	    gabinaryop::bitArrayMutation(*liter_iChrom);
	  }
      }
    \endcode

    \note Use by \cite Casillas:etal:GAclusteringVarK:GA:2003
*/
template <typename T_BITSIZE, typename T_REAL >
inline
void
bitMutation
(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aiobitarray_chrom)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::bitArrayMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(outnput gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: ["
	      << &aiobitarray_chrom << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  
  static std::uniform_int_distribution<uintidx> lstuniformdis_uiMutation0N
     (0,(aiobitarray_chrom.size()-1));
  
  const uintidx aiuintidx_pointMutation(lstuniformdis_uiMutation0N(gmt19937_eng));
    
  aiobitarray_chrom.toggleBit(aiuintidx_pointMutation);
  

#ifdef __VERBOSE_YES
  
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")"
	      << " aiuintidx_pointMutation = " << aiuintidx_pointMutation << '\n';
    aiobitarray_chrom.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
}


/*! \fn void eachBitArrayMutation(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aiobitarray_chrom, const T_REAL airt_probMutation, const FUNCTION function)
    \brief Toggle bits
    \details Each bit of each offspring chromosome mutates with a predefined probability (mutation rate airt_probMutation).
    \param aiobitarray_chrom a gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>
    \param airt_probMutation a real number with mutation rate
    \param function a function generate number with probability distribution
    
    \code{.cpp}

      for ( auto lchrombitarray_iter: lvectorchrom_setO) {
      
         gabinaryop::eachBitArrayMutation
	  (*lchrombitarray_iter,
	   aiinpkkunchevabezdek_inParamGA.getProbMutation(),
	   [&]() 
	   {
	     return uniformdis_real01(gmt19937_eng);
	   }
	  );
      }
      
    \endcode

    \note Use by \cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997 
*/
template <typename T_BITSIZE, typename T_REAL, typename FUNCTION >
void
eachBitArrayMutation
(gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>  &aiobitarray_chrom,
 const T_REAL                          airt_probMutation,
 const FUNCTION                        function
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gabinaryop::eachBitArrayMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(outnput gaencode::ChromosomeBitArray<T_BITSIZE,T_REAL>: [" << &aiobitarray_chrom << ']'
	      << "\nuintidx: aiuintidx_pointMutation = " << airt_probMutation
	      << '\n';
    aiobitarray_chrom.print(std::cout,lpc_labelFunc);
    std::cout << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  for (uintidx li_j = 0; li_j < aiobitarray_chrom.size(); li_j++) {
      if ( function() <= airt_probMutation ) {
	aiobitarray_chrom.toggleBit(li_j);
      }
  }
  
#ifdef __VERBOSE_YES
  
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aiobitarray_chrom.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
}

}  /*END namespace gabinaryop*/

#endif /*GA_BINARY_OPERATOR_HPP*/
