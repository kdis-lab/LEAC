/*! \file ga_generic_operator.hpp
 *
 * \brief genetic operators independent of the data type  
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef GA_GENERIC_OPERATOR_HPP
#define GA_GENERIC_OPERATOR_HPP

#include <random>
#include <iterator>
#include "chromosome_string.hpp"
#include "chromosome_fixedlength.hpp"
#include "chromosome_variablelength.hpp"
#include "ga_function_objective.hpp"
#include "verbose_global.hpp"


extern StdMT19937       gmt19937_eng;

/*! \namespace gagenericop
  \brief Genetic operators independent of the data type
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace gagenericop {

/*! \fn void initializeGenes(ITERATOR iterator_first, const ITERATOR iterator_last, const FUNCTION function)
  \brief Initialize genes
  \details
  \param iterator_first a conteiner iterator where items begin 
  \param iterator_last a conteiner iterator where items end
  \param function a uniform integer _distribution for generate genes
  \code{.cpp}

    for (auto lchromfixleng_iter :lvectorchromfixleng_population) {
      
      gaintegerop::initializeGenes
	(lchromfixleng_iter->begin(),
	 lchromfixleng_iter->end(),
	 [&]() 
	 {
	   return uniformdis_cidx_0K(gmt19937_eng);
	 }
	 );   
    }

  \endcode
 */
template<typename ITERATOR, typename FUNCTION>
inline
void
initializeGenes
(ITERATOR       iterator_first,
 const ITERATOR iterator_last,
 const FUNCTION function
 ) 
{

  for (; iterator_first != iterator_last; ++iterator_first) {
    *iterator_first = function();
  }

}
  
/*! \fn void onePointCrossover(gaencode::ChromFixedLength<T_GENE,T_METRIC> &aochrom_child1, gaencode::ChromFixedLength<T_GENE,T_METRIC> &aochrom_child2, const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent1, const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent2) 
  \brief single-point crossover (Michalewicz, 1992).
  \details
  \param aochrom_child1 a gaencode::ChromFixedLength<T_GENE,T_METRIC>
  \param aochrom_child2 a gaencode::ChromFixedLength<T_GENE,T_METRIC>
  \param aichrom_parent1 a gaencode::ChromFixedLength<T_GENE,T_METRIC>
  \param aichrom_parent2 a ChromFixedLength<T_GENE,T_METRIC> 

  \code{.cpp}

  gaiterator::crossover
    (lvectorchromfixleng_matingPool.begin(),
     lvectorchromfixleng_matingPool.end(),
     lvectorchromfixleng_population.begin(),
     lvectorchromfixleng_population.end(),
     [&](const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* aichrom_parent1,
         const gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>* aichrom_parent2,
         gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>*  aochrom_child1, 
         gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL>*  aochrom_child2
        )
     {

        if ( uniformdis_real01(gmt19937_eng) < 
             aiinpcgaprobfixedk_inParamGA.getProbCrossover()  ) {
	   
           gagenericop::onePointCrossover
             (*aochrom_child1,
              *aochrom_child2,
              *aichrom_parent1,
              *aichrom_parent2
             );

        } //if  Crossover
        else {
          *aochrom_child1 = *aichrom_parent1;
          *aochrom_child2 = *aichrom_parent2;
        }
     }
   );
    
  \endcode
 */
template < typename T_GENE,
	   typename T_METRIC
	   >
void
onePointCrossover
(gaencode::ChromFixedLength<T_GENE,T_METRIC>       &aochrom_child1,
 gaencode::ChromFixedLength<T_GENE,T_METRIC>       &aochrom_child2, 
 const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent1, 
 const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent2
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gatplop::onePointCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "(\noutput gaencode::ChromFixedLength<T_GENE>: aochrom_child1[" 
	      << &aochrom_child1 << "]\n"
	      << "output gaencode::ChromFixedLength<T_GENE>: aochrom_child2[" 
	      << &aochrom_child2 << "]\n";
    std::ostringstream lostrstream_labelparent1;
    lostrstream_labelparent1 << lpc_labelFunc << ":input parent1";
    aichrom_parent1.print(std::cout,lostrstream_labelparent1.str().c_str());
    std::cout << '\n';
    
    std::ostringstream lostrstream_labelparent2;
    lostrstream_labelparent2 << lpc_labelFunc << ":input parent2";
    aichrom_parent2.print(std::cout,lostrstream_labelparent2.str().c_str());
    std::cout << '\n';
  }
#endif /*__VERBOSE_YES*/

  
  static std::uniform_int_distribution<uintidx> uniformdis_uiCrossover1l_1
    (1,gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize()-1);

  const uintidx  lui_randPositionGene(uniformdis_uiCrossover1l_1(gmt19937_eng));
    
  interfacesse::copy
    (aochrom_child1.getString(), 
     aichrom_parent1.getString(), 
     lui_randPositionGene
     );
  interfacesse::copy
    (aochrom_child2.getString(), 
     aichrom_parent2.getString(), 
     lui_randPositionGene
     );
  interfacesse::copy
    (aochrom_child1.getString() + lui_randPositionGene, 
     aichrom_parent2.getString() + lui_randPositionGene, 
     aichrom_parent2.getStringSize() - lui_randPositionGene
     );
  interfacesse::copy
    (aochrom_child2.getString() + lui_randPositionGene, 
     aichrom_parent1.getString() + lui_randPositionGene, 
     aichrom_parent1.getStringSize() - lui_randPositionGene
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
 	      << "\noutput uintidx: lui_randPositionGene = " 
	      <<  lui_randPositionGene << '\n';
    aochrom_child1.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aochrom_child2.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*gabinaryop::onePointCrossove*/


/*! \fn void onePointCrossover(gaencode::ChromVariableLength<T_GENE,T_METRIC> &aiochrom_child1, gaencode::ChromVariableLength<T_GENE,T_METRIC> &aiochrom_child2) 
    \brief One point length variable crossover \cite He:Tan:GAclusteringVarK:TGCA:2012 
    \details
    \param aiochrom_child1 a ChromVariableLength<T_GENE,T_METRIC> 
    \param aiochrom_child2 a ChromVariableLength<T_GENE,T_METRIC>
 */
template < typename T_GENE,
	   typename T_METRIC
	   >
void
onePointCrossover
(gaencode::ChromVariableLength<T_GENE,T_METRIC> &aiochrom_child1,
 gaencode::ChromVariableLength<T_GENE,T_METRIC> &aiochrom_child2
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gatplop::onePointLengthVarCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output gaencode::ChromVariableLength<T_GENE>: aiochrom_child1[" 
	      << &aiochrom_child1 << "]\n"
	      << "\t output gaencode::ChromVariableLength<T_GENE>: aiochrom_child2[" 
	      << &aiochrom_child2 << "]\n";

    std::ostringstream lostrstream_labelChrom1;
    lostrstream_labelChrom1
      << lpc_labelFunc;
     
    aiochrom_child1.print(std::cout,lostrstream_labelChrom1.str().c_str());
    std::cout << '\n';
    
    std::ostringstream lostrstream_labelChrom2;
    lostrstream_labelChrom2
      << lpc_labelFunc;
     
    aiochrom_child2.print(std::cout,lostrstream_labelChrom2.str().c_str());
    std::cout << std::endl;
    
    std::cout << ")"
	      << std::endl; 
  }
#endif //__VERBOSE_YES

 uintidx lui_minStringLength =
    std::min(aiochrom_child1.getStringSize(),aiochrom_child2.getStringSize());
 std::uniform_int_distribution<uintidx> uniformdis_ui0N
   (0,lui_minStringLength-1);
 uintidx lui_positionGene = uniformdis_ui0N(gmt19937_eng);

 T_GENE* larrar_genes = new T_GENE[lui_minStringLength - lui_positionGene];
   
  interfacesse::copy
    (larrar_genes, 
     aiochrom_child1.getString()  + lui_positionGene,
     lui_minStringLength - lui_positionGene
     );
  interfacesse::copy
    (aiochrom_child1.getString()  + lui_positionGene,
     aiochrom_child2.getString() +  lui_positionGene,
     lui_minStringLength - lui_positionGene
     );

  interfacesse::copy
    (aiochrom_child2.getString() + lui_positionGene,
     larrar_genes, 
     lui_minStringLength - lui_positionGene
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelChrom1;
    lostrstream_labelChrom1
      << lpc_labelFunc
      << ":lui_positionGene = " << lui_positionGene;
    aiochrom_child1.print(std::cout,lostrstream_labelChrom1.str().c_str());
    std::cout << '\n';
    
    std::ostringstream lostrstream_labelChrom2;
    lostrstream_labelChrom2
      << lpc_labelFunc
      << ":lui_positionGene = " << lui_positionGene;
    aiochrom_child2.print(std::cout,lostrstream_labelChrom2.str().c_str());
    std::cout << std::endl;
	
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  delete [] larrar_genes;

}

/*PATH-BASED CROSSOVER:---------------------------------------------------------
 */

/*! \fn uintidx GAGRdist(const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_vector1, const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_vector2)
  \brief GAGRdist Definition 2. \cite Chang:etal:GAclustering:GAGR:2009
  \details
  \param aichrom_vector1 a gaencode::ChromFixedLength
  \param aichrom_vector2 a gaencode::ChromFixedLength
*/
template <typename T_GENE, typename T_METRIC >
uintidx   
GAGRdist
(const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_vector1, 
 const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_vector2
 )
{
  uintidx lol_sumDist;

#ifdef __VERBOSE_YES 
  const char* lpc_labelFunc = "gagenericop::GAGRdist";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "(\noutput gaencode::ChromFixedLength<>: aichrom_vector1[" 
	      << &aichrom_vector1 << "]\n"
	      << "output gaencode::ChromFixedLength<>: aichrom_vector2[" 
	      << &aichrom_vector2 << "]\n"
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  const T_GENE* lstr_vector1 = aichrom_vector1.getString();
  const T_GENE* lstr_vector2 = aichrom_vector2.getString();
  const T_GENE* lstr_end
    (lstr_vector1 + gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize());
  
  lol_sumDist = 0;
  while ( lstr_vector1 != lstr_end ) { 
    lol_sumDist += (*lstr_vector1 == *lstr_vector2)?0:1;
    lstr_vector1++;
    lstr_vector2++;
  }

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " lol_sumDist = " << lol_sumDist
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lol_sumDist;
}

/*! \fn void GAGRx1x2bar(T_GENE *aoaT_xr, uintidx  aist_xrLength, T_GENE *aiaT_x1, T_GENE *aiaT_x2)
  \brief GAGRx1x2bar Definition 3. \cite Chang:etal:GAclustering:GAGR:2009
*/
template <typename T_GENE >
void 
GAGRx1x2bar
(T_GENE *aoaT_xr, 
 uintidx  aist_xrLength, 
 T_GENE *aiaT_x1, 
 T_GENE *aiaT_x2
 )
{
  for (uintidx luintidx_i = 0; luintidx_i < aist_xrLength; luintidx_i++) {
    aoaT_xr[luintidx_i] =  (aiaT_x1[luintidx_i] > aiaT_x2[luintidx_i])?aiaT_x1[luintidx_i]: aiaT_x2[luintidx_i];
  }
}

/*! \fn void GAGRx1x2up (T_GENE *aoaT_xr, uintidx aist_xrLength, T_GENE *aiaT_x1, T_GENE  *aiaT_x2)
  \brief GAGRx1x2up Definition 4. \cite Chang:etal:GAclustering:GAGR:2009
*/
template <typename T_GENE >
void 
GAGRx1x2up
(T_GENE    *aoaT_xr, 
 uintidx   aist_xrLength, 
 T_GENE    *aiaT_x1, 
 T_GENE    *aiaT_x2
 )
{
  uintidx list_j;

  list_j = 0;
  array_setArray(aoaT_xr,aiaT_x1,aist_xrLength);
  while (list_j < aist_xrLength && aiaT_x1[list_j] >= aiaT_x2[list_j]) { 
    list_j++;
  }
  if ( list_j < aist_xrLength) {
    aoaT_xr[list_j] = aiaT_x2[list_j];
  }
}

/*! \fn void GAGRx1x2down (T_GENE *aoaT_xr, uintidx aist_xrLength, T_GENE *aiaT_x1, T_GENE *aiaT_x2)
  \brief GAGRx1x2down Definition 5. \cite Chang:etal:GAclustering:GAGR:2009
*/
template <typename T_GENE >
void 
GAGRx1x2down
(T_GENE    *aoaT_xr, 
 uintidx   aist_xrLength, 
 T_GENE    *aiaT_x1, 
 T_GENE    *aiaT_x2
 )
{
  uintidx ll_k;

  arrayt_setArray(aoaT_xr, aiaT_x1, aist_xrLength);
  --aoaT_xr;
  --aiaT_x1;
  --aiaT_x2;
  ll_k = aist_xrLength -1;
  while ((ll_k > 0) && (aiaT_x1[ll_k] < aiaT_x2[ll_k])) {
    ll_k--;
  }
  if ( 0 < ll_k && ll_k <= aist_xrLength) 
    aoaT_xr[ll_k] = aiaT_x2[ll_k];
}


/*! \fn uintidx GAGRcrossoverpathM1n (T_GENE *aioaf_m1n, uintidx aioaf_m1nLength, uintidx aii_initialIdx, T_GENE *aiaf_m1m2Bar)
    \brief  GAGRcrossoverpathM1n \cite Chang:etal:GAclustering:GAGR:2009
 */
template <typename T_GENE >
uintidx 
GAGRcrossoverpathM1n
(T_GENE *aioaf_m1n, 
 uintidx   aioaf_m1nLength,
 uintidx   aii_initialIdx, 
 T_GENE *aiaf_m1m2Bar
 )
{
  while (aii_initialIdx < aioaf_m1nLength && aioaf_m1n[aii_initialIdx]
	 >= aiaf_m1m2Bar[aii_initialIdx]) { 
    aii_initialIdx++;
  }
  if ( aii_initialIdx < aioaf_m1nLength) {
    aioaf_m1n[aii_initialIdx] = aiaf_m1m2Bar[aii_initialIdx];
  }
  return (aii_initialIdx + 1);
}

/*! \fn uintidx GAGRcrossoverpathM2n(T_GENE  *aioaf_m2n, uintidx aii_initialIdx, T_GENE  *aiaf_m2)
  \brief GAGRcrossoverpathM2n \cite Chang:etal:GAclustering:GAGR:2009
*/  
template <typename T_GENE >
uintidx 
GAGRcrossoverpathM2n
(T_GENE  *aioaf_m2n, 
 uintidx aii_initialIdx, 
 T_GENE  *aiaf_m2
 )
{
  /* Parameter adjustments */
  --aioaf_m2n;
  --aiaf_m2;
  while ((aii_initialIdx > 0) && (aioaf_m2n[aii_initialIdx] <= aiaf_m2[aii_initialIdx])) {
    aii_initialIdx--;
  }
  if ( 0 < aii_initialIdx ) 
    aioaf_m2n[aii_initialIdx] = aiaf_m2[aii_initialIdx];

  return aii_initialIdx;
}


/*! \fn void  pathCrossover(gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_child1, gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_child2, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parent1, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parent2, uintidx aiuintidx_GAGRdistM1M2, gafuncobj::GAFunctionObjective<T_GENE,T_METRIC> &aigafo_gaFunctionObjective)
  \brief pathCrossover  \cite Chang:etal:GAclustering:GAGR:2009
  \details Path crossover
  \param aochrom_child1 get chromosome best 
  \param aochrom_child2 get chromosome random
  \param aichrom_parent1 a parent chromosome  
  \param aichrom_parent2 a parent chromosome
  \param aiuintidx_GAGRdistM1M2 the distance between aichrom_parent1 and aichrom_parent1
  \param aigafo_gaFunctionObjective an object to obtain the objective function when the operator is applied
 */
template <typename T_GENE,
	  typename T_METRIC
	   >
void  
pathCrossover
(gaencode::ChromFixedLength<T_GENE,T_METRIC>       &aochrom_child1, 
 gaencode::ChromFixedLength<T_GENE,T_METRIC>       &aochrom_child2,
 const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent1,
 const gaencode::ChromFixedLength<T_GENE,T_METRIC> &aichrom_parent2,
 const uintidx                                     aiuintidx_GAGRdistM1M2,
 gafuncobj::GAFunctionObjective<T_GENE,T_METRIC>   &aigafo_gaFunctionObjective
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gagenericop::GAGRpathCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelChromParent1;
    lostrstream_labelChromParent1
      << "PARENT1:"
      << lpc_labelFunc;	       
    aichrom_parent1.print(std::cout,lostrstream_labelChromParent1.str().c_str());
    
    std::cout << '\n';
    std::ostringstream lostrstream_labelChromParent2;
    lostrstream_labelChromParent2 
      << "PARENT2:"
      << lpc_labelFunc;	       
    aichrom_parent2.print(std::cout,lostrstream_labelChromParent2.str().c_str());
    
    std::cout << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  

  std::uniform_int_distribution<uintidx> uniformdis_uintidx0GAGRdistM1M2
		(0,aiuintidx_GAGRdistM1M2);

  uintidx    lit_randGAGRdistM1M2(uniformdis_uintidx0GAGRdistM1M2(gmt19937_eng));
  
  uintidx    lst_j;
  T_METRIC   lfitT_fitnessM11;
  
  T_GENE *larrayT_m11 =  
    new T_GENE[gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize()];  
  T_GENE *larrayT_m2 = 
    new T_GENE[gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize()];
  T_GENE *larrayT_m1m2b = 
    new T_GENE[gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize()];

  
  interfacesse::copy
    (larrayT_m11,
     aichrom_parent1.getString(),
     gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize()
     );

  interfacesse::copy
    (larrayT_m2,
     aichrom_parent2.getString(),
     gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize()
     );
  
  uintidx lst_idxPath = 0;
   
  gagenericop::GAGRx1x2bar
    (larrayT_m1m2b, 
     aochrom_child1.getStringSize(),
     larrayT_m11, 
     larrayT_m2
     );
   
#ifdef __VERBOSE_YES
  uintidx      lst_pathi = 1;  
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "\naiT_GAGRdistM1M2 = ";
   
    inout::containerprint
      (larrayT_m11,
       larrayT_m11+gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize(),
       std::cout,
       "larrayT_m11 = ",
       ','
       );
    std::cout << '\n';
    
    inout::containerprint
      (larrayT_m2,
       larrayT_m2+gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize(),
       std::cout,
       "larrayT_m2 = ",
       ','
       );
    std::cout << '\n';

    inout::containerprint
      (larrayT_m1m2b,
       larrayT_m1m2b+gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize(),
       std::cout,
       "larrayT_m1m2b = "
       );
    std::cout << '\n';
    
    std::cout << "\nlit_randGAGRdistM1M2 = " << lit_randGAGRdistM1M2
	      << "\nlst_idxPath = "  << lst_idxPath << ' ';
    inout::containerprint
      (larrayT_m11,
       larrayT_m11+gaencode::ChromFixedLength<T_GENE,T_METRIC>::stcgetStringSize(),
       std::cout,
       "M11 = ",
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  if ( lit_randGAGRdistM1M2 == lst_idxPath ) { 
    aochrom_child2.setString(larrayT_m11);
    aochrom_child2.setFitness(aichrom_parent1.getFitness());
    aochrom_child2.setObjetiveFunc(aichrom_parent1.getObjetiveFunc());
    aochrom_child2.setValidString(aichrom_parent1.getValidString());
  }
  lst_j       = 0;
  lst_idxPath = 1;
  while ( (lst_j = gagenericop::GAGRcrossoverpathM1n
	   (larrayT_m11, aochrom_child1.getStringSize(), lst_j, larrayT_m1m2b)
	   ) 
	  <=  aochrom_child1.getStringSize() ) {

    std::pair<T_METRIC,bool> lpair_objFuncSSEM11 = 
      aigafo_gaFunctionObjective.getObjetiveFunc(larrayT_m11);
    lfitT_fitnessM11 = aigafo_gaFunctionObjective.getFitness(lpair_objFuncSSEM11.first);

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<"\nlst_idxPath = " << lst_idxPath 
		<< " M1" << lst_pathi+1 
		<< " = (M1" << lst_pathi;
      inout::containerprint
	(larrayT_m11,
	 larrayT_m11+aochrom_child1.getStringSize(),
	 std::cout,
	 " ^  M1M2b) =",
	 ','
	 );
      std::cout << "\n\tlT_m11Fitness = "  <<  lfitT_fitnessM11
		<< " lT_m11ObjetiveFunc = " << lpair_objFuncSSEM11.first; 
    }
    --geiinparam_verbose;
    lst_pathi= 2;
#endif //__VERBOSE_YES
  
    if ( aochrom_child1.getFitness() < lfitT_fitnessM11 ) { 
      aochrom_child1.setString(larrayT_m11);
      aochrom_child1.setFitness(lfitT_fitnessM11 ); 
      aochrom_child1.setObjetiveFunc( lpair_objFuncSSEM11.first ); 
      aochrom_child1.setValidString( lpair_objFuncSSEM11.second ); 
    }
    if ( lit_randGAGRdistM1M2 == lst_idxPath ) {
      aochrom_child2.setString( larrayT_m11 );
      aochrom_child2.setFitness( lfitT_fitnessM11 ); 
      aochrom_child2.setObjetiveFunc( lpair_objFuncSSEM11.first );
      aochrom_child2.setValidString( lpair_objFuncSSEM11.second );
    }
    lst_idxPath++;
  }
  lst_j = aochrom_child1.getStringSize() - 1;

  while ( (lst_j = gagenericop::GAGRcrossoverpathM2n(larrayT_m11,lst_j,larrayT_m2)) > 0) {
    
    std::pair<T_METRIC,bool>  lpair_objFuncSSEM11 = 
      aigafo_gaFunctionObjective.getObjetiveFunc(larrayT_m11);
    lfitT_fitnessM11 = aigafo_gaFunctionObjective.getFitness(lpair_objFuncSSEM11.first);
    
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout <<"\nidxPath = " << lst_idxPath 
		<< " M2" << lst_pathi 
		<< " = (M2" << lst_pathi-1;  
      inout::containerprint
	(larrayT_m11,
	 larrayT_m11+aochrom_child1.getStringSize(),
	 std::cout,
	 " v M2) = ",
	 ','
	 );
      std::cout << "\n\tlT_m11Fitness = "  <<  lfitT_fitnessM11
		<< " lT_m11ObjetiveFunc = " << lpair_objFuncSSEM11.first
		<< '\n';
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    if ( aochrom_child1.getFitness() < lfitT_fitnessM11 ) { 
      aochrom_child1.setString(larrayT_m11);
      aochrom_child1.setFitness(lfitT_fitnessM11 ); 
      aochrom_child1.setObjetiveFunc( lpair_objFuncSSEM11.first ); 
      aochrom_child1.setValidString( lpair_objFuncSSEM11.second ); 
    }
    if ( lit_randGAGRdistM1M2 == lst_idxPath ) {
      aochrom_child2.setString(larrayT_m11);
      aochrom_child2.setFitness(lfitT_fitnessM11 ); 
      aochrom_child2.setObjetiveFunc( lpair_objFuncSSEM11.first );
      aochrom_child2.setValidString( lpair_objFuncSSEM11.second );
    }
    
    lst_idxPath++;
  }
  
#ifdef __VERBOSE_YES
   if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelChromChild1;
    lostrstream_labelChromChild1 << "CHILD1:"
				 << lpc_labelFunc;
    aochrom_child1.print(std::cout,lostrstream_labelChromChild1.str().c_str());
    std::cout << '\n';
    std::ostringstream lostrstream_labelChromChild2;
    lostrstream_labelChromChild2 << "CHILD2:"
				 << lpc_labelFunc;	       
    aochrom_child2.print(std::cout,lostrstream_labelChromChild2.str().c_str());
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  delete[] larrayT_m1m2b;
  delete[] larrayT_m2;
  delete[] larrayT_m11;
}
  
}  /*END namespace gagenericop*/

#endif /*GA_GENERIC_OPERATOR_HPP*/
