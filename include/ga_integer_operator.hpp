/*! \file ga_integer_operator.hpp
 *
 * \brief genetic integer operators
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef GA_INTEGER_OPERATOR_HPP
#define GA_INTEGER_OPERATOR_HPP

#include <algorithm>    // std::find
#include <iterator>
#include <random>
#include "chromosome_fixedlength.hpp"
#include "linear_algebra_level1.hpp"
#include "vector_utils.hpp"
#include "probability_selection.hpp"
#include "clustering_operator_centroids.hpp"
#include "partition_label.hpp"
#include "dist.hpp"
#include "probability_distribution.hpp"
#include "probability_selection.hpp"

#include "verbose_global.hpp"

extern StdMT19937 gmt19937_eng;

/*! \namespace gaintegerop
  \brief Genetic integer operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaintegerop {



/*! \fn void labelKeep(gaencode::ChromosomeString<T_INTEGER,T_METRIC> &aochrom_relabel, std::vector<uintidx> &aivectorcidx_clustersKeep)
  \brief Label keep
  \details 
  \param aochrom_relabel a gaencode::ChromosomeString<T_INTEGER,T_METRIC>
  \param aivectorcidx_clustersKeep  a std::vector<uintidx> with label to keep
 */
template <typename T_INTEGER,
	  typename T_METRIC //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  >
void 
labelKeep
(gaencode::ChromosomeString<T_INTEGER,T_METRIC> &aochrom_relabel,
 std::vector<uintidx>                 &aivectorcidx_clustersKeep
 ) 
{
  T_INTEGER *lstr_chrom = aochrom_relabel.getString();
  T_INTEGER lcidx_numClustersKeep = 
    (T_INTEGER) aivectorcidx_clustersKeep.size();

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaintegerop::labelKeep";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output Chromosome&: aochrom_relabel[" << &aochrom_relabel << ']'
	      << "\ninput  vector<>&: aivectorcidx_clustersKeep["  
	      << &aivectorcidx_clustersKeep << "]: "
	      << aivectorcidx_clustersKeep
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::vector<T_INTEGER> lvectorinteger_label;

  lvectorinteger_label.reserve(aivectorcidx_clustersKeep.size());
  std::for_each
    (aivectorcidx_clustersKeep.begin(),
     aivectorcidx_clustersKeep.end(),
     [&](uintidx literuintidx_labelK) 
     {
       lvectorinteger_label.push_back((T_INTEGER) literuintidx_labelK);
     }
     ); 

  for (T_INTEGER lcidx_k= 0; lcidx_k < lcidx_numClustersKeep; lcidx_k++)  {
    if ( lvectorinteger_label[lcidx_k] != lcidx_k) {
      for (uintidx luintidx_i = 0; luintidx_i < aochrom_relabel.getStringSize(); luintidx_i++) { 
	if ( lstr_chrom[luintidx_i] == lvectorinteger_label[lcidx_k] )
	  lstr_chrom[luintidx_i] = lcidx_k;
      }
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelFunc;
    lostrstream_labelFunc << lpc_labelFunc; // "Chromosome";  
    aochrom_relabel.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

} /*END labelKeep*/



/*! \fn void oneChangelabel(gaencode::ChromosomeString<T_INTEGER,T_METRIC>  &aochrom_relabel, T_INTEGER  aigene_Cs, T_INTEGER  aigene_Cj)
  \brief One change label
  \details
  \param aochrom_relabel a gaencode::ChromosomeString<T_INTEGER,T_METRIC>
  \param aigene_Cs a integer number
  \param aigene_Cj a integer number 
 */
template <typename T_INTEGER,
	  typename T_METRIC //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  >
void 
oneChangelabel
(gaencode::ChromosomeString<T_INTEGER,T_METRIC>  &aochrom_relabel,
 T_INTEGER                             aigene_Cs,
 T_INTEGER                             aigene_Cj
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaintegerop::oneChangelabel";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output Chromosome&: aochrom_relabel[" << &aochrom_relabel << ']'
	      << "\n input  aigene_Cs = " << aigene_Cs 
	      << " aigene_Cj = " << aigene_Cj << " (Cs --> Cj)" 
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  T_INTEGER *lstr_chrom = aochrom_relabel.getString();
  for (uintidx luintidx_i = 0; luintidx_i < aochrom_relabel.getStringSize(); luintidx_i++) { 
    if ( lstr_chrom[luintidx_i] == aigene_Cs )
      lstr_chrom[luintidx_i] = aigene_Cj;
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:"
				<< geverbosepc_labelstep << ':' << lpc_labelFunc << ':' << geverboseui_idproc
				<< ":Chromosome"
				<< aochrom_relabel;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*--------------------------------------------------------------------------------
  CROSSOVER
  -------------------------------------------------------------------------------*/

/*recombinationD_MX
  

  In general, many recombination operators, in-
  cluding D_MX, when applied to two parent
  strings (P_1, P_2), produce two child strings
  (C_1, C_2):
            D_MX
  P_1, P_2  ---->  C_1, C_2
*/
  
/*! \fn void recombinationD_MX(gaencode::ChromosomeString<uintidx,T_REAL> &aochrom_child1, gaencode::ChromosomeString<uintidx,T_REAL> &aochrom_child2, const gaencode::ChromosomeString<uintidx,T_REAL> &aichrom_parent1, const gaencode::ChromosomeString<uintidx,T_REAL> &aichrom_parent2, const uintidx aiidxinst_numInstances, const T_REAL airt_mixMutationProb)
  \brief recombination D_MX
  \details \cite Lucasius:etal:GAclusteringMedoid:GCA:1993
  \param aochrom_child1 a gaencode::ChromosomeString<uintidx,T_REAL>
  \param aochrom_child2 a gaencode::ChromosomeString<uintidx,T_REAL>
  \param aichrom_parent1 a gaencode::ChromosomeString<uintidx,T_REAL>
  \param aichrom_parent2 a gaencode::ChromosomeString<uintidx,T_REAL>
  \param aiidxinst_numInstances a numeber with numero of instances 
  \param airt_mixMutationProb a real number mutation probability mix
 */
template <typename T_REAL>
void
recombinationD_MX
(gaencode::ChromosomeString<uintidx,T_REAL>       &aochrom_child1,
 gaencode::ChromosomeString<uintidx,T_REAL>       &aochrom_child2, 
 const gaencode::ChromosomeString<uintidx,T_REAL> &aichrom_parent1, 
 const gaencode::ChromosomeString<uintidx,T_REAL> &aichrom_parent2,
 const uintidx                                    aiidxinst_numInstances, 
 const T_REAL                                     airt_mixMutationProb
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaintegerop::recombinationD_MX"; 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax  ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output gaencode::ChromosomeString: aochrom_child1[" << &aochrom_child1 << ']'
	      << "\n output gaencode::ChromosomeString: aochrom_child2[" << &aochrom_child2 << ']'
	      << "\n input  gaencode::ChromosomeString: aichrom_parent1[" << &aichrom_parent1 << ']'
	      << "\n input  gaencode::ChromosomeString: aichrom_parent2[" << &aichrom_parent2 << ']'
	      << "\n input  uintidx: aiidxinst_numInstances = "  
	      << aiidxinst_numInstances
	      << "\n\t input  T_REAL: airt_mixMutationProb = " 
	      << airt_mixMutationProb
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  uintidx luintidx_sizeQ = aichrom_parent1.getStringSize() + aichrom_parent2.getStringSize();
  std::vector<uintidx> lvectoridxinst_stringQ(luintidx_sizeQ);
  
  /*eg.:  using P, = 2 3 7 and P2 = 4 8 2 (hence k = 3) 
    for illustration purposes:

    (1) Mix P, and P2:
    (a) Append copies of P, and P2 to obtain, say, Q:

    Q = 2_1, 3_1, 7_1, 4_2, 8_2, 2_2
  */
  interfacesse::copy
    (lvectoridxinst_stringQ.data(),
     aichrom_parent1.getString(),
     aichrom_parent1.getStringSize()
     );
  
  interfacesse::copy
    (&lvectoridxinst_stringQ[aichrom_parent1.getStringSize()],
     aichrom_parent2.getString(),
     aichrom_parent2.getStringSize()
     );
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax  ) {
    std::ostringstream lostrstream_labelQ;
    lostrstream_labelQ <<  "<(a) APPEND COPIES OF P1 AND P2 TO OBTAIN SAY, Q:" << lpc_labelFunc << ':';
    inout::containerprint
      (lvectoridxinst_stringQ.begin(),
       lvectoridxinst_stringQ.end(),
       std::cout,
       lostrstream_labelQ.str().c_str(),
       ','
       );
    std::cout << "\n";
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  /*(b) Randomly scramble the elements, e.g.:

    Q = 4_2, 2_2, 2_1, 8_2, 7_1, 3_1
   */
  std::shuffle
    (lvectoridxinst_stringQ.begin(),
     lvectoridxinst_stringQ.end(),
     gmt19937_eng
     );
  
  /*(2) Add new material; that is, apply the follow-
    ing built-in mix mutation: With a predetermined
    probability P_m.mix, replace each of the first k
    consecutive elements in Q by a copy of an ele-
    ment indicated randomly, but never more than
    once, in the source set, e.g.:

    Q = 5, 2_2, 7, 8_2, 7_1, 3_1

    Note that in this example, 2 of the k trials
    succeeded (at the first and third position).
   */
  std::uniform_int_distribution<uintidx> uniformdis_ui0M
    (0,aiidxinst_numInstances-1);
  for ( uintidx luintidx_i = 0; luintidx_i < aichrom_parent1.getStringSize(); luintidx_i++) {
    if ( uniformdis_real01(gmt19937_eng) < airt_mixMutationProb ) {
      lvectoridxinst_stringQ[luintidx_i] = uniformdis_ui0M(gmt19937_eng);
      
#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "(2) ADD NEW MATERIAL, Q[" << luintidx_i << "] = " 
		  << lvectoridxinst_stringQ[luintidx_i] << '\n'; 
      }
     --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

    }
  } /*END for luintidx_i*/

  /*(3) Randomly scramble the elements again,
    e.g.:

    Q = 2_2, 7_1, 7, 3_1, 5, 8_2
   */
  std::shuffle
    (lvectoridxinst_stringQ.begin(),
     lvectoridxinst_stringQ.end(),
     gmt19937_eng
     );
  
  /*(4) Build C, by copying k elements from Q,
    starting at the leftmost element and going ele-
    men&vise to the right, subject to the condition
    that elements that are already in C, are skipped:
    
    C_1 = 2, 7, 3

  */
  uintidx luintidx_k = 1;
  uintidx luintidx_l = 1;
  uintidx *larrayidxinstT_c1 = 
    aochrom_child1.getString();

  larrayidxinstT_c1[0] = lvectoridxinst_stringQ[0];
  while ( (luintidx_k < aochrom_child1.getStringSize() ) && (luintidx_l < luintidx_sizeQ)  ) {
    while ( (luintidx_l < luintidx_sizeQ) &&
	    (std::find
	     (larrayidxinstT_c1,
	      &larrayidxinstT_c1[luintidx_k], 
	      lvectoridxinst_stringQ[luintidx_l]) != &larrayidxinstT_c1[luintidx_k] ) 
	    ) 
      ++luintidx_l;
    if ( (luintidx_l < luintidx_sizeQ) )
      larrayidxinstT_c1[luintidx_k++] = lvectoridxinst_stringQ[luintidx_l++];
  }
  if ( luintidx_k < aochrom_child1.getStringSize() ) {
    interfacesse::copy
      (larrayidxinstT_c1,
       aichrom_parent1.getString(),
       aochrom_child1.getStringSize()
       );
  }

  /*(5)  Build C_2 by copying k elements from Q,
    starting at the rightmost element and going ele-
    mentwise to the left, subject to the condition that
    elements that are already in C_2 are skipped:

    C_2 = 8, 5, 3
  */
  luintidx_k = 1;
  luintidx_l =  luintidx_sizeQ - 1;
  uintidx *larrayidxinstT_c2 = 
    aochrom_child2.getString();

  larrayidxinstT_c2[0] = lvectoridxinst_stringQ[luintidx_l--];
  while ( (luintidx_k < aochrom_child2.getStringSize()) && (luintidx_l > 0)  ) {
    while ( (luintidx_l > 0) 
	    && 
	    (std::find
	     (larrayidxinstT_c2,
	      &larrayidxinstT_c2[luintidx_k], 
	      lvectoridxinst_stringQ[luintidx_l]) != &larrayidxinstT_c2[luintidx_k] )
	    )
      --luintidx_l;
    if (luintidx_l > 0)
      larrayidxinstT_c2[luintidx_k++] = lvectoridxinst_stringQ[luintidx_l--];
  }
  if ( luintidx_l == 0 && luintidx_k < aochrom_child2.getStringSize() ) 
    if ( (std::find
	     (larrayidxinstT_c2,
	      &larrayidxinstT_c2[luintidx_k], 
	      lvectoridxinst_stringQ[0]) != &larrayidxinstT_c2[luintidx_k] ) 
	 )
      larrayidxinstT_c2[luintidx_k++] = lvectoridxinst_stringQ[luintidx_l];
  if ( luintidx_k < aochrom_child2.getStringSize() ) {
    interfacesse::copy
      (larrayidxinstT_c2,
       aichrom_parent2.getString(),
       aochrom_child2.getStringSize()
       );
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
#endif /*__VERBOSE_YES*/

} /*recombinationD_MX*/

/*--------------------------------------------------------------------------------
  MUTATION
  -------------------------------------------------------------------------------*/

/*geneticopie_mutationLucasius:
  \cite{Lucasius:etal:GAclusteringMedoid:GCA:1993}

  Mutation. Besides the built-in mutation, a point
  mutation operator, D_PM, is independently ap-
  plied to all strings in the population: one element
  in a child string C is selected randomly with a
  predetermined probability, Pm.point, and, upon
  success, replaced by a copy of an element indi-
  cated randomly in the complementary subset to
  produce C':

     D_PM
  C ----> C'

  where, for instance, C = 3, 6, 8 and C' = 3, 4, 8. (In-
  cidentally, D_PM is the equivalent of one itera-
  tion in D_TM, trade mutation, described in
  [57,581.)

*/


  
/*! \fn void mutation(gaencode::ChromosomeString<T_INTEGER,T_METRIC> &aochrom_offspring, const uintidx aiui_randPositionGene, FUNCTION function)
  \brief Mutate a integer gene 
  \details
  \param aochrom_offspring a gaencode::ChromosomeString<T_INTEGER,T_METRIC>
  \param aiui_randPositionGene a unsigned integer
  \param function a function evaluate a gene
 */
template < typename T_INTEGER,
	   typename T_METRIC,
	   typename FUNCTION
	   >
std::pair<uintidx,T_INTEGER>
mutation 
(gaencode::ChromFixedLength<T_INTEGER,T_METRIC>  &aochrom_offspring,
 //const uintidx                                   aiui_randPositionGene,
 FUNCTION                                        function
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaintegerop::mutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output Chromosome<uintidx>: aochrom_offspring[" << &aochrom_offspring << "]\n"
      //<< " input  uintidx&: aiui_randPositionGene = " << aiui_randPositionGene 
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  static std::uniform_int_distribution<uintidx> uniformdis_uiMutation0N
    (0,gaencode::ChromFixedLength<T_INTEGER,T_METRIC>::stcgetStringSize()-1);

  const uintidx lui_randPositionGene(uniformdis_uiMutation0N(gmt19937_eng));
  std::pair<uintidx,T_INTEGER> lopair_previousGene
    (lui_randPositionGene,aochrom_offspring.getGene(lui_randPositionGene));
  
  
  //T_INTEGER lt_previousGene = aochrom_offspring.getGene(lui_randPositionGene);  
  T_INTEGER lt_newGene     = function(lopair_previousGene.second);

  aochrom_offspring.setGene(lui_randPositionGene,lt_newGene);
  
#ifdef __VERBOSE_YES

  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelOffspring;
    lostrstream_labelOffspring << "lui_randPositionGene = " << lopair_previousGene.first
			       << ", lcidx_previousGene = " << lopair_previousGene.second << " by " << lt_newGene;
    aochrom_offspring.print(std::cout,lostrstream_labelOffspring.str().c_str());	 
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lopair_previousGene;
  //return louintidx_positionGene;
}



/*! \fn bool  isNotValid(gaencode::ChromosomeString<T_INTEGER,T_METRIC>  &aichrom_chromosome, const T_INTEGER aiit_lowerRangeGenes, const T_INTEGER aiit_upperRangeGenes)
  \brief Check if not valid string
  \details Return true if the string has integers in the range aiit_lowerRangeGenes  to aiit_upperRangeGenes otherwise false
  \param aichrom_chromosome a gaencode::ChromosomeString<T_INTEGER,T_METRIC>
  \param aiit_lowerRangeGenes a integer lower value gene
  \param aiit_upperRangeGenes a integer upper value gene
 */
template < typename T_INTEGER,
	   typename T_METRIC
	   > 
bool 
isNotValid
(gaencode::ChromosomeString<T_INTEGER,T_METRIC>  &aichrom_chromosome,
 const T_INTEGER                                 aiit_lowerRangeGenes, 
 const T_INTEGER                                 aiit_upperRangeGenes
 )
{
  for (T_INTEGER li_l = aiit_lowerRangeGenes; li_l < aiit_upperRangeGenes; li_l++) {
    T_INTEGER *lptinteger_find =
      std::find
      (aichrom_chromosome.getString(), 
       aichrom_chromosome.getString() + aichrom_chromosome.getStringSize(),
       li_l
       );
      if ( lptinteger_find == (aichrom_chromosome.getString() + aichrom_chromosome.getStringSize()) )
	return true;
  } 
  return false;
}

/*! \fn void mutationgka (gaencode::ChromFixedLength<T_CLUSTERIDX,T_REAL> &aiochrom_mutate, const T_REAL airt_probMutation, mat::MatrixRow<T_FEATURE> &aomatrixt_centroids, mat::MatrixRow<T_FEATURE_SUM> &aomatrixt_sumInstancesCluster, std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstancesInClusterK, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist)
    \brief mutation GKA 
    \details mutation GKA based on \cite Krishna:Murty:GAClustering:GKA:1999
    \param aiochrom_mutate a gaencode::ChromFixedLength to mutate
    \param airt_probMutation a real number 
    \param aomatrixt_centroids a mat::MatrixRow space to work and store the centroids associated with the chromosome
    \param aomatrixt_sumInstancesCluster a mat::MatrixRow space to work and store the sum of instances per cluster
    \param aovectort_numInstancesInClusterK a std::vector space to work to store the number of instances per cluster
    \param aiiterator_instfirst a input iterator of the instances
    \param aiiterator_instlast a  const input iterator of the instances
    \param aifunc2p_dist an object of type dist::Dist to calculate distances
 */
template < typename T_FEATURE, 
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,  //-1, 0, 1, .., N
	   typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_REAL,
	   typename INPUT_ITERATOR
	   >
void
mutationgka
(gaencode::ChromFixedLength
 <T_CLUSTERIDX,T_REAL>              &aiochrom_mutate,
 const T_REAL                       airt_probMutation,      
 mat::MatrixRow<T_FEATURE>          &aomatrixt_centroids,
 mat::MatrixRow<T_FEATURE_SUM>      &aomatrixt_sumInstancesCluster,
 std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstancesInClusterK,
 INPUT_ITERATOR                     aiiterator_instfirst,
 const INPUT_ITERATOR               aiiterator_instlast,
 const dist::Dist<T_REAL,T_FEATURE> &aifunc2p_dist
 ) 
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaintegerop::mutationgka";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")"
	      << "\n(output ChromFixedLength&: aiochrom_mutate[" << &aiochrom_mutate << "]\n"
	      << "airt_probMutation = " << airt_probMutation << '\n'			
	      << "mat::MatrixRow<T_FEATURE>  aomatrixt_centroids"
	      <<  aomatrixt_centroids.getNumRows() << "x"
	      <<  aomatrixt_centroids.getNumColumns()
	      << "[" <<  &aomatrixt_centroids << "]\n"
	      << "MatrixRow<T_FEATURE_SUM> aomatrixt_sumInstancesCluster"
	      <<  aomatrixt_sumInstancesCluster.getNumRows() << "x"
	      <<  aomatrixt_sumInstancesCluster.getNumColumns()
	      << "[" <<  &aomatrixt_sumInstancesCluster << "]\n"
	      << "aovectort_numInstancesInClusterK"
	      << "[" << &aovectort_numInstancesInClusterK << "]\n"
	      << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << &aiiterator_instlast << "]\n"
	      << "input  dist::Dist<T_METRIC,T_FEATURE> &aifunc2p_dist[" 
	      << &aifunc2p_dist << ']'
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES 

  std::vector<T_REAL> lvector_distINSTiCLUSTER1k
    ( aomatrixt_centroids.getNumRows() );
  
  static std::uniform_real_distribution<T_REAL> lsuniformdis_real01(0.0,1.0);
  const INPUT_ITERATOR liiterator_instfirst  = aiiterator_instfirst;
  
  for (T_CLUSTERIDX *liter_gene = aiochrom_mutate.begin();
       aiiterator_instfirst != aiiterator_instlast;
       ++aiiterator_instfirst, liter_gene++ ) {

    if ( lsuniformdis_real01(gmt19937_eng) < airt_probMutation ) 
      { /*IF BEGIN PROBABILITY*/
	       
	/*Calculate cluster centers,cj's corresponding to sw:*/
	partition::PartitionLabel
	  <T_CLUSTERIDX>
	  lpartition_clusters
	  (aiochrom_mutate.getString(),
	   aiochrom_mutate.getStringSize(),
	   (T_CLUSTERIDX) aomatrixt_centroids.getNumRows()
	   );

	clusteringop::getCentroids
	  (aomatrixt_centroids,
	   aomatrixt_sumInstancesCluster,
	   aovectort_numInstancesInClusterK,
	   lpartition_clusters,
	   liiterator_instfirst,
	   aiiterator_instlast
	   );

	const T_FEATURE* linst_inter =
	  ((data::Instance<T_FEATURE>*) *aiiterator_instfirst)->getFeatures();
	       
	for (T_CLUSTERIDX li_j = 0; 
	     li_j < lpartition_clusters.getNumCluster();
	     li_j++) 
	  {
	    lvector_distINSTiCLUSTER1k[li_j] =
	      aifunc2p_dist
	      (aomatrixt_centroids.getRow(li_j),
	       linst_inter,
	       data::Instance<T_FEATURE>::getNumDimensions()
	       );
	  }
	/*d_sw(i) > 0 (Not singleton cluster) Incorrect is not = 0*/
	if ( lvector_distINSTiCLUSTER1k.at(*liter_gene) 
	     > T_REAL(0.0) )
	  {
	    auto&&  lvectorrt_probDistDistK = 
	      prob::getDistGKA
	      (lvector_distINSTiCLUSTER1k,
	       (T_REAL) 1.01 //Cm >= 1, tobe strictly greater than 1
	       );
	    T_CLUSTERIDX lmcidx_newAllen = 
	      gaselect::getIdxRouletteWheel
	      (lvectorrt_probDistDistK,
	       T_CLUSTERIDX(0)
	       );
	    *liter_gene = lmcidx_newAllen;
	  }
      } /*IF END PROBABILITY*/
  } //END FOR
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelFunc;
    lostrstream_labelFunc << lpc_labelFunc;
    aiochrom_mutate.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
    std::cout << std::endl;   
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES  
}

  
} /*END namespace gaintegerop*/

  
#endif /*GA_INTEGER_OPERATOR_HPP*/
