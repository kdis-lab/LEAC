/*! \file ga_real_operator.hpp
 *
 * \brief genetic real operators
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef GA_REAL_OPERATOR_HPP
#define GA_REAL_OPERATOR_HPP

#include <random>
#include "chromosome_string.hpp"
#include "linear_algebra_level1.hpp"

#include "verbose_global.hpp"


extern std::mt19937       gmt19937_eng;

/*! \namespace garealop
  \brief Genetic real operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace garealop {

  
/*! \fn void void heuristicCrossover(gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_childX, gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_childY, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentX, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentY
 )
  \brief Heuristic crossover \cite Chang:etal:GAclustering:GAGR:2009
  \details if \f$x\f$ is better than \f$y\f$ in term of the fitness \f$x' = x + r(x - y)\f$ and \f$y' =  x\f$ where \f$r = U(0,1)\f$, \f$U(0,1)\f$ is a uniform distribution on interval \f$[0,1]\f$ with mutation probability
  \param aochrom_childX a gaencode::ChromosomeString
  \param aochrom_childY a gaencode::ChromosomeString
  \param aichrom_parentX a gaencode::ChromosomeString
  \param aichrom_parentY a gaencode::ChromosomeString 
 */  
template <typename T_GENE,
	   typename T_METRIC
	   >
void
heuristicCrossover
(gaencode::ChromosomeString<T_GENE,T_METRIC>       &aochrom_childX,
 gaencode::ChromosomeString<T_GENE,T_METRIC>       &aochrom_childY,
 const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentX,
 const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentY
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::heuristicCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output gaencode::ChromosomeString: &aochrom_childX["
	      << &aochrom_childX << "]\n"
	      << " output gaencode::ChromosomeString: &aochrom_childY["
      	      << &aochrom_childY << "]\n"
	      << " input  gaencode::ChromosomeString: &aichrom_parentX[" << &aichrom_parentX << "]"
	      << "\tfitness = " << aichrom_parentX.getFitness() << '\n'
	      << " input  gaencode::ChromosomeString: &aichrom_parentY[" << &aichrom_parentY << "]"
	      << "\tfitness = " << aichrom_parentY.getFitness()
	      << "\n)"
	      << std::endl; 
  }
#endif //__VERBOSE_YES

  static std::uniform_real_distribution<T_METRIC> lsuniformdis_real01(0.0,1.0);

  T_METRIC lrt_r(lsuniformdis_real01(gmt19937_eng));
  
  if ( aichrom_parentX.getFitness() > aichrom_parentY.getFitness() ) { 
    /*x' = x
     */
    interfacesse::copy
      (aochrom_childX.getString(),
       aichrom_parentX.getString(), 
       aochrom_childX.getStringSize()
       );
    
    /*x' = x' + r(x' - y)
     */
  
    interfacesse::aysxpy
      (aochrom_childX.getString(),
       lrt_r,
       aichrom_parentY.getString(),
       aochrom_childX.getStringSize()
       );

    /*y' =  x
     */
    interfacesse::copy
      (aochrom_childY.getString(),
       aichrom_parentX.getString(), 
       aochrom_childY.getStringSize()
       );
   
  }
  else {

    /*x' = x
     */
    interfacesse::copy
      (aochrom_childX.getString(),
       aichrom_parentY.getString(), 
       aichrom_parentY.getStringSize()
       );
    
    /*x' = x' + r(x' - y)
     */
  
    interfacesse::aysxpy
      (aochrom_childX.getString(),
       lrt_r,
       aichrom_parentX.getString(),
       aochrom_childX.getStringSize()
       );

    /*y' =  x
     */
    interfacesse::copy
      (aochrom_childY.getString(),
       aichrom_parentY.getString(), 
       aochrom_childY.getStringSize()
       );
  }

  aochrom_childX.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aochrom_childX.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')' << " r = " << lrt_r << '\n' ;
    aochrom_childX.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aochrom_childY.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

}


/*! \fn void void averageCrossover(gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_childX, gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_childY, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentX, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentY, const T_METRIC airt_lambda                
 )
  \brief Heuristic crossover \cite Chang:etal:GAclustering:GAGR:2009
  \details if \f$x\f$ is better than \f$y\f$ in term of the fitness \f$x' = x + r(x - y)\f$ with mutation probability.
  \param aochrom_childX  a gaencode::ChromosomeString
  \param aochrom_childY  a gaencode::ChromosomeString
  \param aichrom_parentX a gaencode::ChromosomeString
  \param aichrom_parentY a gaencode::ChromosomeString 
  \param airt_lambda a real number 
 */  
template <typename T_GENE,
	   typename T_METRIC
	   >
void
averageCrossover
(gaencode::ChromosomeString<T_GENE,T_METRIC>       &aochrom_childX,
 gaencode::ChromosomeString<T_GENE,T_METRIC>       &aochrom_childY,
 const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentX,
 const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentY,
 const T_METRIC                                    airt_lambda                
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::averageCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output gaencode::ChromosomeString: &aochrom_childX["
	      << &aochrom_childX << "]\n"
              << "\n(output gaencode::ChromosomeString: &aochrom_childY["
	      << &aochrom_childY << "]\n"
	      << " input  gaencode::ChromosomeString: &aichrom_parentX[" << &aichrom_parentX << "]"
	      << "\tfitness = " << aichrom_parentX.getFitness() << '\n'
	      << " input  gaencode::ChromosomeString: &aichrom_parentY[" << &aichrom_parentY << "]"
	      << "\tfitness = " << aichrom_parentY.getFitness()
	      << "\n)"
	      << std::endl; 
  }
#endif //__VERBOSE_YES

  /*x' = y
   */
  interfacesse::copy
    (aochrom_childX.getString(),
     aichrom_parentY.getString(), 
     aochrom_childX.getStringSize()
     );
    
  /*x'  = y  - r(y -  x)
   */
  interfacesse::aysxpy
    (aochrom_childX.getString(),
     -airt_lambda,
     aichrom_parentX.getString(),
     aochrom_childX.getStringSize()
     );

  /*y' = x
   */
  interfacesse::copy
    (aochrom_childY.getString(),
     aichrom_parentX.getString(), 
     aochrom_childY.getStringSize()
     );
    
  /*y'  = x - r(x - y)
   */
  interfacesse::aysxpy
    (aochrom_childY.getString(),
     -airt_lambda,
     aichrom_parentY.getString(),
     aochrom_childY.getStringSize()
     );

  aochrom_childX.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aochrom_childX.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());

  aochrom_childY.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aochrom_childY.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aochrom_childX.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aochrom_childY.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

}


/*! \fn void void BLXalphaCrossover(gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_childX, gaencode::ChromosomeString<T_GENE,T_METRIC> &aochrom_childY, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentX, const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentY, const T_METRIC airt_alpha                
 )
  \brief BLXalphaCrossover
  \details 
  \param aochrom_childX a gaencode::ChromosomeString<T_GENE,T_METRIC>
  \param aochrom_childY a gaencode::ChromosomeString<T_GENE,T_METRIC>
  \param aichrom_parentX a gaencode::ChromosomeString<T_GENE,T_METRIC>
  \param aichrom_parentY a gaencode::ChromosomeString<T_GENE,T_METRIC>  
  \param airt_alpha a real number 
 */
template <typename T_GENE,
	   typename T_METRIC
	   >
void
BLXalphaCrossover
(gaencode::ChromosomeString<T_GENE,T_METRIC>       &aochrom_childX,
 gaencode::ChromosomeString<T_GENE,T_METRIC>       &aochrom_childY,
 const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentX,
 const gaencode::ChromosomeString<T_GENE,T_METRIC> &aichrom_parentY,
 const T_METRIC                                    airt_alpha                
 )
{  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "gaclusteringop::BLXalphaCrossover";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output gaencode::ChromosomeString: &aochrom_childX["
	      << &aochrom_childX << "]\n"
              << "\n(output gaencode::ChromosomeString: &aochrom_childY["
	      << &aochrom_childY << "]\n"
	      << " input  gaencode::ChromosomeString: &aichrom_parentX[" << &aichrom_parentX << "]"
	      << "\tfitness = " << aichrom_parentX.getFitness() << '\n'
	      << " input  gaencode::ChromosomeString: &aichrom_parentY[" << &aichrom_parentY << "]"
	      << "\tfitness = " << aichrom_parentY.getFitness()
	      << "\n)"
	      << std::endl; 
  }
#endif //__VERBOSE_YES
  /*
  const T_GENE* lstr_parentX; // =  aochrom_childX.getString(),
  const T_GENE* lstr_parentY; // =  aochrom_childY.getString(),
  T_GENE* lstr_childX; // =  aochrom_childX.getString(),
  T_GENE* lstr_childY; // =  aochrom_childY.getString(),
  */
  
  const T_GENE* lstr_parentX = aichrom_parentX.getString();
  const T_GENE* lstr_parentY = aichrom_parentY.getString();
  T_GENE* lstr_childX  = aochrom_childX.getString();
  T_GENE* lstr_childY  = aochrom_childY.getString();

  const T_GENE* lstr_end(lstr_parentX + aichrom_parentX.getStringSize());
    
  //for (uintidx lui_i = 0; lui_i <  aichrom_parentX.getStringSize(); lui_i++) {
  while ( lstr_parentX != lstr_end ) { 
    //const T_GENE lgeneX = aichrom_parentX.getGene(lui_i);
    //const T_GENE lgeneY = aichrom_parentY.getGene(lui_i);
    T_GENE lt_cmax = std::max(*lstr_parentX ,*lstr_parentY);
    T_GENE lt_cmin = std::min(*lstr_parentX ,*lstr_parentY);
    T_GENE lt_I = lt_cmax - lt_cmin;
    std::uniform_real_distribution<T_GENE>
      uniformdis_I(lt_cmin-lt_I*airt_alpha,lt_cmax+lt_I*airt_alpha);
    *lstr_childX = uniformdis_I(gmt19937_eng);
    *lstr_childY = uniformdis_I(gmt19937_eng);
    //aochrom_childX.setGene(lui_i,uniformdis_I(gmt19937_eng));
    //aochrom_childY.setGene(lui_i,uniformdis_I(gmt19937_eng));
    lstr_parentX++;
    lstr_parentY++;
    lstr_childX++;
    lstr_childY++;
  }
 
  aochrom_childX.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aochrom_childX.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());

  aochrom_childY.setFitness(-std::numeric_limits<T_METRIC>::max());  
  aochrom_childY.setObjetiveFunc(std::numeric_limits<T_METRIC>::max());

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aochrom_childX.print(std::cout,lpc_labelFunc);
    std::cout << '\n';
    aochrom_childY.print(std::cout,lpc_labelFunc);
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

}

   
/*! \fn void randomMutation(gaencode::ChromFixedLength<T_GENE,T_REAL>  &aiochromfixlen_toMutate)
  \brief randomMutation \cite Maulik:Bandyopadhyay:GAclustering:GAS:2000  
  \details Only one \f$g_i\f$ gene will be mutated. A number \f$\delta\f$ in the range [0,1] is generated with uniform distribution. If the value at that position is \f$g_i\f$, then after mutation it becomes
    \param aiochromfixlen_toMutate a gaencode::ChromFixedLength

    \f[ 
    \cases{
     g_i \cdot (1 \pm 2 \sigma ), &\quad $g_i \neq 0$ \cr 
     \pm 2 \sigma,                &\quad $g_i = 0$. \cr
     }
    \f]

    The ‘+’ or ‘−’ sign occurs with equal probability.

  \param aiochromfixlen_toMutate a gaencode::ChromFixedLength  chromosome to mutate 
 */  
template <typename T_GENE, 
          typename T_REAL
	  > 
void 
randomMutation(gaencode::ChromFixedLength<T_GENE,T_REAL>  &aiochromfixlen_toMutate)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "garealop::randomMutation";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n(output gaencode::ChromosomeString<>&: aiochromfixlen_toMutate["
	      << &aiochromfixlen_toMutate << "]"
	      << '\n';
  }
#endif //__VERBOSE_YES

  static std::uniform_int_distribution<uintidx> lsuniformdis_ui0KD
    (0,gaencode::ChromFixedLength<T_GENE,T_REAL>::stcgetStringSize()-1);
  static std::uniform_real_distribution<T_REAL> lsuniformdis_real01(0.0, 1.0);

  const uintidx  liui_randPositionGene(lsuniformdis_ui0KD(gmt19937_eng));

  T_GENE lt_gene = aiochromfixlen_toMutate.getGene(liui_randPositionGene);
  
#if  DATATYPE_CENTROIDS_ROUND == 0
  const T_GENE  lrt_sign = (lsuniformdis_real01(gmt19937_eng) < 0.5)?T_GENE(1):T_GENE(-1);
  const T_GENE  lrt_delta = T_GENE(lsuniformdis_real01(gmt19937_eng));
  
  const T_GENE lt_newGene = (lt_gene != 0.0)?
    lt_gene  + lrt_sign * 2.0 * lrt_delta  * lt_gene
    :lt_gene + lrt_sign * 2.0 * lrt_delta;

  aiochromfixlen_toMutate.setGene(liui_randPositionGene,lt_newGene);
  
#else
  const DATATYPE_REAL lrt_sign = (lsuniformdis_real01(gmt19937_eng) < 0.5)?
    DATATYPE_REAL(1):DATATYPE_REAL(-1);
  const DATATYPE_REAL lrt_delta = DATATYPE_REAL(lsuniformdis_real01(gmt19937_eng));
  const T_GENE lt_newGene =  (lt_gene != 0)?
    T_GENE( std::round((DATATYPE_REAL) lt_gene + lrt_sign * 2.0 * lrt_delta  * lt_gene))
    :T_GENE(std::round((DATATYPE_REAL) lt_gene + lrt_sign * 2.0 * lrt_delta));

  aiochromfixlen_toMutate.setGene(liui_randPositionGene,lt_newGene);
  
#endif //DATATYPE_CENTROIDS_ROUND

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelFunc;
    lostrstream_labelFunc
      << lpc_labelFunc
      << ",gene position," << liui_randPositionGene
      << "v," << lt_gene
      << ",sign," << lrt_sign * 2
      << ",delta ," << lrt_delta 
      << ",mutation gene,"
      << aiochromfixlen_toMutate.getGene(liui_randPositionGene)
      << ",";
    aiochromfixlen_toMutate.print(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
} /*gaencode::randomMutation*/



  
/*! \fn void muhlenbeinMutation(gaencode::ChromosomeString<T_GENE,T_REAL>  &aiochromstr_toMutate, T_GENE aifeact_ai, T_GENE aifeact_bi)
  \brief Muhlenbein mutation
  \details 
  \param aiochromstr_toMutate a gaencode::ChromosomeString
  \param aifeact_ai lower range limit 
  \param aifeact_bi upper range limit 
*/
template <typename T_GENE, 
          typename T_REAL
	  > 
void
muhlenbeinMutation
(gaencode::ChromosomeString<T_GENE,T_REAL>  &aiochromstr_toMutate,
 T_GENE                                     aifeact_ai,
 T_GENE                                     aifeact_bi
)
{
  static std::uniform_real_distribution<T_REAL> uniformdis_real01(0, 1);
  // 0    1    2      3     4         5         6         7
  static T_GENE larray_pow2[]
    = { 1.0, 0.5, 0.25, 0.125, 6.25e-2, 3.125e-2, 1.5625e-2, 7.8125e-3,
	 // 8             9           10            11          12
	3.90625e-3, 1.953125e-3, 9.765625e-4, 4.8828125e-4, 2.44140625e-4,
	// 13        14         15
	1.220703125e-4, 6.103515625e-5, 3.0517578125e-5};
  
   
  const T_GENE lrt_rangi = 0.1*(aifeact_bi-aifeact_ai);
  T_GENE *lrt_gene =  aiochromstr_toMutate.getString();
  const T_GENE *lrt_gene_end(lrt_gene + aiochromstr_toMutate.getStringSize());

  while ( lrt_gene  != lrt_gene_end ) {
    T_GENE lrt_gamma = 0.0;
    for ( uintidx lui_i = 0; lui_i < 16; lui_i++){
      T_GENE lrt_alpha =  (uniformdis_real01(gmt19937_eng) <= 0.0625)?1.0:0.0;
      lrt_gamma += lrt_alpha * larray_pow2[lui_i];
    }

    T_GENE lrt_sign =  (uniformdis_real01(gmt19937_eng) < 0.5)?-1.0:1.0;

    *lrt_gene +=  lrt_sign * lrt_rangi * lrt_gamma;

    ++lrt_gene; 
  }				
}

  
} /*END namespace garealop*/



#endif /*GA_REAL_OPERATOR_HPP*/
