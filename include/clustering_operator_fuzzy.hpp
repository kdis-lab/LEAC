/*! \file clustering_operator_fuzzy.hpp
 *
 * \brief clustering operator fuzzy  
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __CLUSTERING_OPERATOR_FUZZY_HPP
#define __CLUSTERING_OPERATOR_FUZZY_HPP

#include "random_ext.hpp"
#include "matrix_operation.hpp"
#include "verbose_global.hpp"


extern StdMT19937 gmt19937_eng;

/*! \namespace clusteringop
  \brief Clustering operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace clusteringop {


/*! \fn void randomInitialize(mat::MatrixRow<T_REAL>  &aomatrixt_fuzzyPartition)
    \brief Initializes randomly a fuzzy partition
    \details
    \param aomatrixt_fuzzyPartition a random fuzzy partition
 */
template < typename T_REAL >
void
randomInitialize
(mat::MatrixRow<T_REAL>  &aomatrixt_fuzzyPartition) 
{
  T_REAL  li_anc;
  T_REAL  lT_sum;
  T_REAL  lT_uki;
  T_REAL  lT_random;

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::randomInitialize";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
    	      << "\n(output Matrix: aomatrixt_fuzzyPartition[" << &aomatrixt_fuzzyPartition << "]\n"
	      << ")\n";
  }
#endif /*__VERBOSE_YES*/

  static std::uniform_real_distribution<T_REAL> lsuniformdis_real01(0.0,1.0);
  
  //i \in {1,2,..,n}
  for (uintidx lui_i = 0; lui_i < aomatrixt_fuzzyPartition.getNumColumns(); lui_i++) { 
    lT_sum = 1.0;
    //j \in {1,2,..,k}
    for (uintidx luintidx_j = 0; luintidx_j < (aomatrixt_fuzzyPartition.getNumRows()-1); luintidx_j++) {
      lT_random = lsuniformdis_real01(gmt19937_eng); 
      lT_random /= 2.0;
      li_anc = T_REAL( aomatrixt_fuzzyPartition.getNumRows() - (luintidx_j+1));
      lT_uki = lT_sum * (1.0 - std::pow(lT_random, 1.0 / li_anc));
      lT_sum -= lT_uki; 
      aomatrixt_fuzzyPartition(luintidx_j,lui_i) = lT_uki;
    } /*END for luintidx_j*/
    aomatrixt_fuzzyPartition(aomatrixt_fuzzyPartition.getNumRows()-1, lui_i) =  lT_sum;
  } /*END for lui_i*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelCrisp;
    lostrstream_labelCrisp << "<CRISPMATRIX:" << lpc_labelFunc;
    aomatrixt_fuzzyPartition.print(std::cout,lostrstream_labelCrisp.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
}
  
/*! \fn  mat::MatrixRow<T_REAL> fuzzyPartition(const mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, T_REAL airt_m, const dist::Dist<T_REAL,T_FEATURE>  &aifunc2p_dist)
    \brief Get a \f$ U_{x \times n } = [u_{ji}] \f$ fuzzy c-partitions based on \cite Bezdek:ClusterAnalysis:FCM:1974 \cite Bezdek:etal:ClusterAnalysis:FCM:1984
    \details 
   \f[
    u_{ji} = \left( \sum_{j=1}^k \left( \frac{ D(\mu_j'-x_i) }{ D(\mu_j-x_i) } \right)^{\frac{2}{m-1}} \right)^{-1}
   \f]
    \param aimatrixt_centroids a mat::MatrixRow with \f$ \mu_j \f$ centroids of the clusters
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param airt_m a real number with the  \f$ m \f$  weighting exponent
    \param aifunc2p_dist an object of type dist::Dist square distance of preference
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE,
           typename T_REAL
	   >
mat::MatrixRow<T_REAL>
fuzzyPartition
(const mat::MatrixRow<T_FEATURE>       &aimatrixt_centroids,
 INPUT_ITERATOR                        aiiterator_instfirst,
 const INPUT_ITERATOR                  aiiterator_instlast,
 T_REAL                                airt_m,
 const dist::Dist<T_REAL,T_FEATURE>    &aifunc2p_distSquare
 )
{
  /*Change 2 by 1 of the original formula
    considering that the distance is square
  */
  
  const T_REAL lrt_p = (1.0 /(airt_m - 1.0));
 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::fuzzyPartition";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n"
              << "(input  mat::MatrixRow<T_FEATURE>&: aimatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "input const aiiterator_instlast[" << *aiiterator_instlast << "]\n"
	      << "input  T_REAL airt_m = " << airt_m << '\n'
	      << "lrt_p = "  <<  lrt_p << '\n'
	      << "aifunc2p_distSquare\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  const uintidx  lui_numInstances = uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast));

  mat::MatrixRow<T_REAL>
    aomatrixt_u
    (aimatrixt_centroids.getNumRows(),
     lui_numInstances
     );
    
  T_REAL lrt_Dik;
  T_REAL lrt_Djk;
  T_REAL lrt_DikU;
  T_REAL lrt_sumDjkU;
  

  /*j \in {1, 2, ...,k}*/
  for ( uintidx lui_j = 0; lui_j < aimatrixt_centroids.getNumRows(); lui_j++)  {
    
    T_REAL*   lmrt_uRow_j     = aomatrixt_u.getRow(lui_j);
    const T_FEATURE*  lmrt_centroid_j = aimatrixt_centroids.getRow(lui_j); 

    /*i \in {1, 2, ...,n}*/
    uintidx lui_i = 0;
    for ( INPUT_ITERATOR liiterator_instfirst = aiiterator_instfirst;
	 liiterator_instfirst != aiiterator_instlast;
	 ++liiterator_instfirst, lui_i++)
      {

      T_FEATURE* linst_features_i =
	((data::Instance<T_FEATURE>*) *liiterator_instfirst)->getFeatures();
      
      lrt_Dik = aifunc2p_distSquare
	(lmrt_centroid_j,
	 linst_features_i,
	 data::Instance<T_FEATURE>::getNumDimensions()
	 );
  
      lrt_DikU = 1.0 / std::pow(lrt_Dik, lrt_p);
      lrt_sumDjkU = 0.0;
      /*l \in {1, 2, ...,k}*/
      for ( uintidx lui_l = 0; lui_l < aimatrixt_centroids.getNumRows(); lui_l++)  { 
	const T_FEATURE*  lmrT_centroid_l = aimatrixt_centroids.getRow(lui_l); 
	
        lrt_Djk =
	  aifunc2p_distSquare
	  (lmrT_centroid_l,
	   linst_features_i,
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );
	lrt_sumDjkU += 1.0 / std::pow(lrt_Djk, lrt_p);
      } /*END for lui_l*/
      lmrt_uRow_j[lui_i] = lrt_DikU / lrt_sumDjkU;
    } /*END for lui_i*/
  } /*END for lui_j*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelU;
    lostrstream_labelU << "<U:" << lpc_labelFunc;
    aomatrixt_u.print(std::cout,lostrstream_labelU.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return aomatrixt_u;

}

/*getCentroids
  \cite{Bezdek:ClusterAnalysis:FCM:1974}
  \cite{Bezdek:etal:ClusterAnalysis:FCM:1984}
*/
template < typename T_REAL,
	   typename T_FEATURE
	   >
bool 
getCentroids
(mat::MatrixRow<T_FEATURE>  &aomatrixt_centroids, 
 mat::MatrixRow<T_REAL>     &aimatrixt_fuzzyPartition,
 mat::MatrixRow<T_FEATURE>  &aimatrixt_instances, 
 T_REAL                     ait_m,
 mat::MatrixRow<T_REAL>     &aimatrixrowt_work
 )
{
  bool           lob_isValidCentroids;           
 
  lob_isValidCentroids = true;

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "getCentroids:  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t(output bool: lob_isValidCentroids = " 
	      << lob_isValidCentroids << '\n'
	      << "\t output mat::MatrixRow<T_FEATURE>&: aomatrixt_centroids[" 
	      << &aomatrixt_centroids << "]\n"
	      << "\t input  mat::MatrixRow<T_REAL>& aimatrixt_fuzzyPartition["
	      <<  &aimatrixt_fuzzyPartition << "]\n"
	      << "\t input  mat::MatrixRow<T_FEATURE>&: aimatrixt_instances[" 
	      << &aimatrixt_instances << "]\n"
	      << "\t input  T_REAL ait_m = " << ait_m << '\n'
	      << "\t input  mat::MatrixRow<T_REAL>& aimatrixrowt_work[" 
	      <<  &aimatrixrowt_work << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES*/
  
  /* Work_ij = U_ij^m */

 T_REAL* larray_U =  aimatrixt_fuzzyPartition.toArray();
 uintidx luintidx_l = 0;
  std::for_each
    (aimatrixrowt_work.toArray(), 
     aimatrixrowt_work.toArray() + aimatrixrowt_work.getNumElems(),
     [&] (T_REAL &liter_wij) 
     { 
       liter_wij = std::pow(larray_U[luintidx_l++],ait_m); 
     } 
     );

  
#if  DATATYPE_CENTROIDS_ROUND == 0

  /* V =  Work  * Y  (BEST PERFORMANSE T_FEATURE == T_REAL)*/

  interfacesse::gemm
    (aomatrixt_centroids, /*CENTROIDS k x d*/
     aimatrixrowt_work,   /* k x Inst */
     aimatrixt_instances  /* Inst * d */
     );

  for ( uintidx lui_i = 0; lui_i < aomatrixt_centroids.getNumRows(); lui_i++) {
 
    T_REAL lrt_sumWi  = 
      interfacesse::sum 
      (aimatrixrowt_work.getRow(lui_i),
       aimatrixrowt_work.getNumColumns()
       );

    if ( lrt_sumWi == T_REAL(0.0) ) {
      lob_isValidCentroids = false;

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( (geiinparam_verbose <= geiinparam_verboseMax) && lob_isValidCentroids == false ) {
	std::cout << "\nERROR: getCentroids: "
		  << "solution is not valid"  
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES */

    }
    else {
      T_REAL lrt_sumWiInv = 1.0 / lrt_sumWi;  
      interfacesse::scal
	(aomatrixt_centroids.getRow(lui_i),
	 lrt_sumWiInv,
	 aomatrixt_centroids.getNumColumns()
	 );
    }
  } /*END for lui_i */

#else /* DATATYPE_CENTROIDS_ROUND == 0*/
  
  mat::MatrixRow<T_REAL>         *lmatrixPFloatT_tmpV;
 
  lmatrixPFloatT_tmpV = new mat::MatrixRow<T_REAL>(aomatrixt_centroids.getNumRows(), aomatrixt_centroids.getNumColumns() );

  /*V =  Work  * Y  
  */
  interfacesse::gemm
    (*lmatrixPFloatT_tmpV,
     T_REAL(1.0),
     mat::CblasNoTrans,
     aimatrixrowt_work, 
     mat::CblasNoTrans,
     aimatrixt_instances, 
     T_REAL(0.0)
     );

  for ( uintidx lui_i = 0; lui_i < aomatrixt_centroids.getNumRows(); lui_i++) { 
    T_REAL lrt_sumWi  = 
      interfaceclapack_sum 
      (aimatrixrowt_work.getRow(lui_i),
       aimatrixrowt_work.getNumColumns()
       );

    if ( lrt_sumWi == T_REAL(0.0) ) {
      lob_isValidCentroids = false;

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( (geiinparam_verbose <= geiinparam_verboseMax) && lob_isValidCentroids == false ) {
	std::cout << "\nERROR: getCentroids: "
		  << "solution is not valid"  
		  << std::endl;
      }
      --geiinparam_verbose;
#endif /*__VERBOSE_YES */

    }
    else {
      T_REAL lrt_sumWiInv = 1.0 / lrt_sumWi;  
      interfaceclapack_scal
	(lmatrixPFloatT_tmpV->getRow(lui_i),
	 lrt_sumWiInv,
	 lmatrixPFloatT_tmpV->getNumColumns()
	 );

       T_REAL* larrayPFloat_V =  lmatrixPFloatT_tmpV->toArray();
       uintidx luintidx_ll = 0;
       std::for_each
	 (aomatrixt_centroids.toArray(), 
	  aomatrixt_centroids.toArray() + aomatrixt_centroids.getNumElems(),
	  [&] (T_FEATURE &liter_vij) 
	  { 
	    liter_vij = (T_FEATURE) std::round(larrayPFloat_V[luintidx_ll++]); 
	  } 
	  ); 
   
    }
  } /*END for luintidx_l*/

  delete lmatrixPFloatT_tmpV;

#endif /*DATATYPE_CENTROIDS_ROUND*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "getCentroids: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output bool: lob_isValidCentroids = " << lob_isValidCentroids << '\n'
	      << "output mat::MatrixRow<T_FEATURE>&: aomatrixt_centroids[" << &aomatrixt_centroids << "]\n";
    aomatrixt_centroids.print();
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/


  return lob_isValidCentroids; 
}



/*fcmbezdek_updateMembership:
  Update membership
 
  \cite{Bezdek:ClusterAnalysis:FCM:1974}
  \cite{Bezdek:etal:ClusterAnalysis:FCM:1984}

  new ahora la funcion es fuzzyPartition
*/

template < typename T_REAL,
	   typename T_FEATURE
	 >
void
reassignCluster
(mat::MatrixRow<T_REAL>       &aomatrixt_fuzzyPartition,
 mat::MatrixRow<T_FEATURE>    &aimatrixt_instances, 
 mat::MatrixRow<T_FEATURE>    &aimatrixt_centroids,
 T_REAL                       ait_m,
 dist::Dist<T_REAL,T_FEATURE> &aifunc2p_squaredDist
 )
{  

  T_REAL lT_Dik;
  T_REAL lT_Djk;
  T_REAL  lT_DikU;
  T_REAL lT_sumDjkU;
  //T_REAL lT_p;
  uintidx lst_numDimensions;


#ifdef __VERBOSE_YES
   const char* lpc_labelFunc = "clusteringop::reassignCluster";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ')'
  	      << "\t(output mat::MatrixRow<T_REAL>&: aomatrixt_fuzzyPartition[" 
	      << &aomatrixt_fuzzyPartition << "]\n"
              << "\t input  mat::MatrixRow<T_FEATURE>&: aimatrixt_instances[" 
	      << &aimatrixt_instances << "]\n"
	      << "\t input  mat::MatrixRow<T_FEATURE>&: aomatrixt_centroids[" 
	      << &aimatrixt_centroids << "]\n"
	      << "\t input  T_REAL ait_m = " << ait_m << '\n'
	      << "\t aifunc2p_squaredDist\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES*/

  /*Change 1 by 2 of the original formula
   considering that the distance is square
  */
  
  const T_REAL lrt_p = (1.0 /(ait_m - 1.0)); 
  /*i \in {1, 2, ...,k}*/
  lst_numDimensions = aimatrixt_centroids.getNumColumns(); 
  for ( uintidx lui_i = 0; lui_i < aimatrixt_centroids.getNumRows(); lui_i++)  { 

    T_REAL*     lmrT_uRow_i     = aomatrixt_fuzzyPartition.getRow(lui_i);
    T_FEATURE*  lmrT_centroid_i = aimatrixt_centroids.getRow(lui_i); 

    /*j \in {1, 2, ...,n}*/
    for ( uintidx luintidx_j = 0; luintidx_j < aimatrixt_instances.getNumRows(); luintidx_j++)  { 
      T_FEATURE* lmrT_instance_j = aimatrixt_instances.getRow(luintidx_j); 
      lT_Dik = aifunc2p_squaredDist(lmrT_centroid_i,lmrT_instance_j,lst_numDimensions);
      lT_DikU = 1.0 / std::pow(lT_Dik, lrt_p);
      lT_sumDjkU = 0.0;
      
      /*l \in {1, 2, ...,k}*/
      for ( uintidx luintidx_l = 0; luintidx_l < aimatrixt_centroids.getNumRows(); luintidx_l++)  { 
	T_FEATURE*  lmrT_centroid_l = aimatrixt_centroids.getRow(luintidx_l); 
	//lT_Djk = normt_squareAdistance<T_NORM_CENTROIDINSTANCES,T_FEATURE,T_CENTROIDS,T_WEIGHT_MATRIX>
        //(lmrT_instance_j, lmrT_centroid_l, aimatrixrowt_weightMatrixA);
        lT_Djk = aifunc2p_squaredDist(lmrT_centroid_l,lmrT_instance_j,lst_numDimensions);
	lT_sumDjkU += 1.0 / std::pow(lT_Djk, lrt_p);
      } /*END for luintidx_l*/
      lmrT_uRow_i[luintidx_j] = lT_DikU / lT_sumDjkU;
      //matrix_setElem(aomatrixt_fuzzyPartition, lui_i, luintidx_j, lT_Dik / lT_sumDjk);
    } /*END for luintidx_j*/
  } /*END for lui_i*/

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "output mat::MatrixRow<T_REAL>& aomatrixt_fuzzyPartition["
	      << &aomatrixt_fuzzyPartition << "]\n"
	      << aomatrixt_fuzzyPartition;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
}
  
} /*END namespace clusteringop*/
  
#endif /*__CLUSTERING_OPERATOR_FUZZY_HPP*/
