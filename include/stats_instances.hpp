/*! \file stats_instances.hpp
 *
 * \brief Statistics for instances
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef  __STATS_INSTANCES_HPP
#define  __STATS_INSTANCES_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <utility>      // std::pair, std::make_pair
#include "container_instance.hpp"
#include "linear_algebra_level1.hpp"
#include "linear_algebra_level2.hpp"
#include "matrix_operation.hpp"
#include  "common.hpp"

#include "verbose_global.hpp"

/*! \namespace stats
  \brief Model for transforming the data structure of the instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/


namespace stats {

   
/*! \fn void meanVector(T_FEATURE *aoarrayt_mean, const T_FEATURE_NUM  ait_feactureNum, const T_FEATURE_SUM  *aiarrayt_feactureSum, const uintidx aiu_idxVector = 0)
   \brief ￼Calculate average of a vector
   \details For a sum of feactures, calculate the average 
   \param aoarrayt_mean an array to store the average result
   \param ait_feactureNum an integer with the number of feactures
   \param aiarrayt_feactureSum an array with the sum of feactures
   \param aiu_idxVector an identifier of the vector number, for control
*/
template < typename T_FEATURE, 
	   typename T_FEATURE_SUM,
	   typename T_FEATURE_NUM
	   >
void 
meanVector
(T_FEATURE            *aoarrayt_mean,
 const T_FEATURE_NUM  ait_feactureNum,
 const T_FEATURE_SUM  *aiarrayt_feactureSum,
 const uintidx        aiu_idxVector  = 0
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::meanVector";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(output T_FEATURE *aoarrayt_mean[" << aoarrayt_mean << "]\n"
              << " input  const T_FEATURE_NUM  ait_feactureNum = " 
	      << ait_feactureNum << "\n"
	      << " input  T_FEATURE_SUM *aiarrayt_feactureSum[" 
	      << aiarrayt_feactureSum << "]\n";
    
    std::ostringstream lostrstream_labelSum;
    lostrstream_labelSum << "<SUMFEACTURES:" << lpc_labelFunc;
    
    inout::containerprint
      (aiarrayt_feactureSum,
       aiarrayt_feactureSum + data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelSum.str().c_str(),
       ','
       );
    std::cout << "\ninput  const uintidx aiu_idxVector = "
	      << aiu_idxVector << "\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/


  interfacesse::copya
    (aoarrayt_mean,
     T_FEATURE(0),
     data::Instance<T_FEATURE>::getNumDimensions()
     );

  interfacesse::axpyInv
    (aoarrayt_mean,
     ait_feactureNum,
     aiarrayt_feactureSum,
     data::Instance<T_FEATURE>::getNumDimensions()
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";

    std::ostringstream lostrstream_labelMean;
    lostrstream_labelMean << "<MEANFEACTURES:" << lpc_labelFunc;
    
    inout::containerprint
      (aoarrayt_mean,
       aoarrayt_mean + data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelMean.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}


/*! \fn void sumFeactures (T_FEATURE_SUM *aoarrayt_sumFeactures, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const T_FEATURE ait_initialValue)
   \brief Add the instances by attributes
   \details
   \param aoarrayt_sumFeactures an output vector with the sum of the instances 
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
   \param ait_initialValue value with which the output vector is initialized
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE_SUM,
	   typename T_FEATURE
	   >
void 
sumFeactures
(T_FEATURE_SUM          *aoarrayt_sumFeactures,
 INPUT_ITERATOR         aiiterator_instfirst,
 const INPUT_ITERATOR   aiiterator_instlast,
 const T_FEATURE        ait_initialValue 
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::sumFeactures";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  interfacesse::copya
    (aoarrayt_sumFeactures,
     T_FEATURE_SUM(ait_initialValue),
     data::Instance<T_FEATURE>::getNumDimensions()
     );

  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {

    const data::Instance<T_FEATURE>* linst_inter =
      (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    
    interfacesse::axpy
      (aoarrayt_sumFeactures,
       T_FEATURE(1),
       linst_inter->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
  }
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    
    std::ostringstream lostrstream_labelSum;
    lostrstream_labelSum << "<SUMFEACTURES:" << lpc_labelFunc;
      
    inout::containerprint
      (aoarrayt_sumFeactures,
       aoarrayt_sumFeactures+data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelSum.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES  
}


/*! \fn uintidx sumFeacturesFreq (T_FEATURE_SUM *aoarrayt_sumFeactures, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const T_FEATURE ait_initialValue, const FUNCTIONFREQ func_freq
 )
   \brief Add the instances with frequency by attributes
   \details 
   \param aoarrayt_sumFeactures an output vector with the sum of the instances 
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
   \param ait_initialValue value with which the output vector is initialized
   \param func_freq a function that gets the frequency of the instance that is adding up
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE_SUM,
	   typename T_FEATURE,
	   typename FUNCTIONFREQ
	   >
uintidx
sumFeacturesFreq
(T_FEATURE_SUM         *aoarrayt_sumFeactures,
 INPUT_ITERATOR        aiiterator_instfirst,
 const INPUT_ITERATOR  aiiterator_instlast,
 const T_FEATURE       ait_initialValue,
 const FUNCTIONFREQ    func_freq
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::sumFeacturesFreq";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  typedef decltype(func_freq(*aiiterator_instfirst)) FreqType;

  uintidx luintidx_numInstances = 0;
  
  interfacesse::copya
    (aoarrayt_sumFeactures,
     T_FEATURE_SUM(ait_initialValue),
     data::Instance<T_FEATURE>::getNumDimensions()
     );

  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    
    const data::Instance<T_FEATURE>*
      linst_inter =  static_cast<data::Instance<T_FEATURE>*>(*aiiterator_instfirst);
    
    FreqType lit_freqInst = func_freq(*aiiterator_instfirst);
    
    luintidx_numInstances += lit_freqInst;
    interfacesse::axpy
      (aoarrayt_sumFeactures,
#if  DATATYPE_CENTROIDS_ROUND == 0
       (T_FEATURE_SUM) lit_freqInst,
#else
       lit_freqInst,
#endif //DATATYPE_CENTROIDS_ROUND
       linst_inter->getFeatures(),
       data::Instance<T_FEATURE>::getNumDimensions()
       );
  }
  

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
	      << " luintidx_numInstances = " <<  luintidx_numInstances << '\n';
    
    std::ostringstream lostrstream_labelSum;
    lostrstream_labelSum << "<SUMFEACTURES:" << lpc_labelFunc;
      
    inout::containerprint
      (aoarrayt_sumFeactures,
       aoarrayt_sumFeactures+ data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelSum.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return luintidx_numInstances;
}


/*! \fn std::vector<T_FEATURE> meanRow (const mat::MatrixRow<T_FEATURE>& aimatrixrowt_instances)
   \brief  Obtiene un vector con la media de una matriz  por columnas
   \details Get a vector with the mean of one matrix per row
   \param aimatrixrowt_instances a matrix with the data to average
   \return vector with the mean by row
*/
template < typename T_FEATURE, 
	   typename T_FEATURE_SUM
	   >
std::vector<T_FEATURE>
meanRow
(const mat::MatrixRow<T_FEATURE>& aimatrixrowt_instances)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::meanRow";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input const mat::MatrixRow<T_FEATURE>& aimatrixrowt_instances["
	      << &aimatrixrowt_instances << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  std::vector<T_FEATURE>     
    lovectort_mean
    (aimatrixrowt_instances.getNumRows());

  for (uintidx luintidx_i = 0; luintidx_i < aimatrixrowt_instances.getNumRows(); ++luintidx_i) {
    
    T_FEATURE_SUM lt_sumRow = 
    interfacesse::sum
      (aimatrixrowt_instances.getRow(luintidx_i),
       aimatrixrowt_instances.getNumColumns()
       );
    
#if  DATATYPE_CENTROIDS_ROUND == 0
    
    lovectort_mean[luintidx_i] = lt_sumRow / (T_FEATURE) aimatrixrowt_instances.getNumColumns();

#else
    
    DATATYPE_REAL lt_tmp =
      (DATATYPE_REAL) lt_sumRow / (DATATYPE_REAL) aimatrixrowt_instances.getNumColums();
    lovectort_mean[luintidx_i] = (T_FEATURE) std::round(lt_tmp);
    
#endif /*DATATYPE_CENTROIDS_ROUND*/
    
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMeanFeatures;
    lostrstream_labelMeanFeatures
      << "<MEANROWFEATURES:" << lpc_labelFunc
      << "lovectort_mean[" << &lovectort_mean << "]\n";
    inout::containerprint
      (lovectort_mean.begin(),
       lovectort_mean.end(),
       std::cout,
       lostrstream_labelMeanFeatures.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lovectort_mean; 
}

/*! \fn void sumFeacturesSQ (T_FEATURE_SUM *aoarrayt_sumFeacturesSQ, T_FEATURE *aiarrayt_meanFeactures, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR  aiiterator_instlast)
   \brief Gets the sum of the difference of the instances with the mean and the result square
   \details 
   \param aoarrayt_sumFeacturesSQ an output vector with the sum
   \param aiarrayt_meanFeactures an input vector with the mean of the instances
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE_SUM,
	   typename T_FEATURE
	   >
void 
sumFeacturesSQ
(T_FEATURE_SUM         *aoarrayt_sumFeacturesSQ,
 T_FEATURE             *aiarrayt_meanFeactures,
 INPUT_ITERATOR        aiiterator_instfirst,
 const INPUT_ITERATOR  aiiterator_instlast
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::sumFeacturesSQ";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
  interfacesse::copya
    (aoarrayt_sumFeacturesSQ,
     T_FEATURE_SUM(0),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
 
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    data::Instance<T_FEATURE>* linst_inter =
      (data::Instance<T_FEATURE>*)(*aiiterator_instfirst);
    const T_FEATURE  *larrayT_feacturesi = linst_inter->getFeatures();
    for ( uintidx  li_j = 0; li_j < data::Instance<T_FEATURE>::getNumDimensions(); li_j++) {
      T_FEATURE lt_tmp =  larrayT_feacturesi[li_j] - aiarrayt_meanFeactures[li_j];
      aoarrayt_sumFeacturesSQ[li_j] += T_FEATURE_SUM(lt_tmp * lt_tmp);
    }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    
    std::ostringstream lostrstream_labelSum;
    lostrstream_labelSum << "<SUMSQFEACTURES:" << lpc_labelFunc;
      
    inout::containerprint
      (aoarrayt_sumFeacturesSQ,
       aoarrayt_sumFeacturesSQ+data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelSum.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES 

}


/*! \fn void sumFeacturesFreqSQ(T_FEATURE_SUM *aoarrayt_sumFeacturesSQ, T_FEATURE *aiarrayt_meanFeactures, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR  aiiterator_instlast, const FUNCTIONFREQ func_freq)
   \brief Gets the sum of the difference of the instances with the mean and the result square. For instances frequently
   \details 
   \param aoarrayt_sumFeacturesSQ an output vector with the sum
   \param aiarrayt_meanFeactures an input vector with the mean of the instances
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
   \param func_freq a function that gets the frequency of the instance that is adding up
*/
template < typename INPUT_ITERATOR,
	   typename T_FEATURE_SUM,
	   typename T_FEATURE,
	   typename FUNCTIONFREQ
	   >
void 
sumFeacturesFreqSQ
(T_FEATURE_SUM         *aoarrayt_sumFeacturesSQ,
 T_FEATURE             *aiarrayt_meanFeactures,
 INPUT_ITERATOR        aiiterator_instfirst,
 const INPUT_ITERATOR  aiiterator_instlast,
 const FUNCTIONFREQ    func_freq
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::sumFeacturesFreqSQ";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  interfacesse::copya
    (aoarrayt_sumFeacturesSQ,
     T_FEATURE_SUM(0),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
 
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    data::Instance<T_FEATURE>* linst_inter =
      (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
    const T_FEATURE  *larrayT_feacturesi =  linst_inter->getFeatures();
    const T_FEATURE lit_freqInst = (T_FEATURE) func_freq(*aiiterator_instfirst);
    for ( uintidx li_j = 0; li_j < data::Instance<T_FEATURE>::getNumDimensions(); li_j++)
      {
      T_FEATURE lt_tmp =
	lit_freqInst * larrayT_feacturesi[li_j] - aiarrayt_meanFeactures[li_j];
      aoarrayt_sumFeacturesSQ[li_j] += T_FEATURE_SUM(lt_tmp * lt_tmp);
    }
  }
    
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    
    std::ostringstream lostrstream_labelSum;
    lostrstream_labelSum << "<SUMSQFEACTURES:" << lpc_labelFunc;
      
    inout::containerprint
      (aoarrayt_sumFeacturesSQ,
       aoarrayt_sumFeacturesSQ+data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelSum.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES 

}


/*! \fn void varTodesvstd(T_FEATURE* aioarrayt_varTodesvstd)
   \brief Obtain the standard deviation given the variance
   \details
   \param aioarrayt_varTodesvstd an input and output vector with the variances of the attributes of the instances
*/
template < typename T_FEATURE > 
void
varTodesvstd
(T_FEATURE* aioarrayt_varTodesvstd)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::varTodesvstd";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
    	      << "(input  aioarrayt_varTodesvstd[" <<  aioarrayt_varTodesvstd << "]\n"
	      << ")"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  
std::transform
    (aioarrayt_varTodesvstd,
     aioarrayt_varTodesvstd+data::Instance<T_FEATURE>::getNumDimensions(),
     aioarrayt_varTodesvstd,
     [](T_FEATURE lt_dim)
     {
       
#if  DATATYPE_CENTROIDS_ROUND == 0

       return std::sqrt(lt_dim);

#else

       DATATYPE_REAL lt_tmp = 
	 (DATATYPE_REAL) lt_dim;
       lt_tmp = std::sqrt(lt_tmp);
    
       return  (T_FEATURE) std::round(lt_tmp);

#endif //DATATYPE_CENTROIDS_ROUND

     }
     );

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelDesvStdFeatures;
    lostrstream_labelDesvStdFeatures << "<DESVSTDINSTANCE:" << lpc_labelFunc;
    inout::containerprint
      (aioarrayt_varTodesvstd,
       aioarrayt_varTodesvstd+data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelDesvStdFeatures.str().c_str(),
       ','
       );
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
}



/*! \fn void minFeatures (T_FEATURE *aoarrayt_minFeactures, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR  aiiterator_instlast)
   \brief Obtain the minimum attributes of the instances
   \details 
   \param aoarrayt_minFeactures an output vector with the minimum attributes
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE
	   >
void 
minFeatures
(T_FEATURE             *aoarrayt_minFeactures,
 INPUT_ITERATOR        aiiterator_instfirst,
 const INPUT_ITERATOR  aiiterator_instlast
 )
{ 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::minFeatures";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 

  data::Instance<T_FEATURE>* linst_inter =
    (data::Instance<T_FEATURE>*)(*aiiterator_instfirst);
  interfacesse::copy
    (aoarrayt_minFeactures,
     linst_inter->getFeatures(),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
  ++aiiterator_instfirst;
  
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    linst_inter =  (data::Instance<T_FEATURE>*)(*aiiterator_instfirst);
    const T_FEATURE *larrayT_feacturesi = linst_inter->getFeatures();
    for (uintidx li_j = 0; li_j < data::Instance<T_FEATURE>::getNumDimensions(); li_j++)
      {
	if ( larrayT_feacturesi[li_j] < aoarrayt_minFeactures[li_j] )
	  aoarrayt_minFeactures[li_j] = larrayT_feacturesi[li_j];
      }
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelminFeatures;
    lostrstream_labelminFeatures << "<MINFEATURES:" << lpc_labelFunc;

     inout::containerprint
      (aoarrayt_minFeactures,
       aoarrayt_minFeactures+data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelminFeatures.str().c_str(),
       ','
       );
    
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

/*! \fn  void maxFeatures(T_FEATURE *aoarrayt_maxFeactures, INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR  aiiterator_instlast )
   \brief Obtain the maximum attributes of the instances
   \details 
   \param aoarrayt_maxFeactures an output vector with the maximum attributes
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
*/
template < typename INPUT_ITERATOR,
           typename T_FEATURE
	   >
void 
maxFeatures
(T_FEATURE             *aoarrayt_maxFeactures,
 INPUT_ITERATOR        aiiterator_instfirst,
 const INPUT_ITERATOR  aiiterator_instlast
 )
{ 
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::maxFeatures";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << ")"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/ 

  data::Instance<T_FEATURE>* linst_inter =
    (data::Instance<T_FEATURE>*) *aiiterator_instfirst;
  interfacesse::copy
    (aoarrayt_maxFeactures,
     linst_inter->getFeatures(),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
  ++aiiterator_instfirst;
  
  for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
    linst_inter =  (data::Instance<T_FEATURE>*)(*aiiterator_instfirst);
    const T_FEATURE *larrayT_feacturesi = linst_inter->getFeatures();
    for (uintidx li_j = 0; li_j < data::Instance<T_FEATURE>::getNumDimensions() ; li_j++)
      {
      if ( larrayT_feacturesi[li_j] > aoarrayt_maxFeactures[li_j] )
	aoarrayt_maxFeactures[li_j] = larrayT_feacturesi[li_j];
    }
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelmaxFeatures;
    lostrstream_labelmaxFeatures << "<MAXFEATURES:" << lpc_labelFunc;

    inout::containerprint
      (aoarrayt_maxFeactures,
       aoarrayt_maxFeactures+data::Instance<T_FEATURE>::getNumDimensions(),
       std::cout,
       lostrstream_labelmaxFeatures.str().c_str(),
       ','
       );
   
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
}


/*! \fn mat::MatrixRow<T_FEATURE> matrixVarcovar(INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast)
   \brief  Get covariance matrix, only for T_FEATURE real
   \details
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
*/
template < typename INPUT_ITERATOR>
auto 
matrixVarcovar
(INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR         aiiterator_instlast
 ) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::matrixVarcovar";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast["
	      << *aiiterator_instlast << "]\n"
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  
  mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>       
    lomatrixrowt_varcovar
    (data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions(),
     data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions()
     );

  mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>&& lmatrixrowt_y =
    data::toMatrixRowTrans
    (aiiterator_instfirst,
     aiiterator_instlast,
     data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getHomogeneousCoord()
     );
  
  std::vector<decltype((*aiiterator_instfirst)->getAttribute(0))>&&
    lvector_meanFeactures =
    meanRow<decltype((*aiiterator_instfirst)->getAttribute(0)),decltype((*aiiterator_instfirst)->getAttribute(0))>
    (lmatrixrowt_y);

  for (uintidx luintidx_i = 0; luintidx_i < lmatrixrowt_y.getNumRows(); ++luintidx_i) {

    interfacesse::transy
      (lmatrixrowt_y.getRow(luintidx_i),
       -lvector_meanFeactures[luintidx_i],
       lmatrixrowt_y.getNumColumns()
       );
  }
  
  /* Y' * Y
   */
  const decltype((*aiiterator_instfirst)->getAttribute(0))
    lirt_alpha =
    decltype((*aiiterator_instfirst)->getAttribute(0))(1.0) /
    (decltype((*aiiterator_instfirst)->getAttribute(0))) ( lmatrixrowt_y.getNumColumns()-1);
  interfacesse::gemm
    (lomatrixrowt_varcovar,
     lmatrixrowt_y,
     lmatrixrowt_y,
     mat::CblasNoTrans,
     mat::CblasTrans,
     lirt_alpha,
     0.0
     );
 
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelCovarMatrix;
    lostrstream_labelCovarMatrix
      << "<COVARMATRIX:" << lpc_labelFunc
      << ":lomatrixrowt_varcovar[" << &lomatrixrowt_varcovar << ']';
    lomatrixrowt_varcovar.print(std::cout,lostrstream_labelCovarMatrix.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lomatrixrowt_varcovar;

}


/*! \fn mat::MatrixRow<T_FEATURE> getMatrixDiagonal(T_FEATURE* aiarrayt_varianceFeactures)
  \brief Get weight matrix diagonal for calculate dist induced  \cite Bezdek:ClusterAnalysis:FCM:1974 \cite Bezdek:etal:ClusterAnalysis:FCM:1984 \cite Bezdek:etal:GAclustering:GA:1994
  \details 
  \param aiarrayt_varianceFeactures a array with variance feactures
 */
template <typename T_FEATURE>
mat::MatrixRow<T_FEATURE> 
getMatrixDiagonal
(T_FEATURE* aiarrayt_varianceFeactures)
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::getMatrixDiagonal"; 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input aiarrayt_varianceFeactures[" << aiarrayt_varianceFeactures << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  mat::MatrixRow<T_FEATURE> lomatrix_distWeight
    (data::Instance<T_FEATURE>::getNumDimensions(),
     data::Instance<T_FEATURE>::getNumDimensions()
     );
  lomatrix_distWeight.initialize();
 
  for (uintidx lst_i = 0; lst_i < lomatrix_distWeight.getNumRows(); lst_i++) {

    if ( aiarrayt_varianceFeactures[lst_i]  == 0.0 ) {
      std::ostringstream lostrstream_errormsj;
      lostrstream_errormsj
	<< "\n\nstats::getMatrixDiagonal: "
	<< "the atributo "
	<< lst_i + 1
	<< " of the instances has zero variance,\nit is not possible to calculate 1/var\n\n";
      throw std::runtime_error(lostrstream_errormsj.str());
    }
    lomatrix_distWeight(lst_i, lst_i) = 1.0 / aiarrayt_varianceFeactures[lst_i];
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelWeightMatrix;
    lostrstream_labelWeightMatrix
      << "<WEIGHTMATRIX:" << lpc_labelFunc
      << ":lomatrix_distWeight[" << &lomatrix_distWeight << ']';
    lomatrix_distWeight.print(std::cout,lostrstream_labelWeightMatrix.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
  return lomatrix_distWeight;

}


/*! \fn mat::MatrixRow<T_FEATURE> getMatrixMahalonobis(INPUT_ITERATOR  aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
  \brief Get weight matrix Mahalonobis for calculate dist induced  \cite Bezdek:ClusterAnalysis:FCM:1974 \cite Bezdek:etal:ClusterAnalysis:FCM:1984 \cite Bezdek:etal:GAclustering:GA:1994
  \details 
  \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
  \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
 */
template < typename INPUT_ITERATOR>
auto 
getMatrixMahalonobis
(INPUT_ITERATOR       aiiterator_instfirst,
 const INPUT_ITERATOR aiiterator_instlast
 ) -> mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "stats::getMatrixMahalonobis"; 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "(input INPUT_ITERATOR: aiiterator_instfirst[" << &aiiterator_instfirst
              << " input INPUT_ITERATOR: aiiterator_instlast[" << &aiiterator_instlast
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  mat::MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>&& lomatrix_distWeight =
    stats::matrixVarcovar 
    (aiiterator_instfirst,
     aiiterator_instlast
     );

  if ( mat::inverse(lomatrix_distWeight)  > 0) {
    throw std::runtime_error
      ("\n\nstats::getMatrixMahalonobis: "
       "interfacelapack_tgetrf failed to compute inverse matrix,\nsome attribute of the instances is dependent\n\n"
       );
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelWeightMatrix;
    lostrstream_labelWeightMatrix
      << "<WEIGHTMATRIX:" << lpc_labelFunc
      << ":lomatrix_distWeight[" << &lomatrix_distWeight << ']';
    lomatrix_distWeight.print(std::cout,lostrstream_labelWeightMatrix.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lomatrix_distWeight;

}

/*! void transform (INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, const mat::MatrixRow<T_FEATURE> &aimatrixrow_transform)
   \brief Make a transformation to the instances
   \details 
   \f[
    x_{i}^{'} = A \cdot x_{i}
   \f]
   \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
   \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
   \param aimatrixrow_transform a transformation matrix
*/
template < typename INPUT_ITERATOR,
	     typename T_FEATURE
	     > 
void
transform
(INPUT_ITERATOR                   aiiterator_instfirst,
const INPUT_ITERATOR              aiiterator_instlast,
const mat::MatrixRow<T_FEATURE>   &aimatrixrow_transform
)
{
#ifdef __VERBOSE_YES
  INPUT_ITERATOR  literator_instfirst =  aiiterator_instfirst;
   ++geiinparam_verbose;
   if ( geiinparam_verbose <= geiinparam_verboseMax ) {
     std::cout << "stats::transform  IN"
	       << '(' << geiinparam_verbose << ')'
	       << "(input INPUT_ITERATOR aiiterator_instfirst["
	       << *aiiterator_instfirst << "]\n"
	       << " input INPUT_ITERATOR aiiterator_instlast["
	       << *aiiterator_instlast << "]\n"
	       << "\n\t input  MatrixColumn<T_FEATURE>: &aimatrixrow_transform[" 
	       << &aimatrixrow_transform
	       << "\n\t)"
	       << std::endl;
   }
#endif // __VERBOSE_YES

   if (data::Instance<T_FEATURE>::getHomogeneousCoord() == false)
     throw  std::invalid_argument("stats::transform: instances to be homogeneous coord");
   
   T_FEATURE *larrayt_instFeactureTrans
     = new T_FEATURE[data::Instance<T_FEATURE>::getNumDimensions()];
 
   for (; aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {

     T_FEATURE* linstfeacture_inter =
       ((data::Instance<T_FEATURE>*) *aiiterator_instfirst)->getFeatures();
     
	interfacesse::gemv
	  (larrayt_instFeactureTrans,
	   aimatrixrow_transform,
	   linstfeacture_inter
	   );
	
	interfacesse::copy
	  (linstfeacture_inter,
	   larrayt_instFeactureTrans,
	   data::Instance<T_FEATURE>::getNumDimensions()
	   );     
   }
     

#ifdef __VERBOSE_YES 
   if ( geiinparam_verbose <= geiinparam_verboseMax ) {
     std::cout << "stats::transform OUT"
	       << '(' << geiinparam_verbose << ')';  
     ++geiinparam_verbose;
     if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       for (; literator_instfirst != aiiterator_instlast; ++literator_instfirst) {
	 data::Instance<T_FEATURE>* linst_inter =
	   (data::Instance<T_FEATURE>*) *literator_instfirst;
	 std::cout << *linst_inter << '\n';
       }
     }
     --geiinparam_verbose;
     std::cout << std::endl;
   }
   --geiinparam_verbose;
#endif //__VERBOSE_YES

   delete[] larrayt_instFeactureTrans;
}

template < typename T_IDX,
	   typename T
	 >
void rank(std::vector<std::pair<T_IDX,T> > &aiovectorpair_rank)
{
  const uintidx lui_n = (uintidx) aiovectorpair_rank.size();
  uintidx lui_j = 1, lui_ji, lui_jt;
  T       lt_rank;

  std::sort
    (aiovectorpair_rank.begin(),
     aiovectorpair_rank.end(),
     [](const std::pair<T_IDX,T> &lpair_a,
	const std::pair<T_IDX,T> &lpair_b
	)
     {
       return lpair_a.second < lpair_b.second;
     }
     );
  
  while ( lui_j < lui_n ) {
    if ( aiovectorpair_rank[lui_j].second != aiovectorpair_rank[lui_j-1].second) {
      aiovectorpair_rank[lui_j-1].second = lui_j;
      ++lui_j;
    } else {
      for (lui_jt=lui_j+1;
	   lui_jt <= lui_n && aiovectorpair_rank[lui_jt-1].second
	     == aiovectorpair_rank[lui_j-1].second;
	   lui_jt++);
      lt_rank = 0.5*(lui_j+lui_jt-1);
      for (lui_ji=lui_j;lui_ji<=(lui_jt-1);lui_ji++)
	aiovectorpair_rank[lui_ji-1].second = lt_rank;
      //T lt_t = lui_jt-lui_j;
      lui_j=lui_jt;
    }
  } 
  if( lui_j== lui_n) aiovectorpair_rank[lui_n-1].second=lui_n;

  std::sort
    (aiovectorpair_rank.begin(),
     aiovectorpair_rank.end(),
     [](const std::pair<T_IDX,T> &lpair_a,
	const std::pair<T_IDX,T> &lpair_b
	)
     {
       return lpair_a.first < lpair_b.first;
     }
     );  
}
  
} /*END namespace stats
   */

#endif /*__STATS_INSTANCES_HPP*/
