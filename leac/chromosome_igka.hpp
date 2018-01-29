/*! \file chromosome_igka.hpp
 * 
 * \brief chromosome IGKA \cite Lu:etal:GAclusteringLabel:IGKA:2004 and \cite Lu:etal:GAclusteringLabel:FGKA:2004  
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CHROMOSOME_IGKA_HPP
#define CHROMOSOME_IGKA_HPP

#include <utility>      // std::move
#include <vector>
#include "matrix.hpp"
#include "chromosome_fixedlength.hpp"
#include "instance.hpp"
#include "partition_label.hpp"
#include "linear_algebra_level1.hpp"
#include "clustering_operator_centroids.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {


/*! \class ChromosomeIGKA
  \brief Chromosome IGKA
*/
template <class T_CLUSTERIDX,
	  class T_METRIC,
	  class T_FEATURE,
	  class T_FEATURE_SUM,
	  class T_INSTANCES_CLUSTER_K
	  >
class ChromosomeIGKA: 
  public gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>
{
public:
  ChromosomeIGKA
  (const uintidx   aiuintidx_numClusterK
   ) 
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>()
    , _matrixrowt_centroids
    (aiuintidx_numClusterK,data::Instance<T_FEATURE>::getNumDimensions())
    , _matrixrowt_sumInstancesCluster
    (aiuintidx_numClusterK,data::Instance<T_FEATURE>::getNumDimensions())
    , _vectorT_numInstancesClusterK(aiuintidx_numClusterK)
    , _arrayT_WCV(new T_METRIC[aiuintidx_numClusterK])
    , _matrixrowt_distINSTiCLUSTE1k
    (gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize(),
     aiuintidx_numClusterK
     )
    , _arrayT_distSumINSTiCLUSTE1k
    (new T_METRIC[gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize()])
    , _arrayT_distMaxINSTiCLUSTEk
    (new T_METRIC[gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize()])
    , _arrayT_distMinINSTiCLUSTEk
    (new T_METRIC[gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize()])
    , _arraycidx_idxMinINSTiCLUSTEk
    (new T_CLUSTERIDX[gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize()])
    , _arraycidx_idxMaxINSTiCLUSTEk
    (new T_CLUSTERIDX[gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize()])
    , _matrixrowt_sumInstancesClusterDelta
    (aiuintidx_numClusterK,
     data::Instance<T_FEATURE>::getNumDimensions()
     )
    , _vectorT_numInstancesClusterKDelta(aiuintidx_numClusterK,0)
    , _vectorb_changeCluster(aiuintidx_numClusterK,false)

    , _cidx_numClusterNotNull(0)
    , _t_TWCVDelta(0)
    , _b_selected(false)
  {     
  }

  //copy constructor
  ChromosomeIGKA
  (const ChromosomeIGKA
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
    T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K>                   &aich_chromosome
   )
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>(aich_chromosome)
    , _matrixrowt_centroids(aich_chromosome._matrixrowt_centroids)
    , _matrixrowt_sumInstancesCluster(aich_chromosome._matrixrowt_sumInstancesCluster)
    , _vectorT_numInstancesClusterK(aich_chromosome._vectorT_numInstancesClusterK)

    , _arrayT_WCV(new T_METRIC[aich_chromosome._matrixrowt_centroids.getNumRows()])
    , _matrixrowt_distINSTiCLUSTE1k(aich_chromosome._matrixrowt_distINSTiCLUSTE1k)

    , _arrayT_distSumINSTiCLUSTE1k(new T_METRIC[aich_chromosome.getStringSize()])
    , _arrayT_distMaxINSTiCLUSTEk(new T_METRIC[aich_chromosome.getStringSize()])
    , _arrayT_distMinINSTiCLUSTEk(new T_METRIC[aich_chromosome.getStringSize()])
    , _arraycidx_idxMinINSTiCLUSTEk(new T_CLUSTERIDX[aich_chromosome.getStringSize()])
    , _arraycidx_idxMaxINSTiCLUSTEk(new T_CLUSTERIDX[aich_chromosome.getStringSize()])

    , _matrixrowt_sumInstancesClusterDelta
    (aich_chromosome._matrixrowt_sumInstancesClusterDelta)
    , _vectorT_numInstancesClusterKDelta
      (aich_chromosome._vectorT_numInstancesClusterKDelta)
    , _vectorb_changeCluster(aich_chromosome._vectorb_changeCluster)

    , _cidx_numClusterNotNull(aich_chromosome._cidx_numClusterNotNull)
    , _t_TWCVDelta(aich_chromosome._t_TWCVDelta)
    , _b_selected(aich_chromosome._b_selected)
  {
    
    interfacesse::copy
      (_arrayT_WCV,
       aich_chromosome._arrayT_WCV,
       _matrixrowt_centroids.getNumRows()
       );

    interfacesse::copy
      (_arrayT_distSumINSTiCLUSTE1k,
       aich_chromosome._arrayT_distSumINSTiCLUSTE1k,
       aich_chromosome.getStringSize()
       );

    interfacesse::copy
      (_arrayT_distMaxINSTiCLUSTEk,
       aich_chromosome._arrayT_distMaxINSTiCLUSTEk,
       aich_chromosome.getStringSize()
       );

    interfacesse::copy
      (_arrayT_distMinINSTiCLUSTEk,
       aich_chromosome._arrayT_distMinINSTiCLUSTEk,
       aich_chromosome.getStringSize()
       );

    interfacesse::copy
      (_arraycidx_idxMinINSTiCLUSTEk,
       aich_chromosome._arraycidx_idxMinINSTiCLUSTEk,
       aich_chromosome.getStringSize()
       );

    interfacesse::copy
      (_arraycidx_idxMaxINSTiCLUSTEk,
       aich_chromosome._arraycidx_idxMaxINSTiCLUSTEk,
       aich_chromosome.getStringSize()
       );
  }
  
  //move constructor
  ChromosomeIGKA
  (ChromosomeIGKA
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
    T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K>                 &&aich_chromosome
   )
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>(aich_chromosome)
    , _matrixrowt_centroids
      (std::move(aich_chromosome._matrixrowt_centroids))
    , _matrixrowt_sumInstancesCluster
      (std::move(aich_chromosome._matrixrowt_sumInstancesCluster))
    , _vectorT_numInstancesClusterK
      (std::move((aich_chromosome._vectorT_numInstancesClusterK)))

    , _arrayT_WCV(aich_chromosome._arrayT_WCV)
    , _matrixrowt_distINSTiCLUSTE1k
      (std::move(aich_chromosome._matrixrowt_distINSTiCLUSTE1k))
    , _arrayT_distSumINSTiCLUSTE1k(aich_chromosome._arrayT_distSumINSTiCLUSTE1k) 
    , _arrayT_distMaxINSTiCLUSTEk(aich_chromosome._arrayT_distMaxINSTiCLUSTEk) 
    , _arrayT_distMinINSTiCLUSTEk(aich_chromosome._arrayT_distMinINSTiCLUSTEk)
    , _arraycidx_idxMinINSTiCLUSTEk(aich_chromosome._arraycidx_idxMinINSTiCLUSTEk)
    , _arraycidx_idxMaxINSTiCLUSTEk(aich_chromosome._arraycidx_idxMaxINSTiCLUSTEk)

    , _matrixrowt_sumInstancesClusterDelta
      (std::move(aich_chromosome._matrixrowt_sumInstancesClusterDelta))
    , _vectorT_numInstancesClusterKDelta
      (std::move(aich_chromosome._vectorT_numInstancesClusterKDelta))
    , _vectorb_changeCluster(std::move(aich_chromosome._vectorb_changeCluster))

    , _cidx_numClusterNotNull(aich_chromosome._cidx_numClusterNotNull)
    , _t_TWCVDelta(aich_chromosome._t_TWCVDelta)
    , _b_selected(aich_chromosome._b_selected)
  {
    aich_chromosome._arrayT_WCV = NULL;
    aich_chromosome._arrayT_distSumINSTiCLUSTE1k = NULL;
    aich_chromosome._arrayT_distMaxINSTiCLUSTEk = NULL; 
    aich_chromosome._arrayT_distMinINSTiCLUSTEk = NULL;
    aich_chromosome._arraycidx_idxMinINSTiCLUSTEk = NULL;
    aich_chromosome._arraycidx_idxMaxINSTiCLUSTEk = NULL;
  }

  virtual ~ChromosomeIGKA() 
  {
   
    if (_arrayT_WCV != NULL)
      delete[] _arrayT_WCV;
    if (_arrayT_distSumINSTiCLUSTE1k != NULL)
      delete[] _arrayT_distSumINSTiCLUSTE1k;
    if (_arrayT_distMaxINSTiCLUSTEk != NULL)
      delete[] _arrayT_distMaxINSTiCLUSTEk;
    if (_arrayT_distMinINSTiCLUSTEk != NULL)
      delete _arrayT_distMinINSTiCLUSTEk;
    if (_arraycidx_idxMinINSTiCLUSTEk != NULL)
      delete[] _arraycidx_idxMinINSTiCLUSTEk;
    if (_arraycidx_idxMaxINSTiCLUSTEk != NULL)
      delete[] _arraycidx_idxMaxINSTiCLUSTEk;
   
  }
  
  void initialize
  (const std::vector<data::Instance<T_FEATURE>* >  &aivectorptinst_instances,
   dist::Dist<T_METRIC,T_FEATURE>                  &aifunc2p_dista
   )
  {

    T_CLUSTERIDX lcidx_numClusterK  = 
      (T_CLUSTERIDX) _matrixrowt_centroids.getNumRows();
    
    partition::PartitionLabel
      <T_CLUSTERIDX>
      lpartition_clusters
      (this->getString(),
       this->getStringSize(),
       lcidx_numClusterK
       );
    
    this->_matrixrowt_sumInstancesCluster.initialize();
    
    interfacesse::copya
      (this->_vectorT_numInstancesClusterK.data(),
       T_INSTANCES_CLUSTER_K(0),
       (uintidx) this->_vectorT_numInstancesClusterK.size()
       );
    
    T_CLUSTERIDX lcidx_numClusterNull =
      clusteringop::getCentroids
      (this->_matrixrowt_centroids,
       this->_matrixrowt_sumInstancesCluster,
       this->_vectorT_numInstancesClusterK,
       lpartition_clusters,
       aivectorptinst_instances.begin(),
       aivectorptinst_instances.end()
       );
    
    this->_cidx_numClusterNotNull = 
      this->_matrixrowt_centroids.getNumRows() - lcidx_numClusterNull;
    this->setValidString
      ((uintidx) this->_cidx_numClusterNotNull != _matrixrowt_centroids.getNumRows() );

    /*CALCULATE DISTANCES
     */
    interfacesse::copya
      (_arrayT_WCV,
       T_METRIC(0),
       _matrixrowt_centroids.getNumRows()
       );
    
    for ( uintidx luintidx_i = 0; 
	  luintidx_i < aivectorptinst_instances.size(); 
	  luintidx_i++) 
      { /*BEGIN FOR DIST*/
	data::Instance<T_FEATURE>* liter_iInstance =
	  aivectorptinst_instances.at(luintidx_i);
	this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] = 0.0; 
	this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] = std::numeric_limits<T_METRIC>::max();
	this->_arrayT_distSumINSTiCLUSTE1k[luintidx_i] = 0.0;
	for (  T_CLUSTERIDX li_k = 0; li_k < lcidx_numClusterK; li_k++) {
	  T_METRIC lt_distINSTiCLUSTERk =
	    aifunc2p_dista
	    (this->_matrixrowt_centroids.getRow(li_k),
	     liter_iInstance->getFeatures(),
	     liter_iInstance->getNumDimensions() 
	     );
	  this->_arrayT_distSumINSTiCLUSTE1k[luintidx_i] += lt_distINSTiCLUSTERk;

	  if ( lt_distINSTiCLUSTERk > this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] ) {
	    this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] = lt_distINSTiCLUSTERk;
	    this->_arraycidx_idxMaxINSTiCLUSTEk[luintidx_i] = li_k;
	  }
	  if ( lt_distINSTiCLUSTERk < this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] ) {
	    this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] = lt_distINSTiCLUSTERk;
	    this->_arraycidx_idxMinINSTiCLUSTEk[luintidx_i] = li_k;
	  }

	  _matrixrowt_distINSTiCLUSTE1k(luintidx_i,li_k) = lt_distINSTiCLUSTERk;
	}
	_arrayT_WCV[this->getGene(luintidx_i)] +=
	  (T_METRIC) _matrixrowt_distINSTiCLUSTE1k(luintidx_i,this->getGene(luintidx_i));
      } /*END FOR DIST*/
    T_METRIC lT_TWCV = 
      interfacesse::sum
      (_arrayT_WCV,
       _matrixrowt_centroids.getNumRows()
       );
    this->setObjetiveFunc(lT_TWCV);
    
    _matrixrowt_sumInstancesClusterDelta.initialize();
    
  }
  
   inline 
  mat::MatrixRow<T_FEATURE>&  getCentroids()
  {
    return this->_matrixrowt_centroids;
  }

  inline 
  T_METRIC  getDistINSTiCLUSTE1k
  (uintidx aiuintidx_idxInstance, T_CLUSTERIDX aicidx_clusterK)  
  {
    return this->_matrixrowt_distINSTiCLUSTE1k(aiuintidx_idxInstance,aicidx_clusterK);  
  }

  inline 
  T_METRIC  getDistINSTiCLUSTE1k
  (uintidx aiuintidx_idxInstance, uintidx aiuintidx_clusterK)  
  {
    return this->_matrixrowt_distINSTiCLUSTE1k(aiuintidx_idxInstance,aiuintidx_clusterK);  
  }

  inline
  T_METRIC getDistSumINSTiCLUSTE1k(uintidx aiuintidx_idxInstance)
  {
    return this->_arrayT_distSumINSTiCLUSTE1k[aiuintidx_idxInstance]; 
  }

  inline
  T_METRIC getDistMaxINSTiCLUSTEk(uintidx aiuintidx_idxInstance)
  {
    return this->_arrayT_distMaxINSTiCLUSTEk[aiuintidx_idxInstance]; 
  }

  inline 
  T_CLUSTERIDX  getIdxMinINSTiCLUSTEk(uintidx aiuintidx_idxInstance)  
  { 
    return this->_arraycidx_idxMinINSTiCLUSTEk[aiuintidx_idxInstance];
  }
  
  inline T_CLUSTERIDX  getNumClusterNotNull()
  {
    return this->_cidx_numClusterNotNull;
  }

  inline T_METRIC getLegalityRation() const
  {
    return  (T_METRIC) this->_cidx_numClusterNotNull 
      / (T_METRIC) this->_matrixrowt_centroids.getNumRows();
  }

  inline bool getSelected() 
  {
    return this->_b_selected;
  }

  inline void setSelected(bool aib_selected) 
  {
    this->_b_selected = aib_selected;
  }
  
  void accumulateUpdate
  (T_CLUSTERIDX           aicidx_newAllen,
   T_CLUSTERIDX           aicidx_oldAllen,
   data::Instance<T_FEATURE>*  aiptinst_instance   
   )
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeIGKA::incrementalUpdate";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeIGKA<>: this[" <<  this
		<< "]\n";
	std::ostringstream lostrstream_labelChromDelt;
	lostrstream_labelChromDelt
	  << "<SUMINSTANCESCLUSTERDELTA:"
	  << lpc_labelFunc;
	_matrixrowt_sumInstancesClusterDelta.print
	  (std::cout,lostrstream_labelChromDelt.str().c_str(),',',';');
	std::cout <<  "\n)"
		  << std::endl;
    }
#endif //__VERBOSE_YES
    
    _vectorb_changeCluster.at(aicidx_oldAllen) = true;
    _vectorb_changeCluster.at(aicidx_newAllen) = true;
    --_vectorT_numInstancesClusterKDelta[aicidx_oldAllen];
    ++_vectorT_numInstancesClusterKDelta[aicidx_newAllen];

    interfacesse::axpy
      (_matrixrowt_sumInstancesClusterDelta.getRow(aicidx_oldAllen),
       T_FEATURE(-1),
       aiptinst_instance->getFeatures(), 
       aiptinst_instance->getNumDimensions()
       );

    interfacesse::axpy
      (_matrixrowt_sumInstancesClusterDelta.getRow(aicidx_newAllen),
       T_FEATURE(1),
       aiptinst_instance->getFeatures(), 
       aiptinst_instance->getNumDimensions()
       );

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelChromDelt;
	lostrstream_labelChromDelt
	  << "<SUMINSTANCESCLUSTERDELTA:"
	  << lpc_labelFunc;
	_matrixrowt_sumInstancesClusterDelta.print
	  (std::cout,lostrstream_labelChromDelt.str().c_str(),',',';');
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }
  
  void incrementalUpdate
  (std::vector
   <data::Instance<T_FEATURE>* >    &aivectorptinst_instances,
   dist::Dist<T_METRIC,T_FEATURE>   &aifunc2p_dista
   ) 
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeIGKA::incrementalUpdate";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeIGKA<>: this[" <<  this
		<< "]\n";
      ChromosomeIGKA
	<T_CLUSTERIDX,
	 T_METRIC,
	 T_FEATURE,
	 T_FEATURE_SUM,
	 T_INSTANCES_CLUSTER_K>
	::print();
	std::cout << '\n';
	std::ostringstream lostrstream_labelChromDelt;
	lostrstream_labelChromDelt
	  << "<SUMINSTANCESCLUSTERDELTA:"
	  << lpc_labelFunc;
	_matrixrowt_sumInstancesClusterDelta.print
	  (std::cout,lostrstream_labelChromDelt.str().c_str(),',',';');
	std::cout <<  "\n)"
		  << std::endl;
    }
#endif //__VERBOSE_YES
    
   const  T_CLUSTERIDX lcidx_numClusterK  = 
      (T_CLUSTERIDX) _matrixrowt_centroids.getNumRows();

    _t_TWCVDelta = 0;

    for (T_CLUSTERIDX lcidx_k = 0; lcidx_k < lcidx_numClusterK; lcidx_k++) {

      if ( _vectorb_changeCluster[lcidx_k] != true )  continue;

      T_INSTANCES_CLUSTER_K lick_oldNumInstancesClusterK = 
	_vectorT_numInstancesClusterK[lcidx_k];

      _vectorT_numInstancesClusterK[lcidx_k] +=
	_vectorT_numInstancesClusterKDelta[lcidx_k];
      
      T_FEATURE_SUM* larrayrowt_sumInstCluster_k 
	 =  _matrixrowt_sumInstancesCluster.getRow(lcidx_k);
       T_FEATURE_SUM* larrayrowt_sumInstClusterDelta_k 
	 =  _matrixrowt_sumInstancesClusterDelta.getRow(lcidx_k);
       uintidx luintidx_i = 0;
       std::for_each
	 (larrayrowt_sumInstCluster_k, 
	  larrayrowt_sumInstCluster_k + _matrixrowt_sumInstancesCluster.getNumColumns(),
	  [&] (T_FEATURE_SUM &liter_sumInstCluster_ki) 
	  { 
	    liter_sumInstCluster_ki += larrayrowt_sumInstClusterDelta_k[luintidx_i]; 
	    larrayrowt_sumInstClusterDelta_k[luintidx_i++] = (T_FEATURE_SUM) 0;
	  } 
	  );
       
       /*MISING axpy long a long*/
       
      if ( _vectorT_numInstancesClusterK[lcidx_k] > 0 ) {
	stats::meanVector
	  (_matrixrowt_centroids.getRow(lcidx_k),
	   _vectorT_numInstancesClusterK[lcidx_k],
	   _matrixrowt_sumInstancesCluster.getRow(lcidx_k),
	   lcidx_k
	   );
      }
      
      _vectorb_changeCluster[lcidx_k] = false;
      _vectorT_numInstancesClusterKDelta[lcidx_k] = 0;

      if ( lick_oldNumInstancesClusterK == 0 &&  _vectorT_numInstancesClusterK[lcidx_k] > 0 )
	_cidx_numClusterNotNull++;
      else if (lick_oldNumInstancesClusterK > 0 && _vectorT_numInstancesClusterK[lcidx_k] == 0 )
	_cidx_numClusterNotNull--;

      this->setValidString
	(this->_cidx_numClusterNotNull != lcidx_numClusterK );
      
      T_METRIC lT_oldWCVK = _arrayT_WCV[lcidx_k];
      _arrayT_WCV[lcidx_k] = 0;
      T_FEATURE* larrayrowt_centroid =  _matrixrowt_centroids.getRow(lcidx_k);
      
      for ( uintidx luintidx_i = 0; 
	    luintidx_i < aivectorptinst_instances.size(); 
	    luintidx_i++) 
	{ /*BEGIN FOR DIST*/
	  if ( this->getGene(luintidx_i) == lcidx_k ) {
	    data::Instance<T_FEATURE>* liter_iInstance = 
	      aivectorptinst_instances.at(luintidx_i);     
	    _arrayT_WCV[lcidx_k] += 
	      aifunc2p_dista
	      (larrayrowt_centroid,
	       liter_iInstance->getFeatures(),
	       liter_iInstance->getNumDimensions()
	       );
	  }
	}
      this->_t_TWCVDelta +=  _arrayT_WCV[lcidx_k] - lT_oldWCVK;

      for ( uintidx luintidx_i = 0; 
	    luintidx_i < aivectorptinst_instances.size(); 
	    luintidx_i++) 
	{ /*BEGIN FOR DIST*/
	  T_METRIC lt_oldDistINSTiCLUSTERk = 
	    _matrixrowt_distINSTiCLUSTE1k(luintidx_i,lcidx_k);
	  data::Instance<T_FEATURE>* liter_iInstance = 
	    aivectorptinst_instances.at(luintidx_i);
	  T_METRIC lt_distINSTiCLUSTERk =
	    aifunc2p_dista
	    (larrayrowt_centroid, 
	     liter_iInstance->getFeatures(),
	     liter_iInstance->getNumDimensions() 
	     );
	  _matrixrowt_distINSTiCLUSTE1k(luintidx_i,lcidx_k) = lt_distINSTiCLUSTERk;

	  if ( lt_distINSTiCLUSTERk > this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] ) {
	    this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] = lt_distINSTiCLUSTERk;
	    this->_arraycidx_idxMaxINSTiCLUSTEk[luintidx_i] = lcidx_k;
	  }
	  if ( lt_distINSTiCLUSTERk < this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] ) {
	    this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] = lt_distINSTiCLUSTERk;
	    this->_arraycidx_idxMinINSTiCLUSTEk[luintidx_i] = lcidx_k;
	  }
	    
	  if (lt_distINSTiCLUSTERk < lt_oldDistINSTiCLUSTERk  &&  
	      this->_arraycidx_idxMaxINSTiCLUSTEk[luintidx_i] == lcidx_k) {

	    this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] = 0.0;
 
	    for (uintidx luintidx_l = 0; luintidx_l < this->_matrixrowt_centroids.getNumRows(); luintidx_l++) {
	      if (this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] < _matrixrowt_distINSTiCLUSTE1k(luintidx_i,luintidx_l)) {
		this->_arrayT_distMaxINSTiCLUSTEk[luintidx_i] = _matrixrowt_distINSTiCLUSTE1k(luintidx_i,luintidx_l);
		this->_arraycidx_idxMaxINSTiCLUSTEk[luintidx_i] = luintidx_l;
	      }
	    }
	  }

	  if (lt_distINSTiCLUSTERk > lt_oldDistINSTiCLUSTERk  &&  
	      this->_arraycidx_idxMinINSTiCLUSTEk[luintidx_i] == lcidx_k) {
	    this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] = std::numeric_limits<T_METRIC>::max();
	    for (uintidx luintidx_l = 0; luintidx_l < this->_matrixrowt_centroids.getNumRows(); luintidx_l++) {
	      if (this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] > _matrixrowt_distINSTiCLUSTE1k(luintidx_i,luintidx_l)) {
		this->_arrayT_distMinINSTiCLUSTEk[luintidx_i] = _matrixrowt_distINSTiCLUSTE1k(luintidx_i,luintidx_l);
		this->_arraycidx_idxMinINSTiCLUSTEk[luintidx_i] = luintidx_l;
	      }
	    }
	  }
	  this->_arrayT_distSumINSTiCLUSTE1k[lcidx_k] += 
	    lt_distINSTiCLUSTERk - lt_oldDistINSTiCLUSTERk;
	}
    } /*END FOR K*/ 
    T_METRIC lT_TWCV =  this->getObjetiveFunc() + this->_t_TWCVDelta; 
    this->setObjetiveFunc(lT_TWCV);

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
       ChromosomeIGKA
	<T_CLUSTERIDX,
	 T_METRIC,
	 T_FEATURE,
	 T_FEATURE_SUM,
	 T_INSTANCES_CLUSTER_K>
	::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  } /*END accumulateUpdate*/

  
  ChromosomeIGKA
  <T_CLUSTERIDX,
   T_METRIC,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>& 
  operator=
  (const ChromosomeIGKA
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
    T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K> &aichrom_b)
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeIGKA::operator:copy=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeIGKA<>: &&aichrom_b[" <<  &aichrom_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichrom_b ) {
      gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(aichrom_b);

      this->_matrixrowt_centroids = aichrom_b._matrixrowt_centroids;
      this->_matrixrowt_sumInstancesCluster = aichrom_b._matrixrowt_sumInstancesCluster;
      this->_vectorT_numInstancesClusterK = aichrom_b._vectorT_numInstancesClusterK;
      this->_matrixrowt_distINSTiCLUSTE1k = aichrom_b._matrixrowt_distINSTiCLUSTE1k;

      interfacesse::copy
      (_arrayT_WCV,
       aichrom_b._arrayT_WCV,
       this->_matrixrowt_centroids.getNumRows() 
       );
      
    interfacesse::copy
      (_arrayT_distSumINSTiCLUSTE1k,
       aichrom_b._arrayT_distSumINSTiCLUSTE1k,
       aichrom_b.getStringSize()
       );

    interfacesse::copy
      (_arrayT_distMaxINSTiCLUSTEk,
       aichrom_b._arrayT_distMaxINSTiCLUSTEk,
       aichrom_b.getStringSize()
       );

    interfacesse::copy
      (_arrayT_distMinINSTiCLUSTEk,
       aichrom_b._arrayT_distMinINSTiCLUSTEk,
       aichrom_b.getStringSize()
       );

    interfacesse::copy
      (_arraycidx_idxMinINSTiCLUSTEk,
       aichrom_b._arraycidx_idxMinINSTiCLUSTEk,
       aichrom_b.getStringSize()
       );

    interfacesse::copy
      (_arraycidx_idxMaxINSTiCLUSTEk,
       aichrom_b._arraycidx_idxMaxINSTiCLUSTEk,
       aichrom_b.getStringSize()
       );

    this->_matrixrowt_sumInstancesClusterDelta = 
      aichrom_b._matrixrowt_sumInstancesClusterDelta;
    _vectorT_numInstancesClusterKDelta = aichrom_b._vectorT_numInstancesClusterKDelta;
    _vectorb_changeCluster = aichrom_b._vectorb_changeCluster;

    this->_cidx_numClusterNotNull = aichrom_b._cidx_numClusterNotNull;
    this->_t_TWCVDelta = aichrom_b._t_TWCVDelta;
    this->_b_selected = aichrom_b._b_selected;
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeIGKA
	<T_CLUSTERIDX,
	 T_METRIC,
	 T_FEATURE,
	 T_FEATURE_SUM,
	 T_INSTANCES_CLUSTER_K>
	::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    
    return *this;
  }

  ChromosomeIGKA
  <T_CLUSTERIDX,
   T_METRIC,
   T_FEATURE,
   T_FEATURE_SUM,
   T_INSTANCES_CLUSTER_K>&
  operator=
  (ChromosomeIGKA
   <T_CLUSTERIDX,
    T_METRIC,
    T_FEATURE,
    T_FEATURE_SUM,
    T_INSTANCES_CLUSTER_K>
   &&aichrom_b)
  {
    if ( this != &aichrom_b ) {
      gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(aichrom_b);
      this->_matrixrowt_centroids = std::move(aichrom_b._matrixrowt_centroids);
      this->_matrixrowt_sumInstancesCluster = std::move(aichrom_b._matrixrowt_sumInstancesCluster);
      this->_arrayT_numInstancesClusterK = std::move(aichrom_b._arrayT_numInstancesClusterK);
      
      this->_arrayT_WCV = aichrom_b._arrayT_WCV;
      this->_matrixrowt_distINSTiCLUSTE1k = std::move(aichrom_b._matrixrowt_distINSTiCLUSTE1k);
      this->_arrayT_distSumINSTiCLUSTE1k  =    aichrom_b._arrayT_distSumINSTiCLUSTE1k;
      this->_arrayT_distMaxINSTiCLUSTEk   =   aichrom_b._arrayT_distMaxINSTiCLUSTEk; 
      this->_arrayT_distMinINSTiCLUSTEk   =    aichrom_b._arrayT_distMinINSTiCLUSTEk;
      this->_arraycidx_idxMinINSTiCLUSTEk =     aichrom_b._arraycidx_idxMinINSTiCLUSTEk;
      this->_arraycidx_idxMaxINSTiCLUSTEk =     aichrom_b._arraycidx_idxMaxINSTiCLUSTEk;

      this->_matrixrowt_sumInstancesClusterDelta = std::move(aichrom_b._matrixrowt_sumInstancesClusterDelta);
      this->_vectorT_numInstancesClusterKDelta = std::move(aichrom_b._vectorT_numInstancesClusterKDelta);
      this->_vectorb_changeCluster = std::move(aichrom_b._vectorb_changeCluster);

      this->_cidx_numClusterNotNull = aichrom_b._cidx_numClusterNotNull;
      this->_t_TWCVDelta = aichrom_b._t_TWCVDelta;
      
      this->_b_selected = aichrom_b._b_selected;
    
      aichrom_b._arrayT_WCV = NULL;
      aichrom_b._lmatrixT_distINSTiCLUSTE1k = NULL;
      aichrom_b._arrayT_distSumINSTiCLUSTE1k = NULL;
      aichrom_b._arrayT_distMaxINSTiCLUSTEk = NULL; 
      aichrom_b._arrayT_distMinINSTiCLUSTEk = NULL;
      aichrom_b._arraycidx_idxMinINSTiCLUSTEk = NULL;
      aichrom_b._arraycidx_idxMaxINSTiCLUSTEk = NULL;
      
 
    }

    return *this;
  }

  virtual void  print
  (std::ostream &os=std::cout,
   const char* aipc_label = "",
   const char aic_delimCoef=',',
   const char aic_delimRow=';'
   ) const
  {

    std::ostringstream lostrstream_labelChrom;
    lostrstream_labelChrom
      << aipc_label
      <<  ":e" << aic_delimCoef << this->getLegalityRation();
      
    gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::
      print(os,lostrstream_labelChrom.str().c_str(),aic_delimCoef,aic_delimRow);
    os << '\n';
    std::ostringstream lostrstream_labelChromCent;
    lostrstream_labelChromCent
      << "<CENTROIDS:CHROMOSOMEFGKA:"
      << aipc_label;
    _matrixrowt_centroids.print(os,lostrstream_labelChromCent.str().c_str(),aic_delimCoef,aic_delimRow);

  }
  
protected:

  mat::MatrixRow<T_FEATURE>     _matrixrowt_centroids;
  mat::MatrixRow<T_FEATURE_SUM> _matrixrowt_sumInstancesCluster;
  std::vector
  <T_INSTANCES_CLUSTER_K> _vectorT_numInstancesClusterK;
  
  T_METRIC                *_arrayT_WCV;
  mat::MatrixRow<T_METRIC>          _matrixrowt_distINSTiCLUSTE1k;
  T_METRIC                  *_arrayT_distSumINSTiCLUSTE1k; 
  T_METRIC                  *_arrayT_distMaxINSTiCLUSTEk; 
  T_METRIC                  *_arrayT_distMinINSTiCLUSTEk;
  T_CLUSTERIDX     *_arraycidx_idxMinINSTiCLUSTEk;
  T_CLUSTERIDX     *_arraycidx_idxMaxINSTiCLUSTEk;

  mat::MatrixRow<T_FEATURE_SUM> _matrixrowt_sumInstancesClusterDelta;
  std::vector
  <T_INSTANCES_CLUSTER_K> _vectorT_numInstancesClusterKDelta;
  std::vector<bool>       _vectorb_changeCluster;

  T_CLUSTERIDX     _cidx_numClusterNotNull;
  T_METRIC                _t_TWCVDelta;
  bool                    _b_selected; /*used*/
}; /*END  ChromosomeIGKA */

} /*END namespace gaencode*/

#endif  /* CHROMOSOME_IGKA_HPP */
