/*! \file clustering_operator_crispmatrix.hpp
 *
 * \brief clustering operator crispmatrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CRISPMATRIX_CLUSTERING_HPP
#define CRISPMATRIX_CLUSTERING_HPP

#include <iostream>
#include <vector>
#include "crisp_matrix.hpp"
#include "instance_class.hpp"
#include "dist.hpp"
#include "nearestinstance_operator.hpp"
#include "clustering_operator_centroids.hpp"
#include "bitmatrix_matrix_operator.hpp"

#include "verbose_global.hpp"

extern std::mt19937       gmt19937_eng;

/*! \namespace clusteringop
  \brief Clustering operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace clusteringop {

/*! \fn void initialize(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>& aiocrispmatrix_partition)
    \brief Initializes a crisp matrix randomly
    \details
    \param aiocrispmatrix_partition a mat::CrispMatrix 
*/
template < typename T_BITSIZE,
	   typename T_CLUSTERIDX
	   >
void
randomInitialize
(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>& aiocrispmatrix_partition)
{
  mat::BitArray<T_BITSIZE> larraybit_rowBitOn(aiocrispmatrix_partition.getNumRows());
  
  std::uniform_int_distribution<uintidx> luniformdis_ui0K
    (0,aiocrispmatrix_partition.getNumRows()-1);
    
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::initialize";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output mat::CrispMatrix<T_BITSIZE>: aiocrispmatrix_partition["
	      << &aiocrispmatrix_partition << "]\n"
	      << ')'
	      << std::endl;
  }
#endif //__VERBOSE_YES

  aiocrispmatrix_partition.initialize();
  larraybit_rowBitOn.initialize();
    
  for (uintidx li_j = 0; li_j < aiocrispmatrix_partition.getNumColumns(); li_j++) {
    
    uintidx  luintidx_idxRow = luniformdis_ui0K(gmt19937_eng);
   
    aiocrispmatrix_partition.setBit(luintidx_idxRow,li_j);
    larraybit_rowBitOn.setBit(luintidx_idxRow);
  }
 
  //ALL ROWS WHEN LESS HAVE 1
  while ( larraybit_rowBitOn.getNumBitOn() != aiocrispmatrix_partition.getNumRows() )
    {
      aiocrispmatrix_partition.initialize();
      larraybit_rowBitOn.initialize();
    
      for (uintidx li_j = 0; li_j < aiocrispmatrix_partition.getNumColumns(); li_j++) {
	uintidx  luintidx_idxRow = luniformdis_ui0K(gmt19937_eng); 
	
	aiocrispmatrix_partition.setBit(luintidx_idxRow,li_j);
	larraybit_rowBitOn.setBit(luintidx_idxRow);
      }
     
    }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelCrispMatrix;
    lostrstream_labelCrispMatrix << "<BITCRISPMATRIX:"
				    << geverbosepc_labelstep
				    << ':' << lpc_labelFunc
				    << "aiocrispmatrix_partition["
				    << geverboseui_idproc << &aiocrispmatrix_partition << ']';
    aiocrispmatrix_partition.print(std::cout,lostrstream_labelCrispMatrix.str().c_str(),',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

}

/*! \fn void initializeOptimal(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>& aiocrispmatrix_partition)
    \brief Initializes the crisp matrix with the instance class as its cluster \cite Bezdek:etal:GAclustering:GA:1994
    \details
    \param aiocrispmatrix_partition a mat::CrispMatrix
    \param aiiterator_instfirst an InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast an InputIterator to the final positions of the sequence of instances
    \param function_instClass a function that gets the index of the class to which it belongs
*/
template < typename T_BITSIZE,
	   typename T_CLUSTERIDX,
	   typename INPUT_ITERATOR,
	   typename FUNCTION_INSTCLASS
	   >
void initializeOptimal
(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aiocrispmatrix_partition,
 INPUT_ITERATOR                           aiiterator_instfirst,
 const INPUT_ITERATOR                     aiiterator_instlast,
 FUNCTION_INSTCLASS                       function_instClass 
 )
{
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "clusteringop::initializeOptimal:  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t(output CrispMatrix: aiocrispmatrix_partition["
	      << &aiocrispmatrix_partition << "]\n"
	      << "\t input aiiterator_instfirst[" << *aiiterator_instfirst << "]\n"
	      << "\t input aiiterator_instlast[" <<  *aiiterator_instlast << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES*/
 
  aiocrispmatrix_partition.initialize();

  for (uintidx li_j = 0;
       aiiterator_instfirst != aiiterator_instlast;
       ++aiiterator_instfirst, li_j++)
    {
    aiocrispmatrix_partition.setBit
      (function_instClass(*aiiterator_instfirst),li_j);
    }
  /*for ( uintidx li_j = 0; li_j < aivectorptinst_instances.size(); li_j++) {
    data::InstanceClass
	<T_FEATURE,
	 T_INSTANCES_CLUSTER_K,
	 T_CLUSTERIDX>
	*linst_instancesClassU = 
      (data::InstanceClass
	  <T_FEATURE,
	   T_INSTANCES_CLUSTER_K,
	   T_CLUSTERIDX>*)  
	 aivectorptinst_instances.at(li_j); 
       T_CLUSTERIDX lcidx_instClusterK = 
	 linst_instancesClassU->getClassIdx();
       aiocrispmatrix_partition.setBit
	(lcidx_instClusterK,li_j);
	}*/

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "clusteringop::initializeOptimal: OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output mat::CrispMatrix: aiocrispmatrix_partition["
	      << &aiocrispmatrix_partition << "]\n";
    aiocrispmatrix_partition.print();
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

/*! \fn T_CLUSTERIDX getCentroids(mat::MatrixRow<T_FEATURE> &aomatrixt_centroids, mat::MatrixRow<T_FEATURE_SUM> &aomatrixt_sumInstancesCluster, std::vector<T_INSTANCES_CLUSTER_K> &aovectort_numInstancesInClusterK, mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aicrispmatrix_partition, mat::MatrixRow<T_FEATURE> &aimatrixt_instances)
  \brief Calculate centroids 
  \details Calculate centroids from the Crisp Matrix partition. Returns the sum of instances and number of instances per cluster. Also the number and null clusters
  \param aomatrixt_centroids an output mat::MatrixRow with centroids
  \param aomatrixt_sumInstancesCluster a mat::MatrixRow with the sum of instances per cluster
  \param aovectort_numInstancesInClusterK a std::vector with the number of instances per cluster
  \param aicrispmatrix_partition an input mat::CrispMatrix partition  of instances in clusters
  \param aimatrixt_instances an input mat::MatrixRow with the instances
*/
template < typename T_FEATURE,
	   typename T_FEATURE_SUM, 
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_BITSIZE,
	   typename T_CLUSTERIDX  
	   >
T_CLUSTERIDX 
getCentroids
(mat::MatrixRow<T_FEATURE>                &aomatrixt_centroids,
 mat::MatrixRow<T_FEATURE_SUM>          &aomatrixt_sumInstancesCluster,
 std::vector<T_INSTANCES_CLUSTER_K>       &aovectort_numInstancesInClusterK,
 mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aicrispmatrix_partition,
 mat::MatrixRow<T_FEATURE>                &aimatrixt_instances 
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::getCentroids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
  	      << "(output mat::MatrixRow<T_FEATURE>: aomatrixt_centroids["
	      << &aomatrixt_centroids << "]\n"
	      << " input  mat::CrispMatrix: aicrispmatrix_partition["
	      <<  &aicrispmatrix_partition << "]\n"
	      << " input  mat::MatrixRow<T_INSTANCES>&: aimatrixt_instances["
	      << &aimatrixt_instances << "]\n"
	      << ")\n";
  }
#endif /*__VERBOSE_YES*/  

  aomatrixt_sumInstancesCluster.initialize(T_FEATURE(0.0));
  
  interfacesse::copya
    (aovectort_numInstancesInClusterK.data(),
     T_INSTANCES_CLUSTER_K(0),
     aomatrixt_sumInstancesCluster.getNumRows()
     );
    
  /*\sum_{k=1}^{n} w_{ik} X_k 
   */
  mat::bitgemm
    (aomatrixt_sumInstancesCluster,
     aicrispmatrix_partition, 
     aimatrixt_instances
     );
  
  /*   \sum_{k=1}^{n} w_{ik}   */
  mat::getRowsNumBitOn
    (aovectort_numInstancesInClusterK,
     aicrispmatrix_partition
     );
  
  /*Compute centroides v_i */
  T_CLUSTERIDX locidx_numClusterNull;
 
   meanCentroids
    (locidx_numClusterNull,
     aomatrixt_centroids, 
     aomatrixt_sumInstancesCluster, 
     aovectort_numInstancesInClusterK
     );

#ifdef __VERBOSE_YES
   if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
	      << "locidx_numClusterNull = " << locidx_numClusterNull << '\n';
    std::ostringstream lostrstream_labelCentroids;
    lostrstream_labelCentroids << "<CENTROIDS MEAN:" << lpc_labelFunc;
    aomatrixt_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
    if ( locidx_numClusterNull > 0 ) {
	std::cout << "\nERROR: solution is not valid, cluster without instances";
    }
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return locidx_numClusterNull;
}


/*! \fn void getPartition(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aobcrispmatrix_partition, mat::MatrixRow<T_FEATURE> &aimatrixt_instances, mat::MatrixRow<T_FEATURE> &aimatrixt_centroids, dist::Dist<T_DIST,T_FEATURE> &aifunc2p_dist)
\brief Get partition
  \details Calculate a partition crisp matrix partition with the nearest neighbor rule.
  \param aobcrispmatrix_partition an output mat::CrispMatrix partition  of instances in clusters.
  \param aimatrixt_instances an input mat::MatrixRow with the instances
  \param aimatrixt_centroids an input mat::MatrixRow with centroids
  \param aifunc2p_dist an object of type dist::Dist to calculate distances 
*/
template <typename T_FEATURE, 
	  typename T_CLUSTERIDX, //-1, 0, 1, .., K
	  typename T_BITSIZE,
	  typename T_DIST // _CENTROIDINSTANCES
	  >
void
getPartition
(mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> &aobcrispmatrix_partition,
 mat::MatrixRow<T_FEATURE>                &aimatrixt_instances,
 mat::MatrixRow<T_FEATURE>                &aimatrixt_centroids,
 dist::Dist<T_DIST,T_FEATURE>             &aifunc2p_dist
) 
{
  T_DIST  lT_distMinCentInst;

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::getPartition";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output mat::CrispMatrix: aobcrispmatrix_partition["
	      << &aobcrispmatrix_partition << "]\n"
	      << "\t input  mat::MatrixRow<T_INSTANCES>: aimatrixt_instances[" << &aimatrixt_instances << "]\n"
	      << "\t input  mat::MatrixRow<T_CENTROIDS>: aimatrixt_centroids[" <<  &aimatrixt_centroids << "]\n"
	      << "\t input  aifunc2p_dista[" << &aifunc2p_dist << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  for (uintidx lui_i = 0; lui_i < aimatrixt_instances.getNumRows(); lui_i++) {
    T_FEATURE         *larrayrowt_instances = aimatrixt_instances.getRow(lui_i); 
    T_CLUSTERIDX lmgidx_nearestCluster = 
      nearest::NN
      <T_CLUSTERIDX,T_FEATURE,T_DIST>
      (lT_distMinCentInst,
       aimatrixt_centroids,
       larrayrowt_instances,
       aifunc2p_dist
       );  
    aobcrispmatrix_partition.setMember(lui_i, lmgidx_nearestCluster);
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    aobcrispmatrix_partition.print(std::cout,lpc_labelFunc,',',';');
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

}

} /*END namespace crispmatrix*/


#endif /*CRISPMATRIX_CLUSTERING_HPP*/
