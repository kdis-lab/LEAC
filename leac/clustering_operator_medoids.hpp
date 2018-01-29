/*! \file clustering_operator_medoids.hpp
 *
 * \brief clustering operator medoids
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CLUSTERING_OPERATOR_MEDOIDS_HPP
#define __CLUSTERING_OPERATOR_MEDOIDS_HPP

#include <utility>      // std::pair
#include "dist_matrix_dissimilarity.hpp"
#include "instance.hpp"
#include "partition_linked.hpp"


#include "verbose_global.hpp"

/*! \namespace clusteringop
  \brief Clustering operators
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace clusteringop {
  

/*nearestRepresentative:
   
   Asigna cada instancia a un grupo mas cercano y calcula la suma 
   de las distancias de cada Medoid mj con sus instancias. utilizando
   PartitionLinked
   
   \cite{Sheng:Xiaohui:GAclusteringMedoid:HKA:2004}

 */
template < typename T_DIST,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
std::tuple< T_DIST, ds::PartitionLinked<T_CLUSTERIDX> >
nearestRepresentative
(uintidx                         *aiarrayuiidx_medoids,
 T_CLUSTERIDX                    aicidx_numClusterK,
 const mat::MatrixTriang<T_DIST> &aimatrixtriagt_dissimilarity
 ) 
{ 
#ifdef __VERBOSE_YES
  const char* lpc_label = "clusteringop::clusterNearestRepresentative:";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_label
              << "  IN(" << geiinparam_verbose << ")\n"
	      << "\t input T_INSTANCES_IDX: *aiarrayuiidx_medoids ["
	      << aiarrayuiidx_medoids << "]\n"
	      << "\t input mat::MatrixTriang: aimatrixtriagt_dissimilarity["  
	      <<  &aimatrixtriagt_dissimilarity << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES 

  T_DIST                     lot_cost = T_DIST(0);

  //INITIALIZE PARTITION LINK
  ds::PartitionLinked<T_CLUSTERIDX>
    lopartlink_partition
    (aimatrixtriagt_dissimilarity.getNumRows(),
     (uintidx) aicidx_numClusterK
     );
  
  T_DIST  lrt_distMinMedoidsInst;
  for (uintidx luintidx_i = 0; luintidx_i < aimatrixtriagt_dissimilarity.getNumRows(); luintidx_i++) {
    // FIND MEDOIDS NEAREST 
    T_CLUSTERIDX lmgidx_j = 
      nearest::medoidsNN
      <T_CLUSTERIDX,T_DIST>
      (lrt_distMinMedoidsInst,
       luintidx_i,
       aiarrayuiidx_medoids,
       aicidx_numClusterK,
       aimatrixtriagt_dissimilarity
       );

    lopartlink_partition.addInstanceToCluster(lmgidx_j,luintidx_i);
  
    lot_cost += (T_DIST) lrt_distMinMedoidsInst;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
     std::cout << lpc_label
	      << " OUT(" << geiinparam_verbose << ")\n"
	      << "\noutput PartitionLinked: [" << &lopartlink_partition << "]\n";
  
    lopartlink_partition.print
      (std::cout,
       lpc_label
       );
    std::cout  << "output lot_cost: " << lot_cost 
	       << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_tuple
    (  lot_cost
     , lopartlink_partition 
     );
  
}


/*clusterNearestRepresentative:
  Asigna cada instancia a un grupo mas cercano 
  y calcula la suma de las distancias de cada Medoid mj con sus
  instancias. Utilizando un arreglo de pertenencia.

  \cite{Kaufman:Rousseeuw:Book:ClusterAnalysis:1990}

*/

template < typename T_DIST,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
std::pair<T_DIST,uintidx>
nearestRepresentative
(T_CLUSTERIDX                     *aoarraycidx_memberShip,
 const uintidx                    *aiarrayuiidx_medoids,
 const T_CLUSTERIDX               aicidx_numKMedoids,
 const mat::MatrixTriang<T_DIST>  &aimatrixtriagrt_dissimilarity
 ) 
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::nearestRepresentative";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\t(output T_CLUSTERIDX*: aoarraycidx_memberShip[" 
	      << aoarraycidx_memberShip << "]\n"
              << "\t  input uintidx* aiarrayuiidx_medoids[" <<  aiarrayuiidx_medoids << "]\n"
	      << "\t  input T_CLUSTERIDX: aicidx_numKMedoids = "
	      << aicidx_numKMedoids << '\n'
	      << "\t  input  mat::MatrixTriang: aimatrixtriagrt_dissimilarity["  
	      <<  &aimatrixtriagrt_dissimilarity << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES */

  T_DIST   lot_cost = T_DIST(0);
  uintidx  louintidx_threshold = 0;
  T_DIST  lrt_distMinMedoidsInst;
  
  for (uintidx lui = 0; lui < aimatrixtriagrt_dissimilarity.getNumRows(); lui++) {
    /* FIND MEDOIDS NEAREST */
    T_CLUSTERIDX lmgidx_j = 
      nearest::medoidsNN
      (lrt_distMinMedoidsInst,
       lui,
       aiarrayuiidx_medoids,
       aicidx_numKMedoids,
       aimatrixtriagrt_dissimilarity
       );
    if (aoarraycidx_memberShip[lui] != lmgidx_j) {
      ++louintidx_threshold;
      aoarraycidx_memberShip[lui] = lmgidx_j;
    }
    lot_cost += (T_DIST) lrt_distMinMedoidsInst;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
              << " lot_cost = " << lot_cost
              << " louintidx_threshold = " << louintidx_threshold
	      << '\n';
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:" << geverbosepc_labelstep
				<< ':' << lpc_labelFunc
				<< ":aoarraycidx_memberShip<>["
				<< aoarraycidx_memberShip<< ']';
    inout::containerprint
      (aoarraycidx_memberShip,
       aoarraycidx_memberShip + aimatrixtriagrt_dissimilarity.getNumRows(),
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout << std::endl;	      
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return  std::make_pair(lot_cost,louintidx_threshold);

}


/*computeTotalCostSwap:
 */
template < typename T_DIST,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
T_DIST
computeTotalCostSwap
(const uintidx                    aiuiidx_medoidRand, /*IS A INDEX TO INSTANCE*/ 
 const T_CLUSTERIDX               *aiarraycidx_memberShip,
 uintidx                          *aoarrayuiidx_medoids,
 const T_CLUSTERIDX               aicidx_numClusterK,
 const mat::MatrixTriang<T_DIST>  &aimatrixtriagrt_dissimilarity
 ) 
{
  T_DIST               lot_cost;
  T_DIST               lrt_minMedoidsRandToInstance_i;
  T_DIST               lrt_minCurrentClusterKToInstance_i;
  T_CLUSTERIDX  lcidx_clusterKNonRepre; 
  uintidx              luintidx_previousMedoids;

  lot_cost = T_DIST(0);
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "pamkmedoid_computeTotalCostSwap: IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t(output lot_cost = " << lot_cost << '\n'
	      << "\t input T_INSTANCES_IDX: aiuiidx_medoidRand = " 
	      << aiuiidx_medoidRand << '\n'
	      << "\t input T_CLUSTERIDX*: aiarraycidx_memberShip[" 
	      << aiarraycidx_memberShip << "]\n"
	      << "\t input T_INSTANCES_IDX*: aoarrayuiidx_medoids["
	      << aoarrayuiidx_medoids << "]\n"
	      << "\t input T_CLUSTERIDX: aicidx_numClusterK = " 
	      << aicidx_numClusterK << '\n'
	      << "\t input  mat::MatrixTriang: aimatrixtriagrt_dissimilarity["  
	      <<  &aimatrixtriagrt_dissimilarity << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES */

  lcidx_clusterKNonRepre = aiarraycidx_memberShip[ aiuiidx_medoidRand ];
  luintidx_previousMedoids = aoarrayuiidx_medoids[lcidx_clusterKNonRepre];
  aoarrayuiidx_medoids[lcidx_clusterKNonRepre] = aiuiidx_medoidRand;

  for (uintidx lui_i = 0; lui_i < aimatrixtriagrt_dissimilarity.getNumRows(); lui_i++) {
    if ( lcidx_clusterKNonRepre == aiarraycidx_memberShip[lui_i] ) {
      /* FIND DISTANCE NEAREST CONSY */
      nearest::medoidsNN 
	(lrt_minCurrentClusterKToInstance_i,
	 lui_i,
	 aoarrayuiidx_medoids,
	 aicidx_numClusterK,
	 aimatrixtriagrt_dissimilarity
	 );
    }
    else {
      uintidx luintidx_clusterKcurrentInsti = 
	aoarrayuiidx_medoids[aiarraycidx_memberShip[lui_i]];
 
      lrt_minCurrentClusterKToInstance_i = 
	( lui_i > luintidx_clusterKcurrentInsti)?
	aimatrixtriagrt_dissimilarity( lui_i, luintidx_clusterKcurrentInsti ):
	aimatrixtriagrt_dissimilarity( luintidx_clusterKcurrentInsti, lui_i );
      
      lrt_minMedoidsRandToInstance_i = (lui_i > aiuiidx_medoidRand)?
	aimatrixtriagrt_dissimilarity( lui_i, aiuiidx_medoidRand ):
	aimatrixtriagrt_dissimilarity( aiuiidx_medoidRand, lui_i );

      if ( lrt_minMedoidsRandToInstance_i < lrt_minCurrentClusterKToInstance_i )
	lrt_minCurrentClusterKToInstance_i = lrt_minMedoidsRandToInstance_i;
    }
    lot_cost += (T_DIST) lrt_minCurrentClusterKToInstance_i;
  }
  
  aoarrayuiidx_medoids[lcidx_clusterKNonRepre] = luintidx_previousMedoids;

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "pamkmedoid_computeTotalCostSwap: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << " lot_cost = " << lot_cost << '\n';
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lot_cost;

}

  
template  < typename T_DIST,
	    typename T_CLUSTERIDX //-1, 0, 1, .., K
	   >
T_DIST
computeCostInstanceClusterJ
(uintidx                     aiuintidx_medoids, 
 T_CLUSTERIDX         aicidx_Cj,
  ds::PartitionLinked        /*Relation each object in X to the cluster Cj*/
 <T_CLUSTERIDX>       &aipartlink_partition,           
 mat::MatrixTriang
 <T_DIST> &aimatrixtriagt_dissimilarity
 )
{
  T_DIST  loT_cost;
  T_DIST  lrt_instancesij;
  ds::IteratorPartitionLinked<T_CLUSTERIDX>  literpart_j(&aipartlink_partition);
  

  loT_cost = T_DIST(0);
  for (literpart_j.begin(aicidx_Cj);literpart_j.end();literpart_j.next() ) {
      lrt_instancesij = (aiuintidx_medoids>literpart_j.getValue())?
	aimatrixtriagt_dissimilarity
	(aiuintidx_medoids,literpart_j.getValue()):
	aimatrixtriagt_dissimilarity
	(literpart_j.getValue(),aiuintidx_medoids);
      loT_cost += (T_DIST) lrt_instancesij;
    }

  return loT_cost;
}

/*medoidsCenterPnearestNeighbors:
  \cite{Sheng:Xiaohui:GAclusteringMedoid:HKA:2004}
*/
template  <typename T_DIST,
           typename T_CLUSTERIDX //-1, 0, 1, .., K
	   >
void
centerPnearestNeighbors
(uintidx                            *aioarrayuintidx_medoids, 
 const uintidx                      aiuintidx_nearestNeighborsP,
 ds::PartitionLinked<T_CLUSTERIDX>  &aipartlink_partition, /*Relation each object in X to the cluster Cj*/             
 mat::MatrixTriang<T_DIST>          &aimatrixtriagt_dissimilarity
 )
{
  /*BUSCAL EL MEDOID QUE ESTE MAS EN EL CENTRO PARA CADA GRUPO*/

  std::vector<std::pair<uintidx,T_DIST> > lvectorpairT_subsetCjNearestMj;
  std::pair<uintidx,T_DIST> lpair_idxInstanceDistMj
    (UINTIDX_NIL, //DATATYPE_INSTANCE_IDX_NULL,
     std::numeric_limits<T_DIST>::max()
     );
  
  ds::IteratorPartitionLinked<T_CLUSTERIDX>   literpart_j(&aipartlink_partition);
  T_DIST                                      lrt_minSumij;
  uintidx                                     lidxinst_medoidsCenteri;
  mat::BitArray<unsigned long>   
    lbitarray_medoidEvaluate
    (aimatrixtriagt_dissimilarity.getNumRows());
  
  bool lb_medoidChange;

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::medoidsCenterPnearestNeighbors";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "\n(output uintidx*: aioarrayuintidx_medoids[" 
	      << aioarrayuintidx_medoids << ']'
              << "\n input  uintidx: aiuintidx_nearestNeighborsP = " << aiuintidx_nearestNeighborsP
	      << "\n input  ds::PartitionLinked: aipartlink_partition " << &aipartlink_partition
	      << "\n input  MatrixTriang: aiinstf_instances[" <<  &aimatrixtriagt_dissimilarity
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  lbitarray_medoidEvaluate.initialize();
  
  lvectorpairT_subsetCjNearestMj.reserve(aiuintidx_nearestNeighborsP);
  for (uintidx lui_i = 0; lui_i < aiuintidx_nearestNeighborsP; lui_i++) {
    lvectorpairT_subsetCjNearestMj.push_back
      (std::make_pair
       (UINTIDX_NIL,
       std::numeric_limits<T_DIST>::max()
       )
       );
  }

  T_CLUSTERIDX lcidx_numClusterK = aipartlink_partition.getNumPartitions();
  
  for ( T_CLUSTERIDX lcidx_Cj = 0; 
	lcidx_Cj < lcidx_numClusterK; 
	lcidx_Cj++
	) 
    {
      
      do {

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax  ) {
	  std::cout << "Medoids C_" << lcidx_Cj << ": " << aioarrayuintidx_medoids[lcidx_Cj] << '\n';
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
	lb_medoidChange = false;
	//NEAREST NEIGHBOORS P
	lidxinst_medoidsCenteri = aioarrayuintidx_medoids[lcidx_Cj];
	lrt_minSumij = 
	  computeCostInstanceClusterJ
	  (lidxinst_medoidsCenteri,
	   lcidx_Cj,
	   aipartlink_partition, 
	   aimatrixtriagt_dissimilarity
	   );
	lbitarray_medoidEvaluate.setBit(lidxinst_medoidsCenteri);
   
	uintidx lui_numCsubset = 0; /*T_INSTANCES_IDX*/
	for (typename std::vector<std::pair<uintidx,T_DIST> >::iterator
	       literpair_cjNearestMj = lvectorpairT_subsetCjNearestMj.begin();
	     literpair_cjNearestMj != lvectorpairT_subsetCjNearestMj.end();
	     ++literpair_cjNearestMj)
	  {
	    literpair_cjNearestMj->first = UINTIDX_NIL;
	    // DATATYPE_INSTANCE_IDX_NULL; 
	    literpair_cjNearestMj->second = 
	      std::numeric_limits<T_DIST>::max(); 
	  }
	/*SEARCH Csubset p Nearest neighbors*/
	for ( literpart_j.begin(lcidx_Cj); literpart_j.end(); literpart_j.next() ) {

	  if ( lbitarray_medoidEvaluate.getBit(literpart_j.getValue()) ) {
	    lpair_idxInstanceDistMj.first  = literpart_j.getValue();
	    lpair_idxInstanceDistMj.second = 
	      (aioarrayuintidx_medoids[lcidx_Cj]>literpart_j.getValue())?
	      aimatrixtriagt_dissimilarity
	      (aioarrayuintidx_medoids[lcidx_Cj],literpart_j.getValue()):
	      aimatrixtriagt_dissimilarity
	      (literpart_j.getValue(),aioarrayuintidx_medoids[lcidx_Cj]);
	    uintidx lui_iterCsubset = 0;
	    while (lui_iterCsubset < aiuintidx_nearestNeighborsP &&  lpair_idxInstanceDistMj.second 
		   > lvectorpairT_subsetCjNearestMj[lui_iterCsubset].second )
	      ++lui_iterCsubset;
	    if ( lui_iterCsubset < aiuintidx_nearestNeighborsP ) {
	      for (uintidx lui_i = aiuintidx_nearestNeighborsP-1; lui_i > lui_iterCsubset; lui_i--) {
		std::swap
		  (lvectorpairT_subsetCjNearestMj[lui_i],
		   lvectorpairT_subsetCjNearestMj[lui_i-1]
		   );
	      }
	      std::swap
		(lvectorpairT_subsetCjNearestMj[lui_iterCsubset],
		 lpair_idxInstanceDistMj
		 );
	      ++lui_numCsubset;
	    }
	  }
	} /*  FOR SEARCH Csubset p Nearest neighbors*/

	if ( lui_numCsubset >  aiuintidx_nearestNeighborsP )
	  lui_numCsubset =  aiuintidx_nearestNeighborsP;
	
#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "Csubset p Nearest neighbors = [";
	  if ( lui_numCsubset > 0 ) {
	    for ( uintidx lui_i = 0; lui_i < lui_numCsubset -1; lui_i++ ) {
	      std::cout << lvectorpairT_subsetCjNearestMj[lui_i].first << ":" 
			<< lvectorpairT_subsetCjNearestMj[lui_i].second  
			<< "\t";
	    }
	    std::cout << lvectorpairT_subsetCjNearestMj[lui_numCsubset -1].first << ":" 
		      << lvectorpairT_subsetCjNearestMj[lui_numCsubset -1].second;
	  }  
	  std::cout << "]\n";
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      
	for ( uintidx lui_i = 0; lui_i < lui_numCsubset; lui_i++ ) {

	  T_DIST lrt_sumij =
	    computeCostInstanceClusterJ
	    (lvectorpairT_subsetCjNearestMj[lui_i].first,
	     lcidx_Cj,
	     aipartlink_partition, 
	     aimatrixtriagt_dissimilarity
	     );
	 
#ifdef __VERBOSE_YES
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << "SumDist = " << lrt_sumij << "\n";
	  }
	  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

	  if  ( lrt_sumij < lrt_minSumij ) {
	    lidxinst_medoidsCenteri = lvectorpairT_subsetCjNearestMj[lui_i].first;
	    lrt_minSumij    = lrt_sumij;
	  } 
	} /* For i*/
	
	if ( aioarrayuintidx_medoids[lcidx_Cj] != lidxinst_medoidsCenteri ) {
	  lb_medoidChange = true;
	  aioarrayuintidx_medoids[lcidx_Cj] = lidxinst_medoidsCenteri;
       }

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "MinSumDist: " << lrt_minSumij 
		    << " lidxinst_medoidsCenteri: " << lidxinst_medoidsCenteri  << '\n';
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      } while ( lb_medoidChange );

    } /*FOR  K */


#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMedoids;
    lostrstream_labelMedoids << "<MEDOIDS:"
			     << lpc_labelFunc;
    inout::containerprint
      (aioarrayuintidx_medoids,
       aioarrayuintidx_medoids + (uintidx)lcidx_numClusterK,
       std::cout,
       lostrstream_labelMedoids.str().c_str(),
       ','
       ); 
     std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
  
}


/*! \fn void updateMedoids(uintidx *aoarrayuiidx_medoids, T_CLUSTERIDX aicidx_numClusterK, uintidx aiuiidx_nearestNeighborsP, mat::MatrixTriang<T_DIST>  &aimatrixtriagt_dissimilarity)
    \brief  Update medoids \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004
    \details For each cluster \f$C_j\f$ finds the most representative object
    -# Assign each object in \f$x_i\f$ to the cluster \f$C_j\f$ with the closest medoid
    -# For each cluster \f$C_j\f$, repeat until the medoid does not change 
         - Choose a subset \f$C_{subset}\f$ in Cj the corresponds to \f$m_j\f$ and its \f$p\f$ nearest neighbors of \f$m_j\f$.
         - Calulate the new medoid
            \f[
                m_{j}^* = { \hbox{arg min} \atop { x_i \in C_{subset}} } { \sum \atop { x_i' \in C_j} } \Vert x_i - x_i' \Vert 
            \f] 
         - if \f$m_j\f$ is different from \f$m_{j}^*\f$ replace with the new medoid
    -# Repat step 1 and 2 until \f$k\f$ medoids do not change    
    \param aoarrayuiidx_medoids an array with the indexes of the representative objects initially
    \param aicidx_numClusterK cluster number
    \param aiuiidx_nearestNeighborsP neatest neighbors 
    \param aimatrixtriagt_dissimilarity a distance matrix
 */
template < typename T_CLUSTERIDX,    //-1, 0, 1, .., K
           typename T_DIST
	   >
void updateMedoids
(uintidx                    *aoarrayuiidx_medoids,
 T_CLUSTERIDX               aicidx_numClusterK,
 uintidx                    aiuiidx_nearestNeighborsP,
 mat::MatrixTriang<T_DIST>  &aimatrixtriagt_dissimilarity
 )
{
  /*STEP 1. FIX THE NUMBER OF CLUSTER K AND THE NUMBER OF NEAREST NEIGHBORS P
   */

#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::updateMedoids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
 	      << "(output uintidx: [" << aoarrayuiidx_medoids << "]\n"
	      << " input T_CLUSTERIDX: aicidx_numClusterK = "
	      << aicidx_numClusterK << "]\n"
	      << " input  MatrixTriang: aimatrixtriagt_dissimilarity["  
	      <<  &aimatrixtriagt_dissimilarity << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES */

  
  /*STEP 2. RANDOMY
   */

  /*STEP 3. ASSIGN EACH OBJECT IN X TO THE CLUSTER Cj
   */

  T_DIST                    lt_cost = T_DIST(0);
  ds::PartitionLinked  
    <T_CLUSTERIDX>           lpartlink_partition;
  
  std::tie(lt_cost, lpartlink_partition) =
    nearestRepresentative
    (aoarrayuiidx_medoids,
     aicidx_numClusterK,
     aimatrixtriagt_dissimilarity
     );

  /*STEP 4.
   */
  centerPnearestNeighbors
    (aoarrayuiidx_medoids, 
     aiuiidx_nearestNeighborsP,
     lpartlink_partition,
     aimatrixtriagt_dissimilarity
     );
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMedoids;
    lostrstream_labelMedoids << "<MEDOIDS:"
			     << lpc_labelFunc;
    inout::containerprint
      (aoarrayuiidx_medoids,
       aoarrayuiidx_medoids + (uintidx)aicidx_numClusterK,
       std::cout,
       lostrstream_labelMedoids.str().c_str(),
       ','
       ); 
     std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
 
}

} /*END namespace medoids*/

#endif /*__CLUSTERING_OPERATOR_MEDOIDS_HPP*/
