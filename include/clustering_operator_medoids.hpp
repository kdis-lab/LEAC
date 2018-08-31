/*! \file clustering_operator_medoids.hpp
 *
 * \brief clustering operator medoids
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CLUSTERING_OPERATOR_MEDOIDS_HPP__
#define __CLUSTERING_OPERATOR_MEDOIDS_HPP__

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
  

/*! \fn std::tuple< T_DIST, ds::PartitionLinked<T_CLUSTERIDX> > nearestRepresentative (T_INSTANCEIDX *aioarrayui_idxMedoids, T_CLUSTERIDX aicidx_numKMedoids, const mat::MatrixTriang<T_DIST> &aimattriag_dissimilarity) 
   \brief nearestRepresentative
   \details Assign each instance to a closer group and calculate the sum of the distances of each Medoid mj with its instances
   \param aioarrayui_idxMedoids an array with the indexes of the representative objects initially
   \param aicidx_numKMedoids cluster number
   \param aimattriag_dissimilarity a distance matrix
 */
template < typename T_INSTANCEIDX,
           typename T_DIST,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
std::tuple< T_DIST, ds::PartitionLinked<T_CLUSTERIDX> >
nearestRepresentative
(T_INSTANCEIDX                   *aioarrayui_idxMedoids,
 T_CLUSTERIDX                    aicidx_numKMedoids,
 const mat::MatrixTriang<T_DIST> &aimattriag_dissimilarity
 ) 
{ 
#ifdef __VERBOSE_YES
  const char* lpc_label = "clusteringop::nearestRepresentative:";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_label
              << "  IN(" << geiinparam_verbose << ")\n"
	      << "\t input T_INSTANCES_IDX: *aioarrayui_idxMedoids ["
	      << aioarrayui_idxMedoids << "]\n"
	      << "\t input mat::MatrixTriang: aimattriag_dissimilarity["  
	      <<  &aimattriag_dissimilarity << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES 

  T_DIST               lor_cost = T_DIST(0);
  const T_INSTANCEIDX  lui_numInstances = T_INSTANCEIDX(aimattriag_dissimilarity.getNumRows());
  //Initialize PartitionLink
  ds::PartitionLinked<T_CLUSTERIDX>
    lopartlink_partition
     (lui_numInstances, 
     (T_INSTANCEIDX) aicidx_numKMedoids
     );
  
  T_DIST  lrt_distMinMedoidsInst;
  
  for (T_INSTANCEIDX luintidx_i = 0; luintidx_i < lui_numInstances; luintidx_i++) {
    //Find medoids nearest
    T_CLUSTERIDX lcidx_mj = 
      nearest::medoidsNN
      <T_CLUSTERIDX,T_DIST>
      (lrt_distMinMedoidsInst,
       luintidx_i,
       aioarrayui_idxMedoids,
       aicidx_numKMedoids,
       aimattriag_dissimilarity
       );

    lopartlink_partition.addInstanceToCluster(lcidx_mj,luintidx_i);
  
    lor_cost += (T_DIST) lrt_distMinMedoidsInst;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
     std::cout << lpc_label
	      << " OUT(" << geiinparam_verbose << ")\n"
	      << "output PartitionLinked: [" << &lopartlink_partition << "]\n";
    lopartlink_partition.print
      (std::cout,
       lpc_label
       );
    std::cout  << "\noutput lor_cost: " << lor_cost 
	       << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return std::make_tuple
    (  lor_cost
     , lopartlink_partition 
     );
  
}

/*! \fn std::pair<T_DIST,T_INSTANCEIDX> nearestRepresentative (T_CLUSTERIDX *aoarraycidx_memberShip, const T_INSTANCEIDX *aioarrayui_idxMedoids, const T_CLUSTERIDX aicidx_numKMedoids, const mat::MatrixTriang<T_DIST>  &aimattriag_dissimilarity)
   \brief nearestRepresentative \cite Kaufman:Rousseeuw:Book:ClusterAnalysis:1990
   \details Assign each instance to a closer group and calculate the sum of the distances of each Medoid mj with its instances
   \param aoarraycidx_memberShip an array of out with the labels belong to the cluster 
   \param aioarrayui_idxMedoids an array of input with the indexes of the representative objects
   \param aicidx_numKMedoids cluster number
   \param aimattriag_dissimilarity a distance matrix
 */
template < typename T_INSTANCEIDX,
           typename T_DIST,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
std::pair<T_DIST,T_INSTANCEIDX>
nearestRepresentative
(T_CLUSTERIDX                     *aoarraycidx_memberShip,
 const T_INSTANCEIDX              *aioarrayui_idxMedoids,
 const T_CLUSTERIDX               aicidx_numKMedoids,
 const mat::MatrixTriang<T_DIST>  &aimattriag_dissimilarity
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
              << "\t  input T_INSTANCEIDX* aioarrayui_idxMedoids[" 
	      <<  aioarrayui_idxMedoids << "]\n"
	      << "\t  input T_CLUSTERIDX: aicidx_numKMedoids = "
	      << aicidx_numKMedoids << '\n'
	      << "\t  input  mat::MatrixTriang: aimattriag_dissimilarity["  
	      <<  &aimattriag_dissimilarity << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES */

  T_DIST              lor_cost = T_DIST(0);
  T_INSTANCEIDX       louintidx_threshold = 0;  
  T_DIST              lrt_distMinMedoidsInst;
  const T_INSTANCEIDX lui_numInstances = 
    T_INSTANCEIDX(aimattriag_dissimilarity.getNumRows());

  for (T_INSTANCEIDX lui = 0; lui < lui_numInstances; lui++) {
    /* Find medoids nearest */
    T_CLUSTERIDX lcidx_mj = 
      nearest::medoidsNN
      (lrt_distMinMedoidsInst,
       lui,
       aioarrayui_idxMedoids,
       aicidx_numKMedoids,
       aimattriag_dissimilarity
       );
    if (aoarraycidx_memberShip[lui] != lcidx_mj) {
      ++louintidx_threshold;
      aoarraycidx_memberShip[lui] = lcidx_mj;
    }
    lor_cost += (T_DIST) lrt_distMinMedoidsInst;
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ')'
              << " lor_cost = " << lor_cost
              << " louintidx_threshold = " << louintidx_threshold
	      << '\n';
    std::ostringstream lostrstream_labelMemberShip;
    lostrstream_labelMemberShip << "<MEMBERCLUSTER:" << geverbosepc_labelstep
				<< ':' << lpc_labelFunc
				<< ":aoarraycidx_memberShip<>["
				<< aoarraycidx_memberShip<< ']';
    inout::containerprint
      (aoarraycidx_memberShip,
       aoarraycidx_memberShip + aimattriag_dissimilarity.getNumRows(),
       std::cout,
       lostrstream_labelMemberShip.str().c_str(),
       ','
       );
    std::cout << std::endl;	      
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return  std::make_pair(lor_cost,louintidx_threshold);

}

/*! \fn T_DIST computeTotalCostSwap (const T_INSTANCEIDX aiuiidx_medoidRand, const T_CLUSTERIDX *aiarraycidx_memberShip, T_INSTANCEIDX *aoarrayuiidx_medoids, const T_CLUSTERIDX aicidx_numKMedoids, const mat::MatrixTriang<T_DIST>  &aimattriag_dissimilarity)
   \brief computeTotalCostSwap
   \details calculates the cost of changing a randomly selected medoid to a defined one
   \param aiuiidx_medoidRand an index of a randomly selected medoid to investigate if the cost improves
   \param aiarraycidx_memberShip an aioarraycidx_memberShip an array of indexes belonging to previously assigned or UNKNOWN_CLUSTER_IDX
   \param aoarrayuiidx_medoids an array of output with the indexes of the representative objects
   \param aicidx_numKMedoids cluster number
   \param aimattriag_dissimilarity a distance matrix
 */
template < typename T_INSTANCEIDX,
           typename T_DIST,
	   typename T_CLUSTERIDX    //-1, 0, 1, .., K
	   >
T_DIST
computeTotalCostSwap
(const T_INSTANCEIDX              aiuiidx_medoidRand, /*IS A INDEX TO INSTANCE*/ 
 const T_CLUSTERIDX               *aiarraycidx_memberShip,
 T_INSTANCEIDX                    *aoarrayuiidx_medoids,
 const T_CLUSTERIDX               aicidx_numKMedoids,
 const mat::MatrixTriang<T_DIST>  &aimattriag_dissimilarity
 ) 
{
  T_DIST          lor_cost;
  T_DIST          lrt_minMedoidsRandToInstance_i;
  T_DIST          lrt_minCurrentClusterKToInstance_i;
  T_CLUSTERIDX    lcidx_clusterKNonRepre; 
  T_INSTANCEIDX   luintidx_previousMedoids;

  lor_cost = T_DIST(0);
  
#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "pamkmedoid_computeTotalCostSwap: IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t(output lor_cost = " << lor_cost << '\n'
	      << "\t input T_INSTANCES_IDX: aiuiidx_medoidRand = " 
	      << aiuiidx_medoidRand << '\n'
	      << "\t input T_CLUSTERIDX*: aiarraycidx_memberShip[" 
	      << aiarraycidx_memberShip << "]\n"
	      << "\t input T_INSTANCES_IDX*: aoarrayuiidx_medoids["
	      << aoarrayuiidx_medoids << "]\n"
	      << "\t input T_CLUSTERIDX: aicidx_numKMedoids = " 
	      << aicidx_numKMedoids << '\n'
	      << "\t input  mat::MatrixTriang: aimattriag_dissimilarity["  
	      <<  &aimattriag_dissimilarity << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES */

  lcidx_clusterKNonRepre = aiarraycidx_memberShip[ aiuiidx_medoidRand ];
  luintidx_previousMedoids = aoarrayuiidx_medoids[lcidx_clusterKNonRepre];
  aoarrayuiidx_medoids[lcidx_clusterKNonRepre] = aiuiidx_medoidRand;

  for (T_INSTANCEIDX lui_i = 0; lui_i < aimattriag_dissimilarity.getNumRows(); lui_i++) {
    if ( lcidx_clusterKNonRepre == aiarraycidx_memberShip[lui_i] ) {
      /* FIND DISTANCE NEAREST CONSY */
      nearest::medoidsNN 
	(lrt_minCurrentClusterKToInstance_i,
	 lui_i,
	 aoarrayuiidx_medoids,
	 aicidx_numKMedoids,
	 aimattriag_dissimilarity
	 );
    }
    else {
      T_INSTANCEIDX luintidx_clusterKcurrentInsti = 
	aoarrayuiidx_medoids[aiarraycidx_memberShip[lui_i]];
 
      lrt_minCurrentClusterKToInstance_i = 
	( lui_i > luintidx_clusterKcurrentInsti)?
	aimattriag_dissimilarity( lui_i, luintidx_clusterKcurrentInsti ):
	aimattriag_dissimilarity( luintidx_clusterKcurrentInsti, lui_i );
      
      lrt_minMedoidsRandToInstance_i = (lui_i > aiuiidx_medoidRand)?
	aimattriag_dissimilarity( lui_i, aiuiidx_medoidRand ):
	aimattriag_dissimilarity( aiuiidx_medoidRand, lui_i );

      if ( lrt_minMedoidsRandToInstance_i < lrt_minCurrentClusterKToInstance_i )
	lrt_minCurrentClusterKToInstance_i = lrt_minMedoidsRandToInstance_i;
    }
    lor_cost += (T_DIST) lrt_minCurrentClusterKToInstance_i;
  }
  
  aoarrayuiidx_medoids[lcidx_clusterKNonRepre] = luintidx_previousMedoids;

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "pamkmedoid_computeTotalCostSwap: OUT"
	      << '(' << geiinparam_verbose << ')'
	      << " lor_cost = " << lor_cost << '\n';
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lor_cost;

}

/*! \fn T_DIST computeCostInstanceClusterJ (T_INSTANCEIDX aiuintidx_medoids, T_CLUSTERIDX aicidx_Cj, ds::PartitionLinked<T_CLUSTERIDX> &aipartlink_partition, mat::MatrixTriang<T_DIST>  &aimattriag_dissimilarity)
   \brief computeCostInstanceClusterJ
   \details calculates the cost or sum of Euclidean distances from the most representative to all the objects belonging to the cluster j
   \param aiuintidx_medoids an index of a medoid to investigate the cost
   \param aicidx_Cj an integer with the cluster index j
   \param aipartlink_partition a partition::Partition of objects in clusters
   \param aimattriag_dissimilarity a distance matrix
 */  
template  < typename T_INSTANCEIDX,
            typename T_DIST,
	    typename T_CLUSTERIDX //-1, 0, 1, .., K
	   >
T_DIST
computeCostInstanceClusterJ
(T_INSTANCEIDX              aiuintidx_medoids, 
 T_CLUSTERIDX               aicidx_Cj,
  ds::PartitionLinked       /*Relation each object in X to the cluster Cj*/
 <T_CLUSTERIDX>             &aipartlink_partition,           
 mat::MatrixTriang<T_DIST>  &aimattriag_dissimilarity
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::computeCostInstanceClusterJ";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n" 
	      << "As more respresentatico: " << aiuintidx_medoids 
	      << " in C_" << aicidx_Cj
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  T_DIST  lor_costCj;
  T_DIST  lrt_instancesij;
  ds::IteratorPartitionLinked<T_CLUSTERIDX>  literpart_j(&aipartlink_partition);
  
  lor_costCj = T_DIST(0);
  for (literpart_j.begin(aicidx_Cj);literpart_j.end();literpart_j.next() ) {
      lrt_instancesij = (aiuintidx_medoids>literpart_j.getValue())?
	aimattriag_dissimilarity
	(aiuintidx_medoids,literpart_j.getValue()):
	aimattriag_dissimilarity
	(literpart_j.getValue(),aiuintidx_medoids);
      lor_costCj += (T_DIST) lrt_instancesij;
    }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      <<": OUT(" << geiinparam_verbose << ')' 
	      << " lor_costCj = " << lor_costCj << '\n';
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lor_costCj;
}


/*! \fn T_DIST updateMedoids(T_INSTANCEIDX  *aiui_idxInstMedoids, const T_CLUSTERIDX aicidx_numKMedoids, const T_INSTANCEIDX aiu_nearestNeighborsP, mat::MatrixTriang<T_DIST> &aimattriag_dissimilarity)
    \brief  Update medoids, local search heuristic  \cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004
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
    \param aioarrayui_idxMedoids an array with the indexes of the representative objects initially
    \param aicidx_numKMedoids cluster number
    \param aiu_nearestNeighborsP size of subset P of search neighbors 
    \param aimattriag_dissimilarity a distance matrix
 */
template  <typename T_INSTANCEIDX,
           typename T_DIST,
           typename T_CLUSTERIDX //-1, 0, 1, .., K
	   >
T_DIST
updateMedoids
(T_INSTANCEIDX                      *aioarrayui_idxMedoids, 
 const T_CLUSTERIDX                 aicidx_numKMedoids,
 const T_INSTANCEIDX                aiu_nearestNeighborsP,
 mat::MatrixTriang<T_DIST>          &aimattriag_dissimilarity
 )
{
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "clusteringop::updateMedoids";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMedoids;
    lostrstream_labelMedoids << "<MEDOIDS:"
			     << lpc_labelFunc;
    inout::containerprint
      (aioarrayui_idxMedoids,
       aioarrayui_idxMedoids + (T_INSTANCEIDX) aicidx_numKMedoids,
       std::cout,
       lostrstream_labelMedoids.str().c_str(),
       ','
       ); 
    std::cout << "\n input  T_INSTANCEIDX: aiu_nearestNeighborsP = " 
	      << aiu_nearestNeighborsP
	      << "\n input  MatrixTriang: aiinstf_instances[" <<  &aimattriag_dissimilarity
	      << "\n)"
	      << std::endl;
    T_DIST lt_objetiveFunc = 
      um::SSEMedoid
      (aioarrayui_idxMedoids, 
       aicidx_numKMedoids,
       aimattriag_dissimilarity
       );
    std::cout << "lt_objetiveFunc: " << lt_objetiveFunc << std::endl;
  }
#endif /*__VERBOSE_YES*/

  /*Vector to store the Mj closest to Cj
   */
  std::vector<T_INSTANCEIDX> lvectorui_subsetCjNearestMj(aiu_nearestNeighborsP);
  mat::BitArray<unsigned long>   
    lbitarray_medoidEvaluate
    (aimattriag_dissimilarity.getNumRows());
  T_DIST lor_cost = T_DIST(0);
  bool lb_medoidChange;

  do {

    lb_medoidChange = false;

    /* Assign each object in X to the cluster Cj with the
      closest medoid under Euclidean distance metric.
    */
    ds::PartitionLinked  <T_CLUSTERIDX>  lpartlink_partition;
    std::tie(lor_cost, lpartlink_partition) =
      nearestRepresentative
      (aioarrayui_idxMedoids,
       aicidx_numKMedoids,
       aimattriag_dissimilarity
       );
    lbitarray_medoidEvaluate.initialize();
    
    ds::IteratorPartitionLinked<T_CLUSTERIDX>   literpart_j(&lpartlink_partition);

     /*Update k medoids. For j=1 to the number de cluster k do
      */
    for ( T_CLUSTERIDX lcidx_Cj = 0; lcidx_Cj < aicidx_numKMedoids;  lcidx_Cj++ ) 
      {

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax  ) {
	  std::cout << "Search for Medoids: " << aioarrayui_idxMedoids[lcidx_Cj] 
                    << " in C_" << lcidx_Cj << '\n';
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/
      

	/*Which have no been evaluated
	  before current iteration
	 */
	lbitarray_medoidEvaluate.setBit(aioarrayui_idxMedoids[lcidx_Cj]);

	/*a) Within the cluster Cj, choose a subset Cs
	  that corresponds to m, and its p nearest
	  neighbors (which have no been evaluated
	  before current iteration) of m,.
	*/
	for (auto& xkinSubsetCj: lvectorui_subsetCjNearestMj) {
	  /*SEARCH Csubset p Nearest neighbors*/
	  xkinSubsetCj = UINTIDX_NIL;
	  T_DIST lt_distMinXkinSubsetCj = std::numeric_limits<T_DIST>::max();
	  for ( literpart_j.begin(lcidx_Cj); literpart_j.end(); literpart_j.next() ) {
	    if ( !lbitarray_medoidEvaluate.getBit(literpart_j.getValue()) ) {
	      T_DIST lt_distmjtoP = 
		(aioarrayui_idxMedoids[lcidx_Cj]>literpart_j.getValue())?
		aimattriag_dissimilarity
		(aioarrayui_idxMedoids[lcidx_Cj],literpart_j.getValue()):
		aimattriag_dissimilarity
		(literpart_j.getValue(),aioarrayui_idxMedoids[lcidx_Cj]);
	      if ( lt_distmjtoP < lt_distMinXkinSubsetCj ) {
		xkinSubsetCj  = literpart_j.getValue();
		lt_distMinXkinSubsetCj = lt_distmjtoP;
		
	      }
	    }
	  }
	  if (xkinSubsetCj != UINTIDX_NIL) {
	    lbitarray_medoidEvaluate.setBit(xkinSubsetCj);
	  }
	}

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {

	  std::ostringstream lostrstream_subsetCjNearestMj;
	  lostrstream_subsetCjNearestMj << "<SUBSETCJNEARESTMJ:"
					<< lpc_labelFunc;
	  inout::containerprint
	    (lvectorui_subsetCjNearestMj.begin(),
	     lvectorui_subsetCjNearestMj.end(),
	     std::cout,
	     lostrstream_subsetCjNearestMj.str().c_str(),
	     ','
	     ); 
	  std::cout << std::endl;
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

	/*b) Calculate cost for the new medoid
	 */
	T_INSTANCEIDX lui_idxInstNewMj = aioarrayui_idxMedoids[lcidx_Cj];
	T_DIST lrt_MjSumD = 
	  computeCostInstanceClusterJ
	  (lui_idxInstNewMj,
	   lcidx_Cj,
	   lpartlink_partition, 
	   aimattriag_dissimilarity
	   );
	for (auto& xkinSubsetCj: lvectorui_subsetCjNearestMj) {
	  if ( xkinSubsetCj != UINTIDX_NIL ) {
	    T_DIST lrt_xkinSubsetCjSumD = 
	      computeCostInstanceClusterJ
	      (xkinSubsetCj,
	       lcidx_Cj,
	       lpartlink_partition, 
	       aimattriag_dissimilarity
	       );
	    if ( lrt_xkinSubsetCjSumD < lrt_MjSumD ) {
	      lui_idxInstNewMj = xkinSubsetCj;
	      lrt_MjSumD = lrt_xkinSubsetCjSumD;
	    }
	  }
	}

#ifdef __VERBOSE_YES
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << "lrt_MjSumD: " << lrt_MjSumD
		    << " aioarrayui_idxMedoids[lcidx_Cj]: " << aioarrayui_idxMedoids[lcidx_Cj]  
	            << " lui_idxInstNewMj: " << lui_idxInstNewMj  
		    << '\n';
	}
	--geiinparam_verbose;
#endif /*__VERBOSE_YES*/

	/*It replaces the new Mj found with lower cost, because of the previous
	 */
	if ( lui_idxInstNewMj != aioarrayui_idxMedoids[lcidx_Cj]  ) {
	  aioarrayui_idxMedoids[lcidx_Cj] =  lui_idxInstNewMj; 
	  lb_medoidChange = true;
	}

      } //For Cj

  } while (lb_medoidChange);

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n";
    std::ostringstream lostrstream_labelMedoids;
    lostrstream_labelMedoids << "<MEDOIDS:"
			     << lpc_labelFunc;
    inout::containerprint
      (aioarrayui_idxMedoids,
       aioarrayui_idxMedoids + (T_INSTANCEIDX)aicidx_numKMedoids,
       std::cout,
       lostrstream_labelMedoids.str().c_str(),
       ','
       ); 
    std::cout << '\n';
    std::cout <<  "lor_cost = " << lor_cost;
    std::cout << std::endl;

    T_DIST lt_objetiveFunc = 
      um::SSEMedoid
      (aioarrayui_idxMedoids, 
       aicidx_numKMedoids,
       aimattriag_dissimilarity
       );
    std::cout << "lt_objetiveFunc: " << lt_objetiveFunc << std::endl;

  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lor_cost;
}

} /*END namespace medoids*/

#endif /*__CLUSTERING_OPERATOR_MEDOIDS_HPP__*/
