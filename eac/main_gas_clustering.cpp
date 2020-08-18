/*! \file main_gas_clustering.cpp
 *
 *   MAIN GET METRIC SATURDAY 30 JUNE
 * \brief Main program of Evolutionary Algorithms
 *
 * \details This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#include <iostream>

/*Fixed-K --- Encode label
 */
#ifdef  ALG_GA_CLUSTERING_LABELBASED_MURTHY_AND_CHOWDHURY_1996
#include "datatype_instance_real.hpp"
#include "gaclustering_fklabel.hpp"
#endif /*ALG_GA_CLUSTERING_LABELBASED_MURTHY_AND_CHOWDHURY_1996*/

#ifdef ALG_GKA_FKLABEL_KRISHNA_AND_MURTY_1999
#include "datatype_instance_real.hpp"
#include "gka_fklabel.hpp"
#endif /*ALG_GKA_FKLABEL_KRISHNA_AND_MURTY_1999*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004
#include "datatype_instance_real.hpp"
#include "igka_fklabel.hpp"
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/

#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004
#include "datatype_instance_real.hpp"
#include "igka_fklabel.hpp"
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/


/*Fixed-K --- Encode crisp matrix
 */
#ifdef ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994
#include "datatype_instance_real.hpp"
#include "gaclustering_fkcrispmatrix.hpp"
#endif /*ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994*/


/*Fixed-K ---  Encode centroids:
 */
#ifdef  ALG_CBGA_FKCENTROID_FRANTI_ETAL_1997
#ifdef  __INTEGER_INSTANCES__
#include "datatype_instance_integer.hpp"
#else
#include "datatype_instance_real.hpp"
#endif /*__INTEGER_INSTANCES__*/
#include "cbga_fkcentroid.hpp"
#endif /*ALG_CBGA_FKCENTROID_FRANTI_ETAL_1997*/

#ifdef  ALG_GAS_FKCENTROID_MAULIK_BANDYOPADHYAY_2000
#include "datatype_instance_real.hpp"
#include "gas_fkcentroid.hpp"
#endif /*ALG_GAS_FKCENTROID_MAULIK_BANDYOPADHYAY_2000*/

#ifdef  ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002
#include "datatype_instance_real.hpp"
#include "kga_fkcentroid.hpp"
#endif /*ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002*/

#ifdef  ALG_GAGR_FKCENTROID_CHANG_ETAL_2009
#include "datatype_instance_real.hpp"
#include "gagr_fkcentroid.hpp"
#endif /*ALG_GAGR_FKCENTROID_CHANG_ETAL_2009*/


/*Fixed-K --- Encode medoid
 */
#ifdef  ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993
#include "datatype_instance_real.hpp"
#include "clustering_operator_centroids.hpp"
#include "clustering_operator_medoids.hpp"
#include "gca_fkmedoid.hpp"
#endif /*ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993*/

#ifdef ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997
#include "datatype_instance_real.hpp"
#include "gaprototypes_fkmedoid.hpp"
#endif /*ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997*/

#ifdef  ALG_HKA_FKMEDOID_SHENG_LUI2004
#include "datatype_instance_real.hpp"
#include "clustering_operator_centroids.hpp"
#include "clustering_operator_medoids.hpp"
#include "hka_fkmedoid.hpp"
#endif /*ALG_HKA_FKMEDOID_SHENG_LUI2004*/


/*Variable-K --- Encode label
 */
#ifdef ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003
#include "datatype_instance_real.hpp"
#include "cga_vklabel.hpp"
#endif /*ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003*/

#if defined(ALG_EAC_VKLABEL_HRUSCHKA_CAMPELLO_CASTRO_2006) ||	\
  defined(ALG_EACI_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
  defined(ALG_EACII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
  defined(ALG_EACIII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
  defined(ALG_FEAC_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006)
#include "datatype_instance_real.hpp"
#include "feac_vklabel.hpp"
#endif /*ALG_EAC_VKLABEL_HRUSCHKA_CAMPELLO_CASTRO_2006) ||	\
	 ALG_EACI_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
         ALG_EACII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
	 ALG_EACIII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
	 ALG_FEAC_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006)
       */

#if defined(ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012) ||	\
  defined(ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012)
#include "datatype_instance_real.hpp"
#include "gga_vklabel.hpp"
#endif /*ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012 ||	\
	 ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
       */

/* Variable-K --- Encode centroids 
 */
#ifdef ALG_VGA_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2001
#include "datatype_instance_real.hpp"
#include "vga_vkcentroid.hpp"
#endif /*ALG_VGA_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2001*/

#ifdef ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002
#include "datatype_instance_real.hpp"
#include "gcuk_vkcentroid.hpp"
#endif /*ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002*/

#ifdef ALG_TGCA_VKCENTROID_HE_AND_TAN_2012
#include "datatype_instance_real.hpp"
#include "tgca_vkcentroid.hpp"
#endif /*ALG_TGCA_VKCENTROID_HE_AND_TAN_2012*/


/*Variable-K --- Encode other
 */
#ifdef ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001
#include "datatype_instance_real.hpp"
#include "clustering_vksubclusterbinary.hpp"
#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/


#ifdef ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003
#include "datatype_instance_real.hpp"
#include "gaclustering_vktreebinary.hpp"
#endif /*ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003*/

#ifdef __ALG_GAMIX__
#include "datatype_instance_real.hpp"
#include "ga_mix.hpp"
#endif /*__ALG_GAMIX__*/

#ifdef __ALG_GAMIX_SORT__
#include "datatype_instance_real.hpp"
#include "ga_mix_sort.hpp"
#define DATATYPE_INSTANCEIDX uintidx
#endif /*__ALG_GAMIX_SORT__*/

#ifdef __ALG_GADUAL_SORTFITNESSGENES__
#include "datatype_instance_real.hpp"
#include "gadual_sortfitnessgenes.hpp"
#define DATATYPE_INSTANCEIDX uintidx
#endif /*__ALG_GADUAL_SORTFITNESSGENES__*/

#ifdef __ALG_GAMIX_SIMPLE__
#include "datatype_instance_real.hpp"
#include "gamixsimple_fkcentroid.hpp"
#define DATATYPE_INSTANCEIDX uintidx
#endif /*__ALG_GAMIX_SIMPLE__*/

#ifdef __ALG_GAMIXGAGR__
#include "datatype_instance_real.hpp"
#include "gamixgagr_fkcentroid.hpp"
//#define DATATYPE_INSTANCEIDX uintidx
#endif /*__ALG_GAMIXGAGR__*/

#ifdef __ALG_GAMIXGAGR_CFARTHEST__
#include "datatype_instance_real.hpp"
#include "gamixgagrcfarthest_fkcentroid.hpp"
//#define DATATYPE_INSTANCEIDX uintidx
#endif /*__ALG_GAMIXGAGR_CFARTHEST__*/


#ifdef __ALG_GAMIXSQRTN_SIMPLE__
#include "datatype_instance_real.hpp"
#include "gamixsimplesqrtn_fkcentroid.hpp"
#define DATATYPE_INSTANCEIDX uintidx
#endif /*__ALG_GAMIXSQRTN_SIMPLE__*/

#include <leac.hpp>
#include "inparamclustering_getparameter.hpp"
#include "instances_read.hpp"
#include "bar_progress.hpp"

 
/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) 
{

#ifdef __VERBOSE_YES
  const char* lpc_labeMain = "main_gas_clustering";
  geverbosepc_labelstep = lpc_labeMain;
#endif /*__VERBOSE_YES*/
  

#ifdef __ALG_GAMIX__

  inout::InParamPcPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GAMIX",
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(1000);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.001);

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GAMIX__*/


#ifdef __ALG_GAMIX_SORT__

  inout::InParamGACDual
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GAMIX_SORT",
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );

  linparam_ClusteringGA.setNumMaxGenerations(2000);
  linparam_ClusteringGA.setSizePopulation(200);
  linparam_ClusteringGA.setSizeMatingPool(100);
  linparam_ClusteringGA.setAvgSimilarityMinimumPopulation(0.90);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.15);
  linparam_ClusteringGA.setSizePopulationSecondGenetics(20); 

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GAMIX_SORT__*/

#ifdef __ALG_GAMIX_SIMPLE__

  inout::InParamGACDual
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GAMIX_SIMPLE",
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );

  linparam_ClusteringGA.setNumMaxGenerations(2000);
  linparam_ClusteringGA.setSizePopulation(200);
  linparam_ClusteringGA.setSizeMatingPool(100);
  linparam_ClusteringGA.setAvgSimilarityMinimumPopulation(0.90);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.15);
  linparam_ClusteringGA.setSizePopulationSecondGenetics(20); 

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GAMIX_SIMPLE__*/


#ifdef __ALG_GAMIXGAGR__

  inout::InParamGACDual
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GAMIX_GAGR",
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );

  linparam_ClusteringGA.setNumMaxGenerations(1000);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setSizeMatingPool(50);
  linparam_ClusteringGA.setAvgSimilarityMinimumPopulation(-1.0);
  linparam_ClusteringGA.setProbCrossover(-1.0); 
  linparam_ClusteringGA.setProbMutation(0.7);
  linparam_ClusteringGA.setSizePopulationSecondGenetics(-1); 

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GAMIXGAGR__*/


#ifdef __ALG_GAMIXGAGR_CFARTHEST__

  inout::InParamGACDual
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    (
#if defined(__CROSSSAMPLE_TWO_RANDOM__)
     "GAMIX_GAGR_CROSSTWORAND",
#elif defined(__CROSSS_CIRCULAR__)
     "GAMIX_GAGR_CROSSCIRCULAR",
#else //ONE RANDON AND ONE FAR
     "GAMIX_GAGR_CFARTHEST",
#endif
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );

  linparam_ClusteringGA.setNumMaxGenerations(1000);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setSizeMatingPool(50);
  linparam_ClusteringGA.setAvgSimilarityMinimumPopulation(-1.0);
#if defined(__CROSSS_CIRCULAR__)
  linparam_ClusteringGA.setProbCrossover(0.9); 
#else 
  linparam_ClusteringGA.setProbCrossover(-1.0); 
#endif
  linparam_ClusteringGA.setProbMutation(0.7);
  linparam_ClusteringGA.setSizePopulationSecondGenetics(-1); 

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GAMIXGAGR_CFARTHEST__*/


#ifdef __ALG_GAMIXSQRTN_SIMPLE__

  inout::InParamGACDual
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GAMIXSQRTN_SIMPLE",
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );

  linparam_ClusteringGA.setNumMaxGenerations(2000);
  linparam_ClusteringGA.setSizePopulation(200);
  linparam_ClusteringGA.setSizeMatingPool(100);
  linparam_ClusteringGA.setAvgSimilarityMinimumPopulation(0.90);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.15);
  linparam_ClusteringGA.setSizePopulationSecondGenetics(20); 

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GAMIXSQRTN_SIMPLE__*/

#ifdef __ALG_GADUAL_SORTFITNESSGENES__

  inout::InParamGACDual
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GADUAL_FITNESSGENES",
     "2018",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );

  linparam_ClusteringGA.setNumMaxGenerations(2000);
  linparam_ClusteringGA.setSizePopulation(200);
  linparam_ClusteringGA.setSizeMatingPool(100);
  linparam_ClusteringGA.setAvgSimilarityMinimumPopulation(0.90);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.15);
  linparam_ClusteringGA.setSizePopulationSecondGenetics(20); 

    /*OUT:
   */
  inout::OutParamGACDual
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);

#endif /*__ALG_GADUAL_SORTFITNESSGENES__*/


#ifdef ALG_GA_CLUSTERING_LABELBASED_MURTHY_AND_CHOWDHURY_1996   
  /*INPUT: PARAMETER
   */
  inout::InParamPcPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GA",
     "Murthy and Chowdhury 1996", 
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(10000);
  linparam_ClusteringGA.setSizePopulation(6);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.5);

  /*OUT:
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);
  
#endif /*ALG_GA_CLUSTERING_LABELBASED_MURTHY_AND_CHOWDHURY_1996*/


#ifdef  ALG_CBGA_FKCENTROID_FRANTI_ETAL_1997
  /*INPUT: PARAMETER
   */
  inout::InParamCBGA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K,
     DATATYPE_INSTANCE_FREQUENCY
     > 
    linparam_ClusteringGA
    ("CBGA","Pasi Fränti1; Juha Kivijärvi2; Timo Kaukoranta2 and Olli Nevalainen 1997", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ,
     INPARAMCLUSTERING_CBGA_SELECMETH_ELITIST1,
     (int) 0 /*NUMBER OF GLA ITERATIONS*/
     ); 
  /*NUMBER MAXIMUM NUMBER OF GENERATION*/ 
  linparam_ClusteringGA.setNumMaxGenerations(256); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(8);      
  linparam_ClusteringGA.setProbMutation(0.01);

  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::Distortion);
    
#endif /*ALG_CBGA_FKCENTROID_FRANTI_ETAL_1997*/


#ifdef ALG_GKA_FKLABEL_KRISHNA_AND_MURTY_1999
  /*INPUT: PARAMETER
   */
  inout::InParamPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GKA",
     "Krishna and Murty 1999", 
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );
  linparam_ClusteringGA.setNumMaxGenerations(200);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setProbMutation(0.05);

  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSE);

#endif /*ALG_GKA_FKLABEL_KRISHNA_AND_MURTY_1999*/


#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004 
  /*INPUT: PARAMETER
   */
  inout::InParamPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K> 
    linparam_ClusteringGA
    ("IGKA",
     "Lu etal 2004", 
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );
  linparam_ClusteringGA.setNumMaxGenerations(200);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setProbMutation(0.05);

  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSE);
#endif /*ALG_IGKA_FKLABEL_LU_ETAL2004*/

#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004 
  /*INPUT: PARAMETER
   */
  inout::InParamPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("FGKA",
     "Lu etal 2004", 
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );
  linparam_ClusteringGA.setNumMaxGenerations(200);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setProbMutation(0.05);

  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSE);
#endif /*ALG_FGKA_FKLABEL_LU_ETAL2004*/


#ifdef  ALG_GAS_FKCENTROID_MAULIK_BANDYOPADHYAY_2000
  /*INPUT: PARAMETER
   */
  inout::InParamPcPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GAs",
     "Maulik and Bandyopadhyay 2000",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(100);
  linparam_ClusteringGA.setSizePopulation(100);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.001);
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);
    
#endif /*ALG_GAS_FKCENTROID_MAULIK_BANDYOPADHYAY_2000*/

#ifdef  ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002
  /*INPUT: PARAMETER
   */
  inout::InParamPcPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("KGA",
     "Bandyopadhyay and Maulik 2002", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  
  linparam_ClusteringGA.setNumMaxGenerations(1000);
  linparam_ClusteringGA.setSizePopulation(50);
  linparam_ClusteringGA.setProbCrossover(0.8);
  linparam_ClusteringGA.setProbMutation(0.001);
 
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);
    
#endif /*ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002*/

#ifdef ALG_VGA_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2001
  /*INPUT: PARAMETER
   */
  inout::InParamPcPmVk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("VGA",
     "Bandyopadhyay and Maulik 2001", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(1000);
  linparam_ClusteringGA.setSizePopulation(50);  
  linparam_ClusteringGA.setProbCrossover(0.8);  
  linparam_ClusteringGA.setProbMutation(0.05); 
  /*The values Kmin and Kmax are token to 2 and 10, respectively
   */ 
  linparam_ClusteringGA.setNumClusterKMinimum(2);   
  linparam_ClusteringGA.setNumClusterKMaximum(INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED);    

  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamGAC(inout::IndexI);
#endif /*ALG_VGA_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2001*/
  

#ifdef ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002
  /*INPUT: PARAMETER
   */
  inout::InParamPcPmVk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GCUK",
     "Bandyopadhyay and Maulik 2002", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(100);
  linparam_ClusteringGA.setSizePopulation(50);  
  linparam_ClusteringGA.setProbCrossover(0.8);  
  linparam_ClusteringGA.setProbMutation(0.001); 
  /*The values Kmin and Kmax are token to 2 and 10, respectively
   */ 
  linparam_ClusteringGA.setNumClusterKMinimum(2);   
  linparam_ClusteringGA.setNumClusterKMaximum(10);    

  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamGAC(inout::DBindex);
#endif /*ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002*/
  

#ifdef ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003
  /*INPUT: PARAMETER
   */
#ifdef _INITIATES_KMIN_KMAX_POPULATION_
  inout::InParamGenWOChgVk
    <DATATYPE_BITSIZE,
     DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GA_CASILLAS2003","Casillas and Gonzalez and Martinez 2003",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );  
#endif // _INITIATES_KMIN_KMAX_POPULATION_

#ifdef _INITIATES_RANDOM_POPULATION_
  inout::InParamGAClusteringVKTreeBinary
    <DATATYPE_REAL,DATATYPE_CLUSTERIDX> 
    linparam_ClusteringGA
    ("GA_C","Casillas and Gonzalez and Martinez 2003",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );

  linparam_ClusteringGA.setNumClusterKMinimum(INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED);
  linparam_ClusteringGA.setNumClusterKMaximum(INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED);
  
#endif //_INITIATES_RANDOM_POPULATION_

  /*NUMBER MAXIMUM NUMBER OF ITERATION*/ 
  linparam_ClusteringGA.setNumMaxGenerations(0); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(0);
  /*PROBABILITY CROSSOVER*/
  linparam_ClusteringGA.setProbCrossover(0.8);
  /*PROBABILITY MUTATION*/  
  linparam_ClusteringGA.setProbMutation(0.008);
  linparam_ClusteringGA.setNumNotChangeStop(3);     
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamGAC(inout::VRC);
 
#endif /*ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003*/


#ifdef ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001
  /*INPUT: PARAMETER
   */
  inout::InParamSubClusterBinaryVk
    <DATATYPE_REAL,
     DATATYPE_BITSIZE,
     DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K> 
    linparam_ClusteringGA
    (
#if defined _INITIALIZED_ONLY_ONCE
     "CLUSTERINGPLUS",
#else
     "CLUSTERING",
#endif
     "Lin Yu Tseng and Shiueng Bien Yang", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setU(1.4);
  linparam_ClusteringGA.setLambda(0.125);
  linparam_ClusteringGA.setW1(1.0);   
  linparam_ClusteringGA.setW2(3.0);
  /*NUMBER MAXIMUM NUMBER OF ITERATION*/ 
  linparam_ClusteringGA.setNumMaxGenerations(100);
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(50);       
  /*PROBABILITY CROSSOVER*/
  linparam_ClusteringGA.setProbCrossover(0.8); 
  /*PROBABILITY MUTATION*/
  linparam_ClusteringGA.setProbMutation(0.05); 
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamGAC(inout::IntraInterClust);
  

#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/

#ifdef ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003   
  /*INPUT: PARAMETER*/
  inout::InParamWithoutPcPmVk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("CGA",
     "Hruschka and Ebecken 2003", 
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumClusterKMaximum
    (INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED);
  
  linparam_ClusteringGA.setNumMaxGenerations(200);
  
  /*Populations formed by 20 genotypes
   */
  linparam_ClusteringGA.setSizePopulation(20);
  
  /*OUT:
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::Silhouette); 
#endif /*ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003*/
  
#ifdef ALG_EAC_VKLABEL_HRUSCHKA_CAMPELLO_CASTRO_2006
  /*INPUT: PARAMETER
   */
  inout::InParamFEAC
    <DATATYPE_FEATURE,
     DATATYPE_REAL,
     DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("EAC",
     "Hruschka Eduardo R.; Campello; Ricardo J. G. B. and de Castro Leandro N.",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     ); 
  /*NUMBER MAXIMUM NUMBER OF GENERATION*/
  linparam_ClusteringGA.setNumMaxGenerations(500); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(20);    
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSilhouette);
#endif /*ALG_EAC_VKLABEL_HRUSCHKA_CAMPELLO_CASTRO_2006*/


#ifdef ALG_EACI_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
  /*INPUT: PARAMETER
   */
  inout::InParamFEAC
    <DATATYPE_FEATURE,
     DATATYPE_REAL,
     DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K>
    linparam_ClusteringGA
    ("EAC-I",
     "Vinicius S. Alves; Ricardo J. G. B. Campello and Eduardo R. Hruschka",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     ); 

  /*NUMBER MAXIMUM NUMBER OF GENERATION*/
  linparam_ClusteringGA.setNumMaxGenerations(500); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(20);    
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSilhouette);
#endif /*ALG_EACI_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006*/


#ifdef ALG_EACII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
  /*INPUT: PARAMETER
   */
  inout::InParamFEAC
    <DATATYPE_FEATURE,
     DATATYPE_REAL,
     DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K>
    linparam_ClusteringGA
    ("EAC-II",
     "Vinicius S. Alves; Ricardo J. G. B. Campello and Eduardo R. Hruschka",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     ); 

  /*NUMBER MAXIMUM NUMBER OF GENERATION*/
  linparam_ClusteringGA.setNumMaxGenerations(500); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(20);    
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSilhouette);
#endif /*ALG_EACII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006*/


#ifdef ALG_EACIII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
  /*INPUT: PARAMETER
   */
  inout::InParamFEAC
    <DATATYPE_FEATURE,
     DATATYPE_REAL,
     DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K>
    linparam_ClusteringGA
    ("EAC-III",
     "Vinicius S. Alves; Ricardo J. G. B. Campello and Eduardo R. Hruschka",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     ); 

  /*NUMBER MAXIMUM NUMBER OF GENERATION*/
  linparam_ClusteringGA.setNumMaxGenerations(500); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(20);    
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSilhouette);
#endif /*ALG_EACIII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006*/


#ifdef  ALG_FEAC_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
  /*INPUT: PARAMETER
   */
  inout::InParamFEAC
    <DATATYPE_FEATURE,
     DATATYPE_REAL,
     DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    (
#ifdef __FITNESS_SIMPLIFIED_SILHOUETTE__
     "FEAC_SSILHOUETTE",
#endif /*__FITNESS_SIMPLIFIED_SILHOUETTE__*/
#ifdef __FITNESS_RAND_INDEX__
     "FEAC_RANDINDEX",
#endif /*__FITNESS_RAND_INDEX__*/
     "Vinicius S. Alves; Ricardo J. G. B. Campello and Eduardo R. Hruschka",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     ); 
  /*NUMBER MAXIMUM NUMBER OF GENERATION*/ 
  linparam_ClusteringGA.setNumMaxGenerations(500);
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(20);    

  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SSilhouette);
#endif /*ALG_FEAC_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006*/

#ifdef ALG_TGCA_VKCENTROID_HE_AND_TAN_2012
  /*INPUT: PARAMETER
   */
#ifdef __INITIALIZATION_RANDOM_SAMPLING__
  inout::InParamTGCA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL
     > 
    linparam_ClusteringGA
    ("TGCA_INITRANDOM","He Hong and Tan Yonghong 2012",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );
#elif __DELETE_EMPTY_CLUSTER__
  inout::InParamTGCA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL
     > 
    linparam_ClusteringGA
    ("TGCA_DELEMPTYCLUSTER","He Hong and Tan Yonghong 2012",
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );
#else /*PROPOSED BY THE AUTHOR*/
  inout::InParamTGCA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("TGCA","He Hong and Tan Yonghong 2012", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     );
#endif /*__INITIALIZATION RANDOM SAMPLING__*/

  /*NUMBER MAXIMUM NUMBER OF ITERATION*/ 
  linparam_ClusteringGA.setNumMaxGenerations(200);
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(60);      
  linparam_ClusteringGA.setNumSubpopulationsCross(4);
  /*PROBABILITY CROSSOVER*/
  linparam_ClusteringGA.setProbCrossover(0.8);      
  /*PROBABILITY MUTATION:
    the mutation probabilty pm varies with ki and the variation
    tendency of pm can be also divided into two stages.
  */
  linparam_ClusteringGA.setProbMutation             
    (UNKNOWN_CLUSTER_REAL);      
  /*if k is unknown, its value can be
    in the set K, where K=[kmin,kmax], 
    kmin=2 and kmax = roundo(sqrt(M)).
  */
 
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX
     >
    loop_outParamGAC(inout::VRC);
#endif /*ALG_TGCA_VKCENTROID_HE_AND_TAN_2012*/

#ifdef  ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012
  /*INPUT: PARAMETER
   */
  inout::InParamGGA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GGA_DBINDEX",
     "Agustin-Blas L.E. and Salcedo-Sanz S. and Jimenez-Fernandez S. and "
     "Carro-Calvo L. and Del Ser J. and Portilla-Figueras J.A.",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::DBindex);
  
#endif /*ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012*/

#ifdef ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
  /*INPUT: PARAMETER
   */
  inout::InParamGGA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_REAL, 
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GGA_SILHOUETTE",
     "Agustin-Blas L.E. and Salcedo-Sanz S. and Jimenez-Fernandez S. and "
     "Carro-Calvo L. and Del Ser J. and Portilla-Figueras, J.A.",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX
     >
    loop_outParamGAC(inout::Silhouette);
#endif /*ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012*/
  
#ifdef  ALG_GAGR_FKCENTROID_CHANG_ETAL_2009
  /*INPUT: PARAMETER
   */
  inout::InParamAdaptivePcPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GAGR",
     "Dong-Xia Chang and Xian-Da Zhang and Chang-Wen Zheng 2009", 
     inout::CENTROIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ
     ); 
  /*NUMBER MAXIMUM NUMBER OF GENERATION*/ 
  linparam_ClusteringGA.setNumMaxGenerations(50); 
  /*SIZE POPULATION*/
  linparam_ClusteringGA.setSizePopulation(50);    
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>   
    loop_outParamGAC(inout::SSE);    
#endif /*ALG_GAGR_FKCENTROID_CHANG_ETAL_2009*/


#ifdef ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993 
  /*INPUT: PARAMETER
   */
  inout::InParamGCA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GCA",
     "C.B. Lucasius and A.D. Dane and G. Kateman",
     inout::MEDOIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(200);
  linparam_ClusteringGA.setSizePopulation(200);
  /*Mix recombination probability x7*/
  linparam_ClusteringGA.setProbCrossover(0.9);
  /*Point mutation probability    x8*/ 
  linparam_ClusteringGA.setProbMutation(0.4);
  /*Mix mutation probability   x9*/
  linparam_ClusteringGA.setProbMixMutation(0.125);

  /*OUT
   */
  inout::OutParamGAMedoid
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);
  
#endif /*ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993*/

#ifdef ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997
   /*INPUT: PARAMETER
   */
  inout::InParamGAPrototypesFk
    <DATATYPE_BITSIZE,
     DATATYPE_CLUSTERIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("GAPROTOTYPES",
     "L. I. Kuncheva and J. C. Bezdek",
     inout::MEDOIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_INDUCED
     );
  linparam_ClusteringGA.setNumMaxGenerations(500);
  linparam_ClusteringGA.setSizePopulation(20);
  linparam_ClusteringGA.setProbCrossover(0.5);   
  linparam_ClusteringGA.setProbMutation(0.015);     
  linparam_ClusteringGA.setPini(-1.0);
  linparam_ClusteringGA.setAlpha(10.0);
  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::J1);
#endif /*ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997*/
  

#ifdef ALG_HKA_FKMEDOID_SHENG_LUI2004
  /*INPUT: PARAMETER
   */
  inout::InParamHKA
    <DATATYPE_CLUSTERIDX,
     DATATYPE_INSTANCEIDX,
     DATATYPE_REAL,
     DATATYPE_FEATURE,         
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     > 
    linparam_ClusteringGA
    ("HKA",
     "Weiguo Sheng and Xiaohui Liu",
     inout::MEDOIDS,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  linparam_ClusteringGA.setNumMaxGenerations(200);
  linparam_ClusteringGA.setSizePopulation(200);
  /*Mix recombination probability x7*/
  linparam_ClusteringGA.setProbCrossover(0.95);
  /*Point mutation probability    x8*/ 
  linparam_ClusteringGA.setProbMutation(0.02);
  /*Mix mutation probability   x9*/
  linparam_ClusteringGA.setProbMixMutation(0.05);
  linparam_ClusteringGA.setOrderTournament((size_t) 2);    
  linparam_ClusteringGA.setNearestNeighbors((size_t) 3);     
  linparam_ClusteringGA.setProbSearchHeuristic(0.2);

  /*OUT
   */
  inout::OutParamGAMedoid
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamGAC(inout::SED);
#endif /*ALG_HKA_FKMEDOID_SHENG_LUI2004*/

#ifdef ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994
  /*\cite{Sheng:Xiaohui:GAclusteringMedoid:HKA:2004}
   */
  inout::InParamWithoutPcPmFk
    <DATATYPE_CLUSTERIDX,
     DATATYPE_BITSIZE,
     DATATYPE_FEATURE,
     DATATYPE_FEATURE_SUM,
     DATATYPE_INSTANCES_CLUSTER_K
     >
    linparam_ClusteringGA
    ("GA_B",
     "James C. Bezdek and Srinivas Boggavarapu and Lawrence O. Hall and Amine Bensaid",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_INDUCED
     );

  linparam_ClusteringGA.setNumMaxGenerations(2000);
  linparam_ClusteringGA.setSizePopulation(200);
  linparam_ClusteringGA.setSizeMatingPool(100);

  /*OUT
   */
  inout::OutParamGAC
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX> 
    loop_outParamGAC(inout::J1);
#endif /*ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994*/
  
  inout::OutParamClusteringMetric<DATATYPE_REAL> loop_outParamlusteringMetric;

  /*READ PARAMETER*/
  inparamclustering_getParameter(linparam_ClusteringGA, argc, argv);

#ifdef __VERBOSE_YES
  std::cout << std::boolalpha;
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labeMain
	      << ":  IN(" << geiinparam_verbose << ")\n"
	      << "argv: ";
    std::cout << argv[0] << ' '; 
    for (int li_i = 1; li_i < argc;  li_i++ ) {
      std::cout << ' '  << argv[li_i];
    } 
    std::cout << std::endl;
  }
#endif /*__VERBOSE_YES*/ 
  
  /*BAR PRINTING*/
  if ( linparam_ClusteringGA.getProgressBarPrinting() ) 
    inout::barprogress_update
      (0, 
       linparam_ClusteringGA.getTimesRunAlgorithm() * 
       linparam_ClusteringGA.getNumFilesInstance()
       );

  for (uintidx __lst_i = 0; __lst_i < linparam_ClusteringGA.getNumFilesInstance(); __lst_i++) {
  
    linparam_ClusteringGA.setCurrentFileInstance(__lst_i);

    /*READ: INSTANCES 
     */
#ifdef __INSTANCES_WITH_FREQUENCY

    auto lpairvec_dataset =  inout::dataSetReadWithFreq(linparam_ClusteringGA);
 
#else /* INSTANCES WITHOUT FREQUENCY */

    auto lpairvec_dataset =  inout::dataSetRead(linparam_ClusteringGA);
 
#endif /*__INSTANCES_WITH_FREQUENCY*/

    std::pair<std::vector<std::string>,std::string> 
      lpairvecstrstr_instanceDimName = 
      inout::instancesReadDimName
      (linparam_ClusteringGA,
       false
       );
    
    linparam_ClusteringGA.setNumInstances
      ((uintidx) lpairvec_dataset.first.size());
    linparam_ClusteringGA.setNumDimensionsInstances
      (data::Instance<DATATYPE_FEATURE>::getNumDimensions());

    linparam_ClusteringGA.setNumInstancesTest
      ( (uintidx) lpairvec_dataset.second.size());

    /*STATISTICS OF INSTANCES
     */  
    DATATYPE_FEATURE *larray_meanFeactures =
      new DATATYPE_FEATURE[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];
  
    DATATYPE_FEATURE *larray_desvstdFeactures =
      new DATATYPE_FEATURE[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];

#ifdef __INSTANCES_WITH_FREQUENCY
    /*-----------BEGIN STATISTICS OF INSTANCES__INSTANCES_WITH_FREQUENCY
     */
    {
      DATATYPE_FEATURE_SUM* larray_sumFeatureTmp =
	new DATATYPE_FEATURE_SUM[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];

      const uintidx lui_numInstances = 
	stats::sumFeacturesFreq
	(larray_sumFeatureTmp,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_FEATURE(0),
	 [](const  data::Instance<DATATYPE_FEATURE>* linst_iter) ->
	 const DATATYPE_INSTANCE_FREQUENCY
	 {
	   data::InstanceFreq<DATATYPE_FEATURE,DATATYPE_INSTANCE_FREQUENCY>*
	     linstfreq_inst = 
	     (data::InstanceFreq<DATATYPE_FEATURE,DATATYPE_INSTANCE_FREQUENCY>*)
	     linst_iter;
	 
	   return linstfreq_inst->getFrequency();
	 }
	 );
  
      stats::meanVector
	(larray_meanFeactures,
	 lui_numInstances,
	 larray_sumFeatureTmp
	 );

      stats::sumFeacturesFreqSQ
	(larray_sumFeatureTmp,
	 larray_meanFeactures,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 [](const  data::Instance<DATATYPE_FEATURE>* linst_iter)
	 -> const DATATYPE_INSTANCE_FREQUENCY
	 {
	   data::InstanceFreq<DATATYPE_FEATURE,DATATYPE_INSTANCE_FREQUENCY>*
	     linstfreq_inst = 
	     (data::InstanceFreq<DATATYPE_FEATURE,DATATYPE_INSTANCE_FREQUENCY>*)
	     linst_iter;
	   return linstfreq_inst->getFrequency();
	 }
	 );

      stats::meanVector
	(larray_desvstdFeactures,
	 (uintidx) lui_numInstances-1,
	 larray_sumFeatureTmp
	 );

      stats::varTodesvstd(larray_desvstdFeactures);

      delete [] larray_sumFeatureTmp;
    
    }
#else /*BEGIN STATISTICS OF INSTANCES__INSTANCES_WITHOUT_FREQUENCY
       */    
    {
      DATATYPE_FEATURE_SUM* larray_sumFeatureTmp =
	new DATATYPE_FEATURE_SUM[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];

      stats::sumFeactures
	(larray_sumFeatureTmp,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_FEATURE(0)
	 );

      stats::meanVector
	(larray_meanFeactures,
	 (DATATYPE_INSTANCES_CLUSTER_K) lpairvec_dataset.first.size(),
	 larray_sumFeatureTmp
	 );
    
      stats::sumFeacturesSQ
	(larray_sumFeatureTmp,
	 larray_meanFeactures,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );

      stats::meanVector
	(larray_desvstdFeactures,
	 (uintidx) lpairvec_dataset.first.size()-1,
	 larray_sumFeatureTmp
	 );

      stats::varTodesvstd(larray_desvstdFeactures);
    
      delete [] larray_sumFeatureTmp;
      
    }

#endif /*END STATISTICS OF INSTANCES__INSTANCES_WITHOUT_FREQUENCY
	*/

    /*Declaration of a reference to a generic object 
       of type of dist::Dist.
       DATATYPE_REAL is the type of data obtained when 
         calculating the distance.
       DATATYPE_FEATURE is the data type of the dimensions 
         of the instances or objects
     */
    dist::Dist<DATATYPE_REAL,DATATYPE_FEATURE>
      *pfunct2p_distAlg = NULL;
    dist::Dist<DATATYPE_REAL,DATATYPE_FEATURE>
      *pfunct2p_distEuclidean =
      new dist::Euclidean<DATATYPE_REAL,DATATYPE_FEATURE>();
    dist::Dist<DATATYPE_REAL,DATATYPE_FEATURE>
      *pfunct2p_distEuclideanSq =
      new dist::EuclideanSquared<DATATYPE_REAL,DATATYPE_FEATURE>();
    
    /*Create a dist::Dist object based on the parameter provided 
         by the user, which will have a polymorphic behavior
     */
    switch ( linparam_ClusteringGA.getOpDistance() ) {
  
    case INPARAMCLUSTERING_DISTANCE_EUCLIDEAN:
      pfunct2p_distAlg = 
	new dist::Euclidean<DATATYPE_REAL,DATATYPE_FEATURE>();
      break;
    
    case INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ:
      pfunct2p_distAlg = 
	new dist::EuclideanSquared<DATATYPE_REAL,DATATYPE_FEATURE>();
      break;

#if  DATATYPE_CENTROIDS_ROUND == 0
    
    case INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_INDUCED:
      pfunct2p_distAlg = 
	new dist::Induced<DATATYPE_REAL,DATATYPE_FEATURE>
	( mat::getIdentity
	  <DATATYPE_REAL>
	  (data::Instance<DATATYPE_FEATURE>::getNumDimensions()));
      break;
      
    case INPARAMCLUSTERING_DISTANCE_DIAGONAL_INDUCED:
      pfunct2p_distAlg = 
	new dist::Induced<DATATYPE_REAL,DATATYPE_FEATURE>
	( stats::getMatrixDiagonal<DATATYPE_FEATURE>
	  (larray_desvstdFeactures));
      break;
      
    case  INPARAMCLUSTERING_DISTANCE_MAHALONOBIS_INDUCED:      
      pfunct2p_distAlg = 
	new dist::Induced<DATATYPE_REAL,DATATYPE_FEATURE>
	( stats::getMatrixMahalonobis
	  (lpairvec_dataset.first.begin(),
	   lpairvec_dataset.first.end()
	   )
	  );
      break;
      
#endif //DATATYPE_CENTROIDS_ROUND
      
    default:
      throw  std::invalid_argument("main_gas_clustering: undefined norm");
      break;
    }
  
    for (int  li_l = 1; li_l <= linparam_ClusteringGA.getTimesRunAlgorithm(); li_l++) {

      if ( linparam_ClusteringGA.getRandomSeed().size() == 0 ||
	   linparam_ClusteringGA.getTimesRunAlgorithm() != 1 || 
	   linparam_ClusteringGA.getNumFilesInstance() != 1 )  {
	
	std::string lstr_seed_seq =
	  randomext::setSeed();
	linparam_ClusteringGA.setRandomSeed( lstr_seed_seq );
      }
      else {

	randomext::setSeed(linparam_ClusteringGA.getRandomSeed());
       
      }

      loop_outParamGAC.initialize(li_l);
      loop_outParamlusteringMetric.initialize(li_l);
    

      if ( !linparam_ClusteringGA.getPrintMulLine() ) {
       /*BEGIN PRINT PARAMETERS
	 */
	inout::OutFileName   lofn_filename;
	std::ostream& lostream_out = 
	  lofn_filename.openFile(linparam_ClusteringGA.getFileNameTimesRun());
	lostream_out << "_inout" <<  inout::OutFileName::getDelim() << "in" << inout::OutFileName::getDelim();
	linparam_ClusteringGA.print(lostream_out,inout::OutFileName::getDelim());
	lostream_out << std::endl;
	lofn_filename.closeFile();
      
      } /*END PRINT PARAMETERS
	 */

#ifdef ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994

      gaencode::ChromosomeCrispMatrix
	<DATATYPE_BITSIZE,DATATYPE_CLUSTERIDX,DATATYPE_REAL>&&
	lchrom_best =  
	eac::gaclustering_fkcrispmatrix
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
        
#ifdef __VERBOSE_YES
      geverbosepc_labelstep = lpc_labeMain;
#endif /*__VERBOSE_YES*/ 

      /*DECODE CHROMOSOME*/
      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids /*CENTROIDS ASOCIATED TO U_i*/
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE>&& lmatrixT_y =
	data::toMatrixRow
	(lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>
	lmatrixT_sumWX
	(lomatrixrowt_centroids.getNumRows(),
	 lomatrixrowt_centroids.getNumColumns()
	 );
      
      std::vector<DATATYPE_INSTANCES_CLUSTER_K>
	lvectorT_sumWik(lomatrixrowt_centroids.getNumRows());
    
      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lmatrixT_sumWX,
	 lvectorT_sumWik,
	 lchrom_best, 
	 lmatrixT_y 
	 );
    
      partition::PartitionCrispMatrix
	<DATATYPE_BITSIZE,DATATYPE_CLUSTERIDX>
	lpartition_clusters(lchrom_best);

#endif /*ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994*/

#ifdef ALG_GA_CLUSTERING_LABELBASED_MURTHY_AND_CHOWDHURY_1996   

      gaencode::ChromFixedLength<DATATYPE_CLUSTERIDX,DATATYPE_REAL>&& lchrom_best = 
	eac::gaclustering_fklabel
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
    
      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 lchrom_best.getStringSize(),
	 linparam_ClusteringGA.getNumClusterK()
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 DATATYPE_INSTANCES_CLUSTER_K(0)
	 );
      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lomatrixrowt_sumInstCluster,
	 lovectort_numInstClusterK,
	 lpartition_clusters,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );
      
#endif /*ALG_GA_CLUSTERING_LABELBASED_MURTHY_AND_CHOWDHURY_1996*/

      /*\cite{Franti:etal:GAclustering:gafranti:1997}*/
#ifdef  ALG_CBGA_FKCENTROID_FRANTI_ETAL_1997
      
      auto lchrom_best = 
	eac::cbga_fkcentroid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );


      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 lchrom_best.getCodeBook().toArray()
	 );
      
      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getPartition().getMembersShip(),
	 (uintidx) lpairvec_dataset.first.size(),
	 linparam_ClusteringGA.getNumClusterK()
	 );
       
#endif /*ALG_CBGA_FKCENTROID_FRANTI_ETAL_1997*/


      /*\cite{Krishna:Murty:GAClustering:GKA:1999}
       */
#ifdef ALG_GKA_FKLABEL_KRISHNA_AND_MURTY_1999
    
      gaencode::ChromFixedLength<DATATYPE_CLUSTERIDX,DATATYPE_REAL>&& lchrom_best = 
	eac::gka_fklabel 
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
     
      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 lchrom_best.getStringSize(),
	 linparam_ClusteringGA.getNumClusterK()
	 );
    
      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 DATATYPE_INSTANCES_CLUSTER_K(0)
	 );

      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lomatrixrowt_sumInstCluster,
	 lovectort_numInstClusterK,
	 lpartition_clusters,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );
    
#endif /*ALG_GKA_FKLABEL_KRISHNA_AND_MURTY_1999*/

      /* \cite{Lu:etal:GAclusteringLabel:FGKA:2004}
	 \cite{Lu:etal:GAclusteringLabel:IGKA:2004}
      */
#ifdef ALG_FGKA_FKLABEL_LU_ETAL2004 
      
      gaencode::ChromosomeFGKA<DATATYPE_CLUSTERIDX,DATATYPE_REAL>&& lchrom_best = 
	eac::igka_fklabel
	<DATATYPE_CLUSTERIDX, //-1, 0, 1, .., K
	 DATATYPE_REAL,
	 DATATYPE_FEATURE,
	 DATATYPE_FEATURE_SUM,
	 DATATYPE_INSTANCES_CLUSTER_K>
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first,
	 *pfunct2p_distAlg
	 );

      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 lchrom_best.getStringSize(),
	 linparam_ClusteringGA.getNumClusterK()
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) linparam_ClusteringGA.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 DATATYPE_INSTANCES_CLUSTER_K(0)
	 );

      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lomatrixrowt_sumInstCluster,
	 lovectort_numInstClusterK,
	 lpartition_clusters,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );
#endif /* ALG_FGKA_FKLABEL_LU_ETAL2004 */

    
#ifdef ALG_IGKA_FKLABEL_LU_ETAL2004 
      
      gaencode::ChromosomeIGKA
	<DATATYPE_CLUSTERIDX,
	 DATATYPE_REAL,
	 DATATYPE_FEATURE,
	 DATATYPE_FEATURE_SUM,
	 DATATYPE_INSTANCES_CLUSTER_K>
	lchrom_best =
	eac::igka_fklabel
	<DATATYPE_CLUSTERIDX, //-1, 0, 1, .., K
	 DATATYPE_REAL,
	 DATATYPE_FEATURE,
	 DATATYPE_FEATURE_SUM,
	 DATATYPE_INSTANCES_CLUSTER_K>
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first,
	 *pfunct2p_distAlg
	 );
    
      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 lchrom_best.getStringSize(),
	 linparam_ClusteringGA.getNumClusterK()
	 );

      mat::MatrixRow<DATATYPE_FEATURE>& lomatrixrowt_centroids =
	lchrom_best.getCentroids();
      
#endif /* ALG_IGKA_FKLABEL_LU_ETAL2004 */


      /*\cite{Maulik:Bandyopadhyay:GAclustering:GAS:2000}
       */
#ifdef ALG_GAS_FKCENTROID_MAULIK_BANDYOPADHYAY_2000

      gaencode::ChromFixedLength<DATATYPE_FEATURE,DATATYPE_REAL>&& lchrom_best = 
	eac::gas_fkcentroid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
 	
      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 lchrom_best.getString()
	 );

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );

#endif /*ALG_GAS_FKCENTROID_MAULIK_BANDYOPADHYAY_2000*/


#if defined(__ALG_GAMIX__) ||	\
  defined(__ALG_GAMIX_SIMPLE__) || \
  defined(__ALG_GAMIXGAGR__) || \
  defined(__ALG_GAMIXGAGR_CFARTHEST__) || \
  defined(__ALG_GAMIXSQRTN_SIMPLE__) || \
  defined(__ALG_GAMIX_SORT__) || \
  defined(__ALG_GADUAL_SORTFITNESSGENES__)


      auto lchrom_best =
	eac::gamix
	(loop_outParamGAC,
	linparam_ClusteringGA,
	lpairvec_dataset.first.begin(),
	lpairvec_dataset.first.end(),
	*pfunct2p_distAlg
	);

      //std::cout << "LLEGUE AL MAIN" << std::endl;

      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids(lchrom_best.getCentroids());

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	lpairvec_dataset.first.begin(),
	lpairvec_dataset.first.end(),
	DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	*pfunct2p_distAlg
	);

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_medoids
	( (uintidx) linparam_ClusteringGA.getNumClusterK(),
	  data::Instance<DATATYPE_FEATURE>::getNumDimensions()
	  );

      clusteringop::initialize
	(lomatrixrowt_medoids,
	 lpairvec_dataset.first.begin(),
	 lchrom_best.begin()
	 );

      std::pair<DATATYPE_REAL,bool> lpair_sedmedoid =
	um::SSE
	(lomatrixrowt_medoids,
	  lpairvec_dataset.first.begin(),
	  lpairvec_dataset.first.end(),
	  *pfunct2p_distEuclidean
	  );

	 loop_outParamGAC.setMedoidSED(lpair_sedmedoid.first);

	 //std::cout << "VOY A LAS METRICAS" << std::endl;

#endif /*__ALG_GAMIX__*/
    
      /*\cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
       */
#ifdef ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002 
      gaencode::ChromFixedLength<DATATYPE_FEATURE,DATATYPE_REAL>&& lchrom_best = 
	eac::kga_fkcentroid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
 	
      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 lchrom_best.getString()
	 );

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );
    
#endif /*ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002*/


      /*\cite{Hruschka:Ebecken:GAClusteringLabelKVar:CGA:2003}
       */
#ifdef ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003
    
      gaencode::ChromFixedLength<DATATYPE_CLUSTERIDX,DATATYPE_REAL>&& lchrom_best = 
	eac::cga_vklabel
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
    
      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 lchrom_best.getStringSize(),
	 loop_outParamGAC.getNumClusterK()
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) loop_outParamGAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) loop_outParamGAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) loop_outParamGAC.getNumClusterK(),
	 DATATYPE_INSTANCES_CLUSTER_K(0)
	 );

      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lomatrixrowt_sumInstCluster,
	 lovectort_numInstClusterK,
	 lpartition_clusters,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );
#endif /*ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003*/

#if defined(ALG_EAC_VKLABEL_HRUSCHKA_CAMPELLO_CASTRO_2006) ||	\
  defined(ALG_EACI_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
  defined(ALG_EACII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
  defined(ALG_EACIII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006) ||	\
  defined(ALG_FEAC_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006)

#ifdef __FITNESS_RAND_INDEX__
      if ( linparam_ClusteringGA.getClassInstanceColumn() == 0 ) {
          throw  std::invalid_argument("main_gas_clustering: for the fitness rand index, define the class column");
      }
#endif /*__FITNESS_RAND_INDEX__*/      

      auto lchrom_best =
	eac::feca_vklabel
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
      
      mat::MatrixRow<DATATYPE_FEATURE>&&
	lomatrixrowt_centroids =
	lchrom_best.getCentroids().getMatrix();

      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 lchrom_best.getStringSize(),
	 loop_outParamGAC.getNumClusterK()
	 );

#endif /*ALG_EAC_VKLABEL_HRUSCHKA_CAMPELLO_CASTRO_2006) ||	\
	 ALG_EACI_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
	 ALG_EACII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
	 ALG_EACIII_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006
	 ALG_FEAC_VKLABEL_ALVES_CAMPELLO_HRUSCHKA_2006)
       */

      /*\cite{Tseng:Yang:GAclusteringVarK:CLUSTERING:2001}
       */
#ifdef ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001
      
      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids;
      std::vector<DATATYPE_CLUSTERIDX> lovectormmidx_memberShip;

      partition::PartitionDisjSets
	<DATATYPE_CLUSTERIDX>            lmembclassdisjsets_BkinCi;
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> lvectort_numInstCi;
   
      gaencode::ChromosomeBitArray<DATATYPE_BITSIZE,DATATYPE_REAL> lchrom_best;
      mat::MatrixRow<DATATYPE_FEATURE> lmatrixrowt_Vi;
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> lvectort_numInstBi;
      partition::PartitionDisjSets<DATATYPE_CLUSTERIDX> lmembclassdisjsets_Bi;
    
      std::tie(lchrom_best,lmatrixrowt_Vi,lvectort_numInstBi,lmembclassdisjsets_Bi) =
	eac::clustering_vksubclusterbinary
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );

      if ( loop_outParamGAC.getEndingCondition() ) {
      
	std::tie(lomatrixrowt_centroids,lmembclassdisjsets_BkinCi,lvectort_numInstCi) =
	  clusteringop::getClusters
	  (lchrom_best,
	   lmatrixrowt_Vi,
	   lvectort_numInstBi,
	   *pfunct2p_distAlg,
	   (DATATYPE_CLUSTERIDX) 0
	   );

	lovectormmidx_memberShip =
	  clusteringop::bkinCiToMemberShip
	  (lmembclassdisjsets_Bi,
	   lmembclassdisjsets_BkinCi 
	   );
      }

      partition::PartitionLabelVector
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lovectormmidx_memberShip,
	 loop_outParamGAC.getNumClusterK()
	 );

#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/

   
#ifdef ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003

      auto lpair_chrombestPiMST =
	eac::gaclustering_vktreebinary
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );

      auto &lchrom_best = lpair_chrombestPiMST.first;
    
      partition::PartitionDisjSets
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(graph::component
	 (lpair_chrombestPiMST.second,
	  lpair_chrombestPiMST.first
	  )
	 );
     
      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) loop_outParamGAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) loop_outParamGAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) loop_outParamGAC.getNumClusterK(),
	 DATATYPE_INSTANCES_CLUSTER_K(0)
	 );

      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lomatrixrowt_sumInstCluster,
	 lovectort_numInstClusterK,
	 lpartition_clusters,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );

#endif /*ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003*/


#ifdef ALG_VGA_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2001

     gaencode::ChromVariableLength<DATATYPE_FEATURE,DATATYPE_REAL>
	lchrom_best =
	eac::vga_vkcentroid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );

     mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids
	( lchrom_best.getStringSize() / data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	  data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	  lchrom_best.getString()
	  );

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );
      
#endif /*ALG_VGA_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2001*/
      
    
#ifdef ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids =
	eac::gcuk_vkcentroid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
      mat::MatrixRow<DATATYPE_FEATURE>&  lchrom_best = lomatrixrowt_centroids;
    
      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );
#endif /*ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002*/

    
#ifdef ALG_TGCA_VKCENTROID_HE_AND_TAN_2012

      gaencode::ChromVariableLength<DATATYPE_FEATURE,DATATYPE_REAL>
	lchrom_best =
	eac::tgca_vkcentroid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );

      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids
	( lchrom_best.getStringSize() / data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	  data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	  lchrom_best.getString()
	  );

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );

#endif /*ALG_TGCA_VKCENTROID_HE_AND_TAN_2012*/
  
#if defined(ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012) ||	\
  defined(ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012)
   
      gaencode::ChromosomeGGA<DATATYPE_CLUSTERIDX,DATATYPE_REAL>&& lchrom_best =
	eac::gga_vklabel
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
    
      partition::PartitionLabel
	<DATATYPE_CLUSTERIDX>
	lpartition_clusters
	(lchrom_best.getString(),
	 (uintidx) lpairvec_dataset.first.size(),/*String Size */
	 loop_outParamGAC.getNumClusterK()  
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) loop_outParamGAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) loop_outParamGAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) loop_outParamGAC.getNumClusterK(),
	 DATATYPE_INSTANCES_CLUSTER_K(0)
	 );

      clusteringop::getCentroids
	(lomatrixrowt_centroids,
	 lomatrixrowt_sumInstCluster,
	 lovectort_numInstClusterK,
	 lpartition_clusters,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end()
	 );

    
#endif /*ALG_GGA_VKLABEL_DBINDEX_AGUSTIN_ETAL_2012 ||	\
	 ALG_GGA_VKLABEL_SILHOUETTE_AGUSTIN_ETAL_2012
       */

      /*\cite{Chang:etal:GAclustering:GAGR:2009}*/
#ifdef  ALG_GAGR_FKCENTROID_CHANG_ETAL_2009

      gaencode::ChromFixedLength<DATATYPE_FEATURE,DATATYPE_REAL>&& lchrom_best = 
	eac::gagr_fkcentroid 
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );      

      mat::MatrixRow<DATATYPE_FEATURE> lomatrixrowt_centroids
	((uintidx) linparam_ClusteringGA.getNumClusterK(),
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 lchrom_best.getString()
	 );

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );
#endif /*ALG_GAGR_FKCENTROID_CHANG_ETAL_2009*/

      /*\cite{Lucasius:etal:GAclusteringMedoid:GCA:1993}
       */
#ifdef  ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993
    
      /*(0) calculate matrix dissimilarity
       */
      mat::MatrixTriang<DATATYPE_REAL>&& 
	lmatrixtriagT_dissimilarity = 
	dist::getMatrixDissimilarity
	(lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
    
      gaencode::ChromFixedLength<uintidx,DATATYPE_REAL>&& lchrom_best =
	eac::gca_fkmedoid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lmatrixtriagT_dissimilarity
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids /* Medoids */
	( (uintidx) linparam_ClusteringGA.getNumClusterK(),
	  data::Instance<DATATYPE_FEATURE>::getNumDimensions()
	  );

      clusteringop::initialize
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lchrom_best.begin()
	 );
    
      partition::PartitionMedoids
	<uintidx,
	 DATATYPE_CLUSTERIDX,
	 DATATYPE_REAL
	 >
	lpartition_clusters
	(lchrom_best.getString(),
	 linparam_ClusteringGA.getNumClusterK(),
	 lmatrixtriagT_dissimilarity
	 );
    
#endif /*ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993*/

      /*\cite{Sheng:Xiaohui:GAclusteringMedoid:HKA:2004}
       */
#ifdef ALG_HKA_FKMEDOID_SHENG_LUI2004
    
      /*(0) calculate matrix dissimilarity
       */
      mat::MatrixTriang<DATATYPE_REAL>&& 
	lmatrixtriagT_dissimilarity = 
	dist::getMatrixDissimilarity
	(lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );
    
      gaencode::ChromFixedLength<uintidx,DATATYPE_REAL>&& lchrom_best = 
	eac::hka_fkmedoid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lmatrixtriagT_dissimilarity
	 );

      //std::cout << "HKACROMO: " << lchrom_best << std::endl;

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids /* Medoids */
	( (uintidx) linparam_ClusteringGA.getNumClusterK(),
	  data::Instance<DATATYPE_FEATURE>::getNumDimensions()
	  );

      clusteringop::initialize
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lchrom_best.begin()
	 );

      partition::PartitionMedoids
	<uintidx, 
	 DATATYPE_CLUSTERIDX,
	 DATATYPE_REAL
	 >
	lpartition_clusters
	(lchrom_best.getString(),
	 linparam_ClusteringGA.getNumClusterK(),
	 lmatrixtriagT_dissimilarity
	 );

#endif /*ALG_HKA_FKMEDOID_SHENG_LUI2004*/

#ifdef ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997

      gaencode::ChromosomeBitArray<DATATYPE_BITSIZE,DATATYPE_REAL>&& lchrom_best =
	eac::gaprototypes_fkmedoid
	(loop_outParamGAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );

      std::list<uintidx>&& listui_idxInst =
	lchrom_best.getIdxWithBitOn();
      
      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids    /* Medoids */                         
	( (uintidx) listui_idxInst.size(),
	  data::Instance<DATATYPE_FEATURE>::getNumDimensions()
	  );

      clusteringop::initialize
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 listui_idxInst.begin()
	 );

      auto lpartition_clusters = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );

#endif /*ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997*/

      /*EVALUATE MEASURES ------------------------------------------------------
       */
    
      /*MATCHING MATRIX*/
      sm::ConfusionMatchingMatrix<DATATYPE_INSTANCES_CLUSTER_K>
	lmatchmatrix_confusion;
    

      { /*BEGIN EVALUATE MEASURES*/
      
#ifdef __VERBOSE_YES
	geverbosepc_labelstep = "evaluate measures";
	++geiinparam_verbose;
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << lpc_labeMain << ':' << geverbosepc_labelstep
		    << ":  IN(" << geiinparam_verbose << ")\n";
	  std::ostringstream lostrstream_labelCentroids;
	  lostrstream_labelCentroids << "<CENTROIDS:" << geverbosepc_labelstep;
	  lomatrixrowt_centroids.print(std::cout,lostrstream_labelCentroids.str().c_str(),',',';');
	  std::cout << '\n';
	  lpartition_clusters.print();
	  std::cout << std::endl;
	}
#endif /*__VERBOSE_YES*/
    
    
#ifdef __INSTANCES_WITH_FREQUENCY
	/*-----------BEGIN __INSTANCES_WITH_FREQUENCY*/
	
	/*Calculate confusion matrix or matching matrix  
	  for calculate: Rand Index
	*/  
	if ( linparam_ClusteringGA.getClassInstanceColumn() ) {

	  lmatchmatrix_confusion = 
	    sm::getConfusionMatrix
	    (lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     lpartition_clusters,
	     [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	     -> DATATYPE_INSTANCES_CLUSTER_K
	     {
	       data::InstanceClassFreq
		 <DATATYPE_FEATURE,
		  DATATYPE_INSTANCE_FREQUENCY,
		  DATATYPE_INSTANCES_CLUSTER_K,
		  DATATYPE_CLUSTERIDX>
		 *linstClassFreq_iter =
		 (data::InstanceClassFreq
		  <DATATYPE_FEATURE,
		  DATATYPE_INSTANCE_FREQUENCY,
		  DATATYPE_INSTANCES_CLUSTER_K,
		  DATATYPE_CLUSTERIDX>*)
		 aiinst_iter;
	       return  (DATATYPE_INSTANCES_CLUSTER_K) linstClassFreq_iter->getFrequency(); 
	     },
	     [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	     -> DATATYPE_CLUSTERIDX
	     {

	       data::InstanceClassFreq
		 <DATATYPE_FEATURE,
		  DATATYPE_INSTANCE_FREQUENCY,
		  DATATYPE_INSTANCES_CLUSTER_K,
		  DATATYPE_CLUSTERIDX>
		 *linstClassFreq_iter =
		 (data::InstanceClassFreq
		  <DATATYPE_FEATURE,
		  DATATYPE_INSTANCE_FREQUENCY,
		  DATATYPE_INSTANCES_CLUSTER_K,
		  DATATYPE_CLUSTERIDX>*)
		 aiinst_iter;
	       
	       return linstClassFreq_iter->getClassIdx();
	       
	     }
	     );	  
  
	}
	else {

	  lmatchmatrix_confusion =
	    sm::getMatchingMatrix
	    (lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     lpartition_clusters,
	     [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	     -> DATATYPE_INSTANCES_CLUSTER_K
	     {
	       data::InstanceFreq
		 <DATATYPE_FEATURE,
		  DATATYPE_INSTANCE_FREQUENCY
		  >
		 *lptinstfreq_iter =
		 (data::InstanceFreq
		  <DATATYPE_FEATURE,
		  DATATYPE_INSTANCE_FREQUENCY
		  >*)
		 aiinst_iter;

	       return  (DATATYPE_INSTANCES_CLUSTER_K) lptinstfreq_iter->getFrequency(); 
	     }
	     );  
	}
	/*------------END __INSTANCES_WITH_FREQUENCY
	 */

	
#else /*------------BEGIN __INSTANCES_WITHOUT_FREQUENCY
       */

	//CHECAR EL ALGORITMO TGCA POR QI TERMINA CON 1
	//if (lpartition_clusters.getNumCluster() > 1 ) {
	  if ( linparam_ClusteringGA.getClassInstanceColumn() )  {  

	    lmatchmatrix_confusion = 
	      sm::getConfusionMatrix
	      (lpairvec_dataset.first.begin(),
	       lpairvec_dataset.first.end(),
	       lpartition_clusters,
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	       -> DATATYPE_INSTANCES_CLUSTER_K
	       {
		 return  DATATYPE_INSTANCES_CLUSTER_K(1);
	       },
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	       -> DATATYPE_CLUSTERIDX
	       {
		 data::InstanceClass
		   <DATATYPE_FEATURE,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>
		   *linstclass_iter = 
		   (data::InstanceClass
		    <DATATYPE_FEATURE,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>*)
		   aiinst_iter;
	       
		 return linstclass_iter->getClassIdx();
	       
	       }
	       );	  
  
	  }
	  else {

	    lmatchmatrix_confusion =
	      sm::getMatchingMatrix
	      (lpairvec_dataset.first.begin(),
	       lpairvec_dataset.first.end(),
	       lpartition_clusters,
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	       -> DATATYPE_INSTANCES_CLUSTER_K
	       {
		 return  DATATYPE_INSTANCES_CLUSTER_K(1);
	       }
	       );
	  
	  }

	  //} //IF CUANDO TERMINA CON 1

#endif /*-----------END __INSTANCES_WITHOUT_FREQUENCY
	*/

#ifdef __VERBOSE_YES
	if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	  std::cout << lpc_labeMain << ':' << geverbosepc_labelstep
		    << ": OUT(" << geiinparam_verbose << ')'
		    << std::endl;
	}
	--geiinparam_verbose;
	geverbosepc_labelstep = lpc_labeMain;
#endif //__VERBOSE_YES//

      } //END EVALUATE MEASURES//


      /*BEGIN TEST METRIC-------------------------------------------------------
       */
    
      /*MATCHING MATRIX TEST*/
      sm::ConfusionMatchingMatrix<DATATYPE_INSTANCES_CLUSTER_K>
	lmatchmatrix_confusionTest;

      /* Make centroids partition
       */
      auto lpartition_clustersTest = 
	partition::makePartition
	(lomatrixrowt_centroids,
	 lpairvec_dataset.second.begin(),
	 lpairvec_dataset.second.end(),
	 DATATYPE_CLUSTERIDX(lomatrixrowt_centroids.getNumRows()),
	 *pfunct2p_distAlg
	 );
      
      { /*BEGIN  METRIC TEST*/

	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	
#ifdef __VERBOSE_YES
	  geverbosepc_labelstep = "evaluate measures test";
	  ++geiinparam_verbose;
	  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	    std::cout << lpc_labeMain << ':' << geverbosepc_labelstep
		      << ":  IN(" << geiinparam_verbose << ')'
		      << std::endl;
	  }  
#endif /*__VERBOSE_YES*/

	
#ifdef __INSTANCES_WITH_FREQUENCY

	  /*Calculate confusion matrix or matching matrix  
	    for calculate: Rand Index
	  */
	  if ( linparam_ClusteringGA.getClassInstanceColumn() ) {

	    lmatchmatrix_confusionTest =
	      sm::getConfusionMatrix
	      (lpairvec_dataset.second.begin(),
	       lpairvec_dataset.second.end(),
	       lpartition_clustersTest,
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	       -> DATATYPE_INSTANCES_CLUSTER_K
	       {
		 data::InstanceClassFreq
		   <DATATYPE_FEATURE,
		    DATATYPE_INSTANCE_FREQUENCY,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>
		   *linstClassFreq_iter =
		   (data::InstanceClassFreq
		    <DATATYPE_FEATURE,
		    DATATYPE_INSTANCE_FREQUENCY,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>*)
		   aiinst_iter;
		 return  (DATATYPE_INSTANCES_CLUSTER_K) linstClassFreq_iter->getFrequency(); 
	       },
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter ) -> DATATYPE_CLUSTERIDX
	       {

		 data::InstanceClassFreq
		   <DATATYPE_FEATURE,
		    DATATYPE_INSTANCE_FREQUENCY,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>
		   *linstClassFreq_iter =
		   (data::InstanceClassFreq
		    <DATATYPE_FEATURE,
		    DATATYPE_INSTANCE_FREQUENCY,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>*)
		   aiinst_iter;
	       
		 return linstClassFreq_iter->getClassIdx();
	       
	       }
	       );	  
  
	  }
	  else {

	    lmatchmatrix_confusionTest =
	      sm::getMatchingMatrix
	      (lpairvec_dataset.second.begin(),
	       lpairvec_dataset.second.end(),
	       lpartition_clustersTest,
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter ) -> DATATYPE_INSTANCES_CLUSTER_K
	       {
		 data::InstanceFreq
		   <DATATYPE_FEATURE,
		    DATATYPE_INSTANCE_FREQUENCY
		    >
		   *lptinstfreq_iter =
		   (data::InstanceFreq
		    <DATATYPE_FEATURE,
		    DATATYPE_INSTANCE_FREQUENCY
		    >*)
		   aiinst_iter;

		 return  (DATATYPE_INSTANCES_CLUSTER_K) lptinstfreq_iter->getFrequency(); 
	       }
	       );  
	  }

#else /*__INSTANCES_WITHOUT_FREQUENCY FOR TEST*/
      
	  if ( linparam_ClusteringGA.getClassInstanceColumn() )  {  

	    lmatchmatrix_confusionTest = 
	      sm::getConfusionMatrix
	      (lpairvec_dataset.second.begin(),
	       lpairvec_dataset.second.end(),
	       lpartition_clustersTest,
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	       -> DATATYPE_INSTANCES_CLUSTER_K
	       {
		 return  DATATYPE_INSTANCES_CLUSTER_K(1);
	       },
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter )
	       -> DATATYPE_CLUSTERIDX
	       {
		 data::InstanceClass
		   <DATATYPE_FEATURE,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>
		   *linstclass_iter = 
		   (data::InstanceClass
		    <DATATYPE_FEATURE,
		    DATATYPE_INSTANCES_CLUSTER_K,
		    DATATYPE_CLUSTERIDX>*)
		   aiinst_iter;
	       
		 return linstclass_iter->getClassIdx();
	       
	       }
	       );	  
  
	  }
	  else {

	    lmatchmatrix_confusionTest =
	      sm::getMatchingMatrix
	      (lpairvec_dataset.second.begin(),
	       lpairvec_dataset.second.end(),
	       lpartition_clustersTest,
	       [](const data::Instance<DATATYPE_FEATURE>* aiinst_iter ) -> DATATYPE_INSTANCES_CLUSTER_K
	       {
		 return  DATATYPE_INSTANCES_CLUSTER_K(1);
	       }
	       );
	  
	  }
	
#endif /*__INSTANCES_WITH_FREQUENCY*/
	
	  { /*BEGIN  TEST METRIC RUNNING
	     */

	    std::vector<DATATYPE_INSTANCES_CLUSTER_K>&& lvectort_numInstClusterKTest =
	      lmatchmatrix_confusionTest.getInstClusterK();
     
	    auto li_matrixPartitionClusterNullTest =
	      std::count_if
	      (lvectort_numInstClusterKTest.begin(),
	       lvectort_numInstClusterKTest.end(),
	       [](DATATYPE_INSTANCES_CLUSTER_K liter_numClusterK)
	       {return liter_numClusterK == 0;}
	       );

	    /*CHECK CONDITION END
	     */
	    loop_outParamGAC.setNumClusterKTest
	      (DATATYPE_CLUSTERIDX(lvectort_numInstClusterKTest.size())
	       - DATATYPE_CLUSTERIDX(li_matrixPartitionClusterNullTest)
	       );
	    if ( li_matrixPartitionClusterNullTest != 0 )  {
	      loop_outParamGAC.setEndingConditionTest(false);
	    }
	    else  {
	      loop_outParamGAC.setEndingConditionTest(true);
	    }
#ifdef __VERBOSE_YES
	    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	      std::cout << lpc_labeMain << ':' << geverbosepc_labelstep
			<< ": OUT(" << geiinparam_verbose << ')'
			<< std::endl;
	    }
	    --geiinparam_verbose;
	    geverbosepc_labelstep = lpc_labeMain;
#endif //__VERBOSE_YES//

	  } /*BEGIN  TEST METRIC RUNNING
	     */       
	} /*IF*/

      } /*MODIFICAR END  METRIC TEST*/
    
      /*MODIFICAR END TEST METRIC-------------------------------------------------------
       */

      /*OUT: PRINT------------------------------------------------------------
       */
    
      { /*BEGIN PRINT PARAMETERS*/
       
	ds::PartitionLinkedNumInst
	  <DATATYPE_CLUSTERIDX,DATATYPE_INSTANCES_CLUSTER_K>&&
	  lpartlinknuminst_memberShip =
	  ds::getPartitionLinkedNumInst
	  (lpairvec_dataset.first.begin(),
	   lpairvec_dataset.first.end(),
	   lpartition_clusters,
	   [&](data::Instance<DATATYPE_FEATURE>* liter_inst)
	   {
	     return DATATYPE_INSTANCES_CLUSTER_K(1);
	   }
	   );
	
	/* SED (0)
	   */
	  std::pair<DATATYPE_REAL,bool> lpair_sed =
	    um::SSE
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     *pfunct2p_distEuclidean
	     );
	    
	  loop_outParamlusteringMetric.setMetricFunc
	    (inout::SED,
	     lpair_sed.first
	     );

	/* SSE (1)
	   */
	  std::pair<DATATYPE_REAL,bool> lpair_sse =
	     um::SSE
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     *pfunct2p_distEuclideanSq
	     );

	  loop_outParamlusteringMetric.setMetricFunc
	    (inout::SSE,
	     lpair_sse.first
	     );

	  /* Distortion (2)
	   */
	  loop_outParamlusteringMetric.setMetricFunc
	    (inout::Distortion,
	     (lmatchmatrix_confusion.getNumObjetos()  * 
	      data::Instance<DATATYPE_FEATURE>::getNumDimensions() != 0)? 
	      lpair_sse.first / 
	     ((DATATYPE_REAL) lmatchmatrix_confusion.getNumObjetos() * 
	      (DATATYPE_REAL) data::Instance<DATATYPE_FEATURE>::getNumDimensions())
	     :measuare_undefSSE(DATATYPE_REAL)
	     );

	  /*J1 (3)
	   */
	   loop_outParamlusteringMetric.setMetricFunc
	    (inout::J1,
	     lpair_sse.first
	     );
	
	/*
	  RECUPERAR ESTA PARA CUANDO SE QUIERA INSTANCIAS CON FRECUENCIA
	  loop_outParamlusteringMetric.setMetricFunc
	  (inout::IndexI,
	  um::indexI
	  (lmetrict_e1,
	  lmetrict_SSE, //lt_SSw
	  lmetrict_maxDistCjCjp,//lt_DkMax,
	  (uintidx) loop_outParamlusteringMetric.getNumClusterK()
	  ) 
	  );
	*/

	/*CS measure (5) OK
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::CSmeasure,
	   um::CSmeasure
	   (lpairvec_dataset.first.begin(),
	    lomatrixrowt_centroids,
	    lpartlinknuminst_memberShip,
	    *pfunct2p_distEuclidean
	    )
	   );

	/*DunnIndex (6) OK
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::DunnIndex,
	   um::DunnIndex
	   (lpairvec_dataset.first.begin(),
	    lpartlinknuminst_memberShip,
	    *pfunct2p_distEuclidean
	    )
	   );
	
	/*SDunnIndex (7)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::SDunnIndex,
	   um::simplifiedDunnIndex
	   (lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpartlinknuminst_memberShip,
	    *pfunct2p_distEuclidean
	    )
           );
	/*	
#ifdef __VERBOSE_YES
	geiinparam_verbose = __VERBOSE_CHECK;
#endif //__VERBOSE_YES
	*/

	/*Silhouette (8) OK
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Silhouette,
	   um::silhouette
	   (lpairvec_dataset.first.begin(),
	    lpartlinknuminst_memberShip,
	    *pfunct2p_distEuclidean
	    )
	   );
	
	/*
#ifdef __VERBOSE_YES
	mat::MatrixTriang<DATATYPE_REAL>&& 
	  lmatrixtriagT_dissimilarityCheck = 
	  clusteringop::getMatrixDissimilarity
	  (lpairvec_dataset.first.begin(),
	   lpairvec_dataset.first.end(),
	   *pfunct2p_distEuclidean
	   );
	  
	std::cout << "Other um::silhouette: " <<
	  um::silhouette
	  (lmatrixtriagT_dissimilarityCheck,
	   lpartlinknuminst_memberShip
	   );
	  
	geiinparam_verbose =  __VERBOSE_NOTCHECK;
#endif //__VERBOSE_YES	 
	*/
		
	{ /*SSilhouette (9) OK
	   */	  
	  std::vector<DATATYPE_INSTANCES_CLUSTER_K>&& lvectort_numInstClusterK =
	  lmatchmatrix_confusion.getInstClusterK();

	  std::vector<DATATYPE_REAL>&&  lvectort_partialSilhouette =
	    um::simplifiedSilhouette
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     lpartition_clusters,
	     lvectort_numInstClusterK,
	     *pfunct2p_distEuclidean
	     );

	  DATATYPE_REAL lmetrict_sumPartialSilhouette = 
	    interfacesse::sum
	    (lvectort_partialSilhouette.data(),
	     (uintidx) lvectort_partialSilhouette.size()
	     );

	  loop_outParamlusteringMetric.setMetricFunc
	    (inout::SSilhouette,
	     (lvectort_partialSilhouette.size() != 0)?
	     lmetrict_sumPartialSilhouette /(DATATYPE_REAL) lvectort_partialSilhouette.size():
	     measuare_undefSilhouette(DATATYPE_REAL)
	     );
	} /*END SSilhouette (9)*/

	/*DBindex (10) OK 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::DBindex,
	   um::dbindex
	   (lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    lpartition_clusters,
	    *pfunct2p_distEuclidean
	    )
	   );

	/*Variance Ratio Criterion (11) OK
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::VRC,
	   um::VRC
	   (lomatrixrowt_centroids,
            lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),   
	    lpartition_clusters,
	    *pfunct2p_distEuclideanSq
	    )
	   );  

	/*WB-index (12)*/

	loop_outParamlusteringMetric.setMetricFunc
	  (inout::WBIndex,
	   um::WBIndex
	   (lomatrixrowt_centroids,
            lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),   
	    lpartition_clusters,
	    *pfunct2p_distEuclideanSq
	    )
	   );

	/*SSB (13)*/
	{
	   std::vector<DATATYPE_INSTANCES_CLUSTER_K>&& lvectort_numInstClusterK =
	     lmatchmatrix_confusion.getInstClusterK();

	   loop_outParamlusteringMetric.setMetricFunc
	     (inout::SSB,
	      um::ssb
	      (lomatrixrowt_centroids,
	       larray_meanFeactures,
	       lvectort_numInstClusterK,
	       *pfunct2p_distEuclideanSq
	       )
	      );
	}

	/*Score function (14)
	  The higher the value of the SF , the more suitable 
	  the number of clusters
	  \cite{Sandro:Benny:IanF:ClusteringMeasure:2007}
	*/
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::ScoreFunction, 
	   um::scoreFunction
	   (lomatrixrowt_centroids,
            lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),   
	    lpartition_clusters,
	    *pfunct2p_distEuclidean
	    )
	   );
	
	/*Rand Index (15): 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::RandIndex,
	   sm::randIndex
	   <DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>
	   (lmatchmatrix_confusion)
	   );

	/*Purity (16)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Purity,
	   sm::purity
	   <DATATYPE_REAL,
	   DATATYPE_INSTANCES_CLUSTER_K>
	   (lmatchmatrix_confusion)
	   );

	/*F-measure (17)
         */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Fmeasure,
	   sm::fmeasure
	   <DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>
	   (lmatchmatrix_confusion)
	   );
	
	/*Jaccard index (18):  
	    Jaccard index is quite similar to Rand index, but
	    this measure does not consider the number of correct assign-
	    ments when two elements are assigned to different clusters.
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::JaccardIndex,
	   sm::jaccardIndex
	   <DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K> //lmatchmatrix_confusion.getAgreementObjects(),
	   (lmatchmatrix_confusion)
	   );
	
        /*Precision (19): 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::precision,
	   sm::precision
	   <DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>
	   (lmatchmatrix_confusion)
	   );
	
        /*Recall (20):                                                                                                
         */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::recall,
	   sm::recall
	   <DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>
	   (lmatchmatrix_confusion)
	   );

	/*Misclassified (21)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Misclassified,
	   (lomatrixrowt_centroids.getNumRows() > 0)?
	   (DATATYPE_REAL)lmatchmatrix_confusion.getMisclassified()
	   : (DATATYPE_REAL) lmatchmatrix_confusion.getNumObjetos()
	   );

	/*Pairs a (22) 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Pairs_a,
	   lmatchmatrix_confusion.getSomeClassUSomeClusterV() 
	   );

	/*Pairs b (23) 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Pairs_b,
	   lmatchmatrix_confusion.getSomeClassUDiffClusterV() 
	   );

	/*Pairs c (24) 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Pairs_c,
	   lmatchmatrix_confusion.getDiffClassUSomeClusterV()
	   );

	/*Pairs d (25) 
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Pairs_d,
	   lmatchmatrix_confusion.getDiffClassUDiffClusterV()
	   );

	loop_outParamlusteringMetric.setPercentageSensitivity
	  (sm::getStrSenSpe
	   (lmatchmatrix_confusion, 
	    true,
	    DATATYPE_REAL(100.0))
	   );

	loop_outParamlusteringMetric.setPercentageSpecificity
	  (sm::getStrSenSpe
	   (lmatchmatrix_confusion,
	    false,
	    DATATYPE_REAL(100.0))
	   );
	

	/* Fuzzy measures 26-29
	 */
	mat::MatrixRow<DATATYPE_REAL>&& lmatrixt_u =
	  clusteringop::fuzzyPartition
	  (lomatrixrowt_centroids,
	   lpairvec_dataset.first.begin(),
	   lpairvec_dataset.first.end(),
	   2.0,
	   *pfunct2p_distEuclideanSq
	   );

	/*Jm (26)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Jm,
	   um::jm
	   (lmatrixt_u,
	    lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    2.0,
	    *pfunct2p_distEuclideanSq
	    )
	   );


	/*INDEX I (27) OK
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::IndexI,
	   um::indexI
	   (lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    lpartition_clusters,
	    *pfunct2p_distEuclidean
	    )
	   );

	/*INDEX I (26) OK
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::IndexI_cFuzzy,
	   um::indexI
	   (lmatrixt_u,
	    lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    *pfunct2p_distEuclidean
	    )
	   );

	/*Xie-Beni index (27)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::XieBeniIndex,
	   um::xb
	   (lmatrixt_u,
	    lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    *pfunct2p_distEuclideanSq
	    )
	   );

	loop_outParamlusteringMetric.setMetricFunc
	  (inout::XieBeniIndex_crisp,
	   um::xb
	   (lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    lpartition_clusters,
	    *pfunct2p_distEuclideanSq
	    )
	   );


	/*Entropy (28)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Entropy,
	   um::entropy(lmatrixt_u)
	   );

	/*Partition coefficient (29)
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::PartitionCoefficient,
	   um::partitioncoefficient(lmatrixt_u)
	   );
      
	/*Overlap
	 */
	loop_outParamlusteringMetric.setMetricFunc
	  (inout::Overlap,
	   um::overlap
	   (lomatrixrowt_centroids,
	    lpairvec_dataset.first.begin(),
	    lpairvec_dataset.first.end(),
	    lpartlinknuminst_memberShip,
	    *pfunct2p_distEuclidean
	    )
	   );
    
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {

	  ds::PartitionLinkedNumInst
	    <DATATYPE_CLUSTERIDX,DATATYPE_INSTANCES_CLUSTER_K>&&
	    lpartlinknuminst_memberShipTest =
	    ds::getPartitionLinkedNumInst
	    (lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),
	     lpartition_clustersTest,
	     [&](data::Instance<DATATYPE_FEATURE>* liter_inst)
	     {
	       return DATATYPE_INSTANCES_CLUSTER_K(1);
	     }
	     );


	  /* SED (0)
	   */
	  std::pair<DATATYPE_REAL,bool> lpair_sedTest =
	    um::SSE
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),
	     *pfunct2p_distEuclidean
	     );
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::SED,
	     lpair_sedTest.first
	     );

	  /* SSE (1)
	   */
	  std::pair<DATATYPE_REAL,bool> lpair_sseTest =
	    um::SSE
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),
	     *pfunct2p_distEuclideanSq
	     );
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::SSE,
	     lpair_sseTest.first
	     );

	  /* Distortion (2)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Distortion,
	     (lmatchmatrix_confusionTest.getNumObjetos()  * data::Instance<DATATYPE_FEATURE>::getNumDimensions() != 0)? 
	      lpair_sseTest.first / 
	     ((DATATYPE_REAL) lmatchmatrix_confusionTest.getNumObjetos()* 
	      (DATATYPE_REAL) data::Instance<DATATYPE_FEATURE>::getNumDimensions())
	     :measuare_undefSSE(DATATYPE_REAL)
	     );

	  /*J1 (3)
	   */
	   loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::J1,
	     lpair_sseTest.first
	     );

	  /*CS measure (5) OK
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::CSmeasure,
	     um::CSmeasure
	     (lpairvec_dataset.second.begin(),
	      lomatrixrowt_centroids,
	      lpartlinknuminst_memberShipTest,
	      *pfunct2p_distEuclidean
	      )
	     );

	  /*DunnIndex (6) OK
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::DunnIndex,
	     um::DunnIndex
	     (lpairvec_dataset.second.begin(),
	      lpartlinknuminst_memberShipTest,
	      *pfunct2p_distEuclidean,
	      true
	      )
	     );
	  
	  /*SDunnIndex (7)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::SDunnIndex,
	     um::simplifiedDunnIndex
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpartlinknuminst_memberShipTest,
	      *pfunct2p_distEuclidean,
	      true
	      )
	     );
	 
	  /*Silhouette (8) OK
	   */  
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Silhouette,
	     um::silhouette
	     (lpairvec_dataset.second.begin(),
	      lpartlinknuminst_memberShipTest,
	      *pfunct2p_distEuclidean
	      )
	     );
	  
	  {/*SSilhouette (9) OK
	    */
	    std::vector<DATATYPE_INSTANCES_CLUSTER_K>&& lvectort_numInstClusterKTest =
	      lmatchmatrix_confusionTest.getInstClusterK();

	    std::vector<DATATYPE_REAL>&&  lvectort_partialSilhouetteTest =
	      um::simplifiedSilhouette
	      (lomatrixrowt_centroids,
	       lpairvec_dataset.second.begin(),
	       lpairvec_dataset.second.end(),
	       lpartition_clustersTest,
	       lvectort_numInstClusterKTest,
	       *pfunct2p_distEuclidean
	       );

	    DATATYPE_REAL lmetrict_partialSilhouetteTest = 
	      interfacesse::sum
	      (lvectort_partialSilhouetteTest.data(),
	       (uintidx) lvectort_partialSilhouetteTest.size()
	       );

	    loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::SSilhouette,
	     (lvectort_partialSilhouetteTest.size() != 0)?
	     lmetrict_partialSilhouetteTest /(DATATYPE_REAL) lvectort_partialSilhouetteTest.size():
	     measuare_undefSilhouette(DATATYPE_REAL)
	     );
	  }

	  /*DBindex (10) OK 
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::DBindex,
	     um::dbindex
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      lpartition_clustersTest,
	      *pfunct2p_distEuclidean
	      )
	     );

	  /*Variance Ratio Criterion (11) OK
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::VRC,
	     um::VRCreeval
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),   
	      lpartition_clustersTest,
	      *pfunct2p_distEuclideanSq,
	      true && !lpartlinknuminst_memberShip.haveNullCluster()
	      )
	     );

	  /*WB-index (12)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::WBIndex,
	     um::WBIndexreeval
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),  
	      lpartition_clustersTest,
	      *pfunct2p_distEuclideanSq
	      )
	     );

	  /*SSB (13)*/
	  {
	    std::vector<DATATYPE_INSTANCES_CLUSTER_K>&& lvectort_numInstClusterKTest =
	      lmatchmatrix_confusionTest.getInstClusterK();

	    loop_outParamlusteringMetric.setMetricFuncTest
	      (inout::SSB,
	       um::ssbreeval
	       (lomatrixrowt_centroids,
		lpairvec_dataset.second.begin(),
		lpairvec_dataset.second.end(),
		lvectort_numInstClusterKTest,
		*pfunct2p_distEuclideanSq
		)
	       );
	    
	  }
	 
	  /*Score function (14)
	    The higher the value of the SF , the more suitable 
	    the number of clusters
	    \cite{Sandro:Benny:IanF:ClusteringMeasure:2007}
	  */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::ScoreFunction,
	     um::scoreFunctionreeval
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),   
	      lpartition_clustersTest,
	      *pfunct2p_distEuclidean
	      )
	     );

	  /*Rand Index (15) OK
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::RandIndex,
	     sm::randIndex<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>
	     (lmatchmatrix_confusionTest)
	     );

	  /*Purity (16)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Purity,
	     sm::purity
	     <DATATYPE_REAL,
	     DATATYPE_INSTANCES_CLUSTER_K>
	     (lmatchmatrix_confusionTest)
	     );

	  /*F-measure (17):
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Fmeasure,
	     sm::fmeasure<DATATYPE_REAL,DATATYPE_INSTANCES_CLUSTER_K>
	     (lmatchmatrix_confusionTest)
	     );

	  /*Jaccard index (18):  
	    Jaccard index is quite similar to Rand index, but
	    this measure does not consider the number of correct assign-
	    ments when two elements are assigned to different clusters.
	  */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::JaccardIndex,
	     sm::jaccardIndex<DATATYPE_REAL,DATATYPE_INSTANCES_CLUSTER_K>
	     (lmatchmatrix_confusionTest)
	     );
	  
	  /*Precision (19): 
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::precision,
	     sm::precision<DATATYPE_REAL,DATATYPE_INSTANCES_CLUSTER_K>
	     (lmatchmatrix_confusionTest)
	     );
	  
	  /*Recall (20):
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::recall,
	     sm::recall<DATATYPE_REAL,DATATYPE_INSTANCES_CLUSTER_K>
	     (lmatchmatrix_confusionTest)
	     );

	  /*Misclassified (21)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Misclassified,
	     (lomatrixrowt_centroids.getNumRows() > 0)? 
	     (DATATYPE_REAL)lmatchmatrix_confusionTest.getMisclassified()
	     : (DATATYPE_REAL) lmatchmatrix_confusionTest.getNumObjetos() 
	     );

	  /*Pairs a (22)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Pairs_a,
	     lmatchmatrix_confusionTest.getSomeClassUSomeClusterV()
	     );

	  /*Pairs b (23)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Pairs_b,
	     lmatchmatrix_confusionTest.getSomeClassUDiffClusterV()
	     );

	  /*Pairs c (24)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Pairs_c,
	     lmatchmatrix_confusionTest.getDiffClassUSomeClusterV()
	     );

	  /*Pairs d (25)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Pairs_d,
	     lmatchmatrix_confusionTest.getDiffClassUDiffClusterV()
	     );
	 
	  loop_outParamlusteringMetric.setPercentageSensitivityTest
	    (sm::getStrSenSpe
	     (lmatchmatrix_confusionTest,
	      true,
	      DATATYPE_REAL(100.0))
	     );

	  loop_outParamlusteringMetric.setPercentageSpecificityTest
	    (sm::getStrSenSpe
	     (lmatchmatrix_confusionTest,
	      false,
	      DATATYPE_REAL(100.0))
	     );

	  /* Fuzzy measures 26-29
	   */
	  mat::MatrixRow<DATATYPE_REAL>&& lmatrixt_uTest =
	    clusteringop::fuzzyPartition
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),
	     2.0,
	     *pfunct2p_distEuclideanSq
	     );


	  /*Jm (26)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Jm,
	     um::jm
	     (lmatrixt_uTest,
	      lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      2.0,
	      *pfunct2p_distEuclideanSq
	      )
	     );


	  /*INDEX I (4) OK
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::IndexI,
	     um::indexIreeval
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      lpartition_clustersTest,
	      *pfunct2p_distEuclidean
	      ) 
	     );

	   loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::IndexI_cFuzzy,
	     um::indexIreeval
	     (lmatrixt_uTest,
	      lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      *pfunct2p_distEuclidean
	      ) 
	     );

	  /*Xie-Beni index (27)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::XieBeniIndex,
	     um::xb
	     (lmatrixt_uTest,
	      lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      *pfunct2p_distEuclideanSq
	      )
	     );

	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::XieBeniIndex_crisp,
	     um::xb
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      lpartition_clustersTest,
	      *pfunct2p_distEuclideanSq
	      )
	     );
	  
	  /*Entropy (28)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Entropy,
	     um::entropy(lmatrixt_uTest)
	     );

	  /*Partition coefficient (29)
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::PartitionCoefficient,
	     um::partitioncoefficient(lmatrixt_uTest)
	     );

	  /*Overlap 
	   */
	  loop_outParamlusteringMetric.setMetricFuncTest
	    (inout::Overlap,
	     um::overlap
	     (lomatrixrowt_centroids,
	      lpairvec_dataset.second.begin(),
	      lpairvec_dataset.second.end(),
	      lpartlinknuminst_memberShipTest,
	      *pfunct2p_distEuclidean
	      )
	     );
	 
	} //IF TEST
      
	
	 inout::OutFileName   lofn_filename;
	 std::ostream& lostream_outparam = 
	   lofn_filename.openFile(linparam_ClusteringGA.getFileNameTimesRun());
	 lostream_outparam << std::boolalpha;

	 
	 if ( linparam_ClusteringGA.getPrintMulLine() ) {

	   lostream_outparam
	     << "\nIN:\n"
	     << "    Algorithmo name: "
	     << linparam_ClusteringGA.getAlgorithmoName() << '\n'
	     << "           Based on: "
	     << linparam_ClusteringGA.getAlgorithmoAuthor() << '\n'
	     << "        Metric used: "
	     << loop_outParamGAC.getNameUsedObjetiveFunc()
	     << "\n\n"
	  
	     << "           Data set: "
	     << linparam_ClusteringGA.getCurrentFileInstance() << '\n'
	     << "Number of instances: "
	     << linparam_ClusteringGA.getNumInstances() << '\n'
	     << "         Dimensions: "
	     << linparam_ClusteringGA.getNumDimensionsInstances() << '\n';

	   if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	     lostream_outparam
	       << '\n'
	       << "      Data set test: "
	       << linparam_ClusteringGA.getCurrentFileInstanceTest() << '\n'
	       << "Number of instances: "
	       << linparam_ClusteringGA.getNumInstancesTest() << '\n';
	   }
	   lostream_outparam
	     << '\n'
	     << "       Random seed: "
	     <<  linparam_ClusteringGA.getRandomSeed() << '\n'
	     << '\n';

	   lostream_outparam << "OUT:\n\n";

	   lchrom_best.print(lostream_outparam,"CROMOSOME:BEST",',',';');
		      
	   lostream_outparam << "\n\n";

	   lostream_outparam << "       Cluster number (K): "
			     << loop_outParamGAC.getNumClusterK() << "\n"
			     <<  std::right << std::setw(25)
			     <<  loop_outParamGAC.getNameUsedObjetiveFunc() << ": "
			     << loop_outParamGAC.getObjetiveFuncRun() << "\n";

	   if ( loop_outParamGAC.getUsedObjetiveFunc() != inout::SED ) {
	     /*std::pair<DATATYPE_REAL,bool> lpair_sed =
	       um::SSE
	       (lomatrixrowt_centroids,
		lpairvec_dataset.first.begin(),
		lpairvec_dataset.first.end(),
		*pfunct2p_distEuclidean
		);*/
	     lostream_outparam << "                      SED: "
			       << loop_outParamlusteringMetric.getObjetiveFunc(inout::SED);
	       //		       << lpair_sed.first;
	     /*if ( lpair_sed.second == false ) 
	       lostream_outparam << " Has group without objects";
	     */
	     lostream_outparam << '\n';
	   }

	   if ( loop_outParamGAC.getUsedObjetiveFunc() != inout::DBindex ) {
	     lostream_outparam
	       << "                 DB-index: "
	       << loop_outParamlusteringMetric.getObjetiveFunc(inout::DBindex)
	       /*<<
	       um::dbindex
	       (lomatrixrowt_centroids,
		lpairvec_dataset.first.begin(),
		lpairvec_dataset.first.end(),
		lpartition_clusters,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';
	   }

	   if ( loop_outParamGAC.getUsedObjetiveFunc() != inout::Silhouette ) {
	     lostream_outparam
	       << "               Silhouette: "	
	       << loop_outParamlusteringMetric.getObjetiveFunc(inout::Silhouette)

	       /*um::silhouette
	       (lpairvec_dataset.first.begin(),
		lpartlinknuminst_memberShip,
		*pfunct2p_distEuclidean
		)
	       */
	       << '\n';
	   }

	   if ( loop_outParamGAC.getUsedObjetiveFunc() != inout::VRC ) {
	     lostream_outparam
	       << "                      VRC: "	
	       << loop_outParamlusteringMetric.getObjetiveFunc(inout::VRC)
	       /*um::VRC
	       (lomatrixrowt_centroids,
		lpairvec_dataset.first.begin(),
		lpairvec_dataset.first.end(),   
		lpartition_clusters,
		*pfunct2p_distEuclideanSq
		)*/
	       << '\n';
	   }
      
	   if ( loop_outParamGAC.getUsedObjetiveFunc() != inout::CSmeasure ) {
	     lostream_outparam
	       << "               CS measure: "
	       << loop_outParamlusteringMetric.getObjetiveFunc(inout::CSmeasure)
	       /* um::CSmeasure
	       (lpairvec_dataset.first.begin(),
		lomatrixrowt_centroids,
		lpartlinknuminst_memberShip,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';
	   }
      
	   if ( loop_outParamGAC.getUsedObjetiveFunc() != inout::DunnIndex ) {
	     lostream_outparam
	       << "             Dunn's index: "
	       << loop_outParamlusteringMetric.getObjetiveFunc(inout::DunnIndex)
	       /* um::DunnIndex
	       (lpairvec_dataset.first.begin(),
		lpartlinknuminst_memberShip,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';
	   }
    
	   if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {

	     /*ds::PartitionLinkedNumInst
	       <DATATYPE_CLUSTERIDX,DATATYPE_INSTANCES_CLUSTER_K>&&
	       lpartlinknuminst_memberShipTest =
	       ds::getPartitionLinkedNumInst
	       (lpairvec_dataset.second.begin(),
		lpairvec_dataset.second.end(),
		lpartition_clusters,
		[&](data::Instance<DATATYPE_FEATURE>* liter_inst)
		{
		  return DATATYPE_INSTANCES_CLUSTER_K(1);
		}
		);
	     std::pair<DATATYPE_REAL,bool> lpair_sedTest =
	       um::SSE
	       (lomatrixrowt_centroids,
		lpairvec_dataset.second.begin(),
		lpairvec_dataset.second.end(),
		*pfunct2p_distEuclidean
		);
	     */
	     lostream_outparam << '\n';
	     lostream_outparam
	       << "            Test data SED: "
	       << loop_outParamlusteringMetric.getObjetiveFuncTest(inout::SED);
	       /*  lpair_sedTest.first;
	     if ( lpair_sedTest.second == false ) 
	     lostream_outparam << " Has group without objects";*/
	     lostream_outparam << '\n';
	     
	     lostream_outparam
	       << "       Test data DB-index: "
	       << loop_outParamlusteringMetric.getObjetiveFuncTest(inout::DBindex)
	       /*um::dbindex
	       (lomatrixrowt_centroids,
		lpairvec_dataset.second.begin(),
		lpairvec_dataset.second.end(),
		lpartition_clusters,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';

	     lostream_outparam
	       << "     Test data Silhouette: "	
	       << loop_outParamlusteringMetric.getObjetiveFuncTest(inout::Silhouette)
	       /* um::silhouette
	       (lpairvec_dataset.second.begin(),
		lpartlinknuminst_memberShipTest,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';

	     lostream_outparam
	       << "            Test data VRC: "	
	       << loop_outParamlusteringMetric.getObjetiveFuncTest(inout::VRC)
	       /*um::VRC
	       (lomatrixrowt_centroids,
		lpairvec_dataset.second.begin(),
		lpairvec_dataset.second.end(),   
		lpartition_clusters,
		*pfunct2p_distEuclideanSq
		)*/
	       << '\n';

	     lostream_outparam
	       << "     Test data CS measure: "
	       << loop_outParamlusteringMetric.getObjetiveFuncTest(inout::CSmeasure)
	       /* um::CSmeasure
	       (lpairvec_dataset.second.begin(),
		lomatrixrowt_centroids,
		lpartlinknuminst_memberShipTest,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';
     
      
	     lostream_outparam
	       << "   Test data Dunn's index: "
	       << loop_outParamlusteringMetric.getObjetiveFuncTest(inout::DunnIndex)
	       /* um::DunnIndex
	       (lpairvec_dataset.second.begin(),
		lpartlinknuminst_memberShipTest,
		*pfunct2p_distEuclidean
		)*/
	       << '\n';
	   }
      
	   lostream_outparam << "     Execution time (seg): "
			     << loop_outParamGAC.getAlgorithmRunTime() << '\n'
			     << "Generations find the best: "
			     << loop_outParamGAC.getIterationGetsBest() << '\n';
    
	   lostream_outparam << std::endl;
	 }
	
	 else { /*FOR RUNTIME TEST*/

	   lostream_outparam.precision(COMMON_COUT_PRECISION);
	 

	 lostream_outparam << "_inout" <<  inout::OutFileName::getDelim() << "out" << inout::OutFileName::getDelim();
	 linparam_ClusteringGA.print(lostream_outparam,inout::OutFileName::getDelim());
	 //lostream_outparam <<  inout::OutFileName::getDelim();
	 loop_outParamGAC.print(lostream_outparam,inout::OutFileName::getDelim());
	 //lostream_outparam <<  inout::OutFileName::getDelim();
	 loop_outParamlusteringMetric.print(lostream_outparam,inout::OutFileName::getDelim());
	 lostream_outparam << std::endl;
	 

	   }
	 lofn_filename.closeFile();


      } /*END PRINT PARAMETERS*/
    
    
      /*BEGIN PRINT CENTROIDS*/
      if ( linparam_ClusteringGA.getOutFileCentroids() != NULL )  { 

	inout::OutFileName lofn_filenamecentroids;
	std::ostream&      lostream_outcentroids = 
	  lofn_filenamecentroids.openFile(linparam_ClusteringGA.getOutFileCentroids());

	if (  linparam_ClusteringGA.getPrintCentroidsFormat() ) {

	  std::ostringstream lostrstream_labelNameCol;
	
	  /*WITH HEADER*/
	  if (linparam_ClusteringGA.getHaveHeaderFileInstance()) {
	    for (uintidx li_l = 0; li_l < lpairvecstrstr_instanceDimName.first.size()-1; li_l++) {
	      lostrstream_labelNameCol
		<< lpairvecstrstr_instanceDimName.first.at(li_l)
		<< lofn_filenamecentroids.getDelim();  
	    }
	    lostrstream_labelNameCol 
	      << lpairvecstrstr_instanceDimName.first.at(lpairvecstrstr_instanceDimName.first.size()-1); 
	  }
	
	  lomatrixrowt_centroids.r8mat_print
	    ("Centroids:",
	     lostrstream_labelNameCol.str().c_str(),
	     inout::OutFileName::getDelim(),
	     5,
	     lostream_outcentroids
	     );
	
	  lostream_outcentroids << std::endl;
	
	}
      
	else {
		
	  std::ostringstream lostrstream_labelCentroids;

	  lostream_outcentroids.precision(COMMON_COUT_PRECISION);
	  lostream_outcentroids << std::boolalpha;

#ifdef __MEDOIDS_ALGORITHM
	  lostrstream_labelCentroids << "<MEDOIDS"; 
#else 	   
	  lostrstream_labelCentroids << "<CENTROIDS"; 
#endif  /* __MEDOIDS_ALGORITHM */  

	  lostrstream_labelCentroids
	    << lofn_filenamecentroids.getDelim() <<  "_algorithmo" 
	    << lofn_filenamecentroids.getDelim()
	    << linparam_ClusteringGA.getAlgorithmoName();
	  lostrstream_labelCentroids
	    << lofn_filenamecentroids.getDelim() << "_author" 
	    << lofn_filenamecentroids.getDelim()
	    << linparam_ClusteringGA.getAlgorithmoAuthor();
	  
	  lostrstream_labelCentroids
	    << lofn_filenamecentroids.getDelim() 
	    << "_runnig date" << lofn_filenamecentroids.getDelim() 
	    << linparam_ClusteringGA.getRunningTimeId();
	  lostrstream_labelCentroids
	    << lofn_filenamecentroids.getDelim() 
	    << "_number run" << lofn_filenamecentroids.getDelim() 
	    << loop_outParamGAC.getNumRunningAlgorithm();
	  lostrstream_labelCentroids
	    << lofn_filenamecentroids.getDelim() 
	    << "_times run"  << lofn_filenamecentroids.getDelim() 
	    << linparam_ClusteringGA.getTimesRunAlgorithm();
	  
	  lostrstream_labelCentroids
	    << lofn_filenamecentroids.getDelim()
	    << "_dataset_training"
	    << lofn_filenamecentroids.getDelim()
	    << linparam_ClusteringGA.getCurrentFileInstance();
      
	  if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	
	    lostrstream_labelCentroids
	      << lofn_filenamecentroids.getDelim()
	      << "_dataset_testing"
	      << lofn_filenamecentroids.getDelim()
	      << linparam_ClusteringGA.getCurrentFileInstanceTest();
	  }

	  /*WITH HEADER*/
	  if (linparam_ClusteringGA.getHaveHeaderFileInstance()) {
	    lostrstream_labelCentroids << lofn_filenamecentroids.getDelim();
	    for (uintidx li_l = 0; li_l < lpairvecstrstr_instanceDimName.first.size()-1; li_l++) {
	      lostrstream_labelCentroids
		<< lpairvecstrstr_instanceDimName.first.at(li_l)
		<< lofn_filenamecentroids.getDelim();  
	    }
	    lostrstream_labelCentroids 
	      << lpairvecstrstr_instanceDimName.first.at(lpairvecstrstr_instanceDimName.first.size()-1); 
	  }
		  
	  lomatrixrowt_centroids.print
	    (lostream_outcentroids,
	     lostrstream_labelCentroids.str().c_str(),
	     ',',
	     ';'
	     );
	
	  lostream_outcentroids << std::endl;

#ifdef ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001

	  std::ostringstream lostrstream_labelVi;
	  lostrstream_labelVi << "<Vi";

	  lostrstream_labelVi << lofn_filenamecentroids.getDelim() <<  "_algorithmo" 
			      << lofn_filenamecentroids.getDelim()
			      << linparam_ClusteringGA.getAlgorithmoName();
	  lostrstream_labelVi << lofn_filenamecentroids.getDelim() << "_author" 
			      << lofn_filenamecentroids.getDelim()
			      << linparam_ClusteringGA.getAlgorithmoAuthor();
	  
	  lostrstream_labelVi << lofn_filenamecentroids.getDelim() 
			      << "_runnig date" << lofn_filenamecentroids.getDelim() 
			      << linparam_ClusteringGA.getRunningTimeId();
	  lostrstream_labelVi << lofn_filenamecentroids.getDelim() 
			      << "_number run" << lofn_filenamecentroids.getDelim() 
			      << loop_outParamGAC.getNumRunningAlgorithm();
	  lostrstream_labelVi << lofn_filenamecentroids.getDelim() 
			      << "_times run"  << lofn_filenamecentroids.getDelim() 
			      << linparam_ClusteringGA.getTimesRunAlgorithm();
	  
	  lostrstream_labelVi << lofn_filenamecentroids.getDelim()
			      << "_dataset_training"
			      << lofn_filenamecentroids.getDelim()
			      << linparam_ClusteringGA.getCurrentFileInstance();
	
	  lmatrixrowt_Vi.print
	    (lostream_outcentroids,
	     lostrstream_labelVi.str().c_str(),
	     ',',
	     ';'
	     );
	  lostream_outcentroids << std::endl;
	
#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/
	
	}

	lofn_filenamecentroids.closeFile();

      } /*END PRINT CENTROIDES*/

      /*BEGIN PRINT MEMBERSHIP*/
      if ( linparam_ClusteringGA.getOutFileMemberShip() != NULL )  { 

	inout::OutFileName lofn_filenamemembership; 
	std::ostream&      lostream_outcmembership = 
	  lofn_filenamemembership.openFile(linparam_ClusteringGA.getOutFileMemberShip());

	std::ostringstream lostrstream_labelMemberTraining;
	lostrstream_labelMemberTraining << "<MEMBERCLUSTERTRAINING";

	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() << "_algorithmo" 
	  << lofn_filenamemembership.getDelim()
	  << linparam_ClusteringGA.getAlgorithmoName();
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() << "_author" 
	  << lofn_filenamemembership.getDelim()
	  << linparam_ClusteringGA.getAlgorithmoAuthor(); 
    
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_runnig date" << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_number run" << lofn_filenamemembership.getDelim() 
	  << loop_outParamGAC.getNumRunningAlgorithm();
	  
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_runnig date" << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_number run" << lofn_filenamemembership.getDelim() 
	  << loop_outParamGAC.getNumRunningAlgorithm();
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_times run"  << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getTimesRunAlgorithm();
	  
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim()
	  << "_dataset_training"
	  << lofn_filenamemembership.getDelim()
	  << linparam_ClusteringGA.getCurrentFileInstance();
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	  lostrstream_labelMemberTraining
	    << lofn_filenamemembership.getDelim()
	    << "_dataset_testing"
	    << lofn_filenamemembership.getDelim()
	    << linparam_ClusteringGA.getCurrentFileInstanceTest();
	  
	}
	  
	lpartition_clusters.print
	  (lostream_outcmembership,
	   lostrstream_labelMemberTraining.str().c_str(),
	   ','
	   );

#ifdef ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001

	std::ostringstream lostrstream_labelMemberBi;
	lostrstream_labelMemberBi << "<Bi"
				  << "TRAINING";

	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() << "_algorithmo" 
	  << lofn_filenamemembership.getDelim()
	  << linparam_ClusteringGA.getAlgorithmoName();
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() << "_author" 
	  << lofn_filenamemembership.getDelim()
	  << linparam_ClusteringGA.getAlgorithmoAuthor(); 
    
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_runnig date" << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_number run" << lofn_filenamemembership.getDelim() 
	  << loop_outParamGAC.getNumRunningAlgorithm();
	  
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_runnig date" << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_number run" << lofn_filenamemembership.getDelim() 
	  << loop_outParamGAC.getNumRunningAlgorithm();
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_times run"  << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getTimesRunAlgorithm();
	  
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim()
	  << "_dataset_training"
	  << lofn_filenamemembership.getDelim()
	  << linparam_ClusteringGA.getCurrentFileInstance();
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	  lostrstream_labelMemberBi
	    << lofn_filenamemembership.getDelim()
	    << "_dataset_testing"
	    << lofn_filenamemembership.getDelim()
	    << linparam_ClusteringGA.getCurrentFileInstanceTest();
	  
	}

	lmembclassdisjsets_Bi.print
	  (lostream_outcmembership,
	   lostrstream_labelMemberBi.str().c_str(),
	   ','
	   );
      

#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/
      
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {

	  std::ostringstream lostrstream_labelMemberTest;
	  lostrstream_labelMemberTest << "<MEMBERCLUSTERTEST";

	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim() << "_algorithmo" 
	    << lofn_filenamemembership.getDelim() << linparam_ClusteringGA.getAlgorithmoName();
	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim() << "_author" 
	    << lofn_filenamemembership.getDelim() << linparam_ClusteringGA.getAlgorithmoAuthor();
	  
	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim() 
	    << "_runnig date" << lofn_filenamemembership.getDelim() 
	    << linparam_ClusteringGA.getRunningTimeId();
	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim() 
	    << "_number run" << lofn_filenamemembership.getDelim() 
	    << loop_outParamGAC.getNumRunningAlgorithm();
	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim() 
	    << "_times run"  << lofn_filenamemembership.getDelim() 
	    << linparam_ClusteringGA.getTimesRunAlgorithm();
	  
	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim()
	    << "_dataset_training"
	    << lofn_filenamemembership.getDelim()
	    << linparam_ClusteringGA.getCurrentFileInstance();
	 
	  lostrstream_labelMemberTest
	    << lofn_filenamemembership.getDelim()
	    << "_dataset_testing"
	    << lofn_filenamemembership.getDelim()
	    << linparam_ClusteringGA.getCurrentFileInstanceTest();
	  lpartition_clustersTest.print
	    (lostream_outcmembership,
	     lostrstream_labelMemberTest.str().c_str(),
	     ','
	     );
	 
	}
	  	  
	lofn_filenamemembership.closeFile();
      } /*END PRINT MEMBERSHIP*/

      /*BEGIN PRINT TABLE OF PARTITION*/
      if ( linparam_ClusteringGA.getOutFilePartitionsTable() != NULL )  {
      
	inout::OutFileName  lofn_filenametablepartition;
	std::ostream& lostream_outtablepartition = 
	  lofn_filenametablepartition.openFile( linparam_ClusteringGA.getOutFilePartitionsTable()  );

	if (  linparam_ClusteringGA.getPrintTableFormat() ) {

	  lmatchmatrix_confusion.r8mat_print
	    ("Partition table:",
	     5,
	     lostream_outtablepartition
	     );
	
	}
      
	else {
	
	  std::ostringstream lostrstream_labelPartitionsTable;
	  lostrstream_labelPartitionsTable //lostream_out
	    << "<PARTITIONTABLE" //artitions table" 
	    << lofn_filenametablepartition.getDelim() <<  "_algorithmo" 
	    << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getAlgorithmoName()
	    << lofn_filenametablepartition.getDelim() << "_author" 
	    << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getAlgorithmoAuthor()
	    << lofn_filenametablepartition.getDelim() << "runnig date"
	    << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getRunningTimeId()
	    << lofn_filenametablepartition.getDelim() << "number run"
	    << lofn_filenametablepartition.getDelim()
	    << loop_outParamGAC.getNumRunningAlgorithm()
	    << lofn_filenametablepartition.getDelim() << "times run"
	    << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getTimesRunAlgorithm();
	
	  lmatchmatrix_confusion.print
	    (lostream_outtablepartition,
	     lostrstream_labelPartitionsTable.str().c_str(),
	     ',',
	     ';'
	     );
	  lostream_outtablepartition << std::endl;
	}
        /*lostream_outparam << "\n";
	lostream_outparam << "               Rand index: "
			  << sm::randIndex<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusion);

	lostream_outparam << '\n';
	lostream_outparam << "                   Purity: "
			  << sm::purity<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusion);

	lostream_outparam << '\n';
	lostream_outparam << "                Precision: "
			  << sm::precision<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusion);

	lostream_outparam << '\n';
	lostream_outparam << "                   Recall: "
			  << sm::recall<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusion);
	
	lostream_outparam << '\n' << std::endl;
	*/
	  	  
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {

	  if (  linparam_ClusteringGA.getPrintTableFormat() ) {

	    lmatchmatrix_confusionTest.r8mat_print
	      ("Partition table test:",
	       5,
	       lostream_outtablepartition
	       );
	
	  }
      
	  else {
	  
	    std::ostringstream lostrstream_labelPartitionsTableTest;
	    lostrstream_labelPartitionsTableTest //lostream_out
	      << "<PARTITIONTABLETEST" //artitions table" 
	      << lofn_filenametablepartition.getDelim() <<  "_algorithmo" 
	      << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getAlgorithmoName()
	      << lofn_filenametablepartition.getDelim() << "_author" 
	      << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getAlgorithmoAuthor()
	      << lofn_filenametablepartition.getDelim() << "runnig date"
	      << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getRunningTimeId()
	      << lofn_filenametablepartition.getDelim() << "number run"
	      << lofn_filenametablepartition.getDelim() << loop_outParamGAC.getNumRunningAlgorithm()
	      << lofn_filenametablepartition.getDelim() << "times run"
	      << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getTimesRunAlgorithm();
	
	    lmatchmatrix_confusionTest.print
	      (lostream_outtablepartition,
	       lostrstream_labelPartitionsTableTest.str().c_str(),
	       ',',
	       ';'
	       );

	    lostream_outtablepartition << std::endl;
	  }

	  /* lostream_outparam << '\n';
	  lostream_outparam << "    Test data Rand index: "
			    << sm::randIndex<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusionTest);

	  lostream_outparam << '\n';
	  lostream_outparam << "        Test data Purity: "
			    << sm::purity<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusionTest);

	  lostream_outparam << '\n';
	  lostream_outparam << "     Test data Precision: "
			    << sm::precision<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusionTest);

	  lostream_outparam << '\n';
	  lostream_outparam << "        Test data Recall: "
			    << sm::recall<DATATYPE_REAL, DATATYPE_INSTANCES_CLUSTER_K>(lmatchmatrix_confusionTest);
	
	  lostream_outparam << std::endl;
	  */
	} /* Test*/
      
	lofn_filenametablepartition.closeFile();
      } /*END PRINT TABLE OF PARTITION*/

      //lofn_filenameparam.closeFile();

      /*BEGIN PRINT GRAPH*/
      if ( linparam_ClusteringGA.getOutFileGraph() != NULL )  { 

	inout::OutFileName   lofn_filenamegraph;
	std::ostream&        lostream_outgraph = 
	  lofn_filenamegraph.openFile(linparam_ClusteringGA.getOutFileGraph());

	std::ostringstream lostrstream_labelGraph;
	lostream_outgraph << "<GRAPH";

	lostrstream_labelGraph
	  << lofn_filenamegraph.getDelim() <<  "_algorithmo" 
	  << lofn_filenamegraph.getDelim()
	  << linparam_ClusteringGA.getAlgorithmoName();
	lostrstream_labelGraph
	  << lofn_filenamegraph.getDelim() << "_author" 
	  << lofn_filenamegraph.getDelim()
	  << linparam_ClusteringGA.getAlgorithmoAuthor();
	  
	lostrstream_labelGraph
	  << lofn_filenamegraph.getDelim() 
	  << "_runnig date" << lofn_filenamegraph.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelGraph
	  << lofn_filenamegraph.getDelim() 
	  << "_number run" << lofn_filenamegraph.getDelim() 
	  << loop_outParamGAC.getNumRunningAlgorithm();
	lostrstream_labelGraph
	  << lofn_filenamegraph.getDelim() 
	  << "_times run"  << lofn_filenamegraph.getDelim() 
	  << linparam_ClusteringGA.getTimesRunAlgorithm();
	  
	lostrstream_labelGraph
	  << lofn_filenamegraph.getDelim()
	  << "_dataset_training"
	  << lofn_filenamegraph.getDelim()
	  << linparam_ClusteringGA.getCurrentFileInstance();
      
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	
	  lostrstream_labelGraph
	    << lofn_filenamegraph.getDelim()
	    << "_dataset_testing"
	    << lofn_filenamegraph.getDelim()
	    << linparam_ClusteringGA.getCurrentFileInstanceTest();
	}

#ifdef ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003
      
	auto &&lvectorlist_tmpgraphpi =
	  graph::fromPiToGraph 
	  (lpair_chrombestPiMST.second,
	   lpair_chrombestPiMST.first
	   );
    
	inout::vectorlistprint
	  (lvectorlist_tmpgraphpi,
	   lostream_outgraph,
	   lostrstream_labelGraph.str().c_str(),
	   ';'
	   );
      
	lostream_outgraph  << std::endl;

#endif /*ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003*/
	
	lofn_filenamegraph.closeFile();
      
      } /*BEGIN PRINT GRAPH*/

    
      /*BAR PRINTING*/
      if ( linparam_ClusteringGA.getProgressBarPrinting() ) 
	inout::barprogress_update
	  (__lst_i * linparam_ClusteringGA.getTimesRunAlgorithm()+ li_l, 
	   linparam_ClusteringGA.getTimesRunAlgorithm() * 
	   linparam_ClusteringGA.getNumFilesInstance()
	   );

    } /*FOR TIMES RUN ALGORITHM*/

    /*FREE MEMORY
      DELETE INSTANCES
    */
    delete [] larray_meanFeactures;
    delete [] larray_desvstdFeactures;
    
    for ( auto  liter_instance: lpairvec_dataset.first ) 
      delete liter_instance;

    for ( auto  liter_instanceTest: lpairvec_dataset.second ) 
      delete liter_instanceTest;

    delete pfunct2p_distAlg;
    delete pfunct2p_distEuclidean;
    delete pfunct2p_distEuclideanSq;
  
  } /*END for NumFilesInstance */

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labeMain 
	      << ": OUT(" << geiinparam_verbose << ')'
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    
  return 0;
 
} /*END MAIN*/

