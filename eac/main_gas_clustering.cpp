/*! \file main_gas_clustering.cpp
 *
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

#include <leac.hpp>
#include "inparamclustering_getparameter.hpp"
#include "instances_read.hpp"
//#include "instances_classfrequency_read.hpp"
#include "bar_progress.hpp"

 
/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) 
{

#ifdef __VERBOSE_YES
  const char* lpc_labeMain = "main_gas_clustering";
  geverbosepc_labelstep = lpc_labeMain;
#endif /*__VERBOSE_YES*/
  
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSE);
  
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

  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::Distortion);
    
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

  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::TWCV);

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

  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::TWCV);
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

  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::TWCV);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSE);
    
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSE);
    
#endif /*ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002*/


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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamEAC(inout::DBindex);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamEAC(inout::VRC);
 
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>      loop_outParamEAC(inout::IntraInterClust);
  

#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/

#ifdef ALG_CGA_VKLABEL_HRUSCHKA_EBECKEN_2003   
  /*INPUT: PARAMETER*/
  inout::InParamPcPmVk
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
  
  linparam_ClusteringGA.setNumMaxGenerations(10000);
  
  /*Populations formed by 20 genotypes
   */
  linparam_ClusteringGA.setSizePopulation(20);
  linparam_ClusteringGA.setProbCrossover(0.5);
  linparam_ClusteringGA.setProbMutation(0.25);
  
  /*OUT:
   */
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::Silhouette); 
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSilhouette);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSilhouette);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSilhouette);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSilhouette);
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
    ("FEAC",
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSilhouette);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX
     >
    loop_outParamEAC(inout::VRC);
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
     "Carro-Calvo L. and Del Ser J. and Portilla-Figueras, J.A.",
     inout::LABEL,
     INPARAMCLUSTERING_DISTANCE_EUCLIDEAN
     );
  /*OUT
   */
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::DBindex);
  
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX
     >
    loop_outParamEAC(inout::Silhouette);
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>   
    loop_outParamEAC(inout::TWCV);    
#endif /*ALG_GAGR_FKCENTROID_CHANG_ETAL_2009*/


#ifdef ALG_GCA_FKMEDOID_LUCASIUS_ETAL1993 
  /*INPUT: PARAMETER
   */
  inout::InParamGCA
    <DATATYPE_CLUSTERIDX,
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSE);
  
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::J1);
#endif /*ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997*/
  

#ifdef ALG_HKA_FKMEDOID_SHENG_LUI2004
  /*INPUT: PARAMETER
   */
  inout::InParamHKA
    <DATATYPE_CLUSTERIDX,
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
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX>
    loop_outParamEAC(inout::SSE);
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

  linparam_ClusteringGA.setNumMaxGenerations(200);
  linparam_ClusteringGA.setSizePopulation(200);
  linparam_ClusteringGA.setSizeMatingPool(100);

  /*OUT
   */
  inout::OutParamEAClustering
    <DATATYPE_REAL,
     DATATYPE_CLUSTERIDX> 
    loop_outParamEAC(inout::J1);
#endif /*ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994*/
  
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

    std::vector<std::string> 
      lvectorstr_instanceDimName = 
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

    /*DISTANCE
     */
    dist::Dist<DATATYPE_REAL,DATATYPE_FEATURE>
      *pfunct2p_distAlg = NULL;
    dist::Dist<DATATYPE_REAL,DATATYPE_FEATURE>
      *pfunct2p_distEuclidean =
      new dist::Euclidean<DATATYPE_REAL,DATATYPE_FEATURE>();
    dist::Dist<DATATYPE_REAL,DATATYPE_FEATURE>
      *pfunct2p_distEuclideanSq =
      new dist::EuclideanSquared<DATATYPE_REAL,DATATYPE_FEATURE>();
    
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

      loop_outParamEAC.initialize(li_l);
    
      { /*BEGIN PRINT PARAMETERS
	 */
	inout::OutFileName   lofn_filename;
	std::ostream& lostream_out = 
	  lofn_filename.openFile(linparam_ClusteringGA.getFileNameTimesRun());
	lostream_out
	  << "\nIN:\n"
	  << "    Algorithmo name: "
	  << linparam_ClusteringGA.getAlgorithmoName() << '\n'
	  << "           Based on: "
	  << linparam_ClusteringGA.getAlgorithmoAuthor() << '\n'
	  << "        Metric used: "
	  << loop_outParamEAC.getNameUsedObjetiveFunc()
	  << "\n\n"
	  
	  << "           Data set: "
	  << linparam_ClusteringGA.getCurrentFileInstance() << '\n'
	  << "Number of instances: "
	  << linparam_ClusteringGA.getNumInstances() << '\n'
	  << "         Dimensions: "
	  << linparam_ClusteringGA.getNumDimensionsInstances() << '\n';

	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {
	  lostream_out
	    << '\n'
	    << "      Data set test: "
	    << linparam_ClusteringGA.getCurrentFileInstanceTest() << '\n'
	    << "Number of instances: "
	    << linparam_ClusteringGA.getNumInstancesTest() << '\n';
	}
	lostream_out
	  << '\n'
	  << "       Random seed: "
	  <<  linparam_ClusteringGA.getRandomSeed() << '\n'
	  << '\n';
	lofn_filename.closeFile();
      
      } /*END PRINT PARAMETERS
	 */
    

#ifdef ALG_GACLUSTERING_FKCRISPMATRIX_BEZDEK_ETAL_1994

      gaencode::ChromosomeCrispMatrix
	<DATATYPE_BITSIZE,DATATYPE_CLUSTERIDX,DATATYPE_REAL>&&
	lchrom_best =  
	eac::gaclustering_fkcrispmatrix
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
    
      /*\cite{Bandyopadhyay:Maulik:GAclustering:KGA:2002}
       */
#ifdef ALG_KGA_FKCENTROID_BANDYOPADHYAY_MAULIK_2002 
      gaencode::ChromFixedLength<DATATYPE_FEATURE,DATATYPE_REAL>&& lchrom_best = 
	eac::kga_fkcentroid
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	 loop_outParamEAC.getNumClusterK()
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) loop_outParamEAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) loop_outParamEAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) loop_outParamEAC.getNumClusterK(),
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

      auto lchrom_best =
	eac::feca_vklabel
	(loop_outParamEAC,
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
	 loop_outParamEAC.getNumClusterK()
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
	(loop_outParamEAC,
	 linparam_ClusteringGA,
	 lpairvec_dataset.first.begin(),
	 lpairvec_dataset.first.end(),
	 *pfunct2p_distAlg
	 );

      if ( loop_outParamEAC.getEndingCondition() ) {
      
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
	 loop_outParamEAC.getNumClusterK()
	 );

#endif /*ALG_CLUSTERING_VKSUBCLUSTERBINARY_TSENG_YANG_2001*/

   
#ifdef ALG_GA_CLUSTERING_VKTREEBINARY_CASILLAS_GONZALEZ_MARTINEZ_2003

      auto lpair_chrombestPiMST =
	eac::gaclustering_vktreebinary
	(loop_outParamEAC,
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
	((uintidx) loop_outParamEAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) loop_outParamEAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) loop_outParamEAC.getNumClusterK(),
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

    
#ifdef ALG_GCUK_VKCENTROID_BANDYOPADHYAY_AND_MAULIK_2002

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids =
	eac::gcuk_vkcentroid
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	 loop_outParamEAC.getNumClusterK()  
	 );

      mat::MatrixRow<DATATYPE_FEATURE> 
	lomatrixrowt_centroids
	((uintidx) loop_outParamEAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions() 
	 );

      mat::MatrixRow<DATATYPE_FEATURE_SUM>       
	lomatrixrowt_sumInstCluster
	((uintidx) loop_outParamEAC.getNumClusterK(), 
	 data::Instance<DATATYPE_FEATURE>::getNumDimensions(),
	 DATATYPE_FEATURE_SUM(0)
	 );
	
      std::vector<DATATYPE_INSTANCES_CLUSTER_K> 
	lovectort_numInstClusterK
	((uintidx) loop_outParamEAC.getNumClusterK(),
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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
	(loop_outParamEAC,
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

#endif /*ALG_HKA_FKMEDOID_SHENG_LUI2004*/

#ifdef ALG_GAPROTOTYPES_FKMEDOID_KUNCHEVA_BEZDEK_1997

      gaencode::ChromosomeBitArray<DATATYPE_BITSIZE,DATATYPE_REAL>&& lchrom_best =
	eac::gaprototypes_fkmedoid
	(loop_outParamEAC,
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

	if (lpartition_clusters.getNumCluster() > 1 ) {
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

	} 

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
	    loop_outParamEAC.setNumClusterKTest
	      (DATATYPE_CLUSTERIDX(lvectort_numInstClusterKTest.size())
	       - DATATYPE_CLUSTERIDX(li_matrixPartitionClusterNullTest)
	       );
	    if ( li_matrixPartitionClusterNullTest != 0 )  {
	      loop_outParamEAC.setEndingConditionTest(false);
	    }
	    else  {
	      loop_outParamEAC.setEndingConditionTest(true);
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

      } /*END  METRIC TEST*/
    
      /*END TEST METRIC-------------------------------------------------------
       */

      /*OUT: PRINT------------------------------------------------------------
       */

      inout::OutFileName lofn_filenameparam;
      std::ostream& lostream_outparam = 
	lofn_filenameparam.openFile(linparam_ClusteringGA.getFileNameTimesRun());

    
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
      
	lostream_outparam << "OUT:\n\n";

	lchrom_best.print(lostream_outparam,"CROMOSOME:BEST",',',';');
		      
	lostream_outparam << "\n\n";

	lostream_outparam << "       Cluster number (K): "
			  << loop_outParamEAC.getNumClusterK() << "\n"
			  <<  std::right << std::setw(25)
			  <<  loop_outParamEAC.getNameUsedObjetiveFunc() << ": "
			  << loop_outParamEAC.getObjetiveFuncRun() << "\n";

	if ( loop_outParamEAC.getUsedObjetiveFunc() != inout::SSE ) {
	  std::pair<DATATYPE_REAL,bool> lpair_sse =
	    um::SSE
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     *pfunct2p_distEuclidean
	     );
	  lostream_outparam << "                      SSE: "
			    << lpair_sse.first;
	  if ( lpair_sse.second == false ) 
	    lostream_outparam << " Has group without objects";
	  lostream_outparam << '\n';
	}

	if ( loop_outParamEAC.getUsedObjetiveFunc() != inout::DBindex ) {
	  lostream_outparam
	    << "                 DB-index: "
	    <<
	    um::dbindex
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),
	     lpartition_clusters,
	     *pfunct2p_distEuclidean
	     )
	    << '\n';
	}

	if ( loop_outParamEAC.getUsedObjetiveFunc() != inout::Silhouette ) {
	  lostream_outparam
	    << "               Silhouette: "	
	    <<
	    um::silhouette
	    (lpairvec_dataset.first.begin(),
	     lpartlinknuminst_memberShip,
	     lpartlinknuminst_memberShip.getVectorNumInstClusterK(),
	     *pfunct2p_distEuclidean
	     )
	    << '\n';
	}

	if ( loop_outParamEAC.getUsedObjetiveFunc() != inout::VRC ) {
	  lostream_outparam
	    << "                      VRC: "	
	    <<
	    um::VRC
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.first.begin(),
	     lpairvec_dataset.first.end(),   
	     lpartition_clusters,
	     *pfunct2p_distEuclideanSq
	     )
	    << '\n';
	}
      
	if ( loop_outParamEAC.getUsedObjetiveFunc() != inout::CSmeasure ) {
	  lostream_outparam
	    << "               CS measure: "
	    <<
	    um::CSmeasure
	    (lpairvec_dataset.first.begin(),
	     lomatrixrowt_centroids,
	     lpartlinknuminst_memberShip,
	     *pfunct2p_distEuclidean
	     )
	    << '\n';
	}
      
	if ( loop_outParamEAC.getUsedObjetiveFunc() != inout::DunnIndex ) {
	  lostream_outparam
	    << "             Dunn's index: "
	    <<
	    um::DunnIndex
	    (lpairvec_dataset.first.begin(),
	     lpartlinknuminst_memberShip,
	     *pfunct2p_distEuclidean
	     )
	    << '\n';
	}
    
	if ( linparam_ClusteringGA.getNumFilesInstanceTest() > 0 ) {

	  ds::PartitionLinkedNumInst
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
	 
	  std::pair<DATATYPE_REAL,bool> lpair_sseTest =
	    um::SSE
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),
	     *pfunct2p_distEuclidean
	     );
	  lostream_outparam << '\n';
	  lostream_outparam
	    << "            Test data SSE: "
	    << lpair_sseTest.first;
	  if ( lpair_sseTest.second == false ) 
	    lostream_outparam << " Has group without objects";
	  lostream_outparam << '\n';

	  lostream_outparam
	    << "       Test data DB-index: "
	    <<
	    um::dbindex
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),
	     lpartition_clusters,
	     *pfunct2p_distEuclidean
	     )
	    << '\n';

	  lostream_outparam
	    << "     Test data Silhouette: "	
	    <<
	    um::silhouette
	    (lpairvec_dataset.second.begin(),
	     lpartlinknuminst_memberShipTest,
	     lpartlinknuminst_memberShipTest.getVectorNumInstClusterK(),
	     *pfunct2p_distEuclidean
	     )
	    << '\n';

	  lostream_outparam
	    << "            Test data VRC: "	
	    <<
	    um::VRC
	    (lomatrixrowt_centroids,
	     lpairvec_dataset.second.begin(),
	     lpairvec_dataset.second.end(),   
	     lpartition_clusters,
	     *pfunct2p_distEuclideanSq
	     )
	    << '\n';

	  lostream_outparam
	    << "     Test data CS measure: "
	    <<
	    um::CSmeasure
	    (lpairvec_dataset.second.begin(),
	     lomatrixrowt_centroids,
	     lpartlinknuminst_memberShipTest,
	     *pfunct2p_distEuclidean
	     )
	    << '\n';
     
      
	  lostream_outparam
	    << "   Test data Dunn's index: "
	    <<
	    um::DunnIndex
	    (lpairvec_dataset.second.begin(),
	     lpartlinknuminst_memberShipTest,
	     *pfunct2p_distEuclidean
	     )
	    << '\n';
	}
      
	lostream_outparam << "     Execution time (seg): "
			  << loop_outParamEAC.getAlgorithmRunTime() << '\n'
			  << "Generations find the best: "
			  << loop_outParamEAC.getIterationGetsBest() << '\n';
    
	lostream_outparam << std::endl;

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
	    for (uintidx li_l = 0; li_l < lvectorstr_instanceDimName.size()-1; li_l++) {
	      lostrstream_labelNameCol
		<< lvectorstr_instanceDimName.at(li_l)
		<< lofn_filenamecentroids.getDelim();  
	    }
	    lostrstream_labelNameCol 
	      << lvectorstr_instanceDimName.at(lvectorstr_instanceDimName.size()-1); 
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
	    << loop_outParamEAC.getNumRunningAlgorithm();
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
	      << linparam_ClusteringGA.getCurrentFileInstance();
	  }

	  /*WITH HEADER*/
	  if (linparam_ClusteringGA.getHaveHeaderFileInstance()) {
	    lostrstream_labelCentroids << lofn_filenamecentroids.getDelim();
	    for (uintidx li_l = 0; li_l < lvectorstr_instanceDimName.size()-1; li_l++) {
	      lostrstream_labelCentroids
		<< lvectorstr_instanceDimName.at(li_l)
		<< lofn_filenamecentroids.getDelim();  
	    }
	    lostrstream_labelCentroids 
	      << lvectorstr_instanceDimName.at(lvectorstr_instanceDimName.size()-1); 
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
			      << loop_outParamEAC.getNumRunningAlgorithm();
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
	  << loop_outParamEAC.getNumRunningAlgorithm();
	  
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_runnig date" << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelMemberTraining
	  << lofn_filenamemembership.getDelim() 
	  << "_number run" << lofn_filenamemembership.getDelim() 
	  << loop_outParamEAC.getNumRunningAlgorithm();
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
	    << linparam_ClusteringGA.getCurrentFileInstance();
	  
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
	  << loop_outParamEAC.getNumRunningAlgorithm();
	  
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_runnig date" << lofn_filenamemembership.getDelim() 
	  << linparam_ClusteringGA.getRunningTimeId();
	lostrstream_labelMemberBi
	  << lofn_filenamemembership.getDelim() 
	  << "_number run" << lofn_filenamemembership.getDelim() 
	  << loop_outParamEAC.getNumRunningAlgorithm();
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
	    << linparam_ClusteringGA.getCurrentFileInstance();
	  
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
	    << loop_outParamEAC.getNumRunningAlgorithm();
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
	    << linparam_ClusteringGA.getCurrentFileInstance();
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
	    << loop_outParamEAC.getNumRunningAlgorithm()
	    << lofn_filenametablepartition.getDelim() << "times run"
	    << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getTimesRunAlgorithm();
	
	  lmatchmatrix_confusion.print
	    (lostream_outtablepartition,
	     lostrstream_labelPartitionsTable.str().c_str(),
	     ',',
	     ';'
	     );	
	}
        lostream_outparam << "\n";
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
	      << lofn_filenametablepartition.getDelim() << loop_outParamEAC.getNumRunningAlgorithm()
	      << lofn_filenametablepartition.getDelim() << "times run"
	      << lofn_filenametablepartition.getDelim() << linparam_ClusteringGA.getTimesRunAlgorithm();
	
	    lmatchmatrix_confusionTest.print
	      (lostream_outtablepartition,
	       lostrstream_labelPartitionsTableTest.str().c_str(),
	       ',',
	       ';'
	       );
	  }

	  lostream_outparam << '\n';
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
	} /* Test*/
      
	lofn_filenametablepartition.closeFile();
      } /*END PRINT TABLE OF PARTITION*/

      lofn_filenameparam.closeFile();

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
	  << loop_outParamEAC.getNumRunningAlgorithm();
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
	    << linparam_ClusteringGA.getCurrentFileInstance();
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

