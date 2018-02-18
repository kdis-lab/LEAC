/*! \file datatype_instance_real.hpp
 *
 * \brief datatype instance real
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef DATA_TYPE_INSTANCE_REAL_HPP
#define DATA_TYPE_INSTANCE_REAL_HPP

/* For better performance, the data type of the instances 
   must be equal to the centroids
 */
#define DATATYPE_FEATURE              double 
#define DATATYPE_FEATURE_SUM          double

/* CENTROIDS ROUND :
   if data for floating point to integers must be  
       DATATYPE_CENTROIDS_ROUND 1
   else data all floating
       DATATYPE_CENTROIDS_ROUND 0
*/
#define DATATYPE_CENTROIDS_ROUND      0

#define DATATYPE_INSTANCES_CLUSTER_K  long

#define DATATYPE_CLUSTERIDX           int    // -1, 0, .., K  

/*When instances consist of occurrences 
  or frequency
 */
#define DATATYPE_INSTANCE_FREQUENCY   int  

/*￼Data type for metrics and distances
 */
#define DATATYPE_REAL                 double

/*Tipo de dato para almacenar bits
*/
#define DATATYPE_BITSIZE unsigned int

#endif /*DATA_TYPE_INSTANCE_REAL_HPP*/

