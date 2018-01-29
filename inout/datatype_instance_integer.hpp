/*! \file datatype_instance_integer.hpp
 *
 * \brief datatype instance integer
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef DATA_TYPE_INSTANCE_INTEGER_HPP
#define DATA_TYPE_INSTANCE_INTEGER_HPP


#define DATATYPE_FEATURE               int
#define DATATYPE_FEATURE_SUM           long

/* CENTROIDS ROUND :
   if data for floating point to integers must be  
       DATATYPE_CENTROIDS_ROUND 1
   else data all floating
       DATATYPE_CENTROIDS_ROUND 0
*/

#define DATATYPE_CENTROIDS_ROUND         1

#define DATATYPE_INSTANCES_CLUSTER_K    long

#define DATATYPE_CLUSTERIDX             int    // -1, 0, .., K  

/*When instances consist of occurrences 
  or frequency
 */
#define DATATYPE_INSTANCE_FREQUENCY      int  

/*￼Data type for metrics and distances
 */
#define DATATYPE_REAL           double

/*Tipo de dato para almacenar bits
*/
#define DATATYPE_BITSIZE unsigned int

#endif /*DATA_TYPE_INSTANCE_INTEGER_HPP*/

