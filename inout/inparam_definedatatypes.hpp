/*! \file inparam_definedatatypes.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_DEFINE_DATA_TYPE_HPP
#define __IN_PARAM_DEFINE_DATA_TYPE_HPP


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamDefineFeatFeatSumInstK
    \brief Define data types for the input parameters of the templates
*/
template <typename T_BITSIZE,
          typename T_FEATURE,         
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K 
	  > 
class InParamDefineBitSizeFeatFeatSumInstK {
public:
  InParamDefineBitSizeFeatFeatSumInstK() {}
  ~InParamDefineBitSizeFeatFeatSumInstK() {}
};
  
/*! \class InParamDefineFeatFeatSumInstK
    \brief Define data types for the input parameters of the templates
*/
template <typename T_FEATURE,         
	  typename T_FEATURE_SUM,
	  typename T_INSTANCES_CLUSTER_K 
	  > 
class InParamDefineFeatFeatSumInstK {
public:
  InParamDefineFeatFeatSumInstK() {}
  ~InParamDefineFeatFeatSumInstK() {}
};

/*! \class InParamDefineFeatInstK
    \brief Define data types for the input parameters of the templates
*/
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K 
	  > 
class InParamDefineFeatInstK {
public:
  InParamDefineFeatInstK() {}
  ~InParamDefineFeatInstK() {}
};

/*! \class InParamDefineFeat
    \brief Define data types for the input parameters of the templates
*/
template <typename T_FEATURE> 
class InParamDefineFeat {
public:
  InParamDefineFeat() {}
  ~InParamDefineFeat() {}
};

} /*END namespace inout
   */

#endif /*__IN_PARAM_DEFINE_DATA_TYPE_HPP*/
