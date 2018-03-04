/*! \file inparam_fixedk.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_FIXEDK_HPP
#define IN_PARAM_FIXEDK_HPP

#include "inparam_clustering.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
#define INPARAMCLUSTERING_DEFAULT_CLUSTERK     3
  
/*! \class InParamFk
  \brief Input parameter for algorithms  with Fixed K 
*/
template <class T_CLUSTERIDX> //-1, 0, 1, .., K
class InParamFk {
public:
  InParamFk() 
    : _idxK_numClusterK(INPARAMCLUSTERING_DEFAULT_CLUSTERK)
  {}
  ~InParamFk() {}

  inline void setNumClusterK(T_CLUSTERIDX aiidxK_numClusterK) 
  {
    this->_idxK_numClusterK = aiidxK_numClusterK;
  }

  inline const T_CLUSTERIDX getNumClusterK() const 
  {
    return this->_idxK_numClusterK;
  }

  virtual void print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    aipf_outFile << aic_separator << "_k"   
		 << aic_separator << this->_idxK_numClusterK; 
  
  }
protected:
  T_CLUSTERIDX _idxK_numClusterK; //-1, 0, 1, .., K
}; /*InParamFk*/

} /*END namespace inout
   */

#endif /*IN_PARAM_FIXEDK_HPP*/
