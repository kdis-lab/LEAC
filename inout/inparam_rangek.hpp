/*! \file inparam_rangek.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_RANGEK_HPP
#define IN_PARAM_RANGEK_HPP


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

#define INPARAMCLUSTERING_DEFAULT_CLUSTERK_MIN 2
#define INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED -1


template <class T_CLUSTERIDX> //-1, 0, 1, .., K
class InParamRangeK {
public:
  InParamRangeK() 
    : _idxK_numClusterKMinimum(INPARAMCLUSTERING_DEFAULT_CLUSTERK_MIN)
    , _idxK_numClusterKMaximum(INPARAMCLUSTERING_DEFAULT_CLUSTERK_UNDEFINED)
  {}
  ~InParamRangeK() {}

  inline void setNumClusterKMinimum(T_CLUSTERIDX aiidxK_numClusterKMinimum) 
  {
    this->_idxK_numClusterKMinimum = aiidxK_numClusterKMinimum;
  }

  inline const T_CLUSTERIDX getNumClusterKMinimum() 
  {
    return this->_idxK_numClusterKMinimum;
  }

  inline void setNumClusterKMaximum(T_CLUSTERIDX aiidxK_numClusterKMaximum) 
  {
    this->_idxK_numClusterKMaximum = aiidxK_numClusterKMaximum;
  }

  inline const T_CLUSTERIDX getNumClusterKMaximum() 
  {
    return this->_idxK_numClusterKMaximum;
  }

  virtual void print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    aipf_outFile << aic_separator << "_k-minimum"   
		 << aic_separator << this->_idxK_numClusterKMinimum
		 << aic_separator << "_k-maximum"   
		 << aic_separator << this->_idxK_numClusterKMaximum;
  }

protected:
  
  T_CLUSTERIDX _idxK_numClusterKMinimum;//-1, 0, 1, .., K,    // st_numClusterK;
  T_CLUSTERIDX _idxK_numClusterKMaximum;//-1, 0, 1, .., K,    // st_numClusterK;
  
}; /*InParamRangeK*/
  

} /*END namespace inout
   */

#endif /*IN_PARAM_RANGEK_HPP*/
