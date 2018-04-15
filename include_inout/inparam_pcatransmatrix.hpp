/*! \file inparam_pcatransmatrix.hpp
 *
 * \brief Definition of PCA program parameters for calculate transformation matrix
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_PCA_TRANS_MATRIX_HPP__
#define __IN_PARAM_PCA_TRANS_MATRIX_HPP__

#include "inparam.hpp"

#define  __INPARAM_PCATRANSMATRIX__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/


namespace  inout {

  
typedef enum {COVARIANCE, CORRELATION} ParamPCAmatrix;
typedef enum {EIGEN, SVD}              ParamPCAmethod;


/*! \class InParamPCAtransmatrix
  \brief Input parameter for PCA transformation matrix 
*/
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
class InParamPCAtransmatrix 
  : public InParamAlgorithmo
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public: 
  InParamPCAtransmatrix
  (std::string        ais_algorithmoName,
   std::string        ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut)
    :  InParamAlgorithmo(ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut)
    ,  InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , _ps_outFileTransMatrix(NULL)
  {}

  ~InParamPCAtransmatrix() {}


  inline void setOutFileTransMatrix(char* aips_fileName) 
  {
    this->_ps_outFileTransMatrix = aips_fileName;
  }

  inline char* getOutFileTransMatrix() 
  {
    return this->_ps_outFileTransMatrix;
  }
  
  virtual void print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamAlgorithmo::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }

protected:
 
  char   *_ps_outFileTransMatrix; 
  
}; /*InParamPCAtransmatrix*/

} /* END namespace inout*/

#endif /*__IN_PARAM_PCA_TRANS_MATRIX_HPP__*/
