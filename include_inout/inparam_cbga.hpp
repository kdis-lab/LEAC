/*! \file inparam_cbga.hpp
 *
 * \brief Definition of CBGA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __IN_PARAM_CBGA_HPP__
#define __IN_PARAM_CBGA_HPP__

#include "inparam_fixedk.hpp"
#include "inparam_pm.hpp"
#include "inparam_readinst.hpp"

#define __INPARAM_CBGA__

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
#define INPARAMCLUSTERING_CBGA_SELECMETH_ROULETTE 0
#define INPARAMCLUSTERING_CBGA_SELECMETH_ELITIST1 1
#define INPARAMCLUSTERING_CBGA_SELECMETH_ELITIST2 2
#define INPARAMCLUSTERING_CBGA_SELECMETH_ZIGZAG   3
#define INPARAMCLUSTERING_CBGA_SELECMETH {"roulette", "elitist1", "elitist2", "zigzag", (char *) NULL }  
  
/*! \class InParamCBGA
  \brief Input parameter for CBGA \cite Franti:etal:GAclustering:gafranti:1997
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_REAL,
	   typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_INSTANCE_FREQUENCY
	   >
class InParamCBGA
  : public InParamPm<T_REAL>
  , public InParamFk<T_CLUSTERIDX>
  , public InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY>
{
public:
  InParamCBGA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm,
   int                aii_opSelectMethod,
   int                aii_numGLAIterations 
   ):
    InParamPm<T_REAL>
    (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut, aii_opNorm)
    , InParamFk<T_CLUSTERIDX>()
    , InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY>()
    , i_opSelectMethod(aii_opSelectMethod)
    , i_numGLAIterations(aii_numGLAIterations) 
  {}

  ~InParamCBGA() {}

  inline void setOpSelectMethod(int aii_opSelectMethod) 
  {
    this->i_opSelectMethod = aii_opSelectMethod;
  }

  inline const  int getOpSelectMethod() const  
  {
    return this->i_opSelectMethod;
  }

  inline void setNumGLAIterations(int aii_numGLAIterations) 
  {
    this->i_numGLAIterations = aii_numGLAIterations;
  }

  inline const int getNumGLAIterations() const 
  {
    return this->i_numGLAIterations;
  }

  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    const char  *las_opSelectMethod[] = INPARAMCLUSTERING_CBGA_SELECMETH;

    InParamPm<T_REAL>::print(aipf_outFile,aic_separator);
    InParamFk<T_CLUSTERIDX>::print(aipf_outFile,aic_separator); 
    InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY>::print(aipf_outFile,aic_separator);
      aipf_outFile << aic_separator << "_selection method"            
		   << aic_separator << las_opSelectMethod[this->i_opSelectMethod]; 
      aipf_outFile << aic_separator << "_number of GLA iterations"
		   << aic_separator << i_numGLAIterations;
  }
  
protected:
  
  int    i_opSelectMethod; 
  int    i_numGLAIterations;

}; 

} /* END namespace inparam
   */

#endif /*__IN_PARAM_CBGA_HPP__*/
