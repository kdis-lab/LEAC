/*! \file inparam_clustering_max_iter.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_CLUSTERING_MAX_ITER_HPP
#define IN_PARAM_CLUSTERING_MAX_ITER_HPP

#include "inparam_clustering.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamClusteringMaxIter
  \brief Input parameter for algorithm with iteration maximum
*/
class InParamClusteringMaxIter:
    public InParamClustering
{
public:
  InParamClusteringMaxIter
  (std::string        ais_algorithmoName,
   std::string        ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opNorm
   ):
    InParamClustering
    (ais_algorithmoName,
     ais_algorithmoAuthor,
     aiato_algTypeOut,
     aii_opNorm) {}

  ~InParamClusteringMaxIter() {}
  
  inline void setNumMaxIter(COMMON_IDOMAIN aiiT_numMaxIter) 
  {
    this->iT_numMaxIter = aiiT_numMaxIter;
  }

  inline COMMON_IDOMAIN getNumMaxIter()
  {
    return this->iT_numMaxIter;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamClustering
      ::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_number maximum iterations" 
		 << aic_separator <<  this->iT_numMaxIter;
  }
protected:
  COMMON_IDOMAIN  iT_numMaxIter;
};


} /*END namespace inout 
   */

#endif /*IN_PARAM_CLUSTERING_MAX_ITER_HPP*/
