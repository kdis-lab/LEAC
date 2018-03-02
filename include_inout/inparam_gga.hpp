/*! \file inparam_gaclustering_gga.hpp
 *
 * \brief Definition of GGA program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_GACLUSTERING_GGA_HPP
#define IN_PARAM_GACLUSTERING_GGA_HPP

#include "inparam_clustering.hpp"
#include "inparam_rangek.hpp"
#include "inparam_readinst.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class InParamGAClusteringGGA
  \brief Input parameter for GGA \cite Agustin:etal:GAclusteringVarK:GGA:2012
*/
template < typename T_CLUSTERIDX, //-1, 0, 1, .., K
	   typename T_REAL,
	   typename T_FEATURE,
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K
	   >
class InParamGAClusteringGGA
  : public InParamClustering
  , public InParamRangeK<T_CLUSTERIDX>
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public:
  InParamGAClusteringGGA
  (const std::string& ais_algorithmoName,
   const std::string& ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int         aii_opNorm
   )
    : InParamClustering
      (ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut,aii_opNorm)
    , InParamRangeK<T_CLUSTERIDX>()
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , _ui_subPopulationSize(20)
    , _ui_numIsland(4)
    , _rt_pe(0.5)
    , _it_numMaxGenerations(100)
    , _rt_pci(0.8)  //Probability Crossover Initial
    , _rt_pcf(0.4)  //Probability Crossover Final
    , _rt_pmi(0.05) //Probability Mutation Initial 
    , _rt_pmf(0.2)  //Probability Mutation Final
    , _rt_pbi(0.1)  //Probability local search
    , _rt_pbf(0.05)  //Probability local search
  {}

  ~InParamGAClusteringGGA() {}

  inline void setSubPopulationSize(uintidx aiui_subPopulationSize) 
  {
    this->_ui_subPopulationSize = aiui_subPopulationSize;
  }

  inline uintidx getSubPopulationSize() 
  {
    return this->_ui_subPopulationSize;
  }

  inline void setNumIsland(uintidx aiui_numIsland) 
  {
    this->_ui_numIsland = aiui_numIsland;
  }

  inline uintidx getNumIsland() 
  {
    return this->_ui_numIsland;
  }

  inline void setPe(T_REAL airt_pe) 
  {
    this->_rt_pe = airt_pe;
  }

  inline T_REAL getPe()
  {
    return this->_rt_pe;
  }

  inline void setNumMaxGenerations(COMMON_IDOMAIN aiT_numMaxGenerations) 
  {
    this->_it_numMaxGenerations = aiT_numMaxGenerations;
  }

  inline COMMON_IDOMAIN getNumMaxGenerations() 
  {
    return this->_it_numMaxGenerations;
  }

  inline void setPci(T_REAL airt_pci) 
  {
    this->_rt_pci = airt_pci;
  }

  inline T_REAL getPci()
  {
    return this->_rt_pci;
  }

  inline void setPcf(T_REAL airt_pcf) 
  {
    this->_rt_pcf = airt_pcf;
  }

  inline T_REAL getPcf()
  {
    return this->_rt_pcf;
  }

  inline void setPmi(T_REAL airt_pmi) 
  {
    this->_rt_pmi = airt_pmi;
  }

  inline T_REAL getPmi()
  {
    return this->_rt_pmi;
  }

  inline void setPmf(T_REAL airt_pmf) 
  {
    this->_rt_pmf = airt_pmf;
  }

  inline T_REAL getPmf()
  {
    return this->_rt_pmf;
  }


  inline void setPbi(T_REAL airt_pbi) 
  {
    this->_rt_pbi = airt_pbi;
  }

  inline T_REAL getPbi()
  {
    return this->_rt_pbi;
  }

  inline void setPbf(T_REAL airt_pbf) 
  {
    this->_rt_pbf = airt_pbf;
  }

  inline T_REAL getPbf()
  {
    return this->_rt_pbf;
  }
  
  virtual void  print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    InParamClustering::print(aipf_outFile,aic_separator);
    InParamRangeK<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_sub-population-size"   
		 << aic_separator << this->_ui_subPopulationSize;
    aipf_outFile << aic_separator << "_number-island"   
		 << aic_separator << this->_ui_numIsland;
    aipf_outFile << aic_separator << "_probability of migration"   
		 << aic_separator << this->_rt_pe;
    aipf_outFile << aic_separator << "_number maximum generations" 
		 << aic_separator << this->_it_numMaxGenerations;
    InParamRangeK<T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_initial probability crossover"   
		 << aic_separator << this->_rt_pci;
    aipf_outFile << aic_separator << "_final probability crossover"   
		 << aic_separator << this->_rt_pcf;
    aipf_outFile << aic_separator << "_initial probability mutation" 
		 << aic_separator << this->_rt_pmi;
    aipf_outFile << aic_separator << "_final probability mutation" 
		 << aic_separator << this->_rt_pmf;

    aipf_outFile << aic_separator << "_initial probability local search"   
		 << aic_separator << this->_rt_pbi;
    aipf_outFile << aic_separator << "_final probability local searchr"   
		 << aic_separator << this->_rt_pbf;
    
  }
  
protected:
  uintidx         _ui_subPopulationSize;
  uintidx         _ui_numIsland;
  T_REAL   _rt_pe; /*probability of migration*/
  COMMON_IDOMAIN _it_numMaxGenerations;
  
  T_REAL  _rt_pci; /*probCrossoverInitial*/
  T_REAL  _rt_pcf; /*probCrossoverFinal*/
  T_REAL  _rt_pmi; /*probMutationInitial*/ 
  T_REAL  _rt_pmf; /*probMutationFinal*/
  T_REAL  _rt_pbi; /*Probability local search*/
  T_REAL  _rt_pbf; /*Probability local search*/
};

} /* END namespace inout
   */

#endif /*IN_PARAM_GACLUSTERING_GGA_HPP*/
