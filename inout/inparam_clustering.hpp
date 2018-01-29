/*! \file inparam_clustering.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_CLUSTERING_HPP
#define IN_PARAM_CLUSTERING_HPP
  

#include <iostream>
#include <fstream>
#include <string>
#include "inparam.hpp"
#include "standardize_variable.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  

#define INPARAMCLUSTERING_MAXEXECUTIONTIME     36000.0

/* FOR CLUSTERK_MAX DEPENDS OF NUMBER INSTANCES SQRT(NUMBER INSTANCES)*/ 
#define INPARAMCLUSTERING_DEFAULT_MAX_ITER     100
#define INPARAMCLUSTERING_DEFAULT_MIN_THRESHOLD 0
#define INPARAMCLUSTERING_DEFAULT_EPSILON       1e-4
#define INPARAMCLUSTERING_DEFAULT_WEIGHTING_EXPONENT 2

#define INPARAMGACLUSTERING_PROBCROOSSOVER_IND -1.0f
#define INPARAMGACLUSTERING_PROBMUTATION_IND   -1.0f

#define INPARAMCLUSTERING_DISTANCE_EUCLIDEAN           0
#define INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_SQ        1
#define INPARAMCLUSTERING_DISTANCE_EUCLIDEAN_INDUCED   2
#define INPARAMCLUSTERING_DISTANCE_DIAGONAL_INDUCED    3    
#define INPARAMCLUSTERING_DISTANCE_MAHALONOBIS_INDUCED 4

#if  DATATYPE_CENTROIDS_ROUND == 0

#define INPARAMCLUSTERING_DISTANCE_TYPE {"euclidean", "euclidean_sq", "euclidean_induced", "diagonal_induced", "mahalonobis_induced", (char *) NULL }

#else

#define INPARAMCLUSTERING_DISTANCE_TYPE {"euclidean", "euclidean_sq", (char *) NULL }
  
#endif /*DATATYPE_CENTROIDS_ROUND*/

  
/*! \class InParamClustering
  \brief Input parameter for clustering algorithm 
*/
/*template < typename T_FEATURE,         
	   typename T_FEATURE_SUM,
	   typename T_INSTANCES_CLUSTER_K 
	   >
*/
class InParamClustering:
  public InParamAlgorithmo
{
public:
  InParamClustering
  (std::string        ais_algorithmoName,
   std::string        ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut,
   int                aii_opDistance)
    :  InParamAlgorithmo(ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut)
    ,  i_opDistance(aii_opDistance)
    ,  ps_outFileCentroids(NULL)
    ,  ps_outFileMemberShip(NULL)
    ,  ps_outFilePartitionsTable(NULL)
    ,  ps_outFileGraph(NULL)
    ,  ps_fileNamePlotStatObjetiveFunc(NULL)
    , _str_randomSeed()
    , _rd_maxExecutiontime(INPARAMCLUSTERING_MAXEXECUTIONTIME)
    , _b_printCentroidsFormat(false)
    , _b_printTableFormat(false)
  {}

  ~InParamClustering() {}
  
  /*FOR ALL ALGORITHS*/
 
  inline char* getFileNamePlotStatObjetiveFunc() const
  { 
    return this->ps_fileNamePlotStatObjetiveFunc; 
  }

  inline bool  getWithPlotStatObjetiveFunc() const
  {
    return (this->ps_fileNamePlotStatObjetiveFunc != NULL);
  }

    inline void setFileNamePlotStatObjetiveFunc(char *aips_nameFile) 
  {
    this->ps_fileNamePlotStatObjetiveFunc = aips_nameFile;
  }

  inline void setOutFileMemberShip(char* aips_fileName) 
  {
    this->ps_outFileMemberShip = aips_fileName;
  }

  inline char* getOutFileMemberShip() 
  {
    return this->ps_outFileMemberShip;
  }

  void setRandomSeed(const std::string& aistr_randomSeed) 
  {
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "InParamClustering::setRandomSeed";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
	      << "const std::string& aistr_randomSeed = " << aistr_randomSeed
	      << ")\n";
    }
#endif /*__VERBOSE_YES*/

    this->_str_randomSeed = aistr_randomSeed;

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc
	      << ": OUT(" << geiinparam_verbose << ")\n"
		 << "_str_randomSeed = " <<  this->_str_randomSeed
		 << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  
  }

  void setRandomSeed(const char* aipc_randomSeed) 
  {
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       std::cout << "random_seed  IN"
		 << '(' << geiinparam_verbose << ')'
		 << "\n\t const char* aipc_randomSeed = " << aipc_randomSeed
		 << "\n\t)\n";
    }
#endif /*__VERBOSE_YES*/

    this->_str_randomSeed = aipc_randomSeed;

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
       std::cout << "random_seed OUT"
		 << '(' << geiinparam_verbose << ')'
		 << "\n _str_randomSeed = " <<  this->_str_randomSeed
		 << std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/  
  }

  inline const std::string& getRandomSeed() const 
  {
    return this->_str_randomSeed;
  }

  inline void setOutFilePartitionsTable(char* aips_fileName) 
  {
    this->ps_outFilePartitionsTable = aips_fileName;
  }

  inline char* getOutFilePartitionsTable() 
  {
    return this->ps_outFilePartitionsTable;
  }

  inline void setOutFileCentroids(char* aips_fileName) 
  {
    this->ps_outFileCentroids = aips_fileName;
  }

  inline char* getOutFileCentroids() 
  {
    return this->ps_outFileCentroids;
  }

  inline void setOpDistance(int aii_opDistance) 
  {
    this->i_opDistance = aii_opDistance;
  }

  inline int getOpDistance() 
  {
    return this->i_opDistance;
  }

  inline double getMaxExecutiontime() 
  {
    return _rd_maxExecutiontime;
  }

  inline void setMaxExecutiontime(double aird_maxExecutiontime) 
  {
    _rd_maxExecutiontime = aird_maxExecutiontime;
  }

  inline bool getPrintCentroidsFormat() 
  {
    return  _b_printCentroidsFormat;
  }

  inline void setPrintCentroidsFormat(bool aib_printCentroidsFormat) 
  {
    _b_printCentroidsFormat = aib_printCentroidsFormat;
  }

  inline bool getPrintTableFormat() 
  {
    return  _b_printTableFormat;
  }

  inline void setPrintTableFormat(bool aib_printTableFormat) 
  {
    _b_printTableFormat = aib_printTableFormat;
  }

  inline void setOutFileGraph(char* aips_fileName) 
  {
    this->ps_outFileGraph = aips_fileName;
  }

  inline char* getOutFileGraph() 
  {
    return this->ps_outFileGraph;
  }
 
  
  virtual void print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    const char  *las_opFuncDistance[] = INPARAMCLUSTERING_DISTANCE_TYPE;

    InParamAlgorithmo::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_norm"                       
		 << aic_separator << las_opFuncDistance[this->i_opDistance];
    aipf_outFile << aic_separator << "_random seed" 
		 << aic_separator << this->_str_randomSeed;
    aipf_outFile << aic_separator << "_maximum execution time"
		 << aic_separator << _rd_maxExecutiontime;
    aipf_outFile << aic_separator << "_print centroids format"
		 << aic_separator << _b_printCentroidsFormat;
    aipf_outFile << aic_separator << "_print table format"
		 << aic_separator << _b_printTableFormat;
    
  }
  
protected:
  
  int                 i_opDistance; 
  char                *ps_outFileCentroids;
  char                *ps_outFileMemberShip;
  char                *ps_outFilePartitionsTable;
  char                *ps_outFileGraph;
  char                *ps_fileNamePlotStatObjetiveFunc;

  std::string         _str_randomSeed;
  double              _rd_maxExecutiontime;

  bool                _b_printCentroidsFormat;
  bool                _b_printTableFormat;
 
}; /*InParamClustering*/

}

#endif /*IN_PARAM_CLUSTERING_HPP*/
