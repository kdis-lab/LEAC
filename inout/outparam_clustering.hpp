/*! \file outparam_clustering.hpp
 *
 * \brief out parameters clustering
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __OUT_PARAM_CLUSTERING_HPP
#define __OUT_PARAM_CLUSTERING_HPP

#include <iostream>
#include <string>
#include <time.h>

#include "outfilename.hpp"
#include "execution_time.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
  
#define OUTPARAMCLUSTERING_METRIC_NaN      -1.0
#define OUTPARAMCLUSTERING_FITNESS_NaN     -1.0
#define OUTPARAMCLUSTERING_INT_NaN         -1

#define OUTPARAMCLUSTERING_FCMPARTCOEF_NaN -1.0
#define OUTPARAMCLUSTERING_FCMENTROPY_NaN  -1.0
#define OUTPARAMCLUSTERING_NaN             -1.0
#define OUTPARAMCLUSTERING_MAX_SIZE_NAME_FILE 50


const char *garray_nameObjetiveFunc[] = {
  
  /*0*/ "SSE",
  /*1*/ "TWCV",
  /*2*/ "Distortion", /*\cite Franti:etal:GAclustering:gafranti:1997*/
  /*3*/ "J1",
 
  /*4*/ "Index I",
  /*5*/ "CS measure",
  /*6*/ "Dunn's index",
  /*7*/ "SDunn's index",
  /*8*/ "Silhouette",
  /*9*/ "Simplified Silhouette",
  /*10*/ "DB-index",
  /*11*/ "Variance Ratio Criterion",
  /*12*/ "Score Function",
  
  /*13*/ "Rand Index",
  /*14*/ "Purity",
  /*15*/ "F-measure",
  /*16*/ "Jaccard Index",
  /*17*/ "precision",
  /*18*/ "recall",
  /*19*/ "Misclassified",
  /*20*/ "pairs a",
  /*21*/ "pairs b",
  /*22*/ "pairs c",
  /*23*/ "pairs d",

  /*24*/ "Jm",
  /*25*/ "Xie-Beni index",
  /*26*/ "Entropy",
  /*27*/ "Partition Coefficient",
  /*28*/ "Intra- and inter-clust",
  
  /*29*/ "NotUse",
  (char *) NULL
};
  
typedef enum {
  
  SSE=0,
  TWCV=1,
  Distortion=2, 
  J1=3,
  
  IndexI=4,
  CSmeasure=5,
  DunnIndex=6,
  SDunnIndex=7,
  Silhouette=8,
  SSilhouette=9,
  DBindex=10,
  VRC=11, 
  ScoreFunction=12,
  
  RandIndex=13,
  Purity=14,
  Fmeasure=15,
  JaccardIndex=16,
  precision=17,
  recall=18,
  Misclassified=19,
  Pairs_a=20,
  Pairs_b=21,
  Pairs_c=22,
  Pairs_d=23,

  Jm=24,
  XieBeniIndex=25,
  Entropy=26,
  PartitionCoefficient=27,
  IntraInterClust=28,
  
  NotUse=29
  
} OutParamNameObjectiveFunc;


/*! \class OutParamAlgorithm
  \brief Output parameters for genetic algorithm
*/
class OutParamAlgorithm  {
public:
  OutParamAlgorithm()
  {
  }
  
  virtual ~OutParamAlgorithm() { }

  void initialize(int aii_numRunAlgorithm)
  {
    this->i_numRunAlgorithm   = aii_numRunAlgorithm;
    this->et_algorithmRunTime = 0.0;
    this->_b_endingCondition  = true;
    this->_b_endingConditionTest  = true;
    this->_s_runtimeMessage   = "ok";
  }

  inline void setAlgorithmRunTime(runtime::ExecutionTime aiet_algorithmRunTime)  
  {
    this->et_algorithmRunTime = aiet_algorithmRunTime;
  }

  inline  runtime::ExecutionTime getAlgorithmRunTime()  
  {
    return this->et_algorithmRunTime;
  }

  inline const int getNumRunningAlgorithm() 
  {
    return this->i_numRunAlgorithm;
  }


  inline void setEndingCondition(const bool aib_endingCondition)
  {
    this->_b_endingCondition = aib_endingCondition;
  }

  inline void setEndingConditionTest(const bool aib_endingCondition)
  {
    this->_b_endingConditionTest = aib_endingCondition;
  }

  inline void setRuntimeMessage(const std::string ais_runtimeMessage) 
  {
    this->_s_runtimeMessage = ais_runtimeMessage;
  }
  
  inline const bool getEndingCondition()
  {
    return this->_b_endingCondition;
  }

   inline const bool getEndingConditionTest()
  {
    return this->_b_endingConditionTest;
  }

  inline const std::string& getRuntimeMessage() 
  {
    return this->_s_runtimeMessage;
  }

  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
   aipf_outFile	 << aic_separator << "_number run" 
		 << aic_separator << this->i_numRunAlgorithm
		 << aic_separator << "_execution time seconds" 
		 << aic_separator << this->et_algorithmRunTime
		 << aic_separator << "_end condition" 
		 << aic_separator << this->_b_endingCondition
		 << aic_separator << "_end condition test" 
		 << aic_separator << this->_b_endingConditionTest
		 << aic_separator << "_run time message" 
		 << aic_separator << this->_s_runtimeMessage;
    
  }
  friend std::ostream& operator<<(std::ostream& os, const OutParamAlgorithm  &aioutparamalg_outParam)
  {
    aioutparamalg_outParam.print(os,OutFileName::getDelim());
    
    return os;
  }

protected:

  int           i_numRunAlgorithm;
  runtime::ExecutionTime et_algorithmRunTime;
  bool          _b_endingCondition;
  bool          _b_endingConditionTest;
  std::string   _s_runtimeMessage;
  
}; /*END class OutParamAlgorithm*/


/*! \class OutParamClustering
  \brief Output parameters for a generic clustering algorithm 
*/
template < typename T_METRIC,
	   typename T_CLUSTERIDX
	   >
class OutParamClustering:
  public OutParamAlgorithm
{
public:
  OutParamClustering(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc):
     OutParamAlgorithm()
    , _enum_usedObjectiveFunc(aienum_usedObjectiveFunc)
    , s_fileNameOutPlotStatObjetiveFunc()
     
     //  , _str_percentageSensitivity()
     //, _str_percentageSensitivityTest()
     //, _str_percentageSpecificity()
     //, _str_percentageSpecificityTest()
     
  {
   
    //uintidx  luintidx_numObjetiveFunc  = this->getNumObjetiveFunc();

    // _arrayrt_metricFunc     = new T_METRIC[luintidx_numObjetiveFunc];
    //_arrayrt_metricFuncTest = new T_METRIC[luintidx_numObjetiveFunc];
  }
 
  virtual ~OutParamClustering()
  {
    //delete[] _arrayrt_metricFunc;
    //delete[] _arrayrt_metricFuncTest;
  }

  /*
  const uintidx getNumObjetiveFunc() const
  {
    uintidx  luintidx_numObjetiveFunc;

    luintidx_numObjetiveFunc = 0;
    while ( garray_nameObjetiveFunc[luintidx_numObjetiveFunc] !=  NULL ) {
      ++luintidx_numObjetiveFunc;
      }
    --luintidx_numObjetiveFunc; //FOR NOT USE METRIC

    return luintidx_numObjetiveFunc;

  }

  */
  
  void initialize(int aii_numRunAlgorithm)
  {
   
    OutParamAlgorithm::initialize(aii_numRunAlgorithm);

    _rt_objetiveFuncRun = OUTPARAMCLUSTERING_METRIC_NaN;
   
    this->_it_iterationGetsBest  = 0; 
    this->et_runTimeGetsBest    = 0.0;
   
    this->_l_totalInvalidOffspring = OUTPARAMCLUSTERING_INT_NaN;
    this->_idxK_numClusterK     = OUTPARAMCLUSTERING_INT_NaN;
    this->_idxK_numClusterKTest = OUTPARAMCLUSTERING_INT_NaN;

    /*uintidx  luintidx_numObjetiveFunc  = this->getNumObjetiveFunc();
    for (uintidx lui_i = 0; lui_i < luintidx_numObjetiveFunc; ++lui_i) {
      this->_arrayrt_metricFunc[lui_i] = OUTPARAMCLUSTERING_METRIC_NaN;
      this->_arrayrt_metricFuncTest[lui_i] = OUTPARAMCLUSTERING_METRIC_NaN;
    }
    _str_percentageSensitivity = "";
    _str_percentageSensitivityTest = "";
    _str_percentageSpecificity = "";
    _str_percentageSpecificityTest = "";
    */
  }
                                                                    
  inline void setRunTimeGetsBest(runtime::ExecutionTime aiet_runTimeGetsBest)  
  {
    this->et_runTimeGetsBest = aiet_runTimeGetsBest;
  }
  
  inline void setIterationGetsBest(COMMON_IDOMAIN aiit_iterationGetsBest) 
  {
    this->_it_iterationGetsBest = aiit_iterationGetsBest;
  }

  inline runtime::ExecutionTime getRunTimeGetsBest()
  {
    return this->et_runTimeGetsBest;
  }

  inline COMMON_IDOMAIN getIterationGetsBest()
  {
    return this->_it_iterationGetsBest;
  }

  void setTotalInvalidOffspring(long ail_totalInvalidOffspring) 
  {
    this->_l_totalInvalidOffspring = ail_totalInvalidOffspring;
  }

  void incTotalInvalidOffspring() 
  {
    ++this->_l_totalInvalidOffspring;
  }

  void sumTotalInvalidOffspring(long ail_invalidOffspring ) 
  {
    this->_l_totalInvalidOffspring += ail_invalidOffspring;
  }

  long getTotalInvalidOffspring() 
  {
    return this->_l_totalInvalidOffspring;
  }

  inline void setMetricFuncRun
  (T_METRIC  ait_objetiveFuncRun)	
  {
    _rt_objetiveFuncRun = ait_objetiveFuncRun;
  }

  inline const  T_METRIC getObjetiveFuncRun() const	
  {
    return _rt_objetiveFuncRun;
  }

  /*
  inline void setMetricFunc
  (OutParamNameObjectiveFunc aienum_usedObjectiveFunc, 
   T_METRIC            ait_objetiveFunc)	
  {
    this->_arrayrt_metricFunc[aienum_usedObjectiveFunc] = ait_objetiveFunc;
  }

  inline void setMetricFuncTest
  (OutParamNameObjectiveFunc aienum_usedObjectiveFunc, 
   T_METRIC            ait_objetiveFuncTest)	
  {
    this->_arrayrt_metricFuncTest[aienum_usedObjectiveFunc] = ait_objetiveFuncTest;
  }
  

  inline T_METRIC getObjetiveFunc
  (OutParamNameObjectiveFunc aienum_metricFunc)	
  {
    return this->_arrayrt_metricFunc[aienum_metricFunc];
  }

  inline T_METRIC getObjetiveFuncTest
  (OutParamNameObjectiveFunc aienum_metricFunc)	
  {
    return this->_arrayrt_metricFuncTest[aienum_metricFunc];
  }
  
  
  inline void setUsedObjetiveFunc(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc)
  {
    this->_enum_usedObjectiveFunc = aienum_usedObjectiveFunc;
  }

  */
  
  inline const OutParamNameObjectiveFunc getUsedObjetiveFunc()
  {
    return this->_enum_usedObjectiveFunc;
  }

  inline const std::string  getNameUsedObjetiveFunc()
  {
    return std::string(garray_nameObjetiveFunc[_enum_usedObjectiveFunc]);
  }

  void  setFileNameOutPlotStatObjetiveFunc
  (const char* aipstr_fileNamePlot,
   const int   aii_timesRunAlgorithm
   )
  {
    char ls_fileNameOutPlotStatObjetiveFunc[OUTPARAMCLUSTERING_MAX_SIZE_NAME_FILE];
    
    sprintf(ls_fileNameOutPlotStatObjetiveFunc,
            "%s_%02d_%02d.dat",
	    aipstr_fileNamePlot, 
            this->i_numRunAlgorithm,
            aii_timesRunAlgorithm
            );

    s_fileNameOutPlotStatObjetiveFunc.assign(ls_fileNameOutPlotStatObjetiveFunc);
  }

  void  setFileNameOutPlotStatObjetiveFunc
   (const char*       aipstr_fileNamePlot,
    const int          aii_timesRunAlgorithm,
   const std::string& aistr_prefix
   )
  {
    char ls_fileNameOutPlotStatObjetiveFunc[OUTPARAMCLUSTERING_MAX_SIZE_NAME_FILE];
    
    sprintf(ls_fileNameOutPlotStatObjetiveFunc,
            "%s_%s_%02d_%02d.dat",
	    aipstr_fileNamePlot,
	    aistr_prefix.c_str(),
            aii_timesRunAlgorithm
            );

    s_fileNameOutPlotStatObjetiveFunc.assign(ls_fileNameOutPlotStatObjetiveFunc);
  }
  
  inline const std::string getFileNameOutPlotStatObjetiveFunc()
  {
    return s_fileNameOutPlotStatObjetiveFunc;
  }

  inline void setNumClusterK(T_CLUSTERIDX aiidxK_numClusterK) 
  {
    this->_idxK_numClusterK = aiidxK_numClusterK;
  }

  inline const T_CLUSTERIDX getNumClusterK() 
  {
    return this->_idxK_numClusterK; 
  }

  inline void setNumClusterKTest(T_CLUSTERIDX aiidxK_numClusterKTest) 
  {
    this->_idxK_numClusterKTest = aiidxK_numClusterKTest;
  }

  inline const T_CLUSTERIDX getNumClusterKTest() 
  {
    return this->_idxK_numClusterKTest; 
  }


  /*
  inline void setPercentageSensitivity(std::string aistr_percentageSensitivity) 
  {
    _str_percentageSensitivity = aistr_percentageSensitivity;
  }

  inline const std::string& getPercentageSensitivity() 
  {
    return _str_percentageSensitivity;
  }

  inline void setPercentageSensitivityTest(std::string aistr_percentageSensitivityTest) 
  {
    _str_percentageSensitivityTest = aistr_percentageSensitivityTest;
  }

  inline const std::string& getPercentageSensitivityTest() 
  {
    return _str_percentageSensitivityTest;
  }

  
  inline void setPercentageSpecificity(std::string aistr_percentageSpecificity) 
  {
    _str_percentageSpecificity = aistr_percentageSpecificity;
  }

  inline const std::string& getPercentageSpecificity() 
  {
    return _str_percentageSpecificity;
  }

  inline void setPercentageSpecificityTest(std::string aistr_percentageSpecificityTest) 
  {
    _str_percentageSpecificityTest = aistr_percentageSpecificityTest;
  }

  inline const std::string& getPercentageSpecificityTest() 
  {
    return _str_percentageSpecificityTest;
  }
  */
  
  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    
    OutParamAlgorithm::print(aipf_outFile,aic_separator);
    aipf_outFile << aic_separator << "_file name plot funtion hist";
    if ( s_fileNameOutPlotStatObjetiveFunc.length() > 0 ) {  
      aipf_outFile << aic_separator
		   << this->s_fileNameOutPlotStatObjetiveFunc
		   << ".plot";
    }
    else {
      aipf_outFile << aic_separator << "none";
    }
    aipf_outFile << aic_separator << "_metric used" 
		 << aic_separator
		 << garray_nameObjetiveFunc[_enum_usedObjectiveFunc];

    /*
    uintidx  luintidx_numObjetiveFunc  = this->getNumObjetiveFunc();
    for (uintidx lui_i = 0; lui_i < luintidx_numObjetiveFunc; ++lui_i) {
      aipf_outFile <<  aic_separator << "_" << garray_nameObjetiveFunc[lui_i]
		   << aic_separator << this->_arrayrt_metricFunc[lui_i];
    }

    for (uintidx lui_i = 0; lui_i < luintidx_numObjetiveFunc; ++lui_i) {
      aipf_outFile << aic_separator << ":" << garray_nameObjetiveFunc[lui_i]
		   << aic_separator << this->_arrayrt_metricFuncTest[lui_i];
    }
    */

    aipf_outFile << aic_separator << "_iterations need" 
		 << aic_separator << this->_it_iterationGetsBest;
    aipf_outFile << aic_separator << "_time seconds need best" 
		 << aic_separator << this->et_runTimeGetsBest;
    aipf_outFile << aic_separator << "_outK"   
		 << aic_separator << this->_idxK_numClusterK;
    aipf_outFile << aic_separator << ":outK"   
		 << aic_separator << this->_idxK_numClusterKTest;
   
    aipf_outFile << aic_separator << "_solution overridden in the run of the algorithm" 
		 << aic_separator << this->_l_totalInvalidOffspring;
    
    //aipf_outFile << aic_separator << "_sensitivity"   
    //		 << aic_separator << _str_percentageSensitivity;
    
    //aipf_outFile << aic_separator << "_specificity"   
    //		 << aic_separator <<  _str_percentageSpecificity;

    //aipf_outFile << aic_separator << ":sensitivity"   
    //		 << aic_separator << _str_percentageSensitivityTest;
    
    //aipf_outFile << aic_separator << ":specificity"   
    //		 << aic_separator <<  _str_percentageSpecificityTest;

     aipf_outFile << aic_separator << "_objetivefuncrun"
		  << aic_separator << _rt_objetiveFuncRun;
   
  }
protected:
  /*sum of squared error (SSE) metric (clustering metric)*/
  /*within group sum of squares (WGSS) or minimum variance payoff 
   J_1:M_c x L_c --> R \cite{Bezdek:etal:GAclustering:GA:1994} */
  OutParamNameObjectiveFunc _enum_usedObjectiveFunc;
 
  T_METRIC                  _rt_objetiveFuncRun;
  //T_METRIC                  *_arrayrt_metricFunc;
  //T_METRIC                  *_arrayrt_metricFuncTest;
  COMMON_IDOMAIN           _it_iterationGetsBest;
  runtime::ExecutionTime    et_runTimeGetsBest;
  T_CLUSTERIDX             _idxK_numClusterK; //-1, 0, 1, .., K, K ALGORITHMS SEEK TO OTHERWISE THAT K EQUAL THE INPUT OUTPUT
  T_CLUSTERIDX             _idxK_numClusterKTest;
  std::string               s_fileNameOutPlotStatObjetiveFunc;  
  
  long                      _l_totalInvalidOffspring;
  
  /*std::string               _str_percentageSensitivity;
  std::string               _str_percentageSensitivityTest;
  std::string               _str_percentageSpecificity;
  std::string               _str_percentageSpecificityTest;
  */
  
}; /*END clas OutParamClustering*/


/*! \class OutParamClusteringMetric
  \brief Output parameters for a generic clustering algorithm 
*/
template < typename T_METRIC
	   //, typename T_CLUSTERIDX
	   >
class OutParamClusteringMetric
    // public OutParamClustering
{
public:
   OutParamClusteringMetric() :
  // (const OutParamNameObjectiveFunc aienum_usedObjectiveFunc):
    //OutParamClustering()
     //, _enum_usedObjectiveFunc(aienum_usedObjectiveFunc)
     //, s_fileNameOutPlotStatObjetiveFunc()
     
      _str_percentageSensitivity()
    , _str_percentageSensitivityTest()
    , _str_percentageSpecificity()
    , _str_percentageSpecificityTest()
     
  {
   
    uintidx  luintidx_numObjetiveFunc  = this->getNumObjetiveFunc();

    _arrayrt_metricFunc     = new T_METRIC[luintidx_numObjetiveFunc];
    _arrayrt_metricFuncTest = new T_METRIC[luintidx_numObjetiveFunc];
    
  }
 
  virtual ~OutParamClusteringMetric()
  {
    delete[] _arrayrt_metricFunc;
    delete[] _arrayrt_metricFuncTest;
  }

  const uintidx getNumObjetiveFunc() const
  {
    uintidx  luintidx_numObjetiveFunc;

    luintidx_numObjetiveFunc = 0;
    while ( garray_nameObjetiveFunc[luintidx_numObjetiveFunc] !=  NULL ) {
      ++luintidx_numObjetiveFunc;
      }
    --luintidx_numObjetiveFunc; //FOR NOT USE METRIC

    return luintidx_numObjetiveFunc;

  }

  void initialize(int aii_numRunAlgorithm)
  {
   
    //OutParamClustering::initialize(aii_numRunAlgorithm);

    uintidx  luintidx_numObjetiveFunc  = this->getNumObjetiveFunc();

    //_rt_objetiveFuncRun = OUTPARAMCLUSTERING_METRIC_NaN;
   
    for (uintidx lui_i = 0; lui_i < luintidx_numObjetiveFunc; ++lui_i) {
      this->_arrayrt_metricFunc[lui_i] = OUTPARAMCLUSTERING_METRIC_NaN;
      this->_arrayrt_metricFuncTest[lui_i] = OUTPARAMCLUSTERING_METRIC_NaN;
    }
   
    /*this->_it_iterationGetsBest  = 0; 
    this->et_runTimeGetsBest    = 0.0;
   
    this->_l_totalInvalidOffspring = 
    this->_idxK_numClusterK   = OUTPARAMCLUSTERING_INT_NaN;
    this->_idxK_numClusterKTest  = OUTPARAMCLUSTERING_INT_NaN;
    */
    _str_percentageSensitivity = "";
    _str_percentageSensitivityTest = "";
    _str_percentageSpecificity = "";
    _str_percentageSpecificityTest = "";
  }

  /*
  inline void setRunTimeGetsBest(runtime::ExecutionTime aiet_runTimeGetsBest)  
  {
    this->et_runTimeGetsBest = aiet_runTimeGetsBest;
  }
  
  inline void setIterationGetsBest(COMMON_IDOMAIN aiit_iterationGetsBest) 
  {
    this->_it_iterationGetsBest = aiit_iterationGetsBest;
  }

  inline runtime::ExecutionTime getRunTimeGetsBest()
  {
    return this->et_runTimeGetsBest;
  }

  inline COMMON_IDOMAIN getIterationGetsBest()
  {
    return this->_it_iterationGetsBest;
  }

  void setTotalInvalidOffspring(long ail_totalInvalidOffspring) 
  {
    this->_l_totalInvalidOffspring = ail_totalInvalidOffspring;
  }

  void incTotalInvalidOffspring() 
  {
    ++this->_l_totalInvalidOffspring;
  }

  void sumTotalInvalidOffspring(long ail_invalidOffspring ) 
  {
    this->_l_totalInvalidOffspring += ail_invalidOffspring;
  }

  long getTotalInvalidOffspring() 
  {
    return this->_l_totalInvalidOffspring;
  }

  inline void setMetricFuncRun
  (T_METRIC  ait_objetiveFuncRun)	
  {
    _rt_objetiveFuncRun = ait_objetiveFuncRun;
  }

  inline const  T_METRIC getObjetiveFuncRun() const	
  {
    return _rt_objetiveFuncRun;
  }
  */
  inline void setMetricFunc
  (OutParamNameObjectiveFunc aienum_usedObjectiveFunc, 
   T_METRIC            ait_objetiveFunc)	
  {
    this->_arrayrt_metricFunc[aienum_usedObjectiveFunc] = ait_objetiveFunc;
  }

  inline void setMetricFuncTest
  (OutParamNameObjectiveFunc aienum_usedObjectiveFunc, 
   T_METRIC            ait_objetiveFuncTest)	
  {
    this->_arrayrt_metricFuncTest[aienum_usedObjectiveFunc] = ait_objetiveFuncTest;
  }

  inline T_METRIC getObjetiveFunc
  (OutParamNameObjectiveFunc aienum_metricFunc)	
  {
    return this->_arrayrt_metricFunc[aienum_metricFunc];
  }

  inline T_METRIC getObjetiveFuncTest
  (OutParamNameObjectiveFunc aienum_metricFunc)	
  {
    return this->_arrayrt_metricFuncTest[aienum_metricFunc];
  }

  inline void setUsedObjetiveFunc(const OutParamNameObjectiveFunc aienum_usedObjectiveFunc)
  {
    this->_enum_usedObjectiveFunc = aienum_usedObjectiveFunc;
  }

  inline const OutParamNameObjectiveFunc getUsedObjetiveFunc()
  {
    return this->_enum_usedObjectiveFunc;
  }

  /*
  void  setFileNameOutPlotStatObjetiveFunc
  (const char* aipstr_fileNamePlot,
   const int   aii_timesRunAlgorithm
   )
  {
    char ls_fileNameOutPlotStatObjetiveFunc[OUTPARAMCLUSTERING_MAX_SIZE_NAME_FILE];
    
    sprintf(ls_fileNameOutPlotStatObjetiveFunc,
            "%s_%02d_%02d.dat",
	    aipstr_fileNamePlot, 
            this->i_numRunAlgorithm,
            aii_timesRunAlgorithm
            );

    s_fileNameOutPlotStatObjetiveFunc.assign(ls_fileNameOutPlotStatObjetiveFunc);
  }

  void  setFileNameOutPlotStatObjetiveFunc
   (const char*       aipstr_fileNamePlot,
    const int          aii_timesRunAlgorithm,
   const std::string& aistr_prefix
   )
  {
    char ls_fileNameOutPlotStatObjetiveFunc[OUTPARAMCLUSTERING_MAX_SIZE_NAME_FILE];
    
    sprintf(ls_fileNameOutPlotStatObjetiveFunc,
            "%s_%s_%02d_%02d.dat",
	    aipstr_fileNamePlot,
	    aistr_prefix.c_str(),
            aii_timesRunAlgorithm
            );

    s_fileNameOutPlotStatObjetiveFunc.assign(ls_fileNameOutPlotStatObjetiveFunc);
  }
  
  inline const std::string getFileNameOutPlotStatObjetiveFunc()
  {
    return s_fileNameOutPlotStatObjetiveFunc;
  }
  */

  /*
  inline void setNumClusterK(T_CLUSTERIDX aiidxK_numClusterK) 
  {
    this->_idxK_numClusterK = aiidxK_numClusterK;
  }

  inline const T_CLUSTERIDX getNumClusterK() 
  {
    return this->_idxK_numClusterK; 
  }

  inline void setNumClusterKTest(T_CLUSTERIDX aiidxK_numClusterKTest) 
  {
    this->_idxK_numClusterKTest = aiidxK_numClusterKTest;
  }

  inline const T_CLUSTERIDX getNumClusterKTest() 
  {
    return this->_idxK_numClusterKTest; 
  }
  */
 
  inline void setPercentageSensitivity(std::string aistr_percentageSensitivity) 
  {
    _str_percentageSensitivity = aistr_percentageSensitivity;
  }

  inline const std::string& getPercentageSensitivity() 
  {
    return _str_percentageSensitivity;
  }

  inline void setPercentageSensitivityTest(std::string aistr_percentageSensitivityTest) 
  {
    _str_percentageSensitivityTest = aistr_percentageSensitivityTest;
  }

  inline const std::string& getPercentageSensitivityTest() 
  {
    return _str_percentageSensitivityTest;
  }

  
  inline void setPercentageSpecificity(std::string aistr_percentageSpecificity) 
  {
    _str_percentageSpecificity = aistr_percentageSpecificity;
  }

  inline const std::string& getPercentageSpecificity() 
  {
    return _str_percentageSpecificity;
  }

  inline void setPercentageSpecificityTest(std::string aistr_percentageSpecificityTest) 
  {
    _str_percentageSpecificityTest = aistr_percentageSpecificityTest;
  }

  inline const std::string& getPercentageSpecificityTest() 
  {
    return _str_percentageSpecificityTest;
  }

  
  virtual void  print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {
    
    // OutParamClustering::print(aipf_outFile,aic_separator);
    /* aipf_outFile << aic_separator << "_file name plot funtion hist";
    if ( s_fileNameOutPlotStatObjetiveFunc.length() > 0 ) {  
      aipf_outFile << aic_separator << this->s_fileNameOutPlotStatObjetiveFunc << ".plot";
    }
    else {
      aipf_outFile << aic_separator << "none";
    }
    aipf_outFile << aic_separator << "_metric used" 
		 << aic_separator << garray_nameObjetiveFunc[_enum_usedObjectiveFunc];
    */
    uintidx  luintidx_numObjetiveFunc  = this->getNumObjetiveFunc();
    for (uintidx lui_i = 0; lui_i < luintidx_numObjetiveFunc; ++lui_i) {
      aipf_outFile <<  aic_separator << "_" << garray_nameObjetiveFunc[lui_i]
		   << aic_separator << this->_arrayrt_metricFunc[lui_i];
    }

    aipf_outFile << aic_separator << "_sensitivity"   
		 << aic_separator << _str_percentageSensitivity;
    
    aipf_outFile << aic_separator << "_specificity"   
		 << aic_separator <<  _str_percentageSpecificity;

    for (uintidx lui_i = 0; lui_i < luintidx_numObjetiveFunc; ++lui_i) {
      aipf_outFile << aic_separator << ":" << garray_nameObjetiveFunc[lui_i]
		   << aic_separator << this->_arrayrt_metricFuncTest[lui_i];
    }

    /* aipf_outFile << aic_separator << "_iterations need" 
		 << aic_separator << this->_it_iterationGetsBest;
    aipf_outFile << aic_separator << "_time seconds need best" 
		 << aic_separator << this->et_runTimeGetsBest;
    aipf_outFile << aic_separator << "_outK"   
		 << aic_separator << this->_idxK_numClusterK;
    aipf_outFile << aic_separator << ":outK"   
		 << aic_separator << this->_idxK_numClusterKTest;
   
    aipf_outFile << aic_separator << "_solution overridden in the run of the algorithm" 
		 << aic_separator << this->_l_totalInvalidOffspring;
    */    

    aipf_outFile << aic_separator << ":sensitivity"   
		 << aic_separator << _str_percentageSensitivityTest;
    
    aipf_outFile << aic_separator << ":specificity"   
		 << aic_separator <<  _str_percentageSpecificityTest;

    /* aipf_outFile << aic_separator << "_objetivefuncrun"
		  << aic_separator << _rt_objetiveFuncRun;
    */
   
  }

  friend std::ostream& operator<<(std::ostream& os, const OutParamClusteringMetric<T_METRIC>  &aiom_outParamMetric)
  {
    aiom_outParamMetric.print(os,OutFileName::getDelim());
    
    return os;
  }
  
protected:
  /*sum of squared error (SSE) metric (clustering metric)*/
  /*within group sum of squares (WGSS) or minimum variance payoff 
   J_1:M_c x L_c --> R \cite{Bezdek:etal:GAclustering:GA:1994} */
  //OutParamNameObjectiveFunc _enum_usedObjectiveFunc;
 
  //T_METRIC                  _rt_objetiveFuncRun;
  T_METRIC                  *_arrayrt_metricFunc;
  T_METRIC                  *_arrayrt_metricFuncTest;
  //COMMON_IDOMAIN           _it_iterationGetsBest;
  //runtime::ExecutionTime    et_runTimeGetsBest;
  /*-1, 0, 1, .., K, K ALGORITHMS SEEK TO OTHERWISE
    THAT K EQUAL THE INPUT OUTPUT
  */
  //T_CLUSTERIDX             _idxK_numClusterK; 
  //T_CLUSTERIDX             _idxK_numClusterKTest;
  //std::string              s_fileNameOutPlotStatObjetiveFunc;  
  
  //long                      _l_totalInvalidOffspring;
  
  std::string               _str_percentageSensitivity;
  std::string               _str_percentageSensitivityTest;
  std::string               _str_percentageSpecificity;
  std::string               _str_percentageSpecificityTest;
  
}; /*END clas OutParamClusteringMetric*/

  
} /*END namespace inout
   */

#endif /*__OUT_PARAM_CLUSTERING_HPP*/
