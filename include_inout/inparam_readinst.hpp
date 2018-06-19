/*! \file inparam_readinst.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __IN_PARAM_READINST_HPP
#define __IN_PARAM_READINST_HPP

#include <string>
#include <vector>
#include "common.hpp"


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

#define INPARAM_FILE_SEPARATOR_LENGTH  20
#define INPARAM_FILE_SEPARATOR_DEFAULT  ","

typedef enum { INPARAM_FORMATINSTANCEFILE_UCI=0, 
	       INPARAM_FORMATINSTANCEFILE_KEEL=1} 
  EnumFormatInstanceFile;

#define INPARAMCLUSTERING_FORMATINSTANCEFILE {"uci", "keel", (char *) NULL }

#define INPARAMCLUSTERING_KEEL_FILE_SEPARATOR ", "

/*! \class InParamReadInst
  \brief Input parameter for read instances
*/
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX 
	  > 
class InParamReadInst {
public:
  InParamReadInst():
    _ui_idxCurrentFileInstance(0),
    _vectorstr_filesInstance(0),
    _vectorstr_filesInstanceTest(0),
    _enum_formatInstanceFile(INPARAM_FORMATINSTANCEFILE_UCI),
    _b_haveHeaderFile(false),
    _s_separatorAttributes(INPARAM_FILE_SEPARATOR_DEFAULT),
    _ps_selectAttributes(),
    _ui_classInstanceColumn(0),
    _ui_clusterInstanceColumn(0),
    _ui_idInstanceColumn(0),
    _ui_instancesFrequencyColumn(0),
    _ui_idMultiInstanceColumn(0),
    _ui_classMultiInstColumn(0),
    _ui_numInstances(0),
    _ui_numDimensionsInstances(0)
  {}
  
  InParamReadInst
  (const char *aips_fileNameInstance,
   EnumFormatInstanceFile aienum_formatInstanceFile,
   const char *ais_separatorAttributes,
   const char *aips_selectAttributes,
   bool aib_haveHeaderFile = false
   ):
    _ui_idxCurrentFileInstance(0),
    _vectorstr_filesInstance(0),
    _vectorstr_filesInstanceTest(0),
    _enum_formatInstanceFile(aienum_formatInstanceFile),
    _b_haveHeaderFile(aib_haveHeaderFile),
    _s_separatorAttributes(ais_separatorAttributes),
    _ps_selectAttributes(aips_selectAttributes),
    _ui_classInstanceColumn(0),
    _ui_clusterInstanceColumn(0),
    _ui_idInstanceColumn(0),
    _ui_instancesFrequencyColumn(0),
    _ui_idMultiInstanceColumn(0),
    _ui_classMultiInstColumn(0),
    _ui_numInstances(0),
    _ui_numDimensionsInstances(0),
    _ui_numInstancesTest(0)
  {
    _vectorstr_filesInstance.push_back( std::string(aips_fileNameInstance) ); 
  }

  ~InParamReadInst() {}

  inline void clearVectorFilesInstance()
  {
    _vectorstr_filesInstance.clear();
    _vectorstr_filesInstanceTest.clear();
  }

  const std::vector<std::string>& getVectorFilesInstance() const
  {
    return _vectorstr_filesInstance;
  }

  inline 
  void
   setVectorFilesInstance
  (const std::vector<std::string> &aivectorstr_filesInstance) 
  {
    _vectorstr_filesInstance = aivectorstr_filesInstance;
  }

  inline 
  uintidx getNumFilesInstance() 
  {
    return (uintidx)  _vectorstr_filesInstance.size();
  }

  inline
  void 
  setVectorFilesInstanceTest
  (const std::vector<std::string> &aivectorstr_filesInstanceTest)
  { 
    _vectorstr_filesInstanceTest = aivectorstr_filesInstanceTest;
  }

  const std::vector<std::string>& getVectorFilesInstanceTest() const
  { 
    return _vectorstr_filesInstanceTest;
  }

  uintidx getNumFilesInstanceTest()
  { 
    return (uintidx) _vectorstr_filesInstanceTest.size();
  }

  inline void setFileInstance(const char *aips_nameFile)
  {
    this->_vectorstr_filesInstance.push_back(std::string(aips_nameFile)); 
  }

  inline void setFileInstance(const std::string aistr_nameFile)
  {
    this->_vectorstr_filesInstance.push_back(aistr_nameFile); 
  }

  inline void setFileInstanceTest(const char *aips_nameFileTest)
  {
    this->_vectorstr_filesInstanceTest.push_back(std::string(aips_nameFileTest)); 
  }

  inline void setFileInstanceTest(const std::string aistr_nameFileTest)
  {
    this->_vectorstr_filesInstanceTest.push_back(aistr_nameFileTest); 
  }

  inline const std::string& getCurrentFileInstance() const
  { 
    return this->_vectorstr_filesInstance.at(this->_ui_idxCurrentFileInstance); 
  }

  inline const std::string& getCurrentFileInstanceTest() const
  { 
    return this->_vectorstr_filesInstanceTest.at(this->_ui_idxCurrentFileInstance); 
  }

  inline void setCurrentFileInstance(uintidx aiui_i)
  {
    this->_ui_idxCurrentFileInstance = aiui_i;
  }

  inline void setSelectAttributes(const char *aips_selectAttributes) 
  {
    this->_ps_selectAttributes = aips_selectAttributes;
  }

  inline const std::string& getSelectAttributes()
  { 
    return this->_ps_selectAttributes; 
  }

  inline void setSeparateAttributes(const char* aipc_separateAttributes)
  {
    this->_s_separatorAttributes.assign(aipc_separateAttributes);
  }

  inline const std::string& getSeparateAttributes()
  { 
    return this->_s_separatorAttributes; 
  }

  inline char getSeparateAttributesChar()
  { 
    return this->_s_separatorAttributes.at(0); 
  }
  
  inline void setHaveHeaderFileInstance(bool aib_bool)
  {
    this->_b_haveHeaderFile = aib_bool;
  }

  inline bool getHaveHeaderFileInstance()
  { 
    return this->_b_haveHeaderFile; 
  }

  inline void setClassInstanceColumn(uintidx aiui_classInstanceColumn) 
  {
    this->_ui_classInstanceColumn = aiui_classInstanceColumn;
  }

  inline uintidx getClassInstanceColumn()
  { 
    return this->_ui_classInstanceColumn; 
  }


  inline void setClusterInstanceColumn(uintidx aiui_clusterInstanceColumn) 
  {
    this-> _ui_clusterInstanceColumn = aiui_clusterInstanceColumn;
  }

  inline uintidx getClusterInstanceColumn()
  { 
    return this->_ui_clusterInstanceColumn; 
  }


  inline void setIDInstanceColumn(uintidx aiui_idInstanceColumn)
  {
    this->_ui_idInstanceColumn = aiui_idInstanceColumn;
  }

  inline uintidx getIDInstanceColumn()
  {
    return this->_ui_idInstanceColumn;
  }

  inline void setInstanceFrequencyColumn(uintidx aiui_instancesFrequencyColumn)
  {
    this->_ui_instancesFrequencyColumn = aiui_instancesFrequencyColumn;
  }

  inline uintidx getInstanceFrequencyColumn()
  {
    return this->_ui_instancesFrequencyColumn;
  }

  inline void setFormatInstanceFile
  (EnumFormatInstanceFile aienum_formatInstanceFile)
  {
    _enum_formatInstanceFile = aienum_formatInstanceFile;
  }

  inline EnumFormatInstanceFile getFormatInstanceFile()
  {
   return _enum_formatInstanceFile;
  }

  inline void setIDMultiInstanceColumn(uintidx aiui_idMultiInstanceColumn)
  {
    this->_ui_idMultiInstanceColumn = aiui_idMultiInstanceColumn;
  }

  inline uintidx getIDMultiInstanceColumn()
  {
    return this->_ui_idMultiInstanceColumn;
  }

  inline void setClassMultiInstColumn(uintidx aiui_classMultiInstColumn)
  {
    this->_ui_classMultiInstColumn = aiui_classMultiInstColumn;
  }

  inline uintidx getClassMultiInstColumn()
  {
    return this->_ui_classMultiInstColumn;
  }

  inline void setNumInstances(uintidx aiui_numInstances) 
  {
    _ui_numInstances = aiui_numInstances;
  }

  inline const uintidx getNumInstances() const
  {
    return _ui_numInstances;
  }


  inline void setNumInstancesTest(uintidx aiui_numInstancesTest) 
  {
    _ui_numInstancesTest = aiui_numInstancesTest;
  }

  inline const uintidx getNumInstancesTest() const
  {
    return _ui_numInstancesTest;
  }
  
  inline void setNumDimensionsInstances(uintidx aiui_numDimensionsInstances) 
  {
    _ui_numDimensionsInstances = aiui_numDimensionsInstances;
  }

  inline const uintidx getNumDimensionsInstances() const 
  {
    return _ui_numDimensionsInstances;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {
    const char  *larray_opFormatFile[] = INPARAMCLUSTERING_FORMATINSTANCEFILE;

    aipf_outFile << aic_separator << "_format file" 
		 << aic_separator << larray_opFormatFile[this->_enum_formatInstanceFile];

    aipf_outFile << aic_separator << "_dataset" 
		 << aic_separator << this->getCurrentFileInstance();
    aipf_outFile << aic_separator << "_n" 
		 << aic_separator << this->getNumInstances();
    aipf_outFile << aic_separator << "_d" 
		 << aic_separator << this->getNumDimensionsInstances();
    aipf_outFile << aic_separator << "_select-attributes" 
		 << aic_separator << _ps_selectAttributes;
    aipf_outFile << aic_separator << "_class-column" 
		 << aic_separator << _ui_classInstanceColumn;
    
    if ( this->_vectorstr_filesInstanceTest.size() > 0 ) {
      aipf_outFile << aic_separator << ":dataset test" 
		   << aic_separator << this->getCurrentFileInstanceTest();
      aipf_outFile << aic_separator << ":n" 
		   << aic_separator << this->getNumInstancesTest();
    }

  }

  
protected:

  uintidx                  _ui_idxCurrentFileInstance;
  std::vector<std::string> _vectorstr_filesInstance;
  std::vector<std::string> _vectorstr_filesInstanceTest;
  EnumFormatInstanceFile  _enum_formatInstanceFile;
  bool         _b_haveHeaderFile;
  std::string  _s_separatorAttributes;
  std::string  _ps_selectAttributes;
  uintidx      _ui_classInstanceColumn;
  uintidx      _ui_clusterInstanceColumn;
  uintidx      _ui_idInstanceColumn;
  uintidx      _ui_instancesFrequencyColumn;
  uintidx      _ui_idMultiInstanceColumn;
  uintidx      _ui_classMultiInstColumn;
  
  uintidx      _ui_numInstances;
  uintidx      _ui_numDimensionsInstances;
  uintidx      _ui_numInstancesTest;

}; /*InParamReadInst*/

/*! \class InParamReadInstFreq
  \brief Input parameter for read instances with frequency
*/
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY
	  > 
class InParamReadInstFreq 
  : public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> {
public:
  InParamReadInstFreq()
    : InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
  {}
  InParamReadInstFreq
  (const char *aips_fileNameInstance,
   EnumFormatInstanceFile aienum_formatInstanceFile,
   const char *ais_separatorAttributes,
   const char *aips_selectAttributes,
   bool aib_haveHeaderFile = false
   ):
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
    (aips_fileNameInstance,
     aienum_formatInstanceFile,
     ais_separatorAttributes,
     aips_selectAttributes,
     aib_haveHeaderFile
     )
  {}
};

} /*END namespace inout 
   */

#endif /*__IN_PARAM_READINST_HPP*/
