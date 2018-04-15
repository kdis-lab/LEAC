/*! \file instances_read.hpp
 *
 * \brief instances read
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef INSTANCES_READ_HPP
#define INSTANCES_READ_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cstring>
#include <utility>      // std::pair, std::make_pair
#include <stddef.h>
#include <string.h>
#include "instance_class.hpp"
#include "instance_frequency.hpp"
#include "instance_class_frequency.hpp"

#include "vector_utils.hpp"
#include "container_out.hpp"
#include "inparam_readinst.hpp"

#include "verbose_global.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
    
/*! \fn std::pair<uintidx,uintidx> readNumInstances(inout::InParamReadInst  &aiipri_inParamReadInst, const bool aib_fileTest = false) 
    \brief Read the number of instances or objects 
    \details 
    \param aiipri_inParamReadInst a inout::InParamReadInst with the necessary parameters to read a data set file
    \param aib_fileTest a bool to specify if the data set is a test
 */ 
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX
	  > 
std::pair<uintidx,uintidx>
readNumInstances
(inout::InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>  &aiipri_inParamReadInst,
 const    bool             aib_fileTest = false
 ) 
{
  const std::string lstr_fileInstance
    ((aib_fileTest == false)?aiipri_inParamReadInst.getCurrentFileInstance():
     aiipri_inParamReadInst.getCurrentFileInstanceTest()
     );   
  std::ifstream lfp_file(lstr_fileInstance);
  if (!lfp_file) {
    std::string lstr_error("readNumInstances: no file input data ");
    lstr_error += lstr_fileInstance;
    throw  std::invalid_argument(lstr_error);
  }

  uintidx        lost_numInstances;  
  uintidx        lost_numDimensions = 0; 
  uintidx        luintidx_countLine;
 
  LineSplit
    lls_lineSplit
    (aiipri_inParamReadInst.getSeparateAttributes(),
     aiipri_inParamReadInst.getSelectAttributes()
     );

#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "readNumInstances  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t input  const std::string &aistr_fileInstance = " 
	      << lstr_fileInstance
	      << "\n\t input  InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  std::string lstr_linedata;
    
  /*COUNT THE NUMBER OF INSTANCES*/
  lost_numInstances = 0; 
  luintidx_countLine = 0;
  /*READ FIRST INSTANCE*/ 
  while ( std::getline(lfp_file, lstr_linedata) )  {
    ++luintidx_countLine;
    if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@'  ||  lstr_linedata.at(0) == '#') ) {
      if ( lls_lineSplit.split(lstr_linedata) > 0 ) {
	++lost_numInstances;
	lost_numDimensions = lls_lineSplit.getNumSelectColumns();
	break;
      }
    }
  }
    
  /*READ THE REST OF THE INSTANCE*/ 
  while ( std::getline(lfp_file, lstr_linedata) )  {
    ++luintidx_countLine;
    if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#' ) ) {
      if (  lls_lineSplit.split(lstr_linedata) > 0 ) {
	++lost_numInstances;
      }
      else {
	std::cerr << "readNumInstances: warning on line " 
		  << luintidx_countLine 
		  << " of the file "
		  << lstr_fileInstance 
		  << " does not take it as an instance\n";
      }
    }
  }

  /*COUNT THE NUMBER OF DIMENSIONS*/
    
  lfp_file.close();
    
  if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) --lost_numInstances;

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "readNumInstances OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output uintidx: lost_numInstances = " << lost_numInstances 
	      << "\tuintidx: lost_numDimensions = " << lost_numDimensions
	      << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return std::make_pair(lost_numInstances,lost_numDimensions);

}

/*! \fn std::vector<std::string>  instancesReadDimName (inout::InParamReadInst &aiipri_inParamReadInst, const bool aib_fileTest = false) 
    \brief Read the name of the dimensions of instances or objects 
    \details 
    \param aiipri_inParamReadInst a inout::InParamReadInst with the necessary parameters to read a data set file
    \param aib_fileTest a bool to specify if the data set is a test
 */
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX 
	  > 
std::vector<std::string>  
instancesReadDimName
(inout::InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiipri_inParamReadInst,
 const    bool             aib_fileTest = false        
 ) 
{
  const std::string lstr_fileInstance
    ((aib_fileTest == false)?aiipri_inParamReadInst.getCurrentFileInstance():
     aiipri_inParamReadInst.getCurrentFileInstanceTest()
     );
  std::ifstream lfp_file(lstr_fileInstance);
  if (!lfp_file) {
    std::string lstr_error("instancesReadDimName: no file input data ");
    lstr_error += lstr_fileInstance;
    throw  std::invalid_argument(lstr_error);
  }
  std::vector<std::string> lovectorstr_dimensionsName;
 
  LineSplit lls_lineSplit;
   
#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadDimName  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t input  InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << "]\n"
	      << "\t)"
	      << std::endl;
  }
#endif //__VERBOSE_YES

  std::string lstr_linedata;
    
  if ( aiipri_inParamReadInst.getFormatInstanceFile()
       == inout::INPARAM_FORMATINSTANCEFILE_KEEL ) {
    std::string lstr_inputs("@inputs");
      
    while ( std::getline(lfp_file, lstr_linedata) )  {
      size_t lstidx_findInput;
      
      if ( (lstidx_findInput = lstr_linedata.find(lstr_inputs)) != std::string::npos ) {
	std::string lstr_header 
	  (lstr_linedata,
	   lstidx_findInput+lstr_inputs.size()+1,
	   lstr_linedata.size()
	   );
	
	lls_lineSplit.setSeparator(aiipri_inParamReadInst.getSeparateAttributes());
	lls_lineSplit.split( lstr_header );
	lls_lineSplit.getVectorString(lovectorstr_dimensionsName);
	break;
      }
    } /*end while*/     
  }
  else if (aiipri_inParamReadInst.getFormatInstanceFile() == inout::INPARAM_FORMATINSTANCEFILE_UCI) {
    //JUMPING FILE COMMENTS
    while ( std::getline(lfp_file, lstr_linedata) )  {
      if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#') )
	break;
    }
    //IF HAVE HEADER READ 
    if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {  
      lls_lineSplit.setSeparator(aiipri_inParamReadInst.getSeparateAttributes());
      lls_lineSplit.getSelectColumns(aiipri_inParamReadInst.getSelectAttributes());
      if ( lls_lineSplit.split( lstr_linedata ) > 0 ) {
	lls_lineSplit.getVectorString(lovectorstr_dimensionsName);
      }
    }
    else {
      for (uintidx list_l = 0; 
	   list_l < lls_lineSplit.getNumSelectColumns(); 
	   list_l++) 
	{
	  std::ostringstream lostrstream_nameDimension;
	  lostrstream_nameDimension << 'x' << list_l;
	  lovectorstr_dimensionsName.push_back(lostrstream_nameDimension.str()); 
	}
    }
  }

  lfp_file.close();
       
#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadDimName OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<std::string>: [" 
	      << &lovectorstr_dimensionsName << ']';
    for (const std::string &listr :lovectorstr_dimensionsName)
      std::cout << ',' << listr;
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectorstr_dimensionsName;  
}

/*! \fn std::vector<data::Instance<T_FEATURE>* > instancesRead (inout::InParamReadInst &aiipri_inParamReadInst, const bool aib_fileTest = false)  
    \brief Read the instances or objects 
    \details 
    \param aiipri_inParamReadInst a inout::InParamReadInst with the necessary parameters to read a data set file
    \param aib_fileTest a bool to specify if the data set is a test
 */
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX 
	  > 
std::vector<data::Instance<T_FEATURE>* > 
instancesRead
(inout::InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX> &aiipri_inParamReadInst,
 const    bool          aib_fileTest = false         
 ) 
{
  const std::string lstr_fileInstance
    ((aib_fileTest == false)?aiipri_inParamReadInst.getCurrentFileInstance():
     aiipri_inParamReadInst.getCurrentFileInstanceTest()
     );
  std::ifstream lfp_file(lstr_fileInstance);
  if (!lfp_file) {
    std::string lstr_error("instancesRead with frequency: no file input data ");
    lstr_error += lstr_fileInstance;
    throw  std::invalid_argument(lstr_error);
  }
  std::vector<data::Instance<T_FEATURE>* >  lovectorptinst_instances;

  LineSplit
    lls_lineSplit
    (aiipri_inParamReadInst.getSeparateAttributes(),
     aiipri_inParamReadInst.getSelectAttributes()
     );

#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesRead  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(output vector<Instance<> >: [" << &lovectorptinst_instances << "]\n"
	      << "n\t input  const std::string &aistr_fileInstance = " 
	      << lstr_fileInstance
	      << "\t input  inout::InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << "]\n"
	      << "\t)\n";
  }
#endif /*__VERBOSE_YES*/
    
  std::pair<uintidx,uintidx> lpair_dimInstance =
    readNumInstances
    (aiipri_inParamReadInst,
     aib_fileTest
     );
  lovectorptinst_instances.reserve(lpair_dimInstance.first); 
  data::Instance<T_FEATURE>::setNumDimensions(lpair_dimInstance.second);
 
  std::string lstr_linedata;
  //JUMPING FILE COMMENTS
  while ( std::getline(lfp_file, lstr_linedata) )  {
    if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#') )
      break;
  }
  /*READ HEADER*/
  if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {
    std::getline(lfp_file, lstr_linedata);
  }
    
  do {
    // BEGIN IF LINE COMMENTS
    if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#')) {
      if ( lls_lineSplit.split( lstr_linedata ) > 0 ) {
	/*GET INSTANCE*/
	data::Instance
	  <T_FEATURE> 
	  *lptinst_new = 
	  new data::Instance
	  <T_FEATURE>();
	  
	lptinst_new->readFeature(lls_lineSplit);
	/*IF BEGIN READ ID*/
	if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { 
	  lptinst_new->setId
	    (lls_lineSplit.getItem( aiipri_inParamReadInst.getIDInstanceColumn() ));
	} /*IF END READ ID*/
	lovectorptinst_instances.push_back(lptinst_new);
      }
    }
  } while ( std::getline(lfp_file, lstr_linedata) );  
 
  lfp_file.close();
        
#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesRead OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<Instance<T_FEATURE>* >: [" 
	      << &lovectorptinst_instances << "]\n";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      for ( auto  literinst_instance: lovectorptinst_instances ) 
	std::cout << *literinst_instance << '\n';
    }
    --geiinparam_verbose;
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lovectorptinst_instances;

} /* instancesRead */


/*! \fn std::vector<data::Instance<T_FEATURE>* > instancesReadWithClass(inout::InParamReadInst   &aiipri_inParamReadInst, const boolaib_fileTest = false)  
    \brief Read the instances or objects with class
    \details 
    \param aiipri_inParamReadInst a inout::InParamReadInst with the necessary parameters to read a data set file
    \param aib_fileTest a bool to specify if the data set is a test
 */
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX 
	  >
std::vector<data::Instance<T_FEATURE>* > 
instancesReadWithClass
(inout::InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>   &aiipri_inParamReadInst,
 const    bool            aib_fileTest = false         
 ) 
{
  const std::string lstr_fileInstance
    ((aib_fileTest == false)?aiipri_inParamReadInst.getCurrentFileInstance():
     aiipri_inParamReadInst.getCurrentFileInstanceTest()
     );
  std::ifstream lfp_file(lstr_fileInstance);
  if (!lfp_file) {
    std::string lstr_error("instancesReadWithClass: no file input data ");
    lstr_error += lstr_fileInstance;
    throw  std::invalid_argument(lstr_error);
  }
 
  std::vector<data::Instance<T_FEATURE>* >  lovectorptinst_instances;

  std::istringstream                      liss_stringstream;
 
  LineSplit
    lls_lineSplit
    (aiipri_inParamReadInst.getSeparateAttributes(),
     aiipri_inParamReadInst.getSelectAttributes()
     );
  
#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadWithClass  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(output vector<InstanceClass<> >: [" << &lovectorptinst_instances << ']'
	      << "\n\t input  const std::string &aistr_fileInstance = " 
	      << lstr_fileInstance
	      << "\n\t input  inout::InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::pair<uintidx,uintidx> lpair_dimInstance =
    readNumInstances
    (aiipri_inParamReadInst,
     aib_fileTest
     );
  lovectorptinst_instances.reserve(lpair_dimInstance.first);
  data::Instance<T_FEATURE>::setNumDimensions(lpair_dimInstance.second);

  std::string lstr_linedata;
  //JUMPING FILE COMMENTS
  while ( std::getline(lfp_file, lstr_linedata) )  {
    if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' || lstr_linedata.at(0) == '#' ) )
      break;
  }
    
  /*READ HEADER*/
  if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {
    std::getline(lfp_file, lstr_linedata);
  }

  /*READ INSTANCES */
    
  do {
    /*BEGIN IF LINE COMMENTS*/
    if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#')) {
      if ( lls_lineSplit.split(lstr_linedata) > 0 ) {

	/*GET INSTANCE*/
	data::InstanceClass
	  <T_FEATURE,
	   T_INSTANCES_CLUSTER_K,
	   T_CLUSTERIDX> 
	  *lptinstclass_new = 
	  new data::InstanceClass
	  <T_FEATURE,
	   T_INSTANCES_CLUSTER_K,
	   T_CLUSTERIDX>
	  ();

	lptinstclass_new->readFeature(lls_lineSplit);
	  
	if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { /*IF BEGIN READ ID*/
	  lptinstclass_new->setId
	    (lls_lineSplit.getItem( aiipri_inParamReadInst.getIDInstanceColumn() ));
	} /*IF END READ ID*/

	if ( aiipri_inParamReadInst.getClassInstanceColumn() ) { /*IF BEGIN CLASS LABEL*/
	   
	  std::string lstr_keyMapClass;

	  if ( aiipri_inParamReadInst.getClusterInstanceColumn() ) {
	    lstr_keyMapClass = 
	      lls_lineSplit.getItem( aiipri_inParamReadInst.getClassInstanceColumn() ) 
	      +"_"
	      + lls_lineSplit.getItem( aiipri_inParamReadInst.getClusterInstanceColumn() ); 
	  } 
	  else {
	    lstr_keyMapClass = 
	      lls_lineSplit.getItem( aiipri_inParamReadInst.getClassInstanceColumn() );
	   
	  }

	  lptinstclass_new->setClassIdx(lstr_keyMapClass);

	}  /*IF END CLASS LABEL*/

	if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { /*IF BEGIN ID*/
	   
	  lptinstclass_new->setId
	    (lls_lineSplit.getItem( aiipri_inParamReadInst.getIDInstanceColumn() ));
	} /*IF END ID*/
	  
	lovectorptinst_instances.push_back(lptinstclass_new);
      }
    }
  } while ( std::getline(lfp_file, lstr_linedata) );  
  
  if ( aiipri_inParamReadInst.getClassInstanceColumn() ) {
    data::InstanceIterfazClass
      <T_INSTANCES_CLUSTER_K,
       T_CLUSTERIDX>
      ::setVectorClassLabel(); 
  }
  
  lfp_file.close();
        
#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadWithClass OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<InstanceClass<T_FEATURE>* >: [" 
	      << &lovectorptinst_instances << "]\n";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      for ( auto  literinst_instance: lovectorptinst_instances ) {
	literinst_instance->print(std::cout,'\t');
	std::cout << '\n';
      }
    }
    --geiinparam_verbose;
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lovectorptinst_instances;

} /*end instancesRead: instances with class*/

/*! \fn std::vector<data::Instance<T_FEATURE>* > instancesReadWithFreq (inout::InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY> &aiipri_inParamReadInst, const bool aib_fileTest)
    \brief Read the instances or objects with frequency
    \details 
    \param aiipri_inParamReadInst a inout::InParamReadInstFreq with the necessary parameters to read a data set file
    \param aib_fileTest a bool to specify if the data set is a test
 */
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY
	  > 
std::vector<data::Instance<T_FEATURE>* > 
instancesReadWithFreq
(inout::InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY> &aiipri_inParamReadInst,
 const    bool            aib_fileTest = false
 ) 
{
  const std::string lstr_fileInstance
    ((aib_fileTest == false)?aiipri_inParamReadInst.getCurrentFileInstance():
     aiipri_inParamReadInst.getCurrentFileInstanceTest()
     );
  std::ifstream lfp_file(lstr_fileInstance);
  if (!lfp_file) {
    std::string lstr_error("instancesReadWithFreq: no file input data ");
    lstr_error += lstr_fileInstance;
    throw  std::invalid_argument(lstr_error);
  }
  std::vector<data::Instance<T_FEATURE>* >  lovectorptinst_instances;
  std::istringstream                  liss_stringstream;
  inout::LineSplit
    lls_lineSplit
    (aiipri_inParamReadInst.getSeparateAttributes(),
     aiipri_inParamReadInst.getSelectAttributes()
     );
  
#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadWithFreq  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(output vector<InstanceFreq<>* >: [" 
	      << &lovectorptinst_instances << ']'
	      << "\n\t input  const std::string &aistr_fileInstance = " 
	      << lstr_fileInstance
	      << "\n\t input  inout::InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::pair<uintidx,uintidx> lpair_dimInstance =
    readNumInstances
    (aiipri_inParamReadInst,
     aib_fileTest
     );
  lovectorptinst_instances.reserve(lpair_dimInstance.first); 
  data::Instance<T_FEATURE>::setNumDimensions(lpair_dimInstance.second);
  
  /*OPEN FILE INSTANCES*/
    std::string lstr_linedata;
    /*JUMPING FILE COMMENTS*/
    while ( std::getline(lfp_file, lstr_linedata) )  {
      if ( (lstr_linedata.size() > 0) &&
	   !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#') )
	break;
    }
    /*READ HEADER*/
    if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {
      std::getline(lfp_file, lstr_linedata);
    }

    /*READ INSTANCES */
    T_INSTANCE_FREQUENCY lt_frequency = 1;

    do {
      /*BEGIN IF LINE COMMENTS*/
      if ( (lstr_linedata.size() > 0) &&
	   !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#')) {
	if ( lls_lineSplit.split(lstr_linedata) > 0 ) {
	  /*GET INSTANCE*/
	  if ( aiipri_inParamReadInst.getInstanceFrequencyColumn() ) 
	    { /*IF BEGIN INSTANCES FREQUENCY*/
	    liss_stringstream.clear();
	    liss_stringstream.str
	      ( lls_lineSplit.getItem( aiipri_inParamReadInst.getInstanceFrequencyColumn() ));
	    liss_stringstream >> lt_frequency; 
	    } /*IF END INSTANCES FREQUENCY*/

	  data::InstanceFreq
	    <T_FEATURE,
	     T_INSTANCE_FREQUENCY
	     > 
	    *lptinst_new = 
	    new data::InstanceFreq
	    <T_FEATURE,
	     T_INSTANCE_FREQUENCY
	     >
	    (lt_frequency);
	  lptinst_new->readFeature(lls_lineSplit);
	  if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { /*IF BEGIN READ ID*/
	    lptinst_new->setId
	      (lls_lineSplit.getItem( aiipri_inParamReadInst.getIDInstanceColumn() ));
	  } /*IF END READ ID*/
	  lovectorptinst_instances.push_back(lptinst_new);
	}
      }
    } while ( std::getline(lfp_file, lstr_linedata) );  
  
    lfp_file.close();
        
#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadWithFreq OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<InstanceFreq<>* >: [" 
	      << &lovectorptinst_instances << "]\n";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      for ( auto  literinst_instance: lovectorptinst_instances ) {
	literinst_instance->print(std::cout,'\t');
	std::cout << '\n';
      }
    }
    --geiinparam_verbose;
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lovectorptinst_instances;

} /*END instancesReadWithFreq
   */


/*! \fn std::vector<data::Instance<T_FEATURE>* > instancesReadWithFreqClass (inout::InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY> &aiipri_inParamReadInst, const bool aib_fileTest)
    \brief Read the instances or objects with frequency
    \details 
    \param aiipri_inParamReadInst a inout::InParamReadInstFreq with the necessary parameters to read a data set file
    \param aib_fileTest a bool to specify if the data set is a test
 */
template <typename T_FEATURE,         
	  typename T_INSTANCES_CLUSTER_K,
	  typename T_CLUSTERIDX,
	  typename T_INSTANCE_FREQUENCY
	  > 
std::vector<data::Instance<T_FEATURE>* > 
instancesReadWithFreqClass
(inout::InParamReadInstFreq<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX,T_INSTANCE_FREQUENCY> &aiipri_inParamReadInst,
 const    bool            aib_fileTest=false        
 ) 
{
  const std::string lstr_fileInstance
    ((aib_fileTest == false)?aiipri_inParamReadInst.getCurrentFileInstance():
     aiipri_inParamReadInst.getCurrentFileInstanceTest()
     );
  std::ifstream lfp_file(lstr_fileInstance);
  if (!lfp_file) {
    std::string lstr_error("instancesReadWithFreqClass: no file input data ");
    lstr_error += lstr_fileInstance;
    throw  std::invalid_argument(lstr_error);
  }
  std::vector<data::Instance<T_FEATURE>* >  lovectorptinst_instances;
  std::istringstream                        liss_stringstream;
  LineSplit
    lls_lineSplit
    (aiipri_inParamReadInst.getSeparateAttributes(),
     aiipri_inParamReadInst.getSelectAttributes()
     );
  
#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadWithFreqClass  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n\t(output vector<InstanceClassFreq<>* >: [" << &lovectorptinst_instances << ']'
	      << "\n\t input  const std::string &aistr_fileInstance = " 
	      << lstr_fileInstance
	      << "\n\t input  inout::InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << ']'
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/

  std::pair<uintidx,uintidx> lpair_dimInstance =
    readNumInstances
    (aiipri_inParamReadInst,
     aib_fileTest
     );
  lovectorptinst_instances.reserve(lpair_dimInstance.first); 
  data::Instance<T_FEATURE>::setNumDimensions(lpair_dimInstance.second);
  /*OPEN FILE INSTANCES */
    std::string lstr_linedata;
    /*JUMPING FILE COMMENTS*/
    while ( std::getline(lfp_file, lstr_linedata) )  {
      if ( (lstr_linedata.size() > 0) &&
	   !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#') )
	break;
    }
    /*READ HEADER*/
    if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {
      std::getline(lfp_file, lstr_linedata);
    }
    /*READ INSTANCES */
    T_INSTANCE_FREQUENCY lt_frequency = 1;
    do {
      /*BEGIN IF LINE COMMENTS #*/
      if ( (lstr_linedata.size() > 0) &&
	   !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#')) {
	if ( lls_lineSplit.split(lstr_linedata) > 0 ) {
	  
	  /*GET INSTANCE*/	  
	  data::InstanceClassFreq
	    <T_FEATURE,
	     T_INSTANCE_FREQUENCY,
	     T_INSTANCES_CLUSTER_K,
	     T_CLUSTERIDX> 
	    *lptinstclass_new = 
	    new data::InstanceClassFreq
	    <T_FEATURE,
	     T_INSTANCE_FREQUENCY,
	     T_INSTANCES_CLUSTER_K,
	     T_CLUSTERIDX>
	    ();
	  lptinstclass_new->readFeature(lls_lineSplit);  

	  if ( aiipri_inParamReadInst.getIDInstanceColumn() ) {
	    /*IF BEGIN READ ID*/
	    lptinstclass_new->setId
	      (lls_lineSplit.getItem( aiipri_inParamReadInst.getIDInstanceColumn() ));
	  } /*IF END READ ID*/

	  if ( aiipri_inParamReadInst.getInstanceFrequencyColumn() ) 
	    { /*IF BEGIN INSTANCES FREQUENCY*/
	    liss_stringstream.clear();
	    liss_stringstream.str
	      ( lls_lineSplit.getItem( aiipri_inParamReadInst.getInstanceFrequencyColumn() ));
	    liss_stringstream >> lt_frequency;
	    lptinstclass_new->setFrequency(lt_frequency);
	    } /*IF END INSTANCES FREQUENCY*/

	  if ( aiipri_inParamReadInst.getClassInstanceColumn() ) {
	    /*IF BEGIN CLASS LABEL*/
	    std::string lstr_keyMapClass;
	    if ( aiipri_inParamReadInst.getClusterInstanceColumn() ) {
	      lstr_keyMapClass = 
		lls_lineSplit.getItem( aiipri_inParamReadInst.getClassInstanceColumn() ) 
		+"_"
		+ lls_lineSplit.getItem( aiipri_inParamReadInst.getClusterInstanceColumn() ); 
	    }	    
	    else {
	      lstr_keyMapClass =
		lls_lineSplit.getItem( aiipri_inParamReadInst.getClassInstanceColumn() );
	    }
	    lptinstclass_new->setClassIdx(lstr_keyMapClass);
       
	  } /* IF END CLASS LABEL*/

	  lovectorptinst_instances.push_back(lptinstclass_new);
	}
      }
    } while ( std::getline(lfp_file, lstr_linedata) );  
  
    if ( aiipri_inParamReadInst.getClassInstanceColumn() ) {
      data::InstanceIterfazClass
	<T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::setVectorClassLabel(); 
    }
    lfp_file.close();
        
#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadWithFreqClass OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<InstanceClassFreq<T_FEATURE>* >: [" 
	      << &lovectorptinst_instances << "]\n";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      for ( auto  literinst_instance: lovectorptinst_instances ) {
	literinst_instance->print(std::cout,'\t');
	std::cout << '\n';
      }
    }
    --geiinparam_verbose;
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lovectorptinst_instances;

} /*END instancesReadWithFreqClass*/



} /*END namespace inout
   */

#endif /*INSTANCES_READ_HPP*/
