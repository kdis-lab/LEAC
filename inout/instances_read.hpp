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
#ifdef __INSTANCES_WITH_FREQUENCY
#include "instance_frequency.hpp"
#include "instance_class_frequency.hpp"
#else
#include "instance.hpp"
#include "instance_class.hpp"
#endif /*__INSTANCES_WITH_FREQUENCY*/

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
    
/*readNumInstances:  
 */ 
std::pair<uintidx,uintidx>
readNumInstances
(inout::InParamReadInst  &aiipri_inParamReadInst,
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


std::vector<std::string>  
instancesReadDimName
(inout::InParamReadInst &aiipri_inParamReadInst,
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
      
      //lls_lineSplit.setSeparator(aiipri_inParamReadInst.getSeparateAttributes());
      //lls_lineSplit.getSelectColumns(aiipri_inParamReadInst.getSelectAttributes());
      while ( std::getline(lfp_file, lstr_linedata) )  {
	size_t lstidx_findInput;
	//std::cout << "instancesReadDimName: " << lstr_linedata << std::endl;
	if ( (lstidx_findInput = lstr_linedata.find(lstr_inputs)) != std::string::npos ) {
	  std::string lstr_header 
	    (lstr_linedata,
	     lstidx_findInput+lstr_inputs.size()+1,
	     lstr_linedata.size()
	     );
	  //std::cout << "instancesReadDimName found: " << lstr_header << std::endl;
	  lls_lineSplit.setSeparator(aiipri_inParamReadInst.getSeparateAttributes());
	  lls_lineSplit.split( lstr_header );
	  lls_lineSplit.getVectorString(lovectorstr_dimensionsName);
	  break;
        }
	/*if ( lstr_linedata.size() > 0 ) && lls_lineSplit.split( lstr_linedata ) > 0 ) {
	  std::cout << "instancesReadDimName: " << lstr_linedata << std::endl;
	  //lls_lineSplit.getVectorString(lovectorstr_dimensionsName);
	  if ( lstr_inputs.compare(lls_lineSplit.getItem((uintidx) 1)) == 0 ) {
	    for (uintidx list_l = 2; 
		 list_l < lls_lineSplit.getNumSelectColumns(); 
		 list_l++) 
	      {  
		lovectorstr_dimensionsName.push_back(lls_lineSplit.getItem(list_l)); 
	      }
	    break;
	  }
	  }*/
      } /*end while*/     
    }
    else if (aiipri_inParamReadInst.getFormatInstanceFile() == inout::INPARAM_FORMATINSTANCEFILE_UCI) {
      //JUMPING FILE COMMENTS
      while ( std::getline(lfp_file, lstr_linedata) )  {
	if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#') )
	  break;
      }
      // SI TIENE ENCABEZADO LEER ENCABEZADO 
      if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {  
	lls_lineSplit.setSeparator(aiipri_inParamReadInst.getSeparateAttributes());
	lls_lineSplit.getSelectColumns(aiipri_inParamReadInst.getSelectAttributes());
	if ( lls_lineSplit.split( lstr_linedata ) > 0 ) {
	  lls_lineSplit.getVectorString(lovectorstr_dimensionsName);
	}
      }
      else {
	//char  ls_nameDimension[20];
	for (uintidx list_l = 0; 
	     list_l < lls_lineSplit.getNumSelectColumns(); 
	     list_l++) 
	  {
	     std::ostringstream lostrstream_nameDimension;
	     lostrstream_nameDimension << 'x' << list_l;
	     //sprintf(ls_nameDimension,"x%zu",list_l);
	    //ls_columnName.assign(ls_linedata);
	     lovectorstr_dimensionsName.push_back(lostrstream_nameDimension.str()); //std::string(ls_nameDimension)); 
	  }
      }
    }

    lfp_file.close();
        
    /*}
  catch (std::ios_base::failure &fail) {
    std::string 
      lstr_error
      ("instancesReadDimName: no file input data ");
    lstr_error += aistr_fileInstance;
    //lstr_error += aiipri_inParamReadInst.getFileNameInstance();
    throw  std::invalid_argument(lstr_error);
    }*/

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadDimName OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<std::string>: [" 
	      << &lovectorstr_dimensionsName << ']'
	      << '\n'
	      << lovectorstr_dimensionsName
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectorstr_dimensionsName;  
}


/*std::vector<std::string>  
instancesReadHeader
(InParamReadInst &aiipri_inParamReadInst) 
{
  std::vector<std::string> lovectorstr_dimensionsName;
  std::ifstream            lfp_file;
 
  LineSplit
    lls_lineSplit
    (aiipri_inParamReadInst.getSeparateAttributes(),
     aiipri_inParamReadInst.getSelectAttributes()
     );

#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadHeader  IN"
	      << '(' << geiinparam_verbose << ")\n"
	      << "\t input  InParamReadInst: &aiipri_inParamReadInst[" 
	      << &aiipri_inParamReadInst << "]\n"
	      << "\t)\n";
  }
#endif //__VERBOSE_YES

  try {
    lfp_file.open(aiipri_inParamReadInst.getFileNameInstance().c_str());
    std::string lstr_linedata;
    //JUMPING FILE COMMENTS
    while ( std::getline(lfp_file, lstr_linedata) )  {
      if ( (lstr_linedata.size() > 0) && (lstr_linedata.at(0) == '@') )
	break;
    }
    // SI TIENE ENCABEZADO LEER ENCABEZADO 
    if ( aiipri_inParamReadInst.getHaveHeaderFileInstance() ) {  
      if ( lls_lineSplit.split( lstr_linedata ) > 0 ) {
	if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { //IF BEGIN ID
	  lovectorstr_dimensionsName.push_back
	    (lls_lineSplit.getItem(aiipri_inParamReadInst.getIDInstanceColumn()));
	}
	lls_lineSplit.getVectorString(lovectorstr_dimensionsName);
	if ( aiipri_inParamReadInst.getClassInstanceColumn() ) {  // IF BEGIN CLASS LABEL
	  lovectorstr_dimensionsName.push_back
	    (lls_lineSplit.getItem(aiipri_inParamReadInst.getClassInstanceColumn()));
	}
	if (aiipri_inParamReadInst.getInstanceFrequencyColumn()) { //IF BEGIN INSTANCES FREQUENCY
	  lovectorstr_dimensionsName.push_back
	    (lls_lineSplit.getItem(aiipri_inParamReadInst.getInstanceFrequencyColumn()));
	}
      }
    }
    else {
      if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { //IF BEGIN ID
	lovectorstr_dimensionsName.push_back(std::string("id")); 
      }
      char  ls_nameDimension[20];
      for (uintidx list_l = 0; list_l < lls_lineSplit.getNumSelectColumns(); list_l++) {  
	sprintf(ls_nameDimension,"x%zu",list_l);
	lovectorstr_dimensionsName.push_back(std::string(ls_nameDimension)); 
      }
      if ( aiipri_inParamReadInst.getClassInstanceColumn() ) { //IF BEGIN ID
	lovectorstr_dimensionsName.push_back(std::string("class")); 
      }
      if (aiipri_inParamReadInst.getInstanceFrequencyColumn()) { //IF BEGIN INSTANCES FREQUENCY
	lovectorstr_dimensionsName.push_back(std::string("frequency")); 
       }
    }
 
    lfp_file.close();
        
  }
  catch (std::ios_base::failure &fail) {
    std::string lstr_error("instancesReadHeader: no file input data ");
    lstr_error += aiipri_inParamReadInst.getFileNameInstance();
    throw  std::invalid_argument(lstr_error);
  }

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "instancesReadHeader OUT"
	      << '(' << geiinparam_verbose << ")\n"
	      << "output std::vector<std::string>: [" 
	      << &lovectorstr_dimensionsName << ']';
    std::cout << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES

  return lovectorstr_dimensionsName;  
  }*/


 
/*instancesRead
  return std::vector<Instance<T_FEATURE>* > 
*/

#ifdef __INSTANCES_WITH_FREQUENCY
template < typename T_FEATURE,
	   typename T_INSTANCE_FREQUENCY
	   >
#else
template < typename T_FEATURE
	   >
#endif /*__INSTANCES_WITH_FREQUENCY*/
std::vector<data::Instance<T_FEATURE>* > 
instancesRead
(inout::InParamReadInst &aiipri_inParamReadInst,
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
  //std::vector<Instance<T_FEATURE,T_INSTANCE_FREQUENCY>* >  lovectorptinst_instances;

  //std::ifstream            lfp_file;
  //lfp_file.exceptions(std::ios::failbit);
  //using namespace std;
  //istringstream   liss_stringstream;
  //FILE            *lfp_file;
  //char            ls_linedata[INSTANCESREAD_TEXTLENG];
  //uintidx          luintidx_numInstance;
  //uintidx          luintidx_numDimensions;
  //std::string     ls_columnName;


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
    
  //uintidx luintidx_numInstance =
  std::pair<uintidx,uintidx> lpair_dimInstance =
    readNumInstances
    (//aiipri_inParamReadInst.getCurrentFileInstance(),
     aiipri_inParamReadInst,
     aib_fileTest
     );
  lovectorptinst_instances.reserve(lpair_dimInstance.first); 
  data::Instance<T_FEATURE>::setNumDimensions(lpair_dimInstance.second);
  
  //try {
  // lfp_file.open(aistr_fileInstance.c_str());
    //lfp_file.open(aiptstr_fileInstance);
    //lfp_file.open(aiipri_inParamReadInst.getFileNameInstance().c_str());
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
  
    /*READ INSTANCES */
    /*Instance
      <T_FEATURE,
       T_INSTANCE_FREQUENCY>
       ::setNumDimensions(0);*/
    
    do {
      //std::cout << "begin Procesa :" << lstr_linedata << std::endl;
      if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#')) {// BEGIN IF LINE COMMENTS
	if ( lls_lineSplit.split( lstr_linedata ) > 0 ) {
	  /*GET INSTANCE*/
	  //if ( Instance<T_FEATURE,T_INSTANCE_FREQUENCY>::getNumDimensions() == (uintidx) 0 )
	  /*if ( Instance<T_FEATURE>::getNumDimensions() == (uintidx) 0 )
	    Instance<T_FEATURE>::setNumDimensions
	      (lls_lineSplit.getNumSelectColumns());
	  */
	  //std::cout << " NumDimensions = " << Instance<T_FEATURE>::getNumDimensions() << std::endl;
#ifdef __INSTANCES_WITH_FREQUENCY
	  //std::cout << "read InstanceFreqOne  se creo\n";
	  data::InstanceFreq<T_FEATURE,T_INSTANCE_FREQUENCY> 
	    *lptinst_new = 
	    new data::InstanceFreq
	    <T_FEATURE,T_INSTANCE_FREQUENCY>();
#else
	  //std::cout << "read Instance  se creo\n";
	  data::Instance
	    <T_FEATURE> 
	    *lptinst_new = 
	    new data::Instance
	    <T_FEATURE>();
#endif /*__INSTANCES_WITH_FREQUENCY*/
	  
	  lptinst_new->readFeature(lls_lineSplit);
	  if ( aiipri_inParamReadInst.getIDInstanceColumn() ) { /*IF BEGIN READ ID*/
	    lptinst_new->setId
	      (lls_lineSplit.getItem( aiipri_inParamReadInst.getIDInstanceColumn() ));
	  } /*IF END READ ID*/
	  lovectorptinst_instances.push_back(lptinst_new);
	}
      }
      //std::cout << "end Procesa :" << lstr_linedata << std::endl;
    } while ( std::getline(lfp_file, lstr_linedata) );  
 
    lfp_file.close();
        
    /*}
  catch (std::ios_base::failure &fail) {
    std::string lstr_error("instancesRead: no file input data ");
    lstr_error += aistr_fileInstance;
    //lstr_error += aiipri_inParamReadInst.getFileNameInstance();
    throw  std::invalid_argument(lstr_error);
  }
    */

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
	//std::cout << literinst_instance->getToString() << '\n';
    }
    --geiinparam_verbose;
    std::cout << std::endl;
    //this->print();
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return lovectorptinst_instances;

} /* instancesRead */


/*instancesRead: instances with class
 */
#ifdef __INSTANCES_WITH_FREQUENCY
template < typename T_FEATURE,
	   typename T_INSTANCE_FREQUENCY,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
#else
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
#endif /*__INSTANCES_WITH_FREQUENCY*/
std::vector<data::Instance<T_FEATURE>* > 
instancesReadWithClass
(inout::InParamReadInst   &aiipri_inParamReadInst,
 const    bool              aib_fileTest = false         
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
      if ( (lstr_linedata.size() > 0) && !(lstr_linedata.at(0) == '@' ||  lstr_linedata.at(0) == '#')) {/*BEGIN IF LINE COMMENTS*/
	if ( lls_lineSplit.split(lstr_linedata) > 0 ) {

	  /*GET INSTANCE*/
	 
#ifdef __INSTANCES_WITH_FREQUENCY

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

#else
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

#endif /*__INSTANCES_WITH_FREQUENCY*/

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


} /*END namespace inout
   */

#endif /*INSTANCES_READ_HPP*/
