/*! \file instances_classfrequency_read.hpp
 *
 * \brief instances class frequency read
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef INSTANCES_CLASSFREQUENCY_READ_HPP
#define INSTANCES_CLASSFREQUENCY_READ_HPP

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
#include "instance_frequency.hpp"
#include "instance_class_frequency.hpp"
#include "vector_utils.hpp"
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
    
/*instancesReadWithFreq: instances with frequency
 */
template < typename T_FEATURE,
	   typename T_INSTANCE_FREQUENCY
	   >
std::vector<data::Instance<T_FEATURE>* > 
instancesReadWithFreq
(inout::InParamReadInst   &aiipri_inParamReadInst,
 const    bool            aib_fileTest         
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

} /*instancesRead: instances with frequency*/



/*instancesReadWithFreqClass: instances with frequency and class
 */
template < typename T_FEATURE,
	   typename T_INSTANCE_FREQUENCY,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	   >
std::vector<data::Instance<T_FEATURE>* > 
instancesReadWithFreqClass
(inout::InParamReadInst   &aiipri_inParamReadInst,
 const    bool            aib_fileTest         
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

} /*instancesReadWithFreqClass: instances with frequency and  class*/


} /*END namespace data 
   */

#endif /*INSTANCES_CLASSFREQUENCY_READ_HPP*/
