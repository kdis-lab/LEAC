/*! \file main_stdvar_milligan_cooper1988.cpp
 *
 * \brief Main program for the standardization of instances \cite milligan:cooper:clusteranalysis:1988
 *
 * \details This file is part of the LEAC.\n\n
 * Implementation of the standardization of variables based on the paper:\n
 * G.W. Milligan and M.C. Cooper. A study of standardization\n 
 * of variables in cluster analysis. J. Classification,\n
 * 5:181–204, 1988.\n
 * <a href="http://www.springerlink.com/content/t588424722r23031">doi:http://www.springerlink.com/content/t588424722r23031</a>\n.
 * \n
  * Library Evolutionary Algorithms for Clustering (LEAC) is a library\n
 * for the implementation of evolutionary and genetic algorithms\n
 * focused on the partition type clustering problem. Based on the\n
 * current standards of the <a href="http://en.cppreference.com">C++</a> language, as well as on Standard\n
 * Template Library <a href="http://en.cppreference.com/w/cpp/container">STL</a> 
 * and also  <a href="http://www.openblas.net/">OpenBLAS</a> to have a better performance.\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#include "inparamclustering_getparameter.hpp"
#include "instances_read.hpp"
#include "datatype_instance_real.hpp"
#include "outfilename.hpp"
#include "file_utils.hpp"
#include "stats_instances.hpp"
#include "linear_algebra_level1.hpp"
#include "linear_algebra_level2.hpp"

#include "verbose.hpp"

/*---< main() >-------------------------------------------------------------
*/
int main(int argc, char **argv) 
{
#ifdef __VERBOSE_YES
   const char* lpc_labeMain = "main_stdvar_milligan_cooper1988.cpp";
   geverbosepc_labelstep = lpc_labeMain;
#endif /*__VERBOSE_YES*/
                    
  inout::InParamStdVar
    listdvar_inParam
    ("STDVAR",
     "Milligan and Cooper 1988",
     inout::OTHER
     );
  listdvar_inParam.setStandardizationVar(STD_VAR_Z1);
  
  /*READ PARAMETER
   */
  inparamclustering_getParameter
    (listdvar_inParam, argc, argv);

#ifdef __VERBOSE_YES
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "main:  IN" 
	      << '(' << geiinparam_verbose << ')'
	      << "\nargv: ";
    std::cout << argv[0] << ' '; 
    for (int li_i = 1; li_i < argc;  li_i++ ) {
      std::cout << ' '  << argv[li_i];
    } 
    std::cout << "\n";
  }
#endif /*__VERBOSE_YES*/  

  
  /*READ: INSTANCES 
   */
  std::vector<data::Instance<DATATYPE_FEATURE>* >  lvectorptinst_instances;

  if ( listdvar_inParam.getClassInstanceColumn() ) {
      
    data::InstanceIterfazClass
      <DATATYPE_INSTANCES_CLUSTER_K,DATATYPE_CLUSTERIDX>
      ::initialize();

    lvectorptinst_instances =
      inout::instancesReadWithClass
      <DATATYPE_FEATURE,
       DATATYPE_INSTANCES_CLUSTER_K,
       DATATYPE_CLUSTERIDX>
      (listdvar_inParam,
       false
       );  
  }
  else { /*Instances*/

    lvectorptinst_instances =
      inout::instancesRead
      <DATATYPE_FEATURE>
      (listdvar_inParam,
       false
       );
    
  }

  std::vector<std::string> 
    lvectorstr_instanceDimName = 
    inout::instancesReadDimName
    (listdvar_inParam,
     false
     );

  /*STATISTICS OF INSTANCES
   */
  DATATYPE_FEATURE_SUM* larray_sumFeature =
    new DATATYPE_FEATURE_SUM[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];

  DATATYPE_FEATURE *larray_meanFeactures =
    new DATATYPE_FEATURE[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];

  DATATYPE_FEATURE_SUM* larray_sumFeatureSQ =
    new DATATYPE_FEATURE_SUM[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];
    
  DATATYPE_FEATURE *larray_desvstdFeactures =
    new DATATYPE_FEATURE[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];
                    
  DATATYPE_FEATURE *larray_minFeactures =
    new DATATYPE_FEATURE[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];

  DATATYPE_FEATURE *larray_maxFeactures =
    new DATATYPE_FEATURE[data::Instance<DATATYPE_FEATURE>::getNumDimensions()];
      
  /*STANDARDIZATION OF VARIABLES
   */
  if ( lvectorptinst_instances.size() > 0 ) {
    
    switch ( listdvar_inParam.getStandardizationVar() ) {
    case STD_VAR_Z0:
      break;
    case STD_VAR_Z1:

      /*SUM INSTANCES
       */
      stats::sumFeactures
	(larray_sumFeature,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 DATATYPE_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_meanFeactures,
	 (uintidx) lvectorptinst_instances.size(),
	 larray_sumFeature
	 );

      stats::sumFeacturesSQ
	(larray_sumFeatureSQ,
	 larray_meanFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );

      stats::meanVector
	(larray_desvstdFeactures,
	 (uintidx) lvectorptinst_instances.size()-1,
	 larray_sumFeatureSQ
	 );

      stats::varTodesvstd(larray_desvstdFeactures);
  
      std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures(); 
	   interfacesse::axpy
	     (lt_feactures,
	      DATATYPE_FEATURE(-1),
	      larray_meanFeactures,
	      data::Instance<DATATYPE_FEATURE>::getNumDimensions()
	      );
	   for (uintidx lui_i = 0;
		lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
		++lui_i)
	     {
	       lt_feactures[lui_i] /= larray_desvstdFeactures[lui_i];
	     }
	 }
	 );			
      break;
    case STD_VAR_Z2:

      /*SUM INSTANCES
       */
      stats::sumFeactures
	(larray_sumFeature,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 DATATYPE_FEATURE(0)
	 );
  
      stats::meanVector
	(larray_meanFeactures,
	 (uintidx) lvectorptinst_instances.size(),
	 larray_sumFeature
	 );

      stats::sumFeacturesSQ
	(larray_sumFeatureSQ,
	 larray_meanFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );

      stats::meanVector
	(larray_desvstdFeactures,
	 (uintidx) lvectorptinst_instances.size()-1,
	 larray_sumFeatureSQ
	 );

      stats::varTodesvstd(larray_desvstdFeactures);
  
      std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   for (uintidx lui_i = 0;
		lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
		++lui_i)
	     {
	       lt_feactures[lui_i] /= larray_desvstdFeactures[lui_i];
	     }
	 }
	 );
      break;
    case STD_VAR_Z3:

      stats::maxFeatures
	(larray_maxFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );
      std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   for (uintidx lui_i = 0;
		lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
		++lui_i)
	     {
	       lt_feactures[lui_i] /= larray_maxFeactures[lui_i];
	     }
	 }
	 );
       
      break;
      
    case STD_VAR_Z4:

      stats::maxFeatures
	(larray_maxFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );

      stats::minFeatures
	(larray_minFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );
      
      std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   for (uintidx lui_i = 0;
		lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
		++lui_i)
	     {
	       lt_feactures[lui_i] /= (larray_maxFeactures[lui_i] - larray_minFeactures[lui_i]);
	     }
	 }
	 );
      
      break;
    case STD_VAR_Z5:
      stats::maxFeatures
	(larray_maxFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );

      stats::minFeatures
	(larray_minFeactures,
	 lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end()
	 );
      
      std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   for (uintidx lui_i = 0;
		lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
		++lui_i)
	     {
	       lt_feactures[lui_i] = (lt_feactures[lui_i] -larray_minFeactures[lui_i]) /
		 (larray_maxFeactures[lui_i] - larray_minFeactures[lui_i]);
	     }
	 }
	 );
      break;
    case STD_VAR_Z6:
       std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   for (uintidx lui_i = 0;
		lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
		++lui_i)
	     {
	       lt_feactures[lui_i] /= larray_sumFeature[lui_i];
	     }
	 }
	 );
      break;
    case STD_VAR_Z7:
    {
      std::vector<std::pair<uintidx,DATATYPE_FEATURE> >
	lvector_idxFeature(lvectorptinst_instances.size());
      uintidx lui_consecutive = 0;
      for ( auto & lpair_idxFeacture: lvector_idxFeature ) {
	lpair_idxFeacture.first = lui_consecutive++;
      }

      for (uintidx lui_i = 0;
	   lui_i < data::Instance<DATATYPE_FEATURE>::getNumDimensions();
	   lui_i++) {

	uintidx lui_j = 0;
        std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   lvector_idxFeature[lui_j++].second =  lt_feactures[lui_i];
	 }
	 );
	
	stats::rank(lvector_idxFeature);

	lui_j = 0;
        std::for_each
	(lvectorptinst_instances.begin(),
	 lvectorptinst_instances.end(),
	 [&](data::Instance<DATATYPE_FEATURE>* linst_iter)
	 {
	   DATATYPE_FEATURE* lt_feactures =  linst_iter->getFeatures();
	   lt_feactures[lui_i] = lvector_idxFeature[lui_j++].second; 
	 }
	 );
      }      
      break;
    }
    default:
      throw  std::invalid_argument("main_gas_clustering: undefined standardization of variables");
      break;
    }
  }

  delete [] larray_meanFeactures;
  delete [] larray_desvstdFeactures;
  delete [] larray_sumFeature;
  delete [] larray_sumFeatureSQ;
  delete [] larray_minFeactures;
  delete [] larray_maxFeactures;
  
  /*PRINT PARAMETERS*/
  inout::OutFileName   lofn_filename;
  std::ostream& lostream_out = 
    lofn_filename.openFile(listdvar_inParam.getFileNameTimesRun());
  lostream_out.precision(COMMON_COUT_PRECISION);
  lostream_out << std::boolalpha;

  if (listdvar_inParam.getHaveHeaderFileInstance()) {
    for (uintidx li_l = 0; li_l < lvectorstr_instanceDimName.size()-1; li_l++) {
      lostream_out << lvectorstr_instanceDimName.at(li_l) << lofn_filename.getDelim();  
    }
    lostream_out
      << lvectorstr_instanceDimName.at(lvectorstr_instanceDimName.size()-1)
      << std::endl;
  }
 
  for ( auto liIter_inst: lvectorptinst_instances ) {
    liIter_inst->print(lostream_out,inout::OutFileName::getDelim());
    lostream_out << '\n';
  }

  lofn_filename.closeFile();

#ifdef __VERBOSE_YES
  geverbosepc_labelstep = lpc_labeMain;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "main: OUT" 
	      << '(' << geiinparam_verbose << ")\n"; 
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
     
  return 0;  
}
