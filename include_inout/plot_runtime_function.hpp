/*! \file plot_runtime_function.hpp
 *
 * \brief plot runtime function
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef PLOT_FUNCTION_HIST_HPP
#define PLOT_FUNCTION_HIST_HPP


#include <vector>
#include <string>

#include "inparam_clustering.hpp"
#include "list_runtime_function.hpp"
#include "outparam_clustering.hpp"

#include "runtime_statistic.hpp"

/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace runtime {
  

/*! \fn void plot_funtionHist(runtime::ListRuntimeFunction<COMMON_IDOMAIN> &ailfhT_listFuntionHist, const inparam::InParamClustering &aiipc_inParamClustering, OutParamClustering<T_METRIC,T_CLUSTERIDX>  &aiopcT_outParamClustering) 
  \brief Display function of metrics
  \details 
  \param ailfhT_listFuntionHist a runtime::ListRuntimeFunction<COMMON_IDOMAIN>
  \param aiipc_inParamClustering a const inparam::InParamClustering
  \param aiopcT_outParamClustering a OutParamClustering<T_METRIC,T_CLUSTERIDX>
*/
template < typename T_METRIC,
	   typename T_CLUSTERIDX
	  >
void 
plot_funtionHist
(runtime::ListRuntimeFunction
 <COMMON_IDOMAIN>                 &ailfhT_listFuntionHist,
 const inout::InParamClustering &aiipc_inParamClustering,
 inout::OutParamClustering
 <T_METRIC,T_CLUSTERIDX>          &aiopcT_outParamClustering
 )
{

  std::ofstream lf_outPlot;
  std::string   ls_fileNameOutPlot(aiipc_inParamClustering.getFileNamePlotStatObjetiveFunc());
  static int    lsi_pointtype  = 2;
  std::vector<RuntimeFunction* >& lv_listPtrRuntimeFunction 
    = ailfhT_listFuntionHist.getListRuntimeFunction();

  ls_fileNameOutPlot.append(".plot");
  std::ifstream lifs_fileExists(ls_fileNameOutPlot.c_str());
  
  lf_outPlot.open(ls_fileNameOutPlot.c_str(), std::ios::out | std::ios::app );
  if (lifs_fileExists) {
    lf_outPlot << ",\\\n";
  }
  else {
    lf_outPlot << "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n";
    lf_outPlot << "set ylabel '" << ailfhT_listFuntionHist.getLabelY() << "'\n";
    lf_outPlot << "set xlabel '" << ailfhT_listFuntionHist.getLabelX() << "'\n";
    lf_outPlot << "set grid\n\n";
    
    lf_outPlot << "plot \\\n";
    
  }
 
  for (uintidx lui_i=0; lui_i < lv_listPtrRuntimeFunction.size(); lui_i++) {  
    lf_outPlot << "'" << aiopcT_outParamClustering.getFileNameOutPlotStatObjetiveFunc() << "'" 
	       << " using  1:" << lui_i + 2 
	       << " with " << aiipc_inParamClustering.getGnuPlotCoreStyles();
    if ( aiipc_inParamClustering.getGnuPlotCoreStyles().compare("linespoints") == 0 ) {
      lf_outPlot << " pointsize 1 lt -1 pt " << lsi_pointtype++;  
    } 
    else if ( aiipc_inParamClustering.getGnuPlotCoreStyles().compare("point") == 0 ) { 
      lf_outPlot << " pointsize 1 lt -1 pt " << lsi_pointtype++;  
    }
    lf_outPlot << " title '" 
	       << lv_listPtrRuntimeFunction.at(lui_i)->getName() 
	       << '_' << aiopcT_outParamClustering.getNumRunningAlgorithm()  
	       << '_' << aiipc_inParamClustering.getTimesRunAlgorithm() << "'";

    if  (lui_i < (lv_listPtrRuntimeFunction.size() - 1) )
      lf_outPlot << ",\\\n";
  }

  lf_outPlot.close();

}


} /*END namespace executiontime 
   */

#endif  /* PLOT_FUNCTION_HIST_HPP */


