/*! \file inparam.cpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#include <iostream>
#include <cstdlib>
#include "inparam.hpp"
#include "getsubopt.h"

char inout::OutFileName::_c_delim = OUTFILENAME_SEPARATOR_DEFAULT;

void inout::InParamAlgorithmo::errorArgument
(char       *apc_command, 
 const char *apc_options, 
 const char *aas_arguments[])
{
  int li_i;

  li_i = 0;
  std::cerr << apc_command <<": ambiguous argument for `--" << apc_options << "'\n";
  std::cerr << "Valid arguments are:\n";
  while ( aas_arguments[li_i] != NULL)  {
    std::cerr << "`" << aas_arguments[li_i] << "' ";
    ++li_i;
    if  (aas_arguments[li_i] != NULL)
      std::cerr << ", ";
  }  
  std::cerr << "\nTry `ls --" << apc_command << "' for more information." << std::endl;
  exit(-1);                           
}


bool inout::InParamAlgorithmo::isYesNo(char  *optarg, char *apc_command, const char *apc_options) 
{
  const char *las_optTokensYesNo[] = {"no", "yes", (char *) NULL };
  char *lps_optsubValue = NULL;
  bool lob_yesNo = true;
  int  li_idxSubOpt;

  if (optarg) {
    if ( (li_idxSubOpt = getsubopt_getsubopt(&optarg, las_optTokensYesNo, &lps_optsubValue)) != -1) {
      lob_yesNo = li_idxSubOpt;
    }
    else {
      this->errorArgument(apc_command, apc_options, las_optTokensYesNo);
    }	
  }
  
  return lob_yesNo;    
}

