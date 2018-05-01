/*! \file inparam.hpp
 *
 * \brief Definition of input parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
 
#ifndef IN_PARAM_HPP
#define IN_PARAM_HPP

#include <cstddef>
#include <time.h>
#include "outfilename.hpp"
#include "getsubopt.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
typedef enum {LABEL, CENTROIDS, CRISPPARTITION, MEDOIDS, DENSITY_BASED,HIERARCHICAL,OTHER} InParam_algTypeOut;

#define INPARAM_SIZE_STRING_ID_TIME 50

/*! \class InParamAlgorithmo
  \brief Input parameter for define propertys of the algorithmos
*/
  class InParamAlgorithmo { 
public:
  InParamAlgorithmo
  (std::string        ais_algorithmoName, 
   std::string        ais_algorithmoAuthor, 
   InParam_algTypeOut aiato_algTypeOut
   ):
    s_algorithmoName(ais_algorithmoName), 
    s_algorithmoAuthor(ais_algorithmoAuthor), 
    ato_algTypeOut(aiato_algTypeOut),
 
    ps_fileNameTimesRun(NULL),
    b_progressBarPrinting(false),
    i_timesRunAlgorithm(1),
    s_gnuplotCoreStyles("linespoints")
 
  {
    time_t     ltt_runTime       = time(NULL);
    struct tm  *ltm_localrunTime = localtime(&ltt_runTime);
    char       ls_runTimeId[INPARAM_SIZE_STRING_ID_TIME];

    sprintf(ls_runTimeId,
	    "%04d%02d%02d:%02d%02d%02d",
	    ltm_localrunTime->tm_year+1900,
	    ltm_localrunTime->tm_mon+1,
	    ltm_localrunTime->tm_mday,
	    ltm_localrunTime->tm_hour,
	    ltm_localrunTime->tm_min,
	    ltm_localrunTime->tm_sec
	    );
    _s_runTimeId.assign(ls_runTimeId);
  }

  ~InParamAlgorithmo() {}

  void errorArgument(char *apc_command, const char *apc_options, const char *aas_arguments[])
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
  std::cerr << "\nTry " << apc_command << "' --help for more information." << std::endl;
  exit(-1);                           
  }
   

  bool isYesNo(char  *optarg, char *apc_command, const char *apc_options)
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
  
  inline void setFileNameTimesRun(char *aips_nameFile) 
  {
    this->ps_fileNameTimesRun = aips_nameFile;
  }


  inline void setProgressBarPrinting(bool aib_bool)   
  {
    this->b_progressBarPrinting = aib_bool;
  }

  inline void setTimesRunAlgorithm(int aii_timesRun)
  {
    this->i_timesRunAlgorithm = aii_timesRun;
  }

  inline void setGnuPlotCoreStyles(const char* aips_nameCoreStyle)
  {
    this->s_gnuplotCoreStyles.assign(aips_nameCoreStyle);
  }

  inline std::string getAlgorithmoAuthor() const
  { 
    return this->s_algorithmoAuthor; 
  }

  inline std::string getAlgorithmoName() const
  { 
    return this->s_algorithmoName; 
  }

  inline char* getFileNameTimesRun() const
  { 
    return this->ps_fileNameTimesRun; 
  }


  inline int getTimesRunAlgorithm() const
  { 
    return this->i_timesRunAlgorithm; 
  }

  inline bool getProgressBarPrinting()
  {
    return this->b_progressBarPrinting;
  }

  inline const std::string getGnuPlotCoreStyles() const
  { 
    return this->s_gnuplotCoreStyles; 
  }

  inline InParam_algTypeOut getAlgorithTypeOut()
  { 
    return this->ato_algTypeOut; 
  }

  inline const std::string& getRunningTimeId()
  {
    return this->_s_runTimeId;
  }

  virtual void print(std::ostream&  aipf_outFile=std::cout, const char aic_separator=',') const
  {  
    aipf_outFile << "_algorithmo" 
		 << aic_separator << this->s_algorithmoName;
    aipf_outFile << aic_separator << "_author" 
		 << aic_separator << this->s_algorithmoAuthor;
 
    aipf_outFile << aic_separator << "_times run" 
		 << aic_separator << this->i_timesRunAlgorithm
		 << aic_separator << "_runnig date"
		 << aic_separator << this->_s_runTimeId;
  }

  
  friend std::ostream& operator<<(std::ostream& os, const InParamAlgorithmo  &aiinparamalg_inParam)
  {
    aiinparamalg_inParam.print(os);
    
    return os;
  }

protected:

  std::string        s_algorithmoName;
  std::string        s_algorithmoAuthor;
  InParam_algTypeOut ato_algTypeOut;
   
  char               *ps_fileNameTimesRun;
  bool               b_progressBarPrinting;
  int                i_timesRunAlgorithm;
  std::string        _s_runTimeId;
  std::string        s_gnuplotCoreStyles;
  
}; /* InParamAlgorithmo */

} /*END namespace inout 
   */

#endif /*IN_PARAM_HPP*/
