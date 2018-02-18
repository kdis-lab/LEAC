/*! \file outfilename.hpp
 *
 * \brief out file name
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __OUT_FILE_NAME_HPP
#define __OUT_FILE_NAME_HPP

#include <iostream>
#include <fstream> //std::ofstream
#include <string.h>


/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {

  
#define OUTFILENAME_STDOUT             "stdout"
#define OUTFILENAME_SEPARATOR_DEFAULT  ','

class OutFileName {
public:
  OutFileName() 
    :  _stdstremamsize_precisionDefault(std::cout.precision())
  {}
  ~OutFileName() {} 

  std::ostream& openFile
  (const char* aiptc_outfileName, 
   const std::ios_base::openmode ai_openmode = std::ios::out | std::ios::app
   ) 
  {
    const char      lps_fileStdout[] = OUTFILENAME_STDOUT;

    if ( (aiptc_outfileName == NULL) || (strcmp(lps_fileStdout,aiptc_outfileName) == 0) ) {
      _ptostream_ostream = &std::cout;
    }
    else {
      _ofstream_filename.open(aiptc_outfileName, ai_openmode );
      _ptostream_ostream = &_ofstream_filename;
    }
  
    return *_ptostream_ostream;

  } 

  void closeFile() 
  {        
    if (_ptostream_ostream != NULL) {
      _ofstream_filename.precision(_stdstremamsize_precisionDefault);
      if ( _ptostream_ostream != &std::cout ) {
	_ofstream_filename.close();
      }
    }
  }

  static void setDelim(const char aic_delim)
  {
    OutFileName::_c_delim = aic_delim;
  }

  static  const char getDelim()  
  {
    return OutFileName::_c_delim;
  }

protected:

  std::streamsize  _stdstremamsize_precisionDefault;
  std::ostream     *_ptostream_ostream;
  std::ofstream    _ofstream_filename;

  static char      _c_delim;

}; /*OutFileName*/


} /*END namespace  outparam
   */


#endif /*__OUT_FILE_NAME_HPP*/
