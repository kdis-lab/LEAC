/*! \file line_split.hpp
 *
 * \brief line split
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __LINE_SPLIT_HPP
#define __LINE_SPLIT_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <string.h>
#include "common.hpp"

#include "verbose_global.hpp"

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
/*! \class LineSplit
  \brief Used for processing files to read
*/
class LineSplit {
public:
  LineSplit()
    : _str_separator()
    , _vectorst_selectColumns(0)
    , _vextorptc_columns(0)
  { }
  
  LineSplit
  (std::string  aistr_separator, 
   std::string  aistr_selectColumns /*If NULL is the maximum*/
   )
    : _str_separator(aistr_separator)
    , _vectorst_selectColumns(0)
    , _vextorptc_columns(0)
  {
    this->getSelectColumns(aistr_selectColumns);
  }

  ~LineSplit()
  {
  }

  inline 
  void setSeparator(const std::string &aistr_separator) 
  {
    _str_separator = aistr_separator;
  }

  void  getSelectColumns(std::string aistr_selectColumns)
  {
    char   *lpc_readItem;
    char   *lpc_selectColumns;
    uintidx luintidx_limInf;
    uintidx luintidx_limSup;
    char   ls_selectColumns[ aistr_selectColumns.length() + 2];
  
    using namespace std;
    istringstream   liss_stringstream;

    _vectorst_selectColumns.clear();
    if ( !aistr_selectColumns.empty() ) {
      strcpy (ls_selectColumns, aistr_selectColumns.c_str());
      lpc_readItem = strtok(ls_selectColumns,",");
      while (lpc_readItem != NULL ) {
	lpc_selectColumns = strchr(lpc_readItem,'-');    
	if (lpc_selectColumns == NULL) {
	  liss_stringstream.clear();
	  liss_stringstream.str(lpc_readItem);
	  liss_stringstream >> luintidx_limInf;
	  _vectorst_selectColumns.push_back( luintidx_limInf );
	}
	else {
	  *lpc_selectColumns = '\0';
	  liss_stringstream.clear();
	  liss_stringstream.str(lpc_readItem);
	  liss_stringstream >> luintidx_limInf;
	  liss_stringstream.clear();
	  liss_stringstream.str(lpc_selectColumns+1);
	  liss_stringstream >> luintidx_limSup;
	  for (uintidx luintidx_i = luintidx_limInf;
	       luintidx_i <= luintidx_limSup; luintidx_i++) {
	    _vectorst_selectColumns.push_back(luintidx_i);
	  }
	}
	lpc_readItem = strtok(NULL, ",");
      }	
    } /*IF !NULL*/
  }

  uintidx split(std::string& str) 
  {
    uintidx luintidx_start, luintidx_end = 0;

    str.erase(str.find_last_not_of("\n\r")+1);

    _vextorptc_columns.clear();
    while (luintidx_end < str.size()) {
      luintidx_start = luintidx_end;
      while (luintidx_start < str.size() &&
	     (_str_separator.find(str[luintidx_start]) != std::string::npos)) {
	luintidx_start += _str_separator.size();  // skip initial whitespace
      }
      luintidx_end = luintidx_start;
      while (luintidx_end < str.size() &&
	     (_str_separator.find(str[luintidx_end]) == std::string::npos)) {

        luintidx_end++; // skip to end of word
      }
      if (luintidx_end-luintidx_start != 0) {  // just ignore zero-length strings.
	_vextorptc_columns.push_back
	  (std::string(str, luintidx_start, luintidx_end-luintidx_start));
      }
    }
    //InitializeConsecutive 0..N
    if (_vectorst_selectColumns.size() == 0) {
      uintidx li_i = 1;
      for ( auto  liter_select: _vextorptc_columns ) 
	_vectorst_selectColumns.push_back(li_i++);
    }
    
    return (uintidx) _vextorptc_columns.size();

  }
  
  uintidx getNumSelectColumns() 
  {
    return (uintidx) _vectorst_selectColumns.size();
  }
    
  inline std::string& getItemSelect(uintidx aist_i)
  {
    return _vextorptc_columns.at(_vectorst_selectColumns.at(aist_i)-1);
  }

  inline std::string& getItem(uintidx aist_i)
  {
    return _vextorptc_columns.at(aist_i-1);
  }

  void  getVectorString(std::vector<std::string> &aovectorstr_item)
  {
    std::string  ls_columnName;
    for (uintidx luintidx_j = 0;
	 luintidx_j < _vectorst_selectColumns.size(); luintidx_j++) {
      ls_columnName.assign
	(_vextorptc_columns.at(_vectorst_selectColumns.at(luintidx_j)-1));
      aovectorstr_item.push_back(ls_columnName);
    }
  }

  inline std::vector<std::string>& getVectorRawString()
  {
    return _vextorptc_columns;
  }
  
protected:

  std::string              _str_separator;
  std::vector<uintidx>     _vectorst_selectColumns;
  std::vector<std::string> _vextorptc_columns;

}; 

} /*END namespace inout
   */


#endif /* __LINE_SPLIT_HPP */
