/*! \file file_utils.hpp
 *
 * \brief file utils
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __FILE_UTILS_HPP
#define __FILE_UTILS_HPP

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  

std::string fileutils_subName
(const std::string aistr_name, 
 const std::string aistr_nameSub
)
{
  std::string lostr_newName;

  size_t lst_posPath = aistr_name.find_last_of("/\\");
  if ( lst_posPath > 0 )
    lostr_newName.append(aistr_name.substr(0,lst_posPath+1));
  
  std::string lstr_filename = aistr_name.substr(lst_posPath+1);

  size_t lst_posPoint = lstr_filename.find_last_of(".");
  if ( lst_posPoint == std::string::npos ) {
    lst_posPoint = lstr_filename.size();
  }
  
  lstr_filename.insert(lst_posPoint,aistr_nameSub);
  
  lostr_newName.append(lstr_filename);

  return lostr_newName;
}

bool fileutils_isdir(const char *aiptstr_dirname)
{
  struct stat lstruct_stat;

  if (aiptstr_dirname == NULL )
    return false;
  stat(aiptstr_dirname,&lstruct_stat);

  return S_ISDIR(lstruct_stat.st_mode)?true:false;
}

std::vector<std::string> 
fileutils_listdir(const char *aiptstr_dirname)
{
  DIR *lptdir_dir;
  struct dirent *lstruct_dirent;
  std::vector<std::string> lovector_directoryFile(0);
  lovector_directoryFile.reserve(10);

  lptdir_dir = opendir(aiptstr_dirname);
  if ( lptdir_dir == NULL ) { 
    throw std::invalid_argument("fileutils_listdir: directory does not exist");
  }
  std::string lstr_dirname(aiptstr_dirname);
  while ( (lstruct_dirent = readdir(lptdir_dir)) != NULL ) {
    std::string lstr_file = std::string(lstruct_dirent->d_name);
    if ( (lstr_file.compare(".") != 0) && (lstr_file.compare("..") != 0) )  
      lovector_directoryFile.push_back(lstr_dirname + "/" +lstr_file);
  }
  
  std::sort(lovector_directoryFile.begin(),lovector_directoryFile.end());
  
  closedir(lptdir_dir);

  return lovector_directoryFile;
}

} /*END namespace inout 
   */

#endif /*__FILE_UTILS_HPP*/
