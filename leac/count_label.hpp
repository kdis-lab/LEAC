/*! \file count_label.hpp
 *
 * \brief count label
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef _COUNT_LABEL_HPP
#define _COUNT_LABEL_HPP 

#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

/*! \namespace data
  \brief Module for the handling of instances or also called objects or points.
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class CountLabel
  \brief Stores the number of occurrences of a label for an attribute in a dataset
*/
template < typename T_IDX_LABEL, 
	   typename T_COUNT_LABEL
	   >
class CountLabel {
public:
  
  CountLabel() :
    _str_label(), 
    _idx_label(UNKNOWN_CLUSTER_IDX), 
    _t_coutlabel(0)  
  {
  }
    
  CountLabel
  (const std::string  &aistr_label, 
   const T_IDX_LABEL  aiidx_label
   ): 
    _str_label(aistr_label), 
    _idx_label(aiidx_label), 
    _t_coutlabel(1)  
  {
  }

  //copy constructor
  CountLabel
  (const CountLabel
   <T_IDX_LABEL,
    T_COUNT_LABEL>    &aicl_classLabel
   ): 
     _str_label(aicl_classLabel._str_label), 
     _idx_label(aicl_classLabel._idx_label), 
     _t_coutlabel(aicl_classLabel._t_coutlabel)  
  {
    
  }
  
  //move constructor 
  CountLabel
  (CountLabel
   <T_IDX_LABEL,
    T_COUNT_LABEL>    &&aicl_classLabel
   ): 
     _str_label(aicl_classLabel._str_label), 
    _idx_label(aicl_classLabel._idx_label), 
    _t_coutlabel(aicl_classLabel._t_coutlabel)  
  {
    aicl_classLabel._idx_label = UNKNOWN_CLUSTER_IDX; 
    aicl_classLabel._t_coutlabel = 0;
  }

  ~CountLabel() 
  {
  }


  CountLabel<T_IDX_LABEL,T_COUNT_LABEL>&
  operator=
  (CountLabel<T_IDX_LABEL,
    T_COUNT_LABEL>    &&aicl_classLabel
   )  
  {
    if ( this != &aicl_classLabel ) {
      this->_str_label = aicl_classLabel._str_label; 
      this->_idx_label = aicl_classLabel._idx_label; 
      this->_t_coutlabel = aicl_classLabel._t_coutlabel;

      aicl_classLabel._idx_label = UNKNOWN_CLUSTER_IDX; 
      aicl_classLabel._t_coutlabel = 0;
    }

    return *this;
    
  }


  CountLabel<T_IDX_LABEL,T_COUNT_LABEL>&
  operator=
  (const CountLabel<T_IDX_LABEL,
    T_COUNT_LABEL>    &aicl_classLabel
   )  
  {
    if ( this != &aicl_classLabel ) {
      this->_str_label = aicl_classLabel._str_label; 
      this->_idx_label = aicl_classLabel._idx_label; 
      this->_t_coutlabel = aicl_classLabel._t_coutlabel;
    }

    return *this;
    
  }


  //inline const std::string& getNameClass()
  inline const std::string& getLabel()
  {
    return this->_str_label;
  }

  inline const T_IDX_LABEL getIdxLabel()
  {
    return this->_idx_label;
  }

  inline T_COUNT_LABEL getNumLabel()
  {
    return this->_t_coutlabel;
  }
 
  void increasesCountLabel() {
    ++this->_t_coutlabel;
  }

  bool operator< 
    (const CountLabel
     <T_IDX_LABEL,
      T_COUNT_LABEL>    &aicl_classLabel) const
  {
    return _idx_label < aicl_classLabel._idx_label;
  }

  void print(std::ostream &aipf_outFile=std::cout, char aic_separator=',')
  {
    aipf_outFile << "(idx" << aic_separator << this->_idx_label << aic_separator 
		 << "labe" << aic_separator << this->_str_label << aic_separator 
		 << "number of label"<< aic_separator << this->_t_coutlabel 
		 << ")";
  }
protected:
  std::string   _str_label;
  T_IDX_LABEL   _idx_label;
  T_COUNT_LABEL _t_coutlabel;
};

    
template < typename T_IDX_LABEL,
	   typename T_COUNT_LABEL
	   >
std::vector
<CountLabel
 <T_IDX_LABEL,
  T_COUNT_LABEL> >  
mapcountlabel_mapToVector
(std::map<const std::string,CountLabel<T_IDX_LABEL,T_COUNT_LABEL> >  &aimap_countLabel) 
{
  std::vector
    <CountLabel
     <T_IDX_LABEL,
      T_COUNT_LABEL> > lovector_instancesClassLabel;
    
  for (typename std::map<const std::string, 
	 CountLabel<T_IDX_LABEL,T_COUNT_LABEL> >::iterator it  = aimap_countLabel.begin(); 
       it != aimap_countLabel.end(); ++it) {
    lovector_instancesClassLabel.push_back(it->second);
  }

  if ( lovector_instancesClassLabel.size() > 0 ) {
    std::sort
      (lovector_instancesClassLabel.begin(), 
       lovector_instancesClassLabel.end() 
       );
  }
  
  return lovector_instancesClassLabel;
}

  
} /*END namespace data 
   */

#endif  /* _COUNT_LABEL_HPP */
