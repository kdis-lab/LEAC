/*! \file instance.hpp
 *
 * \brief Module for the handling of instances or also called objects or points
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef __INSTANCE_HPP
#define __INSTANCE_HPP

#include <cmath>
#include <limits>
#include "line_split.hpp"
#include "outfilename.hpp"
#include "common.hpp"


/*! \namespace data
  \brief Module for the handling of instances or also called objects or points.
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace data {

/*! \class Instance
  \brief Patterns, usually vectors in a multidimensional space
*/
template < class T_FEATURE > 
class Instance
{
public:
 
  Instance()
    : _ptc_id(NULL)
    , _arrayt_feature(new T_FEATURE[_stcuintidx_numDimensions])  
  {
    if (_stcb_homogeneousCoord ) {
    _arrayt_feature[_stcuintidx_numDimensions-1] = T_FEATURE(1);
    }
  }

  Instance(T_FEATURE*  aiarrayt_feature)
    : _ptc_id(NULL)
    , _arrayt_feature(new T_FEATURE[_stcuintidx_numDimensions])  
  {
    interfaceclapack_copy
      (_arrayt_feature,
       aiarrayt_feature,
       _stcuintidx_numDimensions
      );
  }

  //copy constructor
  Instance(const Instance<T_FEATURE> &aiinst_b)
    : _ptc_id(NULL)
    , _arrayt_feature(new T_FEATURE[_stcuintidx_numDimensions])
  {
    if (aiinst_b._ptc_id != NULL ) {
      this->_ptc_id = new char[strlen(aiinst_b._ptc_id) + 1];
      strcpy(this->_ptc_id, aiinst_b._ptc_id);
    }
    interfaceclapack_copy
      (_arrayt_feature,
       aiinst_b._arrayt_feature,
       _stcuintidx_numDimensions
      );
  }

  //move constructor
  Instance(Instance<T_FEATURE> &&aiinst_b)
    : _ptc_id(aiinst_b._ptc_id)
    , _arrayt_feature(aiinst_b._arrayt_feature)
  {
    aiinst_b._ptc_id = NULL;
    aiinst_b._arrayt_feature = NULL;
  }

  virtual ~Instance() 
  { 
    if ( this->_ptc_id != NULL ) {
      delete [] this->_ptc_id;
    }
    if ( this->_arrayt_feature != NULL ) {
      delete[] _arrayt_feature;
    }
  }

  Instance<T_FEATURE>& operator=(const Instance<T_FEATURE> &aiinst_b)
  {
    if( this != &aiinst_b ) {

      if ( this->_ptc_id != NULL ) { 
	delete [] this->_ptc_id;
	this->_ptc_id = NULL;
      }

      if (aiinst_b._ptc_id != NULL ) {
	this->_ptc_id = new char[strlen(aiinst_b._ptc_id) + 1];
	strcpy(this->_ptc_id, aiinst_b._ptc_id);
      }
      interfaceclapack_copy
	(_arrayt_feature,
	 aiinst_b._arrayt_feature,
	 _stcuintidx_numDimensions
	 );      
    }

    return *this;

  }

  Instance<T_FEATURE>& operator=(Instance<T_FEATURE> &&aiinst_b)
  {

    if( this != &aiinst_b ) {
      if ( this->_ptc_id != NULL )  
	delete [] this->_ptc_id;
      delete[] _arrayt_feature;
      
      _ptc_id = aiinst_b._ptc_id;
      _arrayt_feature = aiinst_b._arrayt_feature;
 
      aiinst_b._ptc_id = NULL;
      aiinst_b._arrayt_feature = NULL;

    }

    return *this;
  }

  inline static void setNumDimensions(uintidx aiuintidx_numDimensions) 
  {
    _stcuintidx_numDimensions =
      (_stcb_homogeneousCoord)?(aiuintidx_numDimensions+1):aiuintidx_numDimensions;
  }

  inline static uintidx getNumDimensions()
  {
    return  _stcuintidx_numDimensions;
  } 

  inline static void setHomogeneousCoord(bool aib_homogeneousCoord) 
  {
    _stcb_homogeneousCoord = aib_homogeneousCoord;
  }

  inline static bool getHomogeneousCoord()
  {
    return _stcb_homogeneousCoord;
  }

  inline static T_FEATURE type()
  {
    return T_FEATURE(1);
  } 
  
#if  DATATYPE_CENTROIDS_ROUND == 0

  inline static T_FEATURE getInfiniteFeature()
  {
    return std::numeric_limits<T_FEATURE>::infinity();
  }

  inline static bool isInfiniteFeature(T_FEATURE aiT_feacture)
  {
    return std::isinf(aiT_feacture);
  }

#else

  inline static T_FEATURE getInfiniteFeature()
  {
    return std::numeric_limits<T_FEATURE>::max();
  }

  inline static bool isInfiniteFeature(T_FEATURE aiT_feacture)
  {
    return (aiT_feacture == std::numeric_limits<T_FEATURE>::max());
  }

#endif /*DATATYPE_CENTROIDS_ROUND*/

  inline const T_FEATURE* getFeatures() const
  {
    return _arrayt_feature;
  } 

  inline T_FEATURE* getFeatures()
  {
    return _arrayt_feature;
  } 

  inline const T_FEATURE getAttribute(uintidx aui_idxAttribute) const
  {
    return _arrayt_feature[aui_idxAttribute];
  } 

  void readFeature(inout::LineSplit &ails_linesplit)
  {
    using namespace std;
    istringstream liss_stringstream;
    uintidx        luintidx_numDimensions =
      (_stcb_homogeneousCoord)?_stcuintidx_numDimensions-1:_stcuintidx_numDimensions;
    
    for (uintidx luintidx_j = 0; luintidx_j < luintidx_numDimensions; luintidx_j++) {
      liss_stringstream.clear();
    
        liss_stringstream.str(ails_linesplit.getItemSelect(luintidx_j));
      liss_stringstream >> _arrayt_feature[luintidx_j]; 
    
    }
    
  } 
  
  void setId(const std::string &aistr_id)
  {
    if ( this->_ptc_id == NULL ) {
      this->_ptc_id = new char[aistr_id.size() + 1];
    }
    else {
      if ( strlen(this->_ptc_id) != aistr_id.size() ) {
	delete [] this->_ptc_id;
	this->_ptc_id = new char[aistr_id.size() + 1];
      }
    }
    strcpy(this->_ptc_id,aistr_id.c_str());
  }

  inline const char* getId() const
  {
    return this->_ptc_id;
  }

  friend std::ostream& operator<<(std::ostream& os, const Instance<T_FEATURE> &aiinst_instance)
  {
    aiinst_instance.print(os,inout::OutFileName::getDelim());
    
    return os;
  }

  virtual std::string
  getToString(const char aic_delim='\t')
  {
    std::stringstream lss_instance;
    uintidx            luintidx_numDimensions;
    
    luintidx_numDimensions = _stcuintidx_numDimensions-1;
    if ( _ptc_id != NULL ) 
      lss_instance << _ptc_id << aic_delim;
    for( uintidx luintidx_j = 0; luintidx_j < luintidx_numDimensions; luintidx_j++){
      lss_instance << _arrayt_feature[luintidx_j] << aic_delim;
    } 
    lss_instance << _arrayt_feature[luintidx_numDimensions];

    return lss_instance.str();
  }

  virtual void print(std::ostream &os=std::cout, const char aic_delim='\t') const
  {
    uintidx luintidx_numDimensions = _stcuintidx_numDimensions-1;
    if ( _ptc_id != NULL ) 
      os << _ptc_id << aic_delim;
    for( uintidx luintidx_j = 0; luintidx_j < luintidx_numDimensions; luintidx_j++){
      os << _arrayt_feature[luintidx_j] << aic_delim;
    } 
    os << _arrayt_feature[luintidx_numDimensions];
  } 

protected:
  
  char           *_ptc_id;
  T_FEATURE      *_arrayt_feature;  
  static uintidx _stcuintidx_numDimensions;
  static bool    _stcb_homogeneousCoord;

}; /*Instances*/

  
template <class T_FEATURE>
uintidx Instance<T_FEATURE>::_stcuintidx_numDimensions = 0;

template <class T_FEATURE>
bool Instance<T_FEATURE>::_stcb_homogeneousCoord = false;

  
} /*END namespace data 
   */

#endif /*__INSTANCE_HPP*/
