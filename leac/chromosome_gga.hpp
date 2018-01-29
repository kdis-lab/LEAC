/*! \file chromosome_gga.hpp
 *
 * \brief chromosome GGA \cite Agustin:etal:GAclusteringVarK:GGA:2012
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef __CHROMOSOME_GGA_HPP__
#define __CHROMOSOME_GGA_HPP__

#include <algorithm>  
#include <vector>
#include "chromosome_fixedlength.hpp"

#include "verbose_global.hpp"


/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {
  
  
/*! \class ChromosomeGGA
  \brief Chromosome encoding two parts
  \details  The encoding is carried out by separating each individual two parts: c = [l|g], the first part is the element section, whereas the second part is called the group section of the individual \cite Agustin:etal:GAclusteringVarK:GGA:2012
*/ 
template <class T_CLUSTERIDX,
	  class T_METRIC //DATA TYPE OBJETIVE FUNCTION AND T_FITNESS, 
	  >
class ChromosomeGGA:
  public ChromFixedLength<T_CLUSTERIDX,T_METRIC>
{   
public:
  ChromosomeGGA()
    :  ChromFixedLength<T_CLUSTERIDX,T_METRIC>()
    ,  _uintidx_groupSecCapacity((uintidx) 0)
    ,  _pts_groupSecString(NULL)
    ,  _uintidx_groupSecSize((uintidx) 0)
  {
    
  }

  ChromosomeGGA(const T_METRIC airt_objetiveFunc, const T_METRIC airt_fitness)
    :  ChromFixedLength<T_CLUSTERIDX,T_METRIC>(airt_objetiveFunc, airt_fitness)
    ,  _uintidx_groupSecCapacity((uintidx) 0)
    ,  _pts_groupSecString(NULL)
    ,  _uintidx_groupSecSize((uintidx) 0)
    
  { }
  
  ChromosomeGGA(const uintidx aiuintidx_groupSecSize)
  : ChromFixedLength<T_CLUSTERIDX,T_METRIC>()
  , _uintidx_groupSecCapacity((uintidx) aiuintidx_groupSecSize + 5)
  , _pts_groupSecString(new T_CLUSTERIDX[_uintidx_groupSecCapacity])
  , _uintidx_groupSecSize(aiuintidx_groupSecSize)
  { }
  
  //move constructor 
  ChromosomeGGA(ChromosomeGGA<T_CLUSTERIDX,T_METRIC> &&aichrom_b)
    : ChromFixedLength<T_CLUSTERIDX,T_METRIC>(aichrom_b)
    , _uintidx_groupSecCapacity(aichrom_b._uintidx_groupSecCapacity)
    , _pts_groupSecString(aichrom_b._pts_groupSecString)
    , _uintidx_groupSecSize(aichrom_b._uintidx_groupSecSize)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "Move:ChromosomeGGA::ChromosomeGGA";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeGGA<>: &&aichrom_b[" <<  &aichrom_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeGGA<T_CLUSTERIDX,T_METRIC>::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    aichrom_b._pts_groupSecString = NULL;
    aichrom_b._uintidx_groupSecSize = 0;
    
  }

  //copy constructor
  ChromosomeGGA
  (const ChromosomeGGA<T_CLUSTERIDX,T_METRIC> &aichrom_b)
    : ChromFixedLength<T_CLUSTERIDX,T_METRIC>(aichrom_b)
    , _uintidx_groupSecCapacity(aichrom_b._uintidx_groupSecSize + 5) 
    , _pts_groupSecString(new T_CLUSTERIDX[aichrom_b._uintidx_groupSecSize + 5])
    , _uintidx_groupSecSize(aichrom_b._uintidx_groupSecSize)
  {

      interfacesse::copy
	(this->_pts_groupSecString, aichrom_b._pts_groupSecString, this->_uintidx_groupSecSize);

  }
  
  virtual ~ChromosomeGGA() 
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeGGA::ChromosomeGGA~";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ':' << geverbosepc_labelstep
		<< '[' << geverboseui_idproc << ':' << this << ']'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    if ( this->_pts_groupSecString != NULL ) {
      delete[] _pts_groupSecString; 
    }
  }

  inline static void setElementSize(uintidx aiuintidx_elementSize) 
  {
     ChromFixedLength<T_CLUSTERIDX,T_METRIC>::setStringSize(aiuintidx_elementSize);
 
  }
   
  inline static const uintidx getElementSize()
  {
    return ChromFixedLength<T_CLUSTERIDX,T_METRIC>::stcgetStringSize();
 
  }
  
  void initializeGroupSec()
  {
    T_CLUSTERIDX lcidx_k = 0;
    for(uintidx lui_k = 0; lui_k <  _uintidx_groupSecSize; lui_k++) {
      _pts_groupSecString[lui_k] = lcidx_k;
      ++lcidx_k;
    }
  }
  
  ChromosomeGGA<T_CLUSTERIDX,T_METRIC>& 
  operator=(const ChromosomeGGA<T_CLUSTERIDX,T_METRIC> &aichrom_b)
  { 
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeGGA::operator:copy=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeGGA<>: &aichrom_b[" <<  &aichrom_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

    
    if ( this != &aichrom_b ) { 
      ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(aichrom_b);
      if ( _uintidx_groupSecCapacity <= aichrom_b._uintidx_groupSecSize ) {
	if ( this->_pts_groupSecString != NULL )  delete[] _pts_groupSecString;
	_uintidx_groupSecCapacity = aichrom_b._uintidx_groupSecSize + 5;
	
	_pts_groupSecString = new T_CLUSTERIDX[_uintidx_groupSecCapacity];
      }
      
      _uintidx_groupSecSize = aichrom_b._uintidx_groupSecSize;
      
      interfacesse::copy
	(_pts_groupSecString, 
	 aichrom_b._pts_groupSecString, 
	 _uintidx_groupSecSize
	 );
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      
      ChromosomeGGA<T_CLUSTERIDX,T_METRIC>::print();
      assert
	(((this->_pts_groupSecString!=NULL)?this->getNumClusterK():-1) ==
	 ((aichrom_b._pts_groupSecString!=NULL)?aichrom_b.getNumClusterK():-1) );
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    
    return *this;
  }

  ChromosomeGGA<T_CLUSTERIDX,T_METRIC>& 
  operator=(ChromosomeGGA<T_CLUSTERIDX,T_METRIC> &&aichrom_b)
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeGGA::operator:move=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeGGA<>: &aichrom_b[" <<  &aichrom_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichrom_b ) {
      ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(aichrom_b);
      if ( this->_pts_groupSecString != NULL ) {
	delete[] _pts_groupSecString; 
      }
      this->_uintidx_groupSecCapacity = aichrom_b._uintidx_groupSecCapacity;
      this->_pts_groupSecString = aichrom_b._pts_groupSecString;
      
      this->_uintidx_groupSecSize = aichrom_b._uintidx_groupSecSize;
      
      aichrom_b._uintidx_groupSecCapacity = 0;
      aichrom_b._pts_groupSecString = NULL;
      aichrom_b._uintidx_groupSecSize = 0;
      
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeGGA<T_CLUSTERIDX,T_METRIC>::print();
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    return *this;
  }

  inline const T_CLUSTERIDX getNumClusterK()  const
  {
    assert( _uintidx_groupSecSize < _uintidx_groupSecCapacity ||
	    assert_msg( _uintidx_groupSecSize << "<" << _uintidx_groupSecCapacity));
    return _pts_groupSecString[_uintidx_groupSecSize-1] + 1; 
  }
  
  inline const uintidx getGroupSecSize() const 
  {
    return this->_uintidx_groupSecSize;
  }

  void increaseGroupSecSize()
  {
    ++_uintidx_groupSecSize;
    if ( _uintidx_groupSecCapacity <= _uintidx_groupSecSize ) {
      if ( this->_pts_groupSecString != NULL )  delete[] _pts_groupSecString;
      _uintidx_groupSecCapacity += 5;
      _pts_groupSecString = new T_CLUSTERIDX[_uintidx_groupSecCapacity];
	
      this->initializeGroupSec();
    }
    else {
      _pts_groupSecString[_uintidx_groupSecSize-1] =  _uintidx_groupSecSize-1;
    }
  }

  inline void decrementGroupSecSize(uintidx aiuintidx_decrement)
  {
    _uintidx_groupSecSize -= aiuintidx_decrement;  
  }
  
  inline void setGroupSecString(T_CLUSTERIDX *aips_string) 
  {
    interfacesse::copy
      (this->_pts_groupSecString, 
       aips_string, 
       this->_uintidx_groupSecSize
       );
  }

  inline T_CLUSTERIDX* getGroupSecString() 
  {
    return this->_pts_groupSecString;
  }

  inline const T_CLUSTERIDX* getGroupSecString() const
  {
    return this->_pts_groupSecString;
  }

  void deleteGroupNull(T_CLUSTERIDX aicidx_k)
  {
    T_CLUSTERIDX *larraycidx_iChrom = this->getString();
    const T_CLUSTERIDX *larraycidx_iChromEnd = this->getString() + this->_stcui_stringSize;
   
    while (larraycidx_iChrom != larraycidx_iChromEnd) {
      if ( aicidx_k < *larraycidx_iChrom )
	--(*larraycidx_iChrom);
      ++larraycidx_iChrom;
    }
    this->decrementGroupSecSize(1);
      
  }

  bool  isValid() const 
  {
    bool lob_isValid = true;
    std::vector<uintidx> lvectorui_countgene(this->getNumClusterK(),0);
    
    T_CLUSTERIDX lcidx_k = 0;
    uintidx lui_k = 0;
      while ( (lui_k < _uintidx_groupSecSize) &&  (_pts_groupSecString[lui_k] == lcidx_k ) ) {
	++lcidx_k; ++lui_k;
      }
    if (lui_k < _uintidx_groupSecSize)
      lob_isValid = false;

    if ( lob_isValid ) {

      const T_CLUSTERIDX *larraycidx_iChrom = this->getString();
	    
      for (uintidx lui_i = 0; lui_i <  this->_stcui_stringSize; lui_i++) {
	  ++lvectorui_countgene.at(*larraycidx_iChrom);
	  ++larraycidx_iChrom;
	}

      const auto la_numNullGene =
	std::count_if
	(lvectorui_countgene.begin(),
	 lvectorui_countgene.end(),
	 [] (const T_CLUSTERIDX aiit_num) {return aiit_num == 0;}
	 );

      if ( la_numNullGene != 0)
	lob_isValid = false;

    }

    return lob_isValid;

  }
  
  virtual void  print
  (std::ostream &os=std::cout,
   const char* aipc_label   = "",
   const char aic_delimCoef = ',',
   const char aic_delimRow  = ';'
   ) const 
  {

#if defined(__VERBOSE_YES)
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    os << "CHROMAGUSTIN2012[" << this << ']'
       << ",lengelem,"  << this->getElementSize()
       << ",lenggroup," << _uintidx_groupSecSize
       << ",capacity,"  << _uintidx_groupSecCapacity
       << ",codek," <<  (( _pts_groupSecString!=NULL)?this->getNumClusterK():-1)  
       << '>';
#else
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
#endif

    
    if ( this->_pts_string != NULL ) {
      for(uintidx li_j = 0; li_j < ( this->_stcui_stringSize - 1 ); li_j++) {
	os << this->_pts_string[li_j] << aic_delimCoef;
      }
      os << this->_pts_string[(this->_stcui_stringSize-1)];
    }
        
    if ( this->_pts_groupSecString != NULL ) {
      os << '|';
      for(uintidx li_j = 0; li_j < ( this->_uintidx_groupSecSize - 1 ); li_j++) {
	os << this->_pts_groupSecString[li_j] << aic_delimCoef;
      }
      os << this->_pts_groupSecString[(this->_uintidx_groupSecSize-1)];
    }
  }

protected:

  uintidx      _uintidx_groupSecCapacity;
  T_CLUSTERIDX *_pts_groupSecString;
  uintidx      _uintidx_groupSecSize;
  
}; //End ChromosomeGGA

} /*END namespace gaencode*/

#endif  /*__CHROMOSOME_GGA_HPP__*/


