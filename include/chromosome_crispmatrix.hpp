/*! \file chromosome_crispmatrix.hpp
 *
 * \brief chromosome crispmatrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CHROMOSOME_BIT_CRISP_MATRIX_HPP
#define CHROMOSOME_BIT_CRISP_MATRIX_HPP

#include "chromosome_base.hpp"
#include "crisp_matrix.hpp"
#include "verbose_global.hpp"


/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {

/*! \class ChromosomeCrispMatrix
  \brief Chromosome bit crisp matrix
  \details 
*/ 
template < class T_BITSIZE,
	   class T_CLUSTERIDX,
	   class T_METRIC
	   >
class ChromosomeCrispMatrix
  : public ChromosomeBase<T_METRIC> 
  , public mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX> 
{
public:
  ChromosomeCrispMatrix()
    : ChromosomeBase<T_METRIC>()
    , mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>()
  { }
  
  ChromosomeCrispMatrix
  (const uintidx aiuintidx_numRows, 
   const uintidx aiuintidx_numColumns
   )
    : ChromosomeBase<T_METRIC>()
    , mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>(aiuintidx_numRows,aiuintidx_numColumns)
  { }

    //copy constructor
  ChromosomeCrispMatrix
  (const ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC>& aichrombitcrispmatrix_b)
    : ChromosomeBase<T_METRIC>(aichrombitcrispmatrix_b)
    , mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>(aichrombitcrispmatrix_b)
  { }

  //move constructor
  ChromosomeCrispMatrix
  (ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC> &&aichrombitcrispmatrix_b)
    : ChromosomeBase<T_METRIC>(aichrombitcrispmatrix_b)
    , mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>(aichrombitcrispmatrix_b)
  { }

  virtual ~ChromosomeCrispMatrix()
  {
  }
  
  //move copy
  ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC>& operator=
  (const ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC> &aichrombitcrispmatrix_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeCrispMatrix::operator=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeCrispMatrix<>: &["
		<<  &aichrombitcrispmatrix_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if( this !=  &aichrombitcrispmatrix_b ){
      ChromosomeBase<T_METRIC>::operator=(aichrombitcrispmatrix_b);
      mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>::operator=(aichrombitcrispmatrix_b);
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelFunc;
      lostrstream_labelFunc <<  geverbosepc_labelstep
			    << ':' << lpc_labelFunc;
      ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC>::print
	(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    return *this;
  }
 
  //move asigned
  ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC>& operator=
  (ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC> &&aichrombitcrispmatrix_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeCrispMatrix::operatormove=";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeCrispMatrix<>: &&[" <<  &aichrombitcrispmatrix_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if( this !=  &aichrombitcrispmatrix_b ){
      ChromosomeBase<T_METRIC>::operator=(aichrombitcrispmatrix_b);
      mat::CrispMatrix<T_BITSIZE,T_CLUSTERIDX>::operator=(aichrombitcrispmatrix_b);
    }


#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      std::ostringstream lostrstream_labelFunc;
      lostrstream_labelFunc <<  geverbosepc_labelstep
			    << ':' << lpc_labelFunc;
      ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC>::print
	(std::cout,lostrstream_labelFunc.str().c_str(),',',';');
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    return *this;
    
  }
 
  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label   = "",
   const char aic_delimCoef =',',
   const char aic_delimRow  ='\n'
   ) const
  {
    
#if defined(__VERBOSE_YES)
    std::string  lpc_labelChrom("CHROMOSOMEBITMATRIX");
    ChromosomeBase<T_METRIC>::print(os,aipc_label);
    os << ',';
    mat::BitMatrix<T_BITSIZE>::print(os,lpc_labelChrom.c_str());
#else
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    os << ',';
    mat::BitMatrix<T_BITSIZE>::print(os);
#endif
    
  }

  friend std::ostream& operator<<
  (std::ostream& os, 
   const ChromosomeCrispMatrix<T_BITSIZE,T_CLUSTERIDX,T_METRIC> &aichrombitcrispmatrix_b)
  {
    aichrombitcrispmatrix_b.print(os,""); 
    
    return os;
  }

 
protected:

  
}; /*ChromosomeCrispMatrix*/

} /*END namespace gaencode*/

#endif  /*CHROMOSOME_BIT_CRISP_MATRIX_HPP*/
