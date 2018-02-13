/*! \file chromosome_withrownull.hpp
 *
 * \brief chromosome with row null
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef CHROMOSOME_MATRIXROWNULL_HPP
#define CHROMOSOME_MATRIXROWNULL_HPP

#include "chromosome_base.hpp"
#include "matrix_withrownull.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {

/*! \class ChromosomeMatrixWithRowNull
  \brief Matrix with row null  chromosome 
*/
template < class T_GENE, class T_METRIC >
class ChromosomeMatrixWithRowNull
  : public ChromosomeBase<T_METRIC> 
  , public mat::MatrixWithRowNull<T_GENE> 
{
public:
  ChromosomeMatrixWithRowNull() 
    : ChromosomeBase<T_METRIC>()
    , mat::MatrixWithRowNull<T_GENE>()  
  {
  }
  
  ChromosomeMatrixWithRowNull
  (const uintidx aiuintidx_numRows, 
   const uintidx aiuintidx_numRowsMax, 
   const uintidx aiuintidx_numColumns
   )
    : ChromosomeBase<T_METRIC>()
    , mat::MatrixWithRowNull<T_GENE>
      (aiuintidx_numRows,aiuintidx_numRowsMax,aiuintidx_numColumns)
  {
  }

  //copy constructor
  ChromosomeMatrixWithRowNull
  (const ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &aichromwrn_b)
    : ChromosomeBase<T_METRIC>(aichromwrn_b)
    , mat::MatrixWithRowNull<T_GENE>(aichromwrn_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeMatrixWithRowNull::ChromosomeMatrixWithRowNullCopy";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeMatrixWithRowNull<>: &&aichromwrn_b[" <<  &aichromwrn_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeMatrixWithRowNull<T_GENE,T_METRIC>::print(std::cout,lpc_labelFunc);
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
  }

  //move constructor
  ChromosomeMatrixWithRowNull
  (ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &&aichromwrn_b)
    : ChromosomeBase<T_METRIC>(aichromwrn_b)
    , mat::MatrixWithRowNull<T_GENE>(aichromwrn_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeMatrixWithRowNull::ChromosomeMatrixWithRowNullMove";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeMatrixWithRowNull<>: &&aichromwrn_b[" <<  &aichromwrn_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeMatrixWithRowNull<T_GENE,T_METRIC>::print(std::cout,lpc_labelFunc);
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
  }

  virtual ~ChromosomeMatrixWithRowNull()
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeMatrixWithRowNull::ChromosomeMatrixWithRowNull~";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ':' << geverbosepc_labelstep
		<< '[' << geverboseui_idproc << ':' << this << ']'
		<< std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
  }

  ChromosomeMatrixWithRowNull<T_GENE,T_METRIC>& 
   operator=(const ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &aichromwrn_b)
  {
#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeMatrixWithRowNull::operator=Copy";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeMatrixWithRowNull<>: &&aichromwrn_b[" <<  &aichromwrn_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichromwrn_b  ) {
      ChromosomeBase<T_METRIC>::operator=(aichromwrn_b);
      mat::MatrixWithRowNull<T_GENE>::operator=(aichromwrn_b);
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeMatrixWithRowNull<T_GENE,T_METRIC>::print(std::cout,lpc_labelFunc);
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES

    return *this;
  }

  ChromosomeMatrixWithRowNull<T_GENE,T_METRIC>&
   operator=(ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &&aichromwrn_b)
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "ChromosomeMatrixWithRowNull::operator=Move";
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ":  IN(" << geiinparam_verbose << ")\n"
	        << "(input  ChromosomeMatrixWithRowNull<>: &&aichromwrn_b[" <<  &aichromwrn_b
		<< "]\n)"
		<< std::endl;
    }
#endif //__VERBOSE_YES
    
    if ( this != &aichromwrn_b  ) {
      ChromosomeBase<T_METRIC>::operator=(aichromwrn_b);
      mat::MatrixWithRowNull<T_GENE>::operator=(aichromwrn_b);
    }

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << lpc_labelFunc
		<< ": OUT(" << geiinparam_verbose << ")\n";
      ChromosomeMatrixWithRowNull<T_GENE,T_METRIC>::print(std::cout,lpc_labelFunc);
      std::cout << std::endl;
    }
    --geiinparam_verbose;
#endif //__VERBOSE_YES
    
    return *this;
  }

  
  virtual void  print
  (std::ostream &os=std::cout,
   const char* aipc_label   = "",
   const char aic_delimCoef = ',',
   const char aic_delimRow  = ';'
   ) const 
  {
#if defined(__VERBOSE_YES)
    std::string  lpc_labelChrom(":MATRIXWITHROWNULL");
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    mat::MatrixWithRowNull<T_GENE>::print(os,lpc_labelChrom.c_str(),aic_delimCoef,aic_delimRow);
#else
    ChromosomeBase<T_METRIC>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
    mat::MatrixWithRowNull<T_GENE>::print(os,aipc_label,aic_delimCoef,aic_delimRow);
#endif
    
  }

  friend std::ostream& operator<<
  (std::ostream& os, 
   const ChromosomeMatrixWithRowNull<T_GENE,T_METRIC> &aichromwrn_b
   )
  {
    aichromwrn_b.print(os);
    
    return os;
  }

protected:
 
}; /*ChromosomeMatrixWithRowNull*/

} /*END namespace gaencode*/

#endif  /*CHROMOSOME_MATRIXROWNULL_HPP*/
