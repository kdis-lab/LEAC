/*! \file transformation.hpp
 *
 * \brief Transformation operations with matrices
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef  TRANSFORMATION_HPP
#define  TRANSFORMATION_HPP

#include <vector>
#include "matrix.hpp"
#include "container_instance.hpp"
#include "stats_instances.hpp"


/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/


namespace mat {

  
template < typename T>  
MatrixRow<T>       
translation
(const std::vector<T>& aivector_translation)
{

#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "mat::translation  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n(input const std::vector<T>& aivector_translation[" 
	      << &aivector_translation << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  MatrixRow<T>  lomatrixrow_translation
    ((uintidx)(aivector_translation.size() + 1),
     (uintidx)(aivector_translation.size() + 1)
     );
  lomatrixrow_translation.initialize();
  
  for (uintidx lui_i = 0; lui_i < aivector_translation.size(); ++lui_i) {
    lomatrixrow_translation(lui_i,lui_i) = T(1);
    lomatrixrow_translation(lui_i,aivector_translation.size()) =
      aivector_translation[lui_i];
  }

  lomatrixrow_translation
    (aivector_translation.size(),
     aivector_translation.size()
     ) = T(1);

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "mat::translation OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\nMatrixRow<T>: lomatrixrow_translation[" << &lomatrixrow_translation << "]\n"
	      <<  lomatrixrow_translation
	      <<  std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lomatrixrow_translation;
  
}


template < typename T>  
MatrixRow<T>       
scale
(const std::vector<T>& aivector_scale)
{
  
#ifdef __VERBOSE_YES 
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "mat::scale  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\n(input const std::vector<T>& aivector_scale[" 
	      << &aivector_scale << ']'
	      << "\n)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  MatrixRow<T>  lomatrixrow_scale
    ((uintidx) (aivector_scale.size() + 1),
     (uintidx) (aivector_scale.size() + 1)
     );
  lomatrixrow_scale.initialize();
  
  for (uintidx lui_i = 0; lui_i < aivector_scale.size(); ++lui_i) {
    lomatrixrow_scale(lui_i,lui_i) = aivector_scale[lui_i];
  }

  lomatrixrow_scale
    (aivector_scale.size(),
     aivector_scale.size()
     ) = T(1);

#ifdef __VERBOSE_YES 
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "mat::scale OUT"
	      << '(' << geiinparam_verbose << ')'
	      << "\nMatrixRow<T>: lomatrixrow_scale[" << &lomatrixrow_scale << "]\n"
	      <<  lomatrixrow_scale
	      << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
  return lomatrixrow_scale;
  
}


template<typename T>
void point
(MatrixRow<T>& aimatrixt_trans,
 T             *aiarrayt_point,
 long int      aili_idxCoordX,
 long int      aili_idxCoordY,
 long int      aili_idxCoordZ,
 FILE          *aifile_outCoord
)
{
  double ld_x = 0.0;  
  double ld_y = 0.0;  
  double ld_z = 0.0;  

  if (-1 < aili_idxCoordX  )  {
    ld_x = (double)
      interfacesse::dot
      (aimatrixt_trans.getRow(aili_idxCoordX),
       aiarrayt_point,
       aimatrixt_trans.getNumColumns()
       );
  }

  if ( -1 < aili_idxCoordY )  {
    ld_y = (double)
      interfacesse::dot
      (aimatrixt_trans.getRow(aili_idxCoordY),
       aiarrayt_point,
       aimatrixt_trans.getNumColumns()
       );
  }

  if ( -1 < aili_idxCoordZ )  {
    ld_z = (double)
      interfacesse::dot
      (aimatrixt_trans.getRow(aili_idxCoordZ),
       aiarrayt_point,
       aimatrixt_trans.getNumColumns()
       );
  }

  fprintf(aifile_outCoord, "%.18e,%.18e,%.18e\n", ld_x, ld_y, ld_z);
}

template<typename T>
void point
(MatrixRow<T>& aimatrixt_trans,
 T             *aiarrayt_point,
 const char    *aiastr_label,
 long int      aili_idxCoordX,
 long int      aili_idxCoordY,
 long int      aili_idxCoordZ,
 FILE          *aifile_outCoord
)
{
  double ld_x = 0.0; 
  double ld_y = 0.0; 
  double ld_z = 0.0; 

  if (-1 < aili_idxCoordX  )  {
    ld_x = (double)
      interfacesse::dot
      (aimatrixt_trans.getRow(aili_idxCoordX),
       aiarrayt_point,
       aimatrixt_trans.getNumColumns()
       );
  }

  if ( -1 < aili_idxCoordY )  {
    ld_y = (double)
      interfacesse::dot
      (aimatrixt_trans.getRow(aili_idxCoordY),
       aiarrayt_point,
       aimatrixt_trans.getNumColumns()
       );
  }

  if ( -1 < aili_idxCoordZ )  {
    ld_z = (double)
      interfacesse::dot
      (aimatrixt_trans.getRow(aili_idxCoordZ),
       aiarrayt_point,
       aimatrixt_trans.getNumColumns()
       );
  }

  fprintf(aifile_outCoord, "%s,%.18e,%.18e,%.18e\n",aiastr_label,ld_x, ld_y, ld_z);
}

/*getMatrixTransformPCA  
  el dominio debe ser solo para reales 
  Por que cuando se escala los enteros no se puede
 */
template < typename INPUT_ITERATOR>
auto 
getMatrixTransformPCA
(INPUT_ITERATOR               aiiterator_instfirst,
 const INPUT_ITERATOR         aiiterator_instlast
 ) -> MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> 
{
  const uintidx lui_numInstances
    (uintidx(std::distance(aiiterator_instfirst,aiiterator_instlast)));
  
#ifdef __VERBOSE_YES
  const char* lpc_labelFunc = "mat::getMatrixTransformPCA";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
              << ":  IN(" << geiinparam_verbose << ")\n"
              << "\n(input INPUT_ITERATOR aiiterator_instfirst["
	      << *aiiterator_instfirst << "]\n"
	      << " input INPUT_ITERATOR aiiterator_instlast[" << *aiiterator_instlast << "]"
	      << "\n\t)"
	      << std::endl;
  }
#endif /*__VERBOSE_YES*/
  
  uintidx lui_numRows =
    (data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getHomogeneousCoord())?
    data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions() -1:
    data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getNumDimensions();

  MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>&&
    lmatrixrowt_instancesTrans =
    data::toMatrixRowTrans
    (aiiterator_instfirst,
     aiiterator_instlast,
     data::Instance<decltype((*aiiterator_instfirst)->getAttribute(0))>::getHomogeneousCoord()
     );

  /*Calculate Mean of Instances for columna
   */
  std::vector<decltype((*aiiterator_instfirst)->getAttribute(0))>&&
    lvector_meanFeactures =
    stats::meanRow<decltype((*aiiterator_instfirst)->getAttribute(0)),
    decltype((*aiiterator_instfirst)->getAttribute(0))>
    (lmatrixrowt_instancesTrans);
  
  /*Calculate StdDesv of Instances for Columna
   */
  std::vector<decltype((*aiiterator_instfirst)->getAttribute(0))> lvector_desvstdFeactures
    (lui_numRows);
  
  for (uintidx luintidx_i = 0; luintidx_i < lui_numRows; ++luintidx_i) {

    interfacesse::transy
      (lmatrixrowt_instancesTrans.getRow(luintidx_i),
       -lvector_meanFeactures[luintidx_i],
       lui_numInstances
       );
     
    decltype((*aiiterator_instfirst)->getAttribute(0)) lt_sumDim2 =
      interfacesse::dot
      (lmatrixrowt_instancesTrans.getRow(luintidx_i),
       lmatrixrowt_instancesTrans.getRow(luintidx_i),
       lui_numInstances
       );

    lvector_desvstdFeactures[luintidx_i] =
      std::sqrt( lt_sumDim2 / (decltype((*aiiterator_instfirst)->getAttribute(0))) (lui_numInstances -1));
    
      }

       for (uintidx luintidx_i = 0; luintidx_i < lui_numRows; ++luintidx_i) {
	 lvector_meanFeactures[luintidx_i] = -lvector_meanFeactures[luintidx_i];
	 lvector_desvstdFeactures[luintidx_i] =
	   (lvector_desvstdFeactures[luintidx_i] ==
	    (decltype((*aiiterator_instfirst)->getAttribute(0))) 0.0)?1.0:
	 decltype((*aiiterator_instfirst)->getAttribute(0))(1) / lvector_desvstdFeactures[luintidx_i];
       }
  
       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>   lmatrixrow_transMeans =
       translation(lvector_meanFeactures);
  
       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>   lmatrixrow_scaleStdDesv =
       scale(lvector_desvstdFeactures);
  
       /*Matrix of Varcovar for Instances
	*/
  
       /* Normalize instances  DesvStd
	*/
       for (uintidx luintidx_i = 0; luintidx_i < lui_numRows; ++luintidx_i) {
      
	 interfacesse::scal
	   (lmatrixrowt_instancesTrans.getRow(luintidx_i),
	    lvector_desvstdFeactures[luintidx_i],
	    lui_numInstances
	    );
       }

       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>       
	 lmatrixrowt_varcovar
	 (lui_numRows,
	  lui_numRows
	  );
  
       /* Y' * Y
	*/
       for(uintidx luintidx_i = 0; luintidx_i < lmatrixrowt_varcovar.getNumRows(); luintidx_i++) {
	 for(uintidx luintidx_j = 0; luintidx_j <= luintidx_i; luintidx_j++) {
	   decltype((*aiiterator_instfirst)->getAttribute(0)) lt_dot = 
	     interfacesse::dot
	     (lmatrixrowt_instancesTrans.getRow(luintidx_i),
	      lmatrixrowt_instancesTrans.getRow(luintidx_j),
	      lmatrixrowt_instancesTrans.getNumColumns()
	      );
	   lmatrixrowt_varcovar(luintidx_i,luintidx_j) =
	     lt_dot / (decltype((*aiiterator_instfirst)->getAttribute(0))) (lmatrixrowt_instancesTrans.getNumColumns()-1);

	 }
       }
       lmatrixrowt_varcovar.cpyLoweToUpper();

       /*EigenVal
	*/
       std::vector<decltype((*aiiterator_instfirst)->getAttribute(0))>  lvectort_eigenVal;
       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>    lmatrixcolumnt_eigenVec;
       std::tie(lvectort_eigenVal,lmatrixcolumnt_eigenVec) =
	 syevd(lmatrixrowt_varcovar);

       
       std::vector<std::pair<uintidx,decltype((*aiiterator_instfirst)->getAttribute(0))> >
	 lvectorpair_idxEigenValSort;
       lvectorpair_idxEigenValSort.reserve(lvectort_eigenVal.size());
       for (uintidx luintidx_i = 0; luintidx_i < lvectort_eigenVal.size(); ++luintidx_i)  {
	 lvectorpair_idxEigenValSort.emplace_back(luintidx_i,lvectort_eigenVal[luintidx_i]);
       }

       std::sort
	 (lvectorpair_idxEigenValSort.begin(), lvectorpair_idxEigenValSort.end(), 
	  [](const std::pair<uintidx,decltype((*aiiterator_instfirst)->getAttribute(0))> &left,
	     const std::pair<uintidx,decltype((*aiiterator_instfirst)->getAttribute(0))> &right) 
	  {
	    return left.second > right.second; 
	  }
	  );
  
       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))> lmatrixrowt_pcaLoadings
	 (lui_numRows + 1,
	  lui_numRows + 1
	  );

       lmatrixrowt_pcaLoadings.initialize();

       lmatrixcolumnt_eigenVec = getTranspose(lmatrixcolumnt_eigenVec);
       
       for (uintidx luintidx_i = 0; luintidx_i < lvectorpair_idxEigenValSort.size(); ++luintidx_i)  {
	 //Copy Row    
	 lmatrixrowt_pcaLoadings.copyRow
	   (luintidx_i,
	    lmatrixcolumnt_eigenVec.getRow(lvectorpair_idxEigenValSort[luintidx_i].first),
	    lui_numRows
	    );
       }

       lmatrixrowt_pcaLoadings
	 (lui_numRows,
	  lui_numRows
	  ) = (decltype((*aiiterator_instfirst)->getAttribute(0))) 1;

       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>  lomatrixrow_tmp
	 (lui_numRows + 1,
	  lui_numRows + 1
	  );
  
       interfacesse::gemm
	 (lomatrixrow_tmp,
	  lmatrixrow_scaleStdDesv,
	  lmatrixrow_transMeans
	  );

       MatrixRow<decltype((*aiiterator_instfirst)->getAttribute(0))>  lomatrixrow_transPCA
	 (lui_numRows + 1,
	  lui_numRows + 1
	  );

       interfacesse::gemm
	 (lomatrixrow_transPCA,
	  lmatrixrowt_pcaLoadings,
	  lomatrixrow_tmp
	  );

#ifdef __VERBOSE_YES
       if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	 std::cout << lpc_labelFunc
		   << ": OUT(" << geiinparam_verbose << ")\n";
	 std::ostringstream lostrstream_labelPCA;
	 lostrstream_labelPCA
	   << "<PCAMATRIX:" << lpc_labelFunc
	   << ":lomatrixrow_transPCA[" << &lomatrixrow_transPCA << ']';
	 lomatrixrow_transPCA.print(std::cout,lostrstream_labelPCA.str().c_str(),',',';');
	 std::cout << std::endl;
       }
       --geiinparam_verbose;
#endif //__VERBOSE_YES

  
       return lomatrixrow_transPCA;

}




} /*END namespace trans*/


#endif /*TRANSFORMATION_HPP*/
