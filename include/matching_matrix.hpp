/*! \file matching_matrix.hpp
 *
 * \brief matching matrix
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#ifndef MATCHING_MATRIX_HPP
#define MATCHING_MATRIX_HPP

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include "matrix.hpp"
#include "partition.hpp"

#include "verbose_global.hpp"

/*! \namespace sm
  \brief  Supervised measures 
  \details For a given algorithm's are measures that evaluate the quality of the cluster for previously partition of instances in clusters
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace sm {

/*! \class PartitionMatrix
  \brief Partition Matrix also known as matching matrix is similar to the matrix of confusion, but is used in the unsupervised learning
*/
template < class T_INSTANCES_CLUSTER_K>
class PartitionMatrix: public mat::MatrixRow<T_INSTANCES_CLUSTER_K> 
{
public:
  PartitionMatrix()
    : mat::MatrixRow<T_INSTANCES_CLUSTER_K>()
    , _vectorstr_labelClusterV(0)
  {}
    
  PartitionMatrix
  (const uintidx                    aist_numClusterK,
   const std::vector<std::string>&  aivectorstr_labelClusterV)
   : mat::MatrixRow<T_INSTANCES_CLUSTER_K>(aivectorstr_labelClusterV.size()+1,aist_numClusterK+1)
   , _vectorstr_labelClusterV(aivectorstr_labelClusterV)
  {}

  PartitionMatrix(PartitionMatrix<T_INSTANCES_CLUSTER_K> &&B)
    : mat::MatrixRow<T_INSTANCES_CLUSTER_K>(B)
    , _vectorstr_labelClusterV(B._vectorstr_labelClusterV)
  {}

  PartitionMatrix(const PartitionMatrix<T_INSTANCES_CLUSTER_K> &B)
    : mat::MatrixRow<T_INSTANCES_CLUSTER_K>(B)
    , _vectorstr_labelClusterV(B._vectorstr_labelClusterV)
  {}

  virtual ~PartitionMatrix() { } 

  PartitionMatrix<T_INSTANCES_CLUSTER_K>&
  operator=(const PartitionMatrix<T_INSTANCES_CLUSTER_K> &B)
  {
    if ( this != &B) {
      mat::MatrixRow<T_INSTANCES_CLUSTER_K>::operator=(B);
      _vectorstr_labelClusterV = B._vectorstr_labelClusterV;
    }
    return *this;
  }
  
  PartitionMatrix<T_INSTANCES_CLUSTER_K>
  operator=(PartitionMatrix<T_INSTANCES_CLUSTER_K> &&B)
  {
    if ( this != &B) {
      mat::MatrixRow<T_INSTANCES_CLUSTER_K>::operator=(B);
      _vectorstr_labelClusterV = B._vectorstr_labelClusterV;
    }
    return *this;
  }
  
  void sumPartitionTable()
  {
    T_INSTANCES_CLUSTER_K lul_sumTmp; 
    uintidx list_numClassU = this->getNumRows()-1;
    uintidx list_numClusterV = this->getNumColumns()-1;

    /*SUM ROWS*/
    for (uintidx li_i = 0; li_i < list_numClassU; li_i++) {
      lul_sumTmp = T_INSTANCES_CLUSTER_K(0);      
      for (uintidx li_j = 0; li_j < list_numClusterV; li_j++) { 
	lul_sumTmp += (*this)(li_i,li_j);		
      }
      (*this)(li_i, list_numClusterV) =  lul_sumTmp;
    }
  
    /*SUM COLUMS*/
    for (uintidx li_j = 0; li_j <= list_numClusterV; li_j++) {
      lul_sumTmp = T_INSTANCES_CLUSTER_K(0); 
      for (uintidx li_i = 0; li_i < list_numClassU; li_i++){  
	lul_sumTmp += (*this)(li_i,li_j);
      }
      (*this)(list_numClassU, li_j) =  lul_sumTmp;
    }
  }

  inline 
  const T_INSTANCES_CLUSTER_K 
  getNumObjetos() const
  {
    return (*this)
      (this->getNumRows()-1,
       this->getNumColumns()-1
       );
  }

  const virtual T_INSTANCES_CLUSTER_K getSomeClassUSomeClusterV() const
  {
    return T_INSTANCES_CLUSTER_K(-1);
  }

  const virtual T_INSTANCES_CLUSTER_K getDiffClassUDiffClusterV() const
  {
    return T_INSTANCES_CLUSTER_K(-1);
  }

  
 
  /*getAgreementObjects:
    Used to calculate the rand index defined by
    \cite humber:arabie:clusteranalysis:randindex:1988
  */ 
  const virtual T_INSTANCES_CLUSTER_K getAgreementObjects() const
  {
    return T_INSTANCES_CLUSTER_K(-1);
  }


  const virtual T_INSTANCES_CLUSTER_K getDiffClassUSomeClusterV() const
  {
    return T_INSTANCES_CLUSTER_K(-1);
  }

  

  const T_INSTANCES_CLUSTER_K getMisclassified() const
  {

#ifdef __VERBOSE_YES
    const char* lpc_labelFunc = "PartitionMatrix::getMisclassified:  IN";
  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << lpc_labelFunc 
	      << ":  IN(" << geiinparam_verbose << ')'
	      << "\n( input PartitionMatrix: this[" << this << ']'
	      << "\n)"
	      << std::endl;
  }
#endif //__VERBOSE_YES
  

  uintidx luintidx_numClassU = this->getNumRows()-1;
  uintidx luintidx_numClusterV = this->getNumColumns()-1;

  T_INSTANCES_CLUSTER_K loit_misclassified = 
    (*this)(luintidx_numClassU,luintidx_numClusterV);
  
  if ( luintidx_numClusterV != 0 ) {

     
    std::vector<uintidx> lvectort_idxMax(luintidx_numClusterV,0);
    std::vector<T_INSTANCES_CLUSTER_K> lovectort_sumInstCluster(luintidx_numClassU,0);
    
    for (uintidx li_i = 1; li_i < luintidx_numClassU; li_i++) {
      for (uintidx li_j = 0; li_j < luintidx_numClusterV; li_j++) {
        if  ( (*this)(li_i,li_j) > (*this)(lvectort_idxMax[li_j],li_j) ) {
	  lvectort_idxMax[li_j] = li_i;
	}
      }
    }

    for (uintidx li_j = 0; li_j < luintidx_numClusterV; li_j++) {

      lovectort_sumInstCluster[lvectort_idxMax[li_j]] += (*this)(lvectort_idxMax[li_j],li_j);
    
    }

    loit_misclassified = T_INSTANCES_CLUSTER_K(0);
    
    for (uintidx li_i = 0; li_i < luintidx_numClassU; li_i++) {
      loit_misclassified +=
	(*this)(li_i,luintidx_numClusterV) - lovectort_sumInstCluster.at(li_i); 
    }    
  }
  
#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
          std::cout << lpc_labelFunc
		    << ": OUT(" << geiinparam_verbose << ')' << " loit_misclassified = " << loit_misclassified
		    << std::endl;
  }
  --geiinparam_verbose;
#endif //__VERBOSE_YES
  
    return loit_misclassified;
    
  }
  
  void setLabelClusterV(std::vector<std::string> aivectorstr_labelClusterV)
  {
    _vectorstr_labelClusterV = aivectorstr_labelClusterV;
  }

  const std::vector<std::string>& getlabelClass() const
  {
    return _vectorstr_labelClusterV;
  }

  /*getInstClusterK:
      Return number instances by cluster
    
      precondition: before run method sumPartitionTable()
      postcontition:
   */ 
  std::vector<T_INSTANCES_CLUSTER_K> getInstClusterK()
  {
#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionMatrix::getInstClusterK:  IN"
		<< '(' << geiinparam_verbose << ')'
		<< "\n\t( input PartitionMatrix: this[" << this << "]\n"
		<< "\t)"
		<< std::endl;
    }
#endif /*__VERBOSE_YES*/

    T_INSTANCES_CLUSTER_K *lpt_rowSum = this->getRow( this->getNumRows() -1);
    std::vector<T_INSTANCES_CLUSTER_K>  lovectort_numInstClusterK
      (lpt_rowSum, lpt_rowSum + this->getNumColumns()-1);

#ifdef __VERBOSE_YES
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout << "PartitionMatrix::getInstClusterK: OUT"
		<< '(' << geiinparam_verbose << ')'
		<< "\n\tstd::vector<>  lovectort_numInstClusterK[" << this << "]\n"
		<< lovectort_numInstClusterK
		<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
       
    return lovectort_numInstClusterK;
  }

  virtual void  print
  (std::ostream &os=std::cout,
   const char *aipc_label   = "",
   const char aic_delimCoef = ',',
   const char aic_delimRow  = ';'
   ) const
  {
   
    std::ostringstream lostrstream_labelMatchingMatrix;
    lostrstream_labelMatchingMatrix <<  aipc_label;
    
    lostrstream_labelMatchingMatrix
      << aic_delimCoef << "_cluster" << aic_delimRow;
    for (uintidx li_j = 0; li_j < this->getNumColumns()-1; li_j++) {
      lostrstream_labelMatchingMatrix  << "C_" << (li_j+1) << aic_delimRow;
    }
    lostrstream_labelMatchingMatrix << "sums";

    lostrstream_labelMatchingMatrix
      << aic_delimCoef << "_class" << aic_delimRow;
    for (uintidx li_j = 0; li_j < this->getNumRows()-1; li_j++) {
      lostrstream_labelMatchingMatrix
	<< _vectorstr_labelClusterV.at(li_j)
	<< aic_delimRow;
    }
    lostrstream_labelMatchingMatrix << "sums";
    
    mat::MatrixRow<T_INSTANCES_CLUSTER_K>::print
      (os,
       lostrstream_labelMatchingMatrix.str().c_str(),
       aic_delimCoef,
       aic_delimRow
       );
  }

 
  void r8mat_print
  (const char   *aipc_label   = "",
   const int    aii_numColumns = 5,
   std::ostream &os=std::cout
  ) const
//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in row-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, T A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
  {
    //# define INCX 5

 
    const int m = (int) this->getNumRows();
    const int n = (int) this->getNumColumns();
    
    //T_FEATURE a[];
    int ilo = 1;
    int jlo = 1;
    int ihi = (int) this->getNumRows();
    int jhi = (int) this->getNumColumns();
  
  
    int i;
    int i2hi;
    int i2lo;
    int j;
    int j2hi;
    int j2lo;

    os << "\n";
    os << aipc_label << "\n";

    size_t lst_maxLable = 7;
    for (auto &lstr_i:_vectorstr_labelClusterV ) {
      if ( lstr_i.length() >  lst_maxLable ) {
	lst_maxLable = lstr_i.length();
      }
    }
     
    if ( m <= 0 || n <= 0 )
      {
	os << "\n";
	os << "  (None)\n";
	return;
      }
    //
    //  Print the columns of the matrix, in strips of 5.
    //
    for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + aii_numColumns )
      {
	j2hi = j2lo + aii_numColumns - 1;
	if ( n < j2hi )
	  {
	    j2hi = n;
	  }
	if ( jhi < j2hi )
	  {
	    j2hi = jhi;
	  }
	os << "\n";
	//
	//  For each column J in the current range...
	//
	//  Write the header.
	//
	os << std::right << std::setw(lst_maxLable+2) << "Cluster: ";
	for ( j = j2lo; j <= j2hi; j++ )
	  {
	    if ( j == n )
	      os << "sum";
	    else
	      os << std::left <<std::setw(7) << j - 1 << "       ";
	  }
	/*if ( j == n )
	  os << "sum";
	else
	  os << j;
	*/
	//os << "\n";
	os << "\nClass\n";
	os << "\n";
	//
	//  Determine the range of the rows in this strip.
	//
	if ( 1 < ilo )
	  {
	    i2lo = ilo;
	  }
	else
	  {
	    i2lo = 1;
	  }
	if ( ihi < m )
	  {
	    i2hi = ihi;
	  }
	else
	  {
	    i2hi = m;
	  }

	
	for ( i = i2lo; i <= i2hi; i++ )
	  {
	    //
	    //  Print out (up to) 5 entries in row I, that lie in the current strip.
	    // _vectorstr_labelClusterV.at(li_j)
	    if ( i < i2hi )
	      os << std::right << std::setw(lst_maxLable) << _vectorstr_labelClusterV.at(i - 1) << ": ";
	    else
	      os << std::right << std::setw(lst_maxLable+2) <<  "sum: ";
	    for ( j = j2lo; j <= j2hi; j++ )
	      {
		os <<  std::left  << std::setw(12) << this->_arrayt_data[j-1+(i-1)*n] << "  ";
		//os << std::setw(12) << this->_arrayt_data[i-1+(j-1)*m] << "  ";
	      }
	    os << "\n";
	  }
      }

    return;
    //# undef INCX
    
  }

protected:

  std::vector<std::string>  _vectorstr_labelClusterV;
  
}; /*PartitionMatrixRow*/


/*getSensitivity  
  matchingmatrix_getPercSuccessesClass
 */
template < class T_METRIC,
	   class T_INSTANCES_CLUSTER_K
	   >
std::vector<T_METRIC>
getSensitivity
(const PartitionMatrix
 <T_INSTANCES_CLUSTER_K>& aipartmatrix_a,
 const T_METRIC           airt_percentage = 100.0
 )
{  
  uintidx luintidx_numClassU = aipartmatrix_a.getNumRows()-1;
  uintidx luintidx_numClusterV = aipartmatrix_a.getNumColumns()-1;
  std::vector<T_METRIC> lovectort_perSuccessClass(luintidx_numClassU,T_METRIC(0.0));

  if ( luintidx_numClusterV != 0 ) {
     
    std::vector<uintidx> lvectort_idxMax(luintidx_numClusterV,0);
    std::vector<T_INSTANCES_CLUSTER_K> lovectort_sumInstCluster(luintidx_numClassU,0);
    std::vector<T_INSTANCES_CLUSTER_K> lovectort_totalInstCluster(luintidx_numClassU,0);

    for (uintidx li_i = 1; li_i < luintidx_numClassU; li_i++) {
      for (uintidx li_j = 0; li_j < luintidx_numClusterV; li_j++) {
        if  ( aipartmatrix_a(li_i,li_j) > aipartmatrix_a(lvectort_idxMax[li_j],li_j) ) {
	  lvectort_idxMax[li_j] = li_i;
	}
      }
    }

    for (uintidx li_j = 0; li_j < luintidx_numClusterV; li_j++) {
      lovectort_sumInstCluster[lvectort_idxMax[li_j]] += aipartmatrix_a(lvectort_idxMax[li_j],li_j);
      lovectort_totalInstCluster[lvectort_idxMax[li_j]] += aipartmatrix_a(luintidx_numClassU,li_j);
    }

    for (uintidx li_i = 0; li_i < luintidx_numClassU; li_i++) {
      if ( lovectort_totalInstCluster.at(li_i) != 0 ) 
      lovectort_perSuccessClass.at( li_i) =
	airt_percentage * lovectort_sumInstCluster.at(li_i) /lovectort_totalInstCluster.at(li_i);
    }
    
  }
  
  return lovectort_perSuccessClass;

}


template < class T_METRIC,
	   class T_INSTANCES_CLUSTER_K
	   >
std::string
getStrSenSpe
(const PartitionMatrix
 <T_INSTANCES_CLUSTER_K>& aipartmatrix_a,
 const bool               aib_sensitivity = true,
 const T_METRIC           airt_percentage = 100.0
 )
{
  std::ostringstream lostrstream_strSensitivity;
  uintidx luintidx_numClassU   = aipartmatrix_a.getNumRows()-1;
  uintidx luintidx_numClusterV = aipartmatrix_a.getNumColumns()-1;

  if ( luintidx_numClusterV > 0 ) {
     
    std::vector<T_METRIC> lvectort_perSensitivity(luintidx_numClassU,T_METRIC(0.0));
  
    std::vector<uintidx> lvectort_idxMax(luintidx_numClusterV,0);
    std::vector<T_INSTANCES_CLUSTER_K> lovectort_sumInstCluster(luintidx_numClassU,0);
    std::vector<T_INSTANCES_CLUSTER_K> lovectort_totalInstCluster(luintidx_numClassU,0);

    for (uintidx li_i = 1; li_i < luintidx_numClassU; li_i++) {
      for (uintidx li_j = 0; li_j < luintidx_numClusterV; li_j++) {
        if  ( aipartmatrix_a(li_i,li_j) > aipartmatrix_a(lvectort_idxMax[li_j],li_j) ) {
	  lvectort_idxMax[li_j] = li_i;
	}
      }
    }

    for (uintidx li_j = 0; li_j < luintidx_numClusterV; li_j++) {
      lovectort_sumInstCluster[lvectort_idxMax[li_j]] += aipartmatrix_a(lvectort_idxMax[li_j],li_j);
      lovectort_totalInstCluster[lvectort_idxMax[li_j]] += aipartmatrix_a(luintidx_numClassU,li_j);
    }

    for (uintidx li_i = 0; li_i < luintidx_numClassU; li_i++) {
      if ( lovectort_totalInstCluster.at(li_i) != 0 ) { 
	lvectort_perSensitivity.at( li_i ) = (aib_sensitivity)?
	  airt_percentage * lovectort_sumInstCluster.at(li_i) /
	  lovectort_totalInstCluster.at(li_i):
	  airt_percentage * lovectort_sumInstCluster.at(li_i) /
	  aipartmatrix_a(li_i,luintidx_numClusterV);
      }
    }

    const std::vector<std::string>&  lvectorstr_labelClass = 
      aipartmatrix_a.getlabelClass();

    lostrstream_strSensitivity.precision(COMMON_COUT_PRECISION);
    
    lostrstream_strSensitivity
      <<  lvectorstr_labelClass[0]
      << ';'
      << lvectort_perSensitivity[0];
    
    for ( size_t _lst_ci = 1; _lst_ci < lvectorstr_labelClass.size(); _lst_ci++) {
      lostrstream_strSensitivity
	<< ';' 
	<< lvectorstr_labelClass[_lst_ci] 
	<< ';' 
	<< lvectort_perSensitivity[_lst_ci];
    }
   
  }
  
  return lostrstream_strSensitivity.str();

}



/*! \class ConfusionMatchingMatrix
  \brief  Confusion Matching Matrix also known as an error matrix, is a specific table layout that allows visualization of the performance of an algorithm, typically a supervised learning
*/
template < class T_INSTANCES_CLUSTER_K>
class ConfusionMatchingMatrix: public 
PartitionMatrix<T_INSTANCES_CLUSTER_K> 
{
public:

  ConfusionMatchingMatrix()
    : PartitionMatrix<T_INSTANCES_CLUSTER_K>()
  {}
  
  ConfusionMatchingMatrix
  (const uintidx                   aist_numClusterK,
   const std::vector<std::string>& aivectorstr_labelClusterV)
    :PartitionMatrix<T_INSTANCES_CLUSTER_K>(aist_numClusterK,aivectorstr_labelClusterV)
  {}
  
  virtual ~ConfusionMatchingMatrix() {}

  //(i) = a
  const virtual T_INSTANCES_CLUSTER_K getSomeClassUSomeClusterV() const
  {
    uintidx                lui_numClassU = this->getNumRows()-1;
    uintidx                lui_numClusterV = this->getNumColumns()-1;
    
     T_INSTANCES_CLUSTER_K lit_sumnij2 = 0;
     for (uintidx lui_i = 0; lui_i < lui_numClassU; lui_i++) {
       for (uintidx lui_j = 0; lui_j < lui_numClusterV; lui_j++) {
	 T_INSTANCES_CLUSTER_K lit_nij = (*this)(lui_i,lui_j); 
	 lit_sumnij2 += lit_nij * ( lit_nij -1);
      }
    }  
     return lit_sumnij2 / 2; 
  }


  //(ii) = d
  const virtual T_INSTANCES_CLUSTER_K getDiffClassUDiffClusterV() const
  {
 
    T_INSTANCES_CLUSTER_K lit_sumniDot = _getSumClassniDot2();
 
    T_INSTANCES_CLUSTER_K lit_sumnDotj = _getSumClusternDotj2();
      
    T_INSTANCES_CLUSTER_K  lit_sumnij2 = _getSumnij2();
    
    T_INSTANCES_CLUSTER_K lit_numObj = this->getNumObjetos();
    
    T_INSTANCES_CLUSTER_K loit_diffClassUDiffClusterV = 
      (lit_numObj * lit_numObj +  lit_sumnij2 - (lit_sumniDot + lit_sumnDotj)) / 2;
    
    return loit_diffClassUDiffClusterV;
    
  }

  /*getAgreementObjects:
    Used to calculate the rand index defined by
    \cite{humber:arabie:clusteranalysis:randindex:1988}

    (i) + (ii) = a + d = A
  */ 
  const virtual T_INSTANCES_CLUSTER_K getAgreementObjects() const
  {
   
    T_INSTANCES_CLUSTER_K lit_sumniDot = _getSumClassniDot2();
   
    T_INSTANCES_CLUSTER_K lit_sumnDotj = _getSumClusternDotj2();
   
    T_INSTANCES_CLUSTER_K  lit_sumnij2 = _getSumnij2();

    T_INSTANCES_CLUSTER_K lit_numObj = this->getNumObjetos();
    
    T_INSTANCES_CLUSTER_K loit_agreementObjects = 
      lit_numObj *( lit_numObj-1) / 2 +  lit_sumnij2 - (lit_sumniDot + lit_sumnDotj) / 2;
    
    
    return loit_agreementObjects;

  }

  /*c (iii) objects in the pair are placed in different classes in U(class) and in
    the same class in V(cluster);
    U is class
    V cluster
  */
  const virtual T_INSTANCES_CLUSTER_K getDiffClassUSomeClusterV() const
  {
    /*  T_INSTANCES_CLUSTER_K lit_sumnDotj = _getSumClusternDotj2();
    T_INSTANCES_CLUSTER_K lit_sumnij2  = _getSumnij2();
    */
    return (_getSumClusternDotj2()-_getSumnij2()) /2;
    
  }

  /*b (iv) objects in the pair are placed in the same class in U and in
    different classes in V.
   */
  const virtual T_INSTANCES_CLUSTER_K getSomeClassUDiffClusterV() const
  {
    return ( _getSumClassniDot2() - _getSumnij2() ) /2;
  }

  
protected:


  const T_INSTANCES_CLUSTER_K _getSumClassniDot2() const
  {
    uintidx  lui_numClassU = this->getNumRows()-1;
    uintidx  lui_numClusterV = this->getNumColumns()-1;    

    T_INSTANCES_CLUSTER_K lit_sumniDot = 0;
    for (uintidx lui_i = 0; lui_i < lui_numClassU; lui_i++) {
      T_INSTANCES_CLUSTER_K lit_niDot = (*this)(lui_i,lui_numClusterV);
      lit_sumniDot += lit_niDot * lit_niDot;
    }

    return lit_sumniDot;
    
  }


  const T_INSTANCES_CLUSTER_K _getSumClusternDotj2() const
  {
    uintidx  lui_numClassU   = this->getNumRows()-1;
    uintidx  lui_numClusterV = this->getNumColumns()-1;    

    T_INSTANCES_CLUSTER_K lit_sumnDotj = 0;
    for (uintidx lui_j = 0; lui_j < lui_numClusterV; lui_j++) {
      T_INSTANCES_CLUSTER_K lit_nDotj = (*this)(lui_numClassU,lui_j);
      lit_sumnDotj += lit_nDotj * lit_nDotj;
    }

    return lit_sumnDotj;
  }  
  
  const T_INSTANCES_CLUSTER_K _getSumnij2() const
  {
    uintidx  lui_numClassU   = this->getNumRows()-1;
    uintidx   lui_numClusterV = this->getNumColumns()-1;

    T_INSTANCES_CLUSTER_K lit_sumnij2 = 0;
    for (uintidx lui_i = 0; lui_i < lui_numClassU; lui_i++) {
      for (uintidx lui_j = 0; lui_j < lui_numClusterV; lui_j++) {
	T_INSTANCES_CLUSTER_K lit_nij = (*this)(lui_i,lui_j); 
	lit_sumnij2 += lit_nij * lit_nij;
      }
    }
    return lit_sumnij2;
  }
  
}; /*End class ConfusionMatchingMatrix*/


/*! \fn sm::ConfusionMatchingMatrix getConfusionMatrix(INPUT_ITERATOR aiiterator_instfirst, const INPUT_ITERATOR aiiterator_instlast, partition::Partition<T_CLUSTERIDX>  &aipartition_clusters, const FUNCINSTFREQUENCY func_instfrequency, const FUNCINSTCLASS func_instclass) 
    \brief Get confusion matching matrix
    \details 
    \param aiiterator_instfirst a InputIterator to the initial positions of the sequence of instances
    \param aiiterator_instlast a InputIterator to the final positions of the sequence of instances
    \param aipartition_clusters a partition::Partition of instances in clusters
    \param func_instfrequency a function that returns the frequency of the instance
    \param func_instclass a function that returns the class of the instance
 */
template < typename INPUT_ITERATOR,
	   typename T_CLUSTERIDX,
	   typename FUNCINSTFREQUENCY,
	   typename FUNCINSTCLASS
	   >
auto 
getConfusionMatrix
(INPUT_ITERATOR                             aiiterator_instfirst,
 const INPUT_ITERATOR                       aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>  &aipartition_clusters, 
 const FUNCINSTFREQUENCY                    func_instfrequency,
 const FUNCINSTCLASS                        func_instclass
 ) -> ConfusionMatchingMatrix<decltype(func_instfrequency(*aiiterator_instfirst))>
{
  typedef decltype(func_instfrequency(*aiiterator_instfirst)) ResultType;

  const T_CLUSTERIDX   laicidx_numClusterK =  aipartition_clusters.getNumCluster();
  ConfusionMatchingMatrix<ResultType>
    lomatchmatrix_confusion
    ((uintidx) laicidx_numClusterK,
     data::InstanceIterfazClass
     <ResultType,
     T_CLUSTERIDX>
     ::getLabelClass((uintidx) laicidx_numClusterK )
     );

  lomatchmatrix_confusion.initialize();

  //uintidx li_i = 0;
  for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
 
    T_CLUSTERIDX lcidx_xinK =
      aipartition_clusters.next();//aipartition_clusters.getClusterIdx(li_i);
    if ( 0 <= lcidx_xinK  && lcidx_xinK < laicidx_numClusterK ) {
	  
      lomatchmatrix_confusion
	((uintidx) func_instclass(*aiiterator_instfirst), 
	 (uintidx) lcidx_xinK
	 ) += func_instfrequency(*aiiterator_instfirst);  

    }
    //++li_i; 
  }

  lomatchmatrix_confusion.sumPartitionTable();
  
  return lomatchmatrix_confusion;
} /*matchingmatrix_getConfusionMatrix
   */


template < typename  INPUT_ITERATOR,
	   typename  T_CLUSTERIDX,
	   typename  FUNCINSTFREQUENCY
	   >
auto 
getMatchingMatrix
(INPUT_ITERATOR                             aiiterator_instfirst,
 const INPUT_ITERATOR                       aiiterator_instlast,
 partition::Partition<T_CLUSTERIDX>  &aipartition_clusters, 
 const FUNCINSTFREQUENCY                    func_instfrequency
 ) -> ConfusionMatchingMatrix<decltype(func_instfrequency(*aiiterator_instfirst))>
{
  typedef decltype(func_instfrequency(*aiiterator_instfirst)) ResultType;
  
  const T_CLUSTERIDX   laicidx_numClusterK =  aipartition_clusters.getNumCluster();

  ConfusionMatchingMatrix<ResultType>
    lomatchmatrix_confusion
    ((uintidx) laicidx_numClusterK,
     data::InstanceIterfazClass
     <ResultType,
     T_CLUSTERIDX>
     ::getLabelClass((uintidx) laicidx_numClusterK)
     );
 
  lomatchmatrix_confusion.initialize();

  for (aipartition_clusters.begin(); aiiterator_instfirst != aiiterator_instlast; ++aiiterator_instfirst) {
 
    T_CLUSTERIDX lcidx_xinK =
      aipartition_clusters.next();
    if ( 0 <= lcidx_xinK  && lcidx_xinK < laicidx_numClusterK ) {
	  
       lomatchmatrix_confusion
	((uintidx) lcidx_xinK, 
	 (uintidx) lcidx_xinK
	 ) += func_instfrequency(*aiiterator_instfirst); 

    }
    
  }

  lomatchmatrix_confusion.sumPartitionTable();
  
  return lomatchmatrix_confusion;
  
} /*getConfusionMatrix
   */

} /*END namespace sm clustering measure*/


#endif /*MATCHING_MATRIX_HPP*/
