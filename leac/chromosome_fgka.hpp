/*! \file chromosome_fgka.hpp
 *
 * \brief chromosome FGKA \cite Lu:etal:GAclusteringLabel:FGKA:2004
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef CHROMOSOME_FGKA_HPP
#define CHROMOSOME_FGKA_HPP

#include "chromosome_fixedlength.hpp"

/*! \namespace gaencode
  \brief Encode chromosome
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace gaencode {

/*! \class ChromosomeFGKA
  \brief Chromosome fixed length string and legality ration
*/    
template <class T_CLUSTERIDX,
	  class T_METRIC
	  >
class ChromosomeFGKA: 
  public gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>
{
public:
  ChromosomeFGKA()
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>()
    , _cidx_numClusterNotNull(0)
    , _b_selected(false)
  {
  }

  //copy constructor
  ChromosomeFGKA
  (const ChromosomeFGKA<T_CLUSTERIDX,T_METRIC> &aich_chromosome)
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>(aich_chromosome)
    , _cidx_numClusterNotNull(aich_chromosome._cidx_numClusterNotNull)
  {
  }

  //move constructor
  ChromosomeFGKA
  (ChromosomeFGKA<T_CLUSTERIDX,T_METRIC> &&aich_chromosome)
    : gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>(aich_chromosome)
    , _cidx_numClusterNotNull(aich_chromosome._cidx_numClusterNotNull)
  {
  }
    
  ~ChromosomeFGKA()
  {
  }

  ChromosomeFGKA<T_CLUSTERIDX,T_METRIC>&
  operator=(const ChromosomeFGKA<T_CLUSTERIDX,T_METRIC>& aichrom_b)
  {
    if ( this != &aichrom_b ) {
      gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(aichrom_b);
      this->_cidx_numClusterNotNull = aichrom_b._cidx_numClusterNotNull;
    }
    
    return *this;
  }

  ChromosomeFGKA<T_CLUSTERIDX,T_METRIC>&
  operator=(ChromosomeFGKA<T_CLUSTERIDX,T_METRIC>&& aichrom_b)
  {
    if ( this != &aichrom_b ) {
      gaencode::ChromFixedLength<T_CLUSTERIDX,T_METRIC>::operator=(aichrom_b);
      this->_cidx_numClusterNotNull = aichrom_b._cidx_numClusterNotNull;
    }
    
    return *this;
  }

  inline static void setNumClusterK(const T_CLUSTERIDX aicidx_numClusterK) 
  {
    _stcidx_numClusterK = aicidx_numClusterK;
  }

  inline static uintidx getNumClusterK()
  {
    return  _stcidx_numClusterK;
  }
  
  inline T_CLUSTERIDX  getNumClusterNotNull()
  {
    return this->_cidx_numClusterNotNull;
  }

  inline void setNumClusterNotNull(T_CLUSTERIDX aicidx_numClusterNotNull)
  {
    this->_cidx_numClusterNotNull =
      aicidx_numClusterNotNull;
  }

  inline T_METRIC getLegalityRation() const
  {
    return  (T_METRIC) this->_cidx_numClusterNotNull 
      / (T_METRIC) this->_stcidx_numClusterK;
  }

  inline bool getSelected() 
  {
    return this->_b_selected;
  }

  inline void setSelected(bool aib_selected) 
  {
    this->_b_selected = aib_selected;
  }
  
protected:
  T_CLUSTERIDX         _cidx_numClusterNotNull;
  bool                        _b_selected; /*used*/
  static T_CLUSTERIDX  _stcidx_numClusterK;
}; /*END  ChromosomeFGKA */

template <class T_CLUSTERIDX,
	  class T_METRIC
	  >
T_CLUSTERIDX ChromosomeFGKA<T_CLUSTERIDX,T_METRIC>::_stcidx_numClusterK = 0;

} /*END namespace gaencode*/
  
#endif  /* CHROMOSOME_FGKA_HPP */
