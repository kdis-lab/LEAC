/*! \file interval_positiveintegers.hpp
 *
 * \brief interval positive integers
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef INTERVAL_POSITIVE_INTEGER_HPP
#define INTERVAL_POSITIVE_INTEGER_HPP


/*! \namespace runtime
  \brief Module for obtaining run-time statistics
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace runtime {


/*! \class IntervalPositiveInteger
  \brief Interval of positive integer
*/
template < typename T_INTEGERDATATYPE >
class IntervalPositiveInteger {
public:
  IntervalPositiveInteger(T_INTEGERDATATYPE aiui_sizeInterval) :
    ui_upperBound(0), ui_sizeInterval(aiui_sizeInterval) {}

  ~IntervalPositiveInteger() {}

  inline void initialize() { 
    this->ui_upperBound = 0;
  }

  T_INTEGERDATATYPE getLowerBound()
  {  
    return (this->ui_upperBound < this->ui_sizeInterval)?0
      :this->ui_upperBound - this->ui_sizeInterval;
  }

  T_INTEGERDATATYPE getUpperBound() 
  {
    return this->ui_upperBound;
  }

  T_INTEGERDATATYPE upperBoundToEvaluate() {
    return this->ui_upperBound > this->ui_sizeInterval? this->ui_upperBound
      :this->ui_sizeInterval;
  }

  inline void increaseUpperBound() { 
    ++this->ui_upperBound;
  }

  inline void setSizeInterval(T_INTEGERDATATYPE aiui_sizeInterval)   
  {
    this->ui_sizeInterval = aiui_sizeInterval;
  }

  inline T_INTEGERDATATYPE getSizeInterval() 
  {
    return this->ui_sizeInterval;
  }

private:
  T_INTEGERDATATYPE   ui_upperBound;    /*MAXIMO VALOR EVALUADO*/
  T_INTEGERDATATYPE   ui_sizeInterval;  
};


} /*END namespace executiontime 
   */


#endif  /*INTERVAL_POSITIVE_INTEGER_HPP*/


