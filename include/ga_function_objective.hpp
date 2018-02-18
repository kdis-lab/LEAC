/*! \file ga_function_objective.hpp
 *
 * \brief definition of the abstract class function objective
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef GA_FUNCTION_OBJECTIVE_HPP
#define GA_FUNCTION_OBJECTIVE_HPP

/*! \namespace gafuncobj
  \brief Definition of objective function
  \details 
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace gafuncobj {

/*! \class GAFunctionObjective
  \brief Definition of the abstract class GAFunctionObjective. It is used to evaluate an objective function, in a virtual way in some other function.
*/  
template < class T_GENE,       //DATATYPE OF CHROMOSOME
           class T_METRIC
	   >
class GAFunctionObjective {
public:
  GAFunctionObjective()  
  {
  }

  ~GAFunctionObjective() {}

  virtual std::pair<T_METRIC,bool> getObjetiveFunc(T_GENE*) = 0;

  virtual T_METRIC getFitness(T_METRIC) = 0;

protected:

};

}  /*END namespace gafuncobj*/
  
#endif /*GA_FUNCTION_OBJECTIVE_HPP*/
