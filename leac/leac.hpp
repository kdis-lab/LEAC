/*! \file leac.hpp
 *
 * \brief  LEAC
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

//Headers Domain of the problem

#include <matrix.hpp>

#include <clustering_operator_centroids.hpp>
#include <clustering_operator_crispmatrix.hpp>
#include <clustering_operator_fuzzy.hpp>
#include <clustering_operator_medoids.hpp>
#include <clustering_operator_hierarchical.hpp>
#include <clustering_operator_label.hpp>

#include <partition_centroids.hpp>
#include <partition_disjsets.hpp>
#include <partition_fuzzy.hpp>
#include <partition_label.hpp>
#include <partition_labelvector.hpp>
#include <partition_medoids.hpp>
#include <partition_crispmatrix.hpp>

#include <matching_matrix.hpp>
#include <supervised_measures.hpp>

#include <unsupervised_measures.hpp>

//Headers GA and EA

#include <probability_distribution.hpp>
#include <probability_selection.hpp>

#include <chromosome_crispmatrix.hpp>
#include <chromosome_fixedlength.hpp>

//Headers Extension of chromosomes

#include <chromosome_fgka.hpp>
#include <chromosome_igka.hpp>

#include <ga_binary_operator.hpp>
#include <ga_integer_operator.hpp>
#include <ga_generic_operator.hpp>
#include <ga_clustering_operator.hpp>
#include <ga_real_operator.hpp>
#include <ga_iterator.hpp>

#include <ga_function_objective.hpp>
#include <ga_function_objective_sse.hpp>

#include <dist_euclidean.hpp>

//Headers Other tools

#include <random_ext.hpp>
#include <verbose.hpp> //Headers To check correct, use the directive -D __VERBOSE_YES

#include <graph_utils.hpp>
