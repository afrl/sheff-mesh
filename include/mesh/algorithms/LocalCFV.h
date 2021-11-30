#pragma once

#include "mesh/Mesh.h"
#include "mesh/relations/ElementAdjacencyList.h"

#include <unordered_set>

namespace mesh
{
struct LocalCFVResult
{
    /// Probability of failure for each element.  Elements not assigned to
    /// a critical failure volume region will have a value of 0.
    std::vector<double> probability_by_element;
    /// The critical failure volume for each element.  Elements not assigned to
    /// a critical failure volume region will have a value of 0.
    std::vector<double> cfv_by_element;
    /// The index of the critical failure volume region for each element.  Elements not assigned to
    /// a critical failure volume region will have a value of -1.
    std::vector<int> cfv_region_by_element;
};

/// @brief TODO document
LocalCFVResult calculate_local_cfv(
  const Mesh& mesh,
  const ElementAdjacencyList& adjacency_list,
  const utility::Slice<const double, std::size_t>& failure_index_by_element,
  const utility::Slice<const double, std::size_t>& volume_by_element,
  const double& min_cfv,
  const double& V_0,
  const double& alpha,
  const double& stopping_probability,
  const double& local_search_stopping_criteria_max_volume,
  const double& local_search_stopping_criteria_failure_index_ratio);
}