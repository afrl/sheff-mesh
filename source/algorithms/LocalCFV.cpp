#include "mesh/algorithms/LocalCFV.h"

namespace mesh
{
namespace detail
{
struct CFVSearchState
{
    struct CFVRegionSearchState
    {
        // Each of the vectors will be same length indicating the volume, most recent element added,
        // and the lower bound probability per a local sampling for the corresponding CFV region.
        std::vector<double> volume_by_local_sample;
        std::vector<LocalIndex> element_index_by_sample;
        std::vector<double> failure_probability_by_local_sample;

        /// After the local search, identify_cfv sets this value, which is the index of the sample
        /// for which the probability lower bound is maximized and that is larger than the minimum
        /// CFV allowed (a physics and statistics based limit).
        std::size_t sample_index_corresponding_to_cfv;
    };

    double curr_cumulative_volume = 0.0;
    std::unordered_set<LocalIndex> local_visited_elements;
    std::unordered_set<LocalIndex> local_adjacent_elements;

    std::vector<CFVRegionSearchState> cfv_region_states;
};

void add_local_element(const LocalIndex& i,
                       const ElementAdjacencyList& adjacency_list,
                       const utility::Slice<const double, std::size_t>& volume_by_element,
                       CFVSearchState& search_state)
{
    search_state.curr_cumulative_volume += volume_by_element[i];


    search_state.cfv_region_states.back().element_index_by_sample.push_back(i);
    search_state.cfv_region_states.back().volume_by_local_sample.push_back(
      search_state.curr_cumulative_volume);

    search_state.local_visited_elements.insert(i);
    search_state.local_adjacent_elements.erase(i);

    for (auto j : adjacency_list.get_adjacent_element_indices(i, ElementAdjacencyType::Face))
    {
        if (search_state.local_visited_elements.count(j) == 0)
        {
            search_state.local_adjacent_elements.insert(j);
        }
    }
}

void identify_possible_elements_for_next_cfv_region(
  const LocalIndex& starting_element,
  const ElementAdjacencyList& adjacency_list,
  const utility::Slice<const double, std::size_t>& failure_index_by_element,
  const utility::Slice<const double, std::size_t>& volume_by_element,
  const double& min_cfv,
  const double local_searh_stopping_criteria_max_volume,
  const double local_searh_stopping_criteria_n_neg_grad_samples,
  CFVSearchState& search_state)
{
    search_state.cfv_region_states.push_back(CFVSearchState::CFVRegionSearchState());

    add_local_element(starting_element, adjacency_list, volume_by_element, search_state);

    bool search_stopping_criteria_met = false;
    int n_neg_grad_samples            = 0;
    while (!search_stopping_criteria_met)
    {
        LocalIndex adj_elem_id_with_max_failure_index = 0;
        double max_value                              = -1.0;
        for (auto& i : search_state.local_adjacent_elements)
        {
            if (failure_index_by_element[i] > max_value)
            {
                adj_elem_id_with_max_failure_index = i;
                max_value                          = failure_index_by_element[i];
            }
        }

        add_local_element(adj_elem_id_with_max_failure_index,
                          adjacency_list,
                          volume_by_element,
                          search_state);

        auto last_probability_grad =
          (*search_state.cfv_region_states.back().failure_probability_by_local_sample.end() - 1) -
          (*search_state.cfv_region_states.back().failure_probability_by_local_sample.end() - 2);
        if (last_probability_grad > 0)
        {
            n_neg_grad_samples = 0;
        }
        else
        {
            n_neg_grad_samples++;
        }


        // Stop when the volume reaches the maximum cutoff or when enough samples are found with the
        // negative gradient
        search_stopping_criteria_met =
          search_state.curr_cumulative_volume > min_cfv &&
          (search_state.curr_cumulative_volume > local_searh_stopping_criteria_max_volume ||
           n_neg_grad_samples > local_searh_stopping_criteria_n_neg_grad_samples);
    }
}

void identify_cfv(const utility::Slice<const double, std::size_t>& failure_index_by_element,
                  const double& min_cfv,
                  const double& V_0,
                  const double& alpha,
                  CFVSearchState::CFVRegionSearchState& local_search_state)
{
    for (auto i : utility::range(local_search_state.volume_by_local_sample.size()))
    {
        const double A = (-local_search_state.volume_by_local_sample[i] / V_0);
        const double B =
          pow(failure_index_by_element[local_search_state.element_index_by_sample[i]], alpha);
        const double probability =
          1.0 -
          exp((-local_search_state.volume_by_local_sample[i] / V_0) *
              pow(failure_index_by_element[local_search_state.element_index_by_sample[i]], alpha));
        local_search_state.failure_probability_by_local_sample.push_back(probability);
    }

    auto search_start       = std::lower_bound(local_search_state.volume_by_local_sample.begin(),
                                         local_search_state.volume_by_local_sample.end(),
                                         min_cfv);
    auto search_start_index = search_start - local_search_state.volume_by_local_sample.begin();

    auto iter = std::max_element(local_search_state.failure_probability_by_local_sample.begin() +
                                   search_start_index,
                                 local_search_state.failure_probability_by_local_sample.end());
    auto cfv_sample_index = iter - local_search_state.failure_probability_by_local_sample.begin();
    local_search_state.sample_index_corresponding_to_cfv = cfv_sample_index;
}
}



LocalCFVResult calculate_local_cfv(
  const Mesh& mesh,
  const ElementAdjacencyList& adjacency_list,
  const utility::Slice<const double, std::size_t>& failure_index_by_element,
  const utility::Slice<const double, std::size_t>& volume_by_element,
  const double& min_cfv,
  const double& V_0,
  const double& alpha,
  const double& stopping_probability,
  const double& local_searh_stopping_criteria_max_volume,
  const double& local_searh_stopping_criteria_n_neg_grad_samples)
{
    LocalCFVResult result;
    result.probability_by_element.resize(mesh.elements().size(), 0.0);
    result.cfv_by_element.resize(mesh.elements().size(), 0.0);
    result.cfv_region_by_element.resize(mesh.elements().size(), -1);

    std::vector<mesh::LocalIndex> element_indices_sorted_by_failure_index(mesh.elements().size());
    std::iota(element_indices_sorted_by_failure_index.begin(),
              element_indices_sorted_by_failure_index.end(),
              0);
    std::sort(
      element_indices_sorted_by_failure_index.begin(),
      element_indices_sorted_by_failure_index.end(),
      [&failure_index_by_element](const mesh::LocalIndex& lhs, const mesh::LocalIndex& rhs) -> bool
      { return failure_index_by_element[lhs] > failure_index_by_element[rhs]; });


    int curr_region_index        = 0;
    auto curr_sorted_start_index = 0;
    double last_probability      = 1.0;
    detail::CFVSearchState search_state;
    // Determine a new CFV region until the probability of failure falls below the stopping
    // probability
    while (last_probability > stopping_probability)
    {
        // Move the next highest stressed element that is not already set to start a local
        // volume calculation
        while (failure_index_by_element
                 [element_indices_sorted_by_failure_index[curr_sorted_start_index]])
        {
            curr_sorted_start_index++;
        }

        mesh::LocalIndex starting_element_index =
          element_indices_sorted_by_failure_index[curr_sorted_start_index];

        detail::identify_possible_elements_for_next_cfv_region(
          starting_element_index,
          adjacency_list,
          failure_index_by_element,
          volume_by_element,
          min_cfv,
          local_searh_stopping_criteria_max_volume,
          local_searh_stopping_criteria_n_neg_grad_samples,
          search_state);

        auto local_search_state = search_state.cfv_region_states.back();
        detail::identify_cfv(failure_index_by_element, min_cfv, V_0, alpha, local_search_state);

        const auto cfv_sample_index = local_search_state.sample_index_corresponding_to_cfv;
        const auto& critical_failure_probability =
          local_search_state.failure_probability_by_local_sample[cfv_sample_index];
        const auto& critical_failure_volume =
          local_search_state.volume_by_local_sample[cfv_sample_index];

        for (auto i : utility::range(cfv_sample_index + 1))
        {
            auto& element_id = local_search_state.element_index_by_sample[i];
            if (result.cfv_by_element[element_id] == -1)
            {
                result.probability_by_element[element_id] = critical_failure_probability;
                result.cfv_by_element[element_id]         = critical_failure_volume;
                result.cfv_region_by_element[element_id]  = curr_region_index;
            }
        }

        // std::cout << "CFV: " << critical_failure_volume
        //          << "  probability: " << critical_failure_probability << std::endl;

        curr_region_index++;
        last_probability = critical_failure_volume;
    }

    return result;
}
}