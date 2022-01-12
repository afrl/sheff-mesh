#include "mesh/algorithms/LocalCFV.h"

//#include "mesh/algorithms/IO.h"

#include <fstream>

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
        // The index of the element added to each sample
        std::vector<LocalIndex> element_index_by_sample;
        // The failure index of the element added to each sample (not the minimum used for the lower
        // bound)
        std::vector<double> new_failure_index_by_local_sample;
        std::vector<double> min_failure_index_by_local_sample;
        std::vector<double> vol_avg_failure_index_by_local_sample;
        // Failure probability of the entire sampled volume failing assuming it is at a stress of
        // the lower bound.
        std::vector<double> failure_probability_by_local_sample;
        // Failure probability based on vol. avg stress instead of lower bound
        std::vector<double> alt_failure_probability_by_local_sample;

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
                       const ComponentIndex& component_id,
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

    for (auto e : adjacency_list.get_adjacent_elements(i, ElementAdjacencyType::Face))
    {
        if (e.component_id() == component_id &&
            search_state.local_visited_elements.count(e.local_index()) == 0)
        {
            search_state.local_adjacent_elements.insert(e.local_index());
        }
    }
}

void identify_possible_elements_for_next_cfv_region(
  const LocalIndex& starting_element,
  const Mesh& mesh,
  const ElementAdjacencyList& adjacency_list,
  const utility::Slice<const double, std::size_t>& failure_index_by_element,
  const utility::Slice<const double, std::size_t>& volume_by_element,
  const double& min_cfv,
  const double local_search_stopping_criteria_max_volume,
  const double local_search_stopping_criteria_failure_index_ratio,
  CFVSearchState& search_state)
{
    search_state.curr_cumulative_volume = 0;
    search_state.local_adjacent_elements.clear();
    // search_state.local_visited_elements.clear();
    search_state.cfv_region_states.push_back(CFVSearchState::CFVRegionSearchState());

    add_local_element(starting_element,
                      mesh.elements()[starting_element].component_id(),
                      adjacency_list,
                      volume_by_element,
                      search_state);
    auto& starting_failure_index = failure_index_by_element[starting_element];

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
                          mesh.elements()[adj_elem_id_with_max_failure_index].component_id(),
                          adjacency_list,
                          volume_by_element,
                          search_state);

        // Stop when the volume reaches the maximum cutoff or the ratio of the new element's failure
        // index to the failure index of the starting element is smaller than the given tolerance.
        const auto failure_index_ratio =
          failure_index_by_element[adj_elem_id_with_max_failure_index] / starting_failure_index;
        search_stopping_criteria_met =
          search_state.curr_cumulative_volume > min_cfv &&
          (search_state.curr_cumulative_volume > local_search_stopping_criteria_max_volume ||
           failure_index_ratio < local_search_stopping_criteria_failure_index_ratio);
    }
}

void identify_cfv(const utility::Slice<const double, std::size_t>& failure_index_by_element,
                  const double& min_cfv,
                  const double& V_0,
                  const double& alpha,
                  CFVSearchState::CFVRegionSearchState& local_search_state)
{
    double curr_min_failure_index        = 1e9;
    double curr_failure_index_integrated = 0.0;
    double curr_volume_integrated        = 0.0;
    for (auto i : utility::range(local_search_state.volume_by_local_sample.size()))
    {
        // Equation: 1 - exp(A * B)
        const double new_fi =
          failure_index_by_element[local_search_state.element_index_by_sample[i]];
        if (new_fi < curr_min_failure_index)
        {
            curr_min_failure_index = new_fi;
        }
        curr_failure_index_integrated += new_fi * local_search_state.volume_by_local_sample[i];
        curr_volume_integrated += local_search_state.volume_by_local_sample[i];

        const double A                   = (-local_search_state.volume_by_local_sample[i] / V_0);
        const double B                   = pow(curr_min_failure_index, alpha);
        const double failure_probability = 1.0 - exp(A * B);

        const double B_alt = pow(curr_failure_index_integrated / curr_volume_integrated, alpha);
        const double failure_probability_alt = 1.0 - exp(A * B_alt);

        local_search_state.new_failure_index_by_local_sample.push_back(
          failure_index_by_element[local_search_state.element_index_by_sample[i]]);
        local_search_state.min_failure_index_by_local_sample.push_back(curr_min_failure_index);
        local_search_state.vol_avg_failure_index_by_local_sample.push_back(
          curr_failure_index_integrated / curr_volume_integrated);
        local_search_state.failure_probability_by_local_sample.push_back(failure_probability);
        local_search_state.alt_failure_probability_by_local_sample.push_back(
          failure_probability_alt);
    }

    // Start the search for the max probability at the minimum CFV volume allowed
    auto search_start       = std::lower_bound(local_search_state.volume_by_local_sample.begin(),
                                         local_search_state.volume_by_local_sample.end(),
                                         min_cfv);
    auto search_start_index = search_start - local_search_state.volume_by_local_sample.begin();

    // End the search at the last failure index found?
    // const auto& last_failure_index = local_search_state.min_failure_index_by_local_sample.back();

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
  const double& local_search_stopping_criteria_max_volume,
  const double& local_search_stopping_criteria_failure_index_ratio)
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


    int curr_region_index         = 0;
    auto curr_sorted_start_index  = 0;
    int n_regions_below_threshold = 0;
    detail::CFVSearchState search_state;
    // Determine a new CFV region until the probability of failure falls below the stopping
    // probability
    while (n_regions_below_threshold < 4)
    {
        // Move the next highest stressed element that is not already set to start a local
        // volume calculation
        while (curr_sorted_start_index < element_indices_sorted_by_failure_index.size() &&
               result.cfv_region_by_element
                   [element_indices_sorted_by_failure_index[curr_sorted_start_index]] != -1)
        {
            curr_sorted_start_index++;
        }

        if (curr_sorted_start_index == element_indices_sorted_by_failure_index.size() ||
            failure_index_by_element
                [element_indices_sorted_by_failure_index[curr_sorted_start_index]] < 1.0)
        {
            break;
        }

        mesh::LocalIndex starting_element_index =
          element_indices_sorted_by_failure_index[curr_sorted_start_index];

        detail::identify_possible_elements_for_next_cfv_region(
          starting_element_index,
          mesh,
          adjacency_list,
          failure_index_by_element,
          volume_by_element,
          min_cfv,
          local_search_stopping_criteria_max_volume,
          local_search_stopping_criteria_failure_index_ratio,
          search_state);

        auto& local_search_state = search_state.cfv_region_states.back();
        detail::identify_cfv(failure_index_by_element, min_cfv, V_0, alpha, local_search_state);

        /*
        // Debug: write sequence of volumes
        if (curr_region_index == 3)
        {
            mesh::Mesh curr_v;
            std::vector<double> sigma_lb;
            int spacing = local_search_state.volume_by_local_sample.size() / 240;
            std::cout << "spacing = " << spacing << std::endl;
            int counter = 0;
            for (auto i : utility::range(local_search_state.volume_by_local_sample.size()))
            {
                auto& new_element = local_search_state.element_index_by_sample[i];
                curr_v.copy(mesh, new_element);

                if (i % spacing == 0)
                {
                    sigma_lb.resize(i + 1);
                    auto& fi = local_search_state.min_failure_index_by_local_sample[i];
                    std::fill(sigma_lb.begin(), sigma_lb.end(), fi * 2.41e9);
                    mesh::write_vtk(curr_v,
                                    "region_" + std::to_string(curr_region_index) + "_v2_" +
                                      std::to_string(counter) + ".vtk",
                                    false,
                                    {},
                                    sigma_lb);
                    counter++;
                }
            }
            std::cout << "Finished debug output" << std::endl;
        }
        */

        const auto cfv_sample_index = local_search_state.sample_index_corresponding_to_cfv;
        const auto& critical_failure_probability =
          local_search_state.failure_probability_by_local_sample[cfv_sample_index];
        const auto& critical_failure_volume =
          local_search_state.volume_by_local_sample[cfv_sample_index];


        if (critical_failure_probability > stopping_probability)
        {
            for (auto i : utility::range(cfv_sample_index + 1))
            {
                auto& element_id = local_search_state.element_index_by_sample[i];
                if (result.cfv_region_by_element[element_id] == -1)
                {
                    result.probability_by_element[element_id] = critical_failure_probability;
                    result.cfv_by_element[element_id]         = critical_failure_volume;
                    result.cfv_region_by_element[element_id]  = curr_region_index;
                }
            }

            curr_region_index++;
        }
        else
        {
            search_state.cfv_region_states.pop_back();
            curr_sorted_start_index++;
        }
    }

    for (auto i : utility::indices(search_state.cfv_region_states))
    {
        auto& state       = search_state.cfv_region_states[i];
        auto component_id = mesh.elements()[state.element_index_by_sample.front()].component_id();
        std::ofstream out("region_" + std::to_string(i) + "_component_" +
                          std::to_string(component_id) + ".csv");
        out
          << "sample_id,new_element_id,new_element_failure_index,sample_volume,min_failure_index,"
             "vol_avg_failure_index,failure_probability,alt_failure_probability,cfv_sample_index\n";
        for (auto j : utility::indices(state.element_index_by_sample))
        {
            out << j << "," << state.element_index_by_sample[j] << ","
                << state.new_failure_index_by_local_sample[j] << ","
                << state.volume_by_local_sample[j] << ","
                << state.min_failure_index_by_local_sample[j] << ","
                << state.vol_avg_failure_index_by_local_sample[j] << ","
                << state.failure_probability_by_local_sample[j] << ","
                << state.alt_failure_probability_by_local_sample[j] << ","
                << state.sample_index_corresponding_to_cfv << "\n";
        }
        out.close();
    }

    return result;
}
}