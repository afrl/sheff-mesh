#include "mesh/common/Selection.h"

namespace mesh
{
void mesh::ElementSelection::run_selection(const Mesh& mesh)
{
    this->selected_indices_.clear();

    if (this->non_local_criteria_.empty())
    {
        for (auto& element : mesh.elements())
        {
            bool meets_criteria = true;
            for (auto& criterion : this->local_criteria_)
            {
                if (!criterion(element))
                {
                    meets_criteria = false;
                    break;
                }
            }
            if (meets_criteria)
            {
                this->selected_indices_.insert(element.local_index());
            }
        }
    }
    else
    {
        // Get the intersection of all non-local criteria
        std::set<LocalIndex> non_local_criteria_intersection;
        bool first_criteria = true;
        for (auto& criterion : this->non_local_criteria_)
        {
            if (first_criteria)
            {
                criterion(mesh, non_local_criteria_intersection);
                first_criteria = false;
            }
            else
            {
                std::set<LocalIndex> criteria_set;
                criterion(mesh, criteria_set);

                // Take the intersection of criteria_set and the current intersection of previous
                // criteria (non_local_criteria_intersection) and put the result in
                // tmp_intersection.
                std::set<LocalIndex> tmp_intersection;
                std::set_intersection(non_local_criteria_intersection.begin(),
                                      non_local_criteria_intersection.end(),
                                      criteria_set.begin(),
                                      criteria_set.end(),
                                      std::inserter(tmp_intersection, tmp_intersection.begin()));

                // We don't need tmp_intersection so just swap with non_local_criteria_intersection
                // rather than copy
                non_local_criteria_intersection.swap(tmp_intersection);
            }
        }

        // Now run the local criteria over the intersection of the non-local criteria to get
        // the final set.
        for (auto& element_index : non_local_criteria_intersection)
        {
            bool meets_criteria = true;
            for (auto& criterion : this->local_criteria_)
            {
                if (!criterion(mesh.elements()[element_index]))
                {
                    meets_criteria = false;
                    break;
                }
            }
            if (meets_criteria)
            {
                this->selected_indices_.insert(element_index);
            }
        }
    }
}

void NodeSelection::run_selection(const Mesh& mesh)
{
    this->selected_indices_.clear();

    if (this->non_local_criteria_.empty())
    {
        for (auto& node : mesh.nodes())
        {
            bool meets_criteria = true;
            for (auto& criterion : this->local_criteria_)
            {
                if (!criterion(node))
                {
                    meets_criteria = false;
                    break;
                }
            }
            if (meets_criteria)
            {
                this->selected_indices_.insert(node.local_index());
            }
        }
    }
    else
    {
        // Get the intersection of all non-local criteria
        std::set<LocalIndex> non_local_criteria_intersection;
        bool first_criteria = true;
        for (auto& criterion : this->non_local_criteria_)
        {
            if (first_criteria)
            {
                criterion(mesh, non_local_criteria_intersection);
                first_criteria = false;
            }
            else
            {
                std::set<LocalIndex> criteria_set;
                criterion(mesh, criteria_set);

                // Take the intersection of criteria_set and the current intersection of
                // previous criteria (non_local_criteria_intersection) and put the result in
                // tmp_intersection.
                std::set<LocalIndex> tmp_intersection;
                std::set_intersection(non_local_criteria_intersection.begin(),
                                      non_local_criteria_intersection.end(),
                                      criteria_set.begin(),
                                      criteria_set.end(),
                                      std::inserter(tmp_intersection, tmp_intersection.begin()));

                // We don't need tmp_intersection so just swap with
                // non_local_criteria_intersection rather than copy
                non_local_criteria_intersection.swap(tmp_intersection);
            }
        }

        // Now run the local criteria over the intersection of the non-local criteria to get
        // the final set.
        for (auto& node_index : non_local_criteria_intersection)
        {
            bool meets_criteria = true;
            for (auto& criterion : this->local_criteria_)
            {
                if (!criterion(mesh.nodes()[node_index]))
                {
                    meets_criteria = false;
                    break;
                }
            }
            if (meets_criteria)
            {
                this->selected_indices_.insert(node_index);
            }
        }
    }
}
}