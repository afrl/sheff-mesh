#include "mesh/element/ElementContainer.h"

#include "mesh/common/Permutation.h"
#include "mesh/element/ElementAccessor.h"
#include "mesh/element/ElementIterator.h"
#include "mesh/node/NodeAccessor.h"
#include "mesh/node/NodeContainer.h"
#include "mesh/node/NodeIterator.h"

namespace mesh
{
ElementContainer::ElementContainer(NodeContainer& nodes)
  : nodes_(nodes)
{
    connectivity_offsets_ = { 0 };
}

bool ElementContainer::operator==(const ElementContainer& other) const
{
    return global_indices_ == other.global_indices_ && types_ == other.types_ &&
           connectivities_ == other.connectivities_ &&
           connectivity_offsets_ == other.connectivity_offsets_ &&
           parent_partitions_ == other.parent_partitions_ &&
           component_ids_ == other.component_ids_ &&
           groups_ == other.groups_ &&
           lcs_ == other.lcs_;
}

ElementIterator ElementContainer::begin()
{
    return ElementIterator(nodes_, *this);
}

ElementConstIterator ElementContainer::begin() const
{
    return ElementConstIterator(nodes_, *this);
}

ElementConstIterator ElementContainer::end() const
{
    return ElementConstIterator(nodes_, *this, static_cast<LocalIndex>(global_indices_.size()));
}

ElementAccessor ElementContainer::operator[](LocalIndex i)
{
    return ElementAccessor(nodes_, *this, i);
}

ElementConstAccessor ElementContainer::operator[](LocalIndex i) const
{
    return ElementConstAccessor(nodes_, *this, i);
}

ElementAccessor ElementContainer::back()
{
    return ElementAccessor(nodes_, *this, size() - 1);
}

ElementConstAccessor ElementContainer::back() const
{
    return ElementConstAccessor(nodes_, *this, size() - 1);
}

LocalIndex ElementContainer::size() const
{
    return static_cast<LocalIndex>(global_indices_.size());
}

ElementGroup& ElementContainer::group(const std::string group_name)
{
    for (auto& group : groups_)
    {
        if (group.name == group_name)
        {
            return group;
        }
    }
    groups_.push_back({ group_name, {} });
    return groups_.back();
}

const ElementGroup& ElementContainer::group(const std::string group_name) const
{
    for (auto& group : groups_)
    {
        if (group.name == group_name)
        {
            return group;
        }
    }
    throw _UTIL_MSG_EXCEPT_C(
      "The group '" + group_name +
      "' was not found and the container is const, so a new one could not be created.");
}

std::vector<ElementGroup>& ElementContainer::groups()
{
    return groups_;
}

const std::vector<ElementGroup>& ElementContainer::groups() const
{
    return groups_;
}

void ElementContainer::add(ElementType type, const std::vector<LocalIndex>& connectivity)
{
    add(type, utility::slice(connectivity));
}

void ElementContainer::add(ElementType type,
                           const std::vector<LocalIndex>& connectivity,
                           GlobalIndex i,
                           PartitionIndex parent_partition,
                           ComponentIndex component_index)
{
    add(type, utility::slice(connectivity), i, parent_partition, component_index);
}

void ElementContainer::add(ElementType type, const utility::Slice<const LocalIndex>& connectivity)
{
    global_indices_.push_back(static_cast<GlobalIndex>(global_indices_.size()));
    types_.push_back(type);
    connectivities_.insert(std::end(connectivities_),
                           std::begin(connectivity),
                           std::end(connectivity));
    connectivity_offsets_.push_back(connectivities_.size());
    parent_partitions_.push_back(0);
    component_ids_.push_back(0);
    lcs_.push_back(QuaternionContainer());
}

void ElementContainer::add(ElementType type,
                           const utility::Slice<const LocalIndex>& connectivity,
                           GlobalIndex i,
                           PartitionIndex parent_partition,
                           ComponentIndex component_index)
{
    global_indices_.push_back(i);
    types_.push_back(type);
    connectivities_.insert(std::end(connectivities_),
                           std::begin(connectivity),
                           std::end(connectivity));
    connectivity_offsets_.push_back(connectivities_.size());
    parent_partitions_.push_back(parent_partition);
    component_ids_.push_back(component_index);
    lcs_.push_back(QuaternionContainer());
}

ElementAccessor ElementContainer::add_copy(const ElementConstAccessor& other_element)
{
    global_indices_.push_back(other_element.global_index());
    types_.push_back(other_element.type());
    connectivities_.insert(std::end(connectivities_),
                           std::begin(other_element.connectivity()),
                           std::end(other_element.connectivity()));
    connectivity_offsets_.push_back(connectivities_.size());
    lcs_.push_back(QuaternionContainer());
    lcs_.back().copy(other_element.local_coord_sys());
    parent_partitions_.push_back(other_element.parent_partition());
    component_ids_.push_back(other_element.component_id());

    return back();
}

void ElementContainer::replace_connectivity(const LocalIndex& element_index,
                                            const utility::Slice<LocalIndex>& connectivity)
{
    auto conn_start         = connectivity_offsets_[element_index];
    auto existing_conn_size = connectivity_offsets_[element_index + 1] - conn_start;
    auto n_to_copy          = std::min(connectivity.size(), existing_conn_size);

    // Copy over the common size
    for (auto i : utility::range(n_to_copy))
    {
        connectivities_[conn_start + i] = connectivity[i];
    }
    if (connectivity.size() > existing_conn_size)
    {
        connectivities_.insert(connectivities_.begin() + conn_start + n_to_copy,
                               connectivity.begin() + n_to_copy,
                               connectivity.end());

        auto diff = connectivity.size() - n_to_copy;
        for (auto i = element_index + 1; i < connectivity_offsets_.size(); i++)
        {
            connectivity_offsets_[i] += diff;
        }
    }
    else if (existing_conn_size > connectivity.size())
    {
        auto diff = existing_conn_size - n_to_copy;
        connectivities_.erase(connectivities_.begin() + conn_start + n_to_copy,
                              connectivities_.begin() + conn_start + n_to_copy + diff);

        for (auto i = element_index + 1; i < connectivity_offsets_.size(); i++)
        {
            connectivity_offsets_[i] -= diff;
        }
    }
}

void ElementContainer::remove_stable(LocalIndex i)
{
    _UTIL_ASSERT_C(groups_.size() == 0, "TODO: implement removal of element from groups.");

    global_indices_.erase(global_indices_.begin() + i);
    types_.erase(types_.begin() + i);
    parent_partitions_.erase(parent_partitions_.begin() + i);
    component_ids_.erase(component_ids_.begin() + i);
    lcs_.erase(lcs_.begin() + i);

    auto offset_start = connectivity_offsets_[i];
    auto offset_end   = connectivity_offsets_[i + 1];
    auto conn_size    = offset_end - offset_start;
    connectivities_.erase(connectivities_.begin() + offset_start,
                          connectivities_.begin() + offset_end);
    connectivity_offsets_.erase(connectivity_offsets_.begin() + i);
    for (auto j = i; j < connectivity_offsets_.size(); j++)
    {
        connectivity_offsets_[j] -= conn_size;
    }
}

void ElementContainer::remove_stable(std::vector<LocalIndex>& v)
{
    for (auto& i : v)
    {
        remove_stable(i);
    }
}

void ElementContainer::swap(ElementContainer& other)
{
    global_indices_.swap(other.global_indices_);
    types_.swap(other.types_);
    connectivities_.swap(other.connectivities_);
    connectivity_offsets_.swap(other.connectivity_offsets_);
    parent_partitions_.swap(other.parent_partitions_);
    component_ids_.swap(other.component_ids_);
    groups_.swap(other.groups_);
    lcs_.swap(other.lcs_);
}
void ElementContainer::permute(const std::vector<LocalIndex>& backward_permutation_indices)
{
    detail::apply_permutation(global_indices_, backward_permutation_indices);
    detail::apply_permutation(types_, backward_permutation_indices);
    detail::apply_permutation(parent_partitions_, backward_permutation_indices);
    detail::apply_permutation(component_ids_, backward_permutation_indices);

    // Connectivities have to be permuted specially
    std::vector<LocalIndex> new_connectivities;
    std::vector<std::size_t> new_connectivity_offsets;
    new_connectivities.reserve(connectivities_.size());
    new_connectivity_offsets.reserve(connectivity_offsets_.size());
    new_connectivity_offsets.push_back(0);
    for (auto new_index : utility::indices(backward_permutation_indices))
    {
        const auto& old_index = backward_permutation_indices[new_index];
        for (auto i = connectivity_offsets_[old_index]; i < connectivity_offsets_[old_index + 1];
             i++)
        {
            new_connectivities.push_back(connectivities_[i]);
        }
        new_connectivity_offsets.push_back(new_connectivities.size());
    }
    connectivities_.swap(new_connectivities);
    connectivity_offsets_.swap(new_connectivity_offsets);

    if (groups_.size() > 0)
    {
        throw _UTIL_MSG_EXCEPT_C(
          "TODO implement permutation of element indices within element groups.");
    }

    if (!lcs_.empty())
    {
        detail::apply_permutation(lcs_, backward_permutation_indices);
    }
}

void ElementContainer::permute_connectivity(
  const std::vector<LocalIndex>& node_forward_permutation_indices)
{
    for (auto& c : connectivities_)
    {
        c = node_forward_permutation_indices[c];
    }
}

void ElementContainer::clear()
{
    global_indices_.clear();
    types_.clear();
    connectivities_.clear();
    connectivity_offsets_ = { 0 };
    parent_partitions_.clear();
    component_ids_.clear();
    groups_.clear();
    lcs_.clear();
}
}