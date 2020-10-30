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
           parent_partitions_ == other.parent_partitions_ && groups_ == other.groups_ &&
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
                           GlobalIndex i)
{
    add(type, utility::slice(connectivity), i);
}

void ElementContainer::add(ElementType type,
                           const std::vector<LocalIndex>& connectivity,
                           GlobalIndex i,
                           PartitionIndex parent_partition)
{
    add(type, utility::slice(connectivity), i, parent_partition);
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
    lcs_.push_back(QuaternionContainer());
}

void ElementContainer::add(ElementType type,
                           const utility::Slice<const LocalIndex>& connectivity,
                           GlobalIndex i)
{
    global_indices_.push_back(i);
    types_.push_back(type);
    connectivities_.insert(std::end(connectivities_),
                           std::begin(connectivity),
                           std::end(connectivity));
    connectivity_offsets_.push_back(connectivities_.size());
    parent_partitions_.push_back(0);
    lcs_.push_back(QuaternionContainer());
}

void ElementContainer::add(ElementType type,
                           const utility::Slice<const LocalIndex>& connectivity,
                           GlobalIndex i,
                           PartitionIndex parent_partition)
{
    global_indices_.push_back(i);
    types_.push_back(type);
    connectivities_.insert(std::end(connectivities_),
                           std::begin(connectivity),
                           std::end(connectivity));
    connectivity_offsets_.push_back(connectivities_.size());
    parent_partitions_.push_back(parent_partition);
    lcs_.push_back(QuaternionContainer());
}

ElementAccessor ElementContainer::add_copy(ElementConstAccessor& other_element)
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

    return back();
}

void ElementContainer::swap(ElementContainer& other)
{
    global_indices_.swap(other.global_indices_);
    types_.swap(other.types_);
    connectivities_.swap(other.connectivities_);
    connectivity_offsets_.swap(other.connectivity_offsets_);
    parent_partitions_.swap(other.parent_partitions_);
    groups_.swap(other.groups_);
    lcs_.swap(other.lcs_);
}
void ElementContainer::permute(const std::vector<LocalIndex>& backward_permutation_indices)
{
    detail::apply_permutation(global_indices_, backward_permutation_indices);
    detail::apply_permutation(types_, backward_permutation_indices);
    detail::apply_permutation(parent_partitions_, backward_permutation_indices);

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
        throw _UTIL_MSG_EXCEPT_C("TODO implement permutation of element indices within element groups.");
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
}