#include "mesh/node/NodeContainer.h"

#include "mesh/common/Permutation.h"
#include "mesh/node/NodeAccessor.h"
#include "mesh/node/NodeIterator.h"

namespace mesh
{
NodeContainer::NodeContainer() {}

bool NodeContainer::operator==(const NodeContainer& other) const
{
    return global_indices_ == other.global_indices_ && points_ == other.points_ &&
           parent_partitions_ == other.parent_partitions_ && groups_ == other.groups_;
}

NodeIterator NodeContainer::begin()
{
    return NodeIterator(*this);
}

NodeConstIterator NodeContainer::begin() const
{
    return NodeConstIterator(*this);
}

NodeConstIterator NodeContainer::end() const
{
    return NodeConstIterator(*this, static_cast<LocalIndex>(global_indices_.size()));
}

NodeAccessor NodeContainer::operator[](LocalIndex i)
{
    return NodeAccessor(*this, i);
}

NodeConstAccessor NodeContainer::operator[](LocalIndex i) const
{
    return NodeConstAccessor(*this, i);
}

LocalIndex NodeContainer::size() const
{
    return static_cast<LocalIndex>(global_indices_.size());
}

NodeGroup& NodeContainer::group(const std::string group_name)
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

const NodeGroup& NodeContainer::group(const std::string group_name) const
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

std::vector<NodeGroup>& NodeContainer::groups()
{
    return groups_;
}

const std::vector<NodeGroup>& NodeContainer::groups() const
{
    return groups_;
}

void NodeContainer::add(const Point& p)
{
    auto next_index = static_cast<GlobalIndex>(global_indices_.size());
    global_indices_.push_back(next_index);
    points_.push_back(p);
    parent_partitions_.push_back(0);
}

void NodeContainer::add(const Point& p, GlobalIndex i)
{
    global_indices_.push_back(i);
    points_.push_back(p);
    parent_partitions_.push_back(0);
}

void NodeContainer::add(const Point& p, GlobalIndex i, PartitionIndex parent_partition)
{
    global_indices_.push_back(i);
    points_.push_back(p);
    parent_partitions_.push_back(parent_partition);
}

void NodeContainer::remove(LocalIndex i)
{
    global_indices_[i] = global_indices_.back();
    global_indices_.pop_back();
    points_[i] = points_.back();
    points_.pop_back();
    parent_partitions_[i] = parent_partitions_.back();
    parent_partitions_.pop_back();
}

void NodeContainer::remove_stable(LocalIndex i)
{
    global_indices_.erase(global_indices_.begin() + i);
    points_.erase(points_.begin() + i);
    parent_partitions_.erase(parent_partitions_.begin() + i);
}

void NodeContainer::swap(NodeContainer& other)
{
    global_indices_.swap(other.global_indices_);
    points_.swap(other.points_);
    parent_partitions_.swap(other.parent_partitions_);
    groups_.swap(other.groups_);
}

void NodeContainer::permute(const std::vector<LocalIndex>& backward_permutation_indices)
{
    detail::apply_permutation(global_indices_, backward_permutation_indices);
    detail::apply_permutation(points_, backward_permutation_indices);
    detail::apply_permutation(parent_partitions_, backward_permutation_indices);
    if (groups_.size() > 0)
    {
        throw _UTIL_MSG_EXCEPT_C("TODO implement permutation of node indices within node groups.");
    }
}

void NodeContainer::clear()
{
    global_indices_.clear();
    points_.clear();
    parent_partitions_.clear();
    groups_.clear();
}
}