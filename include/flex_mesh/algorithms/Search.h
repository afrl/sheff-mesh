#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include <utility>
#include <vector>

namespace flex_mesh
{
	template<typename ObjectType>
	class SelectionCriteria
	{
	public:
		SelectionCriteria() {}
		virtual ~SelectionCriteria() {}

		virtual const bool eval_criterion(const ObjectType& obj) const = 0;
		virtual std::unique_ptr<SelectionCriteria> clone() const = 0;
	};

	/// The most general criteria that can be used (a lambda function).
	template<typename ObjectType>
	class LambdaCriteria : public SelectionCriteria<ObjectType>
	{
	public:
		LambdaCriteria() = delete;
		LambdaCriteria(std::function<bool(const ObjectType& obj)> criteria)
			: criteria_lambda_(criteria) {}
		virtual ~LambdaCriteria() {}

		virtual const bool eval_criterion(const ObjectType& obj) const override
		{
			return criteria_lambda_(obj);
		}
		virtual std::unique_ptr<SelectionCriteria> clone() const override
		{
			return std::make_unique<LambdaCriteria>(criteria_lambda_);
		}

	protected:
		std::function<bool(const ObjectType& obj)> criteria_lambda_;
	};

	/// A class defining a selection of nodes/elements that can be used to get the objects
	/// fitting the selection from the mesh.  This still needs significant development.
	template<typename MeshType, typename ObjectType>
	struct Selection
	{
	public:
		using ElementTypename = typename MeshType::ElementTypename;
		Selection() {}
		Selection(const Selection& other)
		{
			for (auto& c : other.criteria_)
			{
				criteria_.push_back(c.clone());
			}
		}
		Selection(Selection&& other)
		{
			selected_objects_.swap(other.selected_objects_);
			criteria_.swap(other.criteria_);
		}
		virtual ~Selection() {}

		void add_criteria(std::unique_ptr<SelectionCriteria<ObjectType>> criteria)
		{
			criteria_.push_back(std::move(criteria));
		}
		virtual void run_selection(const MeshType& mesh) = 0;
		const std::set<ObjectType *>& get_selection() const
		{
			return selected_objects_;
		}

	protected:
		std::set<ObjectType *> selected_objects_;
		std::vector<std::unique_ptr<SelectionCriteria<ObjectType>>> criteria_;
	};

	/// A class defining a selection of elements that can be used to get the elements
	/// fitting the selection from the mesh.  This still needs significant development.
	template<typename MeshTypename, typename ElementTypename>
	struct ElementSelection : public Selection<MeshTypename, ElementTypename>
	{
	public:
		virtual ~ElementSelection() {}

		virtual void run_selection(const MeshTypename& mesh) override
		{
			this->selected_objects_.clear();
			for (auto& e_ptr : mesh.elements())
			{
				bool meets_criteria = true;
				for (auto& criterion : criteria_)
				{
					if (!criterion->eval_criterion(*e_ptr))
					{
						meets_criteria = false;
						break;
					}
				}
				if (meets_criteria)
				{
					this->selected_objects_.insert(e_ptr.get());
				}
			}
		}
	};

	template<typename MeshTypename, typename NodeTypename>
	struct NodeSelection : public Selection<MeshTypename, NodeTypename>
	{
	public:
		virtual ~NodeSelection() {}

		virtual void run_selection(const MeshTypename& mesh) override
		{
			this->selected_objects_.clear();
			for (auto& n_ptr : mesh.nodes())
			{
				bool meets_criteria = true;
				for (auto& criterion : criteria_)
				{
					if (!criterion->eval_criterion(*n_ptr))
					{
						meets_criteria = false;
						break;
					}
				}
				if (meets_criteria)
				{
					this->selected_objects_.insert(n_ptr.get());
				}
			}
		}
	};

	/// Given a container of element pointers, this method will find all of the
	/// elements with the given component id.
	/// Returns a pair of vectors.  The first in the pair is the indices relative to the
	/// given elements container, and the second in the pair is the pointers.
	template<typename ElementContainer>
	static std::pair<std::vector<LocalIndex>, std::vector<typename ElementContainer::value_type>> find_elements_in_component(
		const ElementContainer& elements,
		const ComponentIndex component_id);

	/// Given a container of node pointers, this method will find all of the
	/// nodes that appear in an node with the given component id.
	/// Returns a pair of vectors.  The first in the pair is the indices relative to the
	/// given nodes container, and the second in the pair is the pointers.
	template<typename NodePtr, typename ElementContainer>
	static std::pair<std::vector<LocalIndex>, std::vector<NodePtr>> find_nodes_in_component(
		const ElementContainer& elements,
		const ComponentIndex component_id);

	/// Returns the minimum coordinate (minimum x, minimum, y, minimum z) based
	/// on a collectino of node pointers.
	template<typename NodeContainer>
	static std::array<FloatType, 3> min_coord(const NodeContainer& nodes);

	/// Returns the maximum coordinate (maximum x, maximum, y, maximum z) based
	/// on a collectino of node pointers.
	template<typename NodeContainer>
	static std::array<FloatType, 3> max_coord(const NodeContainer& nodes);

	/// Returns the nodes and elements respectively for the given partition in the mesh
	template<typename MeshType>
	static std::pair<std::set<typename MeshType::NodeTypename *>, std::set<typename MeshType::ElementTypename *>> get_partition_nodes_and_elements(const MeshType& mesh, const PartitionIndex current_partition = 0);

	/// Returns the number of partitions (assuming they are contiguous 0-based when assigned)
	template<typename MeshType>
	PartitionIndex get_num_partitions(const MeshType& mesh);

	/// Sort the nodes and elements by the parent partition id (ascending)
	template<typename MeshType>
	void sort_by_partition(const MeshType& mesh);

	/////////////////
	// Definitions //
	/////////////////

	template<typename ElementContainer>
	std::pair<std::vector<LocalIndex>, std::vector<typename ElementContainer::value_type>> find_elements_in_component(const ElementContainer & elements, const ComponentIndex component_id)
	{
		using ElementPtr = typename ElementContainer::value_type;
		LocalIndex count = 0;
		std::vector<LocalIndex> relative_indices;
		std::vector<ElementPtr> elements_in_component;
		for (auto& e_ptr : elements) {
			if (e_ptr->get_component_id() == component_id) {
				relative_indices.push_back(count);
				elements_in_component.push_back(e_ptr);
			}
			count++;
		}
		return { relative_indices, elements_in_component };
	}

	template<typename NodePtr, typename ElementContainer>
	std::pair<std::vector<LocalIndex>, std::vector<NodePtr>> find_nodes_in_component(const ElementContainer & elements, const ComponentIndex component_id)
	{
		std::set<std::pair<LocalIndex, NodePtr>> nodes_in_component_set;
		LocalIndex count = 0;
		for (auto& e_ptr : elements) {
			if (e_ptr->get_component_id() == component_id) {
				for (auto& n_ptr : *e_ptr) {
					nodes_in_component_set.insert({ count, n_ptr });
				}
			}
			count++;
		}
		std::vector<LocalIndex> relative_indices;
		std::vector<NodePtr> nodes_in_component;
		for (auto& pair : nodes_in_component_set) {
			relative_indices.push_back(pair.first);
			nodes_in_component.push_back(pair.second);
		}
		return { relative_indices, nodes_in_component };
	}

	template<typename NodeContainer>
	std::array<FloatType, 3> min_coord(const NodeContainer & nodes)
	{
		Point min = nodes.front()->get_coord();
		for (auto& n_ptr : nodes) {
			for (auto i = 0; i < 3; i++) {
				if ((*n_ptr)[i] < min[i]) {
					min[i] = (*n_ptr)[i];
				}
			}
		}
		return min;
	}
	template<typename NodeContainer>
	std::array<FloatType, 3> max_coord(const NodeContainer & nodes)
	{
		Point max = nodes.front()->get_coord();
		for (auto& n_ptr : nodes) {
			for (auto i = 0; i < 3; i++) {
				if ((*n_ptr)[i] > max[i]) {
					max[i] = (*n_ptr)[i];
				}
			}
		}
		return max;
	}

	template<typename MeshType>
	std::pair<std::set<typename MeshType::NodeTypename*>, std::set<typename MeshType::ElementTypename*>> get_partition_nodes_and_elements(const MeshType & mesh, const PartitionIndex current_partition)
	{
		std::set<typename MeshType::NodeTypename *> nodes_in_part;
		std::set<typename MeshType::ElementTypename *> elements_in_part;
		// Loop over the elements
		for (auto& e_ptr : mesh.elements()) {
			// If one of the nodes is owned by this partition,
			// then all the nodes in the element to the partition need to be in
			// the partition's mesh file, and the element needs to be in the
			// partition's mesh file as well (regardless if the element is owned or not).
			for (auto& n_ptr : *e_ptr) {
				if (n_ptr->get_parent_partition() == current_partition) {
					elements_in_part.insert(e_ptr.get());
					for (auto& n_ptr2 : *e_ptr)
						nodes_in_part.insert(n_ptr2);
					break;
				}
			}
		} // End creating sets

		return{ nodes_in_part, elements_in_part };
	}

	template<typename MeshType>
	PartitionIndex get_num_partitions(const MeshType & mesh)
	{
		// Get the number of partitions first
		std::set<PartitionIndex> partitions;
		for (auto& n_ptr : mesh.nodes())
			partitions.insert(n_ptr->get_parent_partition());
		// Should we spend the time to ensure that the partitions are contiguous (no gaps)
		return static_cast<PartitionIndex>(partitions.size());
	}

	template<typename MeshType>
	void sort_by_partition(const MeshType & mesh)
	{
		using NodeUniquePtr = std::unique_ptr<typename MeshType::NodeTypename>;
		using ElementUniquePtr = std::unique_ptr<typename MeshType::ElementTypename>;

		std::sort(mesh.nodes().begin(), mesh.nodes().end(),
			[](const NodeUniquePtr& n1, const NodeUniquePtr& n2) -> bool {
			return n1->get_parent_partition() < n2->get_parent_partition();
		});

		std::sort(mesh.nodes().begin(), mesh.nodes().end(),
			[](const ElementUniquePtr& e1, const ElementUniquePtr& e2) -> bool {
			return e1->get_parent_partition() < e2->get_parent_partition();
		});
	}
}