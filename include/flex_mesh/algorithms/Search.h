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
	/// Given a container of element pointers, this method will find all of the
	/// elements with the given component id.
	template<typename ElementContainer>
	static std::vector<typename ElementContainer::value_type> find_elements_in_component(
		const ElementContainer& elements,
		const ComponentIndex component_id);

	/// Given a container of element pointers, this method will find all of the
	/// nodes that appear in an element with the given component id.
	template<typename NodePtr, typename ElementContainer>
	static std::vector<NodePtr> find_nodes_in_component(
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

	/////////////////
	// Definitions //
	/////////////////

	template<typename ElementContainer>
	std::vector<typename ElementContainer::value_type> find_elements_in_component(const ElementContainer & elements, const ComponentIndex component_id)
	{
		using ElementPtr = ElementContainer::value_type;
		std::vector<ElementPtr> elements_in_component;
		for (auto& e_ptr : elements) {
			if (e_ptr->get_component_id() == component_id) {
				elements_in_component.push_back(e_ptr);
			}
		}
		return elements_in_component;
	}

	template<typename NodePtr, typename ElementContainer>
	std::vector<NodePtr> find_nodes_in_component(const ElementContainer & elements, const ComponentIndex component_id)
	{
		std::set<NodePtr> nodes_in_component_set;
		for (auto& e_ptr : elements) {
			if (e_ptr->get_component_id() == component_id) {
				for (auto& n_ptr : *e_ptr) {
					nodes_in_component_set.insert(n_ptr);
				}
			}
		}
		std::vector<NodePtr> nodes_in_component(nodes_in_component_set.begin(), nodes_in_component_set.end());
		return nodes_in_component;
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
}