#pragma once

#include "utility/Exceptions.h"

#include "Element.h"
#include "PolyhedronHelpers.h"

namespace flex_mesh
{

	template<typename NodeType>
	class PolyhedronElement : public Element<NodeType>
	{
	public:

		using Parent = Element<NodeType>;
		using ConnSize = typename Parent::ConnSize;

		//////////////////
		// Constructors //
		//////////////////

		PolyhedronElement() {};
		PolyhedronElement(const PolyhedronElement&) = default;
		virtual ~PolyhedronElement() {};

		/// Initializes the polyhedron element provided the connecitivty of the faces.
		/// Note: nodes should be repeated in the connectivity of each face as needed.
		void initialize(const std::vector<Connectivity<NodeType, ConnSize>>& polyhedron_faces);
		/// Throws an error to prevent use
		virtual void initialize(const ElementType type,
			const std::vector<NodeType *> & connectivity,
			const LCSAssignment lcs_assignment = LCSAssignment::None) override;

		///////////////////
		// Interpolation //
		///////////////////

		/// Resizes and fills the Eigen matrix with the shape function values
		virtual void get_shape_function_values(
			const Point& cartesian_coordinate,
			ElementVector& shape_function_values) const override;

		/// Resizes and fills the Eigen matrix with the shape function derivatives
		/// Note: derivatives are in terms of parammetric coordinates (dPhi/dxi)
		virtual void get_shape_function_derivatives(
			const Point& cartesian_coordinate,
			ElementMatrix3D& shape_function_derivs) const override;

	protected:
		/// TODO document
		utility::DynamicArray<utility::DynamicArray<ConnSize, ConnSize>, ConnSize> polyhedron_faces_;
		utility::DynamicArray<Vertex<const NodeType>, ConnSize> graph_;
	};

	/////////////////
	// Definitions //
	/////////////////


	template<typename NodeType>
	inline void PolyhedronElement<NodeType>::initialize(const std::vector<Connectivity<NodeType, ConnSize>>& polyhedron_faces)
	{
		type_ = ElementType::POLYHEDRON;

		std::set<NodeType*> unique_nodes;
		for (auto& face : polyhedron_faces)
			for (auto& n_ptr : face)
				unique_nodes.insert(n_ptr);

		Parent::connectivity_ = unique_nodes;
		std::map<NodeType*, ConnSize> node_to_conn_index;
		for (ConnSize i = 0; i < Parent::connectivity_.size(); i++)
			node_to_conn_index[Parent::connectivity_[i]] = i;

		// Create face structure
		for (auto i = 0; i < polyhedron_faces.size(); i++) {
			if (polyhedron_faces[i].size() != 3)
				throw _UTIL_MSG_EXCEPT("Currently, only polyhedrons with triangular faces are accepted. TODO implement triangularization of faces.");

			polyhedron_faces_.push_back(DynamicArray<ConnSize, ConnSize>(polyhedron_faces[i].size()));
			for (auto j = 0; j < polyhedron_faces[i].size(); j++)
				polyhedron_faces_[i][j] = node_to_conn_index[polyhedron_faces[i][j]];
		}

		// Build the graph
		// Note: we have half edges for the faces
		// In other words, when we are looping over the edges for the faces, the connecitivty will
		// will cause each edge to be encountered twice, but each time the connectivity of the edge
		// will be reversed.  So each node will have the neighboring nodes updated properly.
		graph_.resize(Parent::connectivity_.size());
		for (auto& face : polyhedron_faces_) {
			for (auto i = 0; i < face.size() - 1; i++)
				graph_[face[i]].neighbors.push_back(Parent::connectivity_[face[i + 1]]);
			graph_[face.back()].neighbors.push_back(Parent::connectivity_[face.front()]);
		}

		// DEBUG Print the graph
		/*
		for (auto node_index = 0; node_index < graph_.size(); node_index++) {
			std::cout << "- " << node_index << std::endl;
			for (auto& n : graph_[node_index].neighbors)
				std::cout << "|- " << n->get_global_number() << std::endl;
		}
		*/
	}

	template<typename NodeType>
	inline void PolyhedronElement<NodeType>::initialize(const ElementType type, const std::vector<NodeType*>& connectivity, const LCSAssignment lcs_assignment)
	{
		throw _UTIL_MSG_EXCEPT("Polyhedral elements require the other initialize method.");
	}

	template<typename NodeType>
	inline void PolyhedronElement<NodeType>::get_shape_function_values(const Point & p, ElementVector & shape_function_values) const
	{
		// Following http://graphics.stanford.edu/~wicke/publications/FEM/wbg-afmocp-07.pdf
		// Step 1: Compute weight ratios (Eqn. 7)
		//  - 1.1: Compute constants c for each neighbor (2 constants per a node per a neighbor) (Eqn 8)
		//         Note: the constant is transitive: c_(a,b) = c_(b,a)
		//         Therefore a matrix can be formed
		//  - 1.2: Compute volume of tets formed by node with its neighbors (4 constants per a node per a neighbor)
		//         Note: the volume of tets are fully transitive V_(a,b,c,d) = V_(a,b,d,c) = V(b,a,c,d) = ...
		//         We can use a hash of the combination as an identifier to only calculate it once
		// Step 2: Normalize weight ratios to get shape functions (Eqn 9)
		// Limitations:
		// - Only for convex polyhedrons
		// - Bad shape functions near faces (good near center)

		// Compute weight ratios
		shape_function_values.resize(Parent::connectivity_.size());
		const Vector3D x(p[0], p[1], p[2]);
		ShapeFunctionHelper<NodeType> cv_helper;
		double sum_weights = 0;
		for (auto node_index = 0; node_index < Parent::connectivity_.size(); node_index++) {
			const auto& n_i = *Parent::connectivity_[node_index];
			const Vertex<const NodeType>& node_as_vertex = graph_[node_index];

			shape_function_values[node_index] = 0;
			for (ConnSize j = 0; j < node_as_vertex.neighbors.size(); j++) {
				const auto& n_j = *node_as_vertex.get_neighbor(j);
				const auto& n_jp1 = *node_as_vertex.get_neighbor(j + 1);
				const auto& n_jm1 = *node_as_vertex.get_neighbor(j - 1);

				const double c_j_jp1 = cv_helper.compute_c(n_j, n_jp1, x);
				const double c_i_j = cv_helper.compute_c(n_i, n_j, x);
				const double V_jm1_jp1_j = cv_helper.compute_v(n_jm1, n_jp1, n_j, x);
				const double V_i_j_jp1 = cv_helper.compute_v(n_i, n_j, n_jp1, x);
				const double V_i_jm1_j = cv_helper.compute_v(n_i, n_jm1, n_j, x);

				shape_function_values[node_index] += c_j_jp1 / V_i_j_jp1 + (c_i_j * V_jm1_jp1_j) / (V_i_jm1_j * V_i_j_jp1);
			}
			sum_weights += shape_function_values[node_index];
		}
		// Normalize weights to ensure partition of unity (PUM)
		shape_function_values /= sum_weights;
	}

	template<typename NodeType>
	inline void PolyhedronElement<NodeType>::get_shape_function_derivatives(const Point & cartesian_coordinate, ElementMatrix3D & shape_function_derivs) const
	{
		shape_function_derivs.resize(Parent::connectivity_.size(), 3);
		/*
		// Compute weight ratios and gradient of weight ratios
		std::vector<double> weights();
		const Vector3D x(cartesian_coordinate[0], cartesian_coordinate[1], cartesian_coordinate[2]);
		ShapeFunctionHelper<NodeType> cv_helper;
		for (auto node_index = 0; node_index < Parent::connectivity_.size(); node_index++) {
			const auto& n_i = *Parent::connectivity_[node_index];
			const Vertex<const NodeType>& node_as_vertex = graph_[node_index];

			weights[node_index] = 0;
			for (ConnSize j = 0; j < node_as_vertex.neighbors.size(); j++) {
				const auto& n_j = *node_as_vertex.get_neighbor(j);
				const auto& n_jp1 = *node_as_vertex.get_neighbor(j + 1);
				const auto& n_jm1 = *node_as_vertex.get_neighbor(j - 1);

				const double c_j_jp1 = cv_helper.compute_c(n_j, n_jp1, x);
				const double c_i_j = cv_helper.compute_c(n_i, n_j, x);
				const double V_jm1_jp1_j = cv_helper.compute_v(n_jm1, n_jp1, n_j, x);
				const double V_i_j_jp1 = cv_helper.compute_v(n_i, n_j, n_jp1, x);
				const double V_i_jm1_j = cv_helper.compute_v(n_i, n_jm1, n_j, x);

				weights[node_index] += c_j_jp1 / V_i_j_jp1 + (c_i_j * V_jm1_jp1_j) / (V_i_jm1_j * V_i_j_jp1);
			}
		}

		double sum_weights = 0;
		for (auto& w : weights)
			sum_weights += w;
		*/
	}

};