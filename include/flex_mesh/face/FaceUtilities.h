#pragma once

#include "flex_mesh/face/Face.h"

#include <map>
#include <set>

namespace flex_mesh
{
	/// Determine if a matching face (same centroid) is in a set of other faces
	template<typename FaceType>
	static bool is_face_in_set(const std::set<FaceType*, ConnHashComparison<FaceType>>& face_set);


	/// Get a set of all nodes that are common to the two faces
	template<typename FaceType, typename NodeType>
	static std::set<NodeType> get_common_nodes(const FaceType& a, const FaceType& b);

	/// Get a map of all modes that are coincident (but are not the same node)
	template<typename FaceType, typename NodeType>
	static std::map<NodeType*, NodeType*> get_coincident_nodes(const FaceType& a, const FaceType& b);

	/// Get pairs of coincident nodes, whether or not they are the same node
	template<typename FaceType, typename NodeType>
	static std::map<NodeType*, NodeType*> get_node_pairs(const FaceType& a, const FaceType& b);

	/// Return true if this face and OtherFace have all the same nodes
	bool check_matched(const FacePtr &OtherFace) const;

	/////////////////
	// Definitions //
	/////////////////

	template<typename FaceType>
	bool is_face_in_set(const std::set<FaceType*, ConnHashComparison<FaceType>>& face_set)
	{
		return face_set.find(this->shared_from_this()) != face_set.end();
	}

	template<typename FaceType, typename NodeType>
	std::set<NodeType*> get_common_nodes(const FaceType & a, const FaceType & b)
	{
		std::set<NodeType*> a_nodes(a.begin(), a.end());
		std::set<NodeType*> b_nodes(b.begin(), b.end());
		std::set<NodeType*> result;

		std::set_intersection(a_nodes.begin(), a_nodes.end(),
			b_nodes.begin(), b_nodes.end(),
			std::inserter(result, result.end()));
		return result;
	}

	template<typename FaceType, typename NodeType>
	std::map<NodeType*, NodeType*> get_coincident_nodes(const FaceType & a, const FaceType & b)
	{
		double old_tol = NodeType::get_compare_tolerance();
		NodeType::set_compare_tolerance(1e-15); // Set a really tight tolerance
		std::set<NodeType*, ObjectComparison<NodeType>> b_nodes(b.begin(), b.end());
		std::map<NodeType*, NodeType*> result;
		for (auto a_node_ptr : a) {
			auto b_find = b_nodes.find(a_node_ptr);
			if (b_find != b_nodes.end()) {
				result[a_node_ptr] = *b_find;
			}
		}
		NodeType::set_compare_tolerance(old_tol); // Restore tolerance
		return result;
	}

	template<typename FaceType, typename NodeType>
	std::map<NodeType*, NodeType*> get_node_pairs(const FaceType & a, const FaceType & b)
	{
		old_tol = NodeType::get_compare_tolerance();
		NodeType::set_compare_tolerance(1e-15); // Set a really tight tolerance
		std::map<NodeType*, NodeType*> node_pairs;
		std::set<NodeType*, ObjectComparison<NodeType>> a_nodes(a.begin(), a.end());
		for (auto b_node_ptr : b) {
			auto a_coincident_node_find = a_nodes.find(b_node_ptr);
			if (a_coincident_node_find == a_nodes.end())
				throw _FM_MSG_EXCEPT("Trying to get coincident node pairs for faces which don't have all coincident nodes.");
			else
				node_pairs[*a_coincident_node_find] = b_node_ptr;
		}
		return node_pairs;
	}

}