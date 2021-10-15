#include "mesh/kd_tree/KDTree.h"

namespace mesh
{
	NodeContainerKDTree::NodeContainerKDTree(const NodeContainer& nodes)
		: adapter_(nodes)
		, index_(3 /*dim*/, adapter_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */))
	{
		index_.buildIndex();
	}

	NodeConstAccessor NodeContainerKDTree::closest_node(const Point& p) const
	{
		const size_t num_results = 1;

		size_t ret_index;
		double out_dist_sqr;
		nanoflann::KNNResultSet<double> result_set(num_results);
		result_set.init(&ret_index, &out_dist_sqr);

		index_.findNeighbors(result_set, &p[0], nanoflann::SearchParams(10));

		return adapter_.nodes_[static_cast<LocalIndex>(ret_index)];
	}

	IteratorRange<NodeSubrangeConstIterator> NodeContainerKDTree::closest_nodes(const Point& p,
		const int n) const
	{
		const size_t num_results = n;

		std::vector<size_t> ret_indices(num_results);
		std::vector<double> out_dists_sqr(num_results);
		nanoflann::KNNResultSet<double> result_set(num_results);
		result_set.init(&ret_indices[0], &out_dists_sqr[0]);

		index_.findNeighbors(result_set, &p[0], nanoflann::SearchParams(10));

		std::vector<mesh::LocalIndex> cast_ret_indices(num_results);
		for (auto i : utility::indices(ret_indices))
		{
			cast_ret_indices[i] = static_cast<LocalIndex>(ret_indices[i]);
		}

		return {
			NodeSubrangeConstIterator(adapter_.nodes_,
									  cast_ret_indices.data(),
									  static_cast<LocalIndex>(num_results)), // begin
			NodeSubrangeConstIterator(adapter_.nodes_,
									  cast_ret_indices.data(),
									  static_cast<LocalIndex>(num_results),
									  static_cast<LocalIndex>(num_results)) // end
		};
	}

	IteratorRange<NodeSubrangeConstIterator> NodeContainerKDTree::in_radius(
		const Point& p,
		const double& search_radius) const
	{
		std::vector<std::pair<size_t, double>> ret_matches;

		nanoflann::SearchParams params;
		// params.sorted = false;

		// Note: the documentation for radiusSearch seems to be incorrect.  Upon testing
		// the parameter "radius" seems to actually be the square radius.  Internally,
		// nanoflann is using L2-norm as the metric since we used the L2_Simple_Adaptor
		// (which is of course the square of the distance).
		const size_t n_matches =
			index_.radiusSearch(&p[0], search_radius * search_radius, ret_matches, params);

		std::vector<mesh::LocalIndex> cast_ret_indices(n_matches);
		for (auto i : utility::range(n_matches))
		{
			cast_ret_indices[i] = static_cast<LocalIndex>(ret_matches[i].first);
		}

		return {
			NodeSubrangeConstIterator(adapter_.nodes_,
									  cast_ret_indices.data(),
									  static_cast<LocalIndex>(n_matches)), // begin
			NodeSubrangeConstIterator(adapter_.nodes_,
									  cast_ret_indices.data(),
									  static_cast<LocalIndex>(n_matches),
									  static_cast<LocalIndex>(n_matches)) // end
		};
	}
}