#pragma once

#include "flex_mesh/node/Node.h"

namespace flex_mesh
{
	/// Returns the mean value of the coordinates from a list of nodes
	template<typename NodePointerContainer>
	static Point get_mean_coordinate_of_nodes(const NodePointerContainer& nodes);
	/// Returns the distance between two nodes
	static double distance(const Node& a, const Node& b);
	static double distance(const Point& a, const Point& b);

	/////////////////
	// Definitions //
	/////////////////

	template<typename NodePointerContainer>
	Point get_mean_coordinate_of_nodes(const NodePointerContainer& nodes)
	{
		Point mean = { 0, 0, 0 };
		const FloatType num_nodes = static_cast<FloatType>(nodes.size());
		for (auto& n_ptr : nodes)
			for (auto j = 0; j < mean.size(); j++)
				mean[j] += (*n_ptr)[j] / num_nodes;
		return mean;
	}

	double distance(const Node & a, const Node & b)
	{
		return (a.as_eigen_vec3d() - b.as_eigen_vec3d()).norm();
	}

	double distance(const Point & a, const Point & b)
	{
		double sum_sq = 0;
		for (auto i = 0; i < a.size(); i++)
			sum_sq += (a[i] - b[i]) * (a[i] - b[i]);
		return sqrt(sum_sq);
	}
}