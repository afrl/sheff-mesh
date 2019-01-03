#pragma once

#include "Defines.h"

#include <vector>

namespace flex_mesh
{

	/// A templated struct to use with STL containers that store a shared_ptr of something
	/// and allows the container to sort the shared_ptr's using the objects in the shared_ptr
	/// comparison operator (<).
	/// Example use: set<shared_ptr<NodeTopology>, ObjectComparison<NodeTopology>> will
	/// sort the set and maintain uniqueness based on the NodeTopology::<, which uses a
	/// tolerance based comparison.  This will maintain that the set only has one node per
	/// a location.  In constrast, set<shared_ptr<NodeTopology>> will keep unique node
	/// objects based on address (two nodes are free to have the same coordinate within
	/// the set).
	template<class T>
	struct ObjectComparison {
		bool operator() (const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const
		{
			return *lhs < *rhs;
		}
		bool operator() (const T* lhs, const T* rhs) const
		{
			return *lhs < *rhs;
		}
	};

	/// Basically, just like ObjectComparison but add comparison by number if they are the same.
	template<class T>
	struct ObjectComparisonWithNumber {
		bool operator() (const T* lhs, const T* rhs) const
		{
			if (*lhs < *rhs)
				return true;
			else if (*rhs < *lhs == false) // Share the same coordinate
				return lhs->get_local_number() < rhs->get_local_number();
			else
				return false;
		}
	};

	/// Similar to ObjectComparison, but the sorted connectivity is hashed.  Consequently,
	/// if used with a set, then it will result in an entry for each object without duplicates
	/// that have the same connectivity (even if the order is different).
	/// This is fantastic for removing matched faces and getting the outer set, even when CZ's
	/// are used (which can have the same centroid but you want the one with the unique connectivity
	/// for MPC's).
	template<typename T>
	struct ConnHashComparison {
		// The less than operator
		bool operator() (const T* lhs, const T* rhs) const
		{
			return lhs->get_conn_hash_id() < rhs->get_conn_hash_id();
		}
		bool operator() (const T& lhs, const T& rhs) const
		{
			return lhs.get_conn_hash_id() < rhs.get_conn_hash_id();
		}
	};

	/// The operator of this struct will hash the values in the vector
	/// and combine the hashes.  This is very similar to Boost's hash_combine.
	/// Should be used with things like vector<double>
	template<class T>
	struct VectorValueHasher {
		static std::size_t hash(std::vector<T> const& vec) {
			std::size_t seed = vec.size();
			for (auto& i : vec) {
				seed ^= std::hash<std::size_t>{}((std::size_t) i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};

	/// The operator of this struct will hash the address of the pointers in the vector
	/// and combine the hashes.  This is very similar to Boost's hash_combine.
	/// Should be used with things like vector<shared_ptr<NodeTopology>>
	template<class T>
	struct VectorPtrHasher {
		static std::size_t hash(std::vector<T*> const& vec) {
			std::size_t seed = vec.size();
			for (auto& i : vec) {
				seed ^= std::hash<std::size_t>{}((std::size_t) &(*i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};

	/// Returns if a < b while taking into account a numerical tolerance
	static bool compare_points(const Point& a, const Point& b, const double tolerance);
	static Eigen::Map<const Vector3D> point_as_vector(const Point& p);

	/////////////////
	// Definitions //
	/////////////////

	bool compare_points(const Point & a, const Point & b, const double tolerance)
	{
		return (a[0] < b[0] - tolerance)
			|| ((fabs(a[0] - b[0]) < tolerance) && a[1] < b[1] - tolerance)
			|| ((fabs(a[0] - b[0]) < tolerance) && (fabs(a[1] - b[1]) < tolerance) && (a[2] < b[2] - tolerance));
	}

	Eigen::Map<const Vector3D> point_as_vector(const Point & p)
	{
		return Eigen::Map<const Vector3D>(&p[0], p.size());
	}
}
