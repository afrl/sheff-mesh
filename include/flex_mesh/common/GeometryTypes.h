#pragma once

#include "Exceptions.h"

#include <vector>

namespace flex_mesh
{
	////////////////
	// Edge Types //
	////////////////

	/// Describes the different types of edges in terms of the geometry/topology.
	/// Note: edges are a bit weird in 2D or 1D.  In 1D/2D the edge will be same
	///       as a 1D/2D face basically.
	enum class EdgeType : char {
		UNKNGEOM = 0, ///< Used if a conversion fails
		E2N = 10, ///< 2-node edge
		E3N = 11, ///< 3-node edge
	};

	/// Returns the number of nodes that the geometry type will have.
	static unsigned int get_num_nodes_for_type(const EdgeType& geom_type);

	////////////////
	// Face Types //
	////////////////

	/// Describes the different types of faces in terms of the geometry/topology.
	/// Note: faces are a bit weird in 1D.  In 1D, a face will consist of a single node.
	enum class FaceType : char {
		UNKNGEOM = 0, ///< Used if a conversion fails
		// 1D types 10-19
		L1D1N = 11, ///< 1-D Linear 1-node
		// 2D types 20-24 is linear and 25-29 is quadratic
		L2D2N = 21, ///< 2-D 2-node face (linear)
		L2D3N = 22, ///< 2-D 3-node face (quadratic)
		// Linear 3D 30-39
		T3D3N = 31, ///< triangular 3-node face (linear)
		Q3D4N = 32, ///< quadrilateral 4-node face (linear)
		// Quadratic 3D 41-49
		T3D6N = 41, ///< triangular 6-node face (quadratic)
		Q3D8N = 42, ///< quadrilateral 8-node face (quadratic)
	};

	/// Returns the number of nodes that the geometry type will have.
	static unsigned int get_num_nodes_for_type(const FaceType& geom_type);
	/// Given an face geometry type, this function will fill the vectors:
	/// - types: with the types of edges
	/// - conn: with connectivity of each face relative to the face connectivity
	/// \attention Will overwrite anything in the vectors!
	static void get_edge_info(const FaceType& geom_type, std::vector<EdgeType>& types, std::vector<std::vector<unsigned int>>& conn);
	/// Given a geometry type, this function will return the linear version
	/// if one exists.  Otherwise it will return UNKNGEOM.
	static FaceType get_linearized_type(const FaceType& geom_type);
	/// Returns the permuation vector to go from the quadratic connectivity
	/// to the linear connectivity.  If the type is already linear, then
	/// an empty vector is returned.  Will throw for unrecognized types.
	static std::vector<unsigned int> get_linearized_permutation(const FaceType& geom_type);
	/// Returns a bool indicating if the geometry type is linear or not
	static bool is_type_linear(const FaceType& geom_type);

	///////////////////
	// Element Types //
	///////////////////

	/// Dev note: if you add new ones, do not change the number of old types so that old
	/// mesh files remain compatible.
	/// IMPORTANT NOTE: if you add a new element type that you want to be able to vizualize
	/// using ParaView, you need to update the fea_io library, in any related functions.
	/// Additionally, you need to update the ParaViewReader project file vtkFEAIOReader.cxx methods
	/// get_to_vtk_cell_map and get_to_vtk_reordering.
	// TODO check names for ambiguity
	enum class ElementType : char {
		UNKNGEOM = 0, ///< Used if a conversion fails
		// 1D types 10-19
		DEF1D = 10, ///< 1-D Default
		L1D2N = 12, ///< 1-D Linear 2-node linear element
		L1D3N = 13, ///< 1-D Linear 3-node linear element
		// 2D types 20-24 is linear and 25-29 is quadratic
		DEF2D = 20, ///< 2-D Default
		T2D3N = 21, ///< triangular 2D 3-node element (linear)
		Q2D4N = 22, ///< quadrilateral 4-node linear element
		C2D4N = 24, ///< line 4-node cohesive element (linear)
		T2D6N = 25, ///< triangular 2D 6-node element (quadratic)
		Q2D8N = 26, ///< quadrilateral 8-node quadratic element (quadratic)
		C2D6N = 28, ///< line 6-node cohesive element (quadratic)
		// Linear 3D 30-39
		DEF3D = 30, ///< 3-D Default
		T3D4N = 31, ///< 3D tetrahedral 4-node  (linear) element
		P3D5N = 35, ///< 3D pyramid 5-node  (linear) element
		W3D6N = 32, ///< 3D wedge 6-node  (linear) element
		H3D8N = 33, ///< 3D hexahedral 8-node (linear) element
		Q3D4N = 36, ///< quadrilateral shell element  (linear)
		C3D8N = 38, ///< quadrilateral 8-node cohesive element (linear)
		POLYHEDRON = 39, ///< general 3D N-sided polyhedron
		// Quadratic 3D 41-49
		T3D10N = 41, ///< 3D tetrahedral 10-node  (quadratic) element
		P3D13N = 45, ///< 3D pyramid 13-node (quadratic) element
		W3D15N = 42, ///< 3D wedge 15-node  (quadratic) element
		H3D20N = 43, ///< 3D hexahedral 20-node (quadratic) element
		Q3D8N = 46, ///< quadrilateral shell element  (quadratic)
		C3D16N = 48, ///< quadrilateral 16-node cohesive element (quadratic)
	};

	/// Returns the number of nodes that the geometry type will have.
	static unsigned int get_num_nodes_for_type(const ElementType& geom_type);

	/// Given a geometry type, this function will return the linear version
	/// if one exists.  Otherwise it will return UNKNGEOM.
	ElementType get_linearized_type(const ElementType& geom_type);
	/// Returns the permuation vector to go from the quadratic connectivity
	/// to the linear connectivity.  If the type is already linear, then
	/// an empty vector is returned.  Will throw for unrecognized types.
	std::vector<unsigned int> get_linearized_permutation(const ElementType& geom_type);
	/// Returns a bool indicating if the geometry type is linear or not
	bool is_type_linear(const ElementType& geom_type);

	/// Returns the corresponding quadratic type for a given linear type
	ElementType get_quadratic_type(const ElementType& geom_type);
	/// For a given type of linear element, this will return the connectivity index
	/// map from linear to quadratic for the nodes that already exist (corner nodes)
	/// Example for T2D3N, the vector will be:
	/// [0]: 0
	/// [1]: 2
	/// [2]: 4
	/// which means that the first node in the linear connectivity belongs at index 0
	/// in the quadratic connectivity, etc.
	std::vector<unsigned int> get_quadratic_indices_of_linear_nodes(const ElementType& geom_type);
	/// Given an element geometry type, this function will fill the vectors:
	/// - quadraticIndices: the indicides of hte midside nodes in the quadratic connectivity vector
	/// - linearNodePairIndices: the corresponding pair the midside node lies on in terms of the linear
	///   connectivity vector.
	/// \attention Will overwrite anything in the vectors!
	void get_quadratic_indices_of_midside_nodes(
		const ElementType& geom_type,
		std::vector<unsigned int>& quadratic_indices,
		std::vector<std::pair<unsigned int, unsigned int>>& linear_node_pair_indices);

	/// Given an element geometry type, this function will fill the vectors:
	/// - types: with the types of faces
	/// - conn: with connectivity of each face relative to the element connectivity
	/// \attention Will overwrite anything in the vectors!
	void get_face_info(const ElementType& geom_type, const unsigned int num_nodes_in_element, std::vector<FaceType>& types, std::vector<std::vector<unsigned int>>& conn);
	/// Given an element geometry type, this function will fill the vectors:
	/// - types: with the types of edges
	/// - conn: with connectivity of each element relative to the element connectivity
	/// \attention Will overwrite anything in the vectors!
	void get_edge_info(const ElementType & geom_type, const unsigned int num_nodes_in_element, std::vector<EdgeType>& types, std::vector<std::vector<unsigned int>>& conn);


	/////////////
	// Helpers //
	/////////////

	/// A helper used in get_edge_info and get_face_info
	static std::vector<std::vector<unsigned int>> connectivity_helper(int num_faces, int* nodes_per_face, int* face_nodes);

	//////////////////
	// Definitiions //
	//////////////////

	unsigned int flex_mesh::get_num_nodes_for_type(const EdgeType & geom_type)
	{
		if (geom_type == EdgeType::E2N) return 2;
		else if (geom_type == EdgeType::E3N) return 3;
		else throw _FM_MSG_EXCEPT("Edge type not supported");
	}

	unsigned int flex_mesh::get_num_nodes_for_type(const FaceType & geom_type)
	{
		if (geom_type == FaceType::L1D1N) return 1;
		else if (geom_type == FaceType::L2D2N) return 2;
		else if (geom_type == FaceType::L2D3N) return 3;
		else if (geom_type == FaceType::T3D3N) return 3;
		else if (geom_type == FaceType::Q3D4N) return 4;
		else if (geom_type == FaceType::T3D6N) return 6;
		else if (geom_type == FaceType::Q3D8N) return 8;
		else throw _FM_MSG_EXCEPT("Face type not supported");
	}

	void flex_mesh::get_edge_info(const FaceType & geom_type, std::vector<EdgeType>& types, std::vector<std::vector<unsigned int>>& conn)
	{
		conn.clear();
		// 2D Connectivities
		if (geom_type == FaceType::L2D2N) {
			types = std::vector<EdgeType>(1, EdgeType::E2N);
			conn.push_back({ 0, 1 });
		}
		else if (geom_type == FaceType::L2D3N) {
			types = std::vector<EdgeType>(1, EdgeType::E3N);
			conn.push_back({ 0, 1, 2 });
		}
		// 3D Connectivities
		else if (geom_type == FaceType::T3D3N) {
			int NumEdges = 3;
			int NodesPerEdge[] = { 2, 2, 2 };
			int EdgeIndices[] = { 0, 1,
				1, 2,
				2, 0 };
			types = std::vector<EdgeType>(NumEdges, EdgeType::E2N);
			conn = connectivity_helper(NumEdges, NodesPerEdge, EdgeIndices);
		}
		else if (geom_type == FaceType::Q3D4N) {
			int NumEdges = 4;
			int NodesPerEdge[] = { 2, 2, 2, 2 };
			int EdgeIndices[] = { 0, 1,
				1, 2,
				2, 3,
				3, 0 };
			types = std::vector<EdgeType>(NumEdges, EdgeType::E2N);
			conn = connectivity_helper(NumEdges, NodesPerEdge, EdgeIndices);
		}
		else if (geom_type == FaceType::T3D6N) {
			int NumEdges = 3;
			int NodesPerEdge[] = { 3, 3, 3 };
			int EdgeIndices[] = { 0, 1, 2,
				2, 3, 4,
				4, 5, 0 };
			types = std::vector<EdgeType>(NumEdges, EdgeType::E3N);
			conn = connectivity_helper(NumEdges, NodesPerEdge, EdgeIndices);
		}
		else if (geom_type == FaceType::Q3D8N) {
			int NumEdges = 4;
			int NodesPerEdge[] = { 3, 3, 3, 3 };
			int EdgeIndices[] = { 0, 1, 2,
				2, 3, 4,
				4, 5, 6,
				6, 7, 0 };
			types = std::vector<EdgeType>(NumEdges, EdgeType::E3N);
			conn = connectivity_helper(NumEdges, NodesPerEdge, EdgeIndices);
		}
		else
			throw _FM_MSG_EXCEPT("Edge connectivity not implemented.");
	}

	FaceType flex_mesh::get_linearized_type(const FaceType & geom_type)
	{
		// 2D
		if (geom_type == FaceType::L2D3N) return FaceType::L2D2N;
		// 3D
		else if (geom_type == FaceType::Q3D8N) return FaceType::Q3D4N;
		else if (geom_type == FaceType::T3D6N) return FaceType::T3D3N;
		else return FaceType::UNKNGEOM;
	}

	std::vector<unsigned int> flex_mesh::get_linearized_permutation(const FaceType & geom_type)
	{
		// We are going to cheat a bit here, since the permutations will act just like elements
		// 2D
		if (geom_type == FaceType::L2D3N) return get_linearized_permutation(ElementType::L1D3N);
		// 3D
		else if (geom_type == FaceType::Q3D8N) return get_linearized_permutation(ElementType::Q2D8N);
		else if (geom_type == FaceType::T3D6N) return get_linearized_permutation(ElementType::T2D6N);
		else throw _FM_MSG_EXCEPT("Unsupported face type.");
	}

	bool flex_mesh::is_type_linear(const FaceType & geom_type)
	{
		return geom_type == FaceType::L1D1N ||
			geom_type == FaceType::L2D2N ||
			geom_type == FaceType::Q3D4N ||
			geom_type == FaceType::T3D3N;
	}

	std::vector<std::vector<unsigned int>> flex_mesh::connectivity_helper(int num_faces, int* nodes_per_face, int* face_nodes)
	{
		std::vector<std::vector<unsigned int>> faces_connectivities(num_faces);
		int FaceNodesIndex = 0;
		for (int face = 0; face < num_faces; face++) {
			faces_connectivities[face].resize(nodes_per_face[face]);
			for (int node = 0; node < nodes_per_face[face]; node++) {
				faces_connectivities[face][node] = face_nodes[FaceNodesIndex];
				FaceNodesIndex++;
			}
		}
		return faces_connectivities;
	}

	unsigned int flex_mesh::get_num_nodes_for_type(const ElementType & geom_type)
	{
		// 1D
		if (geom_type == ElementType::L1D2N) return 2;
		else if (geom_type == ElementType::L1D3N) return 3;
		// 2D
		else if (geom_type == ElementType::T2D3N) return 3;
		else if (geom_type == ElementType::Q2D4N) return 4;
		else if (geom_type == ElementType::C2D4N) return 4;
		else if (geom_type == ElementType::T2D6N) return 6;
		else if (geom_type == ElementType::Q2D8N) return 8;
		else if (geom_type == ElementType::C2D6N) return 6;
		// 3D
		else if (geom_type == ElementType::T3D4N) return 4;
		else if (geom_type == ElementType::P3D5N) return 5;
		else if (geom_type == ElementType::W3D6N) return 6;
		else if (geom_type == ElementType::H3D8N) return 8;
		else if (geom_type == ElementType::Q3D4N) return 4;
		else if (geom_type == ElementType::C3D8N) return 8;
		else if (geom_type == ElementType::T3D10N) return 10;
		else if (geom_type == ElementType::P3D13N) return 13;
		else if (geom_type == ElementType::W3D15N) return 15;
		else if (geom_type == ElementType::H3D20N) return 20;
		else if (geom_type == ElementType::Q3D8N) return 8;
		else if (geom_type == ElementType::C3D16N) return 16;
		// No match
		else return 0;
	}
};