#pragma once

#include "utility/Exceptions.h"
#include "Defines.h"

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
	static ConnIndex get_num_nodes_for_type(const EdgeType& geom_type);

    /// A simple struct to hold information about an edge that is part of a face/element.
    struct EdgeInfo
    {
        EdgeType type;
        std::vector<ConnIndex> connectivity;
    };

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
	static ConnIndex get_num_nodes_for_type(const FaceType& geom_type);
	/// Given an face geometry type, this function will fill the vectors:
	/// - types: with the types of edges
	/// - conn: with connectivity of each face relative to the face connectivity
	/// \attention Will overwrite anything in the vectors!
	static const std::vector<EdgeInfo> get_edge_info(const FaceType& geom_type);
	/// Given a geometry type, this function will return the linear version
	/// if one exists.  Otherwise it will return UNKNGEOM.
	static FaceType get_linearized_type(const FaceType& geom_type);
	/// Returns the permuation vector to go from the quadratic connectivity
	/// to the linear connectivity.  If the type is already linear, then
	/// an empty vector is returned.  Will throw for unrecognized types.
	static std::vector<ConnIndex> get_linearized_permutation(const FaceType& geom_type);
	/// Returns a bool indicating if the geometry type is linear or not
	static bool is_type_linear(const FaceType& geom_type);

    /// A simple struct to hold information about a face that is part
    /// of a element.
    struct FaceInfo
    {
        FaceType type;
        std::vector<ConnIndex> edges;
        std::vector<ConnIndex> connectivity;
    };

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
	enum class ElementType : ConnIndex {
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
		C3D16N = 48 ///< quadrilateral 16-node cohesive element (quadratic)
	};
	/// The maximum value that a ElementType enumeration can have.  This is only
	/// used for classes like Mapping, which use the underlying type as an index.
	static constexpr ConnIndex MAX_ELEMENT_TYPE_VALUE = 48;

	/// Returns the number of nodes that the geometry type will have.
	static ConnIndex get_num_nodes_for_type(const ElementType& geom_type);

    /// This method will use the number of dimensions and number of nodes to infer
    /// the type of solid element.
    static ElementType infer_type(unsigned int dimensionality, ConnIndex num_nodes);

	/// Given a geometry type, this function will return the linear version
	/// if one exists.  Otherwise it will return UNKNGEOM.
	ElementType get_linearized_type(const ElementType& geom_type);
	/// Returns the permuation vector to go from the quadratic connectivity
	/// to the linear connectivity.  If the type is already linear, then
	/// an empty vector is returned.  Will throw for unrecognized types.
    static std::vector<ConnIndex> get_linearized_permutation(const ElementType& geom_type) { return std::vector<ConnIndex>(); }
	/// Returns a bool indicating if the geometry type is linear or not
	static bool is_type_linear(const ElementType& geom_type);

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
	std::vector<ConnIndex> get_quadratic_indices_of_linear_nodes(const ElementType& geom_type);
	/// Given an element geometry type, this function will fill the vectors:
	/// - quadraticIndices: the indicides of hte midside nodes in the quadratic connectivity vector
	/// - linearNodePairIndices: the corresponding pair the midside node lies on in terms of the linear
	///   connectivity vector.
	/// \attention Will overwrite anything in the vectors!
	void get_quadratic_indices_of_midside_nodes(
		const ElementType& geom_type,
		std::vector<ConnIndex>& quadratic_indices,
		std::vector<std::pair<ConnIndex, ConnIndex>>& linear_node_pair_indices);

    /// Given an element geometry type, returns the edges of the element
    static const std::vector<EdgeInfo> get_edge_info(const ElementType& geom_type);
    /// Given an element geometry type, returns the faces of the element
    static const std::vector<FaceInfo> get_face_info(const ElementType& geom_type);

	/////////////
	// Helpers //
	/////////////

	/// A helper used in get_edge_info and get_face_info
	static std::vector<std::vector<ConnIndex>> connectivity_helper(
        const ConnIndex num_faces,
        const ConnIndex* nodes_per_face,
        const ConnIndex* face_nodes);

	/////////////////
	// VTK Helpers //
	/////////////////

	/// Returns the integer used to indicate VTK cell types for a given ElementType
	static int get_vtk_element_type(const ElementType& geom_type);
	/// Returns the permuation vector to go from the internal connectivity to the
	/// connectivity used in VTK.  All types are not mappable to a VTK element,
	/// in those cases, it will throw.
	static std::vector<ConnIndex> get_vtk_permutation(const ElementType& geom_type);

	//////////////////
	// Definitiions //
	//////////////////

    ConnIndex get_num_nodes_for_type(const EdgeType & geom_type)
	{
		if (geom_type == EdgeType::E2N) return 2;
		else if (geom_type == EdgeType::E3N) return 3;
		else throw _UTIL_MSG_EXCEPT("Edge type not supported");
	}

    ConnIndex get_num_nodes_for_type(const FaceType & geom_type)
	{
		if (geom_type == FaceType::L1D1N) return 1;
		else if (geom_type == FaceType::L2D2N) return 2;
		else if (geom_type == FaceType::L2D3N) return 3;
		else if (geom_type == FaceType::T3D3N) return 3;
		else if (geom_type == FaceType::Q3D4N) return 4;
		else if (geom_type == FaceType::T3D6N) return 6;
		else if (geom_type == FaceType::Q3D8N) return 8;
		else throw _UTIL_MSG_EXCEPT("Face type not supported");
	}

    const std::vector<EdgeInfo> get_edge_info(const FaceType& geom_type)
    {
        switch (geom_type)
        {
        // 2D Connectivities
        case FaceType::L2D2N:
            return {
                { EdgeType::E2N, {0, 1} }
            };
        case FaceType::L2D3N:
            return {
                { EdgeType::E3N, {0, 1, 2} }
            };
        // 3D Connectivities
        case FaceType::T3D3N:
            return {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 0} }
            };
        case FaceType::Q3D4N:
            return {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 3} },
                { EdgeType::E2N, {3, 0} }
            };
        case FaceType::T3D6N:
            return {
                { EdgeType::E3N, {0, 1, 2} },
                { EdgeType::E3N, {2, 3, 4} },
                { EdgeType::E3N, {4, 5, 0} }
            };
        case FaceType::Q3D8N:
            return {
                { EdgeType::E3N, {0, 1, 2} },
                { EdgeType::E3N, {2, 3, 4} },
                { EdgeType::E3N, {4, 5, 6} },
                { EdgeType::E3N, {6, 7, 0} }
            };
        default:
            throw _UTIL_MSG_EXCEPT("Edge connectivity not implemented.");
        }
    }

    FaceType get_linearized_type(const FaceType & geom_type)
	{
		// 2D
		if (geom_type == FaceType::L2D3N) return FaceType::L2D2N;
		// 3D
		else if (geom_type == FaceType::Q3D8N) return FaceType::Q3D4N;
		else if (geom_type == FaceType::T3D6N) return FaceType::T3D3N;
		else return FaceType::UNKNGEOM;
	}

	std::vector<ConnIndex> get_linearized_permutation(const FaceType & geom_type)
	{
		// We are going to cheat a bit here, since the permutations will act just like elements
		// 2D
		if (geom_type == FaceType::L2D3N) return get_linearized_permutation(ElementType::L1D3N);
		// 3D
		else if (geom_type == FaceType::Q3D8N) return get_linearized_permutation(ElementType::Q2D8N);
		else if (geom_type == FaceType::T3D6N) return get_linearized_permutation(ElementType::T2D6N);
		else throw _UTIL_MSG_EXCEPT("Unsupported face type.");
	}

	bool is_type_linear(const FaceType & geom_type)
	{
		return geom_type == FaceType::L1D1N ||
			geom_type == FaceType::L2D2N ||
			geom_type == FaceType::Q3D4N ||
			geom_type == FaceType::T3D3N;
	}

    const std::vector<EdgeInfo> get_edge_info(const ElementType& geom_type)
    {
        std::vector<EdgeInfo> edge_info;
        switch (geom_type)
        {
        // 2D Connectivities
        case ElementType::T2D3N:
            // Linear triangle
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 0} }
            };
            break;
        case ElementType::T2D6N:
            // Quadratic triangle
            edge_info = {
                { EdgeType::E3N, {0, 1, 2} },
                { EdgeType::E3N, {2, 3, 4} },
                { EdgeType::E3N, {4, 5, 0} }
            };
            break;
        case ElementType::Q2D4N:
            // Linear quadralateral
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 3} },
                { EdgeType::E2N, {3, 0} }
            };
            break;
        case ElementType::Q2D8N:
            // Quadratric quadralateral
            edge_info = {
                { EdgeType::E3N, {0, 1, 2} },
                { EdgeType::E3N, {2, 3, 4} },
                { EdgeType::E3N, {4, 5, 6} },
                { EdgeType::E3N, {6, 7, 0} }
            };
            break;
        case ElementType::C2D4N:
            // Linear 2D cohesive element
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {2, 3} }
            };
            break;
        case ElementType::C2D6N:
            // Quadratic 2D cohesive element
            edge_info = {
                { EdgeType::E2N, {0, 1, 2} },
                { EdgeType::E2N, {5, 4, 3} }
            };
            break;
        // 3D Connectivities
        case ElementType::T3D4N:
            // Linear tetrahedral
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {0, 2} },
                { EdgeType::E2N, {0, 3} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {1, 3} },
                { EdgeType::E2N, {2, 3} }
            };
            break;
        case ElementType::P3D5N:
            // Linear pyramid
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 3} },
                { EdgeType::E2N, {3, 0} },
                { EdgeType::E2N, {0, 4} },
                { EdgeType::E2N, {1, 4} },
                { EdgeType::E2N, {2, 4} },
                { EdgeType::E2N, {3, 4} }
            };
            break;
        case ElementType::W3D6N:
            // Linear wedge
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 0} },
                { EdgeType::E2N, {0, 3} },
                { EdgeType::E2N, {1, 4} },
                { EdgeType::E2N, {2, 5} },
                { EdgeType::E2N, {3, 4} },
                { EdgeType::E2N, {4, 5} },
                { EdgeType::E2N, {5, 3} }
            };
            break;
        case ElementType::H3D8N:
            // Linear hexahedral
            edge_info = {
                { EdgeType::E2N, {0, 1} },
                { EdgeType::E2N, {1, 2} },
                { EdgeType::E2N, {2, 3} },
                { EdgeType::E2N, {3, 0} },
                { EdgeType::E2N, {0, 4} },
                { EdgeType::E2N, {1, 5} },
                { EdgeType::E2N, {2, 6} },
                { EdgeType::E2N, {3, 7} },
                { EdgeType::E2N, {4, 5} },
                { EdgeType::E2N, {5, 6} },
                { EdgeType::E2N, {6, 7} },
                { EdgeType::E2N, {7, 4} }
            };
            break;
        case ElementType::T3D10N:
            // Quadratic tetrahedral
            edge_info = {
                { EdgeType::E3N, {0, 4, 1} },
                { EdgeType::E3N, {0, 6, 2} },
                { EdgeType::E3N, {0, 7, 3} },
                { EdgeType::E3N, {1, 5, 2} },
                { EdgeType::E3N, {1, 8, 3} },
                { EdgeType::E3N, {2, 9, 3} }
            };
            break;
        case ElementType::P3D13N:
            // Quadratic Pyramid
            edge_info = {
                { EdgeType::E3N, {0, 5, 1} },
                { EdgeType::E3N, {1, 6, 2} },
                { EdgeType::E3N, {2, 7, 3} },
                { EdgeType::E3N, {3, 8, 0} },
                { EdgeType::E3N, {0, 9, 4} },
                { EdgeType::E3N, {1, 10, 4} },
                { EdgeType::E3N, {2, 11, 4} },
                { EdgeType::E3N, {3, 12, 4} },
            };
            break;
        case ElementType::W3D15N:
            // Quadratic wedge
            edge_info = {
                { EdgeType::E3N, {0, 1, 2} },
                { EdgeType::E3N, {2, 3, 4} },
                { EdgeType::E3N, {4, 5, 6} },
                { EdgeType::E3N, {0, 6, 9} },
                { EdgeType::E3N, {2, 7, 11} },
                { EdgeType::E3N, {4, 8, 13} },
                { EdgeType::E3N, {9, 10, 11} },
                { EdgeType::E3N, {11, 12, 13} },
                { EdgeType::E3N, {13, 14, 9} }
            };
            break;
        case ElementType::H3D20N:
            // Quadratic hexahedral
            edge_info = {
                { EdgeType::E3N, {0, 1, 2} },
                { EdgeType::E3N, {2, 3, 4} },
                { EdgeType::E3N, {4, 5, 6} },
                { EdgeType::E3N, {6, 7, 0} },
                { EdgeType::E3N, {0, 8, 12} },
                { EdgeType::E3N, {2, 9, 14} },
                { EdgeType::E3N, {4, 10, 16} },
                { EdgeType::E3N, {6, 11, 18} },
                { EdgeType::E3N, {12, 13, 14} },
                { EdgeType::E3N, {14, 15, 16} },
                { EdgeType::E3N, {16, 17, 18} },
                { EdgeType::E3N, {18, 19, 12} }
            };
            break;
        default:
            throw _UTIL_MSG_EXCEPT("Edge connectivity not implemented.");
        }
        return edge_info;
    }

    const std::vector<FaceInfo> get_face_info(const ElementType& geom_type)
    {
        // All the different face definitions
        // Faces defined CCW from outside element and starting at a corner

        // Faces for solid elements
        // 2D Connectivities
        switch (geom_type)
        {
        case ElementType::T2D3N:
        case ElementType::T2D6N:
        case ElementType::Q2D4N:
        case ElementType::Q2D8N:
        case ElementType::C2D4N:
        case ElementType::C2D6N:
            throw _UTIL_MSG_EXCEPT("Faces do not exist for 2D elements, only edges.");
        // 3D Connectivitives
        case ElementType::T3D4N:
            // Linear tetrahedral
            return {
                { FaceType::T3D3N, {0, 3, 1}, {0, 2, 1} },
                { FaceType::T3D3N, {2, 5, 1}, {0, 3, 2} },
                { FaceType::T3D3N, {3, 5, 4}, {1, 2, 3} },
                { FaceType::T3D3N, {0, 4, 2}, {0, 1, 3} }
            };
        case ElementType::P3D5N:
            // Linear pyramid
            return {
                { FaceType::Q3D4N, {0, 1, 2, 3}, {0, 3, 2, 1} },
                { FaceType::T3D3N, {0, 5, 4}, {0, 1, 4} },
                { FaceType::T3D3N, {1, 6, 5}, {1, 2, 4} },
                { FaceType::T3D3N, {2, 7, 6}, {2, 3, 4} },
                { FaceType::T3D3N, {3, 4, 7}, {3, 0, 4} }
            };
        case ElementType::W3D6N:
            // Linear wedge
            return {
                { FaceType::T3D3N, {0, 2, 1}, {0, 2, 1} },
                { FaceType::T3D3N, {6, 7, 8}, {3, 4, 5} },
                { FaceType::Q3D4N, {0, 4, 6, 3}, {0, 1, 4, 3} },
                { FaceType::Q3D4N, {3, 8, 5, 2}, {0, 3, 5, 2} },
                { FaceType::Q3D4N, {1, 5, 7, 4}, {1, 2, 5, 4} }
            };
        case ElementType::H3D8N:
            // Linear hexahedral
            return {
                { FaceType::Q3D4N, {3, 2, 1, 0}, {0, 3, 2, 1} },
                { FaceType::Q3D4N, {4, 11, 7, 3}, {0, 4, 7, 3} },
                { FaceType::Q3D4N, {0, 5, 8, 4}, {0, 1, 5, 4} },
                { FaceType::Q3D4N, {8, 9, 10, 11}, {4, 5, 6, 7} },
                { FaceType::Q3D4N, {1, 6, 9, 5}, {1, 2, 6, 5} },
                { FaceType::Q3D4N, {2, 7, 10, 6}, {2, 3, 7, 6} }
            };
        case ElementType::T3D10N:
            // Quadratic tetrahedral
            return {
                { FaceType::T3D6N, {0, 3, 1}, {0, 6, 2, 5, 1, 4} },
                { FaceType::T3D6N, {2, 5, 1}, {0, 7, 3, 9, 2, 6} },
                { FaceType::T3D6N, {3, 5, 4}, {1, 5, 2, 9, 3, 8} },
                { FaceType::T3D6N, {0, 4, 2}, {0, 4, 1, 8, 3, 7} }
            };
        case ElementType::P3D13N:
            // quadratic pyramid
            return {
                { FaceType::Q3D8N, {0, 1, 2, 3}, {0,8,3,7,2,6,1,5} },
                { FaceType::T3D6N, {0, 5, 4}, {0,5,1,10,4,9} },
                { FaceType::T3D6N, {1, 6, 5}, {1,6,2,11,4,10} },
                { FaceType::T3D6N, {2, 7, 6}, {2,7,3,12,4,11} },
                { FaceType::T3D6N, {3, 4, 7}, {3,8,0,9,4,12} }
            };
        case ElementType::W3D15N:
            // Quadratic wedge
            return {
                { FaceType::T3D6N, {0, 2, 1}, {0,5,4,3,2,1} },
                { FaceType::T3D6N, {6, 7, 8}, {9,10,11,12,13,14} },
                { FaceType::Q3D8N, {0, 4, 6, 3}, {0,1,2,7,11,10,9,6} },
                { FaceType::Q3D8N, {3, 8, 5, 2}, {0,6,9,14,13,8,4,5} },
                { FaceType::Q3D8N, {1, 5, 7, 4}, {2,3,4,8,13,12,11,7} }
            };
        case ElementType::H3D20N:
            // Quadratic hexahedral
            return {
                { FaceType::Q3D8N, {3, 2, 1, 0}, {0,7,6,5,4,3,2,1} },
                { FaceType::Q3D8N, {4, 11, 7, 3}, {0,8,12,19,18,11,6,7} },
                { FaceType::Q3D8N, {0, 5, 8, 4}, {0,1,2,9,14,13,12,8} },
                { FaceType::Q3D8N, {8, 9, 10, 11}, {12,13,14,15,16,17,18,19} },
                { FaceType::Q3D8N, {1, 6, 9, 5}, {2,3,4,10,16,15,14,9} },
                { FaceType::Q3D8N, {2, 7, 10, 6}, {4,5,6,11,18,17,16,10} }
            };
        case ElementType::C3D8N:
            // Linear CZ
            //int FN_3_8[] = { 3,2,1,0,
            //    4,5,6,7 };
        case ElementType::C3D16N:
            // Quadratic CZ
            //int FN_3_16[] = { 7,6,5,4,3,2,1,0,
            //    8,9,10,11,12,13,14,15 };
        default:
            throw _UTIL_MSG_EXCEPT("Face connectivity not implemented.");
        }
    }

	std::vector<std::vector<ConnIndex>> connectivity_helper(
        const ConnIndex num_faces, const ConnIndex* nodes_per_face,
        const ConnIndex* face_nodes)
	{
		std::vector<std::vector<ConnIndex>> faces_connectivities(num_faces);
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

	int get_vtk_element_type(const ElementType & geom_type)
	{
		// 1D
		if (geom_type == ElementType::L1D2N)
			return 3;
		else if (geom_type == ElementType::L1D3N)
			return 21;
		// 2D
		else if (geom_type == ElementType::T2D3N)
			return 5;
		else if (geom_type == ElementType::Q2D4N)
			return 9;
		else if (geom_type == ElementType::T2D6N)
			return 22;
		else if (geom_type == ElementType::Q2D8N)
			return 23;
		//else if (geom_type == ElementType::POLYGON)
		//	return 7;
		// 3D
		else if (geom_type == ElementType::T3D4N)
			return 10;
		else if (geom_type == ElementType::P3D5N)
			return 14;
		else if (geom_type == ElementType::W3D6N)
			return 13;
		else if (geom_type == ElementType::H3D8N)
			return 12;
		else if (geom_type == ElementType::T3D10N)
			return 24;
		else if (geom_type == ElementType::P3D13N)
			return 27;
		else if (geom_type == ElementType::W3D15N)
			return 26;
		else if (geom_type == ElementType::H3D20N)
			return 25;
		//else if (geom_type == ElementType::POLYHEDRON)
		//	return ?;
		else
			throw _UTIL_MSG_EXCEPT("VTK type conversion not implemented.");
	}

	std::vector<ConnIndex> get_vtk_permutation(const ElementType & geom_type)
	{
		auto num_nodes = get_num_nodes_for_type(geom_type);
		std::vector<ConnIndex> permutation(num_nodes);
		// Fill permutation with 0 - num_nodes-1
		std::iota(permutation.begin(), permutation.end(), 0);

		// 1D
		if (geom_type == ElementType::L1D2N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::L1D3N)
			throw _UTIL_MSG_EXCEPT("L1D2N type not implemented.");
		// 2D
		else if (geom_type == ElementType::T2D3N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::Q2D4N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::T2D6N)
			return{ 0,2,4,1,3,5 }; // Different
		else if (geom_type == ElementType::Q2D8N)
			return{ 0,2,4,6,1,3,5,7 }; // Different
		//else if (geom_type == ElementType::POLYGON)
		//	return permutation; // Unchanged
		// 3D
		else if (geom_type == ElementType::T3D4N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::P3D5N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::W3D6N)
			return{ 0, 2, 1, 3, 5, 4 }; // Different
		else if (geom_type == ElementType::H3D8N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::T3D10N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::P3D13N)
			return permutation; // Unchanged
		else if (geom_type == ElementType::W3D15N)
			return{ 0, 4, 2, 9, 13, 11, 5, 3, 1, 14, 12, 10, 6, 8, 7 }; // Different
		else if (geom_type == ElementType::H3D20N)
			return{ 0, 2, 4, 6, 12, 14, 16, 18, 1, 3, 5, 7, 13, 15, 17, 19, 8, 9, 10, 11 }; // Different
		else
			return{};
	}

    ConnIndex get_num_nodes_for_type(const ElementType & geom_type)
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

    ElementType infer_type(unsigned int dimensionality, ConnIndex num_nodes)
    {
        
        switch (dimensionality)
        {
        case 1:
            switch (num_nodes)
            {
            case 2: return ElementType::L1D2N;
            case 3: return ElementType::L1D3N;
            default:
                throw _UTIL_MSG_EXCEPT("Cannot infer element type for 1D with " + std::to_string(num_nodes) + " nodes.");
            }
        case 2:
            switch (num_nodes)
            {
            case 3: return ElementType::T2D3N;
            case 4: return ElementType::Q2D4N;
            case 6: return ElementType::T2D6N;
            case 8: return ElementType::Q2D8N;
            default:
                throw _UTIL_MSG_EXCEPT("Cannot infer element type for 2D with " + std::to_string(num_nodes) + " nodes.");
            }
        case 3:
            switch (num_nodes)
            {
            case 4: return ElementType::T3D4N;
            case 5: return ElementType::P3D5N;
            case 6: return ElementType::W3D6N;
            case 8: return ElementType::H3D8N;
            case 10: return ElementType::T3D10N;
            case 13: return ElementType::P3D13N;
            case 15: return ElementType::W3D15N;
            case 20: return ElementType::H3D20N;
            default:
                throw _UTIL_MSG_EXCEPT("Cannot infer element type for 3D with " + std::to_string(num_nodes) + " nodes.");
            }
        default:
            throw _UTIL_MSG_EXCEPT("Invalid dimensionality.");
        }
    }

	bool is_type_linear(const ElementType & geom_type)
	{
		// 1D
		if (geom_type == ElementType::L1D2N) return true;
		// 2D
		else if (geom_type == ElementType::T2D3N) return true;
		else if (geom_type == ElementType::Q2D4N) return true;
		else if (geom_type == ElementType::C2D4N) return true;
		// 3D
		else if (geom_type == ElementType::T3D4N) return true;
		else if (geom_type == ElementType::P3D5N) return true;
		else if (geom_type == ElementType::W3D6N) return true;
		else if (geom_type == ElementType::H3D8N) return true;
		else if (geom_type == ElementType::Q3D4N) return true;
		else if (geom_type == ElementType::C3D8N) return true;

		// No match
		return false;
	}
};
