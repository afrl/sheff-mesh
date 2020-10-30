#include "pfecDefines.h"

using namespace PFECNativeInterface;

std::map<int, std::string> PFECNativeInterface::get_elem_geom_name_map()
{
    std::map<int, std::string> nameMap;
    // 1D types 10-19
    nameMap[10] = "DEF1D"; // 1-D Default
    nameMap[12] = "L1D2N"; // 1-D Linear 2-node linear element
    nameMap[13] = "L1D3N"; // 1-D Linear 3-node linear element
                           // 2D types 20-24 is linear and 25-29 is quadratic
    nameMap[20] = "DEF2D"; // 2-D Default
    nameMap[21] = "T2D3N"; // triangular 2D 3-node element (linear)
    nameMap[22] = "Q2D4N"; // quadrilateral 4-node linear element
    nameMap[24] = "C2D4N"; // line 4-node cohesive element (linear)
    nameMap[25] = "T2D6N"; // triangular 2D 6-node element (quadratic)
    nameMap[26] = "Q2D8N"; // quadrilateral 8-node quadratic element (quadratic)
    nameMap[28] = "C2D6N"; // line 6-node cohesive element (quadratic)
                           // Linear 3D 30-39
    nameMap[30] = "DEF3D"; // 3-D Default
    nameMap[31] = "T3D4N"; // 3D tetrahedral 4-node  (linear) element
    nameMap[35] = "P3D5N"; // 3D tetrahedral 4-node  (linear) element
    nameMap[32] = "W3D6N"; // 3D wedge 6-node  (linear) element
    nameMap[33] = "H3D8N"; // 3D hexahedral 8-node (linear) element
    nameMap[36] = "Q3D4N"; // quadrilateral shell element  (linear)
    nameMap[38] = "C3D8N"; // quadrilateral 8-node cohesive element (linear)
                           // Quadratic 3D 41-49
    nameMap[41] = "T3D10N"; // 3D tetrahedral 10-node  (quadratic) element
    nameMap[45] = "P3D13N"; // 3D tetrahedral 10-node  (quadratic) element
    nameMap[42] = " W3D15N"; // 3D wedge 15-node  (quadratic) element
    nameMap[43] = "H3D20N"; // 3D hexahedral 20-node (quadratic) element
    nameMap[46] = "Q3D8N"; // quadrilateral shell element  (quadratic)
    nameMap[48] = "C3D16N"; // quadrilateral 16-node cohesive element (quadratic)
    return nameMap;
}

std::map<int, std::string> PFECNativeInterface::get_data_type_int_name_map()
{
    std::map<int, std::string> nameMap;
    nameMap[1] = "Nodal";
    nameMap[2] = "NodalByElement";
    nameMap[4] = "Element";
    nameMap[7] = "QuadPoint";
	nameMap[9] = "Orientation";
    return nameMap;
}

std::map<DataType, std::string> PFECNativeInterface::get_data_type_name_map()
{
	std::map<DataType, std::string> nameMap;
	nameMap[DataType::Nodal] = "Nodal";
	nameMap[DataType::NodalByElement] = "NodalByElement";
	nameMap[DataType::Element] = "Element";
	nameMap[DataType::QuadPoint] = "QuadPoint";
	nameMap[DataType::Orientation] = "Orientation";
	return nameMap;
}

std::map<std::string, DataType> PFECNativeInterface::get_name_data_type_map()
{
	std::map<std::string, DataType> nameMap;
	nameMap["Nodal"] = DataType::Nodal;
	nameMap["NodalByElement"] = DataType::NodalByElement;
	nameMap["Element"] = DataType::Element;
	nameMap["QuadPoint"] = DataType::QuadPoint;
	nameMap["Orientation"] = DataType::Orientation;
	return nameMap;
}

std::map<ElemGeomType, PFEC_INT> PFECNativeInterface::get_conn_num_node_map()
{
    std::map<ElemGeomType, PFEC_INT> connMap;
    // 1D types
    connMap[DEF1D] = 1;
    connMap[L1D2N] = 2;
    connMap[L1D3N] = 3;
    // 2D types
    connMap[DEF2D] = 1;
    connMap[T2D3N] = 3;
    connMap[Q2D4N] = 4;
    connMap[C2D4N] = 4;
    connMap[T2D6N] = 6;
    connMap[Q2D8N] = 8;
    connMap[C2D6N] = 6;
    // 3D types
    connMap[DEF3D] = 1;
    connMap[T3D4N] = 4;
    connMap[P3D5N] = 5;
    connMap[W3D6N] = 6;
    connMap[H3D8N] = 8;
    connMap[Q3D4N] = 4;
    connMap[C3D8N] = 8;
    connMap[T3D10N] = 10;
    connMap[P3D13N] = 13;
    connMap[W3D15N] = 15;
    connMap[H3D20N] = 20;
    connMap[Q3D8N] = 8;
    connMap[C3D16N] = 16;
    return connMap;
}