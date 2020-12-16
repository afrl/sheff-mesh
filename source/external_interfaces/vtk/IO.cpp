#include "mesh/external_interfaces/vtk/IO.h"

#include "mesh/external_interfaces/vtk/TypeConversions.h"

#include <fstream>

namespace mesh
{
void write_vtk(Mesh& mesh,
               const std::string& file_path,
               const bool write_lcs,
               const std::vector<double>& point_data,
               const std::vector<double>& element_data,
               const std::vector<unsigned int>& point_attributes,
               const std::vector<unsigned int>& element_attributes)
{
    std::ofstream fout(file_path.c_str());
    if (!fout)
    {
        throw _UTIL_MSG_EXCEPT("Could not create file " + file_path);
    }

    fout << "# vtk DataFile Version 2.0" << std::endl;
    fout << "Mesh exported from mesh module of SHEFF framework" << std::endl;
    fout << "ASCII" << std::endl;
    fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
    fout << "POINTS " << mesh.nodes().size() << " double" << std::endl;

    for (auto node : mesh.nodes())
    {
        fout << node.x() << " " << node.y() << " " << node.z() << std::endl;
    }

    int total_num_cell_conn_entries = 0;
    for (auto& element : mesh.elements())
    {
        total_num_cell_conn_entries++;
        total_num_cell_conn_entries += static_cast<int>(element.connectivity().size());
    }

    fout << "CELLS " << mesh.elements().size() << " " << total_num_cell_conn_entries << std::endl;

    for (auto& element : mesh.elements())
    {
        auto connectivity = element.connectivity();
        auto vtk_permutation =
          get_vtk_permutation(element.type(), static_cast<ConnIndex>(connectivity.size()));
        if (vtk_permutation.size() > 0)
        {
            // Note: connectivity may have uint8_t type, which will be interpretted as a char
            // by ostream, resulting in an ASCII character instead of the integer. "+0" forces
            // it to be interpretted as an integer.
            fout << connectivity.size() + 0;
            for (unsigned int i = 0; i < connectivity.size(); ++i)
            {
                fout << " " << connectivity[vtk_permutation[i]];
            }
        }
        else
        {
            // How do we handle polygons and polyhedrons?
            throw _UTIL_MSG_EXCEPT("Element type unsupported.");
        }
        fout << std::endl;
    }

    fout << "CELL_TYPES " << mesh.elements().size() << std::endl;

    for (auto& element : mesh.elements())
    {
        fout << get_vtk_element_type(element.type()) << std::endl;
    }

    // Write data
    bool wrote_point_data_header = false;

    // Note: the cell data has to all be together in the VTK file
    fout << "CELL_DATA " << mesh.elements().size() << std::endl;

    /* Note: I think component ID's have been phased out, so commenting out for now.
    fout << "SCALARS ComponentIDs int 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (auto element : mesh.elements())
    {
        fout << e_ptr->get_component_id() << std::endl;
    }
    */

    // Write orientation data
    if (write_lcs)
    {
        fout << "VECTORS LocalMaterialX float\n";
        for (auto& element : mesh.elements())
        {
            auto axis_in_global = element.local_coord_sys()[0]._transformVector(Vector3D::Unit(0));
            for (auto i = 0; i < axis_in_global.size(); i++)
                fout << axis_in_global(i) << " ";
            fout << std::endl;
        }
        fout << "VECTORS LocalMaterialY float\n";
        for (auto& element : mesh.elements())
        {
            auto axis_in_global = element.local_coord_sys()[0]._transformVector(Vector3D::Unit(1));
            for (auto i = 0; i < axis_in_global.size(); i++)
                fout << axis_in_global(i) << " ";
            fout << std::endl;
        }
        fout << "VECTORS LocalMaterialZ float\n";
        for (auto& element : mesh.elements())
        {
            auto axis_in_global = element.local_coord_sys()[0]._transformVector(Vector3D::Unit(2));
            for (auto i = 0; i < axis_in_global.size(); i++)
                fout << axis_in_global(i) << " ";
            fout << std::endl;
        }
    }

    fout << "SCALARS GlobalIDs int 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (auto& element : mesh.elements())
    {
        fout << element.global_index() << std::endl;
    }

    fout << "SCALARS PartitionIDs int 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (auto& element : mesh.elements())
    {
        fout << element.parent_partition() << std::endl;
    }

    fout << "SCALARS ComponentIDs int 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (auto& element : mesh.elements())
    {
        fout << element.component_id() << std::endl;
    }

    if (!element_data.empty())
    {
        if (element_data.size() != mesh.elements().size())
        {
            throw _UTIL_MSG_EXCEPT(
              "The size of the given cell data did not match the number of elements.");
        }

        fout << "SCALARS Data float 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& value : element_data)
        {
            fout << value << std::endl;
        }
    }

    if (!element_attributes.empty())
    {
        if (element_attributes.size() != mesh.elements().size())
        {
            throw _UTIL_MSG_EXCEPT(
              "The size of the given cell attributes did not match the number of elmeents.");
        }

        fout << "SCALARS Attributes unsigned_int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& value : element_attributes)
        {
            fout << value << std::endl;
        }
    }

    if (!wrote_point_data_header)
    {
        fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
        wrote_point_data_header = true;
    }
    fout << "SCALARS GlobalIDs int 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (auto node : mesh.nodes())
    {
        fout << node.global_index() << std::endl;
    }

    fout << "SCALARS PartitionIDs int 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (auto node : mesh.nodes())
    {
        fout << node.parent_partition() << std::endl;
    }

    if (!point_data.empty())
    {
        if (point_data.size() != mesh.nodes().size())
        {
            throw _UTIL_MSG_EXCEPT(
              "The size of the given point data did not match the number of nodes.");
        }

        if (!wrote_point_data_header)
        {
            fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
            wrote_point_data_header = true;
        }
        fout << "SCALARS Data float 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& value : point_data)
        {
            fout << value << std::endl;
        }
    }

    if (!point_attributes.empty())
    {
        if (point_attributes.size() != mesh.nodes().size())
        {
            throw _UTIL_MSG_EXCEPT(
              "The size of the given point attributes did not match the number of nodes.");
        }

        if (!wrote_point_data_header)
        {
            fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
            wrote_point_data_header = true;
        }
        fout << "SCALARS Attributes unsigned_int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& value : point_attributes)
        {
            fout << value << std::endl;
        }
    }

    fout.close();
}
}