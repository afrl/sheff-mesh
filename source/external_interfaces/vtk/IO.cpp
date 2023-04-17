#include "mesh/external_interfaces/vtk/IO.h"

#include "mesh/external_interfaces/vtk/TypeConversions.h"

#include <fstream>

namespace mesh
{
void write_vtk(const Mesh& mesh,
               const std::string& file_path,
               const bool write_lcs,
               const std::vector<double>& point_data,
               const std::vector<double>& element_data,
               const std::vector<unsigned int>& point_attributes,
               const std::vector<unsigned int>& element_attributes,
               const simple_math::GeneralMatrix& point_vector_data,
               const simple_math::GeneralMatrix& element_vector_data)
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

    // Note: the cell data has to all be together in the VTK file
    fout << "CELL_DATA " << mesh.elements().size() << std::endl;

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
              "The size of the given cell attributes did not match the number of elements.");
        }

        fout << "SCALARS Attributes unsigned_int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& value : element_attributes)
        {
            fout << value << std::endl;
        }
    }

    if (element_vector_data.size() > 0)
    {
        _UTIL_ASSERT(
          element_vector_data.rows() == mesh.elements().size(),
          "The size of the given cell vector data did not match the number of elements.");

        auto n_components = element_vector_data.cols();
        fout << "SCALARS VectorData float " + std::to_string(n_components) + "\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto i : utility::range(element_vector_data.rows()))
        {
            for (auto j : utility::range(element_vector_data.cols()))
            {
                fout << element_vector_data(i, j) << " ";
            }
            fout << std::endl;
        }
    }


    // Start of writing point data
    fout << "POINT_DATA " << mesh.nodes().size() << std::endl;

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

        fout << "SCALARS Attributes unsigned_int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& value : point_attributes)
        {
            fout << value << std::endl;
        }
    }

    if (point_vector_data.size() > 0)
    {
        _UTIL_ASSERT(point_vector_data.rows() == mesh.nodes().size(),
                     "The size of the given point vector data did not match the number of nodes.");

        auto n_components = point_vector_data.cols();
        fout << "SCALARS VectorData float " + std::to_string(n_components) + "\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto i : utility::range(point_vector_data.rows()))
        {
            for (auto j : utility::range(point_vector_data.cols()))
            {
                fout << point_vector_data(i, j) << " ";
            }
            fout << std::endl;
        }
    }

    fout.close();
}

void write_vtk(const Mesh& mesh, const std::string& file_path, const MeshAssociatedData& data)
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

    // Note: the cell data has to all be together in the VTK file
    fout << "CELL_DATA " << mesh.elements().size() << std::endl;

    // Write orientation data
    bool has_local_coordinate_systems = false;
    for (auto& element : mesh.elements())
    {
        has_local_coordinate_systems |=
          element.local_coord_sys().base_type() != LCSAssignmentType::None;
    }
    if (has_local_coordinate_systems)
    {
        fout << "VECTORS LocalX float\n";
        for (auto& element : mesh.elements())
        {
            auto axis_in_global = element.local_coord_sys()[0]._transformVector(Vector3D::Unit(0));
            for (auto i = 0; i < axis_in_global.size(); i++)
                fout << axis_in_global(i) << " ";
            fout << std::endl;
        }
        fout << "VECTORS LocalY float\n";
        for (auto& element : mesh.elements())
        {
            auto axis_in_global = element.local_coord_sys()[0]._transformVector(Vector3D::Unit(1));
            for (auto i = 0; i < axis_in_global.size(); i++)
                fout << axis_in_global(i) << " ";
            fout << std::endl;
        }
        fout << "VECTORS LocalZ float\n";
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

    if (!data.empty())
    {
        for (auto& data_field : data.get<DataFieldType::ELEMENT, double>())
        {
            _UTIL_ASSERT(data_field.data().cols() == data_field.component_names.size(),
                         "The # of columns for the cell data did not match the number of names.");

            fout << "SCALARS " << data_field.name << " float " << data_field.component_names.size()
                 << "\n";
            fout << "LOOKUP_TABLE default\n";
            _UTIL_ASSERT(data_field.data().rows() == mesh.elements().size(),
                         "The size of the given cell data did not match the number of elements.");
            auto& data_field_matrix = data_field.data();
            for (auto row : utility::range(data_field_matrix.rows()))
            {
                for (auto col : utility::range(data_field_matrix.cols()))
                {
                    fout << data_field_matrix(row, col) << " ";
                }
                fout << std::endl;
            }
        }

        for (auto& data_field : data.get<DataFieldType::ELEMENT, uint32_t>())
        {
            _UTIL_ASSERT(data_field.data().cols() == data_field.component_names.size(),
                         "The # of columns for the cell data did not match the number of names.");

            fout << "SCALARS " << data_field.name << " unsigned_int " << data_field.component_names.size()
                 << "\n";
            fout << "LOOKUP_TABLE default\n";
            _UTIL_ASSERT(data_field.data().rows() == mesh.elements().size(),
                         "The size of the given cell data did not match the number of elements.");
            auto& data_field_matrix = data_field.data();
            for (auto row : utility::range(data_field_matrix.rows()))
            {
                for (auto col : utility::range(data_field_matrix.cols()))
                {
                    fout << data_field_matrix(row, col) << " ";
                }
                fout << std::endl;
            }
        }
    }

    // Start of writing point data
    fout << "POINT_DATA " << mesh.nodes().size() << std::endl;

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

    if (!data.empty())
    {
        for (auto& data_field : data.get<DataFieldType::NODE, double>())
        {
            _UTIL_ASSERT(data_field.data().cols() == data_field.component_names.size(),
                         "The # of columns for the point data did not match the number of names.");

            fout << "SCALARS " << data_field.name << " float " << data_field.component_names.size()
                 << "\n";
            fout << "LOOKUP_TABLE default\n";
            _UTIL_ASSERT(data_field.data().rows() == mesh.nodes().size(),
                         "The size of the given point data did not match the number of nodes.");
            auto& data_field_matrix = data_field.data();
            for (auto row : utility::range(data_field_matrix.rows()))
            {
                for (auto col : utility::range(data_field_matrix.cols()))
                {
                    fout << data_field_matrix(row, col) << " ";
                }
                fout << std::endl;
            }
        }

        for (auto& data_field : data.get<DataFieldType::NODE, uint32_t>())
        {
            _UTIL_ASSERT(data_field.data().cols() == data_field.component_names.size(),
                         "The # of columns for the point data did not match the number of names.");

            fout << "SCALARS " << data_field.name << " unsigned_int "
                 << data_field.component_names.size() << "\n";
            fout << "LOOKUP_TABLE default\n";
            _UTIL_ASSERT(data_field.data().rows() == mesh.nodes().size(),
                         "The size of the given cell data did not match the number of nodes.");
            auto& data_field_matrix = data_field.data();
            for (auto row : utility::range(data_field_matrix.rows()))
            {
                for (auto col : utility::range(data_field_matrix.cols()))
                {
                    fout << data_field_matrix(row, col) << " ";
                }
                fout << std::endl;
            }
        }
    }

    fout.close();
}
}