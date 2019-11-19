#pragma once

#include "utility/Exceptions.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/Range.h"

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include "flex_mesh/algorithms/Search.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace flex_mesh
{
	/// Write a VTK-type unstructured mesh
	template<typename MeshType>
	static void write_vtk(MeshType& mesh,
		const std::string& file_path,
		const bool write_partition_ids = true,
		const bool write_lcs = true,
        const std::vector<double> & point_data = std::vector<double>());

    template<typename MeshType>
    class AbaqusInterface
    {
    public:
        using NodePtr = typename MeshType::NodeTypename*;
        using ElementPtr = typename MeshType::ElementTypename*;

        static void read(MeshType& mesh, const std::string& file_path);

    private:
        struct AbaqusElementConversion
        {
            ElementType flex_mesh_type;

            std::vector<NodePtr> convert_abaqus_to_flex_connectivity(const std::vector<NodePtr> & connectivity)
            {
                return connectivity;
            }
        };

        static AbaqusElementConversion abaqus_to_flex_type(const std::string & abaqus_type);
    };


	/////////////////
	// Definitions //
	/////////////////

	template<typename MeshType>
	inline void write_vtk(MeshType & mesh, const std::string & file_path, const bool write_partition_ids, const bool write_lcs, const std::vector<double> & point_data)
	{
		std::ofstream fout(file_path.c_str());
		if (!fout)
			throw _UTIL_MSG_EXCEPT("Could not create file " + file_path);

		//mesh.update_global_numbering();

		fout << "# vtk DataFile Version 2.0" << std::endl;
		fout << "Mesh from BetaMesh utility" << std::endl;
		fout << "ASCII" << std::endl;
		fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
		fout << "POINTS " << mesh.nodes().size() << " double" << std::endl;

        LocalIndex curr_local_index = 0;
		for (auto& n_ptr : mesh.nodes())
        {
            if (n_ptr->get_local_number() != curr_local_index)
            {
                throw _UTIL_MSG_EXCEPT("The local numbers of the nodes are not consecutive.  Renumber the nodes (local) before writing to VTK or the connectivity will be messed up.");
            }
			fout << n_ptr->get_coord()[0] << " " << n_ptr->get_coord()[1] << " " << n_ptr->get_coord()[2] << "" << std::endl;
            curr_local_index++;
		}

		int total_num_cell_conn_entries = 0;
		for (auto& e_ptr : mesh.elements())
        {
			total_num_cell_conn_entries++;
			total_num_cell_conn_entries += e_ptr->get_expanded_local_conn().size();
		}

		fout << "CELLS " << mesh.elements().size() << " " << total_num_cell_conn_entries << std::endl;

		for (auto& e_ptr : mesh.elements())
        {
            auto& local_numbered_conn = e_ptr->get_expanded_local_conn();
			auto vtk_permutation = get_vtk_permutation(
                e_ptr->get_type(), local_numbered_conn.size());
			if (vtk_permutation.size() > 0)
            {
				// Note: connectivity may have uint8_t type, which will be interpretted as a char
				// by ostream, resulting in an ASCII character instead of the integer. "+0" forces
				// it to be interpretted as an integer.
				fout << local_numbered_conn.size() + 0;
				for (unsigned int i = 0; i < local_numbered_conn.size(); ++i)
                {
					fout << " " << local_numbered_conn[vtk_permutation[i]];
				}
			}
			else {
				// How do we handle polygons and polyhedrons?
                throw _UTIL_MSG_EXCEPT("Element type unsupported.");
			}
			fout << std::endl;
		}

		fout << "CELL_TYPES " << mesh.elements().size() << std::endl;

		for (auto& e_ptr : mesh.elements())
        {
			fout << get_vtk_element_type(e_ptr->get_type()) << std::endl;
		}

        // Write data
        bool wrote_point_data_header = false;

		// Write the component id's as cell data
		// Note: the cell data has to all be together
		fout << "CELL_DATA " << mesh.elements().size() << std::endl;
		fout << "SCALARS ComponentIDs int 1\n";
		fout << "LOOKUP_TABLE default\n";
		for (auto& e_ptr : mesh.elements()) {
			fout << e_ptr->get_component_id() << std::endl;
		}

		// Write orientation data
		if (write_lcs) {
			fout << "VECTORS LocalMaterialX float\n";
			for (auto& e_ptr : mesh.elements()) {
				Vector3D axis_in_global = e_ptr->get_lcs(0)._transformVector(Vector3D::Unit(0));
				for (auto i = 0; i < axis_in_global.size(); i++) fout << axis_in_global(i) << " ";
				fout << std::endl;
			}
			fout << "VECTORS LocalMaterialY float\n";
			for (auto& e_ptr : mesh.elements()) {
				Vector3D axis_in_global = e_ptr->get_lcs(0)._transformVector(Vector3D::Unit(1));
				for (auto i = 0; i < axis_in_global.size(); i++) fout << axis_in_global(i) << " ";
				fout << std::endl;
			}
			fout << "VECTORS LocalMaterialZ float\n";
			for (auto& e_ptr : mesh.elements()) {
				Vector3D axis_in_global = e_ptr->get_lcs(0)._transformVector(Vector3D::Unit(2));
				for (auto i = 0; i < axis_in_global.size(); i++) fout << axis_in_global(i) << " ";
				fout << std::endl;
			}
		}

        fout << "SCALARS GlobalIDs int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& e_ptr : mesh.elements())
        {
            fout << e_ptr->get_global_number() << std::endl;
        }

		// Write the partition id's if asked for
		if (write_partition_ids) {
			fout << "SCALARS PartitionIDs int 1\n";
			fout << "LOOKUP_TABLE default\n";
			for (auto& e_ptr : mesh.elements()) {
				fout << e_ptr->get_parent_partition() << std::endl;
			}

            if (!wrote_point_data_header)
            {
                fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
                wrote_point_data_header = true;
            }
			fout << "SCALARS PartitionIDs int 1\n";
			fout << "LOOKUP_TABLE default\n";
			for (auto& n_ptr : mesh.nodes()) {
				fout << n_ptr->get_parent_partition() << std::endl;
			}
		}

        if (!wrote_point_data_header)
        {
            fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
            wrote_point_data_header = true;
        }
        fout << "SCALARS GlobalIDs int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto& n_ptr : mesh.nodes())
        {
            fout << n_ptr->get_global_number() << std::endl;
        }

        if (!point_data.empty())
        {
            if (point_data.size() != mesh.nodes().size())
            {
                throw _UTIL_MSG_EXCEPT("The size of the given point data did not match the number of nodes.");
            }

            if (!wrote_point_data_header)
            {
                fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
                wrote_point_data_header = true;
            }
            fout << "SCALARS ExtraData float 1\n";
            fout << "LOOKUP_TABLE default\n";
            for (auto& value : point_data)
            {
                fout << value << std::endl;
            }
        }

		fout.close();
	}

    template<typename MeshType>
    void AbaqusInterface<MeshType>::read(MeshType & mesh, const std::string & input_file_path)
    {
        std::ifstream input_stream(input_file_path);

        GlobalIndex global_node_count = 0;
        GlobalIndex global_element_count = 0;

        auto line = utility::skip_to_line_that_contains("*Node", input_stream);

        // Read nodes
        std::vector<typename MeshType::NodeTypename *> instance_node_ptrs;
        getline(input_stream, line);
        do
        {
            // Replace all commas with spaces so we have only whitespace delimiters
            std::replace(line.begin(), line.end(), ',', ' ');
            auto tokens = utility::tokenize(line);

            simple_math::Point p = { atof(tokens[1].c_str()), atof(tokens[2].c_str()), atof(tokens[3].c_str()) };

            auto new_node = mesh.add_node_by_coord(p);
            new_node->set_local_number(global_node_count);
            new_node->set_global_number(global_node_count);
            global_node_count++;

            instance_node_ptrs.push_back(new_node);

            getline(input_stream, line);
        } while (line.find("*Element") == std::string::npos);

        // Get the type of elements
        auto tokens = utility::tokenize(line);
        auto abaqus_element_type = tokens[1];
        utility::remove_substring(abaqus_element_type, "type=");
        auto abaqus_element_conversion = abaqus_to_flex_type(abaqus_element_type);

        // Read elements
        getline(input_stream, line);
        do
        {
            // Replace all commas with spaces so we have only whitespace delimiters
            std::replace(line.begin(), line.end(), ',', ' ');
            auto tokens = utility::tokenize(line);

            auto new_element = mesh.elements().add();
            new_element->set_local_number(global_element_count);
            new_element->set_global_number(global_element_count);
            global_element_count++;

            std::vector<typename MeshType::NodeTypename *> connectivity;
            for (auto i : utility::range<std::size_t>(1, tokens.size()))
            {
                // Adjust for 0-based
                auto node_instance_index = atoi(tokens[i].c_str()) - 1;
                connectivity.push_back(instance_node_ptrs[node_instance_index]);
            }
            new_element->initialize(abaqus_element_conversion.flex_mesh_type, connectivity);

            getline(input_stream, line);
        } while (line.find("*") == std::string::npos);

        // TODO read multiple instances

        input_stream.close();
    }

    template<typename MeshType>
    inline typename AbaqusInterface<MeshType>::AbaqusElementConversion AbaqusInterface<MeshType>::abaqus_to_flex_type(const std::string & abaqus_type)
    {
        if (utility::istr_is_equal(abaqus_type, "C3D8R"))
        {
            return { ElementType::H3D8N };
        }
        throw _UTIL_MSG_EXCEPT("Conversion not written for Abaqus type: " + abaqus_type);
    }

}