#include "mesh/external_interfaces/bsam_classic/IO.h"

#include "mesh/external_interfaces/abaqus/OptionsBlock.h"
#include "mesh/external_interfaces/abaqus/TypeConversions.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/StringUtils.h"

namespace mesh
{
static bool read_next_cluster(std::ifstream& input_stream,
                              LocalIndex& n_nodes_before_prev_block,
                              LocalIndex& n_elements_before_prev_block,
                              Mesh& mesh)
{
    std::string line;
    while (!input_stream.eof())
    {
        // Note: TODO add support for sections
        OptionBlock block(input_stream, line);
        if (block.block_type == OptionBlock::BlockType::NODE)
        {
            n_nodes_before_prev_block = mesh.nodes().size();

            // Read data lines
            getline(input_stream, line);
            while (line.find("*") == std::string::npos && line != "")
            {
                // Replace all commas with spaces so we have only whitespace delimiters
                std::replace(line.begin(), line.end(), ',', ' ');
                auto tokens = utility::tokenize(line);

                simple_math::Point p = { atof(tokens[1].c_str()),
                                         atof(tokens[2].c_str()),
                                         atof(tokens[3].c_str()) };

                mesh.nodes().add(p);

                getline(input_stream, line);
            }
        }
        else if (block.block_type == OptionBlock::BlockType::ELEMENT)
        {
            assert(block.keywords.count("TYPE") == 1);

            const auto& bsam_element_type = block.keywords["TYPE"];
            // TODO - we might need a specialized conversion for X3D8 and LC3D8
            auto element_type = convert_abaqus_type(bsam_element_type);

            n_elements_before_prev_block = mesh.elements().size();

            // Read data lines
            getline(input_stream, line);
            while (line.find("*") == std::string::npos && line != "")
            {
                // Replace all commas with spaces so we have only whitespace delimiters
                std::replace(line.begin(), line.end(), ',', ' ');
                auto tokens = utility::tokenize(line);

                std::vector<LocalIndex> connectivity;
                for (auto i : utility::range<std::size_t>(1, tokens.size()))
                {
                    // Adjust for 0-based and for the fact that the node numbering resets
                    // for each entity/instance
                    auto node_instance_index =
                      atoi(tokens[i].c_str()) - 1 + n_nodes_before_prev_block;
                    connectivity.push_back(node_instance_index);
                }

                mesh.elements().add(element_type, connectivity);
                getline(input_stream, line);
            }
        }
        else if (block.block_type == OptionBlock::BlockType::NODESET)
        {
            assert(block.keywords.count("NSET") == 1);

            auto node_set_name = block.keywords["NSET"];
            auto& node_group   = mesh.nodes().group(node_set_name);

            if (block.has_option("GENERATE"))
            {
                getline(input_stream, line);
                // Replace all commas with spaces so we have only whitespace delimiters
                std::replace(line.begin(), line.end(), ',', ' ');
                auto tokens = utility::tokenize(line);
                assert(tokens.size() == 3);
                int start = atoi(tokens[0].c_str());
                int end = atoi(tokens[1].c_str());
                int increment = atoi(tokens[2].c_str());
                for (auto i = start; i <= end; i += increment)
                {
                    // Adjust for 0-based and for the fact that the node numbering resets
                    // for each entity/instance
                    auto instance_index = i - 1 + n_nodes_before_prev_block;
                    node_group.node_indices.push_back(instance_index);
                }
                getline(input_stream, line);
            }
            else
            {
                // Read data lines
                getline(input_stream, line);
                while (line.find("*") == std::string::npos && line != "")
                {
                    // Replace all commas with spaces so we have only whitespace delimiters
                    std::replace(line.begin(), line.end(), ',', ' ');
                    auto tokens = utility::tokenize(line);

                    for (auto& i : tokens)
                    {
                        // Adjust for 0-based and for the fact that the node numbering resets
                        // for each entity/instance
                        auto instance_index = atoi(i.c_str()) - 1 + n_nodes_before_prev_block;
                        node_group.node_indices.push_back(instance_index);
                    }
                    getline(input_stream, line);
                }
            }
        }
        else if (block.block_type == OptionBlock::BlockType::ELEMENTSET)
        {
            assert(block.keywords.count("ELSET") == 1);

            auto element_set_name = block.keywords["ELSET"];
            auto& element_group   = mesh.elements().group(element_set_name);

            // Read data lines
            getline(input_stream, line);
            while (line.find("*") == std::string::npos && line != "")
            {
                // Replace all commas with spaces so we have only whitespace delimiters
                std::replace(line.begin(), line.end(), ',', ' ');
                auto tokens = utility::tokenize(line);

                for (auto& i : tokens)
                {
                    // Adjust for 0-based and for the fact that the element numbering resets
                    // for each entity/instance
                    auto instance_index = atoi(i.c_str()) - 1 + n_elements_before_prev_block;
                    element_group.element_indices.push_back(instance_index);
                }
                getline(input_stream, line);
            }
        }
        else if (block.block_type == OptionBlock::BlockType::STOP)
        {
            return true;
        }
        // Note for orientations:
        // * name=ORI means node-based
        // * name=ORI-ELE means element-based
        // * node/element #, fiber direction, 2nd vector (through-thickness), fiber Vf
        // (placeholder)
        else
        {
            // Move to next line
            getline(input_stream, line);
        }

        // Move on to next block
        line = utility::skip_to_line_that_contains("*", input_stream, line);
    }
    return false;
}

Mesh read_bsam_classic(const std::string& file_path)
{
    Mesh mesh;

    std::ifstream input_stream(file_path);
    utility::skip_to_line_that_contains("APPROXIMATION", input_stream);

    LocalIndex n_nodes_before_prev_block    = 0;
    LocalIndex n_elements_before_prev_block = 0;

    while (read_next_cluster(input_stream,
                             n_nodes_before_prev_block,
                             n_elements_before_prev_block,
                             mesh))
    {
    }

    input_stream.close();

    return mesh;
}

std::vector<Mesh> mesh::read_bsam_classic_by_cluster(const std::string& file_path)
{
    std::vector<Mesh> meshes;

    std::ifstream input_stream(file_path);
    utility::skip_to_line_that_contains("APPROXIMATION", input_stream);

    LocalIndex n_nodes_before_prev_block    = 0;
    LocalIndex n_elements_before_prev_block = 0;
    bool read_successful                    = true;

    while (read_successful)
    {
        Mesh mesh;
        n_nodes_before_prev_block    = 0;
        n_elements_before_prev_block = 0;
        read_successful              = read_next_cluster(input_stream,
                                            n_nodes_before_prev_block,
                                            n_elements_before_prev_block,
                                            mesh);
        if (read_successful)
        {
            meshes.push_back(std::move(mesh));
        }
    }

    input_stream.close();

    return meshes;
}
}