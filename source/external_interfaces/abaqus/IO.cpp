#include "mesh/external_interfaces/abaqus/IO.h"

#include "mesh/external_interfaces/abaqus/OptionsBlock.h"
#include "mesh/external_interfaces/abaqus/TypeConversions.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/StringUtils.h"

namespace mesh
{
Mesh read_abaqus(const std::string& file_path)
{
    Mesh mesh;

    std::ifstream input_stream(file_path);

    GlobalIndex global_node_count    = 0;
    GlobalIndex global_element_count = 0;

    auto line = utility::skip_to_line_that_contains("*", input_stream);

    while (!input_stream.eof())
    {
        OptionBlock block(input_stream, line);
        if (block.block_type == OptionBlock::BlockType::NODE)
        {
            // Read data lines
            utility::get_next_line(input_stream, line);
            while (line.find("*") == std::string::npos && line.size() > 0)
            {
                // Replace all commas with spaces so we have only whitespace delimiters
                std::replace(line.begin(), line.end(), ',', ' ');
                auto tokens = utility::tokenize(line);

                simple_math::Point p = { atof(tokens[1].c_str()),
                                         atof(tokens[2].c_str()),
                                         atof(tokens[3].c_str()) };

                mesh.nodes().add(p);

                utility::get_next_line(input_stream, line);
            }
        }
        else if (block.block_type == OptionBlock::BlockType::ELEMENT)
        {
            assert(block.keywords.count("TYPE") == 1);

            auto abaqus_element_type = block.keywords["TYPE"];
            auto element_type        = convert_abaqus_type(abaqus_element_type);

            // Read data lines
            utility::get_next_line(input_stream, line);
            while (line.find("*") == std::string::npos && line != "")
            {
                // Replace all commas with spaces so we have only whitespace delimiters
                std::replace(line.begin(), line.end(), ',', ' ');
                auto tokens = utility::tokenize(line);

                std::vector<LocalIndex> connectivity;
                for (auto i : utility::range<std::size_t>(1, tokens.size()))
                {
                    // Adjust for 0-based
                    auto node_instance_index = atoi(tokens[i].c_str()) - 1;
                    connectivity.push_back(node_instance_index);
                }

                mesh.elements().add(element_type, connectivity);
                utility::get_next_line(input_stream, line);
            }
        }
        else
        {
            // Move on to next block
            line = utility::skip_to_line_that_contains("*", input_stream);
        }
    }

    input_stream.close();

    return mesh;
}
}