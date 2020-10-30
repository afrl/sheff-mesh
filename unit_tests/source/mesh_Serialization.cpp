#include "mesh/Mesh.h"
#include "mesh/algorithms/Creation.h"
#include "utility/serialization/Serialization.h"

#include <iostream>
#include <sstream>

using namespace mesh;

int main()
{
    Mesh linear_grid = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });
    
    // Write the mesh in a JSON format to cout
    {
        utility::JSONSerializer to_json_cout(std::cout);
        to_json_cout(utility::named_ref("mesh", linear_grid));
        // Writes to output stream on destruction
    }
    
    // Count how many bytes would be needed to write to binary
    utility::SerialByteCounter counter;
    counter(linear_grid);
    std::cout << "# of bytes: " << counter.n_bytes() << std::endl;

    // Write in binary format to a string stream
    std::stringstream ss;
    utility::OutputSerializer to_binary_ss(ss);
    to_binary_ss(linear_grid);
    std::cout << "Start binary representation:\n" << ss.str() << "\nEnd binary representation" << std::endl;

    // Deserialize a mesh from the binary representation
    utility::InputDeserializer input_from_ss(ss);
    Mesh linear_grid_deserialized;
    input_from_ss(linear_grid_deserialized);

    assert(linear_grid.nodes().size() == linear_grid_deserialized.nodes().size());
    assert(linear_grid.elements().size() == linear_grid_deserialized.elements().size());

    return 0;
}