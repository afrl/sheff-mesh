#include "mesh/data/MeshAssociatedData.h"
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "utility/serialization/Serialization.h"

#include <assert.h>

using namespace mesh;

int main()
{
    Mesh grid = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

    assert(grid.nodes().size() == 15);
    assert(grid.elements().size() == 8);

    simple_math::DynamicDoubleMatrix node_data    = simple_math::DynamicDoubleMatrix::Random(15, 3);
    simple_math::DynamicDoubleMatrix element_data = simple_math::DynamicDoubleMatrix::Random(8, 3);

    std::vector<double> node_data2 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

    MeshAssociatedData data;

    data.add<DataFieldType::NODE, double>(
      DataField<double>::copy("rand_node_data", { "x", "y", "z" }, node_data));

    data.add<DataFieldType::ELEMENT, double>(
      DataField<double>::ref("rand_element_data", { "x", "y", "z" }, element_data));

    data.add<DataFieldType::NODE, double>(DataField<double>::ref("mapped_node_data", node_data2));

    _UTIL_ASSERT2(data.get<DataFieldType::NODE, double>().size() == 2);
    _UTIL_ASSERT2(data.get<DataFieldType::ELEMENT, double>().size() == 1);

    node_data = simple_math::DynamicDoubleMatrix::Random(15, 3);
    _UTIL_ASSERT2(data.get<DataFieldType::NODE, double>("rand_node_data").data() != node_data);

    element_data = simple_math::DynamicDoubleMatrix::Random(8, 3);
    _UTIL_ASSERT2(data.get<DataFieldType::ELEMENT, double>("rand_element_data").data() ==
                  element_data);

    std::ofstream json_out("data.json");
    utility::JSONSerializer s(json_out);
    s(utility::named_ref("mesh_associated_data", data));

    utility::save("mesh_associated_data_dump", utility::named_ref("mesh_associated_data", data));

    MeshAssociatedData data_copy;
    utility::load("mesh_associated_data_dump", 0, { "mesh_associated_data" }, data_copy);

    assert(data == data_copy);

    return 0;
}