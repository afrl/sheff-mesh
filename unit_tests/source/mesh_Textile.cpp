#include "NaturalSort.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Search.h"
#include "mesh/algorithms/Textile.h"
#include "unit_test_helper.h"
#include "utility/Filesystem.h"
#include "utility/StringUtils.h"

#include <fstream>

using namespace mesh;

void exe_code()
{
    auto tow_ctw_filepaths = get_files_with_extension("./", ".ctw", false);

    Mesh combined_mesh;
    std::ofstream tetwild_tags_file("tags.txt");
    int curr_tow_index = 0;
    for (auto& tow_ctw_filepath : tow_ctw_filepaths)
    {
        auto tow = read_vtms_ctw(tow_ctw_filepath);
        write_vtk(tow.mesh, fs::path(tow_ctw_filepath).stem().string() + "_poly.vtk");
        triangulate_tow_surface(tow);
        write_to_ply(tow.mesh, fs::path(tow_ctw_filepath).stem().string() + ".ply");
        combined_mesh.copy(tow.mesh);

        for (auto& t : tow.mesh.elements())
        {
            tetwild_tags_file << curr_tow_index << std::endl;
        }
        curr_tow_index++;
    }

    /* Add matrix
    auto extent = find_extent(combined_mesh);
    double matrix_factor = 0.8;
    double d             = (1.0 - matrix_factor) * (extent.second[1] - extent.first[1]) / 2.0;
    extent.first[1] += d;
    extent.second[1] -= d;
    auto box_mesh = create_triangular_matrix_box(extent.first, extent.second);
    combined_mesh.copy(box_mesh);

    for (auto& t : box_mesh.elements())
    {
        tetwild_tags_file << curr_tow_index << std::endl;
    }
    */

    write_to_ply(combined_mesh, "combined.ply");
    tetwild_tags_file.close();
}

int main()
{
    auto tow_1 = read_vtms_ctw(get_test_data_filepath("tow_1.ctw"));
    auto tow_2 = read_vtms_ctw(get_test_data_filepath("tow_2.ctw"));
    write_vtk(tow_1.mesh, "tow_1_pretri.vtk");
    write_vtk(tow_2.mesh, "tow_2_pretri.vtk");
    triangulate_tow_surface(tow_1);
    triangulate_tow_surface(tow_2);
    write_to_ply(tow_1.mesh, "tow_1.ply");
    write_to_ply(tow_2.mesh, "tow_2.ply");

    auto extent = find_extent({ find_extent(tow_1.mesh), find_extent(tow_2.mesh) });

    double matrix_factor = 0.8;
    double d             = (1.0 - matrix_factor) * (extent.second[1] - extent.first[1]) / 2.0;
    extent.first[1] += d;
    extent.second[1] -= d;
    auto box_mesh = create_triangular_matrix_box(extent.first, extent.second);
    write_to_ply(box_mesh, "matrix.ply");

    Mesh combined_mesh;
    combined_mesh.copy(tow_1.mesh);
    combined_mesh.copy(tow_2.mesh);
    combined_mesh.copy(box_mesh);
    write_to_ply(combined_mesh, "combined.ply");

    std::ofstream tetwild_tags_file("tags.txt");
    for (auto& t : tow_1.mesh.elements())
    {
        tetwild_tags_file << "0" << std::endl;
    }
    for (auto& t : tow_2.mesh.elements())
    {
        tetwild_tags_file << "1" << std::endl;
    }
    for (auto& t : box_mesh.elements())
    {
        tetwild_tags_file << "2" << std::endl;
    }
    tetwild_tags_file.close();

    return 0;
}