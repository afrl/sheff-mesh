#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "mesh/Mesh.h"

using namespace mesh;

int main()
{
    Mesh mesh;
    mesh.nodes().add({ 0, 0, 0 });
    mesh.nodes().add({ 1, 0, 0 });
    mesh.nodes().add({ 1, 1, 0 });
    mesh.nodes().add({ 0, 1, 0 });
    mesh.elements().add(ElementType::Q2D4N, { 0, 1, 2, 3 });

	ElementAccessor a = mesh.elements()[0];

	assert(mesh.nodes()[0].global_index() == 0);
	assert(mesh.nodes()[0].coord() == Eigen::Vector3d({ 0, 0, 0 }));
    assert(mesh.nodes()[3].global_index() == 3);
    assert(mesh.nodes()[3].coord() == Eigen::Vector3d({ 0, 1, 0 }));

	assert(mesh.elements()[0].type() == ElementType::Q2D4N);
    assert(mesh.elements()[0].node(0).global_index() == 0);
    assert(mesh.elements()[0].node(1).global_index() == 1);
    assert(mesh.elements()[0].node(2).global_index() == 2);
    assert(mesh.elements()[0].node(3).global_index() == 3);

	return 0;
}