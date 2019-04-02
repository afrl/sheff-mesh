#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "flex_mesh/common/Common.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/algorithms/IO.h"

using namespace flex_mesh;

int main()
{
	initialize_flex_mesh(1e-14);

	using NodeTyped = Node;
	using ElementTyped = Element<NodeTyped>;
	using MeshType = Mesh<NodeTyped, ElementTyped>;
	
	MeshType mesh;
	AbaqusInterface<MeshType>::read(mesh, "BarWithHole.inp");

	assert(mesh.nodes().size() == 5724);
	assert(mesh.elements().size() == 3596);
	
	FEAIOInterface<MeshType>::write(mesh, "./", "BarWithHole", "BarWithHole");

	return 0;
}