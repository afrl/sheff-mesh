#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Connectivity.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/node/DistribNode.h"
#include "flex_mesh/node/NodeUtilities.h"
#include "flex_mesh/edge/Edge.h"
#include "flex_mesh/face/Face.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/element/DistribElement.h"
#include "flex_mesh/element/MeshGenElement.h"
#include "flex_mesh/element/PolyhedronElement.h"
#include "flex_mesh/mesh/Mesh.h"

using namespace flex_mesh;

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
	template<typename F, typename ...Args>
	static typename TimeT::rep execution(F&& func, Args&&... args)
	{
		auto start = std::chrono::steady_clock::now();
		std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
		auto duration = std::chrono::duration_cast<TimeT>
			(std::chrono::steady_clock::now() - start);
		return duration.count();
}
};


#if defined(_WIN32) || defined(_WIN64)
#include <process.h>
#define GET_PID ::_getpid()
#else
#include <unistd.h>
#define GET_PID getpid()
#endif

#define RUN_TEST(x) run_test(x, #x);

void run_test(const std::function<bool()>& test_function, const std::string function_name)
{
	bool test_return = test_function();
	if (test_return)
		std::cout << "Test " << function_name << " passed." << std::endl;
	else
		std::cout << "Test " << function_name << " failed." << std::endl;
}

bool test_node()
{
	using NodeTyped = DistribNode<Node>;

	NodeTyped a({ 1, 0, 0 });
	assert(a.get_global_number() == 0);
	assert(a.get_local_number() == 0);
	assert(a.get_coord() == Point({ 1, 0, 0 }));
	assert(a.get_parent_partition() == 0);

	NodeTyped b(1, 1, { 0, 0, 1 }, 1);
	assert(b.get_global_number() == 1);
	assert(b.get_local_number() == 1);
	assert(b.get_coord() == Point({ 0, 0, 1 }));
	assert(b.get_parent_partition() == 1);

	assert(a > b);
	assert(b < a);
	assert(!(a == b));
	assert(NodeTyped::num_based_compare(a, b));
	assert(NodeTyped::part_based_compare(a, b));

	Node* test1 = &a;
	Node** test2 = new Node*;
	*test2 = &a;

	std::array<NodeTyped *, 2> nodes = { &a, &b };
	Point mean = get_mean_coordinate_of_nodes(nodes);
	assert(distance(mean, { 0.5, 0, 0.5 }) < 1e-16);
	assert(fabs(distance(a, b) - sqrt(2)) < 1e-16);

	return true;
}

bool test_connectivity()
{
	using NodeTyped = DistribNode<Node>;

	NodeTyped na(0, 0, { 0, 0, 0 });
	NodeTyped nb(1, 1, { 0.5, 0, 0 });
	std::vector<NodeTyped*> c_vec;

	{
		Connectivity<NodeTyped, uint16_t> c({ &na, &nb });
		assert(*c.begin() == &na);
		assert(c.size() == 2);
		for (auto n_ptr : c)
			assert(n_ptr->get_parent_partition() == 0);
		assert(c[0]->get_global_number() == 0);
		assert(c[1]->get_global_number() == 1);

		c_vec = c.copy_as_vector();
		assert(c_vec[0] == c[0]);
		assert(&c_vec[0] != c.begin());
	}

	assert(c_vec[0]->get_global_number() == 0);
	assert(c_vec[1]->get_global_number() == 1);

	return true;
}

bool test_edge()
{
	using NodeTyped = DistribNode<Node>;
	using EdgeTyped = Edge<DistribNode<Node>>;

	NodeTyped na(0, 0, { 0, 0, 0 });
	NodeTyped nb(1, 1, { 0.5, 0, 0 });
	EdgeTyped a(EdgeType::E2N, { &na, &nb });
	assert(a.get_type() == EdgeType::E2N);
	for (auto n_ptr : a)
		assert(n_ptr->get_parent_partition() == 0);
	// Technically, this could be an invalid assert, since the hash of the node addresses
	// actually be exactly 0, but the changes are astronomically small.
	assert(a.get_conn_hash_id() > 0);
	assert(fabs(a.get_length() - 0.5) < 1e-16);
	Point mean = a.get_mean_coordinate();
	assert(distance(mean, { 0.25, 0, 0 }) < 1e-16);

	NodeTyped nc(2, 2, { 1, 0, 0 });
	a.swap_nodes(&nb, &nc);
	assert(fabs(a.get_length() - 1) < 1e-16);

	EdgeTyped b(EdgeType::E2N, { &na, &nb });
	assert(b < a);

	return true;
}

bool test_face()
{
	using NodeTyped = DistribNode<Node>;
	using FaceTyped = Face<DistribNode<Node>>;

	NodeTyped na(0, 0, { 0, 0, 0 });
	NodeTyped nb(1, 1, { 1, 0, 0 });
	NodeTyped nc(2, 2, { 1, 1, 0 });
	NodeTyped nd(2, 2, { 0, 1, 0 });

	FaceTyped a(FaceType::Q3D4N, { &na, &nb, &nc, &nd });
	assert(a.get_type() == FaceType::Q3D4N);
	for (auto n_ptr : a)
		assert(n_ptr->get_parent_partition() == 0);
	// Technically, this could be an invalid assert, since the hash of the node addresses
	// actually be exactly 0, but the changes are astronomically small.
	assert(a.get_conn_hash_id() > 0);

	Point mean = a.get_mean_coordinate();
	assert(distance(mean, { 0.5, 0.5, 0 }) < 1e-16);

	auto edges = a.get_edges();
	assert(edges.size() == 4);

	auto outward_norm = a.get_outward_normal();
	assert((outward_norm - Vector3D(0, 0, 1)).norm() < 1e-16);
	double projection = a.calc_signed_norm_projection({ 0.5, 0.5, 1 });
	assert(fabs(projection - 1.0) < 1e-16);

	return true;
}

bool test_element()
{
	using NodeTyped = DistribNode<Node>;
	using ElementTyped = MeshGenElement<NodeTyped, DistribElement<NodeTyped, Element<NodeTyped>>>;

	NodeTyped na(0, 0, { 0, 0, 0 });
	NodeTyped nb(1, 1, { 1, 0, 0 });
	NodeTyped nc(2, 2, { 1, 1, 0 });
	NodeTyped nd(2, 2, { 0, 1, 0 });

	ElementTyped a;
	a.initialize(ElementType::Q2D4N, { &na, &nb, &nc, &nd });
	a.set_component_id(1);
	assert(a.get_component_id() == 1);
	a.set_local_number(1);
	assert(a.get_local_number() == 1);
	a.set_global_number(1);
	assert(a.get_global_number() == 1);
	a.set_parent_partition(1);
	assert(a.get_parent_partition() == 1);

	assert(a.contains_node(na));
	assert(a.get_conn_hash_id() > 0);
	Point mean = a.get_mean_coordinate();
	assert(distance(mean, { 0.5, 0.5, 0 }) < 1e-16);

	assert(a.get_orientations().size() == 0);
	a.allocate_lcs(LCSAssignment::None);
	assert(a.get_orientations().size() == 0);
	a.allocate_lcs(LCSAssignment::Element);
	assert(a.get_orientations().size() == 1);
	a.allocate_lcs(LCSAssignment::Component);
	assert(a.get_orientations().size() == 1);
	a.allocate_lcs(LCSAssignment::Node);
	assert(a.get_orientations().size() == 4);

	ElementVector S;
	ElementMatrix3D dS;
	a.get_shape_function_values({ 0, 0, 0 }, S);
	a.get_shape_function_derivatives({ 0, 0, 0 }, dS);

	return true;
}

std::array<std::array<FloatType, 3>, 4> get_tetrahedron_2nd_order_gauss()
{
	std::array<std::array<FloatType, 3>, 4> gauss_points;
	const double val1 = (5.0 - sqrt(5.0)) / 20;
	const double val2 = (5.0 + 3 * sqrt(5.0)) / 20;

	gauss_points[0][0] = val1;
	gauss_points[0][1] = val1;
	gauss_points[0][2] = val1;
	//weight[orderIndex][0] = 1.0 / 24;

	gauss_points[1][0] = val2;
	gauss_points[1][1] = val1;
	gauss_points[1][2] = val1;
	//weight[orderIndex][1] = 1.0 / 24;

	gauss_points[2][0] = val1;
	gauss_points[2][1] = val2;
	gauss_points[2][2] = val1;
	//weight[orderIndex][2] = 1.0 / 24;

	gauss_points[3][0] = val1;
	gauss_points[3][1] = val1;
	gauss_points[3][2] = val2;
	//weight[orderIndex][3] = 1.0 / 24;

	return gauss_points;
}

bool test_polyhedron()
{
	using NodeTyped = DistribNode<Node>;

	NodeTyped n0(0, 0, { 0, 0, 0 });
	NodeTyped n1(1, 1, { 1, 0, 0 });
	NodeTyped n2(2, 2, { 0, 1, 0 });
	NodeTyped n3(3, 3, { 0, 0, 1 });

	// Basically a 4-node tetrahedral
	PolyhedronElement<NodeTyped> poly;
	poly.initialize(ElementType::POLYHEDRON, {
		{&n0, &n2, &n1},
		{&n0, &n1, &n3},
		{&n2, &n0, &n3},
		{&n1, &n2, &n3} });
	Element<NodeTyped> tet;
	tet.initialize(ElementType::T3D4N, { &n0, &n1, &n2, &n3 });

	std::map<LocalIndex, LocalIndex> poly_to_tet_node_map;
	for (auto i = 0; i < poly.get_conn().size(); i++)
		poly_to_tet_node_map[i] = poly.get_conn()[i]->get_global_number();

	ElementVector S_tet;
	ElementVector S_poly;
	ElementMatrix3D tet_node_coordinates;
	Vector3D cartesian_point;

	auto compare_shape = [&poly_to_tet_node_map](const ElementVector& S_tet, const ElementVector& S_poly) -> bool {
		for (auto i = 0; i < S_tet.size(); i++)
			if (fabs(S_poly[i] - S_tet[poly_to_tet_node_map[i]]) > 1e-12)
				return false;
		return true;
	};

	// Test center of tetrahedral (shape functions should all be equal and 0.25)
	Point barycentric_center = { 0.25, 0.25, 0.25 };
	tet.get_shape_function_values(barycentric_center, S_tet);
	for (auto i = 0; i < S_tet.size(); i++)
		assert(fabs(S_tet[i] - 0.25) < 1e-16);
	tet.get_node_coordinates(tet_node_coordinates);
	cartesian_point = tet_node_coordinates.transpose() * S_tet;
	poly.get_shape_function_values({ cartesian_point[0], cartesian_point[1], cartesian_point[2] }, S_poly);
	assert(compare_shape(S_tet, S_poly));

	auto gauss_points = get_tetrahedron_2nd_order_gauss();
	for (auto& gauss_point : gauss_points) {
		tet.get_shape_function_values(gauss_point, S_tet);
		tet.get_node_coordinates(tet_node_coordinates);
		cartesian_point = tet_node_coordinates.transpose() * S_tet;
		poly.get_shape_function_values({ cartesian_point[0], cartesian_point[1], cartesian_point[2] }, S_poly);
		assert(compare_shape(S_tet, S_poly));
	}

	/*
	// Benchmarks
	const int num_calls = 1000000;
	auto bench_tet = [&tet, &barycentric_center, &S_tet, &num_calls]() {
		for (auto i = 0; i < num_calls; i++)
			tet.get_shape_function_values(barycentric_center, S_tet);
	};
	std::cout << "Time for repeated tetrahedron shape function calls "
		<< measure<>::execution(bench_tet)
		<< " milliseconds." << std::endl;

	auto bench_poly = [&poly, &barycentric_center, &S_poly, &num_calls]() {
		for (auto i = 0; i < num_calls; i++)
			poly.get_shape_function_values(barycentric_center, S_poly);
	};
	std::cout << "Time for repeated polyhedron shape function calls "
		<< measure<>::execution(bench_poly)
		<< " milliseconds." << std::endl;
		*/

	return true;
}

bool test_mesh()
{
	using NodeTyped = DistribNode<Node>;
	using ElementTyped = MeshGenElement<NodeTyped, DistribElement<NodeTyped, Element<NodeTyped>>>;
	using MeshType = Mesh<NodeTyped, ElementTyped>;

	MeshType mesh;
	NodeTyped * na = mesh.add_node_by_coord({ 0, 0, 0 });
	NodeTyped * nb = mesh.add_node_by_coord({ 1, 0, 0 });
	NodeTyped * nc = mesh.add_node_by_coord({ 1, 1, 0 });
	NodeTyped * nd = mesh.add_node_by_coord({ 0, 1, 0 });
	ElementTyped * ea = mesh.add_element();
	ea->initialize(ElementType::Q3D4N, { na, nb, nc, nd });
	mesh.update_global_numbering();
	mesh.sort_by_partition();

	// TODO: add asserts

	return true;
}

int main()
{
	initialize_flex_mesh(1e-14);

	RUN_TEST(test_node);
	RUN_TEST(test_connectivity);
	RUN_TEST(test_edge);
	RUN_TEST(test_face);
	RUN_TEST(test_element);
	RUN_TEST(test_polyhedron);
	RUN_TEST(test_mesh);

	//system("pause");
	return 0;
}