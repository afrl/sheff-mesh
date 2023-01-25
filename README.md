## Mesh {#page_mesh_module}

### Overview

```
Important Note:

Within the engineering finite element analysis community, a mesh is traditionally composed of
nodes and elements.  In the applied mathematics community, a mesh is traditionally composed
of points and cells.  For this framework, we adopt the engineering semanics.
```

The primary data structure for the module is mesh::Mesh, which is a very simple structure.  The data structures that store information about the nodes (mesh::NodeContainer) and elements (mesh::ElementContainer) are also quite simple.  Any complexity lies in the interfaces to access information about a node/element, iterators, and the algorithms that operate on a mesh.  Both mesh::NodeContainer and mesh::ElementContainer take the philosophy of storing each piece of data for a node or element in an array that includes the entries for all nodes/elements in the container.  For example, the global number of the nodes are stored as contiguous array, while the coordinates of the nodes are stored in a second contiguous array.  This is similar to how VTK chooses to store data for unstructured grids.  However, since the data is not grouped for each individual node and element, a convenient interface to access information about the nodes and elements is needed.  The most important classes for accessing information are:
* mesh::ElementAccessor - provides an interface to get information about a single element.  It does not store the data, but rather internally holds a reference to a mesh::ElementContainer.
* mesh::ElementConstAccessor - const equivalent of mesh::ElementAccessor that does not allow mutable access to data in the underlying mesh::ElementContainer.
* mesh::ElementIterator - an iterator that sequentially moves through the elements in a mesh::ElementContainer instance. Dereferencing the iterator yields a mesh::ElementAccessor object.
* mesh::ElementConstIterator - const equivalent of mesh::ElementIterator.
* mesh::NodeAccessor - provides an interface to get information about a single node.  Like mesh::ElementAccessor, it internally holds a reference to a mesh::NodeContainer instead of storing data.
* mesh::NodeConstIterator - const equivalent of mesh::NodeAccessor.
* mesh::NodeIterator - an iterator that sequentially moves through the nodes in a mesh::NodeContainer instance. Dereferencing the iterator yields a mesh::NodeAccessor object.
* mesh::NodeConstIterator - const equivalent of mesh::NodeIterator.
* mesh::NodeSubrangeIterator - an iterator that moves over a particular subset of the nodes in a mesh::NodeContainer instance. This is very handy for iterating over the nodes that appear in the connectivity of an element for example.  Dereferencing the iterator yields a mesh::NodeAccessor object.
* mesh::NodeSubrangeConstIterator - const equivalent of mesh::NodeSubrangeIterator.

### Examples

When using the mesh::Mesh, one typically will either
1. start with an empty mesh and add nodes and elements to it or
2. use a built-in algorithm to create a mesh::Mesh instance.

Examples of both of these cases and how data can be accessed in a mesh are shown in the following sections.

#### Creating a new mesh
```cpp
// Relevant includes:
###include "mesh/Mesh.h"
// End includes

mesh::Mesh m;

// Adds nodes given a coordinate.  Each new node takes the next highest local and global number.
m.nodes().add({ 0, 0, 0 });
m.nodes().add({ 1, 0, 0 });
m.nodes().add({ 1, 1, 0 });
m.nodes().add({ 0, 1, 0 });

// Adds a 2D 4-node quadrilateral element with the second parameter indicating the connectivity.
m.elements().add(mesh::ElementType::Q2D4N, { 0, 1, 2, 3 });
```

#### Creating a mesh via a built-in algorithm
```cpp
// Relevant includes:
###include "mesh/Mesh.h"
###include "mesh/algorithms/Creation.h"
###include "mesh/algorithms/IO.h"
// End includes

// Creates a mesh from an Abaqus input script
// Note: There are quite a few features of Abaqus scripts that may not be supported yet.
//       This will be an ongoing process to add support for features as needed.
mesh::Mesh mesh_from_abaqus = mesh::read_abaqus("BarWithHole.inp");

// Creates a 2x4 2D uniform grid of 4-node quadrilateral elements
// The arguments specify the corners of the elements, so each is 1x1 in size.
mesh::Mesh grid = mesh::make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });
```

#### Accessing data in a mesh
```cpp
// Relevant includes:
###include "mesh/Mesh.h"
// End includes

mesh::Mesh m;
// Assume something has built the mesh ...

// Note: node will have the type mesh::NodeAccessor or mesh::NodeConstAccessor
// depending on whether the mesh is const or not.
for (auto node : m.nodes())
{
    // See mesh::NodeAccessorBase for the full list of methods.
    // Note: mesh::NodeAccessor and mesh::NodeConstAccessor are just non-const
    //       const specialization of mesh::NodeAccessorBase.

    // Example of interface:
    node.local_index(); // returns a mesh::LocalIndex
    node.global_index(); // returns a mesh::GlobalIndex
    node.x(); // returns double
    node.y(); // returns double
    node.z(); // returns double
    node.coord(); // returns wrapper for Eigen::Vector3d
    node.parent_rank(); // returns a RankIndex
}

// Note: element will have the type mesh::ElementAccessor or mesh::ElementConstAccessor
// depending on whether the mesh is const or not.
for (auto element : m.elements())
{
    // See mesh::ElementAccessorBase for the full list of methods.
    // Note: mesh::ElementAccessor and mesh::ElementConstAccessor are just non-const
    //       const specialization of mesh::ElementAccessorBase.

    // Example of interface:
    element.local_index(); // returns a mesh::LocalIndex
    element.global_index(); // returns a mesh::GlobalIndex
    element.type(); // returns a mesh::ElementType
    element.connectivity(); // returns a slice of the local indices of the nodes in the element
    // The connectivity can be used in a ranged-based loop:
    for (auto& node_index : element.connectivity()) {}
    // The nodes can be accessed with a ranged-based loop:
    // Note: "node" is either mesh::NodeAccessor or mesh::NodeConstAccessor depending
    // on whether the element is const or not.
    for (auto& node : element.nodes()) {}
    element.node(0); // returns an accessor for the first node in the connectivity
    element.parent_rank(); // returns a RankIndex
    element.local_coord_sys(); // returns a mesh::QuaternionContainer
    // Note: the local coordinate systems is the nominal local coordinate system for the element,
    //       but if an element is enriched (see the finite_element module) then each field could
    //       either refer to this local coordinate system or have their own separate one.
}
```

#### Partitioning the mesh
To perform a finite element analysis on a very large mesh, a high-performance computing system is typically needed, which is composed of multiple computation nodes.  The distributed nature of the hardware requires a distributed analysis.  To accommodate this, it is useful to partition the mesh.  Some terminology related to HPC analyses and mesh partitioning are:
* Rank: an independent process running on a computation node (there can be multiple per computation node).  This is a concept integral to MPI programs.
* Partition: a piece of the mesh that has nothing to do with the model.  It only is meant to accommodate the distributed software.  In this framework, we want one partition per a rank for the purpsoe of analyses.  For post-processing/visualization, multiple partitions can be handled by a single rank (due to the lower computational cost compare to the analysis).
* Parent partition: the partition that "owns" an element or node.  For elements, the parent partition is important to the analysis, since each element should only be assembled into the global matrix once.  For nodes, the parent partition is only relevant for the initial DoF enumeration and post-processing.  For nodes and elements, each can only have one parent partition.
* Ghost elements: elements in the mesh for a partition that are not owned by the partition.  For an elastic analysis, ghost elements are not actually needed, but XFEM analyses and other analyses requiring non-local information require ghost elements.

To partition a mesh with METIS:
```cpp
// Relevant includes:
###include "mesh/algorithms/Creation.h"
###include "mesh/algorithms/IO.h"
###include "mesh/algorithms/Partitioning.h"
// End includes

std::vector<double> coordinates(5);
std::iota(coordinates.begin(), coordinates.end(), 0.0);
mesh::Mesh mesh = mesh::make_linear_grid(coordinates, coordinates, coordinates);

// Partitions the mesh into 3 partitions in such a way to approximately minimize the communcation
// needed during the FEA analysis.
mesh::partition_metis(mesh, 3);

// Extracts a mesh for partition 0 that includes all of the owned elements/nodes and a layer of
// ghost elements.
auto mesh_partition_0 = extract_partition_mesh(mesh, 0);
```

#### Reading/writing the mesh

##### Write to VTK
```cpp
// Relevant includes:
###include "mesh/Mesh.h"
###include "mesh/algorithms/IO.h"
// End includes

mesh::Mesh m;
// Assume something has built the mesh ...
mesh::write_vtk(m, "mesh.vtk");
```

##### Reading/writing to fea_io format
```cpp
// Relevant includes:
###include "mesh/Mesh.h"
###include "mesh/algorithms/IO.h"
###include "mesh/algorithms/Partitioning.h"
// End includes

mesh::Mesh m;
// Assume something has built the mesh ...

// Writes the mesh to mesh_no_meta.m.json/mesh_no_meta.m.bin files without creating a meta file
// The second paramter is the base directory.
mesh::write_fea_io(m, ".", "mesh_no_meta");

// Reads partition 0 directly from the mesh_no_meta.m.json/mesh_no_meta.m.bin file and fills a mesh.
auto m_no_meta = mesh::read_fea_io("./mesh_no_meta", 0);

// Writes the mesh to mesh.m.json/mesh.m.bin files and creates a meta file (mesh.x.json).
// The last parameter is the meta file's name.
// Note: also creates a mesh.g.json/mesh.g.bin files if groups exist within the mesh.
mesh::write_fea_io(m, ".", "mesh", "mesh");

// Reads partition 0 of the mesh correspding to timestep 0 directly from the
// mesh_no_meta.m.json/mesh_no_meta.m.bin file and fills a mesh.
auto m_with_meta = mesh::read_fea_io_from_meta("./mesh.x.json", 0, 0);

// If m is parititioned and distributed across multiple ranks, then you can use
mesh::write_fea_io_distributed(/*args*/);
// to write the mesh in parallel.
```

##### Reading PFEC format
The PFEC (partitioned finite element compressed) file format was developed by Keith Ballard during his doctoral work.  The format has been replaced by serialized and fea_io formats (each accomplish a different goal), but several tools still exist for creating fiber/matrix, multiscale laminated composite, and 3D textile composite models that only write to the PFEC format.  You can read those files in using:
```cpp
// Relevant includes:
###include "mesh/Mesh.h"
###include "mesh/algorithms/IO.h"
// End includes

// Reads partition 0 from the mesh for the first time step of the meta file "example.pfecx"
// Note: there is not the option to read data or meshes from other time steps.
// Note: "domains" get converted into groups with the naming conventions "Domain0", "Domain1", etc.
// Note: orientation data does populate the local coordinate systems in the resulting mesh.
auto m = mesh::read_pfec("example.pfecx", 0);
```

##### Reading Abaqus .inp files
Abaqus input scripts have quite a large set of commands that can be used within them.  Currently, we have a reader that can read simple Abaqus input scripts, but we intend to extend this functionality as needed going forward.  So feel free to contribute or suggest improvements to this algorithm especially.

```cpp
// Relevant includes:
###include "mesh/Mesh.h"
###include "mesh/algorithms/IO.h"
// End includes

auto m = mesh::read_abaqus("example.inp");
```

##### Serialization and deserailization
The utility module provides serialization functionality (refer to documentation in that module for more info), but the mesh module supports serialization for all classes that store data (and a few that do not for convenience to write information to JSON for debugging).  Serialization is a great way to:
* Store the mesh in a compact binary format (though reading the file without the mesh module in a 3rd party code would be difficult)
* Tranfer the mesh as data over a network (deserializing it on the client side)
* Write the mesh or parts of the mesh to JSON for debugging

```cpp
// Relevant includes:
###include "mesh/Mesh.h"
###include "mesh/algorithms/Creation.h"
###include "utility/serialization/Serialization.h"
###include <iostream>
###include <sstream>
// End includes

// Create some mesh for an example
mesh::Mesh linear_grid = mesh::make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

// Write the mesh in a JSON format to cout
{
    utility::JSONSerializer to_json_cout(std::cout);
    to_json_cout(linear_grid);
    // Writes to output stream on destruction, since JSON is a format that cannot
    // be written line by line.
}

// Count how many bytes would be needed to write to binary
utility::SerialByteCounter counter;
counter(linear_grid);
counter.n_bytes(); // returns the # of bytes the binary representation would require

// Write in binary format to a string stream
std::stringstream ss;
utility::OutputSerializer to_binary_ss(ss);
to_binary_ss(linear_grid);
// ss contains the binary representation now

// Deserialize a mesh from the binary representation
utility::InputDeserializer input_from_ss(ss);
Mesh linear_grid_deserialized;
input_from_ss(linear_grid_deserialized);
// Note: linear_grid_deserialized will be the same as linear_grid
```

### Third-Party Mesh Libraries

* @subpage nanoflann A fast K-D tree data structure