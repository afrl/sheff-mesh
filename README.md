# SHEFF Mesh (`sheff-mesh`)

`sheff-mesh` provides unstructured 1D/2D/3D finite element mesh data structures, isoparametric basis functions, spatial search trees, METIS domain partitioning, mesh generation/refinement algorithms, and multi-format I/O for the **Scalable Hierarchical Enrichment Finite-element Framework (SHEFF)**.

## Capabilities

- **Cache-Friendly Mesh Data Structures (`mesh/Mesh.h`, `mesh/node/`, `mesh/element/`, `mesh/face/`, `mesh/surface/`)**:
  - Structure-of-arrays storage in `mesh::NodeContainer` and `mesh::ElementContainer` paired with lightweight proxy accessors and iterators (`NodeAccessor`, `ElementAccessor`, `FaceAccessor`, `SurfaceElementAccessor`, `NodeSubrangeIterator`, `ElementSubrangeIterator`).
  - Topological relations (`ElementAdjacencyList`, `NodeToElementMap`, `ElementToComponentMap`), named `NodeGroup` / `ElementGroup` sets, per-element quaternion local coordinate systems (`QuaternionContainer`), and typed field attachments (`MeshAssociatedData`, `DataFields`).
- **Isoparametric Basis Functions (`mesh/basis/`)**:
  - Linear and quadratic shape functions, parametric derivatives, and isoparametric mappings (`BasisFunctions`, `ElementBasis`, `BasisFunctionsCache`, `IsoparametricMappingUtils`) across lines, triangles, quadrilaterals, tetrahedra, wedges, pyramids, hexahedra, shells, and cohesive elements.
- **Mesh Generation, Manipulation & Partitioning (`mesh/algorithms/`)**:
  - Structured grid creation (`make_linear_grid`), extrusion (`Extrusion.h`), uniform/element refinement and linear-to-quadratic conversion (`Refinement.h`), fiber/matrix microstructure meshing (`FiberMatrix.h`), component connection pairing (`Connections.h`), component interface extraction (`Search.h`), TetGen/Triangle helpers (`TetgenHelpers.h`, `Polygons.h`), K-D tree queries (`kd_tree/KDTree.h`), and METIS partitioning with ghost-layer extraction (`Partitioning.h`).
- **Multi-Format Mesh I/O & Serialization (`mesh/algorithms/IO.h`, `mesh/external_interfaces/`)**:
  - Native parallel/serial `fea_io` (`.m.json`/`.m.bin`/`.x.json`), Abaqus (`.inp`), VTK (`.vtk`), legacy BSAM classic, PFEC (`.pfecx`), and binary/JSON `utility` serialization.

## Project Structure

- `include/mesh/`: Public headers covering `node/`, `element/`, `face/`, `surface/`, `basis/`, `relations/`, `data/`, `kd_tree/`, `algorithms/`, and `external_interfaces/`.
- `source/`: Compiled mesh containers, basis functions, partitioning/refinement algorithms, and format readers/writers.
- `unit_tests/`: Unit tests for mesh construction, basis evaluation, METIS partitioning, search, and file I/O.
- **CMake Package / Target**: `Mesh` / `Mesh::mesh` (static library).

## Dependencies

- **Build System**: [`sheff-build-tools`](https://github.com/afrl/sheff-build-tools) (`SHEFF_TOOLS`), CMake 3.17+, C++17/C++20 compiler.
- **First-Party SHEFF Modules**: [`sheff-fea-io`](https://github.com/afrl/sheff-fea-io) (`FeaIo::fea_io`), [`sheff-math`](https://github.com/afrl/sheff-math) (`SimpleMath::simple_math`), [`sheff-utility`](https://github.com/afrl/sheff-utility) (`Utility::utility`).
- **Second-Party SHEFF Modules**: [`sheff-second-party-mesh`](https://github.com/afrl/sheff-second-party-mesh) (`SecondPartyMesh::second_party_mesh`), [`sheff-second-party-metis`](https://github.com/afrl/sheff-second-party-metis) (`SecondPartyMetis::second_party_metis`).
- **Third-Party Libraries**: [`nlohmann/json`](https://github.com/nlohmann/json), [`RapidJSON`](https://github.com/Tencent/rapidjson).

## Building

When building as part of [SPAWC](https://github.com/afrl/spawc), `sheff-mesh` is configured and built automatically via `make.py`.

To build and test standalone (with `SHEFF_TOOLS` configured):

```bash
python3 make.py --config --build --third-party
python3 make.py --config --build --make-tests --run-tests
```

## Usage Examples

> **Terminology Note**: Following engineering finite element conventions, a `mesh::Mesh` is composed of *nodes* (`mesh::NodeContainer`) and *elements* (`mesh::ElementContainer`).

### Creating and Traversing a Mesh

```cpp
#include "mesh/Mesh.h"
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"

mesh::Mesh m;
m.nodes().add({ 0.0, 0.0, 0.0 });
m.nodes().add({ 1.0, 0.0, 0.0 });
m.nodes().add({ 1.0, 1.0, 0.0 });
m.nodes().add({ 0.0, 1.0, 0.0 });
m.elements().add(mesh::ElementType::Q2D4N, { 0, 1, 2, 3 });

// Or generate structured / imported meshes:
mesh::Mesh grid = mesh::make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });
mesh::Mesh imported = mesh::read_abaqus("BarWithHole.inp");

for (auto element : grid.elements()) {
    auto type = element.type();
    for (auto node : element.nodes()) {
        auto coord = node.coord();
    }
}
```

### Partitioning & Writing Distributed Meshes

```cpp
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Partitioning.h"

mesh::Mesh m = mesh::make_linear_grid({ 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4 });
mesh::partition_metis(m, 3);
auto part0 = mesh::extract_partition_mesh(m, 0);

mesh::write_vtk(part0, "part0.vtk");
mesh::write_fea_io(m, ".", "mesh", "mesh");
```

## Public Release Information

Distribution Statement A. Approved for public release: distribution is unlimited. Case #: AFRL-2026-4513
