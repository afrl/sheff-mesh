#include "mesh/Assembly.h"

namespace mesh
{
using DigitalChain = mesh::Mesh;
using Tow          = Assembly<DigitalChain>;
using Tows         = Assembly<Tow>;

using Mold  = mesh::Mesh;
using Molds = Assembly<Mold>;

//using VTMSGeometry = Assembly<std::variant<Tows, Molds>>;

// Then make a PartIterator?
//  - Need a check to see if the PartType == unused
//  - Use std::get_if
//  - Serialization?

template<typename AssemblyType>
class PartIterator
{
  public:
    PartIterator(const AssemblyType& assembly)
      : assembly_(assembly)
    {
    }

    std::vector<const mesh::Mesh*> meshes()
    {
        return get_meshes(assembly_);
    }

    template<typename T>
    std::vector<const mesh::Mesh*> get_meshes(const T& t)
        requires(std::is_same_v<typename T::part_type::is_assembly, asssembly_tag>)
    {
        // Each object in this assembly is another assembly, so loop over those
        std::vector<const mesh::Mesh*> mesh_ptrs;
        for (auto& subassembly : t.child_objects_)
        {
            auto submesh_ptrs = get_meshes(subassembly);
            mesh_ptrs.insert(mesh_ptrs.end(), submesh_ptrs.begin(), submesh_ptrs.end());
        }
        return mesh_ptrs;
    }
    
    template<typename T>
    std::vector<const mesh::Mesh*> get_meshes(const T& t)
        requires(std::is_same_v<typename T::part_type, mesh::Mesh>)
    {
        std::vector<const mesh::Mesh*> mesh_ptrs;
        // Each object in this assembly is a part
        for (auto& part : t.child_objects_)
        {
            mesh_ptrs.push_back(&part);
        }
        return mesh_ptrs;
    }
    
    const AssemblyType& assembly_;
};
}