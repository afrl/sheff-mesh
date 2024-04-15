#pragma once

#include "Mesh.h"

#include <variant>

namespace mesh
{
using PartIndex = uint32_t;

class AssemblyObject
{
};
/*
namespace detail
{
struct unused
{
};
}

// This concept tests whether the type is equivalent to detail::unused
template<typename T>
concept is_used = requires(T t) { std::is_same_v<T, detail::unused>; };

// Forward declaration
template<typename PartType, typename SubAssemblyType>
class Assembly;

template<typename PartType = detail::unused, typename SubAssemblyType = detail::unused>
class Assembly
{
  public:
    using part_type        = PartType;
    using subassembly_type = SubAssemblyType;

    PartType& add(PartType&& part)
        requires is_used<part_type>
    {
        child_objects_.push_back(std::move(part));
        return std::get<PartType>(child_objects_.back());
    }

    void add(SubAssemblyType&& subassembly)
        requires(is_used<subassembly_type> && is_used<typename subassembly_type::part_type>)
    {
        child_objects_.push_back(std::move(subassembly));
    }

    std::vector<std::variant<SubAssemblyType, PartType>> child_objects_;
};
*/

struct asssembly_tag
{
};

template<typename PartType>
class Assembly
{
  public:
    using part_type   = PartType;
    using is_assembly = asssembly_tag;

    PartType& add(PartType&& part)
    {
        child_objects_.push_back(std::move(part));
        return child_objects_.back();
    }

    PartType* begin()
    {
        return child_objects_.begin();
    }
    PartType* end()
    {
        return child_objects_.end();
    }
    const PartType* begin() const
    {
        return child_objects_.data();
    }
    const PartType* end() const
    {
        return child_objects_.data() + child_objects_.size();
    }
    const PartType* cbegin() const
    {
        return child_objects_.begin();
    }
    const PartType* cend() const
    {
        return child_objects_.end();
    }

    std::vector<PartType> child_objects_;
};


/// @brief Holds a collection of "parts" or other "subassemblies"
///
/// Note: a "part" will refer to an independently meshed piece of an "assembly" or "subassembly".
///
/// @tparam PartType - a type that must have a method `Mesh& PartType::mesh()` and
///         `const Mesh& PartType::mesh() const`
template<typename PartType>
class SubAssembly : AssemblyObject
{
  public:
    SubAssembly();
    SubAssembly(const SubAssembly& other);
    SubAssembly(SubAssembly&& other) noexcept;

    bool operator==(const SubAssembly& other) const;

    //////////////////////////////////////////////////////
    // Accessors/iterators

    PartType& operator[](const PartIndex& i);
    const PartType& operator[](const PartIndex& i) const;
    PartIndex size() const;

    PartType* begin()
    {
        return parts_.begin();
    }
    PartType* end()
    {
        return parts_.end();
    }
    const PartType* begin() const;
    const PartType* end() const;
    const PartType* cbegin() const;
    const PartType* cend() const;

    //////////////////////////////////////////////////////
    // Mutators

    /// @brief Adds a part with an empty mesh to the end of the assembly and returns a reference.
    PartType& add_part();


  private:
    std::vector<PartType> parts_;
};
}