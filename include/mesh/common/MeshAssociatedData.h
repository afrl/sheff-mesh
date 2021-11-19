#pragma once

#include "simple_math/EigenDefines.h"
#include "utility/StringUtils.h"

#include <string>
#include <vector>

namespace mesh
{
struct DataField
{
    std::string name;
    std::vector<std::string> component_names;
    simple_math::GeneralMatrix data;
};

inline DataField to_data_field(const std::string& name, const std::vector<double>& data)
{
    return { name, { "" }, Eigen::Map<const simple_math::GeneralVector>(data.data(), data.size(), 1) };
}

/// @brief A struct that holds data associated with a mesh, such as data at nodes, elements,
/// quadrature points, etc.
struct MeshAssociatedData
{
    /// Holds data fields with a 1-to-1 relationship with nodes.
    std::vector<DataField> node_data;
    /// Holds data fields with a 1-to-1 relationship with elements.
    std::vector<DataField> element_data;

    /// Holds data fields for all nodes within each element.  Note: node will typically have
    /// multiple values since they often appear in multiple elements.
    std::vector<DataField> node_by_element_data;
    /// Element offsets for node_by_element_data.
    std::vector<std::size_t> node_by_element_offsets;

    /// Holds data fields for all nodes within each component.  Note: a node may have multiple
    /// values if it appears in multiple components.
    std::vector<DataField> node_by_component_data;
    /// Component offsets for node_by_component_data.
    std::vector<std::size_t> node_by_component_offsets;

    /// Holds data fields with a 1-to-1 relationship with quadrature points within elements.
    /// Note: quad_point_offsets indicates which quadrature points belong to each element.
    std::vector<DataField> quad_point_data;
    /// Offsets for quadrature points in each element.  For example, the quad points belonging to
    /// element i start at index quad_point_offsets[i] and stops at index quad_point_offsets[i + 1]
    /// (exclusive)
    std::vector<std::size_t> quad_point_offsets;

    //////////////////////////////////////////////////////////////////////////////////////
    // Accessor Methods

    /// @brief Returns true if all data is empty.
    inline bool empty() const
    {
        return node_data.empty() && element_data.empty() && node_by_element_data.empty() &&
               node_by_component_data.empty() && quad_point_data.empty();
    }

    /// @brief Finds a node data field by its name and returns a const reference, but
    /// throws if it doesn't exist.
    inline const DataField& get_node_data(const std::string& field_name) const
    {
        for (auto& d : node_data)
        {
            if (utility::istr_is_equal(d.name, field_name))
            {
                return d;
            }
        }
        throw _UTIL_MSG_EXCEPT_C("Node data field '" + field_name + "' does not exist.");
    }

    /// @brief Finds a element data field by its name and returns a const reference, but
    /// throws if it doesn't exist.
    inline const DataField& get_element_data(const std::string& field_name) const
    {
        for (auto& d : element_data)
        {
            if (utility::istr_is_equal(d.name, field_name))
            {
                return d;
            }
        }
        throw _UTIL_MSG_EXCEPT_C("Element data field '" + field_name + "' does not exist.");
    }

    /// @brief Finds a node by element data field by its name and returns a const reference, but
    /// throws if it doesn't exist.
    inline const DataField& get_node_by_element_data(const std::string& field_name) const
    {
        for (auto& d : node_by_element_data)
        {
            if (utility::istr_is_equal(d.name, field_name))
            {
                return d;
            }
        }
        throw _UTIL_MSG_EXCEPT_C("Node by element data field '" + field_name + "' does not exist.");
    }

    /// @brief Finds a quadrature point data field by its name and returns a const reference, but
    /// throws if it doesn't exist.
    inline const DataField& get_quad_point_data(const std::string& field_name) const
    {
        for (auto& d : quad_point_data)
        {
            if (utility::istr_is_equal(d.name, field_name))
            {
                return d;
            }
        }
        throw _UTIL_MSG_EXCEPT_C("Quad point data field '" + field_name + "' does not exist.");
    }
};
}