#pragma once

#include "DataFields.h"
#include "utility/Slice.h"
#include "utility/StringUtils.h"
#include "utility/Traits.h"
#include "utility/serialization/NamedValue.h"

#include <assert.h>
#include <string>
#include <vector>

#define ADD_DATA_HELPERS_DECL(TYPE_NAME, TYPE, LEVEL)                                             \
    /** @brief Helper method to add TYPE data at LEVEL level */                                   \
    /** */                                                                                        \
    /** Note: same as calling the templated add and this only refers to the data (does not */     \
    /**       make a copy) so make sure it stays in scope */                                      \
    void add_##LEVEL##_##TYPE_NAME(const std::string& name,                                       \
                                   const std::vector<std::string>& component_names,               \
                                     const simple_math::DynamicMatrix<TYPE>& data);               \
                                                                                                  \
    /** @brief Helper method to add TYPE data at LEVEL level */                                   \
    /** */                                                                                        \
    /** Note: same as calling the templated add and this only refers to the data (does not */     \
    /**       make a copy) so make sure it stays in scope */                                      \
    void add_##LEVEL##_##TYPE_NAME(const std::string& name, const std::vector<TYPE>& data);       \
                                                                                                  \
    /** @brief Helper method to add TYPE data at LEVEL level */                                   \
    /** */                                                                                        \
    /** Note: data is initialized to zero's and will own the memory internally */                 \
    void add_##LEVEL##_##TYPE_NAME(const std::string& name, const int n_rows);                    \
                                                                                                  \
    /** @brief Helper method to add TYPE data at LEVEL level */                                   \
    /** */                                                                                        \
    /** Note: data is initialized to zero's and will own the memory internally */                 \
    void add_##LEVEL##_##TYPE_NAME(const std::string& name,                                       \
                                   const std::vector<std::string>& component_names,               \
                                   const int n_rows,                                              \
                                    const int n_cols);                                            \
                                                                                                  \
    /** @brief Helper method to add TYPE data at LEVEL level */                                   \
    /** */                                                                                        \
    /** Note: moves the data to be owned by this object */                                        \
    void add_##LEVEL##_##TYPE_NAME(const std::string& name,                                       \
                                   const std::vector<std::string>& component_names,               \
                                   simple_math::DynamicMatrix<TYPE>&& data);                      \

namespace mesh
{
/// @brief Helper class for MeshAssociatedData. Holds data for a particular data type (uint32_t,
/// double, etc.).
template<typename T>
struct MeshAssociatedDataType
{
    /// Holds data fields with a 1-to-1 relationship with nodes.
    std::vector<DataField<T>> node_data_fields;

    /// Holds data fields with a 1-to-1 relationship with elements.
    std::vector<DataField<T>> element_data_fields;

    /// Holds data fields for all nodes within each element.  Note: node will typically have
    /// multiple values since they often appear in multiple elements.
    std::vector<DataField<T>> node_by_element_data_fields;

    /// Holds data fields for all nodes within each component.  Note: a node may have multiple
    /// values if it appears in multiple components.
    std::vector<DataField<T>> node_by_component_data_fields;

    /// Holds data fields with a 1-to-1 relationship with quadrature points within elements.
    /// Note: quad_point_offsets indicates which quadrature points belong to each element.
    std::vector<DataField<T>> quad_point_data_fields;

    /// Holds data fields with a 1-to-1 relationship with nodes in a node group.
    std::vector<GroupDataField<T>> node_group_data_fields;

    /// Holds data fields with a 1-to-1 relationship with elements in a element group.
    std::vector<GroupDataField<T>> element_group_data_fields;


    MeshAssociatedDataType()                              = default;
    MeshAssociatedDataType(MeshAssociatedDataType&&)      = default;
    MeshAssociatedDataType(const MeshAssociatedDataType&) = delete;

    //////////////////////////////////////////////////////////////////////////////////////
    // Query Methods

    /// @brief Returns true if all data is empty.
    inline bool empty() const;

    /// @brief Tests if this instance the exact same data as another (does not test if they are the
    /// same instance)
    inline bool operator==(const MeshAssociatedDataType<T>& other) const;

    //////////////////////////////////////////////////////////////////////////////////////
    // Serialization Methods

    /// @brief Provides serialization support
    template<typename SerializerType>
    void serialize(SerializerType& s) const;

    /// @brief Provides deserialization support
    template<typename SerializerType>
    void deserialize(SerializerType& s);
};

/// @brief A struct that holds data associated with a mesh, such as data at nodes, elements,
/// quadrature points, etc.
struct MeshAssociatedData
{
    /// Element offsets for node_by_element_data.
    std::vector<std::size_t> node_by_element_offsets;

    /// Component offsets for node_by_component_data.
    std::vector<std::size_t> node_by_component_offsets;

    /// Offsets for quadrature points in each element.  For example, the quad points belonging to
    /// element i start at index quad_point_offsets[i] and stops at index quad_point_offsets[i + 1]
    /// (exclusive)
    std::vector<std::size_t> quad_point_offsets;

    /// Integral data for the mesh
    MeshAssociatedDataType<int16_t> int16_data;
    MeshAssociatedDataType<uint16_t> uint16_data;
    MeshAssociatedDataType<int32_t> int32_data;
    MeshAssociatedDataType<uint32_t> uint32_data;
    MeshAssociatedDataType<int64_t> int64_data;
    MeshAssociatedDataType<uint64_t> uint64_data;

    /// Floating-point data for the mesh
    MeshAssociatedDataType<float> float_data;
    MeshAssociatedDataType<double> double_data;


    //////////////////////////////////////////////////////////////////////////////////////
    // Mutator Methods

    /// @brief Adds a new data field collection of data for a mesh.
    /// @tparam DFT - enumeration for the type of data field to add (node, element, etc.)
    /// @tparam T - the type of the underlying data (supports: double and uint32_t)
    /// @param data - name and data for the new field
    ///
    /// This is a templated kernel to support each type of data that can be added.
    template<DataFieldType DFT, typename T>
    void add(DataField<T>&& data);

    ADD_DATA_HELPERS_DECL(int16, int16_t, node)
    ADD_DATA_HELPERS_DECL(uint16, uint16_t, node)
    ADD_DATA_HELPERS_DECL(int32, int32_t, node)
    ADD_DATA_HELPERS_DECL(uint32, uint32_t, node)
    ADD_DATA_HELPERS_DECL(int64, int64_t, node)
    ADD_DATA_HELPERS_DECL(uint64, uint64_t, node)
    ADD_DATA_HELPERS_DECL(float, float, node)
    ADD_DATA_HELPERS_DECL(double, double, node)

    ADD_DATA_HELPERS_DECL(int16, int16_t, element)
    ADD_DATA_HELPERS_DECL(uint16, uint16_t, element)
    ADD_DATA_HELPERS_DECL(int32, int32_t, element)
    ADD_DATA_HELPERS_DECL(uint32, uint32_t, element)
    ADD_DATA_HELPERS_DECL(int64, int64_t, element)
    ADD_DATA_HELPERS_DECL(uint64, uint64_t, element)
    ADD_DATA_HELPERS_DECL(float, float, element)
    ADD_DATA_HELPERS_DECL(double, double, element)

    ADD_DATA_HELPERS_DECL(int16, int16_t, node_by_element)
    ADD_DATA_HELPERS_DECL(uint16, uint16_t, node_by_element)
    ADD_DATA_HELPERS_DECL(int32, int32_t, node_by_element)
    ADD_DATA_HELPERS_DECL(uint32, uint32_t, node_by_element)
    ADD_DATA_HELPERS_DECL(int64, int64_t, node_by_element)
    ADD_DATA_HELPERS_DECL(uint64, uint64_t, node_by_element)
    ADD_DATA_HELPERS_DECL(float, float, node_by_element)
    ADD_DATA_HELPERS_DECL(double, double, node_by_element)

    ADD_DATA_HELPERS_DECL(int16, int16_t, node_by_component)
    ADD_DATA_HELPERS_DECL(uint16, uint16_t, node_by_component)
    ADD_DATA_HELPERS_DECL(int32, int32_t, node_by_component)
    ADD_DATA_HELPERS_DECL(uint32, uint32_t, node_by_component)
    ADD_DATA_HELPERS_DECL(int64, int64_t, node_by_component)
    ADD_DATA_HELPERS_DECL(uint64, uint64_t, node_by_component)
    ADD_DATA_HELPERS_DECL(float, float, node_by_component)
    ADD_DATA_HELPERS_DECL(double, double, node_by_component)

    ADD_DATA_HELPERS_DECL(int16, int16_t, quad_point)
    ADD_DATA_HELPERS_DECL(uint16, uint16_t, quad_point)
    ADD_DATA_HELPERS_DECL(int32, int32_t, quad_point)
    ADD_DATA_HELPERS_DECL(uint32, uint32_t, quad_point)
    ADD_DATA_HELPERS_DECL(int64, int64_t, quad_point)
    ADD_DATA_HELPERS_DECL(uint64, uint64_t, quad_point)
    ADD_DATA_HELPERS_DECL(float, float, quad_point)
    ADD_DATA_HELPERS_DECL(double, double, quad_point)

    //////////////////////////////////////////////////////////////////////////////////////
    // Query Methods

    /// @brief Returns true if all data is empty.
    inline bool empty() const;

    /// @brief Tests if this instance the exact same data as another (does not test if they are the
    /// same instance)
    bool operator==(const MeshAssociatedData& other) const = default;

    //////////////////////////////////////////////////////////////////////////////////////
    // Accessor Methods

    /// @brief Access all data of the templated type as const.
    /// @tparam DFT - enumeration for the type of data field to add (node, element, etc.)
    /// @tparam T - the type of the underlying data (supports: double and uint32_t)
    template<DataFieldType DFT, typename T>
    const std::vector<DataField<T>>& get() const;

    /// @brief Access all data of the templated type as non-const.
    /// @tparam DFT - enumeration for the type of data field to add (node, element, etc.)
    /// @tparam T - the type of the underlying data (supports: double and uint32_t)
    template<DataFieldType DFT, typename T>
    std::vector<DataField<T>>& get()
    {
        return const_cast<std::vector<DataField<T>>&>(std::as_const(*this).get<DFT, T>());
    }

    /// @brief Access data field by name as const for the templated type.
    /// @tparam DFT - enumeration for the type of data field to add (node, element, etc.)
    /// @tparam T - the type of the underlying data (supports: double and uint32_t)
    template<DataFieldType DFT, typename T>
    const DataField<T>& get(const std::string& field_name) const;

    /// @brief Access data field by name as non-const for the templated type.
    /// @tparam DFT - enumeration for the type of data field to add (node, element, etc.)
    /// @tparam T - the type of the underlying data (supports: double and uint32_t)
    template<DataFieldType DFT, typename T>
    auto& get(const std::string& field_name)
    {
        for (auto& d : get<DFT, T>())
        {
            if (utility::istr_is_equal(d.name, field_name))
            {
                return d;
            }
        }
        throw _UTIL_MSG_EXCEPT_C("Data field '" + field_name + "' does not exist for given type.");
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Serialization Methods

    /// @brief Provides serialization support
    template<typename SerializerType>
    void serialize(SerializerType& s) const;

    /// @brief Provides deserialization support
    template<typename SerializerType>
    void deserialize(SerializerType& s);
};


//////////////////////////////////////////////////////////////////////////////////////////////////
// Inline and templated method definitions

#ifdef WIN32
#pragma region MeshAssociatedDataType - inline and templated method definitions
#endif

template<typename T>
inline bool MeshAssociatedDataType<T>::empty() const
{
    return node_data_fields.empty() && element_data_fields.empty() &&
           node_by_element_data_fields.empty() && node_by_component_data_fields.empty() &&
           quad_point_data_fields.empty() && node_group_data_fields.empty() &&
           element_group_data_fields.empty();
}

template<typename T>
inline bool MeshAssociatedDataType<T>::operator==(const MeshAssociatedDataType<T>& other) const
{
    return node_data_fields == other.node_data_fields &&
           element_data_fields == other.element_data_fields &&
           node_by_element_data_fields == other.node_by_element_data_fields &&
           node_by_component_data_fields == other.node_by_component_data_fields &&
           quad_point_data_fields == other.quad_point_data_fields &&
           node_group_data_fields == other.node_group_data_fields &&
           element_group_data_fields == other.element_group_data_fields;
}

template<typename T>
template<typename SerializerType>
void MeshAssociatedDataType<T>::serialize(SerializerType& s) const
{
    s(utility::named_value("n_node_data_fields", node_data_fields.size()),
      utility::named_value("node_data_fields", utility::slice_c(node_data_fields)),
      utility::named_value("n_element_data_fields", element_data_fields.size()),
      utility::named_value("element_data_fields", utility::slice_c(element_data_fields)),
      utility::named_value("n_node_by_element_data_fields", node_by_element_data_fields.size()),
      utility::named_value("node_by_element_data_fields",
                           utility::slice_c(node_by_element_data_fields)),
      utility::named_value("n_node_by_component_data_fields", node_by_component_data_fields.size()),
      utility::named_value("node_by_component_data_fields",
                           utility::slice_c(node_by_component_data_fields)),
      utility::named_value("n_quad_point_data_fields", quad_point_data_fields.size()),
      utility::named_value("quad_point_data_fields", utility::slice_c(quad_point_data_fields)),
      utility::named_value("n_node_group_data_fields", node_group_data_fields.size()),
      utility::named_value("node_group_data_fields", utility::slice_c(node_group_data_fields)),
      utility::named_value("n_element_group_data_fields", element_group_data_fields.size()),
      utility::named_value("element_group_data_fields",
                           utility::slice_c(element_group_data_fields)));
}

template<typename T>
template<typename SerializerType>
void MeshAssociatedDataType<T>::deserialize(SerializerType& s)
{
    std::size_t n_fields;

    s(n_fields);
    node_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        node_data_fields.push_back(DataField<T>::empty());
    s(utility::slice(node_data_fields));

    s(n_fields);
    element_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        element_data_fields.push_back(DataField<T>::empty());
    s(utility::slice(element_data_fields));

    s(n_fields);
    node_by_element_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        node_by_element_data_fields.push_back(DataField<T>::empty());
    s(utility::slice(node_by_element_data_fields));

    s(n_fields);
    node_by_component_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        node_by_component_data_fields.push_back(DataField<T>::empty());
    s(utility::slice(node_by_component_data_fields));

    s(n_fields);
    quad_point_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        quad_point_data_fields.push_back(DataField<T>::empty());
    s(utility::slice(quad_point_data_fields));

    s(n_fields);
    node_group_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        node_group_data_fields.push_back(GroupDataField<T>::empty());
    s(utility::slice(node_group_data_fields));

    s(n_fields);
    element_group_data_fields.reserve(n_fields);
    for (auto i : utility::range(n_fields))
        element_group_data_fields.push_back(GroupDataField<T>::empty());
    s(utility::slice(element_group_data_fields));
}

#ifdef WIN32
#pragma endregion
#endif


#ifdef WIN32
#pragma region MeshAssociatedData - inline and templated method definitions
#endif

template<DataFieldType DFT, typename T>
void add(const DataField<T>& data)
{
    throw _UTIL_MSG_EXCEPT("Specialization for given types is not implemented.");
}

template<DataFieldType DFT, typename T>
inline const std::vector<DataField<T>>& MeshAssociatedData::get() const
{
    throw _UTIL_MSG_EXCEPT("Specialization for given types is not implemented.");
}
// clang-format off
#define ADD_DATA_HELPERS_DEF(TYPE_NAME, TYPE, LEVEL, LEVEL_ENUM)                                   \
    template<>                                                                                     \
    inline void MeshAssociatedData::add<DataFieldType::LEVEL_ENUM, TYPE>(                          \
      DataField<TYPE> && data)                                                                     \
    {                                                                                              \
        TYPE_NAME##_data.LEVEL##_data_fields.push_back(std::move(data));                           \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    inline const std::vector<DataField<TYPE>>&                                                     \
    MeshAssociatedData::get<DataFieldType::LEVEL_ENUM, TYPE>() const                               \
    {                                                                                              \
        return TYPE_NAME##_data.LEVEL##_data_fields;                                               \
    }                                                                                              \
                                                                                                   \
    inline void MeshAssociatedData::add_##LEVEL##_##TYPE_NAME(                                     \
      const std::string& name,                                                                     \
      const std::vector<std::string>& component_names,                                             \
      const simple_math::DynamicMatrix<TYPE>& data)                                                \
    {                                                                                              \
        add<DataFieldType::LEVEL_ENUM, TYPE>(                                                      \
          DataField<TYPE>::cref(name, component_names, data));                                     \
    }                                                                                              \
                                                                                                   \
    inline void MeshAssociatedData::add_##LEVEL##_##TYPE_NAME(const std::string& name,             \
                                                                const std::vector<TYPE>& data)     \
    {                                                                                              \
        add<DataFieldType::LEVEL_ENUM, TYPE>(DataField<TYPE>::cref(name, data));                   \
    }                                                                                              \
                                                                                                   \
    inline void MeshAssociatedData::add_##LEVEL##_##TYPE_NAME(const std::string& name,             \
                                                                const int n_rows)                  \
    {                                                                                              \
        add<DataFieldType::LEVEL_ENUM, TYPE>(                                                      \
          DataField<TYPE>::copy(name, std::vector<TYPE>(n_rows, 0.0)));                            \
    }                                                                                              \
                                                                                                   \
    inline void MeshAssociatedData::add_##LEVEL##_##TYPE_NAME(                                     \
      const std::string& name,                                                                     \
      const std::vector<std::string>& component_names,                                             \
      const int n_rows,                                                                            \
      const int n_cols)                                                                            \
    {                                                                                              \
        add<DataFieldType::LEVEL_ENUM, TYPE>(                                                      \
          DataField<TYPE>::move(name,                                                              \
                                    component_names,                                               \
                                    simple_math::DynamicMatrix<TYPE>::Zero(n_rows, n_cols)));      \
    }                                                                                              \
                                                                                                   \
    inline void MeshAssociatedData::add_##LEVEL##_##TYPE_NAME(                                     \
      const std::string& name,                                                                     \
      const std::vector<std::string>& component_names,                                             \
      simple_math::DynamicMatrix<TYPE>&& data)                                                     \
    {                                                                                              \
        add<DataFieldType::LEVEL_ENUM, TYPE>(                                                      \
          DataField<TYPE>::move(name, component_names, std::move(data)));                          \
    }

// clang-format on

ADD_DATA_HELPERS_DEF(int16, int16_t, node, NODE)
ADD_DATA_HELPERS_DEF(uint16, uint16_t, node, NODE)
ADD_DATA_HELPERS_DEF(int32, int32_t, node, NODE)
ADD_DATA_HELPERS_DEF(uint32, uint32_t, node, NODE)
ADD_DATA_HELPERS_DEF(int64, int64_t, node, NODE)
ADD_DATA_HELPERS_DEF(uint64, uint64_t, node, NODE)
ADD_DATA_HELPERS_DEF(float, float, node, NODE)
ADD_DATA_HELPERS_DEF(double, double, node, NODE)

ADD_DATA_HELPERS_DEF(int16, int16_t, element, ELEMENT)
ADD_DATA_HELPERS_DEF(uint16, uint16_t, element, ELEMENT)
ADD_DATA_HELPERS_DEF(int32, int32_t, element, ELEMENT)
ADD_DATA_HELPERS_DEF(uint32, uint32_t, element, ELEMENT)
ADD_DATA_HELPERS_DEF(int64, int64_t, element, ELEMENT)
ADD_DATA_HELPERS_DEF(uint64, uint64_t, element, ELEMENT)
ADD_DATA_HELPERS_DEF(float, float, element, ELEMENT)
ADD_DATA_HELPERS_DEF(double, double, element, ELEMENT)

ADD_DATA_HELPERS_DEF(int16, int16_t, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(uint16, uint16_t, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(int32, int32_t, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(uint32, uint32_t, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(int64, int64_t, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(uint64, uint64_t, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(float, float, node_by_element, NODE_BY_ELEMENT)
ADD_DATA_HELPERS_DEF(double, double, node_by_element, NODE_BY_ELEMENT)

ADD_DATA_HELPERS_DEF(int16, int16_t, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(uint16, uint16_t, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(int32, int32_t, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(uint32, uint32_t, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(int64, int64_t, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(uint64, uint64_t, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(float, float, node_by_component, NODE_BY_COMPONENT)
ADD_DATA_HELPERS_DEF(double, double, node_by_component, NODE_BY_COMPONENT)

ADD_DATA_HELPERS_DEF(int16, int16_t, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(uint16, uint16_t, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(int32, int32_t, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(uint32, uint32_t, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(int64, int64_t, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(uint64, uint64_t, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(float, float, quad_point, QUAD)
ADD_DATA_HELPERS_DEF(double, double, quad_point, QUAD)

inline bool MeshAssociatedData::empty() const
{
    return uint16_data.empty() && uint32_data.empty() && uint64_data.empty() &&
           float_data.empty() && double_data.empty();
}

template<DataFieldType DFT, typename T>
inline const DataField<T>& MeshAssociatedData::get(const std::string& field_name) const
{
    for (auto& d : get<DFT, T>())
    {
        if (utility::istr_is_equal(d.name, field_name))
        {
            return d;
        }
    }
    throw _UTIL_MSG_EXCEPT_C("Data field '" + field_name + "' does not exist for given type.");
}

template<typename SerializerType>
void MeshAssociatedData::serialize(SerializerType& s) const
{
    s(utility::named_value("n_node_by_element_offsets", node_by_element_offsets.size()),
      utility::named_value("node_by_element_offsets", utility::slice_c(node_by_element_offsets)),
      utility::named_value("n_node_by_component_offsets", node_by_component_offsets.size()),
      utility::named_value("node_by_component_offsets",
                           utility::slice_c(node_by_component_offsets)),
      utility::named_value("n_quad_point_offsets", quad_point_offsets.size()),
      utility::named_value("quad_point_offsets", utility::slice_c(quad_point_offsets)),
      utility::named_ref("int16_data", int16_data),
      utility::named_ref("uint16_data", uint16_data),
      utility::named_ref("int32_data", int32_data),
      utility::named_ref("uint32_data", uint32_data),
      utility::named_ref("int64_data", int64_data),
      utility::named_ref("uint64_data", uint64_data),
      utility::named_ref("float_data", float_data),
      utility::named_ref("double_data", double_data));
}

template<typename SerializerType>
void MeshAssociatedData::deserialize(SerializerType& s)
{
    std::size_t n_offsets;

    s(n_offsets);
    node_by_element_offsets.resize(n_offsets);
    s(utility::slice(node_by_element_offsets));

    s(n_offsets);
    node_by_component_offsets.resize(n_offsets);
    s(utility::slice(node_by_component_offsets));

    s(n_offsets);
    quad_point_offsets.resize(n_offsets);
    s(utility::slice(quad_point_offsets));

    s(int16_data,
      uint16_data,
      int32_data,
      uint32_data,
      int64_data,
      uint64_data,
      float_data,
      double_data);
}

#ifdef WIN32
#pragma endregion
#endif
}

#undef ADD_DATA_HELPERS_DECL
#undef ADD_DATA_HELPERS_DEF