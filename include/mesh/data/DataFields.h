#pragma once

#include "EigenMatrixHandle.h"
#include "Enums.h"
#include "utility/Slice.h"
#include "utility/serialization/NamedValue.h"

#include <memory>
#include <string>
#include <vector>

namespace mesh
{
/// @brief Holds the actual data and meta data for a data field.
template<typename T>
struct DataField
{
    /// @brief A friendly name for the data field.
    std::string name;
    /// @brief Names of the columns of data.
    std::vector<std::string> component_names;
    /// @brief Data laid out as a matrix where rows correspond to the values, while the columns
    /// correspond to the components.
    std::unique_ptr<EigenMatrixHandleBase<T>> data_handle;

    inline auto& data() const
    {
        return data_handle->get_const();
    }

    inline auto data_mutable()
    {
        return data_handle->get_mutable();
    }

    inline bool operator==(const DataField<T>& other) const
    {
        return name == other.name && component_names == other.component_names &&
               data() == other.data();
    }


    //////////////////////////////////////////////////////////////////////////////////////
    // Serialization Methods

    /// @brief Provides serialization support
    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("name", name),
          utility::named_value("n_component_names", component_names.size()),
          utility::named_value("component_names", utility::slice(component_names)),
          utility::named_value("n_data_rows", data().rows()),
          utility::named_value("data",
                               utility::slice(data().data(), data().rows() * data().cols())));
    }

    /// @brief Provides deserialization support
    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(name);

        std::size_t n_component_names;
        s(n_component_names);

        component_names.resize(n_component_names);
        s(utility::slice(component_names));

        Eigen::Index n_data_rows;
        s(n_data_rows);

        simple_math::DynamicMatrix<T> data =
          simple_math::DynamicMatrix<T>::Zero(n_data_rows, n_component_names);
        s(utility::slice(data.data(), data.rows() * data.cols()));

        data_handle = std::make_unique<EigenMatrixOwningHandle<T>>(data);
    }

    ///////////////////////////////////////////
    // Factories that own the underlying data

    static DataField<T> empty()
    {
        return { "",
                 {},
                 std::make_unique<EigenMatrixOwningHandle<T>>(simple_math::DynamicMatrix<T>()) };
    }

    static DataField<T> copy(const std::string& name, const std::vector<T>& data)
    {
        return { name,
                 { "" },
                 std::make_unique<EigenMatrixOwningHandle<T>>(
                   Eigen::Map<const simple_math::DynamicMatrix<T>>(data.data(), data.size(), 1)) };
    }
    static DataField<T> copy(const std::string& name,
                             const std::vector<std::string>& component_names,
                             const T*& data,
                             const std::size_t& n_rows,
                             const std::size_t& n_cols)
    {
        return { name,
                 component_names,
                 std::make_unique<EigenMatrixOwningHandle<T>>(
                   Eigen::Map<const simple_math::DynamicMatrix<T>>(data, n_rows, n_cols)) };
    }
    static DataField<T> copy(const std::string& name,
                             const std::vector<std::string>& component_names,
                             const simple_math::DynamicMatrix<T>& data)
    {
        EigenMatrixOwningHandle<T> handle(data);
        return { name, component_names, std::make_unique<EigenMatrixOwningHandle<T>>(data) };
    }

    static DataField<T> move(const std::string& name,
                             const std::vector<std::string>& component_names,
                             simple_math::DynamicMatrix<T>&& data)
    {
        return { name, component_names, std::make_unique<EigenMatrixOwningHandle<T>>(data) };
    }

    ///////////////////////////////////////////
    // Factories to create references to data

    // Const references
    static DataField<T> cref(const std::string& name, const std::vector<T>& data)
    {
        return { name,
                 { "" },
                 std::make_unique<EigenMatrixConstMapHandle<T>>(data.data(), data.size(), 1) };
    }
    static DataField<T> cref(const std::string& name, const T* data, const std::size_t size)
    {
        return { name, { "" }, std::make_unique<EigenMatrixConstMapHandle<T>>(data, size, 1) };
    }
    static DataField<T> cref(const std::string& name, const T* data, const Eigen::Index size)
    {
        return { name, { "" }, std::make_unique<EigenMatrixConstMapHandle<T>>(data, size, 1) };
    }
    static DataField<T> cref(const std::string& name, const simple_math::DynamicVector<T>& data)
    {
        return DataField<T>::cref(name, data.data(), data.rows());
    }
    static DataField<T> cref(const std::string& name,
                             const std::vector<std::string>& component_names,
                             const T*& data,
                             const std::size_t& n_rows,
                             const std::size_t& n_cols)
    {
        return { name,
                 component_names,
                 std::make_unique<EigenMatrixConstMapHandle<T>>(data, n_rows, n_cols) };
    }
    static DataField<T> cref(const std::string& name,
                             const std::vector<std::string>& component_names,
                             const simple_math::DynamicMatrix<T>& data)
    {
        return { name, component_names, std::make_unique<EigenMatrixConstHandle<T>>(data) };
    }

    // Mutable references
    static DataField<T> ref(const std::string& name, std::vector<T>& data)
    {
        return { name,
                 { "" },
                 std::make_unique<EigenMatrixMutableMapHandle<T>>(data.data(),
                                                                  (Eigen::Index)data.size(),
                                                                  1) };
    }
    static DataField<T> ref(const std::string& name, T* data, const std::size_t& size)
    {
        return { name, { "" }, std::make_unique<EigenMatrixMutableMapHandle<T>>(data, size, 1) };
    }
    static DataField<T> ref(const std::string& name, simple_math::DynamicVector<T>& data)
    {
        return { name,
                 { "" },
                 std::make_unique<EigenMatrixMutableMapHandle<T>>(data.data(), data.rows(), 1) };
    }
    static DataField<T> ref(const std::string& name,
                            const std::vector<std::string>& component_names,
                            T* data,
                            const std::size_t& n_rows,
                            const std::size_t& n_cols)
    {
        return { name,
                 component_names,
                 std::make_unique<EigenMatrixMutableMapHandle<T>>(data,
                                                                  (Eigen::Index)n_rows,
                                                                  (Eigen::Index)n_cols) };
    }
    static DataField<T> ref(const std::string& name,
                            const std::vector<std::string>& component_names,
                            simple_math::DynamicMatrix<T>& data)
    {
        return { name, component_names, std::make_unique<EigenMatrixMutableHandle<T>>(data) };
    }
};

/// @brief Contains a data field for node groups or element groups.
template<typename T>
struct GroupDataField
{
    std::string group_name;
    DataField<T> data_field;

    inline bool operator==(const GroupDataField<T>& other) const
    {
        return group_name == other.group_name && data_field == other.data_field;
    }

    static GroupDataField<T> empty()
    {
        return { "", DataField<T>::empty() };
    }


    //////////////////////////////////////////////////////////////////////////////////////
    // Serialization Methods

    /// @brief Provides serialization support
    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("group_name", group_name),
          utility::named_ref("data_field", data_field));
    }

    /// @brief Provides deserialization support
    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(group_name);
        s(data_field);
    }
};

}