#pragma once

#include "simple_math/EigenDefines.h"

template<typename T>
class EigenMatrixHandleBase
{
  public:
    virtual ~EigenMatrixHandleBase() {}

    inline Eigen::Index cols() const
    {
        return ref_to_ref_.cols();
    }
    inline Eigen::Index rows() const
    {
        return ref_to_ref_.rows();
    }
    inline const T* data() const
    {
        return ref_to_ref_.data();
    }
    inline const Eigen::Ref<const simple_math::DynamicMatrix<T>>& get_const() const
    {
        return ref_to_ref_;
    }

    virtual Eigen::Ref<simple_math::DynamicMatrix<T>> get_mutable()
    {
        throw _UTIL_MSG_EXCEPT_C(
          "Non-const access to data is not allowed for this type of handle.");
    }

  protected:
    EigenMatrixHandleBase(const Eigen::Ref<const simple_math::DynamicMatrix<T>>& ref)
      : ref_to_ref_(ref)
    {
    }

    const Eigen::Ref<const simple_math::DynamicMatrix<T>>& ref_to_ref_;
};

/// @brief Provides a handle that does not own the underlying Eigen matrix and prevents modication
/// to the underlying data.
///
/// In almost every case, a reference to a particular Eigen class or Eigen::Ref can be used.
/// However, this class provides the very specific use case when something must hold onto a matrix
/// that it either owns or does not own for an indefinite time, while hiding the detail of
/// ownership.
template<typename T>
class EigenMatrixConstHandle : public EigenMatrixHandleBase<T>
{
  public:
    EigenMatrixConstHandle(const Eigen::Ref<const simple_math::DynamicMatrix<T>>& data)
      : EigenMatrixHandleBase<T>(const_ref_)
      , const_ref_(data)
    {
    }
    virtual ~EigenMatrixConstHandle() {}

  protected:
    Eigen::Ref<const simple_math::DynamicMatrix<T>> const_ref_;
};

/// @brief Provides a handle that does not own the underlying Eigen matrix but allows modication
/// to the underlying data.
template<typename T>
class EigenMatrixMutableHandle : public EigenMatrixConstHandle<T>
{
  public:
    EigenMatrixMutableHandle(simple_math::DynamicMatrix<T>& data)
      : EigenMatrixConstHandle<T>(data)
      , ref_(data)
    {
    }
    virtual ~EigenMatrixMutableHandle() {}

    virtual Eigen::Ref<simple_math::DynamicMatrix<T>> get_mutable() override
    {
        return ref_;
    }

  protected:
    Eigen::Ref<simple_math::DynamicMatrix<T>> ref_;
};

/// @brief A handle that does not own the underlying data, prevents modication to the data, and
/// wraps the data with an Eigen matrix.
template<typename T>
class EigenMatrixConstMapHandle : public EigenMatrixHandleBase<T>
{
  public:
    /// @param data - row major data, for example: EigenMatrixConstMapHandle([ 0, 1, 2, 3, 4, 5 ],
    ///               2, 3) results in the matrix [ [0, 1, 2], [3, 4, 5] ]
    /// @param n_rows - number of rows
    /// @param n_cols - number of columns
    EigenMatrixConstMapHandle(const T* data, const Eigen::Index& n_rows, const Eigen::Index& n_cols)
      : EigenMatrixHandleBase<T>(const_ref_)
      , data_map_(data, n_rows, n_cols)
      , const_ref_(data_map_)
    {
    }
    virtual ~EigenMatrixConstMapHandle() {}

  protected:
    Eigen::Map<const simple_math::DynamicMatrix<T>> data_map_;
    Eigen::Ref<const simple_math::DynamicMatrix<T>> const_ref_;
};

/// @brief A handle that does not own the underlying data, prevents modication to the data, and
/// wraps the data with an Eigen matrix.
template<typename T>
class EigenMatrixMutableMapHandle : public EigenMatrixHandleBase<T>
{
  public:
    /// @param data - row major data, for example: EigenMatrixConstMapHandle([ 0, 1, 2, 3, 4, 5 ],
    ///               2, 3) results in the matrix [ [0, 1, 2], [3, 4, 5] ]
    /// @param n_rows - number of rows
    /// @param n_cols - number of columns
    EigenMatrixMutableMapHandle(T* data, const Eigen::Index& n_rows, const Eigen::Index& n_cols)
      : EigenMatrixHandleBase<T>(const_ref_)
      , data_map_(data, n_rows, n_cols)
      , const_ref_(data_map_)
    {
    }
    virtual ~EigenMatrixMutableMapHandle() {}

    virtual Eigen::Ref<simple_math::DynamicMatrix<T>> get_mutable() override
    {
        return data_map_;
    }

  protected:
    Eigen::Map<simple_math::DynamicMatrix<T>> data_map_;
    Eigen::Ref<const simple_math::DynamicMatrix<T>> const_ref_;
};

/// @brief A handle that owns the underlying Eigen matrix and allows modification.
template<typename T>
class EigenMatrixOwningHandle : public EigenMatrixHandleBase<T>
{
  public:
    EigenMatrixOwningHandle(const simple_math::DynamicMatrix<T>& data)
      : EigenMatrixHandleBase<T>(const_ref_)
      , data_(data)
      , const_ref_(data_)
    {
    }
    EigenMatrixOwningHandle(simple_math::DynamicMatrix<T>&& data)
      : EigenMatrixHandleBase<T>(const_ref_)
      , data_(std::move(data))
      , const_ref_(data_)
    {
    }
    virtual ~EigenMatrixOwningHandle() {}

    virtual Eigen::Ref<simple_math::DynamicMatrix<T>> get_mutable() override
    {
        return data_;
    }

  protected:
    simple_math::DynamicMatrix<T> data_;
    Eigen::Ref<const simple_math::DynamicMatrix<T>> const_ref_;
};