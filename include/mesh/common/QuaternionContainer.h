#pragma once

#include "mesh/Defines.h"
#include "simple_math/EigenDefines.h"
#include "utility/Range.h"
#include "utility/Slice.h"
#include "utility/serialization/NamedValue.h"

#include <assert.h>

namespace mesh
{
/// @brief An enum describing what level the local coordinate system will be assigned for a
/// component.
///
/// This could be:
/// - None: no LCS for the entire component
/// - Component: one LCS for the entire component
/// - Element: one LCS for each element in the component
/// - Node: one LCS for each node in the component
enum class LCSAssignmentType : char
{
    None,
    Component,
    Element,
    Node,
    ComponentRef,
    ElementRef,
    NodeRef
};

using Quaternion = simple_math::Quaternion;

/// @brief Stores or refers to the quaternions that define the local coordinate system (LCS)
/// throughout a component, element, or node.
///
/// Depending on the LCS assignment type, there might be one quaternion, one per a node, or just a
/// reference to another QuaternionContainer's data.  There is expected to be one
/// QuaternionContainer object per an element (in the mesh object) and one per a field (within the
/// FEA framework).
class QuaternionContainer
{
  public:
    /// Default constructor creates an assignment type of None
    QuaternionContainer();
    /// Constructor to create a reference to another QuaternionContainer \see refer_to()
    QuaternionContainer(const QuaternionContainer& other);
    /// Move constructor (parameter is put in a valid but empty state afterwards)
    QuaternionContainer(QuaternionContainer&& other) noexcept;
    /// Constructor to copy the quaternions from a container.
    /// Be sure that the type matches the number of quaternions given.
    template<typename T>
    QuaternionContainer(const LCSAssignmentType& type, const T& quaternions);

    void operator=(QuaternionContainer&& other) noexcept;
    bool operator==(const QuaternionContainer& other) const;

    /// Initializes the container to a given assignment type and set all quaternions
    /// to be an identity rotation.
    void initialize(LCSAssignmentType new_type, ConnIndex n_nodes = 0);

    /// Sets all quaternions to be an identity rotation.
    void set_identity();

    /// Sets this QuaternionContainer to refer to another QuaternionContainer.
    /// The other object will retain ownership of the quaternion data,
    /// but this object will now point to the data owned by the owner.
    /// Obviously, this creates a problem if a quaternion is requested
    /// after the object being referred to is deleted, so be sure the
    /// algorithms and structures are designed to avoid such a situation.
    void refer_to(const QuaternionContainer& other);

    /// Copies all data from another QuaternionContainer to this container.  This does not create a
    /// reference.  Important note: if other QuaternionContainer objects already refer to this, then
    /// calling copy will invalidate the other referring QuaternionContainer's.
    void copy(const QuaternionContainer& other);

    /// Begin iterator, exists to enable range loops:
    /// \code{cpp}
    /// QuaternionContainer A;
    /// // initialize A...
    /// for (auto& quaternion : A)
    /// {
    ///     // do something with quaternion ...
    /// }
    /// \endcode
    Quaternion* begin() const;

    /// End iterator, exists to enable range loops, \see QuaternionContainer::begin()
    Quaternion* end() const;

    /// Number of quaternions in the container
    ConnIndex size() const;

    /// Access a particular quaternion in the container given an index.
    /// The returned value depends on the LCS assignment type:
    /// - None: an identity quanternion will be returned (i is ignored).
    /// - Component and Element: the first quanternion will be returned (i is ignored).
    /// - Node: the requested index of the container will be returned
    ///
    /// Note: this does not return the quaternion by reference.
    Quaternion operator[](const ConnIndex i) const;

    /// Yeilds the assignment type (dropping any reference specialization).
    /// For example, a NodeRef will return Node.
    LCSAssignmentType base_type() const;

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("type", type_),
          utility::named_ref("n_quaternions", n_quaternions_),
          utility::named_value("quaternions",
                               utility::Slice<const Quaternion>(quaternions_, n_quaternions_)));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(type_, n_quaternions_);
        quaternions_ = new Quaternion[n_quaternions_];
        s(utility::Slice<Quaternion>(quaternions_, n_quaternions_));
    }

  private:
    LCSAssignmentType type_;
    Quaternion* quaternions_;
    ConnIndex n_quaternions_;

    /// Delete the data if it is initialized and not a reference type.
    /// Otherwise, it does nothing.
    void cleanup_if_needed();
};


//////////////////////////////////////////////////////////////////////////////////////////////
// Template Definitions

template<typename T>
inline QuaternionContainer::QuaternionContainer(const LCSAssignmentType& type, const T& quaternions)
{
    switch (type)
    {
        case LCSAssignmentType::None: assert(quaternions.size() == 0); break;
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef:
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef: assert(quaternions.size() == 1); break;
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef: assert(quaternions.size() > 1); break;
    }

    n_quaternions_ = quaternions.size();
    quaternions_   = new Quaternion[n_quaternions_];
    for (auto i : utility::range(n_quaternions_))
    {
        quaternions_[i] = quaternions[i];
    }
}

}
