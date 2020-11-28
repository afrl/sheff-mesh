#include "mesh/common/QuaternionContainer.h"

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////
// Static helper free-functions

/// Yeilds the equivalent assignment type to the one given that is based on a reference,
/// regardless of whether the given type is based on a reference or not.
static LCSAssignmentType get_ref_equivalent(const LCSAssignmentType& type)
{
    switch (type)
    {
        case LCSAssignmentType::None: return LCSAssignmentType::None;
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef: return LCSAssignmentType::ComponentRef;
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef: return LCSAssignmentType::ElementRef;
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef: return LCSAssignmentType::NodeRef;
        default: throw _UTIL_MSG_EXCEPT("The type is not recognized.");
    }
}

/// Yeilds the equivalent assignment type to the one given that is NOT based on a
/// reference, regardless of whether the given type is based on a reference or not.
static LCSAssignmentType get_non_ref_equivalent(const LCSAssignmentType& type)
{
    switch (type)
    {
        case LCSAssignmentType::None: return LCSAssignmentType::None;
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef: return LCSAssignmentType::Component;
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef: return LCSAssignmentType::Element;
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef: return LCSAssignmentType::Node;
        default: throw _UTIL_MSG_EXCEPT("The type is not recognized.");
    }
}

/// Determines if an assignment type that is based on a reference.
static bool is_ref(const LCSAssignmentType& type)
{
    switch (type)
    {
        case LCSAssignmentType::None:
        case LCSAssignmentType::Component:
        case LCSAssignmentType::Element:
        case LCSAssignmentType::Node: return false;
        default: return true;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////
// QuaternionContainer methods

QuaternionContainer::QuaternionContainer()
  : quaternions_(nullptr)
  , n_quaternions_(0)
  , type_(LCSAssignmentType::None)
{
}

QuaternionContainer::QuaternionContainer(const QuaternionContainer& other)
  : quaternions_(nullptr)
  , n_quaternions_(0)
  , type_(LCSAssignmentType::None)
{
    refer_to(other);
}

QuaternionContainer::QuaternionContainer(QuaternionContainer&& other) noexcept
  : quaternions_(std::move(other.quaternions_))
  , n_quaternions_(std::move(other.n_quaternions_))
  , type_(std::move(other.type_))
{
}

void QuaternionContainer::operator=(QuaternionContainer&& other) noexcept
{
    quaternions_   = std::move(other.quaternions_);
    n_quaternions_ = std::move(other.n_quaternions_);
    type_          = std::move(other.type_);
}

bool QuaternionContainer::operator==(const QuaternionContainer& other) const
{
    bool comparison = type_ == other.type_ && n_quaternions_ == other.n_quaternions_;
    for (auto i : utility::range(n_quaternions_))
    {
        comparison &= quaternions_[i].isApprox(other.quaternions_[i]);
    }
    return comparison;
}

void mesh::QuaternionContainer::initialize(LCSAssignmentType new_type, ConnIndex n_nodes)
{
    if (new_type != type_)
    {
        cleanup_if_needed();

        type_ = new_type;
        switch (type_)
        {
            case LCSAssignmentType::None: n_quaternions_ = 0; break;
            case LCSAssignmentType::Component:
            case LCSAssignmentType::ComponentRef:
            case LCSAssignmentType::Element:
            case LCSAssignmentType::ElementRef: n_quaternions_ = 1; break;
            case LCSAssignmentType::Node:
            case LCSAssignmentType::NodeRef: assert(n_nodes > 0); n_quaternions_ = n_nodes;
        }

        quaternions_ = new Quaternion[n_quaternions_];
    }

    for (auto& quaternion : *this)
    {
        quaternion.setIdentity();
    }
}

void mesh::QuaternionContainer::set_identity()
{
    for (auto& quaternion : *this)
    {
        quaternion.setIdentity();
    }
}

void mesh::QuaternionContainer::refer_to(const QuaternionContainer& other)
{
    cleanup_if_needed();
    type_          = get_ref_equivalent(other.type_);
    quaternions_   = other.quaternions_;
    n_quaternions_ = other.n_quaternions_;
}

void QuaternionContainer::copy(const QuaternionContainer& other)
{
    cleanup_if_needed();
    type_ = get_non_ref_equivalent(other.type_);
    initialize(type_, other.size());
    for (auto i : utility::range(n_quaternions_))
    {
        quaternions_[i] = other.quaternions_[i];
    }
}

Quaternion* mesh::QuaternionContainer::begin() const
{
    return quaternions_;
}

Quaternion* mesh::QuaternionContainer::end() const
{
    return &quaternions_[n_quaternions_];
}

ConnIndex mesh::QuaternionContainer::size() const
{
    return n_quaternions_;
}

Quaternion mesh::QuaternionContainer::operator[](const ConnIndex i) const
{
    switch (type_)
    {
        case LCSAssignmentType::None: return Quaternion::Identity();
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef:
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef: assert(n_quaternions_ == 1); return quaternions_[0];
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef:
            assert(i >= 0 && i < n_quaternions_);
            return quaternions_[i];
        default: throw _UTIL_MSG_EXCEPT("The type is not recognized.");
    }
}

void QuaternionContainer::set(const ConnIndex i, const Quaternion& q)
{
    switch (type_)
    {
        case LCSAssignmentType::None: return;
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef:
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef:
            assert(n_quaternions_ == 1);
            quaternions_[0] = q;
            return;
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef:
            assert(i >= 0 && i < n_quaternions_);
            quaternions_[i] = q;
            return;
        default: throw _UTIL_MSG_EXCEPT("The type is not recognized.");
    }
}

LCSAssignmentType mesh::QuaternionContainer::base_type() const
{
    return get_non_ref_equivalent(type_);
}

void mesh::QuaternionContainer::cleanup_if_needed()
{
    if (quaternions_ != nullptr && !is_ref(type_))
    {
        quaternions_   = new Quaternion[n_quaternions_];
        n_quaternions_ = 0;
    }
}

}