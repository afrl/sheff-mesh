#pragma once

#include "flex_mesh/common/Defines.h"
#include "simple_math/EigenDefines.h"
#include "utility/Range.h"

#include <assert.h>

namespace flex_mesh
{
	/// An enum describing what level the local coordinate system will be
	/// assigned for a component.  This could be:
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

    /// Yeilds the equivalent assignment type to the one given that is based on a reference,
    /// regardless of whether the given type is based on a reference or not.
    static LCSAssignmentType get_ref_equivalent(const LCSAssignmentType& type)
    {
        switch (type)
        {
        case LCSAssignmentType::None:
            return LCSAssignmentType::None;
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef:
            return LCSAssignmentType::ComponentRef;
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef:
            return LCSAssignmentType::ElementRef;
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef:
            return LCSAssignmentType::NodeRef;
        default:
            throw _UTIL_MSG_EXCEPT("The type is not recognized.");
        }
    }

    /// Yeilds the equivalent assignment type to the one given that is NOT based on a
    /// reference, regardless of whether the given type is based on a reference or not.
    static LCSAssignmentType get_non_ref_equivalent(const LCSAssignmentType& type)
    {
        switch (type)
        {
        case LCSAssignmentType::None:
            return LCSAssignmentType::None;
        case LCSAssignmentType::Component:
        case LCSAssignmentType::ComponentRef:
            return LCSAssignmentType::Component;
        case LCSAssignmentType::Element:
        case LCSAssignmentType::ElementRef:
            return LCSAssignmentType::Element;
        case LCSAssignmentType::Node:
        case LCSAssignmentType::NodeRef:
            return LCSAssignmentType::Node;
        default:
            throw _UTIL_MSG_EXCEPT("The type is not recognized.");
        }
    }

    /// Determines if an assignment type that is based on a reference.
    static const bool is_ref(const LCSAssignmentType& type)
    {
        switch (type)
        {
        case LCSAssignmentType::None:
        case LCSAssignmentType::Component:
        case LCSAssignmentType::Element:
        case LCSAssignmentType::Node:
            return false;
        default:
            return true;
        }
    }

    using Quaternion = simple_math::Quaternion;

    /// Stores or refers to the quaternions that define the local coordinate system (LCS)
    /// throughout a component, element, or node.  Depending on the LCS assignment type,
    /// there might be one quaternion, one per a node, or just a reference to another
    /// QuaternionContainer's data.  There is expected to be one QuaternionContainer object
    /// per an element (in the mesh object) and one per a field (within the FEA framework).
	class QuaternionContainer
	{
	public:
        /// Default constructor creates an assignment type of None
        QuaternionContainer()
            : quaternions_(nullptr),
              n_quaternions_(0),
              type_(LCSAssignmentType::None)
        {
        }

        /// Constructor to copy the quaternions from a container.
        /// Be sure that the type matches the number of quaternions given.
		template<typename T>
        QuaternionContainer(
            const LCSAssignmentType& type,
            const T& quaternions)
        {
            switch (type)
            {
            case LCSAssignmentType::None:
                assert(quaternions.size() == 0);
                break;
            case LCSAssignmentType::Component:
            case LCSAssignmentType::ComponentRef:
            case LCSAssignmentType::Element:
            case LCSAssignmentType::ElementRef:
                assert(quaternions.size() == 1);
                break;
            case LCSAssignmentType::Node:
            case LCSAssignmentType::NodeRef:
                assert(quaternions.size() > 1);
                break;
            }

            n_quaterions_ = quaternions.size();
            quaternions_ = new Quaternion[n_quaternions_];
            for (auto i : utility::range_from_zero(n_quaterions_))
            {
                quaternions_[i] = quaternions[i];
            }
        }

        /// Constructor to create a reference to another QuaternionContainer \see refer_to()
        QuaternionContainer(const QuaternionContainer& other)
            : quaternions_(nullptr),
              n_quaternions_(0),
              type_(LCSAssignmentType::None)
        {
            refer_to(other);
        }

        QuaternionContainer(QuaternionContainer&& other) noexcept
            : quaternions_(std::move(other.quaternions_)),
              n_quaternions_(std::move(other.n_quaternions_)),
              type_(std::move(other.type_)) {}

        /// Initializes the container to a given assignment type and set all quaternions
        /// to be an identity rotation.
        void initialize(
            LCSAssignmentType new_type,
            ConnIndex n_nodes = 0)
        {
            if (new_type != type_)
            {
                cleanup_if_needed();

                type_ = new_type;
                switch (type_)
                {
                case LCSAssignmentType::None:
                    n_quaternions_ = 0;
                    break;
                case LCSAssignmentType::Component:
                case LCSAssignmentType::ComponentRef:
                case LCSAssignmentType::Element:
                case LCSAssignmentType::ElementRef:
                    n_quaternions_ = 1;
                    break;
                case LCSAssignmentType::Node:
                case LCSAssignmentType::NodeRef:
                    assert(n_nodes > 0);
                    n_quaternions_ = n_nodes;
                }

                quaternions_ = new Quaternion[n_quaternions_];
            }

            for (auto& quaternion : *this)
            {
                quaternion.setIdentity();
            }
        }

        /// Sets all quaternions to be an identity rotation.
        void set_identity()
        {
            for (auto& quaternion : *this)
            {
                quaternion.setIdentity();
            }
        }

        /// Sets this QuaternionContainer to refer to another QuaternionContainer.
        /// The other object will retain ownership of the quaternion data,
        /// but this object will now point to the data owned by the owner.
        /// Obviously, this creates a problem if a quaternion is requested
        /// after the object being referred to is deleted, so be sure the
        /// algorithms and structures are designed to avoid such a situation.
        void refer_to(const QuaternionContainer& other)
        {
            cleanup_if_needed();
            type_ = get_ref_equivalent(other.type_);
            quaternions_ = other.quaternions_;
            n_quaternions_ = other.n_quaternions_;
        }

        /// Begin iterator, exists to enable range loops:
        /// \code{cpp}
        /// QuaternionContainer A;
        /// // initialize A...
        /// for (auto& quaternion : A)
        /// {
        ///     // do something with quaternion ...
        /// }
        /// \endcode
        Quaternion* const begin() const
        {
            return quaternions_;
        }

        /// End iterator, exists to enable range loops, \see QuaternionContainer::begin()
        Quaternion* const end() const
        {
            return &quaternions_[n_quaternions_];
        }

        ConnIndex size() const
        {
            return n_quaternions_;
        }

        /// Access a particular quaternion in the container given an index.
        /// The returned value depends on the LCS assignment type:
        /// - None: an identity quanternion will be returned (i is ignored).
        /// - Component and Element: the first quanternion will be returned (i is ignored).
        /// - Node: the requested index of the container will be returned
        ///
        /// Note: this does not return the quaternion by reference.
        Quaternion operator[](const ConnIndex i) const
        {
            switch (type_)
            {
            case LCSAssignmentType::None:
                return Quaternion::Identity();
            case LCSAssignmentType::Component:
            case LCSAssignmentType::ComponentRef:
            case LCSAssignmentType::Element:
            case LCSAssignmentType::ElementRef:
                assert(n_quaternions_ == 1);
                return quaternions_[0];
            case LCSAssignmentType::Node:
            case LCSAssignmentType::NodeRef:
                assert(i >= 0 && i < n_quaternions_);
                return quaternions_[i];
            default:
                throw _UTIL_MSG_EXCEPT("The type is not recognized.");
            }
        }

        /// Yeilds the assignment type (dropping any reference specialization).
        /// For example, a NodeRef will return Node.
        LCSAssignmentType base_type() const
        {
            return get_non_ref_equivalent(type_);
        }

    private:
        LCSAssignmentType type_;
        Quaternion* quaternions_;
        ConnIndex n_quaternions_;

        /// Delete the data if it is initialized and not a reference type.
        /// Otherwise, it does nothing.
        void cleanup_if_needed()
        {
            if (quaternions_ != nullptr && !is_ref(type_))
            {
                quaternions_ = new Quaternion[n_quaternions_];
                n_quaternions_ = 0;
            }
        }
	};
}