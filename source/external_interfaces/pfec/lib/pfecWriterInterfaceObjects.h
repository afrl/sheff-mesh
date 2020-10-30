#pragma once

namespace PFECNativeInterface {

    // Inspired by http://stackoverflow.com/questions/470835/a-c-iterator-adapter-which-wraps-and-hides-an-inner-iterator-and-converts-the

    /// The MeshWriter allows you to use this interface to write nodes instead of allocating vectors of data
    /// This is far more efficient if you have a custom node class in some container.
    /// If the container is not STL, then you will need to define an iterator inheriting from the std::iterator
    /// Create an interface what whatever node class you may have, and overload all of the methods
    /// Note: do not delete obj
    template<typename InputIterator>
    class NodeIteratorAdapter {
    public:
		NodeIteratorAdapter() {};
		NodeIteratorAdapter(const InputIterator &it) {
			wrapped_iter_ = it;
        }
        virtual PFEC_IND get_global_ind() = 0;
		virtual PFEC_INT get_owning_partition() = 0;
        virtual PFEC_INT get_coord_len() = 0;
        virtual const double* get_coords() = 0;
		NodeIteratorAdapter& operator++ () {
			wrapped_iter_++;
			return *this;
		}
		bool operator==(NodeIteratorAdapter const & other) const {
			return wrapped_iter_ == other.wrapped_iter_;
		}
		bool operator!=(NodeIteratorAdapter const & other) const {
			return wrapped_iter_ != other.wrapped_iter_;
		}
	protected:
		InputIterator wrapped_iter_;
    };

    template<typename InputIterator>
    class ElementIteratorAdapter {
    public:
		ElementIteratorAdapter() {};
		ElementIteratorAdapter(const InputIterator &it) {
			wrapped_iter_ = it;
        }
        virtual PFEC_IND get_global_ind() = 0;
		virtual PFEC_INT get_owning_partition() = 0;
		virtual PFEC_INT get_domain_id() = 0;
        virtual ElemGeomType get_geom_type() = 0;
        virtual PFEC_INT get_local_conn_len() = 0;
        // Connectivity defined by local partition node numbering (not global)
        virtual PFEC_INT * get_local_conn() = 0;

		ElementIteratorAdapter& operator++ () {
			wrapped_iter_++;
			return *this;
		}
		bool operator==(ElementIteratorAdapter const & other) const {
			return wrapped_iter_ == other.wrapped_iter_;
		}
		bool operator!=(ElementIteratorAdapter const & other) const {
			return wrapped_iter_ != other.wrapped_iter_;
		}
	protected:
		InputIterator wrapped_iter_;
    };
};