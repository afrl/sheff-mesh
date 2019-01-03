#pragma once

#include <functional>

/// This provides a set of type trait helpers that is helpful for determining information
/// about classes at compile time, such as if a class contains a method.  These type trait
/// helper than enable the ability to turn on and off methods of a templated class depending
/// on whether the templated types fit some criteria.  A common situation would be to enable
/// one version of a method if the templated type has some method and use a different method
/// if the templated type does not have the method.  The internals of this can be very confusing.
/// The comments will aim to explain to some degree, but focus on the comments describing how
/// to use such traits.
namespace flex_mesh_traits
{
	/// This class provides a way to test if a type contains a specific method.  The macro
	/// generates the class for a member name.  Note: that this version does not check for
	/// a specific signature.  If the class has a member named <MEMBER_NAME> then the value
	/// of the has_<MEMBER_NAME> will be true.
	/// For example:
	///		GENERATE_HAS_MEMBER(size)
	///		std::cout << boolalpha << has_size(std::vector<int>)::value << std::endl;
	/// would print "true".
	///
	/// Internally:
	/// - This strategy relies on multiple inheritance, priority of members, and SFINAE
	/// - Note: SFINAE means substitution failure is not an error.
	/// - has_<MEMBER_NAME>::Derived inherits from both T and Fallback
	/// - Fallback is a struct that contains a member named <MEMBER_NAME>
	/// - If has_<MEMBER_NAME>::Derived::<MEMBER_NAME> results in the
	///   has_<MEMBER_NAME>::Fallback::<MEMBER_NAME> then we know that T does not have a member
	///   called <MEMBER_NAME>, since T::<MEMBER_NAME> would have been given priority.
	/// - has_<MEMBER_NAME>::deduce tests if the has_<MEMBER_NAME>::Derived::<MEMBER_NAME>
	///   matches the template signature with Fallback via template deduction and returns a false_type
	/// - For all other cases, deduce will return std::true_type
	/// - The type of deduce is then used to know if the member exists.
	#define GENERATE_HAS_MEMBER(MEMBER_NAME)							\
	template <typename T>												\
	class has_##MEMBER_NAME												\
	{																	\
		struct Fallback													\
		{																\
			void get_local_number() {}									\
		};																\
		struct Derived : public T, public Fallback {};					\
		template <typename T, T t>  class Helper {};					\
		template <typename U>											\
		static std::false_type deduce(U*,								\
			Helper<void (Fallback::*)(), &U::get_local_number>* = 0);	\
		static std::true_type deduce(...);								\
	public:																\
		using type = decltype(deduce((Derived*)(0)));					\
		static constexpr bool value = type::value;						\
	};																	\

	//GENERATE_HAS_MEMBER(get_local_number);
	//GENERATE_HAS_MEMBER(set_local_number);

	#undef GENERATE_HAS_MEMBER
}
