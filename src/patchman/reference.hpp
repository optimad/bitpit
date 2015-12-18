//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_REFERENCE_HPP__
#define __PATCHMAN_REFERENCE_HPP__

/*! \file */

namespace pman {

/*!
	\ingroup PatchMan
	@{
*/

/*!
	@brief Class that emulates a reference to an element of type T.

	@details
	Reference is a class that emulates a reference to an element of type T.
	Unlike a classic reference, a Reference object can be copied (it is
	both copy-constructible and copy-assignable) and can referece to a
	special value that is considered a null element.

	@tparam T is the type of the referred element
*/


template<class T>
class Reference : public std::reference_wrapper<T>
{

public:

	/*!
		Creates a new reference to the specified value.

		If no value is passed, a default value is used.

		\param value is an lvalue reference, whose reference is stored
		             in the object.
	*/
	Reference(T &value = m_null) noexcept
		: std::reference_wrapper<T>(value)
	{

	}

	/*!
		Copy constructor.

		\param x is a Reference object of the same type (i.e.,
		             with the same template parameter, T).
	*/
	Reference(const Reference &x) noexcept
		: std::reference_wrapper<T>(x)
	{

	}

	/*!
		Constructs the null element in place.

		\param args the arguments forwarded to construct the element
	*/
	template <class... Args>
	static void define_null(Args&&... args)
	{
		m_null = T(std::forward<Args>(args)...);
	}

	/*!
		Checks if the reference is null

		\result true if the reference is null, false otherwise.
	*/
	bool is_null()
	{
		return ((*this).get() == m_null);
	}

	/*!
		Nullifies the the reference
	*/
	void nullify()
	{
		(*this) = m_null;
	}

private:
	static T m_null;

	Reference(T &&value) = delete;

};

template<class T>
T Reference<T>::m_null = T();

/*!
	@}
*/

}

#endif
