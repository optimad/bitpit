/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __BITPIT_PIERCED_ITERATOR_HPP__
#define __BITPIT_PIERCED_ITERATOR_HPP__

#define  __PI_REFERENCE__ typename PiercedIterator<value_t, id_t, value_no_cv_t>::reference
#define  __PI_POINTER__   typename PiercedIterator<value_t, id_t, value_no_cv_t>::pointer

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <memory>
#include <type_traits>
#include <vector>

namespace bitpit{

template<typename value_t, typename id_t>
class PiercedVector;

/*!
	\ingroup containers

	@brief Iterator for the class PiercedVector

	Usage: use <tt>PiercedVector<value_t, id_t>::iterator</tt>
	to declare an iterator for a pierced vector, use
	<tt>PiercedVector<vaue_t, id_t>::const_iterator</tt> to declare
	a const iterator for a pierced vector.

	@tparam value_t The type of the elements stored in the vector
	@tparam id_t The type of the ids to associate to the elements
*/
template<typename value_t, typename id_t = long,
         typename value_no_cv_t = typename std::remove_cv<value_t>::type>
class PiercedIterator
	: public std::iterator<std::forward_iterator_tag, value_no_cv_t, std::ptrdiff_t, value_t*, value_t&>
{
	static_assert(std::is_integral<id_t>::value, "Signed integer required for id.");
	static_assert(std::numeric_limits<id_t>::is_signed, "Signed integer required for id.");

	 friend class PiercedIterator<value_no_cv_t, id_t>;

private:
	/*!
		Container.
	*/
	template<typename PV_value_t, typename PV_id_t>
	using Container = PiercedVector<PV_value_t, PV_id_t>;

	/*
		Container type

		When building a const_iterator the pointer to the container has to
		be declared const.
	*/
	typedef
		typename std::conditional<std::is_const<value_t>::value,
			const Container<value_no_cv_t, id_t>,
			Container<value_no_cv_t, id_t>
		>::type

		container_t;

public:
	/*!
		Constant container
	*/
	typedef container_t container_type;

	/*!
		Type of data stored in the container
	*/
	typedef value_t value_type;

	/*!
		Type of ids stored in the container
	*/
	typedef id_t id_type;

	// Friendships
	template<typename PV_value_t, typename PV_id_t>
	friend class PiercedVector;

	// Constructors
	PiercedIterator();

	// General methods
	void swap(PiercedIterator& other) noexcept;

	container_type & getContainer() const;

	// Methos to extract information on the current element
	id_t getId(const id_t &fallback = -1) const noexcept;
	std::size_t getRawIndex() const noexcept;

	// Operators
	PiercedIterator& operator++();
	PiercedIterator operator++(int);

	__PI_REFERENCE__ operator*() const;
	__PI_POINTER__ operator->() const;

	template<typename U = value_t, typename U_no_cv = value_no_cv_t,
             typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type = 0>
	operator PiercedIterator<const U_no_cv, id_t>() const;

	/*!
		Two-way comparison.
	*/
	template<typename other_value_t, typename other_id_t = long,
         typename other_value_no_cv_t = typename std::remove_cv<value_t>::type>
	bool operator==(const PiercedIterator<other_value_t, other_id_t, other_value_no_cv_t>& rhs) const
	{
		return (m_container == rhs.m_container) && (m_pos == rhs.m_pos);
	}

	/*!
		Two-way comparison.
	*/
	template<typename other_value_t, typename other_id_t = long,
         typename other_value_no_cv_t = typename std::remove_cv<value_t>::type>
	bool operator!=(const PiercedIterator<other_value_t, other_id_t, other_value_no_cv_t>& rhs) const
	{
		return (m_container != rhs.m_container) || (m_pos != rhs.m_pos);
	}

private:
	/*!
		Internal pointer to the container.
	*/
	container_t *m_container;

	/*!
		Position inside the container.
	*/
	size_t m_pos;

	// Constructors
	explicit PiercedIterator(container_t *container, const size_t &pos);

};

}

// Include the implementation
#include "piercedIterator.tpp"

#endif
