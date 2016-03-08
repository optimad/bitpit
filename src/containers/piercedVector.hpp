/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __BITPIT_PIERCED_VECTOR_TPP__
#define __BITPIT_PIERCED_VECTOR_TPP__

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <deque>
#include <iostream>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

namespace bitpit{

template<typename value_t, typename id_t>
class PiercedVector;

/*!
	\ingroup containers

	@brief Iterator for the class PiercedVector

	Usage: Use <tt>PiercedVector<value_t>::iterator</tt> to declare an iterator
	for a pierced vector, use <tt>PiercedVector<Type>::const_iterator</tt> to
	declare a const iterator for a pierced vector.

	@tparam value_t The type of the objects stored in the vector
*/
template<typename value_t, typename id_t = long,
         typename value_no_cv_t = typename std::remove_cv<value_t>::type,
		 typename id_no_cv_t = typename std::remove_cv<id_t>::type>
class PiercedIterator
	: public std::iterator<std::forward_iterator_tag, value_no_cv_t, std::ptrdiff_t, value_t*, value_t&>
{
	static_assert(std::is_integral<id_t>::value, "Signed integer required for id.");
	static_assert(std::numeric_limits<id_t>::is_signed, "Signed integer required for id.");

private:
	/*!
		Container.
	*/
	template<typename PV_value_t, typename PV_id_t>
	using Container = PiercedVector<PV_value_t, PV_id_t>;

public:
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

	// Operators
	PiercedIterator& operator++();
	PiercedIterator operator++(int);

	value_t & operator*() const;
	value_t * operator->() const;
	operator PiercedIterator<const value_t, const id_t>() const;

	/*!
		Two-way comparison.
	*/
	template<typename other_value_t, typename other_id_t = long,
         typename other_value_no_cv_t = typename std::remove_cv<value_t>::type,
		 typename other_id_no_cv_t = typename std::remove_cv<id_t>::type>
	bool operator==(const PiercedIterator<other_value_t, other_id_t>& rhs) const
	{
		return (m_container == rhs.m_container) && (m_pos == rhs.m_pos);
	}

	/*!
		Two-way comparison.
	*/
	template<typename other_value_t, typename other_id_t = long,
         typename other_value_no_cv_t = typename std::remove_cv<value_t>::type,
		 typename other_id_no_cv_t = typename std::remove_cv<id_t>::type>
	bool operator!=(const PiercedIterator<other_value_t, other_id_t>& rhs) const
	{
		return (m_container != rhs.m_container) || (m_pos != rhs.m_pos);
	}

private:
	/*!
		Internal pointer to the container.
	*/
	Container<value_no_cv_t, id_no_cv_t> *m_container;

	/*!
		Position inside the container.
	*/
	size_t m_pos;

	// Constructors
	explicit PiercedIterator(Container<value_no_cv_t, id_no_cv_t> *container, const size_t &pos);
	explicit PiercedIterator(const Container<value_no_cv_t, id_no_cv_t> *container, const size_t &pos);

};

/*!
	\ingroup containers

	@brief Metafunction for generating of a pierced vector.

	@details
	Usage: Use <tt>PiercedVector<value_t></tt> to declare a pierced vector.

	PiercedVector can work only with objects that are identified by a
	unique id.

	Internally all the holes are stored in a single vector. The first part
	of this vector contains the "regular" holes, whereas the last part
	contains the "pending" holes. The space reserved to the pending holes
	is fixed. When this space if full, the 'holes_flush' function will be
	called and all pending holes will be converted to regular holes.
	New positions for inserting new elements will be searched first among
	the pending holes and then among the regular holes.

	@tparam value_t The type of the objects stored in the vector
	@tparam id_t The type of the ids to associate to the objects
*/
template<typename value_t, typename id_t = long>
class PiercedVector
{
	static_assert(std::is_integral<id_t>::value, "Signed integer required for id.");
	static_assert(std::numeric_limits<id_t>::is_signed, "Signed integer required for id.");

private:
	/*!
		Maximum number of pending deletes before the changes are flushed.
	*/
	static const std::size_t MAX_PENDING_HOLES;

public:
	// Friendships
	template<typename PI_value_t, typename PI_id_t, typename PI_value_no_cv_t, typename PI_id_no_cv_t>
	friend class PiercedIterator;

	/*!
		Type of data stored in the container
	*/
	typedef value_t value_type;

	/*!
		Type of ids stored in the container
	*/
	typedef id_t id_type;

	/*!
		Iterator for the pierced array.
	*/
	typedef PiercedIterator<value_t> iterator;

	/*!
		Constant iterator for the pierced array.
	*/
	typedef PiercedIterator<const value_t, const long> const_iterator;

	/*!
		Iterator for the pierced array raw container.
	*/
	typedef typename std::vector<value_t>::iterator raw_iterator;

	/*!
		Constant iterator for the pierced array raw container.
	*/
	typedef typename std::vector<value_t>::const_iterator raw_const_iterator;

	/*!
		Functional for compare the position of two elements
	*/
	struct position_less
	{
		position_less(PiercedVector<value_t> &vector)
		{
			m_vector = &vector;
		}

		bool operator()(const id_t &id_1, const id_t &id_2) const
		{
			return m_vector->get_pos_from_id(id_1) < m_vector->get_pos_from_id(id_2);
		}

		PiercedVector<value_t> *m_vector;
	};

	/*!
		Functional for compare the position of two elements
	*/
	struct position_greater
	{
		position_greater(PiercedVector<value_t> &vector)
		{
			m_vector = &vector;
		}

		bool operator()(const id_t &id_1, const id_t &id_2) const
		{
			return m_vector->get_pos_from_id(id_1) > m_vector->get_pos_from_id(id_2);
		}

		PiercedVector<value_t> *m_vector;
	};

	// Contructors
	PiercedVector();
	PiercedVector(std::size_t n);

	// Methods that modify the contents of the container
	iterator push_back(const id_t &id, value_t &&value);

	iterator reclaim(const id_t &id);
	iterator reclaim_after(const id_t &referenceId, const id_t &id);
	iterator reclaim_back(const id_t &id);
	iterator reclaim_before(const id_t &referenceId, const id_t &id);

	iterator move_after(const id_t &referenceId, const id_t &id, bool delayed = false);
	iterator move_before(const id_t &referenceId, const id_t &id, bool delayed = false);

	iterator insert(const id_t &id, value_t &&value);
	iterator insert_after(const id_t &referenceId, const id_t &id, value_t &&value);
	iterator insert_before(const id_t &referenceId, const id_t &id, value_t &&value);

	iterator replace(id_t id, value_t &&value);

	void update_id(const id_t &currentId, const id_t &updatedId);

	template<typename... Args>
	iterator emplace(const id_t &id, Args&&... args);
	template<typename... Args>
	iterator emplace_after(const id_t &referenceId, const id_t &id, Args&&... args);
	template<typename... Args>
	void emplace_back(const id_t &id, Args&&... args);
	template<typename... Args>
	iterator emplace_before(const id_t &referenceId, const id_t &id, Args&&... args);

	template<typename... Args>
	iterator emreplace(id_t id, Args&&... args);

	iterator erase(id_t id, bool delayed = false);

	void pop_back();

	void swap(const id_t &id_first, const id_t &id_second);

	// Methods that modify the container as a whole
	void clear(bool release = true);
	void flush();
	void reserve(std::size_t n);
	void resize(std::size_t n);
	void sort();
	void squeeze();
	void swap(PiercedVector& x) noexcept;

	// Methods that extract information on the container
	std::size_t capacity();
	bool contiguous() const;
	void dump();
	bool empty() const;
	bool is_iterator_slow();
	std::size_t max_size() const;
	std::size_t size() const;

	// Methods that extract information on the contents of the container
	bool exists(id_t id);
	std::size_t extract_flat_index(id_t id) const;

	std::vector<id_t> get_ids(bool ordered = true);
	id_t get_size_marker(const size_t &targetSize, const id_t &fallback = -1);

	// Methods that extract the contents of the container
	value_t * data() noexcept;

	value_t & back();
	const value_t & back() const;

	value_t & front();
	const value_t & front() const;

	value_t & at(const id_t &id);
	const value_t & at(const id_t &id) const;

	value_t & raw_at(const std::size_t &pos);
	const value_t & raw_at(const std::size_t &pos) const;
	std::size_t raw_index(id_t id) const;

	const value_t & operator[](const id_t &id) const;
	value_t & operator[](const id_t &id);

	// Iterators
	iterator get_iterator(const id_t &id) noexcept;
	const_iterator get_const_iterator(const id_t &id) const noexcept;

	iterator begin() noexcept;
	iterator end() noexcept;
	const_iterator begin() const noexcept;
	const_iterator end() const noexcept;
	const_iterator cbegin() const noexcept;
	const_iterator cend() const noexcept;

	raw_iterator raw_begin() noexcept;
	raw_iterator raw_end() noexcept;
	raw_const_iterator raw_begin() const noexcept;
	raw_const_iterator raw_end() const noexcept;
	raw_const_iterator raw_cbegin() const noexcept;
	raw_const_iterator raw_cend() const noexcept;

private:
	/*!
		Hasher for the id map.

		Since the id are uniques, the hasher can be a function that
		takes the id and cast it to a size_t.

		The hasher is defined as a struct, because a struct can be
		passed as an object into metafunctions (meaning that the type
		deduction for the template paramenters can take place, and
		also meaning that inlining is easier for the compiler). A bare
		function would have to be passed as a function pointer.
		To transform a function template into a function pointer,
		the template would have to be manually instantiated (with a
		perhaps unknown type argument).

	*/
	struct PiercedHasher {
		/*!
			Function call operator that casts the specified
			value to a size_t.

			\tparam U type of the value
			\param value is the value to be casted
			\result Returns the value casted to a size_t.
		*/
		template<typename U>
		constexpr std::size_t operator()(U&& value) const noexcept
		{
			return static_cast<std::size_t>(std::forward<U>(value));
		}
	};

	/*!
		Container used for storing holes
	*/
	typedef std::vector<std::size_t> hole_container;

	/*!
		Hole iterator
	*/
	typedef hole_container::iterator hole_iterator;

	/*!
		Vector that will hold the elements.
	*/
	std::vector<value_t>m_v;

	/*!
		Vector that will hold the ids.
	*/
	std::vector<id_t> m_ids;

	/*!
		Container that will hold a list of the holes present in
		the piecrecd vector.
	*/
	hole_container m_holes;

	/*!
		Iterator pointing to the first regular hole
	*/
	hole_iterator m_holes_regular_begin;

	/*!
		Iterator pointing to the last regular hole
	*/
	hole_iterator m_holes_regular_end;

	/*!
		Iterator pointing to the first pending hole
	*/
	hole_iterator m_holes_pending_begin;

	/*!
		Iterator pointing to the last pending hole
	*/
	hole_iterator m_holes_pending_end;

	/*!
		Tracks if the regular holes are sorted
	*/
	bool m_holes_regular_sorted;

	/*!
		Tracks if the pending holes are sorted
	*/
	bool m_holes_pending_sorted;

	/*!
		Map that links the id of the elements and their position
		inside the internal vector.
	*/
	std::unordered_map<id_t, std::size_t, PiercedHasher> m_pos;

	/*!
		Position of the first element in the internal vector.
	*/
	std::size_t m_first_pos;

	/*!
		Position of the last element in the internal vector.
	*/
	std::size_t m_last_pos;

	/*!
		Position of the first dirty element.

		After the first dirty position the id of the holes can not be
		properly defined, meaning that the iterator can take longer to
		iterate through the elements.
	*/
	std::size_t m_first_dirty_pos;

	/*!
		Compares the id of the elements in the specified position.

		\param pos_x is the position to the first element to compare
		\param y is the position to the second element to compare
		\result Returns true if the element x has an id lower than the element
		y, false otherwise. Negative ids are special ids and are considered
		higher than positive ids.
	*/
	struct id_less
	{
		const std::vector<id_t> &m_ids;

		id_less(const std::vector<id_t> &ids)
			: m_ids(ids)
		{
		}

	    inline bool operator() (const std::size_t &pos_x, const std::size_t &pos_y)
	    {
			id_t id_x = m_ids[pos_x];
			id_t id_y = m_ids[pos_y];

		    if (id_x >= 0 && id_y < 0) {
			    return true;
		    } else if (id_x < 0 && id_y >= 0) {
			    return false;
		    } else if (id_x >= 0) {
			    return (id_x < id_y);
		    } else {
			    return (id_x > id_y);
		    }
	    }
	};

	iterator get_iterator_from_pos(const std::size_t &pos) noexcept;
	const_iterator get_const_iterator_from_pos(const std::size_t &pos) const noexcept;

	std::size_t fill_pos(const std::size_t &pos, const id_t &id);
	std::size_t fill_pos_append(const id_t &id);
	std::size_t fill_pos_insert(const std::size_t &pos, const id_t &id);
	std::size_t fill_pos_head(const id_t &id);
	std::size_t fill_pos_tail(const id_t &id);
	std::size_t fill_pos_after(const std::size_t &referencePos, const id_t &id);
	std::size_t fill_pos_before(const std::size_t &referencePos, const id_t &id);

	void pierce_pos(const std::size_t &pos, bool flush = true);

	void holes_clear(bool release = true);
	std::size_t holes_count();
	std::size_t holes_count_pending();
	std::size_t holes_count_regular();
	void holes_flush();
	void holes_clear_pending();
	void holes_clear_pending(const long &offset, const long &nRegulars);
	void holes_sort_pending();
	void holes_sort_regular();

	std::size_t find_prev_used_pos(std::size_t pos);
	std::size_t find_next_used_pos(std::size_t pos);
	bool is_pos_empty(std::size_t pos);
	std::size_t get_pos_from_id(id_t id) const;
    void set_pos_id(const std::size_t &pos, const id_t &id);
	void update_empty_pos_id(const std::size_t &pos, const std::size_t &nextUsedPos);
	void update_first_used_pos(const std::size_t &updated_first_pos);
	void update_last_used_pos(const std::size_t &updated_last_pos);

	std::size_t storage_size() const;
	void storage_resize(size_t n);

    template<typename order_t>
    void reorder_vector(std::vector<size_t>& order, std::vector<order_t>& v, const size_t &size);

};

}

// Include the implementation
#include "piercedVector.tpp"

#endif
