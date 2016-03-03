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

// To check if the provided template argument implements the needed methods,
// the following Stackoverflow questions are used:
//
// http://stackoverflow.com/questions/257288/is-it-possible-to-write-a-c-template-to-check-for-a-functions-existence
// http://stackoverflow.com/questions/16976720/how-to-i-restrict-a-template-class-to-certain-types

namespace utils {

/*!
	\ingroup containerUtils

	Checks if the template parameter has a method called "get_id".
*/
template <typename T>
class has_get_id
{
    typedef char true_type;
    typedef long false_type;

    template <typename C> static true_type test(decltype(&C::get_id));
    template <typename C> static false_type test(...);

public:
    enum { value = (sizeof(test<T>(0)) == sizeof(true_type)) };
};

/*!
	\ingroup containerUtils

	Checks if the template parameter has a method called "set_id".
*/
template <typename T>
class has_set_id
{
    typedef char true_type;
    typedef long false_type;

    template <typename C> static true_type test(decltype(&C::set_id));
    template <typename C> static false_type test(...);

public:
    enum { value = (sizeof(test<T>(0)) == sizeof(true_type)) };
};

}

/*!
	\ingroup containers

	@brief Iterator for the class PiercedVector

	@details
	PiercedIterator can work only with objects that are identified by a
	unique id. The objects must implement a function, called get_id,
	that returns the id of the object. The ids have to be positive, but
	the type that defines the id must be a signed integer (negative
	id are used by PiercedVector to store special information in the
	elements).

	Usage: Use <tt>PiercedVector<T>::iterator</tt> to declare an iterator
	for a pierced vector, use <tt>PiercedVector<Type>::const_iterator</tt> to
	declare a const iterator for a pierced vector.

	@tparam T The type of the objects stored in the vector
*/

template<class T, class unqualified_T = typename std::remove_cv<T>::type>
class PiercedIterator
	: public std::iterator<std::forward_iterator_tag, unqualified_T,
				std::ptrdiff_t, T*, T&>
{
	// PiercedIterator can work only with calsses that are identified by a
	// unique id (i.e., classes that implements get_id)
	static_assert(utils::has_get_id<T>::value, "Provided class does not implement get_id");

private:
	/*!
		Iterator for the internal vector that holds the elements in
		the pierced array.
	*/
	typedef typename std::vector<unqualified_T>::iterator BaseIterator;

	/*!
		Const iterator for the internal vector that holds the elements in
		the pierced array.
	*/
	typedef typename std::vector<unqualified_T>::const_iterator BaseConstIterator;

	/*!
		Type id_type is the type of the ids.

		It is automatically defined as the type returned by the
		get_id function .
	*/
	typedef decltype(std::declval<unqualified_T>().get_id()) id_type;

	/*!
		Special id value that identifies every dummy element past
		the end of the pierced vector.
	*/
	static const id_type SENTINEL_ID;

	/*!
		Internal iterator.
	*/
	T *m_itr;

	/*!
		Creates a new iterator and initializes it with the
		specified pointer.
	*/
	explicit PiercedIterator(T *value)
		: m_itr(value)
	{
	}

public:

	/*!
		Creates a new uninitialized iterator
	*/
	PiercedIterator()
	{
	}

	/*!
		Creates a new iterator and initializes it with the position of
		the base iterator recevied in input.
	*/
	explicit PiercedIterator(BaseIterator iterator)
		: m_itr(&(*iterator))
	{
	}

	/*!
		Creates a new iterator and initializes it with the position of
		the const base iterator recevied in input.
	*/
	explicit PiercedIterator(BaseConstIterator iterator)
		: m_itr(&(*iterator))
	{
	}

	/*!
		Exchanges the values of the current iterator and
		the iterator recevied as argument.

		\param other the iterator to exchange values with
	*/
	void swap(PiercedIterator& other) noexcept
	{
		using std::swap;
		swap(m_itr, other.m_itr);
	}

	/*!
		Pre-increment operator.
	*/
	PiercedIterator& operator++ ()
	{
		size_t delta = 1;
		while (true) {
			m_itr += delta;

			id_type id = m_itr->get_id();
			if (id == SENTINEL_ID || id >= 0) {
				return *this;
			} else {
				delta = - id;
			}
		}

		assert(false);
	}

	/*!
		Post-increment operator.
	*/
	PiercedIterator operator++ (int)
	{
		PiercedIterator tmp(m_itr);

		++(*this);

		return tmp;
	}

	/*!
		Two-way comparison.
	*/
	template<class other_T, class unqualified_other_T = typename std::remove_cv<other_T>::type>
	bool operator == (const PiercedIterator<other_T>& rhs) const
	{
		return m_itr == rhs.m_itr;
	}

	/*!
		Two-way comparison.
	*/
	template<class other_T, class unqualified_other_T = typename std::remove_cv<other_T>::type>
	bool operator != (const PiercedIterator<other_T>& rhs) const
	{
		return m_itr != rhs.m_itr;
	}

	/*!
		Deference operator.

		\result A reference to the element currently pointed to by the
		        iterator.
	*/
	T& operator* () const
	{
		return *m_itr;
	}

	/*!
		Deference operator.

		\result A reference to the element currently pointed to by the
		        iterator.
	*/
	T* operator-> () const
	{
		return m_itr;
	}

	/*!
		Assignment operator.

		\param iterator is the base type iterator that holds the
		                position to be set
		\result The updated iterator.
	*/
	PiercedIterator & operator= (BaseIterator iterator)
	{
		m_itr = &(*iterator);

		return *this;
	}

	/*!
		Converts the iterator to a const_iterator.
	*/
	operator PiercedIterator<const T>() const
	{
		return PiercedIterator<const T>(m_itr);
	}
};

// Definition of static constants of PiercedIterator
template<class T, class unqualified_T>
const typename PiercedIterator<T, unqualified_T>::id_type
	PiercedIterator<T, unqualified_T>::SENTINEL_ID = std::numeric_limits<id_type>::min();

/*!
	\ingroup containers

	@brief Metafunction for generating of a pierced vector.

	@details
	Usage: Use <tt>PiercedVector<T></tt> to declare a pierced vector.

	PiercedVector can work only with objects that are identified by a
	unique id. The objects must implement a function, called get_id,
	that returns the id of the object and a function, called set_id,
	that sets the id of the objects. The ids have to be positive, but
	the type that defines the id must be a signed integer (negative
	id are used by PiercedVector to store special information in the
	elements).

	Internally all the holes are stored in a single vector. The first part
	of this vector contains the "regular" holes, whereas the last part
	contains the "pending" holes. The space reserved to the pending holes
	is fixed. When this space if full, the 'holes_flush' function will be
	called and all pending holes will be converted to regular holes.
	New positions for inserting new elements will be searched first among
	the pending holes and then among the regular holes.

	@tparam T The type of the objects stored in the vector
*/

template <class T>
class PiercedVector
{
	// PiercedVector can work only with calsses that are identified by a
	// unique id (i.e., classes that implements set_id and get_id)
	static_assert(utils::has_get_id<T>::value, "Provided class does not implement get_id");
	static_assert(utils::has_set_id<T>::value, "Provided class does not implement set_id");

private:
	/*!
		Type size_type is an unsigned integral type.
	*/
	typedef std::size_t size_type;

	/*!
		Member type value_type is the type of the elements in the
		container, defined as an alias of the first class template
		parameter (T).
	*/
	typedef T value_type;

	/*!
		Type id_type is the type of the ids.

		It is automatically defined as the type returned by the
		get_id function .
	*/
	typedef decltype(std::declval<T>().get_id()) id_type;

	/*!
		Special id value that identifies every dummy element past
		the end of the pierced vector.
	*/
	static const id_type SENTINEL_ID;

	/*!
		At the end of the piecred vector, after all stored elements,
		there should always be at least one sentinel dummy element.
	*/
	static const size_type REQUIRED_SENTINEL_COUNT;

	/*!
		Maximum number of pending deletes before the changes are flushed.
	*/
	static const size_type MAX_PENDING_HOLES;

	/*!
		Number of usable positions in the vector.

		The vector must contain at least some sentinel values at
		its end. Moreover, since some hole values have a special
		meaning, the number of usable position should be decreased
		to match the maximum number of holes that is possible to
		store.
	*/
	static const size_type USABLE_POS_COUNT;

public:

	/*!
		Iterator for the pierced array.
	*/
	typedef PiercedIterator<value_type> iterator;

	/*!
		Constant iterator for the pierced array.
	*/
	typedef PiercedIterator<const value_type> const_iterator;

	/*!
		Iterator for the pierced array raw container.
	*/
	typedef typename std::vector<T>::iterator raw_iterator;

	/*!
		Constant iterator for the pierced array raw container.
	*/
	typedef typename std::vector<T>::const_iterator raw_const_iterator;

	/*!
		Constructs an empty pierced vector with no elements.
	*/
	PiercedVector()
	{
		clear();
	}

	/*!
		Constructs a pierced vector with a capacity at least enough
		to contain n elements.

		\param n the minimum capacity requested for the vector
	*/
	PiercedVector(size_type n)
	{
		clear();

		reserve(n);
	}

	/*!
		Returns a reference to the element with the specified id. If
		there is no element with the specified id, an exception is
		thrown.

		\param id the id of the element
		\result A reference to the element with the specified id.
	*/
	value_type & at(const id_type &id)
	{
		return (*this)[id];
	}

	/*!
		Returns a constant reference to the element with the specified
		id. If there is no element with the specified id, an exception
		is thrown.

		\param id the id of the element
		\result A constant reference to the element with the specified
		        id.
	*/
	const value_type & at(const id_type &id) const
	{
		return (*this)[id];
	}

	/*!
		Returns a reference to the last element of the vector. If
		the vector is empty, an exception is thrown.

		\result A reference to the last element of the vector.
	*/
	value_type & back()
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_last_pos];
	}

	/*!
		Returns a constant reference to the last element of the vector.
		If the vector is empty, an exception is thrown.

		\result A constant reference to the last element of the vector.
	*/
	const value_type & back() const
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_last_pos];
	}

	/*!
		Returns an iterator pointing to the first element in the
		vector.

		\result An iterator pointing to the first element in the vector.
	*/
	iterator begin()
	{
		if (empty()) {
			return end();
		}

		return iterator(raw_begin() + m_first_pos);
	}

	/*!
		Returns a constant iterator pointing to the first element
		in the vector.

		\result A constant iterator pointing to the first element in
		the vector.
	*/
	const_iterator begin() const noexcept
	{
		return cbegin();
	}

	/*!
		Returns the size of the storage space currently allocated
		for the vector, expressed in terms of elements.

		\result The size of the currently allocated storage capacity
		        in the vector, measured in terms of the number elements
		        it can hold.
	*/
	size_type capacity()
	{
		return m_v.capacity() - REQUIRED_SENTINEL_COUNT;
	}

	/*!
		Returns an conts_iterator pointing to the first element in the
		vector.

		\result A const_iterator pointing to the first element in
		        the vector.
	*/
	const_iterator cbegin() const noexcept
	{
		if (empty()) {
			return cend();
		}

		return const_iterator(m_v.cbegin() + m_first_pos);
	}

	/*!
		Returns an const_iterator referring to the past-the-end element
		in the vector.

		\result A const_iterator referring to the past-the-end element
		        in the vector.
	*/
	const_iterator cend() const noexcept
	{
		return const_iterator(m_v.cbegin() + m_last_pos + 1);
	}

	/*!
		Removes all elements from the vector (which are destroyed),
		leaving the container with a size of 0.

		\param release if it's true the memory hold by the container will
		be released, otherwise the container will be cleared but its
		memory will not be relased
	*/
	void clear(bool release = true)
	{
		// Clear storage
		m_v.clear();
		if (release) {
			std::vector<value_type>().swap(m_v);
		}
		storage_resize(0);

		// Reset first and last counters
		m_first_pos = 0;
		m_last_pos  = 0;

		// Clear holes
		holes_clear(release);

		// Clear position map
		m_pos.clear();
		std::unordered_map<id_type, size_type, PiercedHasher>().swap(m_pos);

		// There are no dirty positions
		m_first_dirty_pos = m_last_pos + 1;
	}

	/*!
		Returns whether the vector is contiguous (i.e. whether it contains
		no holes).

		\result true if the container is contiguous, false otherwise.
	*/
	bool contiguous() const
	{
		return m_holes.empty();
	}

	/*!
		Returns a direct pointer to the memory array used internally
		by the vector to store its owned elements.

		\result A pointer to the first element in the array used
		        internally by the vector.

	*/
	value_type * data() noexcept
	{
		return m_v.data();
	}

	/*!
		Dumps to screen the internal data.
	*/
	void dump()
	{
		std::cout << "----------------[ DUMP ]----------------" << std::endl;

		std::cout << std::endl;
		std::cout << " size: " << size() << std::endl;

		std::cout << " m_holes_regular_begin: " << std::distance(m_holes.begin(), m_holes_regular_begin) << std::endl;
		std::cout << " m_holes_regular_end  : " << std::distance(m_holes.begin(), m_holes_regular_end) << std::endl;

		std::cout << std::endl;
		std::cout << " Regular holes: " << std::endl;
		for (auto k = m_holes_regular_begin; k < m_holes_regular_end; ++k) {
			std::cout << *k << std::endl;
		}

		std::cout << std::endl;
		std::cout << " m_holes_pending_begin: " << std::distance(m_holes.begin(), m_holes_pending_begin) << std::endl;
		std::cout << " m_holes_pending_end  : " << std::distance(m_holes.begin(), m_holes_pending_end) << std::endl;

		std::cout << std::endl;
		std::cout << " Pending holes" << std::endl;
		for (auto k = m_holes_pending_begin; k < m_holes_pending_end; ++k) {
			std::cout << *k << std::endl;
		}

		std::cout << std::endl;
		std::cout << " m_first_pos: " << m_first_pos << std::endl;
		std::cout << " m_last_pos: " <<  m_last_pos << std::endl;
		std::cout << " Stored ids: " << std::endl;
		for (size_t k = 0; k <= m_last_pos; ++k) {
			id_type id = m_v[k].get_id();
			std::cout << id;
			if (exists(id)) {
				std::cout << " ( " << get_pos_from_id(id) << ")";
			} else {
				std::cout << " ( negative id!)";
			}
			std::cout << std::endl;
		}

		std::cout << "----------------------------------------" << std::endl;
	}

	/*!
		The container is extended by inserting a new element. This
		new element is constructed in place using args as the
		arguments for its construction.

		\param args the arguments forwarded to construct the new element
		\result An iterator that points to the the newly inserted
		        element.
	*/
	template <class... Args>
	iterator emplace(Args&&... args)
	{
		size_type pos = fill_pos_head();

		return _emplace(pos, std::forward<Args>(args)...);
	}

	/*!
		The container is extended by inserting a new element. The element
		will have a position that is between the element with the
		specified reference id and the end of the container.

		\param referenceId is the id of the element after which the
		new element will be inserted
		\param args the arguments forwarded to construct the new element
		\result An iterator that points to the newly inserted element.
	*/
	template <class... Args>
	iterator emplace_after(const id_type &referenceId, Args&&... args)
	{
		size_type pos = fill_pos_after(get_pos_from_id(referenceId));

		return _emplace(pos, std::forward<Args>(args)...);
	}

	/*!
		Inserts a new element at the end of the vector, right after
		its current last element. This new element is constructed
		in place using args as the arguments for its construction.

		\param args the arguments forwarded to construct the new element
	*/
	template <class... Args>
	void emplace_back(Args&&... args)
	{
		size_type pos = fill_pos_append();

		_emplace(pos, std::forward<Args>(args)...);
	}

	/*!
		The container is extended by inserting a new element. This new
		element is constructed in place using args as the arguments for
		its construction. The element will have a position that is between
		the beginning of the container and the element with the specified
		reference id.

		\param referenceId is the id of the element before which the
		new element will be inserted
		\param args the arguments forwarded to construct the new element
		\result An iterator that points to the newly inserted element.
	*/
	template <class... Args>
	iterator emplace_before(const id_type &referenceId, Args&&... args)
	{
		size_type pos = fill_pos_before(get_pos_from_id(referenceId));

		return _emplace(pos, std::forward<Args>(args)...);
	}

	/*!
		Returns whether the vector is empty (i.e. whether its size
		is 0).

		\result true if the container size is 0, false otherwise.
	*/
	bool empty() const
	{
		return m_pos.empty();
	}

	/*!
		The element with the specified id is replaced with a new element.
		This new element is constructed in place using args as the
		arguments for its construction.

		\param id is the id of the element that will be replaced
		\param args the arguments forwarded to construct the new element
		\result An iterator that points to the newly inserted element.
	*/
	template <class... Args>
	iterator emreplace(id_type id, Args&&... args)
	{
		// Position
		size_t pos = m_pos.at(id);

		// Id of the element that is currently occupying the position
		id_type id_prev = m_v[pos].get_id();

		// Replace the element
		m_v[pos] = T(std::forward<Args>(args)...);

		// Update the map
		if (id != id_prev) {
			unlink_id(id_prev);
			link_id(id, pos);
		}

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Returns an iterator referring to the past-the-end element
		in the vector.

		\result An iterator referring to the past-the-end element
		        in the vector.
	*/
	iterator end()
	{
		return iterator(raw_begin() + m_last_pos + 1);
	}

	/*!
		Returns a constant iterator referring to the past-the-end
		element in the vector.

		\result A constant iterator referring to the past-the-end
		element in the vector.
	*/
	const_iterator end() const noexcept
	{
		return cend();
	}

	/*!
		Removes from the vector the element with the specified id.
		If the id does not exists the function throws an
		exception.

		Element is not deleted from the internal vector, instead its
		id is changed to mark the position as empty and allow the
		container to reuse that position.

		\param id the id of the element to erase
		\param delayed if true the deletion of the element will
		be delayed until a flush is called
		\result An iterator pointing to the new location of the
		        element that followed the element erased by the
		        function call. This is the container end if the
		        operation erased the last element in the sequence.
	*/
	iterator erase(id_type id, bool delayed = false)
	{
		return _erase(get_pos_from_id(id), delayed);
	}

	/*!
		Checks if a given id exists in the vector.

		\param id the id to look for
		\result Returns true is the given id exists in the vector,
		        otherwise it returns false.
	*/
	bool exists(id_type id)
	{
		return (m_pos.count(id) != 0);
	}

	/*!
		Gets the flat index of the element with the specified id.

		A flat id is the id associated to a numbering scheme that starts
		from the element in the first position of the container and is
		incremented by one for each element in the container. The first
		element will have a flat id equal to 0, the last element will
		have a flat id equal to (nElements - 1).

		If there is no element with the specified id, an exception is
		thrown.

		\param id the id of the element for witch the flat id is requested
		\result The flat index of the element with the specified id.
	*/
	size_type extract_flat_index(id_type id) const
	{
		size_t pos = get_pos_from_id(id);
		size_t flat = pos - m_first_pos;
		if (m_holes.size() > 0) {
			auto holeBound = lower_bound(m_holes.cbegin(), m_holes.cend(), pos);
			size_t nHolesBefore = std::distance(m_holes.cbegin(), holeBound);

			flat -= nHolesBefore;
		}

		return flat;
	}

	/*!
		Flush all pending changes.
	*/
	void flush()
	{
		// Flush pending holes
		holes_flush();
	}

	/*!
		Returns a reference to the first element of the vector. If
		the vector is empty, an exception is thrown.

		\result A reference to the first element of the vector.
	*/
	value_type & front()
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_first_pos];
	}

	/*!
		Returns a constant reference to the first element of the vector.
		If the vector is empty, an exception is thrown.

		\result A constant reference to the first element of the vector.
	*/
	const value_type & front() const
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_first_pos];
	}

	/*!
		Gets a vector containing the ids of the elements stored in
		the vector.

		\param ordered if is true the ids will be sorted in ascending
		\              order, otherwise the ids will be in random
		               order.
		\result A vector with the id of the elements in stored in the
		        vector.
	*/
	std::vector<id_type> get_ids(bool ordered = true)
	{
		std::vector<id_type> ids;
		ids.reserve(size());

		if (ordered) {
			for(auto item : m_pos) {
				typename std::vector<id_type>::iterator itr = lower_bound(ids.begin(), ids.end(), item.first);
				ids.insert(itr, item.first);
			}
		} else {
			for(auto const &value : *this) {
				ids.push_back(value.get_id());
			}
		}

		return ids;
	}

	/*!
		Returns the id of the elmement before which there is the requested
		number of other elements. If this element does not exist the
		fallback value will be returned.

		\param targetSize is the number of elements that needs to be
		contained before the marker
		\param fallback is the fallback value to be returned if the
		marker cannot be found
		\return The id of the elmement before which there is the requested
		number of other elements. If this element does not exist the
		fallback value will be returned.
	*/
	id_type get_size_marker(const size_t &targetSize, const id_type &fallback = SENTINEL_ID)
	{
		// If the size is zero, we return the first element, if the target
		// size is equal to the size minus one we return the last element,
		// if the target size is greater or equal the current container size
		// we return the fallback value.
		if (targetSize == 0) {
			return m_v[m_first_pos].get_id();
		} else if (targetSize == (size() - 1)) {
			return m_v[m_last_pos].get_id();
		} else if (targetSize >= size() || targetSize <= 0) {
			return fallback;
		}

		// Sort the holes
		holes_sort_regular();
		holes_sort_pending();

		// Iterate to find the position before wihch there is the
		// requeste number of element.
		hole_iterator regular_hole_itr = m_holes_regular_end;
		hole_iterator pending_hole_itr = m_holes_pending_end;

		size_type nEmpties  = 0;
		size_type markerPos = targetSize;
		while (true) {
			if (is_pos_empty(markerPos)) {
				markerPos = find_next_used_pos(markerPos - 1);
			}

			// Count the number of holes and pending deletes before the
			// current marker position
			if (regular_hole_itr != m_holes_regular_begin) {
				hole_iterator itr_previous = regular_hole_itr;
				regular_hole_itr = std::upper_bound(m_holes_regular_begin, regular_hole_itr, markerPos, std::greater<size_type>());
				nEmpties += std::distance(regular_hole_itr, itr_previous);
			}

			if (pending_hole_itr != m_holes_pending_begin) {
				hole_iterator itr_previous = pending_hole_itr;
				pending_hole_itr = std::upper_bound(m_holes_pending_begin, pending_hole_itr, markerPos, std::greater<size_type>());
				nEmpties += std::distance(pending_hole_itr, itr_previous);
			}

			// Get the marker size
			//
			// If we have reached the target size we can exit, otherwise
			// we update the marker and we continue iterating
			size_type markerSize = markerPos - nEmpties;
			if (markerSize == targetSize) {
				break;
			} else {
				markerPos += targetSize - markerSize;
			}
		}

		return m_v[markerPos].get_id();
	}

	/*!
		The container is extended by inserting a new element.

		\param value is the value to be copied (or moved) to the
		            inserted elements.
		\result An iterator that points to the newly inserted element.
	*/
	iterator insert(value_type &&value)
	{
		size_type pos = fill_pos_head();

		return _insert(pos, std::move(value));
	}

	/*!
		The container is extended by inserting a new element. The element
		will have a position that is between the element with the
		specified reference id and the end of the container.

		\param value is the value to be copied (or moved) to the
		inserted element
		\param referenceId is the id of the element after which the
		new element will be inserted
		\result An iterator that points to the newly inserted element.
	*/
	iterator insert_after(const id_type &referenceId, value_type &&value)
	{
		size_type pos = fill_pos_after(get_pos_from_id(referenceId));

		return _insert(pos, std::move(value));
	}

	/*!
		The container is extended by inserting a new element. The element
		will have a position that is between the beginning of the
		container and the element with the specified reference id.

		\param value is the value to be copied (or moved) to the
		inserted element
		\param referenceId is the id of the element before which the
		new element will be inserted
		\result An iterator that points to the newly inserted element.
	*/
	iterator insert_before(const id_type &referenceId, value_type &&value)
	{
		size_type pos = fill_pos_before(get_pos_from_id(referenceId));

		return _insert(pos, std::move(value));
	}


	/*!
		Checks if the container is in a state that can slow down the iterator.

		If there are dirty positions in the vector, the ids associated to
		those positions will not point directly to a non-empty element.
		This means that the iterator will require a loop to reach the next
		non-empty position and this can slow down the loop. Calling the
		'flush' function allow to recover the best performances.

		\result Return true if the container is in a state that can slow down
		the iterator, false otherwise.
	*/
	bool is_iterator_slow()
	{
		return (m_first_dirty_pos <= m_last_pos);
	}

	/*!
		Returns the maximum number of elements that the vector can hold.

		This is the maximum potential size the container can reach due
		to known system or library implementation limitations, but the
		container is by no means guaranteed to be able to reach that
		size: it can still fail to allocate storage at any point before
		that size is reached.
	*/
	size_type max_size() const
	{
		return USABLE_POS_COUNT;
	}

	/*!
		Move the specified element after the element with the given
		reference id.

		\param referenceId is the id of the element after which the
		new element will be moved
		\param id is the id of the element that will be moved
		\param delayed if true some changes can remain in a pending state
		until a flush is called
		\result An iterator that points to the moved element.
	*/
	iterator move_after(const id_type &referenceId, const id_type &id, bool delayed = false)
	{
		size_type updatedPos = fill_pos_after(get_pos_from_id(referenceId));
		size_type currentPos = get_pos_from_id(id);

		return _move(currentPos, updatedPos, delayed);
	}

	/*!
		Move the specified element before the element with the given
		reference id.

		\param referenceId is the id of the element before which the
		new element will be moved
		\param id is the id of the element that will be moved
		\param delayed if true some changes can remain in a pending state
		until a flush is called
		\result An iterator that points to the moved element.
	*/
	iterator move_before(const id_type &referenceId, const id_type &id, bool delayed = false)
	{
		size_type updatedPos = fill_pos_before(get_pos_from_id(referenceId));
		size_type currentPos = get_pos_from_id(id);

		return _move(currentPos, updatedPos, delayed);
	}

	/*!
		Removes the last element in the vector, effectively reducing
		the container size by one.

		Element is not deleted from the internal vector, instead its
		id is changed to mark the position as empty and allow the
		container to reuse that position.
	*/
	void pop_back()
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		_erase(m_last_pos, false);
	}

	/*!
		Adds a new element at the end of the vector, after its current
		last element.

		The content of value is copied (or moved) to the new element.

		\param value the value to be copied (or moved) to the new
		             element
		\result An iterator that points to the newly inserted element.
	*/
	iterator push_back(value_type &&value)
	{
		size_type pos = fill_pos_append();

		return _insert(pos, std::move(value));
	}

	/*!
		Returns a reference to the element at the specified position.

		\param pos the position of the element
		\result A reference to the element in the specified position.
	*/
	value_type & raw_at(const size_type &pos)
	{
		return m_v[pos];
	}

	/*!
		Returns a constant reference to the element at the specified
		position.

		\param pos the position of the element
		\result A constant reference to the element in the specified
		        position.
	*/
	const value_type & raw_at(const size_type &pos) const
	{
		return m_v[pos];
	}

	/*!
		Returns an iterator pointing to the first element in the
		raw container.

		\result An iterator pointing to the first element in the raw
		        container.
	*/
	raw_iterator raw_begin()
	{
		return m_v.begin();
	}

	/*!
		Returns a constant iterator pointing to the first element
		in the raw container.

		\result A constant iterator pointing to the first element in
		the raw container.
	*/
	raw_const_iterator raw_begin() const noexcept
	{
		return raw_cbegin();
	}

	/*!
		Returns an conts_iterator pointing to the first element in the
		raw container.

		\result A const_iterator pointing to the first element in
		        the raw container.
	*/
	raw_const_iterator raw_cbegin() const noexcept
	{
		return m_v.cbegin();
	}

	/*!
		Returns an const_iterator referring to the past-the-end element
		in raw container.

		\result A const_iterator referring to the past-the-end element
		        in raw container.
	*/
	raw_const_iterator raw_cend() const noexcept
	{
		return m_v.cend();
	}

	/*!
		Returns an iterator referring to the past-the-end element
		in the raw container.

		\result An iterator referring to the past-the-end element
		        in the raw container.
	*/
	raw_iterator raw_end()
	{
		return m_v.end();
	}

	/*!
		Returns a constant iterator referring to the past-the-end
		element in the raw container.

		\result A constant iterator referring to the past-the-end
		element in the raw container.
	*/
	raw_const_iterator raw_end() const noexcept
	{
		return raw_cend();
	}

	/*!
		Gets the row index of the element with the specified id.

		If there is no element with the specified id, an exception is
		thrown.

		\param id the id of the element for witch the raw id is requested
		\result The row index of the element with the specified id.
	*/
	size_type raw_index(id_type id) const
	{
		return get_pos_from_id(id);
	}

	/*!
		Gets an element from a the first position marked as empty and
		assignes to it the specified id. Except for setting the id,
		the element is not modified. Therefore it will still contain
		the data of the element that was previously occupying the
		position or it will be empty if there was no empty position
		and a new element has been created.

		\param id is the id that will be assigned to the element
	*/
	iterator reclaim(const id_type &id)
	{
		size_type pos = fill_pos_head();

		return _reclaim(pos, id);
	}

	/*!
		Gets an element marked as empty and assignes to it the specified
		id. The element will have a position that is between the element
		with the specified reference id and the end of the container.
		Except for setting the id, the element is not modified. Therefore
		it will still contain the data of the element that was previously
		occupying the position or it will be empty if there was no empty
		position and a new element has been created.
		The container is extended by inserting a new element. The element
		will be inserted .

		\param referenceId is the id of the element after which an
		empty position will be reclaimed
		\param id is the id that will be assigned to the element
		\result An iterator that points to the newly inserted element.
	*/
	iterator reclaim_after(const id_type &referenceId, const id_type &id)
	{
		size_type pos = fill_pos_after(get_pos_from_id(referenceId));

		return _reclaim(pos, id);
	}

	/*!
		Gets an element from the first position marked as empty
		past the last element assignes to it the specified id.
		Except for setting the id, the element is not modified.
		Therefore it will still contain the data of the element
		that was previously occupying the position or it will be
		empty if there was no empty position and a new element
		has been created.

		\param id is the id that will be assigned to the element
	*/
	iterator reclaim_back(const id_type &id)
	{
		size_type pos = fill_pos_append();

		return _reclaim(pos, id);
	}

	/*!
		Gets an element marked as empty and assignes to it the specified
		id. The element will have a position that is between the begin
		of the container and the specified reference id.
		Except for setting the id, the element is not modified. Therefore
		it will still contain the data of the element that was previously
		occupying the position or it will be empty if there was no empty
		position and a new element has been created.
		The container is extended by inserting a new element. The element
		will be inserted .

		\param referenceId is the id of the element before which an
		empty position will be reclaimed
		\param id is the id that will be assigned to the element
		\result An iterator that points to the newly inserted element.
	*/
	iterator reclaim_before(const id_type &referenceId, const id_type &id)
	{
		size_type pos = fill_pos_before(get_pos_from_id(referenceId));

		return _reclaim(pos, id);
	}

	/*!
		The element with the specified id is replaced with a new element.

		\param id is the id of the element that will be replaced
		\param value is the value to be moved to the inserted elements.
		\result An iterator that points to the newly inserted element.
	*/
	iterator replace(id_type id, value_type &&value)
	{
		// Position
		size_t pos = m_pos.at(id);

		// Id of the element that is currently occupying the position
		id_type id_prev = m_v[pos].get_id();

		// Replace the element
		m_v[pos] = std::move(value);

		// Update the map
		if (id != id_prev) {
			unlink_id(id_prev);
			link_id(id, pos);
		}

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Requests that the vector capacity be at least enough to contain
		n elements.

		If n is greater than the current vector capacity, the function
		causes the container to reallocate its storage increasing its
		capacity to n (or greater).

		In all other cases, the function call does not cause a
		reallocation and the vector capacity is not affected.

		\param n the minimum capacity requested for the vector
	*/
	void reserve(size_type n)
	{
		m_v.reserve(n + REQUIRED_SENTINEL_COUNT);
	}

	/*!
		 Resizes the container so that it contains n elements.

		If n is smaller than the current container size, the content
		is reduced to its first n elements, removing those beyond
		(and destroying them).

		If n is greater than the current container size, space is
		reserved in the storage to allow the container to reach the
		requested size.

		If n is also greater than the current container capacity, an
		automatic reallocation of the allocated storage space takes
		place.

		Notice that this function changes the actual content of the
		container by erasing elements from it.

		\param n is the new container size, expressed in number of
		elements.
	*/
	void resize(size_type n)
	{
		// If the size of the vector is already the requested size
		// there is nothing to do.
		if (n == size()) {
			return;
		}

		// A request for a size equal to 0 is equivalent to a clear.
		if (n == 0) {
			clear();
			return;
		}

		// If the requested size is greater that the current size we
		// may need to reserve space in the storage to allow the
		// container to reach the requested size.
		if (n > size()) {
			reserve(n);
			return;
		}

		// If the requested size is smaller that the current size
		// we need to perform a real resize.

		// Flush holes
		holes_flush();

		// Find the id of the last element
		id_type last_id = get_size_marker(n - 1);

		// Find the updated last position
		size_type updated_last_pos = get_pos_from_id(last_id);

		// Delete all ids of the elements beyond the updated position
		// of the last element
		iterator itr(raw_begin() + updated_last_pos);
		itr++;
		while (itr != end()) {
			unlink_id(itr->get_id());
			itr++;
		}

		// Update the last position
		update_last_used_pos(updated_last_pos);
	}

	/*!
		Returns the number of elements in the vector.

		This is the number of actual objects held in the vector,
		which is not necessarily equal to its storage capacity.

		\result The number of elements in the container.
	*/
	size_type size() const
	{
		return m_pos.size();
	}


	/*!
		Sorts the elements of the vector in ascending id order.
	*/
	void sort()
	{
		// Squeeze the container
		squeeze();

		// Sort the elements of the vector
		std::sort(m_v.begin(), m_v.end() - 1, less_than_id());

		// Update positions of the ids
		for (size_type pos = 0; pos <= m_last_pos; pos++) {
			link_id(m_v[pos].get_id(), pos, false);
		}
	}

	/*!
		Requests the container to compact the elements and reduce
		its capacity to fit its size.

		The request is non-binding, and the function can leave
		the vector with a capacity greater than its size.

		This may cause a reallocation, but has no effect on the vector
		size and cannot alter its elements.
	*/
	void squeeze()
	{
		// Flush changes
		flush();

		// Compact the vector
		size_type nHoles = holes_count();
		if (nHoles != 0) {
			// Move the elements
			size_type firstPosToUpdate;
			if (m_first_pos == 0) {
				firstPosToUpdate = *(m_holes_regular_end - 1);
			} else {
				firstPosToUpdate = 0;
			}

			size_type offset = 0;
			for (size_type pos = firstPosToUpdate; pos <= m_last_pos; pos++) {
				if (offset < nHoles && *(m_holes_regular_end - offset - 1) == pos) {
					++offset;
					continue;
				}

				id_type id = m_v[pos].get_id();
				size_type updatedPos = pos - offset;

				m_v[updatedPos] = std::move(m_v[pos]);
				link_id(id, updatedPos, false);
			}

			// Clear the holes
			holes_clear();

			// Reset first and last counters
			update_first_used_pos(0);
			update_last_used_pos(size() - 1);
		}

		// Shrink to fit
		m_v.shrink_to_fit();
	}

	/*!
		Exchanges the content of the container by the content of x,
		which is another vector object of the same type. Sizes may
		differ.

		After the call to this member function, the elements in this
		container are those which were in x before the call, and the
		elements of x are those which were in this. All iterators,
		references and pointers remain valid for the swapped objects.

		\param x Another vector container of the same type (i.e.,
		         instantiated with the same template parameters, T and
		         Alloc) whose content is swapped with that of this
		         container.
	*/
	void swap(PiercedVector& x) noexcept
	{
		std::swap(x.m_first_pos, m_first_pos);
		std::swap(x.m_last_pos, m_last_pos);
		std::swap(x.m_first_dirty_pos, m_first_dirty_pos);
		std::swap(x.m_v, m_v);
		std::swap(x.m_holes, m_holes);
		std::swap(x.m_holes_regular_begin, m_holes_regular_begin);
		std::swap(x.m_holes_regular_end, m_holes_regular_end);
		std::swap(x.m_holes_regular_sorted, m_holes_regular_sorted);
		std::swap(x.m_holes_pending_begin, m_holes_pending_begin);
		std::swap(x.m_holes_pending_end, m_holes_pending_end);
		std::swap(x.m_holes_pending_sorted, m_holes_pending_sorted);
		std::swap(x.m_pos, m_pos);
	}

	/*!
		Swap the elements with the specified id.

		\param id_first is the id of the first element to be swapped
		\param id_second is the id of the second element to be swapped
	*/
	void swap(const id_type &id_first, const id_type &id_second)
	{
		// Positions
		size_t pos_first  = m_pos.at(id_first);
		size_t pos_second = m_pos.at(id_second);

		// Swap the elements
		T tmp = std::move(m_v[pos_first]);
		m_v[pos_first]  = std::move(m_v[pos_second]);
		m_v[pos_second] = std::move(tmp);

		// Relink the ids
		link_id(id_first, pos_second, false);
		link_id(id_second, pos_first, false);
	}

	/*!
		Updates the id of the specified element.

		\param currentId is the current id of the element
		\param updatedId is the new id of the element
	*/
	void update_id(const id_type &currentId, const id_type &updatedId)
	{
		const size_t pos = get_pos_from_id(currentId);
		m_v[pos].set_id(updatedId);
		link_id(updatedId, pos, false);
	}

	/*!
		Returns a constant reference to the element with the
		specified id. If there is no element with the specified id,
		an exception is thrown.

		\param id the id of the element
		\result A constant reference to the element with the specified
		        id.
	*/
	const value_type & operator[](const id_type &id) const
	{
		size_type pos = get_pos_from_id(id);

		return m_v[pos];
	}

	/*!
		Returns a reference to the element with the
		specified id. If there is no element with the specified id,
		an exception is thrown.

		\param id the id of the element
		\result A reference to the element with the specified id.
	*/
	value_type & operator[](const id_type &id)
	{
		size_type pos = get_pos_from_id(id);

		return m_v[pos];
	}

	/*!
		Functional for compare the position of two elements
	*/
	struct position_less
	{
		position_less(PiercedVector<T> &vector)
		{
			m_vector = &vector;
		}

		bool operator()(const id_type &id_1, const id_type &id_2) const
		{
			return m_vector->get_pos_from_id(id_1) < m_vector->get_pos_from_id(id_2);
		}

		PiercedVector<T> *m_vector;
	};


	/*!
		Functional for compare the position of two elements
	*/
	struct position_greater
	{
		position_greater(PiercedVector<T> &vector)
		{
			m_vector = &vector;
		}

		bool operator()(const id_type &id_1, const id_type &id_2) const
		{
			return m_vector->get_pos_from_id(id_1) > m_vector->get_pos_from_id(id_2);
		}

		PiercedVector<T> *m_vector;
	};

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
	typedef std::vector<size_type> hole_container;

	/*!
		Hole iterator
	*/
	typedef hole_container::iterator hole_iterator;

	/*!
		Vector that will hold the elements.
	*/
	std::vector<value_type>m_v;

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
	std::unordered_map<id_type, size_type, PiercedHasher> m_pos;

	/*!
		Position of the first element in the internal vector.
	*/
	size_type m_first_pos;

	/*!
		Position of the last element in the internal vector.
	*/
	size_type m_last_pos;

	/*!
		Position of the first dirty element.

		After the first dirty position the id of the holes can not be
		properly defined, meaning that the iterator can take longer to
		iterate through the elements.
	*/
	size_type m_first_dirty_pos;

	/*!
		The container is extended by inserting a new element in
		the specified hole. This new element is constructed in
		place using args as the arguments for its construction.

		\param pos is the position where the new element will be inserted
		\param args the arguments forwarded to construct the new element
		\result An iterator that points to the newly inserted element.
	*/
	template <class... Args>
	iterator _emplace(const size_t &pos, Args&&... args)
	{
		// Insert the element
		m_v[pos] = T(std::forward<Args>(args)...);

		// Add the id to the map
		link_id(m_v[pos].get_id(), pos);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Removes from the vector the element at the specified position.

		Element is not deleted from the internal vector, instead its
		id is changed to mark the position as empty and allow the
		container to reuse that position.

		\param pos the position of the element to erase
		\result An iterator pointing to the new location of the
		        element that followed the element erased by the
		        function call. This is the container end if the
		        operation erased the last element in the sequence.
	*/
	iterator _erase(size_type pos, bool delayed = false)
	{
		// Delete id from map
		unlink_id(m_v[pos].get_id());

		// Push the position into the list of holes
		pierce_pos(pos, !delayed);

		// Return the iterator to the element following the one erased
		iterator itr;
		if (empty() || pos >= m_last_pos) {
			itr = end();
		} else {
			itr = raw_begin() + find_next_used_pos(pos);
		}

		return itr;
	}

	/*!
		The container is extended by inserting a new element in
		the specified hole.

		\param pos is the position where the new element will be inserted
		\param value is the value to be copied (or moved) to the
		inserted elements.
		\result An iterator that points to the the newly inserted
		element.
	*/
	iterator _insert(const size_t &pos, value_type &&value)
	{
		// Insert the element
		m_v[pos] = std::move(value);

		// Add the id to the map
		link_id(m_v[pos].get_id(), pos);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Move the element in the specified position.

		\param currentPos is the current position of the element
		\param updatedPos is the new position of the element
		\param delayed if true some changes can remain in a pending state
		until a flush is called
		\result An iterator that points to the moved element.
	*/
	iterator _move(const size_t &currentPos, const size_t &updatedPos, bool delayed = false)
	{
		// Move the element
		//
		// Current position is reset to avoid leaving elements in an
		// inconsistent state
		m_v[updatedPos] = std::move(m_v[currentPos]);

		// Update the map
		link_id(m_v[updatedPos].get_id(), updatedPos, false);

		// Pierce the position
		pierce_pos(currentPos, !delayed);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + updatedPos;

		return itr;
	}

	/*!
		Gets an element from a position marked as empty and assignes
		to it the specified id. Except for setting the id, the element
		is not modified. Therefore it will still contain
		the data of the element that was previously occupying the
		position or it will be empty if there was no empty position
		and a new element has been created.

		\param pos is the position where the new element will be inserted
		\param id is the id that will be assigned to the element
		\result An iterator that points to the the reclaimed element.
	*/
	iterator _reclaim(const size_t &pos, const id_type &id)
	{
		// Set the id of the element
		m_v[pos].set_id(id);

		// Add the id to the map
		link_id(id, pos);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Gets the position in the storage vector of the element with the
		specified id.

		If there is no element with the specified id, an exception is
		thrown.

		\param id the id to look for
		\result The position in the storage vector of the element with
		the specified id.
	*/
	size_type get_pos_from_id(id_type id) const
	{
		return m_pos.at(id);
	}

	/*!
		Gets a position for storing a new element.

		The position as alyaws appeneded to the end of the container.

		\result A position for storing a new element.
	*/
	size_type fill_pos_append()
	{
		// Extend the container
		size_type updated_last_pos = m_last_pos;
		if (!empty()) {
			++updated_last_pos;
		}
		update_last_used_pos(updated_last_pos);

		return m_last_pos;
	}

	/*!
		Gets a position for storing a new element.

		The specified position is made available.

		\param pos is the position the will be make available
		\result A position for storing a new element.
	*/
	size_type fill_pos_specific(const size_type &pos)
	{
		// Extend the container
		//
		// If the requester position is beyond the positions of the vector
		// just return the newly appended position
		size_type appendedPos = fill_pos_append();
		if (pos >= m_last_pos) {
			return appendedPos;
		}

		// Shift the elements after the reference position
		for (size_t i = m_last_pos; i > pos; --i) {
			m_v[i] = std::move(m_v[i - 1]);

			id_type id = m_v[i].get_id();
			if (id > 0) {
				link_id(id, i, false);
			}
		}

		// Reset the, now empty, element
		//
		// We need to avoid that this element could be in an inconsistent state.
		reset_pos(pos);

		// Update the regular holes
		if (m_holes_regular_begin != m_holes_regular_end) {
			hole_iterator change_begin = m_holes_regular_begin;
			hole_iterator change_end   = upper_bound(m_holes_regular_begin, m_holes_regular_end, pos, std::greater<size_type>());
			for (auto itr = change_begin; itr != change_end; itr++) {
				(*itr)++;
			}
		}

		// Update the pending holes
		if (m_holes_pending_begin != m_holes_pending_end) {
			hole_iterator change_begin = m_holes_pending_begin;
			hole_iterator change_end   = upper_bound(m_holes_pending_begin, m_holes_pending_end, pos, std::greater<size_type>());
			for (auto itr = change_begin; itr != change_end; itr++) {
				(*itr)++;
			}
		}

		return pos;
	}

	/*!
		Gets a position for storing a new element.

		The first available position starting from the head of the container
		is returned.

		\result A position for storing a new element.
	*/
	size_type fill_pos_head()
	{
		// If there are holes we can fill a hole.
		long nRegulars = holes_count_regular();
		long nPendings = holes_count_pending();
		long nHoles    = nRegulars + nPendings;
		if (nHoles != 0) {
			// Sort the holes
			if (m_holes_pending_begin != m_holes_pending_end) {
				holes_sort_pending();
			} else {
				holes_sort_regular();
			}

			// The last element of the hole's container is the hole we need to
			// use if we are filling from the head.
			size_type pos = m_holes.back();
			m_holes.pop_back();

			// Update the iterators
			if (nHoles == 1) {
				holes_clear();
			} else {
				if (nRegulars >= 1 || nPendings == 1) {
					m_holes_regular_end   = m_holes.end();
					m_holes_pending_begin = m_holes_regular_end;
				}
				m_holes_pending_end = m_holes.end();
			}

			// Update first position counter
			//
			// If the vector contains a hole, this means that is not empty
			// and that the hole is before the last element, therefore
			// only the first position counter may have changed.
			if (pos < m_first_pos) {
				update_first_used_pos(pos);
			}

			// Return the position filled
			return pos;
		}

		// There are no holes nor pending delete: use an append fill.
		return fill_pos_append();
	}

	/*!
		Gets a position for storing a new element.

		The first available position starting from the tail of the container
		is returned.

		\result A position for storing a new element.
	*/
	size_type fill_pos_tail()
	{
		// If there are holes we can fill a hole.
		long nRegulars = holes_count_regular();
		long nPendings = holes_count_pending();
		long nHoles    = nRegulars + nPendings;
		if (nHoles != 0) {
			// First search among pendings and the among regulars
			size_type pos;
			if (m_holes_pending_begin != m_holes_pending_end) {
				holes_sort_pending();
				pos = *m_holes_pending_begin;
			} else {
				holes_sort_regular();
				pos = *m_holes_regular_begin;
			}

			// Update the iterators
			if (nHoles == 1) {
				holes_clear();
			} else if (nPendings == 1) {
				holes_clear_pending();
			} else if (nPendings > 1) {
				++m_holes_pending_begin;
			} else {
				++m_holes_regular_begin;
			}

			// Update first position counter
			//
			// It is not possible that a hole is past the last element of the
			// vector. We should olny consider that case where a holes is
			// below the first element.
			if (pos < m_first_pos) {
				update_first_used_pos(pos);
			}

			// If previos element is a hole, its id need to be udated
			if (pos > 0 && is_pos_empty(pos - 1)) {
				update_empty_pos_id(pos - 1, pos);
			}

			// Return the position filled
			return pos;
		}

		// There are no holes nor pending delete: use an append fill.
		return fill_pos_append();
	}

	/*!
		Gets a position for storing a new element.

		The first available position after the specified position is
		returned.

		\param referencePos is the position of the element after which the
		new available position will be searched for
		\result A position for storing a new element.
	*/
	size_type fill_pos_after(const size_type &referencePos)
	{
		// Check if we can fill a hole
		//
		// The last hole should hava a position higher than the reference
		// position
		if (holes_count() != 0) {
			// First search among pendings and the among regulars
			size_type lastHole;
			if (m_holes_pending_begin != m_holes_pending_end) {
				holes_sort_pending();
				lastHole = *(m_holes_pending_begin);
			} else {
				holes_sort_regular();
				lastHole = *(m_holes_regular_begin);
			}

			if (lastHole > referencePos) {
				return fill_pos_tail();
			}
		}

		// We have to append the element at the end of the vector
		return fill_pos_append();
	}

	/*!
		Gets a position for storing a new element.

		The first available position before the specified position is
		returned.

		\param referencePos is the position of the element before which the
		new available position will be searched for
		\result A position for storing a new element.
	*/
	size_type fill_pos_before(const size_type &referencePos)
	{
		// Check if we can fill a hole
		//
		// The first hole available should be in a position lower than the
		// reference position.
		if (holes_count() != 0) {
			// First search among pendings and the among regulars
			size_type firstHole;
			if (m_holes_pending_begin != m_holes_pending_end) {
				holes_sort_pending();
				firstHole = *(m_holes_pending_end - 1);
			} else {
				holes_sort_regular();
				firstHole = *(m_holes_regular_end - 1);
			}

			if (firstHole < referencePos) {
				return fill_pos_head();
			}
		}

		// We have to insert the element at the specified position
		return fill_pos_specific(referencePos);
	}

	/*!
		Returns the first non-empty position before the specified
		starting position.

		If the starting position is the first posistion, an
		exception is thrown.

		\param pos starting position
		\result The firt non-empty position before the starting
		        position.
	*/
	size_type find_prev_used_pos(size_type pos)
	{
		size_type prev_pos = pos;
		while (true) {
			if (prev_pos == m_first_pos) {
				throw std::out_of_range ("Already in the firts position");
			}
			prev_pos--;

			id_type prev_id = m_v[prev_pos].get_id();
			if (prev_id >= 0) {
				return prev_pos;
			}
		}
	}

	/*!
		Returns the first non-empty position after the specified
		starting position.

		If the starting position is the last posistion, an
		exception is thrown.

		\param pos starting position
		\result The firt non-empty position after the starting
		        position.
	*/
	size_type find_next_used_pos(size_type pos)
	{
		size_type next_pos   = pos;
		size_type next_delta = 1;
		while (true) {
			if (next_pos == m_last_pos) {
				throw std::out_of_range ("Already in the last position");
			}
			next_pos += next_delta;

			id_type next_id = m_v[next_pos].get_id();
			if (next_id >= 0) {
				return next_pos;
			} else {
				next_delta = - next_id;
			}
		}
	}

	/*!
		Clear the list of available holes.
	*/
	void holes_clear(bool release = true)
	{
		m_holes.clear();
		if (release) {
			hole_container().swap(m_holes);
		}
		holes_clear_pending(0, 0);
		m_holes_regular_sorted = true;
	}

	/*!
		Count the available holes.

		\result The number of available holes.
	*/
	size_type holes_count()
	{
		return holes_count_pending() + holes_count_regular();
	}

	/*!
		Count the pending holes.

		\result The number of pending holes.
	*/
	size_type holes_count_pending()
	{
		return std::distance(m_holes_pending_begin, m_holes_pending_end);
	}

	/*!
		Count the regular holes.

		\result The number of regular holes.
	*/
	size_type holes_count_regular()
	{
		return std::distance(m_holes_regular_begin, m_holes_regular_end);
	}

	/*!
		Flushes the list of available holes.

		All the pending hole are converted to regular holes and new
		space is reserved for future pending holes.
	*/
	void holes_flush()
	{
		// If there are no pending holes there is nothing to do
		if (m_holes_pending_begin == m_holes_pending_end) {
			return;
		}

		// Convert pending holes to regular ones
		for (auto itr = m_holes_pending_begin; itr != m_holes_pending_end; ++itr) {
			const size_type &pos = *itr;

			// Update the id of the element in the specified position
			update_empty_pos_id(pos);

			// Move the pending holes into the list of regular holes
			//
			// If there is space available at the beginning of the holes, try
			// using pending holes to fill that gap.
			if (m_holes_regular_begin != m_holes.begin()) {
				--m_holes_regular_begin;
				*m_holes_regular_begin = pos;

				// Regular holes are no more sorted
				if (m_holes_regular_sorted) {
					m_holes_regular_sorted = false;
				}
			} else {
				if (itr != m_holes_regular_end) {
					*m_holes_regular_end = pos;
				}
				++m_holes_regular_end;
				++m_holes_pending_begin;

				// Check if regular holes are still sorted
				if (m_holes_regular_sorted) {
					size_type nRegulars = holes_count_regular();
					if (nRegulars > 1 && (*(m_holes_regular_end - 1) > *(m_holes_regular_end - 2))) {
						m_holes_regular_sorted = false;
					}
				}
			}
		}

		// Move the holes at the beginning of the vector
		//
		// The iterators will be updated when clearing the pending holes.
		size_type nRegulars = holes_count_regular();
		if (m_holes_regular_begin != m_holes.begin()) {
			size_type offset = std::distance(m_holes.begin(), m_holes_regular_begin);
			for (size_type k = 0; k < nRegulars; ++k) {
				m_holes[k] = m_holes[k + offset];
			}
		}

		// Resize the vector
		holes_clear_pending(0, nRegulars);

		// There are no more dirty positions
		m_first_dirty_pos = m_last_pos + 1;
	}

	/*!
		Reset pending holes
	*/
	void holes_clear_pending()
	{
		long offset    = std::distance(m_holes.begin(), m_holes_regular_begin);
		long nRegulars = holes_count_regular();

		holes_clear_pending(offset, nRegulars);
	}

	/*!
		Reset pending holes

		\param offset is the distance between the first regulat hole and the
		begin of the hole's container
		\param nRegulars is the number of regulars holes in the hole's
		container
	*/
	void holes_clear_pending(const long &offset, const long &nRegulars)
	{
		m_holes.reserve(offset + nRegulars + MAX_PENDING_HOLES);
		m_holes.resize(offset + nRegulars);

		m_holes_regular_begin = m_holes.begin() + offset;
		m_holes_regular_end   = m_holes_regular_begin + nRegulars;
		m_holes_pending_begin = m_holes_regular_end;
		m_holes_pending_end   = m_holes_pending_begin;
	}

	/*!
		Sort the list of pending holes in descendent order
	*/
	void holes_sort_pending()
	{
		if (m_holes_pending_sorted) {
			return;
		}

		std::sort(m_holes_pending_begin, m_holes_pending_end, std::greater<size_type>());
		m_holes_pending_sorted = true;
	}

	/*!
		Sort the list of regular holes in descendent order
	*/
	void holes_sort_regular()
	{
		if (m_holes_regular_sorted) {
			return;
		}

		std::sort(m_holes_regular_begin, m_holes_regular_end, std::greater<size_type>());
		m_holes_regular_sorted = true;
	}

	/*!
		Returns if the specified position is empty.

		A position is considered empty if the element in that
		position has an id less than 0.

		\param pos the position to check
		\result true is the position is empty, false otherwise.
	*/
	bool is_pos_empty(size_type pos)
	{
		return (m_v[pos].get_id() < 0);
	}

	/*!
		Updates the position associated with the specified id.

		\param id is the id of the element
		\param pos is the position that will be associated with
		the id
		\param checkUnique specifies if a check of the uniqueness
		of the id will be performed
	*/
	void link_id(const id_type id, const size_t pos, bool checkUnique = true)
	{
		// Check uniqueness of the id
		if (checkUnique && exists(id)) {
			throw std::out_of_range ("Duplicate id");
		}

		// Ids needs to be positive
		if (id < 0) {
			throw std::out_of_range ("Negative id");
		}

		// Add id to the map
		m_pos[id] = pos;
	}

	/*!
		Marks a position as empty.

		The position is inserted in the list of holes. If the list of pending
		holes is full, a flush is called before adding the hole. This means
		that the hole is always added as a pending hole.

		\param hole is the position of the new hole
	*/
	void pierce_pos(const size_type &pos, bool flush = true)
	{
		// If removing the last position, there is no need to add the
		// position to the holes, it's enough to update the last position
		// counter or clear the container if this was the last hole.
		if (pos == m_last_pos) {
			if (empty()) {
				clear();
			} else {
				size_type updated_last_pos = find_prev_used_pos(m_last_pos);
				update_last_used_pos(updated_last_pos);
			}
			return;
		}

		// Reset the position
		reset_pos(pos);

		// If removing the first position, update the counter
		if (pos == m_first_pos) {
			size_type updated_first_pos = find_next_used_pos(m_first_pos);
			update_first_used_pos(updated_first_pos);
		}

		// If the list of pending holes is full, flush the holes.
		if (m_holes.size() == m_holes.capacity()) {
			holes_flush();
		}

		// Add the hole at the end of the pending holes
		m_holes.push_back(pos);
		m_holes_pending_end = m_holes.end();

		// Check if pending holes are sorted
		if (m_holes_pending_sorted) {
			size_t nPendings = holes_count_pending();
			if (nPendings > 1 && (*(m_holes_pending_end - 1) > *(m_holes_pending_end - 2))) {
				m_holes_pending_sorted = false;
			}
		}

		// Flush
		if (flush) {
			holes_flush();
		}
	}

	/*!
		Reset the element in the specified position.

		\param pos is the position that will be reset
	*/
	void reset_pos(const size_type &pos)
	{
		m_v[pos] = T();
		update_empty_pos_id(pos, false);
		m_first_dirty_pos = std::min(pos, m_first_dirty_pos);
	}

	/*!
		Returns the size of the storage expressed in terms of
		elements that the storage contains (including the holes).

		\result The size of the storage.
	*/
	size_type storage_size() const
	{
		return m_v.size() - REQUIRED_SENTINEL_COUNT;
	}

	/*!
		Resize the storage.

		In order for the iterator to correctly identify the last
		non-empty position, the container needs to store, after all
		the elements, at least one sentinel. A sentinel element is a
		dummy element with the special id SENTINEL_ID. All elements
		after the last non-empty position are sentinel elements.

		\param n is the new container size, expressed in number of
		elements.
	*/
	void storage_resize(size_t n)
	{
		size_t previous_raw_size = m_v.size();
		if (n == previous_raw_size - REQUIRED_SENTINEL_COUNT + 1) {
			m_v.emplace_back();
			T &sentinel = m_v.back();
			sentinel.set_id(SENTINEL_ID);
		} else {
			m_v.resize(n + REQUIRED_SENTINEL_COUNT);

			size_t current_raw_size = m_v.size();
			for (size_t k = std::min(n, previous_raw_size); k < current_raw_size; ++k) {
				m_v[k].set_id(SENTINEL_ID);
			}
		}
	}

	/*!
		Updates the id of the element in the specified position to make
		it an empty element. If needed, updates also of the id of the
		element before the specified position. The function will figure
		out the position of the next non-empty element.

		\param pos is the position to update
		\param recursive controls if the
	*/
	void update_empty_pos_id(const size_type &pos, bool recursive = true)
	{
		// Position of the next non-empty element
		size_type nextUsedPos;
		if (pos >= m_last_pos) {
			nextUsedPos = m_last_pos;
		} else {
			nextUsedPos = find_next_used_pos(pos);
		}

		// Update the id
		update_empty_pos_id(pos, nextUsedPos, recursive);
	}

	/*!
		Updates the id of the element in the specified position to make
		it an empty element. If needed, updates also of the id of the
		element before the specified position.

		The id of an empty element contains the distance, measured in
		number of elements, between the current element and the next
		non-empty element (the distance is negative). The id of an
		element past the last non-empty position is set to the special
		value SENTINEL_ID.

		\param pos is the position to update
		\param nextUsedPos is the position of the next non-empty element
	*/
	void update_empty_pos_id(const size_type &pos, const size_type &nextUsedPos, bool recursive = true)
	{
		// Id of the element
		id_type id;
		if (nextUsedPos <= pos) {
			id = SENTINEL_ID;
		} else {
			id = pos - nextUsedPos;
		}

		// Update the id of the element in the current position
		m_v[pos].set_id(id);

		if (!recursive) {
			return;
		}

		// Update the id of the elements in previous positions
		if (pos > 0) {
			size_type prevPos = pos - 1;
			while (is_pos_empty(prevPos)) {
				if (id != SENTINEL_ID) {
					id--;
				}
				m_v[prevPos].set_id(id);

				if (prevPos > 0) {
					--prevPos;
				} else {
					break;
				}
			}
		}
	}

	/*!
		Update the first used position.
	*/
	void update_first_used_pos(const size_type &updated_first_pos)
	{
		m_first_pos = updated_first_pos;
	}

	/*!
		Update the last used position.
	*/
	void update_last_used_pos(const size_type &updated_last_pos)
	{
		// Hole needs to be updated only if last position has been decrease
		bool update_holes = (holes_count() > 0) && (updated_last_pos < m_last_pos);

		// Update the last position
		m_last_pos = updated_last_pos;

		// Resize the vector
		storage_resize(m_last_pos + 1);

		// If we don't need to update the holes we can exit now
		if (!update_holes) {
			return;
		}

		// Remove regular holes beyond the updated last position
		holes_sort_regular();
		m_holes_regular_begin = std::lower_bound(m_holes_regular_begin, m_holes_regular_end, m_last_pos, std::greater<size_type>());
		if (m_holes_regular_begin == m_holes_regular_end) {
			m_holes_regular_begin = m_holes.begin();
			m_holes_regular_end   = m_holes_regular_begin;
		}

		// Remove pending holes beyond the updated last position
		holes_sort_pending();
		m_holes_pending_begin = std::lower_bound(m_holes_pending_begin, m_holes_pending_end, m_last_pos, std::greater<size_type>());
		if (m_holes_pending_begin == m_holes_pending_end) {
			m_holes_pending_begin = m_holes_regular_end;
			m_holes_pending_end   = m_holes_pending_begin;
		}

		// Resize the hole's container
		if (m_holes_pending_end != m_holes.end()) {
			m_holes.resize(std::distance(m_holes.begin(), m_holes_pending_end));
		}
	}

	/*!
		Removes the specified id from the map.

		\param id is the id that will be removed from the list
	*/
	void unlink_id(const id_type id)
	{
		m_pos.erase(id);
	}

	/*!
		Compares the id of the specified values.

		\param x first values to compare
		\param y second values to compare
		\result Returns true if the x has an id lower that y, false
		        otherwise. Negative ids are special ids and are
		        considered higher than positive ids.
	*/
	struct less_than_id
	{
	    inline bool operator() (const T &x, const T &y)
	    {
			id_type id_x = x.get_id();
			id_type id_y = y.get_id();

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
};

// Definition of static constants of PiercedVector
template<class T>
const typename PiercedVector<T>::id_type
	PiercedVector<T>::SENTINEL_ID = std::numeric_limits<id_type>::min();

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::MAX_PENDING_HOLES = 16384;

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::REQUIRED_SENTINEL_COUNT = 1;

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::USABLE_POS_COUNT = std::numeric_limits<size_type>::max() - REQUIRED_SENTINEL_COUNT;

}

#endif
