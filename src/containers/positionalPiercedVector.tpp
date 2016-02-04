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
#ifndef __BITPIT_POSITIONAL_PIERCED_VECTOR_TPP__
#define __BITPIT_POSITIONAL_PIERCED_VECTOR_TPP__

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

#include "piercedVector.tpp"


namespace bitpit{

/*!
	@ingroup containers

	@brief Iterator for the class PositionalPiercedVector

	@details
	PositionalPiercedIterator can work only with objects that are identified by a
	unique id. The objects must implement a function, called get_id,
	that returns the id of the object. The ids have to be positive, but
	the type that defines the id must be a signed integer (negative
	id are used by PositionalPiercedVector to store special information in the
	elements).

	Usage: Use <tt>PositionalPiercedVector<T>::iterator</tt> to declare an iterator
	for a pierced vector, use <tt>PositionalPiercedVector<Type>::const_iterator</tt> to
	declare a const iterator for a pierced vector.

	@tparam T The type of the objects stored in the vector
*/

template<class T, class unqualified_T = typename std::remove_cv<T>::type>
class PositionalPiercedIterator
	: public std::iterator<std::forward_iterator_tag, unqualified_T,
				std::ptrdiff_t, T*, T&>
{
	// PositionalPiercedIterator can work only with calsses that are identified by a
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
	explicit PositionalPiercedIterator(T *value)
		: m_itr(value)
	{
	}

public:

	/*!
		Creates a new uninitialized iterator
	*/
	PositionalPiercedIterator()
	{
	}

	/*!
		Creates a new iterator and initializes it with the position of
		the base iterator recevied in input.
	*/
	explicit PositionalPiercedIterator(BaseIterator iterator)
		: m_itr(&(*iterator))
	{
	}

	/*!
		Creates a new iterator and initializes it with the position of
		the const base iterator recevied in input.
	*/
	explicit PositionalPiercedIterator(BaseConstIterator iterator)
		: m_itr(&(*iterator))
	{
	}

	/*!
		Exchanges the values of the current iterator and
		the iterator recevied as argument.

		\param other the iterator to exchange values with
	*/
	void swap(PositionalPiercedIterator& other) noexcept
	{
		using std::swap;
		swap(m_itr, other.m_itr);
	}

	/*!
		Pre-increment operator.
	*/
	PositionalPiercedIterator& operator++ ()
	{
		m_itr++;

		id_type id = m_itr->get_id();
		if (id != SENTINEL_ID && id < 0) {
			m_itr += - m_itr->get_id();
		}

		return *this;
	}

	/*!
		Post-increment operator.
	*/
	PositionalPiercedIterator operator++ (int)
	{
		PositionalPiercedIterator tmp(m_itr);

		++(*this);

		return tmp;
	}

	/*!
		Two-way comparison.
	*/
	template<class other_T, class unqualified_other_T = typename std::remove_cv<other_T>::type>
	bool operator == (const PositionalPiercedIterator<other_T>& rhs) const
	{
		return m_itr == rhs.m_itr;
	}

	/*!
		Two-way comparison.
	*/
	template<class other_T, class unqualified_other_T = typename std::remove_cv<other_T>::type>
	bool operator != (const PositionalPiercedIterator<other_T>& rhs) const
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
	PositionalPiercedIterator & operator= (BaseIterator iterator)
	{
		m_itr = &(*iterator);

		return *this;
	}

	/*!
		Converts the iterator to a const_iterator.
	*/
	operator PositionalPiercedIterator<const T>() const
	{
		return PositionalPiercedIterator<const T>(m_itr);
	}
};

// Definition of static constants of PositionalPiercedIterator
template<class T, class unqualified_T>
const typename PositionalPiercedIterator<T, unqualified_T>::id_type
	PositionalPiercedIterator<T, unqualified_T>::SENTINEL_ID = std::numeric_limits<id_type>::min();

/*!
	@ingroup containers

	@brief Metafunction for generating of a pierced vector.

	@details
	Usage: Use <tt>PositionalPiercedVector<T></tt> to declare a positional
	pierced vector. This is a vector that can handle holes and that
	automatically assigns to the contained elements an id which is
	the position of the element in the container.

	PositionalPiercedVector can work only with objects that are identified by a
	unique id. The objects must implement a function, called get_id,
	that returns the id of the object and a function, called set_id,
	that sets the id of the objects. The ids have to be positive, but
	the type that defines the id must be a signed integer (negative
	id are used by PositionalPiercedVector to store special information in the
	elements).

	@tparam T The type of the objects stored in the vector
*/

template <class T>
class PositionalPiercedVector
{
	// PositionalPiercedVector can work only with calsses that are identified by a
	// unique id (i.e., classes that implements set_id and get_id)
	static_assert(utils::has_get_id<T>::value, "Provided class does not implement get_id");
	static_assert(utils::has_set_id<T>::value, "Provided class does not implement set_id");

private:
	enum FillType {
		FILL_APPEND,
		FILL_FIRST
	};

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
	typedef PositionalPiercedIterator<value_type> iterator;

	/*!
		Constant iterator for the pierced array.
	*/
	typedef PositionalPiercedIterator<const value_type> const_iterator;

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
	PositionalPiercedVector()
	{
		clear();
	}

	/*!
		Move constructor.

		\param other vector with the new contents
	*/
	PositionalPiercedVector(PositionalPiercedVector&& other) = default;

	/*!
		Constructs a pierced vector with a capacity at least enough
		to contain n elements.

		\param n the minimum capacity requested for the vector
	*/
	PositionalPiercedVector(size_type n)
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
	*/
	void clear()
	{
		// Clear storage
		m_v.clear();
		std::vector<value_type>().swap(m_v);
		storage_resize(0);

		// Reset first and last counters
		m_first_pos = 0;
		m_last_pos  = 0;

		// Clear holes
		m_holes.clear();
		std::deque<size_type>().swap(m_holes);

		// Clear pending changes
		m_pending_deletes.clear();
		std::deque<size_type>().swap(m_pending_deletes);

		// Clear dirty flag
		m_dirty = false;
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
		return _emplace(FILL_FIRST, std::forward<Args>(args)...);
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
		_emplace(FILL_APPEND, std::forward<Args>(args)...);
	}

	/*!
		Returns whether the vector is empty (i.e. whether its size
		is 0).

		\result true if the container size is 0, false otherwise.
	*/
	bool empty() const
	{
		return (size() == 0);
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
		Flush all pending changes.
	*/
	void flush()
	{
		if (!m_dirty) {
			return;
		}

		// Flush pending deletes
		auto pendingIter = m_pending_deletes.rbegin();
		while (pendingIter != m_pending_deletes.rend()) {
			pierce_pos(*pendingIter);
			++pendingIter;
		}
		std::deque<size_type>().swap(m_pending_deletes);

		// Done
		m_dirty = false;
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
		return _erase(id, delayed);
	}

	/*!
		Checks if a given id exists in the vector.

		\param id the id to look for
		\result Returns true is the given id exists in the vector,
		        otherwise it returns false.
	*/
	bool exists(id_type id)
	{
		return !is_pos_empty(id);
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
			for (auto const &value : *this) {
				const id_type id = value.get_id();
				typename std::vector<id_type>::iterator itr = lower_bound(ids.begin(), ids.end(), id);
				ids.insert(itr, id);
			}
		} else {
			for (auto const &value : *this) {
				ids.push_back(value.get_id());
			}
		}

		return ids;
	}

	/*!
		The container is extended by inserting a new element.

		\param value is the value to be copied (or moved) to the
		            inserted elements.
	*/
	iterator insert(value_type &&value)
	{
		return _insert(FILL_FIRST, std::move(value));
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

		while (!m_pending_deletes.empty()) {
			if (m_pending_deletes.back() < m_last_pos) {
				break;
			}

			size_type pos = pending_deletes_pop_back();
			if (pos != m_last_pos) {
				_erase(pos, false);
			}
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
		return _insert(FILL_APPEND, std::move(value));
	}

	/*!
		Returns a reference to the element at the specified position.

		\param pos the position of the element
		\result A reference to the element in the specified position.
	*/
	value_type & raw_at(const id_type &pos)
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
	const value_type & raw_at(const id_type &pos) const
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
		Gets the row index of the element with the specified id.

		If there is no element with the specified id, an exception is
		thrown.

		\param id the id of the element for witch the raw id is requested
		\result The row index of the element with the specified id.
	*/
	size_type raw_index(id_type id) const
	{
		return id;
	}

	/*!
		Gets an element from a the first position marked as empty.
		The element is not modified. Therefore it will still contain
		the data of the element that was previously occupying the
		position or it will be empty if there was no empty position
		and a new element has been created.
	*/
	iterator reclaim()
	{
		return _reclaim(FILL_FIRST);
	}

	/*!
		Gets an element from the first position marked as empty
		past the last element. The element is not modified.
		Therefore it will still contain the data of the element
		that was previously occupying the position or it will be
		empty if there was no empty position and a new element
		has been created.
	*/
	iterator reclaim_back()
	{
		return _reclaim(FILL_APPEND);
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
		size_t pos = id;
		size_t flat = pos - m_first_pos;
		if (m_holes.size() > 0) {
			auto holeBound = lower_bound(m_holes.cbegin(), m_holes.cend(), pos);
			size_t nHolesBefore = std::distance(m_holes.cbegin(), holeBound);

			flat -= nHolesBefore;
		}

		return flat;
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
		size_t pos = id;

		// Replace the element
		m_v[pos] = T(std::forward<Args>(args)...);

		// Set the id
		m_v[pos].set_id(pos);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
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
		size_t pos = id;

		// Replace the element
		m_v[pos] = std::move(value);

		// Set the id
		m_v[pos].set_id(pos);

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

		If n is greater than the current container size, the content
		is expanded by inserting at the end as many dummy elements as
		needed to reach a size of n. 

		If n is also greater than the current container capacity, an
		automatic reallocation of the allocated storage space takes
		place.

		Notice that this function changes the actual content of the
		container by inserting or erasing elements from it.

		\param n is the new container size, expressed in number of
		         elements.
	*/
	void resize (size_type n)
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
		// may need to resize the storage to reach the requested size.
		if (n > size()) {
			size_t previousStorageSize = storage_size();
			if (n < previousStorageSize) {
				return;
			}

			storage_resize(n);
			return;
		}

		// If the requested size is smaller that the current size
		// we need to perform a real resize.

		// Find the updated position of the last element
		size_type updated_last_pos = n - 1;
		if (!m_holes.empty()) {
			size_type nHoles = 0;
			std::deque<size_type>::iterator it_begin = m_holes.begin();
			std::deque<size_type>::iterator it_end   = m_holes.end();
			while (true) {
				std::deque<size_type>::iterator it_hole = upper_bound(it_begin, it_end, updated_last_pos);
				nHoles = std::distance(it_begin, it_hole);
				if (nHoles == 0) {
					break;
				}

				updated_last_pos += nHoles;

				it_begin = it_hole;
			}
		}

		// Delete all holes above the updated last position
		holes_delete_after(updated_last_pos);

		// Resize the vector
		storage_resize(updated_last_pos + 1);

		// Update the position of the last element
		m_last_pos = updated_last_pos;
	}

	/*!
		Returns the number of elements in the vector.

		This is the number of actual objects held in the vector,
		which is not necessarily equal to its storage capacity.

		\result The number of elements in the container.
	*/
	size_type size() const
	{
		return m_last_pos - m_first_pos + 1 - m_holes.size();
	}

	/*!
		Requests the container to compact the elements and reduce
		its capacity to fit its size.

		The request is non-binding, and the function can leave
		the vector with a capacity greater than its size.

		This may cause a reallocation, but has no effect on the vector
		size and cannot alter its elements.

		\return A mapper between the old and the new element ids.
	*/
	std::vector<id_type> squeeze()
	{
		// Flush changes
		flush();

		// Get the number of the elements
		size_type nElements = size();

		// Compact the vector
		std::vector<id_type> mapper(nElements);
		if (!m_holes.empty()) {
			// Move the elements
			size_type offset = 0;
			size_type nHoles = m_holes.size();
			for (size_type pos = 0; pos <= m_last_pos; pos++) {
				if (offset < nHoles && m_holes[offset] == pos) {
					mapper[pos] = SENTINEL_ID;
					offset++;
					continue;
				}

				if (offset == 0) {
					mapper[pos] = pos;
					continue;
				}

				id_type id = m_v[pos].get_id();
				size_type updatedPos = pos - offset;

				m_v[updatedPos] = std::move(m_v[pos]);
				m_v[updatedPos].set_id(updatedPos);

				mapper[pos] = updatedPos;
			}

			// Reset first and last counters
			m_first_pos = 0;
			m_last_pos  = nElements - 1;

			// There are no more holes
			std::deque<size_type>().swap(m_holes);
		} else {
			for (size_type pos = 0; pos <= nElements; pos++) {
				mapper[pos] = pos;
			}
		}

		// Resize
		storage_resize(size());
		m_v.shrink_to_fit();

		return mapper;
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
	void swap(PositionalPiercedVector& x) noexcept
	{
		std::swap(x.m_first_pos, m_first_pos);
		std::swap(x.m_last_pos, m_last_pos);
		std::swap(x.m_v, m_v);
		std::swap(x.m_holes, m_holes);
		std::swap(x.m_pending_deletes, m_pending_deletes);
	}

	/*!
		Assigns new contents to the vector, replacing its current
		contents.

		\param other vector with the new contents
	*/
	PositionalPiercedVector& operator=(PositionalPiercedVector&& other) = default;

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
		size_type pos = id;

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
		size_type pos = id;

		return m_v[pos];
	}

	/*!
		Functional for compare the position of two elements
	*/
	struct position_less
	{
		position_less(PositionalPiercedVector<T> &vector)
		{
			m_vector = &vector;
		}

		bool operator()(const id_type &id_1, const id_type &id_2) const
		{
			return id_1 < id_2;
		}

		PositionalPiercedVector<T> *m_vector;
	};


	/*!
		Functional for compare the position of two elements
	*/
	struct position_greater
	{
		position_greater(PositionalPiercedVector<T> &vector)
		{
			m_vector = &vector;
		}

		bool operator()(const id_type &id_1, const id_type &id_2) const
		{
			return id_1 > id_2;
		}

		PositionalPiercedVector<T> *m_vector;
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
		Vector that will hold the elements.
	*/
	std::vector<value_type>m_v;

	/*!
		Container that will hold a list of the holes present in
		the piecrecd vector.
	*/
	std::deque<size_type> m_holes;

	/*!
		Tracks if the container in a dirty status.
	*/
	bool m_dirty;

	/*!
		Container that will hold a list of the pending deletes.
	*/
	std::deque<size_type> m_pending_deletes;

	/*!
		Position of the first element in the internal vector.
	*/
	size_type m_first_pos;

	/*!
		Position of the last element in the internal vector.
	*/
	size_type m_last_pos;

	/*!
		The container is extended by inserting a new element in
		the specified hole. This new element is constructed in
		place using args as the arguments for its construction.

		\param fillType is the fill-pattern that will be used to
		identify the position
		\param args the arguments forwarded to construct the new element
		\result An iterator that points to the newly inserted element.
	*/
	template <class... Args>
	iterator _emplace(FillType fillType, Args&&... args)
	{
		// Position of the element
		size_type pos = fill_pos(fillType);

		// Insert the element
		m_v[pos] = T(std::forward<Args>(args)...);

		// Set the id
		m_v[pos].set_id(pos);

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
		// Free the position
		if (delayed) {
			pending_deletes_add(pos);
		} else {
			pierce_pos(pos);
		}

		// Return the iterator to the element following the one erased
		iterator itr;
		if (empty() || pos >= m_last_pos) {
			itr = end();
		} else {
			itr = raw_begin() + next_used_pos(pos);
		}

		return itr;
	}

	/*!
		The container is extended by inserting a new element in
		the specified hole.

		\param fillType is the fill-pattern that will be used to
		identify the position
		\param value is the value to be copied (or moved) to the
		            inserted elements.
		\result An iterator that points to the the newly inserted
		        element.

	*/
	iterator _insert(FillType fillType, value_type &&value)
	{
		// Position of the element
		size_type pos = fill_pos(fillType);

		// Insert the element
		m_v[pos] = std::move(value);

		// Set the id
		m_v[pos].set_id(pos);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Gets an element from a position marked as empty. The element
		is not modified. Therefore it will still contain
		the data of the element that was previously occupying the
		position or it will be empty if there was no empty position
		and a new element has been created.

		\param fillType is the fill-pattern that will be used to
		identify the position
		\result An iterator that points to the the reclaimed element.
	*/
	iterator _reclaim(FillType fillType)
	{
		// Position of the element
		size_type pos = fill_pos(fillType);

		// Set the id of the element
		m_v[pos].set_id(pos);

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
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
		    if (x.get_id() >= 0 && y.get_id() < 0) {
			    return true;
		    } else if (x.get_id() < 0 && y.get_id() >= 0) {
			    return false;
		    } else if (x.get_id() >= 0) {
			    return (x.get_id() < y.get_id());
		    } else {
			    return (x.get_id() > y.get_id());
		    }
	    }
	};

	/*!
		Gets a position in which store an element.

		\param fillType is the fill-pattern that will be used to
		identify the position
	*/
	size_type fill_pos(FillType fillType)
	{
		// If the container is empty or if there are no holes nor
		// pending deletes, elements can only be appened at the
		// end of the vector.
		if (fillType == FILL_FIRST) {
			if (empty() || (m_holes.empty() && m_pending_deletes.empty())) {
				fillType = FILL_APPEND;
			}
		}

		// Find the position
		size_type pos;
		if (fillType == FILL_APPEND) {
			if (!empty()) {
				m_last_pos++;
			}
			storage_resize(m_last_pos + 1);

			pos = m_last_pos;
			if (!m_pending_deletes.empty()) {
				pending_deletes_delete(pos);
			}
		} else if (fillType == FILL_FIRST) {
			if (m_pending_deletes.empty()) {
				// Get first hole
				pos = holes_pop();

				// Update first and last counters
				if (m_last_pos < pos) {
					m_last_pos = pos;
				}

				if (m_first_pos > pos) {
					m_first_pos = pos;
				}

				// If previos element is a hole, its id need to be udated
				if (pos > 0 && is_pos_empty(pos - 1)) {
					update_empty_pos_id(pos - 1);
				}
			} else {
				pos = pending_deletes_pop_back();
			}

			assert(pos < m_v.size() - 1);
		}

		return pos;
	}

	/*!
		Add a position to the holes.

		The container that hold the hole list has to be always
		kept in ascending order.

		\param pos the position to be added to the holes
	*/
	void holes_add(size_type pos)
	{
		std::deque<size_type>::iterator itr = lower_bound(m_holes.begin(), m_holes.end(), pos);
		m_holes.insert(itr, pos);
	}

	/*!
		Removes a hole.

		\param hole is the hole to be removed
		\result Returns true if the hole was found and successfuly
		removed, otherwise it returns false.
	*/
	bool holes_delete(size_type hole)
	{
		std::deque<size_type>::iterator itr = lower_bound(m_holes.begin(), m_holes.end(), hole);
		if (itr == m_holes.end()) {
			return false;
		}

		m_holes.erase(itr);
		return true;
	}

	/*!
		Deletes all holes after a specified position.

		\param pos the position after wich all holes have to be
		           deleted
	*/
	void holes_delete_after(size_type pos)
	{
		if (m_holes.empty()) {
			return;
		} else if (m_holes.back() <= pos) {
			return;
		}

		std::deque<size_type>::iterator itr = upper_bound(m_holes.begin(), m_holes.end(), pos);
		m_holes.erase(itr, m_holes.end());
	}

	/*!
		Gets the position associated with the first hole and deletes
		that hole from the list.

		\result The position associated with the first hole.
	*/
	size_type holes_pop()
	{
		size_type pos = m_holes.front();
		m_holes.pop_front();

		return pos;
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
		Returns the first non-empty position after the specified
		starting position.

		If the starting position is the last posistion, an
		exception is thrown.

		\param pos starting position
		\result The firt non-empty position after the starting
		        position.
	*/
	size_type next_used_pos(size_type pos)
	{
		if (pos == m_last_pos) {
			throw std::out_of_range ("Already in the last position");
		}

		size_type next_pos = pos + 1;
		id_type next_id    = m_v[next_pos].get_id();
		if (next_id >= 0) {
			return next_pos;
		} else {
			return next_pos - next_id;
		}
	}

	/*!
		Add a position to the pending deletes list.

		The container that hold the pending deletes has to be always
		kept in ascending order.

		\param pos the position to be added to the pending deletes list
	*/
	void pending_deletes_add(size_type pos)
	{
		m_dirty = true;

		std::deque<size_type>::iterator itr = lower_bound(m_pending_deletes.begin(), m_pending_deletes.end(), pos);
		m_pending_deletes.insert(itr, pos);
	}

	/*!
		Gets the position associated with the first pending delete
                and deletes that pending delete from the list.

		\result The position associated with the first pending delete.
	*/
	size_t pending_deletes_pop_front()
	{
		size_t pos = m_pending_deletes.front();
		m_pending_deletes.pop_front();

		return pos;
	}

	/*!
		Gets the position associated with the first pending delete
                and deletes that pending delete from the list.

		\result The position associated with the first pending delete.
	*/
	size_t pending_deletes_pop_back()
	{
		size_t pos = m_pending_deletes.back();
		m_pending_deletes.pop_back();

		return pos;
	}

	/*!
		Removes a pending delete.

		\param pending_delete is the pending delete to be removed
		\result Returns true if the pending delete was found and
		successfuly removed, false otherwise.
	*/
	bool pending_deletes_delete(size_type pending_delete)
	{
		std::deque<size_type>::iterator itr = lower_bound(m_pending_deletes.begin(), m_pending_deletes.end(), pending_delete);
		if (itr == m_pending_deletes.end()) {
			return false;
		}

		m_pending_deletes.erase(itr);
		return true;
	}

	/*!
		Mark the position as empty.

		\param pos the position to be marked as empty
	*/
	void pierce_pos(size_type pos)
	{
		// Update first and last counters
		if (empty()) {
			m_last_pos  = 0;
			m_first_pos = 0;
		} else {
			if (m_last_pos == pos) {
				m_last_pos = prev_used_pos(pos);
			}

			if (m_first_pos == pos) {
				m_first_pos = next_used_pos(pos);
			}
		}

		// Update id of the empty element
		update_empty_pos_id(pos);

		// Hole
		//
		// An empty position is considered a hole, only if it's before
		// the last used position. All holes after the last used
		// position need to be removed.
		if (pos < m_last_pos) {
			holes_add(pos);
		} else {
			holes_delete_after(m_last_pos);
		}
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
	size_type prev_used_pos(size_type pos)
	{
		if (pos == m_first_pos) {
			throw std::out_of_range ("Already in the firts position");
		}

		size_type prev_pos = pos - 1;
		id_type prev_id = m_v[prev_pos].get_id();

		if (prev_id >= 0) {
			return prev_pos;
		} else {
			return prev_used_pos(prev_pos);
		}
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
		m_v.resize(n + REQUIRED_SENTINEL_COUNT);

		size_t current_raw_size = m_v.size();
		for (size_t k = std::min(n, previous_raw_size); k < current_raw_size; ++k) {
			m_v[k].set_id(SENTINEL_ID);
		}
	}

	/*!
		Updates the id of the element in the specified position to make
		it an empty element. If needed, updates also of the id of the
		element before the specified position.

		The id of an empty element contains the distance, measured in
		number of elements, between the current element the next
		non-empty element (the distance is negative). The id of an
		element past the last non-empty position is set to the special
		value SENTINEL_ID.

		\param pos the specified position
	*/
	void update_empty_pos_id(size_type pos)
	{
		// Position of the next non-empty element
		id_type distanceFromNonEmpty;
		if (empty() || pos > m_last_pos) {
			distanceFromNonEmpty = SENTINEL_ID;
		} else if (is_pos_empty(pos + 1)) {
			distanceFromNonEmpty = m_v[pos + 1].get_id() - 1;
		} else {
			distanceFromNonEmpty = - 1;
		}

		// Update the id of the element in the current position
		m_v[pos].set_id(distanceFromNonEmpty);

		// Update the id of the elements the previous positions
		if (pos > m_first_pos) {
			size_type prevPos = pos - 1;
			while (prevPos >= m_first_pos && is_pos_empty(prevPos)) {
				m_v[prevPos].set_id(distanceFromNonEmpty - (pos - prevPos));

				--prevPos;
			}
		}
	}

};

// Definition of static constants of PositionalPiercedVector
template<class T>
const typename PositionalPiercedVector<T>::id_type
	PositionalPiercedVector<T>::SENTINEL_ID = std::numeric_limits<id_type>::min();

template<class T>
const typename PositionalPiercedVector<T>::size_type
	PositionalPiercedVector<T>::REQUIRED_SENTINEL_COUNT = 1;

template<class T>
const typename PositionalPiercedVector<T>::size_type
	PositionalPiercedVector<T>::USABLE_POS_COUNT = std::numeric_limits<size_type>::max() - REQUIRED_SENTINEL_COUNT;

}

#endif
