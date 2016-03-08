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

/*!
	Creates a new uninitialized iterator
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::PiercedIterator()
	: m_container(nullptr), m_pos(0)
{
}

/*!
	Exchanges the values of the current iterator and
	the iterator recevied as argument.

	\param other the iterator to exchange values with
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
void PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::swap(PiercedIterator& other) noexcept
{
	std::swap(m_container, other.m_container);
	std::swap(m_pos, other.m_pos);
}

/*!
	Pre-increment operator.
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
PiercedIterator<T, id_type, T_no_cv, id_type_no_cv> & PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::operator++()
{
	size_t delta = 1;
	while (true) {
		m_pos += delta;
		if (m_pos > m_container->m_last_pos) {
			m_pos = m_container->m_last_pos + 1;
			return *this;
		}

		id_type id = m_container->m_ids[m_pos];
		if (id >= 0) {
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
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
PiercedIterator<T, id_type, T_no_cv, id_type_no_cv> PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::operator++(int)
{
	PiercedIterator tmp(m_container, m_pos);

	++(*this);

	return tmp;
}

/*!
	Deference operator.

	\result A reference to the element currently pointed to by the
			iterator.
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
T & PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::operator*() const
{
	return m_container->m_v[m_pos];
}

/*!
	Deference operator.

	\result A reference to the element currently pointed to by the
			iterator.
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
T * PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::operator->() const
{
	return &(m_container->m_v[m_pos]);
}

/*!
	Converts the iterator to a const_iterator.
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::operator PiercedIterator<const T, const id_type>() const
{
	return PiercedIterator<const T, const id_type>(m_container, m_pos);
}

/*!
	Creates a new iterator and initializes it with the position of
	the const base iterator recevied in input.
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::PiercedIterator(Container<T_no_cv, id_type_no_cv> *container, const size_t &pos)
	: m_container(const_cast<Container<T_no_cv, id_type_no_cv> *>(container)), m_pos(pos)
{
}

/*!
	Creates a new iterator and initializes it with the position of
	the const base iterator recevied in input.
*/
template<typename T, typename id_type, typename T_no_cv, typename id_type_no_cv>
PiercedIterator<T, id_type, T_no_cv, id_type_no_cv>::PiercedIterator(const Container<T_no_cv, id_type_no_cv> *container, const size_t &pos)
	: m_container(const_cast<Container<T_no_cv, id_type_no_cv> *>(container)), m_pos(pos)
{
}

// Definition of static constants of PiercedVector
template<typename T, typename id_type>
const std::size_t
	PiercedVector<T, id_type>::MAX_PENDING_HOLES = 16384;

/*!
	Constructs an empty pierced vector with no elements.
*/
template<typename T, typename id_type>
PiercedVector<T, id_type>::PiercedVector()
{
	clear();
}

/*!
	Constructs a pierced vector with a capacity at least enough
	to contain n elements.

	\param n the minimum capacity requested for the vector
*/
template<typename T, typename id_type>
PiercedVector<T, id_type>::PiercedVector(std::size_t n)
{
	clear();

	reserve(n);
}

/*!
	Adds a new element at the end of the vector, after its current
	last element.

	The content of value is copied (or moved) to the new element.

	\param value the value to be copied (or moved) to the new
					element
	\result An iterator that points to the newly inserted element.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::push_back(const id_type &id, value_type &&value)
{
	// Fille the position
	std::size_t pos = fill_pos_append(id);

	// Insert the element
	m_v[pos] = std::move(value);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::reclaim(const id_type &id)
{
	std::size_t pos = fill_pos_head(id);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::reclaim_after(const id_type &referenceId, const id_type &id)
{
	std::size_t pos = fill_pos_after(get_pos_from_id(referenceId), id);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::reclaim_back(const id_type &id)
{
	std::size_t pos = fill_pos_append(id);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::reclaim_before(const id_type &referenceId, const id_type &id)
{
	std::size_t pos = fill_pos_before(get_pos_from_id(referenceId), id);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::move_after(const id_type &referenceId, const id_type &id, bool delayed)
{
	// Save the element
	std::size_t initialPos = get_pos_from_id(id);
	T temp = std::move(m_v[initialPos]);

	// Pierce the position
	pierce_pos(initialPos, !delayed);

	// Insert the element in the updated position
	return insert_after(referenceId, id, std::move(temp));
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::move_before(const id_type &referenceId, const id_type &id, bool delayed)
{
	// Save the element
	std::size_t initialPos = get_pos_from_id(id);
	T temp = std::move(m_v[initialPos]);

	// Pierce the position
	pierce_pos(initialPos, !delayed);

	// Insert the element in the updated position
	return insert_before(referenceId, id, std::move(temp));
}

/*!
	The container is extended by inserting a new element.

	\param value is the value to be copied (or moved) to the
				inserted elements.
	\result An iterator that points to the newly inserted element.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::insert(const id_type &id, value_type &&value)
{
	// Fill the position
	std::size_t pos = fill_pos_head(id);

	// Insert the element
	m_v[pos] = std::move(value);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::insert_after(const id_type &referenceId, const id_type &id, value_type &&value)
{
	// Fill the position
	std::size_t pos = fill_pos_after(get_pos_from_id(referenceId), id);

	// Insert the element
	m_v[pos] = std::move(value);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::insert_before(const id_type &referenceId, const id_type &id, value_type &&value)
{
	// Fill the position
	std::size_t pos = fill_pos_before(get_pos_from_id(referenceId), id);

	// Insert the element
	m_v[pos] = std::move(value);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
}

/*!
	The element with the specified id is replaced with a new element.

	\param id is the id of the element that will be replaced
	\param value is the value to be moved to the inserted elements.
	\result An iterator that points to the newly inserted element.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::replace(id_type id, value_type &&value)
{
	// Position
	size_t pos = m_pos.at(id);

	// Replace the element
	m_v[pos] = std::move(value);

	// Update the map
	id_type id_prev = m_ids[pos];
	if (id != id_prev) {
		set_pos_id(pos, id);
	}

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
}

/*!
	Updates the id of the specified element.

	\param currentId is the current id of the element
	\param updatedId is the new id of the element
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::update_id(const id_type &currentId, const id_type &updatedId)
{
	set_pos_id(get_pos_from_id(currentId), updatedId);
}

/*!
	The container is extended by inserting a new element. This
	new element is constructed in place using args as the
	arguments for its construction.

	\param args the arguments forwarded to construct the new element
	\result An iterator that points to the the newly inserted
			element.
*/
template<typename T, typename id_type>
template<class... Args>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::emplace(const id_type &id, Args&&... args)
{
	// Fill the position
	std::size_t pos = fill_pos_head(id);

	// Insert the element
	m_v[pos] = T(std::forward<Args>(args)...);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
template<class... Args>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::emplace_after(const id_type &referenceId, const id_type &id, Args&&... args)
{
	// Fill the position
	std::size_t pos = fill_pos_after(get_pos_from_id(referenceId), id);

	// Insert the element
	m_v[pos] = T(std::forward<Args>(args)...);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
}

/*!
	Inserts a new element at the end of the vector, right after
	its current last element. This new element is constructed
	in place using args as the arguments for its construction.

	\param args the arguments forwarded to construct the new element
*/
template<typename T, typename id_type>
template<class... Args>
void PiercedVector<T, id_type>::emplace_back(const id_type &id, Args&&... args)
{
	// Fill the position
	std::size_t pos = fill_pos_append(id);

	// Insert the element
	m_v[pos] = T(std::forward<Args>(args)...);
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
template<typename T, typename id_type>
template<class... Args>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::emplace_before(const id_type &referenceId, const id_type &id, Args&&... args)
{
	// Fill the position
	std::size_t pos = fill_pos_before(get_pos_from_id(referenceId), id);

	// Insert the element
	m_v[pos] = T(std::forward<Args>(args)...);

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
}

/*!
	The element with the specified id is replaced with a new element.
	This new element is constructed in place using args as the
	arguments for its construction.

	\param id is the id of the element that will be replaced
	\param args the arguments forwarded to construct the new element
	\result An iterator that points to the newly inserted element.
*/
template<typename T, typename id_type>
template<class... Args>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::emreplace(id_type id, Args&&... args)
{
	// Get the position of the element
	size_t pos = m_pos.at(id);

	// Replace the element
	m_v[pos] = T(std::forward<Args>(args)...);

	// Update the map
	id_type id_prev = m_ids[pos];
	if (id != id_prev) {
		set_pos_id(pos, id);
	}

	// Return the iterator that points to the element
	return get_iterator_from_pos(pos);
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
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::erase(id_type id, bool delayed)
{
	// Position
	size_t pos = m_pos.at(id);

	// Pierce the position
	pierce_pos(pos, !delayed);

	// Return the iterator to the next element
	if (empty() || pos >= m_last_pos) {
		return end();
	} else {
		return get_iterator_from_pos(find_next_used_pos(pos));
	}
}

/*!
	Removes the last element in the vector, effectively reducing
	the container size by one.

	Element is not deleted from the internal vector, instead its
	id is changed to mark the position as empty and allow the
	container to reuse that position.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::pop_back()
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	} else if (size() == 1) {
		clear();
	} else {
		std::size_t updated_last_pos = find_prev_used_pos(m_last_pos);
		update_last_used_pos(updated_last_pos);
	}
}

/*!
	Swap the elements with the specified id.

	\param id_first is the id of the first element to be swapped
	\param id_second is the id of the second element to be swapped
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::swap(const id_type &id_first, const id_type &id_second)
{
	// Positions
	size_t pos_first  = m_pos.at(id_first);
	size_t pos_second = m_pos.at(id_second);

	// Swap the elements
	T tmp = std::move(m_v[pos_first]);
	m_v[pos_first]  = std::move(m_v[pos_second]);
	m_v[pos_second] = std::move(tmp);

	// Relink the ids
	set_pos_id(pos_second, id_first);
	set_pos_id(pos_first, id_second);
}

/*!
	Removes all elements from the vector (which are destroyed),
	leaving the container with a size of 0.

	\param release if it's true the memory hold by the container will
	be released, otherwise the container will be cleared but its
	memory will not be relased
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::clear(bool release)
{
	// Clear storage
	m_ids.clear();
	m_v.clear();
	if (release) {
		std::vector<id_type>().swap(m_ids);
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
	std::unordered_map<id_type, std::size_t, PiercedHasher>().swap(m_pos);

	// There are no dirty positions
	m_first_dirty_pos = m_last_pos + 1;
}

/*!
	Flush all pending changes.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::flush()
{
	// Flush pending holes
	holes_flush();
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
template<typename T, typename id_type>
void PiercedVector<T, id_type>::reserve(std::size_t n)
{
	m_ids.reserve(n);
	m_v.reserve(n);
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
template<typename T, typename id_type>
void PiercedVector<T, id_type>::resize(std::size_t n)
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
	std::size_t updated_last_pos = get_pos_from_id(last_id);

	// Update the last position
	update_last_used_pos(updated_last_pos);
}

/*!
	Sorts the elements of the vector in ascending id order.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::sort()
{
	// Squeeze the container
	squeeze();

	// The container has been squeezed, there are no holes
	size_t containerSize = size();

	// Evaluates the permutations
	//
	// The permutation for ids and values are the same, however the reored
	// function will destroy the permutation on output, so two different
	// permutations are needed.
	std::vector<std::size_t> id_permutation;
	id_permutation.resize(containerSize);
	std::iota(id_permutation.begin(), id_permutation.end(), 0);
	std::sort(id_permutation.begin(), id_permutation.end(), id_less(m_ids));

	std::vector<std::size_t> value_permutation(id_permutation);

	// Sort the container
	reorder_vector<id_type>(id_permutation, m_ids, containerSize);
	reorder_vector<T>(value_permutation, m_v, containerSize);
}

/*!
	Requests the container to compact the elements and reduce
	its capacity to fit its size.

	The request is non-binding, and the function can leave
	the vector with a capacity greater than its size.

	This may cause a reallocation, but has no effect on the vector
	size and cannot alter its elements.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::squeeze()
{
	// Flush changes
	flush();

	// Compact the vector
	std::size_t nHoles = holes_count();
	if (nHoles != 0) {
		// Move the elements
		std::size_t firstPosToUpdate;
		if (m_first_pos == 0) {
			firstPosToUpdate = *(m_holes_regular_end - 1);
		} else {
			firstPosToUpdate = 0;
		}

		std::size_t offset = 0;
		for (std::size_t pos = firstPosToUpdate; pos <= m_last_pos; pos++) {
			if (offset < nHoles && *(m_holes_regular_end - offset - 1) == pos) {
				++offset;
				continue;
			}

			id_type id = m_ids[pos];
			std::size_t updatedPos = pos - offset;

			set_pos_id(updatedPos, id);
			m_v[updatedPos] = std::move(m_v[pos]);
		}

		// Clear the holes
		holes_clear();

		// Reset first and last counters
		update_first_used_pos(0);
		update_last_used_pos(size() - 1);
	}

	// Shrink to fit
	m_ids.shrink_to_fit();
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
template<typename T, typename id_type>
void PiercedVector<T, id_type>::swap(PiercedVector& x) noexcept
{
	std::swap(x.m_first_pos, m_first_pos);
	std::swap(x.m_last_pos, m_last_pos);
	std::swap(x.m_first_dirty_pos, m_first_dirty_pos);
	std::swap(x.m_ids, m_ids);
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
	Returns the size of the storage space currently allocated
	for the vector, expressed in terms of elements.

	\result The size of the currently allocated storage capacity
			in the vector, measured in terms of the number elements
			it can hold.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::capacity()
{
	return m_v.capacity();
}

/*!
	Returns whether the vector is contiguous (i.e. whether it contains
	no holes).

	\result true if the container is contiguous, false otherwise.
*/
template<typename T, typename id_type>
bool PiercedVector<T, id_type>::contiguous() const
{
	return m_holes.empty();
}

/*!
	Dumps to screen the internal data.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::dump()
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
		std::cout << m_ids[k] << std::endl;
	}

	std::cout << std::endl;
	std::cout << " Poistion map: " << std::endl;
	for (auto itr = m_pos.cbegin(); itr != m_pos.cend(); ++itr) {
		std::cout << itr->first << " -> " << itr->second << std::endl;
	}

	std::cout << "----------------------------------------" << std::endl;
}

/*!
	Returns whether the vector is empty (i.e. whether its size
	is 0).

	\result true if the container size is 0, false otherwise.
*/
template<typename T, typename id_type>
bool PiercedVector<T, id_type>::empty() const
{
	return m_pos.empty();
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
template<typename T, typename id_type>
bool PiercedVector<T, id_type>::is_iterator_slow()
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
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::max_size() const
{
	return m_v.max_size();
}

/*!
	Returns the number of elements in the vector.

	This is the number of actual objects held in the vector,
	which is not necessarily equal to its storage capacity.

	\result The number of elements in the container.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::size() const
{
	return m_pos.size();
}

/*!
	Checks if a given id exists in the vector.

	\param id the id to look for
	\result Returns true is the given id exists in the vector,
			otherwise it returns false.
*/
template<typename T, typename id_type>
bool PiercedVector<T, id_type>::exists(id_type id)
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
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::extract_flat_index(id_type id) const
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
	Gets a vector containing the ids of the elements stored in
	the vector.

	\param ordered if is true the ids will be sorted in ascending
	\              order, otherwise the ids will be in random
					order.
	\result A vector with the id of the elements in stored in the
			vector.
*/
template<typename T, typename id_type>
std::vector<id_type> PiercedVector<T, id_type>::get_ids(bool ordered)
{
	// Initialize the vector wth the ids
	std::vector<id_type> ids;
	ids.reserve(size());

	// Extract the ids
	size_t pos = m_first_pos;

	ids.push_back(m_ids[pos]);
	do {
		pos = find_next_used_pos(pos);
		ids.push_back(m_ids[pos]);
	} while (pos != m_last_pos);

	// Sort the ids
	if (ordered) {
		std::sort(ids.begin(), ids.end());
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
template<typename T, typename id_type>
id_type PiercedVector<T, id_type>::get_size_marker(const size_t &targetSize, const id_type &fallback)
{
	// If the size is zero, we return the first element, if the target
	// size is equal to the size minus one we return the last element,
	// if the target size is greater or equal the current container size
	// we return the fallback value.
	if (targetSize == 0) {
		return m_ids[m_first_pos];
	} else if (targetSize == (size() - 1)) {
		return m_ids[m_last_pos];
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

	std::size_t nEmpties  = 0;
	std::size_t markerPos = targetSize;
	while (true) {
		if (is_pos_empty(markerPos)) {
			markerPos = find_next_used_pos(markerPos - 1);
		}

		// Count the number of holes and pending deletes before the
		// current marker position
		if (regular_hole_itr != m_holes_regular_begin) {
			hole_iterator itr_previous = regular_hole_itr;
			regular_hole_itr = std::upper_bound(m_holes_regular_begin, regular_hole_itr, markerPos, std::greater<std::size_t>());
			nEmpties += std::distance(regular_hole_itr, itr_previous);
		}

		if (pending_hole_itr != m_holes_pending_begin) {
			hole_iterator itr_previous = pending_hole_itr;
			pending_hole_itr = std::upper_bound(m_holes_pending_begin, pending_hole_itr, markerPos, std::greater<std::size_t>());
			nEmpties += std::distance(pending_hole_itr, itr_previous);
		}

		// Get the marker size
		//
		// If we have reached the target size we can exit, otherwise
		// we update the marker and we continue iterating
		std::size_t markerSize = markerPos - nEmpties;
		if (markerSize == targetSize) {
			break;
		} else {
			markerPos += targetSize - markerSize;
		}
	}

	return m_ids[markerPos];
}

/*!
	Returns a direct pointer to the memory array used internally
	by the vector to store its owned elements.

	\result A pointer to the first element in the array used
			internally by the vector.

*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::value_type * PiercedVector<T, id_type>::data() noexcept
{
	return m_v.data();
}

/*!
	Returns a reference to the last element of the vector. If
	the vector is empty, an exception is thrown.

	\result A reference to the last element of the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::back()
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
template<typename T, typename id_type>
const typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::back() const
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	}

	return m_v[m_last_pos];
}

/*!
	Returns a reference to the first element of the vector. If
	the vector is empty, an exception is thrown.

	\result A reference to the first element of the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::front()
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
template<typename T, typename id_type>
const typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::front() const
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	}

	return m_v[m_first_pos];
}

/*!
	Returns a reference to the element with the specified id. If
	there is no element with the specified id, an exception is
	thrown.

	\param id the id of the element
	\result A reference to the element with the specified id.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::at(const id_type &id)
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
template<typename T, typename id_type>
const typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::at(const id_type &id) const
{
	return (*this)[id];
}

/*!
	Returns a reference to the element at the specified position.

	\param pos the position of the element
	\result A reference to the element in the specified position.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::raw_at(const std::size_t &pos)
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
template<typename T, typename id_type>
const typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::raw_at(const std::size_t &pos) const
{
	return m_v[pos];
}

/*!
	Gets the row index of the element with the specified id.

	If there is no element with the specified id, an exception is
	thrown.

	\param id the id of the element for witch the raw id is requested
	\result The row index of the element with the specified id.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::raw_index(id_type id) const
{
	return get_pos_from_id(id);
}

/*!
	Returns a constant reference to the element with the
	specified id. If there is no element with the specified id,
	an exception is thrown.

	\param id the id of the element
	\result A constant reference to the element with the specified
			id.
*/
template<typename T, typename id_type>
const typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::operator[](const id_type &id) const
{
	std::size_t pos = get_pos_from_id(id);

	return m_v[pos];
}

/*!
	Returns a reference to the element with the
	specified id. If there is no element with the specified id,
	an exception is thrown.

	\param id the id of the element
	\result A reference to the element with the specified id.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::value_type & PiercedVector<T, id_type>::operator[](const id_type &id)
{
	std::size_t pos = get_pos_from_id(id);

	return m_v[pos];
}

/*!
	Gets an iterator pointing to the specified element.

	\param id is the id of the specified iterator.
	\result An iterator pointing to the specified element.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::get_iterator(const id_type &id) noexcept
{
	const size_t pos = get_pos_from_id(id);

	return get_iterator_from_pos(pos);
}

/*!
	Gets a constant iterator pointing to the specified element.

	\param id is the id of the specified iterator.
	\result A constant iterator pointing to the specified element.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::const_iterator PiercedVector<T, id_type>::get_const_iterator(const id_type &id) const noexcept
{
	const size_t pos = get_pos_from_id(id);

	return get_const_iterator_from_pos(pos);
}

/*!
	Returns an iterator pointing to the first element in the
	vector.

	\result An iterator pointing to the first element in the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::begin()
{
	return get_iterator_from_pos(m_first_pos);
}

/*!
	Returns an iterator referring to the past-the-end element
	in the vector.

	\result An iterator referring to the past-the-end element
			in the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::end()
{
	return get_iterator_from_pos(m_last_pos + 1);
}

/*!
	Returns a constant iterator pointing to the first element
	in the vector.

	\result A constant iterator pointing to the first element in
	the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::const_iterator PiercedVector<T, id_type>::begin() const noexcept
{
	return cbegin();
}

/*!
	Returns a constant iterator referring to the past-the-end
	element in the vector.

	\result A constant iterator referring to the past-the-end
	element in the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::const_iterator PiercedVector<T, id_type>::end() const noexcept
{
	return cend();
}

/*!
	Returns an conts_iterator pointing to the first element in the
	vector.

	\result A const_iterator pointing to the first element in
			the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::const_iterator PiercedVector<T, id_type>::cbegin() const noexcept
{
	return get_const_iterator_from_pos(m_first_pos);
}

/*!
	Returns an const_iterator referring to the past-the-end element
	in the vector.

	\result A const_iterator referring to the past-the-end element
			in the vector.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::const_iterator PiercedVector<T, id_type>::cend() const noexcept
{
	return get_const_iterator_from_pos(m_last_pos + 1);
}

/*!
	Returns an iterator pointing to the first element in the
	raw container.

	\result An iterator pointing to the first element in the raw
			container.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::raw_iterator PiercedVector<T, id_type>::raw_begin()
{
	return m_v.begin();
}

/*!
	Returns an iterator referring to the past-the-end element
	in the raw container.

	\result An iterator referring to the past-the-end element
			in the raw container.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::raw_iterator PiercedVector<T, id_type>::raw_end()
{
	return m_v.end();
}

/*!
	Returns a constant iterator pointing to the first element
	in the raw container.

	\result A constant iterator pointing to the first element in
	the raw container.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::raw_const_iterator PiercedVector<T, id_type>::raw_begin() const noexcept
{
	return raw_cbegin();
}

/*!
	Returns a constant iterator referring to the past-the-end
	element in the raw container.

	\result A constant iterator referring to the past-the-end
	element in the raw container.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::raw_const_iterator PiercedVector<T, id_type>::raw_end() const noexcept
{
	return raw_cend();
}

/*!
	Returns an conts_iterator pointing to the first element in the
	raw container.

	\result A const_iterator pointing to the first element in
			the raw container.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::raw_const_iterator PiercedVector<T, id_type>::raw_cbegin() const noexcept
{
	return m_v.cbegin();
}

/*!
	Returns an const_iterator referring to the past-the-end element
	in raw container.

	\result A const_iterator referring to the past-the-end element
			in raw container.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::raw_const_iterator PiercedVector<T, id_type>::raw_cend() const noexcept
{
	return m_v.cend();
}

/*!
	Gets an iterator pointing to the element in the specified position.

	\param pos is the position of the element
	\result An iterator pointing to the element in the specified position.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::iterator PiercedVector<T, id_type>::get_iterator_from_pos(const std::size_t &pos) noexcept
{
	if (empty() || pos > m_last_pos) {
		return iterator(this, m_last_pos + 1);
	} else {
		return iterator(this, pos);
	}
}

/*!
	Gets a constant iterator pointing to the element in the specified
	position.

	\param pos is the position of the element
	\result A constant iterator pointing to the element in the specified
	position.
*/
template<typename T, typename id_type>
typename PiercedVector<T, id_type>::const_iterator PiercedVector<T, id_type>::get_const_iterator_from_pos(const std::size_t &pos) const noexcept
{
	if (empty() || pos > m_last_pos) {
		return const_iterator(this, m_last_pos + 1);
	} else {
		return const_iterator(this, pos);
	}
}

/*!
	Fills the specified position with a given id.

	\param pos is the position to fill
	\param id is the id that eill be associated to the position
	\result The position that has bill filled.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos(const std::size_t &pos, const id_type &id)
{
	// Ids needs to be positive
	if (id < 0) {
		throw std::out_of_range ("Negative id");
	}

	// Handle duplicate ids
	if (exists(id)) {
		throw std::out_of_range ("Duplicate id");
	}

	// Associate an id to the position
	set_pos_id(pos, id);

	return pos;
}

/*!
	Gets a position for storing a new element.

	The position as alyaws appeneded to the end of the container.

	\result A position for storing a new element.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos_append(const id_type &id)
{
	return fill_pos_insert(storage_size(), id);
}

/*!
	Gets a position for storing a new element.

	The specified position is made available.

	\param pos is the position the will be make available
	\result A position for storing a new element.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos_insert(const std::size_t &pos, const id_type &id)
{
	assert(pos <= m_last_pos + 1);
	if (pos > m_last_pos + 1) {
		throw std::out_of_range ("Unable to insert elements past the last position");
	}

	// Extend the container
	update_last_used_pos(storage_size());

	// Make room for the new element
	if (pos < m_last_pos) {
		// Shift the elements after the requested position
		for (size_t i = m_last_pos; i > pos; --i) {
			id_type id = m_ids[i - 1];
			if (id > 0) {
				set_pos_id(i, id);
				m_v[i] = std::move(m_v[i - 1]);
			} else {
				m_ids[i] = m_ids[i - 1];
			}
		}

		// Reset the, now empty, element
		m_v[pos] = T();

		// Update the regular holes
		if (m_holes_regular_begin != m_holes_regular_end) {
			hole_iterator change_begin = m_holes_regular_begin;
			hole_iterator change_end   = upper_bound(m_holes_regular_begin, m_holes_regular_end, pos, std::greater<std::size_t>());
			for (auto itr = change_begin; itr != change_end; itr++) {
				(*itr)++;
			}
		}

		// Update the pending holes
		if (m_holes_pending_begin != m_holes_pending_end) {
			hole_iterator change_begin = m_holes_pending_begin;
			hole_iterator change_end   = upper_bound(m_holes_pending_begin, m_holes_pending_end, pos, std::greater<std::size_t>());
			for (auto itr = change_begin; itr != change_end; itr++) {
				(*itr)++;
			}
		}
	}

	// Fill the position
	fill_pos(pos, id);

	return pos;
}

/*!
	Gets a position for storing a new element.

	The first available position starting from the head of the container
	is returned.

	\result A position for storing a new element.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos_head(const id_type &id)
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
		std::size_t pos = m_holes.back();
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

		// Fill the position
		return fill_pos(pos, id);
	}

	// There are no holes nor pending delete: use an append fill.
	return fill_pos_append(id);
}

/*!
	Gets a position for storing a new element.

	The first available position starting from the tail of the container
	is returned.

	\result A position for storing a new element.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos_tail(const id_type &id)
{
	// If there are holes we can fill a hole.
	long nRegulars = holes_count_regular();
	long nPendings = holes_count_pending();
	long nHoles    = nRegulars + nPendings;
	if (nHoles != 0) {
		// First search among pendings and the among regulars
		std::size_t pos;
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

		// If previos element is a hole, its id and the ids of the
		// contigous holes need to be udated
		if (pos > 0) {
			std::size_t nextUsedPos = pos;
			std::size_t prevPos = pos - 1;
			while (is_pos_empty(prevPos)) {
				update_empty_pos_id(prevPos, nextUsedPos);
				if (prevPos > 0) {
					--prevPos;
				} else {
					break;
				}
			}
		}

		// Fill the position
		return fill_pos(pos, id);
	}

	// There are no holes nor pending delete: use an append fill.
	return fill_pos_append(id);
}

/*!
	Gets a position for storing a new element.

	The first available position after the specified position is
	returned.

	\param referencePos is the position of the element after which the
	new available position will be searched for
	\result A position for storing a new element.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos_after(const std::size_t &referencePos, const id_type &id)
{
	// Check if we can fill a hole
	//
	// The last hole should hava a position higher than the reference
	// position
	if (holes_count() != 0) {
		// First search among pendings and the among regulars
		std::size_t lastHole;
		if (m_holes_pending_begin != m_holes_pending_end) {
			holes_sort_pending();
			lastHole = *(m_holes_pending_begin);
		} else {
			holes_sort_regular();
			lastHole = *(m_holes_regular_begin);
		}

		if (lastHole > referencePos) {
			return fill_pos_tail(id);
		}
	}

	// We have to append the element at the end of the vector
	return fill_pos_append(id);
}

/*!
	Gets a position for storing a new element.

	The first available position before the specified position is
	returned.

	\param referencePos is the position of the element before which the
	new available position will be searched for
	\result A position for storing a new element.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::fill_pos_before(const std::size_t &referencePos, const id_type &id)
{
	// Check if we can fill a hole
	//
	// The first hole available should be in a position lower than the
	// reference position.
	if (holes_count() != 0) {
		// First search among pendings and the among regulars
		std::size_t firstHole;
		if (m_holes_pending_begin != m_holes_pending_end) {
			holes_sort_pending();
			firstHole = *(m_holes_pending_end - 1);
		} else {
			holes_sort_regular();
			firstHole = *(m_holes_regular_end - 1);
		}

		if (firstHole < referencePos) {
			return fill_pos_head(id);
		}
	}

	// We have to insert the element at the specified position
	return fill_pos_insert(referencePos, id);
}

/*!
	Marks a position as empty.

	The position is inserted in the list of holes. If the list of pending
	holes is full, a flush is called before adding the hole. This means
	that the hole is always added as a pending hole.

	\param hole is the position of the new hole
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::pierce_pos(const std::size_t &pos, bool flush)
{
	// If removing the last position, there is no need to add the
	// position to the holes, it's enough to update the last position
	// counter or clear the container if this was the last hole.
	if (pos == m_last_pos) {
		if (size() == 1) {
			clear();
		} else {
			std::size_t updated_last_pos = find_prev_used_pos(m_last_pos);
			update_last_used_pos(updated_last_pos);
		}
		return;
	}

	// Remove the id from the map
	id_type id = m_ids[pos];
	m_pos.erase(id);

	// Reset the element
	m_v[pos] = T();

	// Reset the position
	size_t nextUsedPos = find_next_used_pos(pos);
	update_empty_pos_id(pos, nextUsedPos);
	m_first_dirty_pos = std::min(pos, m_first_dirty_pos);

	// If removing the first position, update the counter
	if (pos == m_first_pos) {
		std::size_t updated_first_pos = find_next_used_pos(m_first_pos);
		update_first_used_pos(updated_first_pos);
	}

	// If the list of pending holes is full, flush the holes.
	if (m_holes.size() == m_holes.capacity()) {
		holes_flush();
	}

	// Add the hole at the end of the pending holes
	m_holes.push_back(pos);
	m_holes_pending_end = m_holes.end();

	// Check if pending holes are still sorted
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
	Clear the list of available holes.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::holes_clear(bool release)
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
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::holes_count()
{
	return holes_count_pending() + holes_count_regular();
}

/*!
	Count the pending holes.

	\result The number of pending holes.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::holes_count_pending()
{
	return std::distance(m_holes_pending_begin, m_holes_pending_end);
}

/*!
	Count the regular holes.

	\result The number of regular holes.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::holes_count_regular()
{
	return std::distance(m_holes_regular_begin, m_holes_regular_end);
}

/*!
	Flushes the list of available holes.

	All the pending hole are converted to regular holes and new
	space is reserved for future pending holes.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::holes_flush()
{
	// If there are no pending holes there is nothing to do
	if (m_holes_pending_begin == m_holes_pending_end) {
		return;
	}

	// Update the id of the empty elements
	//
	// The list of pending holes is sorted, in this way we can iterate
	// from the last pending hole in the container to the first one.
	// We start updating the id of the last pending hole, updating also
	// the id of the contiguous holes before that one. Then we advance
	// to the next hole, skipping the positions that have already been
	// updated.
	holes_sort_pending();

	auto itr = m_holes_pending_begin;
	size_t pos = m_last_pos + 1;
	do {
		if (*itr >= pos) {
			itr++;
			continue;
		}

		pos = *itr;
		size_t next_used_pos = find_next_used_pos(pos);
		do {
			update_empty_pos_id(pos, next_used_pos);
			if (pos > 0) {
				pos--;
			} else {
				break;
			}
		} while (is_pos_empty(pos));
	} while (pos > 0 && itr != m_holes_pending_end);

	// Move the pending holes into the list of regular holes
	for (auto itr = m_holes_pending_begin; itr != m_holes_pending_end; ++itr) {
		const std::size_t &pos = *itr;

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
				std::size_t nRegulars = holes_count_regular();
				if (nRegulars > 1 && (*(m_holes_regular_end - 1) > *(m_holes_regular_end - 2))) {
					m_holes_regular_sorted = false;
				}
			}
		}
	}

	// Move the holes at the beginning of the vector
	//
	// The iterators will be updated when clearing the pending holes.
	std::size_t nRegulars = holes_count_regular();
	if (m_holes_regular_begin != m_holes.begin()) {
		std::size_t offset = std::distance(m_holes.begin(), m_holes_regular_begin);
		for (std::size_t k = 0; k < nRegulars; ++k) {
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
template<typename T, typename id_type>
void PiercedVector<T, id_type>::holes_clear_pending()
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
template<typename T, typename id_type>
void PiercedVector<T, id_type>::holes_clear_pending(const long &offset, const long &nRegulars)
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
template<typename T, typename id_type>
void PiercedVector<T, id_type>::holes_sort_pending()
{
	if (m_holes_pending_sorted) {
		return;
	}

	std::sort(m_holes_pending_begin, m_holes_pending_end, std::greater<std::size_t>());
	m_holes_pending_sorted = true;
}

/*!
	Sort the list of regular holes in descendent order
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::holes_sort_regular()
{
	if (m_holes_regular_sorted) {
		return;
	}

	std::sort(m_holes_regular_begin, m_holes_regular_end, std::greater<std::size_t>());
	m_holes_regular_sorted = true;
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
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::find_prev_used_pos(std::size_t pos)
{
	std::size_t prev_pos = pos;
	while (true) {
		if (prev_pos == m_first_pos) {
			throw std::out_of_range ("Already in the firts position");
		}
		prev_pos--;

		id_type prev_id = m_ids[prev_pos];
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
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::find_next_used_pos(std::size_t pos)
{
	std::size_t next_pos   = pos;
	std::size_t next_delta = 1;
	while (true) {
		if (next_pos == m_last_pos) {
			throw std::out_of_range ("Already in the last position");
		}
		next_pos += next_delta;

		id_type next_id = m_ids[next_pos];
		if (next_id >= 0) {
			return next_pos;
		} else {
			next_delta = - next_id;
		}
	}
}

/*!
	Returns if the specified position is empty.

	A position is considered empty if the element in that
	position has an id less than 0.

	\param pos the position to check
	\result true is the position is empty, false otherwise.
*/
template<typename T, typename id_type>
bool PiercedVector<T, id_type>::is_pos_empty(std::size_t pos)
{
	return (m_ids[pos] < 0);
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
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::get_pos_from_id(id_type id) const
{
	return m_pos.at(id);
}

/*!
	Associate a position to the specified id.

	\param pos is the position to associate
	\param id is the id that will be associated to the position
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::set_pos_id(const std::size_t &pos, const id_type &id)
{
	m_ids[pos] = id;
	m_pos[id]  = pos;
}

/*!
	Updates the id of the specified position to mark it as empty element.

	The id of an empty element contains the distance, measured in
	number of elements, between the current element and the next
	non-empty element (the distance is negative).

	The next used position has to be greater than the position to
	update. Otherwise, undefined behavior occurs.

	\param pos is the position to update
	\param nextUsedPos is the position of the next non-empty element
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::update_empty_pos_id(const std::size_t &pos, const std::size_t &nextUsedPos)
{
	assert(nextUsedPos > pos);

	m_ids[pos] = pos - nextUsedPos;
}

/*!
	Update the first used position.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::update_first_used_pos(const std::size_t &updated_first_pos)
{
	m_first_pos = updated_first_pos;
}

/*!
	Update the last used position.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::update_last_used_pos(const std::size_t &updated_last_pos)
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
	m_holes_regular_begin = std::lower_bound(m_holes_regular_begin, m_holes_regular_end, m_last_pos, std::greater<std::size_t>());
	if (m_holes_regular_begin == m_holes_regular_end) {
		m_holes_regular_begin = m_holes.begin();
		m_holes_regular_end   = m_holes_regular_begin;
	}

	// Remove pending holes beyond the updated last position
	holes_sort_pending();
	m_holes_pending_begin = std::lower_bound(m_holes_pending_begin, m_holes_pending_end, m_last_pos, std::greater<std::size_t>());
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
	Returns the size of the storage expressed in terms of
	elements that the storage contains (including the holes).

	\result The size of the storage.
*/
template<typename T, typename id_type>
std::size_t PiercedVector<T, id_type>::storage_size() const
{
	return m_v.size();
}

/*!
	Resize the storage.

	\param n is the new container size, expressed in number of
	elements.
*/
template<typename T, typename id_type>
void PiercedVector<T, id_type>::storage_resize(size_t n)
{
	std::size_t initialSize = storage_size();
	if (n == initialSize + 1) {
		m_v.emplace_back();
		m_ids.emplace_back();
	} else {
		// Delete the ids of the elements that will be removed
		for (std::size_t pos = n; pos < initialSize; ++pos) {
			id_type id = m_ids[pos];
			if (id >= 0) {
				m_pos.erase(m_ids[pos]);
			}
		}

		// Resize the internal vectors
		m_ids.resize(n);
		m_v.resize(n);
	}
}

/*!
	Order a vector according to a reordering vector.

	\tparam order_t is the type of data that needs to be rodered
	\param v is a reference to the reording vector
	\param v is an iterator tj the vector that will be reordered
*/
template<typename T, typename id_type>
template<typename order_t>
void PiercedVector<T, id_type>::reorder_vector(std::vector<size_t>& order, std::vector<order_t>& v, const size_t &size)
{
	for (size_t i = 0; i < size; i++) {
		size_t j;
		while (i != (j = order[i])) {
			size_t k = order[j];

			order_t temp = std::move(v[j]);
			v[j] = std::move(v[k]);
			v[k] = std::move(temp);

			std::swap(order[i], order[j]);
		}
	}
}

}
