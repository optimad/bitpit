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

#ifndef __BITPIT_PIERCED_VECTOR_TPP__
#define __BITPIT_PIERCED_VECTOR_TPP__

namespace bitpit {

// Definition of static constants of PiercedVector
template<typename value_t, typename id_t>
const std::size_t
	PiercedVector<value_t, id_t>::MAX_PENDING_HOLES = 16384;

/*!
	Constructs an empty pierced vector with no elements.
*/
template<typename value_t, typename id_t>
PiercedVector<value_t, id_t>::PiercedVector()
{
	clear();
}

/*!
	Constructs a pierced vector with a capacity at least enough
	to contain n elements.

	\param n the minimum capacity requested for the vector
*/
template<typename value_t, typename id_t>
PiercedVector<value_t, id_t>::PiercedVector(std::size_t n)
{
	clear();

	reserve(n);
}

/*!
	Adds a new element at the end of the vector, after its current
	last element.

	The content of value is copied (or moved) to the new element.

	\param id is the id that will be assigned to the element
	\param value the value to be copied (or moved) to the new
					element
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::pushBack(const id_t &id, value_t &&value)
{
	// Fill the position
	StoragePosition storagePosition = fillAppend(id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = std::move(value);
	} else {
		m_v.insert(m_v.begin() + pos, std::move(value));
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaim(const id_t &id)
{
	// Fill the position
	StoragePosition storagePosition = fillHead(id);

	// If necessary insert an empty element
	size_t pos = storagePosition.pos;
	if (!storagePosition.ready) {
		m_v.emplace(m_v.begin() + pos);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimAfter(const id_t &referenceId, const id_t &id)
{
	StoragePosition storagePosition = fillAfter(getPosFromId(referenceId), id);

	// If necessary insert an empty element
	size_t pos = storagePosition.pos;
	if (!storagePosition.ready) {
		m_v.emplace(m_v.begin() + pos);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimBack(const id_t &id)
{
	StoragePosition storagePosition = fillAppend(id);

	// If necessary insert an empty element
	size_t pos = storagePosition.pos;
	if (!storagePosition.ready) {
		m_v.emplace(m_v.begin() + pos);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimBefore(const id_t &referenceId, const id_t &id)
{
	StoragePosition storagePosition = fillBefore(getPosFromId(referenceId), id);

	// If necessary insert an empty element
	size_t pos = storagePosition.pos;
	if (!storagePosition.ready) {
		m_v.emplace(m_v.begin() + pos);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::moveAfter(const id_t &referenceId, const id_t &id, bool delayed)
{
	// Save the element
	std::size_t initialPos = getPosFromId(id);
	value_t temp = std::move(m_v[initialPos]);

	// Pierce the position
	pierce(initialPos, !delayed);

	// Insert the element in the updated position
	return insertAfter(referenceId, id, std::move(temp));
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::moveBefore(const id_t &referenceId, const id_t &id, bool delayed)
{
	// Save the element
	std::size_t initialPos = getPosFromId(id);
	value_t temp = std::move(m_v[initialPos]);

	// Pierce the position
	pierce(initialPos, !delayed);

	// Insert the element in the updated position
	return insertBefore(referenceId, id, std::move(temp));
}

/*!
	The container is extended by inserting a new element.

	\param id is the id that will be associated to the element
	\param value is the value to be copied (or moved) to the
				inserted elements.
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insert(const id_t &id, const value_t &value)
{
	// Fill the position
	StoragePosition storagePosition = fillHead(id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value;
	} else {
		m_v.insert(m_v.begin() + pos, value);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	The container is extended by inserting a new element. The element
	will have a position that is between the element with the
	specified reference id and the end of the container.

	\param referenceId is the id of the element after which the
	new element will be inserted
	\param id is the id that will be associated to the element
	\param value is the value to be copied (or moved) to the
	inserted element
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insertAfter(const id_t &referenceId, const id_t &id, const value_t &value)
{
	// Fill the position
	StoragePosition storagePosition = fillAfter(getPosFromId(referenceId), id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value;
	} else {
		m_v.insert(m_v.begin() + pos, value);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	The container is extended by inserting a new element. The element
	will have a position that is between the beginning of the
	container and the element with the specified reference id.

	\param referenceId is the id of the element before which the
	new element will be inserted
	\param id is the id that will be associated to the element
	\param value is the value to be copied (or moved) to the
	inserted element
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insertBefore(const id_t &referenceId, const id_t &id, const value_t &value)
{
	// Fill the position
	StoragePosition storagePosition = fillBefore(getPosFromId(referenceId), id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value;
	} else {
		m_v.insert(m_v.begin() + pos, value);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	The element with the specified id is replaced with a new element.

	\param id is the id of the element that will be replaced
	\param value is the value to be moved to the inserted elements.
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::replace(id_t id, value_t &&value)
{
	// Position
	size_t pos = m_pos.at(id);

	// Replace the element
	m_v[pos] = std::move(value);

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	Updates the id of the specified element.

	\param currentId is the current id of the element
	\param updatedId is the new id of the element
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::updateId(const id_t &currentId, const id_t &updatedId)
{
	// Validate the id
	validateId(updatedId);

	// Update the id
	setPosId(getPosFromId(currentId), updatedId);
	m_pos.erase(currentId);
}

/*!
	The container is extended by inserting a new element. This
	new element is constructed in place using args as the
	arguments for its construction.

	\param id is the id that will be associated to the element
	\param args the arguments forwarded to construct the new element
	\result An iterator that points to the the newly inserted
			element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplace(const id_t &id, Args&&... args)
{
	// Fill the position
	StoragePosition storagePosition = fillHead(id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value_t(std::forward<Args>(args)...);
	} else {
		m_v.emplace(m_v.begin() + pos, std::forward<Args>(args)...);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	The container is extended by inserting a new element. The element
	will have a position that is between the element with the
	specified reference id and the end of the container.

	\param referenceId is the id of the element after which the
	new element will be inserted
	\param id is the id that will be associated to the element
	\param args the arguments forwarded to construct the new element
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplaceAfter(const id_t &referenceId, const id_t &id, Args&&... args)
{
	// Fill the position
	StoragePosition storagePosition = fillAfter(getPosFromId(referenceId), id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value_t(std::forward<Args>(args)...);
	} else {
		m_v.emplace(m_v.begin() + pos, std::forward<Args>(args)...);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	Inserts a new element at the end of the vector, right after
	its current last element. This new element is constructed
	in place using args as the arguments for its construction.

	\param id is the id that will be associated to the element
	\param args the arguments forwarded to construct the new element
*/
template<typename value_t, typename id_t>
template<typename... Args>
void PiercedVector<value_t, id_t>::emplaceBack(const id_t &id, Args&&... args)
{
	// Fill the position
	StoragePosition storagePosition = fillAppend(id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value_t(std::forward<Args>(args)...);
	} else {
		m_v.emplace(m_v.begin() + pos, std::forward<Args>(args)...);
	}
}

/*!
	The container is extended by inserting a new element. This new
	element is constructed in place using args as the arguments for
	its construction. The element will have a position that is between
	the beginning of the container and the element with the specified
	reference id.

	\param referenceId is the id of the element before which the
	new element will be inserted
	\param id is the id that will be associated to the element
	\param args the arguments forwarded to construct the new element
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplaceBefore(const id_t &referenceId, const id_t &id, Args&&... args)
{
	// Fill the position
	StoragePosition storagePosition = fillBefore(getPosFromId(referenceId), id);

	// Insert the element
	size_t pos = storagePosition.pos;
	if (storagePosition.ready) {
		m_v[pos] = value_t(std::forward<Args>(args)...);
	} else {
		m_v.emplace(m_v.begin() + pos, std::forward<Args>(args)...);
	}

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
}

/*!
	The element with the specified id is replaced with a new element.
	This new element is constructed in place using args as the
	arguments for its construction.

	\param id is the id of the element that will be replaced
	\param args the arguments forwarded to construct the new element
	\result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emreplace(id_t id, Args&&... args)
{
	// Get the position of the element
	size_t pos = m_pos.at(id);

	// Replace the element
	m_v[pos] = value_t(std::forward<Args>(args)...);

	// Return the iterator that points to the element
	return getIteratorFromPos(pos);
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
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::erase(id_t id, bool delayed)
{
	// Position
	size_t pos = m_pos.at(id);

	// Pierce the position
	pierce(pos, !delayed);

	// Return the iterator to the next element
	if (pos + 1 < m_end_pos) {
		return getIteratorFromPos(findNextUsedPos(pos));
	} else {
		return end();
	}
}

/*!
	Removes the last element in the vector, effectively reducing
	the container size by one.

	Element is not deleted from the internal vector, instead its
	id is changed to mark the position as empty and allow the
	container to reuse that position.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::popBack()
{
	if (empty()) {
		throw std::out_of_range("Vector is empty");
	} else if (size() == 1) {
		clear();
	} else {
		std::size_t last_used_pos = findPrevUsedPos(m_end_pos - 1);
		storageShrink(last_used_pos + 1);
	}
}

/*!
	Swap the elements with the specified id.

	\param id_first is the id of the first element to be swapped
	\param id_second is the id of the second element to be swapped
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::swap(const id_t &id_first, const id_t &id_second)
{
	// Positions
	size_t pos_first  = m_pos.at(id_first);
	size_t pos_second = m_pos.at(id_second);

	// Swap the positions
	swapPosIds(pos_first, id_first, pos_second, id_second);

	// Swap the values
	std::swap(m_v[pos_first], m_v[pos_second]);
}

/*!
	Removes all elements from the vector (which are destroyed),
	leaving the container with a size of 0.

	\param release if it's true the memory hold by the container will
	be released, otherwise the container will be cleared but its
	memory will not be relased
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::clear(bool release)
{
	// Clear positions
	m_ids.clear();
	m_pos.clear();
	if (release) {
		std::vector<id_t>().swap(m_ids);
		std::unordered_map<id_t, std::size_t, PiercedHasher>().swap(m_pos);
	}

	// Clear values
	m_v.clear();
	if (release) {
		std::vector<value_t>().swap(m_v);
	}

	// Reset begin and end
	setBeginPos(0);
	setEndPos(0);

	// Clear holes
	holesClear(release);

	// There are no dirty positions
	m_dirty_begin_pos = m_end_pos;
}

/*!
	Flush all pending changes.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::flush()
{
	// Flush pending holes
	holesFlush();
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
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::reserve(std::size_t n)
{
	m_ids.reserve(n);
	m_pos.reserve(n);
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
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::resize(std::size_t n)
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
	holesFlush();

	// Find the id of the last element
	id_t last_stored_id = getSizeMarker(n - 1);

	// Find the last position
	std::size_t last_used_pos = getPosFromId(last_stored_id);

	// Shrink the storage
	storageShrink(last_used_pos + 1);
}

/*!
	Sorts the elements of the vector in ascending id order.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::sort()
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
	for (size_t i = 0; i < containerSize; ++i) {
		id_permutation[i] = i;
	}
	std::sort(id_permutation.begin(), id_permutation.end(), idLess(m_ids));

	std::vector<std::size_t> value_permutation(id_permutation);

	// Sort the container
	utils::reorderVector<id_t>(id_permutation, m_ids, containerSize);
	utils::reorderVector<value_t>(value_permutation, m_v, containerSize);
}

/*!
	Requests the container to compact the elements and reduce
	its capacity to fit its size.

	The request is non-binding, and the function can leave
	the vector with a capacity greater than its size.

	This may cause a reallocation, but has no effect on the vector
	size and cannot alter its elements.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::squeeze()
{
	// Flush changes
	flush();

	// Compact the vector
	std::size_t nHoles = holesCount();
	if (nHoles != 0) {
		// Move the elements
		std::size_t firstPosToUpdate;
		if (m_begin_pos == 0) {
			firstPosToUpdate = *(m_holes_regular_end - 1);
		} else {
			firstPosToUpdate = 0;
		}

		std::size_t offset = 0;
		for (std::size_t pos = firstPosToUpdate; pos < m_end_pos; pos++) {
			if (offset < nHoles && *(m_holes_regular_end - offset - 1) == pos) {
				++offset;
				continue;
			}

			id_t id = m_ids[pos];
			std::size_t updatedPos = pos - offset;

			setPosId(updatedPos, id);
			setPosEmptyId(pos, pos + 1);
			m_v[updatedPos] = std::move(m_v[pos]);
		}

		// Clear the holes
		holesClear();

		// Reset begin and end
		setBeginPos(0);
		setEndPos(size());

		// Shrink the container
		storageShrink(size(), true);
	}

	// Shrink to fit
	shrinkToFit();
}

/*!
	Requests the container to reduce its capacity to fit its size. This
	method will NOT compact the elements, leaving the existing holes
	unaltered.

	The request is non-binding, and the function can leave
	the vector with a capacity greater than its size.

	This may cause a reallocation, but has no effect on the vector
	size and cannot alter its elements not the holes.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::shrinkToFit()
{
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
				instantiated with the same template parameters, value_t and
				Alloc) whose content is swapped with that of this
				container.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::swap(PiercedVector& x) noexcept
{
	std::swap(x.m_begin_pos, m_begin_pos);
	std::swap(x.m_end_pos, m_end_pos);
	std::swap(x.m_dirty_begin_pos, m_dirty_begin_pos);
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::capacity() const
{
	return m_v.capacity();
}

/*!
	Returns whether the vector is contiguous (i.e. whether it contains
	no holes).

	\result true if the container is contiguous, false otherwise.
*/
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::contiguous() const
{
	return (holesCount() == 0);
}

/*!
	Dumps to screen the internal data.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::dump()
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
	std::cout << " m_begin_pos: " << m_begin_pos << std::endl;
	std::cout << " m_end_pos: " <<  m_end_pos << std::endl;
	std::cout << " Stored ids: " << std::endl;
	if (m_ids.size() > 0) {
		for (size_t k = 0; k < m_end_pos; ++k) {
			std::cout << m_ids[k] << std::endl;
		}
	} else {
		std::cout << "None" << std::endl;
	}

	std::cout << std::endl;
	std::cout << " Poistion map: " << std::endl;
	if (m_pos.size() > 0) {
		for (auto itr = m_pos.cbegin(); itr != m_pos.cend(); ++itr) {
			std::cout << itr->first << " -> " << itr->second << std::endl;
		}
	} else {
		std::cout << "None" << std::endl;
	}

	std::cout << "----------------------------------------" << std::endl;
}

/*!
	Returns whether the vector is empty (i.e. whether its size
	is 0).

	\result true if the container size is 0, false otherwise.
*/
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::empty() const
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
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::isIteratorSlow()
{
	return (m_dirty_begin_pos < m_end_pos);
}

/*!
	Returns the maximum number of elements that the vector can hold.

	This is the maximum potential size the container can reach due
	to known system or library implementation limitations, but the
	container is by no means guaranteed to be able to reach that
	size: it can still fail to allocate storage at any point before
	that size is reached.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::maxSize() const
{
	return m_v.max_size();
}

/*!
	Returns the number of elements in the vector.

	This is the number of actual objects held in the vector,
	which is not necessarily equal to its storage capacity.

	\result The number of elements in the container.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::size() const
{
	return m_pos.size();
}

/*!
	Checks if a given id exists in the vector.

	\param id the id to look for
	\result Returns true is the given id exists in the vector,
			otherwise it returns false.
*/
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::exists(id_t id) const
{
	return (m_pos.count(id) != 0);
}

/*!
	Returns a constant iterator to the first element with the specified id.
	If no such element is found, the function returns a constant iterator
	referring to the past-the-end element of the container.

	\param id the id to look for
	\result Returns a constant iterator to the first element with the
	specified id. If no such element is found, the function returns a
	constant iterator referring to the past-the-end element of the
	container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::find(id_t id) const
{
	if (!exists(id)) {
		return cend();
	}

	return getConstIterator(id);
}

/*!
	Returns an iterator to the first element with the specified id.
	If no such element is found, the function returns a constant iterator
	referring to the past-the-end element of the container.

	\param id the id to look for
	\result Returns an iterator to the first element with the specified id.
	If no such element is found, the function returns an iterator referring
	to the past-the-end element of the container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::find(id_t id)
{
	if (!exists(id)) {
		return end();
	}

	return getIterator(id);
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::evalFlatIndex(id_t id)
{
	size_t pos  = getPosFromId(id);

	// Initialize flat id with the position of the element
	size_t flat = pos;

	// Subtract pending holes before position
	if (holesCountPending() > 0) {
		holesSortPending();
		auto hole_itr = std::upper_bound(m_holes_pending_begin, m_holes_pending_end, pos, std::greater<std::size_t>());
		size_t nHolesBefore = std::distance(hole_itr, m_holes_pending_end);

		flat -= nHolesBefore;
	}

	// Subtract regular holes before position
	if (holesCountRegular() > 0) {
		holesSortRegular();
		auto hole_itr = std::upper_bound(m_holes_regular_begin, m_holes_regular_end, pos, std::greater<std::size_t>());
		size_t nHolesBefore = std::distance(hole_itr, m_holes_regular_end);

		flat -= nHolesBefore;
	}

	// Done
	return flat;
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::getRawIndex(id_t id) const
{
	return getPosFromId(id);
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
template<typename value_t, typename id_t>
std::vector<id_t> PiercedVector<value_t, id_t>::getIds(bool ordered) const
{
	size_t nIds = size();

	// Initialize the vector
	std::vector<id_t> ids;
	if (nIds == 0) {
		return ids;
	}

	// Resize the vector
	ids.resize(nIds);

	// Extract the ids
	size_t n   = 0;
	size_t pos = m_begin_pos;
	while (true) {
		ids[n] = m_ids[pos];
		if (n == nIds - 1) {
			break;
		}

		n++;
		pos = findNextUsedPos(pos);
	}

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
template<typename value_t, typename id_t>
id_t PiercedVector<value_t, id_t>::getSizeMarker(const size_t &targetSize, const id_t &fallback)
{
	// If the size is zero, we return the first element, if the target
	// size is equal to the size minus one we return the last element,
	// if the target size is greater or equal the current container size
	// we return the fallback value.
	if (targetSize >= size()) {
		return fallback;
	} else if (targetSize == 0) {
		return m_ids[m_begin_pos];
	} else if (targetSize == (size() - 1)) {
		return m_ids[m_end_pos - 1];
	}

	// Sort the holes
	holesSortRegular();
	holesSortPending();

	// Iterate to find the position before wihch there is the
	// requeste number of element.
	hole_iterator regular_hole_itr = m_holes_regular_end;
	hole_iterator pending_hole_itr = m_holes_pending_end;

	std::size_t nEmpties  = 0;
	std::size_t markerPos = targetSize;
	while (true) {
		if (isPosEmpty(markerPos)) {
			markerPos = findNextUsedPos(markerPos - 1);
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
template<typename value_t, typename id_t>
__PV_POINTER__ PiercedVector<value_t, id_t>::data() noexcept
{
	return m_v.data();
}

/*!
	Returns a reference to the last element of the vector. If
	the vector is empty, an exception is thrown.

	\result A reference to the last element of the vector.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::back()
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	}

	return m_v[m_end_pos - 1];
}

/*!
	Returns a constant reference to the last element of the vector.
	If the vector is empty, an exception is thrown.

	\result A constant reference to the last element of the vector.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::back() const
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	}

	return m_v[m_end_pos - 1];
}

/*!
	Returns a reference to the first element of the vector. If
	the vector is empty, an exception is thrown.

	\result A reference to the first element of the vector.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::front()
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	}

	return m_v[m_begin_pos];
}

/*!
	Returns a constant reference to the first element of the vector.
	If the vector is empty, an exception is thrown.

	\result A constant reference to the first element of the vector.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::front() const
{
	if (empty()) {
		throw std::out_of_range ("Vector is empty");
	}

	return m_v[m_begin_pos];
}

/*!
	Returns a reference to the element with the specified id. If
	there is no element with the specified id, an exception is
	thrown.

	\param id the id of the element
	\result A reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::at(const id_t &id)
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
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::at(const id_t &id) const
{
	return (*this)[id];
}

/*!
	Returns a reference to the element at the specified position.

	\param pos the position of the element
	\result A reference to the element in the specified position.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::rawAt(const std::size_t &pos)
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
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::rawAt(const std::size_t &pos) const
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::rawIndex(id_t id) const
{
	return getPosFromId(id);
}

/*!
	Returns a constant reference to the element with the
	specified id. If there is no element with the specified id,
	an exception is thrown.

	\param id the id of the element
	\result A constant reference to the element with the specified
			id.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::operator[](const id_t &id) const
{
	std::size_t pos = getPosFromId(id);

	return m_v[pos];
}

/*!
	Returns a reference to the element with the
	specified id. If there is no element with the specified id,
	an exception is thrown.

	\param id the id of the element
	\result A reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::operator[](const id_t &id)
{
	std::size_t pos = getPosFromId(id);

	return m_v[pos];
}

/*!
	Gets an iterator pointing to the specified element.

	\param id is the id of the specified iterator.
	\result An iterator pointing to the specified element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::getIterator(const id_t &id) noexcept
{
	const size_t pos = getPosFromId(id);

	return getIteratorFromPos(pos);
}

/*!
	Gets a constant iterator pointing to the specified element.

	\param id is the id of the specified iterator.
	\result A constant iterator pointing to the specified element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::getConstIterator(const id_t &id) const noexcept
{
	const size_t pos = getPosFromId(id);

	return getConstIteratorFromPos(pos);
}

/*!
	Gets an iterator pointing to the specified position.

	\param id is the id of the specified iterator.
	\result An iterator pointing to the specified position.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::getIteratorFromRawIndex(const std::size_t &rawIndex) noexcept
{
	return getIteratorFromPos(rawIndex);
}

/*!
	Gets a constant iterator pointing to the specified position.

	\param id is the id of the specified iterator.
	\result A constant iterator pointing to the specified position.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::getConstIteratorFromRawIndex(const std::size_t &rawIndex) const noexcept
{
	return getConstIteratorFromPos(rawIndex);
}

/*!
	Returns an iterator pointing to the first element in the
	vector.

	\result An iterator pointing to the first element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::begin() noexcept
{
	return getIteratorFromPos(m_begin_pos);
}

/*!
	Returns an iterator referring to the past-the-end element
	in the vector.

	\result An iterator referring to the past-the-end element
			in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::end() noexcept
{
	return getIteratorFromPos(m_end_pos);
}

/*!
	Returns a constant iterator pointing to the first element
	in the vector.

	\result A constant iterator pointing to the first element in
	the vector.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::begin() const noexcept
{
	return cbegin();
}

/*!
	Returns a constant iterator referring to the past-the-end
	element in the vector.

	\result A constant iterator referring to the past-the-end
	element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::end() const noexcept
{
	return cend();
}

/*!
	Returns an conts_iterator pointing to the first element in the
	vector.

	\result A const_iterator pointing to the first element in
			the vector.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::cbegin() const noexcept
{
	return getConstIteratorFromPos(m_begin_pos);
}

/*!
	Returns an const_iterator referring to the past-the-end element
	in the vector.

	\result A const_iterator referring to the past-the-end element
			in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::cend() const noexcept
{
	return getConstIteratorFromPos(m_end_pos);
}

/*!
	Returns an iterator pointing to the first element in the
	raw container.

	\result An iterator pointing to the first element in the raw
			container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::raw_iterator PiercedVector<value_t, id_t>::rawBegin() noexcept
{
	return m_v.begin();
}

/*!
	Returns an iterator referring to the past-the-end element
	in the raw container.

	\result An iterator referring to the past-the-end element
			in the raw container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::raw_iterator PiercedVector<value_t, id_t>::rawEnd() noexcept
{
	return m_v.end();
}

/*!
	Returns a constant iterator pointing to the first element
	in the raw container.

	\result A constant iterator pointing to the first element in
	the raw container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::raw_const_iterator PiercedVector<value_t, id_t>::rawBegin() const noexcept
{
	return rawCbegin();
}

/*!
	Returns a constant iterator referring to the past-the-end
	element in the raw container.

	\result A constant iterator referring to the past-the-end
	element in the raw container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::raw_const_iterator PiercedVector<value_t, id_t>::rawEnd() const noexcept
{
	return rawCend();
}

/*!
	Returns an conts_iterator pointing to the first element in the
	raw container.

	\result A const_iterator pointing to the first element in
			the raw container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::raw_const_iterator PiercedVector<value_t, id_t>::rawCbegin() const noexcept
{
	return m_v.cbegin();
}

/*!
	Returns an const_iterator referring to the past-the-end element
	in raw container.

	\result A const_iterator referring to the past-the-end element
			in raw container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::raw_const_iterator PiercedVector<value_t, id_t>::rawCend() const noexcept
{
	return m_v.cend();
}

/*!
	Gets an iterator pointing to the element in the specified position.

	\param pos is the position of the element
	\result An iterator pointing to the element in the specified position.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::getIteratorFromPos(const std::size_t &pos) noexcept
{
	return iterator(this, pos);
}

/*!
	Gets a constant iterator pointing to the element in the specified
	position.

	\param pos is the position of the element
	\result A constant iterator pointing to the element in the specified
	position.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::getConstIteratorFromPos(const std::size_t &pos) const noexcept
{
	return const_iterator(this, pos);
}

/*!
	Fills a position and assigns to it the specified id.

	The position filled will be the first position available.

	\param id is the id that will be associated to the position
	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillHead(const id_t &id)
{
	// Check if there is a pending hole that can be filled
	if (m_holes_pending_begin != m_holes_pending_end) {
		// Sort pending holes
		holesSortPending();

		// The last hole is the one with the lowest position
		return fillHole(m_holes_pending_end - 1, id);
	}

	// Check if there is a regular hole that can be filled
	if (m_holes_regular_begin != m_holes_regular_end) {
		// Sort regular holes
		holesSortRegular();

		// The last hole is the one with the lowest position
		return fillHole(m_holes_regular_end - 1, id);
	}

	// There are no holes that can be filled
	return fillAppend(id);
}

/*!
	Fills a position and assigns to it the specified id.

	The position filled will be the first position available.

	\param id is the id that will be associated to the position
	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillTail(const id_t &id)
{
	// Check if there is a pending hole that can be filled
	if (m_holes_pending_begin != m_holes_pending_end) {
		// Sort pending holes
		holesSortPending();

		// The first hole is the one with the highest position
		return fillHole(m_holes_pending_begin, id);
	}

	// Check if there is a regular hole that can be filled
	if (m_holes_regular_begin != m_holes_regular_end) {
		// Sort regular holes
		holesSortRegular();

		// The first hole is the one with the lowest position
		return fillHole(m_holes_regular_begin, id);
	}

	// There are no holes that can be filled
	return fillAppend(id);
}

/*!
	Fills a position and assigns to it the specified id.

	The position filled will be before the specified reference position.

	\param referencePos is the position of the element before which the
	new available position will be searched for
	\param id is the id that will be associated to the position
	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillBefore(const std::size_t &referencePos, const id_t &id)
{
	// Check if there is a pending hole that can be filled
	if (m_holes_pending_begin != m_holes_pending_end) {
		// Sort pending holes
		holesSortPending();

		// The last hole is the one with the lowest position, if the
		// position of this hole is greater than the reference position
		// it is possible to use it.
		hole_iterator holeItr = m_holes_pending_end - 1;
		if (*holeItr < referencePos) {
			return fillHole(holeItr, id);
		}
	}

	// Check if there is a regular hole that can be filled
	if (m_holes_regular_begin != m_holes_regular_end) {
		// Sort regular holes
		holesSortRegular();

		// The last hole is the one with the lowest position, if the
		// position of this hole is greater than the reference position
		// it is possible to use it.
		hole_iterator holeItr = m_holes_regular_end - 1;
		if (*holeItr < referencePos) {
			return fillHole(holeItr, id);
		}
	}

	// There are no holes that can be filled
	return fillInsert(referencePos, id);
}

/*!
	Fills a position and assigns to it the specified id.

	The position filled will be after the specified reference position.

	\param referencePos is the position of the element after which the
	new available position will be searched for
	\param id is the id that will be associated to the position
	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillAfter(const std::size_t &referencePos, const id_t &id)
{
	// Check if there is a pending hole that can be filled
	if (m_holes_pending_begin != m_holes_pending_end) {
		// Sort pending holes
		holesSortPending();

		// The first hole is the one with the highest position, if the
		// position of this hole is greater than the reference position
		// it is possible to use it.
		hole_iterator holeItr = m_holes_pending_begin;
		if (*holeItr > referencePos) {
			return fillHole(holeItr, id);
		}
	}

	// Check if there is a regular hole that can be filled
	if (m_holes_regular_begin != m_holes_regular_end) {
		// Sort regular holes
		holesSortRegular();

		// The first hole is the one with the highest position, if the
		// position of this hole is greater than the reference position
		// it is possible to use it.
		hole_iterator holeItr = m_holes_regular_begin;
		if (*holeItr > referencePos) {
			return fillHole(holeItr, id);
		}
	}

	// There are no holes that can be filled
	return fillAppend(id);
}

/*!
	Fills the specified hole with the given id.

	\param pos is the position to fill
	\param id is the id that will be associated to the position
	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillHole(const hole_iterator &holeItr, const id_t &id)
{
	// Validate the id
	validateId(id);

	// Get the position of the hole
	std::size_t pos = *holeItr;

	// Fill the hole
	setPosId(pos, id);

	// Remove the hole from the container
	if (holeItr >= m_holes_pending_begin) {
		int nPendings = holesCountPending();
		if (nPendings > 1) {
			if (holeItr == (m_holes_pending_end - 1)) {
				--m_holes_pending_end;
			} else if (holeItr == m_holes_pending_begin) {
				++m_holes_pending_begin;
			} else {
				throw std::out_of_range ("Only holes at the beginning or at the end of the container can be filled");
			}
		} else {
			holesClearPending();
		}
	} else {
		int nRegulars = holesCountRegular();
		if (nRegulars > 1) {
			if (holeItr == (m_holes_regular_end - 1)) {
				--m_holes_regular_end;
			} else if (holeItr == m_holes_regular_begin) {
				++m_holes_regular_begin;
			} else {
				throw std::out_of_range ("Only holes at the beginning or at the end of the container can be filled");
			}
		} else {
			holesClearRegular();
		}
	}

	// Update the begin position
	//
	// There are no holes past the end of the vector, this means that only
	// the begin position may need an update.
	if (pos < m_begin_pos) {
		setBeginPos(pos);
	}

	// Update the position of the empty elements before the current one
	//
	// It is necessary to ensure that the empty elements just before the
	// current one are correctly set, otherwise it will not be possible
	// to iterate through the container.
	size_t previousPos = pos;
	while (previousPos > 0) {
		--previousPos;
		if (!isPosEmpty(previousPos)) {
			break;
		}

		id_t previousId = m_ids[previousPos];
		if (previousId == -1 || (size_t) std::abs(previousId) == (pos - previousPos)) {
			break;
		}

		setPosEmptyId(previousPos, pos);
	}

	// Return the filled location
	return StoragePosition(pos, true);
}

/*!
	Fills a position and assigns to it the specified id.

	A new position is created at the end of the container.

	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillAppend(const id_t &id)
{
	// Check if the id is valid
	validateId(id);

	// Add the id
	m_ids.emplace_back();
	id_t &storedId = m_ids.back();
	storedId = id;

	// Update last used position
	setEndPos(m_ids.size());

	// Update the id map
	m_pos[id] = m_end_pos - 1;

	// Return the filled location
	return StoragePosition(m_end_pos - 1, false);
}

/*!
	Fills a position and assigns to it the specified id.

	A new postion is created at the specified position.

	\param pos is the position at which the id will be inserted
	\param id is the id that will be associated to the position
	\result The position that has been filled.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::StoragePosition PiercedVector<value_t, id_t>::fillInsert(const std::size_t &pos, const id_t &id)
{
	// If the position at which we want to insert an element is a hole there
	// was an error somewhere. Before inserting a new position the hole needs
	// to be filled first.
	if (isPosEmpty(pos)) {
		throw std::out_of_range ("Before inserting a new position the hole needs to be filled first");
	}

	// We cannot insert elements past the last position
	if (pos >= m_end_pos) {
		throw std::out_of_range ("Unable to insert elements past the last position");
	}

	// Check if the id is valid
	validateId(id);

	// Add the id
	m_ids.emplace(m_ids.begin() + pos, id);

	// Update last used position
	setEndPos(m_ids.size());

	// Update the id map
	for (size_t i = pos + 1; i < m_end_pos; ++i) {
		id_t id_i = m_ids[i];
		if (id_i >= 0) {
			m_pos[id_i] = i;
		}
	}
	m_pos[id] = pos;

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

	// Return the filled location
	return StoragePosition(pos, false);
}

/*!
	Marks a position as empty.

	The position is inserted in the list of holes. If the list of pending
	holes is full, a flush is called before adding the hole. This means
	that the hole is always added as a pending hole.

	\param pos is the position of the new hole
	\param flush controls if the holes will be flush after piercing the
	the position
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::pierce(const std::size_t &pos, bool flush)
{
	// If removing the last position, there is no need to add the
	// position to the holes, it's enough to update the last position
	// counter or clear the container if this was the last hole.
	if (pos + 1 == m_end_pos) {
		if (size() == 1) {
			clear();
		} else {
			std::size_t last_used_pos = findPrevUsedPos(m_end_pos - 1);
			storageShrink(last_used_pos + 1);
		}
		return;
	}

	// Remove the id from the map
	id_t id = m_ids[pos];
	m_pos.erase(id);

	// Reset the element
	m_v[pos] = value_t();

	// Reset the position
	size_t nextUsedPos = findNextUsedPos(pos);
	setPosEmptyId(pos, nextUsedPos);
	m_dirty_begin_pos = std::min(pos, m_dirty_begin_pos);

	// If removing the first position, update the counter
	if (pos == m_begin_pos) {
		std::size_t begin = findNextUsedPos(m_begin_pos);
		setBeginPos(begin);
	}

	// If the list of pending holes is full, flush the holes.
	if (m_holes_pending_end == m_holes.end()) {
		holesFlush();
	}

	// Add the hole at the end of the pending holes
	*m_holes_pending_end = pos;
	m_holes_pending_end++;

	// Check if pending holes are still sorted
	if (m_holes_pending_sorted) {
		size_t nPendings = holesCountPending();
		if (nPendings > 1 && (*(m_holes_pending_end - 1) > *(m_holes_pending_end - 2))) {
			m_holes_pending_sorted = false;
		}
	}

	// Flush
	if (flush) {
		holesFlush();
	}
}

/*!
	Clear the list of available holes.

	\param relase if set to true the memory previously hold by holes'
	container will be released
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesClear(bool release)
{
	// Clear the vector
	holesResize(0, 0, 0, release);

	// There are no holes, therefore all holes are sorted
	m_holes_regular_sorted = true;
	m_holes_pending_sorted = true;
}

/*!
	Clear regular holes.

	\param relase if set to true the memory previously hold by holes'
	container will be released
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesClearRegular(bool release)
{
	// Release the memory
	if (release) {
		m_holes.shrink_to_fit();
	}

	// Reset regulr holes iterators
	m_holes_regular_begin = m_holes.begin();
	m_holes_regular_end   = m_holes_regular_begin;

	// There are no holes, therefore all holes are sorted
	m_holes_regular_sorted = true;
}

/*!
	Clear pending holes

	\param relase if set to true the memory previously hold by holes'
	container will be released
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesClearPending(bool release)
{
	// Clear section of the container associated with the pending holes
	long offset    = std::distance(m_holes.begin(), m_holes_regular_begin);
	long nRegulars = holesCountRegular();
	holesResize(offset, nRegulars, 0, release);

	// There are no pending holes, therefore pending holes are sorted
	m_holes_pending_sorted = true;
}

/*!
	Resize the container of the pending holes

	\param offset is the distance between the first regular hole and the
	begin of the hole's container
	\param nRegulars is the number of regulars holes
	\param nPendings  the number of pending holes
	\param relase if set to true the memory previously hold by holes'
	container will be released
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesResize(size_t offset, size_t nRegulars, size_t nPendings, bool release)
{
	if (release) {
		m_holes.shrink_to_fit();
	}

	m_holes.resize(offset + nRegulars + MAX_PENDING_HOLES);

	m_holes_regular_begin = m_holes.begin() + offset;
	m_holes_regular_end   = m_holes_regular_begin + nRegulars;
	m_holes_pending_begin = m_holes_regular_end;
	m_holes_pending_end   = m_holes_pending_begin + nPendings;
}

/*!
	Count the available holes.

	\result The number of available holes.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::holesCount() const
{
	return holesCountPending() + holesCountRegular();
}

/*!
	Count the pending holes.

	\result The number of pending holes.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::holesCountPending() const
{
	return std::distance(m_holes_pending_begin, m_holes_pending_end);
}

/*!
	Count the regular holes.

	\result The number of regular holes.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::holesCountRegular() const
{
	return std::distance(m_holes_regular_begin, m_holes_regular_end);
}

/*!
	Flushes the list of available holes.

	All the pending hole are converted to regular holes and new
	space is reserved for future pending holes.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesFlush()
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
	holesSortPending();

	auto itr = m_holes_pending_begin;
	size_t pos = m_end_pos;
	do {
		if (*itr >= pos) {
			itr++;
			continue;
		}

		pos = *itr;
		size_t next_used_pos = findNextUsedPos(pos);
		do {
			setPosEmptyId(pos, next_used_pos);
			if (pos > 0) {
				pos--;
			} else {
				break;
			}
		} while (isPosEmpty(pos));
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
				std::size_t nRegulars = holesCountRegular();
				if (nRegulars > 1 && (*(m_holes_regular_end - 1) > *(m_holes_regular_end - 2))) {
					m_holes_regular_sorted = false;
				}
			}
		}
	}

	// Move the holes at the beginning of the vector
	std::size_t nRegulars = holesCountRegular();
	if (nRegulars != 0 && m_holes_regular_begin != m_holes.begin()) {
		std::size_t offset = std::distance(m_holes.begin(), m_holes_regular_begin);
		for (std::size_t k = 0; k < nRegulars; ++k) {
			m_holes[k] = m_holes[k + offset];
		}

		m_holes_regular_begin = m_holes.begin();
		m_holes_regular_end   = m_holes_regular_begin + nRegulars;
	}

	// Resize the vector
	holesClearPending();

	// There are no more dirty positions
	m_dirty_begin_pos = m_end_pos;
}


/*!
	Sort the list of pending holes in descendent order
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesSortPending()
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
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::holesSortRegular()
{
	if (m_holes_regular_sorted) {
		return;
	}

	std::sort(m_holes_regular_begin, m_holes_regular_end, std::greater<std::size_t>());
	m_holes_regular_sorted = true;
}

/*!
	Validates the specified id.

	\param id is the id to validate
	\result Return true if the id is valid, otherwise it throws an exception.
*/
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::validateId(const id_t &id)
{
	// Ids needs to be positive
	if (id < 0) {
		throw std::out_of_range ("Negative id");
	}

	// Handle duplicate ids
	if (exists(id)) {
		throw std::out_of_range ("Duplicate id");
	}

	// Id is valid
	return true;
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::findPrevUsedPos(std::size_t pos) const
{
	std::size_t prev_pos = pos;
	while (true) {
		if (prev_pos == m_begin_pos) {
			throw std::out_of_range ("Already in the firts position");
		}
		prev_pos--;

		id_t prev_id = m_ids[prev_pos];
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::findNextUsedPos(std::size_t pos) const
{
	std::size_t next_pos   = pos;
	std::size_t next_delta = 1;
	while (true) {
		if (next_pos + 1 == m_end_pos) {
			throw std::out_of_range ("Already in the last position");
		}
		next_pos += next_delta;

		id_t next_id = m_ids[next_pos];
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
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::isPosEmpty(std::size_t pos) const
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
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::getPosFromId(id_t id) const
{
	return m_pos.at(id);
}

/*!
	Associate a position to the specified id.

	\param pos is the position to associate
	\param id is the id that will be associated to the position
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::setPosId(const std::size_t &pos, const id_t &id)
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
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::setPosEmptyId(const std::size_t &pos, const std::size_t &nextUsedPos)
{
	assert(nextUsedPos > pos);

	m_ids[pos] = pos - nextUsedPos;
}

/*!
	Swaps two positions.

	\param pos_1 is the first position to swap
	\param id_1 is the id associated to the first position
	\param pos_2 is the second position to swap
	\param id_2 is the id associated to the second position
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::swapPosIds(const std::size_t &pos_1, const id_t &id_1, const std::size_t &pos_2, const id_t &id_2)
{
	std::swap(m_ids[pos_1], m_ids[pos_2]);
	std::swap(m_pos[id_1], m_pos[id_2]);
}

/*!
	Set the first used position.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::setBeginPos(const std::size_t &pos)
{
	m_begin_pos = pos;
}

/*!
	Set the last used position.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::setEndPos(const std::size_t &pos)
{
	m_end_pos = pos;
}

/*!
	Shrink the container so that it contains n raw position.

	\param n is the new container storage size, expressed in number of raw
	positions.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::storageShrink(std::size_t n, bool force)
{
	size_t initialStorageSize = storageSize();

	// We can only shrink the container
	if (n > initialStorageSize) {
		throw std::out_of_range ("The container can only be shrunk");
	}

	// Check if we actually need to shrink the container
	if (n == initialStorageSize && !force) {
		return;
	}

	// When the new last position is before the first one this is equivalent
	// to a clear
	if (n < (m_begin_pos + 1)) {
		clear();
		return;
	}

	// Delete the ids of the elements that will be removed
	for (std::size_t pos = n; pos < initialStorageSize; ++pos) {
		id_t id = m_ids[pos];
		if (id >= 0) {
			m_pos.erase(id);
		}
	}

	// Resize the internal vectors
	m_ids.resize(n);
	m_v.resize(n);

	// Update the last position
	setEndPos(n);

	// If we don't need to update the holes we can exit now
	if (holesCount() == 0) {
		return;
	}

	// Remove regular holes beyond the updated last position
	holesSortRegular();
	m_holes_regular_begin = std::lower_bound(m_holes_regular_begin, m_holes_regular_end, m_end_pos - 1, std::greater<std::size_t>());
	if (m_holes_regular_begin == m_holes_regular_end) {
		m_holes_regular_begin = m_holes.begin();
		m_holes_regular_end   = m_holes_regular_begin;
	}

	// Remove pending holes beyond the updated last position
	holesSortPending();
	m_holes_pending_begin = std::lower_bound(m_holes_pending_begin, m_holes_pending_end, m_end_pos - 1, std::greater<std::size_t>());
	if (m_holes_pending_begin == m_holes_pending_end) {
		m_holes_pending_begin = m_holes_regular_end;
		m_holes_pending_end   = m_holes_pending_begin;
	}
}

/*!
	Returns the size of the storage expressed in terms of
	elements that the storage contains (including the holes).

	\result The size of the storage.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::storageSize() const
{
	return m_v.size();
}

}

#endif
