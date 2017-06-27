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

/**
* Constructs an empty pierced vector with no elements.
*
* For increase the performances, the synchronization of the internal storage
* is handled outside the kernel.
*/
template<typename value_t, typename id_t>
PiercedVector<value_t, id_t>::PiercedVector()
    : PiercedKernel<id_t>(),
      PiercedStorage<value_t, id_t>(1, static_cast<PiercedKernel<id_t> *>(this), PiercedKernel<id_t>::SYNC_MODE_DISABLED)
{
}

/**
* Constructs a pierced vector with a capacity at least enough
* to contain n elements.
*
* For increase the performances, the synchronization of the internal storage
* is handled outside the kernel.
*
* \param n the minimum capacity requested for the container
*/
template<typename value_t, typename id_t>
PiercedVector<value_t, id_t>::PiercedVector(std::size_t n)
    : PiercedKernel<id_t>(n),
      PiercedStorage<value_t, id_t>(1, static_cast<PiercedKernel<id_t> *>(this), PiercedKernel<id_t>::SYNC_MODE_DISABLED)
{
}

/**
* Gets an element from a the first position marked as empty and
* assignes to it the specified id. Except for setting the id,
* the element is not modified. Therefore it will still contain
* the data of the element that was previously occupying the
* position or it will be empty if there was no empty position
* and a new element has been created.
*
* \param id is the id that will be assigned to the element
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaim(id_t id)
{
    // Fill a position
    FillAction reclaimAction = PiercedKernel<id_t>::fillHead(id);

    // Reclaim the element associated to the position
    return reclaimValue(reclaimAction);
}

/**
* Gets an element marked as empty and assignes to it the specified
* id. The element will have a position that is between the element
* with the specified reference id and the end of the container.
* Except for setting the id, the element is not modified. Therefore
* it will still contain the data of the element that was previously
* occupying the position or it will be empty if there was no empty
* position and a new element has been created.
* The container is extended by inserting a new element. The element
* will be inserted .
*
* \param referenceId is the id of the element after which an
* empty position will be reclaimed
* \param id is the id that will be assigned to the element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimAfter(const id_t &referenceId, id_t id)
{
    // Fill a position
    FillAction reclaimAction = PiercedKernel<id_t>::fillAfter(referenceId, id);

    // Reclaim the element associated to the position
    return reclaimValue(reclaimAction);
}

/**
* Gets an element from the first position marked as empty
* past the last element assignes to it the specified id.
* Except for setting the id, the element is not modified.
* Therefore it will still contain the data of the element
* that was previously occupying the position or it will be
* empty if there was no empty position and a new element
* has been created.
*
* \param id is the id that will be assigned to the element
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimBack(id_t id)
{
    // Fill a position
    FillAction reclaimAction = PiercedKernel<id_t>::fillAppend(id);

    // Reclaim the element associated to the position
    return reclaimValue(reclaimAction);
}

/**
* Gets an element marked as empty and assignes to it the specified
* id. The element will have a position that is between the begin
* of the container and the specified reference id.
* Except for setting the id, the element is not modified. Therefore
* it will still contain the data of the element that was previously
* occupying the position or it will be empty if there was no empty
* position and a new element has been created.
* The container is extended by inserting a new element. The element
* will be inserted .
*
* \param referenceId is the id of the element before which an
* empty position will be reclaimed
* \param id is the id that will be assigned to the element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimBefore(const id_t &referenceId, id_t id)
{
    // Fill a position
    FillAction reclaimAction = PiercedKernel<id_t>::fillBefore(referenceId, id);

    // Reclaim the element associated to the position
    return reclaimValue(reclaimAction);
}

/**
* Move the specified element before the element with the given
* reference id.
*
* \param referenceId is the id of the element before which the
* new element will be moved
* \param id is the id of the element that will be moved
* \param delayed if true some changes can remain in a pending state
* until a flush is called
* \result An iterator that points to the moved element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::moveBefore(const id_t &referenceId, id_t id, bool delayed)
{
    // Update the position of the element in the kernel
    MoveAction moveAction = PiercedKernel<id_t>::moveBefore(referenceId, id, !delayed);

    // Move the element to the new position
    return moveValue(moveAction);
}

/**
* Move the specified element after the element with the given
* reference id.
*
* \param referenceId is the id of the element after which the
* new element will be moved
* \param id is the id of the element that will be moved
* \param delayed if true some changes can remain in a pending state
* until a flush is called
* \result An iterator that points to the moved element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::moveAfter(const id_t &referenceId, id_t id, bool delayed)
{
    // Update the position of the element in the kernel
    MoveAction moveAction = PiercedKernel<id_t>::moveAfter(referenceId, id, !delayed);

    // Move the element to the new position
    return moveValue(moveAction);
}

/**
* The container is extended by inserting a new element.
*
* \param id is the id that will be associated to the element
* \param value is the value to be copied (or moved) to the
*             inserted elements.
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insert(id_t id, const value_t &value)
{
    // Fill a position
    FillAction insertAction = PiercedKernel<id_t>::fillHead(id);

    // Insert the new value
    return insertValue(insertAction, value);
}

/**
* The container is extended by inserting a new element. The element
* will have a position that is between the element with the
* specified reference id and the end of the container.
*
* \param referenceId is the id of the element after which the
* new element will be inserted
* \param id is the id that will be associated to the element
* \param value is the value to be copied (or moved) to the
* inserted element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insertAfter(const id_t &referenceId, id_t id, const value_t &value)
{
    // Fill a position
    FillAction insertAction = PiercedKernel<id_t>::fillAfter(referenceId, id);

    // Insert the new value
    return insertValue(insertAction, value);
}

/**
* The container is extended by inserting a new element. The element
* will have a position that is between the beginning of the
* container and the element with the specified reference id.
*
* \param referenceId is the id of the element before which the
* new element will be inserted
* \param id is the id that will be associated to the element
* \param value is the value to be copied (or moved) to the
* inserted element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insertBefore(const id_t &referenceId, id_t id, const value_t &value)
{
    // Fill a position
    FillAction insertAction = PiercedKernel<id_t>::fillBefore(referenceId, id);

    // Insert the new value
    return insertValue(insertAction, value);
}

/**
* The element with the specified id is replaced with a new element.
*
* \param id is the id of the element that will be replaced
* \param value is the value to be moved to the inserted elements.
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::replace(id_t id, value_t &&value)
{
    // Position
    std::size_t pos = PiercedKernel<id_t>::getPos(id);

    // Replace the value
    PiercedStorage<value_t, id_t>::rawSet(pos, std::move(value));

    // Return the iterator that points to the element
    return PiercedStorage<value_t, id_t>::getIteratorFromPos(pos);
}

/**
* Adds a new element at the end of the container, after its current last element.
*
* The content of value is copied (or moved) to the new element.
*
* \param id is the id that will be assigned to the element
* \param value the value to be copied (or moved) to the new element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::pushBack(id_t id, const value_t &value)
{
    // Fill a position
    FillAction insertAction = PiercedKernel<id_t>::fillAppend(id);

    // Insert the new value
    return insertValue(insertAction, value);
}

/**
* The container is extended by inserting a new element. This new element is
* constructed in place using args as the arguments for its construction.
*
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
* \result An iterator that points to the the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplace(id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedKernel<id_t>::fillHead(id);

    // Create the new value in-place
    return emplaceValue(emplaceAction, std::forward<Args>(args)...);
}

/**
* The container is extended by inserting a new element. The element
* will have a position that is between the element with the
* specified reference id and the end of the container.
*
* \param referenceId is the id of the element after which the
* new element will be inserted
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplaceAfter(const id_t &referenceId, id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedKernel<id_t>::fillAfter(referenceId, id);

    // Create the new value in-place
    return emplaceValue(emplaceAction, std::forward<Args>(args)...);
}

/**
* Inserts a new element at the end of the container, right after
* its current last element. This new element is constructed
* in place using args as the arguments for its construction.
*
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
*/
template<typename value_t, typename id_t>
template<typename... Args>
void PiercedVector<value_t, id_t>::emplaceBack(id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedKernel<id_t>::fillAppend(id);

    // Create the new value in-place
    emplaceValue(emplaceAction, std::forward<Args>(args)...);
}

/**
* The container is extended by inserting a new element. This new
* element is constructed in place using args as the arguments for
* its construction. The element will have a position that is between
* the beginning of the container and the element with the specified
* reference id.
*
* \param referenceId is the id of the element before which the
* new element will be inserted
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplaceBefore(const id_t &referenceId, id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedKernel<id_t>::fillBefore(referenceId, id);

    // Create the new value in-place
    return emplaceValue(emplaceAction, std::forward<Args>(args)...);
}

/**
* The element with the specified id is replaced with a new element.
* This new element is constructed in place using args as the
* arguments for its construction.
*
* \param id is the id of the element that will be replaced
* \param args are the arguments forwarded to construct the new element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emreplace(id_t id, Args&&... args)
{
    // Position
    std::size_t pos = PiercedKernel<id_t>::getPos(id);

    // Replace the value
    PiercedStorage<value_t, id_t>::rawEmreplace(pos, std::forward<Args>(args)...);

    // Return the iterator that points to the element
    return PiercedStorage<value_t, id_t>::getIteratorFromPos(pos);
}

/**
* Removes from the container the element with the specified id. If the id does
* not exists the function throws an exception.
*
* Element is overwritten with an empty element and the id associated to its
* position is updated to mark the position as empty and allow the container
* to reuse that position.
*
* \param id is the id of the element to erase
* \param delayed if true some changes can remain in a pending state
* until a flush is called
* \result An iterator pointing to the new location of the element that followed
* the element erased by the function call. This is the container end if the
* operation erased the last element in the sequence.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::erase(id_t id, bool delayed)
{
    // Erase the position
    EraseAction eraseAction = PiercedKernel<id_t>::erase(id, !delayed);

    // Erase the value
    return eraseValue(eraseAction);
}

/**
* Removes the last element in the container, effectively reducing
* the container size by one.
*
* Element is not deleted from the internal vector, instead its
* id is changed to mark the position as empty and allow the
* container to reuse that position.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::popBack()
{
    // Erase the position
    EraseAction eraseAction = PiercedKernel<id_t>::popBack();

    // Erase the value
    eraseValue(eraseAction);
}

/**
* Swap the elements with the specified id.
*
* \param id_first is the id of the first element to be swapped
* \param id_second is the id of the second element to be swapped
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::swap(id_t id_first, id_t id_second)
{
    // Update the kernel
    SwapAction swapAction = PiercedKernel<id_t>::swap(id_first, id_second);

    // Update the storage
    swapValues(swapAction);
}

/**
* Removes all elements from the container (which are destroyed), leaving the
* container with a size of 0.
*
* \param release if it's true the memory hold by the container will be
* released, otherwise the container will be cleared but its memory will
* not be relased
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::clear(bool release)
{
    // Update the kernel
    PiercedKernel<id_t>::clear(release);

    // Update the storage
    PiercedStorage<value_t, id_t>::rawClear(release);
}

/**
* Requests that the container capacity be at least enough to contain n elements.
*
* If n is greater than the current container capacity, the function causes the
* container to reallocate its storage increasing its capacity to n (or
* greater).
*
* In all other cases, the function call does not cause a reallocation and the
* container capacity is not affected.
*
* \param n the minimum capacity requested for the container
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::reserve(std::size_t n)
{
    // Update the kernel
    PiercedKernel<id_t>::reserve(n);

    // Update the storgae
    PiercedStorage<value_t, id_t>::rawReserve(n);
}

/**
* Resizes the container so that it contains n elements.
*
* If n is smaller than the current container size, the content is reduced to
* its first n elements, removing those beyond (and destroying them).
*
* If n is greater than the current container size, space is reserved in the
* storage to allow the container to reach the requested size.
*
* If n is also greater than the current container capacity, an automatic
* reallocation of the allocated storage space takes place.
*
* Notice that this function changes the actual content of the container by
* erasing elements from it.
*
* \param n is the new container size, expressed in number of elements.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::resize(std::size_t n)
{
    // Update the kernel
    PiercedKernel<id_t>::resize(n);

    // Update the storage
    PiercedStorage<value_t, id_t>::rawResize(PiercedKernel<id_t>::rawSize());
}

/**
* Sorts the elements of the container in ascending id order.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::sort()
{
    // Update the kernel
    std::vector<std::size_t> permutations = PiercedKernel<id_t>::sort();

    // Update the storage
    PiercedStorage<value_t, id_t>::rawReorder(permutations);
}

/**
* Requests the container to compact the elements and reduce its capacity to
* fit its size.
*
* The request is non-binding, and the function can leave the container with a
* capacity greater than its size.
*
* This may cause a reallocation, but has no effect on the container size and
* cannot alter its elements.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::squeeze()
{
    // Update the kernel
    std::vector<std::size_t> permutations = PiercedKernel<id_t>::squeeze();

    // Update the storage
    PiercedStorage<value_t, id_t>::rawReorder(permutations);
    PiercedStorage<value_t, id_t>::rawResize(PiercedKernel<id_t>::rawSize());
    PiercedStorage<value_t, id_t>::rawShrinkToFit();
}

/**
* Requests the container to reduce its capacity to fit its size. This method
* will NOT compact the elements, leaving the existing holes unaltered.
*
* The request is non-binding, and the function can leave the container with a
* capacity greater than its size.
*
* This may cause a reallocation, but has no effect on the container size and
* cannot alter its elements not the holes.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::shrinkToFit()
{
    // Update the kernel
    PiercedKernel<id_t>::shrinkToFit();

    // Update the storage
    PiercedStorage<value_t, id_t>::rawShrinkToFit();
}

/**
* Exchanges the content of the container by the content of x, which is another
* container of the same type. Sizes may differ.
*
* After the call to this member function, the elements in this container are
* those which were in x before the call, and the elements of x are those
* which were in this. All iterators, references and pointers remain valid
* for the swapped objects.
*
* \param x is another container of the same type (i.e., instantiated with
* the same template parameters) whose content is swapped with that of
* this container.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::swap(PiercedVector &x) noexcept
{
    // Swap kernel data
    PiercedKernel<id_t>::swap(x);

    // Swap storage data
    PiercedStorage<value_t, id_t>::swap(x, false);
}

/**
* Get the kernel of the vector.
*
* \result The kernel of the vector.
*/
template<typename value_t, typename id_t>
const PiercedKernel<id_t> & PiercedVector<value_t, id_t>::getKernel() const
{
    return *this;
}

/**
* Get the storage of the vector.
*
* \result The storage of the vector.
*/
template<typename value_t, typename id_t>
const PiercedStorage<value_t, id_t> & PiercedVector<value_t, id_t>::getStorage() const
{
    return *this;
}

/**
* Dumps to screen the internal data.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::dump()
{
    PiercedKernel<id_t>::dump();
}

/*!
* Gets the raw index associated to the element with the specified id.
*
* \return The raw index associated to the element with the specified id.
*/
template<typename value_t, typename id_t>
std::size_t PiercedVector<value_t, id_t>::rawIndex(id_t id) const
{
    return PiercedKernel<id_t>::getRawIndex(id);
}


/**
* Checks if a given id exists in the container.
*
* \param id is the id to look for
* \result Returns true is the given id exists in the container, otherwise it
* returns false.
*/
template<typename value_t, typename id_t>
bool PiercedVector<value_t, id_t>::exists(id_t id) const
{
    return PiercedKernel<id_t>::contains(id);
}

/**
* Returns a reference to the first element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A reference to the first element of the container.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::front()
{
    if (PiercedKernel<id_t>::empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(PiercedKernel<id_t>::front());
}

/**
* Returns a constant reference to the first element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A constant reference to the first element of the container.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::front() const
{
    if (PiercedKernel<id_t>::empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(PiercedKernel<id_t>::front());
}

/**
* Returns a reference to the last element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A reference to the last element of the container.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::back()
{
    if (PiercedKernel<id_t>::empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(PiercedKernel<id_t>::back());
}

/**
* Returns a constant reference to the last element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A constant reference to the last element of the container.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::back() const
{
    if (PiercedKernel<id_t>::empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(PiercedKernel<id_t>::back());
}

/**
* Returns a constant iterator to the first element with the specified id.
*
* If no such element is found, the function returns a constant iterator
* referring to the past-the-end element of the container.
*
* \param id is the id to look for
* \result Returns a constant iterator to the first element with the
* specified id. If no such element is found, the function returns a
* constant iterator referring to the past-the-end element of the
* container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::const_iterator PiercedVector<value_t, id_t>::find(id_t id) const
{
    if (!PiercedKernel<id_t>::contains(id)) {
        return cend();
    }

    return getConstIterator(id);
}

/**
* Returns an iterator to the first element with the specified id.
*
* If no such element is found, the function returns a constant iterator
* referring to the past-the-end element of the container.
*
* \param id is the id to look for
* \result Returns an iterator to the first element with the specified id.
* If no such element is found, the function returns an iterator referring
* to the past-the-end element of the container.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::find(id_t id)
{
    if (!PiercedKernel<id_t>::contains(id)) {
        return end();
    }

    return getIterator(id);
}

/**
* Returns a reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::at(id_t id)
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Returns a constant reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A constant reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::at(id_t id) const
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Returns a reference to the element at the specified position.
*
* \param pos the position of the element
* \result A reference to the element in the specified position.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::rawAt(std::size_t pos)
{
    return PiercedStorage<value_t, id_t>::rawAt(pos, 0);
}

/**
* Returns a constant reference to the element at the specified position.
*
* \param pos the position of the element
* \result A constant reference to the element in the specified position.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::rawAt(std::size_t pos) const
{
    return PiercedStorage<value_t, id_t>::rawAt(pos, 0);
}

/**
* Returns a constant reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A constant reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PV_CONST_REFERENCE__ PiercedVector<value_t, id_t>::operator[](id_t id) const
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Returns a reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PV_REFERENCE__ PiercedVector<value_t, id_t>::operator[](id_t id)
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Register the specified storage
*
* \param storage is the storage that will be registered
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename value_t, typename id_t>
template<typename data_t>
void PiercedVector<value_t, id_t>::registerStorage(PiercedStorage<data_t, id_t> *storage, PiercedSyncMaster::SyncMode syncMode)
{
    storage->setKernel(this, syncMode);
}

/**
* Unregister the specified storage
*
* \param storage is the storage that will be unregistered
*/
template<typename value_t, typename id_t>
template<typename data_t>
void PiercedVector<value_t, id_t>::unregisterStorage(PiercedStorage<data_t, id_t> *storage)
{
    storage->unsetKernel();
}

/**
* Check if te specified storage is registered.
*
* \param storage is the storage to check
*/
template<typename value_t, typename id_t>
template<typename data_t>
bool PiercedVector<value_t, id_t>::isStorageRegistered(const PiercedStorage<data_t, id_t> *storage) const
{
    return PiercedKernel<id_t>::isSlaveRegistered(storage);
}

/**
* Get the synchronization mode for the specified storage
*
* \param storage is the storage for which the synchronization mode is requested
* \result The synchronization mode of the storage
*/
template<typename value_t, typename id_t>
template<typename data_t>
typename PiercedSyncMaster::SyncMode PiercedVector<value_t, id_t>::getStorageSyncMode(const PiercedStorage<data_t, id_t> *storage) const
{
    return PiercedKernel<id_t>::getSlaveSyncMode(storage);
}

/**
* Reclaim the specified element.
*
* \param action is the fill action that defines how to reclaim the element
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::reclaimValue(const FillAction &action)
{
    std::size_t pos = action.info[PiercedSyncAction::INFO_POS];
    switch (static_cast<typename FillAction::FillActionType>(action.type)) {

    case FillAction::TYPE_OVERWRITE:
    {
        // Nothing to do
        break;
    }

    case FillAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawEmplace(pos);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawEmplaceBack();
        pos = PiercedKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedStorage<value_t, id_t>::getIteratorFromPos(pos);
}

/**
* Insert the specified element.
*
* \param action is the fill action that defines how to insert the element
* \param value is the value that will be assigned to the element
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::insertValue(const FillAction &action, const value_t &value)
{
    std::size_t pos = action.info[PiercedSyncAction::INFO_POS];
    switch (static_cast<typename FillAction::FillActionType>(action.type)) {

    case FillAction::TYPE_OVERWRITE:
    {
        PiercedStorage<value_t, id_t>::rawSet(pos, value);
        break;
    }

    case FillAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawInsert(pos, 1, value);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawPushBack(value);
        pos = PiercedKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedStorage<value_t, id_t>::getIteratorFromPos(pos);
}

/**
* Insert the specified element constructing it in-place.
*
* \param action is the fill action that defines how to insert the element
* \param args are the arguments forwarded to the elements' construct when
* synchronizing the action
*/
template<typename value_t, typename id_t>
template<typename... Args>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emplaceValue(const FillAction &action, Args&&... args)
{
    std::size_t pos = action.info[PiercedSyncAction::INFO_POS];
    switch (static_cast<typename FillAction::FillActionType>(action.type)) {

    case FillAction::TYPE_OVERWRITE:
    {
        PiercedStorage<value_t, id_t>::rawEmreplace(pos, std::forward<Args>(args)...);
        break;
    }

    case FillAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawEmplace(pos, std::forward<Args>(args)...);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawEmplaceBack(std::forward<Args>(args)...);
        pos = PiercedKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedStorage<value_t, id_t>::getIteratorFromPos(pos);
}

/**
* Move the specified positions of the storage.
*
* \param action is the move action that defines how to move the element
* \result An iterator pointing to the new position of the element.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::moveValue(const MoveAction &action)
{
    std::size_t posNew = action.info[PiercedSyncAction::INFO_POS_FIRST];
    std::size_t posOld = action.info[PiercedSyncAction::INFO_POS_SECOND];
    switch (static_cast<typename MoveAction::MoveActionType>(action.type)) {

    case MoveAction::TYPE_MOVE_OVERWRITE:
    {
        PiercedStorage<value_t, id_t>::rawSet(posNew, std::move(PiercedStorage<value_t, id_t>::rawAt(posOld)));
        break;
    }

    case MoveAction::TYPE_MOVE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawInsert(posNew, 1, std::move(PiercedStorage<value_t, id_t>::rawAt(posOld)));
        break;
    }

    case MoveAction::TYPE_MOVE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedStorage<value_t, id_t>::rawPushBack(std::move(PiercedStorage<value_t, id_t>::rawAt(posOld)));
        posNew = PiercedKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Clear the old position
    PiercedStorage<value_t, id_t>::rawEmreplace(posOld);

    // Return the iterator to the new position
    return PiercedStorage<value_t, id_t>::getIteratorFromPos(posNew);
}

/**
* Erases the specified positions of the storage.
*
* The position of the storage will be overwritten with an empty element.
*
* \result An iterator pointing to the new location of the element that
* followed the last element erased.
*/
template<typename value_t, typename id_t>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::eraseValue(const EraseAction &action)
{
    switch (static_cast<typename EraseAction::EraseActionType>(action.type)) {

    case EraseAction::TYPE_PIERCE:
    {
        // Nothing to do. To improve performance the element will not be
        // cleared.
        break;
    }

    case EraseAction::TYPE_SHRINK:
    {
        std::size_t rawSize = action.info[PiercedSyncAction::INFO_SIZE];
        PiercedStorage<value_t, id_t>::rawResize(rawSize);
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the next element
    std::size_t nextPos = action.info[PiercedSyncAction::INFO_POS_NEXT];

    return PiercedStorage<value_t, id_t>::getIteratorFromPos(nextPos);
}

/**
* Swap the specified elements.
*
* \param action is the swap action that defines how to swap the element
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::swapValues(const SwapAction &action)
{
    switch (static_cast<typename SwapAction::SwapActionType>(action.type)) {

    case SwapAction::TYPE_SWAP:
    {
        std::size_t posFirst  = action.info[PiercedSyncAction::INFO_POS_FIRST];
        std::size_t posSecond = action.info[PiercedSyncAction::INFO_POS_SECOND];
        PiercedStorage<value_t, id_t>::rawSwap(posFirst, posSecond);
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }
}

}

#endif
