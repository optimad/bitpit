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
    : PiercedVectorKernel<id_t>(),
      PiercedVectorStorage<value_t, id_t>(1, this, PiercedVectorKernel<id_t>::SYNC_MODE_DISABLED)
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
    : PiercedVectorKernel<id_t>(n),
      PiercedVectorStorage<value_t, id_t>(1, this, PiercedVectorKernel<id_t>::SYNC_MODE_DISABLED)
{
}

/**
* Copy constructor
*
* \param x is another container of the same type (i.e., instantiated with
* the same template parameters) whose content is copied in this container.
*/
template<typename value_t, typename id_t>
PiercedVector<value_t, id_t>::PiercedVector(const PiercedVector<value_t, id_t> &x)
    : PiercedVectorKernel<id_t>(x),
      PiercedVectorStorage<value_t, id_t>(x, this, x.getSyncMode())
{
}

/**
* Copy assignment operator.
*
* \param x is another container of the same type (i.e., instantiated with
* the same template parameters) whose content is copied in this container.
*/
template<typename value_t, typename id_t>
PiercedVector<value_t, id_t> & PiercedVector<value_t, id_t>::operator=(PiercedVector<value_t, id_t> x)
{
    this->swap(x);

    return *this;
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
    FillAction reclaimAction = PiercedVectorKernel<id_t>::fillHead(id);

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
    FillAction reclaimAction = PiercedVectorKernel<id_t>::fillAfter(referenceId, id);

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
    FillAction reclaimAction = PiercedVectorKernel<id_t>::fillAppend(id);

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
    FillAction reclaimAction = PiercedVectorKernel<id_t>::fillBefore(referenceId, id);

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
    MoveAction moveAction = PiercedVectorKernel<id_t>::moveBefore(referenceId, id, !delayed);

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
    MoveAction moveAction = PiercedVectorKernel<id_t>::moveAfter(referenceId, id, !delayed);

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
    FillAction insertAction = PiercedVectorKernel<id_t>::fillHead(id);

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
    FillAction insertAction = PiercedVectorKernel<id_t>::fillAfter(referenceId, id);

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
    FillAction insertAction = PiercedVectorKernel<id_t>::fillBefore(referenceId, id);

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
    std::size_t pos = PiercedVectorKernel<id_t>::getPos(id);

    // Replace the value
    PiercedVectorStorage<value_t, id_t>::rawSet(pos, std::move(value));

    // Return the iterator that points to the element
    return PiercedVectorStorage<value_t, id_t>::rawFind(pos);
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
    FillAction insertAction = PiercedVectorKernel<id_t>::fillAppend(id);

    // Insert the new value
    return insertValue(insertAction, value);
}

/**
* The container is extended by inserting a new element. If the element can
* reuse an existing position that position will be initialize using args
* as the argument for its initialization otherwise a new element will be
* created in-place using args as the arguments for its construction.
*
* This function is only enabled if the object stored in the container has
* an initialization function with a signature like "void initialize(Args...)".
*
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct or initialize the
* new element
* \result An iterator that points to the the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type *>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emreclaim(id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillHead(id);

    // Create the new value in-place
    return emreclaimValue<Args...>(emplaceAction, std::forward<Args>(args)...);
}

/**
* The container is extended by inserting a new element. The element will have
* a position that is between the element with the specified reference id and
* the end of the container. If the element can reuse an existing position
* that position will be initialize using args as the argument for its
* initialization otherwise a new element will be created in-place using args
* as the arguments for its construction.
*
* This function is only enabled if the object stored in the container has
* an initialization function with a signature like "void initialize(Args...)".
*
* \param referenceId is the id of the element after which the
* new element will be inserted
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type *>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emreclaimAfter(const id_t &referenceId, id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillAfter(referenceId, id);

    // Create the new value in-place
    return emreclaimValue(emplaceAction, std::forward<Args>(args)...);
}

/**
* Inserts a new element at the end of the container, right after its current
* last element. If the element can reuse an existing position that position
* will be initialize using args as the argument for its initialization
* otherwise a new element will be created in-place using args as the arguments
* for its construction.
*
* This function is only enabled if the object stored in the container has
* an initialization function with a signature like "void initialize(Args...)".
*
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type *>
void PiercedVector<value_t, id_t>::emreclaimBack(id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillAppend(id);

    // Create the new value in-place
    emreclaimValue(emplaceAction, std::forward<Args>(args)...);
}

/**
* The container is extended by inserting a new element. The element will have
* a position that is between the beginning of the container and the element
* with the specified reference id. If the element can reuse an existing
* position that position will be initialize using args as the argument for its
* initialization otherwise a new element will be created in-place using args
* as the arguments for its construction.
*
* This function is only enabled if the object stored in the container has
* an initialization function with a signature like "void initialize(Args...)".
*
* \param referenceId is the id of the element before which the
* new element will be inserted
* \param id is the id that will be associated to the element
* \param args are the arguments forwarded to construct the new element
* \result An iterator that points to the newly inserted element.
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type *>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emreclaimBefore(const id_t &referenceId, id_t id, Args&&... args)
{
    // Fill a position
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillBefore(referenceId, id);

    // Create the new value in-place
    return emreclaimValue(emplaceAction, std::forward<Args>(args)...);
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
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillHead(id);

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
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillAfter(referenceId, id);

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
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillAppend(id);

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
    FillAction emplaceAction = PiercedVectorKernel<id_t>::fillBefore(referenceId, id);

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
    std::size_t pos = PiercedVectorKernel<id_t>::getPos(id);

    // Replace the value
    PiercedVectorStorage<value_t, id_t>::rawEmreplace(pos, std::forward<Args>(args)...);

    // Return the iterator that points to the element
    return PiercedVectorStorage<value_t, id_t>::rawFind(pos);
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
    EraseAction eraseAction = PiercedVectorKernel<id_t>::erase(id, !delayed);

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
    EraseAction eraseAction = PiercedVectorKernel<id_t>::popBack();

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
    SwapAction swapAction = PiercedVectorKernel<id_t>::swap(id_first, id_second);

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
    ClearAction clearAction = PiercedVectorKernel<id_t>::clear(release);

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(clearAction);
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
    ReserveAction reserveAction = PiercedVectorKernel<id_t>::reserve(n);

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(reserveAction);
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
    ResizeAction resizeAction = PiercedVectorKernel<id_t>::resize(n);

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(resizeAction);
}

/**
* Sorts the elements of the container in ascending id order.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::sort()
{
    // Update the kernel
    SortAction sortAction = PiercedVectorKernel<id_t>::sort();

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(sortAction);
}

/**
* Sorts the container after the element with the reference id in ascending
* id order.
*
* \param referenceId is the id of the element after which the container will
* be sorted
* \param inclusive if true the reference element will be sorted, otherwise
* the sorting will stop at the element following the reference
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::sortAfter(id_t referenceId, bool inclusive)
{
    // Update the kernel
    SortAction sortAction = PiercedVectorKernel<id_t>::sortAfter(referenceId, inclusive);

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(sortAction);
}

/**
* Sorts the container before the element with the reference id in ascending
* id order.
*
* \param referenceId is the id of the element before which the container will
* be sorted
* \param inclusive if true the reference element will be sorted, otherwise
* the sorting will stop at the element preceding the reference
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::sortBefore(id_t referenceId, bool inclusive)
{
    // Update the kernel
    SortAction sortAction = PiercedVectorKernel<id_t>::sortBefore(referenceId, inclusive);

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(sortAction);
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
    SqueezeAction squeezeAction = PiercedVectorKernel<id_t>::squeeze();

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(squeezeAction);
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
    ShrinkToFitAction shrinkToFitAction = PiercedVectorKernel<id_t>::shrinkToFit();

    // Update the storage
    PiercedVectorStorage<value_t, id_t>::commitSyncAction(shrinkToFitAction);
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
    // The swap will swap also the slave-master information. This is not what
    // we want, therefore the two pierced storage will be unregistered and the
    // registered again after the swap. When the kernel is unset the storage
    // can't be clear, otherwise its contents will be lost.
    PiercedVectorStorage<value_t, id_t>::detachKernel();
    x.PiercedVectorStorage<value_t, id_t>::detachKernel();

    // Swap kernel data
    PiercedVectorKernel<id_t>::swap(x);

    // Swap storage data
    PiercedVectorStorage<value_t, id_t>::swap(x);

    // Re-register the storages
    //
    // The function that sets a dynamic kernel may throw an exception if the
    // kernel is already set or if we are trying to set a null kernel. Here
    // neither of the two cases can happen, because the kernel has been
    // previously cleared and the kernel we are trying to set is not null.
    try {
        PiercedVectorStorage<value_t, id_t>::setDynamicKernel(this, PiercedVectorKernel<id_t>::SYNC_MODE_DISABLED);
        x.PiercedVectorStorage<value_t, id_t>::setDynamicKernel(&x, PiercedVectorKernel<id_t>::SYNC_MODE_DISABLED);
    } catch (const std::runtime_error &exception) {
        assert(false && "Error while swapping the PiercedVector!");
    }
}

/**
* Get a constant reference to the kernel of the vector.
*
* \result A constant reference to the kernel of the vector.
*/
template<typename value_t, typename id_t>
const PiercedVectorKernel<id_t> & PiercedVector<value_t, id_t>::getKernel() const
{
    return *this;
}

/**
* Get a constant reference to the storage of the vector.
*
* \result A constant reference to the storage of the vector.
*/
template<typename value_t, typename id_t>
const PiercedVectorStorage<value_t, id_t> & PiercedVector<value_t, id_t>::getStorage() const
{
    return *this;
}

/**
* Dumps to screen the internal data.
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::dump() const
{
    PiercedVectorKernel<id_t>::dump();
}

/**
* Restore the container.
*
* \param stream is the stream data should be read from
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<PiercedVectorStorage<T, id_t>::has_restore()>::type *>
void PiercedVector<value_t, id_t>::restore(std::istream &stream)
{
    restoreKernel(stream);

    PiercedVectorStorage<T, id_t>::restore(stream);
}

/**
* Dump the container.
*
* \param stream is the stream data should be written to
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<PiercedVectorStorage<T, id_t>::has_dump()>::type *>
void PiercedVector<value_t, id_t>::dump(std::ostream &stream) const
{
    dumpKernel(stream);

    PiercedVectorStorage<T, id_t>::dump(stream);
}

/**
* Restore the kernel's container.
*
* \param stream is the stream data should be read from
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::restoreKernel(std::istream &stream)
{
    PiercedVectorKernel<id_t>::restore(stream);

    PiercedVectorStorage<value_t, id_t>::rawResize(PiercedVectorKernel<id_t>::rawSize());
}

/**
* Dump the kernel's container.
*
* \param stream is the stream data should be written to
*/
template<typename value_t, typename id_t>
void PiercedVector<value_t, id_t>::dumpKernel(std::ostream &stream) const
{
    PiercedVectorKernel<id_t>::dump(stream);
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
        PiercedVectorStorage<value_t, id_t>::rawEmplace(pos);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawEmplaceBack();
        pos = PiercedVectorKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedVectorStorage<value_t, id_t>::rawFind(pos);
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
        PiercedVectorStorage<value_t, id_t>::rawSet(pos, value);
        break;
    }

    case FillAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawInsert(pos, 1, value);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawPushBack(value);
        pos = PiercedVectorKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedVectorStorage<value_t, id_t>::rawFind(pos);
}

/**
* Insert the specified element constructing it in-place.
*
* This function is only enabled if the object stored in the container has
* an initialization function with a signature like "void initialize(Args...)".
*
* \param action is the fill action that defines how to insert the element
* \param args are the arguments forwarded to the elements' construct when
* synchronizing the action
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type *>
typename PiercedVector<value_t, id_t>::iterator PiercedVector<value_t, id_t>::emreclaimValue(const FillAction &action, Args&&... args)
{
    std::size_t pos = action.info[PiercedSyncAction::INFO_POS];
    switch (static_cast<typename FillAction::FillActionType>(action.type)) {

    case FillAction::TYPE_OVERWRITE:
    {
        PiercedVectorStorage<value_t, id_t>::rawInitialize(pos, std::forward<Args>(args)...);
        break;
    }

    case FillAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawEmplace(pos, std::forward<Args>(args)...);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawEmplaceBack(std::forward<Args>(args)...);
        pos = PiercedVectorKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedVectorStorage<value_t, id_t>::rawFind(pos);
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
        PiercedVectorStorage<value_t, id_t>::rawEmreplace(pos, std::forward<Args>(args)...);
        break;
    }

    case FillAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawEmplace(pos, std::forward<Args>(args)...);
        break;
    }

    case FillAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawEmplaceBack(std::forward<Args>(args)...);
        pos = PiercedVectorKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Return the iterator to the position where the element was inserted
    return PiercedVectorStorage<value_t, id_t>::rawFind(pos);
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
        PiercedVectorStorage<value_t, id_t>::rawSet(posNew, std::move(PiercedVectorStorage<value_t, id_t>::rawAt(posOld)));
        break;
    }

    case MoveAction::TYPE_MOVE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawInsert(posNew, 1, std::move(PiercedVectorStorage<value_t, id_t>::rawAt(posOld)));
        break;
    }

    case MoveAction::TYPE_MOVE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        PiercedVectorStorage<value_t, id_t>::rawPushBack(std::move(PiercedVectorStorage<value_t, id_t>::rawAt(posOld)));
        posNew = PiercedVectorKernel<id_t>::getLastUsedPos();
        break;
    }

    default:
    {
        BITPIT_UNREACHABLE("This action is not handled");
        break;
    }

    }

    // Clear the old position
    PiercedVectorStorage<value_t, id_t>::rawEmreplace(posOld);

    // Return the iterator to the new position
    return PiercedVectorStorage<value_t, id_t>::rawFind(posNew);
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
        PiercedVectorStorage<value_t, id_t>::rawResize(rawSize);
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

    return PiercedVectorStorage<value_t, id_t>::rawFind(nextPos);
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
        PiercedVectorStorage<value_t, id_t>::rawSwap(posFirst, posSecond);
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
