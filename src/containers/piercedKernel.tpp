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

#ifndef __BITPIT_PIERCED_KERNEL_TPP__
#define __BITPIT_PIERCED_KERNEL_TPP__

namespace bitpit {

// Definition of static constants of PiercedKernel
template<typename id_t>
const std::size_t PiercedKernel<id_t>::MAX_PENDING_HOLES = 16384;

/**
* Constructs an empty pierced kernel with no elements.
*/
template<typename id_t>
PiercedKernel<id_t>::PiercedKernel()
    : BasePiercedKernel()
{
    clear();
}

/**
* Constructs a pierced kernel with a capacity at least enough
* to contain n elements.
*
* \param n the minimum capacity requested for the kernel
*/
template<typename id_t>
PiercedKernel<id_t>::PiercedKernel(std::size_t n)
    : BasePiercedKernel()
{
    clear();

    reserve(n);
}

/**
* Updates the id of the specified element.
*
* \param currentId is the current id of the element
* \param updatedId is the new id of the element
*/
template<typename id_t>
void PiercedKernel<id_t>::updateId(const id_t &currentId, const id_t &updatedId)
{
    // Validate the id
    validateId(updatedId);

    // Update the id
    setPosId(getPos(currentId), updatedId);
    m_pos.erase(currentId);
}

/**
* Removes all elements from the kernel.
*
* \param release if it's true the memory hold by the kernel will be released,
* otherwise the kernel will be cleared but its memory will not be relased
*/
template<typename id_t>
void PiercedKernel<id_t>::clear(bool release)
{
    // Clear positions
    m_ids.clear();
    m_pos.clear();
    if (release) {
        std::vector<id_t>().swap(m_ids);
        std::unordered_map<id_t, std::size_t, PiercedHasher>().swap(m_pos);
    }

    // Reset begin and end
    setBeginPos(0);
    setEndPos(0);

    // Clear holes
    holesClear(release);

    // There are no dirty positions
    m_dirty_begin_pos = m_end_pos;

    // Update the storage
    PiercedSyncAction syncAction(PiercedSyncAction::TYPE_CLEAR);
    processSyncAction(syncAction);
}

/**
* Flush all pending changes.
*/
template<typename id_t>
void PiercedKernel<id_t>::flush()
{
    // Flush pending holes
    holesFlush();
}

/**
* Requests that the kernel capacity be at least enough to contain n elements.
*
* If n is greater than the current kernel capacity, the function causes the
* kernel to reallocate its data structured increasing its capacity to n
* (or greater).
*
* In all other cases, the function call does not cause a reallocation and
* the kernel capacity is not affected.
*
* \param n the minimum capacity requested for the kernel
*/
template<typename id_t>
void PiercedKernel<id_t>::reserve(std::size_t n)
{
    m_ids.reserve(n);
    m_pos.reserve(n);
}

/**
* Resizes the kernel so that it contains n elements.
*
* If n is smaller than the current kernel size, the content is reduced to its
* first n elements, removing those beyond.
*
* If n is greater than the current kernel size, space is reserved to allow
* the kernel to reach the requested size.
*
* If n is also greater than the current kernel capacity, an automatic
* reallocation of the allocated kernel space takes place.
*
* Notice that this function changes the actual content of the kernel by
* erasing elements from it.
*
* \param n is the new kernel size, expressed in number of elements.
*/
template<typename id_t>
void PiercedKernel<id_t>::resize(std::size_t n)
{
    // If the size of the kernel is already the requested size there is
    // nothing to do.
    if (n == size()) {
        return;
    }

    // A request for a size equal to 0 is equivalent to a clear.
    if (n == 0) {
        clear();
        return;
    }

    // If the requested size is greater that the current size we may need to
    // reserve space to allow the kernel to reach the requested size.
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
    std::size_t last_used_pos = getPos(last_stored_id);

    // Shrink the kernel
    shrink(last_used_pos + 1);

    // Update the storage
    PiercedSyncAction syncAction(PiercedSyncAction::TYPE_RESIZE);
    syncAction.info[PiercedSyncAction::INFO_SIZE] = rawSize();
}

/**
* Sorts the elements of the kernel in ascending id order.
*/
template<typename id_t>
std::vector<std::size_t> PiercedKernel<id_t>::sort()
{
    // Squeeze the kernel
    std::vector<std::size_t> squeezePermutations = squeeze();

    // The kernel has been squeezed, there are no holes
    std::size_t nElements = size();

    // Evaluates the sort permutations
    std::vector<std::size_t> sortPermutations(nElements);
    for (std::size_t i = 0; i < nElements; ++i) {
        sortPermutations[i] = i;
    }
    std::sort(sortPermutations.begin(), sortPermutations.end(), idLess(m_ids));

    // Combine squeeze and sort permutations
    std::vector<std::size_t> permutations(nElements);
    for (std::size_t i = 0; i < nElements; ++i) {
        permutations[i] = squeezePermutations[sortPermutations[i]];
    }

    // Sort the ids
    //
    // NOTE: the reored function will destroy the permutation on output
    utils::reorderVector<id_t>(sortPermutations, m_ids, nElements);

    // Update storage
    PiercedSyncAction permutationAction(PiercedSyncAction::TYPE_REORDER);
    permutationAction.importData(sortPermutations);
    processSyncAction(permutationAction);

    // Return the permutations
    return permutations;
}

/**
* Requests the kernel to compact the elements and reduce its capacity to
* fit its size.
*
* The request is non-binding, and the function can leave the kernel with a
* capacity greater than its size.
*
* This may cause a reallocation, but has no effect on the kernel size and
* cannot alter its elements.
*/
template<typename id_t>
std::vector<std::size_t> PiercedKernel<id_t>::squeeze()
{
    std::vector<std::size_t> permutations;
    permutations.resize(rawSize());

    // Flush changes
    flush();

    // Compact the kernel
    std::size_t nHoles = holesCount();
    if (nHoles != 0) {
        // Move the elements
        std::size_t firstPosToUpdate;
        if (m_begin_pos == 0) {
            firstPosToUpdate = *(m_holes_regular_end - 1);
        } else {
            firstPosToUpdate = 0;
        }

        for (std::size_t pos = 0; pos < firstPosToUpdate; ++pos) {
            permutations[pos] = pos;
        }

        std::size_t offset = 0;
        for (std::size_t pos = firstPosToUpdate; pos < m_end_pos; pos++) {
            if (offset < nHoles && *(m_holes_regular_end - offset - 1) == pos) {
                permutations[m_end_pos - 1 - offset] = pos;
                ++offset;
                continue;
            }

            id_t id = m_ids[pos];
            std::size_t updatedPos = pos - offset;

            setPosId(updatedPos, id);
            setPosEmptyId(pos, pos + 1);
            permutations[updatedPos] = pos;
        }

        // Clear the holes
        holesClear();

        // Reset begin and end
        setBeginPos(0);
        setEndPos(size());

        // Update storage
        PiercedSyncAction permutationAction(PiercedSyncAction::TYPE_REORDER);
        permutationAction.importData(permutations);
        processSyncAction(permutationAction);

        // Shrink the kernel
        shrink(size(), true);
    } else {
        for (std::size_t pos = 0; pos < m_end_pos; ++pos) {
            permutations[pos] = pos;
        }
    }

    // Shrink to fit
    shrinkToFit();

    return permutations;
}

/**
* Requests the kernel to reduce its capacity to fit its size. This method
* will NOT compact the elements, leaving the existing holes unaltered.
*
* The request is non-binding, and the function can leave the kernel with a
* capacity greater than its size.
*
* This may cause a reallocation, but has no effect on the kernel size and
* cannot alter its elements not the holes.
*/
template<typename id_t>
void PiercedKernel<id_t>::shrinkToFit()
{
    m_ids.shrink_to_fit();

    // Update the storage
    PiercedSyncAction syncAction(PiercedSyncAction::TYPE_SHRINK_TO_FIT);
    processSyncAction(syncAction);
}

/**
* Exchanges the content of the kernel by the content of x, which is another
* kernel object of the same type. Sizes may differ.
*
* After the call to this member function, the elements in this kernel are
* those which were in x before the call, and the elements of x are those
* which were in this. All iterators, references and pointers remain valid
* for the swapped objects.
*
* \param x Another kernel of the same type (i.e., instantiated with the same
* template parameters) whose content is swapped with that of this kernel.
*/
template<typename id_t>
void PiercedKernel<id_t>::swap(PiercedKernel &x) noexcept
{
    PiercedSyncMaster::swap(x);

    std::swap(x.m_begin_pos, m_begin_pos);
    std::swap(x.m_end_pos, m_end_pos);
    std::swap(x.m_dirty_begin_pos, m_dirty_begin_pos);
    std::swap(x.m_ids, m_ids);
    std::swap(x.m_pos, m_pos);
    std::swap(x.m_holes, m_holes);
    std::swap(x.m_holes_regular_begin, m_holes_regular_begin);
    std::swap(x.m_holes_regular_end, m_holes_regular_end);
    std::swap(x.m_holes_regular_sorted, m_holes_regular_sorted);
    std::swap(x.m_holes_pending_begin, m_holes_pending_begin);
    std::swap(x.m_holes_pending_end, m_holes_pending_end);
    std::swap(x.m_holes_pending_sorted, m_holes_pending_sorted);
}

/**
* Dumps to screen the internal data.
*/
template<typename id_t>
void PiercedKernel<id_t>::dump() const
{
    std::cout << "----------------[ DUMP ]----------------" << std::endl;

    std::cout << std::endl;
    std::cout << " size: " << size() << std::endl;

    std::cout << " m_holes_regular_begin: " << std::distance<holes_const_iterator>(m_holes.begin(), m_holes_regular_begin) << std::endl;
    std::cout << " m_holes_regular_end  : " << std::distance<holes_const_iterator>(m_holes.begin(), m_holes_regular_end) << std::endl;

    std::cout << std::endl;
    std::cout << " Regular holes: " << std::endl;
    for (auto k = m_holes_regular_begin; k < m_holes_regular_end; ++k) {
        std::cout << *k << std::endl;
    }

    std::cout << std::endl;
    std::cout << " m_holes_pending_begin: " << std::distance<holes_const_iterator>(m_holes.begin(), m_holes_pending_begin) << std::endl;
    std::cout << " m_holes_pending_end  : " << std::distance<holes_const_iterator>(m_holes.begin(), m_holes_pending_end) << std::endl;

    std::cout << std::endl;
    std::cout << " Pending holes" << std::endl;
    for (auto k = m_holes_pending_begin; k < m_holes_pending_end; ++k) {
        std::cout << *k << std::endl;
    }

    std::cout << std::endl;
    std::cout << " m_begin_pos: " << m_begin_pos << std::endl;
    std::cout << " m_end_pos: " <<  m_end_pos << std::endl;
    std::cout << " Stored ids: " << std::endl;
    if (rawSize() > 0) {
        for (std::size_t k = 0; k < m_end_pos; ++k) {
            std::cout << m_ids[k] << std::endl;
        }
    } else {
        std::cout << "None" << std::endl;
    }

    std::cout << std::endl;
    std::cout << " Poistion map: " << std::endl;
    if (size() > 0) {
        for (auto itr = m_pos.cbegin(); itr != m_pos.cend(); ++itr) {
            std::cout << itr->first << " -> " << itr->second << std::endl;
        }
    } else {
        std::cout << "None" << std::endl;
    }

    std::cout << "----------------------------------------" << std::endl;
}

/**
* Check the integrity of the kernel.
*
* \param stream is the stream data should be written to
*/
template<typename id_t>
void PiercedKernel<id_t>::checkIntegrity() const
{
    for (const auto &entry : m_pos) {
        id_t id = entry.first;
        std::size_t pos = entry.second;
        if (m_ids[pos] != id) {
            std::cout << " Position " << pos << " should contain the element with id " << id << std::endl;
            std::cout << " but it contains the element with id " << m_ids[pos] << std::endl;
            throw std::runtime_error("Integrity check error");
        }
    }

    for (std::size_t pos = m_begin_pos; pos < m_end_pos; ++pos) {
        id_t id = m_ids[pos];
        if (id >= 0) {
            if (m_pos.count(id) == 0) {
                std::cout << " Position " << pos << " contains the element with id " << id << std::endl;
                std::cout << " but the position map doesn't contain that element" << std::endl;
                throw std::runtime_error("Integrity check error");
            } else if (m_pos.at(id) != pos) {
                std::cout << " Position " << pos << " contains the element with id " << id << std::endl;
                std::cout << " but the position map claims that the element is at position " << m_pos.at(id) << std::endl;
                throw std::runtime_error("Integrity check error");
            }
        }
    }
}

/**
* Returns whether the kernel is contiguous (i.e. whether it contains
* no holes).
*
* \result true if the kernel is contiguous, false otherwise.
*/
template<typename id_t>
bool PiercedKernel<id_t>::contiguous() const
{
    return (holesCount() == 0);
}

/**
* Returns whether the kernel is empty (i.e. whether its size is 0).
*
* \result true if the kernel size is 0, false otherwise.
*/
template<typename id_t>
bool PiercedKernel<id_t>::empty() const
{
    return m_pos.empty();
}

/**
* Checks if the kernel is in a state that can slow down the iterator.
*
* If there are dirty positions in the kernel, the ids associated to
* those positions will not point directly to a non-empty element.
* This means that the iterator will require a loop to reach the next
* non-empty position and this can slow down the loop. Calling the
* 'flush' function allow to recover the best performances.
*
* \result Return true if the kernel is in a state that can slow down
* the iterator, false otherwise.
*/
template<typename id_t>
bool PiercedKernel<id_t>::isIteratorSlow()
{
    return (m_dirty_begin_pos < m_end_pos);
}

/**
* Returns the maximum number of elements that the kernel can hold.
*
* This is the maximum potential size the kernel can reach due to known system
* or library implementation limitations, but the kernel is by no means
* guaranteed to be able to reach that size: it can still fail to allocate
* space at any point before that size is reached.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::maxSize() const
{
    return m_ids.max_size();
}

/**
* Returns the number of elements in the kernel.
*
* This is the number of actual objects held in the kernel, which is not
* necessarily equal to its capacity.
*
* \result The number of elements in the kernel.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::size() const
{
    return m_pos.size();
}

/**
* Returns the number of raw positions in the kernel.
*
* This is the number of raw positions in the kernel, which is not necessarily
* equal the number of actual objects held in the kernel nor to its capacity.
*
* \result The number of raw positions in the kernel.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::rawSize() const
{
    return m_ids.size();
}

/**
* Returns the size of the space currently allocated for the kernel, expressed
* in terms of elements.
*
* \result The size of the currently allocated capacity in the kernel, measured
* in terms of the number elements it can hold.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::capacity() const
{
    return m_ids.capacity();
}

/**
* Checks if the kernel contains the specified id.
*
* \param id is the id to look for
* \result Returns true if the kernel contains the specified id, otherwise it
* returns false.
*/
template<typename id_t>
bool PiercedKernel<id_t>::contains(id_t id) const
{
    return (m_pos.count(id) != 0);
}

/*!
* Gets the raw index associated to the element with the specified id.
*
* \return The raw index associated to the element with the specified id.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::getRawIndex(id_t id) const
{
    return getPos(id);
}

/**
* Gets the flat index of the element with the specified id.
*
* A flat id is the id associated to a numbering scheme that starts
* from the element in the first position of the kernel and is
* incremented by one for each element in the kernel. The first
* element will have a flat id equal to 0, the last element will
* have a flat id equal to (nElements - 1).
*
* If there is no element with the specified id, an exception is
* thrown.
*
* \param id is the id of the element for witch the flat id is requested
* \result The flat index of the element with the specified id.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::evalFlatIndex(id_t id)
{
    std::size_t pos  = getPos(id);

    // Initialize flat id with the position of the element
    std::size_t flat = pos;

    // Subtract pending holes before position
    if (holesCountPending() > 0) {
        holesSortPending();
        auto holes_itr = std::upper_bound(m_holes_pending_begin, m_holes_pending_end, pos, std::greater<std::size_t>());
        std::size_t nHolesBefore = std::distance(holes_itr, m_holes_pending_end);

        flat -= nHolesBefore;
    }

    // Subtract regular holes before position
    if (holesCountRegular() > 0) {
        holesSortRegular();
        auto holes_itr = std::upper_bound(m_holes_regular_begin, m_holes_regular_end, pos, std::greater<std::size_t>());
        std::size_t nHolesBefore = std::distance(holes_itr, m_holes_regular_end);

        flat -= nHolesBefore;
    }

    // Done
    return flat;
}

/**
* Gets a vector containing the ids of the elements stored in the kernel.
*
* \param ordered if is true the ids will be sorted in ascending order,
* otherwise the ids will be in random order
* \result A vector with the id of the elements in stored in the kernel.
*/
template<typename id_t>
std::vector<id_t> PiercedKernel<id_t>::getIds(bool ordered) const
{
    std::size_t nIds = size();

    // Initialize the vector
    std::vector<id_t> ids;
    if (nIds == 0) {
        return ids;
    }

    // Resize the vector
    ids.resize(nIds);

    // Extract the ids
    std::size_t n   = 0;
    std::size_t pos = m_begin_pos;
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

/**
* Returns the id of the elmement before which there is the requested number
* of other elements. If this element does not exist the fallback value will
* be returned.
*
* \param targetSize is the number of elements that needs to be contained
* before the marker
* \param fallback is the fallback value to be returned if the marker cannot
* be found
* \return The id of the elmement before which there is the requested number
* of other elements. If this element does not exist the fallback value will
* be returned.
*/
template<typename id_t>
id_t PiercedKernel<id_t>::getSizeMarker(std::size_t targetSize, const id_t &fallback)
{
    // If the size is zero, we return the first element, if the target
    // size is equal to the size minus one we return the last element,
    // if the target size is greater or equal the current kernel size
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
    holes_iterator regular_holes_itr = m_holes_regular_end;
    holes_iterator pending_holes_itr = m_holes_pending_end;

    std::size_t nEmpties  = 0;
    std::size_t markerPos = targetSize;
    while (true) {
        if (isPosEmpty(markerPos)) {
            markerPos = findNextUsedPos(markerPos - 1);
        }

        // Count the number of holes and pending deletes before the
        // current marker position
        if (regular_holes_itr != m_holes_regular_begin) {
            holes_iterator itr_previous = regular_holes_itr;
            regular_holes_itr = std::upper_bound(m_holes_regular_begin, regular_holes_itr, markerPos, std::greater<std::size_t>());
            nEmpties += std::distance(regular_holes_itr, itr_previous);
        }

        if (pending_holes_itr != m_holes_pending_begin) {
            holes_iterator itr_previous = pending_holes_itr;
            pending_holes_itr = std::upper_bound(m_holes_pending_begin, pending_holes_itr, markerPos, std::greater<std::size_t>());
            nEmpties += std::distance(pending_holes_itr, itr_previous);
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

/**
* Gets the position of the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element for witch the raw id is requested
* \result The position of the element with the specified id.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::getPos(id_t id) const
{
    return m_pos.at(id);
}

/**
* Gets the position of the first element in the container.
*
* If there is no element with the specified id, an exception is thrown.
*
* \result The position of the first element in the container.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::getFirstUsedPos() const
{
    if (empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return m_begin_pos;
}

/**
* Gets the position of the first element in the container.
*
* If there is no element with the specified id, an exception is thrown.
*
* \result The position of the first element in the container.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::getLastUsedPos() const
{
    if (empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return m_end_pos - 1;
}


/**
* Fills a position and assigns to it the specified id.
*
* The position filled will be the first position available.
*
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillHead(id_t id)
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

/**
* Fills a position and assigns to it the specified id.
*
* The position filled will be the first position available.
*
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillTail(id_t id)
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

/**
* Fills a position and assigns to it the specified id.
*
* The position filled will be after the specified reference position.
*
* \param referenceId is the id of the element before which the new available
* position will be searched for
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillAfter(id_t referenceId, id_t id)
{
    // Get the reference position
    std::size_t referencePos = getPos(referenceId);

    // Check if there is a pending hole that can be filled
    if (m_holes_pending_begin != m_holes_pending_end) {
        // Sort pending holes
        holesSortPending();

        // The first hole is the one with the highest position, if the
        // position of this hole is greater than the reference position
        // it is possible to use it.
        holes_iterator holeItr = m_holes_pending_begin;
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
        holes_iterator holeItr = m_holes_regular_begin;
        if (*holeItr > referencePos) {
            return fillHole(holeItr, id);
        }
    }

    // There are no holes that can be filled
    return fillAppend(id);
}

/**
* Fills a position and assigns to it the specified id.
*
* The position filled will be before the specified reference position.
*
* \param referenceId is the id of the element before which the new available
* position will be searched for
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillBefore(id_t referenceId, id_t id)
{
    // Get the reference position
    std::size_t referencePos = getPos(referenceId);

    // Check if there is a pending hole that can be filled
    if (m_holes_pending_begin != m_holes_pending_end) {
        // Sort pending holes
        holesSortPending();

        // The last hole is the one with the lowest position, if the
        // position of this hole is greater than the reference position
        // it is possible to use it.
        holes_iterator holeItr = m_holes_pending_end - 1;
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
        holes_iterator holeItr = m_holes_regular_end - 1;
        if (*holeItr < referencePos) {
            return fillHole(holeItr, id);
        }
    }

    // There are no holes that can be filled
    return fillInsert(referencePos, id);
}

/**
* Fills a position and assigns to it the specified id.
*
* A new position is created at the end of the kernel.
*
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillAppend(id_t id)
{
    // Check if the id is valid
    validateId(id);

    // Add the id
    m_ids.emplace_back();
    id_t &storedId = m_ids.back();
    storedId = id;

    // Update last used position
    setEndPos(rawSize());

    // Update the id map
    m_pos[id] = m_end_pos - 1;

    // Update the storage
    FillAction fillAction(FillAction::TYPE_APPEND);
    fillAction.info[PiercedSyncAction::INFO_POS] = m_end_pos - 1;
    processSyncAction(fillAction);

    // Return the fill action
    return fillAction;
}

/**
* Fills the specified hole with the given id.
*
* \param pos is the position to fill
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillHole(const holes_iterator &holeItr, id_t id)
{
    // Validate the id
    validateId(id);

    // Get the position of the hole
    std::size_t pos = *holeItr;

    // Fill the hole
    setPosId(pos, id);

    // Remove the holes from the kernel
    if (holeItr >= m_holes_pending_begin) {
        int nPendings = holesCountPending();
        if (nPendings > 1) {
            if (holeItr == (m_holes_pending_end - 1)) {
                --m_holes_pending_end;
            } else if (holeItr == m_holes_pending_begin) {
                ++m_holes_pending_begin;
            } else {
                throw std::out_of_range("Only holes at the beginning or at the end of the kernel can be filled");
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
                throw std::out_of_range("Only holes at the beginning or at the end of the kernel can be filled");
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
    // to iterate through the kernel.
    std::size_t previousPos = pos;
    while (previousPos > 0) {
        --previousPos;
        if (!isPosEmpty(previousPos)) {
            break;
        }

        id_t previousId = m_ids[previousPos];
        if (previousId == -1 || (std::size_t) std::abs(previousId) == (pos - previousPos)) {
            break;
        }

        setPosEmptyId(previousPos, pos);
    }

    // Update the storage
    FillAction fillAction(FillAction::TYPE_OVERWRITE);
    fillAction.info[PiercedSyncAction::INFO_POS] = pos;
    processSyncAction(fillAction);

    // Return the fill action
    return fillAction;
}

/**
* Fills a position and assigns to it the specified id.
*
* A new postion is created at the specified position.
*
* \param pos is the position at which the id will be inserted
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillInsert(std::size_t pos, id_t id)
{
    // If the position at which we want to insert an element is a hole there
    // was an error somewhere. Before inserting a new position the hole needs
    // to be filled first.
    if (isPosEmpty(pos)) {
        throw std::out_of_range("Before inserting a new position the hole needs to be filled first");
    }

    // We cannot insert elements past the last position
    if (pos >= m_end_pos) {
        throw std::out_of_range("Unable to insert elements past the last position");
    }

    // Check if the id is valid
    validateId(id);

    // Add the id
    m_ids.emplace(m_ids.begin() + pos, id);

    // Update last used position
    setEndPos(rawSize());

    // Update the id map
    for (std::size_t i = pos + 1; i < m_end_pos; ++i) {
        id_t id_i = m_ids[i];
        if (id_i >= 0) {
            m_pos[id_i] = i;
        }
    }
    m_pos[id] = pos;

    // Update the regular holes
    if (m_holes_regular_begin != m_holes_regular_end) {
        holes_iterator change_begin = m_holes_regular_begin;
        holes_iterator change_end   = upper_bound(m_holes_regular_begin, m_holes_regular_end, pos, std::greater<std::size_t>());
        for (auto itr = change_begin; itr != change_end; itr++) {
            (*itr)++;
        }
    }

    // Update the pending holes
    if (m_holes_pending_begin != m_holes_pending_end) {
        holes_iterator change_begin = m_holes_pending_begin;
        holes_iterator change_end   = upper_bound(m_holes_pending_begin, m_holes_pending_end, pos, std::greater<std::size_t>());
        for (auto itr = change_begin; itr != change_end; itr++) {
            (*itr)++;
        }
    }

    // Update the storage
    FillAction fillAction(FillAction::TYPE_INSERT);
    fillAction.info[PiercedSyncAction::INFO_POS] = pos;
    processSyncAction(fillAction);

    // Return the fill action
    return fillAction;
}

/**
* Move the specified id after the element with the reference id.
*
* \param referenceId is the id of the element after which the specified id
* will be moved
* \param id is the id that will be moved
* \param flush controls if the holes will be flush after updating the kernel
* \result The synchronization action associated with the move.
*/
template<typename id_t>
typename PiercedKernel<id_t>::MoveAction PiercedKernel<id_t>::moveAfter(id_t referenceId, id_t id, bool flush)
{
    assert(referenceId != id);

    // Disables the synchronization
    bool syncEnabled = isSyncEnabled();
    if (syncEnabled) {
        setSyncEnabled(false);
    }

    // Pierce the position
    std::size_t initialPos = getPos(id);
    pierce(initialPos, flush);

    // Insert the element in the updated position
    FillAction fillAction = fillAfter(referenceId, id);

    // Re-enables the synchronization
    if (syncEnabled) {
        setSyncEnabled(true);
    }

    // Update the storage
    typename MoveAction::MoveActionType moveActionType;
    switch (static_cast<typename FillAction::FillActionType>(fillAction.type)) {

    case FillAction::TYPE_APPEND:
        moveActionType = MoveAction::TYPE_APPEND;
        break;

    case FillAction::TYPE_INSERT:
        moveActionType = MoveAction::TYPE_INSERT;
        break;

    case FillAction::TYPE_OVERWRITE:
        moveActionType = MoveAction::TYPE_OVERWRITE;
        break;

    default:
        BITPIT_UNREACHABLE("This action is not handled");
        moveActionType = MoveAction::TYPE_UNDEFINED;
        break;

    }

    MoveAction moveAction(moveActionType);
    moveAction.info[PiercedSyncAction::INFO_POS_FIRST]  = initialPos;
    moveAction.info[PiercedSyncAction::INFO_POS_SECOND] = fillAction.info[PiercedSyncAction::INFO_POS];
    processSyncAction(moveAction);

    // Return the move action
    return moveAction;
}

/**
* Move the specified id before the element with the reference id.
*
* \param referenceId is the id of the element after which the specified id
* will be moved
* \param id is the id that will be moved
* \param flush controls if the holes will be flush after updating the kernel
* \result The synchronization action associated with the move.
*/
template<typename id_t>
typename PiercedKernel<id_t>::MoveAction PiercedKernel<id_t>::moveBefore(id_t referenceId, id_t id, bool flush)
{
    assert(referenceId != id);

    // Disables the synchronization
    bool syncEnabled = isSyncEnabled();
    if (syncEnabled) {
        setSyncEnabled(false);
    }

    // Pierce the position
    std::size_t initialPos = getPos(id);
    pierce(initialPos, flush);

    // Insert the element in the updated position
    FillAction fillAction = fillBefore(referenceId, id);

    // Re-enables the synchronization
    if (syncEnabled) {
        setSyncEnabled(true);
    }

    // Update the storage
    typename MoveAction::MoveActionType moveActionType;
    switch (static_cast<typename FillAction::FillActionType>(fillAction.type)) {

    case FillAction::TYPE_APPEND:
        moveActionType = MoveAction::TYPE_APPEND;
        break;

    case FillAction::TYPE_INSERT:
        moveActionType = MoveAction::TYPE_INSERT;
        break;

    case FillAction::TYPE_OVERWRITE:
        moveActionType = MoveAction::TYPE_OVERWRITE;
        break;

    default:
        BITPIT_UNREACHABLE("This action is not handled");
        moveActionType = MoveAction::TYPE_UNDEFINED;
        break;

    }

    MoveAction moveAction(moveActionType);
    moveAction.info[PiercedSyncAction::INFO_POS_FIRST]  = initialPos;
    moveAction.info[PiercedSyncAction::INFO_POS_SECOND] = fillAction.info[PiercedSyncAction::INFO_POS];
    processSyncAction(moveAction);

    // Return the move action
    return moveAction;
}


/**
* Swap the elements with the specified id.
*
* \param id_first is the id of the first element to be swapped
* \param id_second is the id of the second element to be swapped
*/
template<typename id_t>
typename PiercedKernel<id_t>::SwapAction PiercedKernel<id_t>::swap(id_t id_first, id_t id_second)
{
    // Positions
    std::size_t pos_first  = m_pos.at(id_first);
    std::size_t pos_second = m_pos.at(id_second);

    // Swap the positions
    swapPosIds(pos_first, id_first, pos_second, id_second);

    // Update the storage
    SwapAction swapAction(SwapAction::TYPE_SWAP);
    swapAction.info[PiercedSyncAction::INFO_POS_FIRST]  = pos_first;
    swapAction.info[PiercedSyncAction::INFO_POS_SECOND] = pos_second;
    processSyncAction(swapAction);

    // Return the swap action
    return swapAction;
}

/**
* Removes from the kernel the element with the specified id. If the id does
* not exists the function throws an exception.
*
* The position of the element is marked as empty, this tells the kernel that
* the position is no longer in use.
*
* \param id is the id of the element to erase
* \param flush controls if the holes will be flush after updating the kernel
* \result The synchronization action associated with the erase.
*/
template<typename id_t>
typename PiercedKernel<id_t>::EraseAction PiercedKernel<id_t>::erase(id_t id, bool flush)
{
    // Position
    std::size_t pos = m_pos.at(id);

    // Pierce the position
    pierce(pos, flush);

    // Update the storage
    if (pos + 1 < m_end_pos) {
        EraseAction pierceAction(EraseAction::TYPE_PIERCE);
        pierceAction.info[PiercedSyncAction::INFO_POS]      = pos;
        pierceAction.info[PiercedSyncAction::INFO_POS_NEXT] = findNextUsedPos(pos);
        processSyncAction(pierceAction);

        return pierceAction;
    } else {
        EraseAction resizeAction(EraseAction::TYPE_SHRINK);
        resizeAction.info[PiercedSyncAction::INFO_SIZE] = rawSize();
        processSyncAction(resizeAction);

        return resizeAction;
    }
}

/*!
    Removes the last element in the vector, effectively reducing
    the container size by one.

    Element is not deleted from the internal vector, instead its
    id is changed to mark the position as empty and allow the
    container to reuse that position.
*/
template<typename id_t>
typename PiercedKernel<id_t>::EraseAction PiercedKernel<id_t>::popBack()
{
    // Erase the last element
    if (empty()) {
        throw std::out_of_range("Vector is empty");
    } else if (size() == 1) {
        clear();
    } else {
        std::size_t last_used_pos = findPrevUsedPos(m_end_pos - 1);
        shrink(last_used_pos + 1);
    }

    // Return the erase action
    EraseAction eraseAction(EraseAction::TYPE_SHRINK);
    eraseAction.info[PiercedSyncAction::INFO_SIZE] = rawSize();
    processSyncAction(eraseAction);

    return eraseAction;
}

/**
* Marks a position as empty.
*
* The position is inserted in the list of holes. If the list of pending
* holes is full, a flush is called before adding the hole. This means
* that the hole is always added as a pending hole.
*
* \param pos is the position of the new hole
* \param flush controls if the holes will be flush after piercing the
* the position
*/
template<typename id_t>
void PiercedKernel<id_t>::pierce(std::size_t pos, bool flush)
{
    // If removing the last position, there is no need to add the
    // position to the holes, it's enough to update the last position
    // counter or clear the kernel if this was the last hole.
    if (pos + 1 == m_end_pos) {
        if (size() == 1) {
            clear();
        } else {
            std::size_t last_used_pos = findPrevUsedPos(m_end_pos - 1);
            shrink(last_used_pos + 1);
        }
        return;
    }

    // Remove the id from the map
    id_t id = m_ids[pos];
    m_pos.erase(id);

    // Reset the position
    std::size_t nextUsedPos = findNextUsedPos(pos);
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
        std::size_t nPendings = holesCountPending();
        if (nPendings > 1 && (*(m_holes_pending_end - 1) > *(m_holes_pending_end - 2))) {
            m_holes_pending_sorted = false;
        }
    }

    // Flush
    if (flush) {
        holesFlush();
    }
}

/**
* Clear the list of available holes.
*
* \param relase if set to true the memory previously hold by holes'
* kernel will be released
*/
template<typename id_t>
void PiercedKernel<id_t>::holesClear(bool release)
{
    // Clear the holes
    holesResize(0, 0, 0, release);

    // There are no holes, therefore all holes are sorted
    m_holes_regular_sorted = true;
    m_holes_pending_sorted = true;
}

/**
* Clear regular holes.
*
* \param relase if set to true the memory previously hold by holes'
* kernel will be released
*/
template<typename id_t>
void PiercedKernel<id_t>::holesClearRegular(bool release)
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

/**
* Clear pending holes
*
* \param relase if set to true the memory previously hold by holes'
* kernel will be released
*/
template<typename id_t>
void PiercedKernel<id_t>::holesClearPending(bool release)
{
    // Clear section of the kernel associated with the pending holes
    long offset    = std::distance(m_holes.begin(), m_holes_regular_begin);
    long nRegulars = holesCountRegular();
    holesResize(offset, nRegulars, 0, release);

    // There are no pending holes, therefore pending holes are sorted
    m_holes_pending_sorted = true;
}

/**
* Resize the kernel of the pending holes
*
* \param offset is the distance between the first regular hole and the
* begin of the hole's kernel
* \param nRegulars is the number of regulars holes
* \param nPendings  the number of pending holes
* \param relase if set to true the memory previously hold by holes'
* kernel will be released
*/
template<typename id_t>
void PiercedKernel<id_t>::holesResize(std::size_t offset, std::size_t nRegulars, std::size_t nPendings, bool release)
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

/**
* Count the available holes.
*
* \result The number of available holes.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::holesCount() const
{
    return holesCountPending() + holesCountRegular();
}

/**
* Count the pending holes.
*
* \result The number of pending holes.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::holesCountPending() const
{
    return std::distance(m_holes_pending_begin, m_holes_pending_end);
}

/**
* Count the regular holes.
*
* \result The number of regular holes.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::holesCountRegular() const
{
    return std::distance(m_holes_regular_begin, m_holes_regular_end);
}

/**
* Flushes the list of available holes.
*
* All the pending hole are converted to regular holes and new
* space is reserved for future pending holes.
*/
template<typename id_t>
void PiercedKernel<id_t>::holesFlush()
{
    // If there are no pending holes there is nothing to do
    if (m_holes_pending_begin == m_holes_pending_end) {
        return;
    }

    // Update the id of the empty elements
    //
    // The list of pending holes is sorted, in this way we can iterate
    // from the last pending hole in the kernel to the first one.
    // We start updating the id of the last pending hole, updating also
    // the id of the contiguous holes before that one. Then we advance
    // to the next hole, skipping the positions that have already been
    // updated.
    holesSortPending();

    auto itr = m_holes_pending_begin;
    std::size_t pos = m_end_pos;
    do {
        if (*itr >= pos) {
            itr++;
            continue;
        }

        pos = *itr;
        std::size_t next_used_pos = findNextUsedPos(pos);
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
        std::size_t pos = *itr;

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


/**
* Sort the list of pending holes in descendent order
*/
template<typename id_t>
void PiercedKernel<id_t>::holesSortPending()
{
    if (m_holes_pending_sorted) {
        return;
    }

    std::sort(m_holes_pending_begin, m_holes_pending_end, std::greater<std::size_t>());
    m_holes_pending_sorted = true;
}

/**
* Sort the list of regular holes in descendent order
*/
template<typename id_t>
void PiercedKernel<id_t>::holesSortRegular()
{
    if (m_holes_regular_sorted) {
        return;
    }

    std::sort(m_holes_regular_begin, m_holes_regular_end, std::greater<std::size_t>());
    m_holes_regular_sorted = true;
}

/**
* Validates the specified id.
*
* \param id is the id to validate
* \result Return true if the id is valid, otherwise it throws an exception.
*/
template<typename id_t>
bool PiercedKernel<id_t>::validateId(id_t id)
{
    // Ids needs to be positive
    if (id < 0) {
        throw std::out_of_range("Negative id");
    }

    // Handle duplicate ids
    if (contains(id)) {
        throw std::out_of_range("Duplicate id");
    }

    // Id is valid
    return true;
}

/**
* Returns the position of the last element in the kernel. If the kernel
* is empty, an exception is thrown.
*
* \result the position of the last element in the kernel.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::back() const
{
    if (empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return (m_end_pos - 1);
}

/**
* Returns the position of the first element in the kernel. If the kernel
* is empty, an exception is thrown.
*
* \result the position of the first element in the kernel.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::front() const
{
    if (empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return m_begin_pos;
}

/**
* Returns the first non-empty position before the specified starting position.
*
* If the starting position is the first posistion, an exception is thrown.
*
* \param pos starting position
* \result The firt non-empty position before the starting position.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::findPrevUsedPos(std::size_t pos) const
{
    std::size_t prev_pos = pos;
    while (true) {
        if (prev_pos == m_begin_pos) {
            throw std::out_of_range("Already in the firts position");
        }
        prev_pos--;

        id_t prev_id = m_ids[prev_pos];
        if (prev_id >= 0) {
            return prev_pos;
        }
    }
}

/**
* Returns the first non-empty position after the specified starting position.
*
* If the starting position is the last posistion, an exception is thrown.
*
* \param pos starting position
* \result The firt non-empty position after the starting position.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::findNextUsedPos(std::size_t pos) const
{
    std::size_t next_pos   = pos;
    std::size_t next_delta = 1;
    while (true) {
        if (next_pos + 1 == m_end_pos) {
            throw std::out_of_range("Already in the last position");
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

/**
* Returns if the specified position is empty.
*
* A position is considered empty if the element in that position has an id
* less than 0.
*
* \param pos the position to check
* \result true is the position is empty, false otherwise.
*/
template<typename id_t>
bool PiercedKernel<id_t>::isPosEmpty(std::size_t pos) const
{
    return (m_ids[pos] < 0);
}

/**
* Associate a position to the specified id.
*
* \param pos is the position to associate
* \param id is the id that will be associated to the position
*/
template<typename id_t>
void PiercedKernel<id_t>::setPosId(std::size_t pos, id_t id)
{
    m_ids[pos] = id;
    m_pos[id]  = pos;
}

/**
* Updates the id of the specified position to mark it as empty element.
*
* The id of an empty element contains the distance, measured in
* number of elements, between the current element and the next
* non-empty element (the distance is negative).
*
* The next used position has to be greater than the position to
* update. Otherwise, undefined behavior occurs.
*
* \param pos is the position to update
* \param nextUsedPos is the position of the next non-empty element
*/
template<typename id_t>
void PiercedKernel<id_t>::setPosEmptyId(std::size_t pos, std::size_t nextUsedPos)
{
    assert(nextUsedPos > pos);

    m_ids[pos] = pos - nextUsedPos;
}

/**
* Swaps two positions.
*
* \param pos_1 is the first position to swap
* \param id_1 is the id associated to the first position
* \param pos_2 is the second position to swap
* \param id_2 is the id associated to the second position
*/
template<typename id_t>
void PiercedKernel<id_t>::swapPosIds(std::size_t pos_1, id_t id_1, std::size_t pos_2, id_t id_2)
{
    std::swap(m_ids[pos_1], m_ids[pos_2]);
    std::swap(m_pos[id_1], m_pos[id_2]);
}

/**
* Set the first used position.
*/
template<typename id_t>
void PiercedKernel<id_t>::setBeginPos(std::size_t pos)
{
    m_begin_pos = pos;
}

/**
* Set the last used position.
*/
template<typename id_t>
void PiercedKernel<id_t>::setEndPos(std::size_t pos)
{
    m_end_pos = pos;
}

/**
* Shrink the kernel so that it contains n raw positions.
*
* \param n is the new kernel size, expressed in number of raw positions.
*/
template<typename id_t>
void PiercedKernel<id_t>::shrink(std::size_t n, bool force)
{
    // We can only shrink the kernel
    if (n > m_end_pos) {
        throw std::out_of_range("The kernel can only be shrunk");
    }

    // Check if we actually need to shrink the kernel
    if (n == m_end_pos && !force) {
        return;
    }

    // When the new last position is before the first one this is equivalent
    // to a clear
    if (n < (m_begin_pos + 1)) {
        clear();
        return;
    }

    // Delete the ids of the elements that will be removed
    for (std::size_t pos = n; pos < m_end_pos; ++pos) {
        id_t id = m_ids[pos];
        if (id >= 0) {
            m_pos.erase(id);
        }
    }

    // Resize the internal vectors
    m_ids.resize(n);

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

    // Update the storage
    PiercedSyncAction syncAction(PiercedSyncAction::TYPE_RESIZE);
    syncAction.info[PiercedSyncAction::INFO_SIZE] = n;
    processSyncAction(syncAction);
}

/**
* Register the specified storage
*
* \param storage is the storage that will be registered
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename id_t>
void PiercedKernel<id_t>::registerStorage(PiercedSyncSlave *storage, PiercedSyncMaster::SyncMode syncMode)
{
    registerSlave(storage, syncMode);
}

/**
* Unregister the specified storage->
*
* \param storage is the storage that will be unregistered
*/
template<typename id_t>
void PiercedKernel<id_t>::unregisterStorage(const PiercedSyncSlave *storage)
{
    unregisterSlave(storage);
}

/**
* Check if te specified storage is registered.
*
* \param storage is the storage to check
*/
template<typename id_t>
bool PiercedKernel<id_t>::isStorageRegistered(const PiercedSyncSlave *storage) const
{
    return isSlaveRegistered(storage);
}

/**
* Get the synchronization mode for the specified storage
*
* \param storage is the storage for which the synchronization mode is requested
* \result The synchronization mode of the storage->
*/
template<typename id_t>
typename PiercedSyncMaster::SyncMode PiercedKernel<id_t>::getStorageSyncMode(const PiercedSyncSlave *storage) const
{
    return getSlaveSyncMode(storage);
}

/**
* Restore the vector.
*
* \param stream is the stream data should be read from
*/
template<typename id_t>
void PiercedKernel<id_t>::restore(std::istream &stream)
{
    // Ids data
    std::size_t nIds;
    utils::binary::read(stream, nIds);
    m_pos.reserve(nIds);
    for (std::size_t n = 0; n < nIds; ++n) {
        std::size_t id;
        utils::binary::read(stream, id);

        std::size_t pos;
        utils::binary::read(stream, pos);

        m_pos.insert({id, pos});
    }

    // Postions data
    std::size_t nPositions;
    utils::binary::read(stream, nPositions);
    m_ids.resize(nPositions);
    for (std::size_t n = 0; n < nPositions; ++n) {
        std::size_t id;
        utils::binary::read(stream, id);

        m_ids[n] = id;
    }

    utils::binary::read(stream, m_begin_pos);
    utils::binary::read(stream, m_end_pos);
    utils::binary::read(stream, m_dirty_begin_pos);

    // Holes data
    std::size_t nHoles = size();
    utils::binary::read(stream, nHoles);
    m_holes.resize(nHoles);
    for (std::size_t n = 0; n < nHoles; ++n) {
        std::size_t pos;
        utils::binary::read(stream, pos);

        m_holes[n] = pos;
    }

    std::size_t distance;
    utils::binary::read(stream, distance);
    m_holes_regular_begin = m_holes.begin() + distance;
    utils::binary::read(stream, distance);
    m_holes_regular_end = m_holes.begin() + distance;
    utils::binary::read(stream, distance);
    m_holes_pending_begin = m_holes.begin() + distance;
    utils::binary::read(stream, distance);
    m_holes_pending_end = m_holes.begin() + distance;
    utils::binary::read(stream, m_holes_regular_sorted);
    utils::binary::read(stream, m_holes_pending_sorted);

    // Synchronization data
    PiercedSyncMaster::restore(stream);
}

/**
* Dump the vector.
*
* \param stream is the stream data should be written to
*/
template<typename id_t>
void PiercedKernel<id_t>::dump(std::ostream &stream) const
{
    // Ids data
    std::size_t nIds = size();
    utils::binary::write(stream, nIds);
    for (const auto &entry : m_pos) {
        utils::binary::write(stream, entry.first);
        utils::binary::write(stream, entry.second);
    }

    // Postions data
    std::size_t nPositions = size();
    utils::binary::write(stream, nPositions);
    for (std::size_t pos : m_ids) {
        utils::binary::write(stream, pos);
    }

    utils::binary::write(stream, m_begin_pos);
    utils::binary::write(stream, m_end_pos);
    utils::binary::write(stream, m_dirty_begin_pos);

    // Holes data
    std::size_t nHoles = size();
    utils::binary::write(stream, nHoles);
    for (std::size_t hole : m_holes) {
        utils::binary::write(stream, hole);
    }

    utils::binary::write(stream, std::distance<holes_const_iterator>(m_holes.begin(), m_holes_regular_begin));
    utils::binary::write(stream, std::distance<holes_const_iterator>(m_holes.begin(), m_holes_regular_end));
    utils::binary::write(stream, std::distance<holes_const_iterator>(m_holes.begin(), m_holes_pending_begin));
    utils::binary::write(stream, std::distance<holes_const_iterator>(m_holes.begin(), m_holes_pending_end));
    utils::binary::write(stream, m_holes_regular_sorted);
    utils::binary::write(stream, m_holes_pending_sorted);

    // Synchronization data
    PiercedSyncMaster::dump(stream);
}

}

#endif
