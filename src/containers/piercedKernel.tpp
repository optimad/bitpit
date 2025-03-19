/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
* Destructor
*/
template<typename id_t>
PiercedKernel<id_t>::~PiercedKernel()
{
    for (PiercedStorageSyncSlave<id_t> *storage : getStorages()) {
        storage->unsetKernel();
    }
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
* \param release if set to true, the containers that hold kernel data are
* requested to release all unneeded memory (it is a non-binding request)
*/
template<typename id_t>
typename PiercedKernel<id_t>::ClearAction PiercedKernel<id_t>::clear(bool release)
{
    ClearAction syncAction = _clear(release);

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
}


/**
* Removes all elements from the kernel.
*
* The function will NOT process the sync action.
*
* \param release if set to true, the containers that hold holes data are
* requested to release all unneeded memory (it is a non-binding request)
*/
template<typename id_t>
typename PiercedKernel<id_t>::ClearAction PiercedKernel<id_t>::_clear(bool release)
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

    // Generate the sync action
    ClearAction syncAction;

    return syncAction;
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
typename PiercedKernel<id_t>::ReserveAction PiercedKernel<id_t>::reserve(std::size_t n)
{
    ReserveAction syncAction = _reserve(n);

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
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
* The function will NOT process the sync action.
*
* \param n the minimum capacity requested for the kernel
*/
template<typename id_t>
typename PiercedKernel<id_t>::ReserveAction PiercedKernel<id_t>::_reserve(std::size_t n)
{
    // Update the kernel
    m_ids.reserve(n);
    m_pos.reserve(n);

    // Generate the sync action
    ReserveAction syncAction;
    syncAction.info[PiercedSyncAction::INFO_SIZE] = n;

    return syncAction;
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
typename PiercedKernel<id_t>::ResizeAction PiercedKernel<id_t>::resize(std::size_t n)
{
    ResizeAction syncAction = _resize(n);

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
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
* The function will NOT process the sync action.
*
* \param n is the new kernel size, expressed in number of elements.
*/
template<typename id_t>
typename PiercedKernel<id_t>::ResizeAction PiercedKernel<id_t>::_resize(std::size_t n)
{
    // If the size of the kernel is already the requested size there is
    // nothing to do.
    if (n == size()) {
        ResizeAction syncAction(ResizeAction::TYPE_NOOP);

        return syncAction;
    }

    // A request for a size equal to 0 is equivalent to a clear.
    if (n == 0) {
        _clear();

        ResizeAction syncAction(ResizeAction::TYPE_CLEAR);

        return syncAction;
    }

    // If the requested size is greater that the current size we may need to
    // reserve space to allow the kernel to reach the requested size.
    if (n > size()) {
        _reserve(n);

        ResizeAction syncAction(ResizeAction::TYPE_RESERVE);
        syncAction.info[PiercedSyncAction::INFO_SIZE] = rawSize();

        return syncAction;
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
    rawShrink(last_used_pos + 1);

    // Generate the sync action
    ResizeAction syncAction(ResizeAction::TYPE_RESIZE);
    syncAction.info[PiercedSyncAction::INFO_SIZE] = rawSize();

    return syncAction;
}

/**
* Sorts the elements of the kernel in ascending id order.
*/
template<typename id_t>
typename PiercedKernel<id_t>::SortAction PiercedKernel<id_t>::sort()
{
    SortAction syncAction = _sort(m_begin_pos, m_end_pos);

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
}

/**
* Sorts the kernel after the element with the reference id in ascending
* id order.
*
* \param referenceId is the id of the element after which the kernel will
* be sorted
* \param inclusive if true the reference element will be sorted, otherwise
* the sorting will stop at the element following the reference
*/
template<typename id_t>
typename PiercedKernel<id_t>::SortAction PiercedKernel<id_t>::sortAfter(id_t referenceId, bool inclusive)
{
    // Get the reference position
    std::size_t referencePos = getPos(referenceId);
    if (!inclusive) {
        referencePos++;
    }

    // Sorte the
    SortAction syncAction = _sort(referencePos, m_end_pos);

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
}

/**
* Sorts the kernel before the element with the reference id in ascending
* id order.
*
* \param referenceId is the id of the element before which the kernel will
* be sorted
* \param inclusive if true the reference element will be sorted, otherwise
* the sorting will stop at the element preceding the reference
*/
template<typename id_t>
typename PiercedKernel<id_t>::SortAction PiercedKernel<id_t>::sortBefore(id_t referenceId, bool inclusive)
{
    // Get the reference position
    std::size_t referencePos = getPos(referenceId);
    if (inclusive) {
        referencePos++;
    }

    // Sort the container
    SortAction syncAction = _sort(m_begin_pos, referencePos);

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
}

/**
* Sorts the elements of the kernel in ascending id order.
*
* The function will NOT process the sync action.
*
* \param beginPos is the first position that will be sorted
* \param endPos is the position past the last element that will be sorted
*/
template<typename id_t>
typename PiercedKernel<id_t>::SortAction PiercedKernel<id_t>::_sort(std::size_t beginPos, std::size_t endPos)
{
    // Squeeze the kernel
    //
    // After the squeeze there will be no holes in the kernel.
    SqueezeAction squeezeAction = _squeeze();
    const std::vector<std::size_t> &squeezePermutations = *(squeezeAction.data);
    bool squeezed = (squeezeAction.data.get() != nullptr);

    // Get updated kernel size
    std::size_t updatedKernelRawSize = rawSize();

    // Update the sort range
    if (squeezed) {
        std::size_t previousBeginPos = beginPos;
        std::size_t previousEndPos   = endPos;
        for (std::size_t i = m_begin_pos; i < m_end_pos; ++i) {
            std::size_t previousPos = squeezePermutations[i];
            if (previousPos == previousBeginPos) {
                beginPos = i;
            }
            if ((previousPos + 1) == previousEndPos) {
                endPos = i + 1;
                break;
            }
        }
    }

    // Evaluates the sort permutations
    std::vector<std::size_t> sortPermutations(updatedKernelRawSize);
    for (std::size_t i = 0; i < updatedKernelRawSize; ++i) {
        sortPermutations[i] = i;
    }

    std::sort(sortPermutations.begin() + beginPos, sortPermutations.begin() + endPos, idLess(m_ids));

    // Create the sync action
    SortAction syncAction;
    syncAction.info[PiercedSyncAction::INFO_SIZE] = updatedKernelRawSize;
    if (squeezed) {
        syncAction.importData(std::vector<std::size_t>(squeezePermutations));
        std::vector<std::size_t> &permutations = *(syncAction.data);
        for (std::size_t i = beginPos; i < endPos; ++i) {
            permutations[i] = squeezePermutations[sortPermutations[i]];
        }
    } else {
        syncAction.importData(std::vector<std::size_t>(sortPermutations));
    }

    // Sort the ids
    //
    // NOTE: the reorder function will destroy the permutation vector
    utils::reorderVector<id_t>(sortPermutations, m_ids, updatedKernelRawSize);

    // Update the positions
    m_pos.clear();
    for (std::size_t i = 0; i < updatedKernelRawSize; ++i) {
        m_pos[m_ids[i]] = i;
    }

    // Return the permutations
    return syncAction;
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
typename PiercedKernel<id_t>::SqueezeAction PiercedKernel<id_t>::squeeze()
{
    SqueezeAction syncAction = _squeeze();

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
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
*
* The function will NOT process the sync action.
*/
template<typename id_t>
typename PiercedKernel<id_t>::SqueezeAction PiercedKernel<id_t>::_squeeze()
{
    // Flush changes
    flush();

    // Sort the holes
    //
    // After being flushed, the container will contain only regular holes.
    holesSortRegular();

    // Get kernel size
    std::size_t kernelSize    = size();
    std::size_t kernelRawSize = rawSize();

    // Initialize the sync action
    SqueezeAction syncAction;
    syncAction.info[PiercedSyncAction::INFO_SIZE] = kernelSize;

    // Compact the kernel
    std::size_t nHoles = holesCount();
    if (nHoles != 0) {
        // Initialize sync data
        syncAction.importData(std::vector<std::size_t>(kernelRawSize));
        std::vector<std::size_t> &permutations = *(syncAction.data);

        // Move the elements
        std::size_t firstPosToUpdate;
        if (m_begin_pos == 0) {
            firstPosToUpdate = m_holes[m_holes_regular_end - 1];
        } else {
            firstPosToUpdate = 0;
        }

        for (std::size_t pos = 0; pos < firstPosToUpdate; ++pos) {
            permutations[pos] = pos;
        }

        std::size_t offset = 0;
        for (std::size_t pos = firstPosToUpdate; pos < m_end_pos; pos++) {
            if (offset < nHoles && m_holes[m_holes_regular_end - offset - 1] == pos) {
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

        for (std::size_t pos = m_end_pos; pos < rawSize(); ++pos) {
            permutations[pos] = pos;
        }

        // Clear the holes
        holesClear(true);

        // Reset begin position
        setBeginPos(0);

        // Shrink the kernel
        rawShrink(kernelSize);
    }

    // Shrink to fit
    _shrinkToFit();

    return syncAction;
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
typename PiercedKernel<id_t>::ShrinkToFitAction PiercedKernel<id_t>::shrinkToFit()
{
    ShrinkToFitAction syncAction = _shrinkToFit();

    // Update the storage
    processSyncAction(syncAction);

    return syncAction;
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
*
* The function will NOT process the sync action.
*/
template<typename id_t>
typename PiercedKernel<id_t>::ShrinkToFitAction PiercedKernel<id_t>::_shrinkToFit()
{
    // Update the kernel
    m_ids.shrink_to_fit();

    // Generate the sync action
    ShrinkToFitAction syncAction;

    return syncAction;
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
* \param other is another kernel of the same type (i.e., instantiated with the
* same template parameters) whose content is swapped with that of this kernel.
*/
template<typename id_t>
void PiercedKernel<id_t>::swap(PiercedKernel &other) noexcept
{
    PiercedSyncMaster::swap(other);

    std::swap(other.m_begin_pos, m_begin_pos);
    std::swap(other.m_end_pos, m_end_pos);
    std::swap(other.m_dirty_begin_pos, m_dirty_begin_pos);
    std::swap(other.m_ids, m_ids);
    std::swap(other.m_pos, m_pos);
    std::swap(other.m_holes, m_holes);
    std::swap(other.m_holes_regular_begin, m_holes_regular_begin);
    std::swap(other.m_holes_regular_end, m_holes_regular_end);
    std::swap(other.m_holes_regular_sorted, m_holes_regular_sorted);
    std::swap(other.m_holes_pending_begin, m_holes_pending_begin);
    std::swap(other.m_holes_pending_end, m_holes_pending_end);
    std::swap(other.m_holes_pending_sorted, m_holes_pending_sorted);
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
* Dumps to screen the internal data.
*/
template<typename id_t>
void PiercedKernel<id_t>::dump() const
{
    std::cout << "----------------[ DUMP ]----------------" << std::endl;

    std::cout << std::endl;
    std::cout << " size: " << size() << std::endl;

    std::cout << " m_holes_regular_begin: " << m_holes_regular_begin << std::endl;
    std::cout << " m_holes_regular_end  : " << m_holes_regular_end << std::endl;

    std::cout << std::endl;
    std::cout << " Regular holes: " << std::endl;
    for (std::size_t k = m_holes_regular_begin; k < m_holes_regular_end; ++k) {
        std::cout << m_holes[k] << std::endl;
    }

    std::cout << std::endl;
    std::cout << " m_holes_pending_begin: " << m_holes_pending_begin << std::endl;
    std::cout << " m_holes_pending_end  : " << m_holes_pending_end << std::endl;

    std::cout << std::endl;
    std::cout << " Pending holes" << std::endl;
    for (std::size_t k = m_holes_pending_begin; k < m_holes_pending_end; ++k) {
        std::cout << m_holes[k] << std::endl;
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
*/
template<typename id_t>
void PiercedKernel<id_t>::checkIntegrity() const
{
    // Check if the elements and their position match
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

    // Check if the number of elements visited by the iterator matches the size of the kernel
    std::size_t n = 0;
    for (auto it = begin(); it != end(); ++it) {
        ++n;
    }

    if (n != size()) {
        std::cout << " The number of elements visited by the iterator doesn't match the size of the kernel" << std::endl;
        std::cout << " The iterator has visited " << n << "elements" << std::endl;
        std::cout << " The size of the kernel is " << size() << std::endl;
        throw std::runtime_error("Integrity check error");
    }

    // Check the consistency of the ids of elements associated with empty positions
    //
    // The id of an empty element contains the distance, measured in number of elements, between
    // the current element and the next element the iterator should jump to (the distance is
    // negative). If the id of an element is lower than -1, the ids of the following elements
    // should be consecutive until they reach -1.
    if (!empty()) {
        std::size_t pos = 0;
        while (pos < m_ids.size() - 1) {
            for (std::size_t i = pos; i < m_ids.size(); ++i) {
                // Skip non-empty positions
                if (m_ids[i] >= 0) {
                    pos = i + 1;
                    break;
                }

                // Skip empty positions with a distance from a non-empty element equal to one
                if (m_ids[i + 1] >= 0 && m_ids[i] == - 1) {
                    pos = i + 2;
                    break;
                }

                // Check if ids are consecutive
                if ((m_ids[i] - m_ids[i + 1] != -1) && m_ids[i] != -1) {
                    std::cout << " The ids associated with empty positions are not consecutive" << std::endl;
                    std::cout << " Position " << i << " is associated with id " << m_ids[i] << std::endl;
                    std::cout << " Position " << (i + 1) << " is associated with id " << m_ids[i + 1] << std::endl;
                    throw std::runtime_error("Integrity check error");
                }
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

/**
* Searches the container for elements with the specified id and returns the
* number of matches.
*
* Because all ids in a pierced kernel are unique, the function can only return
* 1 (if the element is found) or zero (otherwise).
*
* \param id is the id to look for
* \result Returns 1 if the pierced kernel contains an element wit the specified
* id, or zero otherwise.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::count(id_t id) const
{
    return m_pos.count(id);
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
std::size_t PiercedKernel<id_t>::evalFlatIndex(id_t id) const
{
    std::size_t pos  = getPos(id);

    // Initialize flat id with the position of the element
    std::size_t flat = pos;

    // Subtract pending holes before position
    if (holesCountPending() > 0) {
        holes_const_iterator pending_begin_itr = m_holes.cbegin() + m_holes_pending_begin;
        holes_const_iterator pending_end_itr   = m_holes.cbegin() + m_holes_pending_end;

        std::size_t nHolesBefore = 0;
        for (auto itr = pending_begin_itr; itr != pending_end_itr; ++itr) {
            if (*itr <= pos) {
                ++nHolesBefore;
            }
        }

        flat -= nHolesBefore;
    }

    // Subtract regular holes before position
    if (holesCountRegular() > 0) {
        holes_const_iterator regular_begin_itr = m_holes.cbegin() + m_holes_regular_begin;
        holes_const_iterator regular_end_itr   = m_holes.cbegin() + m_holes_regular_end;

        std::size_t nHolesBefore = 0;
        for (auto itr = regular_begin_itr; itr != regular_end_itr; ++itr) {
            if (*itr <= pos) {
                ++nHolesBefore;
            }
        }

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
    holes_iterator regular_begin_itr = m_holes.begin() + m_holes_regular_begin;
    holes_iterator regular_hole_itr  = m_holes.begin() + m_holes_regular_end;

    holes_iterator pending_begin_itr = m_holes.begin() + m_holes_pending_begin;
    holes_iterator pending_hole_itr  = m_holes.begin() + m_holes_pending_end;

    std::size_t nEmpties  = 0;
    std::size_t markerPos = targetSize;
    while (true) {
        if (isPosEmpty(markerPos)) {
            markerPos = findNextUsedPos(markerPos - 1);
        }

        // Count the number of holes and pending deletes before the
        // current marker position
        if (regular_hole_itr != regular_begin_itr) {
            holes_iterator itr_previous = regular_hole_itr;
            regular_hole_itr = std::upper_bound(regular_begin_itr, regular_hole_itr, markerPos, std::greater<std::size_t>());
            nEmpties += std::distance(regular_hole_itr, itr_previous);
        }

        if (pending_hole_itr != pending_begin_itr) {
            holes_iterator itr_previous = pending_hole_itr;
            pending_hole_itr = std::upper_bound(pending_begin_itr, pending_hole_itr, markerPos, std::greater<std::size_t>());
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

/**
* Gets a constant iterator pointing to the specified element.
*
* \param id is the id
* \result A constant iterator pointing to the specified element.
*/
template<typename id_t>
typename PiercedKernel<id_t>::const_iterator PiercedKernel<id_t>::find(const id_t &id) const noexcept
{
    auto iterator = m_pos.find(id);
    if (iterator != m_pos.end()) {
        return rawFind(iterator->second);
    } else {
        return end();
    }
}

/**
* Gets a constant iterator pointing to the specified position.
*
* \param pos is the position
* \result A constant iterator pointing to the specified position.
*/
template<typename id_t>
typename PiercedKernel<id_t>::const_iterator PiercedKernel<id_t>::rawFind(std::size_t pos) const noexcept
{
    return const_iterator(this, pos);
}

/*!
* Returns a constant iterator pointing to the first element in the vector.
*
* \result A constant iterator pointing to the first element in the vector.
*/
template<typename id_t>
typename PiercedKernel<id_t>::const_iterator PiercedKernel<id_t>::begin() const noexcept
{
    return cbegin();
}

/*!
* Returns a constant iterator referring to the past-the-end element in the
* vector.
*
* \result A constant iterator referring to the past-the-end element in the
* vector.
*/
template<typename id_t>
typename PiercedKernel<id_t>::const_iterator PiercedKernel<id_t>::end() const noexcept
{
    return cend();
}

/*!
* Returns an conts_iterator pointing to the first element in the vector.
*
* \result A const_iterator pointing to the first element in the vector.
*/
template<typename id_t>
typename PiercedKernel<id_t>::const_iterator PiercedKernel<id_t>::cbegin() const noexcept
{
    return rawFind(m_begin_pos);
}

/*!
* Returns an const_iterator referring to the past-the-end element in the
* vector.
*
* \result A const_iterator referring to the past-the-end element in the vector.
*/
template<typename id_t>
typename PiercedKernel<id_t>::const_iterator PiercedKernel<id_t>::cend() const noexcept
{
    return rawFind(m_end_pos);
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
* Gets the position of the last element in the container.
*
* If there is no element with the specified id, an exception is thrown.
*
* \result The position of the last element in the container.
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
        std::size_t hole = m_holes_pending_begin;
        if (m_holes[hole] > referencePos) {
            return fillHole(hole, id);
        }
    }

    // Check if there is a regular hole that can be filled
    if (m_holes_regular_begin != m_holes_regular_end) {
        // Sort regular holes
        holesSortRegular();

        // The first hole is the one with the highest position, if the
        // position of this hole is greater than the reference position
        // it is possible to use it.
        std::size_t hole = m_holes_regular_begin;
        if (m_holes[hole] > referencePos) {
            return fillHole(hole, id);
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
        std::size_t hole = m_holes_pending_end - 1;
        if (m_holes[hole] < referencePos) {
            return fillHole(hole, id);
        }
    }

    // Check if there is a regular hole that can be filled
    if (m_holes_regular_begin != m_holes_regular_end) {
        // Sort regular holes
        holesSortRegular();

        // The last hole is the one with the lowest position, if the
        // position of this hole is greater than the reference position
        // it is possible to use it.
        std::size_t hole = m_holes_regular_end - 1;
        if (m_holes[hole] < referencePos) {
            return fillHole(hole, id);
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
    m_ids.emplace_back(id);

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
* \param hole is the position of the hole to fill
* \param id is the id that will be associated to the position
* \result The synchronization action associated with the fill.
*/
template<typename id_t>
typename PiercedKernel<id_t>::FillAction PiercedKernel<id_t>::fillHole(std::size_t hole, id_t id)
{
    // Validate the id
    validateId(id);

    // Get the position of the hole
    std::size_t pos = m_holes[hole];

    // Fill the hole
    setPosId(pos, id);

    // Remove the holes from the kernel
    if (hole >= m_holes_pending_begin) {
        int nPendings = holesCountPending();
        if (nPendings > 1) {
            if (hole == (m_holes_pending_end - 1)) {
                --m_holes_pending_end;
            } else if (hole == m_holes_pending_begin) {
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
            if (hole == (m_holes_regular_end - 1)) {
                --m_holes_regular_end;
            } else if (hole == m_holes_regular_begin) {
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
        holes_iterator regular_begin_itr = m_holes.begin() + m_holes_regular_begin;
        holes_iterator regular_end_itr   = m_holes.begin() + m_holes_regular_end;

        holes_iterator update_begin_itr = regular_begin_itr;
        holes_iterator update_end_itr   = upper_bound(regular_begin_itr, regular_end_itr, pos, std::greater<std::size_t>());
        for (auto itr = update_begin_itr; itr != update_end_itr; itr++) {
            (*itr)++;
        }
    }

    // Update the pending holes
    if (m_holes_pending_begin != m_holes_pending_end) {
        holes_iterator pending_begin_itr = m_holes.begin() + m_holes_pending_begin;
        holes_iterator pending_end_itr   = m_holes.begin() + m_holes_pending_end;

        holes_iterator update_begin_itr = pending_begin_itr;
        holes_iterator update_end_itr   = upper_bound(pending_begin_itr, pending_end_itr, pos, std::greater<std::size_t>());
        for (auto itr = update_begin_itr; itr != update_end_itr; itr++) {
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
    // The kernel cannt be empty
    if (empty()) {
        throw std::out_of_range("Vector is empty");
    }

    // Erase the last element
    std::size_t updatedRawSize;
    if (size() == 1) {
        updatedRawSize = 0;
    } else {
        updatedRawSize = findPrevUsedPos(m_end_pos - 1) + 1;
    }

    rawShrink(updatedRawSize);

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
        std::size_t updatedRawSize;
        if (size() == 1) {
            updatedRawSize = 0;
        } else {
            updatedRawSize = findPrevUsedPos(m_end_pos - 1) + 1;
        }

        rawShrink(updatedRawSize);

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
    if (m_holes_pending_end == m_holes.size()) {
        holesFlush();
    }

    // Add the hole at the end of the pending holes
    m_holes[m_holes_pending_end] = pos;
    m_holes_pending_end++;

    // Check if pending holes are still sorted
    if (m_holes_pending_sorted) {
        std::size_t nPendings = holesCountPending();
        if (nPendings > 1 && (m_holes[m_holes_pending_end - 1] > m_holes[m_holes_pending_end - 2])) {
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
* \param release if set to true, the containers that hold holes data are
* requested to release all unneeded memory (it is a non-binding request)
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
* \param release if set to true, the containers that hold holes data are
* requested to release all unneeded memory (it is a non-binding request)
*/
template<typename id_t>
void PiercedKernel<id_t>::holesClearRegular(bool release)
{
    // Release the memory
    if (release) {
        m_holes.shrink_to_fit();
    }

    // Reset regulr holes iterators
    m_holes_regular_begin = 0;
    m_holes_regular_end   = m_holes_regular_begin;

    // There are no holes, therefore all holes are sorted
    m_holes_regular_sorted = true;
}

/**
* Clear pending holes
*
* \param release if set to true, the containers that hold holes data are
* requested to release all unneeded memory (it is a non-binding request)
*/
template<typename id_t>
void PiercedKernel<id_t>::holesClearPending(bool release)
{
    // Clear section of the kernel associated with the pending holes
    long offset    = m_holes_regular_begin;
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
* \param release if set to true, the containers that hold holes data are
* requested to release all unneeded memory (it is a non-binding request)
*/
template<typename id_t>
void PiercedKernel<id_t>::holesResize(std::size_t offset, std::size_t nRegulars, std::size_t nPendings, bool release)
{
    m_holes.resize(offset + nRegulars + MAX_PENDING_HOLES);
    if (release) {
        m_holes.shrink_to_fit();
    }

    m_holes_regular_begin = offset;
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
    return (m_holes_pending_end - m_holes_pending_begin);
}

/**
* Count the regular holes.
*
* \result The number of regular holes.
*/
template<typename id_t>
std::size_t PiercedKernel<id_t>::holesCountRegular() const
{
    return (m_holes_regular_end - m_holes_regular_begin);
}

/**
* Flushes the list of available holes.
*
* All the pending hole are converted to regular holes and new
* space is reserved for future pending holes.
*
* This function will not sort the regular holes.
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

    std::size_t hole = m_holes_pending_begin;
    std::size_t pos = m_end_pos;
    do {
        if (m_holes[hole] >= pos) {
            hole++;
            continue;
        }

        pos = m_holes[hole];
        std::size_t next_used_pos = findNextUsedPos(pos);
        do {
            setPosEmptyId(pos, next_used_pos);
            if (pos > 0) {
                pos--;
            } else {
                break;
            }
        } while (isPosEmpty(pos));
    } while (pos > 0 && hole != m_holes_pending_end);

    // Move the pending holes into the list of regular holes
    for (std::size_t hole = m_holes_pending_begin; hole < m_holes_pending_end; ++hole) {
        std::size_t pos = m_holes[hole];

        // If there is space available at the beginning of the holes, try
        // using pending holes to fill that gap.
        if (m_holes_regular_begin != 0) {
            --m_holes_regular_begin;
            m_holes[m_holes_regular_begin] = pos;

            // Regular holes are no more sorted
            if (m_holes_regular_sorted) {
                m_holes_regular_sorted = false;
            }
        } else {
            if (hole != m_holes_regular_end) {
                m_holes[m_holes_regular_end] = pos;
            }
            ++m_holes_regular_end;
            ++m_holes_pending_begin;

            // Check if regular holes are still sorted
            if (m_holes_regular_sorted) {
                std::size_t nRegulars = holesCountRegular();
                if (nRegulars > 1 && (m_holes[m_holes_regular_end - 1] > m_holes[m_holes_regular_end - 2])) {
                    m_holes_regular_sorted = false;
                }
            }
        }
    }

    // Move the holes at the beginning of the vector
    std::size_t nRegulars = holesCountRegular();
    if (nRegulars != 0 && m_holes_regular_begin != 0) {
        for (std::size_t k = 0; k < nRegulars; ++k) {
            m_holes[k] = m_holes[m_holes_regular_begin + k];
        }

        m_holes_regular_begin = 0;
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

    holes_iterator pending_begin_itr = m_holes.begin() + m_holes_pending_begin;
    holes_iterator pending_end_itr   = m_holes.begin() + m_holes_pending_end;
    std::sort(pending_begin_itr, pending_end_itr, std::greater<std::size_t>());
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

    holes_iterator regular_begin_itr = m_holes.begin() + m_holes_regular_begin;
    holes_iterator regular_end_itr   = m_holes.begin() + m_holes_regular_end;
    std::sort(regular_begin_itr, regular_end_itr, std::greater<std::size_t>());
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
* Get the storages assocaited with the kernel.
*
* \result The storages assocaited with the kernel.
*/
template<typename id_t>
std::vector<PiercedStorageSyncSlave<id_t> *> PiercedKernel<id_t>::getStorages()
{
    std::size_t storageIdx = 0;
    std::vector<PiercedStorageSyncSlave<id_t> *> storges(m_slaves.size());
    for (auto &slaveEntry : m_slaves) {
        PiercedSyncSlave *slave = slaveEntry.first;
        PiercedStorageSyncSlave<id_t> *piercedStorage = dynamic_cast<PiercedStorageSyncSlave<id_t> *>(slave);
        if (piercedStorage) {
            storges[storageIdx] = piercedStorage;
            ++storageIdx;
        }
    }

    return storges;
}

/**
* Get the storages assocaited with the kernel.
*
* \result The storages assocaited with the kernel.
*/
template<typename id_t>
std::vector<const PiercedStorageSyncSlave<id_t> *> PiercedKernel<id_t>::getStorages() const
{
    std::size_t storageIdx = 0;
    std::vector<PiercedStorageSyncSlave<id_t> *> storges(m_slaves.size());
    for (const auto &slaveEntry : m_slaves) {
        const PiercedSyncSlave *slave = slaveEntry.first;
        const PiercedStorageSyncSlave<id_t> *piercedStorage = dynamic_cast<const PiercedStorageSyncSlave<id_t> *>(slave);
        if (piercedStorage) {
            storges[storageIdx] = piercedStorage;
            ++storageIdx;

        }
    }

    return storges;
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
* The id of an empty element contains the distance, measured in number of elements, between
* the current element and the next element the iterator should jump to (the distance is
* negative).
*
* The next used position has to be greater than the position to update. Otherwise, undefined
* behavior occurs.
*
* \param pos is the position to update
* \param nextUsedPos is the position of the next non-empty element
*/
template<typename id_t>
void PiercedKernel<id_t>::setPosEmptyId(std::size_t pos, std::size_t nextUsedPos)
{
    assert(nextUsedPos > pos);

    // Early return if we are updating the last position
    //
    // When updating the last position, the id can only be set to point to the next position
    // (i.e., set equal to -1).
    if (pos == (m_end_pos - 1)) {
        m_ids[pos] = id_t(-1);
        return;
    }

    // Set the id of the position
    //
    // The id of an empty element contains the distance, measured in number of elements, between
    // the current element and the next element the iterator should jump to (the distance is
    // negative). If the next position is non-empty or a regular hole, the value of the id will
    // be set in such a way to point to the next non-empty element (this will make the iterator
    // as fast as possible). If the next position is a pending hole, we need to set the id to
    // -1, otherwise we may break the iterator. For example,
    //
    //      ID    1  2  3 -1  4 -1 -1  5
    //     POS    0  1  2  3  4  5  6  7
    //
    // when setting the id associated with position 4, it should be set to -1 and not to -3,
    // otherwise if the position 5 will be filled, the iterator will not able to reach it (it
    // will reach position 4 and will wrongly skip to position 7).
    id_t offset = id_t(pos - nextUsedPos);
    if (m_ids[pos + 1] == offset + 1) {
        m_ids[pos] = offset;
    } else {
        m_ids[pos] = id_t(-1);
    }
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
void PiercedKernel<id_t>::rawShrink(std::size_t n)
{
    // We can only shrink the kernel
    if (n > m_end_pos) {
        throw std::out_of_range("The kernel can only be shrunk");
    }

    // Check if we actually need to shrink the kernel
    if (n == m_end_pos) {
        return;
    }

    // When the new last position is before the first one this is equivalent
    // to a clear
    if (n < (m_begin_pos + 1)) {
        _clear();
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

    // Update the holes
    if (holesCount() != 0) {
        // Remove regular holes beyond the updated last position
        holesSortRegular();
        holes_iterator regular_begin_itr = m_holes.begin() + m_holes_regular_begin;
        holes_iterator regular_end_itr   = m_holes.begin() + m_holes_regular_end;
        regular_begin_itr = std::lower_bound(regular_begin_itr, regular_end_itr, m_end_pos - 1, std::greater<std::size_t>());
        m_holes_regular_begin = std::distance(m_holes.begin(), regular_begin_itr);
        if (m_holes_regular_begin == m_holes_regular_end) {
            m_holes_regular_begin = 0;
            m_holes_regular_end   = m_holes_regular_begin;
        }

        // Remove pending holes beyond the updated last position
        holesSortPending();
        holes_iterator pending_begin_itr = m_holes.begin() + m_holes_pending_begin;
        holes_iterator pending_end_itr   = m_holes.begin() + m_holes_pending_end;
        pending_begin_itr = std::lower_bound(pending_begin_itr, pending_end_itr, m_end_pos - 1, std::greater<std::size_t>());
        m_holes_pending_begin = std::distance(m_holes.begin(), pending_begin_itr);
        if (m_holes_pending_begin == m_holes_pending_end) {
            m_holes_pending_begin = m_holes_regular_end;
            m_holes_pending_end   = m_holes_pending_begin;
        }
    }
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
        id_t id;
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
    std::size_t nHoles;
    utils::binary::read(stream, nHoles);
    m_holes.resize(nHoles);
    for (std::size_t n = 0; n < nHoles; ++n) {
        std::size_t pos;
        utils::binary::read(stream, pos);

        m_holes[n] = pos;
    }

    utils::binary::read(stream, m_holes_regular_begin);
    utils::binary::read(stream, m_holes_regular_end);
    utils::binary::read(stream, m_holes_pending_begin);
    utils::binary::read(stream, m_holes_pending_end);
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
    std::size_t nIds = m_pos.size();
    utils::binary::write(stream, nIds);
    for (const auto &entry : m_pos) {
        utils::binary::write(stream, entry.first);
        utils::binary::write(stream, entry.second);
    }

    // Postions data
    std::size_t nPositions = m_ids.size();
    utils::binary::write(stream, nPositions);
    for (std::size_t pos : m_ids) {
        utils::binary::write(stream, pos);
    }

    utils::binary::write(stream, m_begin_pos);
    utils::binary::write(stream, m_end_pos);
    utils::binary::write(stream, m_dirty_begin_pos);

    // Holes data
    std::size_t nHoles = m_holes.size();
    utils::binary::write(stream, nHoles);
    for (std::size_t hole : m_holes) {
        utils::binary::write(stream, hole);
    }

    utils::binary::write(stream, m_holes_regular_begin);
    utils::binary::write(stream, m_holes_regular_end);
    utils::binary::write(stream, m_holes_pending_begin);
    utils::binary::write(stream, m_holes_pending_end);
    utils::binary::write(stream, m_holes_regular_sorted);
    utils::binary::write(stream, m_holes_pending_sorted);

    // Synchronization data
    PiercedSyncMaster::dump(stream);
}

}

#endif
