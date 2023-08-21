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

#define __BITPIT_LEVELSET_CACHED_SRC__

# include "levelSetCachedObject.hpp"

# include "bitpit_patchkernel.hpp"

namespace bitpit {

// Explicit instantization
template class LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>;
template class LevelSetNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>;
template class LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>;

template class LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
template class LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>>;
template class LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>;

template class LevelSetCachedObject<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
template class LevelSetCachedObject<LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>>;
template class LevelSetCachedObject<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>;

/*!
 * Constructor
 *
 * \param kernel is the container associated with the storage manager
 */
LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::LevelSetNarrowBandCache(Kernel *kernel)
    : LevelSetExternalPiercedStorageManager(kernel, KERNEL_SYNC_MODE_AUTOMATIC),
      LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>()
{
    m_values    = addStorage<double>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_JOURNALED);
    m_gradients = addStorage<std::array<double, 3>>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_JOURNALED);

    m_narrowBandFlag = addStorage<char>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_JOURNALED);
    m_narrowBandFlag->fill(0);
}

/*!
 * Insert a kernel entry with the specified id.
 *
 * This function cannot be used.
 *
 * \param id is id of the entry
 * \param sync controls if the storages will be synched
 * \result A kernel iterator pointing to the newly created entry.
 */
LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::KernelIterator LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::insert(long id, bool sync)
{
    BITPIT_UNUSED(sync);

    KernelIterator iterator = m_kernel->find(id);
    std::size_t rawId = iterator.getRawIndex();

    m_narrowBandFlag->rawAt(rawId) = 1;

    return iterator;
}

/*!
 * Erase from the kernel the entry with the specified id.
 *
 * This function cannot be used.
 *
 * \param id is id of the cell
 * \param sync controls if the kernel will be synched
 */
void LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::erase(long id, bool sync)
{
    BITPIT_UNUSED(sync);

    KernelIterator iterator = m_kernel->find(id);
    std::size_t rawId = iterator.getRawIndex();

    m_narrowBandFlag->rawAt(rawId) = 0;
}

/*!
 * Checks if the kernel contains an entry with the specified id.
 *
 * \param id is id of the entry
 * \result Return true if kernel contains an entry with the specified id,
 * false otherwise.
 */
bool LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::contains(long id) const
{
    KernelIterator iterator = m_kernel->find(id);
    std::size_t rawId = iterator.getRawIndex();

    return (m_narrowBandFlag->rawAt(rawId) == 1);
}

/*!
 * Get a kernel iterator for the entry with the specified id.
 *
 * \param id is id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::KernelIterator LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::find(long id) const
{
    KernelIterator iterator = m_kernel->find(id);
    std::size_t rawId = iterator.getRawIndex();
    if (m_narrowBandFlag->rawAt(rawId) == 0) {
        return m_kernel->end();
    }

    return iterator;
}

/*!
 * Get a kernel iterator for the entry with the specified raw id.
 *
 * \param rawId is raw id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::KernelIterator LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::rawFind(std::size_t rawId) const
{
    if (m_narrowBandFlag->rawAt(rawId) == 0) {
        return m_kernel->end();
    }

     return m_kernel->rawFind(rawId);
}

/*!
 * Get a reference to the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result A reference to the levelset value of the specified entry.
 */
double & LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::getValue(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_values->rawAt(rawId);
}

/*!
 * Get the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset value of the specified entry.
 */
double LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::getValue(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_values->rawAt(rawId);
}

/*!
 * Get a reference to the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result A reference to the levelset gradient of the specified entry.
 */
std::array<double, 3> & LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::getGradient(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_gradients->rawAt(rawId);
}

/*!
 * Get the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset gradient of the specified entry.
 */
const std::array<double, 3> & LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::getGradient(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_gradients->rawAt(rawId);
}

/*!
 * Clear the kernel of the storage manager.
 */
void LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::clearKernel()
{
    LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>::clearKernel();

    m_narrowBandFlag->fill(0);
}

/*!
 * Dump the kernel of the storage manager.
 *
 * \param stream is the output stream
 */
void LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::dumpKernel(std::ostream &stream)
{
    LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>::dumpKernel(stream);

    m_narrowBandFlag->dump(stream);
}

/*!
 * Restore the kernel of the storage manager.
 *
 * \param stream is the input stream
 */
void LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::restoreKernel(std::istream &stream)
{
    LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>::restore(stream);

    m_narrowBandFlag->restore(stream);
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>::swap(LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager> &other) noexcept
{
    LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>::swap(other);

    std::swap(other.m_narrowBandFlag, m_narrowBandFlag);
}

/*!
 * Constructor
 */
LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>::LevelSetNarrowBandCache()
    : LevelSetInternalPiercedStorageManager(), LevelSetNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>()
{
    // It is faster to use a concurrent synchronization because items will be added/removed
    // to the kernel one at the time.
    m_values    = addStorage<double>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_CONCURRENT);
    m_gradients = addStorage<std::array<double, 3>>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_CONCURRENT);
}

/*!
 * Get the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset value of the specified entry.
 */
double LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>::getValue(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_values->rawAt(rawId);
}

/*!
 * Get a reference to the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result A reference to the levelset value of the specified entry.
 */
double & LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>::getValue(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_values->rawAt(rawId);
}

/*!
 * Get the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset gradient of the specified entry.
 */
const std::array<double, 3> & LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>::getGradient(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_gradients->rawAt(rawId);
}

/*!
 * Get a reference to the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result A reference to the levelset gradient of the specified entry.
 */
std::array<double, 3> & LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>::getGradient(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_gradients->rawAt(rawId);
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>::swap(LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager> &other) noexcept
{
    LevelSetNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>::swap(other);
}

/*!
 * Constructor
 *
 * \param nItems are the maximum number of items the cache will hold
 */
LevelSetNarrowBandCache<LevelSetDirectStorageManager>::LevelSetNarrowBandCache(std::size_t nItems)
    : LevelSetDirectStorageManager(nItems), LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>()
{
    m_values    = addStorage<double>(getStorageCount());
    m_gradients = addStorage<std::array<double, 3>>(getStorageCount());

    m_narrowBandFlag = addStorage<char>(getStorageCount());
    std::fill(m_narrowBandFlag->begin(), m_narrowBandFlag->end(), 0);
}

/*!
 * Insert a kernel entry with the specified id.
 *
 * This function cannot be used.
 *
 * \param id is id of the entry
 * \param sync controls if the storages will be synched
 * \result A kernel iterator pointing to the newly created entry.
 */
LevelSetNarrowBandCache<LevelSetDirectStorageManager>::KernelIterator LevelSetNarrowBandCache<LevelSetDirectStorageManager>::insert(long id, bool sync)
{
    BITPIT_UNUSED(sync);

    (*m_narrowBandFlag)[id] = 1;

    return id;
}

/*!
 * Erase from the kernel the entry with the specified id.
 *
 * This function cannot be used.
 *
 * \param id is id of the cell
 * \param sync controls if the kernel will be synched
 */
void LevelSetNarrowBandCache<LevelSetDirectStorageManager>::erase(long id, bool sync)
{
    BITPIT_UNUSED(sync);

    (*m_narrowBandFlag)[id] = 0;
}

/*!
 * Checks if the kernel contains an entry with the specified id.
 *
 * \param id is id of the entry
 * \result Return true if kernel contains an entry with the specified id,
 * false otherwise.
 */
bool LevelSetNarrowBandCache<LevelSetDirectStorageManager>::contains(long id) const
{
    return ((*m_narrowBandFlag)[id] == 1);
}

/*!
 * Get a kernel iterator for the entry with the specified id.
 *
 * \param id is id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetNarrowBandCache<LevelSetDirectStorageManager>::KernelIterator LevelSetNarrowBandCache<LevelSetDirectStorageManager>::find(long id) const
{
    if ((*m_narrowBandFlag)[id] == 0) {
        return m_nItems;
    }

    return id;
}

/*!
 * Get a kernel iterator for the entry with the specified raw id.
 *
 * \param rawId is raw id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetNarrowBandCache<LevelSetDirectStorageManager>::KernelIterator LevelSetNarrowBandCache<LevelSetDirectStorageManager>::rawFind(std::size_t rawId) const
{
    if ((*m_narrowBandFlag)[rawId] == 0) {
        return m_nItems;
    }

    return rawId;
}

/*!
 * Get the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset value of the specified entry.
 */
double LevelSetNarrowBandCache<LevelSetDirectStorageManager>::getValue(const KernelIterator &itr) const
{
    return (*m_values)[itr];
}

/*!
 * Get a reference to the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result A reference to the levelset value of the specified entry.
 */
double & LevelSetNarrowBandCache<LevelSetDirectStorageManager>::getValue(const KernelIterator &itr)
{
    return (*m_values)[itr];
}

/*!
 * Get the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset gradient of the specified entry.
 */
const std::array<double, 3> & LevelSetNarrowBandCache<LevelSetDirectStorageManager>::getGradient(const KernelIterator &itr) const
{
    return (*m_gradients)[itr];
}

/*!
 * Get a reference to the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result A reference to the levelset gradient of the specified entry.
 */
std::array<double, 3> & LevelSetNarrowBandCache<LevelSetDirectStorageManager>::getGradient(const KernelIterator &itr)
{
    return (*m_gradients)[itr];
}

/*!
 * Clear the kernel of the storage manager.
 */
void LevelSetNarrowBandCache<LevelSetDirectStorageManager>::clearKernel()
{
    LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>::clearKernel();

    std::fill(m_narrowBandFlag->begin(), m_narrowBandFlag->end(), 0);
}

/*!
 * Dump the kernel of the storage manager.
 *
 * \param stream is the output stream
 */
void LevelSetNarrowBandCache<LevelSetDirectStorageManager>::dumpKernel(std::ostream &stream)
{
    LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>::dumpKernel(stream);

    for (char flag : *m_narrowBandFlag) {
        utils::binary::write(stream, flag);
    }
}

/*!
 * Restore the kernel of the storage manager.
 *
 * \param stream is the input stream
 */
void LevelSetNarrowBandCache<LevelSetDirectStorageManager>::restoreKernel(std::istream &stream)
{
    LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>::restore(stream);

    for (char &flag : *m_narrowBandFlag) {
        utils::binary::read(stream, flag);
    }
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetNarrowBandCache<LevelSetDirectStorageManager>::swap(LevelSetNarrowBandCache<LevelSetDirectStorageManager> &other) noexcept
{
    LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>::swap(other);

    std::swap(other.m_narrowBandFlag, m_narrowBandFlag);
}

/*!
 * Create the narrow band cache.
 *
 * @param object is the levelset object for which the ache will be created
 */
std::shared_ptr<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>> LevelSetNarrowBandCacheFactory<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>::create(LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>> *object)
{
    VolumeKernel *mesh = object->getKernel()->getMesh();
    PiercedVector<Cell, long> &cells = mesh->getCells();

    return std::shared_ptr<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>(new LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>(&cells));
}

/*!
 * Create the narrow band cache.
 *
 * @param object is the levelset object for which the ache will be created
 */
std::shared_ptr<LevelSetNarrowBandCache<LevelSetDirectStorageManager>> LevelSetNarrowBandCacheFactory<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>::create(LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetDirectStorageManager>> *object)
{
    const VolumeKernel *mesh = object->getKernel()->getMesh();
    const std::size_t nCells = mesh->getCellCount();

    return std::shared_ptr<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>(new LevelSetNarrowBandCache<LevelSetDirectStorageManager>(nCells));
}

}
