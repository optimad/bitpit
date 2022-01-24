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

#include "levelSetStorage.hpp"

namespace bitpit {

/*!
 * \ingroup levelset
 * \interface LevelSetContainerWrapper
 * \brief Is the base class for defining container wrappers.
 */

/*!
 * Constructor.
 *
 * \param container is the container
 */
LevelSetContainerWrapper::LevelSetContainerWrapper(void *container)
    : m_container(container)
{
}

/*!
 * Get a pointer to the container.
 *
 * \result A pointer to the container.
 */
void * LevelSetContainerWrapper::getContainer()
{
    return m_container;
}

/*!
 * Get a constant pointer to the container.
 *
 * \result A constant pointer to the container.
 */
const void * LevelSetContainerWrapper::getContainer() const
{
    return m_container;
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetContainerWrapper::swap(LevelSetContainerWrapper &other) noexcept
{
    std::swap(other.m_container, m_container);
}

/*!
 * \ingroup levelset
 * \interface LevelSetExternalPiercedStorageManager
 * \brief Is the template class for defining levelset storages managers
 * whose kernel is an external PiercedKernel (where external means that
 * the kernel is defined outside the manager).
 */

/*!
 * Constructor.
 */
LevelSetExternalPiercedStorageManager::LevelSetExternalPiercedStorageManager()
    : LevelSetExternalPiercedStorageManager(nullptr)
{
}

/*!
 * Constructor.
 *
 * \param kernel is the container associated with the storage manager
 */
LevelSetExternalPiercedStorageManager::LevelSetExternalPiercedStorageManager(Kernel *kernel)
    : LevelSetStorageManager<PiercedKernel<long>>(kernel)
{
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
LevelSetExternalPiercedStorageManager::KernelIterator LevelSetExternalPiercedStorageManager::insert(long id, bool sync)
{
    BITPIT_UNUSED(id);
    BITPIT_UNUSED(sync);

    return find(id);
}

/*!
 * Erase from the kernel the entry with the specified id.
 *
 * This function cannot be used.
 *
 * \param id is id of the cell
 * \param sync controls if the kernel will be synched
 */
void LevelSetExternalPiercedStorageManager::erase(long id, bool sync)
{
    BITPIT_UNUSED(id);
    BITPIT_UNUSED(sync);

    // Nothing to do
}

/*!
 * Checks if the kernel contains an entry with the specified id.
 *
 * \param id is id of the entry
 * \result Return true if kernel contains an entry with the specified id,
 * false otherwise.
 */
bool LevelSetExternalPiercedStorageManager::contains(long id) const
{
    return m_kernel->contains(id);
}

/*!
 * Get a kernel iterator for the entry with the specified id.
 *
 * \param id is id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetExternalPiercedStorageManager::KernelIterator LevelSetExternalPiercedStorageManager::find(long id) const
{
    return m_kernel->find(id);
}

/*!
 * Get a kernel iterator for the entry with the specified raw id.
 *
 * \param rawId is raw id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetExternalPiercedStorageManager::KernelIterator LevelSetExternalPiercedStorageManager::rawFind(std::size_t rawId) const
{
    return m_kernel->rawFind(rawId);
}

/*!
 * Get a kernel iterator pointing to the first entry.
 *
 * \result A kernel iterator pointing to the first entry.
 */
LevelSetExternalPiercedStorageManager::KernelIterator LevelSetExternalPiercedStorageManager::begin() const
{
    return m_kernel->cbegin();
}

/*!
 * Get a kernel iterator referring to the past-the-end entry.
 *
 * \result A kernel iterator referring to the past-the-end entry.
 */
LevelSetExternalPiercedStorageManager::KernelIterator LevelSetExternalPiercedStorageManager::end() const
{
    return m_kernel->cend();
}

/*!
 * Synchronize the storages with the kernel.
 *
 * This function cannot be used.
 */
void LevelSetExternalPiercedStorageManager::syncStorages()
{
    // Nothing to do
}

/*!
 * Clear the kernel of the storage manager.
 */
void LevelSetExternalPiercedStorageManager::clearKernel()
{
    // Nothing to do
}

/*!
 * Dump the kernel of the storage manager.
 *
 * \param stream is the output stream
 */
void LevelSetExternalPiercedStorageManager::dumpKernel(std::ostream &stream)
{
    BITPIT_UNUSED(stream);

    // Nothing to do
}

/*!
 * Restore the kernel of the storage manager.
 *
 * \param stream is the input stream
 */
void LevelSetExternalPiercedStorageManager::restoreKernel(std::istream &stream)
{
    BITPIT_UNUSED(stream);

    // Nothing to do
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetExternalPiercedStorageManager::swap(LevelSetExternalPiercedStorageManager &other) noexcept
{
    LevelSetStorageManager<PiercedKernel<long>>::swap(other);
}

/*!
 * \ingroup levelset
 * \interface LevelSetInternalPiercedStorageManager
 * \brief Is the template class for defining levelset storages managers
 * whose kernel is an internal PiercedKernel (where internal means that
 * the kernel is owned by the manager).
 */

/*!
 * Constructor.
 */
LevelSetInternalPiercedStorageManager::LevelSetInternalPiercedStorageManager()
    : LevelSetExternalPiercedStorageManager(&m_internalKernel)
{
}

/*!
 * Insert a kernel entry with the specified id.
 *
 * \param id is id of the entry
 * \param sync controls if the storages will be synched
 * \result A kernel iterator pointing to the newly created entry.
 */
LevelSetInternalPiercedStorageManager::KernelIterator LevelSetInternalPiercedStorageManager::insert(long id, bool sync)
{
    Kernel::FillAction action = m_internalKernel.fillHead(id);
    if (sync) {
        syncStorages();
    }

    return m_internalKernel.rawFind(action.info[PiercedSyncAction::INFO_POS]);
}

/*!
 * Erase from the kernel the entry with the specified id.
 *
 * \param id is id of the cell
 * \param sync controls if the kernel will be synched
 */
void LevelSetInternalPiercedStorageManager::erase(long id, bool sync)
{
    m_internalKernel.erase(id, sync);
    if (sync) {
        syncStorages();
    }
}

/*!
 * Synchronize the storages with the kernel.
 */
void LevelSetInternalPiercedStorageManager::syncStorages()
{
    m_internalKernel.flush();
    m_internalKernel.sync();
}

/*!
 * Clear the kernel of the storage manager.
 */
void LevelSetInternalPiercedStorageManager::clearKernel()
{
    m_internalKernel.clear();
    m_internalKernel.flush();
    m_internalKernel.sync();
}

/*!
 * Dump the kernel of the storage manager.
 *
 * \param stream is the output stream
 */
void LevelSetInternalPiercedStorageManager::dumpKernel(std::ostream &stream)
{
    m_internalKernel.dump(stream);
}

/*!
 * Restore the kernel of the storage manager.
 *
 * \param stream is the input stream
 */
void LevelSetInternalPiercedStorageManager::restoreKernel(std::istream &stream)
{
    m_internalKernel.restore(stream);
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetInternalPiercedStorageManager::swap(LevelSetInternalPiercedStorageManager &other) noexcept
{
    LevelSetExternalPiercedStorageManager::swap(other);

    m_internalKernel.swap(other.m_internalKernel);
}

/*!
 * \ingroup levelset
 * \interface LevelSetDirectStorageManager
 * \brief Is the template class for defining levelset storages managers that
 * don't need a kernel.
 *
 * The kernel of this storage manager is a dummy kernel and contains the number
 * of items that each storage will contain.
 */

/*!
 * Constructor.
 */
LevelSetDirectStorageManager::LevelSetDirectStorageManager()
    : LevelSetDirectStorageManager(0)
{
}

/*!
 * Constructor.
 *
 * \param nItems is the number of items each storage will contain
 */
LevelSetDirectStorageManager::LevelSetDirectStorageManager(std::size_t nItems)
    : LevelSetStorageManager<std::size_t, std::size_t>(nullptr), m_nItems(nItems)
{
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
LevelSetDirectStorageManager::KernelIterator LevelSetDirectStorageManager::insert(long id, bool sync)
{
    BITPIT_UNUSED(sync);

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
void LevelSetDirectStorageManager::erase(long id, bool sync)
{
    BITPIT_UNUSED(id);
    BITPIT_UNUSED(sync);

    // Nothing to do
}

/*!
 * Checks if the kernel contains an entry with the specified id.
 *
 * \param id is id of the entry
 * \result Return true if kernel contains an entry with the specified id,
 * false otherwise.
 */
bool LevelSetDirectStorageManager::contains(long id) const
{
    BITPIT_UNUSED(id);

    return true;
}

/*!
 * Get a kernel iterator for the entry with the specified id.
 *
 * \param id is id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetDirectStorageManager::KernelIterator LevelSetDirectStorageManager::find(long id) const
{
    return id;
}

/*!
 * Get a kernel iterator for the entry with the specified raw id.
 *
 * \param rawId is raw id of the entry
 * \result A kernel iterator pointing to the specified entry.
 */
LevelSetDirectStorageManager::KernelIterator LevelSetDirectStorageManager::rawFind(std::size_t rawId) const
{
    return rawId;
}

/*!
 * Get a kernel iterator pointing to the first entry.
 *
 * \result A kernel iterator pointing to the first entry.
 */
LevelSetDirectStorageManager::KernelIterator LevelSetDirectStorageManager::begin() const
{
    return 0;
}

/*!
 * Get a kernel iterator referring to the past-the-end entry.
 *
 * \result A kernel iterator referring to the past-the-end entry.
 */
LevelSetDirectStorageManager::KernelIterator LevelSetDirectStorageManager::end() const
{
    return m_nItems;
}

/*!
 * Synchronize the storages with the kernel.
 *
 * This function cannot be used.
 */
void LevelSetDirectStorageManager::syncStorages()
{
    // Nothing to do
}

/*!
 * Clear the kernel of the storage manager.
 */
void LevelSetDirectStorageManager::clearKernel()
{
    // Nothing to do
}

/*!
 * Dump the kernel of the storage manager.
 *
 * \param stream is the output stream
 */
void LevelSetDirectStorageManager::dumpKernel(std::ostream &stream)
{
    utils::binary::write(stream, m_nItems);
}

/*!
 * Restore the kernel of the storage manager.
 *
 * \param stream is the input stream
 */
void LevelSetDirectStorageManager::restoreKernel(std::istream &stream)
{
    utils::binary::read(stream, m_nItems);
}

/*!
 * Exchanges the content of the storage manager with the content the specified
 * other storage manager.
 *
 * \param other is another storage manager whose content is swapped with that
 * of this storage manager
 */
void LevelSetDirectStorageManager::swap(LevelSetDirectStorageManager &other) noexcept
{
    LevelSetStorageManager<std::size_t, std::size_t>::swap(other);

    std::swap(other.m_nItems, m_nItems);
}

}
