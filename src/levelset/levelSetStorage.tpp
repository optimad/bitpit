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

#ifndef __BITPIT_LEVELSET_STORAGE_TPP__
#define __BITPIT_LEVELSET_STORAGE_TPP__

namespace bitpit {

/*!
 * \ingroup levelset
 * \interface LevelSetPointerContainerWrapper
 * \brief Is the base class for defining container wrappers around pointers.
 */

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename container_t>
void LevelSetPointerContainerWrapper<container_t>::swap(LevelSetPointerContainerWrapper<container_t> &other) noexcept
{
    LevelSetContainerWrapper::swap(other);
}

/*!
 * \ingroup levelset
 * \interface LevelSetSharedContainerWrapper
 * \brief Is the base class for defining container wrappers around shared
 * pointers.
 */

/*!
 * Constructor.
 *
 * \param container is the container
 */
template<typename container_t>
LevelSetSharedContainerWrapper<container_t>::LevelSetSharedContainerWrapper(const std::shared_ptr<container_t> &container)
    : LevelSetContainerWrapper(container.get()), m_sharedContainer(container)
{
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename container_t>
void LevelSetSharedContainerWrapper<container_t>::swap(LevelSetSharedContainerWrapper<container_t> &other) noexcept
{
    LevelSetContainerWrapper::swap(other);

    std::swap(other.m_sharedContainer, m_sharedContainer);
}

/*!
 * \ingroup levelset
 * \interface LevelSetUniqueContainerWrapper
 * \brief Is the base class for defining container wrappers around shared
 * pointers.
 */

/*!
 * Constructor.
 *
 * \param container is the container
 */
template<typename container_t>
LevelSetUniqueContainerWrapper<container_t>::LevelSetUniqueContainerWrapper(std::unique_ptr<container_t> &&container)
    : LevelSetContainerWrapper(container.get()), m_uniqueContainer(std::move(container))
{
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename container_t>
void LevelSetUniqueContainerWrapper<container_t>::swap(LevelSetUniqueContainerWrapper<container_t> &other) noexcept
{
    LevelSetContainerWrapper::swap(other);

    std::swap(other.m_uniqueContainer, m_uniqueContainer);
}

/*!
 * \ingroup levelset
 * \interface LevelSetBaseStorage
 * \brief Is the base class for defining levelset storages.
 */

/*!
 * Constructor.
 *
 * \param containerWrapper is the container wrapper
 */
template<typename kernel_iterator_t>
LevelSetBaseStorage<kernel_iterator_t>::LevelSetBaseStorage(LevelSetContainerWrapper *containerWrapper)
    : m_containerWrapper(containerWrapper)
{
}

/*!
 * Get a pointer to the container.
 *
 * \result A pointer to the container.
 */
template<typename kernel_iterator_t>
void * LevelSetBaseStorage<kernel_iterator_t>::getContainer()
{
    return m_containerWrapper->getContainer();
}

/*!
 * Get a constant pointer to the container.
 *
 * \result A constant pointer to the container.
 */
template<typename kernel_iterator_t>
const void * LevelSetBaseStorage<kernel_iterator_t>::getContainer() const
{
    return m_containerWrapper->getContainer();
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename kernel_iterator_t>
void LevelSetBaseStorage<kernel_iterator_t>::swap(LevelSetBaseStorage<kernel_iterator_t> &other) noexcept
{
    std::swap(other.m_containerWrapper, m_containerWrapper);
}

/*!
 * \ingroup levelset
 * \interface LevelSetStorage
 * \brief Is the template class for defining levelset storages.
 */

/*!
 * Constructor.
 *
 * \param containerWrapper is the container wrapper
 */
template<typename kernel_iterator_t, typename container_t>
LevelSetStorage<kernel_iterator_t, container_t>::LevelSetStorage(LevelSetContainerWrapper *containerWrapper)
    : LevelSetBaseStorage<kernel_iterator_t>(containerWrapper)
{
}

/*!
 * Clear the storage.
 */
template<typename kernel_iterator_t, typename container_t>
void LevelSetStorage<kernel_iterator_t, container_t>::clear()
{
    Container *container = static_cast<Container *>(this->getContainer());

    container->clear();
}

/*!
 * Dump the storage.
 *
 * \param stream is the output stream
 */
template<typename kernel_iterator_t, typename container_t>
void LevelSetStorage<kernel_iterator_t, container_t>::dump(std::ostream &stream)
{
    const Container *container = static_cast<const Container *>(this->getContainer());

    for (const auto &item : *container) {
        utils::binary::write(stream, item);
    }
}

/*!
 * Restore the storage.
 *
 * \param stream is the input stream
 */
template<typename kernel_iterator_t, typename container_t>
void LevelSetStorage<kernel_iterator_t, container_t>::restore(std::istream &stream)
{
    Container *container = static_cast<Container *>(this->getContainer());

    for (auto &item : *container) {
        utils::binary::read(stream, item);
    }
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename kernel_iterator_t, typename container_t>
void LevelSetStorage<kernel_iterator_t, container_t>::swap(LevelSetStorage<kernel_iterator_t, container_t> &other) noexcept
{
    LevelSetBaseStorage<kernel_iterator_t>::swap(other);
}

#if BITPIT_ENABLE_MPI
/*!
 * Get the size, expressed in bytes, of an item.
 *
 * \result The size, expressed in bytes, of an item.
 */
template<typename kernel_iterator_t, typename container_t>
std::size_t LevelSetStorage<kernel_iterator_t, container_t>::getItemBinarySize() const
{
    return sizeof(typename container_t::value_type);
}

/*!
 * Write narrow band entry into the communication buffer.
 *
 * \param kernelIterator is the kernel iterator pointing to entry associated to
 * the item
 * \param[in,out] buffer buffer for data communication
 */
template<typename kernel_iterator_t, typename container_t>
void LevelSetStorage<kernel_iterator_t, container_t>::writeItem(const KernelIterator &kernelIterator, SendBuffer &buffer)
{
    const Container *container = static_cast<const Container *>(this->getContainer());

    buffer << (*container)[kernelIterator];
}

/*!
 * Read narrow band entry from the communication buffer.
 *
 * \param kernelIterator is the kernel iterator pointing to entry associated to
 * the item
 * \param[in,out] buffer buffer containing data
 */
template<typename kernel_iterator_t, typename container_t>
void LevelSetStorage<kernel_iterator_t, container_t>::readItem(const KernelIterator &kernelIterator, RecvBuffer &buffer)
{
    Container *container = static_cast<Container *>(this->getContainer());

    buffer >> (*container)[kernelIterator];
}
#endif

/*!
 * \ingroup levelset
 * \interface LevelSetPiercedStorage
 * \brief Is the template class for defining levelset storages based on
 * a PiercedStorage class.
 */

/*!
 * Constructor.
 *
 * \param containerWrapper is the container wrapper
 */
template<typename value_t, typename id_t>
LevelSetPiercedStorage<value_t, id_t>::LevelSetPiercedStorage(LevelSetContainerWrapper *containerWrapper)
    : LevelSetBaseStorage<KernelIterator>(containerWrapper)
{
}

/*!
 * Clear the storage.
 */
template<typename value_t, typename id_t>
void LevelSetPiercedStorage<value_t, id_t>::clear()
{
    // Nothing to do
}

/*!
 * Dump the storage.
 *
 * \param stream is the output stream
 */
template<typename value_t, typename id_t>
void LevelSetPiercedStorage<value_t, id_t>::dump(std::ostream &stream)
{
    const Container *container = static_cast<const Container *>(this->getContainer());

    container->dump(stream);
}

/*!
 * Restore the storage.
 *
 * \param stream is the input stream
 */
template<typename value_t, typename id_t>
void LevelSetPiercedStorage<value_t, id_t>::restore(std::istream &stream)
{
    Container *container = static_cast<Container *>(this->getContainer());

    container->restore(stream);
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename value_t, typename id_t>
void LevelSetPiercedStorage<value_t, id_t>::swap(LevelSetPiercedStorage<value_t, id_t> &other) noexcept
{
    LevelSetBaseStorage<KernelIterator>::swap(other);
}

#if BITPIT_ENABLE_MPI
/*!
 * Get the size, expressed in bytes, of an item.
 *
 * \result The size, expressed in bytes, of an item.
 */
template<typename value_t, typename id_t>
std::size_t LevelSetPiercedStorage<value_t, id_t>::getItemBinarySize() const
{
    return sizeof(typename Container::value_type);
}

/*!
 * Write narrow band entry into the communication buffer.
 *
 * \param kernelIterator is the kernel iterator pointing to entry associated to
 * the item
 * \param[in,out] buffer buffer for data communication
 */
template<typename value_t, typename id_t>
void LevelSetPiercedStorage<value_t, id_t>::writeItem(const KernelIterator &kernelIterator, SendBuffer &buffer)
{
    const Container *container = static_cast<const Container *>(this->getContainer());
    std::size_t itemRawId = kernelIterator.getRawIndex();

    buffer << container->rawAt(itemRawId);
}

/*!
 * Read narrow band entry from the communication buffer.
 *
 * \param kernelIterator is the kernel iterator pointing to entry associated to
 * the item
 * \param[in,out] buffer buffer containing data
 */
template<typename value_t, typename id_t>
void LevelSetPiercedStorage<value_t, id_t>::readItem(const KernelIterator &kernelIterator, RecvBuffer &buffer)
{
    Container *container = static_cast<Container *>(this->getContainer());
    std::size_t itemRawId = kernelIterator.getRawIndex();

    buffer >> container->rawAt(itemRawId);
}
#endif

/*!
 * \ingroup levelset
 * \interface LevelSetDirectStorage
 * \brief Is the template class for defining levelset storages that .
 */

/*!
 * Constructor.
 *
 * \param containerWrapper is the container wrapper
 */
template<typename value_t>
LevelSetDirectStorage<value_t>::LevelSetDirectStorage(LevelSetContainerWrapper *containerWrapper)
    : LevelSetBaseStorage<KernelIterator>(containerWrapper)
{
}

/*!
 * Clear the storage.
 */
template<typename value_t>
void LevelSetDirectStorage<value_t>::clear()
{
    Container *container = static_cast<Container *>(this->getContainer());

    container->clear();
}

/*!
 * Dump the storage.
 *
 * \param stream is the output stream
 */
template<typename value_t>
void LevelSetDirectStorage<value_t>::dump(std::ostream &stream)
{
    const Container *container = static_cast<const Container *>(this->getContainer());

    for (const auto &item : *container) {
        utils::binary::write(stream, item);
    }
}

/*!
 * Restore the storage.
 *
 * \param stream is the input stream
 */
template<typename value_t>
void LevelSetDirectStorage<value_t>::restore(std::istream &stream)
{
    Container *container = static_cast<Container *>(this->getContainer());

    for (auto &item : *container) {
        utils::binary::read(stream, item);
    }
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename value_t>
void LevelSetDirectStorage<value_t>::swap(LevelSetDirectStorage<value_t> &other) noexcept
{
    LevelSetBaseStorage<KernelIterator>::swap(other);
}

#if BITPIT_ENABLE_MPI
/*!
 * Get the size, expressed in bytes, of an item.
 *
 * \result The size, expressed in bytes, of an item.
 */
template<typename value_t>
std::size_t LevelSetDirectStorage<value_t>::getItemBinarySize() const
{
    return sizeof(value_t);
}

/*!
 * Write narrow band entry into the communication buffer.
 *
 * \param kernelIterator is the kernel iterator pointing to entry associated to
 * the item
 * \param[in,out] buffer buffer for data communication
 */
template<typename value_t>
void LevelSetDirectStorage<value_t>::writeItem(const KernelIterator &kernelIterator, SendBuffer &buffer)
{
    const Container *container = static_cast<const Container *>(this->getContainer());

    buffer << (*container)[kernelIterator];
}

/*!
 * Read narrow band entry from the communication buffer.
 *
 * \param kernelIterator is the kernel iterator pointing to entry associated to
 * the item
 * \param[in,out] buffer buffer containing data
 */
template<typename value_t>
void LevelSetDirectStorage<value_t>::readItem(const KernelIterator &kernelIterator, RecvBuffer &buffer)
{
    Container *container = static_cast<Container *>(this->getContainer());

    buffer >> (*container)[kernelIterator];
}
#endif

/*!
 * \ingroup levelset
 * \interface LevelSetStorageManager
 * \brief Is the template class for defining levelset storages managers.
 */

/*!
 * Constructor.
 *
 * \param kernel is the container associated with the storage manager
 * \param kernelSyncMode is the synchronization mode of the kernel
 */
template<typename kernel_t, typename kernel_iterator_t>
LevelSetStorageManager<kernel_t, kernel_iterator_t>::LevelSetStorageManager(Kernel *kernel, KernelSyncMode kernelSyncMode)
    : m_dirty(true), m_kernel(kernel), m_kernelSyncMode(kernelSyncMode)
{
}

/*!
 * Check if the storage manager is dirty.
 *
 * \result Returns true if the storage manager is dirty, false otherwise.
 */
template<typename kernel_t, typename kernel_iterator_t>
bool LevelSetStorageManager<kernel_t, kernel_iterator_t>::isDirty() const
{
    return m_dirty;
}

/*!
 * Set the storage manager as dirty.
 *
 * \param dirty if set to true the storage manager will be set as dirty
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::setDirty(bool dirty)
{
    m_dirty = dirty;
}

/*!
 * Get a constant pointer to the kernel associated with manager.
 *
 * \result A constant pointer to the kernel associated with manager.
 */
template<typename kernel_t, typename kernel_iterator_t>
const typename LevelSetStorageManager<kernel_t, kernel_iterator_t>::Kernel * LevelSetStorageManager<kernel_t, kernel_iterator_t>::getKernel() const
{
    return m_kernel;
}

/*!
 * Get the synchronization mode of the kernel.
 *
 * \result The synchronization mode of the kernel.
 */
template<typename kernel_t, typename kernel_iterator_t>
typename LevelSetStorageManager<kernel_t, kernel_iterator_t>::KernelSyncMode LevelSetStorageManager<kernel_t, kernel_iterator_t>::getKernelSyncMode() const
{
    return m_kernelSyncMode;
}

/*!
 * Clear the storage manager.
 *
 * Storages will not be erased, only their contents will be cleared.
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::clear()
{
    // Clear storages
    for (const auto &entry : m_storages) {
        entry.second->clear();
    }

    // Clear kernel
    clearKernel();
}

/*!
 * Count the number of storages in the manager.
 *
 * \result The number of storages in the manager.
 */
template<typename kernel_t, typename kernel_iterator_t>
std::size_t LevelSetStorageManager<kernel_t, kernel_iterator_t>::getStorageCount() const
{
    return m_storages.size();
}

/*!
 * Add a new storage.
 *
 * If a storage with the same id is already in the manager, an exception will
 * be thrown.
 *
 * \param id is the id that will be assigned to the storage
 * \param container is the container associated with the storage
 * \result Returns a pointer to the container.
 */
template<typename kernel_t, typename kernel_iterator_t>
template<typename container_t>
container_t * LevelSetStorageManager<kernel_t, kernel_iterator_t>::addStorage(int id, container_t *container)
{
    m_containers.emplace_back(new LevelSetPointerContainerWrapper<container_t>(container));
    LevelSetContainerWrapper *containerWrapper = m_containers.back().get();

    std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>> storage = std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>>(new LevelSetStorage<kernel_iterator_t, container_t>(containerWrapper));

    return static_cast<container_t *>(addStorage(id, std::move(storage)));
}

/*!
 * Add a new storage.
 *
 * If a storage with the same id is already in the manager, an exception will
 * be thrown.
 *
 * \param id is the id that will be assigned to the storage
 * \param container is the container associated with the storage
 * \result Returns a pointer to the container.
 */
template<typename kernel_t, typename kernel_iterator_t>
template<typename container_t>
container_t * LevelSetStorageManager<kernel_t, kernel_iterator_t>::addStorage(int id, const std::shared_ptr<container_t> &container)
{
    m_containers.emplace_back(new LevelSetSharedContainerWrapper<container_t>(container));
    LevelSetContainerWrapper *containerWrapper = m_containers.back().get();

    std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>> storage = std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>>(new LevelSetStorage<kernel_iterator_t, container_t>(containerWrapper));

    return static_cast<container_t *>(addStorage(id, std::move(storage)));
}

/*!
 * Add a new storage.
 *
 * If a storage with the same id is already in the manager, an exception will
 * be thrown.
 *
 * \param id is the id that will be assigned to the storage
 * \param container is the container associated with the storage
 * \result Returns a pointer to the container.
 */
template<typename kernel_t, typename kernel_iterator_t>
template<typename container_t>
container_t * LevelSetStorageManager<kernel_t, kernel_iterator_t>::addStorage(int id, std::unique_ptr<container_t> &&container)
{
    m_containers.emplace_back(new LevelSetUniqueContainerWrapper<container_t>(std::move(container())));
    LevelSetContainerWrapper *containerWrapper = m_containers.back().get();

    std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>> storage = std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>>(new LevelSetStorage<kernel_iterator_t, container_t>(containerWrapper));

    return static_cast<container_t *>(addStorage(id, std::move(storage)));
}

/*!
 * Add a new storage.
 *
 * If a storage with the same id is already in the manager, an exception will
 * be thrown.
 *
 * \param id is the id that will be assigned to the storage
 * \param storage is the storage
 * \result Returns a pointer to the container.
 */
template<typename kernel_t, typename kernel_iterator_t>
void * LevelSetStorageManager<kernel_t, kernel_iterator_t>::addStorage(int id, std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>> &&storage)
{
    if (m_storages.count(id) != 0) {
        throw std::runtime_error("The manager already contains a storage with the specified id.");
    }

    m_storages.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple(std::move(storage)));

    return getContainer(id);
}

/*!
 * Get a constant pointer to the container associated with the specified
 * storage.
 *
 * \param id is the id of the storage
 * \result A constant pointer to the container associated with the specified
 * storage.
 */
template<typename kernel_t, typename kernel_iterator_t>
template<typename container_t>
container_t * LevelSetStorageManager<kernel_t, kernel_iterator_t>::getContainer(int id)
{
    return static_cast<container_t *>(getContainer(id));
}

/*!
 * Get a constant pointer to the container associated with the specified
 * storage.
 *
 * \param id is the id of the storage
 * \result A constant pointer to the container associated with the specified
 * storage.
 */
template<typename kernel_t, typename kernel_iterator_t>
template<typename container_t>
const container_t * LevelSetStorageManager<kernel_t, kernel_iterator_t>::getContainer(int id) const
{
    return static_cast<const container_t *>(getContainer(id));
}

/*!
 * Get a constant pointer to the container associated with the specified
 * storage.
 *
 * \param id is the id of the storage
 * \result A constant pointer to the container associated with the specified
 * storage.
 */
template<typename kernel_t, typename kernel_iterator_t>
void * LevelSetStorageManager<kernel_t, kernel_iterator_t>::getContainer(int id)
{
    return m_storages.at(id)->getContainer();
}

/*!
 * Get a constant pointer to the container associated with the specified
 * storage.
 *
 * \param id is the id of the storage
 * \result A constant pointer to the container associated with the specified
 * storage.
 */
template<typename kernel_t, typename kernel_iterator_t>
const void * LevelSetStorageManager<kernel_t, kernel_iterator_t>::getContainer(int id) const
{
    return m_storages.at(id)->getContainer();
}

/*!
 * Dump the storage manager.
 *
 * \param stream is the output stream
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::dump(std::ostream &stream)
{
    // Dump manager properties
    utils::binary::write(stream, m_dirty);

    // Dump kernel
    dumpKernel(stream);

    // Dump storages
    for (const auto &entry : m_storages) {
        entry.second->dump(stream);
    }
}

/*!
 * Restore the storage manager.
 *
 * \param stream is the input stream
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::restore(std::istream &stream)
{
    // Restore manager properties
    utils::binary::read(stream, m_dirty);

    // Restore kernel
    restoreKernel(stream);

    // Restore storages
    for (auto &entry : m_storages) {
        entry.second->restore(stream);
    }
}

#if BITPIT_ENABLE_MPI
/*!
 * Write storage manager data into the buffers.
 *
 * The function will receive a list of ids that will be sent and will write
 * into the buffer the items whose id are contained in the specified list.
 *
 * \param ids list of ids that will be send
 * \param[in,out] buffer buffer for data communication
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::write(const std::vector<long> &ids, SendBuffer &buffer)
{
    // Evaluate the size of the buffer
    std::size_t nBufferItems = 0;
    for (long id : ids) {
        if (contains(id)) {
            ++nBufferItems;
        }
    }

    std::size_t bufferSize = buffer.getSize() + sizeof(std::size_t);
    for (const auto &storage : m_storages) {
        bufferSize += nBufferItems * storage.second->getItemBinarySize();
    }
    buffer.setSize(bufferSize);

    // Fill the buffer
    buffer << nBufferItems;
    for (std::size_t k = 0; k < ids.size(); ++k) {
        // Get an iterator pointing to the item
        long id = ids[k];
        KernelIterator itr = find(id);
        if (itr == end()) {
            continue;
        }

        // Write the index of the cell
        buffer << k;

        // Write item data
        for (const auto &storage : m_storages) {
            storage.second->writeItem(itr, buffer);
        }
    }
}

/*!
 * Read storage data from the buffer.
 *
 * The function will receive a list of ids that has been received and will read
 * from the buffer the items of the sotrage.
 *
 * \param ids list of ids that has been received
 * \param[in,out] buffer buffer containing the data
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::read(const std::vector<long> &ids, RecvBuffer &buffer)
{
    std::size_t nReceviedItems;
    buffer >> nReceviedItems;

    for (std::size_t i = 0; i < nReceviedItems; ++i) {
        // Read the id of the item
        std::size_t k;
        buffer >> k;
        long id = ids[k];

        // Create the item
        KernelIterator itr = find(id);
        if (itr == end()) {
            itr = insert(id);
        }

        // Read item data
        for (auto &storage : m_storages) {
            storage.second->readItem(itr, buffer);
        }
    }
}
#endif

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
template<typename kernel_t, typename kernel_iterator_t>
void LevelSetStorageManager<kernel_t, kernel_iterator_t>::swap(LevelSetStorageManager<kernel_t, kernel_iterator_t> &other) noexcept
{
    std::swap(other.m_dirty, m_dirty);
    std::swap(other.m_kernel, m_kernel);
    std::swap(other.m_storages, m_storages);
    std::swap(other.m_containers, m_containers);
}

/*!
 * Add a new storage.
 *
 * If a storage with the same id is already in the manager, an exception will
 * be thrown.
 *
 * \param id is the id that will be assigned to the storage
 * \param nFields are the number of fields
 * \param syncMode is the synchronization mode of the storage
 * \result Returns a pointer to the container.
 */
template<typename value_t>
LevelSetExternalPiercedStorageManager::Storage<value_t> * LevelSetExternalPiercedStorageManager::addStorage(int id, int nFields)
{
    std::unique_ptr<Storage<value_t>> container = std::unique_ptr<Storage<value_t>>(new Storage<value_t>(nFields, m_kernel, m_storageSyncMode));
    m_containers.emplace_back(new LevelSetUniqueContainerWrapper<Storage<value_t>>(std::move(container)));
    LevelSetContainerWrapper *containerWrapper = m_containers.back().get();

    std::unique_ptr<LevelSetBaseStorage<KernelIterator>> storage = std::unique_ptr<LevelSetBaseStorage<KernelIterator>>(new LevelSetPiercedStorage<value_t>(containerWrapper));

    return static_cast<Storage<value_t> *>(LevelSetStorageManager<Kernel, KernelIterator>::addStorage(id, std::move(storage)));
}

/*!
 * Add a new storage.
 *
 * If a storage with the same id is already in the manager, an exception will
 * be thrown.
 *
 * \param id is the id that will be assigned to the storage
 * \result Returns a pointer to the container.
 */
template<typename value_t>
LevelSetDirectStorageManager::Storage<value_t> * LevelSetDirectStorageManager::addStorage(int id)
{
    std::unique_ptr<Storage<value_t>> container = std::unique_ptr<Storage<value_t>>(new Storage<value_t>(m_nItems));
    m_containers.emplace_back(new LevelSetUniqueContainerWrapper<Storage<value_t>>(std::move(container)));
    LevelSetContainerWrapper *containerWrapper = m_containers.back().get();

    std::unique_ptr<LevelSetBaseStorage<KernelIterator>> storage = std::unique_ptr<LevelSetBaseStorage<KernelIterator>>(new LevelSetDirectStorage<value_t>(containerWrapper));

    return static_cast<Storage<value_t> *>(LevelSetStorageManager<Kernel, KernelIterator>::addStorage(id, std::move(storage)));
}

}

#endif 
