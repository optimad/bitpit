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

#ifndef __BITPIT_LEVELSET_STORAGE_HPP__
#define __BITPIT_LEVELSET_STORAGE_HPP__

#include <memory>

#include "bitpit_containers.hpp"
#if BITPIT_ENABLE_MPI
#include "bitpit_communications.hpp"
#endif

namespace bitpit {

class LevelSetContainerWrapper
{

public:
    virtual ~LevelSetContainerWrapper() = default;

    void * getContainer();
    const void * getContainer() const;

    void swap(LevelSetContainerWrapper &other) noexcept;

protected:
    void *m_container;

    LevelSetContainerWrapper(void *container);

};

template<typename container_t>
class LevelSetPointerContainerWrapper : public LevelSetContainerWrapper
{

public:
    using LevelSetContainerWrapper::LevelSetContainerWrapper;

    void swap(LevelSetPointerContainerWrapper<container_t> &other) noexcept;

};

template<typename container_t>
class LevelSetSharedContainerWrapper : public LevelSetContainerWrapper
{

public:
    LevelSetSharedContainerWrapper(const std::shared_ptr<container_t> &container);

    void swap(LevelSetSharedContainerWrapper<container_t> &other) noexcept;

protected:
    std::shared_ptr<container_t> m_sharedContainer;

};

template<typename container_t>
class LevelSetUniqueContainerWrapper : public LevelSetContainerWrapper
{

public:
    LevelSetUniqueContainerWrapper(std::unique_ptr<container_t> &&container);

    void swap(LevelSetUniqueContainerWrapper<container_t> &other) noexcept;

protected:
    std::unique_ptr<container_t> m_uniqueContainer;

};

template<typename kernel_iterator_t>
class LevelSetBaseStorage
{

public:
    typedef kernel_iterator_t KernelIterator;

    virtual ~LevelSetBaseStorage() = default;

    void * getContainer();
    const void * getContainer() const;

    virtual void clear() = 0;

    virtual void dump(std::ostream &stream) = 0;
    virtual void restore(std::istream &stream) = 0;

#if BITPIT_ENABLE_MPI
    virtual std::size_t getItemBinarySize() const = 0;
    virtual void writeItem(const KernelIterator &kernelIterator, SendBuffer &dataBufferbuffer) = 0;
    virtual void readItem(const KernelIterator &kernelIterator, RecvBuffer &dataBufferbuffer) = 0;
#endif

    void swap(LevelSetBaseStorage<kernel_iterator_t> &other) noexcept;

protected:
    LevelSetBaseStorage(LevelSetContainerWrapper *containerWrapper);

private:
    LevelSetContainerWrapper *m_containerWrapper;

};

template<typename kernel_iterator_t, typename container_t>
class LevelSetStorage : public LevelSetBaseStorage<kernel_iterator_t>
{

public:
    typedef container_t Container;
    typedef typename LevelSetBaseStorage<kernel_iterator_t>::KernelIterator KernelIterator;

    LevelSetStorage(LevelSetContainerWrapper *containerWrapper);

    void clear() override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

#if BITPIT_ENABLE_MPI
    std::size_t getItemBinarySize() const override;
    void writeItem(const KernelIterator &kernelIterator, SendBuffer &buffer) override;
    void readItem(const KernelIterator &kernelIterator, RecvBuffer &buffer) override;
#endif

    void swap(LevelSetStorage<kernel_iterator_t, container_t> &other) noexcept;

};

template<typename value_t, typename id_t = long>
class LevelSetPiercedStorage : public LevelSetBaseStorage<typename PiercedKernel<id_t>::const_iterator>
{

public:
    typedef PiercedStorage<value_t, id_t> Container;
    typedef typename LevelSetBaseStorage<typename PiercedKernel<id_t>::const_iterator>::KernelIterator KernelIterator;

    LevelSetPiercedStorage(LevelSetContainerWrapper *containerWrapper);

    void clear() override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

#if BITPIT_ENABLE_MPI
    std::size_t getItemBinarySize() const override;
    void writeItem(const KernelIterator &kernelIterator, SendBuffer &buffer) override;
    void readItem(const KernelIterator &kernelIterator, RecvBuffer &buffer) override;
#endif

    void swap(LevelSetPiercedStorage<value_t, id_t> &other) noexcept;

};

template<typename value_t>
class LevelSetDirectStorage : public LevelSetBaseStorage<std::size_t>
{

public:
    typedef std::vector<value_t> Container;
    typedef std::size_t KernelIterator;

    LevelSetDirectStorage(LevelSetContainerWrapper *containerWrapper);

    void clear() override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

#if BITPIT_ENABLE_MPI
    std::size_t getItemBinarySize() const override;
    void writeItem(const KernelIterator &kernelIterator, SendBuffer &buffer) override;
    void readItem(const KernelIterator &kernelIterator, RecvBuffer &buffer) override;
#endif

    void swap(LevelSetDirectStorage<value_t> &other) noexcept;

};

template<typename kernel_t, typename kernel_iterator_t = typename kernel_t::const_iterator>
class LevelSetStorageManager
{

public:
    /*!
    * Defines the kernel synchronization mode.
    *
    * Possible synchronization modes are:
    *  - manual, in this mode the storage manager is in charge of updating the kernel;
    *  - automatic, in this mode the storage manager doesn't need to updated the kernel,
    *    because its update will be performed automatically.
    */
    enum KernelSyncMode {
        KERNEL_SYNC_MODE_MANUAL,
        KERNEL_SYNC_MODE_AUTOMATIC,
    };

    typedef kernel_t Kernel;
    typedef kernel_iterator_t KernelIterator;

    virtual ~LevelSetStorageManager() = default;

    bool isDirty() const;
    void setDirty(bool dirty);

    const Kernel * getKernel() const;
    KernelSyncMode getKernelSyncMode() const;

    virtual KernelIterator insert(long id, bool sync = true) = 0;
    virtual void erase(long id, bool sync = true) = 0;

    virtual bool contains(long id) const = 0;

    virtual KernelIterator find(long id) const = 0;
    virtual KernelIterator rawFind(std::size_t rawId) const = 0;

    virtual KernelIterator begin() const = 0;
    virtual KernelIterator end() const = 0;

    virtual void clear();

    std::size_t getStorageCount() const;

    template<typename container_t>
    container_t * addStorage(int id, container_t *container);

    template<typename container_t>
    container_t * addStorage(int id, const std::shared_ptr<container_t> &container);

    template<typename container_t>
    container_t * addStorage(int id, std::unique_ptr<container_t> &&container);

    void * addStorage(int id, std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>> &&storage);

    virtual void syncStorages() = 0;

    template<typename container_t>
    container_t * getContainer(int id);

    template<typename container_t>
    const container_t * getContainer(int id) const;

    void * getContainer(int id);

    const void * getContainer(int id) const;

    virtual void dump(std::ostream &stream);
    virtual void restore(std::istream &stream);

#if BITPIT_ENABLE_MPI
    virtual void write(const std::vector<long> &ids, SendBuffer &buffer);
    virtual void read(const std::vector<long> &ids, RecvBuffer &buffer);
#endif

    void swap(LevelSetStorageManager<kernel_t, kernel_iterator_t> &other) noexcept;

protected:
    bool m_dirty; //! Controls if the manager is dirty
    Kernel *m_kernel; //! Kernel associated with the manager
    KernelSyncMode m_kernelSyncMode; //! Kernel synchronization mode
    std::unordered_map<int, std::unique_ptr<LevelSetBaseStorage<kernel_iterator_t>>> m_storages; //! Storages handled by the manager
    std::vector<std::unique_ptr<LevelSetContainerWrapper>> m_containers; //! Containers that are handled by the manager

    LevelSetStorageManager(Kernel *kernel, KernelSyncMode kernelSyncMode);

    virtual void clearKernel() = 0;

    virtual void dumpKernel(std::ostream &stream) = 0;
    virtual void restoreKernel(std::istream &stream) = 0;

};

class LevelSetExternalPiercedStorageManager : public LevelSetStorageManager<PiercedKernel<long>>
{

public:
    template<typename T>
    using Storage = PiercedStorage<T, long>;

    using StorageSyncMode = PiercedSyncMaster::SyncMode;

    KernelIterator insert(long id, bool sync = true) override;
    void erase(long id, bool sync = true) override;

    bool contains(long id) const override;

    KernelIterator find(long id) const override;
    KernelIterator rawFind(std::size_t) const override;

    KernelIterator begin() const override;
    KernelIterator end() const override;

    template<typename value_t>
    Storage<value_t> * addStorage(int id, int nFields);

    void syncStorages() override;

    void swap(LevelSetExternalPiercedStorageManager &other) noexcept;

protected:
    StorageSyncMode m_storageSyncMode; //! Storage synchronization mode

    LevelSetExternalPiercedStorageManager();
    LevelSetExternalPiercedStorageManager(Kernel *kernel, KernelSyncMode kernelSyncMode, StorageSyncMode storageSyncMode);

    void clearKernel() override;

    void dumpKernel(std::ostream &stream) override;
    void restoreKernel(std::istream &stream) override;

};

class LevelSetInternalPiercedStorageManager : public LevelSetExternalPiercedStorageManager
{

public:
    KernelIterator insert(long id, bool sync = true) override;
    void erase(long id, bool sync = true) override;

    void syncStorages() override;

    void swap(LevelSetInternalPiercedStorageManager &other) noexcept;

protected:
    Kernel m_internalKernel; //! Internal pierced kernel

    LevelSetInternalPiercedStorageManager();
    LevelSetInternalPiercedStorageManager(StorageSyncMode storageSyncMode);

    void clearKernel() override;

    void dumpKernel(std::ostream &stream) override;
    void restoreKernel(std::istream &stream) override;

};

class LevelSetDirectStorageManager : public LevelSetStorageManager<std::size_t, std::size_t>
{

public:
    template<typename T>
    using Storage = std::vector<T>;

    KernelIterator insert(long id, bool sync = true) override;
    void erase(long id, bool sync = true) override;

    bool contains(long id) const override;

    KernelIterator find(long id) const override;
    KernelIterator rawFind(std::size_t) const override;

    KernelIterator begin() const override;
    KernelIterator end() const override;

    template<typename value_t>
    Storage<value_t> * addStorage(int id);

    void syncStorages() override;

    void swap(LevelSetDirectStorageManager &other) noexcept;

protected:
    std::size_t m_nItems;

    LevelSetDirectStorageManager();
    LevelSetDirectStorageManager(std::size_t nItems);

    void clearKernel() override;

    void dumpKernel(std::ostream &stream) override;
    void restoreKernel(std::istream &stream) override;

};

}

// Include template implemenetations
#include <levelSetStorage.tpp>

#endif 
