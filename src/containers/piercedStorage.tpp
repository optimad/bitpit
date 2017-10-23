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

#ifndef __BITPIT_CONTAINERS_PIERCED_STORAGE_TPP__
#define __BITPIT_CONTAINERS_PIERCED_STORAGE_TPP__

namespace bitpit {

/**
*   Constructor.
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage()
    : BasePiercedStorage(), m_nFields(1), m_const_kernel(nullptr), m_kernel(nullptr)
{
}

/**
* Constructor.
*
* \param nFields is the number of fields in the storage
* \param kernel is the kernel that will be set
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage(std::size_t nFields)
    : BasePiercedStorage(), m_nFields(nFields), m_const_kernel(nullptr), m_kernel(nullptr)
{
}

/**
* Constructor.
*
* \param nFields is the number of fields in the storage
* \param kernel is the kernel that will be set
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage(std::size_t nFields, const PiercedKernel<id_t> *kernel)
    : BasePiercedStorage(), m_nFields(nFields), m_const_kernel(nullptr), m_kernel(nullptr)
{
    setStaticKernel(kernel);
}

/**
* Constructor.
*
* \param nFields is the number of fields in the storage
* \param kernel is the kernel that will be set
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage(std::size_t nFields, PiercedKernel<id_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
    : BasePiercedStorage(), m_nFields(nFields), m_const_kernel(nullptr), m_kernel(nullptr)
{
    setDynamicKernel(kernel, syncMode);
}

/**
* Constructor.
*
* \param x is another container of the same type (i.e., instantiated with
* the same template parameters) whose content is copied in this container
* \param kernel is the kernel that will be set
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage(const PiercedStorage<value_t, id_t> &x, const PiercedKernel<id_t> *kernel)
    : m_nFields(x.m_nFields), m_fields(x.m_fields),
      m_const_kernel(nullptr), m_kernel(nullptr)
{
    KernelType kernelType = x.getKernelType();
    switch (kernelType) {

    case KERNEL_STATIC:
    {
        if (kernel) {
            setStaticKernel(kernel);
        } else {
            setStaticKernel(x.m_const_kernel);
        }
        break;
    }

    case KERNEL_DYNAMIC:
    {
        throw std::runtime_error("Unable to set a dynamic kernel. The kernel received in input can be only used to set a static kernel");
        break;
    }

    default:
    {
        break;
    }

    }
}

/**
* Constructor.
*
* \param x is another container of the same type (i.e., instantiated with
* the same template parameters) whose content is copied in this container
* \param kernel is the kernel that will be set
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage(const PiercedStorage<value_t, id_t> &x, PiercedKernel<id_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
    : m_nFields(x.m_nFields), m_fields(x.m_fields),
      m_const_kernel(nullptr), m_kernel(nullptr)
{
    KernelType kernelType = x.getKernelType();
    switch (kernelType) {

    case KERNEL_STATIC:
    {
        if (kernel) {
            setStaticKernel(kernel);
        } else {
            setStaticKernel(x.m_const_kernel);
        }
        break;
    }

    case KERNEL_DYNAMIC:
    {
        if (kernel) {
            setDynamicKernel(kernel, syncMode);
        } else {
            setDynamicKernel(x.m_kernel, syncMode);
        }
        break;
    }

    default:
    {
        break;
    }

    }
}

/**
* Copy constructor.
*
* \param x is another container of the same type (i.e., instantiated with
* the same template parameters) whose content is copied in this container
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::PiercedStorage(const PiercedStorage<value_t, id_t> &x)
    : PiercedStorage<value_t, id_t>(x, nullptr)
{
}

/**
* Destructor
*/
template<typename value_t, typename id_t>
PiercedStorage<value_t, id_t>::~PiercedStorage()
{
    detachKernel();
}

/**
* Gets the number of fields in the storage.
*
* \result The number of fields in the storage.
*/
template<typename value_t, typename id_t>
std::size_t PiercedStorage<value_t, id_t>::getFieldCount() const
{
    return m_nFields;
}

/**
* Sets the kernel that will be used by the storage
*
* The storage will NOT be synchronized with the kernel. Every change to the
* kernel can potentially invalidate the link between kernel and storage.
*
* \param kernel is the kernel that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::setStaticKernel(const PiercedKernel<id_t> *kernel)
{
    if (!kernel) {
        throw std::runtime_error("Unable to set the kernel. Provided kernel is not valid.");
    } else if (m_kernel != nullptr) {
        throw std::runtime_error("Unable to set the kernel. The kernel of the storage is already set.");
    }

    m_kernel       = nullptr;
    m_const_kernel = kernel;

    rawResize(m_const_kernel->rawSize());
}

/**
* Sets the kernel that will be used by the storage.
*
* The storage will dynamically synchronized with the kernel.
*
* \param kernel is the kernel that will be set
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::setDynamicKernel(PiercedKernel<id_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
{
    // Set the static kernel
    setStaticKernel(kernel);

    // Register the storage for dynamic synchronization
    m_kernel = kernel;
    m_kernel->registerSlave(this, syncMode);
}

/**
* Unsets the kernel.
*
* \param release if it's true the memory hold by the container will
* be released, otherwise the container will be cleared but its
* memory will not be released
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::unsetKernel(bool release)
{
    // Detach the kernel
    detachKernel();

    // Clear the storage
    rawClear(release);
}

/**
* Detach the kernel without clearing the container.
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::detachKernel()
{
    if (!m_const_kernel) {
        return;
    }

    if (m_kernel) {
        m_kernel->unregisterSlave(this);
        m_kernel = nullptr;
    }

    m_const_kernel = nullptr;
}

/**
* Gets a constant reference to the kernel of the storage
*
* \result A constant reference to the kernel of the storage.
*/
template<typename value_t, typename id_t>
const PiercedKernel<id_t> * PiercedStorage<value_t, id_t>::getKernel() const
{
    return m_const_kernel;
}

/**
* Get the type of kernel set for the storage.
*
* \return The type of kernel set for the storage.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::KernelType PiercedStorage<value_t, id_t>::getKernelType() const
{
    if (m_kernel) {
        return KERNEL_DYNAMIC;
    } else if (m_const_kernel) {
        return KERNEL_STATIC;
    } else {
        return KERNEL_NONE;
    }
}

/**
* Gets the syncronization mode of the storage.
*
* \result The synchronization mode of the storage.
*/
template<typename value_t, typename id_t>
PiercedSyncMaster::SyncMode PiercedStorage<value_t, id_t>::getSyncMode() const
{
    if (getKernelType() == KERNEL_NONE || getKernelType() == KERNEL_STATIC) {
        return PiercedKernel<id_t>::SYNC_MODE_DISABLED;
    } else {
        return m_const_kernel->getSlaveSyncMode(this);
    }
}

/**
* Returns the number of raw positions in the storage.
*
* This is the number of raw positions in the storage, which is not necessarily
* equal the number of actual objects held in the storage nor to its capacity.
* If the storage is not synchronized with its kernel, the raw size of the
* storage may differ from the raw size of the kernel.
*
* \result The number of raw positions in the storage.
*/
template<typename value_t, typename id_t>
std::size_t PiercedStorage<value_t, id_t>::rawSize() const
{
    return (m_fields.size() / m_nFields);
}

/**
* Commit the specified synchronization action.
*
* Is the action requires new elements to be inserted in the storage, empty
* elements will be inserted.
*
* \param action is the synchronization action that will be commited
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::commitSyncAction(const PiercedSyncAction &action)
{
    switch (action.type) {

    case PiercedSyncAction::TYPE_CLEAR:
    {
        rawClear(true);
        break;
    }

    case PiercedSyncAction::TYPE_RESERVE:
    {
        rawReserve(action.info[PiercedSyncAction::INFO_SIZE]);
        break;
    }

    case PiercedSyncAction::TYPE_RESIZE:
    {
        rawResize(action.info[PiercedSyncAction::INFO_SIZE]);
        break;
    }

    case PiercedSyncAction::TYPE_SHRINK_TO_FIT:
    {
        rawShrinkToFit();
        break;
    }

    case PiercedSyncAction::TYPE_REORDER:
    {
        rawReorder(*action.data);
        break;
    }

    case PiercedSyncAction::TYPE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        rawEmplaceBack();
        break;
    }

    case PiercedSyncAction::TYPE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        rawEmplace(action.info[PiercedSyncAction::INFO_POS]);
        break;
    }

    case PiercedSyncAction::TYPE_OVERWRITE:
    {
        break;
    }

    case PiercedSyncAction::TYPE_MOVE_APPEND:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        rawPushBack(std::move(rawAt(action.info[PiercedSyncAction::INFO_POS_FIRST])));
        rawEmreplace(action.info[PiercedSyncAction::INFO_POS_FIRST]);
        break;
    }

    case PiercedSyncAction::TYPE_MOVE_INSERT:
    {
        // Since we are increasing the sotrage by an element at the time
        // calling a reserve will hurt performance badly because this will
        // prevent the automatic reallocation of the storage.
        rawInsert(action.info[PiercedSyncAction::INFO_POS_SECOND], 1, std::move(rawAt(action.info[PiercedSyncAction::INFO_POS_FIRST])));
        rawEmreplace(action.info[PiercedSyncAction::INFO_POS_FIRST]);
        break;
    }

    case PiercedSyncAction::TYPE_MOVE_OVERWRITE:
    {
        rawSet(action.info[PiercedSyncAction::INFO_POS_SECOND], std::move(rawAt(action.info[PiercedSyncAction::INFO_POS_FIRST])));
        rawEmreplace(action.info[PiercedSyncAction::INFO_POS_FIRST]);
        break;
    }

    case PiercedSyncAction::TYPE_PIERCE:
    {
        // Nothing to do. To improve performance the element will not be
        // cleared.
        break;
    }

    case PiercedSyncAction::TYPE_SWAP:
    {
        rawSwap(action.info[PiercedSyncAction::INFO_POS_FIRST], action.info[PiercedSyncAction::INFO_POS_SECOND]);
        break;
    }

    default:
    {
        throw std::runtime_error("Undefined synchronization action.");
        break;
    }

    }
}

/**
* Requests that the storage capacity be at least enough to contain n elements.
*
* \param n is the minimum capacity requested for the vector, expressed
* in number of elements
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawReserve(std::size_t n)
{
    m_fields.reserve(m_nFields * n);
}

/**
* Requests the storage to reduce its capacity to fit its size.
*
* The request is non-binding, and the function can leave the storage with a
* capacity greater than its size.
*
* This may cause a reallocation, but has no effect on the container size and
* cannot alter its elements.
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawShrinkToFit()
{
    m_fields.shrink_to_fit();
}

/**
* Clears the contents of the storage.
*
* \param release if it's true the memory hold by the container will
* be released, otherwise the container will be cleared but its
* memory will not be relased
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawClear(bool release)
{
    if (release) {
        std::vector<value_t>().swap(m_fields);
    } else {
        m_fields.clear();
    }
}

/**
* Erase the specified number of elements from the storage
*
* \param pos is the position of the first element that will be deleted
* \param n is the number of elements that will be deleted
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawErase(std::size_t pos, std::size_t n)
{
    auto itr_begin = m_fields.begin() + pos * m_nFields;
    auto itr_end   = itr_begin + n * m_nFields;

    m_fields.erase(itr_begin, itr_end);
}

/**
* Swaps two elements.
*
* \param pos_first is the position of the first element that will be swapped
* \param pos_second is the position of the second element that will be swapped
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawSwap(std::size_t pos_first, std::size_t pos_second)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        std::swap(m_fields[pos_first * m_nFields + k], m_fields[pos_second * m_nFields + k]);
    }
}

/**
* Reorder the storage according to the specified permutations.
*
* \param permutations are the permutations that wil be applied
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawReorder(const std::vector<std::size_t> &permutations)
{
    std::size_t storageRawSize = rawSize();
    assert(permutations.size() == storageRawSize);

    // Evaluate the permutation of the fields
    std::vector<std::size_t> fieldPermutations(storageRawSize * m_nFields);
    for (std::size_t pos = 0; pos < storageRawSize; ++pos) {
        std::size_t posOffset         = pos * m_nFields;
        std::size_t permutationOffset = permutations[pos] * m_nFields;
        for (std::size_t k = 0; k < m_nFields; ++k) {
            fieldPermutations[posOffset + k] = permutationOffset + k;
        }
    }

    // Sort the fields
    utils::reorderVector<value_t>(fieldPermutations, m_fields, storageRawSize * m_nFields);
}

/**
* Resizes the container so that it contains n elements.
*
* \param size is the new container size, expressed in number of elements
* \param value is the value to be copied (or moved) to the newly created
* elements
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawResize(std::size_t n, const value_t &value)
{
    m_fields.resize(m_nFields * n, value);
}

/**
* Initialize an element.
*
* \param pos is the position of the element to initialize
* \param args are the arguments forwarded to initialize the new element
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedStorage<value_t>::template has_initialize<Args...>()>::type *>
void PiercedStorage<value_t, id_t>::rawInitialize(std::size_t pos, Args&&... args)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        rawInitialize(pos, k, std::forward<Args>(args)...);
    }
}

/**
* Initialize the specified field of an element.
*
* \param pos is the position of the element to initialize
* \param k is the index of the field to initialize
* \param args are the arguments forwarded to initialize the new element
*/
template<typename value_t, typename id_t>
template<typename... Args, typename std::enable_if<PiercedStorage<value_t>::template has_initialize<Args...>()>::type *>
void PiercedStorage<value_t, id_t>::rawInitialize(std::size_t pos, std::size_t k, Args&&... args)
{
    rawAt(pos, k).initialize(std::forward<Args>(args)...);
}

/**
* Insert the specified number of elements in the storage
*
* \param pos is the position where the new elements will be inserted
* \param n is the number of new elements that will be inserted
* \param value is the value to be copied (or moved) to the inserted elements
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawInsert(std::size_t pos, std::size_t n, const value_t &value)
{
    m_fields.insert(m_fields.begin() + pos * m_nFields, n * m_nFields, value);
}

/**
* Adds a new element at the end of the storage, after its current last element.
*
* \param value is the value to be copied (or moved) to the inserted elements
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawPushBack(const value_t &value)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        m_fields.push_back(value);
    }
}

/**
* Create a new element in the storage. This new element is constructed in
* place using args as the arguments for its construction.
*
* \param pos is the position where the new element will be inserted
* \param args are the arguments forwarded to construct the new element
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<!std::is_same<T, bool>::value>::type *, typename... Args>
void PiercedStorage<value_t, id_t>::rawEmplace(std::size_t pos, Args&&... args)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        m_fields.emplace(m_fields.begin() + pos * m_nFields + k, std::forward<Args>(args)...);
    }
}

/**
* Create a new element in the storage. This new element is constructed in
* place using args as the arguments for its construction.
*
* \param pos is the position where the new element will be inserted
* \param value is the value assigned to the new element
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<std::is_same<T, bool>::value>::type *>
void PiercedStorage<value_t, id_t>::rawEmplace(std::size_t pos, bool value)
{
    rawInsert(pos, 1, value);
}

/**
* Adds a new element at the end of the storage, after its current last element.
* This new element is constructed in place using args as the arguments for its
* construction.
*
* \param args are the arguments forwarded to construct the new element
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<!std::is_same<T, bool>::value>::type *, typename... Args>
void PiercedStorage<value_t, id_t>::rawEmplaceBack(Args&&... args)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        m_fields.emplace_back(std::forward<Args>(args)...);
    }
}

/**
* Adds a new element at the end of the storage, after its current last element.
* This new element is constructed in place using args as the arguments for its
* construction.
*
* \param value is the value assigned to the new element
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<std::is_same<T, bool>::value>::type *>
void PiercedStorage<value_t, id_t>::rawEmplaceBack(bool value)
{
    rawPushBack(value);
}

/**
* Replaces the element at the specified position with a new element. This
* new element is constructed in place using args as the arguments for its
* construction.
*
* \param pos is the position where the new element will be inserted
* \param args are the arguments forwarded to construct the new element
*/
template<typename value_t, typename id_t>
template<typename... Args>
void PiercedStorage<value_t, id_t>::rawEmreplace(std::size_t pos, Args&&... args)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        m_fields[pos * m_nFields + k] = value_t(std::forward<Args>(args)...);
    }
}

/**
* Exchanges the content of the storage by the content of x, which is another
* storage object of the same type. Sizes may differ but the number of fields
* has to be the same.
*
* After the call to this member function, the elements in this storage are
* those which were in x before the call, and the elements of x are those
* which were in this. All iterators, references and pointers remain valid
* for the swapped objects.
*
* \param x is another storage of the same type (i.e., instantiated with the
* same template parameters) whose content is swapped with that of this
* storage.
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::swap(PiercedStorage &x) noexcept
{
    // It is only possible to swap two storages with the same number of field.
    // If this condition is not fulfilled we can not continue. However, we
    // cannot throw an exception because the function is declared snoexcept.
    if (x.getFieldCount() != getFieldCount()) {
        std::cout << "It is only possible to swap storages with the same number of fields." << std::endl;
        assert(false);
        exit(EXIT_FAILURE);
    }

    PiercedSyncSlave::swap(x);
    std::swap(x.m_nFields, m_nFields);
    std::swap(x.m_fields, m_fields);
    std::swap(x.m_kernel, m_kernel);
    std::swap(x.m_const_kernel, m_const_kernel);
}

/**
* Assigns the given value to all the elements in the storage.
*
* \param value is the value to be assigned
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::fill(const value_t &value)
{
    std::fill(m_fields.begin(), m_fields.end(), value);
}

/**
* Gets a constant pointer to the memory array used internally by the vector
* to store its owned elements.
*
* \result A constant pointer to the memory array used internally by the
* vector to store its owned elements.
*/
template<typename value_t, typename id_t>
__PS_CONST_POINTER__ PiercedStorage<value_t, id_t>::data() const
{
    return m_fields.data();
}

/**
* Gets a pointer to the memory array used internally by the vector to store
* its owned elements.
*
* \result A pointer to the memory array used internally by the vector to
* store its owned elements.
*/
template<typename value_t, typename id_t>
__PS_POINTER__ PiercedStorage<value_t, id_t>::data()
{
    return m_fields.data();
}

/**
* Gets a pointer to the data of the specified item.
*
* \param id the id of the item
* \param offset is the offset relative to the first field
* \result A pointer to the data of the specfied item.
*/
template<typename value_t, typename id_t>
__PS_POINTER__ PiercedStorage<value_t, id_t>::data(id_t id, std::size_t offset)
{
    std::size_t pos = m_const_kernel->getPos(id);

    return rawData(pos, offset);
}

/**
* Gets a constant pointer to the data of the specified item.
*
* \param id the id of the item
* \param offset is the offset relative to the first field
* \result A constant pointer to the data of the specfied item.
*/
template<typename value_t, typename id_t>
__PS_CONST_POINTER__ PiercedStorage<value_t, id_t>::data(id_t id, std::size_t offset) const
{
    std::size_t pos = m_const_kernel->getPos(id);

    return rawData(pos, offset);
}

/**
* Gets a constant pointer to the data of of the item at the specified raw
* position.
*
* \param pos is the raw position of the item
* \param offset is the offset relative to the first field
* \result A constant pointer to the data of the item at the specified raw
* position.
*/
template<typename value_t, typename id_t>
__PS_POINTER__ PiercedStorage<value_t, id_t>::rawData(std::size_t pos, std::size_t offset)
{
    return (data() + pos * m_nFields + offset);
}

/**
* Gets a pointer to the data of of the item at the specified raw
* position.
*
* \param pos is the raw position of the item
* \param offset is the offset relative to the first field
* \result A pointer to the data of the item at the specified raw
* position.
*/
template<typename value_t, typename id_t>
__PS_CONST_POINTER__ PiercedStorage<value_t, id_t>::rawData(std::size_t pos, std::size_t offset) const
{
    return (data() + pos * m_nFields + offset);
}

/**
* Returns a reference to the first element of the container.
*
* If the container is empty, an exception is thrown.
*
* \param k is the index of the requested field
* \result A reference to the first element of the container.
*/
template<typename value_t, typename id_t>
__PS_REFERENCE__ PiercedStorage<value_t, id_t>::front(std::size_t k)
{
    if (m_const_kernel->empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(m_const_kernel->front(), k);
}

/**
* Returns a constant reference to the first element of the container.
*
* If the container is empty, an exception is thrown.
*
* \param k is the index of the requested field
* \result A constant reference to the first element of the container.
*/
template<typename value_t, typename id_t>
__PS_CONST_REFERENCE__ PiercedStorage<value_t, id_t>::front(std::size_t k) const
{
    if (m_const_kernel->empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(m_const_kernel->front(), k);
}

/**
* Returns a reference to the last element of the container.
*
* If the container is empty, an exception is thrown.
*
* \param k is the index of the requested field
* \result A reference to the last element of the container.
*/
template<typename value_t, typename id_t>
__PS_REFERENCE__ PiercedStorage<value_t, id_t>::back(std::size_t k)
{
    if (m_const_kernel->empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(m_const_kernel->back(), k);
}

/**
* Returns a constant reference to the last element of the container.
*
* If the container is empty, an exception is thrown.
*
* \param k is the index of the requested field
* \result A constant reference to the last element of the container.
*/
template<typename value_t, typename id_t>
__PS_CONST_REFERENCE__ PiercedStorage<value_t, id_t>::back(std::size_t k) const
{
    if (m_const_kernel->empty()) {
        throw std::out_of_range("Vector is empty");
    }

    return rawAt(m_const_kernel->back(), k);
}


/**
* Gets a reference to the requested field of the specfied item.
*
* \param id is the id of the item
* \param k is the index of the requested field
* \result A reference to the requested field of the specfied item.
*/
template<typename value_t, typename id_t>
__PS_REFERENCE__ PiercedStorage<value_t, id_t>::at(id_t id, std::size_t k)
{
    std::size_t pos = m_const_kernel->getPos(id);

    return rawAt(pos, k);
}

/**
* Gets a constant reference to the requested field of the specfied item.
*
* \param id is the id of the item
* \param k is the index of the requested field
* \result A constant reference to the requested field of the specfied item.
*/
template<typename value_t, typename id_t>
__PS_CONST_REFERENCE__ PiercedStorage<value_t, id_t>::at(id_t id, std::size_t k) const
{
    std::size_t pos = m_const_kernel->getPos(id);

    return rawAt(pos, k);
}

/**
* Gets a reference to the first field of the specfied item.
*
* \param id is the id of the item
* \result A reference to the requested field of the specfied item.
*/
template<typename value_t, typename id_t>
__PS_REFERENCE__ PiercedStorage<value_t, id_t>::operator[](id_t id)
{
    return at(id, 0);
}

/**
* Gets a constant reference to the first field of the specfied item.
*
* \param id is the id of the item
* \result A constant reference to the requested field of the specfied item.
*/
template<typename value_t, typename id_t>
__PS_CONST_REFERENCE__ PiercedStorage<value_t, id_t>::operator[](id_t id) const
{
    return at(id, 0);
}

/**
* Copy all the fields of the specified item.
*
* \param id is the id of the item
* \param values is a pointer to the destination
* \param offset is the offset used for setting the fields
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::copy(id_t id, value_t *values) const
{
    std::size_t pos = m_const_kernel->getPos(id);

    rawCopy(pos, getFieldCount(), 0, values);
}

/**
* Set the requested fields of the specified item.
*
* \param id is the id of the item
* \param nFields is the number of fields that will be copied
* \param offset is the offset used for setting the fields
* \param values is a pointer to the destination
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::copy(id_t id, std::size_t nFields, std::size_t offset, value_t *values) const
{
    std::size_t pos = m_const_kernel->getPos(id);

    rawCopy(pos, nFields, offset, values);
}

/**
* Set all the fields of the specified item.
*
* \param id is the id of the item
* \param value is the value that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::set(id_t id, const value_t &value)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        set(id, k, value);
    }
}

/**
* Set the requested field of the specified item.
*
* \param id is the id of the item
* \param k is the index of the requested field
* \param value is the value that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::set(id_t id, std::size_t k, const value_t &value)
{
    std::size_t pos = m_const_kernel->getPos(id);

    rawSet(pos, k, value);
}

/**
* Set all the fields of the specified item.
*
* \param id is the id of the item
* \param values is a pointer to the values that will be set
* \param offset is the offset used for setting the fields
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::set(id_t id, const value_t *values)
{
    std::size_t pos = m_const_kernel->getPos(id);

    rawSet(pos, getFieldCount(), 0, values);
}

/**
* Set the requested fields of the specified item.
*
* \param id is the id of the item
* \param nFields is the number of fields that will be set
* \param offset is the offset used for setting the fields
* \param values is a pointer to the values that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::set(id_t id, std::size_t nFields, std::size_t offset, const value_t *values)
{
    std::size_t pos = m_const_kernel->getPos(id);

    rawSet(pos, nFields, offset, values);
}

/**
* Gets a reference to the requested field of the item at the specified
* raw position.
*
* \param pos is the raw position of the item
* \param k is the index of the requested field
* \result A reference to the requested field of the item at the specified
* raw position.
*/
template<typename value_t, typename id_t>
__PS_REFERENCE__ PiercedStorage<value_t, id_t>::rawAt(std::size_t pos, std::size_t k)
{
    return m_fields[pos * m_nFields + k];
}

/**
* Gets a constant reference to the requested field of the item at the
* specified raw position.
*
* \param pos is the raw position of the item
* \param k is the index of the requested field
* \result A constant reference to the requested field of the item at the
* specified raw position.
*/
template<typename value_t, typename id_t>
__PS_CONST_REFERENCE__ PiercedStorage<value_t, id_t>::rawAt(std::size_t pos, std::size_t k) const
{
    return m_fields[pos * m_nFields + k];
}

/**
* Copy all the fields of the specified item.
*
* \param id is the id of the item
* \param values is a pointer to the destination
* \param offset is the offset used for setting the fields
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawCopy(std::size_t pos, value_t *values) const
{
    rawCopy(pos, getFieldCount(), 0, values);
}

/**
* Set the requested fields of the specified item.
*
* \param id is the id of the item
* \param nFields is the number of fields that will be copied
* \param offset is the offset used for setting the fields
* \param values is a pointer to the destination
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawCopy(std::size_t pos, std::size_t nFields, std::size_t offset, value_t *values) const
{
    nFields = std::max(nFields, getFieldCount() - offset);
    for (std::size_t k = offset; k < (offset + nFields); ++k) {
        values[k] = m_fields[pos * m_nFields + k];
    }
}

/**
* Set all the fields of the item at the specified raw position.
*
* \param pos is the raw position of the item
* \param value is the value that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawSet(std::size_t pos, const value_t &value)
{
    for (std::size_t k = 0; k < m_nFields; ++k) {
        rawSet(pos, k, value);
    }
}

/**
* Set the requested field of the item at the specified raw position.
*
* \param pos is the raw position of the item
* \param k is the index of the requested field
* \param value is the value that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawSet(std::size_t pos, std::size_t k, const value_t &value)
{
    m_fields[pos * m_nFields + k] = value;
}

/**
* Set all the fields of the item at the specified raw position.
*
* \param pos is the raw position of the item
* \param values is a pointer to the values that will be set
* \param offset is the offset used for setting the fields
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawSet(std::size_t pos, const value_t *values)
{
    rawSet(pos, getFieldCount(), 0, values);
}

/**
* Set the requested fields of the item at the specified raw position.
*
* \param pos is the raw position of the item
* \param nFields is the number of fields that will be set
* \param offset is the offset used for setting the fields
* \param values is a pointer to the values that will be set
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::rawSet(std::size_t pos, std::size_t nFields, std::size_t offset, const value_t *values)
{
    nFields = std::max(nFields, getFieldCount() - offset);
    for (std::size_t k = offset; k < (offset + nFields); ++k) {
        m_fields[pos * m_nFields + k] = values[k];
    }
}

/**
* Gets an iterator pointing to the specified element.
*
* \param id is the id of the specified iterator.
* \result An iterator pointing to the specified element.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::iterator PiercedStorage<value_t, id_t>::find(const id_t &id) noexcept
{
    typename PiercedKernel<id_t>::const_iterator iterator = m_const_kernel->find(id);

    return rawFind(iterator.getPos());
}

/**
* Gets a constant iterator pointing to the specified element.
*
* \param id is the id of the specified iterator.
* \result A constant iterator pointing to the specified element.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::const_iterator PiercedStorage<value_t, id_t>::find(const id_t &id) const noexcept
{
    typename PiercedKernel<id_t>::const_iterator iterator = m_const_kernel->find(id);

    return rawFind(iterator.getPos());
}

/**
* Gets an iterator pointing to the specified position.
*
* \param id is the id of the specified iterator.
* \result An iterator pointing to the specified position.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::iterator PiercedStorage<value_t, id_t>::rawFind(std::size_t pos) noexcept
{
    return iterator(this, pos);
}

/**
* Gets a constant iterator pointing to the specified position.
*
* \param id is the id of the specified iterator.
* \result A constant iterator pointing to the specified position.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::const_iterator PiercedStorage<value_t, id_t>::rawFind(std::size_t pos) const noexcept
{
    return const_iterator(this, pos);
}

/*!
* Returns an iterator pointing to the first element in the vector.
*
* \result An iterator pointing to the first element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::iterator PiercedStorage<value_t, id_t>::begin() noexcept
{
    return rawFind(m_const_kernel->m_begin_pos);
}

/*!
* Returns an iterator referring to the past-the-end element in the vector.
*
* \result An iterator referring to the past-the-end element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::iterator PiercedStorage<value_t, id_t>::end() noexcept
{
    return rawFind(m_const_kernel->m_end_pos);
}

/*!
* Returns a constant iterator pointing to the first element in the vector.
*
* \result A constant iterator pointing to the first element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::const_iterator PiercedStorage<value_t, id_t>::begin() const noexcept
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
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::const_iterator PiercedStorage<value_t, id_t>::end() const noexcept
{
    return cend();
}

/*!
* Returns an conts_iterator pointing to the first element in the vector.
*
* \result A const_iterator pointing to the first element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::const_iterator PiercedStorage<value_t, id_t>::cbegin() const noexcept
{
    return rawFind(m_const_kernel->m_begin_pos);
}

/*!
* Returns an const_iterator referring to the past-the-end element in the
* vector.
*
* \result A const_iterator referring to the past-the-end element in the vector.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::const_iterator PiercedStorage<value_t, id_t>::cend() const noexcept
{
    return rawFind(m_const_kernel->m_end_pos);
}

/*!
* Returns an iterator pointing to the first element in the raw container.
*
* \result An iterator pointing to the first element in the raw container.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::raw_iterator PiercedStorage<value_t, id_t>::rawBegin() noexcept
{
    return m_fields.begin();
}

/*!
* Returns an iterator referring to the past-the-end element in the raw
* container.
*
* \result An iterator referring to the past-the-end element in the raw
* container.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::raw_iterator PiercedStorage<value_t, id_t>::rawEnd() noexcept
{
    return m_fields.end();
}

/*!
* Returns a constant iterator pointing to the first element in the raw
* container.
*
* \result A constant iterator pointing to the first element in the raw
* container.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::raw_const_iterator PiercedStorage<value_t, id_t>::rawBegin() const noexcept
{
    return rawCbegin();
}

/*!
* Returns a constant iterator referring to the past-the-end element in the raw
* container.
*
* \result A constant iterator referring to the past-the-end element in the raw
* container.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::raw_const_iterator PiercedStorage<value_t, id_t>::rawEnd() const noexcept
{
    return rawCend();
}

/*!
* Returns an conts_iterator pointing to the first element in the raw container.
*
* \result A const_iterator pointing to the first element in the raw container.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::raw_const_iterator PiercedStorage<value_t, id_t>::rawCbegin() const noexcept
{
    return m_fields.cbegin();
}

/*!
* Returns an const_iterator referring to the past-the-end element in raw
* container.
*
* \result A const_iterator referring to the past-the-end element in raw
* container.
*/
template<typename value_t, typename id_t>
typename PiercedStorage<value_t, id_t>::raw_const_iterator PiercedStorage<value_t, id_t>::rawCend() const noexcept
{
    return m_fields.cend();
}

/**
* Restore the storage.
*
* \param stream is the stream data should be read from
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<std::is_pod<T>::value || PiercedStorage<T, id_t>::has_restore()>::type *>
void PiercedStorage<value_t, id_t>::restore(std::istream &stream)
{
    // Size
    std::size_t nElements;
    utils::binary::read(stream, nElements);
    rawResize(nElements);

    // Fill the storage
    auto begin_itr = m_fields.cbegin() + m_const_kernel->m_begin_pos;
    auto end_itr   = m_fields.cbegin() + m_const_kernel->m_end_pos;

    size_t pos = m_const_kernel->m_begin_pos;
    for (auto itr = begin_itr; itr != end_itr; ++itr) {
        restoreField(stream, m_fields[pos]);
        ++pos;
    }
}

/**
* Restore a field.
*
* \param stream is the stream data should be read from
* \param value on output will contain the restored value
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::restoreField(std::istream &stream, std::vector<bool>::reference value)
{
    bool bool_value;
    utils::binary::read(stream, bool_value);
    value = bool_value;
}

/**
* Restore a field.
*
* \param stream is the stream data should be read from
* \param value on output will contain the restored value
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<std::is_pod<T>::value>::type *>
void PiercedStorage<value_t, id_t>::restoreField(std::istream &stream, T &value)
{
    utils::binary::read(stream, value);
}


/**
* Restore a field.
*
* \param stream is the stream data should be read from
* \param object on output will contain the restored object
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<PiercedStorage<T, id_t>::has_restore()>::type *>
void PiercedStorage<value_t, id_t>::restoreField(std::istream &stream, T &object)
{
    object.restore(stream);
}

/**
* Dump the storage.
*
* \param stream is the stream data should be written to
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<std::is_pod<T>::value || PiercedStorage<T, id_t>::has_dump()>::type *>
void PiercedStorage<value_t, id_t>::dump(std::ostream &stream) const
{
    // Size
    utils::binary::write(stream, rawSize());

    // Fileds
    auto begin_itr = m_fields.cbegin() + m_const_kernel->m_begin_pos;
    auto end_itr   = m_fields.cbegin() + m_const_kernel->m_end_pos;

    size_t pos = m_const_kernel->m_begin_pos;
    for (auto itr = begin_itr; itr != end_itr; ++itr) {
        dumpField(stream, m_fields[pos]);
        ++pos;
    }
}

/**
* Dump a field.
*
* \param stream is the stream data should be written to
* \param value is the value that will be dumped
*/
template<typename value_t, typename id_t>
void PiercedStorage<value_t, id_t>::dumpField(std::ostream &stream, const std::vector<bool>::const_reference &value) const
{
    bool bool_value = value;
    utils::binary::write(stream, bool_value);
}

/**
* Dump a field.
*
* \param stream is the stream data should be written to
* \param value is the value that will be dumped
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<std::is_pod<T>::value>::type *>
void PiercedStorage<value_t, id_t>::dumpField(std::ostream &stream, const T &value) const
{
    utils::binary::write(stream, value);
}

/**
* Dump a field.
*
* \param stream is the stream data should be written to
* \param object is the object that will be dumped
*/
template<typename value_t, typename id_t>
template<typename T, typename std::enable_if<PiercedStorage<T, id_t>::has_dump()>::type *>
void PiercedStorage<value_t, id_t>::dumpField(std::ostream &stream, const T &object) const
{
    object.dump(stream);
}

}

#endif
