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

#ifndef __BITPIT_PROXY_VECTOR_TPP__
#define __BITPIT_PROXY_VECTOR_TPP__

namespace bitpit {

/*!
    Constructor
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t>::ProxyVectorIterator()
    : m_position(nullptr)
{
}

/*!
    Constructor

    This constructor allows to generate a constant iterator from a non
    constnat iterator.

    \param other is the iterator that will be copied
*/
template<typename value_t, typename container_t>
template<typename other_value_t, typename std::enable_if<std::is_const<value_t>::value && !std::is_const<other_value_t>::value && std::is_same<other_value_t, typename std::remove_cv<value_t>::type>::value, int>::type>
ProxyVectorIterator<value_t, container_t>::ProxyVectorIterator(const ProxyVectorIterator<other_value_t, container_t> &other)
    : m_position(other.m_position)
{
}

/*!
    Exchanges the values of the current iterator and
    the iterator recevied as argument.

    \param other is the iterator to exchange values with
*/
template<typename value_t, typename container_t>
void ProxyVectorIterator<value_t, container_t>::swap(ProxyVectorIterator& other) noexcept
{
    std::swap(m_position, other.m_position);
}

/*!
    Pre-increment operator.
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t> & ProxyVectorIterator<value_t, container_t>::operator++()
{
    m_position++;

    return *this;
}

/*!
    Post-increment operator.
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t> ProxyVectorIterator<value_t, container_t>::operator++(int)
{
    ProxyVectorIterator tmp(m_position);

    ++(*this);

    return tmp;
}

/*!
    Pre-decrement operator.
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t> & ProxyVectorIterator<value_t, container_t>::operator--()
{
    m_position--;

    return *this;
}

/*!
    Post-decrement operator.
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t> ProxyVectorIterator<value_t, container_t>::operator--(int)
{
    ProxyVectorIterator tmp(m_position);

    --(*this);

    return tmp;
}

/*!
    Compound assigment operator.

    \param increment is the increment
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t>& ProxyVectorIterator<value_t, container_t>::operator+=(int increment)
{
    m_position += increment;

    return *this;
}

/*!
    Deference operator.

    \result A reference to the element currently pointed to by the iterator.
*/
template<typename value_t, typename container_t>
__PXI_REFERENCE__ ProxyVectorIterator<value_t, container_t>::operator*() const
{
    return *m_position;
}

/*!
    Deference operator.

    \result A reference to the element currently pointed to by the iterator.
*/
template<typename value_t, typename container_t>
__PXI_POINTER__ ProxyVectorIterator<value_t, container_t>::operator->() const
{
    return m_position;
}

/*!
* Copy assignment operator to create a constant iterator from a non-constant
* one.
*
* \param other is the iterator that will be copied
*/
template<typename value_t, typename container_t>
template<typename other_value_t, typename std::enable_if<std::is_const<value_t>::value && !std::is_const<other_value_t>::value && std::is_same<other_value_t, typename std::remove_cv<value_t>::type>::value, int>::type>
ProxyVectorIterator<value_t, container_t> & ProxyVectorIterator<value_t, container_t>::operator=(const ProxyVectorIterator<other_value_t, container_t> &other)
{
    m_position = other.m_position;

    return *this;
}

/*!
    Creates a new iterator and initializes it with the position of
    the const base iterator recevied in input.
*/
template<typename value_t, typename container_t>
ProxyVectorIterator<value_t, container_t>::ProxyVectorIterator(__PXI_POINTER__ position)
    : m_position(position)
{
}

/*!
    Distance operator.

    \param other is the iterator from which the distance will be evaluated
    \result The distance between the specified iterator.
*/
template<typename value_t, typename container_t>
std::size_t ProxyVectorIterator<value_t, container_t>::operator-(const ProxyVectorIterator &other)
{
    return (m_position - other.m_position);
}

/*!
    Constructor.

    \param size is the size of the storage expressed in number of elements.
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::ProxyVectorDummyStorage(std::size_t size)
{
    BITPIT_UNUSED(size);
}

/*!
    Swaps the contents.

    \param other is another storage of the same type
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
void ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::swap(ProxyVectorDummyStorage &other) noexcept
{
    BITPIT_UNUSED(other);
}

/*!
    Return a pointer to the data.

    \result A pointer to the data.
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
typename ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::pointer ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::data()
{
    return nullptr;
}

/*!
    Return a constant pointer to the data.

    \result A constant pointer to the data.
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
typename ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::const_pointer ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::data() const
{
    return nullptr;
}

/*!
    Check if the storage is empty.

    \result Returns true if the storage is empty, false otherwise.
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
bool ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::empty() const
{
    return true;
}

/*!
    Get the size of the storage expressed in number of elements.

    \result The size of the storage, expressed in number of elements.
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
std::size_t ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::size() const
{
    return 0;
}

/*!
    Resize the storage.

    \result The size of the storage expressed in number of elements.
*/
template<typename value_t, typename pointer_t, typename const_pointer_t>
void ProxyVectorDummyStorage<value_t, pointer_t, const_pointer_t>::resize(std::size_t size)
{
    BITPIT_UNUSED(size);
}

/*!
    Memory pool
*/
template<typename value_t, typename container_t, bool thread_safe>
std::vector<std::unique_ptr<container_t>> ProxyVectorStorage<value_t, container_t, thread_safe>::m_containerPool = std::vector<std::unique_ptr<container_t>>();

/*!
    Create a data container.

    \param size is the size of the container, expressed in number of elements.
    \param allowEmpty controls if the container will be created also when the
    requested size is zero
*/
template<typename value_t, typename container_t, bool thread_safe>
std::unique_ptr<container_t> ProxyVectorStorage<value_t, container_t, thread_safe>::createContainer(std::size_t size, bool allowEmpty)
{
    if (size == 0 && !allowEmpty) {
        return std::unique_ptr<container_t>(nullptr);
    }

    if (!thread_safe) {
        if (!m_containerPool.empty()) {
            std::unique_ptr<container_t> container = std::move(m_containerPool.back());
            if (container->size() < size) {
                container->resize(size);
            }

            m_containerPool.resize(m_containerPool.size() - 1);

            return container;
        } else {
            return std::unique_ptr<container_t>(new container_t(size));
        }
    } else {
        return std::unique_ptr<container_t>(new container_t(size));
    }
}

/*!
    Create a data container.

    \param source is the container whose content will be copied in newly
    created container
    \param allowEmpty controls if the container will be created also when the
    source container is empty
*/
template<typename value_t, typename container_t, bool thread_safe>
std::unique_ptr<container_t> ProxyVectorStorage<value_t, container_t, thread_safe>::createContainer(const std::unique_ptr<container_t> &source, bool allowEmpty)
{
    if (!source || source->empty()) {
        if (allowEmpty) {
            return createContainer(0, true);
        } else {
            return std::unique_ptr<container_t>(nullptr);
        }
    }

    if (!thread_safe) {
        if (!m_containerPool.empty()) {
            std::unique_ptr<container_t> container = createContainer(source->size(), false);
            std::copy_n(source->data(), source->size(), container->data());

            return container;
        } else {
            return std::unique_ptr<container_t>(new container_t(*source));
        }
    } else {
        return std::unique_ptr<container_t>(new container_t(*source));
    }
}

/*!
    Delete a data container.
*/
template<typename value_t, typename container_t, bool thread_safe>
void ProxyVectorStorage<value_t, container_t, thread_safe>::destroyContainer(std::unique_ptr<container_t> *container)
{
    if (!(*container)) {
        return;
    }

    if (!thread_safe) {
        if (m_containerPool.size() < MEMORY_POOL_VECTOR_COUNT) {
            if ((*container)->size() > MEMORY_POOL_MAX_CAPACITY) {
                (*container)->resize(MEMORY_POOL_MAX_CAPACITY);
            }

            m_containerPool.emplace_back(std::move(*container));
            container->reset();
        }
    }
}

/*!
    Constructor.

    \param size is the size of the storage expressed in number of elements.
*/
template<typename value_t, typename container_t, bool thread_safe>
ProxyVectorStorage<value_t, container_t, thread_safe>::ProxyVectorStorage(std::size_t size)
    : m_container(createContainer(size, false))
{
}

/*!
    Copy constructor.

    \param x is another storage of the same type (i.e., instantiated with
    the same template parameters) whose content is copied in this container.
*/
template<typename value_t, typename container_t, bool thread_safe>
ProxyVectorStorage<value_t, container_t, thread_safe>::ProxyVectorStorage(const ProxyVectorStorage<value_t, container_t, thread_safe> &other)
    : m_container(createContainer(other.m_container, false))
{
}

/*!
    Copy assigment operator.

    \param x is another storage of the same type (i.e., instantiated with
    the same template parameters) whose content is copied in this container.
*/
template<typename value_t, typename container_t, bool thread_safe>
ProxyVectorStorage<value_t, container_t, thread_safe> & ProxyVectorStorage<value_t, container_t, thread_safe>::operator=(const ProxyVectorStorage<value_t, container_t, thread_safe> &other)
{
    ProxyVectorStorage<value_t, container_t, thread_safe> temporary(other);
    this->swap(temporary);

    return *this;
}

/*!
    Destructor.
*/
template<typename value_t, typename container_t, bool thread_safe>
ProxyVectorStorage<value_t, container_t, thread_safe>::~ProxyVectorStorage()
{
    destroyContainer(&m_container);
}

/*!
    Get a reference to the data container associated with the storage.

    \param forceCreation if set to true and the storage is not associated with
    a container, an empty container will be created
    \result A reference to the data container associated with the storage.
*/
template<typename value_t, typename container_t, bool thread_safe>
container_t * ProxyVectorStorage<value_t, container_t, thread_safe>::container(bool forceCreation)
{
    if (!m_container && forceCreation) {
        m_container = createContainer(0, true);
    }

    return m_container.get();
}

/*!
    Get a reference to the data container associated with the storage.

    \param forceCreation if set to true and the storage is not associated with
    a container, an empty container will be created
    \result A reference to the data container associated with the storage.
*/
template<typename value_t, typename container_t, bool thread_safe>
const container_t * ProxyVectorStorage<value_t, container_t, thread_safe>::container(bool forceCreation) const
{
    if (!m_container && forceCreation) {
        m_container = createContainer(0, true);
    }

    return m_container.get();
}

/*!
    Swaps the contents.

    \param other is another storage of the same type
*/
template<typename value_t, typename container_t, bool thread_safe>
void ProxyVectorStorage<value_t, container_t, thread_safe>::swap(ProxyVectorStorage &other) noexcept
{
    m_container.swap(other.m_container);
}

/*!
    Return a pointer to the data.

    \result A pointer to the data.
*/
template<typename value_t, typename container_t, bool thread_safe>
typename ProxyVectorStorage<value_t, container_t, thread_safe>::pointer ProxyVectorStorage<value_t, container_t, thread_safe>::data()
{
    if (empty()) {
        return nullptr;
    }

    return m_container->data();
}

/*!
    Return a constant pointer to the data.

    \result A constant pointer to the data.
*/
template<typename value_t, typename container_t, bool thread_safe>
typename ProxyVectorStorage<value_t, container_t, thread_safe>::const_pointer ProxyVectorStorage<value_t, container_t, thread_safe>::data() const
{
    if (empty()) {
        return nullptr;
    }

    return m_container->data();
}

/*!
    Check if the storage is empty.

    \result Returns true if the storage is empty, false otherwise.
*/
template<typename value_t, typename container_t, bool thread_safe>
bool ProxyVectorStorage<value_t, container_t, thread_safe>::empty() const
{
    if (!m_container) {
        return true;
    }

    return m_container->empty();
}

/*!
    Get the size of the storage expressed in number of elements.

    \result The size of the storage, expressed in number of elements.
*/
template<typename value_t, typename container_t, bool thread_safe>
std::size_t ProxyVectorStorage<value_t, container_t, thread_safe>::size() const
{
    if (!m_container) {
        return 0;
    }

    return m_container->size();
}

/*!
    Resize the storage.

    \result The size of the storage expressed in number of elements.
*/
template<typename value_t, typename container_t, bool thread_safe>
void ProxyVectorStorage<value_t, container_t, thread_safe>::resize(std::size_t size)
{
    if (size == 0) {
        destroyContainer(&m_container);
        return;

    }

    if (m_container) {
        m_container->resize(size);
    } else {
        m_container = createContainer(size, false);
    }
}

/*!
    Constructor
*/
template<typename value_t, bool thread_safe>
ProxyVector<value_t, thread_safe>::ProxyVector()
    : m_size(0), m_data(nullptr)
{
}

/*!
    Constructor

    The container will create an internal storage that later can be filled
    with data. This is allowed because the container points to constant data,
    i.e., the container is not allowed to change the data it points to. Having
    the data stored internallt or pointing to external data doesn't change the
    behaviour of the container: in both cases it acts as a proxy to some
    constant data.

    \param size is the number elements contained in the data
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
ProxyVector<value_t, thread_safe>::ProxyVector(std::size_t size)
    : ProxyVector<value_t, thread_safe>(INTERNAL_STORAGE, size, size)
{
}

/*!
    Constructor

    The container will create an internal storage that later can be filled
    with data. This is allowed because the container points to constant data,
    i.e., the container is not allowed to change the data it points to. Having
    the data stored internallt or pointing to external data doesn't change the
    behaviour of the container: in both cases it acts as a proxy to some
    constant data.

    \param size is the number elements contained in the data
    \param capacity is the size of the internal storage space expressed in
    number of elements, the capacity of the container cannot be smaller than
    the size of the data, if a smaller capacity is specified the storage will
    be resized using data size
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
ProxyVector<value_t, thread_safe>::ProxyVector(std::size_t size, std::size_t capacity)
    : ProxyVector<value_t, thread_safe>::ProxyVector(INTERNAL_STORAGE, size, capacity)
{
}

/*!
    Constructor

    The container will create an internal storage that later can be filled
    with data. This is allowed because the container points to constant data,
    i.e., the container is not allowed to change the data it points to. Having
    the data stored internallt or pointing to external data doesn't change the
    behaviour of the container: in both cases it acts as a proxy to some
    constant data.

    \param data a pointer to the data
    \param size is the number elements contained in the data
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
ProxyVector<value_t, thread_safe>::ProxyVector(__PXV_POINTER__ data, std::size_t size)
    : ProxyVector<value_t, thread_safe>::ProxyVector(data, size, (data != INTERNAL_STORAGE) ? 0 : size)
{
}

/*!
    Constructor

    If data is set to INTERNAL_STORAGE, the container will create an internal
    storage that later can be filled with data. This is allowed because the
    container points to constant data, i.e., the container is not allowed to
    change the data it points to. Having the data stored internallt or pointing
    to external data doesn't change the behaviour of the container: in both
    cases it acts as a proxy to some constant data.

    \param data a pointer to the data
    \param size is the number elements contained in the data
    \param capacity is the size of the internal storage space expressed in
    number of elements, the capacity of the container cannot be smaller than
    the size of the data, if a smaller capacity is specified the storage will
    be resized using data size
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
ProxyVector<value_t, thread_safe>::ProxyVector(__PXV_POINTER__ data, std::size_t size, std::size_t capacity)
    : m_storage(capacity), m_size(size), m_data((data != INTERNAL_STORAGE) ? data : m_storage.data())
{
}

/*!
    Constructor

    Containers that point to non-constant data cannot use the internal storage.
    This guarantees that all the pointers returned by the container are always
    pointing to the original data (i.e., the container acts as a proxy to the
    original data).

    \param data a pointer to the data
    \param size is the number elements contained in the data
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
ProxyVector<value_t, thread_safe>::ProxyVector(__PXV_POINTER__ data, std::size_t size)
    : m_storage(0), m_size(size), m_data(data)
{
    assert(data != INTERNAL_STORAGE);
}

/*!
    Copy constructor.

    \param other is another container whose content is copied in this container
*/
template<typename value_t, bool thread_safe>
ProxyVector<value_t, thread_safe>::ProxyVector(const ProxyVector &other)
    : m_storage(other.m_storage), m_size(other.m_size), m_data(other.storedData() ? m_storage.data() : other.m_data)
{
}

/*!
    Move constructor.

    We need to explicitly implement the move constructor to workaround a bug
    in gcc, see:

        https://gcc.gnu.org/bugzilla/show_bug.cgi?id=60796
        https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57728

    \param other is another container whose content is moved in this container
*/
template<typename value_t, bool thread_safe>
ProxyVector<value_t, thread_safe>::ProxyVector(ProxyVector &&other)
    : m_storage(std::move(other.m_storage)), m_size(std::move(other.m_size)), m_data(std::move(other.m_data))
{
}

/*!
    Copy assignment operator.

    Assigns new contents to the container, replacing its current contents,
    and modifying its size accordingly.

    \param other is another container whose content is copied in this container
*/
template<typename value_t, bool thread_safe>
ProxyVector<value_t, thread_safe> & ProxyVector<value_t, thread_safe>::operator=(const ProxyVector &other)
{
    if (this != &other) {
        ProxyVector temporary(other);
        temporary.swap(*this);
    }

    return *this;
}

/*!
    Move assignment operator.

    Assigns new contents to the container, replacing its current contents,
    and modifying its size accordingly.

    We need to explicitly implement the move assignment operator to workaround
    a bug in gcc, see:

        https://gcc.gnu.org/bugzilla/show_bug.cgi?id=60796
        https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57728

    \param other is another container whose content is moved in this container
*/
template<typename value_t, bool thread_safe>
ProxyVector<value_t, thread_safe> & ProxyVector<value_t, thread_safe>::operator=(ProxyVector &&other)
{
    if (this != &other) {
        ProxyVector temporary(std::move(other));
        temporary.swap(*this);
    }

    return *this;
}

/*!
    Sets the content of the container.

    If data is set to INTERNAL_STORAGE, the container will create an internal
    storage that later can be filled with data. This is allowed because the
    container points to constant data, i.e., the container is not allowed to
    change the data it points to. Having the data stored internallt or pointing
    to external data doesn't change the behaviour of the container: in both
    cases it acts as a proxy to some constant data.

    \param data a pointer to the data, if the value INTERNAL_STORAGE is
    specified, the proxy will point to the data contained in the internal
    storage
    \param size is the number elements contained in the data
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
void ProxyVector<value_t, thread_safe>::set(__PXV_POINTER__ data, std::size_t size)
{
    std::size_t capacity;
    if (data != INTERNAL_STORAGE) {
        capacity = 0;
    } else {
        capacity = size;
    }
    set(data, size, capacity);
}

/*!
    Sets the content of the container.

    If data is set to INTERNAL_STORAGE, the container will create an internal
    storage that later can be filled with data. This is allowed because the
    container points to constant data, i.e., the container is not allowed to
    change the data it points to. Having the data stored internallt or pointing
    to external data doesn't change the behaviour of the container: in both
    cases it acts as a proxy to some constant data.

    \param data a pointer to the data, if the value INTERNAL_STORAGE is
    specified, the proxy will point to the data contained in the internal
    storage
    \param size is the number elements contained in the data
    \param capacity is the number elements the internal contained should be
    able to contain, the capacity of the container cannot be smaller thatn
    the size of the data, if a smaller capacity is specified the storage
    will be resize using data size
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
void ProxyVector<value_t, thread_safe>::set(__PXV_POINTER__ data, std::size_t size, std::size_t capacity)
{
    m_storage.resize(std::max(size, capacity));

    m_size = size;
    if (data == INTERNAL_STORAGE) {
        if (m_size == 0) {
            m_data = nullptr;
        } else {
            m_data = m_storage.data();
        }
    } else {
        m_data = data;
    }
}

/*!
    Sets the content of the container.

    Containers that point to non-constant data cannot use the internal storage.
    This guarantees that all the pointers returned by the container are always
    pointing to the original data (i.e., the container acts as a proxy to the
    original data).

    \param data a pointer to the data, containers that point to non-constant
    data cannot use the internal storage, hence the value INTERNAL_STORAGE in
    not allowed
    \param size is the number elements contained in the data
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
void ProxyVector<value_t, thread_safe>::set(__PXV_POINTER__ data, std::size_t size)
{
    assert(data != INTERNAL_STORAGE);

    m_storage.resize(0);

    m_size = size;
    m_data = data;
}

/*!
    Returns a direct pointer to the memory of the internal storage.

    If the container is not using the internal storage, a null pointer is
    returned.

    \result A a direct pointer to the memory of the internal storage.
*/
template<typename value_t, bool thread_safe>
__PXV_STORAGE_POINTER__ ProxyVector<value_t, thread_safe>::storedData() noexcept
{
    __PXV_STORAGE_POINTER__ internalData = m_storage.data();
    if (!internalData) {
        return nullptr;
    } else if (internalData != m_data) {
        return nullptr;
    }

    return internalData;
}

/*!
    Returns a constant direct pointer to the memory of the internal storage.

    If the container is not using the internal storage, a null pointer is
    returned.

    \result A a direct pointer to the memory of the internal storage.
*/
template<typename value_t, bool thread_safe>
__PXV_STORAGE_CONST_POINTER__ ProxyVector<value_t, thread_safe>::storedData() const noexcept
{
    __PXV_STORAGE_CONST_POINTER__ internalData = m_storage.data();
    if (!internalData) {
        return nullptr;
    } else if (internalData != m_data) {
        return nullptr;
    }

    return internalData;
}

/*!
    Returns a direct reference to the container associated with the internal
    storage.

    Interacting directly with the container associated with the internal
    storage may leave the proxy vector in an inconsistent state. It's up
    to the caller of this function to guarantee that this will not happen.

    \param forceCreation if set to true and the storage is not associated with
    a container, an empty container will be created
    \result A direct reference to the container associated with the internal
    storage.
*/
template<typename value_t, bool thread_safe>
template<typename U, typename std::enable_if<std::is_const<U>::value, int>::type>
typename ProxyVector<value_t, thread_safe>::container_type * ProxyVector<value_t, thread_safe>::storedDataContainer(bool forceCreation)
{
    return m_storage.container(forceCreation);
}

/*!
    Returns a constant direct reference to the container associated with the
    internal storage.

    Interacting directly with the container associated with the internal
    storage may leave the proxy vector in an inconsistent state. It's up
    to the caller of this function to guarantee that this will not happen.

    \param forceCreation if set to true and the storage is not associated with
    a container, an empty container will be created
    \result A constant direct reference to the container associated with the
    internal storage.
*/
template<typename value_t, bool thread_safe>
template<typename U, typename std::enable_if<std::is_const<U>::value, int>::type>
const typename ProxyVector<value_t, thread_safe>::container_type * ProxyVector<value_t, thread_safe>::storedDataContainer(bool forceCreation) const
{
    return m_storage.container(forceCreation);
}

/*!
    Swaps the content.

    \param other is another container of the same type
*/
template<typename value_t, bool thread_safe>
void ProxyVector<value_t, thread_safe>::swap(ProxyVector &other)
{
    std::swap(m_size, other.m_size);
    std::swap(m_data, other.m_data);

    m_storage.swap(other.m_storage);
}

/*!
    Tests whether two containers are equal.

    \result true if the containers are equal, false otherwise.
*/
template<typename value_t, bool thread_safe>
bool ProxyVector<value_t, thread_safe>::operator==(const ProxyVector& other) const
{
    if (m_size != other.m_size) {
        return false;
    }

    if (!storedData() && !other.storedData()) {
        if (m_data != other.m_data) {
            return false;
        }
    }

    for (std::size_t i = 0; i < m_size; ++i) {
        if (m_data[i] != other.m_data[i]) {
            return false;
        }
    }

    return true;
}

/*!
    Tests whether the container is empty.

    \result true if the container size is 0, false otherwise.
*/
template<typename value_t, bool thread_safe>
bool ProxyVector<value_t, thread_safe>::empty() const
{
    return size() == 0;
}

/*!
    Returns the number of elements in the container

    \result The number of elements in the container.
*/
template<typename value_t, bool thread_safe>
std::size_t ProxyVector<value_t, thread_safe>::size() const
{
    return m_size;
}

/*!
    Returns a direct pointer to the memory where the elments are stored.

    \result A direct pointer to the memory where the elments are stored.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_POINTER__ ProxyVector<value_t, thread_safe>::data() noexcept
{
    return m_data;
}

/*!
    Returns a direct constant pointer to the memory where the elments are
    stored.

    \result A direct constant pointer to the memory where the elments are
    stored.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_POINTER__ ProxyVector<value_t, thread_safe>::data() const noexcept
{
    return m_data;
}

/*!
    Returns a reference to the specified element.

    \param n is the position of the requested element
    \result A reference to the specified element.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t, thread_safe>::operator[](std::size_t n)
{
    return m_data[n];
}

/*!
    Returns a constant reference to the specified element.

    \param n is the position of the requested element
    \result A constant reference to the specified element.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_REFERENCE__ ProxyVector<value_t, thread_safe>::operator[](std::size_t n) const
{
    return m_data[n];
}

/*!
    Returns a reference to the specified element.

    \param n is the position of the requested element
    \result A reference to the specified element.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t, thread_safe>::at(std::size_t n)
{
    return m_data[n];
}

/*!
    Returns a constant reference to the specified element.

    \param n is the position of the requested element
    \result A constant reference to the specified element.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_REFERENCE__ ProxyVector<value_t, thread_safe>::at(std::size_t n) const
{
    return m_data[n];
}

/*!
    Gets a reference to the first element in the container.

    \result A reference to the first element in the container.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t, thread_safe>::front()
{
    return m_data[0];
}

/*!
    Gets a constant reference to the first element in the container.

    \result A constant reference to the first element in the container.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_REFERENCE__ ProxyVector<value_t, thread_safe>::front() const
{
    return m_data[0];
}

/*!
    Gets a reference to the last element in the container.

    \result A reference to the last element in the container.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t, thread_safe>::back()
{
    return m_data[m_size - 1];
}

/*!
    Gets a constant reference to the last element in the container.

    \result A constant reference to the last element in the container.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_REFERENCE__ ProxyVector<value_t, thread_safe>::back() const
{
    return m_data[m_size - 1];
}

/*!
    Returns an iterator pointing to the first element in the container.

    \result An iterator pointing to the first element in the container.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_ITERATOR__ ProxyVector<value_t, thread_safe>::begin()
{
    return iterator(m_data);
}

/*!
    Returns a constant iterator pointing to the first element in the container.

    \result A constant iterator pointing to the first element in the container.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_ITERATOR__ ProxyVector<value_t, thread_safe>::begin() const
{
    return const_iterator(m_data);
}

/*!
    Returns an iterator referring to the past-the-end element in the container.

    \result An iterator referring to the past-the-end element in the container.
*/
template<typename value_t, bool thread_safe>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_ITERATOR__ ProxyVector<value_t, thread_safe>::end()
{
    return iterator(m_data + m_size);
}

/*!
    Returns a constant iterator referring to the past-the-end element in the
    container.

    \result A constant iterator referring to the past-the-end element in the
    container.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_ITERATOR__ ProxyVector<value_t, thread_safe>::end() const
{
    return const_iterator(m_data + m_size);
}

/*!
    Returns a constant iterator pointing to the first element in the container.

    \result A constant iterator pointing to the first element in the container.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_ITERATOR__ ProxyVector<value_t, thread_safe>::cbegin()
{
    return const_iterator(m_data);
}

/*!
    Returns a constant iterator referring to the past-the-end element in the
    container.
git gui

    \result A constant iterator referring to the past-the-end element in the
    container.
*/
template<typename value_t, bool thread_safe>
__PXV_CONST_ITERATOR__ ProxyVector<value_t, thread_safe>::cend()
{
    return const_iterator(m_data + m_size);
}

}

#endif
