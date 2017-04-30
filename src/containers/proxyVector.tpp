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

#ifndef __BITPIT_PROXY_VECTOR_TPP__
#define __BITPIT_PROXY_VECTOR_TPP__

namespace bitpit{

/*!
    Creates a new uninitialized iterator
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv>::ProxyVectorIterator()
    : m_position(nullptr)
{
}

/*!
    Exchanges the values of the current iterator and
    the iterator recevied as argument.

    \param other the iterator to exchange values with
*/
template<typename T, typename T_no_cv>
void ProxyVectorIterator<T, T_no_cv>::swap(ProxyVectorIterator& other) noexcept
{
    std::swap(m_position, other.m_position);
}

/*!
    Pre-increment operator.
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv> & ProxyVectorIterator<T, T_no_cv>::operator++()
{
    m_position++;

    return *this;
}

/*!
    Post-increment operator.
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv> ProxyVectorIterator<T, T_no_cv>::operator++(int)
{
    ProxyVectorIterator tmp(m_position);

    ++(*this);

    return tmp;
}

/*!
    Pre-decrement operator.
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv> & ProxyVectorIterator<T, T_no_cv>::operator--()
{
    m_position--;

    return *this;
}

/*!
    Post-decrement operator.
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv> ProxyVectorIterator<T, T_no_cv>::operator--(int)
{
    ProxyVectorIterator tmp(m_position);

    --(*this);

    return tmp;
}

/*!
    Compound assigment operator.

    \param increment is the increment
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv>& ProxyVectorIterator<T, T_no_cv>::operator+=(int increment)
{
    m_position += increment;

    return *this;
}

/*!
    Deference operator.

    \result A reference to the element currently pointed to by the
            iterator.
*/
template<typename T, typename T_no_cv>
__PXI_REFERENCE__ ProxyVectorIterator<T, T_no_cv>::operator*() const
{
    return *m_position;
}

/*!
    Deference operator.

    \result A reference to the element currently pointed to by the
            iterator.
*/
template<typename T, typename T_no_cv>
__PXI_POINTER__ ProxyVectorIterator<T, T_no_cv>::operator->() const
{
    return m_position;
}

/*!
    Converts the iterator to a const_iterator.
*/
template<typename T, typename T_no_cv>
template<typename U, typename std::enable_if<!std::is_const<U>::value, int>::type>
ProxyVectorIterator<T, T_no_cv>::operator ProxyVectorIterator<const U>() const
{
    return ProxyVectorIterator<const U>(m_position);
}

/*!
    Creates a new iterator and initializes it with the position of
    the const base iterator recevied in input.
*/
template<typename T, typename T_no_cv>
ProxyVectorIterator<T, T_no_cv>::ProxyVectorIterator(T *position)
    : m_position(position)
{
}

/*!
    Distance operator.

    \param other is the iterator from which the distance will be evaluated
    \result The distance between the specified iterator.
*/
template<typename T, typename T_no_cv>
std::size_t ProxyVectorIterator<T, T_no_cv>::operator-(const ProxyVectorIterator &other)
{
    return (m_position - other.m_position);
}

/*!
    Constructor
*/
template<typename T>
ProxyVector<T>::ProxyVector()
    : m_data(nullptr), m_size(0)
{
}

/*!
    Constructor

    \param size is the number elements contained in the data
    \param data a pointer to the data
*/
template<typename T>
ProxyVector<T>::ProxyVector(T *data, std::size_t size)
    : m_data(data), m_size(size)
{
}

/*!
    Constructor

    \param size is the number elements contained in the data
    \param data a pointer to the data
*/
template<typename T>
template<typename U, typename std::enable_if<std::is_const<U>::value, int>::type>
ProxyVector<T>::ProxyVector(std::vector<T_no_cv> &&storage)
    : m_storage(std::unique_ptr<std::vector<T_no_cv>>(new std::vector<T_no_cv>(std::move(storage)))),
      m_data(m_storage->data()), m_size(m_storage->size())
{
}

/*!
    Copy constructor.
*/
template<typename T>
ProxyVector<T>::ProxyVector(const ProxyVector &other)
{
    // Copy the elements
    std::vector<T_no_cv> new_storage(*(other.m_storage));

    // Assign the new memory to the object
    m_storage = std::unique_ptr<std::vector<T_no_cv>>(new std::vector<T_no_cv>(std::move(new_storage)));
    if (other.m_data == other.m_storage->data()) {
        m_data = m_storage->data();
    } else {
        m_data = other.m_data;
    }
    m_size = other.m_size;
}

/*!
    Copy assignment operator.

    Assigns new contents to the container, replacing its current contents,
    and modifying its size accordingly.
*/
template<typename T>
ProxyVector<T> & ProxyVector<T>::operator=(const ProxyVector &other)
{
    if (this != &other) {
        ProxyVector temporary(other);
        temporary.swap(*this);
    }

    return *this;
}

/*!
    Sets the content of the container.

    \param size is the number elements contained in the data
    \param data a pointer to the data
*/
template<typename T>
void ProxyVector<T>::set(T *data, std::size_t size)
{
    m_data = data;
    m_size = size;
}

/*!
    Sets the content of the container.

    \param size is the number elements contained in the data
    \param data a pointer to the data
*/
template<typename T>
template<typename U, typename std::enable_if<std::is_const<U>::value, int>::type>
void ProxyVector<T>::set(std::vector<T_no_cv> &&storage)
{
    m_storage = std::unique_ptr<std::vector<T_no_cv>>(new std::vector<T_no_cv>(std::move(storage)));
    m_data    = m_storage->data();
    m_size    = m_storage->size();
}

/*!
    Clear content.
*/
template<typename T>
void ProxyVector<T>::clear()
{
    m_data = nullptr;
    m_size = 0;

    m_storage.reset();
}

/*!
    Swaps the contents.

    \param other is another container of the same type
*/
template<typename T>
void ProxyVector<T>::swap(ProxyVector &other)
{
    std::swap(m_data, other.m_data);
    std::swap(m_size, other.m_size);

    m_storage.swap(other.m_storage);
}

/*!
    Tests whether two containers are equal.

    \result true if the containers are equal, false otherwise.
*/
template<typename T>
bool ProxyVector<T>::operator==(const ProxyVector& rhs) const
{
    return m_size == rhs.m_size && m_data == rhs.m_data && m_storage == rhs.m_storage;
}

/*!
    Tests whether the container is empty.

    \result true if the container size is 0, false otherwise.
*/
template<typename T>
bool ProxyVector<T>::empty() const
{
    return size() == 0;
}

/*!
    Returns the number of elements in the container

    \result The number of elements in the container.
*/
template<typename T>
std::size_t ProxyVector<T>::size() const
{
    return m_size;
}

/*!
    Returns a direct constant pointer to the memory where the elments are
    stored.

    \result A direct constant pointer to the memory where the elments are
    stored.
*/
template<typename T>
__PXV_CONST_POINTER__ ProxyVector<T>::data() const noexcept
{
    return m_data;
}

/*!
    Returns a direct pointer to the memory where the elments are stored.

    \result A direct pointer to the memory where the elments are stored.
*/
template<typename T>
__PXV_POINTER__ ProxyVector<T>::data() noexcept
{
    return m_data;
}

/*!
    Returns a constant reference to the specified element.

    \param n is the position of the requested element
    \result A constant reference to the specified element.
*/
template<typename T>
__PXV_CONST_REFERENCE__ ProxyVector<T>::operator[](std::size_t n) const
{
    return m_data[n];
}

/*!
    Returns a reference to the specified element.

    \param n is the position of the requested element
    \result A reference to the specified element.
*/
template<typename T>
__PXV_REFERENCE__ ProxyVector<T>::operator[](std::size_t n)
{
    return m_data[n];
}

/*!
    Returns a constant reference to the specified element.

    \param n is the position of the requested element
    \result A constant reference to the specified element.
*/
template<typename T>
__PXV_CONST_REFERENCE__ ProxyVector<T>::at(std::size_t n) const
{
    return m_data[n];
}

/*!
    Returns a reference to the specified element.

    \param n is the position of the requested element
    \result A reference to the specified element.
*/
template<typename T>
__PXV_REFERENCE__ ProxyVector<T>::at(std::size_t n)
{
    return m_data[n];
}

/*!
    Gets a constant reference to the first element in the container.

    \result A constant reference to the first element in the container.
*/
template<typename T>
__PXV_CONST_REFERENCE__ ProxyVector<T>::front() const
{
    return m_data[0];
}

/*!
    Gets a reference to the first element in the container.

    \result A reference to the first element in the container.
*/
template<typename T>
__PXV_REFERENCE__ ProxyVector<T>::front()
{
    return m_data[0];
}

/*!
    Gets a constant reference to the last element in the container.

    \result A constant reference to the last element in the container.
*/
template<typename T>
__PXV_CONST_REFERENCE__ ProxyVector<T>::back() const
{
    return m_data[m_size - 1];
}

/*!
    Gets a reference to the last element in the container.

    \result A reference to the last element in the container.
*/
template<typename T>
__PXV_REFERENCE__ ProxyVector<T>::back()
{
    return m_data[m_size - 1];
}

/*!
    Returns an iterator pointing to the first element in the container.

    \result An iterator pointing to the first element in the container.
*/
template<typename T>
__PXV_ITERATOR__ ProxyVector<T>::begin()
{
    return iterator(m_data);
}

/*!
    Returns an iterator referring to the past-the-end element in the container.

    \result An iterator referring to the past-the-end element in the container.
*/
template<typename T>
__PXV_ITERATOR__ ProxyVector<T>::end()
{
    return iterator(m_data + m_size);
}

/*!
    Returns a constant iterator pointing to the first element in the container.

    \result A constant iterator pointing to the first element in the container.
*/
template<typename T>
__PXV_CONST_ITERATOR__ ProxyVector<T>::begin() const
{
    return const_iterator(m_data);
}

/*!
    Returns a constant iterator referring to the past-the-end element in the
    container.

    \result A constant iterator referring to the past-the-end element in the
    container.
*/
template<typename T>
__PXV_CONST_ITERATOR__ ProxyVector<T>::end() const
{
    return const_iterator(m_data + m_size);
}

/*!
    Returns a constant iterator pointing to the first element in the container.

    \result A constant iterator pointing to the first element in the container.
*/
template<typename T>
__PXV_CONST_ITERATOR__ ProxyVector<T>::cbegin()
{
    return const_iterator(m_data);
}

/*!
    Returns a constant iterator referring to the past-the-end element in the
    container.

    \result A constant iterator referring to the past-the-end element in the
    container.
*/
template<typename T>
__PXV_CONST_ITERATOR__ ProxyVector<T>::cend()
{
    return const_iterator(m_data + m_size);
}

}

#endif
