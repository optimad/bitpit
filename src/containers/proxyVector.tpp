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
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t>::ProxyVectorIterator()
    : m_position(nullptr)
{
}

/*!
    Constructor

    This constructor allows to generate a constant iterator from a non
    constnat iterator.

    \param other is the iterator that will be copied
*/
template<typename value_t, typename value_no_cv_t>
template<typename other_value_t, typename std::enable_if<std::is_const<value_t>::value && !std::is_const<other_value_t>::value && std::is_same<other_value_t, typename std::remove_cv<value_t>::type>::value, int>::type>
ProxyVectorIterator<value_t, value_no_cv_t>::ProxyVectorIterator(const ProxyVectorIterator<other_value_t> &other)
    : m_position(other.m_position)
{
}

/*!
    Exchanges the values of the current iterator and
    the iterator recevied as argument.

    \param other is the iterator to exchange values with
*/
template<typename value_t, typename value_no_cv_t>
void ProxyVectorIterator<value_t, value_no_cv_t>::swap(ProxyVectorIterator& other) noexcept
{
    std::swap(m_position, other.m_position);
}

/*!
    Pre-increment operator.
*/
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t> & ProxyVectorIterator<value_t, value_no_cv_t>::operator++()
{
    m_position++;

    return *this;
}

/*!
    Post-increment operator.
*/
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t> ProxyVectorIterator<value_t, value_no_cv_t>::operator++(int)
{
    ProxyVectorIterator tmp(m_position);

    ++(*this);

    return tmp;
}

/*!
    Pre-decrement operator.
*/
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t> & ProxyVectorIterator<value_t, value_no_cv_t>::operator--()
{
    m_position--;

    return *this;
}

/*!
    Post-decrement operator.
*/
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t> ProxyVectorIterator<value_t, value_no_cv_t>::operator--(int)
{
    ProxyVectorIterator tmp(m_position);

    --(*this);

    return tmp;
}

/*!
    Compound assigment operator.

    \param increment is the increment
*/
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t>& ProxyVectorIterator<value_t, value_no_cv_t>::operator+=(int increment)
{
    m_position += increment;

    return *this;
}

/*!
    Deference operator.

    \result A reference to the element currently pointed to by the iterator.
*/
template<typename value_t, typename value_no_cv_t>
__PXI_REFERENCE__ ProxyVectorIterator<value_t, value_no_cv_t>::operator*() const
{
    return *m_position;
}

/*!
    Deference operator.

    \result A reference to the element currently pointed to by the iterator.
*/
template<typename value_t, typename value_no_cv_t>
__PXI_POINTER__ ProxyVectorIterator<value_t, value_no_cv_t>::operator->() const
{
    return m_position;
}

/*!
* Copy assignment operator to create a constant iterator from a non-constant
* one.
*
* \param other is the iterator that will be copied
*/
template<typename value_t, typename value_no_cv_t>
template<typename other_value_t, typename std::enable_if<std::is_const<value_t>::value && !std::is_const<other_value_t>::value && std::is_same<other_value_t, typename std::remove_cv<value_t>::type>::value, int>::type>
ProxyVectorIterator<value_t, value_no_cv_t> & ProxyVectorIterator<value_t, value_no_cv_t>::operator=(const ProxyVectorIterator<other_value_t> &other)
{
    m_position = other.m_position;

    return *this;
}

/*!
    Creates a new iterator and initializes it with the position of
    the const base iterator recevied in input.
*/
template<typename value_t, typename value_no_cv_t>
ProxyVectorIterator<value_t, value_no_cv_t>::ProxyVectorIterator(__PXI_POINTER__ position)
    : m_position(position)
{
}

/*!
    Distance operator.

    \param other is the iterator from which the distance will be evaluated
    \result The distance between the specified iterator.
*/
template<typename value_t, typename value_no_cv_t>
std::size_t ProxyVectorIterator<value_t, value_no_cv_t>::operator-(const ProxyVectorIterator &other)
{
    return (m_position - other.m_position);
}

/*!
    Constructor
*/
template<typename value_t>
ProxyVector<value_t>::ProxyVector()
    : m_data(nullptr), m_size(0)
{
}

/*!
    Constructor

    \param size is the number elements contained in the data
    \param data a pointer to the data
*/
template<typename value_t>
ProxyVector<value_t>::ProxyVector(value_t *data, std::size_t size)
    : m_data(data), m_size(size)
{
}

/*!
    Constructor

    \param storage is the storage that contains the data
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
ProxyVector<value_t>::ProxyVector(std::vector<value_no_cv_t> &&storage)
    : m_storage(std::unique_ptr<std::vector<value_no_cv_t>>(new std::vector<value_no_cv_t>(std::move(storage)))),
      m_data(m_storage->data()), m_size(m_storage->size())
{
}

/*!
    Copy constructor.
*/
template<typename value_t>
ProxyVector<value_t>::ProxyVector(const ProxyVector &other)
    : m_size(other.m_size)
{
    if (other.m_storage) {
        if (m_storage) {
            m_storage->assign(other.m_storage->begin(), other.m_storage->end());
        } else {
            m_storage = std::unique_ptr<std::vector<value_no_cv_t>>(new std::vector<value_no_cv_t>(*(other.m_storage)));
        }

        if (other.m_data == other.m_storage->data()) {
            m_data = m_storage->data();
        } else {
            m_data = other.m_data;
        }
    } else {
        m_storage.reset();

        m_data = other.m_data;
    }
}

/*!
    Copy assignment operator.

    Assigns new contents to the container, replacing its current contents,
    and modifying its size accordingly.
*/
template<typename value_t>
ProxyVector<value_t> & ProxyVector<value_t>::operator=(const ProxyVector &other)
{
    if (this != &other) {
        ProxyVector temporary(other);
        temporary.swap(*this);
    }

    return *this;
}

/*!
    Sets the content of the container.

    \param data a pointer to the data
    \param size is the number elements contained in the data
*/
template<typename value_t>
void ProxyVector<value_t>::set(value_t *data, std::size_t size)
{
    m_data = data;
    m_size = size;
}

/*!
    Sets the content of the container.

    \param storage is the storage that contains the data
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
typename ProxyVector<value_t>::value_no_cv_t * ProxyVector<value_t>::set(std::vector<value_no_cv_t> &&storage)
{
    m_storage = std::unique_ptr<std::vector<value_no_cv_t>>(new std::vector<value_no_cv_t>(std::move(storage)));
    m_data    = m_storage->data();
    m_size    = m_storage->size();

    return m_storage->data();
}

/*!
    Sets the content of the container.

    \param size is the number of elements the storage should be able to
    contain
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<std::is_const<other_value_t>::value, int>::type>
typename ProxyVector<value_t>::value_no_cv_t * ProxyVector<value_t>::set(std::size_t size)
{
    if (!m_storage) {
        m_storage = std::unique_ptr<std::vector<value_no_cv_t>>(new std::vector<value_no_cv_t>(size));
    } else {
        m_storage->resize(size);
    }
    m_data = m_storage->data();
    m_size = m_storage->size();

    return m_storage->data();
}

/*!
    Clear content.
*/
template<typename value_t>
void ProxyVector<value_t>::clear()
{
    m_data = nullptr;
    m_size = 0;

    m_storage.reset();
}

/*!
    Swaps the content.

    \param other is another container of the same type
*/
template<typename value_t>
void ProxyVector<value_t>::swap(ProxyVector &other)
{
    std::swap(m_data, other.m_data);
    std::swap(m_size, other.m_size);

    m_storage.swap(other.m_storage);
}

/*!
    Tests whether two containers are equal.

    \result true if the containers are equal, false otherwise.
*/
template<typename value_t>
bool ProxyVector<value_t>::operator==(const ProxyVector& other) const
{
    return m_size == other.m_size && m_data == other.m_data && m_storage == other.m_storage;
}

/*!
    Tests whether the container is empty.

    \result true if the container size is 0, false otherwise.
*/
template<typename value_t>
bool ProxyVector<value_t>::empty() const
{
    return size() == 0;
}

/*!
    Returns the number of elements in the container

    \result The number of elements in the container.
*/
template<typename value_t>
std::size_t ProxyVector<value_t>::size() const
{
    return m_size;
}

/*!
    Returns a direct pointer to the memory where the elments are stored.

    \result A direct pointer to the memory where the elments are stored.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_POINTER__ ProxyVector<value_t>::data() noexcept
{
    return m_data;
}

/*!
    Returns a direct constant pointer to the memory where the elments are
    stored.

    \result A direct constant pointer to the memory where the elments are
    stored.
*/
template<typename value_t>
__PXV_CONST_POINTER__ ProxyVector<value_t>::data() const noexcept
{
    return m_data;
}

/*!
    Returns a reference to the specified element.

    \param n is the position of the requested element
    \result A reference to the specified element.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t>::operator[](std::size_t n)
{
    return m_data[n];
}

/*!
    Returns a constant reference to the specified element.

    \param n is the position of the requested element
    \result A constant reference to the specified element.
*/
template<typename value_t>
__PXV_CONST_REFERENCE__ ProxyVector<value_t>::operator[](std::size_t n) const
{
    return m_data[n];
}

/*!
    Returns a reference to the specified element.

    \param n is the position of the requested element
    \result A reference to the specified element.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t>::at(std::size_t n)
{
    return m_data[n];
}

/*!
    Returns a constant reference to the specified element.

    \param n is the position of the requested element
    \result A constant reference to the specified element.
*/
template<typename value_t>
__PXV_CONST_REFERENCE__ ProxyVector<value_t>::at(std::size_t n) const
{
    return m_data[n];
}

/*!
    Gets a reference to the first element in the container.

    \result A reference to the first element in the container.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t>::front()
{
    return m_data[0];
}

/*!
    Gets a constant reference to the first element in the container.

    \result A constant reference to the first element in the container.
*/
template<typename value_t>
__PXV_CONST_REFERENCE__ ProxyVector<value_t>::front() const
{
    return m_data[0];
}

/*!
    Gets a reference to the last element in the container.

    \result A reference to the last element in the container.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_REFERENCE__ ProxyVector<value_t>::back()
{
    return m_data[m_size - 1];
}

/*!
    Gets a constant reference to the last element in the container.

    \result A constant reference to the last element in the container.
*/
template<typename value_t>
__PXV_CONST_REFERENCE__ ProxyVector<value_t>::back() const
{
    return m_data[m_size - 1];
}

/*!
    Returns an iterator pointing to the first element in the container.

    \result An iterator pointing to the first element in the container.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_ITERATOR__ ProxyVector<value_t>::begin()
{
    return iterator(m_data);
}

/*!
    Returns a constant iterator pointing to the first element in the container.

    \result A constant iterator pointing to the first element in the container.
*/
template<typename value_t>
__PXV_CONST_ITERATOR__ ProxyVector<value_t>::begin() const
{
    return const_iterator(m_data);
}

/*!
    Returns an iterator referring to the past-the-end element in the container.

    \result An iterator referring to the past-the-end element in the container.
*/
template<typename value_t>
template<typename other_value_t, typename std::enable_if<!std::is_const<other_value_t>::value, int>::type>
__PXV_ITERATOR__ ProxyVector<value_t>::end()
{
    return iterator(m_data + m_size);
}

/*!
    Returns a constant iterator referring to the past-the-end element in the
    container.

    \result A constant iterator referring to the past-the-end element in the
    container.
*/
template<typename value_t>
__PXV_CONST_ITERATOR__ ProxyVector<value_t>::end() const
{
    return const_iterator(m_data + m_size);
}

/*!
    Returns a constant iterator pointing to the first element in the container.

    \result A constant iterator pointing to the first element in the container.
*/
template<typename value_t>
__PXV_CONST_ITERATOR__ ProxyVector<value_t>::cbegin()
{
    return const_iterator(m_data);
}

/*!
    Returns a constant iterator referring to the past-the-end element in the
    container.

    \result A constant iterator referring to the past-the-end element in the
    container.
*/
template<typename value_t>
__PXV_CONST_ITERATOR__ ProxyVector<value_t>::cend()
{
    return const_iterator(m_data + m_size);
}

}

#endif
