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

#ifndef __BITPIT_FLAT_VECTOR_2D_TPP__
#define __BITPIT_FLAT_VECTOR_2D_TPP__

#include <vector>
#include <cassert>
#include <iostream>
#include <memory>

#include "binary_stream.hpp"

namespace bitpit {

/*!
    Stream operator from class FlatVector2D to communication buffer.
    Stream data from vector to communication buffer

    \param[in] buffer is the output memory stream
    \param[in] vetor is the container to be streamed
    \result Returns the same output stream received in input.
*/
template<class T>
OBinaryStream& operator<<(OBinaryStream &buffer, const FlatVector2D<T> &vector)
{
    buffer << vector.m_index;
    buffer << vector.m_v;

    return buffer;
}

/*!
    Input stream operator from Communication buffer for class FlatVector2D.
    Stream data from communication buffer to vector.

    \param[in] buffer is the input memory stream
    \param[in] vector is the container to be streamed
    \result Returns the same input stream received in input.
*/
template<class T>
IBinaryStream& operator>>(IBinaryStream &buffer, FlatVector2D<T> &vector)
{
    buffer >> vector.m_index;
    buffer >> vector.m_v;

    return buffer;
}

/*!
    Default constructor.
*/
template <class T>
FlatVector2D<T>::FlatVector2D(bool initialize)
    : m_index(initialize ? 1 : 0, 0L)
{
}

/*!
    Creates a new container.

    \param sizes are the sizes of the vectors
    \param value is the value that will be use to initialize the items of
    the vectors
*/
template <class T>
FlatVector2D<T>::FlatVector2D(const std::vector<std::size_t> &sizes, const T &value)
{
    initialize(sizes.size(), sizes.data(), 1, &value, 0);
}

/*!
    Creates a new container.

    \param nVectors is the number of vectors
    \param size is the size of the vectors
    \param value is the value that will be use to initialize the
    items of the vectors
*/
template <class T>
FlatVector2D<T>::FlatVector2D(std::size_t nVectors, std::size_t size, const T &value)
{
    initialize(nVectors, &size, 0, &value, 0);
}

/*!
    Creates a new container.

    \param nVectors is the number of vectors
    \param sizes are the sizes of the vectors
    \param value is the value that will be use to initialize the
    items of the vectors
*/
template <class T>
FlatVector2D<T>::FlatVector2D(std::size_t nVectors, const std::size_t *sizes, const T &value)
{
    initialize(nVectors, sizes, 1, &value, 0);
}

/*!
    Creates a new container.

    \param nVectors is the number of vectors
    \param sizes are the sizes of the vectors
    \param values are the values of the vectors
*/
template <class T>
FlatVector2D<T>::FlatVector2D(std::size_t nVectors, const std::size_t *sizes, const T *values)
{
    initialize(nVectors, sizes, 1, values, 1);
}

/*!
    Creates a new container.

    \param vector2D is a 2D vector that will be used to initialize the
    newly created container
*/
template <class T>
FlatVector2D<T>::FlatVector2D(const std::vector<std::vector<T> > &vector2D)
{
    initialize(vector2D);
}

/*!
    Check if the container has been initialized.

    \result Returns true if the container has been initialized, false otherwise.
*/
template <class T>
bool FlatVector2D<T>::isInitialized() const
{
    return (!m_index.empty());
}

/*!
    Initializes the container.

    \param sizes are the sizes of the vectors
    \param value is the value that will be use to initialize the items of
    the vectors
*/
template <class T>
void FlatVector2D<T>::initialize(const std::vector<std::size_t> &sizes, const T &value)
{
    initialize(sizes.size(), sizes.data(), 1, &value, 1);
}

/*!
    Initializes the container.

    \param nVectors is the number of vectors
    \param size is the size of the vectors
    \param value is the value that will be use to initialize the
    items of the vectors
*/
template <class T>
void FlatVector2D<T>::initialize(std::size_t nVectors, std::size_t size, const T &value)
{
    initialize(nVectors, &size, 0, &value, 0);
}

/*!
    Initializes the container.

    \param nVectors is the number of vectors
    \param sizes are the sizes of the vectors
    \param value is the value that will be use to initialize the
    items of the vectors
*/
template <class T>
void FlatVector2D<T>::initialize(std::size_t nVectors, const std::size_t *sizes, const T &value)
{
    initialize(nVectors, sizes, 1, &value, 0);
}


/*!
    Initializes the container.

    \param nVectors is the number of vectors
    \param sizes are the sizes of the vectors
    \param values are the values of the vectors
*/
template <class T>
void FlatVector2D<T>::initialize(std::size_t nVectors, const std::size_t *sizes, const T *values)
{
    initialize(nVectors, sizes, 1, values, 1);
}

/*!
    Initializes the container.

    \param nVectors is the number of vectors
    \param sizes are the sizes of the vectors
    \param sizesStride is the stride for accessing the sizes
    \param values are the values of each vector
    \param valuesStride is the stride for accessing the values
*/
template <class T>
void FlatVector2D<T>::initialize(std::size_t nVectors,
                                 const std::size_t *sizes, std::size_t sizesStride,
                                 const T *values, std::size_t valuesStride)
{
    std::size_t nItems = 0;
    for (std::size_t i = 0; i < nVectors; ++i) {
        nItems += *(sizes + i * sizesStride);
    }

    // Check if the container need reallocation
    //
    // If the current number of items is different from the number of items
    // the container should contain after the initialization, a reallocation
    // is needed.
    bool reallocateIndex  = (nVectors != size());
    bool reallocateValues = (nItems != getItemCount());

    // Destroy current data structures
    //
    // Index and the storage data will probabily be stored in memory
    // regions contigous to each other. To reduce memory fragmentation
    // it's better to deallocate the container before updating its
    // data structures.
    if (reallocateIndex || reallocateValues) {
        destroy(reallocateIndex, reallocateValues);
    }

    // Initialize the indexes
    if (reallocateIndex) {
        m_index.resize(nVectors + 1);
    }

    m_index[0] = 0;
    for (std::size_t i = 0; i < nVectors; ++i) {
        m_index[i+1] = m_index[i] + *(sizes + i * sizesStride);
    }

    // Initialize the storage
    if (valuesStride == 0) {
        if (reallocateValues) {
            m_v.assign(nItems, *values);
        } else {
            for (std::size_t k = 0; k < nItems; ++k) {
                m_v[k] = *values;
            }
        }
    } else {
        if (reallocateValues) {
            m_v.resize(nItems);
        }

        std::size_t k = 0;
        for (std::size_t i = 0; i < nVectors; ++i) {
            std::size_t subArraySize = *(sizes + i * sizesStride);
            for (std::size_t j = 0; j < subArraySize; ++j) {
                m_v[k++] = *(values + i * valuesStride + j);
            }
        }
    }
}

/*!
    Initializes the container.

    \param vector2D is a 2D vector that will be used to initialize the
    container
*/
template <class T>
void FlatVector2D<T>::initialize(const std::vector<std::vector<T> > &vector2D)
{
    std::size_t nVectors = vector2D.size();

    std::size_t nItems = 0;
    for (std::size_t i = 0; i < nVectors; ++i) {
        nItems += vector2D[i].size();
    }

    // Check if the container need reallocation
    //
    // If the current number of items is different from the number of items
    // the container should contain after the initialization, a reallocation
    // is needed.
    bool reallocateIndex  = (nVectors != size());
    bool reallocateValues = (nItems != getItemCount());

    // Destroy current data structures
    //
    // Index and the storage data will probabily be stored in memory
    // regions contigous to each other. To reduce memory fragmentation
    // it's better to deallocate the container before updating its
    // data structures.
    if (reallocateIndex || reallocateValues) {
        destroy(reallocateIndex, reallocateValues);
    }

    // Initialize the indexes
    if (reallocateIndex) {
        m_index.resize(nVectors + 1);
    }

    m_index[0] = 0;
    for (std::size_t i = 0; i < nVectors; ++i) {
        m_index[i+1] = m_index[i] + vector2D[i].size();
    }

    // Initialize the storage
    if (reallocateValues) {
        m_v.resize(nItems);
    }

    std::size_t k = 0;
    for (std::size_t i = 0; i < nVectors; ++i) {
        std::size_t subArraySize = vector2D[i].size();
        for (std::size_t j = 0; j < subArraySize; ++j) {
            m_v[k++] = vector2D[i][j];
        }
    }
}

/*!
    Initializes the container.

    \param other is antoher container of the same type, whose contents will
    be used to initialize the current container
*/
template <class T>
void FlatVector2D<T>::initialize(const FlatVector2D<T> &other)
{
    m_v.assign(other.m_v.begin(), other.m_v.end());
    m_index.assign(other.m_index.begin(), other.m_index.end());
}

/*!
    Destroy the container.

    After calling this function the container will be non-functional
    until it is re-initialized.
*/
template <class T>
void FlatVector2D<T>::destroy()
{
    destroy(true, true);
}

/*!
    Destroy the container.

    After calling this function the container will be non-functional
    until it is re-initialized.

    \param destroyIndex if true the index data structure will be destoryed
    \param destroyValues if true the values data structure will be destoryed
*/
template <class T>
void FlatVector2D<T>::destroy(bool destroyIndex, bool destroyValues)
{
    if (destroyIndex) {
        m_index.clear();
        m_index.shrink_to_fit();
    }

    if (destroyValues) {
        m_v.clear();
        m_v.shrink_to_fit();
    }
}

/*!
    Requests a change in capacity.

    Requests that the collpased-vector capacity be at least enough to
    contain nVectors vectors and nItems items.

    \param nVectors is the minimum number of vectors that the container
    should be able to contain
    \param nItems is the minimum number of items that the container should
    be able to contain
*/
template <class T>
void FlatVector2D<T>::reserve(std::size_t nVectors, std::size_t nItems)
{
    m_index.reserve(nVectors + 1);
    if (nItems > 0) {
        m_v.reserve(nItems);
    }
}

/*!
    Swaps the contents.

    \param other is another container of the same type
*/
template <class T>
void FlatVector2D<T>::swap(FlatVector2D &other) noexcept
{
    m_index.swap(other.m_index);
    m_v.swap(other.m_v);
}

/*!
    Sets the specified value as the value for all the items in the
    container.

    \param value is the value to fill the container with
*/
template <class T>
void FlatVector2D<T>::fill(T &value)
{
    std::fill(m_v.begin(), m_v.end(), value);
}

/*!
    Tests whether two containers are equal.

    \result true if the containers are equal, false otherwise.
*/
template <class T>
bool FlatVector2D<T>::operator==(const FlatVector2D& rhs) const
{
    return m_index == rhs.m_index && m_v == rhs.m_v;
}

/*!
    Tests whether the container is empty.

    \result true if the container size is 0, false otherwise.
*/
template <class T>
bool FlatVector2D<T>::empty() const
{
    return size() == 0;
}

/*!
    Clears content.

    Removes all items from the container (which are destroyed), leaving
    the container with a size of 0.

    \param release if it's true the memory hold by the container will be
    released, otherwise the container will be cleared but its memory will
    not be relased
*/
template <class T>
void FlatVector2D<T>::clear(bool release)
{
    if (release) {
        std::vector<T>(0).swap(m_v);

        std::vector<size_t>(1, 0L).swap(m_index);
    } else {
        m_v.clear();

        m_index.resize(1);
        m_index[0] = 0;
    }
}

/*!
    Clears the items.

    Removes all items stored in the vectors (which are destroyed), leaving
    each vector with a size of 0.

    \param release if it's true the memory hold by the container will be
    released, otherwise the container will be cleared but its memory will
    not be relased
*/
template <class T>
void FlatVector2D<T>::clearItems(bool release)
{
    std::size_t nVectors = size();
    if (release) {
        std::vector<T>(0).swap(m_v);

        std::vector<size_t>(nVectors + 1, 0L).swap(m_index);
    } else {
        m_v.clear();

        m_index.assign(nVectors + 1, 0L);
    }
}

/*!
    Shrinks to fit

    Requests the container to reduce its capacity to fit its size.
*/
template <class T>
void FlatVector2D<T>::shrinkToFit()
{
    m_v.shrink_to_fit();
    m_index.shrink_to_fit();
}

/*!
    Returns a constant pointer to the first item in the vector used
    internally by the container to store the indices.

    \result A constant pointer to the first item in the vector used
    internally by the container to store the indices.
*/
template <class T>
const std::size_t * FlatVector2D<T>::indices() const noexcept
{
    return m_index.data();
}

/*!
    Returns a constant pointer to the first item in the vector used
    internally by the container to store the indices of the specified
    vector.

    \param i is the index of the vector
    \result A constant pointer to the first item in the vector used
    internally by the container to store the indices of the specified
    vector.
*/
template <class T>
const std::size_t * FlatVector2D<T>::indices(std::size_t i) const noexcept
{
    return (m_index.data() + i);
}

/*!
    Returns a direct pointer to the memory vector used internally by the
    container to store its items.

    \result A pointer to the first item in the vector used
            internally by the container.

*/
template <class T>
T * FlatVector2D<T>::data() noexcept
{
    return m_v.data();
}

/*!
    Returns a direct constant pointer to the memory vector used internally by
    the container to store its items.

    \result A constant pointer to the first item in the vector used
            internally by the container.

*/
template <class T>
const T * FlatVector2D<T>::data() const noexcept
{
    return m_v.data();
}

/*!
    Returns a constant reference to the vector used internally by the
    container to store its items.

    \result A constant reference to the vector used internally by the
    container.

*/
template <class T>
const std::vector<T> & FlatVector2D<T>::vector() const
{
    return m_v;
}

/*!
    Adds an empty vector at the end.

    Adds an empty vector at the end of the container, after its current
    last vector.
*/
template <class T>
void FlatVector2D<T>::pushBack()
{
    pushBack(0);
}

/*!
    Adds a vector with the specified size at the end.

    Adds a vector with the specified size at the end of the vector,
    after its current last item. The content of value is copied
    (or moved) to the new vector.

    \param subArraySize is the size of the vector
    \param value is the value to be copied (or moved) to the new
    item
*/
template <class T>
void FlatVector2D<T>::pushBack(std::size_t subArraySize, const T &value)
{
    std::size_t previousLastIndex = m_index.back();
    m_index.emplace_back(previousLastIndex + subArraySize);

    m_v.resize(m_v.size() + subArraySize, value);
}

/*!
    Adds the specified vector at the end.

    Adds the specified vector at the end of the vector, after its current
    last item.

    \param subArray is the vector that will be added
*/
template <class T>
void FlatVector2D<T>::pushBack(const std::vector<T> &subArray)
{
    pushBack(subArray.size(), subArray.data());
}

/*!
    Adds the specified array at the end.

    Adds the specified array at the end of the vector, after its current
    last item.

    \param subArraySize is the size of the sub array
    \param subArray is a pointer to the sub array will be added
*/
template <class T>
void FlatVector2D<T>::pushBack(std::size_t subArraySize, const T *subArray)
{
    std::size_t previousLastIndex = m_index.back();
    m_index.emplace_back(previousLastIndex + subArraySize);

    m_v.reserve(m_v.size() + subArraySize);
    for (std::size_t j = 0; j < subArraySize; j++) {
        m_v.emplace_back();
        T &storedValue = m_v.back();
        storedValue = subArray[j];
    }
}

/*!
    Adds an item to the last vector.

    Adds an item at the end of to the last vector.

    \param value is the value that will be added
*/
template <class T>
void FlatVector2D<T>::pushBackItem(const T& value)
{
    m_index.back()++;

    m_v.emplace_back(value);
}

/*!
    Adds an item to the last vector.

    Adds an item at the end of to the last vector.

    \param value is the value that will be added
*/
template <class T>
void FlatVector2D<T>::pushBackItem(T &&value)
{
    m_index.back()++;

    m_v.emplace_back(std::move(value));
}

/*!
    Adds an item to the specified vector.

    Adds an item at the end of to the specified last vector.

    \param i is the index of the vector
    \param value is the value that will be added
*/
template <class T>
void FlatVector2D<T>::pushBackItem(std::size_t i, const T &value)
{
    assert(isIndexValid(i));

    m_v.insert(m_v.begin() + m_index[i+1], value);

    std::size_t nIndexes = m_index.size();
    for (std::size_t k = i + 1; k < nIndexes; ++k) {
        m_index[k]++;
    }
}


/*!
    Adds an item to the specified vector.

    Adds an item at the end of to the specified last vector.

    \param i is the index of the vector
    \param value is the value that will be added
*/
template <class T>
void FlatVector2D<T>::pushBackItem(std::size_t i, T &&value)
{
    assert(isIndexValid(i));

    m_v.insert(m_v.begin() + m_index[i+1], std::move(value));

    std::size_t nIndexes = m_index.size();
    for (std::size_t k = i + 1; k < nIndexes; ++k) {
        m_index[k]++;
    }
}

/*!
    Deletes last vector.

    Removes the last vector in the container, effectively reducing the
    container size by one.
*/
template <class T>
void FlatVector2D<T>::popBack()
{
    if (size() == 0) {
        return;
    }

    m_index.pop_back();
    m_v.resize(m_index.back() + 1);
}

/*!
    Deletes last item from last vector.

    Removes the last item from the last vector in the container.
*/
template <class T>
void FlatVector2D<T>::popBackItem()
{
    if (getItemCount(size() - 1) == 0) {
        return;
    }

    m_index.back()--;
    m_v.resize(m_index.back() + 1);
}

/*!
    Deletes last item from specified vector.

    Removes the last item from the specified vector in the container.

    \param i is the index of the vector
*/
template <class T>
void FlatVector2D<T>::popBackItem(std::size_t i)
{
    assert(isIndexValid(i));

    if (getItemCount(i) == 0) {
        return;
    }

    m_v.erase(m_v.begin() + m_index[i+1] - 1);

    std::size_t nIndexes = m_index.size();
    for (std::size_t k = i + 1; k < nIndexes; ++k) {
        m_index[k]--;
    }
}

/*!
    Deletes specified vector.

    Removes from the container the specified vector, effectively reducing
    the container size by one.

    \param i is the index of the vector
*/
template <class T>
void FlatVector2D<T>::erase(std::size_t i)
{
    assert(isIndexValid(i));

    m_v.erase(m_v.begin() + m_index[i], m_v.begin() + m_index[i+1] - 1);
    m_index.erase(m_index.begin() + i + 1);
}

/*!
    Deletes the specified item from a vector.

    \param i is the index of the vector
    \param j is the index of the item that will be removed
*/
template <class T>
void FlatVector2D<T>::eraseItem(std::size_t i, std::size_t j)
{
    assert(isIndexValid(i, j));

    m_v.erase(m_v.begin() + m_index[i] + j);

    std::size_t nIndexes = m_index.size();
    for (std::size_t k = i + 1; k < nIndexes; ++k) {
        m_index[k]--;
    }
}

/*!
    Sets the value of the specified item in a vector.

    \param i is the index of the vector
    \param j is the index of the item that will be removed
    \param value is the value that will be set
*/
template <class T>
void FlatVector2D<T>::setItem(std::size_t i, std::size_t j, const T &value)
{
    assert(isIndexValid(i, j));
    (*this)[i][j] = value;
}

/*!
    Sets the value of the specified item in a vector.

    \param i is the index of the vector
    \param j is the index of the item that will be removed
    \param value is the value that will be set
*/
template <class T>
void FlatVector2D<T>::setItem(std::size_t i, std::size_t j, T &&value)
{
    assert(isIndexValid(i, j));
    (*this)[i][j] = std::move(value);
}

/*!
    Gets a reference of the specified item in a vector.

    \param i is the index of the vector
    \param j is the index of the item that will be removed
    \result A reference to the requested value.
*/
template <class T>
T & FlatVector2D<T>::getItem(std::size_t i, std::size_t j)
{
    assert(isIndexValid(i, j));
    return (*this)[i][j];
}

/*!
    Gets a constant reference of the specified item in a vector.

    \param i is the index of the vector
    \param j is the index of the item that will be removed
    \result A constant reference to the requested value.
*/
template <class T>
const T & FlatVector2D<T>::getItem(std::size_t i, std::size_t j) const
{
    assert(isIndexValid(i, j));
    return (*this)[i][j];
}

/*!
    Gets a constant pointer to the first item of the specified vector.

    \param i is the index of the vector
    \result A constant pointer to the first item of the specified vector.
*/
template <class T>
const T * FlatVector2D<T>::get(std::size_t i) const
{
    assert(!empty());
    assert(isIndexValid(i));
    return (*this)[i];
}

/*!
    Gets a pointer to the first item of the specified vector.

    \param i is the index of the vector
    \result A pointer to the first item of the specified vector.
*/
template <class T>
T * FlatVector2D<T>::get(std::size_t i)
{
    assert(!empty());
    assert(isIndexValid(i));
    return (*this)[i];
}

/*!
    Sets the value of the specified item in a vector.

    \param k is the raw index
    \param value is the value that will be set
*/
template <class T>
void FlatVector2D<T>::rawSetItem(std::size_t k, const T &value)
{
    m_v[k] = value;
}

/*!
    Sets the value of the specified item in a vector.

    \param k is the raw index
    \param value is the value that will be set
*/
template <class T>
void FlatVector2D<T>::rawSetItem(std::size_t k, T &&value)
{
    m_v[k] = std::move(value);
}

/*!
    Gets a reference of the specified item in a vector.

    \param k is the raw index
    \result A reference to the requested value.
*/
template <class T>
T & FlatVector2D<T>::rawGetItem(std::size_t k)
{
    return m_v[k];
}

/*!
    Gets a constant reference of the specified item in a vector.

    \param k is the raw index
    \result A constant reference to the requested value.
*/
template <class T>
const T & FlatVector2D<T>::rawGetItem(std::size_t k) const
{
    return m_v[k];
}

/*!
    Gets a pointer to the first item of the last vector.

    \result A pointer to the first item of the vector.
*/
template <class T>
T * FlatVector2D<T>::back()
{
    return get(size() - 1);
}

/*!
    Gets a pointer to the first item of the first vector.

    \result A pointer to the first item of the vector.
*/
template <class T>
T * FlatVector2D<T>::first()
{
    return get(0);
}

/*!
    Returns the number of vectors in the container

    \result The number of vectors in the container.
*/
template <class T>
std::size_t FlatVector2D<T>::size() const
{
    if (!isInitialized()) {
        return 0;
    }

    return (m_index.size() - 1);
}

/*!
    Returns the size of the storage space currently allocated for
    storing vectors, expressed in terms of items.

    \result The size of the storage space currently allocated for
    storing vectors, expressed in terms of items.
*/
template <class T>
std::size_t FlatVector2D<T>::capacity() const
{
    if (!isInitialized()) {
        return 0;
    }

    return m_index.capacity() - 1;
}

/*!
    Merge the arrays together.
*/
template <class T>
void FlatVector2D<T>::merge()
{
    if (size() == 0) {
        return;
    }

    m_index[1] = m_index.back();
    m_index.resize(2);
}

/*!
    Returns the total size of all the vectors.

    \result The total size of all the vectors.
*/
template <class T>
std::size_t FlatVector2D<T>::getItemCount() const
{
    if (!isInitialized()) {
        return 0;
    }

    return m_v.size();
}

/*!
    Returns the size of the specified vector.

    \param i is the index of the vector
    \result The size of the vector.
*/
template <class T>
std::size_t FlatVector2D<T>::getItemCount(std::size_t i) const
{
    if (!isInitialized()) {
        return 0;
    }

    return m_index[i + 1] - m_index[i];
}

/*!
    Returns the size of the storage space currently allocated for
    storing vectors items, expressed in terms of items.

    \result The size of the storage space currently allocated for
    storing vectors items, expressed in terms of items.
*/
template <class T>
std::size_t FlatVector2D<T>::getItemCapacity() const
{
    return m_v.capacity();
}

/*!
    Returns the buffer size (in bytes) required to store the container.

    \result The buffer size (in bytes) required to store the container.
*/
template <class T>
size_t FlatVector2D<T>::getBinarySize() const
{
    return ((2 + m_index.size())*sizeof(size_t) + m_v.size() * sizeof(T));
}

/*!
    Returns a constant pointer to the first item of the specified vector.

    \param i is the index of the vector
    \result A constant pointer to the first item of the specified vector.
*/
template <class T>
const T* FlatVector2D<T>::operator[](std::size_t i) const
{
    assert(isIndexValid(i));

    if (m_v.empty()) {
        return nullptr;
    }

    std::size_t index = m_index[i];
    return &m_v[index];
}

/*!
    Returns a pointer to the first item of the specified vector.

    \param i is the index of the vector
    \result A pointer to the first item of the specified vector.
*/
template <class T>
T* FlatVector2D<T>::operator[](std::size_t i)
{
    assert(isIndexValid(i));

    if (m_v.empty()) {
        return nullptr;
    }

    std::size_t index = m_index[i];
    return &m_v[index];
}

/*!
    Checks if the specified index is valid.

    \param i is the index of the vector
    \result true if the index is vaid, false otherwise.
*/
template <class T>
bool FlatVector2D<T>::isIndexValid(std::size_t i) const
{
    return (i < size());
}

/*!
    Checks if the specified indexes are valid

    \param i is the index of the vector
    \param j is the index of the item in the vector
    \result true if the indexes are vaid, false otherwise.
*/
template <class T>
bool FlatVector2D<T>::isIndexValid(std::size_t i, std::size_t j) const
{
    if (!isIndexValid(i)) {
        return false;
    }

    return (j < (m_index[i+1] - m_index[i]));
}

}

#endif
