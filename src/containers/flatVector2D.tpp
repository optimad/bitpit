/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __BITPIT_FLAT_VECTOR_2D_HPP__
#define __BITPIT_FLAT_VECTOR_2D_HPP__

#include <vector>
#include <cassert>
#include <iostream>
#include <memory>

#include "binary_stream.hpp"

namespace bitpit{
    template<class T>
    class FlatVector2D;
}

/*!
    Stream operator from class FlatVector2D to communication buffer.
    Stream data from vector to communication buffer

    \param[in] buffer is the output memory stream
    \param[in] vetor is the container to be streamed
    \result Returns the same output stream received in input.
*/
template<class T>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::FlatVector2D<T> &vector)
{
    typename std::vector<T>::const_iterator           it;
    typename std::vector<size_t>::const_iterator      jt;

    buffer << vector.m_index.size() << vector.m_v.size();
    for (jt = vector.m_index.begin(); jt != vector.m_index.end(); ++jt) {
        buffer << *(jt);
    }
    for (it = vector.m_v.begin(); it != vector.m_v.end(); ++it) {
        buffer << *(it);
    }

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
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::FlatVector2D<T> &vector)
{
    size_t                      size_m_v, size_m_index;
    size_t                      i;

    buffer >> size_m_index;
    buffer >> size_m_v;

    vector.m_index.resize(size_m_index, 0);
    vector.m_v.resize(size_m_v);

    for (i = 0; i < size_m_index; ++i) {
        buffer >> vector.m_index[i];
    }

    for (i = 0; i < size_m_v; ++i) {
        buffer >> vector.m_v[i];
    }

    return buffer;
}

namespace bitpit{
/*!
    \ingroup containers
    @{
*/

/*!
    @brief Metafunction for generation of a flattened vector of vectors.

    @details
    Usage: Use <tt>FlatVector2D<Type></tt> to declare a flattened vector of
    vectors.

    @tparam T The type of the objects stored in the vector
*/

template <class T>
class FlatVector2D
{

// Friendship(s)
template<class U>
friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream &buffer, const FlatVector2D<U> &vector);
template<class U>
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream &buffer, FlatVector2D<U> &vector);

// Members and methods
public:

    /*!
        Default constructor.
    */
    FlatVector2D(bool initialize = true)
        : m_index(initialize ? 1 : 0, 0L)
    {
    }

    /*!
        Creates a new container.

        \param sizes are the sizes of the vectors
        \param value is the value that will be use to initialize the items of
        the vectors
    */
    FlatVector2D(const std::vector<int> &sizes, const T &value = T())
    {
        initialize(sizes, value);
    }

    /*!
        Creates a new container.

        \param nVectors is the number of vectors
        \param size is the size of the vectors
        \param value is the value that will be use to initialize the
        items of the vectors
    */
    FlatVector2D(const int &nVectors, const int &size, const T &value = T())
    {
        initialize(nVectors, size, value);
    }

    /*!
        Creates a new container.

        \param vector2D is a 2D vector that will be used to initialize the
        newly created container
    */
    FlatVector2D(const std::vector<std::vector<T> > &vector2D)
    {
        initialize(vector2D);
    }

    /*!
        Copy constructor.
    */
    FlatVector2D(const FlatVector2D &other)
    {
        // Copy the items
        std::vector<T> new_v(other.m_v);
        std::vector<std::size_t> new_index(other.m_index);

        // Assign the new memory to the object
        m_v.swap(new_v);
        m_index.swap(new_index);
    }

    /*!
        Copy assignment operator.

        Assigns new contents to the container, replacing its current contents,
        and modifying its size accordingly.
    */
    FlatVector2D & operator=(const FlatVector2D &other)
    {
        if (this != &other) {
            FlatVector2D temporary(other);
            temporary.swap(*this);
        }

        return *this;
    }

    /*!
        Move assignment operator.

        The move assignment operator "steals" the resources held by the
        argument.
    */
    FlatVector2D & operator=(FlatVector2D &&other) = default;

    /*!
        Initializes the container.

        \param sizes are the sizes of the vectors
        \param value is the value that will be use to initialize the items of
        the vectors
    */
    void initialize(const std::vector<int> &sizes, const T &value = T())
    {
        int nVectors = sizes.size();

        // Destroy the container
        //
        // Index and the storage data will probabily be stored in memory
        // regions contigous to each other. To reduce memory fragmentation
        // it's better to deallocate the container before updating its
        // data structures.
        destroy();

        // Initialize the indexes
        //
        // Clear the vector before the resize to reduce memory fragmentation
        m_index.resize(nVectors + 1);
        m_index.shrink_to_fit();
        for (int i = 0; i < nVectors; ++i) {
            m_index[i+1] = m_index[i] + sizes[i];
        }

        // Initialize the storage
        //
        // Clear the vector before the resize to reduce memory fragmentation
        m_v.assign(m_index[nVectors], value);
        m_v.shrink_to_fit();
    }

    /*!
        Initializes the container.

        \param nVectors is the number of vectors
        \param size is the size of the vectors
        \param value is the value that will be use to initialize the items of
        the vectors
    */
    void initialize(const int &nVectors, const int &size, const T &value = T())
    {
        // Destroy the container
        //
        // Index and the storage data will probabily be stored in memory
        // regions contigous to each other. To reduce memory fragmentation
        // it's better to deallocate the container before updating its
        // data structures.
        destroy();

        // Initialize the indexes
        m_index.resize(nVectors + 1);
        m_index.shrink_to_fit();
        for (int i = 0; i < nVectors; ++i) {
            m_index[i+1] = m_index[i] + size;
        }

        // Initialize the storage
        m_v.assign(m_index[nVectors], value);
        m_v.shrink_to_fit();
    }

    /*!
        Initializes the container.

        \param vector2D is a 2D vector that will be used to initialize the
        container
    */
    void initialize(const std::vector<std::vector<T> > &vector2D)
    {
        int nVectors = vector2D.size();

        // Destroy the container
        //
        // Index and the storage data will probabily be stored in memory
        // regions contigous to each other. To reduce memory fragmentation
        // it's better to deallocate the container before updating its
        // data structures.
        destroy();

        // Initialize the indexes
        m_index.resize(nVectors + 1);
        m_index.shrink_to_fit();
        for (int i = 0; i < nVectors; ++i) {
            m_index[i+1] = m_index[i] + vector2D[i].size();
        }

        // Initialize the storage
        m_v.resize(m_index[nVectors]);
        m_v.shrink_to_fit();

        int k = 0;
        for (int i = 0; i < nVectors; ++i) {
            int subArraySize = vector2D[i].size();
            for (int j = 0; j < subArraySize; ++j) {
                m_v[k++] = vector2D[i][j];
            }
        }
    }

    /*!
        Destroy the container.

        After calling this function the container will be non-functional
        until it is re-initialized.

    */
    void destroy()
    {
        m_index.clear();
        m_index.shrink_to_fit();

        m_v.clear();
        m_v.shrink_to_fit();
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
    void reserve(int nVectors, int nItems = 0)
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
    void swap(FlatVector2D &other)
    {
        m_index.swap(other.m_index);
        m_v.swap(other.m_v);
    }

    /*!
        Sets the specified value as the value for all the items in the
        container.

        \param value is the value to fill the container with
    */
    void fill(T &value)
    {
        std::fill(m_v.begin(), m_v.end(), value);
    }

    /*!
        Tests whether two containers are equal.

        \result true if the containers are equal, false otherwise.
    */
    bool operator==(const FlatVector2D& rhs) const
    {
        return m_index == rhs.m_index && m_v == rhs.m_v;
    }

    /*!
        Tests whether the container is empty.

        \result true if the container size is 0, false otherwise.
    */
    bool empty() const
    {
        return size() == 0;
    }

    /*!
        Clears content.

        Removes all items from the container (which are destroyed), leaving
        the container with a size of 0.
    */
    void clear()
    {
        std::vector<T>(0).swap(m_v);

        std::vector<size_t>(1, 0L).swap(m_index);
    }

    /*!
        Shrinks to fit

        Requests the container to reduce its capacity to fit its size.
    */
    void shrinkToFit()
    {
        m_v.shrink_to_fit();
        m_index.shrink_to_fit();
    }

    /*!
        Returns a direct pointer to the memory vector used internally by the
        container to store its items.

        \result A pointer to the first item in the vector used
                internally by the container.

    */
    T * data() noexcept
    {
        return m_v.data();
    }

    /*!
        Returns a constant reference to the vector used internally by the
        container to store its items.

        \result A constant reference to the vector used internally by the
        container.

    */
    const std::vector<T> & vector() const
    {
        return m_v;
    }

    /*!
        Adds an empty vector at the end.

        Adds an empty vector at the end of the container, after its current
        last vector.
    */
    void pushBack()
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
    void pushBack(const int &subArraySize, const T &value = T())
    {
        std::size_t previousLastIndex = m_index.back();
        m_index.emplace_back();
        std::size_t &lastIndex = m_index.back();
        lastIndex = previousLastIndex + subArraySize;

        m_v.resize(m_v.size() + subArraySize, value);
    }

    /*!
        Adds the specified vector at the end.

        Adds the specified vector at the end of the vector, after its current
        last item.

        \param subArray is the vector that will be added
    */
    void pushBack(const std::vector<T> &subArray)
    {
        int subArraySize = subArray.size();

        std::size_t previousLastIndex = m_index.back();
        m_index.emplace_back();
        std::size_t &lastIndex = m_index.back();
        lastIndex = previousLastIndex + subArraySize;

        m_v.reserve(m_v.size() + subArraySize);
        for (int j = 0; j < subArraySize; j++) {
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
    void pushBackItem(const T& value)
    {
        m_index.back()++;

        m_v.emplace_back();
        T &storedValue = m_v.back();
        storedValue = value;
    }

    /*!
        Adds an item to the specified vector.

        Adds an item at the end of to the specified last vector.

        \param i is the index of the vector
        \param value is the value that will be added
    */
    void pushBackItem(const int &i, const T& value)
    {
        assert(indexValid(i));

        m_v.insert(m_v.begin() + m_index[i+1], value);

        int nIndexes = m_index.size();
        for (int k = i + 1; k < nIndexes; ++k) {
            m_index[k]++;
        }
    }

    /*!
        Deletes last vector.

        Removes the last vector in the container, effectively reducing the
        container size by one.
    */
    void popBack()
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
    void popBackItem()
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
    void popBackItem(const int &i)
    {
        assert(indexValid(i));

        if (getItemCount(i) == 0) {
            return;
        }

        m_v.erase(m_v.begin() + m_index[i+1] - 1);

        int nIndexes = m_index.size();
        for (int k = i + 1; k < nIndexes; ++k) {
            m_index[k]--;
        }
    }

    /*!
        Deletes specified vector.

        Removes from the container the specified vector, effectively reducing
        the container size by one.

        \param i is the index of the vector
    */
    void erase(const int &i)
    {
        assert(indexValid(i));

        m_v.erase(m_v.begin() + m_index[i], m_v.begin() + m_index[i+1] - 1);
        m_index.erase(m_index.begin() + i + 1);
    }

    /*!
        Deletes the specified item from a vector.

        \param i is the index of the vector
        \param j is the index of the item that will be removed
    */
    void eraseItem(const int &i, const int &j)
    {
        assert(indexValid(i, j));

        m_v.erase(m_v.begin() + m_index[i] + j);

        int nIndexes = m_index.size();
        for (int k = i + 1; k < nIndexes; ++k) {
            m_index[k]--;
        }
    }

    /*!
        Sets the value of the specified item in a vector.

        \param i is the index of the vector
        \param j is the index of the item that will be removed
        \param value is the value that will be set
    */
    void setItem(const int &i, const int &j, const T &value)
    {
        assert(indexValid(i, j));
        (*this)[i][j] = value;
    }

    /*!
        Gets a reference of the specified item in a vector.

        \param i is the index of the vector
        \param j is the index of the item that will be removed
        \result A reference to the requested value.
    */
    T & getItem(const int &i, const int &j)
    {
        assert(indexValid(i, j));
        return (*this)[i][j];
    }

    /*!
        Gets a constant reference of the specified item in a vector.

        \param i is the index of the vector
        \param j is the index of the item that will be removed
        \result A constant reference to the requested value.
    */
    const T & getItem(const int &i, const int &j) const
    {
        assert(indexValid(i, j));
        return (*this)[i][j];
    }

    /*!
        Gets a pointer to the first item of the specified vector.

        \param i is the index of the vector
        \result A pointer to the first item of the specified vector.
    */
    const T * get(const int &i) const
    {
        assert(!empty());
        assert(indexValid(i));
        return (*this)[i];
    }

    /*!
        Gets a pointer to the first item of the last vector.

        \result A pointer to the first item of the vector.
    */
    T * back()
    {
        return get(size() - 1);
    }

    /*!
        Gets a pointer to the first item of the first vector.

        \result A pointer to the first item of the vector.
    */
    T * first()
    {
        return get(0);
    }

    /*!
        Returns the number of vectors in the container

        \result The number of vectors in the container.
    */
    int size() const
    {
        return (m_index.size() - 1);
    }

    /*!
        Returns the size of the storage space currently allocated for
        storing vectors, expressed in terms of items.

        \result The size of the storage space currently allocated for
        storing vectors, expressed in terms of items.
    */
    int capacity() const
    {
        return m_index.capacity() - 1;
    }

    /*!
        Returns the total size of all the vectors.

        \result The total size of all the vectors.
    */
    int getItemCount() const
    {
        return m_index[size()];
    }

    /*!
        Returns the size of the specified vector.

        \param i is the index of the vector
        \result The size of the vector.
    */
    int getItemCount(int i) const
    {
        return m_index[i + 1] - m_index[i];
    }

    /*!
        Returns the size of the storage space currently allocated for
        storing vectors items, expressed in terms of items.

        \result The size of the storage space currently allocated for
        storing vectors items, expressed in terms of items.
    */
    int getItemCapacity() const
    {
        return m_v.capacity();
    }

    /*!
        Returns the buffer size (in bytes) required to store the container.

        \result The buffer size (in bytes) required to store the container.
    */
    size_t getBinarySize()
    {
        return ((2 + m_index.size())*sizeof(size_t) + m_v.size() * sizeof(T));
    }

private:
    std::vector<T> m_v;
    std::vector<std::size_t> m_index;

    /*!
        Returns a constant pointer to the first item of the specified vector.

        \param i is the index of the vector
        \result A constant pointer to the first item of the specified vector.
    */
    const T* operator[](const int &i) const
    {
        assert(indexValid(i));

        int index = m_index[i];
        return &m_v[index];
    }

    /*!
        Returns a pointer to the first item of the specified vector.

        \param i is the index of the vector
        \result A pointer to the first item of the specified vector.
    */
    T* operator[](const int &i)
    {
        assert(indexValid(i));

        int index = m_index[i];
        return &m_v[index];
    }

    /*!
        Checks if the specified index is valid.

        \param i is the index of the vector
        \result true if the index is vaid, false otherwise.
    */
    bool indexValid(const int &i) const
    {
        return (i >= 0 && i < size());
    }

    /*!
        Checks if the specified indexes are valid

        \param i is the index of the vector
        \param j is the index of the item in the vector
        \result true if the indexes are vaid, false otherwise.
    */
    bool indexValid(const int &i, const int &j) const
    {
        if (!indexValid(i)) {
            return false;
        }

        return (j >= 0 && j < (int) (m_index[i+1] - m_index[i]));
    }

};

/*!
    @}
*/

}

#endif
