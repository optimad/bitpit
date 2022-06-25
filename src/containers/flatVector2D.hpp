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

template<class T>
OBinaryStream& operator<<(OBinaryStream &buffer, const FlatVector2D<T> &vector);

template<class T>
IBinaryStream& operator>>(IBinaryStream &buffer, FlatVector2D<T> &vector);

/*!
    @ingroup containers

    @brief Metafunction for generation of a flattened vector of vectors.

    @details
    Usage: Use <tt>FlatVector2D<Type></tt> to declare a flattened vector of
    vectors.

    @tparam T The type of the objects stored in the vector
*/

template <class T>
class FlatVector2D
{

template<class U>
friend OBinaryStream& (operator<<) (OBinaryStream &buffer, const FlatVector2D<U> &vector);
template<class U>
friend IBinaryStream& (operator>>) (IBinaryStream &buffer, FlatVector2D<U> &vector);

public:
    FlatVector2D(bool initialize = true);
    FlatVector2D(const std::vector<std::size_t> &sizes, const T &value = T());
    FlatVector2D(std::size_t nVectors, std::size_t size, const T &value = T());
    FlatVector2D(std::size_t nVectors, const std::size_t *sizes, const T &value);
    FlatVector2D(std::size_t nVectors, const std::size_t *sizes, const T *values);
    FlatVector2D(const std::vector<std::vector<T> > &vector2D);
    FlatVector2D(const FlatVector2D &other) = default;
    FlatVector2D(FlatVector2D &&other) = default;

    /*!
        Copy assignment operator.

        Assigns new contents to the container, replacing its current contents,
        and modifying its size accordingly.
    */
    FlatVector2D & operator=(const FlatVector2D &other) = default;

    /*!
        Move assignment operator.

        The move assignment operator "steals" the resources held by the
        argument.
    */
    FlatVector2D & operator=(FlatVector2D &&other) = default;

    bool isInitialized() const;

    void initialize(const std::vector<std::size_t> &sizes, const T &value = T());
    void initialize(std::size_t nVectors, std::size_t size, const T &value = T());
    void initialize(std::size_t nVectors, const std::size_t *sizes, const T &value);
    void initialize(std::size_t nVectors, const std::size_t *sizes, const T *values);
    void initialize(const std::vector<std::vector<T> > &vector2D);
    void initialize(const FlatVector2D<T> &other);

    void destroy();
    void reserve(std::size_t nVectors, std::size_t nItems = 0);
    void swap(FlatVector2D &other) noexcept;
    bool operator==(const FlatVector2D& rhs) const;
    void fill(T &value);
    bool empty() const;

    void clear(bool release = true);
    void clearItems(bool release = true);
    void shrinkToFit();

    const std::size_t * indices() const noexcept;
    const std::size_t * indices(std::size_t i) const noexcept;

    T * data() noexcept;
    const T * data() const noexcept;
    const std::vector<T> & vector() const;

    void pushBack();
    void pushBack(std::size_t subArraySize, const T &value = T());
    void pushBack(const std::vector<T> &subArray);
    void pushBack(std::size_t subArraySize, const T *subArray);
    void pushBackItem(const T& value);
    void pushBackItem(T &&value);
    void pushBackItem(std::size_t i, const T& value);
    void pushBackItem(std::size_t i, T &&value);

    void popBack();
    void popBackItem();
    void popBackItem(std::size_t i);

    void erase(std::size_t i);
    void eraseItem(std::size_t i, std::size_t j);

    void setItem(std::size_t i, std::size_t j, const T &value);
    void setItem(std::size_t i, std::size_t j, T &&value);
    T & getItem(std::size_t i, std::size_t j);
    const T & getItem(std::size_t i, std::size_t j) const;
    const T * get(std::size_t i) const;
    T * get(std::size_t i);

    void rawSetItem(std::size_t k, const T &value);
    void rawSetItem(std::size_t k, T &&value);
    T & rawGetItem(std::size_t k);
    const T & rawGetItem(std::size_t k) const;

    T * back();
    T * first();

    std::size_t size() const;
    std::size_t capacity() const;

    void merge();

    std::size_t getItemCount() const;
    std::size_t getItemCount(std::size_t i) const;
    std::size_t getItemCapacity() const;

    std::size_t getBinarySize() const;

private:
    std::vector<T> m_v;
    std::vector<std::size_t> m_index;

    const T* operator[](std::size_t i) const;
    T* operator[](std::size_t i);
    bool isIndexValid(std::size_t i) const;
    bool isIndexValid(std::size_t i, std::size_t j) const;

    void initialize(std::size_t nVectors, const std::size_t *sizes, std::size_t sizesStride, const T *values, std::size_t valuesStride);

    void destroy(bool destroyIndex, bool destroyValues);

};

}

// Include the implementation
#include "flatVector2D.tpp"

#endif
