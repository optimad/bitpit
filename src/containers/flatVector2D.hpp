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

template<class T>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::FlatVector2D<T> &vector);

template<class T>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::FlatVector2D<T> &vector);

namespace bitpit{

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
friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream &buffer, const FlatVector2D<U> &vector);
template<class U>
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream &buffer, FlatVector2D<U> &vector);

public:
    FlatVector2D(bool initialize = true);
    FlatVector2D(const std::vector<int> &sizes, const T &value = T());
    FlatVector2D(const int &nVectors, const int &size, const T &value = T());
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

    void initialize(const std::vector<int> &sizes, const T &value = T());
    void initialize(const int &nVectors, const int &size, const T &value = T());
    void initialize(const std::vector<std::vector<T> > &vector2D);

    void destroy();
    void reserve(int nVectors, int nItems = 0);
    void swap(FlatVector2D &other) noexcept;
    bool operator==(const FlatVector2D& rhs) const;
    void fill(T &value);
    bool empty() const;

    void clear();
    void shrinkToFit();

    T * data() noexcept;
    const T * data() const noexcept;
    const std::vector<T> & vector() const;

    void pushBack();
    void pushBack(const int &subArraySize, const T &value = T());
    void pushBack(const std::vector<T> &subArray);
    void pushBack(int subArraySize, const T *subArray);
    void pushBackItem(const T& value);
    void pushBackItem(const int &i, const T& value);

    void popBack();
    void popBackItem();
    void popBackItem(const int &i);

    void erase(const int &i);
    void eraseItem(const int &i, const int &j);

    void setItem(const int &i, const int &j, const T &value);
    T & getItem(const int &i, const int &j);
    const T & getItem(const int &i, const int &j) const;
    const T * get(const int &i) const;
    T * get(const int &i);

    void rawSetItem(const int &k, const T &value);
    T & rawGetItem(const int &k);
    const T & rawGetItem(const int &k) const;

    T * back();
    T * first();

    int size() const;
    int capacity() const;

    void merge();

    int getItemCount() const;
    int getItemCount(int i) const;
    int getItemCapacity() const;

    size_t getBinarySize() const;

private:
    std::vector<T> m_v;
    std::vector<std::size_t> m_index;

    const T* operator[](const int &i) const;
    T* operator[](const int &i);
    bool isIndexValid(const int &i) const;
    bool isIndexValid(const int &i, const int &j) const;

    void destroy(bool destroyIndex, bool destroyValues);

};

}

// Include the implementation
#include "flatVector2D.tpp"

#endif
