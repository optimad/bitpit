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

#ifndef __BITPIT_PROXY_VECTOR_HPP__
#define __BITPIT_PROXY_VECTOR_HPP__

#define  __PXI_REFERENCE__ typename ProxyVectorIterator<T, T_no_cv>::reference
#define  __PXI_POINTER__   typename ProxyVectorIterator<T, T_no_cv>::pointer

#define __PXV_REFERENCE__       typename ProxyVector<T>::reference
#define __PXV_CONST_REFERENCE__ typename ProxyVector<T>::const_reference
#define __PXV_POINTER__         typename ProxyVector<T>::pointer
#define __PXV_CONST_POINTER__   typename ProxyVector<T>::const_pointer
#define __PXV_ITERATOR__        typename ProxyVector<T>::iterator
#define __PXV_CONST_ITERATOR__  typename ProxyVector<T>::const_iterator

#include <memory>
#include <vector>

namespace bitpit{

template<typename PXV_T>
class ProxyVector;

/*!
    @ingroup containers
    @brief Iterator for the class ProxyVector

    @tparam T is the type of the objects handled by the ProxyVector
*/
template<typename T, typename T_no_cv = typename std::remove_cv<T>::type>
class ProxyVectorIterator
    : public std::iterator<std::random_access_iterator_tag, T_no_cv, std::ptrdiff_t, T*, T&>
{

template<typename PXV_T>
friend class ProxyVector;

public:
    /*!
        Type of data stored in the container
    */
    typedef T value_type;

    // Constructors
    ProxyVectorIterator();

    // General methods
    void swap(ProxyVectorIterator& other) noexcept;

    // Operators
    ProxyVectorIterator& operator++();
    ProxyVectorIterator operator++(int);

    ProxyVectorIterator& operator--();
    ProxyVectorIterator operator--(int);

    ProxyVectorIterator& operator+=(int increment);

    std::size_t operator-(const ProxyVectorIterator &other);

    __PXI_REFERENCE__ operator*() const;
    __PXI_POINTER__ operator->() const;

    template<typename U = T, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    operator ProxyVectorIterator<const U>() const;

    /*!
        Two-way comparison.
    */
    template<typename U = T, typename U_no_cv = typename std::remove_cv<U>::type>
    bool operator==(const ProxyVectorIterator<U, U_no_cv>& rhs) const
    {
        return (m_position == rhs.m_position);
    }

    /*!
        Two-way comparison.
    */
    template<typename U = T, typename U_no_cv = typename std::remove_cv<U>::type>
    bool operator!=(const ProxyVectorIterator<U, U_no_cv>& rhs) const
    {
        return (m_position != rhs.m_position);
    }

private:
    /*!
        Position inside the container.
    */
    T *m_position;

    // Constructors
    explicit ProxyVectorIterator(T *position);

};

/*!
    @ingroup containers
    @brief Metafunction for generating a list of elements that can be either
    stored in an external vectror or, if the elements are constant, inside
    the container itself.

    @details
    Usage: Use <tt>ProxyVector<Type></tt> to declare a list of elements that
    can be either stored in an external vectror or, if the elements are
    constant, inside the container itself. When the ProxyVector is destroyed,
    the elements of the list will be destroyed only if owned by the ProxyVector o
    bject itself.

    @tparam T is the type of the objects handled by the ProxyVector
*/
template<typename T>
class ProxyVector
{

private:
    typedef typename std::remove_cv<T>::type T_no_cv;

public:
    /*!
        Type of data stored in the container
    */
    typedef T value_type;

    /*!
        Iterator for the pierced array raw container.
    */
    typedef ProxyVectorIterator<T> iterator;

    /*!
        Constant iterator for the pierced array raw container.
    */
    typedef ProxyVectorIterator<const T_no_cv> const_iterator;

    /*!
        Reference
    */
    typedef
        typename std::conditional<std::is_const<T>::value,
            typename std::vector<T_no_cv>::const_reference,
            typename std::vector<T_no_cv>::reference>::type
        reference;

    /*!
        Constant reference
    */
    typedef typename std::vector<T_no_cv>::const_reference const_reference;

    /*!
        Pointer
    */
    typedef
        typename std::conditional<std::is_const<T>::value,
            typename std::vector<T_no_cv>::const_pointer,
            typename std::vector<T_no_cv>::pointer>::type
        pointer;

    /*!
        Constant pointer
    */
    typedef typename std::vector<T_no_cv>::const_pointer const_pointer;

    ProxyVector();
    ProxyVector(T *data, std::size_t size);
    template<typename U = T, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    ProxyVector(std::vector<T_no_cv> &&data);

    ProxyVector(const ProxyVector &other);
    ProxyVector(ProxyVector &&other) = default;

    /*!
        Move assignment operator.

        The move assignment operator "steals" the resources held by the
        argument.
    */
    ProxyVector & operator=(ProxyVector &&other) = default;

    ProxyVector & operator=(const ProxyVector &other);

    void set(T *data, std::size_t size);
    template<typename U = T, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    void set(std::vector<T_no_cv> &&storage);

    void clear();
    void swap(ProxyVector &other);

    bool empty() const;
    std::size_t size() const;
    bool operator==(const ProxyVector& rhs) const;

    __PXV_CONST_POINTER__ data() const noexcept;
    __PXV_POINTER__ data() noexcept;

    __PXV_CONST_REFERENCE__ operator[](std::size_t n) const;
    __PXV_REFERENCE__ operator[](std::size_t n);

    __PXV_CONST_REFERENCE__ at(std::size_t n) const;
    __PXV_REFERENCE__ at(std::size_t n);

    __PXV_CONST_REFERENCE__ front() const;
    __PXV_REFERENCE__ front();

    __PXV_CONST_REFERENCE__ back() const;
    __PXV_REFERENCE__ back();

    __PXV_ITERATOR__ begin();
    __PXV_ITERATOR__ end();

    __PXV_CONST_ITERATOR__ begin() const;
    __PXV_CONST_ITERATOR__ end() const;

    __PXV_CONST_ITERATOR__ cbegin();
    __PXV_CONST_ITERATOR__ cend();

private:
    std::unique_ptr<std::vector<T_no_cv>> m_storage;

    T *m_data;
    std::size_t m_size;

};

// Constant proxy vector
template<typename T>
using ConstProxyVector = ProxyVector<const T>;

}

// Include the implementation
#include "proxyVector.tpp"

// Some commonly used ProxyVectors are instantiated explicitly
namespace bitpit{

extern template class ProxyVector<int>;
extern template class ProxyVector<long>;
extern template class ProxyVector<double>;

}


#endif
