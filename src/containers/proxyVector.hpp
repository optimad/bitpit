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

#ifndef __BITPIT_PROXY_VECTOR_HPP__
#define __BITPIT_PROXY_VECTOR_HPP__

#define  __PXI_POINTER_TYPE__ \
    typename std::conditional<std::is_const<value_t>::value, \
        typename std::vector<value_no_cv_t>::const_pointer, \
        typename std::vector<value_no_cv_t>::pointer>::type

#define  __PXI_REFERENCE_TYPE__ \
    typename std::conditional<std::is_const<value_t>::value, \
        typename std::vector<value_no_cv_t>::const_reference, \
        typename std::vector<value_no_cv_t>::reference>::type

#define  __PXI_REFERENCE__ typename ProxyVectorIterator<value_t, value_no_cv_t>::reference
#define  __PXI_POINTER__   typename ProxyVectorIterator<value_t, value_no_cv_t>::pointer

#define __PXV_REFERENCE__             typename ProxyVector<value_t>::reference
#define __PXV_CONST_REFERENCE__       typename ProxyVector<value_t>::const_reference
#define __PXV_POINTER__               typename ProxyVector<value_t>::pointer
#define __PXV_CONST_POINTER__         typename ProxyVector<value_t>::const_pointer
#define __PXV_STORAGE_POINTER__       typename ProxyVector<value_t>::storage_pointer
#define __PXV_STORAGE_CONST_POINTER__ typename ProxyVector<value_t>::storage_const_pointer
#define __PXV_ITERATOR__              typename ProxyVector<value_t>::iterator
#define __PXV_CONST_ITERATOR__        typename ProxyVector<value_t>::const_iterator

#include <cassert>
#include <memory>
#include <vector>

#include "bitpit_common.hpp"

namespace bitpit {

template<typename PXV_value_t>
class ProxyVector;

/*!
    @ingroup containers
    @brief Iterator for the class ProxyVector

    @tparam value_t is the type of the objects handled by the ProxyVector
    @tparam pointer_t defines a pointer to the type iterated over
    @tparam reference_t defines a reference to the type iterated over
*/
template<typename value_t, typename value_no_cv_t = typename std::remove_cv<value_t>::type>
class ProxyVectorIterator
    : public std::iterator<std::random_access_iterator_tag, value_t, std::ptrdiff_t, __PXI_POINTER_TYPE__, __PXI_REFERENCE_TYPE__>
{

template<typename PXV_value_t>
friend class ProxyVector;

friend class ProxyVectorIterator<typename std::add_const<value_t>::type, value_no_cv_t>;

public:
    /*!
        Iterator category
    */
    typedef std::bidirectional_iterator_tag iterator_category;

    /*!
        Value type
    */
    typedef value_t value_type;

    /*!
        Difference type
    */
    typedef std::ptrdiff_t difference_type;

    /*!
        Pointer type
    */
    typedef __PXI_POINTER_TYPE__ pointer;

    /*!
        Reference type
    */
    typedef __PXI_REFERENCE_TYPE__ reference;

    // Constructors
    ProxyVectorIterator();

    template<typename other_value_t, typename std::enable_if<std::is_const<value_t>::value && !std::is_const<other_value_t>::value && std::is_same<other_value_t, typename std::remove_cv<value_t>::type>::value, int>::type = 0>
    ProxyVectorIterator(const ProxyVectorIterator<other_value_t> &other);

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

    template<typename other_value_t, typename std::enable_if<std::is_const<value_t>::value && !std::is_const<other_value_t>::value && std::is_same<other_value_t, typename std::remove_cv<value_t>::type>::value, int>::type = 0>
    ProxyVectorIterator & operator=(const ProxyVectorIterator<other_value_t> &other);

    /*!
        Two-way comparison.
    */
    bool operator==(const ProxyVectorIterator &other) const
    {
        return (m_position == other.m_position);
    }

    /*!
        Two-way comparison.
    */
    bool operator!=(const ProxyVectorIterator &other) const
    {
        return (m_position != other.m_position);
    }

private:
    /*!
        Position inside the container.
    */
    __PXI_POINTER__ m_position;

    // Constructors
    explicit ProxyVectorIterator(__PXI_POINTER__ position);

};

/*!
    @ingroup containers
    @brief Interface for ProxyVector storages.

    @tparam pointer_t defines a pointer to the data stored
    @tparam const_pointer_t defines a constant pointer to the data stored
*/
template<typename pointer_t, typename const_pointer_t>
class ProxyVectorStorageInterface
{

public:
    virtual ~ProxyVectorStorageInterface() = default;

    virtual bool empty() const = 0;
    virtual std::size_t size() const = 0;

    virtual pointer_t data() = 0;
    virtual const_pointer_t data() const = 0;

    virtual void resize(std::size_t size) = 0;

};

/*!
    @ingroup containers
    @brief Metafunction for generating ProxyVector dummy storages.

    @tparam value_t is the type of the objects handled by the dummy storage
    @tparam pointer_t defines a pointer to the data stored
    @tparam const_pointer_t defines a constant pointer to the data stored
*/
template<typename value_t, typename pointer_t = value_t *, typename const_pointer_t = const value_t *>
class ProxyVectorDummyStorage : public ProxyVectorStorageInterface<pointer_t, const_pointer_t>
{

template<typename PXV_value_t>
friend class ProxyVector;

public:
    typedef pointer_t pointer;
    typedef const_pointer_t const_pointer;

    void swap(ProxyVectorDummyStorage &other) noexcept;

    pointer data() override;
    const_pointer data() const override;

    bool empty() const override;
    std::size_t size() const override;

    void resize(std::size_t size) override;

protected:
    ProxyVectorDummyStorage(std::size_t size = 0);

};

/*!
    @ingroup containers
    @brief Metafunction for generating ProxyVector storages.

    @tparam value_t is the type of the objects handled by the storage
    @tparam pointer_t defines a pointer to the data stored
    @tparam const_pointer_t defines a constant pointer to the data stored
*/
template<typename value_t, typename pointer_t = typename std::vector<value_t>::pointer, typename const_pointer_t = typename std::vector<value_t>::const_pointer>
class ProxyVectorStorage : public ProxyVectorStorageInterface<pointer_t, const_pointer_t>
{

template<typename PXV_value_t>
friend class ProxyVector;

public:
    typedef pointer_t pointer;
    typedef const_pointer_t const_pointer;

    ~ProxyVectorStorage();

    void swap(ProxyVectorStorage &other) noexcept;

    pointer data() override;
    const_pointer data() const override;

    bool empty() const override;
    std::size_t size() const override;

    void resize(std::size_t size) override;

protected:
    ProxyVectorStorage(std::size_t size = 0);
    ProxyVectorStorage(const ProxyVectorStorage &other);
    ProxyVectorStorage(ProxyVectorStorage &&other) = default;

private:
    static const int MEMORY_POOL_VECTOR_COUNT = 10;
    static const int MEMORY_POOL_MAX_CAPACITY = 128;

    static std::vector<std::unique_ptr<std::vector<value_t>>> m_storagePool;

    std::unique_ptr<std::vector<value_t>> createStorage(std::size_t size);
    std::unique_ptr<std::vector<value_t>> createStorage(const std::unique_ptr<std::vector<value_t>> &source);
    void destroyStorage(std::unique_ptr<std::vector<value_t>> *storage);

    std::unique_ptr<std::vector<value_t>> m_storage;

};

/*!
    @ingroup containers
    @brief Metafunction for generating a list of elements that can be either
    stored in an external vectror or, if the elements are constant, inside
    the container itself.

    @details
    Usage: Use <tt>ProxyVector<value_t></tt> to declare a list of elements that
    can be either stored in an external vectror or, if the elements are
    constant, inside the container itself. When the ProxyVector is destroyed,
    the elements of the list will be destroyed only if owned by the ProxyVector o
    bject itself.

    @tparam value_t is the type of the objects handled by the ProxyVector
*/
template<typename value_t>
class ProxyVector
{

private:
    typedef typename std::remove_cv<value_t>::type value_no_cv_t;

public:
    /*!
        Type of data stored in the container
    */
    typedef value_t value_type;

    /*!
        Iterator for the container
    */
    typedef ProxyVectorIterator<value_t, value_no_cv_t> iterator;

    /*!
        Constant iterator for the container
    */
    typedef ProxyVectorIterator<typename std::add_const<value_no_cv_t>::type, value_no_cv_t> const_iterator;

    /*!
        Pointer type
    */
    typedef
        typename std::conditional<std::is_const<value_t>::value,
            typename std::vector<value_no_cv_t>::const_pointer,
            typename std::vector<value_no_cv_t>::pointer>::type
        pointer;

    /*!
        Constant pointer
    */
    typedef typename std::vector<value_no_cv_t>::const_pointer const_pointer;

    /*!
        Storage pointer type
    */
    typedef typename std::vector<value_no_cv_t>::pointer storage_pointer;

    /*!
        Constant storage pointer
    */
    typedef typename std::vector<value_no_cv_t>::const_pointer storage_const_pointer;

    /*!
        Reference type
    */
    typedef
        typename std::conditional<std::is_const<value_t>::value,
            typename std::vector<value_no_cv_t>::const_reference,
            typename std::vector<value_no_cv_t>::reference>::type
        reference;

    /*!
        Constant reference
    */
    typedef typename std::vector<value_no_cv_t>::const_reference const_reference;

    /*!
        Flag to use the internal storage
    */
    static constexpr __PXV_POINTER__ INTERNAL_STORAGE = nullptr;

    ProxyVector();
    template<typename U = value_t, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    ProxyVector(std::size_t size);
    template<typename U = value_t, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    ProxyVector(__PXV_POINTER__ data, std::size_t size);
    template<typename U = value_t, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    ProxyVector(std::size_t size, std::size_t capacity);
    template<typename U = value_t, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    ProxyVector(__PXV_POINTER__ data, std::size_t size, std::size_t capacity);
    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    ProxyVector(__PXV_POINTER__ data, std::size_t size);

    ProxyVector(const ProxyVector &other);
    ProxyVector(ProxyVector &&other) = default;

    ProxyVector & operator=(ProxyVector &&other) = default;
    ProxyVector & operator=(const ProxyVector &other);

    void swap(ProxyVector &other);

    bool empty() const;
    std::size_t size() const;

    bool operator==(const ProxyVector &other) const;

    template<typename U = value_t, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    void set(__PXV_POINTER__ data, std::size_t size);
    template<typename U = value_t, typename std::enable_if<std::is_const<U>::value, int>::type = 0>
    void set(__PXV_POINTER__ data, std::size_t size, std::size_t capacity);
    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    void set(__PXV_POINTER__ data, std::size_t size);

    __PXV_STORAGE_POINTER__ storedData() noexcept;
    __PXV_STORAGE_CONST_POINTER__ storedData() const noexcept;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_POINTER__ data() noexcept;
    __PXV_CONST_POINTER__ data() const noexcept;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_REFERENCE__ operator[](std::size_t n);
    __PXV_CONST_REFERENCE__ operator[](std::size_t n) const;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_REFERENCE__ at(std::size_t n);
    __PXV_CONST_REFERENCE__ at(std::size_t n) const;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_REFERENCE__ front();
    __PXV_CONST_REFERENCE__ front() const;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_REFERENCE__ back();
    __PXV_CONST_REFERENCE__ back() const;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_ITERATOR__ begin();
    __PXV_CONST_ITERATOR__ begin() const;

    template<typename U = value_t, typename std::enable_if<!std::is_const<U>::value, int>::type = 0>
    __PXV_ITERATOR__ end();
    __PXV_CONST_ITERATOR__ end() const;

    __PXV_CONST_ITERATOR__ cbegin();
    __PXV_CONST_ITERATOR__ cend();

private:
    /*!
        Storage type
    */
    typedef
        typename std::conditional<std::is_const<value_t>::value,
            ProxyVectorStorage<value_no_cv_t>,
            ProxyVectorDummyStorage<value_no_cv_t>>::type
        storage_t;

    storage_t m_storage;

    std::size_t m_size;
    __PXV_POINTER__ m_data;

};

// Constant proxy vector
template<typename value_t>
using ConstProxyVector = ProxyVector<const value_t>;

}

// Include the implementation
#include "proxyVector.tpp"

namespace bitpit{

// Some commonly used ProxyVectors are instantiated explicitly
extern template class ProxyVector<int>;
extern template class ProxyVector<long>;
extern template class ProxyVector<double>;

extern template class ProxyVector<const int>;
extern template class ProxyVector<const long>;
extern template class ProxyVector<const double>;

}

#endif
