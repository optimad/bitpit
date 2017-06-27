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

#ifndef __BITPIT_PIERCED_STORAGE_ITERATOR_HPP__
#define __BITPIT_PIERCED_STORAGE_ITERATOR_HPP__

#include <cassert>
#include <iterator>
#include <limits>
#include <type_traits>

#include "piercedKernel.hpp"

#define  __PSI_REFERENCE__ typename PiercedStorageIterator<value_t, id_t, value_no_cv_t>::reference
#define  __PSI_POINTER__   typename PiercedStorageIterator<value_t, id_t, value_no_cv_t>::pointer

namespace bitpit{

template<typename value_t, typename id_t>
class PiercedStorage;

/**
* \ingroup containers
*
* \brief Iterator for the class PiercedStorage
*
* \tparam value_t is the type of elements in the storage
* \tparam id_t is the type of ids associated to the elements
*/
template<typename value_t, typename id_t = long,
         typename value_no_cv_t = typename std::remove_cv<value_t>::type>
class PiercedStorageIterator
    : public std::iterator<std::forward_iterator_tag, value_no_cv_t, std::ptrdiff_t, value_t*, value_t&>
{

friend class PiercedStorageIterator<value_no_cv_t, id_t, value_no_cv_t>;

template<typename PS_value_t, typename PS_id_t>
friend class PiercedStorage;

private:
    /**
    * Storage.
    */
    template<typename PS_value_t, typename PS_id_t>
    using Storage = PiercedStorage<PS_value_t, PS_id_t>;

    /**
    * Storage type
    *
    * When building a const_iterator the pointer to the storage has to be
    * declared const.
    */
    typedef
        typename std::conditional<std::is_const<value_t>::value,
            const Storage<value_no_cv_t, id_t>,
            Storage<value_no_cv_t, id_t>
        >::type

        storage_t;

    /**
    * Kernel type
    */
    typedef typename storage_t::kernel_t kernel_t;

public:
    /*!
    *Constant container
    */
    typedef storage_t storage_type;

    /**
    * Type of values in the storage
    */
    typedef value_t value_type;

    /**
    * Type of ids in the kernel
    */
    typedef id_t id_type;

    // Constructors
    PiercedStorageIterator();

    // General methods
    void swap(PiercedStorageIterator& other) noexcept;

    storage_type & getStorage() const;

    // Methos to extract information on the current element
    id_t getId(const id_t &fallback = -1) const noexcept;
    std::size_t getRawIndex() const noexcept;
    __PSI_REFERENCE__ getValue(std::size_t k = 0) const;

    // Operators
    PiercedStorageIterator& operator++();
    PiercedStorageIterator operator++(int);

    __PSI_REFERENCE__ operator*() const;
    __PSI_POINTER__ operator->() const;

    template<typename U = value_t, typename U_no_cv = value_no_cv_t,
             typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type = 0>
    operator PiercedStorageIterator<const U_no_cv, id_t>() const;

    /**
    * Two-way comparison.
    */
    template<typename other_value_t, typename other_id_t = long,
         typename other_value_no_cv_t = typename std::remove_cv<value_t>::type>
    bool operator==(const PiercedStorageIterator<other_value_t, other_id_t, other_value_no_cv_t>& rhs) const
    {
        return (m_storage == rhs.m_storage) && (m_pos == rhs.m_pos);
    }

    /**
    * Two-way comparison.
    */
    template<typename other_value_t, typename other_id_t = long,
         typename other_value_no_cv_t = typename std::remove_cv<value_t>::type>
    bool operator!=(const PiercedStorageIterator<other_value_t, other_id_t, other_value_no_cv_t>& rhs) const
    {
        return (m_storage != rhs.m_storage) || (m_pos != rhs.m_pos);
    }

private:
    /**
    * Internal pointer to the storage.
    */
    storage_t *m_storage;

    /**
    * Internal pointer to the kernel.
    */
    const kernel_t *m_kernel;

    /**
    * Position inside the kernel.
    */
    std::size_t m_pos;

    // Constructors
    explicit PiercedStorageIterator(storage_t *storage, const std::size_t &pos);

};

}

// Include the implementation
#include "piercedStorageIterator.tpp"

#endif
