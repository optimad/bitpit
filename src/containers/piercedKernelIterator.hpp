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

#ifndef __BITPIT_PIERCED_KERNEL_ITERATOR_HPP__
#define __BITPIT_PIERCED_KERNEL_ITERATOR_HPP__

#include <cassert>
#include <iterator>
#include <limits>
#include <type_traits>

namespace bitpit{

template<typename id_t>
class PiercedKernel;

/**
* \ingroup containers
*
* \brief Iterator for the class PiercedKernel
*
* \tparam value_t is the type of elements in the storage
* \tparam id_t is the type of ids associated to the elements
*/
template<typename id_t = long>
class PiercedKernelIterator
{

template<typename PK_id_t>
friend class PiercedKernel;

template<typename PSI_value_t, typename PSI_id_t, typename PSI_value_no_cv_t>
friend class PiercedStorageIterator;

private:
    /**
    * Kernel.
    */
    template<typename PK_id_t>
    using Kernel = PiercedKernel<PK_id_t>;

    /**
    * Kernel type
    *
    * When building a const_iterator the pointer to the storage has to be
    * declared const.
    */
    typedef Kernel<id_t> kernel_t;

public:
    /*!
    * Iterator category
    */
    typedef std::forward_iterator_tag iterator_category;

    /*!
    * Value type
    */
    typedef id_t value_type;

    /*!
    * Difference type
    */
    typedef std::ptrdiff_t difference_type;

    /*!
    * Pointer type
    */
    typedef const id_t * pointer;

    /*!
    * Reference type
    */
    typedef const id_t & reference;

    /**
    * Kernel type
    */
    typedef kernel_t kernel_type;

    /**
    * Type of ids in the kernel
    */
    typedef id_t id_type;

    // Constructors
    PiercedKernelIterator();

    // General methods
    void swap(PiercedKernelIterator& other) noexcept;

    const kernel_type & getKernel() const;

    // Methos to extract information on the current element
    id_t getId(const id_t &fallback = -1) const noexcept;
    std::size_t getRawIndex() const noexcept;

    // Operators
    PiercedKernelIterator& operator++();
    PiercedKernelIterator operator++(int);

    const id_t & operator*() const;
    const id_t * operator->() const;

    /**
    * Two-way comparison.
    */
    template<typename other_id_t>
    bool operator==(const PiercedKernelIterator<other_id_t>& rhs) const
    {
        return (m_kernel == rhs.m_kernel) && (m_pos == rhs.m_pos);
    }

    /**
    * Two-way comparison.
    */
    template<typename other_id_t>
    bool operator!=(const PiercedKernelIterator<other_id_t>& rhs) const
    {
        return (m_kernel != rhs.m_kernel) || (m_pos != rhs.m_pos);
    }

private:
    /**
    * Internal pointer to the kernel.
    */
    const kernel_t *m_kernel;

    /**
    * Position inside the kernel.
    */
    std::size_t m_pos;

    // Constructors
    explicit PiercedKernelIterator(const kernel_t *kernel, const std::size_t &pos);

};

}

// Include the implementation
#include "piercedKernelIterator.tpp"

#endif
