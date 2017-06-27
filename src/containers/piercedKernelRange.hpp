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

#ifndef __BITPIT_PIERCED_KERNEL_RANGE_HPP__
#define __BITPIT_PIERCED_KERNEL_RANGE_HPP__

namespace bitpit {

/*!
    @brief The PiercedKernelRange allow to iterate using range-based loops over
    a PiercedStorage.
*/
template<typename id_t = long>
class PiercedKernelRange
{

template<typename PK_id_t>
friend class PiercedKernel;

private:
    /**
    * Kernel.
    */
    template<typename PK_id_t>
    using Kernel = PiercedKernel<PK_id_t>;

    /**
    * Kernel type
    */
    typedef Kernel<id_t> kernel_t;

    /*
    * Iterator type
    */
    typedef typename kernel_t::const_iterator const_iterator_t;

public:
    /*! Kernel type */
    typedef kernel_t kernel_type;

    /*! Type of ids stored in the container */
    typedef id_t id_type;

    /*! Type of constant iterator */
    typedef const_iterator_t const_iterator;

    // Constructors
    PiercedKernelRange();
    PiercedKernelRange(const kernel_t *kernel);
    PiercedKernelRange(const kernel_t *kernel, id_t first, id_t last);
    PiercedKernelRange(const const_iterator &begin, const const_iterator &end);

    // General methods
    void swap(PiercedKernelRange &other) noexcept;

    const kernel_type & getKernel() const;

    // Methods to get begin and end
    const_iterator begin() const noexcept;
    const_iterator end() const noexcept;

    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;

    /*!
        Two-way comparison.
    */
    template<typename other_id_t = long>
    bool operator==(const PiercedKernelRange<other_id_t> &rhs) const
    {
        if (m_kernel == rhs.m_kernel) {
            return false;
        }

        if (m_begin_pos == rhs.m_begin_pos) {
            return false;
        }

        if (m_end_pos == rhs.m_end_pos) {
            return false;
        }

        return true;
    }

    /*!
    * Two-way comparison.
    */
    template<typename other_id_t = long>
    bool operator!=(const PiercedKernelRange<other_id_t> &rhs) const
    {
        if (m_kernel != rhs.m_kernel) {
            return true;
        }

        if (m_begin_pos == rhs.m_begin_pos) {
            return true;
        }

        if (m_end_pos != rhs.m_end_pos) {
            return true;
        }

        return false;
    }

protected:
    /*! Container */
    const kernel_t *m_kernel;

    /*! Begin */
    size_t m_begin_pos;

    /*! End */
    size_t m_end_pos;

};

}

// Include the implementation
#include "piercedKernelRange.tpp"

#endif
