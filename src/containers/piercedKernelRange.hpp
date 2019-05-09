/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
    void initialize(const kernel_t *kernel);
    void initialize(const kernel_t *kernel, id_t first, id_t last);
    void initialize(const const_iterator &begin, const const_iterator &end);

    void swap(PiercedKernelRange &other) noexcept;
    std::size_t evalSize() const;

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
        if (m_cbegin == rhs.m_cbegin) {
            return false;
        }

        if (m_cend == rhs.m_cend) {
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
        if (m_cbegin == rhs.m_cbegin) {
            return true;
        }

        if (m_cend != rhs.m_cend) {
            return true;
        }

        return false;
    }

protected:
    /*! Begin */
    const_iterator m_cbegin;

    /*! End */
    const_iterator m_cend;

};

}

// Include the implementation
#include "piercedKernelRange.tpp"

#endif
