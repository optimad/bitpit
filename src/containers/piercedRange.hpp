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

#ifndef __BITPIT_PIERCED_RANGE_HPP__
#define __BITPIT_PIERCED_RANGE_HPP__

#include "piercedVector.hpp"

namespace bitpit {

/*!
    @brief The PiercedRange allow to iterate using range-based loops over
    a PiercedVector.
*/
template<typename value_t, typename id_t = long,
         typename value_no_cv_t = typename std::remove_cv<value_t>::type>
class PiercedRange
{

private:
    /*!
        Container.
    */
    template<typename PV_value_t, typename PV_id_t>
    using Container = PiercedVector<PV_value_t, PV_id_t>;

    /*
        Container type

        When building a const_range the iterator has to be declared const.
    */
    typedef
        typename std::conditional<std::is_const<value_t>::value,
            const Container<value_no_cv_t, id_t>,
            Container<value_no_cv_t, id_t>
        >::type

        container_t;

    /*
        Iterator type

        When building a const_iterator the pointer to the container has to
        be declared const.
    */
    typedef
        typename std::conditional<std::is_const<value_t>::value,
            typename container_t::const_iterator,
            typename container_t::iterator
        >::type

        iterator_t;

    /*
        Const iterator type

        When building a const_iterator the pointer to the container has to
        be declared const.
    */
    typedef typename container_t::const_iterator const_iterator_t;

public:
    /*! Type of container */
    typedef container_t container_type;

    /*! Type of data stored in the container */
    typedef value_t value_type;

    /*! Type of ids stored in the container */
    typedef id_t id_type;

    /*! Type of iterator */
    typedef iterator_t iterator;

    /*! Type of constant iterator */
    typedef const_iterator_t const_iterator;

    // Constructors
    PiercedRange();
    PiercedRange(container_t *container);
    PiercedRange(container_t *container, id_t first, id_t last);
    PiercedRange(iterator begin, iterator end);

    // General methods
    void swap(PiercedRange &other) noexcept;

    // Methods to get begin and end
    template<typename U = value_t, typename U_no_cv = value_no_cv_t,
             typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type = 0>
    iterator begin() noexcept;

    template<typename U = value_t, typename U_no_cv = value_no_cv_t,
             typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type = 0>
    iterator end() noexcept;

    const_iterator begin() const noexcept;
    const_iterator end() const noexcept;

    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;

    /*!
        Two-way comparison.
    */
    template<typename other_value_t, typename other_id_t = long>
    bool operator==(const PiercedRange<other_value_t, other_id_t> &rhs) const
    {
        if (m_container == rhs.m_container) {
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
        Two-way comparison.
    */
    template<typename other_value_t, typename other_id_t = long>
    bool operator!=(const PiercedRange<other_value_t, other_id_t> &rhs) const
    {
        if (m_container != rhs.m_container) {
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

private:
    /*! Container */
    container_t *m_container;

    /*! Begin */
    size_t m_begin_pos;

    /*! End */
    size_t m_end_pos;

};

}

// Include the implementation
#include "piercedRange.tpp"

#endif
