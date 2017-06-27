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

#ifndef __BITPIT_PIERCED_KERNEL_ITERATOR_TPP__
#define __BITPIT_PIERCED_KERNEL_ITERATOR_TPP__

namespace bitpit {

/**
* Creates a new uninitialized iterator
*/
template<typename id_t>
PiercedKernelIterator<id_t>::PiercedKernelIterator()
    : m_kernel(nullptr), m_pos(0)
{
}

/**
* Creates a new iterator and initializes it with the position of the const
* base iterator recevied in input.
*/
template<typename id_t>
PiercedKernelIterator<id_t>::PiercedKernelIterator(const kernel_t *kernel, const std::size_t &pos)
    : m_kernel(kernel), m_pos(pos)
{
}

/**
* Exchanges the values of the current iterator and the iterator recevied as
* argument.
*
* \param other the iterator to exchange values with
*/
template<typename id_t>
void PiercedKernelIterator<id_t>::swap(PiercedKernelIterator& other) noexcept
{
    std::swap(m_kernel, other.m_kernel);
    std::swap(m_pos, other.m_pos);
}

/*!
* Get a constant reference of the kernel associated with the iterator.
*
* \result A constant reference of the kernel associated with the iterator.
*/
template<typename id_t>
const typename PiercedKernelIterator<id_t>::kernel_type & PiercedKernelIterator<id_t>::getKernel() const
{
    return *m_kernel;
}

/**
* Gets the id of the current element.
*
* \return The id of the current element or the fallback value if the iterator
* points to an invalid position.
*/
template<typename id_t>
id_t PiercedKernelIterator<id_t>::getId(const id_t &fallback) const noexcept
{
    id_t id;
    if (m_pos >= m_kernel->m_end_pos) {
        id = fallback;
        return id;
    }

    id = m_kernel->m_ids[m_pos];
    if (id >= 0) {
        return id;
    } else {
        id = fallback;
        return id;
    }
}

/**
* Gets the position of the current element.
*
* \return The position of the current element.
*/
template<typename id_t>
std::size_t PiercedKernelIterator<id_t>::getRawIndex() const noexcept
{
    return m_pos;
}

/**
* Pre-increment operator.
*/
template<typename id_t>
PiercedKernelIterator<id_t> & PiercedKernelIterator<id_t>::operator++()
{
    std::size_t delta = 1;
    while (true) {
        m_pos += delta;
        if (m_pos >= m_kernel->m_end_pos) {
            m_pos = m_kernel->m_end_pos;
            return *this;
        }

        id_t id = m_kernel->m_ids[m_pos];
        if (id >= 0) {
            return *this;
        } else {
            delta = - id;
        }
    }

    assert(false);
}

/**
* Post-increment operator.
*/
template<typename id_t>
PiercedKernelIterator<id_t> PiercedKernelIterator<id_t>::operator++(int)
{
    PiercedKernelIterator tmp(m_kernel, m_pos);

    ++(*this);

    return tmp;
}

/**
* Deference operator.
*
* \result A reference to the element currently pointed to by the iterator.
*/
template<typename id_t>
const id_t & PiercedKernelIterator<id_t>::operator*() const
{
    return m_kernel->m_ids[m_pos];
}

/**
* Deference operator.
*
* \result A reference to the element currently pointed to by the iterator.
*/
template<typename id_t>
const id_t * PiercedKernelIterator<id_t>::operator->() const
{
    return m_kernel->m_ids[m_pos];
}

}

#endif
