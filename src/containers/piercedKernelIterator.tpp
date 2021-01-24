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
PiercedKernelIterator<id_t>::PiercedKernelIterator(const kernel_t *kernel, std::size_t pos)
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
    return getPos();
}

/**
* Gets the position of the current element.
*
* \return The position of the current element.
*/
template<typename id_t>
std::size_t PiercedKernelIterator<id_t>::getPos() const noexcept
{
    return m_pos;
}

/**
* Pre-increment operator.
*
* Increment operator does not check whether it crosses the end() of the
* container. Thus, calling this function if the iterator is already at
* the end of the container results in undefined behavior.
*/
template<typename id_t>
PiercedKernelIterator<id_t> & PiercedKernelIterator<id_t>::operator++()
{
    if (m_pos >= m_kernel->m_end_pos) {
        return *this;
    }

    auto basePosItr    = m_kernel->m_ids.begin();
    auto endPosItr     = basePosItr + m_kernel->m_end_pos;
    auto currentPosItr = basePosItr + m_pos + 1;

    while (currentPosItr < endPosItr) {
        id_t id = *currentPosItr;
        if (id >= 0) {
            break;
        }

        currentPosItr -= id;
    }

    m_pos = currentPosItr - basePosItr;

    return *this;
}

/**
* Post-increment operator.
*
* Increment operator does not check whether it crosses the end() of the
* container. Thus, calling this function if the iterator is already at
* the end of the container results in undefined behavior.
*/
template<typename id_t>
PiercedKernelIterator<id_t> PiercedKernelIterator<id_t>::operator++(int)
{
    PiercedKernelIterator tmp(m_kernel, m_pos);

    ++(*this);

    return tmp;
}

/**
* Pre-decrement operator.
*
* Decrement operator does not check whether it crosses the begin() of the
* container. Thus, calling this function if the iterator is already at
* the begin of the container results in undefined behavior.
*/
template<typename id_t>
PiercedKernelIterator<id_t> & PiercedKernelIterator<id_t>::operator--()
{
    if (m_pos <= m_kernel->m_begin_pos) {
        return *this;
    }

    auto basePosItr    = m_kernel->m_ids.begin();
    auto beginPosItr   = basePosItr + m_kernel->m_begin_pos;
    auto currentPosItr = basePosItr + m_pos - 1;

    while (currentPosItr > beginPosItr) {
        id_t id = *currentPosItr;
        if (id >= 0) {
            break;
        }

        --currentPosItr;
    }

    m_pos = currentPosItr - basePosItr;

    return *this;
}

/**
* Post-decrement operator.
*
* Decrement operator does not check whether it crosses the begin() of the
* container. Thus, calling this function if the iterator is already at
* the begin of the container results in undefined behavior.
*/
template<typename id_t>
PiercedKernelIterator<id_t> PiercedKernelIterator<id_t>::operator--(int)
{
    PiercedKernelIterator tmp(m_kernel, m_pos);

    --(*this);

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
