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

#ifndef __BITPIT_PIERCED_STORAGE_ITERATOR_TPP__
#define __BITPIT_PIERCED_STORAGE_ITERATOR_TPP__

namespace bitpit {

/**
* Creates a new uninitialized iterator
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageIterator<value_t, id_t, value_no_cv_t>::PiercedStorageIterator()
    : m_storage(nullptr), m_kernel(nullptr), m_pos(0)
{
}

/**
* Creates a new iterator and initializes it with the position of the const
* base iterator recevied in input.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageIterator<value_t, id_t, value_no_cv_t>::PiercedStorageIterator(storage_t *storage, const std::size_t &pos)
    : m_storage(storage), m_kernel(&(m_storage->getKernel())), m_pos(pos)
{
}

/**
* Exchanges the values of the current iterator and the iterator recevied as
* argument.
*
* \param other the iterator to exchange values with
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
void PiercedStorageIterator<value_t, id_t, value_no_cv_t>::swap(PiercedStorageIterator& other) noexcept
{
    std::swap(m_storage, other.m_storage);
    std::swap(m_kernel, other.m_kernel);
    std::swap(m_pos, other.m_pos);
}

/*!
* Get a constant reference of the storage associated with the iterator.
*
* \result A constant reference of the storage associated with the iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageIterator<value_t, id_t, value_no_cv_t>::storage_type & PiercedStorageIterator<value_t, id_t, value_no_cv_t>::getStorage() const
{
    return *m_storage;
}

/*!
* Get a constant reference to the kernel iterator.
*
* \result A constant reference to the kernel iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
const PiercedKernelIterator<id_t> & PiercedStorageIterator<value_t, id_t, value_no_cv_t>::getKernelIterator() const
{
    return static_cast<const PiercedKernelIterator<id_t> &>(*this);
}

/**
* Gets the id of the current element.
*
* \return The id of the current element or the fallback value if the iterator
* points to an invalid position.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
id_t PiercedStorageIterator<value_t, id_t, value_no_cv_t>::getId(const id_t &fallback) const noexcept
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
template<typename value_t, typename id_t, typename value_no_cv_t>
std::size_t PiercedStorageIterator<value_t, id_t, value_no_cv_t>::getRawIndex() const noexcept
{
    return m_pos;
}

/**
* Gets the values of the current element.
*
* \param k is the index of the requested field
* \return The values of the current element.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
__PSI_REFERENCE__ PiercedStorageIterator<value_t, id_t, value_no_cv_t>::getValue(std::size_t k) const
{
    if (m_pos >= m_kernel->m_end_pos) {
        throw std::out_of_range("Iterator points to an invalid position.");
    }

    return m_storage->rawAt(m_pos, k);
}

/**
* Pre-increment operator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageIterator<value_t, id_t, value_no_cv_t> & PiercedStorageIterator<value_t, id_t, value_no_cv_t>::operator++()
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
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageIterator<value_t, id_t, value_no_cv_t> PiercedStorageIterator<value_t, id_t, value_no_cv_t>::operator++(int)
{
    PiercedStorageIterator tmp(m_storage, m_pos);

    ++(*this);

    return tmp;
}

/**
* Deference operator.
*
* \result A reference to the element currently pointed to by the iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
__PSI_REFERENCE__ PiercedStorageIterator<value_t, id_t, value_no_cv_t>::operator*() const
{
    return m_storage->rawAt(m_pos, 0);
}

/**
* Deference operator.
*
* \result A reference to the element currently pointed to by the iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
__PSI_POINTER__ PiercedStorageIterator<value_t, id_t, value_no_cv_t>::operator->() const
{
    return m_storage->rawData(m_pos);
}

/**
* Converts the iterator to a const_iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
template<typename U, typename U_no_cv,
         typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type>
PiercedStorageIterator<value_t, id_t, value_no_cv_t>::operator PiercedStorageIterator<const U_no_cv, id_t>() const
{
    return PiercedStorageIterator<const U_no_cv, id_t>(m_storage, m_pos);
}

}

#endif
