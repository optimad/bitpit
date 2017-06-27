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

#ifndef __BITPIT_PIERCED_STORAGE_RANGE_TPP__
#define __BITPIT_PIERCED_STORAGE_RANGE_TPP__

namespace bitpit {

/*!
* Constructor.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange()
    : PiercedKernelRange<id_t>(), m_storage(nullptr)
{
}

/*!
* Constructor.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange(storage_t *storage)
    : PiercedKernelRange<id_t>(&(storage->getKernel())), m_storage(storage)
{
}

/*!
* Constructor.
*
* \param first is the id of the first element in the range
* \param last is the id of the last element in the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange(storage_t *storage, id_t first, id_t last)
    : PiercedKernelRange<id_t>(&(storage->getKernel()), first, last), m_storage(storage)
{
}

/*!
* Constructor.
*
* \param begin is the begin of the range
* \param end is the end of the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange(const iterator &begin, const iterator &end)
    : PiercedKernelRange<id_t>(begin.getKernelIterator(), end.getKernelIterator()), m_storage(&(begin.getStorage()))
{
    if (&(begin.getStorage()) != &(end.getStorage())) {
        throw std::runtime_error("The two iterators belong to different storages");
    }
}

/*!
* Exchanges the values of the current iterator and
* the iterator recevied as argument.
*
* \param other the iterator to exchange values with
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
void PiercedStorageRange<value_t, id_t, value_no_cv_t>::swap(PiercedStorageRange &other) noexcept
{
    PiercedKernelRange<id_t>::swap(other);

    std::swap(m_storage, other.m_storage);
}

/*!
* Get a constant reference of the storage associated with the range.
*
* \result A constant reference of the storage associated with the range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::storage_type & PiercedStorageRange<value_t, id_t, value_no_cv_t>::getStorage() const
{
    return *m_storage;
}

/*!
* Get a constant reference to the kernel range.
*
* \result A constant reference to the kernel range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
const PiercedKernelRange<id_t> & PiercedStorageRange<value_t, id_t, value_no_cv_t>::getKernelRange() const
{
    return static_cast<const PiercedKernelRange<id_t> &>(*this);
}

/*!
* Returns an iterator pointing to the first element in the range.
*
* \result An iterator pointing to the first element in the range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
template<typename U, typename U_no_cv,
         typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::begin() noexcept
{
    return m_storage->getIteratorFromRawIndex(PiercedKernelRange<id_t>::m_begin_pos);
}

/*!
* Returns a constant iterator pointing to the past-the-end element in the
* range.
*
* \result A constant iterator pointing to the past-the-end element in the
* range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
template<typename U, typename U_no_cv,
         typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::end() noexcept
{
    return m_storage->getIteratorFromRawIndex(PiercedKernelRange<id_t>::m_end_pos);
}

/*!
* Returns a constant iterator pointing to the first element in the range.
*
* \result A constant iterator pointing to the first element in the range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::const_iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::begin() const noexcept
{
    return m_storage->getConstIteratorFromRawIndex(PiercedKernelRange<id_t>::m_begin_pos);
}

/*!
* Returns a constant iterator pointing to the past-the-end element in the
* range.
*
* \result A constant iterator pointing to the past-the-end element in the
* range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::const_iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::end() const noexcept
{
    return m_storage->getConstIteratorFromRawIndex(PiercedKernelRange<id_t>::m_end_pos);
}

/*!
* Returns a constant iterator pointing to the first element in the range.
*
* \result A constant iterator pointing to the first element in the range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::const_iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::cbegin() const noexcept
{
    return m_storage->getConstIteratorFromRawIndex(PiercedKernelRange<id_t>::m_begin_pos);
}

/*!
* Returns a constant iterator pointing to the past-the-end element in the
* range.
*
* \result A constant iterator pointing to the past-the-end element in the
* range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::const_iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::cend() const noexcept
{
    return m_storage->getConstIteratorFromRawIndex(PiercedKernelRange<id_t>::m_end_pos);
}

}

#endif
