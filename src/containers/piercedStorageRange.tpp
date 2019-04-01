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
    : PiercedKernelRange<id_t>(),
      m_begin(), m_end()
{
}

/*!
* Constructor.
*
* \param storage is the storage that will be associated to the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange(storage_t *storage)
{
    initialize(storage);
}

/*!
* Constructor.
*
* \param storage is the storage that will be associated to the range
* \param first is the id of the first element in the range
* \param last is the id of the last element in the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange(storage_t *storage, id_t first, id_t last)
{
    initialize(storage, first, last);
}

/*!
* Constructor.
*
* \param begin is the begin of the range
* \param end is the end of the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedStorageRange<value_t, id_t, value_no_cv_t>::PiercedStorageRange(const iterator &begin, const iterator &end)
{
    initialize(begin, end);
}

/*!
* Initialize the range.
*
* \param kernel is the kernel that will be associated to the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
void PiercedStorageRange<value_t, id_t, value_no_cv_t>::initialize(const storage_t *storage)
{
    PiercedKernelRange<id_t>::initialize(&(storage->getKernel()));

    m_begin = storage->begin();
    m_end   = storage->end();
}

/*!
* Initialize the range.
*
* \param kernel is the kernel that will be associated to the range
* \param first is the id of the first element in the range
* \param last is the id of the last element in the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
void PiercedStorageRange<value_t, id_t, value_no_cv_t>::initialize(const storage_t *storage, id_t first, id_t last)
{
    PiercedKernelRange<id_t>::initialize(&(storage->getKernel()), first, last);

    m_begin = storage->find(first);
    m_end   = ++(storage->find(last));
}

/*!
* Initialize the range.
*
* \param begin is the begin of the range
* \param end is the end of the range
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
void PiercedStorageRange<value_t, id_t, value_no_cv_t>::initialize(const iterator &begin, const iterator &end)
{
    if (&(begin.getStorage()) != &(end.getStorage())) {
        throw std::runtime_error("The two iterators belong to different storages");
    }

    PiercedKernelRange<id_t>::initialize(begin.getKernelIterator(), end.getKernelIterator());

    m_begin = begin;
    m_end   = end;
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

    std::swap(m_begin, other.m_begin);
    std::swap(m_end, other.m_end);
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
    return m_begin;
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
    return m_end;
}

/*!
* Returns a constant iterator pointing to the first element in the range.
*
* \result A constant iterator pointing to the first element in the range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::const_iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::begin() const noexcept
{
    return m_begin;
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
    return m_end;
}

/*!
* Returns a constant iterator pointing to the first element in the range.
*
* \result A constant iterator pointing to the first element in the range.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedStorageRange<value_t, id_t, value_no_cv_t>::const_iterator PiercedStorageRange<value_t, id_t, value_no_cv_t>::cbegin() const noexcept
{
    return m_begin;
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
    return m_end;
}

}

#endif
