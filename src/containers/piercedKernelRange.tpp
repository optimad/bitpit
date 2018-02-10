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

#ifndef __BITPIT_PIERCED_KERNEL_RANGE_TPP__
#define __BITPIT_PIERCED_KERNEL_RANGE_TPP__

namespace bitpit {

/*!
* Constructor.
*/
template<typename id_t>
PiercedKernelRange<id_t>::PiercedKernelRange()
    : m_cbegin(),
      m_cend()
{
}

/*!
* Constructor.
*
* \param kernel is the kernel that will be associated to the range
*/
template<typename id_t>
PiercedKernelRange<id_t>::PiercedKernelRange(const kernel_t *kernel)
    : m_cbegin(kernel->cbegin()),
      m_cend(kernel->cend())
{
}

/*!
* Constructor.
*
* \param kernel is the kernel that will be associated to the range
* \param first is the id of the first element in the range
* \param last is the id of the last element in the range
*/
template<typename id_t>
PiercedKernelRange<id_t>::PiercedKernelRange(const kernel_t *kernel, id_t first, id_t last)
    : m_cbegin(kernel->find(first)),
      m_cend(++(kernel->find(last)))
{
}

/*!
* Constructor.
*
* \param begin is the begin of the range
* \param end is the end of the range
*/
template<typename id_t>
PiercedKernelRange<id_t>::PiercedKernelRange(const const_iterator &begin, const const_iterator &end)
    : m_cbegin(begin),
      m_cend(end)
{
    if (&(begin.getKernel()) != &(end.getKernel())) {
        throw std::runtime_error("The two iterators belong to different kernels");
    }
}

/*!
* Exchanges the values of the current range and the range recevied as argument.
*
* \param other the range to exchange values with
*/
template<typename id_t>
void PiercedKernelRange<id_t>::swap(PiercedKernelRange &other) noexcept
{
    std::swap(m_cbegin, other.m_cbegin);
    std::swap(m_cend, other.m_cend);
}

/*!
* Returns a constant iterator pointing to the first element in the range.
*
* \result A constant iterator pointing to the first element in the range.
*/
template<typename id_t>
typename PiercedKernelRange<id_t>::const_iterator PiercedKernelRange<id_t>::begin() const noexcept
{
    return m_cbegin;
}

/*!
* Returns a constant iterator pointing to the past-the-end element in the
* range.
*
* \result A constant iterator pointing to the past-the-end element in the
* range.
*/
template<typename id_t>
typename PiercedKernelRange<id_t>::const_iterator PiercedKernelRange<id_t>::end() const noexcept
{
    return m_cend;
}

/*!
* Returns a constant iterator pointing to the first element in the range.
*
* \result A constant iterator pointing to the first element in the range.
*/
template<typename id_t>
typename PiercedKernelRange<id_t>::const_iterator PiercedKernelRange<id_t>::cbegin() const noexcept
{
    return m_cbegin;
}

/*!
* Returns a constant iterator pointing to the past-the-end element in the
* range.
*
* \result A constant iterator pointing to the past-the-end element in the
* range.
*/
template<typename id_t>
typename PiercedKernelRange<id_t>::const_iterator PiercedKernelRange<id_t>::cend() const noexcept
{
    return m_cend;
}

}

#endif
