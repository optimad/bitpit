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

#ifndef __BITPIT_PIERCED_ITERATOR_TPP__
#define __BITPIT_PIERCED_ITERATOR_TPP__

namespace bitpit {

/*!
	Creates a new uninitialized iterator
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedIterator<value_t, id_t, value_no_cv_t>::PiercedIterator()
	: m_container(nullptr), m_pos(0)
{
}

/*!
	Exchanges the values of the current iterator and
	the iterator recevied as argument.

	\param other the iterator to exchange values with
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
void PiercedIterator<value_t, id_t, value_no_cv_t>::swap(PiercedIterator& other) noexcept
{
	std::swap(m_container, other.m_container);
	std::swap(m_pos, other.m_pos);
}

/*!
	Get a constant reference of the container associated with the iterator.

	\result A constant reference of the container associated with the iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
typename PiercedIterator<value_t, id_t, value_no_cv_t>::container_type & PiercedIterator<value_t, id_t, value_no_cv_t>::getContainer() const
{
	return *m_container;
}

/*!
	Gets the id of the current element.

	\param fallback is the fallback value to be returned if the iterator
	points to an invalid position
	\return The id of the current element or the fallback value if the
	the iterator points to an invalid position.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
id_t PiercedIterator<value_t, id_t, value_no_cv_t>::getId(const id_t &fallback) const noexcept
{
	id_t id;
	if (m_pos >= m_container->m_end_pos) {
		id = fallback;
		return id;
	}

	id = m_container->m_ids[m_pos];
	if (id >= 0) {
		return id;
	} else {
		id = fallback;
		return id;
	}
}

/*!
	Gets the position of the current element.

	\return The position of the current element.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
std::size_t PiercedIterator<value_t, id_t, value_no_cv_t>::getRawIndex() const noexcept
{
	return m_pos;
}

/*!
	Pre-increment operator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedIterator<value_t, id_t, value_no_cv_t> & PiercedIterator<value_t, id_t, value_no_cv_t>::operator++()
{
	size_t delta = 1;
	while (true) {
		m_pos += delta;
		if (m_pos >= m_container->m_end_pos) {
			m_pos = m_container->m_end_pos;
			return *this;
		}

		id_t id = m_container->m_ids[m_pos];
		if (id >= 0) {
			return *this;
		} else {
			delta = - id;
		}
	}

	assert(false);
}

/*!
	Post-increment operator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedIterator<value_t, id_t, value_no_cv_t> PiercedIterator<value_t, id_t, value_no_cv_t>::operator++(int)
{
	PiercedIterator tmp(m_container, m_pos);

	++(*this);

	return tmp;
}

/*!
	Deference operator.

	\result A reference to the element currently pointed to by the
			iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
__PI_REFERENCE__ PiercedIterator<value_t, id_t, value_no_cv_t>::operator*() const
{
	return m_container->m_v[m_pos];
}

/*!
	Deference operator.

	\result A reference to the element currently pointed to by the
			iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
__PI_POINTER__ PiercedIterator<value_t, id_t, value_no_cv_t>::operator->() const
{
	return &(m_container->m_v[m_pos]);
}

/*!
	Converts the iterator to a const_iterator.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
template<typename U, typename U_no_cv,
         typename std::enable_if<std::is_same<U, U_no_cv>::value, int>::type>
PiercedIterator<value_t, id_t, value_no_cv_t>::operator PiercedIterator<const U_no_cv, id_t>() const
{
	return PiercedIterator<const U_no_cv, id_t>(m_container, m_pos);
}

/*!
	Creates a new iterator and initializes it with the position of
	the const base iterator recevied in input.
*/
template<typename value_t, typename id_t, typename value_no_cv_t>
PiercedIterator<value_t, id_t, value_no_cv_t>::PiercedIterator(container_t *container, const size_t &pos)
	: m_container(container), m_pos(pos)
{
}

}

#endif
