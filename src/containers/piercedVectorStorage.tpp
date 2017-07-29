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

#ifndef __BITPIT_PIERCED_VECTOR_STORAGE_TPP__
#define __BITPIT_PIERCED_VECTOR_STORAGE_TPP__

namespace bitpit {

/**
* Returns a reference to the first element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A reference to the first element of the container.
*/
template<typename value_t, typename id_t>
__PVS_REFERENCE__ PiercedVectorStorage<value_t, id_t>::front()
{
    return front(0);
}

/**
* Returns a constant reference to the first element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A constant reference to the first element of the container.
*/
template<typename value_t, typename id_t>
__PVS_CONST_REFERENCE__ PiercedVectorStorage<value_t, id_t>::front() const
{
    return front(0);
}

/**
* Returns a reference to the last element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A reference to the last element of the container.
*/
template<typename value_t, typename id_t>
__PVS_REFERENCE__ PiercedVectorStorage<value_t, id_t>::back()
{
    return back(0);
}

/**
* Returns a constant reference to the last element of the container.
*
* If the container is empty, an exception is thrown.
*
* \result A constant reference to the last element of the container.
*/
template<typename value_t, typename id_t>
__PVS_CONST_REFERENCE__ PiercedVectorStorage<value_t, id_t>::back() const
{
    return back(0);
}

/**
* Returns a reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PVS_REFERENCE__ PiercedVectorStorage<value_t, id_t>::at(id_t id)
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Returns a constant reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A constant reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PVS_CONST_REFERENCE__ PiercedVectorStorage<value_t, id_t>::at(id_t id) const
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Returns a reference to the element at the specified position.
*
* \param pos the position of the element
* \result A reference to the element in the specified position.
*/
template<typename value_t, typename id_t>
__PVS_REFERENCE__ PiercedVectorStorage<value_t, id_t>::rawAt(std::size_t pos)
{
    return PiercedStorage<value_t, id_t>::rawAt(pos, 0);
}

/**
* Returns a constant reference to the element at the specified position.
*
* \param pos the position of the element
* \result A constant reference to the element in the specified position.
*/
template<typename value_t, typename id_t>
__PVS_CONST_REFERENCE__ PiercedVectorStorage<value_t, id_t>::rawAt(std::size_t pos) const
{
    return PiercedStorage<value_t, id_t>::rawAt(pos, 0);
}

/**
* Returns a constant reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A constant reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PVS_CONST_REFERENCE__ PiercedVectorStorage<value_t, id_t>::operator[](id_t id) const
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

/**
* Returns a reference to the element with the specified id.
*
* If there is no element with the specified id, an exception is thrown.
*
* \param id is the id of the element
* \result A reference to the element with the specified id.
*/
template<typename value_t, typename id_t>
__PVS_REFERENCE__ PiercedVectorStorage<value_t, id_t>::operator[](id_t id)
{
    return PiercedStorage<value_t, id_t>::at(id, 0);
}

}

#endif
