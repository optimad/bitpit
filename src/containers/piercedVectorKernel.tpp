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

#ifndef __BITPIT_PIERCED_VECTOR_KERNEL_TPP__
#define __BITPIT_PIERCED_VECTOR_KERNEL_TPP__

namespace bitpit {

/*!
* Gets the raw index associated to the element with the specified id.
*
* \return The raw index associated to the element with the specified id.
*/
template<typename id_t>
std::size_t PiercedVectorKernel<id_t>::rawIndex(id_t id) const
{
    return PiercedKernel<id_t>::getRawIndex(id);
}


/**
* Checks if a given id exists in the container.
*
* \param id is the id to look for
* \result Returns true is the given id exists in the container, otherwise it
* returns false.
*/
template<typename id_t>
bool PiercedVectorKernel<id_t>::exists(id_t id) const
{
    return PiercedKernel<id_t>::contains(id);
}

}

#endif
