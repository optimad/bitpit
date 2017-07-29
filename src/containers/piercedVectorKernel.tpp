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

/**
* Register the specified storage
*
* \param storage is the storage that will be registered
* \param syncMode is the synchronization mode that will be used for the storage
*/
template<typename id_t>
template<typename data_t>
void PiercedVectorKernel<id_t>::registerStorage(PiercedStorage<data_t, id_t> *storage, PiercedSyncMaster::SyncMode syncMode)
{
    storage->setKernel(this, syncMode);
}

/**
* Unregister the specified storage->
*
* \param storage is the storage that will be unregistered
*/
template<typename id_t>
template<typename data_t>
void PiercedVectorKernel<id_t>::unregisterStorage(PiercedStorage<data_t, id_t> *storage)
{
    if (dynamic_cast<const BasePiercedStorage *>(this) == storage) {
        throw std::out_of_range("Internal storage cannot be unregistered.");
    }

    storage->unsetKernel();
}

/**
* Check if te specified storage is registered.
*
* \param storage is the storage to check
*/
template<typename id_t>
template<typename data_t>
bool PiercedVectorKernel<id_t>::isStorageRegistered(const PiercedStorage<data_t, id_t> *storage) const
{
    return PiercedKernel<id_t>::isSlaveRegistered(storage);
}

/**
* Get the synchronization mode for the specified storage
*
* \param storage is the storage for which the synchronization mode is requested
* \result The synchronization mode of the storage->
*/
template<typename id_t>
template<typename data_t>
typename PiercedSyncMaster::SyncMode PiercedVectorKernel<id_t>::getStorageSyncMode(const PiercedStorage<data_t, id_t> *storage) const
{
    return PiercedKernel<id_t>::getSlaveSyncMode(storage);
}

}

#endif
