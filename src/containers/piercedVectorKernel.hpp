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

#ifndef __BITPIT_PIERCED_VECTOR_KERNEL_HPP__
#define __BITPIT_PIERCED_VECTOR_KERNEL_HPP__

#include <cassert>

#include "piercedKernel.hpp"

namespace bitpit {

/**
* \ingroup containers
*
* \brief Base class for the pierced vector kernels.
*/
class BasePiercedVectorKernel {
};

/**
* \ingroup containers
*
* \brief Kernel of the pierced vector.
*/
template<typename id_t = long>
class PiercedVectorKernel : public BasePiercedVectorKernel,
                            public PiercedKernel<id_t> {

public:
    // Methods that extract information on the ids contained in the kernel
    std::size_t rawIndex(id_t id) const;
    bool exists(id_t id) const;

    // Methods for handing the synchronization
    template<typename data_t>
    void registerStorage(PiercedStorage<data_t, id_t> *storage, PiercedSyncMaster::SyncMode syncMode);
    template<typename data_t>
    void unregisterStorage(PiercedStorage<data_t, id_t> *storage);
    template<typename data_t>
    bool isStorageRegistered(const PiercedStorage<data_t, id_t> *storage) const;
    template<typename data_t>
    PiercedSyncMaster::SyncMode getStorageSyncMode(const PiercedStorage<data_t, id_t> *storage) const;

    // Dump and restore
    using PiercedKernel<id_t>::dump;
    using PiercedKernel<id_t>::restore;

protected:
    // Constructors
    using PiercedKernel<id_t>::PiercedKernel;

    // Methods that modify the elements stored in the kernel
    using PiercedKernel<id_t>::fillHead;
    using PiercedKernel<id_t>::fillTail;
    using PiercedKernel<id_t>::fillAfter;
    using PiercedKernel<id_t>::fillBefore;
    using PiercedKernel<id_t>::fillAppend;
    using PiercedKernel<id_t>::fillHole;
    using PiercedKernel<id_t>::fillInsert;

    using PiercedKernel<id_t>::moveAfter;
    using PiercedKernel<id_t>::moveBefore;

    using PiercedKernel<id_t>::swap;

    using PiercedKernel<id_t>::erase;
    using PiercedKernel<id_t>::popBack;

    // Methods that extract information about the elements of the kernel
    using PiercedKernel<id_t>::back;
    using PiercedKernel<id_t>::front;

};

}

// Include the implementation
#include "piercedVectorKernel.tpp"

#endif
