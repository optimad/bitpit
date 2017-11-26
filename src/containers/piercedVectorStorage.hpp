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

#ifndef __BITPIT_PIERCED_VECTOR_STORAGE_HPP__
#define __BITPIT_PIERCED_VECTOR_STORAGE_HPP__

#include <cassert>

#include "piercedStorage.hpp"
#include "piercedVectorKernel.hpp"

#define __PVS_REFERENCE__       typename PiercedVectorStorage<value_t, id_t>::reference
#define __PVS_CONST_REFERENCE__ typename PiercedVectorStorage<value_t, id_t>::const_reference
#define __PVS_POINTER__         typename PiercedVectorStorage<value_t, id_t>::pointer
#define __PVS_CONST_POINTER__   typename PiercedVectorStorage<value_t, id_t>::const_pointer

namespace bitpit {

class BasePiercedVectorStorage {

protected:
    BasePiercedVectorStorage();

};

/**
* \ingroup containers
*
* \brief Kernel of the pierced vector.
*/
template<typename value_t, typename id_t = long>
class PiercedVectorStorage : public BasePiercedVectorStorage,
                             public PiercedStorage<value_t, id_t> {

public:
    // Contructors
    using PiercedStorage<value_t, id_t>::PiercedStorage;

    // Methods that extract the contents of the container
    using PiercedStorage<value_t, id_t>::data;

    __PVS_REFERENCE__ back();
    __PVS_CONST_REFERENCE__ back() const;

    __PVS_REFERENCE__ front();
    __PVS_CONST_REFERENCE__ front() const;

    __PVS_REFERENCE__ at(id_t id);
    __PVS_CONST_REFERENCE__ at(id_t id) const;

    __PVS_REFERENCE__ rawAt(std::size_t pos);
    __PVS_CONST_REFERENCE__ rawAt(std::size_t pos) const;

    __PVS_CONST_REFERENCE__ operator[](id_t id) const;
    __PVS_REFERENCE__ operator[](id_t id);

    using PiercedStorage<value_t, id_t>::find;

    using PiercedStorage<value_t, id_t>::rawFind;

    // Iterators
    using PiercedStorage<value_t, id_t>::begin;
    using PiercedStorage<value_t, id_t>::end;
    using PiercedStorage<value_t, id_t>::cbegin;
    using PiercedStorage<value_t, id_t>::cend;

    using PiercedStorage<value_t, id_t>::rawBegin;
    using PiercedStorage<value_t, id_t>::rawEnd;
    using PiercedStorage<value_t, id_t>::rawCbegin;
    using PiercedStorage<value_t, id_t>::rawCend;

    // Methods for handing the synchronization
    using PiercedStorage<value_t, id_t>::unsetKernel;

    void setStaticKernel(const PiercedVectorKernel<id_t> *kernel);
    void setDynamicKernel(PiercedVectorKernel<id_t> *kernel, PiercedSyncMaster::SyncMode syncMode);
    const PiercedVectorKernel<id_t> * getKernel() const;

    // Dump and restore
    using PiercedStorage<value_t, id_t>::dump;
    using PiercedStorage<value_t, id_t>::restore;

};

}

// Include the implementation
#include "piercedVectorStorage.tpp"

#endif
