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

#ifndef __BITPIT_PIERCED_VECTOR_HPP__
#define __BITPIT_PIERCED_VECTOR_HPP__

#include <cassert>
#include <limits>
#include <unordered_map>
#include <type_traits>
#include <vector>

#include "piercedKernel.hpp"
#include "piercedIterator.hpp"
#include "piercedStorage.hpp"

#define __PV_REFERENCE__       typename PiercedVector<value_t, id_t>::reference
#define __PV_CONST_REFERENCE__ typename PiercedVector<value_t, id_t>::const_reference
#define __PV_POINTER__         typename PiercedVector<value_t, id_t>::pointer
#define __PV_CONST_POINTER__   typename PiercedVector<value_t, id_t>::const_pointer

namespace bitpit {

/**
* \ingroup containers
*
* \brief Base class for the pierced vectors.
*/
class BasePiercedVector {
};

/**
* \ingroup containers
*
* \brief Metafunction for generating a pierced vector.
*
* \details
* Usage: use <tt>PiercedVector<value_t, id_t></tt> to declare a pierced
* vector.
*
* \tparam value_t is the type of the elements stored in the vector
* \tparam id_t is the type of the ids associated to the elements
*/
template<typename value_t, typename id_t = long>
class PiercedVector : public BasePiercedVector,
                      private PiercedKernel<id_t>,
                      private PiercedStorage<value_t, id_t> {

public:
    // Typedefs

    /**
    * Kernel template
    */
    template<typename PK_id_t>
    using Kernel = PiercedKernel<PK_id_t>;

    /**
    * Kernel type
    */
    typedef Kernel<id_t> kernel_type;

    /*!
    * Type of data stored in the container
    */
    typedef typename PiercedStorage<value_t, id_t>::value_type value_type;

    /*!
    * Type of ids stored in the container
    */
    typedef typename PiercedKernel<id_t>::id_type id_type;

    /**
    * Reference
    */
    typedef typename PiercedStorage<value_t, id_t>::reference reference;

    /**
    * Constant reference
    */
    typedef typename PiercedStorage<value_t, id_t>::const_reference const_reference;

    /**
    * Pointer
    */
    typedef typename PiercedStorage<value_t, id_t>::pointer pointer;

    /**
    * Constant pointer
    */
    typedef typename PiercedStorage<value_t, id_t>::const_pointer const_pointer;

    /**
    * Iterator
    */
    typedef typename PiercedStorage<value_t, id_t>::iterator iterator;

    /**
    * Constant iterator
    */
    typedef typename PiercedStorage<value_t, id_t>::const_iterator const_iterator;

    /**
    * Raw iterator
    */
    typedef typename PiercedStorage<value_t, id_t>::raw_iterator raw_iterator;

    /**
    * Raw constant iterator
    */
    typedef typename PiercedStorage<value_t, id_t>::raw_const_iterator raw_const_iterator;

    /**
    * Range
    */
    typedef typename PiercedStorage<value_t, id_t>::range range;

    /**
    * Constant range
    */
    typedef typename PiercedStorage<value_t, id_t>::const_range const_range;

    // Contructors
    PiercedVector();
    PiercedVector(std::size_t n);

    // Methods that modify the contents of the container
    using PiercedKernel<id_t>::updateId;

    iterator reclaim(id_t id);
    iterator reclaimAfter(const id_t &referenceId, id_t id);
    iterator reclaimBack(id_t id);
    iterator reclaimBefore(const id_t &referenceId, id_t id);

    iterator moveAfter(const id_t &referenceId, id_t id, bool delayed = false);
    iterator moveBefore(const id_t &referenceId, id_t id, bool delayed = false);

    iterator insert(id_t id, const value_t &value);
    iterator insertAfter(const id_t &referenceId, id_t id, const value_t &value);
    iterator insertBefore(const id_t &referenceId, id_t id, const value_t &value);

    iterator replace(id_t id, value_t &&value);

    iterator pushBack(id_t id, const value_t &value);

    template<typename... Args>
    iterator emplace(id_t id, Args&&... args);
    template<typename... Args>
    iterator emplaceAfter(const id_t &referenceId, id_t id, Args&&... args);
    template<typename... Args>
    void emplaceBack(id_t id, Args&&... args);
    template<typename... Args>
    iterator emplaceBefore(const id_t &referenceId, id_t id, Args&&... args);

    template<typename... Args>
    iterator emreplace(id_t id, Args&&... args);

    iterator erase(id_t id, bool delayed = false);
    void popBack();

    void swap(id_t id_first, id_t id_second);

    // Methods that modify the container as a whole
    using PiercedKernel<id_t>::flush;

    void clear(bool release = true);
    void reserve(std::size_t n);
    void resize(std::size_t n);
    void sort();
    void squeeze();
    void shrinkToFit();
    void swap(PiercedVector &x) noexcept;

    // Methods that extract information about the container
    using PiercedKernel<id_t>::capacity;
    using PiercedKernel<id_t>::contiguous;
    using PiercedKernel<id_t>::empty;
    using PiercedKernel<id_t>::isIteratorSlow;
    using PiercedKernel<id_t>::maxSize;
    using PiercedKernel<id_t>::size;

    const PiercedKernel<id_t> & getKernel() const;
    const PiercedStorage<value_t, id_t> & getStorage() const;

    void dump();

    // Methods that extract information on the contents of the container
    using PiercedKernel<id_t>::contains;
    using PiercedKernel<id_t>::getRawIndex;
    using PiercedKernel<id_t>::evalFlatIndex;

    std::size_t rawIndex(id_t id) const;
    bool exists(id_t id) const;

    using PiercedKernel<id_t>::getIds;
    using PiercedKernel<id_t>::getSizeMarker;

    // Methods that extract the contents of the container
    using PiercedStorage<value_t, id_t>::data;

    __PV_REFERENCE__ back();
    __PV_CONST_REFERENCE__ back() const;

    __PV_REFERENCE__ front();
    __PV_CONST_REFERENCE__ front() const;

    __PV_REFERENCE__ at(id_t id);
    __PV_CONST_REFERENCE__ at(id_t id) const;

    __PV_REFERENCE__ rawAt(std::size_t pos);
    __PV_CONST_REFERENCE__ rawAt(std::size_t pos) const;

    __PV_CONST_REFERENCE__ operator[](id_t id) const;
    __PV_REFERENCE__ operator[](id_t id);

    const_iterator find(id_t id) const;
    iterator find(id_t id);

    // Iterators
    using PiercedStorage<value_t, id_t>::getIterator;
    using PiercedStorage<value_t, id_t>::getConstIterator;

    using PiercedStorage<value_t, id_t>::getIteratorFromRawIndex;
    using PiercedStorage<value_t, id_t>::getConstIteratorFromRawIndex;

    using PiercedStorage<value_t, id_t>::begin;
    using PiercedStorage<value_t, id_t>::end;
    using PiercedStorage<value_t, id_t>::cbegin;
    using PiercedStorage<value_t, id_t>::cend;

    using PiercedStorage<value_t, id_t>::rawBegin;
    using PiercedStorage<value_t, id_t>::rawEnd;
    using PiercedStorage<value_t, id_t>::rawCbegin;
    using PiercedStorage<value_t, id_t>::rawCend;

    // External storage
    using PiercedStorage<value_t, id_t>::getSyncMaster;

    template<typename data_t>
    void registerStorage(PiercedStorage<data_t, id_t> *storage, PiercedSyncMaster::SyncMode syncMode);
    template<typename data_t>
    void unregisterStorage(PiercedStorage<data_t, id_t> *storage);
    template<typename data_t>
    bool isStorageRegistered(const PiercedStorage<data_t, id_t> *storage) const;
    template<typename data_t>
    PiercedSyncMaster::SyncMode getStorageSyncMode(const PiercedStorage<data_t, id_t> *storage) const;

    using PiercedKernel<id_t>::sync;

private:
    typedef typename PiercedKernel<id_t>::FillAction FillAction;
    typedef typename PiercedKernel<id_t>::MoveAction MoveAction;
    typedef typename PiercedKernel<id_t>::SwapAction SwapAction;
    typedef typename PiercedKernel<id_t>::EraseAction EraseAction;

    iterator reclaimValue(const FillAction &action);
    iterator insertValue(const FillAction &action, const value_t &value);
    template<typename... Args>
    iterator emplaceValue(const FillAction &action, Args&&... args);

    iterator moveValue(const MoveAction &action);

    void swapValues(const SwapAction &action);

    iterator eraseValue(const EraseAction &action);

};

}

// Include the implementation
#include "piercedVector.tpp"

#endif
