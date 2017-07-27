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

#include "piercedVectorKernel.hpp"
#include "piercedVectorStorage.hpp"

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
                      public PiercedVectorKernel<id_t>,
                      public PiercedVectorStorage<value_t, id_t> {

public:
    // Typedefs

    /*!
    * Type of ids stored in the container
    */
    typedef PiercedVectorKernel<id_t> kernel_type;

    /*!
    * Type of ids stored in the container
    */
    typedef typename PiercedVectorKernel<id_t>::id_type id_type;

    /*!
    * Type of data stored in the container
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::value_type value_type;

    /**
    * Iterator
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::iterator iterator;

    /**
    * Constant iterator
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::const_iterator const_iterator;

    /**
    * Raw iterator
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::raw_iterator raw_iterator;

    /**
    * Raw constant iterator
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::raw_const_iterator raw_const_iterator;

    /**
    * Range
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::range range;

    /**
    * Constant range
    */
    typedef typename PiercedVectorStorage<value_t, id_t>::const_range const_range;

    // Contructors
    PiercedVector();
    PiercedVector(std::size_t n);
    PiercedVector(const PiercedVector<value_t, id_t> &x);
    PiercedVector(PiercedVector<value_t, id_t> &&x) = default;

    // Methods that modify the contents of the container
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

    template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type * = nullptr>
    iterator emreclaim(id_t id, Args&&... args);
    template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type * = nullptr>
    iterator emreclaimAfter(const id_t &referenceId, id_t id, Args&&... args);
    template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type * = nullptr>
    void emreclaimBack(id_t id, Args&&... args);
    template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type * = nullptr>
    iterator emreclaimBefore(const id_t &referenceId, id_t id, Args&&... args);

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
    void clear(bool release = true);
    void reserve(std::size_t n);
    void resize(std::size_t n);
    void sort();
    void squeeze();
    void shrinkToFit();
    void swap(PiercedVector &x) noexcept;

    // Methods that extract information about the container
    const PiercedVectorKernel<id_t> & getKernel() const;
    const PiercedVectorStorage<value_t, id_t> & getStorage() const;

    void dump();

    // Methods that extract the contents of the container
    using PiercedVectorStorage<value_t, id_t>::back;
    using PiercedVectorStorage<value_t, id_t>::front;

    using PiercedVectorStorage<value_t, id_t>::at;
    using PiercedVectorStorage<value_t, id_t>::rawAt;

    using PiercedVectorStorage<value_t, id_t>::operator[];

    // Iterators
    using PiercedVectorStorage<value_t, id_t>::begin;
    using PiercedVectorStorage<value_t, id_t>::end;
    using PiercedVectorStorage<value_t, id_t>::cbegin;
    using PiercedVectorStorage<value_t, id_t>::cend;

    using PiercedVectorStorage<value_t, id_t>::rawBegin;
    using PiercedVectorStorage<value_t, id_t>::rawEnd;
    using PiercedVectorStorage<value_t, id_t>::rawCbegin;
    using PiercedVectorStorage<value_t, id_t>::rawCend;

    using PiercedVectorStorage<value_t, id_t>::find;
    using PiercedVectorStorage<value_t, id_t>::rawFind;

    // Dump and restore
    template<typename T = value_t, typename std::enable_if<PiercedVectorStorage<T, id_t>::has_dump_restore>::type * = nullptr>
    void restore(std::istream &stream);

    template<typename T = value_t, typename std::enable_if<PiercedVectorStorage<T, id_t>::has_dump_restore>::type * = nullptr>
    void dump(std::ostream &stream) const;

protected:
    using PiercedVectorStorage<value_t, id_t>::setStaticKernel;
    using PiercedVectorStorage<value_t, id_t>::setDynamicKernel;
    using PiercedVectorStorage<value_t, id_t>::unsetKernel;
    using PiercedVectorStorage<value_t, id_t>::getKernel;
    using PiercedVectorStorage<value_t, id_t>::getKernelType;
    using PiercedVectorStorage<value_t, id_t>::getSyncMode;

private:
    typedef typename PiercedVectorKernel<id_t>::FillAction FillAction;
    typedef typename PiercedVectorKernel<id_t>::MoveAction MoveAction;
    typedef typename PiercedVectorKernel<id_t>::SwapAction SwapAction;
    typedef typename PiercedVectorKernel<id_t>::EraseAction EraseAction;

    iterator reclaimValue(const FillAction &action);
    iterator insertValue(const FillAction &action, const value_t &value);
    template<typename... Args, typename std::enable_if<PiercedVectorStorage<value_t, id_t>::template has_initialize<Args...>()>::type * = nullptr>
    iterator emreclaimValue(const FillAction &action, Args&&... args);
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
