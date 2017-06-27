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

#ifndef __BITPIT_PIERCED_STORAGE_HPP__
#define __BITPIT_PIERCED_STORAGE_HPP__

#include <vector>

#include "piercedIterator.hpp"
#include "piercedRange.hpp"
#include "piercedKernel.hpp"
#include "piercedSync.hpp"

#define __PS_REFERENCE__       typename PiercedStorage<value_t, id_t>::reference
#define __PS_CONST_REFERENCE__ typename PiercedStorage<value_t, id_t>::const_reference
#define __PS_POINTER__         typename PiercedStorage<value_t, id_t>::pointer
#define __PS_CONST_POINTER__   typename PiercedStorage<value_t, id_t>::const_pointer

namespace bitpit {

class BasePiercedStorage : protected virtual PiercedSyncSlave {

protected:
    BasePiercedStorage();

};

/**
* \ingroup containers
*
* \brief Metafunction for generating a pierced storage.
*
* \details
* Usage: use <tt>PiercedStorage<value_t, id_t></tt> to declare a pierced
* storage.
*
* \tparam value_t is the type of the elements stored in the storage
* \tparam id_t is the type of the ids associated to the elements
*/
template<typename value_t, typename id_t = long>
class PiercedStorage : public BasePiercedStorage {

// Friendships
template<typename PI_value_t, typename PI_id_t, typename PI_value_no_cv_t>
friend class PiercedIterator;

private:
    /**
    * Checks if the elements stored in the storage have restore capability
    */
    template <typename T>
    class check_restore
    {

    private:
        template<typename class_t>
        static std::true_type test_restore(decltype(std::declval<class_t>().restore(std::declval<std::ostream &>())) *);

        template<typename class_t>
        static std::false_type test_restore(...);

        template<typename class_t>
        static constexpr bool has_restore()
        {
            return std::is_same<decltype(test_restore<class_t>(nullptr)), std::true_type>();
        };

    public:
        static const bool value = has_restore<T>();

    };

    /**
    * Checks if the elements stored in the storage have dump capability
    */
    template <typename T>
    class check_dump
    {

    private:
        template<typename class_t>
        static std::true_type test_dump(decltype(std::declval<class_t>().dump(std::declval<std::istream &>())) *);

        template<typename class_t>
        static std::false_type test_dump(...);

        template<typename class_t>
        static constexpr bool has_dump()
        {
            return std::is_same<decltype(test_dump<class_t>(nullptr)), std::true_type>();
        };

    public:
        static const bool value = has_dump<T>();

    };

    /**
    * Checks if the elements stored in the storage have initialize capability
    */
    template <typename T, typename Ret, typename... Args>
    class check_initialize
    {

    private:
        template<typename C>
        static constexpr auto test_initialize(C *)
            -> typename std::is_same<decltype(std::declval<C>().initialize(std::declval<Args>()...)), Ret>::type;

        template<typename class_t>
        static constexpr auto test_initialize(...)
            -> std::false_type;

    public:
        static const bool value = std::is_same<decltype(test_initialize<T>(nullptr)), std::true_type>();

    };

public:
    // Template typedef

    /**
    * Kernel template
    */
    template<typename PK_id_t>
    using Kernel = PiercedKernel<PK_id_t>;

    // Typedefs

    /*!
    * Type of data stored in the container
    */
    typedef value_t value_type;

    /*!
    * Type of ids stored in the container
    */
    typedef typename Kernel<id_t>::id_type id_type;

    /**
    * Kernel type
    */
    typedef Kernel<id_t> kernel_t;

    /**
    * Container
    */
    typedef std::vector<value_t> container_t;

    /**
    * Reference
    */
    typedef typename container_t::reference reference;

    /**
    * Constant reference
    */
    typedef typename container_t::const_reference const_reference;

    /**
    * Pointer
    */
    typedef typename container_t::pointer pointer;

    /**
    * Constant pointer
    */
    typedef typename container_t::const_pointer const_pointer;

    /**
    * Iterator
    */
    typedef PiercedIterator<value_t, id_t> iterator;

    /**
    * Constant iterator
    */
    typedef PiercedIterator<const value_t, id_t> const_iterator;

    /**
    * Raw iterator
    */
    typedef typename std::vector<value_t>::iterator raw_iterator;

    /**
    * Raw constant iterator
    */
    typedef typename std::vector<value_t>::const_iterator raw_const_iterator;

    /**
    * Range
    */
    typedef PiercedRange<value_t, id_t> range;

    /**
    * Constant range
    */
    typedef PiercedRange<const value_t, id_t> const_range;

    /**
    * Checks if the storage has the 'restore' capability
    */
    static constexpr bool has_restore()
    {
        return check_restore<value_t>::value;
    };

    /**
    * Checks if the storage has the 'dump' capability
    */
    static constexpr bool has_dump()
    {
        return check_dump<value_t>::value;
    };

    /**
    * Checks if the storage has the 'initialize' capability
    */
    template<typename... Args>
    static constexpr bool has_initialize()
    {
        return check_initialize<value_t, void, Args...>::value;
    };

    // Constructors and initialization
    PiercedStorage();
    PiercedStorage(std::size_t nFields, PiercedKernel<id_t> *kernel = nullptr, PiercedSyncMaster::SyncMode syncMode = PiercedSyncMaster::SYNC_MODE_DISABLED);

    ~PiercedStorage();

    std::size_t getFieldCount() const;

    void setKernel(PiercedKernel<id_t> *kernel, PiercedSyncMaster::SyncMode syncMode);
    void setKernel(const PiercedKernel<id_t> *kernel);
    void unsetKernel();
    const PiercedKernel<id_t> & getKernel() const;

    void setSyncMaster(PiercedSyncMaster *master, PiercedSyncMaster::SyncMode syncMode);
    void unsetSyncMaster();
    PiercedSyncMaster & getSyncMaster();
    const PiercedSyncMaster & getSyncMaster() const;

    // Methods that modify the container as a whole
    void swap(PiercedStorage &x) noexcept;
    void fill(const value_t &value);

    // Methos to access data stored in the container
    __PS_POINTER__ data();
    __PS_CONST_POINTER__ data() const;

    __PS_POINTER__ data(id_t id, std::size_t offset = 0);
    __PS_CONST_POINTER__ data(id_t id, std::size_t offset = 0) const;

    __PS_POINTER__ rawData(std::size_t pos, std::size_t offset = 0);
    __PS_CONST_POINTER__ rawData(std::size_t pos, std::size_t offset = 0) const;

    // Methods for editing the items using their id
    __PS_REFERENCE__ at(id_t id, std::size_t k = 0);
    __PS_CONST_REFERENCE__ at(id_t id, std::size_t k = 0) const;

    __PS_REFERENCE__ operator[](id_t id);
    __PS_CONST_REFERENCE__ operator[](id_t id) const;

    void copy(id_t id, value_t *values) const;
    void copy(id_t id, std::size_t nFields, value_t *values) const;
    void copy(id_t id, std::size_t nFields, std::size_t offset, value_t *values) const;

    void set(id_t id, const value_t &value);
    void set(id_t id, std::size_t k, const value_t &value);
    void set(id_t id, const value_t *values);
    void set(id_t id, std::size_t nFields, const value_t *values);
    void set(id_t id, std::size_t nFields, std::size_t offset, const value_t *values);

    // Methods for editing the items using their position
    __PS_REFERENCE__  rawAt(std::size_t pos, std::size_t offset = 0);
    __PS_CONST_REFERENCE__ rawAt(std::size_t pos, std::size_t offset = 0) const;

    void rawCopy(std::size_t pos, value_t *values) const;
    void rawCopy(std::size_t pos, std::size_t nFields, value_t *values) const;
    void rawCopy(std::size_t pos, std::size_t nFields, std::size_t offset, value_t *values) const;

    void rawSet(std::size_t pos, const value_t &value);
    void rawSet(std::size_t pos, std::size_t k, const value_t &value);
    void rawSet(std::size_t pos, const value_t *values);
    void rawSet(std::size_t pos, std::size_t nFields, const value_t *values);
    void rawSet(std::size_t pos, std::size_t nFields, std::size_t offset, const value_t *values);

    // Iterators
    iterator getIterator(const id_t &id) noexcept;
    const_iterator getConstIterator(const id_t &id) const noexcept;

    iterator getIteratorFromRawIndex(const std::size_t &pos) noexcept;
    const_iterator getConstIteratorFromRawIndex(const std::size_t &pos) const noexcept;

    iterator begin() noexcept;
    iterator end() noexcept;
    const_iterator begin() const noexcept;
    const_iterator end() const noexcept;
    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;

    raw_iterator rawBegin() noexcept;
    raw_iterator rawEnd() noexcept;
    raw_const_iterator rawBegin() const noexcept;
    raw_const_iterator rawEnd() const noexcept;
    raw_const_iterator rawCbegin() const noexcept;
    raw_const_iterator rawCend() const noexcept;

    // Dump and restore
    template<typename T = value_t, typename std::enable_if<std::is_pod<T>::value || PiercedStorage<T, id_t>::has_restore()>::type * = nullptr>
    void restore(std::istream &stream);

    template<typename T = value_t, typename std::enable_if<std::is_pod<T>::value || PiercedStorage<T, id_t>::has_dump()>::type * = nullptr>
    void dump(std::ostream &stream) const;

protected:
    // Methos for getting information on the storage
    std::size_t rawSize() const;

    // Methods that modify the container as a whole
    void swap(PiercedStorage &x, bool swapKernel) noexcept;

    // Methods for synchronizing the storage
    void commitSyncAction(const PiercedSyncAction &action);

    // Methos for updating the storage
    void rawReserve(std::size_t n);
    void rawShrinkToFit();

    void rawClear(bool release);
    void rawErase(std::size_t pos, std::size_t n);

    void rawSwap(std::size_t pos_first, std::size_t pos_second);
    void rawReorder(const std::vector<std::size_t> &permutations);

    void rawResize(std::size_t n, const value_t &value = value_t());

    template<typename... Args, typename std::enable_if<PiercedStorage<value_t>::template has_initialize<Args...>()>::type * = nullptr>
    void rawInitialize(std::size_t pos, Args&&... args);
    template<typename... Args, typename std::enable_if<PiercedStorage<value_t>::template has_initialize<Args...>()>::type * = nullptr>
    void rawInitialize(std::size_t pos, std::size_t k, Args&&... args);

    void rawInsert(std::size_t pos, std::size_t n, const value_t &value);
    void rawPushBack(const value_t &value);
    template<typename T = value_t, typename std::enable_if<!std::is_same<T, bool>::value>::type * = nullptr, typename... Args>
    void rawEmplace(std::size_t pos, Args&&... args);
    template<typename T = value_t, typename std::enable_if<std::is_same<T, bool>::value>::type * = nullptr>
    void rawEmplace(std::size_t pos, bool value = false);
    template<typename T = value_t, typename std::enable_if<!std::is_same<T, bool>::value>::type * = nullptr, typename... Args>
    void rawEmplaceBack(Args&&... args);
    template<typename T = value_t, typename std::enable_if<std::is_same<T, bool>::value>::type * = nullptr>
    void rawEmplaceBack(bool value = false);
    template<typename... Args>
    void rawEmreplace(std::size_t pos, Args&&... args);

    // Iterators
    iterator getIteratorFromPos(const std::size_t &pos) noexcept;
    const_iterator getConstIteratorFromPos(const std::size_t &pos) const noexcept;

private:
    const std::size_t m_nFields;
    container_t m_fields;

    const PiercedKernel<id_t> *m_kernel;
    PiercedSyncMaster *m_master;

    template<typename T = value_t, typename std::enable_if<std::is_pod<T>::value>::type * = nullptr>
    void restoreField(std::istream &stream, T &value);

    template<typename T, typename std::enable_if<PiercedStorage<T, id_t>::has_restore()>::type * = nullptr>
    void restoreField(std::istream &stream, T &object);

    template<typename T = value_t, typename std::enable_if<std::is_pod<T>::value>::type * = nullptr>
    void dumpField(std::ostream &stream, const T &value) const;

    template<typename T, typename std::enable_if<PiercedStorage<T, id_t>::has_dump()>::type * = nullptr>
    void dumpField(std::ostream &stream, const T &object) const;

};

}

// Templates
#include "piercedStorage.tpp"

#endif
