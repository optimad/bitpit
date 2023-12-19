/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

# ifndef __BITPIT_LEVELSET_CACHE_HPP__
# define __BITPIT_LEVELSET_CACHE_HPP__

#include "levelSetCommon.hpp"

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"
#include "bitpit_patchkernel.hpp"

#include <cassert>
#include <functional>
#include <memory>
#include <vector>

namespace bitpit {

class SendBuffer;
class RecvBuffer;

/*!
 * \ingroup levelset
 * \class LevelSetCache
 *
 * \brief The class LevelSetCache is the base class for defining caches.
 */

template<typename key_t>
class LevelSetCache
{
public:
    virtual ~LevelSetCache() = default;

    virtual std::unique_ptr<LevelSetCache<key_t>> clone() const = 0;

    virtual bool contains(const key_t &key) const = 0;

    template<typename Keys>
    std::size_t erase(const Keys &keys);
    virtual std::size_t erase(const key_t &key) = 0;

    virtual void reserve(std::size_t n) = 0;
    virtual void shrink_to_fit() = 0;
    virtual void clear() = 0;

    virtual bool isVolatile() const = 0;

    virtual void dump(std::ostream &stream) = 0;
    virtual void restore(std::istream &stream) = 0;

#if BITPIT_ENABLE_MPI
    virtual std::size_t getEntryBinarySize() const = 0;
    virtual void writeBuffer(const std::vector<key_t> &keys, SendBuffer &buffer) const = 0;
    virtual void readBuffer(const std::vector<key_t> &keys, RecvBuffer &buffer) = 0;
#endif

protected:
    LevelSetCache() = default;
};

/*!
 * \ingroup levelset
 * \class LevelSetValueCacheEntry
 *
 * \brief The class LevelSetValueCacheEntry allows to get read-only access to a value stored in
 * the cache.
 *
 * The value stored in the cache entry may not point directly to the valued contained in the
 * cache. For example, for boolean values, the cache entry will point to an internally defined
 * boolean variable. In this case, the value of the internal variable will be set equal to the
 * value of the cache during the construction of the entry (this means that, if the cached values
 * changes, the entry will not be updated). This is done to make it possible for the cache entry
 * to deal with std::vector containing bool values.
 */

template<typename value_t>
class LevelSetValueCacheBaseEntry
{
public:
    LevelSetValueCacheBaseEntry(bool valid);

    bool isValid() const;

protected:
    const static value_t m_dummyValue;

private:
    bool m_valid;
};

template<typename value_t>
class LevelSetValueCacheEntry : public LevelSetValueCacheBaseEntry<value_t>
{
public:
    LevelSetValueCacheEntry();
    LevelSetValueCacheEntry(const value_t &value);

    const value_t & operator*() const;

private:
    std::reference_wrapper<const value_t> m_value;

};

template<>
struct LevelSetValueCacheEntry<bool> : public LevelSetValueCacheBaseEntry<bool>
{
    LevelSetValueCacheEntry();
    explicit LevelSetValueCacheEntry(const bool &value);
    explicit LevelSetValueCacheEntry(const std::vector<bool>::reference &value);

    bool operator*() const;

private:
    static const std::vector<bool> m_dummyVector;

    const bool &m_value;

    bool m_useVectorValue;
    std::vector<bool>::const_reference m_vectorValue;
};

/*!
 * \ingroup levelset
 * \class LevelSetValueCache
 *
 * \brief The class LevelSetCache is the base class for defining caches that store values.
 */

template<typename key_t, typename value_t>
class LevelSetValueCache : public LevelSetCache<key_t>
{
public:
    typedef LevelSetValueCacheEntry<value_t> Entry;

    virtual Entry findEntry(const key_t &key) const = 0;
    virtual Entry insertEntry(const key_t &key, const value_t &value) = 0;
    virtual Entry insertEntry(const key_t &key, value_t &&value) = 0;

};

/*!
 * \ingroup levelset
 * \class LevelSetContainerBaseCache
 *
 * \brief The class LevelSetContainerBaseCache is the base class for defining caches that stores
 * the values inside a container owned by them.
 */

template<typename key_t, typename container_t, typename value_t = typename container_t::value_type, typename reference_t = value_t &, typename const_reference_t = const value_t &>
class LevelSetContainerBaseCache : public LevelSetValueCache<key_t, value_t>
{
public:
    typedef typename LevelSetValueCache<key_t, value_t>::Entry Entry;

    typedef container_t container_type;
    typedef const key_t &key_type;
    typedef value_t value_type;

    typedef reference_t reference;
    typedef const_reference_t const_reference;

    typedef typename container_t::iterator iterator;
    typedef typename container_t::const_iterator const_iterator;

    iterator begin();
    iterator end();

    const_iterator cbegin() const;
    const_iterator cend() const;

    const_iterator begin() const;
    const_iterator end() const;

    bool contains(const key_t &key) const override;

    Entry findEntry(const key_t &key) const override;
    Entry insertEntry(const key_t &key, const value_t &value) override;
    Entry insertEntry(const key_t &key, value_t &&value) override;

    reference at(const key_t &key);
    const_reference at(const key_t &key) const;

    reference operator[](const key_t &key);
    const_reference operator[](const key_t &key) const;

    virtual iterator find(const key_t &key) = 0;
    virtual const_iterator find(const key_t &key) const = 0;

    virtual iterator insert(const key_t &key, const value_t &value) = 0;
    virtual iterator insert(const key_t &key, value_t &&value) = 0;

#if BITPIT_ENABLE_MPI
    std::size_t getEntryBinarySize() const override;
    void writeBuffer(const std::vector<key_t> &ids, SendBuffer &buffer) const override;
    void readBuffer(const std::vector<key_t> &ids, RecvBuffer &buffer) override;
#endif

protected:
    container_t m_container;

    template<typename... Args>
    LevelSetContainerBaseCache(Args&&... args);

    virtual key_t getKey(const const_iterator &itr) const = 0;
    virtual reference getValue(const iterator &itr) const = 0;
    virtual const_reference getValue(const const_iterator &itr) const = 0;
};

/*!
 * \ingroup levelset
 * \class LevelSetContainerCache
 *
 * \brief The class LevelSetContainerCache is the class for defining caches that stores
 * the values inside a container owned by them.
 */

template<typename key_t, typename container_t>
class LevelSetContainerCache : public LevelSetContainerBaseCache<key_t, container_t>
{
public:
    typedef LevelSetContainerBaseCache<key_t, container_t> Base;

    typedef typename Base::key_type key_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::container_type container_type;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;

};

template<typename key_t, typename value_t>
class LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>> : public LevelSetContainerBaseCache<key_t,
                                                                                                            std::unordered_map<key_t, value_t>,
                                                                                                            typename std::unordered_map<key_t, value_t>::mapped_type>
{
public:
    typedef LevelSetContainerBaseCache<key_t,
                                       std::unordered_map<key_t, value_t>,
                                       typename std::unordered_map<key_t, value_t>::mapped_type>
            Base;

    typedef typename Base::key_type key_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::container_type container_type;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;

    LevelSetContainerCache(std::size_t capacity = 0);

    std::unique_ptr<LevelSetCache<key_t>> clone() const override;

    void reserve(std::size_t n) override;
    void shrink_to_fit() override;
    void clear() override;

    iterator find(const key_t &key) override;
    const_iterator find(const key_t &key) const override;

    iterator insert(const key_t &key, const value_t &value) override;
    iterator insert(const key_t &key, value_t &&value) override;
    std::size_t erase(const key_t &key) override;

    bool isVolatile() const override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

private:
    key_t getKey(const const_iterator &itr) const override;
    reference getValue(const iterator &itr) const override;
    const_reference getValue(const const_iterator &itr) const override;
};


template<typename key_t, typename value_t>
class LevelSetContainerCache<key_t, std::vector<value_t>> : public LevelSetContainerBaseCache<key_t,
                                                                                              std::vector<value_t>,
                                                                                              typename std::vector<value_t>::value_type,
                                                                                              typename std::vector<value_t>::reference,
                                                                                              typename std::vector<value_t>::const_reference>
{
public:
    typedef LevelSetContainerBaseCache<key_t,
                                       std::vector<value_t>,
                                       typename std::vector<value_t>::value_type,
                                       typename std::vector<value_t>::reference,
                                       typename std::vector<value_t>::const_reference>
            Base;

    typedef typename Base::key_type key_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::container_type container_type;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;

    LevelSetContainerCache(std::size_t capacity = 0);

    std::unique_ptr<LevelSetCache<key_t>> clone() const override;

    void reserve(std::size_t n) override;
    void shrink_to_fit() override;
    void clear() override;

    iterator find(const key_t &key) override;
    const_iterator find(const key_t &key) const override;

    iterator insert(const key_t &key, const value_t &value) override;
    iterator insert(const key_t &key, value_t &&value) override;
    std::size_t erase(const key_t &key) override;

    bool isVolatile() const override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

protected:
    key_t getKey(const const_iterator &itr) const override;
    reference getValue(const iterator &itr) const override;
    const_reference getValue(const const_iterator &itr) const override;

private:
    std::vector<bool> m_isCached;
};

template<typename key_t, typename value_t>
class LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>> : public LevelSetContainerBaseCache<key_t,
                                                                                                       PiercedVector<value_t, key_t>,
                                                                                                       typename PiercedVector<value_t, key_t>::value_type,
                                                                                                       typename PiercedVector<value_t, key_t>::reference,
                                                                                                       typename PiercedVector<value_t, key_t>::const_reference>
{
public:
    typedef LevelSetContainerBaseCache<key_t,
                                       PiercedVector<value_t, key_t>,
                                       typename PiercedVector<value_t, key_t>::value_type,
                                       typename PiercedVector<value_t, key_t>::reference,
                                       typename PiercedVector<value_t, key_t>::const_reference>
            Base;

    typedef typename Base::key_type key_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::container_type container_type;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;

    LevelSetContainerCache(std::size_t capacity = 0);

    std::unique_ptr<LevelSetCache<key_t>> clone() const override;

    void reserve(std::size_t n) override;
    void shrink_to_fit() override;
    void clear() override;

    iterator find(const key_t &key) override;
    const_iterator find(const key_t &key) const override;

    iterator insert(const key_t &key, const value_t &value) override;
    iterator insert(const key_t &key, value_t &&value) override;
    std::size_t erase(const key_t &key) override;

    bool isVolatile() const override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

protected:
    key_t getKey(const const_iterator &itr) const override;
    reference getValue(const iterator &itr) const override;
    const_reference getValue(const const_iterator &itr) const override;
};

template<typename key_t, typename value_t>
class LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>> : public LevelSetContainerBaseCache<key_t,
                                                                                                        PiercedStorage<value_t, key_t>,
                                                                                                        typename PiercedStorage<value_t, key_t>::value_type,
                                                                                                        typename PiercedStorage<value_t, key_t>::reference,
                                                                                                        typename PiercedStorage<value_t, key_t>::const_reference>
{
public:
    typedef LevelSetContainerBaseCache<key_t,
                                       PiercedStorage<value_t, key_t>,
                                       typename PiercedStorage<value_t, key_t>::value_type,
                                       typename PiercedStorage<value_t, key_t>::reference,
                                       typename PiercedStorage<value_t, key_t>::const_reference>
            Base;

    typedef typename Base::key_type key_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::container_type container_type;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;

    LevelSetContainerCache(const PiercedKernel<key_t> *m_kernel, PiercedSyncMaster::SyncMode syncMode);

    std::unique_ptr<LevelSetCache<key_t>> clone() const override;

    void reserve(std::size_t n) override;
    void shrink_to_fit() override;
    void clear() override;

    iterator find(const key_t &key) override;
    const_iterator find(const key_t &key) const override;

    iterator insert(const key_t &key, const value_t &value) override;
    iterator insert(const key_t &key, value_t &&value) override;
    std::size_t erase(const key_t &key) override;

    bool isVolatile() const override;

    void dump(std::ostream &stream) override;
    void restore(std::istream &stream) override;

protected:
    key_t getKey(const const_iterator &itr) const override;
    reference getValue(const iterator &itr) const override;
    const_reference getValue(const const_iterator &itr) const override;

private:
    PiercedStorage<bool, key_t> m_isCached;
};

/*!
 * \ingroup levelset
 * \class LevelSetCacheFactory
 *
 * \brief The class LevelSetCacheFactory provides basic functionalities for cache factories.
 */

template<typename key_t>
class LevelSetCacheFactory
{
public:
    virtual ~LevelSetCacheFactory() = default;

    virtual std::unique_ptr<LevelSetCache<key_t>> create() const = 0;

protected:
    LevelSetCacheFactory() = default;
};

/*!
 * \ingroup levelset
 * \class LevelSetValueCacheFactory
 *
 * \brief The class LevelSetValueCacheFactory provides basic functionalities for cache factories.
 */

template<typename key_t, typename value_t>
class LevelSetValueCacheFactory : public LevelSetCacheFactory<key_t>
{
protected:
    LevelSetValueCacheFactory() = default;

};

/*!
 * \ingroup levelset
 * \class LevelSetContainerBaseCacheFactory
 *
 * \brief The class LevelSetContainerBaseCacheFactory provides basic functionalities for cache factories.
 */

template<typename key_t, typename container_t>
class LevelSetContainerBaseCacheFactory : public LevelSetValueCacheFactory<key_t, typename LevelSetContainerCache<key_t, container_t>::value_type>
{
protected:
    LevelSetContainerBaseCacheFactory() = default;

};

/*!
 * \ingroup levelset
 * \class LevelSetContainerCacheFactory
 *
 * \brief The class LevelSetContainerCacheFactory provides basic functionalities for cache factories.
 */

template<typename key_t, typename container_t>
class LevelSetContainerCacheFactory : public LevelSetContainerBaseCacheFactory<key_t, container_t>
{
public:
    typedef LevelSetContainerBaseCacheFactory<key_t, container_t> Base;

    std::unique_ptr<LevelSetCache<key_t>> create() const override;

};

template<typename key_t, typename value_t>
class LevelSetContainerCacheFactory<key_t, std::unordered_map<key_t ,value_t>> : public LevelSetContainerBaseCacheFactory<key_t, std::unordered_map<key_t, value_t>>
{
public:
    typedef LevelSetContainerBaseCacheFactory<key_t, std::unordered_map<key_t, value_t>> Base;

    std::unique_ptr<LevelSetCache<key_t>> create() const override;

};

template<typename key_t, typename value_t>
class LevelSetContainerCacheFactory<key_t, std::vector<value_t>> : public LevelSetContainerBaseCacheFactory<key_t, std::vector<value_t>>
{
public:
    typedef LevelSetContainerBaseCacheFactory<key_t, std::vector<value_t>> Base;

    std::unique_ptr<LevelSetCache<key_t>> create() const override;

};

template<typename key_t, typename value_t>
class LevelSetContainerCacheFactory<key_t, PiercedVector<value_t, key_t>> : public LevelSetContainerBaseCacheFactory<key_t, PiercedVector<value_t, key_t>>
{
public:
    typedef LevelSetContainerBaseCacheFactory<key_t, PiercedVector<value_t, key_t>> Base;

    std::unique_ptr<LevelSetCache<key_t>> create() const override;

};

template<typename key_t, typename value_t>
class LevelSetContainerCacheFactory<key_t, PiercedStorage<value_t, key_t>> : public LevelSetContainerBaseCacheFactory<key_t, PiercedStorage<value_t, key_t>>
{
public:
    typedef LevelSetContainerBaseCacheFactory<key_t, PiercedStorage<value_t, key_t>> Base;

    LevelSetContainerCacheFactory(const PiercedKernel<key_t> *kernel, PiercedSyncMaster::SyncMode syncMode);

    std::unique_ptr<LevelSetCache<key_t>> create() const override;

protected:
    const PiercedKernel<key_t> *m_kernel;
    PiercedSyncMaster::SyncMode m_syncMode;
};

/*!
 * \interface LevelSetCacheCollection
 * \ingroup levelset
 *
 * \brief The class LevelSetCacheCollection allows to store a collection of caches.
 *
 * Creation of the caches is lazy: they will be created the first time they are needed. The
 * creation of the caches is handled by the cache factory defined when the cache was added
 * to the collection. This means that, when iterating the collection, it is necessary to
 * check if the single caches have been allocated.
 */

template<typename key_t>
class LevelSetCacheCollection {

public:
    class Item
    {
    public:
        Item(const std::shared_ptr<LevelSetCacheFactory<key_t>> &factory = nullptr);

        Item(const Item &other);
        Item(Item &&other) = default;

        Item & operator=(const Item &other);
        Item & operator=(Item &&other) = default;

        bool hasFactory() const;

        bool hasCache() const;
        LevelSetCache<key_t> * createCache();
        void destroyCache();

        LevelSetCache<key_t> * getCache(bool allowCreation = true);
        const LevelSetCache<key_t> * getCache(bool allowCreation = true) const;

        template<typename value_t>
        LevelSetValueCache<key_t, value_t> * getCache(bool allowCreation = true);
        template<typename value_t>
        const LevelSetValueCache<key_t, value_t> * getCache(bool allowCreation = true) const;

    private:
        std::shared_ptr<LevelSetCacheFactory<key_t>> m_factory; //!< Factory that will create the caches
        mutable std::unique_ptr<LevelSetCache<key_t>> m_cache; //!< Cache

        LevelSetCache<key_t> * _createCache() const;
    };

    typedef key_t key_type;

    typedef std::vector<Item> Caches;

    typedef typename Caches::iterator iterator;
    typedef typename Caches::const_iterator const_iterator;

    static const std::size_t NULL_CACHE_ID;

    LevelSetCacheCollection() = default;
    LevelSetCacheCollection(const LevelSetCacheCollection &other);
    LevelSetCacheCollection(LevelSetCacheCollection &&other) = default;

    virtual ~LevelSetCacheCollection() = default;

    virtual std::unique_ptr<LevelSetCacheCollection<key_t>> clone() const;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator cbegin() const;
    const_iterator cend() const;

    std::size_t size() const;

    Item & operator[](std::size_t index);
    const Item & operator[](std::size_t index) const;

    Item & at(std::size_t index);
    const Item & at(std::size_t index) const;

    template<typename container_t, typename... Args>
    std::size_t insert(Args&&... args);
    void erase(std::size_t index);

    void clear();

protected:
    Caches m_caches; //!< Caches owned by the collection

};

class ElementCacheCollection : public LevelSetCacheCollection<long>
{
public:
    typedef LevelSetCacheCollection<long> Base;

    typedef Base::key_type key_type;

    typedef LevelSetCache<key_type> Cache;

    template<typename value_t>
    using ValueCache = LevelSetValueCache<key_type, value_t>;

    ElementCacheCollection(const PiercedKernel<key_type> *kernel);

    template<typename container_t, typename... Args, typename std::enable_if<std::is_same<bitpit::PiercedStorage<typename container_t::value_type>, container_t>::value>::type * = nullptr>
    std::size_t insert(Args&&... args);
    template<typename container_t, typename... Args, typename std::enable_if<!std::is_same<bitpit::PiercedStorage<typename container_t::value_type>, container_t>::value>::type * = nullptr>
    std::size_t insert(Args&&... args);

protected:
    const PiercedKernel<long> *m_kernel;
};

}

// Include template implementations
#include "levelSetCache.tpp"

#endif 
