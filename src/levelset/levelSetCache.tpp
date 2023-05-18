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

# ifndef __BITPIT_LEVELSET_CACHE_TPP__
# define __BITPIT_LEVELSET_CACHE_TPP__

namespace bitpit {

/*!
 * Dummy value to be used for invalid entries.
 */
template<typename value_t>
const value_t LevelSetValueCacheBaseEntry<value_t>::m_dummyValue = value_t{};

/*!
 * Constructor.
 *
 * \param valid if set to true, the cache entry is considered valid otherwise the cache entry is
 * \considered invalid. Valid cache entry can be deferenced to access the cached value. Trying to
 * deference an invalid entry results in undefined behaviour.
 */
template<typename value_t>
LevelSetValueCacheBaseEntry<value_t>::LevelSetValueCacheBaseEntry(bool valid)
    : m_valid(valid)
{
}

/*!
 * Check if the cache entry is valid.
 *
 * Valid cache entry can be deferenced to access the cached value. Trying to deference an invalid
 * entry results in undefined behaviour.
 */
template<typename value_t>
bool LevelSetValueCacheBaseEntry<value_t>::isValid() const
{
    return m_valid;
}

/*!
 * Constructor for invalid entries.
 */
template<typename value_t>
LevelSetValueCacheEntry<value_t>::LevelSetValueCacheEntry()
    : LevelSetValueCacheBaseEntry<value_t>(false),
      m_value(std::cref(LevelSetValueCacheBaseEntry<value_t>::m_dummyValue))
{
}

/*!
 * Constructor for valid entries.
 *
 * \param value is the cached value
 */
template<typename value_t>
LevelSetValueCacheEntry<value_t>::LevelSetValueCacheEntry(const value_t &value)
    : LevelSetValueCacheBaseEntry<value_t>(true),
      m_value(value)
{
}

/*!
 * Get the stored value.
 *
 * Trying to deference an invalid entry results in undefined behaviour.
 */
template<typename value_t>
const value_t & LevelSetValueCacheEntry<value_t>::operator*() const
{
    return m_value;
}

/*!
 * Erase the entries associated to the given keys from the specified cache.
 *
 * \param keys are the keys whose entries will be erased
 */
template<typename key_t>
template<typename Keys>
void LevelSetCache<key_t>::erase(const Keys &keys)
{
    for (const key_t &key : keys) {
        erase(key);
    }
}

/*!
 * Constructor.
 *
 * \param args are the arguments that will be used to create the container
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
template<typename... Args>
LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::LevelSetContainerBaseCache(Args&&... args)
    : m_container(std::forward<Args>(args)...)
{
}

/*!
 * Get an iterator pointing to the first entry.
 *
 * \result An iterator pointing to the first entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::iterator LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::begin()
{
    return m_container.begin();
}

/*!
 * Get an iterator referring to the past-the-end entry.
 *
 * \result An iterator referring to the past-the-end entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::iterator LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::end()
{
    return m_container.end();
}

/*!
 * Get a constant iterator pointing to the first entry.
 *
 * \result A constant iterator pointing to the first entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::const_iterator LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::begin() const
{
    return m_container.cbegin();
}

/*!
 * Get a constant iterator referring to the past-the-end entry.
 *
 * \result a constant iterator referring to the past-the-end entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::const_iterator LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::end() const
{
    return m_container.cend();
}

/*!
 * Get a constant iterator pointing to the first entry.
 *
 * \result A constant iterator pointing to the first entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::const_iterator LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::cbegin() const
{
    return m_container.cbegin();
}

/*!
 * Get a constant iterator referring to the past-the-end entry.
 *
 * \result a constant iterator referring to the past-the-end entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::const_iterator LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::cend() const
{
    return m_container.cend();
}

/*!
 * Check if the cache contains an entry for the specified key.
 *
 * \param key is key associated with the entry
 * \result True if the cache contains an entry for the specified key, false otherwise.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
bool LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::contains(const key_t &key) const
{
    return (find(key) != m_container.cend());
}

/*!
 * Find the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an invalid entry will
 * be returned.
 *
 * \param key is key associated with the entry
 * \result The entry associated with the specified key or an invalid entry if the cache doesn't
 * contain an actual entry associated with the key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetValueCache<key_t, value_t>::Entry LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::findEntry(const key_t &key)
{
    iterator itr = find(key);
    if (itr != end()) {
        return Entry(getValue(itr));
    } else {
        return Entry();
    }
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result The entry associated with the key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetValueCache<key_t, value_t>::Entry LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::insertEntry(const key_t &key, const value_t &value)
{
    iterator itr = insert(key, value);

    return Entry(getValue(itr));
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result The entry associated with the key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetValueCache<key_t, value_t>::Entry LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::insertEntry(const key_t &key, value_t &&value)
{
    iterator itr = insert(key, std::move(value));

    return Entry(getValue(itr));
}

/*!
 * Get a reference to the entry associated with the specified key.
 *
 * Attempting to access a non-existent entry results in undefined behavior.
 *
 * \param key is key associated with the entry
 * \result A reference to the entry associated with the specified key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::reference LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::at(const key_t &key)
{
    return *(find(key));
}

/*!
 * Get a constant reference to the entry associated with the specified key.
 *
 * Attempting to access a non-existent entry results in undefined behavior.
 *
 * \param key is key associated with the entry
 * \result A constant reference to the entry associated with the specified key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::const_reference LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::at(const key_t &key) const
{
    return getValue(find(key));
}

/*!
 * Get a reference to the entry associated with the specified key.
 *
 * Attempting to access a non-existent entry results in undefined behavior.
 *
 * \param key is key associated with the entry
 * \result A reference to the entry associated with the specified key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::reference LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::operator[](const key_t &key)
{
    return getValue(find(key));
}

/*!
 * Get a constant reference to the entry associated with the specified key.
 *
 * Attempting to access a non-existent entry results in undefined behavior.
 *
 * \param key is key associated with the entry
 * \result A constant reference to the entry associated with the specified key.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::const_reference LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::operator[](const key_t &key) const
{
    return getValue(find(key));
}

#if BITPIT_ENABLE_MPI
/*!
 * Get the size, expressed in bytes, of an entry.
 *
 * \result The size, expressed in bytes, of an entry.
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
std::size_t LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::getEntryBinarySize() const
{
    return sizeof(bool) + sizeof(value_t);
}

/*!
 * Write the specified entries to the given buffer.
 *
 * \param[in] keys is the list of keys whose data need to be send
 * \param[in,out] buffer is the buffer to write to
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
void LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::writeBuffer(const std::vector<key_t> &keys, SendBuffer &buffer)
{
    for (const key_t &key : keys) {
        const_iterator itr = find(key);
        if (itr != cend()) {
            buffer << true;
            buffer << *itr;
        } else {
            buffer << false;
            buffer << value_t{};
        }
    }
}

/*!
 * Read the specified entries from the given buffer.
 *
 * \param keys is the list of keys whose data need to be received
 * @param[in,out] buffer is the buffer to read from
 */
template<typename key_t, typename container_t, typename value_t, typename reference_t, typename const_reference_t>
void LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::readBuffer(const std::vector<key_t> &keys, RecvBuffer &buffer)
{
    for (const key_t &key : keys) {
        bool isCached;
        buffer >> isCached;

        value_t value;
        buffer >> value;
        if (isCached) {
            insert(key, value);
        }
    }
}
#endif

/*!
 * Constructor.
 *
 * \param capacity is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::LevelSetContainerCache(std::size_t capacity)
    : Base()
{
    Base::m_container.reserve(capacity);
}

/*!
 * Clone the cache.
 *
 * \result The cloned cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::clone() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>(*this));
}

/*!
 * Requests that the cache capacity be at least enough to contain n entries.
 *
 * \param n is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::reserve(std::size_t n)
{
    Base::m_container.reserve(n);
}

/*!
 * Requests the removal of unused capacity.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::shrink_to_fit()
{
    // Nothing to do
}

/*!
 * Clear the cache.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::clear()
{
    Base::m_container.clear();
}

/*!
 * Get an iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result An iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::iterator LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::find(const key_t &key)
{
    return Base::m_container.find(key);
}

/*!
 * Get a constant iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result A constant iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::const_iterator LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::find(const key_t &key) const
{
    return Base::m_container.find(key);
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::iterator LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::insert(const key_t &key, const value_t &value)
{
    auto itr = Base::m_container.find(key);
    if (itr != Base::m_container.end()) {
        itr->second = value;
    } else {
        itr = Base::m_container.insert({key, value}).first;
    }

    return itr;
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::iterator LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::insert(const key_t &key, value_t &&value)
{
    auto itr = Base::m_container.find(key);
    if (itr != Base::m_container.end()) {
        itr->second = std::move(value);
    } else {
        itr = Base::m_container.insert({key, std::move(value)}).first;
    }

    return itr;
}

/*!
 * Erase the cache entry associated with the specified key.
 *
 * \param key is key associated with the entry
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::erase(const key_t &key)
{
    Base::m_container.erase(key);
}

/*!
 * Write the cache to the specified stream.
 *
 * \param stream is the stream to write to
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::dump(std::ostream &stream)
{
    utils::binary::write(stream, Base::m_container.size());
    for (const auto &entry : Base::m_container) {
        utils::binary::write(stream, entry.first);
        utils::binary::write(stream, entry.second);
    }
}

/*!
 * Restore the cache from the specified stream.
 *
 * \param stream is the stream to read from
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::restore(std::istream &stream)
{
    Base::m_container.clear();

    std::size_t nEntries;
    utils::binary::read(stream, nEntries);
    Base::m_container.reserve(nEntries);

    for (std::size_t i = 0; i < nEntries; ++i) {
        key_t key;
        utils::binary::read(stream, key);

        value_t value;
        utils::binary::read(stream, value);

        Base::m_container[key] = std::move(value);
    }
}

/*!
 * Get the key associated to the iterator.
 *
 * \param itr is the iterator
 * \result The key associated to the iterator.
 */
template<typename key_t, typename value_t>
key_t LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::getKey(const const_iterator &itr) const
{
    return itr->first;
}

/*!
 * Get a reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::reference LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::getValue(const iterator &itr) const
{
    return itr->second;
}

/*!
 * Get a constant reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A constant reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::const_reference LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::getValue(const const_iterator &itr) const
{
    return itr->second;
}

/*!
 * Constructor.
 *
 * \param capacity is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
LevelSetContainerCache<key_t, std::vector<value_t>>::LevelSetContainerCache(std::size_t capacity)
    : Base()
{
    Base::m_container.reserve(capacity);
    m_isCached.reserve(capacity);
}

/*!
 * Clone the cache.
 *
 * \result The cloned cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCache<key_t, std::vector<value_t>>::clone() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, std::vector<value_t>>(*this));
}

/*!
 * Requests that the cache capacity be at least enough to contain n entries.
 *
 * \param n is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::vector<value_t>>::reserve(std::size_t n)
{
    Base::m_container.reserve(n);
}

/*!
 * Requests the removal of unused capacity.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::vector<value_t>>::shrink_to_fit()
{
    if (m_isCached.empty()) {
        return;
    }

    std::size_t capacity = m_isCached.size();
    while (capacity > 0 && !m_isCached[capacity - 1]) {
        --capacity;
    }

    Base::m_container.resize(capacity);
    m_isCached.resize(capacity);
}

/*!
 * Clear the cache.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::vector<value_t>>::clear()
{
    std::fill(m_isCached.begin(), m_isCached.end(), false);
}

/*!
 * Get an iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result An iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::vector<value_t>>::iterator LevelSetContainerCache<key_t, std::vector<value_t>>::find(const key_t &key)
{
    if (key < static_cast<key_t>(Base::m_container.size()) && m_isCached[key]) {
        return (Base::m_container.begin() + key);
    } else {
        return Base::m_container.end();
    }
}

/*!
 * Get a constant iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result A constant iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::vector<value_t>>::const_iterator LevelSetContainerCache<key_t, std::vector<value_t>>::find(const key_t &key) const
{
    if (key < static_cast<key_t>(Base::m_container.size()) && m_isCached[key]) {
        return (Base::m_container.cbegin() + key);
    } else {
        return Base::m_container.cend();
    }
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::vector<value_t>>::iterator LevelSetContainerCache<key_t, std::vector<value_t>>::insert(const key_t &key, const value_t &value)
{
    if (key >= static_cast<key_t>(Base::m_container.size())) {
        Base::m_container.resize(key + 1);
        m_isCached.resize(key + 1);
    }

    auto itr = Base::m_container.begin() + key;

    *itr             = value;
    m_isCached[key]  = true;

    return itr;
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::vector<value_t>>::iterator LevelSetContainerCache<key_t, std::vector<value_t>>::insert(const key_t &key, value_t &&value)
{
    if (key >= static_cast<key_t>(Base::m_container.size())) {
        Base::m_container.resize(key + 1);
        m_isCached.resize(key + 1);
    }

    auto itr = Base::m_container.begin() + key;

    *itr             = std::move(value);
    m_isCached[key]  = true;

    return itr;
}

/*!
 * Erase the cache entry associated with the specified key.
 *
 * \param key is key associated with the entry
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::vector<value_t>>::erase(const key_t &key)
{
    m_isCached[key] = false;
}

/*!
 * Write the cache to the specified stream.
 *
 * \param stream is the stream to write to
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::vector<value_t>>::dump(std::ostream &stream)
{
    utils::binary::write(stream, Base::m_container);
    utils::binary::write(stream, m_isCached);
}

/*!
 * Restore the cache from the specified stream.
 *
 * \param stream is the stream to read from
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, std::vector<value_t>>::restore(std::istream &stream)
{
    utils::binary::read(stream, Base::m_container);
    utils::binary::read(stream, m_isCached);
}

/*!
 * Get the key associated to the iterator.
 *
 * \param itr is the iterator
 * \result The key associated to the iterator.
 */
template<typename key_t, typename value_t>
key_t LevelSetContainerCache<key_t, std::vector<value_t>>::getKey(const const_iterator &itr) const
{
    return (itr - Base::m_container.begin());
}

/*!
 * Get a reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::vector<value_t>>::reference LevelSetContainerCache<key_t, std::vector<value_t>>::getValue(const iterator &itr) const
{
    return *itr;
}

/*!
 * Get a constant reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A constant reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, std::vector<value_t>>::const_reference LevelSetContainerCache<key_t, std::vector<value_t>>::getValue(const const_iterator &itr) const
{
    return *itr;
}


/*!
 * Constructor.
 *
 * \param capacity is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::LevelSetContainerCache(std::size_t capacity)
    : Base()
{
    Base::m_container.reserve(capacity);
}

/*!
 * Clone the cache.
 *
 * \result The cloned cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::clone() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>(*this));
}

/*!
 * Requests that the cache capacity be at least enough to contain n entries.
 *
 * \param n is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::reserve(std::size_t n)
{
    Base::m_container.reserve(n);
}

/*!
 * Requests the removal of unused capacity.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::shrink_to_fit()
{
    Base::m_container.squeeze();
}

/*!
 * Clear the cache.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::clear()
{
    Base::m_container.clear();
}

/*!
 * Get an iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result An iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::iterator LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::find(const key_t &key)
{
    return Base::m_container.find(key);
}

/*!
 * Get a constant iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result A constant iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::const_iterator LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::find(const key_t &key) const
{
    return Base::m_container.find(key);
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::iterator LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::insert(const key_t &key, const value_t &value)
{
    auto itr = Base::m_container.find(key);
    if (itr == Base::m_container.end()) {
        itr = Base::m_container.reclaim(key);
    }

    *itr = value;

    return itr;
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::iterator LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::insert(const key_t &key, value_t &&value)
{
    auto itr = Base::m_container.find(key);
    if (itr == Base::m_container.end()) {
        itr = Base::m_container.reclaim(key);
    }

    *itr = std::move(value);

    return itr;
}

/*!
 * Erase the cache entry associated with the specified key.
 *
 * \param key is key associated with the entry
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::erase(const key_t &key)
{
    Base::m_container.erase(key);
}

/*!
 * Write the cache to the specified stream.
 *
 * \param stream is the stream to write to
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::dump(std::ostream &stream)
{
    Base::m_container.dump(stream);
}

/*!
 * Restore the cache from the specified stream.
 *
 * \param stream is the stream to read from
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::restore(std::istream &stream)
{
    Base::m_container.restore(stream);
}

/*!
 * Get the key associated to the iterator.
 *
 * \param itr is the iterator
 * \result The key associated to the iterator.
 */
template<typename key_t, typename value_t>
key_t LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::getKey(const const_iterator &itr) const
{
    return itr.getId();
}

/*!
 * Get a reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::reference LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::getValue(const iterator &itr) const
{
    return *itr;
}

/*!
 * Get a constant reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A constant reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::const_reference LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::getValue(const const_iterator &itr) const
{
    return *itr;
}

/*!
 * Constructor.
 *
 * \param kernel is the kernel associated with the entries
 * \param syncMode is the synchronization mode that will be used for the storage
 */
template<typename key_t, typename value_t>
LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::LevelSetContainerCache(PiercedKernel<key_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
    : Base(1, kernel, syncMode),
      m_isCached(1, kernel, syncMode)
{
    m_isCached.fill(false);
}

/*!
 * Clone the cache.
 *
 * \result The cloned cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::clone() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>(*this));
}

/*!
 * Requests that the cache capacity be at least enough to contain n entries.
 *
 * \param n is the number of entries the cache should be able to contain
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::reserve(std::size_t n)
{
    BITPIT_UNUSED(n);

    // Nothing to do
}

/*!
 * Requests the removal of unused capacity.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::shrink_to_fit()
{
    // Nothing to do
}

/*!
 * Clear the cache.
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::clear()
{
    m_isCached.fill(false);
}

/*!
 * Get an iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result An iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::iterator LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::find(const key_t &key)
{
    iterator itr = Base::m_container.find(key);
    std::size_t rawId = itr.getRawIndex();

    if (m_isCached.rawAt(rawId)) {
        return itr;
    } else {
        return Base::m_container.end();
    }
}

/*!
 * Get a constant iterator pointing to the entry associated with the specified key.
 *
 * If the cache doesn't contain an entry associated with the specified key, an iterator pointing
 * to the container's end will be returned.
 *
 * \param key is key associated with the entry
 * \result A constant iterator pointing to the entry associated with the specified key.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::const_iterator LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::find(const key_t &key) const
{
    const_iterator itr = Base::m_container.find(key);
    std::size_t rawId = itr.getRawIndex();

    if (m_isCached.rawAt(rawId)) {
        return itr;
    } else {
        return Base::m_container.cend();
    }
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::iterator LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::insert(const key_t &key, const value_t &value)
{
    auto itr = Base::m_container.find(key);
    std::size_t rawId = itr.getRawIndex();

    *itr                    = value;
    m_isCached.rawAt(rawId) = true;

    return itr;
}

/*!
 * Insert a cache entry for the specified key.
 *
 * If an entry associated with the same key is already in the cache, the existing cached value
 * will be replaced by the one passed to the function.
 *
 * \param key is key associated with the entry
 * \param value is the value that will be inserted into the cache
 * \result An iterator pointing to the inserted cache entry.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::iterator LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::insert(const key_t &key, value_t &&value)
{
    auto itr = Base::m_container.find(key);
    std::size_t rawId = itr.getRawIndex();

    *itr                    = std::move(value);
    m_isCached.rawAt(rawId) = true;

    return itr;
}

/*!
 * Erase the cache entry associated with the specified key.
 *
 * \param key is key associated with the entry
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::erase(const key_t &key)
{
    auto itr = Base::m_container.find(key);
    std::size_t rawId = itr.getRawIndex();

    m_isCached.rawAt(rawId) = false;
}

/*!
 * Write the cache to the specified stream.
 *
 * \param stream is the stream to write to
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::dump(std::ostream &stream)
{
    Base::m_container.dump(stream);
    m_isCached.dump(stream);
}

/*!
 * Restore the cache from the specified stream.
 *
 * \param stream is the stream to read from
 */
template<typename key_t, typename value_t>
void LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::restore(std::istream &stream)
{
    Base::m_container.restore(stream);
    m_isCached.restore(stream);
}

/*!
 * Get the key associated to the iterator.
 *
 * \param itr is the iterator
 * \result The key associated to the iterator.
 */
template<typename key_t, typename value_t>
key_t LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::getKey(const const_iterator &itr) const
{
    return itr.getId();
}

/*!
 * Get a reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::reference LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::getValue(const iterator &itr) const
{
    return *itr;
}

/*!
 * Get a constant reference to the value pointed by the iterator.
 *
 * \param itr is the iterator
 * \result A constant reference to the value pointed by the iterator.
 */
template<typename key_t, typename value_t>
typename LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::const_reference LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::getValue(const const_iterator &itr) const
{
    return *itr;
}

/*!
 * Create the cache.
 *
 * \result The newly created cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCacheFactory<key_t, std::unordered_map<key_t ,value_t>>::create() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, std::unordered_map<key_t ,value_t>>());
}

/*!
 * Internal function to create the cache.
 *
 * \result The newly created cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCacheFactory<key_t, std::vector<value_t>>::create() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, std::vector<value_t>>());
}

/*!
 * Create a cache of the specified type.
 *
 * \result The newly created cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCacheFactory<key_t, PiercedVector<value_t, key_t>>::create() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>());
}

/*!
 * Constructor.
 *
 * \param kernel is the kernel associated with the entries
 * \param syncMode is the synchronization mode that will be used for the storage
 */
template<typename key_t, typename value_t>
LevelSetContainerCacheFactory<key_t, PiercedStorage<value_t, key_t>>::LevelSetContainerCacheFactory(PiercedKernel<key_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
    : Base(),
      m_kernel(kernel), m_syncMode(syncMode)
{
}

/*!
 * Create a cache of the specified type.
 *
 * \result The newly created cache.
 */
template<typename key_t, typename value_t>
std::unique_ptr<LevelSetCache<key_t>> LevelSetContainerCacheFactory<key_t, PiercedStorage<value_t, key_t>>::create() const
{
    return std::unique_ptr<LevelSetCache<key_t>>(new LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>(m_kernel, m_syncMode));
}

/*!
 * Is the index associated with an invalid cache.
 */
template<typename key_t>
const std::size_t LevelSetCacheCollection<key_t>::NULL_ID = std::numeric_limits<std::size_t>::max();

/*!
 * Copy constructor.
 */
template<typename key_t>
LevelSetCacheCollection<key_t>::LevelSetCacheCollection(const LevelSetCacheCollection &other)
    : m_factories(other.m_factories),
      m_caches(other.m_caches.size())
{
    for (std::size_t i = 0; i < other.m_caches.size(); ++i) {
        m_caches[i] = other.m_caches[i]->clone();
    }
}

/*!
 * Clone the cache collection.
 *
 * \result The cloned cache collection.
 */
template<typename key_t>
std::unique_ptr<LevelSetCacheCollection<key_t>> LevelSetCacheCollection<key_t>::clone() const
{
    return std::unique_ptr<LevelSetCacheCollection<key_t>>(new LevelSetCacheCollection<key_t>(*this));
}

/*!
 * Get an iterator pointing to the first entry.
 *
 * \result An iterator pointing to the first entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::iterator LevelSetCacheCollection<key_t>::begin()
{
    return m_caches.begin();
}

/*!
 * Get an iterator referring to the past-the-end entry.
 *
 * \result An iterator referring to the past-the-end entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::iterator LevelSetCacheCollection<key_t>::end()
{
    return m_caches.end();
}

/*!
 * Get a constant iterator pointing to the first entry.
 *
 * \result A constant iterator pointing to the first entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::const_iterator LevelSetCacheCollection<key_t>::begin() const
{
    return m_caches.cbegin();
}

/*!
 * Get a constant iterator referring to the past-the-end entry.
 *
 * \result a constant iterator referring to the past-the-end entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::const_iterator LevelSetCacheCollection<key_t>::end() const
{
    return m_caches.cend();
}

/*!
 * Get a constant iterator pointing to the first entry.
 *
 * \result A constant iterator pointing to the first entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::const_iterator LevelSetCacheCollection<key_t>::cbegin() const
{
    return m_caches.cbegin();
}

/*!
 * Get a constant iterator referring to the past-the-end entry.
 *
 * \result a constant iterator referring to the past-the-end entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::const_iterator LevelSetCacheCollection<key_t>::cend() const
{
    return m_caches.cend();
}

/*!
 * Insert a new cache.
 *
 * \param args are the arguments that will be used to create the cache factory
 * \result The index associated with the newly created cache.
 */
template<typename id_t>
template<typename container_t, typename... Args>
std::size_t LevelSetCacheCollection<id_t>::insert(Args&&... args)
{
    typedef typename LevelSetContainerCache<id_t, container_t>::value_type value_type;

    // Create the factory
    auto factory = std::shared_ptr<LevelSetValueCacheFactory<id_t, value_type>>(new LevelSetContainerCacheFactory<id_t, container_t>(std::forward<Args>(args)...));

    // Search for an unused index
    //
    // If an index is not associated with a factory, it means that its corresponding cache has
    // been erased and therefore the index can be re-used.
    std::size_t nCaches = m_caches.size();
    for (std::size_t index = 0; index < nCaches; ++index) {
        if (!m_factories[index]) {
            m_factories[index] = factory;

            return index;
        }
    }

    // No indexes can be re-used, a new entry will be added.
    m_factories.emplace_back(factory);
    m_caches.emplace_back();

    return nCaches;
}

/*!
 * Erase the specified cache.
 *
 * \param index is the entry that will be erased.
 */
template<typename key_t>
void LevelSetCacheCollection<key_t>::erase(std::size_t index)
{
    // Early return if the specified index is not valid
    if (index == NULL_ID) {
        return;
    } else if (index >= m_caches.size()) {
        return;
    }

    // To avoid invalidating the index of other caches, we cannot remove the cache and its
    // factory from the collection. A cache will be deleted clearing destroying its storage
    // and its factory. The index associated to the cache will be re-used when a new cache
    // will be inserted.
    m_factories[index].reset();
    m_caches[index].reset();
}

/*!
 * Delete all the registered caches.
 */
template<typename key_t>
void LevelSetCacheCollection<key_t>::clear()
{
    m_factories.clear();
    m_caches.clear();

    m_factories.shrink_to_fit();
    m_caches.shrink_to_fit();
}

/*!
 * Get a pointer to the cache with the specified index.
 *
 * If the requested index is not associated with a cache, a null pointer is returned.
 *
 * \param index is the index of the cache
 * \return A pointer to the cache with the specified index.
 */
template<typename key_t>
LevelSetCache<key_t> * LevelSetCacheCollection<key_t>::at(std::size_t index)
{
    return const_cast<LevelSetCache<key_t> *>(const_cast<const LevelSetCacheCollection<key_t> &>(*this).at(index));
}

/*!
 * Get a constant pointer to the cache with the specified index.
 *
 * If the requested index is not associated with a cache, a null pointer is returned.
 *
 * \param index is the index of the cache
 * \return A constant pointer to the cache with the specified index.
 */
template<typename key_t>
const LevelSetCache<key_t> * LevelSetCacheCollection<key_t>::at(std::size_t index) const
{
    if (index >= m_caches.size()) {
        return nullptr;
    }

    const LevelSetCache<key_t> *cache = m_caches[index].get();
    if (!cache) {
        const std::shared_ptr<LevelSetCacheFactory<key_t>> &cacheFactory = m_factories[index];
        if (!cacheFactory) {
            return nullptr;
        }

        m_caches[index] = cacheFactory->create();
        cache = m_caches[index].get();
    }

    return cache;
}

/*!
 * Get a pointer to the cache with the specified index.
 *
 * If the requested index is not associated with a cache, a null pointer is returned. It's up
 * to the caller of this function to make sure the request value type is the same used when
 * registering the cache.
 *
 * \param index is the index of the cache
 * \return A pointer to the cache with the specified index.
 */
template<typename key_t>
template<typename value_t>
LevelSetValueCache<key_t, value_t> * LevelSetCacheCollection<key_t>::at(std::size_t index)
{
    return const_cast<LevelSetValueCache<key_t, value_t> *>(const_cast<const LevelSetCacheCollection<key_t> &>(*this).at<value_t>(index));
}

/*!
 * Get a constant pointer to the cache with the specified index.
 *
 * If the requested index is not associated with a cache, a null pointer is returned. It's up
 * to the caller of this function to make sure the request value type is the same used when
 * registering the cache.
 *
 * \param index is the index of the cache
 * \return A constant pointer to the cache with the specified index.
 */
template<typename key_t>
template<typename value_t>
const LevelSetValueCache<key_t, value_t> * LevelSetCacheCollection<key_t>::at(std::size_t index) const
{
    const LevelSetCache<key_t> *cache = at(index);
    if (!cache) {
        return nullptr;
    }

    assert(dynamic_cast<const LevelSetValueCache<key_t BITPIT_COMMA value_t> *>(cache));

    return static_cast<const LevelSetValueCache<key_t, value_t> *>(m_caches[index].get());
}

}

#endif 
