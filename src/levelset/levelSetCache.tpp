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
 * Keys not associated with an entry will be skipped without raising any exception.
 *
 * \param keys are the keys whose entries will be erased
 * \result The number of entries that were actually erased.
 */
template<typename key_t>
template<typename Keys>
std::size_t LevelSetCache<key_t>::erase(const Keys &keys)
{
    std::size_t nDeletedEntries = 0;
    for (const key_t &key : keys) {
        nDeletedEntries += erase(key);
    }

    return nDeletedEntries;
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
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::Entry LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::findEntry(const key_t &key) const
{
    const_iterator itr = find(key);
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
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::Entry LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::insertEntry(const key_t &key, const value_t &value)
{
    const_iterator itr = insert(key, value);

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
typename LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::Entry LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::insertEntry(const key_t &key, value_t &&value)
{
    const_iterator itr = insert(key, std::move(value));

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
void LevelSetContainerBaseCache<key_t, container_t, value_t, reference_t, const_reference_t>::writeBuffer(const std::vector<key_t> &keys, SendBuffer &buffer) const
{
    for (const key_t &key : keys) {
        const_iterator itr = find(key);
        if (itr != cend()) {
            buffer << static_cast<unsigned char>(1);
            buffer << getValue(itr);
        } else {
            buffer << static_cast<unsigned char>(0);
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
    value_t value;
    for (const key_t &key : keys) {
        unsigned char isCached;
        buffer >> isCached;
        if (isCached == 1) {
            buffer >> value;

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
 * \result The number of entries that were actually erased.
 */
template<typename key_t, typename value_t>
std::size_t LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::erase(const key_t &key)
{
    return Base::m_container.erase(key);
}

/*!
 * Check if cache content is volatile.
 *
 * When a cache is volatile, its contents might be changed by means undetectable by the cache
 * itself. The cache is still in charge of keeping track of the entries it contains (e.g., it
 * is in charge of erasing the entries that are no longer in the cache), but the values of the
 * entries (and even the existence of the values) might be changed by someone outside the cache.
 * After the values of a volatile cache has been updated outside the cache itself, the cache
 * will be in an inconsistent state until the cache is explicitly updated to reflect the changes
 * already performed on its contents. Accessing the entries whose values have been changed from
 * outside before the cache has been explicitly updated, results in undefined behavior.
 *
 * An example of a volatile cache is a cache in which the entries are stored in a PiercedStorage
 * synchronized with a mesh. When the mesh is adapted, the container that stores the cache entries
 * will be automatically synchronized, but the list of cache entries will remain the same until
 * the cache is explicitly updated. At the time the cache is updated to reflect the changes in the
 * mesh, since the PiercedStorage has already been updated, the entries associated with items that
 * no longer exists in the mesh will no longer have their values accessible and the PiercedStorage
 * will already contain uninitialized values for the newly created items even if they are not yet
 * in the cache. As soon as the mesh is updated, all entries associated to items deleted from the
 * mesh should be considered invalid and cannot be accessed (and the cache will not be aware that
 * these entries have been deleted until it will be explicitly updated). A practical limitation
 * of volatile caches is that it is not possible to communicate the entries that have been send
 * to other processes. When the cache is updated, the values associated with those entries have
 * already been deleted (at the time of the synchronization of the storage).
 *
 * \result Returns true if the cache is volatile, false otherwise.
 */
template<typename key_t, typename value_t>
bool LevelSetContainerCache<key_t, std::unordered_map<key_t, value_t>>::isVolatile() const
{
    return false;
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
 * \result The number of entries that were actually erased.
 */
template<typename key_t, typename value_t>
std::size_t LevelSetContainerCache<key_t, std::vector<value_t>>::erase(const key_t &key)
{
    if (static_cast<std::size_t>(key) >= m_isCached.size()) {
        return 0;
    }

    m_isCached[key] = false;

    return 1;
}

/*!
 * Check if cache content is volatile.
 *
 * When a cache is volatile, its contents might be changed by means undetectable by the cache
 * itself. The cache is still in charge of keeping track of the entries it contains (e.g., it
 * is in charge of erasing the entries that are no longer in the cache), but the values of the
 * entries (and even the existence of the values) might be changed by someone outside the cache.
 * After the values of a volatile cache has been updated outside the cache itself, the cache
 * will be in an inconsistent state until the cache is explicitly updated to reflect the changes
 * already performed on its contents. Accessing the entries whose values have been changed from
 * outside before the cache has been explicitly updated, results in undefined behavior.
 *
 * An example of a volatile cache is a cache in which the entries are stored in a PiercedStorage
 * synchronized with a mesh. When the mesh is adapted, the container that stores the cache entries
 * will be automatically synchronized, but the list of cache entries will remain the same until
 * the cache is explicitly updated. At the time the cache is updated to reflect the changes in the
 * mesh, since the PiercedStorage has already been updated, the entries associated with items that
 * no longer exists in the mesh will no longer have their values accessible and the PiercedStorage
 * will already contain uninitialized values for the newly created items even if they are not yet
 * in the cache. As soon as the mesh is updated, all entries associated to items deleted from the
 * mesh should be considered invalid and cannot be accessed (and the cache will not be aware that
 * these entries have been deleted until it will be explicitly updated). A practical limitation
 * of volatile caches is that it is not possible to communicate the entries that have been send
 * to other processes. When the cache is updated, the values associated with those entries have
 * already been deleted (at the time of the synchronization of the storage).
 *
 * \result Returns true if the cache is volatile, false otherwise.
 */
template<typename key_t, typename value_t>
bool LevelSetContainerCache<key_t, std::vector<value_t>>::isVolatile() const
{
    return false;
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
 * \result The number of entries that were actually erased.
 */
template<typename key_t, typename value_t>
std::size_t LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::erase(const key_t &key)
{
    if (!Base::m_container.contains(key)) {
        return 0;
    }

    Base::m_container.erase(key);

    return 1;
}

/*!
 * Check if cache content is volatile.
 *
 * When a cache is volatile, its contents might be changed by means undetectable by the cache
 * itself. The cache is still in charge of keeping track of the entries it contains (e.g., it
 * is in charge of erasing the entries that are no longer in the cache), but the values of the
 * entries (and even the existence of the values) might be changed by someone outside the cache.
 * After the values of a volatile cache has been updated outside the cache itself, the cache
 * will be in an inconsistent state until the cache is explicitly updated to reflect the changes
 * already performed on its contents. Accessing the entries whose values have been changed from
 * outside before the cache has been explicitly updated, results in undefined behavior.
 *
 * An example of a volatile cache is a cache in which the entries are stored in a PiercedStorage
 * synchronized with a mesh. When the mesh is adapted, the container that stores the cache entries
 * will be automatically synchronized, but the list of cache entries will remain the same until
 * the cache is explicitly updated. At the time the cache is updated to reflect the changes in the
 * mesh, since the PiercedStorage has already been updated, the entries associated with items that
 * no longer exists in the mesh will no longer have their values accessible and the PiercedStorage
 * will already contain uninitialized values for the newly created items even if they are not yet
 * in the cache. As soon as the mesh is updated, all entries associated to items deleted from the
 * mesh should be considered invalid and cannot be accessed (and the cache will not be aware that
 * these entries have been deleted until it will be explicitly updated). A practical limitation
 * of volatile caches is that it is not possible to communicate the entries that have been send
 * to other processes. When the cache is updated, the values associated with those entries have
 * already been deleted (at the time of the synchronization of the storage).
 *
 * \result Returns true if the cache is volatile, false otherwise.
 */
template<typename key_t, typename value_t>
bool LevelSetContainerCache<key_t, PiercedVector<value_t, key_t>>::isVolatile() const
{
    return false;
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
LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::LevelSetContainerCache(const PiercedKernel<key_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
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
 * \result The number of entries that were actually erased.
 */
template<typename key_t, typename value_t>
std::size_t LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::erase(const key_t &key)
{
    // Even if we are calling this function after the storages have been synchronized, we may
    // still need to mark the element associated with the key as non-cached. After synchronizing
    // the storages, the key may be associated with a different element. However, since we are
    // calling this function, no cache have been initialized for this new element. Hence, we need
    // to mark the new element as non-cached. If the key is not associated with a new element,
    // nothing should be done, because the storages don't contain anymore an entry for the key.
    auto itr = Base::m_container.find(key);
    if (itr == Base::m_container.end()) {
        return 0;
    }

    std::size_t rawId = itr.getRawIndex();

    m_isCached.rawAt(rawId) = false;

    return 1;
}

/*!
 * Check if cache content is volatile.
 *
 * When a cache is volatile, its contents might be changed by means undetectable by the cache
 * itself. The cache is still in charge of keeping track of the entries it contains (e.g., it
 * is in charge of erasing the entries that are no longer in the cache), but the values of the
 * entries (and even the existence of the values) might be changed by someone outside the cache.
 * After the values of a volatile cache has been updated outside the cache itself, the cache
 * will be in an inconsistent state until the cache is explicitly updated to reflect the changes
 * already performed on its contents. Accessing the entries whose values have been changed from
 * outside before the cache has been explicitly updated, results in undefined behavior.
 *
 * An example of a volatile cache is a cache in which the entries are stored in a PiercedStorage
 * synchronized with a mesh. When the mesh is adapted, the container that stores the cache entries
 * will be automatically synchronized, but the list of cache entries will remain the same until
 * the cache is explicitly updated. At the time the cache is updated to reflect the changes in the
 * mesh, since the PiercedStorage has already been updated, the entries associated with items that
 * no longer exists in the mesh will no longer have their values accessible and the PiercedStorage
 * will already contain uninitialized values for the newly created items even if they are not yet
 * in the cache. As soon as the mesh is updated, all entries associated to items deleted from the
 * mesh should be considered invalid and cannot be accessed (and the cache will not be aware that
 * these entries have been deleted until it will be explicitly updated). A practical limitation
 * of volatile caches is that it is not possible to communicate the entries that have been send
 * to other processes. When the cache is updated, the values associated with those entries have
 * already been deleted (at the time of the synchronization of the storage).
 *
 * \result Returns true if the cache is volatile, false otherwise.
 */
template<typename key_t, typename value_t>
bool LevelSetContainerCache<key_t, PiercedStorage<value_t, key_t>>::isVolatile() const
{
    return (Base::m_container.getSyncMode() != PiercedSyncMaster::SyncMode::SYNC_MODE_DISABLED);
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
LevelSetContainerCacheFactory<key_t, PiercedStorage<value_t, key_t>>::LevelSetContainerCacheFactory(const PiercedKernel<key_t> *kernel, PiercedSyncMaster::SyncMode syncMode)
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
 * \interface LevelSetCacheCollection::Item
 * \ingroup levelset
 *
 * \brief The class LevelSetCacheCollection::Item defines the items stored in
 * a cache collection. Each item stores the cache and the factory that will be
 * used to create the cache.
 */

/*!
 * Constructor.
 *
 * \param factory is the factory that will be used to create the cache
 */
template<typename key_t>
LevelSetCacheCollection<key_t>::Item::Item(const std::shared_ptr<LevelSetCacheFactory<key_t>> &factory)
    : m_factory(factory),
      m_cache(nullptr)
{
}

/*!
 * Copy constructor.
 *
 * \param other is another item whose content is copied in this item
 */
template<typename key_t>
LevelSetCacheCollection<key_t>::Item::Item(const LevelSetCacheCollection<key_t>::Item &other)
    : m_factory(other.m_factory),
      m_cache(other.m_cache->clone())
{
}

/**
 * Copy assignment operator.
 *
 * \param other is another collection whose content is copied in this collection
 * \return A reference to the item.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::Item & LevelSetCacheCollection<key_t>::Item::operator=(const LevelSetCacheCollection<key_t>::Item &other)
{
    m_factory = other.m_factory;
    m_cache   = other.m_cache->clone();

    return *this;
}

/*!
 * Check if a factory has been set for the item.
 *
 * \result Return true if a factory has been set for the item, false otherwise.
 */
template<typename key_t>
bool LevelSetCacheCollection<key_t>::Item::hasFactory() const
{
    return !!m_factory;
}

/*!
 * Check if a cache has been created for the item.
 *
 * \result Return true if a has has been created for the item, false otherwise.
 */
template<typename key_t>
bool LevelSetCacheCollection<key_t>::Item::hasCache() const
{
    return !!m_cache;
}

/*!
 * Create the cache associated with the item.
 *
 * If the cache has been previously created, the existing cache will be
 * returned. If no factory has been set for the item, a null pointer is
 * returned.
 *
 * \result Return a pointer to the newly created cache.
 */
template<typename key_t>
LevelSetCache<key_t> * LevelSetCacheCollection<key_t>::Item::createCache()
{
    return _createCache();
}

/*!
 * Internal function to create the cache associated with the item.
 *
 * If the cache has been previously created, the existing cache will be
 * returned. If no factory has been set for the item, a null pointer is
 * returned.
 *
 * \result Return a pointer to the newly created cache.
 */
template<typename key_t>
LevelSetCache<key_t> * LevelSetCacheCollection<key_t>::Item::_createCache() const
{
    // Early return if the cache has already been created
    if (hasCache()) {
        return m_cache.get();
    }

    // Create the cache
    if (m_factory) {
        m_cache = m_factory->create();
        return m_cache.get();
    }

    // No cache can be created for the item
    return nullptr;
}

/*!
 * Destroy the cache associated with the item.
 */
template<typename key_t>
void LevelSetCacheCollection<key_t>::Item::destroyCache()
{
    return m_cache.reset();
}

/*!
 * Get a pointer to the cache owned by the item.
 *
 * \param allowCreation controls the behaviour when the cache has not yet
 * been created, if the argument is set to true, the item will create a new
 * cache from scratch, otherwise a null pointer is returned
 * \return A pointer to the cache owned by the item.
 */
template<typename key_t>
LevelSetCache<key_t> * LevelSetCacheCollection<key_t>::Item::getCache(bool allowCreation)
{
    return const_cast<LevelSetCache<key_t> *>(const_cast<const LevelSetCacheCollection<key_t>::Item &>(*this).getCache(allowCreation));
}

/*!
 * Get a constant pointer to the cache owned by the item.
 *
 * \param allowCreation controls the behaviour when the cache has not yet
 * been created, if the argument is set to true, the item will create a new
 * cache from scratch, otherwise a null pointer is returned
 * \return A constant pointer to the cache owned by the item.
 */
template<typename key_t>
const LevelSetCache<key_t> * LevelSetCacheCollection<key_t>::Item::getCache(bool allowCreation) const
{
    // Early return if the cache has already been created
    if (hasCache()) {
        return m_cache.get();
    }

    // Create the cache
    if (allowCreation) {
        return _createCache();
    }

    // No cache can be created for the item
    return nullptr;
}

/*!
 * Get a pointer to the cache owned by the item.
 *
 * \param allowCreation controls the behaviour when the cache has not yet
 * been created, if the argument is set to true, the item will create a new
 * cache from scratch, otherwise a null pointer is returned
 * \return A pointer to the cache owned by the item.
 */
template<typename key_t>
template<typename value_t>
LevelSetValueCache<key_t, value_t> * LevelSetCacheCollection<key_t>::Item::getCache(bool allowCreation)
{
    return const_cast<LevelSetValueCache<key_t, value_t> *>(const_cast<const LevelSetCacheCollection<key_t>::Item &>(*this).getCache<value_t>(allowCreation));
}

/*!
 * Get a constant pointer to the cache owned by the item.
 *
 * \param allowCreation controls the behaviour when the cache has not yet
 * been created, if the argument is set to true, the item will create a new
 * cache from scratch, otherwise a null pointer is returned
 * \return A constant pointer to the cache owned by the item.
 */
template<typename key_t>
template<typename value_t>
const LevelSetValueCache<key_t, value_t> * LevelSetCacheCollection<key_t>::Item::getCache(bool allowCreation) const
{
    assert(!hasCache() || dynamic_cast<const LevelSetValueCache<key_t BITPIT_COMMA value_t> *>(getCache(allowCreation)));

    return static_cast<const LevelSetValueCache<key_t, value_t> *>(getCache(allowCreation));
}

/*!
 * \interface LevelSetCacheCollection
 * \ingroup levelset
 *
 * \brief The class LevelSetCacheCollection allows to store a collection of caches.
 */

/*!
 * Is the index associated with an invalid cache.
 */
template<typename key_t>
const std::size_t LevelSetCacheCollection<key_t>::NULL_CACHE_ID = std::numeric_limits<std::size_t>::max();

/*!
 * Copy constructor.
 *
 * \param other is another collection whose content is copied in this collection
 */
template<typename key_t>
LevelSetCacheCollection<key_t>::LevelSetCacheCollection(const LevelSetCacheCollection &other)
    : m_caches(other.m_caches.size())
{
    for (std::size_t i = 0; i < other.m_caches.size(); ++i) {
        m_caches[i] = other.m_caches[i];
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
 * \result A constant iterator referring to the past-the-end entry.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::const_iterator LevelSetCacheCollection<key_t>::cend() const
{
    return m_caches.cend();
}

/*!
 * Get a size of the collection.
 *
 * \result The size of the collection.
 */
template<typename key_t>
std::size_t LevelSetCacheCollection<key_t>::size() const
{
    return m_caches.size();
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
    // Create the factory
    typedef typename LevelSetContainerCache<id_t, container_t>::value_type value_type;

    auto factory = std::shared_ptr<LevelSetValueCacheFactory<id_t, value_type>>(new LevelSetContainerCacheFactory<id_t, container_t>(std::forward<Args>(args)...));

    // Search for an unused index
    //
    // If an index is not associated with a factory, it means that its corresponding cache has
    // been erased and therefore the index can be re-used.
    std::size_t nCaches = m_caches.size();
    for (std::size_t index = 0; index < nCaches; ++index) {
        if (!m_caches[index].hasFactory()) {
            m_caches[index] = Item(factory);

            return index;
        }
    }

    // No indexes can be re-used, a new entry will be added.
    m_caches.emplace_back(factory);

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
    if (index == NULL_CACHE_ID) {
        return;
    } else if (index >= m_caches.size()) {
        return;
    }

    // To avoid invalidating the index of other caches, we cannot remove the cache and its
    // factory from the collection. A cache will be deleted clearing destroying its storage
    // and its factory. The index associated to the cache will be re-used when a new cache
    // will be inserted.
    m_caches[index] = Item();
}

/*!
 * Delete all the registered caches.
 */
template<typename key_t>
void LevelSetCacheCollection<key_t>::clear()
{
    m_caches.clear();
    m_caches.shrink_to_fit();
}

/*!
 * Get a reference to the cache item with the specified index.
 *
 * A similar member function, Item::at, has the same behavior as this operator
 * function, except that Item::at is bound-checked and signals if the requested
 * position is out of range by throwing an out_of_range exception.
 *
 * \param index is the index of the cache
 * \return A reference to the cache item with the specified index.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::Item & LevelSetCacheCollection<key_t>::operator[](std::size_t index)
{
    return const_cast<LevelSetCacheCollection<key_t>::Item &>(const_cast<const LevelSetCacheCollection<key_t> &>(*this).at(index));
}

/*!
 * Get a constant reference to the cache item with the specified index.
 *
 * A similar member function, Item::at, has the same behavior as this operator
 * function, except that Item::at is bound-checked and signals if the requested
 * position is out of range by throwing an out_of_range exception.
 *
 * \param index is the index of the cache
 * \return A constant reference to the cache item with the specified index.
 */
template<typename key_t>
const typename LevelSetCacheCollection<key_t>::Item & LevelSetCacheCollection<key_t>::operator[](std::size_t index) const
{
    return m_caches[index];
}

/*!
 * Get a reference to the cache item with the specified index.
 *
 * If this specified index is greater than, or equal to, the collection size,
 * an exception of type out_of_range is thrown.
 *
 * \param index is the index of the cache
 * \return A reference to the cache item with the specified index.
 */
template<typename key_t>
typename LevelSetCacheCollection<key_t>::Item & LevelSetCacheCollection<key_t>::at(std::size_t index)
{
    return const_cast<LevelSetCacheCollection<key_t>::Item &>(const_cast<const LevelSetCacheCollection<key_t> &>(*this).at(index));
}

/*!
 * Get a constant reference to the cache item with the specified index.
 *
 * If this specified index is greater than, or equal to, the collection size,
 * an exception of type out_of_range is thrown.
 *
 * \param index is the index of the cache
 * \return A constant reference to the cache item with the specified index.
 */
template<typename key_t>
const typename LevelSetCacheCollection<key_t>::Item & LevelSetCacheCollection<key_t>::at(std::size_t index) const
{
    return m_caches.at(index);
}

/*!
 * Insert a new cache.
 *
 * \param args are the arguments that will be used to create the cache factory
 * \result The index associated with the newly created cache.
 */
template<typename container_t, typename... Args, typename std::enable_if<std::is_same<bitpit::PiercedStorage<typename container_t::value_type>, container_t>::value>::type *>
std::size_t ElementCacheCollection::insert(Args&&... args)
{
    PiercedSyncMaster::SyncMode cacheSyncMode = PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED;

    return LevelSetCacheCollection<long>::insert<container_t>(m_kernel, cacheSyncMode, std::forward<Args>(args)...);
}

/*!
 * Insert a new cache.
 *
 * \param args are the arguments that will be used to create the cache factory
 * \result The index associated with the newly created cache.
 */
template<typename container_t, typename... Args, typename std::enable_if<!std::is_same<bitpit::PiercedStorage<typename container_t::value_type>, container_t>::value>::type *>
std::size_t ElementCacheCollection::insert(Args&&... args)
{
    return Base::insert<container_t, Args...>(std::forward<Args>(args)...);
}

}

#endif 
