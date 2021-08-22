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

# ifndef __BITPIT_LEVELSET_CACHED_TPP__
# define __BITPIT_LEVELSET_CACHED_TPP__

namespace bitpit {

/*!
	@ingroup levelset
	@interface LevelSetNarrowBandCacheBase
	@brief Base class for defining narrow band caches.
*/

/*!
 * Constructor.
 */
template<typename storage_manager_t>
LevelSetNarrowBandCacheBase<storage_manager_t>::LevelSetNarrowBandCacheBase()
    : storage_manager_t()
{
}

/*!
 * Set the specified cache entry.
 *
 * \param itr is an iterator pointing to the cache entry
 * \param value is the levelset value
 * \param gradient is the levelset gradient
 */
template<typename storage_manager_t>
void LevelSetNarrowBandCacheBase<storage_manager_t>::set(const KernelIterator &itr, double value, const std::array<double, 3> &gradient) {

    double &cachedValue = getValue(itr);
    cachedValue = value;

    std::array<double, 3> &cachedGradient = getGradient(itr);
    cachedGradient = gradient;

}

/*!
 * Exchanges the content of the cache with the content the specified other
 * cache.
 *
 * \param other is another cache whose content is swapped with that of this
 * cache
 */
template<typename storage_manager_t>
void LevelSetNarrowBandCacheBase<storage_manager_t>::swap(LevelSetNarrowBandCacheBase &other) noexcept
{
    storage_manager_t::swap(other);

    std::swap(other.m_values, m_values);
    std::swap(other.m_gradients, m_gradients);
}

/*!
	@ingroup levelset
	@interface LevelSetNarrowBandCache
	@brief Narrow band cache.
*/

/*!
 * \ingroup levelset
 * \class LevelSetCachedObjectInterface
 * \brief The class LevelSetCachedObjectInterface allows to define objects
 * that cache narrow band information.
 */

/*!
 * Initialize the cache.
 */
template<typename narrow_band_cache_t>
narrow_band_cache_t * LevelSetCachedObjectInterface<narrow_band_cache_t>::initializeNarrowBandCache()
{
    m_narrowBandCache = LevelSetNarrowBandCacheFactory<narrow_band_cache_t>::create(this);

    return getNarrowBandCache();
}

/*!
 * Get a pointer to the cache.
 *
 * \result A pointer to the cache.
 */
template<typename narrow_band_cache_t>
narrow_band_cache_t * LevelSetCachedObjectInterface<narrow_band_cache_t>::getNarrowBandCache()
{
    return m_narrowBandCache.get();
}

/*!
 * Get a constant pointer to the cache.
 *
 * \result A constant pointer to the cache.
 */
template<typename narrow_band_cache_t>
const narrow_band_cache_t * LevelSetCachedObjectInterface<narrow_band_cache_t>::getNarrowBandCache() const
{
    return m_narrowBandCache.get();
}

/*!
 * Clear the storage.
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObjectInterface<narrow_band_cache_t>::clearNarrowBandCache()
{
    // Narrowband cache
    if (m_narrowBandCache) {
        m_narrowBandCache->clear();
    }
}

/*!
 * Dump the storage.
 *
 * \param stream is the output stream
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObjectInterface<narrow_band_cache_t>::dumpNarrowBandCache(std::ostream &stream)
{
    // Narrowband cache
    if (m_narrowBandCache) {
        m_narrowBandCache->dump( stream );
    }
}

/*!
 * Restore the storage.
 *
 * \param stream is the input stream
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObjectInterface<narrow_band_cache_t>::restoreNarrowBandCache(std::istream &stream)
{
    // Narrowband cache
    if (m_narrowBandCache) {
        m_narrowBandCache->restore( stream );
    }
}

/*!
 * Check if the specified cell lies within the narrow band and hence its
 * levelset is computed exactly.
 *
 * \param[in] id is the cell id
 * \resutl Return true if the cell is in the narrow band, falst otherwise.
 */
template<typename narrow_band_cache_t>
bool LevelSetCachedObjectInterface<narrow_band_cache_t>::isInNarrowBand(long id)const
{
    return m_narrowBandCache->contains(id);
}

/*!
 * Exchanges the content of the object with the content the specified other
 * object.
 *
 * \param other is another object whose content is swapped with that of this
 * object
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObjectInterface<narrow_band_cache_t>::swap(LevelSetCachedObjectInterface &other) noexcept
{
    m_narrowBandCache.swap(other.m_narrowBandCache);
}

/*!
 * \ingroup levelset
 * \class LevelSetNarrowBandCacheFactory
 * \brief The class LevelSetNarrowBandCacheFactory allows to create narrow band
 * cache objects.
 */

/*!
 * Create a narrow band cache.
 *
 * @param ojbect is the levelset object for which the ache will be created
 */
template<typename narrow_band_cache_t>
std::shared_ptr<narrow_band_cache_t> LevelSetNarrowBandCacheFactory<narrow_band_cache_t>::create(LevelSetCachedObjectInterface<narrow_band_cache_t> *object) {

    BITPIT_UNUSED(object);

    return std::shared_ptr<narrow_band_cache_t>(new narrow_band_cache_t());

}

/*!
	@ingroup levelset
	@interface LevelSetCachedObject
	@brief Interface class for all objects which need to store the discrete values of levelset function.
*/

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
template<typename narrow_band_cache_t>
LevelSetCachedObject<narrow_band_cache_t>::LevelSetCachedObject(int id)
    : LevelSetObject(id)
{
}

/*!
 * Sets the kernel for the object
 * @param[in] kernel is the LevelSetKernel
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::setKernel(LevelSetKernel *kernel) {

    LevelSetObject::setKernel(kernel);

    this->initializeNarrowBandCache();

}

/*!
 * Create the storage for sign propagation.
 */
template<typename narrow_band_cache_t>
std::shared_ptr<LevelSetSignStorage> LevelSetCachedObject<narrow_band_cache_t>::createSignStorage() {

    VolumeKernel *mesh = m_kernel->getMesh() ;
    assert(mesh) ;

    return std::shared_ptr<LevelSetSignStorage>(new LevelSetSignStorage(&(mesh->getCells())));

}

/*!
 * Get LevelSetInfo of cell
 * @param[in] i cell idex
 * @return LevelSetInfo of cell
*/
template<typename narrow_band_cache_t>
LevelSetInfo LevelSetCachedObject<narrow_band_cache_t>::getLevelSetInfo( long id)const{

    const narrow_band_cache_t *narrowBandCache = this->getNarrowBandCache();
    typename narrow_band_cache_t::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
    if( narrowBandCacheItr != narrowBandCache->end() ){
        double value = narrowBandCache->getValue(narrowBandCacheItr);
        const std::array<double, 3> &gradient = narrowBandCache->getGradient(narrowBandCacheItr);

        return LevelSetInfo(value, gradient);
    }

    return LevelSetInfo();

} 

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
template<typename narrow_band_cache_t>
double LevelSetCachedObject<narrow_band_cache_t>::getLS( long id)const {

    return getValue(id);

}

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
template<typename narrow_band_cache_t>
double LevelSetCachedObject<narrow_band_cache_t>::getValue( long id)const {

    const narrow_band_cache_t *narrowBandCache = this->getNarrowBandCache();
    typename narrow_band_cache_t::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
    if( narrowBandCacheItr != narrowBandCache->end() ){
        return narrowBandCache->getValue(narrowBandCacheItr);
    }

    return getSign(id) * levelSetDefaults::VALUE;

}

/*!
 * Get the sign of the levelset function
 * @param[in] id cell id
 * @return sign of levelset
 */
template<typename narrow_band_cache_t>
short LevelSetCachedObject<narrow_band_cache_t>::getSign( long id ) const {

    // Check if the sign can be evaluated from narrowband value
    const narrow_band_cache_t *narrowBandCache = this->getNarrowBandCache();
    typename narrow_band_cache_t::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
    if( narrowBandCacheItr != narrowBandCache->end() ){
        double value = narrowBandCache->getValue(narrowBandCacheItr);

        return evalValueSign(value);
    }

    // Check if the sign can be evaluated from the propagation
    if (!isSignStorageDirty()) {
        const LevelSetSignStorage *propagatedSignStorage = getSignStorage();
        LevelSetSignStorage::KernelIterator propagatedSignStorageItr = propagatedSignStorage->find(id);
        LevelSetSignStorage::Sign propagatedSign = propagatedSignStorage->at(propagatedSignStorageItr);
        if (propagatedSign != LevelSetSignStorage::SIGN_UNDEFINED) {
            return static_cast<short>(propagatedSign);
        }
    }

    // Unable to evaluate the sign
    //
    // The sign cannot be evaluated, let's return the defualt sign.
    return levelSetDefaults::SIGN;

}

/*!
 * Get the levelset gradient of cell
 * @param[in] id cell id
 * @return levelset gradient in cell 
 */
template<typename narrow_band_cache_t>
std::array<double,3> LevelSetCachedObject<narrow_band_cache_t>::getGradient(long id) const {

    const narrow_band_cache_t *narrowBandCache = this->getNarrowBandCache();
    typename narrow_band_cache_t::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
    if( narrowBandCacheItr != narrowBandCache->end() ){
        return narrowBandCache->getGradient(narrowBandCacheItr);
    }

    return levelSetDefaults::GRADIENT;

}

/*! 
 * Deletes non-existing items after grid adaption.
 * @param[in] adaptionData are the information about the adaption
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::_clearAfterMeshAdaption( const std::vector<adaption::Info> &adaptionData ){

    // Clear stale narrow band entries
    narrow_band_cache_t *narrowBandCache = this->getNarrowBandCache();
    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        for (long cellId : adaptionInfo.previous) {
            if (narrowBandCache->contains(cellId)) {
                narrowBandCache->erase(cellId, false);
            }
        }
    }

    narrowBandCache->syncStorages();

}

/*! 
 * Clears all levelset information
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::_clear( ){

    // Clear narrow band entries
    this->clearNarrowBandCache();

    // Clear sign propgation storage
    clearSignStorage();
}

/*!
 * Writes LevelSetCachedObject to stream in binary format
 * @param[in] stream output stream
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::_dump( std::ostream &stream ){

    // Narrow band storage
    this->dumpNarrowBandCache( stream );

    // Stored sign
    dumpSignStorage( stream );
}

/*!
 * Reads LevelSetCachedObject from stream in binary format
 * @param[in] stream output stream
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::_restore( std::istream &stream ){

    // Narrow band storage
    this->restoreNarrowBandCache( stream );

    // Stored sign
    restoreSignStorage( stream );
}

#if BITPIT_ENABLE_MPI

/*!
 * Flushing of data to communication buffers for partitioning
 * Sign data is not written into the buffer, because sign storage is kept in
 * sync with the mesh, hence, when this function is called, entries associated
 * with the cells to send as already been deleted.
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::_writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    this->getNarrowBandCache()->write(sendList, dataBuffer);
}

/*!
 * Processing of communication buffer into data structure
 * Sign data is not read from the buffer, because sign storage is kept in
 * sync with the mesh, hence, when the buffer is written, entries associated
 * with the cells to send as already been deleted.
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
template<typename narrow_band_cache_t>
void LevelSetCachedObject<narrow_band_cache_t>::_readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    this->getNarrowBandCache()->read(recvList, dataBuffer);
}

#endif

}

#endif

