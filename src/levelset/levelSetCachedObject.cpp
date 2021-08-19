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

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# include "bitpit_communications.hpp"
# endif

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesianKernel.hpp"
# include "levelSetOctreeKernel.hpp"
# include "levelSetObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetSignPropagator.hpp"
# include "levelSet.hpp"

namespace bitpit {


/*!
	@ingroup levelset
	@interface LevelSetNarrowBandCache
	@brief Base class for defining narrow band caches.
*/

/*!
 * Constructor.
 */
LevelSetNarrowBandCache::LevelSetNarrowBandCache() {

    m_values    = addStorage<double>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_JOURNALED);
    m_gradients = addStorage<std::array<double, 3>>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_JOURNALED);

}

/*!
 * Get the levelset value of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset value of the specified entry.
 */
double LevelSetNarrowBandCache::getValue(const KernelIterator &itr) const {

    std::size_t rawId = itr.getRawIndex();

    return m_values->rawAt(rawId);

}

/*!
 * Get the levelset gradient of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The levelset gradient of the specified entry.
 */
const std::array<double, 3> & LevelSetNarrowBandCache::getGradient(const KernelIterator &itr) const {

    std::size_t rawId = itr.getRawIndex();

    return m_gradients->rawAt(rawId);

}

/*!
 * Set the specified cache entry.
 *
 * \param itr is an iterator pointing to the cache entry
 * \param value is the levelset value
 * \param gradient is the levelset gradient
 */
void LevelSetNarrowBandCache::set(const KernelIterator &itr, double value, const std::array<double, 3> &gradient) {

    std::size_t rawId = itr.getRawIndex();

    m_values->rawAt(rawId)    = value;
    m_gradients->rawAt(rawId) = gradient;

}

/*!
 * Exchanges the content of the cache with the content the specified other
 * cache.
 *
 * \param other is another cache whose content is swapped with that of this
 * cache
 */
void LevelSetNarrowBandCache::swap(LevelSetNarrowBandCache &other) noexcept
{
    LevelSetInternalPiercedStorageManager::swap(other);

    std::swap(other.m_values, m_values);
    std::swap(other.m_gradients, m_gradients);
}

/*!
 * \ingroup levelset
 * \class LevelSetCachedObjectInterface
 * \brief The class LevelSetCachedObjectInterface allows to define objects
 * that cache narrow band information.
 */

/*!
 * Initialize the cache.
 */
LevelSetNarrowBandCache * LevelSetCachedObjectInterface::initializeNarrowBandCache()
{
    m_narrowBandCache = createNarrowBandCache();

    return getNarrowBandCache();
}

/*!
 * Get a pointer to the cache.
 *
 * \result A pointer to the cache.
 */
LevelSetNarrowBandCache * LevelSetCachedObjectInterface::getNarrowBandCache()
{
    return m_narrowBandCache.get();
}

/*!
 * Get a constant pointer to the cache.
 *
 * \result A constant pointer to the cache.
 */
const LevelSetNarrowBandCache * LevelSetCachedObjectInterface::getNarrowBandCache() const
{
    return m_narrowBandCache.get();
}

/*!
 * Clear the storage.
 */
void LevelSetCachedObjectInterface::clearNarrowBandCache()
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
void LevelSetCachedObjectInterface::dumpNarrowBandCache(std::ostream &stream)
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
void LevelSetCachedObjectInterface::restoreNarrowBandCache(std::istream &stream)
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
bool LevelSetCachedObjectInterface::isInNarrowBand(long id)const
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
void LevelSetCachedObjectInterface::swap(LevelSetCachedObjectInterface &other) noexcept
{
    m_narrowBandCache.swap(other.m_narrowBandCache);
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
LevelSetCachedObject::LevelSetCachedObject(int id)
    : LevelSetObject(id)
{
}

/*!
 * Sets the kernel for the object
 * @param[in] kernel is the LevelSetKernel
 */
void LevelSetCachedObject::setKernel(LevelSetKernel *kernel) {

    LevelSetObject::setKernel(kernel);

    initializeNarrowBandCache();

}

/*!
 * Create the storage for the narrow band data.
 */
std::shared_ptr<LevelSetNarrowBandCache> LevelSetCachedObject::createNarrowBandCache() {

    return std::shared_ptr<LevelSetNarrowBandCache>(new LevelSetNarrowBandCache());

}

/*!
 * Create the storage for sign propagation.
 */
std::shared_ptr<LevelSetSignStorage> LevelSetCachedObject::createSignStorage() {

    VolumeKernel *mesh = m_kernelPtr->getMesh() ;
    assert(mesh) ;

    return std::shared_ptr<LevelSetSignStorage>(new LevelSetSignStorage(&(mesh->getCells())));

}

/*!
 * Get LevelSetInfo of cell
 * @param[in] i cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetCachedObject::getLevelSetInfo( long id)const{

    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
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
double LevelSetCachedObject::getLS( long id)const {

    return getValue(id);

}

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetCachedObject::getValue( long id)const {

    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
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
short LevelSetCachedObject::getSign( long id ) const {

    // Check if the sign can be evaluated from narrowband value
    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
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
std::array<double,3> LevelSetCachedObject::getGradient(long id) const {

    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id) ;
    if( narrowBandCacheItr != narrowBandCache->end() ){
        return narrowBandCache->getGradient(narrowBandCacheItr);
    }

    return levelSetDefaults::GRADIENT;

}

/*! 
 * Deletes non-existing items after grid adaption.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSetCachedObject::_clearAfterMeshAdaption( const std::vector<adaption::Info> &adaptionData ){

    // Clear stale narrow band entries
    LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
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
void LevelSetCachedObject::_clear( ){

    // Clear narrow band entries
    clearNarrowBandCache();

    // Clear sign propgation storage
    clearSignStorage();
}

/*!
 * Writes LevelSetCachedObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_dump( std::ostream &stream ){

    // Narrow band storage
    dumpNarrowBandCache( stream );

    // Stored sign
    dumpSignStorage( stream );
}

/*!
 * Reads LevelSetCachedObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_restore( std::istream &stream ){

    // Narrow band storage
    restoreNarrowBandCache( stream );

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
void LevelSetCachedObject::_writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    getNarrowBandCache()->write(sendList, dataBuffer);
}

/*!
 * Processing of communication buffer into data structure
 * Sign data is not read from the buffer, because sign storage is kept in
 * sync with the mesh, hence, when the buffer is written, entries associated
 * with the cells to send as already been deleted.
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetCachedObject::_readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    getNarrowBandCache()->read(recvList, dataBuffer);
}

#endif

}
