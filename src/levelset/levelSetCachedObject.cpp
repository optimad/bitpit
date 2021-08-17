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
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"
# include "levelSetObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetSignPropagator.hpp"
# include "levelSet.hpp"

namespace bitpit {

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
    : LevelSetObject(id),
      m_narrowBandValues(1, &m_narrowBandKernel, PiercedSyncMaster::SYNC_MODE_JOURNALED),
      m_narrowBandGradients(1, &m_narrowBandKernel, PiercedSyncMaster::SYNC_MODE_JOURNALED)
{
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

    auto narrowBandItr = m_narrowBandKernel.find(id) ;
    if( narrowBandItr != m_narrowBandKernel.end() ){
        double value = m_narrowBandValues.rawAt(narrowBandItr.getRawIndex());
        const std::array<double, 3> &gradient = m_narrowBandGradients.rawAt(narrowBandItr.getRawIndex());

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

    auto narrowBandItr = m_narrowBandKernel.find(id) ;
    if( narrowBandItr != m_narrowBandKernel.end() ){
        return m_narrowBandValues.rawAt(narrowBandItr.getRawIndex());
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
    auto narrowBandItr = m_narrowBandKernel.find(id) ;
    if( narrowBandItr != m_narrowBandKernel.end() ){
        double value = m_narrowBandValues.rawAt(narrowBandItr.getRawIndex());

        return evalValueSign(value);
    }

    // Check if the sign can be evaluated from the propagation
    if (!isSignStorageDirty()) {
        LevelSetSignStorage::KernelIterator signStorageItr = getSignStorage()->find(id);
        LevelSetSignStorage::Sign propagatedSign = getSignStorage()->at(signStorageItr);
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

    auto narrowBandItr = m_narrowBandKernel.find(id) ;
    if( narrowBandItr != m_narrowBandKernel.end() ){
        return m_narrowBandGradients.rawAt(narrowBandItr.getRawIndex());
    }

    return levelSetDefaults::GRADIENT;

}

/*! 
 * Deletes non-existing items after grid adaption.
 * @param[in] mapper mapping info
 */
void LevelSetCachedObject::_clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

    // Clear stale narrow band entries
    for (const adaption::Info &adaptionInfo : mapper) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        for (long cellId : adaptionInfo.previous) {
            if (containsNarrowBandEntry(cellId)) {
                eraseNarrowBandEntry(cellId, false);
            }
        }
    }

    syncNarrowBandStorages();


}

/*! 
 * Clears all levelset information
 */
void LevelSetCachedObject::_clear( ){

    // Clear narrow band entries
    m_narrowBandKernel.clear() ;
    syncNarrowBandStorages() ;

    // Clear sign propgation storage
    clearSignStorage();
}

/*!
 * If cell lies within the narrow band and hence its levelset is computed
 * exactly.
 * @param[in] id cell id
 * @return true/false if the centroid is in narrow band
 */
bool LevelSetCachedObject::isInNarrowBand(long id)const{
    return containsNarrowBandEntry(id);
}

/*!
 * Writes LevelSetCachedObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_dump( std::ostream &stream ){

    // Narrow band kernel
    m_narrowBandKernel.dump( stream );

    // Narrow band storages
    m_narrowBandValues.dump( stream );
    m_narrowBandGradients.dump( stream );

    // Stored sign
    dumpSignStorage( stream );
}

/*!
 * Reads LevelSetCachedObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_restore( std::istream &stream ){

    // Narrow band kernel
    m_narrowBandKernel.restore( stream );

    // Narrow band storages
    m_narrowBandValues.restore( stream );
    m_narrowBandGradients.restore( stream );

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

    // Evaluate the size of the buffer
    std::size_t nInfoItems = 0;
    for( long id : sendList){
        if( containsNarrowBandEntry(id)){
            ++nInfoItems ;
        }
    }

    dataBuffer.setSize(dataBuffer.getSize() + sizeof(std::size_t) + nInfoItems* getNarrowBandEntryBinarySize() ) ;

    // Fill the buffer
    dataBuffer << nInfoItems ;

    for( std::size_t k = 0; k < sendList.size(); ++k){
        long id = sendList[k];
        NarrowBandIterator narrowBandItr = getNarrowBandEntryIterator(id) ;
        if( narrowBandItr == m_narrowBandKernel.end() ){
            continue;
        }

        // Write the index of the cell
        dataBuffer << k ;

        // Write the narrowband entry
        writeNarrowBandEntryCommunicationBuffer(narrowBandItr, dataBuffer);
    }
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

    std::size_t nInfoItems ;
    dataBuffer >> nInfoItems ;

    for( std::size_t i=0; i<nInfoItems; ++i){
        // Read the id of the cell
        std::size_t k ;
        dataBuffer >> k ;
        long id = recvList[k] ;

        // Create the narrowband entry
        NarrowBandIterator narrowBandItr = getNarrowBandEntryIterator(id) ;
        if( narrowBandItr == m_narrowBandKernel.end() ){
            narrowBandItr = createNarrowBandEntry(id);
        }

        // Fill the information associated to the narrowband entry
        readNarrowBandEntryCommunicationBuffer(narrowBandItr, dataBuffer);
    }
}

#endif

/*!
 * Create a new narrow band entry.
 *
 * \param id is id of the cell
 * \param sync controls if the kernel will be synched
 * \result An iterator pointing to the newly created narrowband kernel entry.
 */
LevelSetCachedObject::NarrowBandIterator LevelSetCachedObject::createNarrowBandEntry(long id, bool sync) {

    NarrowBandKernel::FillAction action = m_narrowBandKernel.fillHead(id) ;
    if (sync) {
        syncNarrowBandStorages() ;
    }

    return m_narrowBandKernel.rawFind(action.info[PiercedSyncAction::INFO_POS]);

}

/*!
 * Delete a narrow band entry.
 *
 * \param id is id of the cell
 * \param sync controls if the kernel will be synched
 */
void LevelSetCachedObject::eraseNarrowBandEntry(long id, bool sync) {

    m_narrowBandKernel.erase(id, sync) ;
    if (sync) {
        syncNarrowBandStorages() ;
    }

}

/*!
 * Checks if the object contains a narrow band entry for the specified cell.
 *
 * \param id is id of the cell
 * \result Return true is the object contains a narrow band entry for the
 * specified cell, false otherwise.
 */
bool LevelSetCachedObject::containsNarrowBandEntry(long id) const {

    return m_narrowBandKernel.contains(id) ;

}

/*!
 * Get an iterator pointing to the narrowband entry associated with the
 * specified cell.
 *
 * \param id is id of the cell
 * \result An iterator pointing to the narrowband entry associated with the
 * specified cell.
 */
LevelSetCachedObject::NarrowBandIterator LevelSetCachedObject::getNarrowBandEntryIterator(long id) const {

    return m_narrowBandKernel.find(id);

}

/*!
 * Set the information associated with the specified narrow band entry.
 *
 * \param itr is an iterator pointing to the narrowband entry
 * \param value is the levelset value
 * \param gradient is the levelset gradient
 */
void LevelSetCachedObject::setNarrowBandEntry(NarrowBandIterator itr, double value, const std::array<double, 3> &gradient) {

    std::size_t rawId = itr.getRawIndex();

    m_narrowBandValues.rawAt(rawId)    = value;
    m_narrowBandGradients.rawAt(rawId) = gradient;

}

/*!
 * Synchronize narrow band storages with the kernel.
 */
void LevelSetCachedObject::syncNarrowBandStorages() {

    m_narrowBandKernel.flush() ;
    m_narrowBandKernel.sync() ;

}

#if BITPIT_ENABLE_MPI
/*!
 * Get the size, expressed in bytes, of a narrowband entry.
 *
 * \result The size, expressed in bytes, of a narrowband entry.
 */
std::size_t LevelSetCachedObject::getNarrowBandEntryBinarySize() const {

    std::size_t entrySize = sizeof(double) + sizeof(std::array<double, 3>) ;

    return entrySize;

}

/*!
 * Write narrow band entry into the communication buffer.
 *
 * \param narrowBandItr is an iterator pointing to the narrowband entry
 * \param[in,out] dataBuffer is the communication buffer
 */
void LevelSetCachedObject::writeNarrowBandEntryCommunicationBuffer( NarrowBandIterator narrowBandItr, SendBuffer &dataBuffer ){

    std::size_t narrowBandRawId = narrowBandItr.getRawIndex();

    dataBuffer << m_narrowBandValues.rawAt(narrowBandRawId);
    dataBuffer << m_narrowBandGradients.rawAt(narrowBandRawId);

}

/*!
 * Read narrow band entry from the communication buffer.
 *
 * \param narrowBandItr is an iterator pointing to the narrowband entry
 * \param[in,out] dataBuffer is the communication buffer
 */
void LevelSetCachedObject::readNarrowBandEntryCommunicationBuffer( NarrowBandIterator narrowBandItr, RecvBuffer &dataBuffer ){

    std::size_t narrowBandRawId = narrowBandItr.getRawIndex();

    dataBuffer >> m_narrowBandValues.rawAt(narrowBandRawId);
    dataBuffer >> m_narrowBandGradients.rawAt(narrowBandRawId);

}
#endif

}
