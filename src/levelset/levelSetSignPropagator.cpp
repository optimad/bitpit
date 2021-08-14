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

# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetBoundedObject.hpp"
# include "levelSetKernel.hpp"
# include "levelSetObject.hpp"
# include "levelSetSignPropagator.hpp"

namespace bitpit {

/*!
 * \ingroup levelset
 * \class LevelSetSignStorage
 * \brief The class LevelSetSignStorage allows to store the levelset sign
 * on the whole mesh.
 */

const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_UNDEFINED = -2;
const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_NEGATIVE  = -1;
const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_ZERO      =  0;
const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_POSITIVE  =  1;

/*!
 * Constructor.
 */
LevelSetSignStorage::LevelSetSignStorage()
    : m_dirty(true), m_storage(1)
{
}

/*!
 * Check if the stored sign is dirty.
 *
 * \result Returns true if the stored sign is dirty, false
 * otherwise.
 */
bool LevelSetSignStorage::isStoredSignDirty() const
{
    return m_dirty;
}

/*!
 * Set the stored sign as dirty.
 *
 * \param dirty is set to true the stored will be set as dirty
 */
void LevelSetSignStorage::setStoredSignDirty(bool dirty)
{
    m_dirty = dirty;
}

/*!
 * Check if the sign storage has been initialized.
 *
 * \result Return true if the sign storage has been initialized, false
 * otherwise.
 */
bool LevelSetSignStorage::isSignStorageInitialized() const
{
    return (m_storage.getKernel() != nullptr);
}

/*!
 * Initialize the storage.
 *
 * \param cellKernel is the kernel of the cells
 */
void LevelSetSignStorage::initializeSignStorage(bitpit::PiercedKernel<long> *cellKernel)
{
    m_storage.setDynamicKernel(cellKernel, PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
}

/*!
 * Clear the storage.
 *
 * \param release if it's true the memory hold by the container will
 * be released, otherwise the container will be cleared but its
 * memory will not be relased
 */
void LevelSetSignStorage::clearSignStorage(bool release)
{
    m_storage.unsetKernel(release);
    setStoredSignDirty(true);
}

/*!
 * Get the sign of the specified cell.
 *
 * \param id is the id of the cell
 * \result The sign of the specified cell.
 */
LevelSetSignStorage::Sign LevelSetSignStorage::getStoredSign(long id) const
{
    return m_storage.at(id);
}

/*!
 * Get the sign of the specified cell.
 *
 * \param itr is the iterator pointing to the cell
 * \result The sign of the specified cell.
 */
LevelSetSignStorage::Sign LevelSetSignStorage::getStoredSign(const VolumeKernel::CellConstIterator &itr) const
{
    return rawGetStoredSign(itr.getRawIndex());
}

/*!
 * Get the sign of the specified cell.
 *
 * \param rawIndex is the raw index of the cell
 * \result The sign of the specified cell.
 */
LevelSetSignStorage::Sign LevelSetSignStorage::rawGetStoredSign(std::size_t rawIndex) const
{
    return m_storage.rawAt(rawIndex);
}

/*!
 * Set the sign of all the stored cells.
 *
 * \param sign is the sign that will be set
 */
void LevelSetSignStorage::setStoredSign(Sign sign)
{
    m_storage.fill(sign);
}

/*!
 * Set the sign of the specified cell.
 *
 * \param itr is the iterator pointing to the cell
 * \param sign is the sign that will be set
 */
void LevelSetSignStorage::setStoredSign(const VolumeKernel::CellConstIterator &itr, Sign sign)
{
    m_storage.rawAt(itr.getRawIndex()) = sign;
}

/*!
 * Dump storage information.
 *
 * \param stream is the output stream
 */
void LevelSetSignStorage::dumpStoredSign(std::ostream &stream)
{
    utils::binary::write(stream, m_dirty);
    m_storage.dump(stream);
}

/*!
 * Restore storage information.
 *
 * \param stream is the output stream
 */
void LevelSetSignStorage::restoreStoredSign(std::istream &stream)
{
    utils::binary::read(stream, m_dirty);
    m_storage.restore(stream);
}

/*!
    \ingroup levelset
    \class LevelSetSignPropagator
    \brief The class LevelSetSignPropagator allows to propagate the levelset
    sign otuside the narrow band.

    The propagation will start from the cells inside the narrowband and will
    continue

    The sign of the outise the bounding box of all object (external cells) can
    be either positive or negative depending on the orientation of the object.
    There is no need to explicitly propagate the sign into those cells, once
    the sign of the external region is identified, it can be assigned to all
    the cells in the external region. When the propagation reaches the external
    region it can be stopped, the sign of the seed from which the propagation
    has started will be the sign of the external region.
*/

const LevelSetSignPropagator::PropagationState LevelSetSignPropagator::STATE_EXTERNAL = - 1;
const LevelSetSignPropagator::PropagationState LevelSetSignPropagator::STATE_WAITING  =   0;
const LevelSetSignPropagator::PropagationState LevelSetSignPropagator::STATE_REACHED  =   1;

/*!
 * Constructor
 */
LevelSetSignPropagator::LevelSetSignPropagator(VolumeKernel *mesh)
    : m_mesh(mesh)
{
}

/*!
 * Propagate the sign of the signed distance function from narrow band to
 * entire domain.
 *
 * The function assumes that the storage to the propagated sign is not yet
 * initialized. Therefore, the storage will be created and initialized.
 *
 * \param object is the object that whose sign will be propagated
 * \param[in,out] storage is the storage for the propagated sign
 */
void LevelSetSignPropagator::execute(const LevelSetObject *object, LevelSetSignStorage *storage)
{
    // Early return if the stored sign is not dirty
    if (!storage->isStoredSignDirty()) {
        return;
    }

    // Set storage kernel
    storage->initializeSignStorage(&(object->getKernel()->getMesh()->getCells()));

    // Reset stored sign
    storage->setStoredSign(LevelSetSignStorage::SIGN_UNDEFINED);

    // Propagate sign
    propagate(object, storage);

    // The propagated sign is now up-to-date
    storage->setStoredSignDirty(false);
}

/*!
 * Propagate the sign of the signed distance function from narrow band to
 * entire domain.
 *
 * If the storage for the propagated sign has already been initialized, only
 * the entries modified by grid adaptation will be updated. Otherwise the
 * storage will be created and initialized.
 *
 * \param adaptionData are the information about mesh adaption
 * \param object is the object that whose sign will be propagated
 * \param[in,out] storage is the storage for the propagated sign
 */
void LevelSetSignPropagator::execute(const std::vector<adaption::Info> &adaptionData, const LevelSetObject *object, LevelSetSignStorage *storage)
{
    // Early return if the storage has not been initialized
    if (!storage->isSignStorageInitialized()) {
        execute(object, storage);
        return;
    }

    // Early return if the stored sign is not dirty
    if (!storage->isStoredSignDirty()) {
        return;
    }

    // Reset stored sign of new cells
    const PiercedVector<Cell, long> &meshCells = m_mesh->getCells();

    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        for (long cellId : adaptionInfo.current) {
            VolumeKernel::CellConstIterator cellItr = meshCells.find(cellId);
            storage->setStoredSign(cellItr, LevelSetSignStorage::SIGN_UNDEFINED);
        }
    }

    // Propagate sign
    propagate(object, storage);

    // The propagated sign is now up-to-date
    storage->setStoredSignDirty(false);
}

/*!
 * Propagate the sign of the signed distance function from narrow band to
 * entire domain.
 *
 * \param object is the object that whose sign will be propagated
 * \param[in,out] storage is the storage for the propagated sign
 */
void LevelSetSignPropagator::propagate(const LevelSetObject *object, LevelSetSignStorage *storage)
{
    // Initialize propagation information
    initializePropagation(object);

    // Set sign of cells in the narrowband
    //
    // Cells in the narrowband will defined the seed for the propagation.
    VolumeKernel::CellConstIterator cellBegin = m_mesh->cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = m_mesh->cellConstEnd();

    std::vector<std::size_t> rawSeeds;
    for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
        int cellSign = storage->getStoredSign(cellItr);
        if (cellSign == LevelSetSignStorage::SIGN_UNDEFINED) {
            long cellId = cellItr.getId();
            if (object->isInNarrowBand(cellId)) {
                cellSign = object->getSign(cellId);
            }
        }

        if (cellSign != LevelSetSignStorage::SIGN_UNDEFINED) {
            std::size_t cellRawId = cellItr.getRawIndex();
            setSign(cellItr, cellSign, storage);
            rawSeeds.push_back(cellRawId);
        }
    }

    // Use the seeds to propagate the sign
    executeSeedPropagation(rawSeeds, storage);

#if BITPIT_ENABLE_MPI
    // If there are cells with an unknown sign, data communication among
    // ghost cells is needed. However it is only possibly to have cells with
    // an unknown sign for partinioned patches.
    long nGlobalWaiting = m_nWaiting;
    if (m_mesh->isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &nGlobalWaiting, 1, MPI_LONG, MPI_SUM, m_mesh->getCommunicator());
    }

    if (nGlobalWaiting != 0) {
        assert(m_mesh->isPartitioned());
        assert(m_mesh->getProcessorCount() != 1);

        // Initialize the communicator for exchanging the sign of the ghosts
        DataCommunicator dataCommunicator(m_mesh->getCommunicator());

        signed char exchangedSign;
        std::size_t exchangedDataSize = sizeof(exchangedSign);

        // Set the receives
        for (const auto &entry : m_mesh->getGhostCellExchangeTargets()) {
            const int rank = entry.first;
            const auto &list = entry.second;

            dataCommunicator.setRecv(rank, list.size() * exchangedDataSize);
        }

        // Set the sends
        for (const auto &entry : m_mesh->getGhostCellExchangeSources()) {
            const int rank = entry.first;
            auto &list = entry.second;

            dataCommunicator.setSend(rank, list.size() * exchangedDataSize);
        }

        // Communicate sign information among the partitions
        while (nGlobalWaiting != 0) {
            // Start the receives
            for (const auto &entry : m_mesh->getGhostCellExchangeTargets()) {
                const int rank = entry.first;
                dataCommunicator.startRecv(rank);
            }

            // Start the sends
            for (const auto &entry : m_mesh->getGhostCellExchangeSources()) {
                const int rank = entry.first;
                const auto &sendIds = entry.second;
                SendBuffer &buffer = dataCommunicator.getSendBuffer(rank);

                for (long cellId : sendIds) {
                    exchangedSign = LevelSetSignStorage::SIGN_UNDEFINED;
                    if (m_propagationStates.at(cellId) == STATE_REACHED) {
                        exchangedSign = storage->getStoredSign(cellId);
                    }
                    buffer << exchangedSign;
                }

                dataCommunicator.startSend(rank);
            }

            // Receive the sign and propagate the sign
            //
            // If we discover the sign of a ghost, we can use it as a seed.
            rawSeeds.clear();
            int nCompletedRecvs = 0;
            while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
                int rank = dataCommunicator.waitAnyRecv();
                const auto &recvIds = m_mesh->getGhostCellExchangeTargets(rank);
                RecvBuffer &buffer = dataCommunicator.getRecvBuffer(rank);

                // Receive data and detect new seeds
                for (long cellId : recvIds) {
                    buffer >> exchangedSign;
                    if (exchangedSign == LevelSetSignStorage::SIGN_UNDEFINED) {
                        continue;
                    }

                    VolumeKernel::CellConstIterator cellItr = m_mesh->getCells().find(cellId);
                    std::size_t cellRawId = cellItr.getRawIndex();
                    PropagationState cellPropagationState = m_propagationStates.rawAt(cellRawId);
                    if (cellPropagationState == STATE_WAITING) {
                        setSign(cellItr, exchangedSign, storage);
                        rawSeeds.push_back(cellRawId);
                    } else if (cellPropagationState == STATE_REACHED) {
                        assert(object->getSign(cellId) == exchangedSign);
                    }
                }

                ++nCompletedRecvs;
            }

            if (rawSeeds.size() > 0) {
                executeSeedPropagation(rawSeeds, storage);
            }

            // Wait to the sends to finish
            dataCommunicator.waitAllSends();

            // Update the global counter for cells with an unknow sign
            nGlobalWaiting = m_nWaiting;
            MPI_Allreduce(MPI_IN_PLACE, &nGlobalWaiting, 1, MPI_LONG, MPI_SUM, m_mesh->getCommunicator());
        }
    }

    // Communicate the sign of the external region
    //
    // The sign has to be consistent among all the partitions.
    bool exchangeExternalSign;
    if (m_mesh->isPartitioned()) {
        exchangeExternalSign = (m_externalSign != LevelSetSignStorage::SIGN_UNDEFINED);
        MPI_Allreduce(MPI_IN_PLACE, &exchangeExternalSign, 1, MPI_C_BOOL, MPI_LOR, m_mesh->getCommunicator());
    } else {
        exchangeExternalSign = false;
    }

    if (exchangeExternalSign) {
        bool positiveExternalSign = (m_externalSign == 1);
        MPI_Allreduce(MPI_IN_PLACE, &positiveExternalSign, 1, MPI_C_BOOL, MPI_LOR, m_mesh->getCommunicator());

        bool negativeExternalSign = (m_externalSign == -1);
        MPI_Allreduce(MPI_IN_PLACE, &negativeExternalSign, 1, MPI_C_BOOL, MPI_LOR, m_mesh->getCommunicator());

        if (positiveExternalSign && negativeExternalSign) {
            m_externalSign = LevelSetSignStorage::SIGN_UNDEFINED;
        } else if (positiveExternalSign) {
            m_externalSign = 1;
        } else if (negativeExternalSign) {
            m_externalSign = -1;
        } else {
            m_externalSign = LevelSetSignStorage::SIGN_UNDEFINED;
        }
    }
#else
    // Check that the sign has been propagated into all regions
    assert(m_nWaiting == 0);
#endif

    // Assign the sign to the external cells
    if (m_nExternal > 0) {
        // Check if the sign of the external region has been identified
        if (m_externalSign == LevelSetSignStorage::SIGN_UNDEFINED) {
            throw std::runtime_error("Sign of external region not properly identified!");
        }

        // Assign the sign to the cells of the external region
        for (auto cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
            std::size_t cellRawId = cellItr.getRawIndex();
            if (m_propagationStates.rawAt(cellRawId) != STATE_EXTERNAL) {
                continue;
            }

            setSign(cellItr, m_externalSign, storage);
        }
    }

    // Cleanup
    finalizePropagation();
}

/*!
 * Initialize sign propagation
 *
 * \param object is the object that whose sign will be propagated
 */
void LevelSetSignPropagator::initializePropagation(const LevelSetObject *object)
{
    // Initialize propagation state
    m_nWaiting = m_mesh->getCellCount();

    m_propagationStates.setStaticKernel(&(m_mesh->getCells()));
    m_propagationStates.fill(STATE_WAITING);

    // Identify external cells
    //
    // A cell is external if it is completely outside the object bounding box.
    // If the cell may be intersected by the object (i.e., the bounding box of
    // the cell intersects the bounding box of the object), the cell cannot be
    // flagged as external.
    const LevelSetBoundedObject *boundedObject = dynamic_cast<const LevelSetBoundedObject *>(object);

    m_nExternal = 0;
    if (boundedObject) {
        // Get tolerance for distance comparison
        double distanceTolerance = m_mesh->getTol();

        // Evaluate the bounding box of the object
        //
        // The current process may only have the portion of the object needed for
        // evaluating the levelset on the cells of its mesh, therefore we need to
        // evaluate the overall bounding box across all process.
        std::array<double,3> objectBoxMin;
        std::array<double,3> objectBoxMax;
#if BITPIT_ENABLE_MPI
        boundedObject->getGlobalBoundingBox(objectBoxMin, objectBoxMax);
#else
        boundedObject->getBoundingBox(objectBoxMin, objectBoxMax);
#endif

        // Check if the patch intersects the bounding box of the object
        std::array<double,3> patchBoxMin;
        std::array<double,3> patchBoxMax;
        m_mesh->getBoundingBox(patchBoxMin, patchBoxMax);

        bool isPatchIntersected = CGElem::intersectBoxBox(patchBoxMin, patchBoxMax, objectBoxMin, objectBoxMax, 3, distanceTolerance);

        // Detect external cells
        VolumeKernel::CellConstIterator cellBegin = m_mesh->cellConstBegin();
        VolumeKernel::CellConstIterator cellEnd   = m_mesh->cellConstEnd();
        for (VolumeKernel::CellConstIterator itr = cellBegin; itr != cellEnd; ++itr) {
            // Cells inside the narrowband cannot be external
            long cellId = itr.getId();
            if (object->isInNarrowBand(cellId)) {
                continue;
            }

            // Check if the centroid is inside the bounding box
            //
            // Cells with the centroid inside the bounding box of the object
            // cannot be external cells
            if (isPatchIntersected) {
                double geometricTolerance = m_mesh->getTol();
                const std::array<double,3> &cellCentroid = object->getKernel()->computeCellCentroid(cellId);

                bool isCentroidInternal = true;
                for (int i = 0; i < 3; ++i) {
                    if (cellCentroid[i] < objectBoxMin[i] - geometricTolerance || cellCentroid[i] > objectBoxMin[i] + geometricTolerance) {
                        isCentroidInternal = false;
                        break;
                    }
                }

                if (isCentroidInternal) {
                    continue;
                }
            }

            // Check if the cell is inside the bounding box of the object
            std::array<double,3> cellBoxMin;
            std::array<double,3> cellBoxMax;
            m_mesh->evalCellBoundingBox(cellId, &cellBoxMin, &cellBoxMax);

            bool isCellIntersected = CGElem::intersectBoxBox(cellBoxMin, cellBoxMax, objectBoxMin, objectBoxMax, 3, distanceTolerance);
            if (isCellIntersected) {
                continue;
            }

            std::size_t cellRawId = itr.getRawIndex();
            m_propagationStates.rawAt(cellRawId) = STATE_EXTERNAL;
            ++m_nExternal;
        }
        m_nWaiting -= m_nExternal;
    }

    // Initialize the sign of the external region
    m_externalSign = LevelSetSignStorage::SIGN_UNDEFINED;
}

/*!
 * Propagate the sign of the levelset from the specified seeds.
 *
 * Sign will be propagated into both interior and ghost cells of the current
 * process.
 *
 * The sign will NOT be propagated into cells flagged with "EXTERNAL" state
 * (i.e., cells outside the bounding box of all the objects). When propagation
 * reaches the external region, it will be stopped. The sign of the seed from
 * which the propagation has started will define the sign of the external region.
 *
 * \param rawSeeds are the raw ids of the cells that will be used as seeds
 * for the propagation
 * \param[in,out] storage is the storage for the propagated sign
 */
void LevelSetSignPropagator::executeSeedPropagation(const std::vector<std::size_t> &rawSeeds, LevelSetSignStorage *storage)
{
    const PiercedVector<Cell, long> &meshCells = m_mesh->getCells();

    std::vector<std::size_t> rawProcessList;

    std::size_t rawSeedCursor = rawSeeds.size();
    while (rawSeedCursor != 0) {
        // Get a seed
        --rawSeedCursor;
        std::size_t seedRawId = rawSeeds[rawSeedCursor];

        // Get the sign of the seed
        LevelSetSignStorage::Sign seedSign = storage->rawGetStoredSign(seedRawId);
        assert(seedSign >= -1 && seedSign <= 1);

        // Initialize the process list with the seed
        rawProcessList.resize(1);
        rawProcessList[0] = seedRawId;

        // Propagate the sign
        while (!rawProcessList.empty()) {
            std::size_t cellRawId = rawProcessList.back();
            rawProcessList.resize(rawProcessList.size() - 1);

            // Cell information
            VolumeKernel::CellConstIterator cellItr = meshCells.rawFind(cellRawId);

            // Set the sign of the cell
            //
            // We need to set the sign only if it has not already been set
            // (for example, the sign of the seeds is already set).
            LevelSetSignStorage::Sign cellSign = storage->getStoredSign(cellItr);
            if (cellSign == LevelSetSignStorage::SIGN_UNDEFINED) {
                setSign(cellItr, seedSign, storage);
            }

            // Process cell neighbours
            //
            // If a neighbour is waiting for the propagation, add it to the
            // process list. When the propagation reaches an external cell
            // the sign of the seed frow which the propagation started will
            // be the sign of the external region.
            const Cell &cell = *cellItr;
            const long *cellNeighs = cell.getAdjacencies();
            int nCellNeighs = cell.getAdjacencyCount();
            for(int n = 0; n < nCellNeighs; ++n){
                long neighId = cellNeighs[n];
                VolumeKernel::CellConstIterator neighItr = meshCells.find(neighId);
                std::size_t neighRawId = neighItr.getRawIndex();

                PropagationState neighState = m_propagationStates.rawAt(neighRawId);
                if (neighState == STATE_WAITING) {
                    rawProcessList.push_back(neighRawId);
                } else if (neighState == STATE_EXTERNAL) {
                    // If the sign of the external region is unknown it can
                    // be assigned, otherwise check if the current sign is
                    // consistent with the previously evaluated sign.
                    if (m_externalSign == LevelSetSignStorage::SIGN_UNDEFINED) {
                        m_externalSign = seedSign;
                    } else if (m_externalSign != seedSign) {
                        throw std::runtime_error("Mismatch in sign of external region!");
                    }
                }
            }

            // Check if the propagation is complete
            //
            // It can be possible to stop the propagation without processing
            // all the cells in the process list if:
            //  - all cells have been reached by the propagation;
            //  - the sign of the external region have been identified.
            bool emptyWaitingList       = (m_nWaiting == 0);
            bool externalSignIdentified = (m_nExternal == 0) || (m_externalSign != LevelSetSignStorage::SIGN_UNDEFINED);
            if (emptyWaitingList && externalSignIdentified) {
                break;
            }
        }
    }
}

/*!
 * Finalize information used for sign propagation
 */
void LevelSetSignPropagator::finalizePropagation()
{
    // Reset propagation state
    m_propagationStates.unsetKernel(false);
}

/*!
 * Set the sign associated with the cell and update propagation information
 *
 * \param cellId is the id of the cell
 * \param cellSign is the sign associated with the cell
 * \param[out] storage is the storage for the propagated sign
 */
void LevelSetSignPropagator::setSign(const VolumeKernel::CellConstIterator &cellItr, LevelSetSignStorage::Sign cellSign, LevelSetSignStorage *storage)
{
    assert(cellSign >= -1 && cellSign <= 1);

    // Set the sign of the cell
    storage->setStoredSign(cellItr, cellSign);

    // Update the state of the cell
    //
    // If the cell is external, the sign of the eternal region should be
    // updated.
    //
    // If the cell is waiting for the propagation, the counter of the
    // waiting cell list should be updated.
    std::size_t cellRawId = cellItr.getRawIndex();
    PropagationState *cellState = m_propagationStates.rawData(cellRawId);
    if (*cellState == STATE_EXTERNAL) {
        if (m_externalSign == LevelSetSignStorage::SIGN_UNDEFINED) {
            m_externalSign = cellSign;
        } else if (m_externalSign != cellSign) {
            throw std::runtime_error("Mismatch in sign of external region!");
        }
    } else if (*cellState == STATE_WAITING) {
        --m_nWaiting;
    }
    *cellState = STATE_REACHED;
}

}