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

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# include "bitpit_communications.hpp"
# endif
# include <stack>

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
# include "levelSet.hpp"

namespace bitpit {

const int LevelSetCachedObject::PROPAGATION_STATUS_EXTERNAL = - 1;
const int LevelSetCachedObject::PROPAGATION_STATUS_WAITING  =   0;
const int LevelSetCachedObject::PROPAGATION_STATUS_REACHED  =   1;

/*!
	@ingroup levelset
	@interface LevelSetCachedObject
	@brief Interface class for all objects which need to store the discrete values of levelset function.
*/

/*!
 * Destructor
 */
LevelSetCachedObject::~LevelSetCachedObject( ){
}

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetCachedObject::LevelSetCachedObject(int id) : LevelSetObject(id){
}

/*!
 * Get LevelSetInfo of cell
 * @param[in] i cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetCachedObject::getLevelSetInfo( const long &i)const{

    auto itr = m_ls.find(i);
    if ( itr == m_ls.end() ){
        return LevelSetInfo();
    }

    return *itr;

} 

/*!
 * Get the levelset value of cell
 * @param[in] i cell index
 * @return levelset value in cell
 */
double LevelSetCachedObject::getLS( const long &i)const {

    auto itr = m_ls.find(i);
    if ( itr == m_ls.end() ){
        return levelSetDefaults::VALUE;
    }

    return itr->value;

}

/*!
 * Get the levelset gradient of cell
 * @param[in] i cell index
 * @return levelset gradient in cell 
 */
std::array<double,3> LevelSetCachedObject::getGradient(const long &i) const {

    auto itr = m_ls.find(i);
    if ( itr == m_ls.end() ){
        return levelSetDefaults::GRADIENT;
    }

    return itr->gradient;

}

/*! 
 * Deletes non-existing items after grid adaption.
 * @param[in] mapper mapping info
 */
void LevelSetCachedObject::_clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

    for ( auto & map : mapper ){
        if( map.entity == adaption::Entity::ENTITY_CELL ){
            if( map.type == adaption::Type::TYPE_DELETION || 
                map.type == adaption::Type::TYPE_PARTITION_SEND  ||
                map.type == adaption::Type::TYPE_REFINEMENT  ||
                map.type == adaption::Type::TYPE_COARSENING  ){

                for ( auto & parent : map.previous){
                    if( m_ls.exists(parent) ) 
                        m_ls.erase(parent,true) ;
                }
            }
        }
    }

    m_ls.flush() ;

    __clearAfterMeshAdaption( mapper ) ;

    return ;
}

/*! 
 * Deletes non-existing items after grid adaption in derived classes
 * @param[in] mapper mapping info
 */
void LevelSetCachedObject::__clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){
    BITPIT_UNUSED(mapper);
}

/*! 
 * Clears all levelset information
 */
void LevelSetCachedObject::_clear( ){
    m_ls.clear() ;
    __clear() ;
}

/*! 
 * Clears all levelset information in derived classes
 */
void LevelSetCachedObject::__clear( ){
}

/*!
 * Propagate the sign of the signed distance function from narrow band to
 * entire domain.
 */
void LevelSetCachedObject::propagateSign() {

    VolumeKernel const &mesh = *(m_kernelPtr->getMesh()) ;
    const PiercedVector<Cell> &cells = mesh.getCells();

    VolumeKernel::CellConstIterator cellBegin = mesh.cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh.cellConstEnd();

    std::vector<long> seeds;
    seeds.reserve(mesh.getCellCount());

    // Evaluate the bounding box of the object
    std::array<double,3> boxMin;
    std::array<double,3> boxMax;
    getBoundingBox(boxMin, boxMax);

    // Set the initial propagation status of the cells
    //
    // We don't need to set the sign of cells associated to a levelset info,
    // if a cell is associated to a levelset info we know its sign (the value
    // contained in the levelset info may be a dummy value, but the sign is
    // the correct one) and the cell can be used as seeds for propagating the
    // sign to other cells.
    //
    // The sign of the cells outise the bounding box of all objects (external
    // cells) can be either positive or negative depending on the orientation
    // of the objects. There is no need to propagate the sign into those cells,
    // once the sign of the external region is know, it can be assigned to all
    // the cells in the external region. When the propagation reaches the
    // external region it can be stopped, the sign the sign of the seed from
    // which the propagation has started that will be the sign of the external
    // region.
    PiercedStorage<int, long> propagationStatus(1, &cells);
    propagationStatus.fill(PROPAGATION_STATUS_WAITING);

    int externalSign = 0;
    long nExternal   = 0;
    long nWaiting    = mesh.getCellCount();
    for (auto itr = cellBegin; itr != cellEnd; ++itr) {
        long cellId = itr.getId();

        // Process cells associated to a levelset info
        if (m_ls.find(cellId) != m_ls.end()) {
            std::size_t cellRawId = itr.getRawIndex();
            propagationStatus.rawAt(cellRawId) = PROPAGATION_STATUS_REACHED;
            --nWaiting;

            seeds.push_back(cellId);

            continue;
        }

        // Process external cells
        std::array<double, 3> centroid = mesh.evalCellCentroid(cellId);

        bool isExternal = false;
        for (int i = 0; i < 3; ++i) {
            if (centroid[i] < boxMin[i] || centroid[i] > boxMax[i]) {
                isExternal = true;
                break;
            }
        }

        if (isExternal) {
            std::size_t cellRawId = itr.getRawIndex();
            propagationStatus.rawAt(cellRawId) = PROPAGATION_STATUS_EXTERNAL;
            --nWaiting;
            ++nExternal;

            continue;
        }
    }

    // Use the seeds to propagate the sign
    propagateSeedSign(seeds, &propagationStatus, &nWaiting, &externalSign);

#if BITPIT_ENABLE_MPI
    // If there are cells with an unknown sign, data communication among
    // ghost cells is needed. However it is only possibly to have cells with
    // an unknown sign for partinioned patches.
    long nGlobalWaiting = nWaiting;
    if (mesh.isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &nGlobalWaiting, 1, MPI_LONG, MPI_SUM, m_kernelPtr->getCommunicator());
    }

    if (nGlobalWaiting != 0) {
        assert(mesh.isPartitioned());
        assert(mesh.getProcessorCount() != 1);

        // Initialize the communicator for exchanging the sign of the ghosts
        DataCommunicator dataCommunicator(m_kernelPtr->getCommunicator());

        int sign;
        size_t dataSize = sizeof(sign);

        // Set the receives
        for (const auto entry : mesh.getGhostExchangeTargets()) {
            const int rank = entry.first;
            const auto &list = entry.second;

            dataCommunicator.setRecv(rank, list.size() * dataSize);
        }

        // Set the sends
        for (const auto entry : mesh.getGhostExchangeSources()) {
            const int rank = entry.first;
            auto &list = entry.second;

            dataCommunicator.setSend(rank, list.size() * dataSize);
        }

        // Communicate the sign among the partitions
        while (nGlobalWaiting != 0) {
            // Start the receives
            for (const auto entry : mesh.getGhostExchangeTargets()) {
                const int rank = entry.first;
                dataCommunicator.startRecv(rank);
            }

            // Start the sends
            for (const auto entry : mesh.getGhostExchangeSources()) {
                const int rank = entry.first;
                auto &sendIds = entry.second;
                SendBuffer &buffer = dataCommunicator.getSendBuffer(rank);

                for (long cellId : sendIds) {
                    int sign = 0;
                    if (propagationStatus.at(cellId) == PROPAGATION_STATUS_REACHED) {
                        sign = getSign(cellId);
                    }
                    buffer << sign;
                }

                dataCommunicator.startSend(rank);
            }

            // Receive the sign and propagate the sign
            //
            // If we discover the sign of a ghost, we can use it as a seed.
            seeds.clear();
            int nCompletedRecvs = 0;
            while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
                int rank = dataCommunicator.waitAnyRecv();
                const auto &recvIds = mesh.getGhostExchangeTargets(rank);
                RecvBuffer &buffer = dataCommunicator.getRecvBuffer(rank);

                // Receive data and detect new seeds
                for (long cellId : recvIds) {
                    buffer >> sign;

                    std::size_t cellRawId = cells.getRawIndex(cellId);
                    int &cellPropagationStatus = propagationStatus.rawAt(cellRawId);
                    if (cellPropagationStatus == PROPAGATION_STATUS_WAITING) {
                        if (sign != 0) {
                            setSign(cellId, sign);
                            cellPropagationStatus = PROPAGATION_STATUS_REACHED;
                            --nWaiting;
                            seeds.push_back(cellId);
                        }
                    } else if (cellPropagationStatus == PROPAGATION_STATUS_REACHED) {
                        assert(getSign(cellId) == sign);
                    }
                }

                ++nCompletedRecvs;
            }

            if (seeds.size() > 0) {
                propagateSeedSign(seeds, &propagationStatus, &nWaiting, &externalSign);
            }

            // Wait to the sends to finish
            dataCommunicator.waitAllSends();

            // Update the global counter for cells with an unknow sign
            nGlobalWaiting = nWaiting;
            MPI_Allreduce(MPI_IN_PLACE, &nGlobalWaiting, 1, MPI_LONG, MPI_SUM, m_kernelPtr->getCommunicator());
        }
    }

    // Communicate the sign of the external region
    //
    // The sign has to be consistent among all the partitions.
    if (mesh.isPartitioned()) {
        int minExternalSign = externalSign;
        MPI_Allreduce(MPI_IN_PLACE, &minExternalSign, 1, MPI_INT, MPI_MIN, m_kernelPtr->getCommunicator());

        int maxExternalSign = externalSign;
        MPI_Allreduce(MPI_IN_PLACE, &maxExternalSign, 1, MPI_INT, MPI_MAX, m_kernelPtr->getCommunicator());

        if (minExternalSign == maxExternalSign) {
            externalSign = minExternalSign;
        } else if (minExternalSign == 0) {
            externalSign = maxExternalSign;
        } else if (maxExternalSign == 0) {
            externalSign = minExternalSign;
        }
    }
#else
    // Check that the sign has been propagated into all regions
    assert(nWaiting == 0);
#endif

    // Check that the sign of the external region was correctly identified
    if (nExternal != 0 && externalSign == 0) {
        throw std::runtime_error("Sign of external region not properly identified!");
    }

    // Assign the sign to the external cells
    for (auto itr = cellBegin; itr != cellEnd; ++itr) {
        std::size_t cellRawId = itr.getRawIndex();
        if (propagationStatus.rawAt(cellRawId) != PROPAGATION_STATUS_EXTERNAL) {
            continue;
        }

        long cellId = itr.getId();
        setSign(cellId, externalSign);
    }
}

/*!
 * Propagates the sign of the signed distance function from the specified
 * seeds to all reachable cells whose sign has not yet been assigned.
 *
 * The sign will NOT be propagated into cells flagged with the "EXTERNAL"
 * status (i.e., cells outside the bounding box of all the objects). When
 * the propagation reaches the external region it can be stopped, the sign
 * the sign of the seed from which the propagation has startedthat will be
 * the sign of the external region.
 *
 * \param seeds are the seeds to be used for sign propagation
 * \param[in,out] statuses is a flag that defines the propagation status of
 * the cells. On output, this flag will be updated with the new propagation
 * statuses
 * \param[in,out] nWaiting is the number of cells that are waiting for the
 * propagation to reach them. On output, this number will be updated so it's
 * possible to keep track of the cells whose sign is not yet assigned
 * \param[in,out] externalSign is the sign of the external region. If the
 * external sign is not yet defined and the propagation reaches the external
 * region, the parameter will be updated with the sign of the external region
 */
void LevelSetCachedObject::propagateSeedSign(const std::vector<long> &seeds,
                                             PiercedStorage<int, long> *statuses,
                                             long *nWaiting, int *externalSign) {

    VolumeKernel const &mesh = *(m_kernelPtr->getMesh()) ;

    std::vector<long> processList;

    std::size_t seedCursor = seeds.size();
    while (seedCursor != 0) {
        // Get a seed
        --seedCursor;
        long seedId = seeds[seedCursor];

        // Get the sign of the seed
        int seedSign = getSign(seedId);

        // Initialize the process list with the seed
        processList.resize(1);
        processList[0] = seedId;

        // Propagate the sign
        while (!processList.empty()) {
            long cellId = processList.back();
            processList.resize(processList.size() - 1);

            // Process non-seed cells
            bool isSeed = (cellId == seedId);
            if (!isSeed) {
                // Consider only cells waiting for the propagation to reach them
                //
                // The process list is not unique, it can contain duplicate, so
                // we need to make sure not to process a cell multiple times.
                int &cellStatus = statuses->at(cellId);
                if (cellStatus != PROPAGATION_STATUS_WAITING) {
                    continue;
                }

                // Set the sign of the cell
                setSign(cellId, seedSign);
                cellStatus = PROPAGATION_STATUS_REACHED;
                --(*nWaiting);

                // If there are no more waiting cell we can exit
                if (*nWaiting == 0) {
                    return;
                }
            }

            // Process cell neighbours
            //
            // If a neighbour is waiting for the propagation, add it to the
            // process list. When the propagation reaches an external cell
            // the sign of the seed frow which the propagation started will
            // be the sign of the external region.
            const Cell &cell = mesh.getCell(cellId);
            const long *cellNeighs = cell.getAdjacencies() ;
            int nCellNeighs = cell.getAdjacencyCount() ;
            for(int n = 0; n < nCellNeighs; ++n){
                long neighId = cellNeighs[n] ;
                if (neighId < 0) {
                    continue;
                }

                int neighStatus = statuses->at(neighId);
                if (neighStatus == PROPAGATION_STATUS_WAITING) {
                    processList.push_back(neighId);
                } else if (neighStatus == PROPAGATION_STATUS_EXTERNAL) {
                    // If the sign of the external region was unknown it can
                    // be assigned, otherwise check if the current sign is
                    // consistent with the previously evaluated sign.
                    if (*externalSign == 0) {
                        *externalSign = seedSign;
                    } else if (*externalSign != seedSign) {
                        throw std::runtime_error("Mismatch in sign of external region!");
                    }
                }
            }
        }
    }
}

/*!
 * Set the sign of the specified cell.
 *
 * \param id is the id of the cell
 * \param sign is the sign that will be assiged to the cell
 */
void LevelSetCachedObject::setSign(long id, int sign) {

    // The sign has to be set only if it is different from the default sign
    if (sign == levelSetDefaults::SIGN) {
        return;
    }

    // Get the info associated to the id
    PiercedVector<LevelSetInfo>::iterator infoItr = m_ls.find(id);
    if (infoItr == m_ls.end()){
        infoItr = m_ls.emplace(id);
    }

    // Update the sign
    (*infoItr).value = sign * levelSetDefaults::VALUE;
}

/*!
 * Writes LevelSetCachedObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_dump( std::ostream &stream ){

    utils::binary::write(stream, (long) m_ls.size() ) ;
    bitpit::PiercedVector<LevelSetInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    for( infoItr=m_ls.begin(); infoItr!=infoEnd; ++infoItr){
        utils::binary::write(stream, infoItr.getId()) ;
        utils::binary::write(stream, infoItr->value) ;
        utils::binary::write(stream, infoItr->gradient) ;
    }

    __dump(stream) ;
}

/*!
 * Writes information of derived class to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::__dump( std::ostream &stream ){
    BITPIT_UNUSED(stream);
}

/*!
 * Reads LevelSetCachedObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_restore( std::istream &stream ){

    long i, n, id;
    LevelSetInfo cellInfo;

    utils::binary::read(stream, n);

    m_ls.reserve(n);
    for( i=0; i<n; ++i){
        utils::binary::read(stream, id) ;
        utils::binary::read(stream, cellInfo.value) ;
        utils::binary::read(stream, cellInfo.gradient) ;
        m_ls.insert(id, cellInfo) ;
    }

    __restore(stream) ;
}

/*!
 * Reads information of derived class from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::__restore( std::istream &stream ){
    BITPIT_UNUSED(stream);
}

#if BITPIT_ENABLE_MPI

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetCachedObject::_writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    // Evaluate the size of the buffer
    long nItems = 0;
    for( const auto &index : sendList){
        if( m_ls.exists(index)){
            nItems++ ;
        }
    }

    dataBuffer.setSize(dataBuffer.getSize() + sizeof(long) + nItems* (sizeof(long) +4*sizeof(double)) ) ;

    // Fill the buffer
    dataBuffer << nItems ;

    long index = 0;
    for( long id : sendList){
        auto lsInfoItr = m_ls.find(id) ;
        if( lsInfoItr != m_ls.end() ){
            dataBuffer << index ;
            dataBuffer << lsInfoItr->value ;
            dataBuffer << lsInfoItr->gradient ;
        }
        ++index ;
    }

    __writeCommunicationBuffer( sendList, dataBuffer) ;
}

/*!
 * Flushing of data of derived class to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetCachedObject::__writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){
    BITPIT_UNUSED(sendList);
    BITPIT_UNUSED(dataBuffer);
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetCachedObject::_readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    long nItems ;
    dataBuffer >> nItems ;

    for( int i=0; i<nItems; ++i){
        // Determine the id of the element
        long index ;
        dataBuffer >> index ;
        long id = recvList[index] ;

        // Assign the data of the element
        PiercedVector<LevelSetInfo>::iterator infoItr = m_ls.find(id) ;
        if( infoItr == m_ls.end() ){
            infoItr = m_ls.emplace(id) ;
        }

        dataBuffer >> infoItr->value ;
        dataBuffer >> infoItr->gradient ;
    }

    __readCommunicationBuffer( recvList, dataBuffer ) ;
}

/*!
 * Processing of communication buffer of derived class into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetCachedObject::__readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){
    BITPIT_UNUSED(recvList);
    BITPIT_UNUSED(dataBuffer);
}

#endif 

}
