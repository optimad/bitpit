/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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
# endif

# include <stack>
# include "levelSet.hpp"

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

/*!
	@ingroup levelset
	@interface LevelSetObject
	@brief Interface class for all objects with respect to whom the levelset function may be computed.
*/

/*!
 * Destructor
 */
LevelSetObject::~LevelSetObject( ){
};

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetObject::LevelSetObject( int id, bool primary) : m_id(id), m_primary(primary){
};

/*!
 * Get the id 
 * @return id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
};

/*!
 * If the levelset is primary (e.g. of a surface triangulation) or not (e.g. derived by boolean operations between two levelsets)
 * @return if object is primary
 */
bool LevelSetObject::isPrimary( ) const {
    return m_primary ;
};

/*!
 * Returns reference to LevelSetInfo
*/
PiercedVector<LevelSetInfo>& LevelSetObject::getLevelSetInfo(){
    return m_ls ;
} 

/*!
 * Returns reference to LevelSetInfo
*/
LevelSetInfo LevelSetObject::getLevelSetInfo( const long &i)const{
    if( ! m_ls.exists(i) ){
        return (  LevelSetInfo() );
    } else {
        return m_ls[i] ;
    };

} 

/*!
 * Get the Sdf value of the i-th local element of the octree mesh.
 * @param[in] i cell index
 * @return levelset value in cell
 */
double LevelSetObject::getLS( const long &i)const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::VALUE;
    } else {
        return (  m_ls[i].value );
    };

};

/*!
 * Get the Sdf gradient vector of the i-th local element of the octree mesh.
 * @param[in] i cell index
 * @return levelset gradient in cell 
 */
std::array<double,3> LevelSetObject::getGradient(const long &i) const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::GRADIENT;
    } else {
        return (  m_ls[i].gradient );
    };

};


/*!
 * Get the object and part id of projection point
 * @param[in] i cell index
 * @return pair containing object and part id 
 */
int LevelSetObject::getPart(const long &i) const {

    BITPIT_UNUSED(i) ;
    return levelSetDefaults::PART ;

};

/*!
 * Get the sign of the levelset function
 * @param[in] i cell index
 * @return sign of levelset
 */
short LevelSetObject::getSign(const long &i)const{

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::SIGN;
    } else {
        return ( static_cast<short>(sign( m_ls[i].value)) );
    };

};


/*!
 * If cell centroid lies within the narrow band and hence levelset is computet exactly
 * @param[in] i cell index
 * @return true/false if the centroid is in narrow band
 */
bool LevelSetObject::isInNarrowBand(const long &i)const{

    if( ! m_ls.exists(i) ){
        return false;
    } else {
        return ( std::abs( m_ls[i].value) <= m_RSearch );
    };

};

/*!
 * Get the current size of the narrow band.
 * @return size of the current narrow band
 */
double LevelSetObject::getSizeNarrowBand()const{
    return m_RSearch;
};

/*!
 * Manually set the size of the narrow band.
 * @param[in] r size of the narrow band.
 */
void LevelSetObject::setSizeNarrowBand(double r){
    m_RSearch = r;
};

/*! 
 * Deletes non-existing items and items outside the narrow band after grid adaption.
 * @param[in] mapper mapping info
 */
void LevelSetObject::clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

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

    clearAfterMeshAdaptionDerived( mapper ) ;

    return ;
};


/*!
 * Clears data structure after mesh modification
 * @param[in] mapper mapper describing mesh modifications
 */
void LevelSetObject::clearAfterMeshAdaptionDerived( const std::vector<adaption::Info> &mapper ){
    BITPIT_UNUSED(mapper) ;
    return;
}


/*! 
 * Deletes items outside the narrow band after grid adaption.
 * @param[in] newRSearch new size of narrow band
 */
void LevelSetObject::filterOutsideNarrowBand( double newRSearch ){

    PiercedIterator<LevelSetInfo> lsItr = m_ls.begin() ;
    while( lsItr != m_ls.end() ){


        if( std::abs(lsItr->value) > newRSearch ){
            lsItr = m_ls.erase( lsItr.getId(), true );
        } else {
            ++lsItr ;
        }

    };

    m_ls.flush() ;

    filterOutsideNarrowBandDerived(newRSearch) ;

    return ;
};


/*!
 * Clears data structure outside narrow band
 * @param[in] search size of narrow band
 */
void LevelSetObject::filterOutsideNarrowBandDerived( double search ){
    BITPIT_UNUSED(search) ;
    return;
}


/*! 
 * Clears all levelset information
 */
void LevelSetObject::clear( ){
    m_ls.clear() ;
    clearDerived() ;
}

/*! 
 * Clears all levelset information stored in derived class
 */
void LevelSetObject::clearDerived( ){
    return ;
}

/*!
 * Gets the number of support items within the narrow band of cell
 * @param[in] id index of cell
 */
int LevelSetObject::getSupportCount( const long &id )const{
    BITPIT_UNUSED( id) ;
    return 0;
}

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
long LevelSetObject::getSupport( const long &id )const{
    BITPIT_UNUSED(id) ;
    return levelSetDefaults::SUPPORT ;
}

/*!
 * Propagate the sign of the signed distance function from narrow band to entire domain
 */
void LevelSetObject::propagateSign( LevelSetKernel *visitee ) {

    VolumeKernel const &mesh = *(visitee->getMesh()) ;

    // Save the bounding boxes of the object
    std::array<double,3> boxMin;
    std::array<double,3> boxMax;
    getBoundingBox(boxMin, boxMax);

    // Identify the regions
    //
    // A region is a group of cells that share at leas a face. There are three
    // type of regions: external, internal, narrowband, and unknown.
    //
    // An 'external' region is outside all bodies, an 'internal' region is
    // inside of at least one body, a 'narrowband' region contains at least
    // one cell in the narrow band, whereas it is not easly possible to find
    // out the position of an 'unknown' region.
    const int REGION_EXTERNAL   = 0;
    const int REGION_INTERNAL   = 1;
    const int REGION_NARROWBAND = 2;
    const int REGION_UNKNOWN    = 3;

    const std::vector<long> &cellIds = mesh.getCells().getIds();

    int nRegions = 0;
    std::vector<int> regionType;
    std::vector<std::unordered_set<long>> regionCellList;
    std::unordered_set<long> alreadyAssigned;
    for (long id : cellIds) {
        // Skip the cell if is already assigned to a region
        if (alreadyAssigned.count(id) != 0) {
            continue;
        }

        // Creat a new region
        int region = nRegions++;
        regionType.push_back(REGION_EXTERNAL);
        regionCellList.emplace_back();

        // Find all the cells that belongs to the region
        std::stack<long> seeds;
        seeds.push(id);

        while (!seeds.empty()) {
            long seed = seeds.top();
            seeds.pop();

            // Seeds assigned to a region has already been processed
            if (alreadyAssigned.count(seed) != 0) {
                continue;
            }

            // Assign the current region to the seed
            regionCellList[region].insert(seed);

            // Detect the region type
            if (regionType[region] != REGION_NARROWBAND) {
                if (isInNarrowBand(seed)) {
                    regionType[region] = REGION_NARROWBAND;
                } else if (regionType[region] != REGION_UNKNOWN) {
                    std::array<double, 3> centroid = mesh.evalCellCentroid(seed);

                    bool isInside = true;
                    for (int i = 0; i < 3; ++i) {
                        if (centroid[i] < boxMin[i] || centroid[i] > boxMax[i]) {
                            isInside = false;
                            break;
                        }
                    }

                    if (isInside) {
                        regionType[region] = REGION_UNKNOWN;
                    }
                }
            }

            // Add the unassigned neighboors to the seeds
            for (long neigh : mesh.findCellFaceNeighs(seed)) {
                if (alreadyAssigned.count(neigh) == 0) {
                    seeds.push(neigh);
                }
            }

            // The cell has been assigned to a region
            alreadyAssigned.insert(seed);
        }
    }

    // According to the region type propagate the sign.
    //
    // If a region is 'external', all its cells have a positive levelset. If
    // a region is 'narrowband' it is possible to propagate the sign of the
    // narrowband to the whole region. All the cells of an 'unknown' region
    // will have the same sign, but the value hasto be obtain from the
    // neighbouring partitions.
    for (int region = 0; region < nRegions; ++region) {
        if (regionType[region] == REGION_EXTERNAL) {
            assignSign(1, regionCellList[region]);
        } else if (regionType[region] == REGION_NARROWBAND) {
            // The seeds are the cells in the narrowband of the region
            std::stack<long> seeds;
            std::unordered_set<long> alreadyEvaluated;
            for (long id : regionCellList[region]) {
                if (!isInNarrowBand(id)) {
                    continue;
                }

                seeds.push(id);
                alreadyEvaluated.insert(id);
            }

            // Propagate the sign
            while (!seeds.empty()) {
                long seed = seeds.top();
                seeds.pop();

                // Discard seeds with a LS value equal to 0
                double ls = getLS(seed);
                if (utils::DoubleFloatingEqual()(std::abs(ls), (double) 0.)) {
                    continue;
                }

                // Start filling the process list with the neighbours of the
                // seeds that has not yet been evaluated.
                //
                // If a seed is surrounded only by items already evaluated,
                // it can't propagate the sign to anyone.
                std::stack<long> processList;
                for (long neigh : mesh.findCellFaceNeighs(seed)) {
                    if (alreadyEvaluated.count(neigh) == 0) {
                        processList.push(neigh);
                    }
                }

                if (processList.empty()) {
                    continue;
                }

                // Get the sign of the seed
                short seedSign = getSign(seed);

                // Propagate the sign
                while (!processList.empty()) {
                    long id = processList.top();
                    processList.pop();

                    // Get the value associated to the id
                    //
                    // A new value needs to be created only if the sign to
                    // propagate is different from the default sign.
                    PiercedVector<LevelSetInfo>::iterator infoItr = m_ls.find(id) ;
                    if( infoItr == m_ls.end() && seedSign != levelSetDefaults::SIGN ){
                        infoItr = m_ls.reclaim(id) ;
                    }

                    // Update the value
                    if( infoItr != m_ls.end() ) {
                        (*infoItr).value = seedSign * levelSetDefaults::VALUE;
                    }

                    // Add non-evaluated neighs to the process list
                    std::vector<long> neighs = mesh.findCellFaceNeighs( id ) ;
                    for (long neigh : neighs) {
                        if (alreadyEvaluated.count(neigh) == 0) {
                            processList.push(neigh);
                        }
                    }

                    // The cell has been processeed
                    alreadyEvaluated.insert(id);
                }

                // If all cells have been evaluated we can stop the propagation
                if (alreadyEvaluated.size() == regionCellList[region].size()) {
                    break;
                }
            }
        }
    }

#if BITPIT_ENABLE_MPI
    // If the communicator is not set we can exit because the calculation
    // is serial
    if (!visitee->isCommunicatorSet()) {
        return;
    }

    // If there is only one processor we can exit
    int nProcs;
    MPI_Comm_size(visitee->getCommunicator(), &nProcs);
    if (nProcs == 1) {
        return;
    }

    // Initialize the communicator for exchanging the sign among partitions
    DataCommunicator dataCommunicator(visitee->getCommunicator());
    dataCommunicator.setTag(108) ;

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
    while (true) {
        // If all partitions have propgate the sign on all cells we can exit
        bool locallyComplete = true;
        for (int region = 0; region < nRegions; ++region) {
            if (regionType[region] == REGION_UNKNOWN) {
                locallyComplete = false;
                break;
            }
        }

        bool globallyComplete;
        MPI_Allreduce(&locallyComplete, &globallyComplete, 1, MPI_C_BOOL, MPI_LAND, visitee->getCommunicator());
        if (globallyComplete) {
            return;
        }

        // Start the receives
        for (const auto entry : mesh.getGhostExchangeTargets()) {
            const int rank = entry.first;
            dataCommunicator.startRecv(rank);
        }

        // Start the sends
        for (const auto entry : mesh.getGhostExchangeSources()) {
            const int rank = entry.first;
            auto &list = entry.second;
            SendBuffer &buffer = dataCommunicator.getSendBuffer(rank);

            for (long id : list) {
                // Detect the region associated to the ghost
                int region = -1;
                for (int k = 0; k < nRegions; ++k) {
                    if (regionCellList[k].count(id) != 0) {
                        region = k;
                        break;
                    }
                }

                // Detect the sign to communicate.
                if (regionType[region] == REGION_UNKNOWN) {
                    sign = 0;
                } else if (regionType[region] == REGION_EXTERNAL) {
                    sign = 1;
                } else if (regionType[region] == REGION_INTERNAL) {
                    sign = -1;
                } else {
                    sign = getSign(id);
                }

                buffer << sign;
            }

            dataCommunicator.startSend(rank);
        }

        // Check if it is possible to detect the sign of an unkown region
        std::vector<int> regionSign(nRegions, 0);
        for (int region = 0; region < nRegions; ++region) {
            if (regionType[region] == REGION_EXTERNAL) {
                regionSign[region] = 1;
            } else if (regionType[region] == REGION_INTERNAL) {
                regionSign[region] = -1;
            }
        }

        int nCompletedRecvs = 0;
        while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
            int rank = dataCommunicator.waitAnyRecv();
            const auto &list = mesh.getGhostExchangeTargets(rank);
            RecvBuffer &buffer = dataCommunicator.getRecvBuffer(rank);


            for (long id : list) {
                buffer >> sign;
                if (sign == 0) {
                    continue;
                }

                // Detect the region associated to the ghost
                int region = -1;
                for (int k = 0; k < nRegions; ++k) {
                    if (regionCellList[k].count(id) != 0) {
                        region = k;
                        break;
                    }
                }

                // Assign the sign to the region
                if (regionType[region] != REGION_UNKNOWN) {
                    continue;
                } else if (regionSign[region] != 0) {
                    assert(sign == regionSign[region]);
                    continue;
                }

                regionSign[region] = sign;
            }

            ++nCompletedRecvs;
        }

        dataCommunicator.waitAllSends();

        // Set the sign of the fixable unkown regions
        for (int region = 0; region < nRegions; ++region) {
            if (regionType[region] != REGION_UNKNOWN) {
                continue;
            }

            sign = regionSign[region];
            if (sign == 0) {
                continue;
            }

            // Set the type of the region
            if (sign > 0) {
                regionType[region] = REGION_EXTERNAL;
            } else if (sign < 0) {
                regionType[region] = REGION_INTERNAL;
            }

            // Assign the sign to the whole region
            assignSign(sign, regionCellList[region]);
        }
    }
#else
    // Check that the sign has been propagated into all regions
    bool complete = true;
    for (int region = 0; region < nRegions; ++region) {
        if (regionType[region] == REGION_UNKNOWN) {
            complete = false;
            break;
        }
    }

    assert(complete);
#endif

}

/*!
 * Assign the sign to the specified list of cells.
 *
 * \param sign is the sign that will be assiged to the cells
 * \param cells is the list of cells
 */
void LevelSetObject::assignSign(int sign, const std::unordered_set<long> &cells) {

    for (long id : cells) {
        // Get the info associated to the id
        //
        // A new info needs to be created only if the sign to assign is
        // different from the default sign.
        PiercedVector<LevelSetInfo>::iterator infoItr = m_ls.find(id) ;
        if( infoItr == m_ls.end() && sign != levelSetDefaults::SIGN ){
            infoItr = m_ls.reclaim(id) ;
        }

        // Update the sign
        if( infoItr != m_ls.end() ) {
            (*infoItr).value = sign * levelSetDefaults::VALUE;
        }
    }

}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::ostream &stream ){

    bitpit::PiercedVector<LevelSetInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    IO::binary::write(stream, m_id) ;
    IO::binary::write(stream, m_primary) ;
    IO::binary::write(stream, m_RSearch);
    IO::binary::write(stream, (long) m_ls.size() ) ;

    for( infoItr=m_ls.begin(); infoItr!=infoEnd; ++infoItr){
        IO::binary::write(stream, infoItr.getId()) ;
        IO::binary::write(stream, infoItr->value) ;
        IO::binary::write(stream, infoItr->gradient) ;
    };

    _dump(stream) ;
};

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::istream &stream ){

    long i, n, id;
    LevelSetInfo cellInfo;

    IO::binary::read(stream, m_id) ;
    IO::binary::read(stream, m_primary) ;
    IO::binary::read(stream, m_RSearch);
    IO::binary::read(stream, n);

    m_ls.reserve(n);
    for( i=0; i<n; ++i){
        IO::binary::read(stream, id) ;
        IO::binary::read(stream, cellInfo.value) ;
        IO::binary::read(stream, cellInfo.gradient) ;
        m_ls.insert(id, cellInfo) ;
    };

    _restore(stream) ;
};

#if BITPIT_ENABLE_MPI

/*!
 * Checks if MPI communicator is available in underlying mesh.
 * If available, MPI communicator is retreived from mesh (and duplicated if necessary) and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSetObject::assureMPI(LevelSetKernel *visitee){

    return(visitee->assureMPI() ) ;
}

/*!
 * Update of ghost cell;
 */
void LevelSetObject::exchangeGhosts(LevelSetKernel *visitee){

    std::unordered_map<int,std::vector<long>> &sendList =  visitee->getMesh()->getGhostExchangeSources() ;
    std::unordered_map<int,std::vector<long>> &recvList =  visitee->getMesh()->getGhostExchangeTargets() ;

    communicate(visitee,sendList,recvList);

    return ;
}

/*!
 * communicates data structures of kernel and objects.
 * If mapper!=NULL, clearAfterMeshAdaption of kernel and objects will be called between send and receive.
 * @param[in] visitee LevelSetKernel containing the mesh
 * @param[in] sendList list of elements to be send
 * @param[in] recvList list of elements to be received
 * @param[in] mapper mapper containing mesh modifications
 */
void LevelSetObject::communicate( LevelSetKernel *visitee, std::unordered_map<int,std::vector<long>> &sendList, std::unordered_map<int,std::vector<long>> &recvList, std::vector<adaption::Info> const *mapper){

    if( assureMPI(visitee) ){

        MPI_Comm meshComm = visitee->getCommunicator() ;

        DataCommunicator    sizeCommunicator(meshComm) ; 
        DataCommunicator    dataCommunicator(meshComm) ; 

        int offset = 2*getId() ;
        sizeCommunicator.setTag(offset+0) ;
        dataCommunicator.setTag(offset+1) ;

        /* Start receiving of first communications which contain the sizes
         * of the communications of the exchanged data
        */
	    for (const auto entry : sendList) {
            int rank = entry.first;

            sizeCommunicator.setRecv( rank, sizeof(long) ) ;
            sizeCommunicator.startRecv(rank);
        }

        /* Fill the databuffer with its content from the LevelSetObject 
         * base class and specific derived class. Fill the sizebuffer
         * with the size of the databuffer. Start sending both buffers
         */
	    for (const auto entry : sendList) {

	    	// Get the send buffer
	    	int rank = entry.first;

            sizeCommunicator.setSend(rank, sizeof(long) ) ;
            dataCommunicator.setSend(rank, 0 ) ;

            SendBuffer &sizeBuffer = sizeCommunicator.getSendBuffer(rank);
            SendBuffer &dataBuffer = dataCommunicator.getSendBuffer(rank);

            // Write data in databuffer
            writeCommunicationBuffer( entry.second, dataBuffer ) ;

            // Write size of databuffer in sizebuffer
            sizeBuffer << dataBuffer.capacity() ;

            // Start sending
            sizeCommunicator.startSend(rank);
            dataCommunicator.startSend(rank);
	    }

        
        /* Check which size communications have arrived. For those which are
         * available read the size of the databuffer and start receiving.
         * Proceed once all sizeBuffers have arrived.
         */
        int nCompletedRecvs(0);
        long dataSize;

        while (nCompletedRecvs < sizeCommunicator.getRecvCount()) {
            int rank = sizeCommunicator.waitAnyRecv();
            RecvBuffer &sizeBuffer = sizeCommunicator.getRecvBuffer(rank);

            sizeBuffer >> dataSize ;

            dataCommunicator.setRecv(rank,dataSize) ;
            dataCommunicator.startRecv(rank) ;

            ++nCompletedRecvs;
        }
        

        /* In case of mesh adaption, the ids of the moved items must be freed 
         * in order to avoid conflicts when ids are recycled
         */
        if( mapper!=NULL){
            clearAfterMeshAdaption(*mapper);
        }

        /* Check which data communications have arrived. For those which are
         * available start reading the databuffer into the data structure of
         * LevelSetObject and its derived classes.
         */
        nCompletedRecvs = 0;
        while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
            int rank = dataCommunicator.waitAnyRecv();

            RecvBuffer &dataBuffer = dataCommunicator.getRecvBuffer(rank);
            readCommunicationBuffer( recvList[rank], dataBuffer ) ;

            ++nCompletedRecvs;
        }
        

        sizeCommunicator.waitAllSends();
        sizeCommunicator.finalize() ;

        dataCommunicator.waitAllSends();
        dataCommunicator.finalize() ;

    }

    return ;
}


/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetObject::writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    long nItems(0), counter(0) ;


    //determine number of elements to send
    for( const auto &index : sendList){
        if( m_ls.exists(index)){
            nItems++ ;
        }

    }

    dataBuffer << nItems ;
    dataBuffer.setCapacity(dataBuffer.capacity() +nItems* (sizeof(long) +4*sizeof(double)) ) ;


    for( const auto &index : sendList){
        if( m_ls.exists(index)){
            const auto &lsinfo = m_ls[index] ;
            dataBuffer << counter ;
            dataBuffer << lsinfo.value ;
            dataBuffer << lsinfo.gradient ;
            ++nItems ;
        }
        ++counter ;
    }

    writeCommunicationBufferDerived( sendList, dataBuffer) ;

    dataBuffer.squeeze( ) ;

    return;
};

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetObject::writeCommunicationBufferDerived( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    BITPIT_UNUSED(sendList) ;
    BITPIT_UNUSED(dataBuffer) ;
    return;
};

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetObject::readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    long    nItems, index, id;

    dataBuffer >> nItems ;

    for( int i=0; i<nItems; ++i){

        // Determine the id of the element
        dataBuffer >> index ;
        id = recvList[index] ;

        // Assign the data of the element
        PiercedVector<LevelSetInfo>::iterator infoItr ;
        if( !m_ls.exists(id)){
            infoItr = m_ls.reclaim(id) ;
        } else {
            infoItr = m_ls.getIterator(id) ;
        }

        dataBuffer >> infoItr->value ;
        dataBuffer >> infoItr->gradient ;
    }

    readCommunicationBufferDerived( recvList, dataBuffer ) ;

    return;
};

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetObject::readCommunicationBufferDerived( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    BITPIT_UNUSED(recvList) ;
    BITPIT_UNUSED(dataBuffer) ;
    return;
};

#endif 

}
