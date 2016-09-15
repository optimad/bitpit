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
LevelSetObject::LevelSetObject( int id) : m_id(id){
};

/*!
 * Get the id 
 * @return id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
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

    // We don't need to propagate the sign in the narrowband
    //
    // An item is in the narrow band if it has a levelset value that differs
    // from the defualt value.
    std::unordered_set<long> alreadyEvaluated;
    VolumeKernel const &mesh = *(visitee->getMesh()) ;

    PiercedIterator<LevelSetInfo> infoItr = m_ls.begin() ;
    while (infoItr != m_ls.end()) {
        double &value = (*infoItr).value;
        if(!utils::DoubleFloatingEqual()(std::abs(value), levelSetDefaults::VALUE)) {
            alreadyEvaluated.insert(infoItr.getId()) ;
        }
        ++infoItr;
    }

    // If the all cells have the correct value we don't need to progagate the
    // sign
    if (alreadyEvaluated.size() == (size_t) mesh.getCellCount()) {
        return;
    }

    // Define the seed candidates
    //
    // First list cells in the narroband, then all other cells. The cells
    // outisde the narrowband will be used as seeds only if there are regions
    // of the mesh disconnected from the narrow band.
    std::vector<long> seedCandidates;
    seedCandidates.reserve(mesh.getCellCount());

    seedCandidates.assign(alreadyEvaluated.begin(), alreadyEvaluated.end());
    for (const Cell &cell : mesh.getCells()) {
        long cellId = cell.getId();
        if (alreadyEvaluated.count(cellId) == 0) {
            seedCandidates.push_back(cellId);
        }
    }

    // Identify real seeds and propagate the sign
    for (long seed : seedCandidates) {
        // Get the neighbours that still need to be processed
        //
        // If a cell is surrounded only by items already evaluated,
        // this cell can not be uses as a seed.
        std::stack<long> processList;

        int nSeedNeighs;
        const long* seedNeighs;
        std::vector<long> seedNeighList;
        if (mesh.getCells().size() == 0) {
            seedNeighList = mesh.findCellFaceNeighs( seed ) ;
            seedNeighs    = seedNeighList.data() ;
            nSeedNeighs   = seedNeighList.size() ;
        } else {
            const Cell& cell  = mesh.getCell( seed );
            seedNeighs  = cell.getAdjacencies() ;
            nSeedNeighs = cell.getAdjacencyCount() ;
        }

        for ( int n=0; n<nSeedNeighs; ++n ) {
            long neigh = seedNeighs[n] ;
            if(neigh<0){
                continue ;
            }

            if (alreadyEvaluated.count(neigh) == 0) {
                processList.push(neigh);
            }
        }

        if (processList.empty()) {
            continue;
        }

        // Discard seeds with a LS value equal to 0
        //
        // A cell can have a levelset value equal to zero only if it's inside
        // the narrow band, therefore the levelset value returned by getLS
        // is enough to make this check.
        double ls = getLS(seed);
        if( utils::DoubleFloatingEqual()(std::abs(ls), (double) 0.) ) {
            continue;
        }

        // Get the sign of the seed
        short seedSign = ls > 0 ? 1 : -1;

        // Propagate the sign
        while (!processList.empty()) {
            long id = processList.top();
            processList.pop();

            // Get the value associated to the id
            //
            // A new value needs to be created only if the sign to propagate
            // is different from the default sign.
            infoItr = m_ls.find(id) ;
            if( infoItr == m_ls.end() && seedSign != levelSetDefaults::SIGN ){
                infoItr = m_ls.emplace(id) ;
            }

            // Update the value
            if( infoItr != m_ls.end() ) {
                (*infoItr).value = seedSign * levelSetDefaults::VALUE;
            }

            // Add non-evaluated neighs to the process list
            int nNeighs;
            const long* neighs;
            std::vector<long> neighList;
            if (mesh.getCells().size() == 0) {
                neighList = mesh.findCellFaceNeighs( id ) ;
                neighs    = neighList.data() ;
                nNeighs   = neighList.size() ;
            } else {
                const Cell& cell  = mesh.getCell( id );
                neighs  = cell.getAdjacencies() ;
                nNeighs = cell.getAdjacencyCount() ;
            }

            for ( int n=0; n<nNeighs; ++n ) {
                long neigh = neighs[n];
                if(neigh<0){
                    continue;
                }

                if (alreadyEvaluated.count(neigh) == 0) {
                    processList.push(neigh);
                }
            }

            // The item has been processeed
            //
            // If all cells have been evaluated we can stop the propagation
            alreadyEvaluated.insert(id);
        }

        // If all cells have been evaluated we can stop the propagation
        if (alreadyEvaluated.size() == (size_t) mesh.getCellCount()) {
            break;
        }
    }
};


/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::ostream &stream ){

    bitpit::PiercedVector<LevelSetInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    IO::binary::write(stream, m_id) ;
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
