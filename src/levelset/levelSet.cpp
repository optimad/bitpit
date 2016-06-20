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

# include <unordered_set>

# include "bitpit_SA.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
 * @ingroup levelset
 * @class  LevelSetInfo
 *
 * @brief  A public container which includes all information provided by LevelSet
 *
 * LevelSetInfo conatins the following information
 * - distance to closest object
 * - gradient of level set function
 * - the id of the closest object
 * - the patch of closest object which contains the projection point
 * - the id of the surface element closest to the projection point 
 *
 * If a grid point lies within the narrow band of an object, all of these information are available, provided that the clsest object may provide them.
 * On the contrary, if the grid point lies were the levelset has been propagated not all of these information may be available; 
 * in particular the object, patch and segment may be set to the default values.
 *
*/

/*!
 * Default constructor
 */
LevelSetInfo::LevelSetInfo() {
    value = levelSetDefaults::VALUE ;
    gradient = levelSetDefaults::GRADIENT ;
    object = levelSetDefaults::OBJECT ;
    part = levelSetDefaults::PART  ;
    support = levelSetDefaults::SUPPORT ;
}

/*!
 * @ingroup levelset
 * @class  LevelSet
 *
 * @brief  Level Set driver class
 *
 * LevelSet is the main user interface class for computing signed- or unsigned- distance functions on Cartesian or Octree meshes
 * with respect to geometrical objects. The user needs to define the computional kernel by calling setMesh() and the objects which define the zero 
 * levelset via addObject().
 *
 * LevelSet will calculate the exact distance with respect the objects within a narrow band.
 * Outside this narrow band an approximate value will be calculated.
 *
 * The user may set the size of the narrow band explicitly.
 * Alternatively LevelSet will guarantee a least on cell center with exact levelset values across the zero-levelset iso-surface.
 *
 * LevelSet will test if the underlying mesh can provide a MPI communicator.
 * In case LevelSet is parallelized according the underlying mesh partitioning.
*/

/*!
 * Default constructor
 */
LevelSet::LevelSet() {

    m_kernel = NULL ;
    m_object.clear() ;

    m_userRSearch = false ;

    m_signedDF    = true ;
    m_propagateS  = false;
    m_propagateV  = false;

};

/*!
 * Destructor of LevelSet
*/
LevelSet::~LevelSet(){
    if( m_kernel != NULL){
        delete m_kernel ;
    }

    for( auto &object : m_object){
        delete object.second  ;
    }

    m_object.clear() ;

};

/*!
 * Sets the grid on which the levelset function should be computed.
 * Only cartesian and octree patches are supported at this moment.
 * @param[in] mesh computational grid
 */
void LevelSet::setMesh( VolumeKernel* mesh ) {

    if( VolCartesian* cartesian = dynamic_cast<VolCartesian*> (mesh) ){
        m_kernel = new LevelSetCartesian( *cartesian) ;

    } else if( VolOctree* octree = dynamic_cast<VolOctree*> (mesh) ){
        m_kernel = new LevelSetOctree(*octree) ;
    
    } else{
        log::cout() << "Mesh non supported in LevelSet::setMesh()" << std::endl ;
    }; 

    return;
};

/*!
 * Sets the grid on which the levelset function should be computed.
 * @param[in] cartesian cartesian patch
 */
void LevelSet::setMesh( VolCartesian* cartesian ) {

    m_kernel = new LevelSetCartesian( *cartesian) ;

    return;
};

/*!
 * Sets the grid on which the levelset function should be computed.
 * @param[in] octree octree patch
 */
void LevelSet::setMesh( VolOctree* octree ) {

    m_kernel = new LevelSetOctree( *octree) ;

    return;
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfUnstructured> &&segmentation, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, std::move(segmentation)) ;

    return addObject(lsSeg);
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfUnstructured *segmentation, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, segmentation) ;

    return addObject(lsSeg);
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfaceKernel> &&segmentation, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id) ;
    if( SurfUnstructured* surfUnstructured = dynamic_cast<SurfUnstructured*>(segmentation.get()) ){
        segmentation.release();
        std::unique_ptr<SurfUnstructured> surfUnstructuredUPtr = std::unique_ptr<SurfUnstructured>(surfUnstructured) ;
        lsSeg->setSegmentation( std::move(surfUnstructuredUPtr) );
    } else {
        throw std::runtime_error ("Segmentation type not supported");
    }

    return addObject(lsSeg);
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfaceKernel *segmentation, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id) ;
    if( SurfUnstructured* surfUnstructured = dynamic_cast<SurfUnstructured*>(segmentation)){
        lsSeg->setSegmentation( surfUnstructured );
    } else {
        throw std::runtime_error ("Segmentation type not supported");
    }

    return addObject(lsSeg);
};

/*!
 * Clear LevelSet entirely
 */
void LevelSet::clear(){

    if( m_kernel != NULL){
        delete m_kernel ;
    }

    m_object.clear() ;
    return ;
};

/*!
 * Clear all LevelSet Objects but leave kernel unaltered
 */
void LevelSet::clearObject(){
    m_object.clear() ;
    return ;
};

/*!
 * Clear only one LevelSet Object
 * @param[in] id Id of object to be deleted
 */
void LevelSet::clearObject(int id){

    m_object.erase(id) ;
    return ;
};


/*!
 * Adds a generic LevelSetObject
 * @param[in] object generic object
 */
int LevelSet::addObject( LevelSetObject* object ) {


    if( LevelSetSegmentation* segmentation = dynamic_cast<LevelSetSegmentation*>(object) ){
        m_object.insert( {{ object->getId(), new LevelSetSegmentation( *segmentation) }} ) ;
        return object->getId();

    } else {
        return levelSetDefaults::OBJECT;

    };

};

/*!
 * Get all levelset information 
 * @param[in] i index of cell
 * @return levelset information
 */
LevelSetInfo LevelSet::getLevelSetInfo( const long &i)const {
    return( m_kernel->getLevelSetInfo(i) ) ;
};

/*!
 * Get the levelset value of the i-th local cell.
 * @param[in] i index of cell
 * @return value of the i-th cell
 */
double LevelSet::getLS( const long &i)const {
    return( m_kernel->getLS(i) ) ;
};

/*!
 * Get the levelset gradient of the i-th local cell.
 * @param[in] i index of cell
 * @return Array with components of the gradient 
 */
std::array<double,3> LevelSet::getGradient(const long &i) const {
    return( m_kernel->getGradient(i) ) ;
};

/*!
 * Get the id of closest object
 * @param[in] i index of cell
 * @return id of closest object
 */
int LevelSet::getObject(const long &i) const {
    return( m_kernel->getObject(i) ) ;
};

/*!
 * Get the object and part id of the projection point
 * @param[in] i index of cell
 * @return pair containing object and part id
 */
std::pair<int,int> LevelSet::getPart(const long &i) const {

    return (m_kernel->getPart(i)) ;
};

/*!
 * Get the object and support id of the projection point
 * @param[in] i index of cell
 * @return pair containing object and support id
 */
std::pair<int,long> LevelSet::getSupport(const long &i) const {

    return (m_kernel->getSupport(i)) ;
};

/*!
 * Get the sign of the levelset function
 * @param[in] i index of cell
 * @return sign
 */
short LevelSet::getSign(const long &i)const{
    return( m_kernel->getSign(i) ) ;
};

/*!
 * Returns if the cell centroid lies within the narrow band.
 * @param[in] i index of cell
 * @return true/false if the cell centroid is within the narrow band
 */
bool LevelSet::isInNarrowBand(const long &i)const{
    return( m_kernel->isInNarrowBand(i) ) ;
};

/*!
 * Get the current size of the narrow band.
 * @return Physical size of the current narrow band.
 */
double LevelSet::getSizeNarrowBand()const{
    return( m_kernel->getSizeNarrowBand() ) ;
};

/*!
 * Set if the signed or unsigned levelset function should be computed.
 * @param[in] flag true/false for signed /unsigned Level-Set function .
 */
void LevelSet::setSign(bool flag){
    m_signedDF = flag;

};

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation .
 */
void LevelSet::setPropagateSign(bool flag){
    m_propagateS = flag;
};

/*!
 * Set if the levelset value has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation.
 */
void LevelSet::setPropagateValue(bool flag){
    m_propagateV = flag;
};

/*!
 * Manually set the physical size of the narrow band.
 * @param[in] r Size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double r){
    m_userRSearch = true ;
    m_kernel->setSizeNarrowBand(r) ;
};

/*!
 * Computes levelset on given mesh with respect to the objects.
 * This routines needs to be called at least once.
 */
void LevelSet::compute(){

    double RSearch ;

    if( !m_userRSearch){
        RSearch = 0. ;
        for( const auto &visitor : m_object ){
            RSearch = std::max( RSearch, m_kernel->computeSizeNarrowBand( visitor.second ) ) ;
        }

        m_kernel->setSizeNarrowBand(RSearch) ;

    }

    for( const auto &visitor : m_object ){
        visitor.second->computeLSInNarrowBand( m_kernel, RSearch, m_signedDF) ; 
    }


    if( m_propagateS ) m_kernel->propagateSign( m_object ) ;
//    if( propagateV ) propagateValue( ) ;

    return ;
}

/*!
 * Updates the levelset after mesh adaptation.
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::update( const std::vector<adaption::Info> &mapper ){

    // Check the mapper to detect the operations to perform
    bool updateNarrowBand = false;
    std::unordered_map<int,std::vector<long>> sendList, recvList ;

    for( const auto &event : mapper){
        if( event.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        if( event.type == adaption::Type::TYPE_PARTITION_SEND){
            sendList.insert({{event.rank,event.previous}}) ;
        } else if( event.type == adaption::Type::TYPE_PARTITION_RECV){
            recvList.insert({{event.rank,event.current}}) ;
        } else if( event.type != adaption::Type::TYPE_PARTITION_NOTICE){
            updateNarrowBand = true ;
        }
    }

    // Evaluate new narrow band size
    double newRSearch ;
    if (updateNarrowBand) {
        if(m_userRSearch){
            newRSearch = m_kernel->getSizeNarrowBand() ;
        } else {
            newRSearch = m_kernel->updateSizeNarrowBand( mapper, m_object )  ;
        };
    }

    // Clear levelset
    m_kernel->clearAfterMeshMovement(mapper) ;

    // Update narrow band
    if (updateNarrowBand) {
        m_kernel->filterOutsideNarrowBand(newRSearch) ;

        for( const auto &visitor : m_object ){
            visitor.second->updateLSInNarrowBand( m_kernel, mapper, newRSearch, m_signedDF ) ;
        }
    }

#if BITPIT_ENABLE_MPI
    // Parallel communications
    if (sendList.size() > 0 || recvList.size() > 0) {
        communicate(sendList, recvList, &mapper ) ;
    }

    exchangeGhosts() ;
#endif

    // Finish narrow band update
    if (updateNarrowBand) {
        newRSearch = m_kernel->computeSizeNarrowBandFromLS();

        m_kernel->filterOutsideNarrowBand(newRSearch) ;
        m_kernel->setSizeNarrowBand(newRSearch) ;

        for( const auto &visitor : m_object ){
            visitor.second->filterOutsideNarrowBand(m_kernel) ;
        }

        if( m_propagateS ) m_kernel->propagateSign( m_object ) ;
//TODO      if( propagateV ) updatePropagatedValue() ;
    }

    return;

};

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::fstream &stream ){

    bitpit::genericIO::flushBINARY(stream, m_userRSearch);
    bitpit::genericIO::flushBINARY(stream, m_signedDF);
    bitpit::genericIO::flushBINARY(stream, m_propagateS);
    bitpit::genericIO::flushBINARY(stream, m_propagateV);

    m_kernel->dump( stream ) ;

    for( const auto &visitor : m_object ){
        visitor.second->dump( stream ) ;
    }


    return ;
};

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::restore( std::fstream &stream ){

    bitpit::genericIO::absorbBINARY(stream, m_userRSearch);
    bitpit::genericIO::absorbBINARY(stream, m_signedDF);
    bitpit::genericIO::absorbBINARY(stream, m_propagateS);
    bitpit::genericIO::absorbBINARY(stream, m_propagateV);

    m_kernel->restore( stream ) ;

    for( const auto &visitor : m_object ){
        visitor.second->restore( stream ) ;
    }


    return ;
}

#if BITPIT_ENABLE_MPI
/*!
 * Checks if MPI communicator is available in underlying mesh.
 * If available, MPI communicator is retreived from mesh (and duplicated if necessary) and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSet::assureMPI( ){

    return(m_kernel->assureMPI() ) ;
}

/*!
 * Update of ghost cell;
 */
void LevelSet::exchangeGhosts(  ){

    std::unordered_map<int,std::vector<long>> &sendList =  m_kernel->getMesh()->getGhostExchangeSources() ;
    std::unordered_map<int,std::vector<long>> &recvList =  m_kernel->getMesh()->getGhostExchangeTargets() ;

    communicate(sendList,recvList);

    return ;
}

/*!
 * communicates data structures of kernel and objects.
 * If mapper!=NULL, clearAfterMeshMovement of kernel and objects will be called between send and receive.
 * @param[in] sendList list of elements to be send
 * @param[in] recvList list of elements to be received
 * @param[in] mapper mapper containing mesh modifications
 */
void LevelSet::communicate( std::unordered_map<int,std::vector<long>> &sendList, std::unordered_map<int,std::vector<long>> &recvList, std::vector<adaption::Info> const *mapper){

    if( assureMPI() ){

        MPI_Comm meshComm = m_kernel->getCommunicator() ;

        int                 nClasses = 1 + m_object.size() ;
        DataCommunicator    sizeCommunicator(meshComm) ; 
        DataCommunicator    dataCommunicator(meshComm) ; 

        sizeCommunicator.setTag(0) ;
        dataCommunicator.setTag(1) ;

        // start receive of sizes
	    for (const auto entry : sendList) {
            int rank = entry.first;

            sizeCommunicator.setRecv( rank, (1 + nClasses) * sizeof(long) ) ;
            sizeCommunicator.startRecv(rank);

        }

        { // send first number of items and then data itself

	        // Fill the buffer with the given field and start sending the data
	        for (const auto entry : sendList) {

	        	// Get the send buffer
	        	int rank = entry.first;

                sizeCommunicator.setSend(rank, (1 + nClasses) * sizeof(long) ) ;
                dataCommunicator.setSend(rank, 0 ) ;

                SendBuffer &sizeBuffer = sizeCommunicator.getSendBuffer(rank);
                SendBuffer &dataBuffer = dataCommunicator.getSendBuffer(rank);

                m_kernel->writeCommunicationBuffer( entry.second, sizeBuffer, dataBuffer ) ;
                for( const auto &visitor : m_object){
                    visitor.second->writeCommunicationBuffer( entry.second, sizeBuffer, dataBuffer ) ;
                }

                sizeBuffer << dataBuffer.capacity() ;

                // Start the send
                sizeCommunicator.startSend(rank);
                dataCommunicator.startSend(rank);

	        }

        }


        { 
            // as soon as sizes are received start receiving data.
            int nCompletedRecvs(0);

            std::vector<long>   dummy(nClasses) ;
            std::unordered_map<int,std::vector<long> > items ;

            std::vector<long>::iterator itemItr ;

            // receive sizes
            while (nCompletedRecvs < sizeCommunicator.getRecvCount()) {
                int rank = sizeCommunicator.waitAnyRecv();
                RecvBuffer &sizeBuffer = sizeCommunicator.getRecvBuffer(rank);

                items.insert({{rank, dummy}}) ;
                itemItr = items[rank].begin() ;

                sizeBuffer >> *(itemItr) ;

                ++itemItr; 
                for(size_t i = 0; i < m_object.size(); ++i){
                    sizeBuffer >> *(itemItr) ;
                    ++itemItr ;
                };

                long dataSize ;
                sizeBuffer >> dataSize ;

                dataCommunicator.setRecv(rank,dataSize) ;
                dataCommunicator.startRecv(rank) ;

                ++nCompletedRecvs;
            }

            if( mapper!=NULL){
                m_kernel->clearAfterMeshMovement(*mapper);
                for( const auto & visitor : m_object){
                    visitor.second->clearAfterMeshMovement(*mapper);
                }
            }

            //  post-process data from buffer to data within classes
            nCompletedRecvs = 0;
            while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
                int rank = dataCommunicator.waitAnyRecv();

                RecvBuffer &dataBuffer = dataCommunicator.getRecvBuffer(rank);
                itemItr = items[rank].begin();

                m_kernel->readCommunicationBuffer( recvList[rank], *itemItr, dataBuffer ) ;
                ++itemItr ;

                for( const auto &visitor : m_object){
                    visitor.second->readCommunicationBuffer( recvList[rank], *itemItr, dataBuffer ) ;
                    ++itemItr ;
                }

                ++nCompletedRecvs;
            }
        }

        sizeCommunicator.waitAllSends();
        sizeCommunicator.finalize() ;

        dataCommunicator.waitAllSends();
        dataCommunicator.finalize() ;

    }


    return ;
}


#endif 

}

