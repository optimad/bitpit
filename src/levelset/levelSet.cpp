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
# include "communications.hpp"
# endif

# include <unordered_set>

# include "bitpit_SA.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
 * @ingroup levelset
 * @class  LevelSet
 *
 * @brief  Level Set driver class
 *
 * LevelSet is the main user interface class for computing signed- or unsigned- distance functions on Cartesian or Octree meshes
 * with respect to geometrical objects. The user needs to define the computional by calling setMesh() and the objects which define the zero 
 * levelset via addObject().
 *
 * LevelSet will calculate the exact distance with respect within a narrow band.
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
 * Only caretsian and octree patches are supported at this moment.
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
 */
int LevelSet::addObject( SurfUnstructured* segmentation, int id ) {

    if( id == levelSetDefaults::DEFAULT_ID ){
        id = m_object.size() ;
    }

    m_object.insert( {{id, new LevelSetSegmentation(id,segmentation)}} ) ;

    return id;
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 */
int LevelSet::addObject( SurfaceKernel* segmentation, int id ) {

    if( SurfUnstructured* unstruct = dynamic_cast<SurfUnstructured*>(segmentation) ){
        return (addObject(unstruct, id) ) ;

    } else {
        return levelSetDefaults::NULL_ADD;
    };

};

/*!
 * Clear LevelSet entirely
 */
void LevelSet::clear(){

    if( m_kernel != NULL){
        delete m_kernel ;
    }

    std::cout << "cazz" << std::endl ;

    m_object.clear() ;
    return ;
};

/*!
 * Clear all LevelSet Objects
 */
void LevelSet::clearObject(){
    m_object.clear() ;
    return ;
};

/*!
 * Clear LevelSet Object
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
        return levelSetDefaults::NULL_ADD;

    };

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
 * @return Array with components of the Sdf gradient of the i-th local element of the octree mesh.
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
 * Get the id of closest object
 * @param[in] i index of cell
 * @return id of closest object
 */
long LevelSet::getSupport(const long &i) const {

    int id= getObject(i);

    if( id != levelSetDefaults::OBJECT){

        auto objItr = m_object.find(id) ;

        if( objItr!=m_object.end() ){
            return objItr->second->getSupport(i);
        };

    };
    return levelSetDefaults::ELEMENT ;
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


    if( m_propagateS ) m_kernel->propagateSign( m_object[0] ) ; //TODO for several objects
//    if( propagateV ) propagateValue( ) ;

    return ;
}

/*!
 * Updates the levelset after mesh adaptation.
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::update( const std::vector<adaption::Info> &mapper ){

    double  newRSearch ;

    if(m_userRSearch){
        newRSearch = m_kernel->getSizeNarrowBand() ;
    } else {
        newRSearch = m_kernel->updateSizeNarrowBand( mapper )  ;
    };

    m_kernel->clearAfterAdaption(mapper,newRSearch) ;

    for( const auto &visitor : m_object ){
        visitor.second->updateLSInNarrowBand( m_kernel, mapper, newRSearch, m_signedDF ) ;
    }


    if( m_propagateS ) m_kernel->propagateSign( m_object[0] ) ; //TODO several objects
//TODO    if( propagateV ) updatePropagatedValue() ;

    m_kernel->setSizeNarrowBand(newRSearch) ;

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
 * If available MPI communicator is retreived from mesh and duplicated if necessary and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSet::assureMPI( ){

    return(m_kernel->assureMPI() ) ;
}

/*!
 * Frees the MPI communicator.
 */
void LevelSet::finalizeMPI( ){

    m_kernel->finalizeMPI() ;

    for( auto visitor:m_object){
        visitor.second->finalizeMPI() ;
    }

}


/*!
 * Distribution of levelset over available processes after partitioning of mesh
 * @param[in] mapper mapper describing partitioning
 */
void LevelSet::loadBalance( const std::vector<adaption::Info> &mapper ){

    if( assureMPI() ){

        MPI_Comm meshComm = m_kernel->getCommunicator() ;

        int                 nClasses = 1 + m_object.size() ;
        DataCommunicator    sizeCommunicator(meshComm) ; 
        DataCommunicator    dataCommunicator(meshComm) ; 


        // start receive of sizes
        for( const auto &event : mapper){
            if( event.entity == adaption::Entity::ENTITY_CELL){
                if( event.type == adaption::Type::TYPE_PARTITION_RECV){
                    short rank = event.rank;

                    sizeCommunicator.setRecv( rank, 2*sizeof(long)*nClasses ) ;
                    sizeCommunicator.startRecv(rank);

                }
            }
        }

        { // send first number of items and then data itself

            for( const auto &event : mapper){
                if( event.entity == adaption::Entity::ENTITY_CELL){

                    if( event.type == adaption::Type::TYPE_PARTITION_SEND){
                        short rank = event.rank;

                        //determine elements to send
                        sizeCommunicator.setSend(rank, 2*sizeof(long) *nClasses ) ;
                        dataCommunicator.setSend(rank, 0 ) ;

                        //store data in buffer
                        OBinaryStream &sizeBuffer = sizeCommunicator.getSendBuffer(rank);
                        OBinaryStream &dataBuffer = dataCommunicator.getSendBuffer(rank);

                        m_kernel->writeCommunicationBuffer( event.previous, sizeBuffer, dataBuffer ) ;
                        for( const auto &visitor : m_object){
                            visitor.second->writeCommunicationBuffer( event.previous, sizeBuffer, dataBuffer ) ;
                        }

                        // Start the send
                        sizeCommunicator.startSend(rank);
                        dataCommunicator.startSend(rank);
                    }
                }
            }
        }


        { 
            // as soon as sizes are received start receiving data.
            int rank, nCompletedRecvs;
            long data, dataSize ;
            std::vector<long> items(nClasses) ;
            std::vector<long>::iterator  itemItr = items.begin() ;

            // receive sizes
            nCompletedRecvs = 0;
            while (nCompletedRecvs < sizeCommunicator.getRecvCount()) {
                rank = sizeCommunicator.waitAnyRecv();
                IBinaryStream &sizeBuffer = sizeCommunicator.getRecvBuffer(rank);

                sizeBuffer >> *(itemItr) ;
                sizeBuffer >> dataSize ;

                ++itemItr; 
                for( const auto & visitor : m_object){
                    sizeBuffer >> *(itemItr) ;
                    sizeBuffer >> data ;
                    dataSize += data ;
                    ++itemItr ;
                };

                dataCommunicator.setRecv(rank,dataSize) ;
                dataCommunicator.startRecv(rank) ;

                ++nCompletedRecvs;
            }

            //  post-process data from buffer to data within classes
            nCompletedRecvs = 0;
            while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
                rank = dataCommunicator.waitAnyRecv();

                IBinaryStream &dataBuffer = dataCommunicator.getRecvBuffer(rank);
                itemItr = items.begin();

                m_kernel->readCommunicationBuffer( *itemItr, dataBuffer ) ;
                ++itemItr ;

                for( const auto &visitor : m_object){
                    visitor.second->readCommunicationBuffer( *itemItr, dataBuffer ) ;
                    ++itemItr ;
                }

                ++nCompletedRecvs;
            }
        }

        sizeCommunicator.finalize() ;
        dataCommunicator.finalize() ;

    }

    return ;
}

#endif 

}

