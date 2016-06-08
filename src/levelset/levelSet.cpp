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

/*!
 *    \date            10/jul/2014
 *    \authors        Alessandro Alaia
 *    \authors        Haysam Telib
 *    \authors        Edoardo Lombardi
 *    \version        0.1
 *    \copyright        Copyright 2015 Optimad engineering srl. All rights reserved.
 *    \par            License:\n
 *    This version of Class_LevelSet_Stl is released under the LGPL License.
 *
 *    \brief Level Set Manager Class - 3D PABLO Octree specialization
 *
 *    Level Set Stl is a user interface class. One user should (read can...) work only
 *    with this Class and its methods to maintain a signed distance function (Sdf)
 *    computed from a piece-wise linear approximation of a d manifold in a 3D Euclidean
 *    space. Sdf is computed in a narrow band of at least 2 mesh cell centers
 *    around the geometry.
 *    Parallel implementation developed by using the features of PABLO library.
 *
 */

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
    @ingroup levelset
    @class  LevelSetKernel
    @brief  Level Set on Stl Manager Class

    LevelSetKernel is a user interface class. One user should (read can...) work only
    with this class and its methods to maintain a signed distance function (Sdf)
    computed from a piece-wise linear approximation of a d manifold in a 3D Euclidean
    space. Sdf is computed in a narrow band of at least 2 mesh cell centers
    around the geometry.
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

#if BITPIT_ENABLE_MPI
    m_commMPI = MPI_COMM_NULL;
# endif
};

/*!
 * Destructor of LevelSet
*/
LevelSet::~LevelSet(){
    delete m_kernel ;

    for( auto &object : m_object){
        delete object  ;
    }

};

/*!
 */
void LevelSet::setMesh( VolumeKernel* mesh ) {

    if( VolCartesian* cartesian = dynamic_cast<VolCartesian*> (mesh) ){
        m_kernel = new LevelSetCartesian( *cartesian) ;

    } else if( VolOctree* octree = dynamic_cast<VolOctree*> (mesh) ){
        m_kernel = new LevelSetOctree(*octree) ;
    }; 

    return;
};

void LevelSet::setMesh( VolCartesian* cartesian ) {

    m_kernel = new LevelSetCartesian( *cartesian) ;

    return;
};

void LevelSet::setMesh( VolOctree* octree ) {

    m_kernel = new LevelSetOctree( *octree) ;

    return;
};

void LevelSet::addObject( SurfUnstructured* segmentation ) {

    int id = m_object.size() ;
    m_object.push_back( new LevelSetSegmentation(id,segmentation) ) ;

    return;
};

void LevelSet::addObject( LevelSetObject* object ) {


    if( LevelSetSegmentation* segmentation = dynamic_cast<LevelSetSegmentation*>(object) ){
        m_object.push_back( new LevelSetSegmentation(*segmentation) ) ;
    }

    return;
};

/*!
 * Get the Sdf value of the i-th local element of the octree mesh.
 * @param[in] i Local index of target octant.
 * @return Value of the i-th local element of the octree mesh.
 */
double LevelSet::getLS( const long &i)const {
    return( m_kernel->getLS(i) ) ;
};

/*!
 * Get the Sdf gradient vector of the i-th local element of the octree mesh.
 * @param[in] i Local index of target octant.
 * @return Array with components of the Sdf gradient of the i-th local element of the octree mesh.
 */
std::array<double,3> LevelSet::getGradient(const long &i) const {
    return( m_kernel->getGradient(i) ) ;
};

/*!
 * Get the sign of the levelset function
 * @param[in] i Local index of target octant.
 * @return sign
 */
short LevelSet::getSign(const long &i)const{
    return( m_kernel->getSign(i) ) ;
};

/*!
 * Get if the Sdf value of the i-th local element is exactly computed or not.
 * @param[in] i Local index of target octant.
 * @return True/false if the Sdf value is exactly computed (true) or not (false).
 */
bool LevelSet::isInNarrowBand(const long &i)const{
    return( m_kernel->isInNarrowBand(i) ) ;
};

/*!
 * Get the current size of the narrow band.
 * @return Physical size of the current narrow band to guarantee at least one element inside it.
 */
double LevelSet::getSizeNarrowBand()const{
    return( m_kernel->getSizeNarrowBand() ) ;
};

/*!
 * Set if the signed or unsigned LevelSet should be computed.
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

void LevelSet::compute(){

    double RSearch ;

    if( !m_userRSearch){
        RSearch = 0. ;
        for( const auto &visitor : m_object ){
            RSearch = std::max( RSearch, m_kernel->computeSizeNarrowBand( visitor ) ) ;
        }

        m_kernel->setSizeNarrowBand(RSearch) ;

    }

    for( const auto &visitor : m_object ){
        visitor->computeLSInNarrowBand( m_kernel, RSearch, m_signedDF) ; 
    }


    if( m_propagateS ) m_kernel->propagateSign( m_object[0] ) ; //TODO for several objects
//    if( propagateV ) propagateValue( ) ;

    return ;
}

/*!
 * Compute the levelset function 
 */
void LevelSet::update( const std::vector<adaption::Info> &mapper ){

    double  newRSearch ;

    if(m_userRSearch){
        newRSearch = m_kernel->getSizeNarrowBand() ;
    } else {
        newRSearch = std::max( newRSearch, m_kernel->updateSizeNarrowBand( mapper ) ) ;
    };

    m_kernel->clearAfterAdaption(mapper,newRSearch) ;

    for( const auto &visitor : m_object ){
        visitor->updateLSInNarrowBand( m_kernel, mapper, newRSearch, m_signedDF ) ;
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
        visitor->dump( stream ) ;
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
        visitor->restore( stream ) ;
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

    if( m_commMPI == MPI_COMM_NULL){

        MPI_Comm meshComm = m_kernel->getMesh()->getCommunicator() ;

        if( meshComm == MPI_COMM_NULL){
            return false;
        } else {
            MPI_Comm_dup(m_kernel->getMesh()->getCommunicator(), &m_commMPI);
            return true; 
        }
    } else {
        return true ;
    };

}

/*!
 * Frees the MPI communicator.
 */
void LevelSet::finalizeMPI( ){

    if( m_commMPI != MPI_COMM_NULL){
        MPI_Comm_free( &m_commMPI ) ;
    }
}


/*!
 * Repartioning of levelset after partitioning of mesh
 * @param[in] mapper mapper describing partitioning
 */
void LevelSet::loadBalance( const std::vector<adaption::Info> &mapper ){

    if( assureMPI() ){

        int                 nClasses = 1 + m_object.size() ;
        DataCommunicator    sizeCommunicator(m_commMPI) ; 
        DataCommunicator    dataCommunicator(m_commMPI) ; 


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
                        long nItems = event.previous.size() ;

                        sizeCommunicator.setSend(rank, 2*sizeof(long) *nClasses ) ;
                        dataCommunicator.setSend(rank, 0 ) ;

                        //store data in buffer
                        OBinaryStream &sizeBuffer = sizeCommunicator.getSendBuffer(rank);
                        OBinaryStream &dataBuffer = dataCommunicator.getSendBuffer(rank);

                        m_kernel->writeCommunicationBuffer( event.previous, sizeBuffer, dataBuffer ) ;
                        for( const auto &visitor : m_object){
                            visitor->writeCommunicationBuffer( event.previous, sizeBuffer, dataBuffer ) ;
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
            long nItems, data, dataSize ;
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
                    visitor->readCommunicationBuffer( *itemItr, dataBuffer ) ;
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

