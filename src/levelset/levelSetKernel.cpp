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

# include <stack>
# include <unordered_set>

# include "bitpit_SA.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
    @ingroup levelset
    @class  LevelSetKernel
    @brief  Mesh specific implementation to calculate the levelset function

*/

/*!
 * Default constructor.
 */
LevelSetKernel::LevelSetKernel() {

    m_mesh = NULL ;

#if BITPIT_ENABLE_MPI
    m_commMPI = MPI_COMM_NULL;
# endif

};

/*!
 * Constructor
 * @param[in] patch underlying mesh
 */
LevelSetKernel::LevelSetKernel( VolumeKernel *patch): LevelSetKernel() {
    m_mesh = patch ;
};


/*!
 * Destructor of LevelSetKernel
*/
LevelSetKernel::~LevelSetKernel(){
    m_mesh = NULL ;

# if BITPIT_ENABLE_MPI
    freeCommunicator();
# endif

};

/*!
 * Returns reference to LSInfo
*/
PiercedVector<LevelSetKernel::LSInfo>& LevelSetKernel::getLSInfo(){
    return m_ls ;
} 

/*!
 * Returns pointer to underlying mesh.
 * @return pointer to mesh
*/
VolumeKernel* LevelSetKernel::getMesh() const{
    return m_mesh ;
} 


/*!
 * Get the Sdf value of the i-th local element of the octree mesh.
 * @param[in] i cell index
 * @return levelset value in cell
 */
double LevelSetKernel::getLS( const long &i)const {

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
std::array<double,3> LevelSetKernel::getGradient(const long &i) const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::GRADIENT;
    } else {
        return (  m_ls[i].gradient );
    };

};

/*!
 * Get the id of closest object
 * @param[in] i cell index
 * @return id of closest object
 */
int LevelSetKernel::getObject(const long &i) const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::OBJECT;
    } else {
        return (  m_ls[i].object );
    };

};

/*!
 * Get the sign of the levelset function
 * @param[in] i cell index
 * @return sign of levelset
 */
short LevelSetKernel::getSign(const long &i)const{

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
bool LevelSetKernel::isInNarrowBand(const long &i)const{

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
double LevelSetKernel::getSizeNarrowBand()const{
    return m_RSearch;
};

/*!
 * Manually set the size of the narrow band.
 * @param[in] r size of the narrow band.
 */
void LevelSetKernel::setSizeNarrowBand(double r){
    m_RSearch = r;
};

/*!
 * Computes the Level Set Gradient on on cell by first order finite-volume upwind stencil
 * @param[in] I index of cell 
 * @return gradient
 */
std::array<double,3> LevelSetKernel::computeGradientUpwind( const long &I ){

    int                     i, border ;
    std::array<double,3>    normal, gradient = {{0.,0.,0.}};
    long                    owner;

    Cell            &cell = m_mesh->getCell(I) ;
    int                     nI = cell.getInterfaceCount() ;
    const long*             interfaces = cell.getInterfaces() ;
    const long*             neighbours = cell.getAdjacencies() ;

    long                    F, N ;
    double                  value, area;


    for( i=0; i<nI; ++i){
        F = interfaces[i] ;
        N = neighbours[i] ;

        Interface       &interface = m_mesh->getInterface(F) ;

        owner  = interface.getOwner() ;
        border = interface.isBorder() ;

        area   = m_mesh->evalInterfaceArea(F) ;
        normal = m_mesh->evalInterfaceNormal(F) ;

        if( owner != I)
            normal = -1. *normal ;

        value = border* m_ls[I].value + (1-border) *std::min( m_ls[I].value , m_ls[N].value ) ;

        gradient += area *value *normal ;

    };

    gradient = gradient /m_mesh->evalCellVolume(I) ;


    return gradient;

};

/*!
 * Computes the Level Set Gradient on on cell by second order finite-volume central stencil
 * @param[in] I index of cell 
 * @return gradient
 */
std::array<double,3> LevelSetKernel::computeGradientCentral( const long &I ){

    int                     i, border ;
    std::array<double,3>    normal, gradient = {{0.,0.,0.}};
    long                    owner;

    Cell            &cell = m_mesh->getCell(I) ;
    int                     nI = cell.getInterfaceCount() ;
    const long*             interfaces = cell.getInterfaces() ;
    const long*             neighbours = cell.getAdjacencies() ;

    long                    F, N;
    double                  value, area;

    for( i=0; i<nI; ++i){
        F = interfaces[i] ;
        N = neighbours[i] ;

        Interface       &interface = m_mesh->getInterface(F) ;

        owner  = interface.getOwner() ;
        border = interface.isBorder() ;

        area   = m_mesh->evalInterfaceArea(F) ;
        normal = m_mesh->evalInterfaceNormal(F) ;

        if( owner != I)
            normal = -1. *normal ;

        value = border* m_ls[I].value + (1-border) *(m_ls[I].value +m_ls[N].value)/2.  ;

        gradient += area *value *normal ;

    };

    gradient = gradient /m_mesh->evalCellVolume(I) ;


    return gradient;

};

/*!
 * Propagate the sign of the signed distance function from narrow band to entire domain
 */
void LevelSetKernel::propagateSign( std::unordered_map<int,LevelSetObject*> visitors ) {

    // We don't need to propagate the sign in the narrowband
    //
    // An item is in the narrow band if it has a levelset value that differs
    // from the defualt value.
    std::unordered_set<long> alreadyEvaluated;

    PiercedIterator<LSInfo> infoItr = m_ls.begin() ;
    while (infoItr != m_ls.end()) {
        double &value = (*infoItr).value;
        if(!utils::DoubleFloatingEqual()(std::abs(value), levelSetDefaults::VALUE)) {
            alreadyEvaluated.insert(infoItr.getId()) ;
        }
        ++infoItr;
    }

    // If the all cells have the correct value we don't need to progagate the
    // sign
    if (alreadyEvaluated.size() == (size_t) m_mesh->getCellCount()) {
        return;
    }

    // Define the seed candidates
    //
    // First list cells in the narroband, then all other cells. The cells
    // outisde the narrowband will be used as seeds only if there are regions
    // of the mesh disconnected from the narrow band.
    std::vector<long> seedCandidates;
    seedCandidates.reserve(m_mesh->getCellCount());

    seedCandidates.assign(alreadyEvaluated.begin(), alreadyEvaluated.end());
    for (const Cell &cell : m_mesh->getCells()) {
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
        for (long neigh : m_mesh->findCellFaceNeighs( seed )) {
            if (alreadyEvaluated.count(neigh) == 0) {
                processList.push(neigh);
            }
        }

        if (processList.empty()) {
            continue;
        }

        // Discard seeds with a LS value equal to 0
        //
        // If a seed has a value equal to the default value, this means that
        // the cell is outside the narrow band. To have a meaningful levelset
        // value we need to evaulate the levelset from scratch.
        double ls = getLS(seed);
        if(utils::DoubleFloatingEqual()(std::abs(ls), levelSetDefaults::VALUE)) {
            ls = levelSetDefaults::VALUE;
            for (auto entry : visitors) {
                const LevelSetObject *visitor = entry.second ;
                ls = std::min( visitor->evaluateLS(this, seed), ls) ;
            }
        }

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
                infoItr = m_ls.reclaim(id) ;
            }

            // Update the value
            if( infoItr != m_ls.end() ) {
                (*infoItr).value = seedSign * levelSetDefaults::VALUE;
            }

            // Add non-evaluated neighs to the process list
            std::vector<long> neighs = m_mesh->findCellFaceNeighs( id ) ;
            for (long neigh : neighs) {
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
        if (alreadyEvaluated.size() == (size_t) m_mesh->getCellCount()) {
            break;
        }
    }
};

/*!
 * Driver routine for propagtaing levelset value by a fast marching method
 */
void LevelSetKernel::propagateValue( LevelSetObject *visitor ){

    BITPIT_UNUSED(visitor) ;

    // Propagate outwards ------------------------------------------------------- //
    solveEikonal(1.0, 1.0);

    // Propagate inwards -------------------------------------------------------- //
    solveEikonal(1.0, -1.0);


    return ;

};

/*! 
 * Solve the 3D Eikonal equation |grad(u)| = g, using  a fast marching method
 * (Function in the unknown region must be set to the value 1.0e+18)
 * @param[in] g Propagation speed.
 * @param[in] s Velocity sign (+1 --> propagate outwards, -1 --> propagate inwards).
 */
void LevelSetKernel::solveEikonal( double g, double s ){

    long                            N( m_mesh->getCellCount()) ;
    std::unordered_map<long,short>   active ;

    { //FLAG ALIVE, DEAD AND FAR AWAY VERTEXES
        long    m(0), myId ;
        bool    check;

        std::vector<long>               neighs ;
        std::vector<long>::iterator     it, itbeg, itend ;

        for ( const auto &cell : m_mesh->getCells() ){ 

            myId           = cell.getId() ;

            if ( isInNarrowBand(myId) ){ // dead vertex
                active.insert( {{myId, 0}} ) ;

            } else{

                // Loop over neighbors
                check   = false;
                neighs = m_mesh->findCellFaceNeighs(myId) ;

                it    = neighs.begin() ;
                itend = neighs.end() ;

                while( !check &&  it != itend  ){

                    check = (s*m_ls[*it].value >= 0.0) && (abs(m_ls[*it].value) < levelSetDefaults::VALUE) ;
                    ++it ;

                };

                active.insert( {{ myId, 2 - (int) check }} )  ;
                m += (int) check  ;

            }
        } //next i
    }


    { // Construct min heap data structure 
        long                            m(0), I(0), myId, J ;
        double                          value ;

        std::vector<long>               neighs ;
        std::vector<long>::iterator     it, itbeg, itend ;

        std::unordered_map<long,long>   contiguos ;
        std::vector<std::array<int,2>>  map(N), *mapPtr = &map;
        MinPQueue<double, long> heap(m, true, mapPtr);

        for ( const auto &cell : m_mesh->getCells() ){ 

            myId = cell.getId() ;
            contiguos.insert( {{myId,I}} ) ;

            if (active[myId] == 1) {

                // Solve the quadratic form
                value = updateEikonal(s, g, myId, active);

                // Store value into heap
                map[m][0] = I;
                map[I][1] = m; 


                heap.keys[m] = value;
                heap.labels[m] = myId;

                // Update counter
                m++;
            }

            ++I;

        } //next i

        // Build min-heap 
        heap.heap_size = m;
        heap.buildHeap();

        // FAST MARCHING
        while (heap.heap_size > 0) {


            // Extract root
            heap.extract(value, myId);

            LSInfo &lsInfo = m_ls[myId] ;

            // Update level set value
            lsInfo.value = s*updateEikonal(s, g, myId, active);

            // Update m_activeNode to dead node
            active.at(myId) = 0;

            // Upadate far-away neighboors
            neighs = m_mesh->findCellFaceNeighs(myId) ;

            itbeg = neighs.begin() ;
            itend = neighs.end() ;

            for( it=itbeg; it != itend; ++it){

                J = *it;

                // Update the local value
                value = updateEikonal(s, g, J, active);

                if (active.at(J) == 1) { // Update value in min-heap
                    I = contiguos[J] ;
                    heap.modify( map[I][1], value, J );
                }

                else if (active.at(J) == 2) {

                    // Update m_activeNode for neighbor
                    active.at(J) = 1;

                    I = contiguos[J] ;

                    // Insert neighbor into the min heap
                    map[heap.heap_size][0] = I ;
                    map[I][1] = heap.heap_size;

                    heap.insert(value, J);

                };
            };

        } //next item
    }


    return; 
};

/*! 
 * Update scalar field value at mesh vertex on a by locally solving the 3D Eikonal equation.
 * @param[in] s Flag for inwards/outwards propagation (s = -+1).
 * @param[in] g Propagation speed for the 3D Eikonal equation.
 * @param[in] I index of the cartesian cell to be updated.
 * @param[in] active list with dead/ alive/ far information
 * @return Updated value at mesh vertex
 */
double LevelSetKernel::updateEikonal( double s, double g, const long &I, const std::unordered_map<long,short> &active ){

    BITPIT_UNUSED(s) ;
    BITPIT_UNUSED(g) ;
    BITPIT_UNUSED(I) ;
    BITPIT_UNUSED(active) ;

    return levelSetDefaults::VALUE;
};

/*! 
 * Clears all levelset information
 */
void LevelSetKernel::clear( ){
    m_ls.clear() ;
}

/*! 
 * Deletes non-existing items and items outside the narrow band after grid adaption.
 * @param[in] mapper mapping info
 */
void LevelSetKernel::clearAfterMeshMovement( const std::vector<adaption::Info> &mapper ){

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

    return ;
};

/*! 
 * Deletes items outside the narrow band after grid adaption.
 * @param[in] newRSearch new size of narrow band
 */
void LevelSetKernel::filterOutsideNarrowBand( double newRSearch ){

    PiercedIterator<LSInfo> lsItr = m_ls.begin() ;
    while( lsItr != m_ls.end() ){


        if( std::abs(lsItr->value) > newRSearch ){
            lsItr = m_ls.erase( lsItr.getId(), true );
        } else {
            ++lsItr ;
        }

    };

    m_ls.flush() ;

    return ;
};

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetKernel::dump( std::fstream &stream ){

    bitpit::PiercedVector<LSInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    bitpit::genericIO::flushBINARY(stream, m_RSearch);
    bitpit::genericIO::flushBINARY(stream, (long) m_ls.size() ) ;

    for( infoItr=m_ls.begin(); infoItr!=infoEnd; ++infoItr){
        bitpit::genericIO::flushBINARY(stream, infoItr.getId()) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->value) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->gradient) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->object) ;
    };

    return ;
};

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetKernel::restore( std::fstream &stream ){

    long i, n, id;
    LSInfo cellInfo;


    bitpit::genericIO::absorbBINARY(stream, m_RSearch);
    bitpit::genericIO::absorbBINARY(stream, n);

    m_ls.reserve(n);
    for( i=0; i<n; ++i){
        bitpit::genericIO::absorbBINARY(stream, id) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.value) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.gradient) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.object) ;
        m_ls.insert(id, cellInfo) ;
    };

    return ;
};

# if BITPIT_ENABLE_MPI

MPI_Comm LevelSetKernel::getCommunicator() const {
    return m_commMPI;
}

/*!
    Checks if the communicator to be used for parallel communications has
    already been set.

    \result Returns true if the communicator has been set, false otherwise.
*/
bool LevelSetKernel::isCommunicatorSet() const {

    return (getCommunicator() != MPI_COMM_NULL);
}

/*!
    Frees the MPI communicator associated to the patch
*/
void LevelSetKernel::freeCommunicator() {

    if (!isCommunicatorSet()) {
        return;
    }

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled) {
        return;
    }

    MPI_Comm_free(&m_commMPI);
}

/*!
 * Checks if MPI communicator is available in underlying mesh.
 * If available MPI communicator is retreived from mesh and duplicated if necessary and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSetKernel::assureMPI( ){

    if( m_commMPI == MPI_COMM_NULL){

        MPI_Comm meshComm = m_mesh->getCommunicator() ;

        if( meshComm == MPI_COMM_NULL){
            return false;
        } else {
            MPI_Comm_dup(m_mesh->getCommunicator(), &m_commMPI);
            return true; 
        }
    } else {
        return true ;
    };

}

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in/out] sizeBuffer buffer for first communication used to communicate the size of data buffer
 * @param[in/out] dataBuffer buffer for second communication containing data
 */
void LevelSetKernel::writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &sizeBuffer, SendBuffer &dataBuffer ){

    long nItems = sendList.size() ;
    int dataSize = 4*sizeof(double) +sizeof(int) +sizeof(long) ;

    dataBuffer.setCapacity(nItems*dataSize) ;

    //determine elements to send
    long counter(0) ;
    nItems = 0 ;
    for( const auto &index : sendList){
        if( m_ls.exists(index)){
            const auto &lsinfo = m_ls[index] ;
            dataBuffer << counter ;
            dataBuffer << lsinfo.value ;
            dataBuffer << lsinfo.gradient ;
            dataBuffer << lsinfo.object ;
            ++nItems ;
        }
        ++counter ;
    }


    dataBuffer.squeeze( ) ;
    sizeBuffer << nItems ;

    return;
};


/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in] nItems number of items within the buffer
 * @param[in/out] dataBuffer buffer containing the data
 */
void LevelSetKernel::readCommunicationBuffer( const std::vector<long> &recvList, const long &nItems, RecvBuffer &dataBuffer ){

    long    index, id;

    for( int i=0; i<nItems; ++i){
        // Get the id of the element
        dataBuffer >> index ;
        id = recvList[index] ;

        // Assign the data of the element
        PiercedVector<LSInfo>::iterator infoItr ;
        if( !m_ls.exists(id)){
            infoItr = m_ls.reclaim(id) ;
        } else {
            infoItr = m_ls.getIterator(id) ;
        }

        dataBuffer >> infoItr->value ;
        dataBuffer >> infoItr->gradient ;
        dataBuffer >> infoItr->object ;

    }

    return;
};

#endif

}
