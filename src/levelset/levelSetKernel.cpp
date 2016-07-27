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
# include "bitpit_CG.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
    @ingroup levelset
    @interface  LevelSetKernel
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
 * Returns pointer to underlying mesh.
 * @return pointer to mesh
*/
VolumeKernel* LevelSetKernel::getMesh() const{
    return m_mesh ;
} 

/*!
 * Returns reference to LevelSetInfo
*/
PiercedVector<LevelSetInfo>& LevelSetKernel::getLevelSetInfo(){
    return m_ls ;
} 

/*!
 * Returns reference to LevelSetInfo
*/
LevelSetInfo LevelSetKernel::getLevelSetInfo( const long &i)const{
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
int LevelSetKernel::getClosestObject(const long &i) const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::OBJECT;
    } else {
        return (  m_ls[i].object );
    };

};

/*!
 * Get the object and part id of projection point
 * @param[in] i cell index
 * @return pair containing object and part id 
 */
std::pair<int,int> LevelSetKernel::getClosestPart(const long &i) const {

    if( ! m_ls.exists(i) ){
        return ( std::make_pair(levelSetDefaults::OBJECT, levelSetDefaults::PART) ) ;
    } else {
        LevelSetInfo const &ls = m_ls[i] ;
        return (  std::make_pair(ls.object, ls.part) );
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
 * Compute the size of the narrow band using the levelset value
 * @param[in]  mapper mesh modifications
 */
double LevelSetKernel::computeSizeNarrowBandFromLS( ){

    // We need to consider only the cells with a levelset value less than
    // local narrow band (ie. size of the narrowband evalauted using the
    // cell).
    double newRSearch = 0.;
    for (auto itr = m_ls.begin(); itr != m_ls.end(); ++itr) {
        // Discard cells outside the narrow band
        long id = itr.getId() ;
        if( !isInNarrowBand(id) ) {
            continue;
        }

        // Evaluate local search radius
        std::array<double,3> myCenter = computeCellCentroid(id) ;

        const Cell &cell = m_mesh->getCell(id) ;
        const long* neighbours = cell.getAdjacencies() ;
        int N = cell.getAdjacencyCount() ;

        double localRSearch = 0. ;
        for(int n=0; n<N; ++n){
            long neighId = neighbours[n] ;
            if (neighId < 0) {
                continue;
            }

            if( isInNarrowBand(neighId)){
                std::array<double,3> diff = computeCellCentroid(neighId) - myCenter ;
                if( dotProduct(diff, getGradient(id)) > 0){
                    localRSearch = std::max( localRSearch, norm2(diff) );
                }
            }
        } 

        // Discard cells with a levelset greater than the local narrow band
        if ( std::abs( getLS(id) ) > (localRSearch + m_mesh->getTol()) ) {
            continue;
        }

        // Update the levelset
        newRSearch = std::max( localRSearch, newRSearch );
    }

# if BITPIT_ENABLE_MPI
    if( assureMPI() ) {
        double reducedRSearch ;
        MPI_Allreduce( &newRSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, m_commMPI );
        newRSearch = reducedRSearch ;
    }
#endif

    return newRSearch;

};

/*!
 * Clears the geometry cache.
 */
void LevelSetKernel::clearGeometryCache(  ) {

    std::unordered_map<long, std::array<double,3>>().swap( m_cellCentroids ) ;

}

/*!
 * Updates the geometry cache after an adaption.
 */
void LevelSetKernel::updateGeometryCache( const std::vector<adaption::Info> &mapper ) {

    // If there are no cells in the mesh we can just delete all the cache
    if ( m_mesh->getCellCount() == 0) {
        clearGeometryCache();
        return;
    }

    // Remove the previous cells from the cache
    for ( auto & map : mapper ){
        if( map.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for ( auto & previousId : map.previous){
            auto centroidItr = m_cellCentroids.find( previousId ) ;
            if ( centroidItr == m_cellCentroids.end() ) {
                continue ;
            }

            m_cellCentroids.erase( previousId ) ;
        }
    }

}

/*!
 * Computes the centroid of the specfified cell.
 *
 * If the centroid of the cell has been already evaluated, the cached value is
 * returned. Otherwise the cell centroid is evaluated and stored in the cache.
 *
 * @param[in] id is the index of cell
 * @return The centroid of the cell.
 */
const std::array<double,3> & LevelSetKernel::computeCellCentroid( long id ) {

    auto centroidItr = m_cellCentroids.find( id ) ;
    if ( centroidItr == m_cellCentroids.end() ) {
        centroidItr = m_cellCentroids.insert( { id, m_mesh->evalCellCentroid( id ) } ).first ;
    }

    return centroidItr->second;

}

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
        if (F < 0) {
            continue;
        }

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
        if (F < 0) {
            continue;
        }

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
void LevelSetKernel::propagateSign( std::unordered_map<int, std::unique_ptr<LevelSetObject>> &visitors ) {

    // We don't need to propagate the sign in the narrowband
    //
    // An item is in the narrow band if it has a levelset value that differs
    // from the defualt value.
    std::unordered_set<long> alreadyEvaluated;

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

        int nSeedNeighs;
        const long* seedNeighs;
        std::vector<long> seedNeighList;
        if (m_mesh->getCells().size() == 0) {
            seedNeighList = m_mesh->findCellFaceNeighs( seed ) ;
            seedNeighs    = seedNeighList.data() ;
            nSeedNeighs   = seedNeighList.size() ;
        } else {
            Cell& cell  = m_mesh->getCell( seed );
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
        // If a seed has a value equal to the default value, this means that
        // the cell is outside the narrow band. To have a meaningful levelset
        // value we need to evaulate the levelset from scratch.
        double ls = getLS(seed);
        if(utils::DoubleFloatingEqual()(std::abs(ls), levelSetDefaults::VALUE)) {
            ls = levelSetDefaults::VALUE;
            for (auto &entry : visitors) {
                const LevelSetObject *visitor = entry.second.get() ;
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
            if (m_mesh->getCells().size() == 0) {
                neighList = m_mesh->findCellFaceNeighs( id ) ;
                neighs    = neighList.data() ;
                nNeighs   = neighList.size() ;
            } else {
                Cell& cell  = m_mesh->getCell( id );
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

            LevelSetInfo &lsInfo = m_ls[myId] ;

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
void LevelSetKernel::clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

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

    PiercedIterator<LevelSetInfo> lsItr = m_ls.begin() ;
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

    bitpit::PiercedVector<LevelSetInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    bitpit::genericIO::flushBINARY(stream, m_RSearch);
    bitpit::genericIO::flushBINARY(stream, (long) m_ls.size() ) ;

    for( infoItr=m_ls.begin(); infoItr!=infoEnd; ++infoItr){
        bitpit::genericIO::flushBINARY(stream, infoItr.getId()) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->value) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->gradient) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->object) ;
        bitpit::genericIO::flushBINARY(stream, infoItr->part) ;
    };

    return ;
};

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetKernel::restore( std::fstream &stream ){

    long i, n, id;
    LevelSetInfo cellInfo;


    bitpit::genericIO::absorbBINARY(stream, m_RSearch);
    bitpit::genericIO::absorbBINARY(stream, n);

    m_ls.reserve(n);
    for( i=0; i<n; ++i){
        bitpit::genericIO::absorbBINARY(stream, id) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.value) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.gradient) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.object) ;
        bitpit::genericIO::absorbBINARY(stream, cellInfo.part) ;
        m_ls.insert(id, cellInfo) ;
    };

    return ;
};

/*!
 * Checks if the specified cell is inside the given bounding box
 * @param[in] id is the id of the cell
 * @param[in] minPoint is the lower left point of the boungind box
 * @param[in] maxPoint is the upper right point of the boungind box
 * @result True if the specified cell is inside the given bounding box, false
 * otherwise.
 */
double LevelSetKernel::isCellInsideBoundingBox( long id, std::array<double, 3> minPoint, std::array<double, 3> maxPoint ){

    const Cell &cell = m_mesh->getCell(id);
    double tolerance = m_mesh->getTol();

    std::array<double, 3> cellMinPoint;
    std::array<double, 3> cellMaxPoint;

    cellMinPoint.fill(   std::numeric_limits<double>::max() ) ;
    cellMaxPoint.fill( - std::numeric_limits<double>::max() ) ;

    int nVertices = cell.getVertexCount();
    for (int i = 0; i < nVertices; ++i) {
        long vertexId = cell.getVertex(i);
        std::array<double, 3> vertexCoords = m_mesh->getVertexCoords(vertexId);
        for (int d = 0; d < 3; ++d) {
            cellMinPoint[d] = std::min( vertexCoords[d] - tolerance, cellMinPoint[d]) ;
            cellMaxPoint[d] = std::max( vertexCoords[d] + tolerance, cellMaxPoint[d]) ;
        }
    }

    return CGElem::intersectBoxBox(minPoint, maxPoint, cellMinPoint, cellMaxPoint);
};

# if BITPIT_ENABLE_MPI

/*!
 * Returns the MPI communicator stored within LevelSetKernel
 * @return MPI communicator
 */
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

    if (isCommunicatorSet()) {
        return true ;
    }

    MPI_Comm meshComm = m_mesh->getCommunicator() ;
    if( meshComm == MPI_COMM_NULL){
        return false;
    } else {
        MPI_Comm_dup(m_mesh->getCommunicator(), &m_commMPI);
        return true;
    }
}

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] sizeBuffer buffer for first communication used to communicate the size of data buffer
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetKernel::writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &sizeBuffer, SendBuffer &dataBuffer ){

    long nItems = sendList.size() ;
    int dataSize = 4*sizeof(double) +2*sizeof(int) +2*sizeof(long) ;

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
            dataBuffer << lsinfo.part ;
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
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetKernel::readCommunicationBuffer( const std::vector<long> &recvList, const long &nItems, RecvBuffer &dataBuffer ){

    long    index, id;

    for( int i=0; i<nItems; ++i){
        // Get the id of the element
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
        dataBuffer >> infoItr->object ;
        dataBuffer >> infoItr->part ;

    }

    return;
};

#endif

}
