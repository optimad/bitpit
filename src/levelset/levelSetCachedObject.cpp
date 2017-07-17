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
 * Deletes items outside the narrow band after grid adaption.
 * @param[in] newRSearch new size of narrow band
 */
void LevelSetCachedObject::_filterOutsideNarrowBand( double newRSearch ){

    PiercedVector<LevelSetInfo>::iterator lsItr = m_ls.begin() ;
    while( lsItr != m_ls.end() ){


        if( std::abs(lsItr->value) > newRSearch ){
            if(lsItr->value>0){ 
                lsItr = m_ls.erase( lsItr.getId(), true );
            } else {
                lsItr->value = -1.*levelSetDefaults::VALUE;
                lsItr->gradient = levelSetDefaults::GRADIENT;
                ++lsItr ;
            }
        } else {
            ++lsItr ;
        }

    }

    m_ls.flush() ;

    __filterOutsideNarrowBand(newRSearch) ;

}

/*! 
 * Deletes items outside the narrow band after grid adaption in derived classes
 * @param[in] newRSearch new size of narrow band
 */
void LevelSetCachedObject::__filterOutsideNarrowBand( double newRSearch ){
    BITPIT_UNUSED(newRSearch);
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
 * Calculation of the necessary size of the narrow band for a generic LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::computeSizeNarrowBand(){

    double RSearch(levelSetDefaults::VALUE);

    if( LevelSetCartesian* cartesian = dynamic_cast<LevelSetCartesian*>(m_kernelPtr) ){
        RSearch = _computeSizeNarrowBand(cartesian) ;

    } else if( LevelSetOctree* octree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        RSearch = _computeSizeNarrowBand(octree) ;
    }

    return RSearch;
}

/*!
 * Calculation of the necessary size of the narrow band for a cartesian mesh
 * @param[in] visitee Cartesian LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::_computeSizeNarrowBand( LevelSetCartesian *visitee ){

    VolCartesian const *mesh = visitee->getCartesianMesh();
    double RSearch(0.);

    for( int d=0; d<mesh->getDimension(); ++d){
        RSearch = std::max( RSearch, mesh->getSpacing(d) ) ;
    }

    return RSearch;
}

/*!
 * Calculation of the necessary size of the narrow band for an octree mesh
 * @param[in] visitee Octree LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::_computeSizeNarrowBand( LevelSetOctree *visitee ){

    VolOctree& mesh = *(visitee->getOctreeMesh());
    double  RSearch(0.) ;

    double 						size;

    bool                        flagged ;
    std::vector<bool>           nearPoint ; 

    int                         i, level(100) ;
    uint8_t                     j0(0), j1( pow(2,mesh.getDimension())-1) ;

    std::array<double,3>        octrBB0, octrBB1, triBB0, triBB1, C0, C1 ;

    // finest cell in octree
    size = mesh.getTree().getLocalMinSize();

    visitee->getMesh()->getBoundingBox(octrBB0, octrBB1) ;
    getBoundingBox(triBB0, triBB1) ;

    if( CGElem::intersectBoxBox(octrBB0,octrBB1,triBB0,triBB1,C0,C1) ) { //intersect two Bounding Boxes around geometry and local grid

        // snap bounding box to grid and create cartesian grid
        std::array<int,3>    nc ;

        C0 -= size ;
        C1 += size ;

        for( i=0; i<mesh.getDimension(); ++i){
            C0[i] =  octrBB0[i] + size *   (int) ( ( C0[i] - octrBB0[i] ) / size ) ;
            C1[i] =  octrBB0[i] + size * ( (int) ( ( C1[i] - octrBB0[i] ) / size ) +1 ) ;

            nc[i] = round( ( C1[i] - C0[i] ) /size ) ;
        }

        // calculate LS on cartesian mesh and calculate RSearch by finding largest cell throughout flagged cartesian cells
        VolCartesian            cmesh( 0, mesh.getDimension(), C0, C1-C0, nc ) ;
        LevelSet                auxLS ;

        auxLS.setMesh( &cmesh ) ;
        int objectId = auxLS.addObject( std::unique_ptr<LevelSetObject>(clone()) ) ;

        auxLS.setSign(false) ;
        auxLS.compute( ) ;

        const LevelSetObject &obj = auxLS.getObject(objectId);

        std::array<int,3>   i0;
        std::array<int,3>   i1;
        int                 _i, _j, _k, index, twoDAdjust; 

        twoDAdjust = 3-mesh.getDimension() ;

        for( const auto &cell : mesh.getCells() ){

            const long int* conn = cell.getConnect() ;

            C0 = mesh.getVertexCoords(conn[j0]) ;
            C1 = mesh.getVertexCoords(conn[j1]) ;

            i0 = cmesh.locateClosestVertexCartesian(C0);
            i1 = cmesh.locateClosestVertexCartesian(C1);

            i1[2] +=  twoDAdjust ;

            flagged = false ;

            for( _k=i0[2]; _k<i1[2]; ++_k){
                for( _j=i0[1]; _j<i1[1]; ++_j){
                    for( _i=i0[0]; _i<i1[0]; ++_i){

                        index = cmesh.getCellLinearId( _i, _j, _k) ;
                        flagged = flagged || obj.isInNarrowBand(index) ;

                    }
                }
            }
            
            if(flagged) 
                level = std::min( level, mesh.getCellLevel( cell.getId() ) ) ;

        }

        RSearch = visitee->computeRSearchFromLevel( level ) ;

    } //endif intersect

# if BITPIT_ENABLE_MPI
    if( assureMPI() ){
        double reducedRSearch ;
        MPI_Comm meshComm = visitee->getCommunicator() ;
        MPI_Allreduce( &RSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, meshComm );
        RSearch = reducedRSearch ;
    }
#endif

    return RSearch;
}

/*!
 * Update the size of the narrow band after an adaptation of the mesh
 * @param[in] mapper mesh modifications
 * @return size of narrow band
 */
double LevelSetCachedObject::updateSizeNarrowBand(const std::vector<adaption::Info> &mapper){

    double R(levelSetDefaults::VALUE);

    if( LevelSetCartesian *cartesian = dynamic_cast<LevelSetCartesian*>(m_kernelPtr) ){
        R= _updateSizeNarrowBand(cartesian,mapper) ;

    } else if( LevelSetOctree *octree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        R= _updateSizeNarrowBand(octree,mapper) ;
    }

    return R;
}

/*!
 * Update the size of the narrow band after an adaptation of the cartesian mesh
 * @param[in] mapper mesh modifications
 * @param[in] visitee Cartesian LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::_updateSizeNarrowBand(LevelSetCartesian *visitee, const std::vector<adaption::Info> &mapper){

    BITPIT_UNUSED(mapper);

    VolCartesian const &mesh = *(visitee->getCartesianMesh());
    double newRSearch(0.);

    for( int d=0; d<mesh.getDimension(); ++d){
        newRSearch = std::max( newRSearch, mesh.getSpacing(d) ) ;
    }

    return newRSearch;
}

/*!
 * Update the size of the narrow band after an adaptation of the octree mesh
 * @param[in] mapper mesh modifications
 * @param[in] visitee Octree LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::_updateSizeNarrowBand(LevelSetOctree *visitee, const std::vector<adaption::Info> &mapper){

    VolOctree const &mesh = *(visitee->getOctreeMesh());
    double  newRSearch(0.) ;

    bool coarseningInNarrowBand=false ;

    // assumes that LS information is relevant to OLD!!! grid
    // scrrens old narrow band for coarsest elements

    //
    // Get the bounding box of the objects
    std::array<double, 3> objectsMinPoint;
    std::array<double, 3> objectsMaxPoint;
    getBoundingBox( objectsMinPoint, objectsMaxPoint ) ;

    //
    // Cells in the narrow band
    //
    std::unordered_set<long> narrowBandCells;

    // First insert the new cells added to the narrow band
    std::unordered_set<long> removedNBCells;
    for ( auto & info : mapper ){
        if( info.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        bool parentInNarrowBand = false;
        for ( auto & parent : info.previous){
            if( isInNarrowBand(parent) ){
                parentInNarrowBand = true;
                if( info.type == adaption::Type::TYPE_COARSENING ){
                    coarseningInNarrowBand=true ;
                }
                break;
            }
        }

        if (parentInNarrowBand) {
            for ( auto & parent : info.previous){
                removedNBCells.insert(parent);
            }

            for( auto &child : info.current){
                narrowBandCells.insert(child);
            }
        }
    }

    // Now add the cells that were in the narrow band befor and have not been
    // updated
    for (auto cell : mesh.getCells()) {
        long id = cell.getId() ;
        if( removedNBCells.count(id) > 0 || !isInNarrowBand(id) ) {
            continue;
        }

        narrowBandCells.insert(id);
    }

    // Evaluate if the cells in the narrow band are inside the bounding box
    // defined by the objects.
    std::unordered_map<long, bool> isInsideObjectBox;
    for ( long cellId : narrowBandCells ){
        isInsideObjectBox.insert( { cellId, visitee->isCellInsideBoundingBox( cellId, objectsMinPoint, objectsMaxPoint ) } );
    }

    //
    // Get the miminum level in the narrow band
    //
    // We need to consider only the cells that are inside the bounding box
    // defined by the objects, or the cells that have at least a neighbout
    // inside the bounding box defined by the objects.
    std::vector<long> faceNeighs;
    for ( long cellId : narrowBandCells ){
        // Discard cells that are not in the bounding box
        bool discardCell = ! isInsideObjectBox.at(cellId) ;
        if ( discardCell ) {
            faceNeighs.clear();
            mesh.findCellFaceNeighs(cellId, &faceNeighs);
            for ( long neighId : faceNeighs ) {
                if ( narrowBandCells.count(neighId) == 0 ) {
                    continue;
                }

                if (isInsideObjectBox.at(neighId)) {
                    discardCell = false ;
                    break ;
                }
            }
        }

        if (discardCell) {
            continue;
        }

        // Update the level
        newRSearch = std::max( newRSearch, visitee->computeRSearchFromCell(cellId) ) ;

    }

    if(coarseningInNarrowBand==false){
        newRSearch = std::min( newRSearch, getSizeNarrowBand() ) ;
    }

# if BITPIT_ENABLE_MPI
    if( assureMPI() ) {
        double reducedRSearch ;
        MPI_Comm meshComm = visitee->getCommunicator() ;
        MPI_Allreduce( &newRSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, meshComm );
        newRSearch = reducedRSearch ;
    }
#endif

    return newRSearch;

}

/*!
 * Propagate the sign of the signed distance function from narrow band to
 * entire domain.
 */
void LevelSetCachedObject::propagateSign() {

    VolumeKernel const &mesh = *(m_kernelPtr->getMesh()) ;

    VolumeKernel::CellConstIterator cellBegin = mesh.cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh.cellConstEnd();

    PiercedStorage<bool, long> signAssigned(1, &(m_kernelPtr->getMesh()->getCells()));
    signAssigned.fill(false);

    long nUnassigned = mesh.getCellCount();

    std::vector<long> seeds;
    seeds.reserve(mesh.getCellCount());

    // Evaluate the bounding box of the object
    std::array<double,3> boxMin;
    std::array<double,3> boxMax;
    getBoundingBox(boxMin, boxMax);

    // Identify the cells that have a known sign
    //
    // The sign of the cells in the narrowband is known because for those cells
    // we know the value of the levelset. We don't need to set the sign of the
    // cells in the narrowband, instead they will be used as seeds for
    // propagating the sign to other cells.
    //
    // The sign of external cells will always be positive. An external cell is
    // a cell outside the bounding box of the object. Alsot the external cells
    // can be usesd as seeds for sign propagation.
    for (auto itr = cellBegin; itr != cellEnd; ++itr) {
        long cellId = itr.getId();

        // Process cells in the narrowband
        if (isInNarrowBand(cellId)) {
            std::size_t cellRawId = itr.getRawIndex();
            signAssigned.rawAt(cellRawId) = true;
            --nUnassigned;

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
            setSign(cellId, 1);
            std::size_t cellRawId = itr.getRawIndex();
            signAssigned.rawAt(cellRawId) = true;
            --nUnassigned;

            seeds.push_back(cellId);

            continue;
        }
    }

    // Use the seeds to propagate the sign
    propagateSeedSign(seeds, &nUnassigned, &signAssigned);

#if BITPIT_ENABLE_MPI
    // If there are cells with an unknown sign, data communication among
    // ghost cells is needed. At this point, it is only possible to have
    // cells with an unknown sign for It is not possible to have
    long nGlobalUnassigned = nUnassigned;
    if (mesh.isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &nGlobalUnassigned, 1, MPI_LONG, MPI_SUM, m_kernelPtr->getCommunicator());
    }

    if (nGlobalUnassigned != 0) {
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
        while (nGlobalUnassigned != 0) {
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
                    int sign = getSign(cellId);
                    buffer << sign;
                }

                dataCommunicator.startSend(rank);
            }

            // Receive the sign
            //
            // If we discovee the sign of a ghost, we can usi it as a seed.
            int nCompletedRecvs = 0;
            while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
                int rank = dataCommunicator.waitAnyRecv();
                const auto &recvIds = mesh.getGhostExchangeTargets(rank);
                RecvBuffer &buffer = dataCommunicator.getRecvBuffer(rank);

                // Receive data and detect new seeds
                seeds.clear();
                for (long cellId : recvIds) {
                    buffer >> sign;

                    std::size_t cellRawId = mesh.getCells().getRawIndex(cellId);
                    if (!signAssigned.rawAt(cellRawId)) {
                        if (sign != 0) {
                            seeds.push_back(cellId);
                            signAssigned.rawAt(cellRawId) = true;
                        }
                    } else {
                        assert(getSign(cellId) == sign);
                    }
                }

                ++nCompletedRecvs;

                // Propagate the sign of the newly identified seeds
                if (seeds.size() > 0) {
                    propagateSeedSign(seeds, &nUnassigned, &signAssigned);
                }
            }

            // Wait to the sends to finish
            dataCommunicator.waitAllSends();

            // Update the global counter for cells with an unknow sign
            nGlobalUnassigned = nUnassigned;
            MPI_Allreduce(MPI_IN_PLACE, &nGlobalUnassigned, 1, MPI_LONG, MPI_SUM, m_kernelPtr->getCommunicator());
        }
    }
#else
    // Check that the sign has been propagated into all regions
    assert(nUnassigned == 0)
#endif
}

/*!
 * Propagate the sign of the signed distance function from the specified
 * seeds to all reachable cells.
 *
 * \param seeds are the seeds to be used for sign propagation
 * \param[in,out] nUnassigned is the number of cells for which the sign is
 * still unassigned. On output, this number will be updated so it's possible
 * to keep track of the cells whose sign is not yet assigned
 * \param[in,out] signAssigned is a flag that states if the sign of a
 * cell has been assigned. On output, this flag will be updated for
 * all the cells whose sign has been set
 * \result The number of cells for which the sign has been set during this
 * function call.
 */
void LevelSetCachedObject::propagateSeedSign(const std::vector<long> &seeds,
                                             long *nUnassigned, PiercedStorage<bool, long> *signAssigned) {

    VolumeKernel const &mesh = *(m_kernelPtr->getMesh()) ;

    std::vector<long> processList;

    std::size_t seedCursor = seeds.size();
    while (seedCursor != 0) {
        // Get a seed
        --seedCursor;
        long seedId = seeds[seedCursor];

        // Clear process list
        processList.clear();

        // Start filling the process list with the neighbours of the seed
        // that have an unknown sign. If the seed is sourrounded only by
        // cell for which the sign value is know, the seed can't be used
        // to propagate the sign
        const Cell &seed = mesh.getCell(seedId);
        const long *seedNeighs = seed.getAdjacencies() ;
        int nSeedNeighs = seed.getAdjacencyCount() ;
        for(int n = 0; n < nSeedNeighs; ++n){
            long neighId = seedNeighs[n] ;
            if (neighId < 0) {
                continue;
            } else if (signAssigned->at(neighId)) {
                continue;
            }

            processList.push_back(neighId);
        }

        if (processList.empty()) {
            continue;
        }

        // Get the sign of the seed
        int seedSign = getSign(seedId);

        // Propagate the sign
        while (!processList.empty()) {
            long cellId = processList.back();
            processList.resize(processList.size() - 1);

            // Consider only cells with an unknown sign
            //
            // The process list is not unique, it can contain duplcates, so
            // we need to make sure not to process a cell multiple times.
            // If the sign of the cell is know, this means that it has
            // already been processed.
            if (signAssigned->at(cellId)) {
                continue;
            }

            // Set the sign of the cell
            setSign(cellId, seedSign);
            signAssigned->at(cellId) = true;
            --(*nUnassigned);

            // If there are no more unassigned cell we can exit
            if (*nUnassigned == 0) {
                return;
            }

            // Add cells with unknow sign to the process list
            const Cell &cell = mesh.getCell(cellId);
            const long *cellNeighs = cell.getAdjacencies() ;
            int nCellNeighs = cell.getAdjacencyCount() ;
            for(int n = 0; n < nCellNeighs; ++n){
                long neighId = cellNeighs[n] ;
                if (neighId < 0) {
                    continue;
                } else if (signAssigned->at(neighId)) {
                    continue;
                }

                processList.push_back(neighId);
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

    long nItems(0), counter(0) ;


    //determine number of elements to send
    for( const auto &index : sendList){
        if( m_ls.exists(index)){
            nItems++ ;
        }

    }

    dataBuffer << nItems ;
    dataBuffer.setSize(dataBuffer.getSize() +nItems* (sizeof(long) +4*sizeof(double)) ) ;


    for( const auto &index : sendList){
        auto lsInfoItr = m_ls.find(index) ;
        if( lsInfoItr != m_ls.end() ){
            dataBuffer << counter ;
            dataBuffer << lsInfoItr->value ;
            dataBuffer << lsInfoItr->gradient ;
            ++nItems ;
        }
        ++counter ;
    }

    __writeCommunicationBuffer( sendList, dataBuffer) ;

    dataBuffer.squeeze( ) ;
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

    long    nItems, index, id;

    dataBuffer >> nItems ;

    for( int i=0; i<nItems; ++i){

        // Determine the id of the element
        dataBuffer >> index ;
        id = recvList[index] ;

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
