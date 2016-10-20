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
};

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetCachedObject::LevelSetCachedObject(int id) : LevelSetObject(id,true){
};

/*!
 * Get LevelSetInfo of cell
 * @param[in] i cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetCachedObject::getLevelSetInfo( const long &i)const{
    if( ! m_ls.exists(i) ){
        return (  LevelSetInfo() );
    } else {
        return m_ls[i] ;
    };

} 

/*!
 * Get the levelset value of cell
 * @param[in] i cell index
 * @return levelset value in cell
 */
double LevelSetCachedObject::getLS( const long &i)const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::VALUE;
    } else {
        return (  m_ls[i].value );
    };

};

/*!
 * Get the levelset gradient of cell
 * @param[in] i cell index
 * @return levelset gradient in cell 
 */
std::array<double,3> LevelSetCachedObject::getGradient(const long &i) const {

    if( ! m_ls.exists(i) ){
        return levelSetDefaults::GRADIENT;
    } else {
        return (  m_ls[i].gradient );
    };

};

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
};

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

    PiercedIterator<LevelSetInfo> lsItr = m_ls.begin() ;
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

    };

    m_ls.flush() ;

    __filterOutsideNarrowBand(newRSearch) ;

    return ;
};

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
 * @param[in] visitee LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::computeSizeNarrowBand( LevelSetKernel *visitee ){

    double RSearch;

    if( dynamic_cast<LevelSetCartesian*>(visitee) != nullptr ){
        LevelSetCartesian* cartesian = dynamic_cast<LevelSetCartesian*>(visitee) ;
        RSearch = _computeSizeNarrowBand(cartesian) ;

    } else if( dynamic_cast<LevelSetOctree*>(visitee) != nullptr ){
        LevelSetOctree* octree = dynamic_cast<LevelSetOctree*>(visitee) ;
        RSearch = _computeSizeNarrowBand(octree) ;
    }

    return RSearch;
};

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
    };

    return RSearch;
};

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
        };

        // calculate LS on cartesian mesh and calculate RSearch by finding largest cell throughout flagged cartesian cells
        VolCartesian            cmesh( 0, mesh.getDimension(), C0, C1-C0, nc ) ;
        LevelSet                auxLS ;

        auxLS.setMesh( &cmesh ) ;
        int objectId = auxLS.addObject( std::unique_ptr<LevelSetObject>(clone()) ) ;

        auxLS.setSign(false) ;
        auxLS.compute( ) ;

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
                        flagged = flagged || auxLS.isInNarrowBand(index,objectId) ;

                    };
                };
            };
            
            if(flagged) 
                level = std::min( level, mesh.getCellLevel( cell.getId() ) ) ;

        };

        RSearch = visitee->computeRSearchFromLevel( level ) ;

    }; //endif intersect

# if BITPIT_ENABLE_MPI
    if( assureMPI(visitee) ){
        double reducedRSearch ;
        MPI_Comm meshComm = visitee->getCommunicator() ;
        MPI_Allreduce( &RSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, meshComm );
        RSearch = reducedRSearch ;
    }
#endif

    return RSearch;
};

/*!
 * Update the size of the narrow band after an adaptation of the mesh
 * @param[in] mapper mesh modifications
 * @param[in] visitee Cartesian LevelSetKernel
 * @return size of narrow band
 */
double LevelSetCachedObject::updateSizeNarrowBand(LevelSetKernel *visitee, const std::vector<adaption::Info> &mapper){

    double R;

    if( dynamic_cast<LevelSetCartesian*>(visitee) != nullptr ){
        LevelSetCartesian* cartesian = dynamic_cast<LevelSetCartesian*>(visitee) ;
        R= _updateSizeNarrowBand(cartesian,mapper) ;

    } else if( dynamic_cast<LevelSetOctree*>(visitee) != nullptr ){
        LevelSetOctree* octree = dynamic_cast<LevelSetOctree*>(visitee) ;
        R= _updateSizeNarrowBand(octree,mapper) ;
    }

    return R;
};

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
    };

    return newRSearch;
};

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
    for ( long cellId : narrowBandCells ){
        // Discard cells that are not in the bounding box
        bool discardCell = ! isInsideObjectBox.at(cellId) ;
        if ( discardCell ) {
            for ( long neighId : mesh.findCellFaceNeighs(cellId) ) {
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
    if( assureMPI(visitee) ) {
        double reducedRSearch ;
        MPI_Comm meshComm = visitee->getCommunicator() ;
        MPI_Allreduce( &newRSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, meshComm );
        newRSearch = reducedRSearch ;
    }
#endif

    return newRSearch;

};

/*!
 * Propagate the sign of the signed distance function from narrow band to entire domain
 */
void LevelSetCachedObject::propagateSign( LevelSetKernel *visitee ) {

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
void LevelSetCachedObject::assignSign(int sign, const std::unordered_set<long> &cells) {

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
 * Writes LevelSetCachedObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_dump( std::ostream &stream ){

    IO::binary::write(stream, (long) m_ls.size() ) ;
    bitpit::PiercedVector<LevelSetInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    for( infoItr=m_ls.begin(); infoItr!=infoEnd; ++infoItr){
        IO::binary::write(stream, infoItr.getId()) ;
        IO::binary::write(stream, infoItr->value) ;
        IO::binary::write(stream, infoItr->gradient) ;
    };

    __dump(stream) ;
};

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

    IO::binary::read(stream, n);

    m_ls.reserve(n);
    for( i=0; i<n; ++i){
        IO::binary::read(stream, id) ;
        IO::binary::read(stream, cellInfo.value) ;
        IO::binary::read(stream, cellInfo.gradient) ;
        m_ls.insert(id, cellInfo) ;
    };

    __restore(stream) ;
};

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

    __writeCommunicationBuffer( sendList, dataBuffer) ;

    dataBuffer.squeeze( ) ;

    return;
};

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
        PiercedVector<LevelSetInfo>::iterator infoItr ;
        if( !m_ls.exists(id)){
            infoItr = m_ls.reclaim(id) ;
        } else {
            infoItr = m_ls.getIterator(id) ;
        }

        dataBuffer >> infoItr->value ;
        dataBuffer >> infoItr->gradient ;
    }

    __readCommunicationBuffer( recvList, dataBuffer ) ;

    return;
};

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
