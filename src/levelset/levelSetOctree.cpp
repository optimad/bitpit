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

# include "levelSet.hpp"

# include "bitpit_SA.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

/*!
 *  @ingroup    levelset
 *  @class      LevelSetOctree
 *  @brief      Implements LevelSetKernel for octree meshes
 */

/*!
 * Destructor
 */
LevelSetOctree::~LevelSetOctree( ){
    m_octree = NULL ;
};

/*!
 * Constructor
 */
LevelSetOctree::LevelSetOctree(VolOctree & patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_octree = &patch ;
};

VolOctree* LevelSetOctree::getOctreeMesh() const{
    return m_octree ;
}
/*!
 * Calculation of the necessary size of the narrow band around a LevelSetObject in order to guarantee one cell center on each side of the object
 * @param[in] visitor reference object 
 */
double LevelSetOctree::computeSizeNarrowBand( LevelSetObject *visitor ){

    double                      RSearch(0.) ;
    double 						size;

    bool                        flagged ;
    std::vector<bool>           nearPoint ; 

    int                         i, level(100) ;
    uint8_t                     j0(0), j1( pow(2,m_mesh->getDimension())-1) ;

    std::array<double,3>        octrBB0, octrBB1, triBB0, triBB1, C0, C1 ;

    // finest cell in octree
    size = (m_octree->getTree()).getLocalMinSize();

    m_mesh->getBoundingBox(octrBB0, octrBB1) ;
    visitor->getBoundingBox(triBB0, triBB1) ;

    if( CGElem::intersectBoxBox(octrBB0,octrBB1,triBB0,triBB1,C0,C1) ) { //intersect two Bounding Boxes around geometry and local grid

        // snap bounding box to grid and create cartesian grid
        std::array<int,3>    nc ;

        C0 -= size ;
        C1 += size ;

        for( i=0; i<m_octree->getDimension(); ++i){
            C0[i] =  octrBB0[i] + size *   (int) ( ( C0[i] - octrBB0[i] ) / size ) ;
            C1[i] =  octrBB0[i] + size * ( (int) ( ( C1[i] - octrBB0[i] ) / size ) +1 ) ;

            nc[i] = round( ( C1[i] - C0[i] ) /size ) ;
        };

        // calculate LS on cartesian mesh and calculate RSearch by finding largest cell throughout flagged cartesian cells
        VolCartesian            cmesh( 0, m_octree->getDimension(), C0, C1-C0, nc ) ;
        LevelSet                auxLS ;

        auxLS.setMesh( &cmesh ) ;
        int objectId = auxLS.addObject( std::unique_ptr<LevelSetObject>(visitor->clone()) ) ;

        auxLS.setSign(false) ;
        auxLS.compute( ) ;


        std::array<int,3>   i0;
        std::array<int,3>   i1;
        int                 _i, _j, _k, index, twoDAdjust; 

        twoDAdjust = 3-m_octree->getDimension() ;

        for( const auto &cell : m_mesh->getCells() ){

            const long int* conn = cell.getConnect() ;

            C0 = m_mesh->getVertexCoords(conn[j0]) ;
            C1 = m_mesh->getVertexCoords(conn[j1]) ;

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
                level = std::min( level, m_octree->getCellLevel( cell.getId() ) ) ;

        };


        RSearch = computeRSearchFromLevel( level ) ;

    }; //endif intersect

# if BITPIT_ENABLE_MPI
    if( assureMPI() ){
        double reducedRSearch ;
        MPI_Allreduce( &RSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, m_commMPI );
        RSearch = reducedRSearch ;
    }
#endif

    return RSearch;

};

/*!
 * Update the size of the narrow band after an adaptation of the octree mesh
 * @param[in] mapper mesh modifications
 * @param[in] object the new size of the narrow band will be computed for this object
 */
double LevelSetOctree::updateSizeNarrowBand( const std::vector<adaption::Info> &mapper, LevelSetObject *object ){

    bool coarseningInNarrowBand=false ;

    // assumes that LS information is relevant to OLD!!! grid
    // scrrens old narrow band for coarsest elements

    //
    // Get the bounding box of the objects
    std::array<double, 3> objectsMinPoint;
    std::array<double, 3> objectsMaxPoint;

    object->getBoundingBox( objectsMinPoint, objectsMaxPoint ) ;

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
            if( object->isInNarrowBand(parent) ){
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
    auto const & lsInfo = object->getLevelSetInfo() ;
    for (auto itr = lsInfo.begin(); itr != lsInfo.end(); ++itr) {
        long id = itr.getId() ;
        if( removedNBCells.count(id) > 0 || !object->isInNarrowBand(id) ) {
            continue;
        }

        narrowBandCells.insert(id);
    }

    // Evaluate if the cells in the narrow band are inside the bounding box
    // defined by the objects.
    std::unordered_map<long, bool> isInsideObjectBox;
    for ( long cellId : narrowBandCells ){
        isInsideObjectBox.insert( { cellId, isCellInsideBoundingBox( cellId, objectsMinPoint, objectsMaxPoint ) } );
    }

    //
    // Get the miminum level in the narrow band
    //
    // We need to consider only the cells that are inside the bounding box
    // defined by the objects, or the cells that have at least a neighbout
    // inside the bounding box defined by the objects.
    double newRSearch = 0. ;
    for ( long cellId : narrowBandCells ){
        // Discard cells that are not in the bounding box
        bool discardCell = ! isInsideObjectBox.at(cellId) ;
        if ( discardCell ) {
            for ( long neighId : m_mesh->findCellFaceNeighs(cellId) ) {
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
        newRSearch = std::max( newRSearch, computeRSearchFromCell(cellId) ) ;

    }

    if(coarseningInNarrowBand==false){
        newRSearch = std::min( newRSearch, object->getSizeNarrowBand() ) ;
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
 * Compute size of narrow band given a cell.
 * This function assumes that the octree is balanced 2:1
 * @param[in] id is the id of the cell
 */
double LevelSetOctree::computeRSearchFromCell( long id ){

    int level = m_octree->getCellLevel(id) ;

    return  computeRSearchFromLevel( level ) ;

};

/*!
 * Compute size of narrow band given the coarsest element level which is crossed by geometry.
 * This function assumes that the octree is balanced 2:1
 * @param[in] level the level of the coarsest octree
 */
double LevelSetOctree::computeRSearchFromLevel( uint8_t level){

    return  (m_octree->getTree()).levelToSize(level) *sqrt(11.) /2. ;

};

/*!
 * Compute size of smallest octants greater than a given size (typically the size of the narrow band)
 * @param[in] r limit size
 */
double LevelSetOctree::computeSizeFromRSearch( double r){

    PabloUniform &tree = m_octree->getTree() ;

    uint8_t     level ( tree.getLocalMaxDepth() ) ;
    double      size ;

    size = tree.levelToSize(level) ;

    while( size <= r ) {
        level-- ;
        size = tree.levelToSize(level) ;
    };

    return size ;

};

}
