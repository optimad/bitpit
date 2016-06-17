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
        LevelSetObject*         auxSe = visitor->clone() ;

        auxLS.setMesh( &cmesh ) ;
        auxLS.addObject( auxSe ) ;

        auxLS.setSign(false) ;
        auxLS.compute( ) ;
        delete auxSe ;


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
                        flagged = flagged || auxLS.isInNarrowBand(index) ;

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
 * @param[in]  mapper mesh modifications
 */
double LevelSetOctree::updateSizeNarrowBand( const std::vector<adaption::Info> &mapper ){

    double  newRSearch ;
    long    id ;
    int     level(100) ;
    std::vector<bool>    map(mapper.size()) ;
    std::vector<bool>::iterator    mapIt=map.begin()  ;

    PiercedIterator<LSInfo> it=m_ls.begin(), itEnd = m_ls.end() ;

    std::unordered_set<long> nb;

    // assumes that LS information is relevant to OLD!!! grid
    // scrrens old narrow band for coarsest elements
    
    nb.reserve( m_ls.size() ) ;

    while( it!=itEnd ){
        id = it.getId() ;
        if( isInNarrowBand(id) )
            nb.insert(id) ;
        ++it ;
    };


    for ( auto & info : mapper ){

        *mapIt = false ;
        if( info.entity == adaption::Entity::ENTITY_CELL){

            for ( auto & parent : info.previous){

                if( isInNarrowBand(parent) ){
                    *mapIt = true;
                    nb.erase(parent) ;
                };
            };

        }
        ++mapIt ;
    }

    mapIt= map.begin() ;
    for ( auto & info : mapper ){
        if( info.entity == adaption::Entity::ENTITY_CELL){
            if(*mapIt){ //parent was in narrow band
                for( auto &child : info.current){
                    nb.insert(child) ;
                };
            }; //endif parent in narrow band

        }//if on cell

        ++mapIt ;
    };//foreach mesh modification

    for( auto &id : nb){
        level = min( level, m_octree->getCellLevel(id) ) ;
    };

    newRSearch = computeRSearchFromLevel(level) ;

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
