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
	@ingroup    levelset
	@class      LevelSetCartesian
	@brief      Implements partially LevelSet for cartesian meshes
*/

/*!
 * Destructor
 */
LevelSetCartesian::~LevelSetCartesian( ){
    m_cartesian = NULL ;
};

/*!
 * Constructor
 */
LevelSetCartesian::LevelSetCartesian(VolCartesian &patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_cartesian = &patch ;
};

/*!
 * Calculate size of narrow band in order to guarantee one element on each side of geometry
 */
double LevelSetCartesian::computeSizeNarrowBand( LevelSetObject *visitor ){

    BITPIT_UNUSED(visitor) ;

    double RSearch = -1.;

    for( int d=0; d<m_cartesian->getDimension(); ++d){
        RSearch = std::max( RSearch, m_cartesian->getSpacing(d) ) ;
    };

    return RSearch;
};

/*!
 * Update the size of the narrow band after an adaptation of the octree mesh
 * around the linked triangulation.
 */
double LevelSetCartesian::updateSizeNarrowBand( const std::vector<adaption::Info> &mapper ){

    BITPIT_UNUSED(mapper) ;

    double newRSearch = -1. ;

    for( int d=0; d<m_cartesian->getDimension(); ++d){
        newRSearch = std::max( newRSearch, m_cartesian->getSpacing(d) ) ;
    };

    return newRSearch;

};

/*! 
 * Update scalar field value at mesh vertex on a by locally solving the 3D Eikonal equation.
 * @param[in] s Flag for inwards/outwards propagation (s = -+1).
 * @param[in] g Propagation speed for the 3D Eikonal equation.
 * @param[in] I index of the cartesian cell to be updated.
 * @return Updated value at mesh vertex
 */
double LevelSetCartesian::updateEikonal( double s, double g, const long &I ){

    int                d;
    long               J;
    double             h2, delta, value, a(0), b(0), c(0);

    Cell&   cell = m_cartesian->getCell(I) ;

    for( d=0; d<m_cartesian->getDimension(); ++d){ // COMPUTE QUADRATIC FORM COEFFICIENTS FROM UPWIND STENCIL

        value   = levelSetDefaults::VALUE ;

        // Left neighbor
        J   = cell.getAdjacency( 2*d, 0) ;

        LSInfo  &lsInfo = m_ls[J] ;

        if( J >= 0 && lsInfo.active == 0){
            value = std::min(s*lsInfo.value, value);
        };


        // Right neighbor
        J   = cell.getAdjacency( 2*d+1, 0) ;

        if( J >= 0 && lsInfo.active == 0){
            value = std::min(s*lsInfo.value, value);
        };


        // Update coeffs in the quadratic form
        if (value < levelSetDefaults::VALUE) {
            h2 = pow(m_cartesian->getSpacing(d), 2);

            a += 1.0/h2;
            b += -2.0 * value/h2;
            c += std::pow(value, 2)/h2;
        }

    };


    { // SOLVE THE QUADRATIC FORM
        // Quadratic form determinant
        delta = pow(b, 2) - 4.0*a*(c - pow(g, 2));

        // Solution
        value = -(b - sqrt(delta))/(2.0*a);

    }

    return(value); 

};

/*!
 *  @ingroup    levelset
 *  @class      LevelSetOctree
 *  @brief      Implements partially LevelSet for octree meshes
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
 * Initialization of the size of the narrow band around the linked triangulation
 * on the given (and already linked) octree mesh.
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
 * around the linked triangulation.
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
                id = (long) parent;

                if( isInNarrowBand(id) ){
                    *mapIt = true;
                    nb.erase(id) ;
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
                    id = (long) child;
                    nb.insert(id) ;
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
 * Compute size of Narrow Band given a the coarsest element level 
 */
double LevelSetOctree::computeRSearchFromLevel( uint8_t level){

    return  (m_octree->getTree()).levelToSize(level) *sqrt(11.) /2. ;

};

/*!
 * Compute size of Narrow Band given a the coarsest element level 
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
