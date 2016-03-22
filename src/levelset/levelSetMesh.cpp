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
    m_cmesh = NULL ;
};

/*!
 * Constructor
 */
LevelSetCartesian::LevelSetCartesian(VolCartesian &patch ): LevelSet( (static_cast<VolumeKernel*>(&patch)) ){
    m_cmesh = &patch ;
};

/*!
 * Compute the levelset function 
 */
void LevelSetCartesian::compute( LSObject *visitor ){

    if( !m_userRSearch)
        computeSizeNarrowBand( visitor ) ;

    visitor->computeLSInNarrowBand(this) ; 

//    if( propagateS ) propagateSign( visitor ) ;
//    if( propagateV ) propagateValue( visitor ) ;

    return;
};

/*!
 * Calculate size of narrow band in order to guarantee one element on each side of geometry
 */
void LevelSetCartesian::computeSizeNarrowBand( LSObject *visitor ){

    BITPIT_UNUSED(visitor) ;

    RSearch = -1.;

    for( int d=0; d<m_cmesh->getDimension(); ++d){
        RSearch = std::max( RSearch, m_cmesh->getSpacing(d) ) ;
    };

    return ;
};

/*!
 * Update the levelset function 
 */
void LevelSetCartesian::update( LSObject *visitor, std::vector<Adaption::Info> &mapper ){

    BITPIT_UNUSED(mapper) ;

    double  newRSearch ;

    if(m_userRSearch){
        newRSearch = RSearch;
    } else{
        newRSearch = updateSizeNarrowBand( mapper ) ;
    };

    visitor->updateLSInNarrowBand( this, mapper, newRSearch ) ;

//TODO    if( propagateS ) updatePropagatedSign() ;
//TODO    if( propagateV ) updatePropagatedValue() ;

    RSearch = newRSearch ;

    return;

};

/*!
 * Update the size of the narrow band after an adaptation of the octree mesh
 * around the linked triangulation.
 */
double LevelSetCartesian::updateSizeNarrowBand( std::vector<Adaption::Info> &mapper ){

    BITPIT_UNUSED(mapper) ;

    double newRSearch = -1. ;

    for( int d=0; d<m_cmesh->getDimension(); ++d){
        newRSearch = std::max( newRSearch, m_cmesh->getSpacing(d) ) ;
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

    Cell&   cell = m_cmesh->getCell(I) ;

    for( d=0; d<m_cmesh->getDimension(); ++d){ // COMPUTE QUADRATIC FORM COEFFICIENTS FROM UPWIND STENCIL

        value   = levelSetDefaults::VALUE ;

        // Left neighbor
        J   = cell.getAdjacency( 2*d, 0) ;

        LSInfo  &lsInfo = info[J] ;

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
            h2 = pow(m_cmesh->getSpacing(d), 2);

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
    m_omesh = NULL ;
};

/*!
 * Constructor
 */
LevelSetOctree::LevelSetOctree(VolOctree & patch ): LevelSet( (static_cast<VolumeKernel*>(&patch)) ){
    m_omesh = &patch ;
};

/*!
 * Compute the levelset function 
 */
void LevelSetOctree::compute( LSObject *visitor ){


    if( !m_userRSearch)
        computeSizeNarrowBand( visitor ) ;

    visitor->computeLSInNarrowBand(this) ; 

//    if( propagateS ) propagateSign( visitor ) ;
//    if( propagateV ) propagateValue( visitor ) ;

    return;

};

/*!
 * Compute the levelset function 
 */
void LevelSetOctree::update( LSObject *visitor, std::vector<Adaption::Info> &mapper ){

    double  newRSearch ;

    if(m_userRSearch){
        newRSearch = RSearch;
    } else {
        newRSearch = updateSizeNarrowBand( mapper ) ;
    };

    visitor->updateLSInNarrowBand( this, mapper, newRSearch ) ;

//TODO    if( propagateS ) updatePropagatedSign() ;
//TODO    if( propagateV ) updatePropagatedValue() ;

    RSearch = newRSearch ;

    return;

};

/*!
 * Initialization of the size of the narrow band around the linked triangulation
 * on the given (and already linked) octree mesh.
 */
void LevelSetOctree::computeSizeNarrowBand( LSObject *visitor ){

    double 						size;

    bool                        flagged ;
    std::vector<bool>           nearPoint ; 

    int                         i, level(100) ;
    uint8_t                     j0(0), j1( pow(2,m_mesh->getDimension())-1) ;

    std::array<double,3>        octrBB0, octrBB1, triBB0, triBB1, C0, C1 ;

    // finest cell in octree
    size = (m_omesh->getTree()).getLocalMinSize();

    m_mesh->getBoundingBox(octrBB0, octrBB1) ;
    visitor->getBoundingBox(triBB0, triBB1) ;

    if( CGElem::intersectBoxBox(octrBB0,octrBB1,triBB0,triBB1,C0,C1) ) { //intersect two Bounding Boxes around geometry and local grid

        // snap bounding box to grid and create cartesian grid
        std::array<int,3>    nc ;

        C0 -= size ;
        C1 += size ;

        for( i=0; i<m_omesh->getDimension(); ++i){
            C0[i] =  octrBB0[i] + size *   (int) ( ( C0[i] - octrBB0[i] ) / size ) ;
            C1[i] =  octrBB0[i] + size * ( (int) ( ( C1[i] - octrBB0[i] ) / size ) +1 ) ;

            nc[i] = round( ( C1[i] - C0[i] ) /size ) ;
        };

        // calculate LS on cartesian mesh and calculate RSearch by finding largest cell throughout flagged cartesian cells
        VolCartesian            cmesh( 0, m_omesh->getDimension(), C0, C1-C0, nc ) ;
        LevelSetCartesian       auxLS(cmesh) ;
        LSObject*               auxSe = visitor->clone() ;

        auxLS.setSign(false) ;
        auxLS.compute( auxSe ) ;
        delete auxSe ;


        std::array<int,3>   i0;
        std::array<int,3>   i1;
        int                 _i, _j, _k, index; 

        for( const auto &cell : m_mesh->getCells() ){

            const long int* conn = cell.getConnect() ;

            C0 = m_mesh->getVertexCoords(conn[j0]) ;
            C1 = m_mesh->getVertexCoords(conn[j1]) ;

            i0 = cmesh.locateClosestVertexCartesian(C0);
            i1 = cmesh.locateClosestVertexCartesian(C1);

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
                level = std::min( level, m_omesh->getCellLevel( cell.getId() ) ) ;

        };


        RSearch = computeRSearchFromLevel( level ) ;

    }; //endif intersect

    //TODO communicate RSearch
    return ;

};

/*!
 * Update the size of the narrow band after an adaptation of the octree mesh
 * around the linked triangulation.
 */
double LevelSetOctree::updateSizeNarrowBand( std::vector<Adaption::Info> &mapper ){

    double  newRSearch ;
    long    id ;
    int     level(100) ;
    std::vector<bool>    map(mapper.size()) ;
    std::vector<bool>::iterator    mapIt=map.begin()  ;

    PiercedVector<long> nb ;
    PiercedIterator<LSInfo> it=info.begin(), itEnd = info.end() ;

    // ========================================================================== //
    // assumes that LS information is relevant to OLD!!! grid
    // scrrens old narrow band for coarsest elements
    
    nb.reserve( info.size() ) ;

    while( it!=itEnd ){
        id = it.getId() ;
        if( isInNarrowBand(id) )
            nb.insert(id,id) ;
        ++it ;
    };


    for ( auto & info : mapper ){

        *mapIt = false ;
        if( info.entity == Adaption::Entity::ENTITY_CELL){

            for ( auto & parent : info.previous){
                id = (long) parent;

                if( isInNarrowBand(id) ){
                    *mapIt = true;
                    nb.erase(id,true) ;
                };
            };

        }
        ++mapIt ;
    }

    nb.flush() ;

    mapIt= map.begin() ;
    for ( auto & info : mapper ){
        if( info.entity == Adaption::Entity::ENTITY_CELL){
            if(*mapIt){ //parent was in narrow band
                for( auto &child : info.current){
                    id = (long) child;
                    nb.insert(id,id) ;
                };
            }; //endif parent in narrow band

        }//if on cell

        ++mapIt ;
    };//foreach mesh modification

    for( auto &id : nb){
        level = min( level, m_omesh->getCellLevel(id) ) ;
    };

    newRSearch = computeRSearchFromLevel(level) ;

    return newRSearch;

};

/*!
 * Compute size of Narrow Band given a the coarsest element level 
 */
double LevelSetOctree::computeRSearchFromLevel( uint8_t level){

    return  (m_omesh->getTree()).levelToSize(level) *sqrt(11.) /2. ;

};

/*!
 * Compute size of Narrow Band given a the coarsest element level 
 */
int LevelSetOctree::computeLevelFromRSearch( double r){

    PabloUniform &tree = m_omesh->getTree() ;

    uint8_t     level ( tree.getLocalMaxDepth() ) ;
    double      size ;

    r = r /sqrt(11.) *2. ; 

    size = tree.levelToSize(level) ;

    while( size < RSearch - 1.e-8 ) {
        level-- ;
        size = tree.levelToSize(level) ;
    };


    return level;
};

}
