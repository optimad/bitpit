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

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

std::set<long> LevelSetSegmentation::NULL_LIST ;
long LevelSetSegmentation::NULL_ELEMENT=-1 ;

/*!
	@ingroup    levelset
	@class      LSObject
	@brief      LevelSet object.
*/

/*!
 * Destructor
 */
LSObject::~LSObject( ){
};

/*!
 * Constructor
 * @param[in] id Id assigned to object
 */
LSObject::LSObject( int id) : m_id(id){
};

/*!
 * Get the id 
 * @return Id of the object
 */
int LSObject::getId( ) const {
    return m_id ;
};

/*!
	@ingroup    levelset
	@class      LevelSetSegmentation
	@brief      Implements visitor pattern fo segmentated geometries
*/
/*!
 * Constructor of LevelSet_Stl with input parameters.
 * @param[in] id id to be asigned to pierced vector
 * @param[in] list list of simplices
 */
LevelSetSegmentation::SegData::SegData( ) : m_segments(NULL_LIST), m_support(NULL_ELEMENT){
};

/*!
 * Constructor of LevelSet_Stl with input parameters.
 * @param[in] id id to be asigned to pierced vector
 * @param[in] list list of simplices
 */
LevelSetSegmentation::SegData::SegData( const std::set<long> &list) :m_segments(list), m_support(NULL_ELEMENT) {
};

/*!
 * Constructor of LevelSet_Stl with input parameters.
 * @param[in] id id to be asigned to pierced vector
 * @param[in] list list of simplices
 * @param[in] support index of closest simplex
 */
LevelSetSegmentation::SegData::SegData( const std::set<long> &list, const long &support) :m_segments(list), m_support(support){
};

/*!
 * Destructor of LevelSet_Stl.
 */
LevelSetSegmentation::~LevelSetSegmentation() {
    stl = NULL;
};

/*!
 * Constructor of LevelSet_Stl with input parameters.
 * It builds one Sdf object with :
 * @param[in] *Pmesh Pointer to 3D pablo octree mesh to be linked to Sdf.
 * @param[in] *STL Pointer to surface Triangulation to be linked to Sdf.
 */
LevelSetSegmentation::LevelSetSegmentation( int id, SurfUnstructured *STL) :LSObject(id) {

    stl = STL;
    abs_tol = 1.0e-12 ;

};

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetSegmentation::LevelSetSegmentation( const LevelSetSegmentation &other) :LSObject(other.getId() ) {

    stl = other.stl; 
    abs_tol = other.abs_tol ;

};

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetSegmentation* LevelSetSegmentation::clone() const {
    return new LevelSetSegmentation( *this ); 
}

/*!
 * Get the list of simplices wich contain the i-th local element in their narrow band.
 * @param[in] i Local index of target octant.
 * @return set with indices of simplices wich contain the i-th local element in their narrow band.
 */
const std::set<long> & LevelSetSegmentation::getSimplexList(const long &i){

    if( !m_segInfo.exists(i) ){
        return NULL_LIST;
    } else {
        return ( m_segInfo[i].m_segments );
    };

};

/*!
 * Get the index of the closest simplex
 * @param[in] i Local index of target octant.
 * @return index of closest simplex
 */
const long & LevelSetSegmentation::getSupportSimplex(const long &i){

    if( !m_segInfo.exists(i) ){
        return NULL_ELEMENT ;
    } else {
        return ( m_segInfo[i].m_support );
    };
};

/*!
 * Check if cell is in narrowband of any triangle;
 * LevelSetSegmentation::associateSimplexToCell() should have been called prior;
 * @param[in] i Local index of target octant.
 * @return if in narow band
 */
bool LevelSetSegmentation::isInNarrowBand(const long &i){
    return( m_segInfo.exists(i) ) ; 
};

/*!
 * Add a SurfTriPatch object to stl.
 * @param[in] *STL Pointer to surface Triangulation to be added .
 */
std::vector<std::array<double,3>> LevelSetSegmentation::getSimplexVertices( const long &i ) const {

    Cell &cell = stl->getCell(i) ;

    int                                     j, n, N (cell.getVertexCount()) ;
    std::vector<std::array<double,3>>       VS(N) ;

    for( n=0; n<N; ++n){
        j = cell.getVertex(n) ;
        VS[n] = stl->getVertexCoords(j);
    };

    return VS;
};

/*!
 * Update the signed distance function for a given octant with the list of simplices
 * which contain the octant in the narrow band (if et[IDX] is not void) otherwise
 * by using the full triangulation.
 * @param[in] visitee visited mesh 
 * @param[in] search size of narrow band
 * @param[in] filter if triangles should be ereased when outside narrow band (default false)
 */
void LevelSetSegmentation::lsFromSimplex( LevelSet *visitee, const double &search, bool filter){

    VolumeKernel                         &mesh  = *(visitee->m_mesh) ;
    bool                                &signd = visitee->signedDF ;

    long                                id ;
    double                              s, d, value;
    std::array<double,3>                n, xP, P;

    std::set<long>::iterator            it, itend ;
    PiercedIterator<SegData>    segIt, segEnd = m_segInfo.end() ;

    LevelSet::LSInfo                    lsInfo  ;

    for( segIt=m_segInfo.begin(); segIt!=segEnd; ++segIt ){

        id = segIt.getId() ;
        SegData &segInfo = *segIt ;

        lsInfo.value = visitee->getLS( id ) ;

        std::set<long>                      &segs = segInfo.m_segments ;
        long                                &supp = segInfo.m_support ;

        it    = segs.begin();
        itend = segs.end() ;

        P       = mesh.evalCellCentroid(id) ;
        value   = abs( lsInfo.value );

        while( it != itend ){

            infoFromSimplex(P, *it, d, s, xP, n);

            if ( d <= search ){

                if( d<value ) {
                    s       = signd *s + (!signd) *1.;
                    value   = d ;

                    lsInfo.object   = getId();
                    lsInfo.value    = s *d; //TODO check
                    lsInfo.gradient = s *n ;
                    supp            = *it ;
                }

                ++it ;

            } //end if distance

            else {
                if( filter){ 
                    it = segs.erase(it) ; 
                }

                else{
                    ++it ;
                };
            };


        } //end foreach triangle

        if( abs(lsInfo.value) < 1.e17 ){
            if( !visitee->info.exists(id) ){
                visitee->info.insert(id,lsInfo) ;
            } else {
                visitee->info[id]= lsInfo ;
            };
        }

    };// foreach cell

    m_segInfo.flush() ;

    return;

};

/*!
 * Update the signed distance function for a given octant with the list of simplices
 * which contain the octant in the narrow band (if et[IDX] is not void) otherwise
 * by using the full triangulation.
 * @param[in] p coordinates of point
 * @param[in] i index of triangle
 * @param[out] d distance point triangle
 * @param[out] s sign of point wrt to triangle, i.e. according to normal
 * @param[out] x closest point on trinagle
 * @param[out] n normal at closest point on trinagle
 */
void LevelSetSegmentation::infoFromSimplex( const std::array<double,3> &p,const  long &i, double &d, double &s, std::array<double,3> &x, std::array<double,3> &n ) const {

    std::vector<std::array<double,3>>   VS( getSimplexVertices(i) ) ;
    int                                 where;

    d = CGElem::distancePointSimplex( p, VS, x, where);

    if (where == 0) {
        n = stl->evalFacetNormal(i);

    } else if (where > 0) {
        n = stl->evalEdgeNormal(i,where-1) ;

    } else if (where < 0) {
        n = stl->evalVertexNormal(i,-where-1) ;
    }

    s = sign( dotProduct(n, p - x) );

    if( d < abs_tol ){
        n = s *(p - x) /d ; //TODO check
    };

    return ;

};

/*!
 * Finds seed points for given simplex
 * @param[in] VS Simplex
 * @param[out] I indices of seed points
 */
bool LevelSetSegmentation::seedNarrowBand( LevelSetCartesian *visitee, std::vector<std::array<double,3>> &VS, std::vector<int> &I){

    VolumeKernel                 &mesh = *(visitee->m_mesh) ;

    bool                                found(false) ;
    int                                 dim( mesh.getDimension() ) ;
    std::array<double,3>                B0, B1;
    std::vector<std::array<double,3>>   VP ;

    mesh.getBoundingBox(B0, B1) ;

    B0 = B0 - visitee->RSearch ;
    B1 = B1 + visitee->RSearch ;

    I.clear() ;

    for( const auto &P : VS){
        if(  CGElem::intersectPointBox( P, B0, B1, dim ) ) {
            I.push_back( mesh.locatePoint(P) );
            found =  true ;
        };
    }

    if( !found && CGElem::intersectBoxSimplex( B0, B1, VS, VP, dim ) ) {
        for( const auto &P : VP){
            I.push_back( mesh.locatePoint(P) );
            found = true ;
        };
    }

    return found ;
};

/*!
 * Finds a seed point to be used for propgating the levelset sign
 * @param[out] seed index of seed cell
 * @param[out] sign sign at seed cell [+1;-1]
 */
void LevelSetSegmentation::seedSign( LevelSet *visitee, long &seed, double &value) const {


    double                      d, s, lsTmp;
    int                         i, M, N( visitee->m_mesh->getCellCount() );
    std::array<double,3>        P, X, temp;

    seed = -1;

    //------------------------------------------------------------------------------
    // FIND SEED ELEMENT In GRID                                                       //
    i = 0;
    while ( i < N && seed == -1) {
        if (abs(visitee->info[i].value) < 1.0e+17) {
            seed = i;
            value = sign(visitee->info[i].value) ;
        }
        ++i;
    } //next oct

    //--------------------------------------------------------------------------- 
    // if no seed element found, check whole triangulation                                                     
    if( seed == -1 ){

        seed    = int(N/2);
        lsTmp  = 1.0e+18;

        M       = stl->getCellCount();
        P       = visitee->m_mesh->evalCellCentroid(seed) ;

        for( i=0; i<M; ++i){

            infoFromSimplex(P, i, d, s, X, temp);

            d = abs(d) ;
            if ( d < lsTmp) {
                lsTmp = d ;
                value = s ;
            };

        };
    };

    return;


};

/*!
 * Computes axis aligned bounding box of object
 * @param[out] minP minimum point
 * @param[out] maxP maximum point
 */
void LevelSetSegmentation::getBoundingBox( std::array<double,3> &minP, std::array<double,3> &maxP ) const {
    stl->getBoundingBox(minP,maxP) ;
};

/*!
 * Compute the levelset function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. 
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetCartesian *visitee ){

    associateSimplexToCell(visitee) ;
    lsFromSimplex(visitee, visitee->RSearch) ;

    return;
};

/*!
 * Determines the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle).
 */
void LevelSetSegmentation::associateSimplexToCell( LevelSetCartesian *visitee ){

    VolumeKernel                     &mesh = *(visitee->m_mesh) ;
    std::vector<std::array<double,3>>       VS(3);

    std::vector< int >                      stack, temp ;
    std::vector< std::array<double,3> >     cloud ;

    std::vector<double>                     d;
    std::vector<double>::iterator           vit;

    std::vector< std::array<double,3> >     xP ;
    std::vector< int >                      where ;
    std::vector<int>                        flag( mesh.getCellCount(), -1);

    int                                     i, N( stl->getCellCount() );
    double                                  search( visitee->RSearch ) ;

    std::vector<long>                       neighs ;

    PiercedVector<SegData>::iterator data ;


    stack.reserve(128) ;
    temp.reserve(128) ;
    cloud.reserve(128) ;

    xP.resize(1) ;
    where.resize(1) ;

    // --------------------------------------------------------------------------
    // COMPUTE THE SDF VALUE AT EACH MESH POINT                                   //
    //
    for (i = 0; i < N; i++) {

        // Segments vertex ------------------------------------------------------ //
        VS  = getSimplexVertices( i ) ;
        seedNarrowBand( visitee, VS, stack ) ;


        //-----------------------------------------------------------------
        //Propagate from seed
        while (stack.size() > 0) {

            // Extract point from lifo
            cloud.clear() ;
            for( const auto & I : stack){
                cloud.push_back( mesh.evalCellCentroid(I) ) ;  
            };

            d = CGElem::distanceCloudSimplex( cloud, VS, xP, where); 
            vit = d.begin() ;

            for( const auto & I : stack){
                if ( *vit <= search ) {

                    if( m_segInfo.exists(I) ){
                        m_segInfo[I].m_segments.insert(i) ;
                    } else {
                        data = m_segInfo.reclaim(I) ;
                        data->m_segments.insert(i);
                    };


                    neighs  = mesh.findCellFaceNeighs(I) ; 
                    for( const auto &  neigh : neighs){
                        if( flag[neigh] != i) {
                            temp.push_back( neigh) ;
                            flag[neigh] = i ;
                        };
                    }

                } //end if distance

                ++vit ;

            };

            stack.clear() ;
            stack.swap( temp ) ;


        } //end while


    } //end for i

    return;

};

/*!
 * Compute the levelset function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. 
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetCartesian *visitee, std::vector<Adaption::Info> &mapper, double &newRSearch ){

    updateSimplexToCell(visitee, mapper, newRSearch ) ;
    lsFromSimplex(visitee,  newRSearch, true) ;

    return;
};

/*!
 * Update the Sdf of the triangulation after an octree adaptation.
 * Note: Only a single octree adapt with marker (-1,0,1) is permitted.
 */
void LevelSetSegmentation::updateSimplexToCell( LevelSetCartesian *visitee, std::vector<Adaption::Info> &mapper, double &newRSearch ){

    BITPIT_UNUSED(mapper) ;
    BITPIT_UNUSED(newRSearch) ;

    associateSimplexToCell( visitee ) ;

    return;
};

/*!
 * Compute the levelset function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. 
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetOctree *visitee ){

    associateSimplexToCell(visitee) ;
    lsFromSimplex(visitee, visitee->RSearch,true) ;

    return;
};

/*!
 * Compute the signed distance function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. If usersearch = true, RSearch
 * is initialized around the geometry (in order to guarantee at least 2 mesh cell centers
 * around the geometry) and used as imposed size of narrow band.
 */
void LevelSetSegmentation::associateSimplexToCell( LevelSetOctree *visitee){

    VolumeKernel         &mesh = *visitee->m_mesh ;
    int                         dim(mesh.getDimension()) ;
    long                        id, icart;
    int                         i;
    double                      size ;

    std::array<double,3>        C, C0, C1, G0, G1, octrBB0, octrBB1, triBB0, triBB1 ;
    PiercedVector<SegData>::iterator data;

    { // mesh size corresponding to RSearch
        uint8_t    level = visitee->computeLevelFromRSearch( visitee->RSearch ) ;
        size = (visitee->m_omesh->getTree()).levelToSize(level);
    }

    mesh.getBoundingBox(octrBB0,octrBB1) ;
    getBoundingBox( triBB0, triBB1 );

    G0 = triBB0 - visitee->RSearch ;
    G1 = triBB1 + visitee->RSearch ;

    if( CGElem::intersectBoxBox(octrBB0,octrBB1,G0,G1,C0,C1) ) { //intersect two Bounding Boxes around geometry and local grid

        // snap bounding box to grid and create cartesian grid
        std::array<int,3>    nc ;

        for( i=0; i<dim; ++i){
            C0[i] =  octrBB0[i] + size *   (int) ( ( C0[i] - octrBB0[i] ) / size ) ;
            C1[i] =  octrBB0[i] + size * ( (int) ( ( C1[i] - octrBB0[i] ) / size ) +1 ) ;

            nc[i] = round( ( C1[i] - C0[i] ) /size ) ;
        };

        // calculate LS triangle lists on cartesian mesh and map on pablo
        VolCartesian            cmesh( 0, dim, C0, C1-C0, nc) ;

        LevelSetCartesian       auxLS(cmesh) ;
        LevelSetSegmentation    objLS(*this) ;

        auxLS.computeSizeNarrowBand(this) ;
        objLS.associateSimplexToCell(&auxLS) ; 

        for( auto & cell : mesh.getCells() ){
            id = cell.getId() ;
            C  = mesh.evalCellCentroid(id) ;

            if( cmesh.isPointInside(C) ){
                icart = cmesh.locatePoint(C) ;

                if( objLS.isInNarrowBand(icart) ){
                    const std::set<long> &list = objLS.getSimplexList(icart) ;
                    data = m_segInfo.emplace(id, list) ;
                };

            };

        };

    }; //endif intersect

    return;

};

/*!
 * Compute the levelset function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. 
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetOctree *visitee, std::vector<Adaption::Info> &mapper, double &newRSearch ){

    updateSimplexToCell(visitee, mapper, newRSearch ) ; 

    visitee->clearAfterRefinement(mapper) ;

    lsFromSimplex(visitee, newRSearch, true) ;

    return;
};

/*!
 * Update the Sdf of the triangulation after an octree adaptation.
 * Note: Only a single octree adapt with marker (-1,0,1) is permitted.
 */
void LevelSetSegmentation::updateSimplexToCell( LevelSetOctree *visitee, std::vector<Adaption::Info> &mapper, double &newRSearch){

    int         oldLevel, newLevel ;

    oldLevel = visitee->computeLevelFromRSearch( visitee->RSearch ) ;
    newLevel = visitee->computeLevelFromRSearch( newRSearch ) ;

    if( newLevel <= oldLevel ) { //size of narrow band decreased or remained the same -> mapping

        std::unordered_map<long,std::set<long>> oldSegs ;
        std::unordered_map<long,std::set<long>>::iterator oldSegsIt ;

        for ( auto & info : mapper ){
            if( info.entity == Adaption::Entity::ENTITY_CELL ){

                for ( auto & parent : info.previous){ //save old data and delete element
                    if( m_segInfo.exists(parent) ){
                        SegData *seg =  &m_segInfo[parent] ;

                        oldSegs.insert({{ parent, seg->m_segments }}) ;
                        m_segInfo.erase(parent,true) ;
                    }
                }
            }
        }

        m_segInfo.flush() ;

        for ( auto & info : mapper ){ //forall mesh modifications
            if( info.entity == Adaption::Entity::ENTITY_CELL){ //check if changes on cells
                for ( auto & child : info.current){ // forall new elements

                    PiercedVector<SegData>::iterator seg =  m_segInfo.reclaim(child) ;
                    seg->m_segments.clear() ;

                    for ( auto & parent : info.previous){ //take their parents
                        oldSegsIt = oldSegs.find(parent);
                        if( oldSegsIt != oldSegs.end() ) //add their information if any
                            seg->m_segments.insert( oldSegsIt->second.begin(), oldSegsIt->second.end() ) ;
                    }

                }
            }
        }


    }


    else { //size of narrow band increased -> recalculation

        m_segInfo.clear() ;
        associateSimplexToCell(visitee) ; 

    };


    return ;

};

}
