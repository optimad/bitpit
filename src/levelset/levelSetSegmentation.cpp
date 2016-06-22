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

# include <cassert>

# include "levelSet.hpp"

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

/*!
	@ingroup    levelset
	@class      LevelSetSegmentation::SegInfo
	@brief      Implements visitor pattern fo segmentated geometries
*/

/*!
 * Default constructor 
 */
LevelSetSegmentation::SegInfo::SegInfo( ) : m_segments(levelSetDefaults::LIST), m_checked(false){
};

/*!
 * Constructor
 * @param[in] list list of simplices
 */
LevelSetSegmentation::SegInfo::SegInfo( const std::unordered_set<long> &list) :m_segments(list), m_checked(false) {
};

/*!
	@ingroup    levelset
	@class      LevelSetSegmentation
	@brief      Implements visitor pattern fo segmentated geometries
*/

/*!
 * Destructor
 */
LevelSetSegmentation::~LevelSetSegmentation() {
    m_segmentation = NULL;
    m_vertexNormal.clear() ;
};

/*!
 * Constructor
 * @param[in] id identifier of object
 */
LevelSetSegmentation::LevelSetSegmentation( int id) :LevelSetObject(id) {

};

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL unique pointer to surface mesh
 */
LevelSetSegmentation::LevelSetSegmentation( int id, std::unique_ptr<SurfUnstructured> &&STL) :LevelSetObject(id) {

    setSegmentation( std::move(STL) );

};

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL pointer to surface mesh
 */
LevelSetSegmentation::LevelSetSegmentation( int id, SurfUnstructured *STL) :LevelSetObject(id) {

    setSegmentation( STL );

};

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetSegmentation::LevelSetSegmentation( const LevelSetSegmentation &other) : LevelSetObject(other.getId() ) {

    m_segmentation = other.m_segmentation; 
    m_dimension = other.m_dimension ;
    m_vertexNormal = other.m_vertexNormal ;
    if (m_own != nullptr) {
        m_own = unique_ptr<SurfUnstructured>(new SurfUnstructured(*(other.m_own)));
    }
};

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetSegmentation* LevelSetSegmentation::clone() const {
    return new LevelSetSegmentation( *this ); 
}

/*!
 * Set the segmentation
 * @param[in] segmentation unique pointer to surface mesh
 */
void LevelSetSegmentation::setSegmentation( std::unique_ptr<SurfUnstructured> &&segmentation){

    m_own = std::move(segmentation) ;

    setSegmentation( m_own.get() );
}

/*!
 * Set the segmentation
 * @param[in] segmentation pointer to surface mesh
 */
void LevelSetSegmentation::setSegmentation( SurfUnstructured *segmentation){

    std::vector<std::array<double,3>>   vertexNormal ;

    m_segmentation = segmentation;
    m_dimension = m_segmentation->getSpaceDimension() ;

    int  i, nV;
    long segId ;

    for( auto & segment : m_segmentation->getCells() ){

        segId = segment.getId() ;
        nV = segment.getVertexCount() ;

        vertexNormal.resize(nV) ;
        for(i=0; i<nV; ++i){
            vertexNormal[i] = m_segmentation->evalVertexNormal(segId,i) ;
        }

        m_vertexNormal.insert({{segId,vertexNormal}}) ;

    };
}

/*!
 * Get a constant refernce to the segmentation
 * @return constant reference to the segmentation
 */
const SurfUnstructured & LevelSetSegmentation::getSegmentation() const {
    return *m_segmentation ;
};

/*!
 * Get the list of simplices wich contain the cell centroid in their narrow band.
 * @param[in] i cell index
 * @return set with indices of simplices
 */
const std::unordered_set<long> & LevelSetSegmentation::getSimplexList(const long &i) const{

    if( !m_seg.exists(i) ){
        return levelSetDefaults::LIST;
    } else {
        return ( m_seg[i].m_segments );
    };

};

/*!
 * Check if cell is in narrowband of any triangle;
 * @param[in] i cell index
 * @return if in narow band
 */
bool LevelSetSegmentation::isInNarrowBand(const long &i){
    return( m_seg.exists(i) ) ; 
};

/*!
 * Aggregate cell vertex coordinates in one vector
 * @param[in] i cell index
 * @return coordinates of cell vertices
 */
std::vector<std::array<double,3>> LevelSetSegmentation::getSimplexVertices( const long &i ) const {

    Cell &cell = m_segmentation->getCell(i) ;

    int                                     j, n, N (cell.getVertexCount()) ;
    std::vector<std::array<double,3>>       VS(N) ;

    for( n=0; n<N; ++n){
        j = cell.getVertex(n) ;
        VS[n] = m_segmentation->getVertexCoords(j);
    };

    if( N > 3){
        log::cout() << "levelset: only segments and triangles supported in LevelSetSegmentation !!" << std::endl ;
    }

    return VS;
};

/*!
 * Update the levelset function of whole mesh by using associated simplices
 * @param[in] visitee visited mesh 
 * @param[in] search size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 * @param[in] filter if triangles should be ereased when outside narrow band (default false)
 */
void LevelSetSegmentation::lsFromSimplex( LevelSetKernel *visitee, const double &search, const bool & signd, bool filter){

    VolumeKernel                &mesh  = *(visitee->getMesh() ) ;

    long                        id ;
    double                      s, d, value;
    std::array<double,3>        n, xP, P;

    std::unordered_set<long>::iterator  it, itend ;
    PiercedIterator<SegInfo>            segIt ;
    PiercedVector<LevelSetInfo>         &lsInfo = visitee->getLevelSetInfo() ;

    for( segIt=m_seg.begin(); segIt!=m_seg.end(); ++segIt ){

        SegInfo                 &segInfo = *segIt ;

        if( segInfo.m_checked == false){
            segInfo.m_checked = true ;

            std::unordered_set<long>    &segs = segInfo.m_segments ;

            id    = segIt.getId() ;

            it    = segs.begin();
            itend = segs.end() ;

            P = mesh.evalCellCentroid(id) ;

            auto lsInfoItr = lsInfo.find(id) ;
            if( lsInfoItr != lsInfo.end() ){
                value = abs( lsInfoItr->value );
            } else {
                value = 1e18;
            }

            while( it != itend ){

                infoFromSimplex(P, *it, d, s, xP, n);

                if ( d <= search ){

                    if( d<value ) {
                        if (lsInfoItr == lsInfo.end()) {
                            lsInfoItr = lsInfo.reclaim(id) ;
                        }

                        value   = d ;

                        lsInfoItr->object   = getId();
                        lsInfoItr->part     = m_segmentation->getCell(*it).getPID();
                        lsInfoItr->support  = *it ;
                        lsInfoItr->value    = ( signd *s  + (!signd) *1.) *d; 
                        lsInfoItr->gradient = ( signd *1. + (!signd) *s ) *n ;
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

            if( segs.size() == 0 ){
                m_seg.erase(id,true) ;
            };

        };

    };// foreach cell

    m_seg.flush() ;

    return;

};

/*!
 * Computes levelset relevant information at one point with respect to a simplex
 * @param[in] p coordinates of point
 * @param[in] i index of simplex
 * @param[out] d distance point to simplex
 * @param[out] s sign of point wrt to simplex, i.e. according to normal
 * @param[out] x closest point on simplex
 * @param[out] n normal at closest point
 */
void LevelSetSegmentation::infoFromSimplex( const std::array<double,3> &p, const long &i, double &d, double &s, std::array<double,3> &x, std::array<double,3> &n ) const {

    Cell &cell = m_segmentation->getCell(i) ;
    int nV = cell.getVertexCount() ;

    auto itr = m_vertexNormal.find(i) ;
    assert( itr != m_vertexNormal.end() ) ;

    if( nV == 1){
        long id = cell.getVertex(0) ;
        d = norm2( p- m_segmentation->getVertexCoords(id) ) ;
        n.fill(0.) ;

    } else if( nV == 2){
        long id0 = cell.getVertex(0) ;
        long id1 = cell.getVertex(1) ;

        std::array<double,2> lambda ;

        d= CGElem::distancePointSegment( p, m_segmentation->getVertexCoords(id0), m_segmentation->getVertexCoords(id1), x, lambda ) ;
        n  = lambda[0] *itr->second[0] ;
        n += lambda[1] *itr->second[1] ;

        n /= norm2(n) ;

    } else if (nV == 3){
        long id0 = cell.getVertex(0) ;
        long id1 = cell.getVertex(1) ;
        long id2 = cell.getVertex(2) ;

        std::array<double,3> lambda ;

        d= CGElem::distancePointTriangle( p, m_segmentation->getVertexCoords(id0), m_segmentation->getVertexCoords(id1), m_segmentation->getVertexCoords(id2), x, lambda ) ;
        n  = lambda[0] *itr->second[0] ;
        n += lambda[1] *itr->second[1] ;
        n += lambda[2] *itr->second[2] ;

        n /= norm2(n) ;

    } else{
        log::cout() << " simplex not supported in LevelSetSegmentation::infoFromSimplex " << nV << std::endl ;
        
    };

    s = sign( dotProduct(n, p - x) );


    return ;

};

/*!
 * Finds seed points in narrow band within a cartesian mesh for one simplex
 * @param[in] visitee cartesian mesh 
 * @param[in] VS Simplex
 * @param[out] I indices of seed points
 */
bool LevelSetSegmentation::seedNarrowBand( LevelSetCartesian *visitee, std::vector<std::array<double,3>> &VS, std::vector<int> &I){

    VolumeKernel                        &mesh = *(visitee->getMesh()) ;

    bool                                found(false) ;
    int                                 dim( mesh.getDimension() ) ;
    std::array<double,3>                B0, B1;
    std::vector<std::array<double,3>>   VP ;

    mesh.getBoundingBox(B0, B1) ;

    B0 = B0 - visitee->getSizeNarrowBand() ;
    B1 = B1 + visitee->getSizeNarrowBand() ;

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
 * Evaluates the levelset in the specified cell
 * @param[in] visitee pointer to mesh
 * @param[in] id is the cell index
 * @result The value of the levelset.
 */
double LevelSetSegmentation::evaluateLS( LevelSetKernel *visitee, long id) const {

    double                      d, s, ls;
    std::array<double,3>        P, X, temp;

    P  = visitee->getMesh()->evalCellCentroid(id) ;
    ls = levelSetDefaults::VALUE;

    for( auto & segment : m_segmentation->getCells() ){

        infoFromSimplex(P, segment.getId(), d, s, X, temp);

        d = abs(d) ;
        if ( d < ls && !utils::DoubleFloatingEqual()(s, (double) 0.)) {
            ls = d ;
        }

    }

    return ls;

};

/*!
 * Computes axis aligned bounding box of object
 * @param[out] minP minimum point
 * @param[out] maxP maximum point
 */
void LevelSetSegmentation::getBoundingBox( std::array<double,3> &minP, std::array<double,3> &maxP ) const {
    m_segmentation->getBoundingBox(minP,maxP) ;
};

/*!
 * Computes the levelset function within the narrow band
 * @param[in] visitee pointer to mesh
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetKernel *visitee, const double &RSearch, const bool &signd ){

    if( LevelSetCartesian* lsCartesian = dynamic_cast<LevelSetCartesian*>(visitee) ){
        associateSimplexToCell( lsCartesian, RSearch ) ;

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(visitee) ){
        associateSimplexToCell( lsOctree, RSearch ) ;

    };

    lsFromSimplex(visitee, RSearch, signd, true) ;

    return;
};

/*!
 * Updates the levelset function within the narrow band after mesh adaptation.
 * @param[in] visitee pointer to mesh
 * @param[in] mapper information concerning mesh adaption 
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetKernel *visitee, const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd ){

    if( LevelSetCartesian* lsCartesian = dynamic_cast<LevelSetCartesian*>(visitee) ){
        associateSimplexToCell( lsCartesian, RSearch ) ;

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(visitee) ){
        updateSimplexToCell( lsOctree, mapper, RSearch ) ; 

    };

    lsFromSimplex(visitee, RSearch, signd, true) ;

    return;
};

/*!
 * Determines the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle) for cartesian meshes
 * @param[in] visitee pointer to cartesian mesh
 * @param[in] RSearch size of narrow band
 */
void LevelSetSegmentation::associateSimplexToCell( LevelSetCartesian *visitee, const double &RSearch ){

    VolumeKernel                            &mesh = *(visitee->getMesh() ) ;
    std::vector<std::array<double,3>>       VS(3);

    std::vector< int >                      stack, temp ;
    std::vector< std::array<double,3> >     cloud ;

    std::vector<double>                     d;
    std::vector<double>::iterator           vit;

    std::vector< std::array<double,3> >     xP ;
    std::vector< int >                      where ;
    std::vector<int>                        flag( mesh.getCellCount(), -1);

    int                                     i, N( m_segmentation->getCellCount() );

    std::vector<long>                       neighs ;

    PiercedVector<SegInfo>::iterator data ;


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
        seedNarrowBand( visitee, VS, stack ) ; //TODO check if seed is found correctly if segmentation is outside grid but within narrow band


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
                if ( *vit <= RSearch ) {

                    if( m_seg.exists(I) ){
                        m_seg[I].m_segments.insert(i) ;
                    } else {
                        data = m_seg.emplace(I) ;
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
 * Determines the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle) for octree meshes
 * @param[in] visitee pointer to octree mesh
 * @param[in] RSearch size of narrow band
 */
void LevelSetSegmentation::associateSimplexToCell( LevelSetOctree *visitee, const double &RSearch){

    VolumeKernel                &mesh = *(visitee->getMesh()) ;
    int                         dim(mesh.getDimension()) ;
    long                        id, icart;
    int                         i;
    double                      size ;

    std::array<double,3>        C, C0, C1, G0, G1, octrBB0, octrBB1, triBB0, triBB1 ;
    PiercedVector<SegInfo>::iterator data;

    // mesh size corresponding to RSearch
    size = visitee->computeSizeFromRSearch( RSearch ) ;

    mesh.getBoundingBox(octrBB0,octrBB1) ;
    getBoundingBox( triBB0, triBB1 );

    G0 = triBB0 - RSearch ;
    G1 = triBB1 + RSearch ;

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

        double                  localRSearch = auxLS.computeSizeNarrowBand(this) ;

        auxLS.setSizeNarrowBand(localRSearch);
        objLS.associateSimplexToCell(&auxLS, auxLS.getSizeNarrowBand() ) ; 

        for( auto & cell : mesh.getCells() ){
            id = cell.getId() ;
            C  = mesh.evalCellCentroid(id) ;

            if( cmesh.isPointInside(C) ){
                icart = cmesh.locatePoint(C) ;

                if( objLS.isInNarrowBand(icart) ){
                    const std::unordered_set<long> &list = objLS.getSimplexList(icart) ;
                    data = m_seg.emplace(id, list) ;
                };

            };

        };

    }; //endif intersect

    return;

};

/*!
 * Updates the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle) for octree meshes
 * @param[in] visitee pointer to octree mesh
 * @param[in] mapper information concerning mesh adaption 
 * @param[in] newRSearch new size of narrow band
 */
void LevelSetSegmentation::updateSimplexToCell( LevelSetOctree *visitee, const std::vector<adaption::Info> &mapper, const double &newRSearch){

    double      oldSize, newSize ;

    oldSize = visitee->computeSizeFromRSearch( visitee->getSizeNarrowBand() ) ;
    newSize = visitee->computeSizeFromRSearch( newRSearch ) ;

    if( newSize-oldSize <= 1.e-8 ) { //size of narrow band decreased or remained the same -> mapping

        std::vector<std::vector<long>> previousSegments ;

        long nNewElements = 0;
        for ( auto & info : mapper ){
            // Consider only changes on cells
            if( info.entity != adaption::Entity::ENTITY_CELL ){
                continue;
            }

            // Count new elements
            long nCurrentElements = info.current.size() ;
            nNewElements += nCurrentElements ;

            // Save previous data and delete previos elements
            std::vector<long> *parentSegments = nullptr ;
            if ( nCurrentElements > 0 ) {
                previousSegments.emplace_back() ;
                parentSegments = &previousSegments.back() ;
            }

            for ( auto & parent : info.previous){

                PiercedVector<SegInfo>::const_iterator parentSegInfo = m_seg.find(parent) ;
                if( parentSegInfo != m_seg.cend() ){
                    // Add previous segments only if there are current elements
                    // associated with this change.
                    if( parentSegments != nullptr ){
                        auto &segments = parentSegInfo->m_segments ;
                        parentSegments->insert( parentSegments->end(), segments.begin(), segments.end() ) ;
                    }

                    m_seg.erase(parent,true) ;
                }
            }

            // Remove duplicate entries from the list of previous segments
            if( parentSegments != nullptr ){
                std::sort( parentSegments->begin(), parentSegments->end() ) ;
                parentSegments->erase( std::unique(parentSegments->begin(), parentSegments->end()), parentSegments->end() ) ;
            }
        }

        m_seg.flush() ;

        // Update new elements
        if (nNewElements > 0) {
            m_seg.reserve(m_seg.size() + nNewElements) ;

            size_t adaptionIdx = 0;
            for ( auto & info : mapper ){
                // Consider only changes on cells
                if( info.entity != adaption::Entity::ENTITY_CELL ){
                    continue;
                }

                // Skip adaption info with no current elements
                long nCurrentElements = info.current.size() ;
                if (nCurrentElements == 0) {
                    continue;
                }

                // Get the list of parent segments
                std::vector<long> &parentSegments = previousSegments[adaptionIdx] ;

                // Assign the segments of the parents to the childs
                for ( auto & child : info.current){
                    PiercedVector<SegInfo>::iterator childSegInfo = m_seg.emplace(child) ;
                    childSegInfo->m_checked = false;
                    childSegInfo->m_segments.insert( parentSegments.begin(), parentSegments.end() ) ;
                }

                // Increase adaption info counter
                adaptionIdx++;
            }
        }

    } else { //size of narrow band increased -> recalculation

        m_seg.clear() ;
        visitee->clear() ;
        associateSimplexToCell( visitee, visitee->getSizeNarrowBand() ) ; 

    };

};

/*! 
 * Deletes non-existing items 
 * @param[in] mapper mapping info
 */
void LevelSetSegmentation::clearAfterMeshMovement( const std::vector<adaption::Info> &mapper ){

    for ( auto & map : mapper ){
        if( map.entity == adaption::Entity::ENTITY_CELL ){
            if( map.type == adaption::Type::TYPE_DELETION || 
                map.type == adaption::Type::TYPE_PARTITION_SEND  ||
                map.type == adaption::Type::TYPE_REFINEMENT  ||
                map.type == adaption::Type::TYPE_COARSENING  ){

                for ( auto & parent : map.previous){
                    if( m_seg.exists(parent) ) 
                        m_seg.erase(parent,true) ;
                }
            }
        }
    }

    m_seg.flush() ;

    return ;
};

/*!
 * Clears data structure outside narrow band
 * @param[in] visitee LevelSetKernel with narrow band information
 */
void LevelSetSegmentation::filterOutsideNarrowBand( LevelSetKernel *visitee ){

    long id ;

    bitpit::PiercedVector<SegInfo>::iterator segItr ;
    for( segItr = m_seg.begin(); segItr != m_seg.end(); ++segItr){
        id = segItr.getId() ;
        if( ! visitee->isInNarrowBand(id) ){
            m_seg.erase(id,true) ;
        };
    };

    m_seg.flush() ;

    return ;
};

/*!
 * Writes LevelSetSegmentation to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetSegmentation::dumpDerived( std::fstream &stream ){

    int                 s;
    std::vector<long>   temp;

    bitpit::PiercedVector<SegInfo>::iterator segItr, segEnd = m_seg.end() ;

    bitpit::genericIO::flushBINARY( stream, (long) m_seg.size() ) ;

    for( segItr = m_seg.begin(); segItr != segEnd; ++segItr){
        s = segItr->m_segments.size() ;

        temp.resize(s);
        std::copy( segItr->m_segments.begin(), segItr->m_segments.end(), temp.begin() );

        bitpit::genericIO::flushBINARY( stream, segItr.getId() );
        bitpit::genericIO::flushBINARY( stream, s );
        bitpit::genericIO::flushBINARY( stream, temp );
        bitpit::genericIO::flushBINARY( stream, segItr->m_checked );
    }

    return;
};

/*!
 * Reads LevelSetSegmentation from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetSegmentation::restoreDerived( std::fstream &stream ){

    int     s;
    long    i, n, id;
    SegInfo cellData ;
    std::vector<long>   temp;

    bitpit::genericIO::absorbBINARY( stream, n ) ;

    m_seg.reserve(n);

    for( i=0; i<n; ++i){
        bitpit::genericIO::absorbBINARY( stream, id );
        bitpit::genericIO::absorbBINARY( stream, s );

        temp.resize(s) ;
        bitpit::genericIO::absorbBINARY( stream, temp );
        bitpit::genericIO::absorbBINARY( stream, cellData.m_checked );

        std::copy( temp.begin(), temp.end(), std::inserter( cellData.m_segments, cellData.m_segments.end() ) );

        m_seg.insert(id,cellData) ;

    }

    return;
};

# if BITPIT_ENABLE_MPI

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] sizeBuffer buffer for first communication used to communicate the size of data buffer
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetSegmentation::writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &sizeBuffer, SendBuffer &dataBuffer ){

    long nItems = sendList.size(), counter(0) ;
    int dataSize = 10*sizeof(long)  +sizeof(bool) +sizeof(long) +sizeof(int) ;

    dataBuffer.setCapacity(nItems*dataSize) ;

    //determine elements to send
    nItems = 0 ;
    for( const auto &index : sendList){
        if( m_seg.exists(index)){
            const auto &seginfo = m_seg[index] ;
            dataBuffer << counter ;
            dataBuffer << (int) seginfo.m_segments.size() ;
            for( const long & seg : seginfo.m_segments ){
                dataBuffer << seg ;
            };
            dataBuffer << seginfo.m_checked ;
            ++nItems ;
        }

        ++counter ;
    }


    dataBuffer.squeeze() ;
    sizeBuffer << nItems ;

    return;
};

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in] nItems number of items within the buffer
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetSegmentation::readCommunicationBuffer( const std::vector<long> &recvList, const long &nItems, RecvBuffer &dataBuffer ){

    int     s, nSegs ;
    long    index, id, segment ;

    for( int i=0; i<nItems; ++i){
        // Get the id of the element
        dataBuffer >> index ;
        id = recvList[index] ;

        // Assign the data of the element
        PiercedVector<SegInfo>::iterator segItr ;
        if( !m_seg.exists(id)){
            segItr = m_seg.emplace(id) ;
        } else {
            segItr = m_seg.getIterator(id) ;
        }

        dataBuffer >> nSegs ;
        for( s=0; s<nSegs; ++s){
            dataBuffer >> segment ;
            segItr->m_segments.insert(segment) ;
        }

        dataBuffer >> segItr->m_checked ;
    }

    return;
};
# endif

}
