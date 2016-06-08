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

/*!
	@ingroup    levelset
	@class      LevelSetObject
	@brief      LevelSet object.
*/

/*!
 * Destructor
 */
LevelSetObject::~LevelSetObject( ){
};

/*!
 * Constructor
 * @param[in] id Id assigned to object
 */
LevelSetObject::LevelSetObject( int id) : m_id(id){
};

/*!
 * Get the id 
 * @return Id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
};

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::fstream &stream ){
    bitpit::genericIO::flushBINARY(stream, m_id) ;
    dumpDerived(stream) ;
    return;
};

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::fstream &stream ){
    bitpit::genericIO::absorbBINARY(stream, m_id) ;
    restoreDerived(stream) ;
    return;
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
LevelSetSegmentation::SegInfo::SegInfo( ) : m_segments(levelSetDefaults::LIST), m_support(levelSetDefaults::ELEMENT), m_checked(false){
};

/*!
 * Constructor of LevelSet_Stl with input parameters.
 * @param[in] id id to be asigned to pierced vector
 * @param[in] list list of simplices
 */
LevelSetSegmentation::SegInfo::SegInfo( const std::set<long> &list) :m_segments(list), m_support(levelSetDefaults::ELEMENT), m_checked(false) {
};

/*!
 * Constructor of LevelSet_Stl with input parameters.
 * @param[in] id id to be asigned to pierced vector
 * @param[in] list list of simplices
 * @param[in] support index of closest simplex
 */
LevelSetSegmentation::SegInfo::SegInfo( const std::set<long> &list, const long &support) :m_segments(list), m_support(support), m_checked(false){
};

/*!
 * Destructor of LevelSet_Stl.
 */
LevelSetSegmentation::~LevelSetSegmentation() {
    m_segmentation = NULL;
};

/*!
 * Constructor of LevelSet_Stl with input parameters.
 * It builds one Sdf object with :
 * @param[in] *Pmesh Pointer to 3D pablo octree mesh to be linked to Sdf.
 * @param[in] *STL Pointer to surface Triangulation to be linked to Sdf.
 */
LevelSetSegmentation::LevelSetSegmentation( int id, SurfUnstructured *STL) :LevelSetObject(id) {

    m_segmentation = STL;
    m_dimension = m_segmentation->getSpaceDimension() ;

};

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetSegmentation::LevelSetSegmentation( const LevelSetSegmentation &other) :LevelSetObject(other.getId() ) {

    m_segmentation = other.m_segmentation; 
    m_dimension = other.m_dimension ;

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

    if( !m_seg.exists(i) ){
        return levelSetDefaults::LIST;
    } else {
        return ( m_seg[i].m_segments );
    };

};

/*!
 * Get the index of the closest simplex
 * @param[in] i Local index of target octant.
 * @return index of closest simplex
 */
const long & LevelSetSegmentation::getSupportSimplex(const long &i){

    if( !m_seg.exists(i) ){
        return levelSetDefaults::ELEMENT ;
    } else {
        return ( m_seg[i].m_support );
    };
};

/*!
 * Check if cell is in narrowband of any triangle;
 * LevelSetSegmentation::associateSimplexToCell() should have been called prior;
 * @param[in] i Local index of target octant.
 * @return if in narow band
 */
bool LevelSetSegmentation::isInNarrowBand(const long &i){
    return( m_seg.exists(i) ) ; 
};

/*!
 * Add a SurfTriPatch object to stl.
 * @param[in] *STL Pointer to surface Triangulation to be added .
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
 * Update the signed distance function for a given octant with the list of simplices
 * which contain the octant in the narrow band (if et[IDX] is not void) otherwise
 * by using the full triangulation.
 * @param[in] visitee visited mesh 
 * @param[in] search size of narrow band
 * @param[in] filter if triangles should be ereased when outside narrow band (default false)
 */
void LevelSetSegmentation::lsFromSimplex( LevelSetKernel *visitee, const double &search, const bool & signd, bool filter){

    VolumeKernel                &mesh  = *(visitee->getMesh() ) ;

    long                        id ;
    double                      s, d, value;
    std::array<double,3>        n, xP, P;

    std::set<long>::iterator    it, itend ;
    PiercedIterator<SegInfo>    segIt ;
    PiercedVector<LevelSetKernel::LSInfo>       &lsInfo = visitee->getLSInfo() ;

    for( segIt=m_seg.begin(); segIt!=m_seg.end(); ++segIt ){

        SegInfo                 &segInfo = *segIt ;

        if( segInfo.m_checked == false){
            segInfo.m_checked = true ;

            std::set<long>      &segs = segInfo.m_segments ;
            long                &supp = segInfo.m_support ;

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
                        lsInfoItr->value    = ( signd *s  + (!signd) *1.) *d; 
                        lsInfoItr->gradient = ( signd *1. + (!signd) *s ) *n ;
                        supp                = *it ;
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

    //TODO save vertex normals
    if( m_dimension == 2){
        std::array<double,2> lambda ;

        d= CGElem::distancePointSegment( p, VS[0], VS[1], x, lambda ) ;
        n  = lambda[0] *m_segmentation->evalVertexNormal(i,0) ;
        n += lambda[1] *m_segmentation->evalVertexNormal(i,1) ;

    } else {
        std::array<double,3> lambda ;

        d= CGElem::distancePointTriangle( p, VS[0], VS[1], VS[2], x, lambda ) ;
        n  = lambda[0] *m_segmentation->evalVertexNormal(i,0) ;
        n += lambda[1] *m_segmentation->evalVertexNormal(i,1) ;
        n += lambda[2] *m_segmentation->evalVertexNormal(i,2) ;
    };

    s = sign( dotProduct(n, p - x) );


    return ;

};

/*!
 * Finds seed points for given simplex
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
 * Finds a seed point to be used for propgating the levelset sign
 * @param[out] seed index of seed cell
 * @param[out] sign sign at seed cell [+1;-1]
 */
void LevelSetSegmentation::seedSign( LevelSetKernel *visitee, long &seed, double &value) const {

    seed = -1;

    {  // FIND SEED ELEMENT In GRID
        double                      lsTmp;
        const PiercedVector<LevelSetKernel::LSInfo> & lsInfo = visitee->getLSInfo() ;
        PiercedVector<LevelSetKernel::LSInfo>::const_iterator lsItr = lsInfo.begin() , lsEnd=lsInfo.end() ;

        while( lsItr != lsEnd ) {
            lsTmp = (*lsItr).value ;
            if( std::abs(lsTmp) < levelSetDefaults::VALUE ){
                seed = lsItr.getId();
                value = sign(lsTmp) ;
                return ;
            };

            ++lsItr ;
        };
    }

    { // if no seed element found, check whole triangulation

        double                      d, s, lsTmp;
        std::array<double,3>        P, X, temp;

        PiercedIterator<bitpit::Cell>   cellItr = (visitee->getMesh())->getCells().begin() ;

        seed    = (*cellItr).getId() ;
        P       = (visitee->getMesh())->evalCellCentroid(seed) ;


        lsTmp  = levelSetDefaults::VALUE;

        for( auto & segment : m_segmentation->getCells() ){

            infoFromSimplex(P, segment.getId(), d, s, X, temp);

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
    m_segmentation->getBoundingBox(minP,maxP) ;
};

/*!
 * Compute the levelset function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. 
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetKernel *visitee, const double &RSearch, const bool &signd ){

    if( LevelSetCartesian* lsCartesian = dynamic_cast<LevelSetCartesian*>(visitee) ){
        associateSimplexToCell( lsCartesian, RSearch ) ;

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(visitee) ){
        associateSimplexToCell( lsOctree, RSearch ) ;

    };

    lsFromSimplex(visitee, RSearch, signd) ;

    return;
};

/*!
 * Compute the levelset function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. 
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetKernel *visitee, const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd ){

    if( LevelSetCartesian* lsCartesian = dynamic_cast<LevelSetCartesian*>(visitee) ){
        associateSimplexToCell( lsCartesian, RSearch ) ;

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(visitee) ){
        updateSimplexToCell( lsOctree, mapper, RSearch ) ; 

    };

    lsFromSimplex(visitee, RSearch, signd) ;

    return;
};

/*!
 * Determines the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle).
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
                        data = m_seg.reclaim(I) ;
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
 * Compute the signed distance function from a piece-wise linear approximation
 * of a d manifold in a 3D Euclidean space. Level set is computed in narrow band
 * of at least 2 mesh cell centers around the geometry. If usersearch = true, RSearch
 * is initialized around the geometry (in order to guarantee at least 2 mesh cell centers
 * around the geometry) and used as imposed size of narrow band.
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
                    const std::set<long> &list = objLS.getSimplexList(icart) ;
                    data = m_seg.emplace(id, list) ;
                };

            };

        };

    }; //endif intersect

    return;

};

/*!
 * Update the Sdf of the triangulation after an octree adaptation.
 * Note: Only a single octree adapt with marker (-1,0,1) is permitted.
 */
void LevelSetSegmentation::updateSimplexToCell( LevelSetOctree *visitee, const std::vector<adaption::Info> &mapper, const double &newRSearch){

    double      oldSize, newSize ;

    oldSize = visitee->computeSizeFromRSearch( visitee->getSizeNarrowBand() ) ;
    newSize = visitee->computeSizeFromRSearch( newRSearch ) ;

    //TODO does not support coarsening if children were on different processes

    if( newSize-oldSize <= 1.e-8 ) { //size of narrow band decreased or remained the same -> mapping

        { // map segments
            std::unordered_map<long,std::set<long>> oldSegs ;
            std::unordered_map<long,std::set<long>>::iterator oldSegsIt ;

            for ( auto & info : mapper ){
                if( info.entity == adaption::Entity::ENTITY_CELL ){

                    for ( auto & parent : info.previous){ //save old data and delete element
                        if( m_seg.exists(parent) ){
                            SegInfo *seg =  &m_seg[parent] ;

                            oldSegs.insert({{ parent, seg->m_segments }}) ;
                            m_seg.erase(parent,true) ;
                        }
                    }
                }
            }

            m_seg.flush() ;

            for ( auto & info : mapper ){ //forall mesh modifications
                if( info.entity == adaption::Entity::ENTITY_CELL){ //check if changes on cells
                    for ( auto & child : info.current){ // forall new elements

                        PiercedVector<SegInfo>::iterator seg =  m_seg.reclaim(child) ;
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

    } else { //size of narrow band increased -> recalculation

        m_seg.clear() ;
        associateSimplexToCell( visitee, visitee->getSizeNarrowBand() ) ; 

    };


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
        bitpit::genericIO::flushBINARY( stream, segItr->m_support );
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
        bitpit::genericIO::absorbBINARY( stream, cellData.m_support );
        bitpit::genericIO::absorbBINARY( stream, cellData.m_checked );

        std::copy( temp.begin(), temp.end(), std::inserter( cellData.m_segments, cellData.m_segments.end() ) );

        m_seg.insert(id,cellData) ;

    }

    return;
};

# if BITPIT_ENABLE_MPI

/*!
 * Repartioning of levelset after partitioning of mesh
 * @param[in] mapper mapper describing partitioning
 */
void LevelSetSegmentation::writeCommunicationBuffer( const std::vector<long> &previous, OBinaryStream &sizeBuffer, OBinaryStream &dataBuffer ){


    long nItems = previous.size() ;
    int dataSize = 10*sizeof(long)  +sizeof(long) +sizeof(bool) +sizeof(long) +sizeof(int) ;

    //TODO new BITPIT dataBuffer.setCapacity(nItems*dataSize) ;

    //determine elements to send
    nItems = 0 ;
    for( const auto &index : previous){
        if( m_seg.exists(index)){
            const auto &seginfo = m_seg[index] ;
            dataBuffer << index ;
            dataBuffer << (int) seginfo.m_segments.size() ;
            for( const auto & seg : seginfo.m_segments ){
                dataBuffer << seg ;
            };
            dataBuffer << seginfo.m_support ;
            dataBuffer << seginfo.m_checked ;
            ++nItems ;
        }
    }

    //TODO new BITPIT dataBuffer.squeeze() ;
    sizeBuffer << nItems ;
    //TODO new BITPIT sizeBuffer << dataBuffer.capacity() ;

    return;
};

/*!
 * Repartioning of levelset after partitioning of mesh
 * @param[in] mapper mapper describing partitioning
 */
void LevelSetSegmentation::readCommunicationBuffer( const long &nItems, IBinaryStream &dataBuffer ){

    int     s, nSegs ;
    long    index, segment ;
    PiercedVector<SegInfo>::iterator segItr ;

    for( int i=0; i<nItems; ++i){
        dataBuffer >> index ;
        dataBuffer >> nSegs ;

        segItr = m_seg.reclaim(index) ;

        for( s=0; s<nSegs; ++s){
            dataBuffer >> segment ;
            segItr->m_segments.insert(segment) ;
        }

        dataBuffer >> segItr->m_support ;
        dataBuffer >> segItr->m_checked ;
    }

    return;
};
# endif

}
