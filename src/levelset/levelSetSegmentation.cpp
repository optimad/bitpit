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

# include <cassert>

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"

# include "levelSetObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetSegmentation.hpp"

namespace bitpit {

/*!
	@class      LevelSetSegmentation::SegInfo
	@ingroup    levelset
	@brief      Information about the segments
*/

/*!
 * Default constructor 
 */
LevelSetSegmentation::SegInfo::SegInfo( ) {
}

/*!
 * Constructor
 * @param[in] _segments is the list of segment's ids
 * @param[in] _distances is the list of segment's distances
 */
LevelSetSegmentation::SegInfo::SegInfo( const std::vector<long> &_segments, const std::vector<double> &_distances )
    : segments(_segments), distances(_distances) {
}

/*!
	@class      LevelSetSegmentation
	@ingroup    levelset
	@brief      Implements visitor pattern fo segmentated geometries
*/

/*!
 * Destructor
 */
LevelSetSegmentation::~LevelSetSegmentation() {
    m_segmentation = NULL;
    m_vertexNormal.clear() ;
    m_vertexGradient.clear() ;
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] angle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentation::LevelSetSegmentation(int id, double angle) :LevelSetCachedObject(id) {
    setFeatureAngle(angle) ;
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL unique pointer to surface mesh
 * @param[in] angle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentation::LevelSetSegmentation( int id, std::unique_ptr<SurfUnstructured> &&STL, double angle) :LevelSetSegmentation(id,angle) {
    setSegmentation( std::move(STL) );
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL pointer to surface mesh
 * @param[in] angle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentation::LevelSetSegmentation( int id, SurfUnstructured *STL, double angle) :LevelSetSegmentation(id,angle) {
    setSegmentation( STL );
}

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetSegmentation::LevelSetSegmentation( const LevelSetSegmentation &other) : LevelSetSegmentation(other.getId() ) {

    m_segmentation = other.m_segmentation; 
    m_dimension = other.m_dimension ;
    m_featureAngle = other.m_featureAngle ;
    m_vertexNormal = other.m_vertexNormal ;
    m_vertexGradient = other.m_vertexGradient ;
    if (m_own != nullptr) {
        m_own = unique_ptr<SurfUnstructured>(new SurfUnstructured(*(other.m_own)));
    }
}

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

    std::vector<std::array<double,3>>   vertexNormal, vertexGradient ;

    m_segmentation = segmentation;
    m_dimension = m_segmentation->getSpaceDimension() ;

    int  i, nV;
    long segId ;
    double norm, tol = m_segmentation->getTol() ;

    for( auto & segment : m_segmentation->getCells() ){

        segId = segment.getId() ;
        nV = segment.getVertexCount() ;

        vertexNormal.resize(nV) ;
        vertexNormal.shrink_to_fit() ;

        vertexGradient.resize(nV) ;
        vertexGradient.shrink_to_fit() ;

        for(i=0; i<nV; ++i){
            vertexGradient[i] = m_segmentation->evalVertexNormal(segId,i) ;
            vertexNormal[i] = m_segmentation->evalLimitedVertexNormal(segId,i,m_featureAngle) ;
        }

        m_vertexGradient.insert({{segId,vertexGradient}}) ;

        norm = 0. ;
        for(i=0; i<nV; ++i){
            norm += norm2(vertexGradient[i] - vertexNormal[i]) ;
        }

        if( norm >= tol ){
            m_vertexNormal.insert({{segId,vertexNormal}}) ;
        }

    }
}

/*!
 * Get a constant refernce to the segmentation
 * @return constant reference to the segmentation
 */
const SurfUnstructured & LevelSetSegmentation::getSegmentation() const {
    return *m_segmentation ;
}

/*!
 * Set feature angle
 * @param[in] angle feature angle to be used when calculating face normals;
 */
void LevelSetSegmentation::setFeatureAngle( double angle){
    m_featureAngle= angle;
}

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
int LevelSetSegmentation::getPart( const long &id ) const{

    if( m_seg.exists(id)){
        long support = m_seg.at(id).segments.front() ;
        return m_segmentation->getCell(support).getPID();
    } else {
        return levelSetDefaults::PART ;
    }

}

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
long LevelSetSegmentation::getSupport( const long &id ) const{

    if( m_seg.exists(id)){
        return m_seg.at(id).segments.front() ;
    } else {
        return levelSetDefaults::SUPPORT ;
    }

}

/*!
 * Gets the number of support items within the narrow band of cell
 * @param[in] id index of cell
 * @return number of segments in narrow band 
 */
int LevelSetSegmentation::getSupportCount( const long &id ) const{

    if( m_seg.exists(id)){
        return m_seg.at(id).segments.size() ;
    } else {
        return 0 ;
    }

}

/*!
 * Get the list of simplices wich contain the cell centroid in their narrow band.
 * @param[in] i cell index
 * @return set with indices of simplices
 */
const std::vector<long> & LevelSetSegmentation::getSimplexList(const long &i) const{

    if( !m_seg.exists(i) ){
        return levelSetDefaults::LIST;
    } else {
        return ( m_seg[i].segments );
    }

}

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
    }

    if( N > 3){
        log::cout() << "levelset: only segments and triangles supported in LevelSetSegmentation !!" << std::endl ;
    }

    return VS;
}

/*!
 * Get size of support triangle
 * @param[in] i cell index
 * @return charcteristic size of support triangle
 */
double LevelSetSegmentation::getSurfaceFeatureSize( const long &i ) const {

    long support = getSupport(i);
    if (support == levelSetDefaults::SUPPORT) {
        return (- levelSetDefaults::SIZE);
    }

    return getSegmentSize(support);
}

/*!
 * Get the sie of a segment
 * @param[in] id is the id of the segment
 * @return charcteristic size of the segment
 */
double LevelSetSegmentation::getSegmentSize( long id ) const {

    int spaceDimension = m_segmentation->getSpaceDimension();
    if (spaceDimension == 2) {
        return m_segmentation->evalCellArea(id); //TODO check
    } else if (spaceDimension == 3) {
        int dummy;
        return m_segmentation->evalMinEdgeLength(id, dummy);
    }

    return (- levelSetDefaults::SIZE);
}

/*!
 * Get the smallest characterstic size within the triangultaion
 * @return smallest charcteristic size within the triangulation
 */
double LevelSetSegmentation::getMinSurfaceFeatureSize( ) const {

    bool   minimumValid = false;
    double minimumSize  = levelSetDefaults::SIZE;
    for( const Cell &cell : m_segmentation->getCells() ){
        double segmentSize = getSegmentSize(cell.getId());
        if (segmentSize < 0) {
            continue;
        }

        minimumValid = true;
        minimumSize  = std::min(segmentSize, minimumSize);
    }

    if (!minimumValid) {
        minimumSize = - levelSetDefaults::SIZE;
    }

    return minimumSize;
}

/*!
 * Get the largest characterstic size within the triangultaion
 * @return largest charcteristic size within the triangulation
 */
double LevelSetSegmentation::getMaxSurfaceFeatureSize( ) const {

    double maximumSize = - levelSetDefaults::SIZE;
    for( const Cell &cell : m_segmentation->getCells() ){
        double segmentSize = getSegmentSize(cell.getId());
        maximumSize = std::max(segmentSize, maximumSize);
    }

    return maximumSize;
}

/*!
 * Create the segment information for the specified segment-to-cell map
 * @param[in] visitee visited mesh
 * @param[in] search size of narrow band
 * @param[in] segmentToCellMap is the segment-to-cell map that will be
 * processed
 * @return The list of cells associated to the newly created segment
 * information.
 */
std::unordered_set<long> LevelSetSegmentation::createSegmentInfo( LevelSetKernel *visitee,
                                                                  const double &search,
                                                                  SegmentToCellMap &segmentToCellMap ){

    log::cout() << "  Creating segment info... " << std::endl;

    // Count the number of segments that will be processed for each cell
    std::unordered_map<long, int> nSegmentsPerCell;
    for ( auto &entry : segmentToCellMap ) {
        const std::vector<long> &cellList = entry.second ;

        for ( long cell : cellList ) {
            auto countItr = nSegmentsPerCell.find( cell );
            if ( countItr == nSegmentsPerCell.end() ) {
                countItr = nSegmentsPerCell.insert( { cell, 0 } ).first;
            }

            *(countItr)++;
        }
    }

    // Create the needed dat strucures
    std::unordered_set<long> newSegInfo ;
    std::vector<std::array<double,3>> cloud ;

    // Add the segments info
    for ( auto mapItr = segmentToCellMap.begin(); mapItr != segmentToCellMap.end(); ) {
        // Segment data
        long segment = mapItr->first ;
        std::vector<std::array<double,3>> VS = getSimplexVertices( segment ) ;

        // Cell data
        const std::vector<long> &cellList = mapItr->second ;
        size_t cellListCount = cellList.size();

        cloud.resize( cellListCount ) ;
        for ( size_t k = 0; k < cellListCount; ++k ) {
            cloud[k] = visitee->computeCellCentroid( cellList[k] ) ;
        }

        // Eval distances
        std::vector<double> distancesFromSegment = CGElem::distanceCloudSimplex( cloud, VS );
        for ( size_t k = 0; k < cellListCount; ++k ) {
            // Discard segments with a distance greater than the narrow band
            double segmentDistance = distancesFromSegment[k];
            if ( segmentDistance > search ) {
                continue ;
            }

            // Get or create the segment's info
            long cell = cellList[k] ;
            auto segInfoItr = m_seg.find( cell ) ;
            if ( segInfoItr == m_seg.end() ) {
                int segmentCount = nSegmentsPerCell.at(cell);
                segInfoItr = m_seg.emplace( cell );
                segInfoItr->segments.reserve( segmentCount );
                segInfoItr->distances.reserve( segmentCount );

                newSegInfo.insert( cell ) ;
            }

            // Add the segment info
            segInfoItr->segments.push_back(segment) ;
            segInfoItr->distances.push_back(segmentDistance) ;
        }

        mapItr = segmentToCellMap.erase( mapItr );
    }

    // Order the segments from the closes to the farthest from the body
    std::vector<size_t> distanceRank;
    std::vector<size_t> segmentRank;
    for ( long id : newSegInfo ) {
        auto segInfoItr = m_seg.find( id ) ;
        std::vector<long> &segments = segInfoItr->segments;
        size_t nSegments = segments.size();

        // Evaluate ranks
        distanceRank.resize(nSegments);
        for (size_t k = 0; k < nSegments; ++k) {
            distanceRank[k] = k;
        }

        std::vector<double> &distances = segInfoItr->distances;
        std::sort(distanceRank.begin(), distanceRank.end(), DistanceComparator(distances) ) ;

        segmentRank.resize(nSegments);
        for (size_t k = 0; k < nSegments; ++k) {
            segmentRank[k] = distanceRank[k];
        }

        // Order vectors
        utils::reorderVector(distanceRank, distances, nSegments) ;
        distances.shrink_to_fit();

        utils::reorderVector(segmentRank, segments, nSegments) ;
        segments.shrink_to_fit();
    }

    return newSegInfo ;
}

/*!
 * Update the segment list associated to the cells, keeping only the segments
 * with a distance from the body less than the specified narrow band size plus
 * the cell circumcenter.
 * @param[in] search size of narrow band
 */
void LevelSetSegmentation::updateSegmentList( const double &search) {


    log::cout() << "  Updating segment list for cells inside narrow band... " << std::endl;

    bitpit::PiercedVector<SegInfo>::iterator segEnd = m_seg.end();
    for (bitpit::PiercedVector<SegInfo>::iterator segItr = m_seg.begin(); segItr != segEnd; ++segItr) {
        // The filter radius is the narrow band size plus the cell circumcenter.
        long cellId = segItr.getId() ;
        double filterRadius = search + m_kernelPtr->computeCellCircumcircle(cellId);

        // Starting from the farthest segment (the last in the list) we loop
        // backwards until we find the first segment with a distance less
        // that the specified narrow band size.
        //
        // We know that the cell is inside the narrow band, because all the
        // cells outside the narrow band have already been removed. Therefore
        // we need to perform the check up to the second segment (the first one
        // is in the narrow band).
        std::vector<double> &distances = segItr->distances;
        size_t nCurrentSegments = distances.size();

        size_t nSegmentsToKeep = 1;
        for( size_t k = nCurrentSegments - 1; k >= 1; --k) {
            double distance = distances[k];
            if ( distance <= filterRadius ){
                nSegmentsToKeep = k + 1;
                break;
            }
        }

        if (nSegmentsToKeep != nCurrentSegments) {
            distances.resize(nSegmentsToKeep);
            distances.shrink_to_fit();

            std::vector<long> &segments = segItr->segments;
            segments.resize(nSegmentsToKeep);
            segments.shrink_to_fit();
        }
    }
}

/*!
 * Create the levelset info for the specified cell list.
 *
 * This function assumes that the segment information for the specified cells
 * have already been created.
 *
 * @param[in] visitee visited mesh
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 * @param[in] cellList is the list of cells that will be processed
 */
void LevelSetSegmentation::createLevelsetInfo( LevelSetKernel *visitee, const bool & signd,
                                               std::unordered_set<long> &cellList ){


    log::cout() << "  Creating levelset info... " << std::endl;

    for ( long id : cellList ) {
        auto segInfoItr = m_seg.find( id ) ;
        long support = segInfoItr->segments.front();
        const std::array<double,3> &centroid = visitee->computeCellCentroid(id) ;

        double                s, d;
        std::array<double,3>  n, xP;
        infoFromSimplex(centroid, support, d, s, xP, n);

        PiercedVector<LevelSetInfo>::iterator lsInfoItr ;
        if( !m_ls.exists(id)){
            lsInfoItr = m_ls.emplace(id) ;
        } else {
            lsInfoItr = m_ls.getIterator(id) ;
        }

        if( d < std::abs(lsInfoItr->value) ){
            lsInfoItr->value    = ( signd *s  + (!signd) *1.   ) *d ;
            lsInfoItr->gradient = ( signd *1. + (!signd) *s ) *n ;
        }
    }
}

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

    std::array<double,3> g ;

    Cell &cell = m_segmentation->getCell(i) ;
    int nV = cell.getVertexCount() ;

    auto itrNormal = m_vertexNormal.find(i) ;
    auto itrGradient = m_vertexGradient.find(i) ;
    assert( itrGradient != m_vertexGradient.end() ) ;

    if( nV == 1){
        long id = cell.getVertex(0) ;
        d = norm2( p- m_segmentation->getVertexCoords(id) ) ;
        g.fill(0.) ;
        n.fill(0.) ;

    } else if( nV == 2){
        long id0 = cell.getVertex(0) ;
        long id1 = cell.getVertex(1) ;

        std::array<double,2> lambda ;
        int flag ;

        d= CGElem::distancePointSegment( p, m_segmentation->getVertexCoords(id0), m_segmentation->getVertexCoords(id1), x, lambda, flag ) ;

        g = p-x;
        g /= norm2(g);

        n  = lambda[0] *itrGradient->second[0] ;
        n += lambda[1] *itrGradient->second[1] ;
        n /= norm2(n) ;

        g *= sign(dotProduct(g,n));

        if( itrNormal != m_vertexNormal.end() ){
            n  = lambda[0] *itrNormal->second[0] ;
            n += lambda[1] *itrNormal->second[1] ;
            n /= norm2(n) ;

            double kappa ;
            maxval(lambda,kappa);
            kappa = 1. -kappa;

            n *= kappa;
            n += (1.-kappa)*g;
            n /= norm2(n);

        }


    } else if (nV == 3){
        long id0 = cell.getVertex(0) ;
        long id1 = cell.getVertex(1) ;
        long id2 = cell.getVertex(2) ;

        std::array<double,3> lambda ;
        int flag ;

        d= CGElem::distancePointTriangle( p, m_segmentation->getVertexCoords(id0), m_segmentation->getVertexCoords(id1), m_segmentation->getVertexCoords(id2), x, lambda, flag ) ;

        g = p-x;
        g /= norm2(g);

        n  = lambda[0] *itrGradient->second[0] ;
        n += lambda[1] *itrGradient->second[1] ;
        n += lambda[2] *itrGradient->second[2] ;

        g *= sign(dotProduct(g,n));

        if( itrNormal != m_vertexNormal.end() ){
            n  = lambda[0] *itrNormal->second[0] ;
            n += lambda[1] *itrNormal->second[1] ;
            n += lambda[2] *itrNormal->second[2] ;
            n /= norm2(n) ;

            double kappa ;
            maxval(lambda,kappa);
            kappa = 1. -kappa;

            n *= kappa;
            n += (1.-kappa)*g;
            n /= norm2(n);

        }



    } else{
        log::cout() << " simplex not supported in LevelSetSegmentation::infoFromSimplex " << nV << std::endl ;
        
    }

    s = sign( dotProduct(g, p - x) );

}

/*!
 * Finds seed points in narrow band within a cartesian mesh for one simplex
 * @param[in] visitee cartesian mesh 
 * @param[in] VS Simplex
 * @param[out] I indices of seed points
 */
bool LevelSetSegmentation::seedNarrowBand( LevelSetCartesian *visitee, std::vector<std::array<double,3>> &VS, std::vector<int> &I){

    VolCartesian                        &mesh = *(static_cast<VolCartesian*>(visitee->getMesh()));

    bool                                found(false) ;
    int                                 dim( mesh.getDimension() ) ;
    std::array<double,3>                B0, B1;
    std::vector<std::array<double,3>>   VP ;

    mesh.getBoundingBox(B0, B1) ;

    for( int i=0; i<dim; ++i){
        B0[i] -= getSizeNarrowBand() ;
        B1[i] += getSizeNarrowBand() ;
    }

    I.clear() ;

    for( const auto &P : VS){
        if(  CGElem::intersectPointBox( P, B0, B1, dim ) ) {
            I.push_back( mesh.locateClosestCell(P) );
            found =  true ;
        }
    }

    if( !found && CGElem::intersectBoxSimplex( B0, B1, VS, VP, dim ) ) {
        for( const auto &P : VP){
            I.push_back( mesh.locateClosestCell(P) );
            found = true ;
        }
    }

    return found ;
}

/*!
 * Computes axis aligned bounding box of object
 * @param[out] minP minimum point
 * @param[out] maxP maximum point
 */
void LevelSetSegmentation::getBoundingBox( std::array<double,3> &minP, std::array<double,3> &maxP ) const {
    m_segmentation->getBoundingBox(minP,maxP) ;
}

/*!
 * Clear the segmentation and the specified kernel.
 */
void LevelSetSegmentation::__clear( ){

    m_seg.clear() ;
}

/*!
 * Computes the levelset function within the narrow band
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetSegmentation::computeLSInNarrowBand( const double &RSearch, const bool &signd ){

    log::cout() << "Computing levelset within the narrow band... " << std::endl;

    SegmentToCellMap segmentToCellMap;
    if( LevelSetCartesian* lsCartesian = dynamic_cast<LevelSetCartesian*>(m_kernelPtr) ){
        segmentToCellMap = extractSegmentToCellMap( lsCartesian, RSearch ) ;

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        segmentToCellMap = extractSegmentToCellMap( lsOctree, RSearch ) ;

    }

    std::unordered_set<long> addedCells = createSegmentInfo(m_kernelPtr, RSearch, segmentToCellMap) ;

    createLevelsetInfo( m_kernelPtr, signd, addedCells );
}

/*!
 * Updates the levelset function within the narrow band after mesh adaptation.
 * @param[in] mapper information concerning mesh adaption 
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetSegmentation::updateLSInNarrowBand( const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd ){

    // Update is not implemented for Cartesian patches
    if( dynamic_cast<LevelSetCartesian*>(m_kernelPtr) ){
        clear( ) ;
        computeLSInNarrowBand( RSearch, signd ) ;
        return;
    }

    log::cout() << "Updating levelset within the narrow band... " << std::endl;

    // Detect changes in narrow band size
    int narrowBandResizeDirection = 0;
    if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        narrowBandResizeDirection = getNarrowBandResizeDirection( lsOctree, RSearch ) ;
    }

    // If the narrow band size has been increased we can't just update the
    // levelset, we need to rebuild it from scratch.
    if (narrowBandResizeDirection > 0) {
        clear( ) ;
        computeLSInNarrowBand( RSearch, signd ) ;
        return;
    }

    // If the narrow band size has been decreased or is the same as before
    // it is possible to update the segment to cell association
    SegmentToCellMap segmentToCellMap = extractSegmentToCellMap( mapper ) ;

    // Prune previous segment info
    //clearAfterMeshAdaptionDerived( mapper ) ;
    clearAfterMeshAdaption( mapper ) ;

    // Evaluate the levelset for the newly added elements
    std::unordered_set<long> addedCells;
    if (segmentToCellMap.size() != 0 ) {
        addedCells = createSegmentInfo(m_kernelPtr, RSearch, segmentToCellMap) ;

        createLevelsetInfo( m_kernelPtr, signd, addedCells );
    }

}

/*!
 * Determines the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle) for cartesian meshes
 * @param[in] visitee pointer to cartesian mesh
 * @param[in] RSearch size of narrow band
 */
LevelSetSegmentation::SegmentToCellMap LevelSetSegmentation::extractSegmentToCellMap( LevelSetCartesian *visitee, const double &RSearch ){

    VolumeKernel                            &mesh = *(visitee->getMesh() ) ;
    std::vector<std::array<double,3>>       VS(3);

    std::vector< int >                      stack, temp ;
    std::vector< std::array<double,3> >     cloud ;

    std::vector<double>                     d;
    std::vector<double>::iterator           vit;

    std::vector<int>                        flag( mesh.getCellCount(), -1);

    int                                     i, N( m_segmentation->getCellCount() );

    std::vector<long>                       neighs ;

    stack.reserve(128) ;
    temp.reserve(128) ;

    SegmentToCellMap segmentToCellMap;
    segmentToCellMap.reserve(m_segmentation->getCellCount());

    log::cout() << "  Extracting segment-to-cell map from octree patch... " << std::endl;

    // --------------------------------------------------------------------------
    // COMPUTE THE SDF VALUE AT EACH MESH POINT                                   //
    //
    for (i = 0; i < N; i++) {

        std::vector<long> &cellList = segmentToCellMap[i] ;

        // Segments vertex ------------------------------------------------------ //
        VS  = getSimplexVertices( i ) ;
        seedNarrowBand( visitee, VS, stack );


        //-----------------------------------------------------------------
        //Propagate from seed
        size_t stackSize = stack.size();
        while (stackSize > 0) {

            // Extract point from lifo
            cloud.resize(stackSize) ;

            for( size_t k = 0; k < stackSize; ++k) {
                long cell = stack[k];
                cloud[k] = visitee->computeCellCentroid(cell) ;
            }

            d = CGElem::distanceCloudSimplex( cloud, VS); 
            vit = d.begin() ;

            for( const auto & cell : stack){
                if ( *vit <= RSearch ) {

                    cellList.push_back( cell ) ;
                    neighs  = mesh.findCellFaceNeighs(cell) ;

                    for( const auto &  neigh : neighs){
                        if( flag[neigh] != i) {
                            temp.push_back( neigh) ;
                            flag[neigh] = i ;
                        }
                    }

                } //end if distance

                ++vit ;

            }

            stack.clear() ;
            stack.swap( temp ) ;
            stackSize = stack.size() ;


        } //end while

        // The list of cells has to be unique
        std::sort( cellList.begin(), cellList.end() ) ;
        cellList.erase( std::unique(cellList.begin(), cellList.end()), cellList.end() ) ;

    } //end for i

    return segmentToCellMap ;

}

/*!
 * Determines the list of triangles which influence each cell (i.e. cells which are within the narrow band of the triangle) for octree meshes
 * @param[in] visitee pointer to octree mesh
 * @param[in] RSearch size of narrow band
 */
LevelSetSegmentation::SegmentToCellMap LevelSetSegmentation::extractSegmentToCellMap( LevelSetOctree *visitee, const double &RSearch){

    VolumeKernel                &mesh = *(visitee->getMesh()) ;
    int                         dim(mesh.getDimension()) ;
    long                        id, icart;
    int                         i;
    double                      size ;

    std::array<double,3>        C0, C1, octrBB0, octrBB1, triBB0, triBB1 ;
    PiercedVector<SegInfo>::iterator data;

    SegmentToCellMap segmentToCellMap;
    segmentToCellMap.reserve(m_segmentation->getCellCount());

    log::cout() << "  Extracting segment-to-cell map from Octree patch... " << std::endl;

    // mesh size corresponding to RSearch
    size = visitee->computeSizeFromRSearch( RSearch ) ;

    mesh.getBoundingBox(octrBB0,octrBB1) ;
    getBoundingBox( triBB0, triBB1 );

    triBB0 -= RSearch ;
    triBB1 += RSearch ;

    if( CGElem::intersectBoxBox(octrBB0,octrBB1,triBB0,triBB1,C0,C1) ) { //intersect two Bounding Boxes around geometry and local grid

        // snap bounding box to octree grid and create cartesian grid
        std::array<int,3>    nc ;

        for( i=0; i<dim; ++i){
            C0[i] =  octrBB0[i] + size *   (int) ( ( C0[i] - octrBB0[i] ) / size ) ;
            C1[i] =  octrBB0[i] + size * ( (int) ( ( C1[i] - octrBB0[i] ) / size ) +1 ) ;

            nc[i] = round( ( C1[i] - C0[i] ) /size ) ;
        }

        // calculate LS triangle lists on cartesian mesh and map on pablo
        VolCartesian            cmesh( 0, dim, C0, C1-C0, nc) ;

        LevelSetCartesian       auxLS(cmesh) ;
        LevelSetSegmentation    objLS(*this) ;
        objLS.setKernel(&auxLS);

        double localRSearch = (1 + std::sqrt(3.) / 2.) * objLS.computeSizeNarrowBand() ;

        objLS.setSizeNarrowBand(localRSearch);
        SegmentToCellMap auxSegmentToCellMap = extractSegmentToCellMap( &auxLS, localRSearch ) ;
        std::unordered_set<long> inNarrowBandCells = objLS.createSegmentInfo( &auxLS, localRSearch, auxSegmentToCellMap ) ;

        for( auto & cell : mesh.getCells() ){
            id = cell.getId() ;

            const std::array<double,3> &C = visitee->computeCellCentroid(id) ;

            icart = cmesh.locatePoint(C) ;
            if ( icart != Cell::NULL_ID ) {
                if( inNarrowBandCells.count(icart)== 0 ){
                    continue;
                }

                for ( long segment : objLS.getSimplexList(icart) ) {
                    segmentToCellMap[segment].push_back(id);
                }

            }

        }

    } //endif intersect

    return segmentToCellMap;

}

/*!
 * Extract the map that links segments and cells
 * @param[in] mapper information concerning mesh adaption
 * @return The map that links segments and cells
 */
LevelSetSegmentation::SegmentToCellMap LevelSetSegmentation::extractSegmentToCellMap( const std::vector<adaption::Info> &mapper ){

    log::cout() << "  Extracting segment-to-cell map from mapper... " << std::endl;

    // Count the cells that will be associated to the segments
    std::unordered_map<long, size_t> cellsPerSegment;
    for ( const auto &info : mapper ){
        // Consider only changes on cells
        if( info.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        // Get the number of child cells
        long nChildElements = info.current.size() ;
        if (nChildElements == 0) {
            continue;
        }

        // Associate the childs to the list of parent's segments
        for ( const long & parent : info.previous ) {
            PiercedVector<SegInfo>::const_iterator parentSegInfoItr = m_seg.find(parent) ;
            if ( parentSegInfoItr == m_seg.cend() ) {
                continue;
            }

            for ( long segment : parentSegInfoItr->segments ) {
                cellsPerSegment[segment] += nChildElements ;
            }
        }
    }

    // Extract the segment-to-cell map
    SegmentToCellMap segmentToCellMap;
    segmentToCellMap.reserve(cellsPerSegment.size());

    std::unordered_set<long> removeDuplicateList;
    removeDuplicateList.reserve(cellsPerSegment.size());

    for ( const auto &info : mapper ){
        // Consider only changes on cells
        if( info.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        // Get the number of parent cells
        long nParentCells = info.previous.size() ;
        if (nParentCells == 0) {
            continue;
        }

        // Get the number of child cells
        long nChildElements = info.current.size() ;
        if (nChildElements == 0) {
            continue;
        }

        // Associate the childs to the list of parent's segments
        bool possibleDuplicates = (nParentCells > 1);
        for ( const long & parent : info.previous ) {
            PiercedVector<SegInfo>::const_iterator parentSegInfoItr = m_seg.find(parent) ;
            if ( parentSegInfoItr == m_seg.cend() ) {
                continue;
            }

            for ( long segment : parentSegInfoItr->segments ) {
                if ( possibleDuplicates ) {
                    removeDuplicateList.insert( segment ) ;
                }

                std::vector<long> &cellList = segmentToCellMap[segment] ;
                if ( cellList.capacity() == 0 ) {
                    cellList.reserve( cellsPerSegment.at(segment) ) ;
                }
                cellList.insert( cellList.end(), info.current.begin(), info.current.end() ) ;
            }
        }
    }

    // Remove duplicate entries
    for ( long segment : removeDuplicateList ) {
        std::vector<long> &cellList = segmentToCellMap.at(segment) ;
        std::sort( cellList.begin(), cellList.end() ) ;
        cellList.erase( std::unique(cellList.begin(), cellList.end()), cellList.end() ) ;
        cellList.shrink_to_fit();
    }

    return segmentToCellMap;
}

/*!
 * Detects if the requested narrow band size will make the narrow band grow
 * or shrink.
 * @param[in] visitee pointer to octree mesh
 * @param[in] newRSearch new size of narrow band
 * @return Returns 0 if the outer limit of the narrow band will not change,
 * +1 is the narrow band will be growth and -1 is the narrow band will shrink.
 */
int LevelSetSegmentation::getNarrowBandResizeDirection( LevelSetOctree *visitee, const double &newRSearch){

    double oldCellSize = visitee->computeSizeFromRSearch( getSizeNarrowBand() ) ;
    double newCellSize = visitee->computeSizeFromRSearch( newRSearch ) ;

    return sign( newCellSize - oldCellSize );
}

/*!
 * Prune the segment's info removing entries associated to cells that are
 * are not in the mesh anymore
 * @param[in] mapper information concerning mesh adaption
 */
void LevelSetSegmentation::__clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

    log::cout() << "  Clearing segment info... " << std::endl;

    for ( const auto &info : mapper ){
        // Consider only changes on cells
        if( info.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        // Delete only old data that belongs to the current processor
        if (info.type == adaption::Type::TYPE_PARTITION_RECV) {
            continue;
        }

        // Remove info of previous cells
        for ( const long & parent : info.previous ) {
            if ( m_seg.find( parent ) == m_seg.end() ) {
                continue;
            }

            m_seg.erase( parent, true ) ;
        }
    }

    m_seg.flush();
}

/*!
 * Clears data structure outside narrow band
 * @param[in] search size of narrow band
 */
void LevelSetSegmentation::__filterOutsideNarrowBand( double search ){

    long id ;

    bitpit::PiercedVector<SegInfo>::iterator segItr ;
    for( segItr = m_seg.begin(); segItr != m_seg.end(); ++segItr){
        id = segItr.getId() ;
        if( ! isInNarrowBand(id) ){
            m_seg.erase(id,true) ;
        }
    }

    m_seg.flush() ;

    updateSegmentList( search) ;

}

/*!
 * Writes LevelSetSegmentation to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetSegmentation::__dump( std::ostream &stream ){

    utils::binary::write( stream, m_seg.size() ) ;

    bitpit::PiercedVector<SegInfo>::iterator segItr, segEnd = m_seg.end() ;
    for( segItr = m_seg.begin(); segItr != segEnd; ++segItr){
        utils::binary::write( stream, segItr.getId() );
        utils::binary::write( stream, segItr->segments.size() );
        utils::binary::write( stream, segItr->segments );
        utils::binary::write( stream, segItr->distances );
    }
}

/*!
 * Reads LevelSetSegmentation from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetSegmentation::__restore( std::istream &stream ){

    size_t segSize;
    utils::binary::read( stream, segSize ) ;
    m_seg.reserve(segSize);

    for( size_t i=0; i<segSize; ++i){
        long id;
        utils::binary::read( stream, id );

        long nSegments;
        utils::binary::read( stream, nSegments );

        SegInfo cellData ;

        cellData.segments.resize(nSegments) ;
        utils::binary::read( stream, cellData.segments );

        cellData.distances.resize(nSegments) ;
        utils::binary::read( stream, cellData.distances );

        m_seg.insert(id,cellData) ;
    }
}

# if BITPIT_ENABLE_MPI

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetSegmentation::__writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    long nItems(0), counter(0) ;

    //determine number of elements to send
    for( const auto &index : sendList){
        if( m_seg.exists(index)){
            nItems++ ;
            const auto &seginfo = m_seg[index] ;
            counter += seginfo.segments.size() ;
        }

    }

    dataBuffer << nItems ;
    dataBuffer.setCapacity(dataBuffer.capacity() +nItems* (sizeof(long) +sizeof(size_t)) +counter*(sizeof(long)+sizeof(double))) ;

    //determine elements to send
    counter= 0 ;
    for( const auto &index : sendList){
        if( m_seg.exists(index)){
            const auto &seginfo = m_seg[index] ;
            dataBuffer << counter ;
            dataBuffer << (size_t) seginfo.segments.size() ;
            for( const long & seg : seginfo.segments ){
                dataBuffer << seg ;
            }
            for( const double & distance : seginfo.distances ){
                dataBuffer << distance ;
            }
        }

        ++counter ;
    }
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetSegmentation::__readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    long    nItems, index, id ;

    dataBuffer >> nItems ;

    for( int i=0; i<nItems; ++i){

        // Determine the id of the element
        dataBuffer >> index ;
        id = recvList[index] ;

        // Assign the data of the element
        PiercedVector<SegInfo>::iterator segItr ;
        if( !m_seg.exists(id)){
            segItr = m_seg.emplace(id) ;
        } else {
            segItr = m_seg.getIterator(id) ;
        }

        size_t nSegs ;
        dataBuffer >> nSegs ;
        segItr->segments.resize(nSegs) ;
        for( size_t s=0; s<nSegs; ++s){
            dataBuffer >> segItr->segments[s] ;
        }
        segItr->distances.resize(nSegs) ;
        for( size_t s=0; s<nSegs; ++s){
            dataBuffer >> segItr->distances[s] ;
        }
    }
}
# endif

}
