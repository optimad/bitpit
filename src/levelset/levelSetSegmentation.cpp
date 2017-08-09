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
    @class      SegmentationKernel
    @ingroup    levelset
    @brief      Segmentation kernel
*/

/*!
 * Default constructor
 */
SegmentationKernel::SegmentationKernel( ) : m_surface(nullptr), m_featureAngle(0) {
}

/*!
 * Constructor
 */
SegmentationKernel::SegmentationKernel( std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle ) {

    m_ownedSurface = std::shared_ptr<const SurfUnstructured>(surface.release());

    setSurface(m_ownedSurface.get(), featureAngle);
}

/*!
 * Constructor
 */
SegmentationKernel::SegmentationKernel( const SurfUnstructured *surface, double featureAngle ) {

    setSurface(surface, featureAngle);
}

/*!
 * Get feature angle
 * @return feature angle used when calculating face normals;
 */
double SegmentationKernel::getFeatureAngle() const {
    return m_featureAngle;
}

/*!
 * Get segmentation vertex normals
 * @return segmentation vertex normals;
 */
const std::unordered_map<long, std::vector< std::array<double,3>>> & SegmentationKernel::getVertexNormals() const {
    return m_vertexNormals;
}

/*!
 * Get segmentation vertex gradients
 * @return segmentation vertex gradients;
 */
const std::unordered_map<long, std::vector< std::array<double,3>>> & SegmentationKernel::getVertexGradients() const {
    return m_vertexGradients;
}

/*!
 * Get segmentation surface
 * @return segmentation surface;
 */
const SurfUnstructured & SegmentationKernel::getSurface() const {
    return *m_surface;
}

/*!
 * Set the surface
 * @param[in] patch pointer to surface
 */
void SegmentationKernel::setSurface( const SurfUnstructured *surface, double featureAngle){

    std::vector<std::array<double,3>> vertexNormal ;
    std::vector<std::array<double,3>> vertexGradient ;

    m_surface      = surface;
    m_featureAngle = featureAngle;

    double tol = m_surface->getTol() ;
    for( const Cell &segment : m_surface->getCells() ){
        long segmentId = segment.getId() ;
        int nVertices  = segment.getVertexCount() ;

        vertexNormal.resize(nVertices) ;
        vertexGradient.resize(nVertices) ;

        double misalignment = 0. ;
        for( int i = 0; i < nVertices; ++i ){
            vertexGradient[i] = m_surface->evalVertexNormal(segmentId, i) ;
            vertexNormal[i]   = m_surface->evalLimitedVertexNormal(segmentId, i, m_featureAngle) ;

            misalignment += norm2(vertexGradient[i] - vertexNormal[i]) ;
        }

        m_vertexGradients.insert({{segmentId, vertexGradient}}) ;
        if( misalignment >= tol ){
            m_vertexNormals.insert({{segmentId, vertexNormal}}) ;
        }
    }

    // Initialize search tree
    m_searchTreeUPtr = std::unique_ptr<SurfaceSkdTree>(new SurfaceSkdTree(surface));
    m_searchTreeUPtr->build();
}

/*!
 * Get the coordinates of the specified segment's vertices.
 * @param[in] id segmment's id
 * @param[out] coords on output will contain coordinates of the vertices
 */
void SegmentationKernel::getSegmentVertexCoords( long id, std::vector<std::array<double,3>> *coords ) const {

    const Cell &segment = m_surface->getCell(id) ;
    int nVertices = segment.getVertexCount() ;

    coords->resize(nVertices);
    for (int n = 0; n < nVertices; ++n) {
        long vertexId = segment.getVertex(n) ;
        (*coords)[n] = m_surface->getVertexCoords(vertexId);
    }

    if ( nVertices > 3 ) {
        log::cout() << "levelset: only segments and triangles supported in LevelSetSegmentation !!" << std::endl ;
    }
}

/*!
 * Computes levelset relevant information at one point with respect to a segment
 * @param[in] p coordinates of point
 * @param[in] i index of segment
 * @param[out] d distance point to segment
 * @param[out] s sign of point wrt to segment, i.e. according to normal.
 * Care should be taken since the method could return erroneous information when
 * the point p lies on the normal plane. In this case s=0 but due to surface curvature
 * the point may not lie necessary on the surface. This sititaion is easily indentified
 * because the distance != 0.
 * @param[out] x closest point on segment
 * @param[out] n normal at closest point
 */
void SegmentationKernel::getSegmentInfo( const std::array<double,3> &p, const long &i, double &d, double &s, std::array<double,3> &x, std::array<double,3> &n ) const {

    std::array<double,3> g ;

    auto itrNormal = getVertexNormals().find(i) ;
    auto itrGradient = getVertexGradients().find(i) ;
    assert( itrGradient != getVertexGradients().end() ) ;

    const Cell &cell = m_surface->getCell(i) ;
    int nVertices = cell.getVertexCount() ;
    switch (nVertices) {

    case 1:
    {
        long id = cell.getVertex(0) ;
        d = norm2( p- m_surface->getVertexCoords(id) ) ;
        g.fill(0.) ;
        n.fill(0.) ;

        break;
    }

    case 2:
    {
        long id0 = cell.getVertex(0) ;
        long id1 = cell.getVertex(1) ;

        std::array<double,2> lambda ;
        int flag ;

        d= CGElem::distancePointSegment( p, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), x, lambda, flag ) ;

        g = p-x;
        g /= norm2(g);

        n  = lambda[0] *itrGradient->second[0] ;
        n += lambda[1] *itrGradient->second[1] ;
        n /= norm2(n) ;

        g *= sign(dotProduct(g,n));

        if( itrNormal != getVertexNormals().end() ){
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

        break;
    }

    case 3:
    {
        long id0 = cell.getVertex(0) ;
        long id1 = cell.getVertex(1) ;
        long id2 = cell.getVertex(2) ;

        std::array<double,3> lambda ;
        int flag ;

        d= CGElem::distancePointTriangle( p, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), m_surface->getVertexCoords(id2), x, lambda, flag ) ;

        g = p-x;
        g /= norm2(g);

        n  = lambda[0] *itrGradient->second[0] ;
        n += lambda[1] *itrGradient->second[1] ;
        n += lambda[2] *itrGradient->second[2] ;
        n /= norm2(n);

        g *= sign(dotProduct(g,n));

        if( itrNormal != getVertexNormals().end() ){
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

        break;
    }

    default:
    {
        log::cout() << " Segment not supported in SegmentationKernel::getSegmentInfo " << nVertices << std::endl ;

        break;
    }

    }

    s = sign( dotProduct(g, p - x) );

}

/*!
	@struct     LevelSetSegmentation::SegInfo
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
 * @param[in] capacity is the number of segments the data structure should
 * be able to contain without requiring a reallocation of its internal
 * structures
 */
LevelSetSegmentation::SegInfo::SegInfo( std::size_t capacity ) : segments(capacity), distances(capacity) {
    segments.clear();
    distances.clear();
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
 * Constructor
 * @param[in] capacity is the number of segments the data structure should
 * be able to contain without requiring a reallocation of its internal
 * structures
 */
void LevelSetSegmentation::SegInfo::initialize( std::size_t capacity ) {
    segments.clear() ;
    segments.reserve( capacity ) ;

    distances.clear() ;
    distances.reserve( capacity ) ;
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
}

/*!
 * Constructor
 * @param[in] id identifier of object
 */
LevelSetSegmentation::LevelSetSegmentation(int id) : LevelSetCachedObject(id), m_segmentation(nullptr) {
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL unique pointer to surface mesh
 * @param[in] featureAngle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentation::LevelSetSegmentation( int id, std::unique_ptr<const SurfUnstructured> &&STL, double featureAngle) :LevelSetSegmentation(id) {
    setSegmentation( std::move(STL), featureAngle );
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL pointer to surface mesh
 * @param[in] featureAngle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentation::LevelSetSegmentation( int id, const SurfUnstructured *STL, double featureAngle) :LevelSetSegmentation(id) {
    setSegmentation( STL, featureAngle );
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
 * @param[in] patch pointer to surface
 */
void LevelSetSegmentation::setSegmentation( const SurfUnstructured *surface, double featureAngle){

    m_segmentation = std::make_shared<const SegmentationKernel>(surface, featureAngle);
}

/*!
 * Set the segmentation
 * @param[in] patch pointer to surface
 */
void LevelSetSegmentation::setSegmentation( std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle){

    m_segmentation = std::make_shared<const SegmentationKernel>(std::move(surface), featureAngle);
}

/*!
 * Get a constant refernce to the segmentation
 * @return constant reference to the segmentation
 */
const SegmentationKernel & LevelSetSegmentation::getSegmentation() const {
    return *m_segmentation ;
}

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
int LevelSetSegmentation::getPart( const long &id ) const{

    auto itr = m_seg.find(id) ;
    if( itr != m_seg.end() ){
        const SurfUnstructured &m_surface = m_segmentation->getSurface();
        long support = itr->segments.front() ;
        return m_surface.getCell(support).getPID();
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

    auto itr = m_seg.find(id) ;
    if( itr != m_seg.end() ){
        return itr->segments.front() ;
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

    auto itr = m_seg.find(id) ;
    if( itr != m_seg.end() ){
        return itr->segments.size() ;
    } else {
        return 0 ;
    }

}

/*!
 * Get the list of simplices wich contain the cell centroid in their narrow band.
 * @param[in] i cell index
 * @return set with indices of simplices
 */
const std::vector<long> & LevelSetSegmentation::getSimplexList(const long &id) const{

    auto itr = m_seg.find(id) ;
    if( itr != m_seg.end() ){
        return itr->segments ;
    } else {
        return levelSetDefaults::LIST;
    }

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

    const SurfUnstructured &m_surface = m_segmentation->getSurface();

    int spaceDimension = m_surface.getSpaceDimension();
    if (spaceDimension == 2) {
        return m_surface.evalCellArea(id); //TODO check
    } else if (spaceDimension == 3) {
        int dummy;
        return m_surface.evalMinEdgeLength(id, dummy);
    }

    return (- levelSetDefaults::SIZE);
}

/*!
 * Get the smallest characterstic size within the triangultaion
 * @return smallest charcteristic size within the triangulation
 */
double LevelSetSegmentation::getMinSurfaceFeatureSize( ) const {

    const SurfUnstructured &m_surface = m_segmentation->getSurface();

    bool   minimumValid = false;
    double minimumSize  = levelSetDefaults::SIZE;
    for( const Cell &cell : m_surface.getCells() ){
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

    const SurfUnstructured &m_surface = m_segmentation->getSurface();

    double maximumSize = - levelSetDefaults::SIZE;
    for( const Cell &cell : m_surface.getCells() ){
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
                                                                  const SegmentToCellMap &segmentToCellMap ){

    log::cout() << "  Creating segment info... " << std::endl;

    // Count the number of segments that will be processed for each cell
    std::unordered_map<long, int> nSegmentsPerCell;
    for ( auto &entry : segmentToCellMap ) {
        const std::vector<long> &cellList = entry.second ;

        for ( long cell : cellList ) {
            auto countItr = nSegmentsPerCell.find( cell );
            if ( countItr != nSegmentsPerCell.end() ) {
                countItr->second++;
            } else {
                nSegmentsPerCell.emplace( cell, 1 );
            }
        }
    }

    // Create the needed dat strucures
    std::unordered_set<long> newSegInfo ;
    std::vector<std::array<double,3>> cloud ;

    // Add the segments info
    std::vector<std::array<double,3>> VS;
    for ( auto mapItr = segmentToCellMap.begin(); mapItr != segmentToCellMap.end(); ++mapItr) {
        // Segment data
        long segment = mapItr->first ;
        m_segmentation->getSegmentVertexCoords( segment, &VS ) ;

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
                segInfoItr = m_seg.emreclaim( cell, segmentCount );

                newSegInfo.insert( cell ) ;
            }

            // Add the segment info
            segInfoItr->segments.push_back(segment) ;
            segInfoItr->distances.push_back(segmentDistance) ;
        }
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
        utils::reorderVector(segmentRank, segments, nSegments) ;
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
            const double SHRINK_THRESHOLD = 0.5;

            distances.resize(nSegmentsToKeep);
            if (nSegmentsToKeep < SHRINK_THRESHOLD * distances.capacity()) {
                distances.shrink_to_fit();
            }

            std::vector<long> &segments = segItr->segments;
            segments.resize(nSegmentsToKeep);
            if (nSegmentsToKeep < SHRINK_THRESHOLD * segments.capacity()) {
                segments.shrink_to_fit();
            }
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
        m_segmentation->getSegmentInfo(centroid, support, d, s, xP, n);

        PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.find(id) ;
        if( lsInfoItr == m_ls.end() ){
            lsInfoItr = m_ls.emplace(id) ;
        }

        if( d < std::abs(lsInfoItr->value) ){
            lsInfoItr->value    = ( signd *s  + (!signd) *1.   ) *d ;
            lsInfoItr->gradient = ( signd *1. + (!signd) *s ) *n ;
        }
    }
}

/*!
 * Finds seed points in narrow band within a cartesian mesh for one simplex
 * @param[in] visitee cartesian mesh 
 * @param[in] VS Simplex
 * @param[out] I indices of seed points
 */
bool LevelSetSegmentation::seedNarrowBand( LevelSetCartesian *visitee, std::vector<std::array<double,3>> &VS, std::vector<long> &I){

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
    const SurfUnstructured &m_surface = m_segmentation->getSurface();
    m_surface.getBoundingBox(minP,maxP) ;
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
        std::unordered_set<long> addedCells = createSegmentInfo(m_kernelPtr, RSearch, segmentToCellMap) ;
        createLevelsetInfo( m_kernelPtr, signd, addedCells );

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        computeLSInNarrowBand( lsOctree, RSearch, signd ) ;

    }
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

    if( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        log::cout() << "Updating levelset within the narrow band... " << std::endl;
        updateLSInNarrowBand( lsOctree, mapper, RSearch, signd ) ;
        return;
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

    const SurfUnstructured                  &m_surface = m_segmentation->getSurface();

    std::vector<long>                       stack, temp ;
    std::vector< std::array<double,3> >     cloud ;

    std::vector<double>                     d;
    std::vector<double>::iterator           vit;

    std::vector<long>                       flag( mesh.getCellCount(), -1);

    std::vector<long>                       neighs ;

    stack.reserve(128) ;
    temp.reserve(128) ;

    SegmentToCellMap segmentToCellMap;
    segmentToCellMap.reserve(m_surface.getCellCount());

    log::cout() << "  Extracting segment-to-cell map from octree patch... " << std::endl;

    // --------------------------------------------------------------------------
    // COMPUTE THE SDF VALUE AT EACH MESH POINT                                   //
    //
    for (const Cell &segment : m_surface.getCells()) {
        long segmentId = segment.getId();

        std::vector<long> &cellList = segmentToCellMap[segmentId] ;

        // Segments vertex ------------------------------------------------------ //
        m_segmentation->getSegmentVertexCoords( segmentId, &VS ) ;
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
                    neighs.clear();
                    mesh.findCellFaceNeighs(cell, &neighs) ;

                    for( const auto &  neigh : neighs){
                        if( flag[neigh] != segmentId) {
                            temp.push_back( neigh) ;
                            flag[neigh] = segmentId ;
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

    }

    return segmentToCellMap ;

}

/*!
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetOctree *visitee, const double &RSearch, const bool &signd){


    VolumeKernel &mesh = *(visitee->getMesh()) ;

    bool adaptiveSearch(RSearch<0);
    double searchRadius = RSearch;

    long segmentId;
    double distance;
    std::unordered_set<long> addedCells;

    for( const Cell &cell : mesh.getCells() ){
        long cellId = cell.getId();
        std::array<double, 3> cellCentroid = visitee->computeCellCentroid(cellId);

        if(adaptiveSearch){
            double cellSize = mesh.evalCellSize(cellId);
            searchRadius = 0.5*sqrt(3.)*cellSize;
        }

        m_segmentation->m_searchTreeUPtr->findPointClosestCell(cellCentroid, searchRadius, &segmentId, &distance);

        if(segmentId>=0){

            double                s, d;
            std::array<double,3>  n, xP;
            m_segmentation->getSegmentInfo(cellCentroid, segmentId, d, s, xP, n);

            PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.find(cellId) ;
            if( lsInfoItr == m_ls.end() ){
                lsInfoItr = m_ls.emplace(cellId) ;
            }

            if( d < std::abs(lsInfoItr->value) ){
                lsInfoItr->value    = ( signd *s  + (!signd) *1.   ) *d ;
                lsInfoItr->gradient = ( signd *1. + (!signd) *s ) *n ;
            }
        }
    }
}

/*!
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetOctree *visitee, const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd){

    clearAfterMeshAdaption(mapper);

    VolumeKernel &mesh = *(visitee->getMesh()) ;

    bool adaptiveSearch(RSearch<0);
    double searchRadius = RSearch;

    long segmentId;
    double distance;

    for( const auto &event : mapper ){

        if( event.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for( const long &cellId : event.current ){
            std::array<double, 3> cellCentroid = visitee->computeCellCentroid(cellId);

            if(adaptiveSearch){
                double cellSize = mesh.evalCellSize(cellId);
                searchRadius = 0.5*sqrt(3.)*cellSize;
            }

            m_segmentation->m_searchTreeUPtr->findPointClosestCell(cellCentroid, searchRadius, &segmentId, &distance);

            if(segmentId>=0){

                double                s, d;
                std::array<double,3>  n, xP;
                m_segmentation->getSegmentInfo(cellCentroid, segmentId, d, s, xP, n);

                PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.find(cellId) ;
                if( lsInfoItr == m_ls.end() ){
                    lsInfoItr = m_ls.emplace(cellId) ;
                }

                if( d < std::abs(lsInfoItr->value) ){
                    lsInfoItr->value    = ( signd *s  + (!signd) *1.   ) *d ;
                    lsInfoItr->gradient = ( signd *1. + (!signd) *s ) *n ;
                }
            }
        }
    }
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
        auto seginfoItr = m_seg.find(index) ;
        if( seginfoItr != m_seg.end() ){
            nItems++ ;
            counter += seginfoItr->segments.size() ;
        }
    }

    dataBuffer << nItems ;
    dataBuffer.setSize(dataBuffer.getSize() +nItems* (sizeof(long) +sizeof(size_t)) +counter*(sizeof(long)+sizeof(double))) ;

    //determine elements to send
    counter= 0 ;
    for( const auto &index : sendList){
        auto seginfoItr = m_seg.find(index) ;
        if( seginfoItr != m_seg.end() ){
            dataBuffer << counter ;
            dataBuffer << (size_t) seginfoItr->segments.size() ;
            for( const long & seg : seginfoItr->segments ){
                dataBuffer << seg ;
            }
            for( const double & distance : seginfoItr->distances ){
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
        PiercedVector<SegInfo>::iterator segItr = m_seg.find(id) ;
        if( segItr == m_seg.end() ){
            segItr = m_seg.emplace(id) ;
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
