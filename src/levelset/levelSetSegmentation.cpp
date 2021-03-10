/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
 *
 * @param[in,out] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
SegmentationKernel::SegmentationKernel( std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle ) {

    m_ownedSurface = std::move(surface);

    setSurface(m_ownedSurface.get(), featureAngle);
}

/*!
 * Constructor
 *
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
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
 * Get segmentation surface
 * @return segmentation surface;
 */
const SurfUnstructured & SegmentationKernel::getSurface() const {
    return *m_surface;
}

/*!
 * Get search tree
 * @return search tree;
 */
const SurfaceSkdTree & SegmentationKernel::getSearchTree() const {
    return *m_searchTree;
}

/*!
 * Set the surface
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
void SegmentationKernel::setSurface( const SurfUnstructured *surface, double featureAngle){

    // Check if adjacencies are built
    if (surface->getAdjacenciesBuildStrategy() == SurfUnstructured::ADJACENCIES_NONE) {
        throw std::runtime_error ("Segmentation needs adjacencies!") ;
    }

    // Surface information
    m_surface      = surface;
    m_featureAngle = featureAngle;

    // Check if segment is supported
    for( SurfUnstructured::CellConstIterator segmentItr = m_surface->cellConstBegin(); segmentItr != m_surface->cellConstEnd(); ++segmentItr ){
        switch (segmentItr->getType()) {

        case ElementType::VERTEX :
        case ElementType::LINE :
        case ElementType::TRIANGLE :
            break ;

        default:
            throw std::runtime_error ("levelset: only segments and triangles supported in LevelSetSegmentation!") ;
            break ;

        }
    }

    // Segment vertices information
    m_segmentVertexOffset.setStaticKernel(&m_surface->getCells());

    std::size_t nTotalSegmentVertices = 0;
    for( auto itr = m_segmentVertexOffset.begin(); itr != m_segmentVertexOffset.end(); ++itr ){
        *itr = nTotalSegmentVertices;
        nTotalSegmentVertices += m_surface->getCells().rawAt(itr.getRawIndex()).getVertexCount();
    }

    // Normals
    m_segmentNormalsValid.setStaticKernel(&m_surface->getCells());
    m_segmentNormalsValid.fill(false);
    m_segmentNormalsStorage.setStaticKernel(&m_surface->getCells());

    m_unlimitedVertexNormalsValid.setStaticKernel(&m_surface->getVertices());
    m_unlimitedVertexNormalsValid.fill(false);
    m_unlimitedVertexNormalsStorage.setStaticKernel(&m_surface->getVertices());

    m_limitedSegmentVertexNormalValid.resize(nTotalSegmentVertices);

    // Initialize search tree
    m_searchTree = std::unique_ptr<SurfaceSkdTree>(new SurfaceSkdTree(surface));
    m_searchTree->build();
}

/*!
 * Computes levelset relevant information at one point with respect to a segment
 *
 * @param[in] pointCoords coordinates of point
 * @param[in] segmentId index of segment
 * @param[in] signd true is signed distance should be computed
 * @param[out] distance distance point to segment
 * @param[out] gradient levelset gradient
 * @param[out] normal normal at closest point
 */
int SegmentationKernel::getSegmentInfo( const std::array<double,3> &pointCoords, long segmentId, bool signd, double &distance, std::array<double,3> &gradient, std::array<double,3> &normal ) const {

    // Segment information
    SurfUnstructured::CellConstIterator segmentIterator = m_surface->getCellConstIterator(segmentId);
    const Cell &segment = *segmentIterator ;
    ElementType segmentType = segment.getType();
    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds() ;
    int nSegmentVertices = segmentVertexIds.size() ;

    // Projct the point on the surface and evaluate the point-projeciont vector
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    std::array<double,3> pointProjectionVector = pointCoords;
    switch (segmentType) {

    case ElementType::VERTEX :
    {
        long id = segmentVertexIds[0] ;
        pointProjectionVector -= m_surface->getVertexCoords(id);

        break;
    }

    case ElementType::LINE:
    {
        long id0 = segmentVertexIds[0] ;
        long id1 = segmentVertexIds[1] ;
        pointProjectionVector -= CGElem::projectPointSegment( pointCoords, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), lambda);

        break;
    }

    case ElementType::TRIANGLE:
    {
        long id0 = segmentVertexIds[0] ;
        long id1 = segmentVertexIds[1] ;
        long id2 = segmentVertexIds[2] ;
        pointProjectionVector -= CGElem::projectPointTriangle( pointCoords, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), m_surface->getVertexCoords(id2), lambda );

        break;
    }

    default:
    {
        std::runtime_error ("Type of cell not supported.");
        break;
    }

    }

    // Compute surface normal
    normal = computeSurfaceNormal(segmentIterator, lambda);

    // Evaluate distance from surface
    distance = norm2(pointProjectionVector);

    // Check if the point lies on the segmentation
    //
    // If the distance is zero, the point and the projection are coincident,
    // this means that the point lies on the segmentation.
    double distanceTolerance = m_surface->getTol();
    bool pointOnSegmentation = utils::DoubleFloatingEqual()(distance, 0., distanceTolerance, distanceTolerance);

    // Evaluate levelset gradient
    if (!pointOnSegmentation) {
        gradient = pointProjectionVector / distance;
    } else {
        if (signd) {
            gradient = normal;
        } else {
            gradient = {{0., 0., 0.}};
        }
    }

    // Evaluate levelset sign
    //
    // The sign is computed by determining the side of the point with respect
    // to the normal plane. The sign will be zero if the point lies exaclty
    // on the segmentation or on the normal plane. In the latter case the sign
    // must be evaluated taking into account the the curvature of the surface.
    // However, this is not yet implemented.
    std::array<double, 3> pseudoNormal = computePseudoNormal(segmentIterator, lambda);
    int s = sign( dotProduct(pointProjectionVector, pseudoNormal) );
    if (!pointOnSegmentation && s == 0) {
        distance = levelSetDefaults::VALUE;
        gradient = levelSetDefaults::GRADIENT;
        normal   = levelSetDefaults::GRADIENT;

        return 1;
    }

    // Use sign to update levelset information
    //
    // If signed distance are computed, the distance value and gradient
    // need to be changed accordingly. If unsigned distance are computed
    // the orientation of the suraface normal is discarded and in order
    // to agnostic with repect the two sides of the surface
    if (s < 0) {
        distance *= (double) ( signd *s + (!signd) *1);
        gradient *= (double) ( signd *s + (!signd) *1);
        normal   *= (double) ( signd *1 + (!signd) *s);
    }

    return 0;
}

/*!
 * Compute the pseudo-normal at specified point of the given triangle.
 *
 * The algorithm used to evaluate the pseudo-normal depends on the location
 * of the point within the triangle:
 *  - if the point coincides with a vertex, the pseudo-normal is evaluated
 *    as the unlimited normal of the vertex;
 *  - if the point lies on an edge, the pseudo-normal is the average of
 *    the normals of the two segments sharing the edge;
 *  - if the point is inside the segment, the pseudo-normal is evaluated as
 *    the normal of the segment.
 *
 * See "Signed Distance Computation Using the Angle Weighted Pseudo-normal",
 * J. Andreas Bearentzen, Henrik Aanaes, IEEE Transactions on Visualization
 * and Computer Graphics, 2005.
 *
 * To reduce computational times, normals of segments and vertices are cached.
 *
 * @param[in] segmentIterator is an iterator pointing to the segment
 * @param[in] lambda are the barycentric coordinates of the point
 * @return the pseudo-normal at specified point of the given triangle
 */
std::array<double,3> SegmentationKernel::computePseudoNormal( const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const {

    // Early return if the segment is a point
    const Cell &segment = *segmentIterator;
    ElementType segmentType = segment.getType();
    if (segmentType == ElementType::VERTEX) {
        return {{0., 0., 0.}};
    }

    // Evaluate pseudo normal
    int positionFlag;
    if (segmentType == ElementType::LINE) {
        positionFlag = CGElem::convertBarycentricToFlagSegment(lambda, m_surface->getTol());
    } else {
        int nSegmentVertices = segment.getVertexCount();
        positionFlag = CGElem::convertBarycentricToFlagPolygon(nSegmentVertices, lambda, m_surface->getTol());
    }

    std::array<double,3> pseudoNormal;
    if (positionFlag == 0) {
        pseudoNormal = computeSegmentNormal(segmentIterator);
    } else if (positionFlag > 0) {
        int vertex = positionFlag - 1;
        pseudoNormal = computeSegmentVertexNormal(segmentIterator, vertex, false);
    } else if (positionFlag < 0) {
        int edge = (- positionFlag) - 1;
        pseudoNormal = computeSegmentEdgeNormal(segmentIterator, edge);
    }

    return pseudoNormal;
}

/*!
 * Compute the surface-normal at specified point of the given triangle.
 *
 * Surface-normal is evaluated interpolating the unlimited vertex normals at
 * the specified point.
 *
 * To reduce computational times, normals of vertices are cached.
 *
 * @param[in] segmentIterator is an iterator pointing to the segment
 * @param[in] lambda are the barycentric coordinates of the point
 * @return the surface-normal at specified point of the given triangle
 */
std::array<double,3> SegmentationKernel::computeSurfaceNormal( const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const {

    // Early return if the segment is a point
    const Cell &segment = *segmentIterator;
    ElementType segmentType = segment.getType();
    if (segmentType == ElementType::VERTEX) {
        return {{0., 0., 0.}};
    }

    // Evaluate surface normal
    std::size_t nSegmentVertices = segment.getVertexCount();
    std::array<double,3> surfaceNormal = lambda[0] * computeSegmentVertexNormal(segmentIterator, 0, true);
    for (std::size_t i = 1; i < nSegmentVertices; ++i) {
        surfaceNormal += lambda[i] * computeSegmentVertexNormal(segmentIterator, i, true);
    }
    surfaceNormal /= norm2(surfaceNormal);

    return surfaceNormal;
}

/*!
 * Compute the normal of the specified triangle.
 *
 * To reduce computational times, normals of vertices are cached.
 *
 * @param[in] segmentIterator is an iterator pointing to the segment
 * @return the normal of the specified triangle
 */
std::array<double,3> SegmentationKernel::computeSegmentNormal( const SurfUnstructured::CellConstIterator &segmentIterator ) const {

    std::size_t segmentRawId = segmentIterator.getRawIndex();
    std::array<double, 3> *segmentNormal = m_segmentNormalsStorage.rawData(segmentRawId);
    if (!m_segmentNormalsValid.rawAt(segmentRawId)) {
        *segmentNormal = m_surface->evalFacetNormal(segmentIterator->getId());
        m_segmentNormalsValid.rawAt(segmentRawId) = true;
    }

    return *segmentNormal;
}

/*!
 * Compute the normal of the specified triangle's edge.
 *
 * To reduce computational times, normals of vertices are cached.
 *
 * @param[in] segmentIterator is an iterator pointing to the segment
 * @param[in] edge is the local index of the edge
 * @return the normal of the specified triangle's edge
 */
std::array<double,3> SegmentationKernel::computeSegmentEdgeNormal( const SurfUnstructured::CellConstIterator &segmentIterator, int edge ) const {

    long neighId = segmentIterator->getAdjacency(edge);
    SurfUnstructured::CellConstIterator neighIterator = m_surface->getCellConstIterator(neighId);

    std::array<double,3> normal = computeSegmentNormal(segmentIterator);
    normal += computeSegmentNormal(neighIterator);
    normal /= norm2(normal);

    return normal;
}

/*!
 * Compute the normal of the specified triangle's vertex.
 *
 * To reduce computational times, normals of vertices are cached.
 *
 * @param[in] segmentIterator is an iterator pointing to the segment
 * @param[in] vertex is the local index of the vertex
 * @param[in] limited controls is the limited or the unlimited normal will
 * be evaluated
 * @return the normal of the specified triangle's vertex
 */
std::array<double,3> SegmentationKernel::computeSegmentVertexNormal( const SurfUnstructured::CellConstIterator &segmentIterator, int vertex, bool limited ) const {

    // Segment information
    long segmentId = segmentIterator.getId();
    long segmentRawId = segmentIterator.getRawIndex();
    const Cell &segment = *segmentIterator;

    // Update the cache
    //
    // Cache is updated for all the vertices of the segment.
    long vertexId = segment.getVertexId(vertex);
    std::size_t vertexRawId = m_surface->getVertices().getRawIndex(vertexId);
    bool hasUnlimitedNormal = m_unlimitedVertexNormalsValid.rawAt(vertexRawId);

    bool hasLimitedNormal = m_limitedSegmentVertexNormalValid[m_segmentVertexOffset.rawAt(segmentRawId) + vertex];

    if (!hasUnlimitedNormal || !hasLimitedNormal) {
        static std::vector<long> vertexNeighbours;
        vertexNeighbours.clear();
        m_surface->findCellVertexNeighs(segmentId, vertex, &vertexNeighbours);

        std::array<double, 3> limitedVertexNormal;
        std::array<double, 3> unlimitedVertexNormal;
        if (hasUnlimitedNormal) {
            limitedVertexNormal   = m_surface->evalLimitedVertexNormal(segmentId, vertex, vertexNeighbours.size(), vertexNeighbours.data(), m_featureAngle) ;
            unlimitedVertexNormal = m_unlimitedVertexNormalsStorage.rawAt(vertexRawId);
        } else {
            m_surface->evalVertexNormals(segmentId, vertex, vertexNeighbours.size(), vertexNeighbours.data(), m_featureAngle, &unlimitedVertexNormal, &limitedVertexNormal) ;
        }

        // Store vertex limited normal
        //
        // Both limited and unlimited normal are evaluated, however limited
        // normal is only stored if its misalignment with respect to the
        // unlimited normal is greater than a defined tolerance.
        if( !hasLimitedNormal ){
            double misalignment = norm2(unlimitedVertexNormal - limitedVertexNormal) ;
            if( misalignment >= m_surface->getTol() ){
                std::pair<long, int> segmentVertexKey = std::make_pair(segmentId, vertex);
                m_limitedSegmentVertexNormalStorage.insert({segmentVertexKey, std::move(limitedVertexNormal)}) ;
            }
            m_limitedSegmentVertexNormalValid[m_segmentVertexOffset.rawAt(segmentRawId) + vertex] = true;
        }

        // Store vertex unlimited normal
        if ( !hasUnlimitedNormal ) {
            m_unlimitedVertexNormalsStorage.rawAt(vertexRawId) = std::move(unlimitedVertexNormal);
            m_unlimitedVertexNormalsValid.rawAt(vertexRawId) = true ;
        }
    }

    // Get the normal from the cache
    if (limited) {
        std::pair<long, int> segmentVertexKey = std::make_pair(segmentId, vertex);
        auto limitedSegmentVertexNormal = m_limitedSegmentVertexNormalStorage.find(segmentVertexKey);
        if (limitedSegmentVertexNormal != m_limitedSegmentVertexNormalStorage.end()) {
            return limitedSegmentVertexNormal->second;
        }
    }

    return m_unlimitedVertexNormalsStorage.rawAt(vertexRawId);
}

/*!
 *  \class      SurfaceInfo
 *  \ingroup    levelset
 *  \brief      Stores information regarding the projection points, like support element and surface normal
*/

/*!
 * Default constructor
 */
LevelSetSegmentation::SurfaceInfo::SurfaceInfo( ) : support(levelSetDefaults::SUPPORT), normal(levelSetDefaults::GRADIENT) 
{
}

/*!
 * Constructor
 */
LevelSetSegmentation::SurfaceInfo::SurfaceInfo( long index, const std::array<double,3> &vec ) : support(index), normal(vec)
{
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
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
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
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
void LevelSetSegmentation::setSegmentation( const SurfUnstructured *surface, double featureAngle){

    m_segmentation = std::make_shared<const SegmentationKernel>(surface, featureAngle);
}

/*!
 * Set the segmentation
 * @param[in,out] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
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
int LevelSetSegmentation::getPart( long id ) const{

    long supportId = getSupport(id);

    if( supportId != levelSetDefaults::SUPPORT){
        const SurfUnstructured &m_surface = m_segmentation->getSurface();
        return m_surface.getCell(supportId).getPID();
    } else { 
        return levelSetDefaults::PART ;
    }

}

/*!
 * Gets the surface normal at the projection point
 * @param[in] id index of cell
 * @return surface normal
 */
std::array<double,3> LevelSetSegmentation::getNormal( long id ) const{

    auto itr = m_surfaceInfo.find(id) ;
    if( itr != m_surfaceInfo.end() ){
        return itr->normal;
    } else {
        return levelSetDefaults::GRADIENT ;
    }

}


/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
long LevelSetSegmentation::getSupport( long id ) const{

    auto itr = m_surfaceInfo.find(id) ;
    if( itr != m_surfaceInfo.end() ){
        return itr->support;
    } else {
        return levelSetDefaults::SUPPORT ;
    }

}

/*!
 * Get size of support triangle
 * @param[in] i cell index
 * @return charcteristic size of support triangle
 */
double LevelSetSegmentation::getSurfaceFeatureSize( long i ) const {

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

    m_surfaceInfo.clear() ;
}

/*!
 * Computes the levelset function within the narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetSegmentation::computeLSInNarrowBand(bool signd){

    log::cout() << "Computing levelset within the narrow band... " << std::endl;

    if( LevelSetCartesian* lsCartesian = dynamic_cast<LevelSetCartesian*>(m_kernelPtr) ){
        computeLSInNarrowBand( lsCartesian, signd) ;

    } else if ( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        computeLSInNarrowBand( lsOctree, signd) ;

    }
}

/*!
 * Updates the levelset function within the narrow band after mesh adaptation.
 * @param[in] mapper information concerning mesh adaption 
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetSegmentation::updateLSInNarrowBand( const std::vector<adaption::Info> &mapper, bool signd){

    log::cout() << "Updating levelset within the narrow band... " << std::endl;
    if( LevelSetCartesian* lsCartesian= dynamic_cast<LevelSetCartesian*>(m_kernelPtr) ){

        // Update is not implemented for Cartesian patches
        clear( ) ;
        computeLSInNarrowBand( lsCartesian, signd) ;
        return;
    }

    if( LevelSetOctree* lsOctree = dynamic_cast<LevelSetOctree*>(m_kernelPtr) ){
        updateLSInNarrowBand( lsOctree, mapper, signd ) ;
        return;
    }


}

/*!
 * Computes the levelset within the narrow band on an
 * cartesian grid.
 * If the size of the narrow band has been set, the
 * method will compute the levelset values only of those
 * cells within the threshold. 
 * In case the size of the narrow band has not been set,
 * the method will calculate the levelset within a band
 * containing one cell on each side of the surface.
 * @param[in] visitee the octree LevelSetKernel
 * @param[in] signd whether signed distance should be calculated
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetCartesian *visitee, bool signd){

    log::cout() << " Compute levelset on cartesian mesh"  << std::endl;

    // Get mesh information
    VolCartesian &mesh = *(visitee->getCartesianMesh() ) ;
    int meshDimension = mesh.getDimension();

    // Get surface information
    const SurfUnstructured &surface = m_segmentation->getSurface();

    // Define search radius
    double searchRadius = m_narrowBand;
    if (searchRadius < 0.) {
        for (int d = 0; d < meshDimension; ++d) {
            searchRadius = std::max(searchRadius, mesh.getSpacing(d));
        }
    }

    // Define mesh bounding box
    //
    // The bounding box is inflated be the search radius.
    std::array<double,3> meshMinPoint;
    std::array<double,3> meshMaxPoint;
    mesh.getBoundingBox(meshMinPoint, meshMaxPoint) ;
    for (int d = 0; d < meshDimension; ++d) {
        meshMinPoint[d] -= searchRadius;
        meshMaxPoint[d] += searchRadius;
    }

    // Initialize process list
    //
    // Process list is initialized with cells that are certainly inside the
    // narrow band. Those cells are the one that contain the vertices of the
    // segments or the intersection between the segments and the bounding box
    // of the patch.
    std::unordered_set<long> processList;

    std::vector<std::array<double,3>> intersectionPoints;
    std::vector<std::array<double,3>> segmentVertexCoords;
    for (const Cell &segment : surface.getCells()) {
        // Get segment info
        ConstProxyVector<long> segmentVertexIds = segment.getVertexIds();
        std::size_t nSegmentVertices = segmentVertexIds.size();

        // Get segment coordinates
        segmentVertexCoords.resize(nSegmentVertices);
        surface.getVertexCoords(nSegmentVertices, segmentVertexIds.data(), segmentVertexCoords.data());

        // Add to the process list the cells that contain the vertices of the
        // segment or the intersection between the segment and the bounding box
        // of the patch.
        int nInnerVertices = 0;
        for (const std::array<double,3> &vertexPoint : segmentVertexCoords) {
            long cellId = mesh.locatePoint(vertexPoint);
            if (cellId < 0) {
                continue;
            }

            processList.insert(cellId);
            ++nInnerVertices;
        }

        if (nInnerVertices == 0) {
            if (CGElem::intersectBoxPolygon(meshMinPoint, meshMaxPoint, segmentVertexCoords, false, true, true, intersectionPoints, meshDimension)) {
                for (const std::array<double,3> &intersectionPoint : intersectionPoints){
                    long cellId = mesh.locateClosestCell(intersectionPoint);
                    assert(cellId >= 0);
                    processList.insert(cellId);
                }
            }
        }
    }

    // Evaluate the levelset within the narrow band
    //
    // The initial process list is gradually expanded considering all the
    // neighbours with a distance less than the search radius.
    std::unordered_set<long> alreadyProcessed;
    while (!processList.empty()) {
        // Get the cell to process
        long cellId = *(processList.begin());
        processList.erase(cellId);
        alreadyProcessed.insert(cellId);

        // Find segment associated to the cell
        std::array<double,3> cellCentroid = visitee->computeCellCentroid(cellId);

        long segmentId;
        double distance;
        m_segmentation->getSearchTree().findPointClosestCell(cellCentroid, searchRadius, &segmentId, &distance);
        if(segmentId < 0){
            continue;
        }

        // Evaluate levelset information
        std::array<double, 3> gradient;
        std::array<double, 3> normal;
        int error = m_segmentation->getSegmentInfo(cellCentroid, segmentId, signd, distance, gradient, normal);
        if (error) {
            throw std::runtime_error ("Unable to extract the levelset information from segment.");
        }

        PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(cellId) ;
        lsInfoItr->value    = distance;
        lsInfoItr->gradient = gradient;

        PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId);
        infoItr->support = segmentId;
        infoItr->normal = normal;

        // Add cell neighbours to the process list
        const Cell &cell = mesh.getCell(cellId);
        const long *neighbours = cell.getAdjacencies() ;
        int nNeighbours = cell.getAdjacencyCount() ;
        for (int n = 0; n < nNeighbours; ++n) {
            long neighId = neighbours[n];
            if (alreadyProcessed.count(neighId) == 0) {
                processList.insert(neighId);
            }
        }
    }
}

/*!
 * Computes the levelset within the narrow band on an
 * octree grid.
 * If the size of the narrow band has been set, the
 * method will compute the levelset values only of those
 * cells within the threshold. 
 * In case the size of the narrow band has not been set,
 * the method will calculate the levelset within the cells
 * that intersect the surface and within their first neighbours,
 * \param[in] visitee the octree LevelSetKernel
 * \param[in] signd whether signed distance should be calculated
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetOctree *visitee, bool signd){

    VolumeKernel &mesh = *(visitee->getMesh()) ;

    double intersectionFactor = 0.5 * std::sqrt((double) mesh.getDimension());

    std::unordered_set<long> intersectedCells;

    // Evaluate levelset information
    for( const Cell &cell : mesh.getCells() ){

        // Identify the segment associated with the cell
        //
        // The search radius is evaluated as the maximum value between the
        // narroband size and the distance above which the cell will surely
        // not intersect the surface. In this way, cells that intersect the
        // surface are always included in the narrowband, even if their
        // distance from the surface is greater than then narrowband size
        // explicitly set by the user.
        //
        // If no segment is identified the cell is not processed.
        long cellId = cell.getId();
        std::array<double,3> cellCentroid = visitee->computeCellCentroid(cellId);

        double intersectDistance = intersectionFactor * mesh.evalCellSize(cellId);
        double searchRadius = std::max(m_narrowBand, intersectDistance);

        long segmentId;
        double distance;
        m_segmentation->getSearchTree().findPointClosestCell(cellCentroid, searchRadius, &segmentId, &distance);
        if(segmentId < 0){
            continue;
        }

        // Evaluate levelset information
        std::array<double,3> gradient;
        std::array<double,3> normal;
        int error = m_segmentation->getSegmentInfo(cellCentroid, segmentId, signd, distance, gradient, normal);
        if (error) {
            throw std::runtime_error ("Unable to extract the levelset information from segment.");
        }

        PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.find(cellId);
        lsInfoItr = m_ls.emplace(cellId) ;
        lsInfoItr->value    = distance;
        lsInfoItr->gradient = gradient;

        PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId);
        infoItr->support = segmentId;
        infoItr->normal = normal;

        // Update the list of cells that intersects the surface
        //
        // When the narrowband size is not explicitly set, the cell will always
        // intersects the surface because only cells that intersect the surface
        // are considered, otherwise we need to check if the absolute distance
        // associated with the cell is lower than the intersection distance.
        if (m_narrowBand < 0 || intersectDistance < std::abs(lsInfoItr->value)) {
            intersectedCells.insert(cellId);
        }

    }

    // Process the neighbours of the cells that intersect the surface
    //
    // If a cell intersects the surface, we need to evaluate the levelset
    // of all its neigbours.
    for( long cellId : intersectedCells){

        Cell const &cell = mesh.getCell(cellId);
        
        std::array<double,3> cellProjectionPoint = computeProjectionPoint(cellId);

        const long *neighbours = cell.getAdjacencies() ;
        int nNeighbours = cell.getAdjacencyCount() ;
        for (int n = 0; n < nNeighbours; ++n) {
            // Skip the neighbour if it has already been processed
            //
            // The neighbour may already have been processed either because
            // it distance from the segmentation is within the search radius,
            // or because is a neighbour of an intersected cells already
            // processed.
            long neighId = neighbours[n];
            if( m_ls.count(neighId) != 0 ){
                continue;
            }

            // Identify the segment associated with the neighbour
            std::array<double,3> neighCentroid = visitee->computeCellCentroid(neighId);

            double searchRadius = 1.05 * norm2(neighCentroid - cellProjectionPoint);

            long segmentId;
            double distance;
            m_segmentation->getSearchTree().findPointClosestCell(neighCentroid, searchRadius, &segmentId, &distance);
            if (segmentId < 0) {
                assert(false && "Should not pass here");
            }

            // Evaluate negihbour leveset information
            std::array<double,3> gradient;
            std::array<double,3> normal;
            int error = m_segmentation->getSegmentInfo(neighCentroid, segmentId, signd, distance, gradient, normal);
            if (error) {
                throw std::runtime_error ("Unable to extract the levelset information from segment.");
            }

            PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(neighId) ;
            lsInfoItr->value    = distance;
            lsInfoItr->gradient = gradient;

            PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(neighId);
            infoItr->support = segmentId;
            infoItr->normal = normal;
        }
    }
}

/*!
 * Updates the levelset within the narrow band on an
 * octree grid after an grid adaption.
 * If the size of the narrow band has been set, the
 * method will compute the levelset values only of those
 * cells within the threshold. 
 * In case the size of the narrow band has not been set,
 * the method will calculate the levelset within the cells
 * that intersect the surface and within their first neighbours,
 * @param[in] visitee the octree LevelSetKernel
 * @param[in] mapper the adaption mapper
 * @param[in] signd whether signed distance should be calculated
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetOctree *visitee, const std::vector<adaption::Info> &mapper, bool signd){

    VolumeKernel &mesh = *(visitee->getMesh()) ;

    double intersectionFactor = 0.5 * std::sqrt((double) mesh.getDimension());

    std::vector<long> cellsOutsideNarrowband;

    // Evaluate the levelset of the cells
    //
    // When searching for the segment associated to a cell, the search radius
    // is evaluated as the maximum value between the narroband size and the
    // distance above which the cell will surely not intersect the surface.
    for( const auto &event : mapper ){

        if( event.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for( long cellId : event.current ){

            // Identify the segment associated with the cell
            //
            // The search radius is evaluated as the maximum value between the
            // narroband size and the distance above which the cell will surely
            // not intersect the surface. In this way, cells that intersect the
            // surface are always included in the narrowband, even if their
            // distance from the surface is greater than then narrowband size
            // explicitly set by the user.
            //
            // If no segment is identified the cell is not processed.
            std::array<double,3> centroid = visitee->computeCellCentroid(cellId);

            double intersectDistance = intersectionFactor * mesh.evalCellSize(cellId);
            double searchRadius = std::max(m_narrowBand, intersectDistance);

            long segmentId;
            double distance;
            m_segmentation->getSearchTree().findPointClosestCell(centroid, searchRadius, &segmentId, &distance);
            if (segmentId < 0) {
                cellsOutsideNarrowband.push_back(cellId);
                continue;
            }

            // Evaluate levelset information
            std::array<double,3> gradient;
            std::array<double,3> normal;
            int error = m_segmentation->getSegmentInfo(centroid, segmentId, signd, distance, gradient, normal);
            if (error) {
                throw std::runtime_error ("Unable to extract the levelset information from segment.");
            }

            PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(cellId) ;
            lsInfoItr->value    = distance;
            lsInfoItr->gradient = gradient;

            PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId) ;
            infoItr->support = segmentId;
            infoItr->normal = normal;
        }

    }

    // Cells with neighbours that intersect the surface need to be added to
    // the narrowband even if they don't intersect the surface themself or
    // have a distance from the surface greater than the narroband size.
    for( long cellId : cellsOutsideNarrowband){
        const Cell &cell = mesh.getCell(cellId);

        // Consider only cells with a neighbour that intersects the surface
        const long *neighbours = cell.getAdjacencies() ;
        int nNeighbours = cell.getAdjacencyCount() ;

        long intersectedNeighId = Cell::NULL_ID;
        for (int n = 0; n < nNeighbours; ++n) {
            long neighId = neighbours[n];
            if( intersectSurface(neighId,LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) == LevelSetIntersectionStatus::TRUE){
                intersectedNeighId = neighId;
                break;
            }
        }

        if (intersectedNeighId == Cell::NULL_ID) {
            continue;
        }

        // Identify the segment associated with the cell
        std::array<double,3> cellCentroid = visitee->computeCellCentroid(cellId);
        std::array<double,3> neighProjectionPoint = computeProjectionPoint(intersectedNeighId);

        double searchRadius = 1.05 * norm2(cellCentroid - neighProjectionPoint);

        long segmentId;
        double distance;
        m_segmentation->getSearchTree().findPointClosestCell(cellCentroid, searchRadius, &segmentId, &distance);
        if (segmentId < 0) {
            assert(false && "Should not pass here");
            continue;
        }

        // Evaluate levelset information for the cell
        std::array<double,3> gradient;
        std::array<double,3> normal;
        int error = m_segmentation->getSegmentInfo(cellCentroid, segmentId, signd, distance, gradient, normal);
        if (error) {
            throw std::runtime_error ("Unable to extract the levelset information from segment.");
        }

        PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(cellId) ;
        lsInfoItr->value    = distance;
        lsInfoItr->gradient = gradient;

        PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId);
        infoItr->support = segmentId;
        infoItr->normal = normal;
    }
}

/*!
 * Prune the segment's info removing entries associated to cells that are
 * are not in the mesh anymore
 * @param[in] mapper information concerning mesh adaption
 */
void LevelSetSegmentation::__clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

    for ( const auto &info : mapper ){
        // Consider only changes on cells
        if( info.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        // Delete only old data that belongs to the current process
        if (info.type == adaption::Type::TYPE_PARTITION_RECV) {
            continue;
        }

        // Remove info of previous cells
        for ( long  parent : info.previous ) {
            if ( m_surfaceInfo.find( parent ) != m_surfaceInfo.end() ) {
                m_surfaceInfo.erase( parent, true ) ;
            }
        }
    }

    m_surfaceInfo.flush();
}

/*!
 * Writes LevelSetSegmentation to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetSegmentation::__dump( std::ostream &stream ){

    utils::binary::write( stream, m_surfaceInfo.size() ) ;

    bitpit::PiercedVector<SurfaceInfo>::iterator infoItr, infoEnd = m_surfaceInfo.end() ;
    for( infoItr = m_surfaceInfo.begin(); infoItr != infoEnd; ++infoItr){
        utils::binary::write( stream, infoItr.getId() );
        utils::binary::write( stream, infoItr->support );
        utils::binary::write( stream, infoItr->normal );
    }
}

/*!
 * Reads LevelSetSegmentation from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetSegmentation::__restore( std::istream &stream ){

    size_t size;

    long id;
    long support;
    std::array<double,3> normal;

    utils::binary::read( stream, size ) ;
    m_surfaceInfo.reserve(size);

    for( size_t i=0; i<size; ++i){
        utils::binary::read( stream, id );

        utils::binary::read( stream, support );
        utils::binary::read( stream, normal );

        m_surfaceInfo.insert(id, SurfaceInfo(support,normal));
    }
}

# if BITPIT_ENABLE_MPI

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetSegmentation::__writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    // Evaluate the size of the buffer
    long nItems = 0;
    for( const auto &index : sendList){
        auto infoItr = m_surfaceInfo.find(index) ;
        if( infoItr != m_surfaceInfo.end() ){
            nItems++ ;
        }
    }

    dataBuffer.setSize(dataBuffer.getSize() + sizeof(long) + nItems* ( 2*sizeof(long) +3*sizeof(double) ));

    // Fill the buffer
    dataBuffer << nItems ;

    long index = 0 ;
    for( long id : sendList){
        auto infoItr = m_surfaceInfo.find(id) ;
        if( infoItr != m_surfaceInfo.end() ){
            dataBuffer << index ;
            dataBuffer << infoItr->support;
            dataBuffer << infoItr->normal;
        }

        ++index;
    }
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetSegmentation::__readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    long nItems ;
    dataBuffer >> nItems ;

    for( long i=0; i<nItems; ++i){
        // Get the id of the element
        long index;
        dataBuffer >> index;
        long id = recvList[index] ;

        // Assign the data of the element
        auto infoItr = m_surfaceInfo.find(id) ;
        if( infoItr == m_surfaceInfo.end() ){
            infoItr = m_surfaceInfo.emplace(id) ;
        }
        dataBuffer >> infoItr->support;
        dataBuffer >> infoItr->normal;
    }
}
# endif

/*!
 * Computes the LevelSetInfo of a point
 * \param[in] coords coordinates of the point
 * \return the LevelSetInfo
 */
LevelSetInfo LevelSetSegmentation::computeLevelSetInfo(const std::array<double,3> &coords) const {

    long segmentId;
    double distance;
    std::array<double,3> gradient;
    std::array<double,3> normal;

    m_segmentation->getSearchTree().findPointClosestCell(coords, &segmentId, &distance);

    int error = m_segmentation->getSegmentInfo(coords, segmentId, false, distance, gradient, normal);
    if (error) {
        throw std::runtime_error ("Unable to extract the levelset information from segment.");
    }

    return LevelSetInfo(distance,gradient);

}

}
