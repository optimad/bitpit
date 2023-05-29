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

# include "levelSetSegmentationObject.hpp"

namespace bitpit {

/*!
    @class      LevelSetSegmentationKernel
    @ingroup    levelset
    @brief      Segmentation kernel
*/

/*!
 * Default constructor
 */
LevelSetSegmentationKernel::LevelSetSegmentationKernel()
    : m_surface(nullptr),
      m_featureAngle(0)
{
}

/*!
 * Copy constructor
 */
LevelSetSegmentationKernel::LevelSetSegmentationKernel(const LevelSetSegmentationKernel &other)
    : m_surface(other.m_surface),
      m_featureAngle(other.m_featureAngle),
      m_segmentVertexOffset(other.m_segmentVertexOffset),
      m_segmentNormalsValid(other.m_segmentNormalsValid),
      m_segmentNormalsStorage(other.m_segmentNormalsStorage),
      m_unlimitedVertexNormalsValid(other.m_unlimitedVertexNormalsValid),
      m_unlimitedVertexNormalsStorage(other.m_unlimitedVertexNormalsStorage),
      m_limitedSegmentVertexNormalValid(other.m_limitedSegmentVertexNormalValid),
      m_limitedSegmentVertexNormalStorage(other.m_limitedSegmentVertexNormalStorage)
{
    if (other.m_ownedSurface) {
        m_ownedSurface = std::unique_ptr<SurfUnstructured>(new SurfUnstructured(*(other.m_ownedSurface)));
    } else {
        m_ownedSurface = nullptr;
    }

    if (other.m_searchTree) {
        m_searchTree = std::unique_ptr<SurfaceSkdTree>(new SurfaceSkdTree(m_surface));
        m_searchTree->build();
    } else {
        m_searchTree = nullptr;
    }
}

/*!
 * Move constructor
 */
LevelSetSegmentationKernel::LevelSetSegmentationKernel(LevelSetSegmentationKernel &&other)
    : m_surface(std::move(other.m_surface)),
      m_ownedSurface(std::move(other.m_ownedSurface)),
      m_featureAngle(std::move(other.m_featureAngle)),
      m_searchTree(std::move(other.m_searchTree)),
      m_segmentVertexOffset(std::move(other.m_segmentVertexOffset)),
      m_segmentNormalsValid(std::move(other.m_segmentNormalsValid)),
      m_segmentNormalsStorage(std::move(other.m_segmentNormalsStorage)),
      m_unlimitedVertexNormalsValid(std::move(other.m_unlimitedVertexNormalsValid)),
      m_unlimitedVertexNormalsStorage(std::move(other.m_unlimitedVertexNormalsStorage)),
      m_limitedSegmentVertexNormalValid(std::move(other.m_limitedSegmentVertexNormalValid)),
      m_limitedSegmentVertexNormalStorage(std::move(other.m_limitedSegmentVertexNormalStorage))
{
}

/*!
 * Constructor
 *
 * @param[in,out] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
LevelSetSegmentationKernel::LevelSetSegmentationKernel(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle ) {

    setSurface(std::move(surface), featureAngle);
}

/*!
 * Constructor
 *
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
LevelSetSegmentationKernel::LevelSetSegmentationKernel(const SurfUnstructured *surface, double featureAngle ) {

    setSurface(surface, featureAngle);
}

/*!
 * Get feature angle
 * @return feature angle used when calculating face normals;
 */
double LevelSetSegmentationKernel::getFeatureAngle() const {
    return m_featureAngle;
}

/*!
 * Get segmentation surface
 * @return segmentation surface;
 */
const SurfUnstructured & LevelSetSegmentationKernel::getSurface() const {
    return *m_surface;
}

/*!
 * Set the surface
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
void LevelSetSegmentationKernel::setSurface(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle){

    m_ownedSurface = std::move(surface);
    setSurface(m_ownedSurface.get(), featureAngle);
}

/*!
 * Set the surface
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
void LevelSetSegmentationKernel::setSurface(const SurfUnstructured *surface, double featureAngle){

    // Check if adjacencies are built
    if (surface->getAdjacenciesBuildStrategy() == SurfUnstructured::ADJACENCIES_NONE) {
        throw std::runtime_error ("Segmentation needs adjacencies!") ;
    }

    // Surface information
    m_surface      = surface;
    m_featureAngle = featureAngle;

    // Segment vertices information
    m_segmentVertexOffset.setStaticKernel(&m_surface->getCells());

    std::size_t nTotalSegmentVertices = 0;
    for (auto itr = m_segmentVertexOffset.begin(); itr != m_segmentVertexOffset.end(); ++itr) {
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
 * Get search tree
 * @return search tree;
 */
const SurfaceSkdTree & LevelSetSegmentationKernel::getSearchTree() const {
    return *m_searchTree;
}

/*!
 * Evaluate levelset information at the specified point.
 *
 * @param[in] point coordinates of point
 * @param[in] support support associated with the point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @param[out] s if a valid pointer is provided, on output it will contain the sign of the levelset
 * @param[out] value if a valid pointer is provided, on output it will contain the value of the
 * levelset
 * @param[out] value if a valid pointer is provided, on output it will contain the gradient of
 * the levelset
 */
int LevelSetSegmentationKernel::evalLevelsetInfo(const std::array<double,3> &point, long support, bool signedLevelSet,
                                                 double *value, std::array<double,3> *gradient) const
{
    // Segment information
    SurfUnstructured::CellConstIterator segmentIterator = m_surface->getCellConstIterator(support);
    const Cell &segment = *segmentIterator ;
    ElementType segmentType = segment.getType();
    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds() ;
    int nSegmentVertices = segmentVertexIds.size() ;

    // Project the point on the surface and evaluate the point-projection vector
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    std::array<double,3> pointProjectionVector = point;
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
        pointProjectionVector -= CGElem::projectPointSegment(point, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), lambda);

        break;
    }

    case ElementType::TRIANGLE:
    {
        long id0 = segmentVertexIds[0] ;
        long id1 = segmentVertexIds[1] ;
        long id2 = segmentVertexIds[2] ;
        pointProjectionVector -= CGElem::projectPointTriangle(point, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), m_surface->getVertexCoords(id2), lambda );

        break;
    }

    default:
    {
        ConstProxyVector<long> elementVertexIds = m_surface->getFacetOrderedVertexIds(segment);
        BITPIT_CREATE_WORKSPACE(segmentVertexCoors, std::array<double BITPIT_COMMA 3>, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
        m_surface->getVertexCoords(elementVertexIds.size(), elementVertexIds.data(), segmentVertexCoors);
        pointProjectionVector -= CGElem::projectPointPolygon(point, nSegmentVertices, segmentVertexCoors, lambda );

        break;
    }

    }

    // Evaluate distance from surface
    double pointProjectionDistance = norm2(pointProjectionVector);

    // Check if the point lies on the segmentation
    //
    // If the distance is zero, the point and the projection are coincident,
    // this means that the point lies on the segmentation.
    bool pointOnSegmentation = false;
    if (signedLevelSet || gradient) {
        double distanceTolerance = m_surface->getTol();
        pointOnSegmentation = utils::DoubleFloatingEqual()(pointProjectionDistance, 0., distanceTolerance, distanceTolerance);
    }

    // Evaluate point sign
    //
    // The sign is computed by determining the side of the point with respect
    // to the normal plane. The sign will be zero if the point lies exactly
    // on the segmentation or on the normal plane. In the latter case the sign
    // must be evaluated taking into account the the curvature of the surface.
    // However, this is not yet implemented.
    int pointSign;
    if (signedLevelSet) {
        std::array<double, 3> pseudoNormal = computePseudoNormal(segmentIterator, lambda);
        pointSign = sign(dotProduct(pointProjectionVector, pseudoNormal) );
        if (!pointOnSegmentation && pointSign == 0) {
            return 1;
        }
    } else {
        pointSign = 1;
    }

    // Evaluate levelset value
    if (value) {
        *value = pointSign * pointProjectionDistance;
    }

    // Evaluate levelset gradient
    if (gradient) {
        if (!pointOnSegmentation) {
            *gradient = pointProjectionVector / (pointSign * pointProjectionDistance);
        } else {
            *gradient = static_cast<double>(pointSign) * computeSurfaceNormal(segmentIterator, lambda);
        }
    }

    return 0;
}

/*!
 * Get the size of a segment
 * @param[in] segmentId is the id of the segment
 * @return charcteristic size of the segment
 */
double LevelSetSegmentationKernel::getSegmentSize(long segmentId) const {

    int surfaceDimension = m_surface->getDimension();
    if (surfaceDimension == 1) {
        return m_surface->evalCellArea(segmentId); //TODO check
    } else if (surfaceDimension == 2) {
        int dummy;
        return m_surface->evalMinEdgeLength(segmentId, dummy);
    }

    return (- levelSetDefaults::SIZE);
}

/*!
 * Get the size of the smallest segment
 * @return the size of the smallest segment
 */
double LevelSetSegmentationKernel::getMinSegmentSize() const {

    bool   minimumValid = false;
    double minimumSize  = levelSetDefaults::SIZE;
    for (const Cell &cell : m_surface->getCells()) {
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
 * Get the size of the largest segment
 * @return the size of the largest segment
 */
double LevelSetSegmentationKernel::getMaxSegmentSize() const {

    double maximumSize = - levelSetDefaults::SIZE;
    for (const Cell &cell : m_surface->getCells()) {
        double segmentSize = getSegmentSize(cell.getId());
        maximumSize = std::max(segmentSize, maximumSize);
    }

    return maximumSize;
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
std::array<double,3> LevelSetSegmentationKernel::computePseudoNormal(const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const {

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
std::array<double,3> LevelSetSegmentationKernel::computeSurfaceNormal(const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const {

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
std::array<double,3> LevelSetSegmentationKernel::computeSegmentNormal(const SurfUnstructured::CellConstIterator &segmentIterator ) const {

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
std::array<double,3> LevelSetSegmentationKernel::computeSegmentEdgeNormal(const SurfUnstructured::CellConstIterator &segmentIterator, int edge ) const {

    std::array<double,3> normal = computeSegmentNormal(segmentIterator);

    if (segmentIterator->getAdjacencyCount(edge) > 0) {
        long neighId = segmentIterator->getAdjacency(edge);
        SurfUnstructured::CellConstIterator neighIterator = m_surface->getCellConstIterator(neighId);

        normal += computeSegmentNormal(neighIterator);
        normal /= norm2(normal);
    }

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
std::array<double,3> LevelSetSegmentationKernel::computeSegmentVertexNormal(const SurfUnstructured::CellConstIterator &segmentIterator, int vertex, bool limited ) const {

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
        if (!hasLimitedNormal) {
            double misalignment = norm2(unlimitedVertexNormal - limitedVertexNormal) ;
            if (misalignment >= m_surface->getTol()) {
                std::pair<long, int> segmentVertexKey = std::make_pair(segmentId, vertex);
                m_limitedSegmentVertexNormalStorage.insert({segmentVertexKey, std::move(limitedVertexNormal)}) ;
            }
            m_limitedSegmentVertexNormalValid[m_segmentVertexOffset.rawAt(segmentRawId) + vertex] = true;
        }

        // Store vertex unlimited normal
        if (!hasUnlimitedNormal ) {
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
	@class      LevelSetSegmentation
	@ingroup    levelset
	@brief      Implements visitor pattern fo segmentated geometries
*/

/*!
 * Constructor
 * @param[in] id identifier of object
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id)
    : LevelSetObject(id),
      LevelSetSegmentationKernel()
{
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL unique pointer to surface mesh
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id, std::unique_ptr<const SurfUnstructured> &&STL, double featureAngle)
    : LevelSetObject(id),
      LevelSetSegmentationKernel(std::move(STL), featureAngle)
{
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] STL pointer to surface mesh
 * @param[in] featureAngle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id, const SurfUnstructured *STL, double featureAngle)
    : LevelSetObject(id),
      LevelSetSegmentationKernel(STL, featureAngle)
{
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetSegmentationObject * LevelSetSegmentationObject::clone() const {
    return new LevelSetSegmentationObject(*this );
}

/*!
 * Set the cache to be used for the specified field.
 *
 * If an existing cache is already defined for the specified field, it will be destroyed and
 * recreated from scratch.
 *
 * \param field is the field for which cache mode will be set
 * \param cacheMode is the type of cache that will be used for caching field information
 */
void LevelSetSegmentationObject::setFieldCache(LevelSetField field, LevelSetCacheMode cacheMode)
{
    switch(field) {

    case LevelSetField::SUPPORT:
        unregisterFieldCellCache(field);
        registerFieldCellCache<long>(field, cacheMode);
        break;

    default:
        LevelSetObject::setFieldCache(field, cacheMode);

    }
}

/*!
 * Get the list of supported field.
 * @result The list of supported field.
 */
LevelSetFieldset LevelSetSegmentationObject::getSupportedFields() const
{
    LevelSetFieldset supportedFields = LevelSetObject::getSupportedFields();
    supportedFields.insert(LevelSetField::PART);
    supportedFields.insert(LevelSetField::NORMAL);
    supportedFields.insert(LevelSetField::SUPPORT);

    return supportedFields;
}

/*!
 * Get the smallest characteristic size within the triangulation
 * @return smallest characteristic size within the triangulation
 */
double LevelSetSegmentationObject::getMinSurfaceFeatureSize() const {

    return getMinSegmentSize();
}

/*!
 * Get the largest characteristic size within the triangulation.
 * @return largest characteristic size within the triangulation
 */
double LevelSetSegmentationObject::getMaxSurfaceFeatureSize() const {

    return getMaxSegmentSize();
}

/*!
 * Fill cell caches in "narrow band" mode.
 *
 * If the size of the narrow band has been set, the method will fill the caches on the cells
 * that intersect the surface, on all their first neighbours and on the cells with a distance
 * from the surface less than the defined narrow band size.
 *
 * In case the size of the narrow band has not been set, the method will fill the caches on
 * the cells that intersect the surface and on all their first neighbours.
 */
void LevelSetSegmentationObject::fillNarrowBandCellCaches()
{
    log::cout() << "Computing levelset within the narrow band... " << std::endl;

    // Cartesian patches are handled separately
    if (dynamic_cast<LevelSetCartesianKernel*>(m_kernel)) {
        fillCartesianNarrowBandCellCaches() ;
        return;
    }

    // All other patches are handled with the same method.
    LevelSetObject::fillNarrowBandCellCaches();
}

/*!
 * Fill cell caches in "full" mode after a mesh update.
 *
 * If the size of the narrow band has been set, the method will fill the caches on the cells
 * that intersect the surface, on all their first neighbours and on the cells with a distance
 * from the surface less than the defined narrow band size.
 *
 * In case the size of the narrow band has not been set, the method will fill the caches on
 * the cells that intersect the surface and on all their first neighbours.
 *
 * \param adaptionData are the information about the adaption
 */
void LevelSetSegmentationObject::fillNarrowBandCellCaches(const std::vector<adaption::Info> &adaptionData)
{
    log::cout() << "Updating levelset within the narrow band... " << std::endl;

    // Cartesian patches are handled separately
    //
    // Update is not implemented for Cartesian patches, the levelset is cleared and rebuild
    // from scratch.
    if (dynamic_cast<LevelSetCartesianKernel*>(m_kernel)) {
        // Clear caches
        LevelSetFieldset narrowBandCacheFieldset = getCachedFields(LevelSetCacheMode::NARROW_BAND);
        for (LevelSetField field : narrowBandCacheFieldset) {
            CellCache *cache = getFieldCellCache(field);
            cache->clear();
        }

        // Fill the caches from scratch
        fillCartesianNarrowBandCellCaches() ;

        return;
    }

    // All other patches are handled with the same method
    LevelSetObject::fillNarrowBandCellCaches(adaptionData);
}

/*!
 * Fill cell caches in "narrow band" mode.
 *
 * This function is
 *
 * If the size of the narrow band has been set, the method will fill the caches on the cells
 * that intersect the surface, on all their first neighbours and on the cells with a distance
 * from the surface less than the defined narrow band size.
 *
 * In case the size of the narrow band has not been set, the method will fill the caches on
 * the cells that intersect the surface and on all their first neighbours.
 */
void LevelSetSegmentationObject::fillCartesianNarrowBandCellCaches()
{
    // The function needs a Cartesian kernel
    const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<LevelSetCartesianKernel*>(m_kernel);
    if (!dynamic_cast<LevelSetCartesianKernel*>(m_kernel)) {
        throw std::runtime_error("The function needs a Cartesian kernels.");
    }

    // Get field to process
    LevelSetFieldset narrowBandCacheFieldset = getCachedFields(LevelSetCacheMode::NARROW_BAND);
    std::vector<LevelSetField> fieldProcessList(narrowBandCacheFieldset.begin(), narrowBandCacheFieldset.end());
    if (fieldProcessList.empty()) {
        return;
    }

    // Get narrow band cache
    CellValueCache<bool> *narrowBandCache = getCellCache<bool>(m_cellNarrowBandCacheId);

    // Get mesh information
    const VolCartesian &mesh = *(cartesianKernel->getMesh() ) ;
    int meshDimension = mesh.getDimension();
    VolCartesian::MemoryMode meshMemoryMode = mesh.getMemoryMode();

    ElementType meshCellType = mesh.getCellType();
    const ReferenceElementInfo &meshCellTypeInfo = ReferenceElementInfo::getInfo(meshCellType);
    int meshCellFaceCount = meshCellTypeInfo.nFaces;

    std::array<double,3> meshMinPoint;
    std::array<double,3> meshMaxPoint;
    mesh.getBoundingBox(meshMinPoint, meshMaxPoint) ;

    // Get surface information
    const SurfUnstructured &surface = getSurface();

    // Initialize process list
    //
    // Process list is initialized with cells that are certainly inside the
    // narrow band. Those cells are the ones that contain the vertices of the
    // segments or the intersection between the segments and the bounding box
    // of the patch.
    std::unordered_set<long> processList;

    std::vector<std::array<double,3>> intersectionPoints;
    std::vector<std::array<double,3>> segmentVertexCoords;
    for (const Cell &segment : surface.getCells()) {
        // Get segment info
        //
        // Since vertex information will be passed to the CG module, we need to get the
        // vertices in the proper order.
        ConstProxyVector<long> segmentVertexIds = surface.getFacetOrderedVertexIds(segment);
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
                    processList.insert(cellId);
                }
            }
        }
    }

    // Evaluate the levelset within the narrow band
    //
    // The initial process list is gradually expanded considering all the
    // neighbours with a distance less than the search radius.
    std::unordered_set<long> outsideSearchRadius;
    std::unordered_set<long> surfaceNeighbourhood;
    while (!processList.empty()) {
        // Get the cell to process
        long cellId = *(processList.begin());
        processList.erase(processList.begin());

        // Check if the cell is within the narrow band
        //
        // No neighbour check is performed, cells with neighbours that intersect the
        // zero-levelset iso-surface the surface will be processed later.
        //
        // The check should be performed ignoring the narrow band cache, because we are
        // in the process of building that cache.
        double maximumDistance;
        bool cellInsideNarrowBand = _isCellInNarrowBand(cellId, false, &maximumDistance);
        if (!cellInsideNarrowBand) {
            continue;
        }

        // Fill cell caches
        narrowBandCache->insertEntry(cellId, true);
        for (LevelSetField field : fieldProcessList) {
            fillFieldCellCache(cellId, field, maximumDistance);
        }

        // Check if the cell intersects the surface
        //
        // When the narrow band size is not explicitly set, the cell will always
        // intersects the surface because only cells that intersect the surface
        // are considered, otherwise we need to explicitly check if the cell
        // intersects the surface.
        bool isCellIntersected = false;
        if (m_narrowBandSize < 0 || intersectSurface(cellId, LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) == LevelSetIntersectionStatus::TRUE) {
            isCellIntersected = true;
        }

        // Add cell neighbours to the process list
        //
        // Neighbours of cells that intersect the surface needs to be tracked,
        // because they should be added to the narrow band regardless of their
        // levelset value.
        if (meshMemoryMode == VolCartesian::MEMORY_LIGHT) {
            for (int face = 0; face < meshCellFaceCount; ++face) {
                long neighId = mesh.getCellFaceNeighsLinearId(cellId, face);
                if (neighId < 0) {
                    continue;
                } else if (isCellInNarrowBand(neighId)) {
                    continue;
                }

                if (outsideSearchRadius.count(neighId) == 0) {
                    processList.insert(neighId);
                }

                if (isCellIntersected) {
                    surfaceNeighbourhood.insert(neighId);
                }
            }
        } else {
            const Cell &cell = mesh.getCell(cellId);
            const long *neighbours = cell.getAdjacencies() ;
            int nNeighbours = cell.getAdjacencyCount() ;
            for (int n = 0; n < nNeighbours; ++n) {
                long neighId = neighbours[n];
                if (isCellInNarrowBand(neighId)) {
                    continue;
                }

                if (outsideSearchRadius.count(neighId) == 0) {
                    processList.insert(neighId);
                }

                if (isCellIntersected) {
                    surfaceNeighbourhood.insert(neighId);
                }
            }
        }
    }

    // Process the neighbours of the cells that intersect the surface
    //
    // If a cell intersects the surface, all its neighbours shoul be added to
    // the narrow band.
    //
    // The search radius for evaluating the support should be large enugh to
    // make sure the closest segment is always found. At least one neighbour
    // of the cell intersects the surface, this means that the distance from
    // the surface should be less than the distance from the cell centroid (C)
    // to the furthest point on the bounding sphere of the the neghbours (B).
    //
    //            +-------+-------+
    //            |       |       |       cell  : is the processed cell
    //            | C +   |   +----->+B   neigh : is the negihbout that
    //            |       |       |               intersects the surface
    //            +-------+-------+
    //             [cell]   [neigh]
    //
    // If the largest side of the cell is equal to h and the radius of the
    // bounding is equal to R, the maximum distance form the surface is:
    //
    //                        CB = h + R
    //
    // regardless of which neighbour intersects the surface.
    //
    // To avoid apporximation errors we arbitrary increase the search radius
    // by 5%.
    double cellBoundingRadius = cartesianKernel->getCellBoundingRadius();

    double neighbourhoodSearchRadius = 0.;
    for (int d = 0; d < mesh.getDimension(); ++d) {
        neighbourhoodSearchRadius = std::max(mesh.getSpacing(d), neighbourhoodSearchRadius);
    }
    neighbourhoodSearchRadius = 1.05 * (neighbourhoodSearchRadius + cellBoundingRadius);

    for (long cellId : surfaceNeighbourhood) {
        if (isCellInNarrowBand(cellId)) {
            continue;
        }

        // Fill cell caches
        narrowBandCache->insertEntry(cellId, true);
        for (LevelSetField field : fieldProcessList) {
            fillFieldCellCache(cellId, field, neighbourhoodSearchRadius);
        }
    }
}

/*!
 * Fill the specified field cache of the given cell.
 *
 * \param id is the id of the cell whose cache will be filled
 * \param field is the field whose cache will be filled
 * \param searchRadius all the portions of the surface with a distance greater than the search
 * radius will not be considered when evaluating the levelset. Trying to fill the cache cache
 * of a cell whose distance is greater than the search radius results in undefined behaviour.
 * Reducing the search radius could speedup the evaluation of levelset information.
 */
void LevelSetSegmentationObject::fillFieldCellCache(long id, LevelSetField field, double distanceHint)
{
    switch (field) {

    case LevelSetField::SUPPORT:
        evalCellSupport(id, distanceHint);
        break;

    default:
        LevelSetObject::fillFieldCellCache(id, field, distanceHint);

    }
}

/*!
 * Evaluate the part associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetSegmentationObject::evalCellPart(long id) const
{
    long support = evalCellSupport(id);
    const SurfUnstructured &surface = getSurface();
    int part = surface.getCell(support).getPID();

    return part;
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified cell.
 *
 * If unsigned levelset is used, the orientation of the suraface normal is discarded and in order
 * to be agnostic with repect the two sides of the surface
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetSegmentationObject::evalCellNormal(long id, bool signedLevelSet) const
{
    std::array<double, 3> normal = evalCellGradient(id, false);
    if (!signedLevelSet && evalCellSign(id) > 0) {
        normal *= -1.;
    }

    return normal;
}

/*!
 * Evaluate the size of the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The size of the segment closest to the specified cell.
 */
double LevelSetSegmentationObject::evalCellSurfaceFeatureSize(long id) const
{
    long support = evalCellSupport(id);
    double segmentSize = getSegmentSize(support);

    return segmentSize;
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationObject::evalCellSupport(long id) const
{
    // Try fetching the value from the cache
    CellValueCache<long> *cache = getFieldCellCache<long>(LevelSetField::SUPPORT);
    if (cache) {
        typename CellValueCache<long>::Entry cacheEntry = cache->findEntry(id);
        if (cacheEntry.isValid()) {
            return *cacheEntry;
        }
    }

    // Evaluate support
    long support = _evalCellSupport(id);

    // Store support in the cache
    if (cache) {
        cache->insertEntry(id, support);
    }

    return support;
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * The result is cached only if there is a segment within the search range.
 *
 * \param id is the id of the cell
 * \param searchRadius all segment whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationObject::evalCellSupport(long id, double searchRadius) const
{
    // Evaluate support
    //
    // Since we are using a limited search radius, it is not possible to use the support
    // stored in the cache, because it may be outside the search radius.
    long support = _evalCellSupport(id, searchRadius);

    // Store support in the cache
    //
    // Support should be stored only if it is valid.
    if (support >= 0) {
        CellValueCache<long> *cache = getFieldCellCache<long>(LevelSetField::SUPPORT);
        if (cache) {
            cache->insertEntry(id, support);
        }
    }

    return support;
}

/*!
 * Evaluate the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetSegmentationObject::evalPart(const std::array<double,3> &point) const
{
    long support = evalSupport(point);
    const SurfUnstructured &surface = getSurface();
    int part = surface.getCell(support).getPID();

    return part;
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetSegmentationObject::evalNormal(const std::array<double,3> &point, bool signedLevelSet) const
{
    std::array<double, 3> normal = evalGradient(point, false);
    if (!signedLevelSet) {
        if (evalSign(point) < 0) {
            normal *= -1.;
        }
    }

    return normal;
}

/*!
 * Evaluate the size of the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The size of the segment closest to the specified point.
 */
double LevelSetSegmentationObject::evalSurfaceFeatureSize(const std::array<double,3> &point) const
{
    long support = evalSupport(point);
    double segmentSize = getSegmentSize(support);

    return segmentSize;
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationObject::evalSupport(const std::array<double,3> &point) const
{
    return _evalSupport(point, std::numeric_limits<double>::max());
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \param searchRadius all segment whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationObject::evalSupport(const std::array<double,3> &point, double searchRadius) const
{
    return _evalSupport(point, searchRadius);
}

/*!
 * Internal function to check if the specified cell lies within the narrow band.
 *
 * A cell is considered within the narrow band if one of the following conditions hold:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - it has at least one negihbour that intersects the zero-levelset iso-surface and its
 *    sign differs form the sign of the neighbour that intersects the surface.
 *
 * Neighbour check is not reliable if the cell is on the last layer of ghosts.
 *
 * \param[in] id is the cell id
 * \param[in] checkNeighbours is set to true, neighbours are check to detect if the cell
 * should be added to the levelset because it has at least one negihbour that intersects
 * the zero-levelset iso-surface and its sign differs form the sign of the neighbour that
 * intersects the surface
 * \param[out] maximumDistance if a valid pointer is provided and the cell is inside the
 * narrow band, on output will contain a conservative estimate for the distance of the
 * cell from the surface, the distance of the cell from the surface will always be less
 * or equal than the provided estimate
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
bool LevelSetSegmentationObject::_isCellInNarrowBand(long id, bool checkNeighbours, double *maximumDistance) const
{
    // Check if the cell is within the narrow band size or if intersects the surface
    //
    // The check is performed evaluating the support of the cell. The search radius for the check
    // is evaluated as the maximum value between the narrow band size and bounding radius of the
    // cell (which defines the distance above which the cell will surely not intersect the
    // surface). In this way, cells that intersect the surface are always included in the narrow
    // band, even if their distance from the surface is greater than then narrow band size
    // explicitly set by the user.
    double cellBoundingRadius = m_kernel->computeCellBoundingRadius(id);
    double searchRadius       = std::max(m_narrowBandSize, cellBoundingRadius);

    long support = _evalCellSupport(id, searchRadius);
    if (support >= 0) {
        if (maximumDistance) {
            *maximumDistance = searchRadius;
        }

        return true;
    }

    // Process cells with neighbours that intersect the zero-levelset iso-surface
    //
    // Cells with at least a negihbour that intersect the zero-levelset iso-surface need to be
    // added to the narrow band if their sign differs form the sign of the neighbour that
    // intersects the surface.
    if (checkNeighbours) {
        const VolumeKernel *mesh = m_kernel->getMesh();
        bool cellAdjacenciesAvailable = (mesh->getAdjacenciesBuildStrategy() != VolumeKernel::ADJACENCIES_NONE);

        const long *cellNeighs;
        int nCellNeighs;
        if (cellAdjacenciesAvailable) {
            const Cell &cell = mesh->getCell(id);
            cellNeighs = cell.getAdjacencies();
            nCellNeighs = cell.getAdjacencyCount();
        } else {
            std::vector<long> cellNeighStorage;
            mesh->findCellFaceNeighs(id, &cellNeighStorage);
            cellNeighs = cellNeighStorage.data();
            nCellNeighs = cellNeighStorage.size();
        }

        for(int n = 0; n < nCellNeighs; ++n){
            long neighId = cellNeighs[n];

            // Skip neighbours that doesn't intersect the surface
            //
            // The check is performed evaluating the support of the neighbour. The search radius
            // for the check is evaluated as bounding radius of the neighbour (which defines the
            // distance above which the neighbour will surely not intersect the surface).
            double neighSearchRadius = m_kernel->computeCellBoundingRadius(neighId);

            long neighSupport = _evalCellSupport(neighId, neighSearchRadius);
            if (neighSupport < 0) {
                continue;
            }

            // Skip neighbours with the same sign
            long neighSign = evalCellSign(neighId);
            if (neighId == neighSign) {
                continue;
            }

            // Cell is inside the narrow band
            //
            // The cell has a neighbour with opposite sign the intersects the zero-levelset
            // iso-surface.
            if (maximumDistance) {
                std::array<double,3> cellCentroid = m_kernel->computeCellCentroid(id);
                std::array<double,3> neighProjectionPoint = evalCellProjectionPoint(neighId);
                *maximumDistance = norm2(cellCentroid - neighProjectionPoint);
            }

            return true;
        }
    }

    // The cell is not in the narrow band
    if (maximumDistance) {
        *maximumDistance = -1.;
    }

    return false;
}

/*!
 * Evaluate levelset sign at the specified cell.
 *
 * \param id is the id of the cell
 * \result The sign of the levelset at the specified cell.
 */
short LevelSetSegmentationObject::_evalCellSign(long id) const
{
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalSign(centroid);
}

/*!
 * Evaluate levelset value at the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified cell.
 */
double LevelSetSegmentationObject::_evalCellValue(long id, bool signedLevelSet) const
{
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalValue(centroid, signedLevelSet);
}

/*!
 * Evaluate levelset gradient at the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The gradient of the levelset at the specified cell.
 */
std::array<double,3> LevelSetSegmentationObject::_evalCellGradient(long id, bool signedLevelSet) const
{
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalGradient(centroid, signedLevelSet);
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationObject::_evalCellSupport(long id) const
{
    return _evalCellSupport(id, std::numeric_limits<double>::max());
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param searchRadius all segment whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationObject::_evalCellSupport(long id, double searchRadius) const
{
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalSupport(centroid, searchRadius);
}

/*!
 * Evaluate levelset sign at the specified point.
 *
 * \param point are the coordinates of the point
 * \result The sign of the levelset at the specified point.
 */
short LevelSetSegmentationObject::_evalSign(const std::array<double,3> &point) const
{
    return static_cast<short>(sign(_evalValue(point, true)));
}

/*!
 * Evaluate levelset value at the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified point.
 */
double LevelSetSegmentationObject::_evalValue(const std::array<double,3> &point, bool signedLevelSet) const
{
    long support = evalSupport(point);

    double distance;
    int error = evalLevelsetInfo(point, support, signedLevelSet, &distance, nullptr);
    if (error) {
        throw std::runtime_error ("Unable to extract the levelset information from segment.");
    }

    return distance;
}

/*!
 * Evaluate levelset gradient at the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The gradient of the levelset at the specified point.
 */
std::array<double,3> LevelSetSegmentationObject::_evalGradient(const std::array<double,3> &point, bool signedLevelSet) const
{
    long support = evalSupport(point);

    std::array<double,3> gradient;
    int error = evalLevelsetInfo(point, support, signedLevelSet, nullptr, &gradient);
    if (error) {
        throw std::runtime_error ("Unable to extract the levelset information from segment.");
    }

    return gradient;
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationObject::_evalSupport(const std::array<double,3> &point) const
{
    return _evalSupport(point, std::numeric_limits<double>::max());
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \param searchRadius all segment whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationObject::_evalSupport(const std::array<double,3> &point, double searchRadius) const
{
    long closestSegmentId;
    double closestDistance;
    getSearchTree().findPointClosestCell(point, searchRadius, &closestSegmentId, &closestDistance);

    return closestSegmentId;
}

/*!
 * Write the specified field to the given stream.
 *
 * @param[in] field is the field that will be written
 * @param[in] stream output stream
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 */
void LevelSetSegmentationObject::flushField(LevelSetField field, std::fstream &stream, VTKFormat format) const {

    switch(field) {

    case LevelSetField::SUPPORT:
    {
        void (*writeFunctionPtr)(std::fstream &, const long &) = nullptr;

        if (format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<long>;
        } else if (format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<long>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        CellValueCache<long> *cache = getFieldCellCache<long>(LevelSetField::SUPPORT);
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            if (cache) {
                long cellId = cell.getId();
                CellValueCache<long>::Entry cacheEntry = cache->findEntry(cellId);
                if (cacheEntry.isValid()) {
                    (*writeFunctionPtr)(stream, *cacheEntry);
                } else {
                    (*writeFunctionPtr)(stream, levelSetDefaults::SUPPORT);
                }
            } else {
                (*writeFunctionPtr)(stream, levelSetDefaults::SUPPORT);
            }
        }

        break;
    }

    case LevelSetField::PART:
    {
        void (*writeFunctionPtr)(std::fstream &, const int &) = nullptr;

        if (format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<int>;
        } else if (format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<int>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        CellValueCache<long> *cache = getFieldCellCache<long>(LevelSetField::SUPPORT);
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            if (cache) {
                long cellId = cell.getId();
                CellValueCache<long>::Entry cacheEntry = cache->findEntry(cellId);
                if (cacheEntry.isValid()) {
                    (*writeFunctionPtr)(stream, evalCellPart(cellId));
                } else {
                    (*writeFunctionPtr)(stream, levelSetDefaults::PART);
                }
            } else {
                (*writeFunctionPtr)(stream, levelSetDefaults::PART);
            }
        }

        break;
    }

    case LevelSetField::NORMAL:
    {
        void (*writeFunctionPtr)(std::fstream &, const std::array<double,3> &) = nullptr;

        if (format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<std::array<double,3>>;
        } else if (format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<std::array<double,3>>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        CellValueCache<long> *cache = getFieldCellCache<long>(LevelSetField::SUPPORT);
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            if (cache) {
                long cellId = cell.getId();
                CellValueCache<long>::Entry cacheEntry = cache->findEntry(cellId);
                if (cacheEntry.isValid()) {
                    (*writeFunctionPtr)(stream, evalCellNormal(cellId, true));
                } else {
                    (*writeFunctionPtr)(stream, levelSetDefaults::GRADIENT);
                }
            } else {
                (*writeFunctionPtr)(stream, levelSetDefaults::GRADIENT);
            }
        }

        break;
    }

    default:
    {
        LevelSetObject::flushField(field, stream, format);

        break;
    }

    }
}

/*!
 * Get the part associated with the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetSegmentationObject::getPart(long cellId) const
{
    return evalCellPart(cellId);
}

/*!
 * Get the normal of the surface at the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetSegmentationObject::getNormal(long cellId) const
{
    return evalCellNormal(cellId, m_defaultSignedLevelSet);
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationObject::getSupport(long cellId) const
{
    return evalCellSupport(cellId);
}

/*!
 * Get the size of the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The size of the segment closest to the specified cell.
 */
double LevelSetSegmentationObject::getSurfaceFeatureSize(long cellId) const
{
    return evalCellSurfaceFeatureSize(cellId);
}

/*!
 * Get the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetSegmentationObject::getPart(const std::array<double,3> &point) const
{
    return evalPart(point);
}

/*!
 * Get the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetSegmentationObject::getNormal(const std::array<double,3> &point) const
{
    return evalNormal(point, m_defaultSignedLevelSet);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationObject::getSupport(const std::array<double,3> &point) const
{
    return evalSupport(point);
}

/*!
 * Get the size of the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The size of the segment closest to the specified point.
 */
double LevelSetSegmentationObject::getSurfaceFeatureSize(const std::array<double,3> &point) const
{
    return evalSurfaceFeatureSize(point);
}

}
