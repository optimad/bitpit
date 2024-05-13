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

#include "levelSetSegmentationObject.hpp"

#include "bitpit_CG.hpp"
#include "levelSetObject.hpp"

namespace bitpit {

/*!
    @class      LevelSetSegmentationSurfaceInfo
    @ingroup    levelset
    @brief      Segmentation kernel
*/

/*!
 * The default angle that is used to identify sharp edges. If the angle between two segments
 * is bigger than this angle, the enclosed edge is considered as a sharp edge.
 *
 */
const double LevelSetSegmentationSurfaceInfo::DEFAULT_FEATURE_ANGLE = 2. * BITPIT_PI;

/*!
 * Default constructor
 */
LevelSetSegmentationSurfaceInfo::LevelSetSegmentationSurfaceInfo()
    : m_surface(nullptr),
      m_featureAngle(0),
      m_surfaceSmoothing(LevelSetSurfaceSmoothing::LOW_ORDER)
{
}

/*!
 * Copy constructor
 */
LevelSetSegmentationSurfaceInfo::LevelSetSegmentationSurfaceInfo(const LevelSetSegmentationSurfaceInfo &other)
    : m_surface(other.m_surface),
      m_featureAngle(other.m_featureAngle),
      m_surfaceSmoothing(other.m_surfaceSmoothing),
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
 * Constructor
 *
 * @param[in] surfaceSmoothing is the given surface snoothing order
 */
LevelSetSegmentationSurfaceInfo::LevelSetSegmentationSurfaceInfo(LevelSetSurfaceSmoothing surfaceSmoothing)
    : m_surface(nullptr),
      m_featureAngle(0),
      m_surfaceSmoothing(surfaceSmoothing)
{
}

/*!
 * Constructor
 *
 * @param[in,out] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 * @param[in] surfaceSmoothing is the given surface snoothing order
 */
LevelSetSegmentationSurfaceInfo::LevelSetSegmentationSurfaceInfo(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing)
    : m_surfaceSmoothing(surfaceSmoothing)
{
    setSurface(std::move(surface), featureAngle);
}

/*!
 * Constructor
 *
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 * @param[in] surfaceSmoothing is the given surface snoothing order
 */
LevelSetSegmentationSurfaceInfo::LevelSetSegmentationSurfaceInfo(const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing)
{
    setSurface(surface, featureAngle, surfaceSmoothing);
}

/*!
 * Get segmentation surface
 * @return segmentation surface;
 */
const SurfUnstructured & LevelSetSegmentationSurfaceInfo::getSurface() const {
    return *m_surface;
}

/*!
 * Set the surface
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
void LevelSetSegmentationSurfaceInfo::setSurface(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle){

    m_ownedSurface = std::move(surface);
    setSurface(m_ownedSurface.get(), featureAngle);
}

/*!
 * Set the surface
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 * @param[in] surfaceSmoothing is the given surface snoothing order
 */
void LevelSetSegmentationSurfaceInfo::setSurface(const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing){

    // Check if adjacencies are built
    if (surface->getAdjacenciesBuildStrategy() == SurfUnstructured::ADJACENCIES_NONE) {
        throw std::runtime_error ("Segmentation needs adjacencies!") ;
    }

    // Surface information
    m_surface          = surface;
    m_featureAngle     = featureAngle;
    m_surfaceSmoothing = surfaceSmoothing;

    // Segment vertices information
    m_segmentVertexOffset.unsetKernel();
    m_segmentVertexOffset.setStaticKernel(&m_surface->getCells());

    std::size_t nTotalSegmentVertices = 0;
    for (auto itr = m_segmentVertexOffset.begin(); itr != m_segmentVertexOffset.end(); ++itr) {
        *itr = nTotalSegmentVertices;
        nTotalSegmentVertices += m_surface->getCells().rawAt(itr.getRawIndex()).getVertexCount();
    }

    // Normals
    m_segmentNormalsValid.unsetKernel();
    m_segmentNormalsValid.setStaticKernel(&m_surface->getCells());
    m_segmentNormalsValid.fill(false);
    m_segmentNormalsStorage.unsetKernel();
    m_segmentNormalsStorage.setStaticKernel(&m_surface->getCells());

    m_unlimitedVertexNormalsValid.unsetKernel();
    m_unlimitedVertexNormalsValid.setStaticKernel(&m_surface->getVertices());
    m_unlimitedVertexNormalsValid.fill(false);
    m_unlimitedVertexNormalsStorage.unsetKernel();
    m_unlimitedVertexNormalsStorage.setStaticKernel(&m_surface->getVertices());

    m_limitedSegmentVertexNormalValid.assign(nTotalSegmentVertices, false);

    // Initialize search tree
    m_searchTree = std::unique_ptr<SurfaceSkdTree>(new SurfaceSkdTree(surface));
    m_searchTree->build();
}

/*!
 * Set the surface
 * @param[in] surface pointer to surface
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
void LevelSetSegmentationSurfaceInfo::setSurface(const SurfUnstructured *surface, double featureAngle){
    setSurface(surface, featureAngle, LevelSetSurfaceSmoothing::LOW_ORDER);
}

/*!
 * Get search tree
 * @return search tree;
 */
const SurfaceSkdTree & LevelSetSegmentationSurfaceInfo::getSearchTree() const {
    return *m_searchTree;
}

/*!
 * Get feature angle
 * @return feature angle used when calculating face normals.
 */
double LevelSetSegmentationSurfaceInfo::getFeatureAngle() const {
    return m_featureAngle;
}

/*!
 * Get smoothing order (low or high) imposed on surface when calculating the
 * projection point and normal on the surface
 * @return the ssurface smoothing order
 */
bitpit::LevelSetSurfaceSmoothing LevelSetSegmentationSurfaceInfo::getSurfaceSmoothing() const {
    return m_surfaceSmoothing;
}

/*!
 * Evaluate the distance function at the specified point.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[in] signedDistance controls if the signed or unsigned distance will be evaluated
 * @return The distance function at the specified point.
 */
double LevelSetSegmentationSurfaceInfo::evalDistance(const std::array<double, 3> &point,
                                                     const SegmentConstIterator &segmentItr,
                                                     bool signedDistance) const
{
    // Project the point on the surface and evaluate the point-projection vector
    int nSegmentVertices = segmentItr->getVertexCount();
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    std::array<double, 3> pointProjection = evalProjection(point, segmentItr, lambda);
    std::array<double, 3> pointProjectionVector = point - pointProjection;

    // Evaluate unsigned distance
    double unsignedDistance = norm2(pointProjectionVector);
    if (!signedDistance) {
        return unsignedDistance;
    }

    // Signed distance
    //
    // If the sign is null and the point doesn't lie on the segmentation, it lies on the normal
    // plane. This case is not supported, because it would require to evaluate the sign taking
    // into account the the curvature of the surface.
    std::array<double, 3> pseudoNormal = computePseudoNormal(segmentItr, lambda);
    double pointProjectionNormalComponent = dotProduct(pointProjectionVector, pseudoNormal);

    double distanceTolerance = m_surface->getTol();
    if (utils::DoubleFloatingEqual()(pointProjectionNormalComponent, 0., distanceTolerance, distanceTolerance)) {
        bool pointOnSegmentation = utils::DoubleFloatingEqual()(unsignedDistance, 0., distanceTolerance, distanceTolerance);
        if (!pointOnSegmentation) {
            throw std::runtime_error("Unable to evaluate point sign: the point lies on the normal plane!");
        }
    }

    return sign(pointProjectionNormalComponent) * unsignedDistance;
}

/*!
  * Evaluate the distance vector function at the specified point.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @return The distance vector function at the specified point.
 */
std::array<double, 3> LevelSetSegmentationSurfaceInfo::evalDistanceVector(const std::array<double, 3> &point,
                                                                          const SegmentConstIterator &segmentItr) const
{
    int nSegmentVertices = segmentItr->getVertexCount();
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    std::array<double, 3> pointProjection = evalProjection(point, segmentItr, lambda);

    return (point - pointProjection);
}

/*!
 * Evaluate the surface normal at the projection of the specified point.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @return The surface normal at the projection of the specified point.
 */
std::array<double, 3> LevelSetSegmentationSurfaceInfo::evalNormal(const std::array<double, 3> &point,
                                                                  const SegmentConstIterator &segmentItr) const
{
    // Project the point on the surface and evaluate the point-projection vector
    int nSegmentVertices = segmentItr->getVertexCount();
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    evalProjection(point, segmentItr, lambda);

    // Evaluate normal
    return computeSurfaceNormal(segmentItr, lambda);
}

/*!
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is vertical to the normal vectors associated with them.
 *
 * In case the segment is just a vertex, the projection coisides with it.
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] segmentItr is an iterator pointing to the segment on which the surface
 * will be created.
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 * projection point.
 */
void LevelSetSegmentationSurfaceInfo::evalProjectionOnVertex(const std::array<double,3> &point,
                                                             const SegmentConstIterator &segmentItr,
                                                             std::array<double, 3> *projectionPoint,
                                                             std::array<double, 3> *projectionNormal) const
{
    BITPIT_UNUSED(point);

    // Get segment
    const Cell &segment = *segmentItr;
    assert(segment.getType() == ElementType::VERTEX);

    // Get vertex id
    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds();
    long id = segmentVertexIds[0];

    // Compute projection point and normal 
    (*projectionPoint)  = m_surface->getVertexCoords(id);
    (*projectionNormal) = {0., 0., 0.};
}

/*!
 * Evaluate the projection of the given point on the specified linear segment.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 */
void LevelSetSegmentationSurfaceInfo::evalLowOrderProjectionOnLine(const std::array<double, 3> &point,
                                                                   const SegmentConstIterator &segmentItr,
                                                                   std::array<double, 3> *projectionPoint,
                                                                   std::array<double, 3> *projectionNormal) const
{
    const Cell &segment = *segmentItr;
    assert(segment.getType() == ElementType::LINE);

    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds();

    const std::array<double,3> &point0 = m_surface->getVertexCoords(segmentVertexIds[0]);
    const std::array<double,3> &point1 = m_surface->getVertexCoords(segmentVertexIds[1]);

    int nSegmentVertices = segment.getVertexCount();
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    (*projectionPoint) = CGElem::projectPointSegment(point, point0, point1, lambda);

    std::array<double, 3> normal0 = computeSegmentVertexNormal(segmentItr, 0, true);
    std::array<double, 3> normal1 = computeSegmentVertexNormal(segmentItr, 1, true);

    (*projectionNormal) = lambda[0] * normal0 + lambda[1] * normal1;
    (*projectionNormal) /= norm2((*projectionNormal));
}

/*!
 * Evaluate the projection of the given point on the specified triangular segment.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 */
void LevelSetSegmentationSurfaceInfo::evalLowOrderProjectionOnTriangle(const std::array<double, 3> &point,
                                                                       const SegmentConstIterator &segmentItr,
                                                                       std::array<double, 3> *projectionPoint,
                                                                       std::array<double, 3> *projectionNormal) const
{
    const Cell &segment = *segmentItr;
    assert(segment.getType() == ElementType::TRIANGLE);

    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds();

    const std::array<double,3> &point0 = m_surface->getVertexCoords(segmentVertexIds[0]);
    const std::array<double,3> &point1 = m_surface->getVertexCoords(segmentVertexIds[1]);
    const std::array<double,3> &point2 = m_surface->getVertexCoords(segmentVertexIds[2]);

    int nSegmentVertices = segment.getVertexCount();
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    (*projectionPoint) = CGElem::projectPointTriangle(point, point0, point1, point2, lambda);

    std::array<double, 3> normal0 = computeSegmentVertexNormal(segmentItr, 0, true);
    std::array<double, 3> normal1 = computeSegmentVertexNormal(segmentItr, 1, true);
    std::array<double, 3> normal2 = computeSegmentVertexNormal(segmentItr, 2, true);

    (*projectionNormal) = lambda[0] * normal0 + lambda[1] * normal1 + lambda[2] * normal2;
    (*projectionNormal) /= norm2((*projectionNormal));
}

/*!
 * Evaluate the projection of the given point on the specified polygonalsegment.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 */
void LevelSetSegmentationSurfaceInfo::evalLowOrderProjectionOnPolygon(const std::array<double, 3> &point,
                                                                      const SegmentConstIterator &segmentItr,
                                                                      std::array<double, 3> *projectionPoint,
                                                                      std::array<double, 3> *projectionNormal) const
{
    const Cell &segment = *segmentItr;

    ConstProxyVector<long> segmentVertexIds = m_surface->getFacetOrderedVertexIds(segment);

    std::size_t nSegmentVertices = segmentVertexIds.size();
    BITPIT_CREATE_WORKSPACE(segmentVertexCoords, std::array<double BITPIT_COMMA 3>, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    m_surface->getVertexCoords(segmentVertexIds.size(), segmentVertexIds.data(), segmentVertexCoords);

    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    (*projectionPoint) = CGElem::projectPointPolygon(point, nSegmentVertices, segmentVertexCoords, lambda);

    (*projectionNormal) = lambda[0] * computeSegmentVertexNormal(segmentItr, 0, true);
    for (std::size_t i = 1; i < nSegmentVertices; ++i) {
        (*projectionNormal) += lambda[i] * computeSegmentVertexNormal(segmentItr, i, true);
    }
    (*projectionNormal) /= norm2(*projectionNormal);
}

/*!
 * Evaluate the projection of the given point on the specified segment.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 */
void LevelSetSegmentationSurfaceInfo::evalLowOrderProjection(const std::array<double, 3> &point,
                                                             const SegmentConstIterator &segmentItr,
                                                             std::array<double, 3> *projectionPoint,
                                                             std::array<double, 3> *projectionNormal) const
{
    const Cell &segment = *segmentItr;
    ElementType segmentType = segment.getType();
    switch (segmentType) {

    case ElementType::VERTEX:
    {
        evalProjectionOnVertex(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    case ElementType::LINE:
    {
        evalLowOrderProjectionOnLine(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    case ElementType::TRIANGLE:
    {
        evalLowOrderProjectionOnTriangle(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    default:
    {
        evalLowOrderProjectionOnPolygon(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    }
}

/*!
 * Evaluate the projection of the given point on the specified segment.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[out] lambda on output will contain the barycentric coordinates of the projection point
 * @return The coordinates of the projection point.
 */
std::array<double, 3> LevelSetSegmentationSurfaceInfo::evalProjection(const std::array<double, 3> &point,
                                                                      const SegmentConstIterator &segmentItr,
                                                                      double *lambda) const
{
    std::array<double, 3> projectionPoint;
    std::array<double, 3> projectionNormal;
    evalLowOrderProjection(point, segmentItr, &projectionPoint, &projectionNormal);
    BITPIT_UNUSED(projectionNormal);

    const Cell &segment = *segmentItr;
    m_surface->evalBarycentricCoordinates(segment.getId(), point, lambda);
    return projectionPoint;
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
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[in] lambda are the barycentric coordinates of the point
 * @return the pseudo-normal at specified point of the given triangle
 */
std::array<double,3> LevelSetSegmentationSurfaceInfo::computePseudoNormal(const SegmentConstIterator &segmentItr,
                                                                          const double *lambda ) const
{
    // Early return if the segment is a point
    const Cell &segment = *segmentItr;
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
        if (positionFlag > 0) {
            int polygonVertex = positionFlag - 1;
            int elementVertex = m_surface->getFacetOrderedLocalVertex(segment, polygonVertex);

            positionFlag = elementVertex + 1;
        } else if (positionFlag < 0) {
            int polygonEdge = (- positionFlag) - 1;
            int elementEdge = m_surface->getFacetOrderedLocalEdge(segment, polygonEdge);

            positionFlag = - (elementEdge + 1);
        }
    }

    std::array<double,3> pseudoNormal;
    if (positionFlag == 0) {
        pseudoNormal = computeSegmentNormal(segmentItr);
    } else if (positionFlag > 0) {
        int vertex = positionFlag - 1;
        pseudoNormal = computeSegmentVertexNormal(segmentItr, vertex, false);
    } else {
        int edge = (- positionFlag) - 1;
        pseudoNormal = computeSegmentEdgeNormal(segmentItr, edge, false);
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
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[in] lambda are the barycentric coordinates of the point
 * @return the surface-normal at specified point of the given triangle
 */
std::array<double,3> LevelSetSegmentationSurfaceInfo::computeSurfaceNormal(const SegmentConstIterator &segmentItr,
                                                                           const double *lambda ) const
{
    // Early return if the segment is a point
    const Cell &segment = *segmentItr;
    ElementType segmentType = segment.getType();
    if (segmentType == ElementType::VERTEX) {
        return {{0., 0., 0.}};
    }

    // Evaluate surface normal
    std::size_t nSegmentVertices = segment.getVertexCount();
    std::array<double,3> surfaceNormal = lambda[0] * computeSegmentVertexNormal(segmentItr, 0, true);
    for (std::size_t i = 1; i < nSegmentVertices; ++i) {
        surfaceNormal += lambda[i] * computeSegmentVertexNormal(segmentItr, i, true);
    }
    surfaceNormal /= norm2(surfaceNormal);

    return surfaceNormal;
}

/*!
 * Compute the normal of the specified triangle.
 *
 * To reduce computational times, normals of vertices are cached.
 *
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @return the normal of the specified triangle
 */
std::array<double,3> LevelSetSegmentationSurfaceInfo::computeSegmentNormal(const SegmentConstIterator &segmentItr ) const {

    std::size_t segmentRawId = segmentItr.getRawIndex();
    std::array<double, 3> *segmentNormal = m_segmentNormalsStorage.rawData(segmentRawId);
    if (!m_segmentNormalsValid.rawAt(segmentRawId)) {
        *segmentNormal = m_surface->evalFacetNormal(segmentItr->getId());
        m_segmentNormalsValid.rawAt(segmentRawId) = true;
    }

    return *segmentNormal;
}

/*!
 * Compute the normal of the specified triangle's edge.
 *
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[in] edge is the local index of the edge
 * @param[in] limited is a flag controling if the limited or the unlimited normal will
 * be evaluated
 * @return the normal of the specified triangle's edge
 */
std::array<double,3> LevelSetSegmentationSurfaceInfo::computeSegmentEdgeNormal(const SegmentConstIterator &segmentItr, int edge, bool limited ) const {

    long segmentId = segmentItr.getId();

    std::array<double, 3> limitedEdgeNormal;
    std::array<double, 3> unlimitedEdgeNormal;
    m_surface->evalEdgeNormals(segmentId, edge, m_featureAngle, &unlimitedEdgeNormal, &limitedEdgeNormal) ;

    if (limited) {
        return limitedEdgeNormal;
    }

    return unlimitedEdgeNormal;
}

/*!
 * Compute the normal of the specified triangle's vertex.
 *
 * To reduce computational times, normals of vertices are cached.
 *
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[in] vertex is the local index of the vertex
 * @param[in] limited controls is the limited or the unlimited normal will
 * be evaluated
 * @return the normal of the specified triangle's vertex
 */
std::array<double,3> LevelSetSegmentationSurfaceInfo::computeSegmentVertexNormal(const SegmentConstIterator &segmentItr, int vertex, bool limited ) const {

    // Segment information
    long segmentId = segmentItr.getId();
    long segmentRawId = segmentItr.getRawIndex();
    const Cell &segment = *segmentItr;

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
	@class      LevelSetSegmentationBaseObject
	@ingroup    levelset
	@brief      Implements visitor pattern fo segmentated geometries
*/

/*!
 * If the search radius is set equal to the constant AUTOMATIC_SEARCH_RADIUS, the object will try
 * to evaluate the optional search radius for the specified cell. The automatic evaluation of the
 * search radius is possible only for a limited number of cases, when the automatic evaluation
 * cannot be performed, an infinite search radius will be used.
 */
const double LevelSetSegmentationBaseObject::AUTOMATIC_SEARCH_RADIUS = -1;

/*!
 * Create the cache that will be used for storing cell information of the specified field.
 *
 * \param field is the field for which the caches will be registered
 * \param cacheId is the id that will be associated with the cache, if a NULL_ID is specified
 * the cache id will be assigned automatically
 * \result The id associated with the registered cache.
 */
std::size_t LevelSetSegmentationBaseObject::createFieldCellCache(LevelSetField field, std::size_t cacheId)
{
    switch(field) {

    case LevelSetField::SUPPORT:
        return createFieldCellCache<long>(field, cacheId);

    default:
        return LevelSetObject::createFieldCellCache(field, cacheId);

    }
}

/*!
 * Get the list of supported field.
 * @result The list of supported field.
 */
LevelSetFieldset LevelSetSegmentationBaseObject::getSupportedFields() const
{
    LevelSetFieldset supportedFields = LevelSetObject::getSupportedFields();
    supportedFields.push_back(LevelSetField::PART);
    supportedFields.push_back(LevelSetField::NORMAL);
    supportedFields.push_back(LevelSetField::SUPPORT);

    return supportedFields;
}

/*!
 * Fill the specified field cache of the given cell.
 *
 * \param field is the field whose cache will be filled
 * \param id is the id of the cell whose cache will be filled
 */
void LevelSetSegmentationBaseObject::fillFieldCellCache(LevelSetField field, long id)
{
    switch (field) {

    case LevelSetField::SUPPORT:
        evalCellSupport(id);
        break;

    default:
        LevelSetObject::fillFieldCellCache(field, id);

    }
}

/*!
 * Evaluate the surface associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The surface associated with the segment closest to the specified cell.
 */
const SurfUnstructured & LevelSetSegmentationBaseObject::evalCellSurface(long id) const
{
    return _evalCellSurface(id);
}

/*!
 * Evaluate the part associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetSegmentationBaseObject::evalCellPart(long id) const
{
    auto evaluator = [this] (long id)
        {
            return _evalCellPart(id);
        };

    auto fallback = [] (long id)
        {
            BITPIT_UNUSED(id);

            return levelSetDefaults::PART;
        };

    LevelSetField field = LevelSetField::PART;
    int part = evalCellFieldCached<int>(field, id, evaluator, fallback);

    return part;
}

/*!
 * Check if cell intersects the surface.
 *
 * If mode==LevelSetIntersectionMode::FAST_FUZZY the method will compare the levelset
 * value to tangent and bounding radius of a cell. If the value is smaller than the
 * tangent radius LevelSetIntersectionStatus::TRUE is returned, if it is larger than the
 * bounding radius LevelSetIntersectionStatus::FALSE is returned. If it is in-between
 * LevelSetIntersectionStatus::CLOSE is returned.
 *
 * If mode==LevelSetIntersectionMode::FAST_GUARANTEE_TRUE and the levelset value is
 * smaller than the rangent radius LevelSetIntersectionStatus::TRUE is returned,
 * otherwise LevelSetIntersectionStatus::FALSE.
 *
 * If mode==LevelSetIntersectionMode::FAST_GURANTEE_FALSE and the levelset value is
 * larger than the bounding radius LevelSetIntersectionStatus::FALSE is returned,
 * otherwise LevelSetIntersectionStatus::TRUE.
 *
 * If mode==LevelSetIntersectionMode::ACCURATE, the same checks of fuzzy mode are
 * performed, however, in the cases where fuzzy mode would return CLOSE, an additional
 * check on the intersection between the tangent plane at the projection point and the
 * cell is performed. Errors of the method are related to the ratio of surface curvature
 * over cell size.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * cell vertices and has the center in the cell centroid.
 *
 * The tangent sphere is a sphere having the center in the level centroid and tangent
 * to the cell.
 *
 * @param[in] id cell id
 * @param[in] distance is the unsigned distance of the cell centroid from the zero-levelset
 * iso-surface
 * @param[in] mode describes the types of check that should be performed
 * @return indicator regarding intersection
 */
LevelSetIntersectionStatus LevelSetSegmentationBaseObject::_intersectSurface(long id, double distance, LevelSetIntersectionMode mode) const
{
    // Get surface information
    const SurfUnstructured &surface = evalCellSurface(id);

    // Early return if the surface is empty
    if (surface.empty()) {
        return LevelSetIntersectionStatus::FALSE;
    }

    // Evaluate intersection using base class
    return LevelSetObject::_intersectSurface(id, distance, mode);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetSegmentationBaseObject::evalCellNormal(long id, bool signedLevelSet) const
{
    // Evaluate signed normal
    //
    // The normal stored in the cache is unsigned.
    auto evaluator = [this] (long id)
        {
            return _evalCellNormal(id, false);
        };

    auto fallback = [] (long id)
        {
            BITPIT_UNUSED(id);

            return levelSetDefaults::NORMAL;
        };

    LevelSetField field = LevelSetField::NORMAL;
    std::array<double, 3> normal = evalCellFieldCached<std::array<double, 3>>(field, id, evaluator, fallback);

    // Evaluate the normal with the correct signdness
    //
    // If an unsigned evaluation is requested, the orientation of the surface should be discarded
    // and in order to have a normal that is agnostic with respect the two sides of the surface.
    if (signedLevelSet) {
        short cellSign = evalCellSign(id);
        if (cellSign <= 0) {
            normal *= static_cast<double>(cellSign);;
        }
    }

    return normal;
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * The result is cached only if there is a segment within the search range.
 *
 * \param id is the id of the cell
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support. If the search radius is set equal to the
 * constant AUTOMATIC_SEARCH_RADIUS, the object will try to evaluate the optional search radius
 * for the specified cell. The automatic evaluation of the search radius is possible only for
 * a limited number of cases, when the automatic evaluation cannot be performed, an infinite
 * search radius will be used.
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationBaseObject::evalCellSupport(long id, double searchRadius) const
{
    // Evaluate support
    //
    // If we are using a limited search radius, it is not possible to use the support stored
    // in the cache, because it may be outside the search radius.
    auto evaluator = [this, searchRadius] (long id)
        {
            // Evaluate the search radius for support evaluation
            //
            // For the automatic evaluation of the search range, it is possible to use the
            // information about the cell zone:
            //  - cells that are inside the narrow band because their distance from the surface is
            //    less than the narrow band size can use a search radius equal to the narrow band
            //    size
            //  - cells that are inside the narrow band because they intersects the surface can use
            //    a search radius equal to their bounding radius;
            //  - cells that are inside the narrow band because their one of their face neighbours
            //    intersect the surface can use a search radius equal to their bounding radius plus
            //    the bounding diameter of the smallest neighbour that intersects the surface.
            double supportSearchRadius;
            if (searchRadius == AUTOMATIC_SEARCH_RADIUS) {
                LevelSetCellLocation cellLocation = getCellLocation(id);
                if (cellLocation == LevelSetCellLocation::NARROW_BAND_DISTANCE) {
                    supportSearchRadius = this->m_narrowBandSize;
                } else if (cellLocation == LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
                    supportSearchRadius = m_kernel->computeCellBoundingRadius(id);
                } else if (cellLocation == LevelSetCellLocation::NARROW_BAND_NEIGHBOUR) {
                    // Evaluate the bounding diameter of the smallest intersected face neighbour
                    supportSearchRadius = std::numeric_limits<double>::max();
                    auto neighProcessor = [this, &supportSearchRadius](long neighId, int layer) {
                        BITPIT_UNUSED(layer);

                        // Discard neighbours that are not intersected
                        LevelSetCellLocation neighLocation = getCellLocation(neighId);
                        if (neighLocation != LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
                            return false;
                        }

                        // Consider the bounding diameter of the smallest intersected neighbour
                        supportSearchRadius = std::min(2 * m_kernel->computeCellBoundingRadius(neighId), supportSearchRadius);

                        // Continue processing the other neighbours
                        return false;
                    };

                    const VolumeKernel &mesh = *(m_kernel->getMesh()) ;
                    mesh.processCellFaceNeighbours(id, 1, neighProcessor);

                    // Ad the bounding radius of the cell
                    supportSearchRadius += m_kernel->computeCellBoundingRadius(id);
                } else {
                    supportSearchRadius = std::numeric_limits<double>::max();
                }
            } else {
                supportSearchRadius = searchRadius;
            }

            // Evaluate cell support
            return _evalCellSupport(id, supportSearchRadius);
        };

    auto fallback = [] (long id)
        {
            BITPIT_UNUSED(id);

            return levelSetDefaults::SUPPORT;
        };

    LevelSetField field = LevelSetField::SUPPORT;

    long support;
    if (searchRadius < std::numeric_limits<double>::max() && searchRadius != AUTOMATIC_SEARCH_RADIUS) {
        // Evaluate the support from scratch
        support = evalCellField<long>(field, id, evaluator, fallback);

        // Update the cache if the support is valid
        if (support >= 0) {
            fillFieldCellCache(field, id, support);
        }
    } else {
        // Evaluate the support
        support = evalCellFieldCached<long>(field, id, evaluator, fallback);
    }

    return support;
}

/*!
 * Evaluate the surface associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The surface associated with the segment closest to the specified point.
 */
const SurfUnstructured & LevelSetSegmentationBaseObject::evalSurface(const std::array<double,3> &point) const
{
    return _evalSurface(point);
}

/*!
 * Evaluate the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetSegmentationBaseObject::evalPart(const std::array<double,3> &point) const
{
    return _evalPart(point);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetSegmentationBaseObject::evalNormal(const std::array<double,3> &point, bool signedLevelSet) const
{
    return _evalNormal(point, signedLevelSet);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationBaseObject::evalSupport(const std::array<double,3> &point) const
{
    return _evalSupport(point);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationBaseObject::evalSupport(const std::array<double,3> &point, double searchRadius) const
{
    return _evalSupport(point, searchRadius);
}

/*!
 * Evaluate the part associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetSegmentationBaseObject::_evalCellPart(long id) const
{
    long support = evalCellSupport(id);
    const SurfUnstructured &surface = evalCellSurface(id);
    int part = surface.getCell(support).getPID();

    return part;
}

/*!
 * Evaluate the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetSegmentationBaseObject::_evalPart(const std::array<double,3> &point) const
{
    long support = evalSupport(point);
    const SurfUnstructured &surface = evalSurface(point);
    int part = surface.getCell(support).getPID();

    return part;
}

/*!
 * Add the VTK data associated with the specified field.
 *
 * @param[in] field is the field
 * @param[in] objectName is the name that will be associated with the object
 */
void LevelSetSegmentationBaseObject::addVTKOutputData( LevelSetField field, const std::string &objectName)
{
    VTK &vtkWriter = this->m_kernel->getMesh()->getVTK() ;
    std::string name = this->getVTKOutputDataName(field, objectName);

    switch (field) {

        case LevelSetField::SUPPORT:
            vtkWriter.addData<long>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
            break;

        case LevelSetField::PART:
            vtkWriter.addData<int>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
            break;

        case LevelSetField::NORMAL:
            vtkWriter.addData<double>( name, VTKFieldType::VECTOR, VTKLocation::CELL, this);
            break;

        default:
            LevelSetObject::addVTKOutputData(field, objectName);
            break;

    }
}

/*!
 * Get the name that will be used by the VTK writer for the specifed field.
 *
 * @param[in] field is the field
 * @result The name that will be used by the VTK writer for the specifed field.
 */
std::string LevelSetSegmentationBaseObject::getVTKOutputFieldName( LevelSetField field) const
{
    switch (field) {

        case LevelSetField::SUPPORT:
            return "SupportId";

        case LevelSetField::PART:
            return "PartId";

        case LevelSetField::NORMAL:
            return "Normal";

        default:
            return LevelSetObject::getVTKOutputFieldName(field);

    }
}

/*!
 * Write the specified field to the given stream.
 *
 * @param[in] stream output stream
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 * @param[in] field is the field that will be written
 */
void LevelSetSegmentationBaseObject::flushVTKOutputData(std::fstream &stream, VTKFormat format,
                                                        LevelSetField field) const
{
    switch (field) {

    case LevelSetField::SUPPORT:
    {
        auto evaluator = [this] (long id) { return evalCellSupport(id); };
        auto fallback = [] (long id) { BITPIT_UNUSED(id); return levelSetDefaults::SUPPORT; };
        flushVTKOutputData<double>(stream, format, field, evaluator, fallback);
        break;
    }

    case LevelSetField::PART:
    {
        auto evaluator = [this] (long id) { return evalCellPart(id); };
        auto fallback = [] (long id) { BITPIT_UNUSED(id); return levelSetDefaults::PART; };
        flushVTKOutputData<double>(stream, format, field, evaluator, fallback);
        break;
    }

    case LevelSetField::NORMAL:
    {
        auto evaluator = [this] (long id) { return evalCellNormal(id, true); };
        auto fallback = [] (long id) { BITPIT_UNUSED(id); return levelSetDefaults::NORMAL; };
        flushVTKOutputData<double>(stream, format, field, evaluator, fallback);
        break;
    }

    default:
        LevelSetObject::flushVTKOutputData(stream, format, field);
        break;

    }
}

/*!
 * Get the part associated with the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetSegmentationBaseObject::getPart(long cellId) const
{
    return evalCellPart(cellId);
}

/*!
 * Get the normal of the surface at the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetSegmentationBaseObject::getNormal(long cellId) const
{
    return evalCellNormal(cellId, m_defaultSignedLevelSet);
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationBaseObject::getSupport(long cellId) const
{
    return evalCellSupport(cellId);
}

/*!
 * Get the size of the segment closest to the specified cell.
 *
 * \param cellId is the id of the cell
 * \result The size of the segment closest to the specified cell.
 */
double LevelSetSegmentationBaseObject::getSurfaceFeatureSize(long cellId) const
{
    long support = evalCellSupport(cellId);
    if (support < 0) {
        return (- levelSetDefaults::SIZE);
    }

    const SurfUnstructured &surface = evalCellSurface(cellId);

    return surface.evalCellSize(support);
}

/*!
 * Get the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetSegmentationBaseObject::getPart(const std::array<double,3> &point) const
{
    return evalPart(point);
}

/*!
 * Get the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetSegmentationBaseObject::getNormal(const std::array<double,3> &point) const
{
    return evalNormal(point, m_defaultSignedLevelSet);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationBaseObject::getSupport(const std::array<double,3> &point) const
{
    return evalSupport(point);
}

/*!
 * Get the size of the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The size of the segment closest to the specified point.
 */
double LevelSetSegmentationBaseObject::getSurfaceFeatureSize(const std::array<double,3> &point) const
{
    long support = evalSupport(point);
    if (support < 0) {
        return (- levelSetDefaults::SIZE);
    }

    const SurfUnstructured &surface = evalSurface(point);

    return surface.evalCellSize(support);
}

/*!
	@class      LevelSetSegmentationObject
	@ingroup    levelset
	@brief      Implements visitor pattern fo segmentated geometries
*/

/*!
 * Constructor
 * @param[in] id identifier of object
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id)
    : LevelSetSegmentationBaseObject(id),
      m_surfaceInfo(nullptr)
{
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] surface unique pointer to surface mesh
 * @param[in] featureAngle feature angle. If the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id, std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle)
    : LevelSetSegmentationBaseObject(id)
{
    setSurface(std::move(surface), featureAngle);
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] surface pointer to surface mesh
 * @param[in] featureAngle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id, const SurfUnstructured *surface, double featureAngle)
    : LevelSetSegmentationBaseObject(id)
{
    setSurface(surface, featureAngle);
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] surface pointer to surface mesh
 * @param[in] featureAngle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 * @param[in] surfaceSmoothing is the given surface snoothing order
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(int id, const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing)
    : LevelSetSegmentationBaseObject(id)
{
    setSurface(surface, featureAngle, surfaceSmoothing);
}

/*!
 * Copy constructor.
 *
 * \param other is another object whose content is copied in this object
 */
LevelSetSegmentationObject::LevelSetSegmentationObject(const LevelSetSegmentationObject &other)
    : LevelSetSegmentationBaseObject(other)
{
    if (other.m_surfaceInfo) {
        m_surfaceInfo = std::unique_ptr<LevelSetSegmentationSurfaceInfo>(new LevelSetSegmentationSurfaceInfo(*(other.m_surfaceInfo)));
    } else {
        m_surfaceInfo = nullptr;
    }
}

/*!
 * Checks if the object is empty.
 *
 * \result Returns true is the object is empty, false otherwise.
 */
bool LevelSetSegmentationObject::empty() const
{
    return getSurface().empty();
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetSegmentationObject * LevelSetSegmentationObject::clone() const {
    return new LevelSetSegmentationObject(*this );
}

/*!
 * Get segmentation surface
 * @return segmentation surface;
 */
const SurfUnstructured & LevelSetSegmentationObject::getSurface() const {
    return m_surfaceInfo->getSurface();
}

/*!
 * Set the surface
 *
 * Unless explicitly forced, it is not possible to replace an existing surface. Also, when the
 * surface is replaced, the object will not recalculate the levelset on the newly set surface
 * (nor will tell the proxy objects that may depend on the current object to update the
 * levelset values).
 *
 * The feature angle will be set to the defualt value specified by the constant
 * LevelSetSegmentationSurfaceInfo::DEFAULT_FEATURE_ANGLE.
 *
 * @param[in] surface is the surface that will be set
 * @param[in] force controls if it is possible to replace an existing surface.
 */
void LevelSetSegmentationObject::setSurface(std::unique_ptr<const SurfUnstructured> &&surface, bool force){
    setSurface(std::move(surface), LevelSetSegmentationSurfaceInfo::DEFAULT_FEATURE_ANGLE, force);
}

/*!
 * Set the surface
 *
 * Unless explicitly forced, it is not possible to replace an existing surface. Also, when the
 * surface is replaced, the object will not recalculate the levelset on the newly set surface
 * (nor will tell the proxy objects that may depend on the current object to update the
 * levelset values).
 *
 * @param[in] surface is the surface that will be set
 * @param[in] featureAngle is the angle that is used to identify sharp edges. If the angle between
 * two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 * @param[in] force controls if it is possible to replace an existing surface.
 */
void LevelSetSegmentationObject::setSurface(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle, bool force){
    if (m_surfaceInfo) {
        // Check if replacing an existing surface is allowed
        if (!force) {
            throw std::runtime_error ("The surface can only be set once.");
        }

        // Replace the surface
        m_surfaceInfo->setSurface(std::move(surface), featureAngle);
    } else {
        // Set surface
        //
        // Since this is the first time we set the surface, there is no need
        // to clear the caches.
        m_surfaceInfo = std::unique_ptr<LevelSetSegmentationSurfaceInfo>(new LevelSetSegmentationSurfaceInfo(std::move(surface), featureAngle));
    }
}

/*!
 * Set the surface
 *
 * Unless explicitly forced, it is not possible to replace an existing surface. Also, when the
 * surface is replaced, the object will not recalculate the levelset on the newly set surface
 * (nor will tell the proxy objects that may depend on the current object to update the
 * levelset values).
 *
 * The feature angle will be set to the defualt value specified by the constant
 * LevelSetSegmentationSurfaceInfo::DEFAULT_FEATURE_ANGLE.
 *
 * @param[in] surface is the surface that will be set
 * @param[in] force controls if it is possible to replace an existing surface.
 */
void LevelSetSegmentationObject::setSurface(const SurfUnstructured *surface, bool force){
    setSurface(surface, LevelSetSegmentationSurfaceInfo::DEFAULT_FEATURE_ANGLE, force);
}

/*!
 * Set the surface
 *
 * Unless explicitly forced, it is not possible to replace an existing surface. Also, when the
 * surface is replaced, the object will not recalculate the levelset on the newly set surface
 * (nor will tell the proxy objects that may depend on the current object to update the
 * levelset values).
 *
 * @param[in] surface is the surface that will be set
 * @param[in] featureAngle is the angle that is used to identify sharp edges. If the angle between
 * two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 * @param[in] force controls if it is possible to replace an existing surface.
 */
void LevelSetSegmentationObject::setSurface(const SurfUnstructured *surface, double featureAngle, bool force){
    if (m_surfaceInfo) {
        // Check if replacing an existing surface is allowed
        if (!force) {
            throw std::runtime_error ("The surface can only be set once.");
        }

        // Replace the surface
        m_surfaceInfo->setSurface(surface, featureAngle);
    } else {
        // Set surface
        //
        // Since this is the first time we set the surface, there is no need
        // to clear the caches.
        m_surfaceInfo = std::unique_ptr<LevelSetSegmentationSurfaceInfo>(new LevelSetSegmentationSurfaceInfo(surface, featureAngle));
    }
}

/*!
 * Set the surface
 *
 * Unless explicitly forced, it is not possible to replace an existing surface. Also, when the
 * surface is replaced, the object will not recalculate the levelset on the newly set surface
 * (nor will tell the proxy objects that may depend on the current object to update the
 * levelset values).
 *
 * @param[in] surface is the surface that will be set
 * @param[in] featureAngle is the angle that is used to identify sharp edges. If the angle between
 * two segments is bigger than this angle, the enclosed edge is considered as a sharp edge
 * @param[in] surfaceSmoothing is the given surface snoothing order
 * @param[in] force controls if it is possible to replace an existing surface.
 */
void LevelSetSegmentationObject::setSurface(const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing, bool force){
    if (m_surfaceInfo) {
        // Check if replacing an existing surface is allowed
        if (!force) {
            throw std::runtime_error ("The surface can only be set once.");
        }

        // Replace the surface
        m_surfaceInfo->setSurface(surface, featureAngle, surfaceSmoothing);
    } else {
        // Set surface
        //
        // Since this is the first time we set the surface, there is no need
        // to clear the caches.
        m_surfaceInfo = std::unique_ptr<LevelSetSegmentationSurfaceInfo>(new LevelSetSegmentationSurfaceInfo(surface, featureAngle, surfaceSmoothing));
    }
}

/*!
 * Get search tree
 * @return search tree;
 */
const SurfaceSkdTree & LevelSetSegmentationObject::getSearchTree() const {
    return m_surfaceInfo->getSearchTree();
}

/*!
 * Get feature angle
 * @return feature angle used when calculating face normals.
 */
double LevelSetSegmentationObject::getFeatureAngle() const {
    return m_surfaceInfo->getFeatureAngle();
}

/*!
 * Get smoothing order (low or high) imposed on surface when calculating the
 * projection point and normal on the surface
 * @return the ssurface smoothing order
 */
bitpit::LevelSetSurfaceSmoothing LevelSetSegmentationObject::getSurfaceSmoothing() const {
    return m_surfaceInfo->getSurfaceSmoothing();
}

/*!
 * Fill the cache that contains the zone associated to the cells.
 *
 * A cell can be either in the narrow band or in the bulk. It will be considered inside the narrow
 * band if one of the following conditions holds:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 */
void LevelSetSegmentationObject::fillCellLocationCache()
{
    // Cartesian patches are handled separately
    if (dynamic_cast<LevelSetCartesianKernel*>(m_kernel)) {
        fillCartesianCellZoneCache();
        return;
    }

    // All other patches are handled with the base method.
    LevelSetObject::fillCellLocationCache();
}

/*!
 * Fill the cache that contains the zone associated to the cells.
 *
 * A cell can be either in the narrow band or in the bulk. It will be considered inside the narrow
 * band if one of the following conditions holds:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 *
 * \param adaptionData are the information about the mesh update
 */
void LevelSetSegmentationObject::fillCellLocationCache(const std::vector<adaption::Info> &adaptionData)
{
    // Cartesian patches are handled separately
    //
    // Update is not implemented for Cartesian patches, the cells that should be inserted in
    // the cache will be evaluated considering all the mash not just the elements newly added.
    if (dynamic_cast<LevelSetCartesianKernel*>(m_kernel)) {
        fillCartesianCellZoneCache();
        return;
    }

    // All other patches are handled with the base method
    LevelSetObject::fillCellLocationCache(adaptionData);
}

/*!
 * Fill the cache that contains the zone associated to the cells of a Cartesian patch.
 *
 * A cell is considered inside the narrow band if one of the following conditions hold:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 *
 * \result The ids of cells that should be inserted in a cache operating "narrow band" mode.
 */
void LevelSetSegmentationObject::fillCartesianCellZoneCache()
{
    // The function needs a Cartesian kernel
    const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<LevelSetCartesianKernel*>(m_kernel);
    if (!dynamic_cast<LevelSetCartesianKernel*>(m_kernel)) {
        throw std::runtime_error("The function needs a Cartesian kernels.");
    }

    // Get mesh information
    const VolCartesian &mesh = *(cartesianKernel->getMesh() ) ;
    int meshDimension = mesh.getDimension();
    long nCells = mesh.getCellCount();

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

    // Get cache for zone identification
    CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);

    // Cells with an unknown region are in the bulk
    for (long cellId = 0; cellId < nCells; ++cellId) {
        locationCache->insertEntry(cellId, static_cast<char>(LevelSetCellLocation::UNKNOWN));
    }

    // Start filling the list of cells within the narrow band
    //
    // The initial process list is gradually expanded considering all the neighbours
    // inside the narrow band.
    std::vector<long> intersectedCellIds;
    std::vector<bool> alreadyProcessed(nCells, false);
    while (!processList.empty()) {
        // Get the cell to process
        long cellId = *(processList.begin());
        processList.erase(processList.begin());

        // Skip cells already processed
        if (alreadyProcessed[cellId]) {
            continue;
        }
        alreadyProcessed[cellId] = true;

        // Fill location cache for cells geometrically inside the narrow band
        LevelSetCellLocation cellLocation = fillCellGeometricNarrowBandLocationCache(cellId);

        // Track intersected cells
        if (cellLocation == LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
            intersectedCellIds.push_back(cellId);
        }

        // Add unprocessed neighbours of cells inside the narrow band to the process list
        if (cellLocation != LevelSetCellLocation::UNKNOWN) {
            auto neighProcessor = [&processList, &alreadyProcessed](long neighId, int layer) {
                BITPIT_UNUSED(layer);

                // Skip neighbours already processed
                if (alreadyProcessed[neighId]) {
                    return false;
                }

                // Update the process list
                processList.insert(neighId);

                // Continue processing the other neighbours
                return false;
            };

            mesh.processCellFaceNeighbours(cellId, 1, neighProcessor);
        }
    }

    // Identify neighbours of cells that intersect the surface
    for (std::size_t cellId : intersectedCellIds) {
        // Process face neighbours
        auto neighProcessor = [this, &locationCache](long neighId, int layer) {
            BITPIT_UNUSED(layer);

            // Skip neighbours whose region has already been identified
            if (getCellLocation(neighId) != LevelSetCellLocation::UNKNOWN) {
                return false;
            }

            // The neighbour is inside the narrow band
            locationCache->insertEntry(neighId, static_cast<char>(LevelSetCellLocation::NARROW_BAND_NEIGHBOUR));

            // Continue processing the other neighbours
            return false;
        };

        mesh.processCellFaceNeighbours(cellId, 1, neighProcessor);
    }

    // Cells with an unknown region are in the bulk
    for (long cellId = 0; cellId < nCells; ++cellId) {
        CellCacheCollection::ValueCache<char>::Entry locationCacheEntry = locationCache->findEntry(cellId);
        if (*locationCacheEntry == static_cast<char>(LevelSetCellLocation::UNKNOWN)) {
            locationCache->insertEntry(cellId, static_cast<char>(LevelSetCellLocation::BULK));
        }
    }

#if BITPIT_ENABLE_MPI==1
    // Exchange ghost data
    if (mesh.isPartitioned()) {
        std::unique_ptr<DataCommunicator> dataCommunicator = m_kernel->createDataCommunicator();
        startCellCacheExchange(mesh.getGhostCellExchangeSources(), m_cellLocationCacheId, dataCommunicator.get());
        completeCellCacheExchange(mesh.getGhostCellExchangeTargets(), m_cellLocationCacheId, dataCommunicator.get());
    }
#endif
}

/*!
 * Fill location cache for the specified cell if it is geometrically inside the narrow band
 *
 * A cell is geometrically inside the narrow band if its distance from the surface is smaller
 * than the narrow band side or if it intersects the surface.
 *
 * This function may require the evaluation of some levelset fields. To improve performance,
 * it is important to attempt filling the cache of the evaluated fields. It is then up to the
 * caches to decide if the fields can be cached.
 *
 * \param[in] id is the cell id
 * \return The cell location of the cache or LevelSetCellLocation::UNKNOWN if the cell
 * location was not identified.
 */
LevelSetCellLocation LevelSetSegmentationObject::fillCellGeometricNarrowBandLocationCache(long id)
{
    // Early return if the cell is geometrically outside the narrow band
    double searchRadius = std::max(m_kernel->computeCellBoundingRadius(id), m_narrowBandSize);
    long cellSupport = evalCellSupport(id, searchRadius);
    if (cellSupport < 0) {
        return LevelSetCellLocation::UNKNOWN;
    }

    // Evaluate levelset value
    std::array<double,3> cellCentroid = m_kernel->computeCellCentroid(id);

    double cellCacheValue    = _evalValue(cellCentroid, cellSupport, CELL_CACHE_IS_SIGNED);
    double cellUnsigendValue = std::abs(cellCacheValue);

    // Update the cell location cache
    //
    // First we need to check if the cell intersectes the surface, and only if it
    // deosn't we should check if its distance is lower than the narrow band size.
    LevelSetCellLocation cellLocation = LevelSetCellLocation::UNKNOWN;
    if (_intersectSurface(id, cellUnsigendValue, CELL_LOCATION_INTERSECTION_MODE) == LevelSetIntersectionStatus::TRUE) {
        cellLocation = LevelSetCellLocation::NARROW_BAND_INTERSECTED;
    } else if (cellUnsigendValue <= m_narrowBandSize) {
        cellLocation = LevelSetCellLocation::NARROW_BAND_DISTANCE;
    }
    assert(cellLocation != LevelSetCellLocation::UNKNOWN);

    CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);
    locationCache->insertEntry(id, static_cast<char>(cellLocation));

    // Fill the cache of the evaluated fields
    //
    // Now that the cell location has been identified, we can fill the cache of the
    // evaluated fields.
    fillFieldCellCache(LevelSetField::SUPPORT, id, cellSupport);
    fillFieldCellCache(LevelSetField::VALUE, id, cellCacheValue);

    // Return the location
    return cellLocation;
}

/*!
 * Evaluate the surface associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The surface associated with the segment closest to the specified cell.
 */
const SurfUnstructured & LevelSetSegmentationObject::_evalCellSurface(long id) const
{
    BITPIT_UNUSED(id);

    return getSurface();
}

/*!
 * Evaluate levelset sign at the specified cell.
 *
 * \param id is the id of the cell
 * \result The sign of the levelset at the specified cell.
 */
short LevelSetSegmentationObject::_evalCellSign(long id) const
{
    long support = evalCellSupport(id);
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalSign(centroid, support);
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
    long support = evalCellSupport(id);
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalValue(centroid, support, signedLevelSet);
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
    long support = evalCellSupport(id);
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalGradient(centroid, support, signedLevelSet);
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support. If the search radius is set equal to the
 * constant AUTOMATIC_SEARCH_RADIUS, the object will try to evaluate the optional search radius
 * for the specified cell. The automatic evaluation of the search radius is possible only for
 * a limited number of cases, when the automatic evaluation cannot be performed, an infinite
 * search radius will be used.
 * \result The segment closest to the specified cell.
 */
long LevelSetSegmentationObject::_evalCellSupport(long id, double searchRadius) const
{
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalSupport(centroid, searchRadius);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetSegmentationObject::_evalCellNormal(long id, bool signedLevelSet) const
{
    long support = evalCellSupport(id);
    std::array<double,3> centroid = m_kernel->computeCellCentroid(id);

    return _evalNormal(centroid, support, signedLevelSet);
}

/*!
 * Evaluate the surface associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The surface associated with the segment closest to the specified point.
 */
const SurfUnstructured & LevelSetSegmentationObject::_evalSurface(const std::array<double,3> &point) const
{
    BITPIT_UNUSED(point);

    return getSurface();
}

/*!
 * Evaluate levelset sign at the specified point.
 *
 * \param point are the coordinates of the point
 * \result The sign of the levelset at the specified point.
 */
short LevelSetSegmentationObject::_evalSign(const std::array<double,3> &point) const
{
    long support = evalSupport(point);

    return _evalSign(point, support);
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

    return _evalValue(point, support, signedLevelSet);
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

    return _evalGradient(point, support, signedLevelSet);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetSegmentationObject::_evalNormal(const std::array<double,3> &point, bool signedLevelSet) const
{
    long support = evalSupport(point);

    return _evalNormal(point, support, signedLevelSet);
}

/*!
 * Evaluate levelset sign at the specified point.
 *
 * \param point are the coordinates of the point
 * \param support is the the closest segment to the specified point
 * \result The sign of the levelset at the specified point.
 */
short LevelSetSegmentationObject::_evalSign(const std::array<double,3> &point, long support) const
{
    // Throw an error if the support is not valid
    if (support < 0) {
        if (empty()) {
            return levelSetDefaults::SIGN;
        }

        throw std::runtime_error("Unable to evaluate the sign: the support is not valid.");
    }

    LevelSetSegmentationSurfaceInfo::SegmentConstIterator supportItr = getSurface().getCellConstIterator(support);

    return evalValueSign(m_surfaceInfo->evalDistance(point, supportItr, true));
}

/*!
 * Evaluate levelset value at the specified point.
 *
 * \param point are the coordinates of the point
 * \param support is the the closest segment to the specified point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified point.
 */
double LevelSetSegmentationObject::_evalValue(const std::array<double,3> &point, long support,
                                              bool signedLevelSet) const
{
    // Early return if the support is not valid
    //
    // With an invalid support, only the unsigend levelset can be evaluated.
    if (support < 0) {
        if (!signedLevelSet || empty()) {
            return levelSetDefaults::VALUE;
        }

        throw std::runtime_error("With an invalid support, only the unsigend levelset can be evaluated.");
    }

    // Evaluate the distance of the point from the surface
    LevelSetSegmentationSurfaceInfo::SegmentConstIterator supportItr = getSurface().getCellConstIterator(support);
    double distance = m_surfaceInfo->evalDistance(point, supportItr, signedLevelSet);

    // Early return if the point lies on the surface
    if (evalValueSign(distance) == 0) {
        return 0.;
    }

    // Evaluate levelset value
    double value;
    if (signedLevelSet) {
        value = distance;
    } else {
        value = std::abs(distance);
    }

    return value;
}

/*!
 * Evaluate levelset gradient at the specified point.
 *
 * \param point are the coordinates of the point
 * \param support is the the closest segment to the specified point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The gradient of the levelset at the specified point.
 */
std::array<double,3> LevelSetSegmentationObject::_evalGradient(const std::array<double,3> &point, long support,
                                                               bool signedLevelSet) const
{
    // Early return if the support is not valid
    //
    // With an invalid support, only the unsigend levelset can be evaluated.
    if (support < 0) {
        if (!signedLevelSet || empty()) {
            return levelSetDefaults::GRADIENT;
        }

        throw std::runtime_error("With an invalid support, only the unsigend levelset can be evaluated.");
    }

    // Evaluate the distance of the point from the surface
    LevelSetSegmentationSurfaceInfo::SegmentConstIterator supportItr = getSurface().getCellConstIterator(support);
    double distance = m_surfaceInfo->evalDistance(point, supportItr, signedLevelSet);

    // Early return if the point lies on the surface
    if (evalValueSign(distance) == 0) {
        if (signedLevelSet) {
            return m_surfaceInfo->evalNormal(point, supportItr);
        } else {
            return {{0., 0., 0.}};
        }
    }

    // Evaluate levelset gradient
    std::array<double,3> gradient = m_surfaceInfo->evalDistanceVector(point, supportItr) / distance;

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
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The closest segment to the specified point.
 */
long LevelSetSegmentationObject::_evalSupport(const std::array<double,3> &point, double searchRadius) const
{
    long closestSegmentId;
    double closestDistance;
    m_surfaceInfo->getSearchTree().findPointClosestCell(point, searchRadius, &closestSegmentId, &closestDistance);

    return closestSegmentId;
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \param support is the the closest segment to the specified point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetSegmentationObject::_evalNormal(const std::array<double,3> &point, long support,
                                                             bool signedLevelSet) const
{
    // Early return if the support is not valid
    //
    // With an invalid support, only the unsigend levelset can be evaluated.
    if (support < 0) {
        if (!signedLevelSet || empty()) {
            return levelSetDefaults::NORMAL;
        }

        throw std::runtime_error("With an invalid support, only the unsigend levelset can be evaluated.");
    }

    // Evaluate the normal
    //
    // If an unsigned evaluation is requested, the orientation of the surface should be discarded
    // and in order to have a normal that is agnostic with respect the two sides of the surface.
    LevelSetSegmentationSurfaceInfo::SegmentConstIterator supportItr = getSurface().getCellConstIterator(support);
    std::array<double,3> normal = m_surfaceInfo->evalNormal(point, supportItr);
    if (!signedLevelSet) {
        normal *= static_cast<double>(evalSign(point));
    }

    return normal;
}

/*!
 * Get the smallest characteristic size within the triangulation
 * This function is only provided for guarantee backwards compatibility with older versions.
 * It is out of the levelset scope to evaluate the feature size of the surface.
 * @return smallest characteristic size within the triangulation
 */
double LevelSetSegmentationObject::getMinSurfaceFeatureSize() const {

    const SurfUnstructured &surface = getSurface();

    bool   minimumValid = false;
    double minimumSize  = levelSetDefaults::SIZE;
    for (const Cell &cell : surface.getCells()) {
        double segmentSize = surface.evalCellSize(cell.getId());
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
 * Get the largest characteristic size within the triangulation.
 * This function is only provided for guarantee backwards compatibility with older versions.
 * It is out of the levelset scope to evaluate the feature size of the surface.
 * @return largest characteristic size within the triangulation
 */
double LevelSetSegmentationObject::getMaxSurfaceFeatureSize() const {

    const SurfUnstructured &surface = getSurface();

    double maximumSize = - levelSetDefaults::SIZE;
    for (const Cell &cell : surface.getCells()) {
        double segmentSize = surface.evalCellSize(cell.getId());
        maximumSize = std::max(segmentSize, maximumSize);
    }

    return maximumSize;
}

/*!
 * Constructor taking two objects.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] source1 pointer to first source object
 * @param[in] source2 pointer to second source object
 */
LevelSetBooleanObject<LevelSetSegmentationBaseObject>::LevelSetBooleanObject( int id, LevelSetBooleanOperation op, const LevelSetSegmentationBaseObject *source1, const LevelSetSegmentationBaseObject *source2  )
    : LevelSetBooleanBaseObject<LevelSetSegmentationBaseObject>(id, op, source1, source2) {
}

/*!
 * Constructor taking a vector of objects.
 * The boolean operation will be applied recursively on each entry.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] sourceObjects pointers to source objects
 */
LevelSetBooleanObject<LevelSetSegmentationBaseObject>::LevelSetBooleanObject( int id, LevelSetBooleanOperation op, const std::vector<const LevelSetSegmentationBaseObject *> &sourceObjects )
    : LevelSetBooleanBaseObject<LevelSetSegmentationBaseObject>(id, op, sourceObjects) {
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetBooleanObject<LevelSetSegmentationBaseObject> * LevelSetBooleanObject<LevelSetSegmentationBaseObject>::clone() const {
    return new LevelSetBooleanObject<LevelSetSegmentationBaseObject>(*this );
}

/*!
 * Evaluate the source associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The source associated with the segment closest to the specified cell.
 */
const SurfUnstructured & LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalCellSurface(long id) const
{
    return getCellReferenceObject(id)->evalCellSurface(id);
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The segment closest to the specified cell.
 */
long LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalCellSupport(long id, double searchRadius) const
{
    return getCellReferenceObject(id)->evalCellSupport(id, searchRadius);
}

/*!
 * Evaluate the part associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalCellPart(long id) const
{
    return getCellReferenceObject(id)->evalCellPart(id);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalCellNormal(long id, bool signedLevelSet) const
{
    return _evalCellFunction<std::array<double,3>>(id, signedLevelSet, [&id, signedLevelSet] (const LevelSetBooleanResult<LevelSetSegmentationBaseObject> &result)
        {
            const LevelSetSegmentationBaseObject *resultObject = result.getObject();
            if ( !resultObject ) {
                return levelSetDefaults::NORMAL;
            }

            std::array<double,3> normal = resultObject->evalCellNormal(id, signedLevelSet);
            if (signedLevelSet) {
                normal *= static_cast<double>(result.getObjectSign());
            }

            return normal;
        });
}

/*!
 * Evaluate the surface associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The surface associated with the segment closest to the specified point.
 */
const SurfUnstructured & LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalSurface(const std::array<double,3> &point) const
{
    return getReferenceObject(point)->evalSurface(point);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalSupport(const std::array<double,3> &point) const
{
    return getReferenceObject(point)->evalSupport(point);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The closest segment to the specified point.
 */
long LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalSupport(const std::array<double,3> &point, double searchRadius) const
{
    return getReferenceObject(point)->evalSupport(point, searchRadius);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalNormal(const std::array<double,3> &point, bool signedLevelSet) const
{
    return _evalFunction<std::array<double,3>>(point, signedLevelSet, [&point, signedLevelSet] (const LevelSetBooleanResult<LevelSetSegmentationBaseObject> &result)
        {
            const LevelSetSegmentationBaseObject *resultObject = result.getObject();
            if ( !resultObject ) {
                return levelSetDefaults::NORMAL;
            }

            std::array<double,3> normal = resultObject->evalNormal(point, signedLevelSet);
            if (signedLevelSet) {
                normal *= static_cast<double>(result.getObjectSign());
            }

            return normal;
        });
}

/*!
 * Evaluate the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalPart(const std::array<double,3> &point) const
{
    return getReferenceObject(point)->evalPart(point);
}

/*!
 * Constructor.
 *
 * \param[in] id identifier of object
 * \param[in] source pointer to source object
 */
LevelSetComplementObject<LevelSetSegmentationBaseObject>::LevelSetComplementObject(int id, const LevelSetSegmentationBaseObject *source)
    : LevelSetComplementBaseObject<LevelSetSegmentationBaseObject>(id, source)
{
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetComplementObject<LevelSetSegmentationBaseObject> * LevelSetComplementObject<LevelSetSegmentationBaseObject>::clone() const {
    return new LevelSetComplementObject<LevelSetSegmentationBaseObject>(*this );
}

/*!
 * Evaluate the surface associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The surface associated with the segment closest to the specified cell.
 */
const SurfUnstructured & LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalCellSurface(long id) const
{
    return getSourceObject()->evalCellSurface(id);
}

/*!
 * Evaluate the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support. If the search radius is set equal to the
 * constant AUTOMATIC_SEARCH_RADIUS, the object will try to evaluate the optional search radius
 * for the specified cell. The automatic evaluation of the search radius is possible only for
 * a limited number of cases, when the automatic evaluation cannot be performed, an infinite
 * search radius will be used.
 * \result The segment closest to the specified cell.
 */
long LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalCellSupport(long id, double searchRadius) const
{
    return getSourceObject()->evalCellSupport(id, searchRadius);
}

/*!
 * Evaluate the part associated with the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \result The part associated with the segment closest to the specified cell.
 */
int LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalCellPart(long id) const
{
    return getCellReferenceObject(id)->evalCellPart(id);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified cell.
 */
std::array<double,3> LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalCellNormal(long id, bool signedLevelSet) const
{
    std::array<double,3> normal = getSourceObject()->evalCellNormal(id, signedLevelSet);
    if (signedLevelSet) {
        normal *= -1.;
    }

    return normal;
}

/*!
 * Evaluate the surface associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The surface associated with the segment closest to the specified point.
 */
const SurfUnstructured & LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalSurface(const std::array<double,3> &point) const
{
    return getSourceObject()->evalSurface(point);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The closest segment to the specified point.
 */
long LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalSupport(const std::array<double,3> &point) const
{
    return getSourceObject()->evalSupport(point);
}

/*!
 * Evaluate the closest segment to the specified point.
 *
 * \param point are the coordinates of the point
 * \param searchRadius all segments whose distance is greater than the search radius will not
 * be considered for the evaluation of the support
 * \result The closest segment to the specified point.
 */
long LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalSupport(const std::array<double,3> &point, double searchRadius) const
{
    return getSourceObject()->evalSupport(point, searchRadius);
}

/*!
 * Evaluate the part associated with the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \result The part associated with the segment closest to the specified point.
 */
int LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalPart(const std::array<double,3> &point) const
{
    return getSourceObject()->evalPart(point);
}

/*!
 * Evaluate the normal of the surface at the segment closest to the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The normal of the surface at the segment closest to the specified point.
 */
std::array<double,3> LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalNormal(const std::array<double,3> &point, bool signedLevelSet) const
{
    std::array<double,3> normal = getSourceObject()->evalNormal(point, signedLevelSet);
    if (signedLevelSet) {
        normal *= -1.;
    }

    return normal;
}

}
