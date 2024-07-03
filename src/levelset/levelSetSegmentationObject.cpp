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
 * @param[out] the distance function at the specified point
 * @return The distance function at the specified point.
 */
double LevelSetSegmentationSurfaceInfo::evalDistance(const std::array<double, 3> &point,
                                                     const SegmentConstIterator &segmentItr,
                                                     bool signedDistance,
                                                     std::array<double, 3> *distanceVector) const
{
    // Project the point on the surface and evaluate the point-projection vector
    int nSegmentVertices = segmentItr->getVertexCount();
    BITPIT_CREATE_WORKSPACE(lambda, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    std::array<double, 3> pointProjection = evalProjection(point, segmentItr, lambda);
    (*distanceVector) = point - pointProjection;

    // Evaluate unsigned distance
    double unsignedDistance = norm2(*distanceVector);
    if (!signedDistance) {
        return unsignedDistance;
    }

    // Signed distance
    //
    // If the sign is null and the point doesn't lie on the segmentation, it lies on the normal
    // plane. This case is not supported, because it would require to evaluate the sign taking
    // into account the the curvature of the surface.
    std::array<double, 3> pseudoNormal = computePseudoNormal(segmentItr, lambda);
    double pointProjectionNormalComponent = dotProduct(*distanceVector, pseudoNormal);

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
    std::array<double, 3> distanceVector;
    return evalDistance(point, segmentItr, signedDistance, &distanceVector);
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
    std::array<double, 3> projectionPoint;
    std::array<double, 3> projectionNormal;
    evalProjection(point, segmentItr, &projectionPoint, &projectionNormal);

    return projectionNormal;
}

/*!
 * Compute the pseudo-normal at specified point of the given segment.
 *
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[in] lambda are the barycentric coordinates of the point
 * @return the pseudo-normal at specified point of the given segment
 */
std::array<double,3> LevelSetSegmentationSurfaceInfo::evalPseudoNormal(const SegmentConstIterator &segmentItr,
                                                                       const double *lambda ) const
{
    return computePseudoNormal(segmentItr, lambda);
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
 * Evaluate the projection of the given point on a curve created based on
 * the points representing the specified segment. The curve passes from these
 * points and is vertical to the normal vectors associated with them. It's
 * continuous with a continuous normal vector distribution.
 *
 * Consider an arbitrary segment in 2D space defined by nodes "point0" and "poin1",
 * its normal vector "normal_s" and an arbitrary point in space named "point".
 * Let "point_s" be the projection of "point" to the segment
 * and "t[0]" and "t[1] its barycentric coordinates such that:
 * point_s = tau[0] * point0 + tau[1] * point1
 *
 * The constructed curve is described by equation:
 * projectionPoint = point_s + f01 * normal_s
 *
 * Variable "f01" is function of "t[0]" and  guaranties that the curve passes from
 * "point0" and "point1". Thus, f01(0)=0 and f01(1)=0. Also, it guarantees that
 * the curve's tangent computed as
 * d(projectionPoint)/d(t[0]) = point0 - point1 + d(f10)/d(t[0]) * normal_s
 * is vertical to the normals "normal0" and "normal1".
 *
 * Let's introduce variables "lambda10" and "lambda01" defined as the derivative
 * of "f01" for t[0]=0 and t[0]=1, respectively. They are computed as follows:
 *
 * For t[0] = 0: d(projectionPoint)/d(t[0]) * normal1 = 0 <=>
 * lambda10 = ((point1 - point0) * normal1) / (normal_s * normal1)
 * For t[0] = 1: d(projectionPoint)/d(t[0]) * normal0 = 0 <=>
 * lambda01 = ((point0 - point1) * normal0) / (normal_s * normal0)
 *
 * Function "f01" is defined as a polynomial which respects all four conditions.
 * Thus,
 * f01 = t[0] * t[1] * (lambda01 * t[0] + lambda10 * t[1])
 *
 * The normal vector to the curve is defined as the derivative of its
 * tangential vector wrt "t[0]". After some math, its expression is
 * projectionNormal = d(f10)/d(t[0]) * (point1 - point0) + l * l * normal_s,
 * where "l" is the segment length.
 *
 * Since the normals corresponding to the geometry points are equal for two
 * subsequent segments, the developed curve is equipped with a continuous normal
 * vector everywhere.
 *
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] segmentItr is an iterator pointing to the segment on which the surface
 * will be created.
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 * projection point.
 */
void LevelSetSegmentationSurfaceInfo::evalHighOrderProjectionOnLine(const std::array<double,3> &point,
                                                                    const SegmentConstIterator &segmentItr,
                                                                    std::array<double, 3> *projectionPoint,
                                                                    std::array<double, 3> *projectionNormal) const
{
    // Get segment
    const Cell &segment = *segmentItr;
    assert(segment.getType() == ElementType::LINE);

    // Get segment vertices
    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds();

    const std::array<double,3> &point0 = m_surface->getVertexCoords(segmentVertexIds[0]);
    const std::array<double,3> &point1 = m_surface->getVertexCoords(segmentVertexIds[1]);

    // Get projection point on segment
    std::array<double,2> t;
    std::array<double, 3> point_s = CGElem::projectPointSegment(point, point0, point1, &(t[0]));

    // Early return if point_s consides with segment's node
    double distanceTolerance = m_surface->getTol();
    for (int i = 0; i < 2; ++i) {
        if (utils::DoubleFloatingEqual()(t[i], 1., distanceTolerance, distanceTolerance)) {
            (*projectionPoint)  = m_surface->getVertexCoords(segmentVertexIds[i]);
            (*projectionNormal) = computeSegmentVertexNormal(segmentItr, i, true);
            return;
        }
    }

    // Get normal on segment
    std::array<double, 3> facetNormal = computeSegmentNormal(segmentItr);
    std::array<double, 3> normal_s = point - point_s;

    double distance = norm2(normal_s);
    if (utils::DoubleFloatingEqual()(distance, 0., distanceTolerance, distanceTolerance)) {
        normal_s = facetNormal;
    } else {
        normal_s /= distance;
        if (dotProduct(normal_s, facetNormal) < 0.0) {
            normal_s *= -1.0;
        }
    }

    // Get normals on segment's vertices
    std::array<double, 3> normal0 = computeSegmentVertexNormal(segmentItr, 0, true);
    std::array<double, 3> normal1 = computeSegmentVertexNormal(segmentItr, 1, true);

    // Compute lambda coefficients
    std::array<double, 3> edge = point1 - point0;
    double lambda01 = -dotProduct(edge, normal0) / dotProduct(normal_s, normal0); 
    double lambda10 =  dotProduct(edge, normal1) / dotProduct(normal_s, normal1); 

    // Eval polynomial parameterizing the curve
    double product1 = t[0] * t[1];
    double product2 = lambda01 * t[0] + lambda10 * t[1];
    double f        = product1 * product2;
    double df       = product1 * (lambda01 - lambda10) + product2 * (t[1] - t[0]);

    // Eval projection point on curve
    (*projectionPoint) = point_s + f * normal_s;

    // Eval normal vector on curve
    double length = norm2(edge);
    (*projectionNormal)  = df * edge + length * length * normal_s;
    (*projectionNormal) /= norm2((*projectionNormal));

    if (dotProduct(normal_s, (*projectionNormal)) < 0.0) {
       (*projectionNormal) *= -1.0;
    }
}

/*!
 * Evaluate the projection of the given point on a smooth surface created 
 * over a discretized geometry and, referred as global surface. The surface
 * passes from the geometry's points and it's vertical to the normal vectors
 * associated with them. Also, it's continuous and offers normal vector
 * continuity over the different segments of the discretized geometry.
 *
 * In this case the studied segment is a triangle.
 *
 * The developed algorithm does not exist in literature and is described bellow.
 * Consider an arbitrary triangle in 3D space defined by nodes "point0", "poin1"
 * and "point2", its normal vector "normal_s" and an arbitrary point in space
 * named "point". Let "point_s" be the projection of "point" to the triangle
 * and "tau[i]" with i=0,1,2 be its barycentric coordinates such that:
 * point_s = tau[0] * point0 + tau[1] * point1 + tau[2] * point2
 *
 * Let "edgePoint_s[e]" be the projection of "point_s" to the triangle edge "e".
 * It can be described as a function of the two nodes "node0" and "node1"
 * defining the edge as:
 * edgePoint_s[e] = p10 * node0 + p01 * node1,
 * where "p01" and "p10" are the edge barycentric coordinates.
 *
 * Let "edgeNormal_s[e]" be the normal to the triangulated geometry corresponding
 * to "edgePoint_s[e]". This normal remains the same for both triangles sharing
 * the same edge.
 *
 * For each edge "e" a curve is constructed described by equation:
 * edgePoint_p[e] = edgePoint_s[e] + f01 * edgeNormal_s,
 * where "edgePoint_P[e]" is the curve point corresponding to "edgePoint_s[e]".
 *
 * Variable "f01" is function of "p10" and  guaranties that the curve passes from
 * "node0" and "node1". Thus, f01(0)=0 and f01(1)=0. Also, it guararantees that
 * the curve's tangent computed as
 * d(edgePoint_p[e])/d(p10) = node0 - node1 + d(f10)/d(p10) * edgeNormal_s
 * is vertical to the normals "normal0" and "normal1" which correspond to
 * nodes "node0" and "node1". 
 *
 * Let's introduce variables "larmda10" and "lambda01" defined as the derivative
 * of "f01" for p10=0 and p10=1, respectively. They are computed as follows:
 *
 * For p10 = 0: d(edgePoint_p[e])/d(p10) * normal1 = 0 <=>
 * lambda10 = ((node1 - node0) * normal1) / (edgeNormal_s[e] * normal1)
 * For p10 = 1: d(edgePoint_p[e])/d(p10) * normal0 = 0 <=>
 * lambda01 = ((node0 - node1) * normal0) / (edgeNormal_s[e] * normal0)
 *
 * Function "f01" is defined as a polynomial which respects all four conditions.
 * Thus,
 * f01 = p10 * p01 * (lambda01 * p10 + lambda10 * p01)
 *
 * The normal vector to the curve is defined as the derivative of its
 * tangential vector wrt "p10". After some math, its expression is
 * edgeNormal_p[e] = d(f10)/d(p10) * (point1 - point0) + l * l * edgeNormal_s,
 * where "l" is the edge length.
 *
 * The three points "edgePoint_s[e]", e=0,1,2 define a new triangle in
 * space named inscribed triangle.
 * Similarly, the three points "edgePoint_p[e]", e=0,1,2 define a new
 * triangle in space named projection triangle.
 *
 * Point "point_s" can be written by using the barycentric coordinates of
 * the inscribed triangle, namely "t[i]", i=0,1,2:
 * point_s = t[0] * edgePoint_s[0] + t[1] * edgePoint_s[1] + t[2] * edgePoint_s[2]
 * This point is transformed to the projection triangle by the formula:
 * point_m = t[0] * edgePoint_p[0] + t[1] * edgePoint_p[1] + t[2] * edgePoint_p[2]
 * This point is used to compute the projection point on the global surface as
 * projectionPoint = point_m + f * normal_s;
 *
 * Variable "f" is a function of "t[0]" and "t[1]" and it's the expansion of "f01"
 * in 2D. It's responsible for guaranteeing that the projection surface created
 * over the projection triangle will pass from its nodes "edgePoint_p" and will be
 * vertical to the corresponding normals "edgeNormal_p". It's a 2D polynomial
 * defined as
 * f(t[0], t[1]) = sum i=0,2 {
 *                     t[i] * t[j] * (lambda_ij * t[i] + lambda_ji * t[j])
 *               },
 * where "j" is the subsequent vertex of "i", i.e. j=(i+1)%3.
 *
 * The above function satisfies the following conditions:
 * f(0,0) = 0, f(1,0) = 0, f(0,1) = 0
 * meaning that the projection surface passes from the projection triangle's vertices.
 *
 * Also, the two tangents of the projection surface at each projection triangle
 * vertex should be vertical to the corresponding normal. Thus, there
 * are 6 conditions to be satisfied and 6 lambda_ij to be found. If "k" is the third
 * triangle vertex and t[k] = 0, the "lamda_ij" and "lamda_ji" can be found as follows:
 * for t[j] = 0: d(projectionPoint)/d(t[i]) * edgeNormal_p[i] = 0 <=>
 * lambda_ij = (edgePoint_p[i] - edgePoint_p[j] * edgeNormal[i]) / (normal_s * edgeNormal[i])
 * for t[i] = 0: d(projectionPoint)/d(t[j]) * edgeNormal_p[j] = 0 <=>
 * lambda_ji = (edgePoint_p[j] - edgePoint_p[i] * edgeNormal[j]) / (normal_s * edgeNormal[j])
 *
 * The global surface normal to the projectionPoint is computed by considering
 * that "ptojectionPoint" is a function of "tau[0]" and "tau[1]". The unit normal
 * is perpendicular to the two global surface tangent vectors and thus
 * projectionNormal ~ d(projectionPoint)/d(tau[0]) x d(projectionPoint)/d(tau[1])
 * Its direction is corrected to be aligned with "normal_s", hence
 * projectionNormal * normal_s >= 0
 *
 * The above algorithm creates a global surface with the following properties.
 * Firstly, it is clarified that the projection surface is not part of the
 * global surface. This is because for each "point_s" a different inscribed and
 * projection triangle will be created and a different projection surface will
 * be computed. As point_s moves towards an edge "e", the corresponding projection
 * surface will always pass from the "edgePoint_p[e]" and will be always vertical
 * to the "edgeNormal_p[e]". So, when point_s arrives at "e" it will be equal
 * to "edgePoint_p[e]" and it's normal will coincide with the "edgeNormal_p[e]".
 * Thus, the curves constructed for each edge are part of the global smoothed
 * surface over the geometry and the curve normal vectors are respected by the
 * global surface. Since these curves are the same for each triangle neighboring
 * the edge, global surface is continuous and offers a continuous normal over
 * adjacent triangles. Also, since the curves respect the vertices and normals
 * of the discretized geometry, the final surface will respect them as well.
 *
 * Used terminology:
 * - var_s: is any point or vector named "var" positioned on the interior or
 *   boundary of the given triangle represented by segmentItr which is part of
 *   the triangulated geometrical surface
 * - var_p: is any point or vector named "var" positioned on the interior or
 *   boundary of the projection triangle defined above
 * - d0_var : is the derivative of an arbitrary variable "var" wrt the first
 *   barycentric coordinate of the triangle (i.e. "tau[0]")
 * - d1_var : is the derivative of an arbitrary variable "var" wrt the first
 *   barycentric coordinate of the triangle (i.e. "tau[1]")
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] segmentItr is an iterator pointing to the segment on which the surface
 * will be created.
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 * projection point.
 */
void LevelSetSegmentationSurfaceInfo::evalHighOrderProjectionOnTriangle(const std::array<double,3> &point,
                                                                        const SegmentConstIterator &segmentItr,
                                                                        std::array<double, 3> *projectionPoint,
                                                                        std::array<double, 3> *projectionNormal) const
{
    // Get segment
    const Cell &segment = *segmentItr;
    assert(segment.getType() == ElementType::TRIANGLE);

    // Get segment vertices
    ConstProxyVector<long> segmentVertexIds = segment.getVertexIds();

    const std::array<double,3> &point0 = m_surface->getVertexCoords(segmentVertexIds[0]);
    const std::array<double,3> &point1 = m_surface->getVertexCoords(segmentVertexIds[1]);
    const std::array<double,3> &point2 = m_surface->getVertexCoords(segmentVertexIds[2]);

    // Get projection point on segment
    std::array<double,3> tau;
    std::array<double, 3> point_s = CGElem::projectPointTriangle(point, point0, point1, point2, &(tau[0]));
    std::array<double, 3> d0_point_s = point0 - point2;
    std::array<double, 3> d1_point_s = point1 - point2;

    // Early return if point_s consides with segment's node
    double distanceTolerance = m_surface->getTol();
    for (int i = 0; i < 3; ++i) {
        if (utils::DoubleFloatingEqual()(tau[i], 1., distanceTolerance, distanceTolerance)) {
            (*projectionPoint)  = m_surface->getVertexCoords(segmentVertexIds[i]);
            (*projectionNormal) = computeSegmentVertexNormal(segmentItr, i, true);
            return;
        }
    }

    // Get normal on segment
    std::array<double, 3> facetNormal = computeSegmentNormal(segmentItr);
    std::array<double, 3> normal_s = point - point_s;

    double distance = norm2(normal_s);
    if (utils::DoubleFloatingEqual()(distance, 0., distanceTolerance, distanceTolerance)) {
        normal_s = facetNormal;
    } else {
        normal_s /= distance;
        if (dotProduct(normal_s, facetNormal) < 0.0) {
            normal_s *= -1.0;
        }
    }

    // Compute projection point and normal at each edge "e"
    std::array<std::array<double, 3>,3 > edgePoint_s;
    std::array< std::array<double, 3>, 3 > d0_edgePoint_s;
    std::array< std::array<double, 3>, 3 > d1_edgePoint_s;
    std::array< std::array<double, 3>, 3 > edgePoint_p;
    std::array< std::array<double, 3>, 3 > d0_edgePoint_p;
    std::array< std::array<double, 3>, 3 > d1_edgePoint_p;
    std::array< std::array<double, 3>, 3 > edgeNormal_p;
    std::array< std::array<double, 3>, 3 > d0_edgeNormal_p;
    std::array< std::array<double, 3>, 3 > d1_edgeNormal_p;
    for (int e = 0; e < 3; ++e) {
        // Get segment
        const Cell &segment = *segmentItr;

        // Get edge's local and global node ids
        ConstProxyVector<int> edgeLocalVertexIds = segment.getFaceLocalVertexIds(e);
        int localId0 = edgeLocalVertexIds[0]; 
        int localId1 = edgeLocalVertexIds[1]; 
        int id0      = segment.getVertexId(localId0);
        int id1      = segment.getVertexId(localId1);

        // Get vertex information
        const std::array<double,3> &node0 = m_surface->getVertexCoords(id0);
        const std::array<double,3> &node1 = m_surface->getVertexCoords(id1);

        std::array<double, 3> normal0 = computeSegmentVertexNormal(segmentItr, localId0, true);
        std::array<double, 3> normal1 = computeSegmentVertexNormal(segmentItr, localId1, true);

        // Get edge normal
        std::array<double, 3> edgeNormal_s = computeSegmentEdgeNormal(segmentItr, e, true);

        // Compute projection coefficients
        std::array<double, 3> edge = node1 - node0;
        double p01 = dotProduct(point_s - node0, edge) / dotProduct(edge, edge);
        double d0_p01;
        double d1_p01;
        if (p01 > 1.0) {
            p01 = 1.0;
            d0_p01 = 0.0;
            d1_p01 = 0.0;
        } else if (p01 < 0.0) {
            p01 = 0.0;
            d0_p01 = 0.0;
            d1_p01 = 0.0;
        } else {
            d0_p01 = dotProduct(d0_point_s - node0, edge) / dotProduct(edge, edge);
            d1_p01 = dotProduct(d1_point_s - node0, edge) / dotProduct(edge, edge);
        }
        double p10 = 1.0 - p01;
        double d0_p10 = - d0_p01;
        double d1_p10 = - d1_p01;

        // Compute lambda coefficients
        double lambda01 = -dotProduct(edge, normal0) / dotProduct(edgeNormal_s, normal0); 
        double lambda10 =  dotProduct(edge, normal1) / dotProduct(edgeNormal_s, normal1); 

        // Edge projection point
        double product1   = p10 * p01;
        double d_product1 = p01 - p10; // derivative wrt p10

        double product2   = lambda01 * p10 + lambda10 * p01;
        double d_product2 = lambda01 - lambda10;

        double f01     = product1 * product2;
        double d_f01   = d_product1 * product2 + product1 * d_product2;
        double d_d_f01 = 2.0 * (d_product1 * d_product2 - product2);

        edgePoint_s[e] = p10 * node0 + p01 * node1;
        std::array<double, 3> d_edgePoint_s = node0 - node1;

        d0_edgePoint_s[e] = d_edgePoint_s  * d0_p10;
        d1_edgePoint_s[e] = d_edgePoint_s  * d1_p10;

        edgePoint_p[e] = edgePoint_s[e] + f01 * edgeNormal_s;
        std::array<double, 3> d_edgePoint_p = d_edgePoint_s + d_f01 * edgeNormal_s;

        d0_edgePoint_p[e] = d_edgePoint_p * d0_p10;
        d1_edgePoint_p[e] = d_edgePoint_p * d1_p10;

        // Edge surface normal
        double length   = norm2(edge);
        edgeNormal_p[e] = d_f01 * edge + length * length * edgeNormal_s;
        std::array<double, 3> d_edgeNormal_p = d_d_f01 * edge;

        double norm = norm2(edgeNormal_p[e]);
        double d_norm = dotProduct(edgeNormal_p[e], d_edgeNormal_p) / norm;
        edgeNormal_p[e] /= norm;
        d_edgeNormal_p = (d_edgeNormal_p - edgeNormal_p[e] * d_norm) / norm;

        d0_edgeNormal_p[e] = d_edgeNormal_p * d0_p10;
        d1_edgeNormal_p[e] = d_edgeNormal_p * d1_p10;
    }

    // Project surface point to the intermediator polygon
    double sum_t = 0.0;
    double d0_sum_t = 0.0;
    double d1_sum_t = 0.0;
    std::array<double, 3> t;
    std::array<double, 3> d0_t;
    std::array<double, 3> d1_t;
    for (int i = 0; i < 3; i ++) { // loop in vertices
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        std::array<double, 3> v1 = edgePoint_s[j] - point_s;
        std::array<double, 3> d0_v1 = d0_edgePoint_s[j]- d0_point_s;
        std::array<double, 3> d1_v1 = d1_edgePoint_s[j]- d1_point_s;

        std::array<double, 3> v2 = edgePoint_s[k] - point_s;
        std::array<double, 3> d0_v2 = d0_edgePoint_s[k]- d0_point_s;
        std::array<double, 3> d1_v2 = d1_edgePoint_s[k]- d1_point_s;

        std::array<double, 3> normal = crossProduct(v1, v2);
        std::array<double, 3> d0_normal = crossProduct(d0_v1, v2) + crossProduct(v1, d0_v2);
        std::array<double, 3> d1_normal = crossProduct(d1_v1, v2) + crossProduct(v1, d1_v2);

        t[i] = norm2(normal);
        d0_t[i] = dotProduct(normal, d0_normal) / t[i];
        d1_t[i] = dotProduct(normal, d1_normal) / t[i];
        sum_t += t[i];
        d0_sum_t += d0_t[i];
        d1_sum_t += d1_t[i];
    }

    std::array<double, 3> point_m    = {0.0, 0.0, 0.0};
    std::array<double, 3> d0_point_m = {0.0, 0.0, 0.0};
    std::array<double, 3> d1_point_m = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i ++) {
        t[i] /= sum_t;
        d0_t[i] = (d0_t[i] - t[i] * d0_sum_t) / sum_t;
        d1_t[i] = (d1_t[i] - t[i] * d1_sum_t) / sum_t;
        point_m += t[i] * edgePoint_p[i];
        d0_point_m += d0_t[i] * edgePoint_p[i] + t[i] * d0_edgePoint_p[i];
        d1_point_m += d1_t[i] * edgePoint_p[i] + t[i] * d1_edgePoint_p[i];
    }

    // Early return if point_s coinsides with a triangle
    // edge (or a node of the intermediate triangle)
    for (int i = 0; i < 3; ++i) {
        if (utils::DoubleFloatingEqual()(t[i], 1., distanceTolerance, distanceTolerance)) {
            (*projectionPoint)  = edgePoint_p[i];
            (*projectionNormal) = edgeNormal_p[i];
            return;
        }
    }

    // Compute 2D f function
    double f    = 0.0;
    double d0_f = 0.0;
    double d1_f = 0.0;
    for (int e = 0; e < 3; ++e) {
        int i = e;           // previous vertex
        int j = (e + 1) % 3; // next vertex

        // Get vertex information
        std::array<double,3> &node_i    = edgePoint_p[i];
        std::array<double,3> &d0_node_i = d0_edgePoint_p[i];
        std::array<double,3> &d1_node_i = d1_edgePoint_p[i];

        std::array<double,3> &node_j     = edgePoint_p[j];
        std::array<double, 3> &d0_node_j = d0_edgeNormal_p[j];
        std::array<double, 3> &d1_node_j = d1_edgeNormal_p[j];

        std::array<double, 3> &normal_i    = edgeNormal_p[i];
        std::array<double, 3> &d0_normal_i = d0_edgeNormal_p[i];
        std::array<double, 3> &d1_normal_i = d1_edgeNormal_p[i];

        std::array<double, 3> &normal_j    = edgeNormal_p[j];
        std::array<double, 3> &d0_normal_j = d0_edgeNormal_p[j];
        std::array<double, 3> &d1_normal_j = d1_edgeNormal_p[j];

        // Compute lambda parameters
        std::array<double, 3> edge    = node_j - node_i;
        std::array<double, 3> d0_edge = d0_node_j - d0_node_i;
        std::array<double, 3> d1_edge = d1_node_j - d1_node_i;

        double lambda_ij    = -dotProduct(edge, normal_i) / dotProduct(normal_s, normal_i);
        double d0_lambda_ij = -(dotProduct(edge - lambda_ij * normal_s, d0_normal_i) + dotProduct(d0_edge, normal_i)) / dotProduct(normal_s, normal_i);
        double d1_lambda_ij = -(dotProduct(edge - lambda_ij * normal_s, d1_normal_i) + dotProduct(d1_edge, normal_i)) / dotProduct(normal_s, normal_i);

        double lambda_ji    = dotProduct(edge, normal_j) / dotProduct(normal_s, normal_j);
        double d0_lambda_ji = (dotProduct(edge - lambda_ji * normal_s, d0_normal_j) + dotProduct(d0_edge, normal_j)) / dotProduct(normal_s, normal_j);
        double d1_lambda_ji = (dotProduct(edge - lambda_ji * normal_s, d1_normal_j) + dotProduct(d1_edge, normal_j)) / dotProduct(normal_s, normal_j);

        // Compute contribution to the f function
        f += t[i] * t[j] * (lambda_ij * t[i] + lambda_ji * t[j]);
        d0_f += (d0_t[i] * t[j] + t[i] * d0_t[j]) * (lambda_ij * t[i] + lambda_ji * t[j])
              + t[i] * t[j] * (d0_lambda_ij * t[i] + lambda_ij * d0_t[i] + d0_lambda_ji * t[j] + lambda_ji * d0_t[j]);
        d1_f += (d1_t[i] * t[j] + t[i] * d1_t[j]) * (lambda_ij * t[i] + lambda_ji * t[j])
              + t[i] * t[j] * (d1_lambda_ij * t[i] + lambda_ij * d1_t[i] + d1_lambda_ji * t[j] + lambda_ji * d1_t[j]);
    }

    // Compute projection point on curved surface
    (*projectionPoint) = point_m + f * normal_s;

    // Compute projection normal on curved surface
    std::array<double, 3> d0_point_p = d0_point_m + d0_f * normal_s;
    std::array<double, 3> d1_point_p = d1_point_m + d1_f * normal_s;

    (*projectionNormal) = crossProduct(d0_point_p, d1_point_p);
    (*projectionNormal) /= norm2((*projectionNormal));

    if (dotProduct(normal_s, (*projectionNormal)) < 0.0) {
       (*projectionNormal) *= -1.0;
    }
}

/*!
 * Evaluate the projection of the given point on a smooth surface created 
 * over a discretized geometry. The surface passes from the geometry's points
 * and it's vertical to the normal vectors associated with them. Also, it's
 * continuous and offers normal vector continuity over the different segments
 * of the discretized geometry.
 *
 * In this case the segment is a non degenerated polygon. No checks are performed to
 * find out if the specified polygon is degenerate.
 *
 * The developed algorithm does not exist in literature and it's a straightforward
 * expansion of the algorithm developed for triangular segments. The algorithm is
 * described in the documentation of function "evalHighOrderProjectionOnTriangle".
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] segmentItr is an iterator pointing to the segment on which the surface
 * will be created.
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 * projection point. It is computed approximately.
 */
void LevelSetSegmentationSurfaceInfo::evalHighOrderProjectionOnPolygon(const std::array<double,3> &point,
                                                                       const SegmentConstIterator &segmentItr,
                                                                       std::array<double, 3> *projectionPoint,
                                                                       std::array<double, 3> *projectionNormal) const
{
    // Get segment
    const Cell &segment = *segmentItr;

    // Get projection point on segment
    ConstProxyVector<long> segmentVertexIds = m_surface->getFacetOrderedVertexIds(segment);

    std::size_t nSegmentVertices = segmentVertexIds.size();
    BITPIT_CREATE_WORKSPACE(segmentVertexCoords, std::array<double BITPIT_COMMA 3>, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    m_surface->getVertexCoords(segmentVertexIds.size(), segmentVertexIds.data(), segmentVertexCoords);

    BITPIT_CREATE_WORKSPACE(tau, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    std::array<double, 3> point_s = CGElem::projectPointPolygon(point, nSegmentVertices, segmentVertexCoords, tau);

    std::array<double, 3> lowOrderProjectionNormal = tau[0] * computeSegmentVertexNormal(segmentItr, 0, true);
    for (std::size_t i = 1; i < nSegmentVertices; ++i) {
        lowOrderProjectionNormal += tau[i] * computeSegmentVertexNormal(segmentItr, i, true);
    }

    // Early return if point_s consides with segment's node
    double distanceTolerance = m_surface->getTol();
    for (std::size_t i = 0; i < nSegmentVertices; ++i) {
        if (utils::DoubleFloatingEqual()(tau[i], 1., distanceTolerance, distanceTolerance)) {
            (*projectionPoint)  = segmentVertexCoords[i];
            (*projectionNormal) = computeSegmentVertexNormal(segmentItr, i, true);
            return;
        }
    }

    std::array<double, 3> d0_point_s;
    std::array<double, 3> d1_point_s;

    // Define parameterization of point_s
    // Definition: point_s - x_central = t0 * (x_prev - x_central) + t1 * (x_next - x_central)
    // where x_prev, x_central, x_next are subsequent polygon vertices
    // Vertex x_central is the one corresponding to polygon angle closer to 90 degrees
    double cosMin = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < nSegmentVertices; ++i) {
        int i_p = i;                          // previous
        int i_c = (i + 1) % nSegmentVertices; // central
        int i_n = (i + 2) % nSegmentVertices; // next

        const std::array<double,3> &point_p = m_surface->getVertexCoords(segmentVertexIds[i_p]);
        const std::array<double,3> &point_c = m_surface->getVertexCoords(segmentVertexIds[i_c]);
        const std::array<double,3> &point_n = m_surface->getVertexCoords(segmentVertexIds[i_n]);

        std::array<double,3> x_p = point_p - point_c;
        std::array<double,3> x_n = point_n - point_c;

        double cos = std::abs(dotProduct(x_p, x_n)) / (norm2(x_p) * norm2(x_n));
        if (cos < cosMin) {
            cosMin = cos;
            d0_point_s = x_p;
            d1_point_s = x_n;
        }
    }

    // Get normal on segment
    std::array<double, 3> facetNormal = computeSegmentNormal(segmentItr);
    std::array<double, 3> normal_s = point - point_s;
    double distance = norm2(normal_s);
    if (utils::DoubleFloatingEqual()(distance, 0., distanceTolerance, distanceTolerance)) {
        normal_s = facetNormal;
    } else {
        normal_s /= distance;
        if (dotProduct(normal_s, facetNormal) < 0.0) {
            normal_s *= -1.0;
        }
    }

    // Compute projection point and normal at each edge "e"
    int nSegmentEdges = nSegmentVertices;
    BITPIT_CREATE_WORKSPACE(edgePoint_s, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(d0_edgePoint_s, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(d1_edgePoint_s, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(edgePoint_p, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(d0_edgePoint_p, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(d1_edgePoint_p, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(edgeNormal_p, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(d0_edgeNormal_p, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    BITPIT_CREATE_WORKSPACE(d1_edgeNormal_p, std::array<double BITPIT_COMMA 3>, nSegmentEdges, ReferenceElementInfo::MAX_ELEM_EDGES);
    for (int e = 0; e < nSegmentEdges; ++e) {
        // Get edge's local and global node ids
        ConstProxyVector<int> edgeLocalVertexIds = segment.getFaceLocalVertexIds(e);
        int localId0 = edgeLocalVertexIds[0]; 
        int localId1 = edgeLocalVertexIds[1]; 
        int id0      = segment.getVertexId(localId0);
        int id1      = segment.getVertexId(localId1);

        // Get vertex information
        const std::array<double,3> &node0 = m_surface->getVertexCoords(id0);
        const std::array<double,3> &node1 = m_surface->getVertexCoords(id1);

        std::array<double, 3> normal0 = computeSegmentVertexNormal(segmentItr, localId0, true);
        std::array<double, 3> normal1 = computeSegmentVertexNormal(segmentItr, localId1, true);

        // Get edge normal
        std::array<double, 3> edgeNormal_s = computeSegmentEdgeNormal(segmentItr, e, true);

        // Compute projection coefficients
        std::array<double, 3> edge = node1 - node0;
        double p01 = dotProduct(point_s - node0, edge) / dotProduct(edge, edge);
        double d0_p01;
        double d1_p01;
        if (p01 > 1.0) {
            p01 = 1.0;
            d0_p01 = 0.0;
            d1_p01 = 0.0;
        } else if (p01 < 0.0) {
            p01 = 0.0;
            d0_p01 = 0.0;
            d1_p01 = 0.0;
        } else {
            d0_p01 = dotProduct(d0_point_s - node0, edge) / dotProduct(edge, edge);
            d1_p01 = dotProduct(d1_point_s - node0, edge) / dotProduct(edge, edge);
        }
        double p10 = 1.0 - p01;
        double d0_p10 = - d0_p01;
        double d1_p10 = - d1_p01;

        // Compute lambda coefficients
        double lambda01 = -dotProduct(edge, normal0) / dotProduct(edgeNormal_s, normal0); 
        double lambda10 =  dotProduct(edge, normal1) / dotProduct(edgeNormal_s, normal1); 

        // Edge projection point
        double product1 = p10 * p01;
        double d_product1 = p01 - p10; // derivative wrt p10

        double product2 = lambda01 * p10 + lambda10 * p01;
        double d_product2 = lambda01 - lambda10;

        double f01     = product1 * product2;
        double d_f01   = d_product1 * product2 + product1 * d_product2;
        double d_d_f01 = 2.0 * (d_product1 * d_product2 - product2);

        edgePoint_s[e] = p10 * node0 + p01 * node1;
        std::array<double, 3> d_edgePoint_s = node0 - node1;

        d0_edgePoint_s[e] = d_edgePoint_s  * d0_p10;
        d1_edgePoint_s[e] = d_edgePoint_s  * d1_p10;

        edgePoint_p[e] = edgePoint_s[e] + f01 * edgeNormal_s;
        std::array<double, 3> d_edgePoint_p = d_edgePoint_s + d_f01 * edgeNormal_s;

        d0_edgePoint_p[e] = d_edgePoint_p * d0_p10;
        d1_edgePoint_p[e] = d_edgePoint_p * d1_p10;

        // Edge surface normal
        double length   = norm2(edge);
        edgeNormal_p[e] = d_f01 * edge + length * length * edgeNormal_s;
        std::array<double, 3> d_edgeNormal_p = d_d_f01 * edge;

        double norm = norm2(edgeNormal_p[e]);
        double d_norm = dotProduct(edgeNormal_p[e], d_edgeNormal_p) / norm;
        edgeNormal_p[e] /= norm;
        d_edgeNormal_p = (d_edgeNormal_p - edgeNormal_p[e] * d_norm) / norm;

        d0_edgeNormal_p[e] = d_edgeNormal_p * d0_p10;
        d1_edgeNormal_p[e] = d_edgeNormal_p * d1_p10;
    }

    // Early return if point_s coinsides with a polygon's
    // edge (or a node of the intermediate polygon)
    for (int e = 0; e < nSegmentEdges; ++e) {
        double distance = norm2(edgePoint_s[e] - point_s);
        if (utils::DoubleFloatingEqual()(distance, 0.0, distanceTolerance, distanceTolerance)) {
            (*projectionPoint)  = edgePoint_p[e];
            (*projectionNormal) = edgeNormal_p[e];
            return;
        }
    }

    // Project surface point to the intermediator polygon
    // The projection is done by following formula:
    // point_m - point_s = sum i=0,vertices-1 { t[i] * (edgePoint_p[i] - edgePoint_s[i]) }
    // After various numerical experiments it is concluded that a quite smooth interpolation is
    // resulted from the following choice of the interpolation weights "t":
    // t[i] = w[i] / sum j=0,vertices-1 { w[j] }, where
    // w[i] = product j=0,vertices-1, j!=i { d[j]^4 }, where
    // d[i] is the distance between point_s and edgePoint_s[i]
    BITPIT_CREATE_WORKSPACE(d, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    BITPIT_CREATE_WORKSPACE(d0_d, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    BITPIT_CREATE_WORKSPACE(d1_d, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    for (std::size_t i = 0; i < nSegmentVertices; ++i) {
        d[i] = (edgePoint_s[i][0] - point_s[0]) * (edgePoint_s[i][0] - point_s[0])
        + (edgePoint_s[i][1] - point_s[1]) * (edgePoint_s[i][1] - point_s[1])
        + (edgePoint_s[i][2] - point_s[2]) * (edgePoint_s[i][2] - point_s[2]);

        d[i] = d[i] * d[i];

        d0_d[i] = 4.0 * d[i] * dotProduct((edgePoint_s[i] - point_s), (d0_edgePoint_s[i] - d0_point_s));
        d1_d[i] = 4.0 * d[i] * dotProduct((edgePoint_s[i] - point_s), (d1_edgePoint_s[i] - d1_point_s));
    }

    BITPIT_CREATE_WORKSPACE(t, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    BITPIT_CREATE_WORKSPACE(d0_t, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    BITPIT_CREATE_WORKSPACE(d1_t, double, nSegmentVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    double sum = 0.0;
    double d0_sum = 0.0;
    double d1_sum = 0.0;
    for (std::size_t i = 0; i < nSegmentVertices; ++i) {
        t[i]    = 1.0;
        d0_t[i] = 0.0;
        d1_t[i] = 0.0;
        for (std::size_t j = 0; j < nSegmentVertices; ++j) {
            if (i != j) {
                t[i] *= d[j];

                double product = 1.0;
                for (std::size_t k = 0; k < nSegmentVertices; ++k) {
                    if (j != k) {
                        product *= d[k];
                    }
                }
                d0_t[i] += d0_d[j] * product;
                d1_t[i] += d1_d[j] * product;
            }
        }
        sum    += t[i];
        d0_sum += d0_t[i];
        d1_sum += d1_t[i];
    }
    for (std::size_t i = 0; i < nSegmentVertices; ++i) {
        t[i] /= sum;
        d0_t[i] = (d0_t[i] - t[i] * d0_sum) / sum;
        d1_t[i] = (d1_t[i] - t[i] * d1_sum) / sum;
    }

    std::array<double, 3> point_m = {0.0, 0.0, 0.0};
    std::array<double, 3> d0_point_m = {0.0, 0.0, 0.0};
    std::array<double, 3> d1_point_m = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < nSegmentVertices; ++i) {
        point_m += (edgePoint_p[i] - edgePoint_s[i]) * t[i];
        d0_point_m += d0_t[i] * (edgePoint_p[i]- edgePoint_s[i]) + t[i] * (d0_edgePoint_p[i] - d0_edgePoint_s[i]);
        d1_point_m += d1_t[i] * (edgePoint_p[i]- edgePoint_s[i]) + t[i] * (d1_edgePoint_p[i] - d1_edgePoint_s[i]);
    }
    point_m += point_s;
    d0_point_m += d0_point_s;
    d1_point_m += d1_point_s;

    // Compute 2D f function
    double f    = 0.0;
    double d0_f = 0.0;
    double d1_f = 0.0;
    for (int e = 0; e < nSegmentEdges; ++e) {
        int i = e;                       // previous vertex
        int j = (e + 1) % nSegmentEdges; // next vertex

        // Get vertex information
        std::array<double,3> &node_i    = edgePoint_p[i];
        std::array<double,3> &d0_node_i = d0_edgePoint_p[i];
        std::array<double,3> &d1_node_i = d1_edgePoint_p[i];

        std::array<double, 3> &node_j     = edgePoint_p[j];
        std::array<double, 3> &d0_node_j = d0_edgeNormal_p[j];
        std::array<double, 3> &d1_node_j = d1_edgeNormal_p[j];

        std::array<double, 3> &normal_i    = edgeNormal_p[i];
        std::array<double, 3> &d0_normal_i = d0_edgeNormal_p[i];
        std::array<double, 3> &d1_normal_i = d1_edgeNormal_p[i];

        std::array<double, 3> &normal_j    = edgeNormal_p[j];
        std::array<double, 3> &d0_normal_j = d0_edgeNormal_p[j];
        std::array<double, 3> &d1_normal_j = d1_edgeNormal_p[j];

        // Compute lambda parameters
        std::array<double, 3> edge    = node_j - node_i;
        std::array<double, 3> d0_edge = d0_node_j - d0_node_i;
        std::array<double, 3> d1_edge = d1_node_j - d1_node_i;

        double lambda_ij    = -dotProduct(edge, normal_i) / dotProduct(normal_s, normal_i);
        double d0_lambda_ij = -(dotProduct(edge - lambda_ij * normal_s, d0_normal_i) + dotProduct(d0_edge, normal_i)) / dotProduct(normal_s, normal_i);
        double d1_lambda_ij = -(dotProduct(edge - lambda_ij * normal_s, d1_normal_i) + dotProduct(d1_edge, normal_i)) / dotProduct(normal_s, normal_i);

        double lambda_ji    =  dotProduct(edge, normal_j) / dotProduct(normal_s, normal_j);
        double d0_lambda_ji = (dotProduct(edge - lambda_ji * normal_s, d0_normal_j) + dotProduct(d0_edge, normal_j)) / dotProduct(normal_s, normal_j);
        double d1_lambda_ji = (dotProduct(edge - lambda_ji * normal_s, d1_normal_j) + dotProduct(d1_edge, normal_j)) / dotProduct(normal_s, normal_j);

        // Compute contribution to the f function
        f += t[i] * t[j] * (lambda_ij * t[i] + lambda_ji * t[j]);
        d0_f += (d0_t[i] * t[j] + t[i] * d0_t[j]) * (lambda_ij * t[i] + lambda_ji * t[j])
              + t[i] * t[j] * (d0_lambda_ij * t[i] + lambda_ij * d0_t[i] + d0_lambda_ji * t[j] + lambda_ji * d0_t[j]);
        d1_f += (d1_t[i] * t[j] + t[i] * d1_t[j]) * (lambda_ij * t[i] + lambda_ji * t[j])
              + t[i] * t[j] * (d1_lambda_ij * t[i] + lambda_ij * d1_t[i] + d1_lambda_ji * t[j] + lambda_ji * d1_t[j]);
    }

    // Compute projection point on curved surface
    (*projectionPoint) = point_m + f * normal_s;

    // Compute projection normal on curved surface
    std::array<double, 3> d0_point_p = d0_point_m + d0_f * normal_s;
    std::array<double, 3> d1_point_p = d1_point_m + d1_f * normal_s;

    (*projectionNormal) = crossProduct(d0_point_p, d1_point_p);
    (*projectionNormal) /= norm2((*projectionNormal));

    if (dotProduct(normal_s, (*projectionNormal)) < 0.0) {
       (*projectionNormal) *= -1.0;
    }
}

/*!
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is vertical to the normal vectors associated with them.
 *
 * The computed projection point (\vec{x_p}) is not the closest to the given
 * point, but it's defined from the following equation:
 * \vec{x_p} = \vec{x_s} + f * \vec{n_s},
 * where \vec{x_s} and \vec{n_s} is the projection point and boundary normal on the
 * discretized polygonal/plyhedral geometry and f is a function depending on \vec{x_s}
 * which returns a real number.
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] segmentItr is an iterator pointing to the segment on which the surface
 * will be created.
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 * projection point.
 */
void LevelSetSegmentationSurfaceInfo::evalHighOrderProjection(const std::array<double,3> &point,
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
        evalHighOrderProjectionOnLine(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    case ElementType::TRIANGLE:
    {
        evalHighOrderProjectionOnTriangle(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    default:
    {
        evalHighOrderProjectionOnPolygon(point, segmentItr, projectionPoint, projectionNormal);
        return;
    }

    }
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
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the normal to the surface vector on the surface
 */
void LevelSetSegmentationSurfaceInfo::evalProjection(const std::array<double, 3> &point,
                                                     const SegmentConstIterator &segmentItr,
                                                     std::array<double, 3> *projectionPoint,
                                                     std::array<double, 3> *projectionNormal) const
{
    if (m_surfaceSmoothing == LevelSetSurfaceSmoothing::HIGH_ORDER) {
        evalHighOrderProjection(point, segmentItr, projectionPoint, projectionNormal);
    } else {
        evalLowOrderProjection(point, segmentItr, projectionPoint, projectionNormal);
    }
}

/*!
 * Evaluate the projection of the given point on the specified segment.
 *
 * @param[in] point are the coordinates of point
 * @param[in] segmentItr is an iterator pointing to the closest segment
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 */
void LevelSetSegmentationSurfaceInfo::evalProjection(const std::array<double, 3> &point,
                                                     const SegmentConstIterator &segmentItr,
                                                     std::array<double, 3> *projectionPoint) const
{
    std::array<double, 3> projectionNormal;
    evalProjection(point, segmentItr, projectionPoint, &projectionNormal);

    BITPIT_UNUSED(projectionNormal);
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
    evalProjection(point, segmentItr, &projectionPoint);

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
    auto evaluator = [this] (long id) -> std::array<double,3>
        {
            return _evalCellNormal(id, false);
        };

    auto fallback = [] (long id) -> const std::array<double,3> &
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
    // Project the point on the surface and evaluate the point-projection vector
    std::array<double, 3> projectionPoint;
    std::array<double, 3> projectionNormal;
    evalProjection(point, signedLevelSet, &projectionPoint, &projectionNormal);

    return projectionNormal;
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
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is verical to the normal vectors associated with them.
 *
 * If the usage of an  unsigned level set is requested, the orientation of the
 * surface should be discarded in order to have a normal that is agnostic with
 * respect to the two sides of the surface. If the sign is not cached, it
 * should be evaluated from scratch and this will cause the projection
 * evaluation to be performed twice.
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the norrmal to the surface vector on the surface
 * projection point.
 */
void LevelSetSegmentationBaseObject::evalProjection(const std::array<double,3> &point,
                                                    bool signedLevelSet,
                                                    std::array<double, 3> *projectionPoint,
                                                    std::array<double, 3> *projectionNormal) const
{
    _evalProjection(point, signedLevelSet, projectionPoint, projectionNormal);
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
        auto evaluator = [this] (long id) -> long { return evalCellSupport(id); };
        auto fallback = [] (long id) -> long { BITPIT_UNUSED(id); return levelSetDefaults::SUPPORT; };
        flushVTKOutputData<double>(stream, format, field, evaluator, fallback);
        break;
    }

    case LevelSetField::PART:
    {
        auto evaluator = [this] (long id) -> int { return evalCellPart(id); };
        auto fallback = [] (long id) -> int { BITPIT_UNUSED(id); return levelSetDefaults::PART; };
        flushVTKOutputData<double>(stream, format, field, evaluator, fallback);
        break;
    }

    case LevelSetField::NORMAL:
    {
        auto evaluator = [this] (long id) -> std::array<double,3> { return evalCellNormal(id, true); };
        auto fallback = [] (long id) -> const std::array<double,3> & { BITPIT_UNUSED(id); return levelSetDefaults::NORMAL; };
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
    // First we need to check if the cell intersects the surface, and only if it
    // doesn't we  should check if its distance is lower than the narrow band size.
    //
    // When high order smoothing is active, there may be cells whose support is within the search
    // radius, but they are not intersected and their distance is less than the narrow band size.
    // These cells are not geometrically inside the narrow band, they are neighbours of cells
    // geometrically inside the narrow band and as such it's up to the caller of this function to
    // identify their cell location.
    LevelSetCellLocation cellLocation = LevelSetCellLocation::UNKNOWN;
    if (_intersectSurface(id, cellUnsigendValue, CELL_LOCATION_INTERSECTION_MODE) == LevelSetIntersectionStatus::TRUE) {
        cellLocation = LevelSetCellLocation::NARROW_BAND_INTERSECTED;
    } else if (cellUnsigendValue <= m_narrowBandSize) {
        cellLocation = LevelSetCellLocation::NARROW_BAND_DISTANCE;
    }
    assert((getSurfaceSmoothing() == LevelSetSurfaceSmoothing::HIGH_ORDER) || (cellLocation != LevelSetCellLocation::UNKNOWN));

    if (cellLocation != LevelSetCellLocation::UNKNOWN) {
        CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);
        locationCache->insertEntry(id, static_cast<char>(cellLocation));
    }

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

    std::array<double, 3> projectionPoint;
    std::array<double, 3> projectionNormal;
    _evalProjection(centroid, support, signedLevelSet, &projectionPoint, &projectionNormal);

    return projectionNormal;
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
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is verical to the normal vectors associated with them.
 *
 * If the usage of an  unsigned level set is requested, the orientation of the
 * surface should be discarded in order to have a normal that is agnostic with
 * respect to the two sides of the surface. If the sign is not cached, it
 * should be evaluated from scratch and this will cause the projection
 * evaluation to be performed twice.
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the norrmal to the surface vector on the surface
 * projection point.
 */
void LevelSetSegmentationObject::_evalProjection(const std::array<double,3> &point,
                                                 bool signedLevelSet,
                                                 std::array<double, 3> *projectionPoint,
                                                 std::array<double, 3> *projectionNormal) const
{
    // Get closest segment
    long support = evalSupport(point);

    _evalProjection(point, support, signedLevelSet, projectionPoint, projectionNormal);
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
    std::array<double,3> distanceVector;
    double distance = m_surfaceInfo->evalDistance(point, supportItr, signedLevelSet, &distanceVector);

    // Early return if the point lies on the surface
    if (evalValueSign(distance) == 0) {
        if (signedLevelSet) {
            return m_surfaceInfo->evalNormal(point, supportItr);
        } else {
            return {{0., 0., 0.}};
        }
    }

    // Evaluate levelset gradient
    std::array<double,3> gradient = distanceVector / distance;

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
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is verical to the normal vectors associated with them.
 *
 * If the usage of an  unsigned level set is requested, the orientation of the
 * surface should be discarded in order to have a normal that is agnostic with
 * respect to the two sides of the surface. If the sign is not cached, it
 * should be evaluated from scratch and this will cause the projection
 * evaluation to be performed twice.
 *
 * \param[in] point are the coordinates of the given point
 * \param[in] support is the the closest segment to the specified point
 * \param[in] signedLevelSet controls if signed levelset function will be used
 * \param[out] projectionPoint if a valid pointer is provided, contains the coordinates
 * of the projection point on the surface.
 * \param[out] projectionNormal if a valid pointer is provided, contains the normal to
 * the surface on the projection point. If the projection point lies on the zero-level-set surface,
 * a zero vector is returned
 */
void LevelSetSegmentationObject::_evalProjection(const std::array<double,3> &point,
                                                 long support,
                                                 bool signedLevelSet,
                                                 std::array<double, 3> *projectionPoint,
                                                 std::array<double, 3> *projectionNormal) const
{
    // Early return if projection normal and projection point pointers are null
    if (!(projectionNormal) && !(projectionPoint)) {
       return;
    }

    // Early return if the support is not valid
    //
    // With an invalid support, only the unsigend levelset can be evaluated.
    if (support < 0) {
        if (!signedLevelSet || empty()) {
            (*projectionPoint) = levelSetDefaults::POINT;
            (*projectionNormal) = levelSetDefaults::NORMAL;
        }

        throw std::runtime_error("With an invalid support, only the unsigend levelset can be evaluated.");
    }

    // Get closest segment
    LevelSetSegmentationSurfaceInfo::SegmentConstIterator segmentItr = getSurface().getCellConstIterator(support);

    // Eval projection point and normal
    m_surfaceInfo->evalProjection(point, segmentItr, projectionPoint, projectionNormal);

    // Early return if a projection normal null pointer is given
    if (!(projectionNormal)) {
       return;
    }

    // If an unsigned evaluation is requested, the orientation of the surface should be discarded
    // in order to have a normal that is agnostic with respect the two sides of the surface.
    if (!signedLevelSet) {
        short sign = evalSign(point);
        (*projectionNormal) *= static_cast<double>(sign);
    }
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
    return _evalCellFunction<std::array<double,3>>(id, signedLevelSet, [&id, signedLevelSet] (const LevelSetBooleanResult<LevelSetSegmentationBaseObject> &result) -> std::array<double,3>
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
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is verical to the normal vectors associated with them.
 *
 * If the usage of an  unsigned level set is requested, the orientation of the
 * surface should be discarded in order to have a normal that is agnostic with
 * respect to the two sides of the surface. If the sign is not cached, it
 * should be evaluated from scratch and this will cause the projection
 * evaluation to be performed twice.
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the norrmal to the surface vector on the surface
 * projection point.
 */
void LevelSetBooleanObject<LevelSetSegmentationBaseObject>::_evalProjection(const std::array<double,3> &point,
                                                                            bool signedLevelSet,
                                                                            std::array<double, 3> *projectionPoint,
                                                                            std::array<double, 3> *projectionNormal) const
{
    return _evalFunction<void>(point, signedLevelSet, [&point, projectionPoint, projectionNormal, signedLevelSet] (const LevelSetBooleanResult<LevelSetSegmentationBaseObject> &result)
        {
            const LevelSetSegmentationBaseObject *resultObject = result.getObject();
            if ( !resultObject ) {
                (*projectionNormal) = levelSetDefaults::NORMAL;
                (*projectionPoint)  = levelSetDefaults::POINT;
            }

            resultObject->evalProjection(point, signedLevelSet, projectionPoint, projectionNormal);
            if (signedLevelSet) {
                (*projectionNormal) *= static_cast<double>(result.getObjectSign());
            }
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
 * Evaluate the projection of the given point on the surface created based on
 * the points representing the specified segment. The surface passes from these
 * points and is verical to the normal vectors associated with them.
 *
 * If the usage of an  unsigned level set is requested, the orientation of the
 * surface should be discarded in order to have a normal that is agnostic with
 * respect to the two sides of the surface. If the sign is not cached, it
 * should be evaluated from scratch and this will cause the projection
 * evaluation to be performed twice.
 *
 * @param[in] point are the coordinates of the given point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @param[out] projectionPoint The coordinates of the projection point on the surface.
 * @param[out] projectionNormal The coordinates of the norrmal to the surface vector on the surface
 * projection point.
 */
void LevelSetComplementObject<LevelSetSegmentationBaseObject>::_evalProjection(const std::array<double,3> &point,
                                                                               bool signedLevelSet,
                                                                               std::array<double, 3> *projectionPoint,
                                                                               std::array<double, 3> *projectionNormal) const
{
    getSourceObject()->evalProjection(point, signedLevelSet, projectionPoint, projectionNormal);
    if (signedLevelSet) {
        (*projectionNormal) *= -1.;
    }
}

}
