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

# define __BITPIT_LEVELSET_SEGMENTATION_OBJECT_SRC__

# include "levelSetSegmentationObject.hpp"

namespace bitpit {

// Explicit instantization
template class LevelSetSegmentationNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>;
template class LevelSetSegmentationNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>;
template class LevelSetSegmentationNarrowBandCacheBase<LevelSetDirectStorageManager>;

template class LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
template class LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>>;
template class LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>>;

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
SegmentationKernel::SegmentationKernel( std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle ) : m_ownedSurface(std::move(surface)) {

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
    SurfUnstructured::CellConstIterator endItr = m_surface->cellConstEnd();
    for( SurfUnstructured::CellConstIterator segmentItr = m_surface->cellConstBegin(); segmentItr != endItr; ++segmentItr ){
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
 * Constructor
 *
 * \param kernel is the container associated with the storage manager
 */
LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>::LevelSetSegmentationNarrowBandCache(Kernel *kernel)
    : LevelSetExternalPiercedStorageManager(kernel), LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>(kernel), LevelSetSegmentationNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>()
{
    // Cache entries are added and processed one at the time, there is no
    // advantage in using a journaled synchronization.
    m_supportIds     = this->template addStorage<long>(this->getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_CONCURRENT);
    m_surfaceNormals = this->template addStorage<std::array<double, 3>>(this->getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_CONCURRENT);
}

/*!
 * Get a reference to the support id of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The support id of the specified entry.
 */
long & LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>::getSupportId(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_supportIds->rawAt(rawId);
}

/*!
 * Get the support id of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The support id of the specified entry.
 */
long LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>::getSupportId(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_supportIds->rawAt(rawId);
}

/*!
 * Get a reference to the surface normal of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The surface normal of the specified entry.
 */
std::array<double, 3> & LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>::getSurfaceNormal(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_surfaceNormals->rawAt(rawId);
}

/*!
 * Get the surface normal of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The surface normal of the specified entry.
 */
const std::array<double, 3> & LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>::getSurfaceNormal(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_surfaceNormals->rawAt(rawId);
}

/*!
 * Constructor
 */
LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>::LevelSetSegmentationNarrowBandCache()
    : LevelSetInternalPiercedStorageManager(), LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>(), LevelSetSegmentationNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>()
{
    // Cache entries are added and processed one at the time, there is no
    // advantage in using a journaled synchronization.
    m_supportIds     = this->template addStorage<long>(this->getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_CONCURRENT);
    m_surfaceNormals = this->template addStorage<std::array<double, 3>>(this->getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_CONCURRENT);
}

/*!
 * Get a reference to the support id of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The support id of the specified entry.
 */
long & LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>::getSupportId(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_supportIds->rawAt(rawId);
}

/*!
 * Get the support id of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The support id of the specified entry.
 */
long LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>::getSupportId(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_supportIds->rawAt(rawId);
}

/*!
 * Get a reference to the surface normal of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The surface normal of the specified entry.
 */
std::array<double, 3> & LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>::getSurfaceNormal(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_surfaceNormals->rawAt(rawId);
}

/*!
 * Get the surface normal of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The surface normal of the specified entry.
 */
const std::array<double, 3> & LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>::getSurfaceNormal(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_surfaceNormals->rawAt(rawId);
}

/*!
 * Constructor
 *
 * \param nItems are the maximum number of items the cache will hold
 */
LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>::LevelSetSegmentationNarrowBandCache(std::size_t nItems)
    : LevelSetDirectStorageManager(nItems), LevelSetNarrowBandCache<LevelSetDirectStorageManager>(nItems), LevelSetSegmentationNarrowBandCacheBase<LevelSetDirectStorageManager>()
{
    m_supportIds     = this->template addStorage<long>(this->getStorageCount());
    m_surfaceNormals = this->template addStorage<std::array<double, 3>>(this->getStorageCount());
}

/*!
 * Get a reference to the support id of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The support id of the specified entry.
 */
long & LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>::getSupportId(const KernelIterator &itr)
{
    return (*m_supportIds)[itr];
}

/*!
 * Get the support id of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The support id of the specified entry.
 */
long LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>::getSupportId(const KernelIterator &itr) const
{
    return (*m_supportIds)[itr];
}

/*!
 * Get a reference to the surface normal of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The surface normal of the specified entry.
 */
std::array<double, 3> & LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>::getSurfaceNormal(const KernelIterator &itr)
{
    return (*m_surfaceNormals)[itr];
}

/*!
 * Get the surface normal of the specified entry.
 *
 * \param itr is an iterator pointing to the entry
 * \result The surface normal of the specified entry.
 */
const std::array<double, 3> & LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>::getSurfaceNormal(const KernelIterator &itr) const
{
    return (*m_surfaceNormals)[itr];
}

/*!
 * Create the narrow band cache.
 *
 * @param ojbect is the levelset object for which the ache will be created
 */
std::shared_ptr<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>> LevelSetNarrowBandCacheFactory<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>>::create(LevelSetCachedObjectInterface<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>> *object)
{
    VolumeKernel *mesh = object->getKernel()->getMesh();
    PiercedVector<Cell, long> &cells = mesh->getCells();

    return std::shared_ptr<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>>(new LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>(&cells));
}

/*!
 * Create the narrow band cache.
 *
 * @param ojbect is the levelset object for which the ache will be created
 */
std::shared_ptr<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>> LevelSetNarrowBandCacheFactory<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>>::create(LevelSetCachedObjectInterface<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>> *object)
{
    const VolumeKernel *mesh = object->getKernel()->getMesh();
    const std::size_t nCells = mesh->getCellCount();

    return std::shared_ptr<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>>(new LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>(nCells));
}

}
