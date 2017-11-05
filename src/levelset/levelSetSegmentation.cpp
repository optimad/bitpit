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
    const long *segmentConnect = segment.getConnect();

    coords->resize(nVertices);
    for (int n = 0; n < nVertices; ++n) {
        long vertexId = segmentConnect[n] ;
        (*coords)[n] = m_surface->getVertexCoords(vertexId);
    }

    if ( nVertices > 3 ) {
        log::cout() << "levelset: only segments and triangles supported in LevelSetSegmentation !!" << std::endl ;
    }
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
void SegmentationKernel::getSegmentInfo( const std::array<double,3> &pointCoords, const long &segmentId, const bool &signd, double &distance, std::array<double,3> &gradient, std::array<double,3> &normal ) const {


    std::array<double,3> outwards;

    auto itrNormal = getVertexNormals().find(segmentId) ;
    auto itrGradient = getVertexGradients().find(segmentId) ;
    assert( itrGradient != getVertexGradients().end() ) ;

    const Cell &cell = m_surface->getCell(segmentId) ;
    ElementType cellType = cell.getType();
    const long *cellConnect = cell.getConnect();

    std::array<double,3> projectionCoords;

    switch (cellType) {

    case ElementType::VERTEX :
    {
        long id = cellConnect[0] ;

        projectionCoords = m_surface->getVertexCoords(id);


        normal.fill(0.);

        break;
    }

    case ElementType::LINE:
    {
        long id0 = cellConnect[0] ;
        long id1 = cellConnect[1] ;
        std::array<double,2> lambda ;

        projectionCoords = CGElem::projectPointSegment( pointCoords, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), lambda);
        outwards  = lambda[0] *itrGradient->second[0] ;
        outwards += lambda[1] *itrGradient->second[1] ;
        outwards /= norm2(outwards) ;

        if( itrNormal != getVertexNormals().end() ){
            normal  = lambda[0] *itrNormal->second[0] ;
            normal += lambda[1] *itrNormal->second[1] ;
            normal /= norm2(normal) ;

        } else {
            normal = outwards;

        }

        break;
    }

    case ElementType::TRIANGLE:
    {
        long id0 = cellConnect[0] ;
        long id1 = cellConnect[1] ;
        long id2 = cellConnect[2] ;

        std::array<double,3> lambda ;

        projectionCoords = CGElem::projectPointTriangle( pointCoords, m_surface->getVertexCoords(id0), m_surface->getVertexCoords(id1), m_surface->getVertexCoords(id2), lambda );
        outwards  = lambda[0] *itrGradient->second[0] ;
        outwards += lambda[1] *itrGradient->second[1] ;
        outwards += lambda[2] *itrGradient->second[2] ;
        outwards /= norm2(outwards);


        if( itrNormal != getVertexNormals().end() ){
            normal  = lambda[0] *itrNormal->second[0] ;
            normal += lambda[1] *itrNormal->second[1] ;
            normal += lambda[2] *itrNormal->second[2] ;
            normal /= norm2(normal) ;

        } else {
            normal = outwards;

        }

        break;
    }

    default:
    {
        std::runtime_error ("Type of cell not supported.");
        break;
    }

    }

    gradient = pointCoords-projectionCoords;
    distance = norm2(gradient); 
    gradient /= distance;

    // the sign is computed by determining the side of point p
    // with respect to the normal plane 
    double s = sign( dotProduct(gradient, outwards) );

    // if p lies on the normal plane (s=0), but the distance is finite the sign must be evaluated
    // considering the curvature of the surface. Anyhow this situation is not crucial because
    // there should exists another element with smaller distance. The signed is put aribtrarly
    // to positive
    if(utils::DoubleFloatingEqual()(s,0.) && distance>0){
        s = 1.;
    } 

    // If signed distance are computed, the distance value and gradient
    // need to be changed accordingly. If unsigned distance are computed
    // the orientation of the suraface normal is discarded and in order
    // to agnostic with repect the two sides of the surface
    distance *= ( signd *s  + (!signd) *1.);
    gradient *= ( signd *s  + (!signd) *1.);
    normal   *= ( signd *1. + (!signd) *s );

    return;

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
LevelSetSegmentation::SurfaceInfo::SurfaceInfo( long index, std::array<double,3> vec ) : support(index), normal(vec)
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
std::array<double,3> LevelSetSegmentation::getNormal( const long &id ) const{

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
long LevelSetSegmentation::getSupport( const long &id ) const{

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
 * Finds seed points in narrow band within a cartesian mesh for one simplex
 * @param[in] visitee cartesian mesh 
 * @param[in] VS Simplex
 * @param[out] I indices of seed points
 */
bool LevelSetSegmentation::seedNarrowBand( LevelSetCartesian *visitee, std::vector<std::array<double,3>> &VS, double searchRadius, std::vector<long> &I){

    VolCartesian                        &mesh = *(static_cast<VolCartesian*>(visitee->getMesh()));

    bool                                found(false) ;
    int                                 dim( mesh.getDimension() ) ;
    std::array<double,3>                B0, B1;
    std::vector<std::array<double,3>>   VP ;

    mesh.getBoundingBox(B0, B1) ;

    for( int i=0; i<dim; ++i){
        B0[i] -= searchRadius;
        B1[i] += searchRadius;
    }

    I.clear() ;

    for( const auto &P : VS){
        if(  CGElem::intersectPointBox( P, B0, B1, dim ) ) {
            I.push_back( mesh.locateClosestCell(P) );
            found =  true ;
        }
    }

    if( !found && CGElem::intersectBoxPolygon( B0, B1, VS, false, true, true, VP, dim ) ) {
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
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetCartesian *visitee, bool signd){

    VolCartesian &mesh = *(visitee->getCartesianMesh() ) ;
    double searchRadius = m_narrowBand;

    if(searchRadius<0.){
        for( int d=0; d < mesh.getDimension(); ++d){
            searchRadius = std::max( searchRadius, mesh.getSpacing(d) ) ;
        }
    }

    std::vector<std::array<double,3>>       VS(3);

    const SurfUnstructured                  &m_surface = m_segmentation->getSurface();

    std::vector<long>                       stack, temp, neighs, flag( mesh.getCellCount(), -1);

    std::vector< std::array<double,3> >     cloud ;
    std::vector<double>                     cloudDistance;

    double distance;
    std::array<double,3>  gradient, normal;

    stack.reserve(128) ;
    temp.reserve(128) ;

    log::cout() << " Compute levelset on cartesian mesh"  << std::endl;

    for (const Cell &segment : m_surface.getCells()) {
        long segmentId = segment.getId();

        // compute initial seeds, ie the cells where the vertices
        // of the surface element fall in and add them to stack
        m_segmentation->getSegmentVertexCoords( segmentId, &VS ) ;
        seedNarrowBand( visitee, VS, searchRadius, stack );

        ElementType segmentType = segment.getType();


        // propagate from seed
        size_t stackSize = stack.size();
        while (stackSize > 0) {

            // put the cell centroids of the stack into a vector
            // and calculate the distances to the cloud
            cloud.resize(stackSize) ;
            cloudDistance.resize(stackSize);

            for( size_t k = 0; k < stackSize; ++k) {
                long cell = stack[k];
                cloud[k] = visitee->computeCellCentroid(cell) ;
            }

            switch (segmentType) {

            case ElementType::VERTEX :
            {
                for( size_t k=0; k<stackSize; ++k){
                    cloudDistance[k] = norm2( cloud[k]-VS[0]);
                }
                break;
            }

            case ElementType::LINE:
            {
                for( size_t k=0; k<stackSize; ++k){
                    cloudDistance[k] = CGElem::distancePointSegment( cloud[k], VS[0], VS[1]);
                }
                break;
            }

            case ElementType::TRIANGLE:
            {
                cloudDistance = CGElem::distanceCloudTriangle( cloud, VS[0], VS[1], VS[2]); 
                break;
            }

            default:
            {
                std::runtime_error ("Type of cell not supported.");
                break;
            }
            }

            // check each cell of cloud individually
            for( size_t k = 0; k < stackSize; ++k) {

                long &cellId = stack[k];
                double &cellDistance = cloudDistance[k];

                // consider only cells within the search radius
                if ( cellDistance <= searchRadius ) {

                    PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.find(cellId) ;
                    if( lsInfoItr == m_ls.end() ){
                        lsInfoItr = m_ls.emplace(cellId) ;
                    }

                    // check if the computed distance is the closest distance
                    if( cellDistance < std::abs(lsInfoItr->value) ){

                        // compute all necessary information and store them
                        m_segmentation->getSegmentInfo(cloud[k], segmentId, signd, distance, gradient, normal);

                        lsInfoItr->value    = distance;
                        lsInfoItr->gradient = gradient;
                
                        auto infoItr = m_surfaceInfo.find(cellId) ;
                        if( infoItr == m_surfaceInfo.end() ){
                            infoItr = m_surfaceInfo.emplace(cellId) ;
                        }
                        infoItr->support = segmentId;
                        infoItr->normal = normal;
                    }


                    // the new stack is composed of all neighbours
                    // of the old stack. Attention must be paid in 
                    // order not to evaluate the same cell twice
                    neighs.clear();
                    mesh.findCellFaceNeighs(cellId, &neighs) ;
                    for( const auto &  neigh : neighs){
                        if( flag[neigh] != segmentId) {
                            temp.push_back( neigh) ;
                            flag[neigh] = segmentId ;
                        }
                    }

                } //end if distance
            }

            stack.clear() ;
            stack.swap( temp ) ;
            stackSize = stack.size() ;

        } //end while
    }

    return;
}

/*!
 */
void LevelSetSegmentation::computeLSInNarrowBand( LevelSetOctree *visitee, bool signd){

    VolumeKernel &mesh = *(visitee->getMesh()) ;

    bool adaptiveSearch(m_narrowBand<0);
    double searchRadius = m_narrowBand;
    double factor = 0.5 *sqrt( (double) mesh.getDimension() );

    long segmentId;
    double distance;
    std::array<double,3> gradient, normal;
    std::array<double,3> centroid, root;

    std::unordered_set<long> intersects;

    for( const Cell &cell : mesh.getCells() ){

        long cellId = cell.getId();

        centroid = visitee->computeCellCentroid(cellId);

        if(adaptiveSearch){
            double cellSize = mesh.evalCellSize(cellId);
            searchRadius = factor *cellSize;
        }

        m_segmentation->m_searchTreeUPtr->findPointClosestCell(centroid, searchRadius, &segmentId, &distance);

        if(segmentId>=0){

            m_segmentation->getSegmentInfo(centroid, segmentId, signd, distance, gradient, normal);

            PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.find(cellId);
            lsInfoItr = m_ls.emplace(cellId) ;
            lsInfoItr->value    = distance;
            lsInfoItr->gradient = gradient;

            PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId);
            infoItr->support = segmentId;
            infoItr->normal = normal;

            if(adaptiveSearch){
                intersects.insert(cellId);
            }

        }
             
    }

    if(!adaptiveSearch){
        return;
    }
    
    for( const long &cellId : intersects){

        Cell const &cell = mesh.getCell(cellId);
        
        centroid = visitee->computeCellCentroid(cellId);
        root = computeProjectionPoint(cellId);

        const long *neighbours = cell.getAdjacencies() ;
        int nNeighbours = cell.getAdjacencyCount() ;
        for (int n = 0; n < nNeighbours; ++n) {
            long neighId = neighbours[n];

            if (neighId < 0) {
                continue;
            }

            // skip if neigh cell has already been processed
            // either because it is intersected by surface or
            // because it is a neigh to a previous intersect
            if( m_ls.count(neighId) != 0 ){
                continue;
            }

            centroid = visitee->computeCellCentroid(neighId);
            searchRadius =  1.05 *norm2(centroid-root);
            m_segmentation->m_searchTreeUPtr->findPointClosestCell(centroid, searchRadius, &segmentId, &distance);

            if(segmentId>=0){

                m_segmentation->getSegmentInfo(centroid, segmentId, signd, distance, gradient, normal);

                PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(neighId) ;
                lsInfoItr->value    = distance;
                lsInfoItr->gradient = gradient;


                PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(neighId);
                infoItr->support = segmentId;
                infoItr->normal = normal;

            } else {
                assert(false && "Should not pass here");
            }

        }
    }
}

/*!
 */
void LevelSetSegmentation::updateLSInNarrowBand( LevelSetOctree *visitee, const std::vector<adaption::Info> &mapper, bool signd){

    VolumeKernel &mesh = *(visitee->getMesh()) ;

    bool adaptiveSearch(m_narrowBand<0);
    double searchRadius = m_narrowBand;
    double factor = 0.5 *sqrt( (double) mesh.getDimension() );

    long segmentId;
    std::array<double,3> centroid, root;

    double distance;
    std::array<double,3> gradient, normal;

    std::vector<long> unprocessed;

    for( const auto &event : mapper ){

        if( event.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for( const long &cellId : event.current ){
            centroid = visitee->computeCellCentroid(cellId);

            if(adaptiveSearch){
                double cellSize = mesh.evalCellSize(cellId);
                searchRadius =  factor *cellSize;
            }

            m_segmentation->m_searchTreeUPtr->findPointClosestCell(centroid, searchRadius, &segmentId, &distance);

            if(segmentId>=0){

                m_segmentation->getSegmentInfo(centroid, segmentId, signd, distance, gradient, normal);

                PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(cellId) ;
                lsInfoItr->value    = distance;
                lsInfoItr->gradient = gradient;

                PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId) ;
                infoItr->support = segmentId;
                infoItr->normal = normal;

            } else if(adaptiveSearch){
                unprocessed.push_back(cellId);
            }

        }
    }

    if(!adaptiveSearch){
        return;
    }

    for( const long &cellId : unprocessed){

        const Cell &cell = mesh.getCell(cellId);

        const long *neighbours = cell.getAdjacencies() ;
        int nNeighbours = cell.getAdjacencyCount() ;

        int n=0;
        bool iterate = (nNeighbours>0);

        while(iterate){

            long neighId = neighbours[n];

            if(neighId>=0){
                if( intersectSurface(neighId,LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) == LevelSetIntersectionStatus::TRUE){

                    centroid = visitee->computeCellCentroid(cellId);
                    root = computeProjectionPoint(neighId);

                    searchRadius =  1.05 *norm2(centroid-root);
                    m_segmentation->m_searchTreeUPtr->findPointClosestCell(centroid, searchRadius, &segmentId, &distance);

                    if(segmentId>=0){

                        m_segmentation->getSegmentInfo(centroid, segmentId, signd, distance, gradient, normal);

                        PiercedVector<LevelSetInfo>::iterator lsInfoItr = m_ls.emplace(cellId) ;
                        lsInfoItr->value    = distance;
                        lsInfoItr->gradient = gradient;


                        PiercedVector<SurfaceInfo>::iterator infoItr = m_surfaceInfo.emplace(cellId);
                        infoItr->support = segmentId;
                        infoItr->normal = normal;

                    } else {
                        assert(false && "Should not pass here");
                    }
                    iterate = false;
                }
            }

            ++n;
            iterate &= (n<nNeighbours);
        }
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

        // Delete only old data that belongs to the current processor
        if (info.type == adaption::Type::TYPE_PARTITION_RECV) {
            continue;
        }

        // Remove info of previous cells
        for ( const long & parent : info.previous ) {
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

}
