/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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


# if BITPIT_ENABLE_MPI
# include "bitpit_communications.hpp"
# endif

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetObject.hpp"
# include "levelSetKernel.hpp"

namespace bitpit {

/*!
    @interface  LevelSetKernel
    @ingroup levelset
    @brief  Mesh specific implementation to calculate the levelset function

*/

/*!
 * Default constructor.
 */
LevelSetKernel::LevelSetKernel() {
    m_mesh = NULL ;

#if BITPIT_ENABLE_MPI
    m_communicator = MPI_COMM_NULL;
# endif

}

/*!
 * Constructor
 * @param[in] patch underlying mesh
 */
LevelSetKernel::LevelSetKernel( VolumeKernel *patch): LevelSetKernel() {
    m_mesh = patch ;
}

/*!
 * Destructor of LevelSetKernel
*/
LevelSetKernel::~LevelSetKernel(){
    m_mesh = NULL ;

# if BITPIT_ENABLE_MPI
    freeCommunicator();
# endif

}

/*!
 * Returns pointer to underlying mesh.
 * @return pointer to mesh
*/
VolumeKernel* LevelSetKernel::getMesh() const{
    return m_mesh ;
} 

/*!
 * Clears the geometry cache.
 */
void LevelSetKernel::clearGeometryCache(  ) {

    std::unordered_map<long, std::array<double,3>>().swap( m_cellCentroids ) ;

}

/*!
 * Updates the geometry cache after an adaption.
 */
void LevelSetKernel::updateGeometryCache( const std::vector<adaption::Info> &mapper ) {

    // If there are no cells in the mesh we can just delete all the cache
    if ( m_mesh->getCellCount() == 0) {
        clearGeometryCache();
        return;
    }

    // Remove the previous cells from the cache
    for ( auto & map : mapper ){
        if( map.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for ( auto & previousId : map.previous){
            auto centroidItr = m_cellCentroids.find( previousId ) ;
            if ( centroidItr == m_cellCentroids.end() ) {
                continue ;
            }

            m_cellCentroids.erase( previousId ) ;
        }
    }

}

/*!
 * Computes the centroid of the specfified cell.
 *
 * If the centroid of the cell has been already evaluated, the cached value is
 * returned. Otherwise the cell centroid is evaluated and stored in the cache.
 *
 * @param[in] id is the index of cell
 * @return The centroid of the cell.
 */
const std::array<double,3> & LevelSetKernel::computeCellCentroid( long id ) {

    auto centroidItr = m_cellCentroids.find( id ) ;
    if ( centroidItr == m_cellCentroids.end() ) {
        centroidItr = m_cellCentroids.insert( { id, m_mesh->evalCellCentroid( id ) } ).first ;
    }

    return centroidItr->second;

}

/*!
 * Computes the radius of the incircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetKernel::computeCellIncircle( long id ) {

    VolumeKernel *patch = getMesh();
    Cell &cell = patch->getCell(id);

    const long* interfaceIds = cell.getInterfaces();
    int interfaceCount = cell.getInterfaceCount();
    
    std::array<double,3> cellCenter = computeCellCentroid(id);

    double radius = std::numeric_limits<double>::max() ;
    for (int k = 0; k < interfaceCount; ++k) {
        long interfaceId = interfaceIds[k];
        double r = norm2(cellCenter - patch->evalInterfaceCentroid(interfaceId));
        radius = std::min(radius, r);
    }

    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();
    for (int k = 0; k < nCellVertices; ++k) {
        long vertexId = cellVertexIds[k];
        double r = norm2(cellCenter - patch->getVertexCoords(vertexId));
        radius = std::min(radius, r);
    }

    return radius;

}

/*!
 * Computes the radius of the circumcircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetKernel::computeCellCircumcircle( long id ) {

    VolumeKernel *patch = getMesh();
    Cell &cell = patch->getCell(id);

    const long* interfaceIds = cell.getInterfaces();
    int interfaceCount = cell.getInterfaceCount();

    std::array<double,3> cellCenter = computeCellCentroid(id);

    double radius = -std::numeric_limits<double>::max() ;
    for (int k = 0; k < interfaceCount; ++k) {
        long interfaceId = interfaceIds[k];
        double r = norm2(cellCenter - patch->evalInterfaceCentroid(interfaceId));
        radius = std::max(radius, r);
    }

    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();
    for (int k = 0; k < nCellVertices; ++k) {
        long vertexId = cellVertexIds[k];
        double r = norm2(cellCenter - patch->getVertexCoords(vertexId));
        radius = std::max(radius, r);
    }

    return radius;

}

/*!
 * Check if a point lies within the cell
 * @param[in] id is the cell index
 * @param[in] pointCoords are the point coordinates
 * @return true if point is inside, false otherwise
 */
bool LevelSetKernel::isPointInCell(long id, const std::array<double,3> &pointCoords){
    return getMesh()->isPointInside(id,pointCoords);
}

/*!
 * Checks if the specified cell is inside the given bounding box
 * @param[in] id is the id of the cell
 * @param[in] minPoint is the lower left point of the boungind box
 * @param[in] maxPoint is the upper right point of the boungind box
 * @result True if the specified cell is inside the given bounding box, false
 * otherwise.
 */
double LevelSetKernel::isCellInsideBoundingBox( long id, std::array<double, 3> minPoint, std::array<double, 3> maxPoint ){

    const Cell &cell = m_mesh->getCell(id);
    double tolerance = m_mesh->getTol();

    std::array<double, 3> cellMinPoint;
    std::array<double, 3> cellMaxPoint;

    cellMinPoint.fill(   std::numeric_limits<double>::max() ) ;
    cellMaxPoint.fill( - std::numeric_limits<double>::max() ) ;

    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();
    for (int i = 0; i < nCellVertices; ++i) {
        long vertexId = cellVertexIds[i];
        std::array<double, 3> vertexCoords = m_mesh->getVertexCoords(vertexId);
        for (int d = 0; d < 3; ++d) {
            cellMinPoint[d] = std::min( vertexCoords[d] - tolerance, cellMinPoint[d]) ;
            cellMaxPoint[d] = std::max( vertexCoords[d] + tolerance, cellMaxPoint[d]) ;
        }
    }

    return CGElem::intersectBoxBox(minPoint, maxPoint, cellMinPoint, cellMaxPoint);
}

# if BITPIT_ENABLE_MPI

/*!
	Initializes the MPI communicator to be used for parallel communications.
*/
void LevelSetKernel::initializeCommunicator()
{
	// The communicator has to be valid
	MPI_Comm communicator = m_mesh->getCommunicator();
	if (communicator == MPI_COMM_NULL) {
		throw std::runtime_error ("Levelset communicator is not valid");
	}

	// Communication can be set just once
	if (isCommunicatorSet()) {
		int comparison;
		MPI_Comm_compare(communicator, m_communicator, &comparison);
		if (comparison == MPI_IDENT || comparison == MPI_CONGRUENT) {
			return;
		} else {
			throw std::runtime_error ("Levelset communicator cannot be modified.");
		}
	}

	// Create a duplicate of the patch communicator
	MPI_Comm_dup(communicator, &m_communicator);
}

/*!
 * Returns the MPI communicator stored within LevelSetKernel
 * @return MPI communicator
 */
MPI_Comm LevelSetKernel::getCommunicator() const {

    return m_communicator;
}

/*!
    Checks if the communicator to be used for parallel communications has
    already been set.

    \result Returns true if the communicator has been set, false otherwise.
*/
bool LevelSetKernel::isCommunicatorSet() const {

    return (getCommunicator() != MPI_COMM_NULL);
}

/*!
    Frees the MPI communicator associated to the patch
*/
void LevelSetKernel::freeCommunicator() {

    if (!isCommunicatorSet()) {
        return;
    }

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled) {
        return;
    }

    MPI_Comm_free(&m_communicator);
}

#endif

}
