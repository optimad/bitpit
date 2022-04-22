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
VolumeKernel * LevelSetKernel::getMesh() const{
    return m_mesh ;
}

/*!
 * Updates the kernel after an adaption.
 *
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSetKernel::update( const std::vector<adaption::Info> &adaptionData ) {

    BITPIT_UNUSED( adaptionData );

    // Nothing to do
}

/*!
 * Checks if a plane intersects the cell
 * @param[in] id is the index of cell
 * @param[in] root is a point on the plane
 * @param[in] normal is the normal of the plane
 * @param[in] tolerance is the tolerance used for distance comparisons
 * @return true if intersect
 */
bool LevelSetKernel::intersectCellPlane( long id, const std::array<double,3> &root, const std::array<double,3> &normal, double tolerance ) {

    std::array<double,3> minPoint;
    std::array<double,3> maxPoint;
    m_mesh->evalCellBoundingBox(id, &minPoint, &maxPoint);

    int dim = m_mesh->getDimension();
    return CGElem::intersectPlaneBox( root, normal, minPoint, maxPoint, dim, tolerance);
}

/*!
 * Check if a point lies within the cell
 * @param[in] id is the cell index
 * @param[in] pointCoords are the point coordinates
 * @return true if point is inside, false otherwise
 */
bool LevelSetKernel::isPointInCell(long id, const std::array<double,3> &pointCoords) const {
    return m_mesh->isPointInside(id,pointCoords);
}

/*!
 * Checks if the specified cell is inside the given bounding box
 * @param[in] id is the id of the cell
 * @param[in] boxMin is the lower left point of the boungind box
 * @param[in] boxMax is the upper right point of the boungind box
 * @result True if the specified cell is inside the given bounding box, false
 * otherwise.
 */
double LevelSetKernel::isCellInsideBoundingBox( long id, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax ) const {

    double distanceTolerance = m_mesh->getTol();

    std::array<double,3> cellBoxMin;
    std::array<double,3> cellBoxMax;
    m_mesh->evalCellBoundingBox(id, &cellBoxMin, &cellBoxMax);

    return CGElem::intersectBoxBox(cellBoxMin, cellBoxMax, boxMin, boxMax, 3, distanceTolerance);
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

/*!
    Create a data communicator for exchanging data among partitions.

    \result A data communicator for exchanging data among partitions.
*/
std::unique_ptr<DataCommunicator> LevelSetKernel::createDataCommunicator( ) const {

    std::unique_ptr<DataCommunicator> dataCommunicator;
    if (isCommunicatorSet()) {
        dataCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(getCommunicator())) ;
    }

    return dataCommunicator;
}

#endif

/*!
    Create the sign propagator.

    The sign propagator allow to propagate the levelset sign from the narrow
    band to the rest of the domain.

    \result The newlycreated sign propagator.
*/
std::unique_ptr<LevelSetSignPropagator> LevelSetKernel::createSignPropagator( ) const {

    return std::unique_ptr<LevelSetSignPropagator>(new LevelSetSignPropagator(m_mesh)) ;
}

}
