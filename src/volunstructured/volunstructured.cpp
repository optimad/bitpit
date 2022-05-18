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

#include "bitpit_common.hpp"

#include "volunstructured.hpp"

namespace bitpit {

/*!
	\class VolUnstructured
	\ingroup volumepatches

	\brief The VolUnstructured class defines a dummy unstructured volume
	patch.

	VolUnstructured defines a dummy unstructured volume patch. This is just
	a dummy interface, the real implementation will come in a future release.
*/

#if BITPIT_ENABLE_MPI==1
/*!
	Creates an uninitialized partitioned patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
VolUnstructured::VolUnstructured(MPI_Comm communicator, std::size_t haloSize)
	: VolumeKernel(communicator, haloSize, true)
#else
/*!
	Creates an uninitialized serial patch.
*/
VolUnstructured::VolUnstructured()
	: VolumeKernel(true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
VolUnstructured::VolUnstructured(int dimension, MPI_Comm communicator, std::size_t haloSize)
	: VolUnstructured(PatchManager::AUTOMATIC_ID, dimension, communicator, haloSize)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
*/
VolUnstructured::VolUnstructured(int dimension)
	: VolUnstructured(PatchManager::AUTOMATIC_ID, dimension)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
VolUnstructured::VolUnstructured(int id, int dimension, MPI_Comm communicator, std::size_t haloSize)
	: VolumeKernel(id, dimension, communicator, haloSize, true)
#else
/*!
	Creates a patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
*/
VolUnstructured::VolUnstructured(int id, int dimension)
	: VolumeKernel(id, dimension, true)
#endif
{
}

/*!
	Creates a clone of the pach.

	\result A clone of the pach.
*/
std::unique_ptr<PatchKernel> VolUnstructured::clone() const
{
	return std::unique_ptr<VolUnstructured>(new VolUnstructured(*this));
}

/*!
 * Enables or disables expert mode.
 *
 * When expert mode is enabled, it will be possible to change the
 * patch using low level functions (e.g., it will be possible to
 * add individual cells, add vertices, delete cells, ...).
 *
 * \param expert if true, the expert mode will be enabled
 */
void VolUnstructured::setExpert(bool expert)
{
	VolumeKernel::setExpert(expert);
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double VolUnstructured::evalCellVolume(long id) const
{
	const Cell &cell = getCell(id);

	ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
	std::size_t nCellVertices = cellVertexIds.size();
	BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
	getVertexCoords(nCellVertices, cellVertexIds.data(), vertexCoordinates);

	if (isThreeDimensional()) {
		return cell.evalVolume(vertexCoordinates);
	} else {
		return cell.evalArea(vertexCoordinates);
	}
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double VolUnstructured::evalCellSize(long id) const
{
	const Cell &cell = getCell(id);

	ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
	std::size_t nCellVertices = cellVertexIds.size();
	BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
	getVertexCoords(nCellVertices, cellVertexIds.data(), vertexCoordinates);

	return cell.evalSize(vertexCoordinates);
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double VolUnstructured::evalInterfaceArea(long id) const
{
	const Interface &interface = getInterface(id);

	ConstProxyVector<long> interfaceVertexIds = interface.getVertexIds();
	std::size_t nInterfaceVertices = interfaceVertexIds.size();
	BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nInterfaceVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
	getVertexCoords(nInterfaceVertices, interfaceVertexIds.data(), vertexCoordinates);

	if (isThreeDimensional()) {
		return interface.evalArea(vertexCoordinates);
	} else {
		return interface.evalLength(vertexCoordinates);
	}
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> VolUnstructured::evalInterfaceNormal(long id) const
{
	const Interface &interface = getInterface(id);

	ConstProxyVector<long> interfaceVertexIds = interface.getVertexIds();
	std::size_t nInterfaceVertices = interfaceVertexIds.size();
	BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nInterfaceVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
	getVertexCoords(nInterfaceVertices, interfaceVertexIds.data(), vertexCoordinates);

	std::array<double, 3> orientation = {{0., 0., 0.}};
	if (!isThreeDimensional()) {
		long ownerId = interface.getOwner();
		const Cell &owner = getCell(ownerId);

		ConstProxyVector<long> ownerVertexIds = owner.getVertexIds();
		const int nOwnerVertices = ownerVertexIds.size();

		const std::array<double, 3> &V_A = getVertex(ownerVertexIds[0]).getCoords();
		const std::array<double, 3> &V_B = getVertex(ownerVertexIds[1]).getCoords();
		const std::array<double, 3> &V_Z = getVertex(ownerVertexIds[nOwnerVertices - 1]).getCoords();

		orientation = crossProduct(V_B - V_A, V_Z - V_A);
	}

	return interface.evalNormal(vertexCoordinates, orientation);
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int VolUnstructured::_getDumpVersion() const
{
	const int DUMP_VERSION = 2;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void VolUnstructured::_dump(std::ostream &stream) const
{
	// Dump certices
	dumpVertices(stream);

	// Dump cells
	dumpCells(stream);

	// Dump cells
	dumpInterfaces(stream);
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void VolUnstructured::_restore(std::istream &stream)
{
	// Restore certices
	restoreVertices(stream);

	// Restore cells
	restoreCells(stream);

	// Restore cells
	restoreInterfaces(stream);
}

/*!
 * Checks if the specified point is inside the patch.
 *
 * \param[in] point is the point to be checked
 * \result Returns true if the point is inside the patch, false otherwise.
 */
bool VolUnstructured::isPointInside(const std::array<double, 3> &point) const
{
	BITPIT_UNUSED(point);

	throw std::runtime_error ("The function 'isPointInside' is not implemented yet");

	return false;
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the idof the cell
	\param[in] point is the point to be checked
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolUnstructured::isPointInside(long id, const std::array<double, 3> &point) const
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(point);

	throw std::runtime_error ("The function 'isPointInside' is not implemented yet");

	return false;
}

/*!
 * Locates the cell the contains the point.
 *
 * If the point is not inside the patch, the function returns the id of the
 * null element.
 *
 * \param[in] point is the point to be checked
 * \result Returns the linear id of the cell the contains the point. If the
 * point is not inside the patch, the function returns the id of the null
 * element.
 */
long VolUnstructured::locatePoint(const std::array<double, 3> &point) const
{
	BITPIT_UNUSED(point);

	throw std::runtime_error ("The function 'locatePoint' is not implemented yet");

	return Cell::NULL_ID;
}

}
