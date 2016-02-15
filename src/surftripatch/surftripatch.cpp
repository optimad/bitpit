/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#include "surftripatch.hpp"

namespace bitpit {

/*!
	\ingroup surftripatch
	@{
*/

/*!
	\class SurfTriPatch

	\brief The SurfTriPatch class defines an unstructured surface
	triangulation.

	SurfTriPatch defines an unstructured surface triangulation.
*/

/*!
	Creates a new patch.

	\param id is the id of the patch
*/
SurfTriPatch::SurfTriPatch(const int &id)
	: Patch(id, 2, true)
{

}

/*!
	Destroys the patch.
*/
SurfTriPatch::~SurfTriPatch()
{
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
void SurfTriPatch::setExpert(bool expert)
{
	Patch::setExpert(expert);
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double SurfTriPatch::evalCellVolume(const long &id)
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double SurfTriPatch::evalCellSize(const long &id)
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double SurfTriPatch::evalInterfaceArea(const long &id)
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> SurfTriPatch::evalInterfaceNormal(const long &id)
{
	BITPIT_UNUSED(id);

	return {{0., 0., 0.}};
}

/*!
	Updates the patch.

	\result Returns a vector of Adaption::Info that can be used to track
	the changes done during the update.
*/
const std::vector<Adaption::Info> SurfTriPatch::_update(bool trackAdaption)
{
	if (!isDirty()) {
		return std::vector<Adaption::Info>();
	}

	std::cout << ">> Updating surface triangulation mesh\n";

	// Adaption info
	std::vector<Adaption::Info> adaptionData;
	if (trackAdaption) {

	}

	// Done
	return adaptionData;
}

/*!
	Marks a cell for refinement.

	This is a void function since mesh refinement is not implemented
	for SurfTri patches.

	\param id is the id of the cell that needs to be refined
*/
bool SurfTriPatch::_markCellForRefinement(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Marks a cell for coarsening.

	This is a void function since mesh refinement is not implemented
	for SurfTri patches.

	\param id the cell to be refined
*/
bool SurfTriPatch::_markCellForCoarsening(const long &id)
{
	BITPIT_UNUSED(id);

	return false;
}

/*!
	Enables cell balancing.

	This is a void function since mesh refinement is not implemented
	for SurfTri patches.

	\param id is the id of the cell
	\param enabled defines if enable the balancing for the specified cell
*/
bool SurfTriPatch::_enableCellBalancing(const long &id, bool enabled)
{
	BITPIT_UNUSED(id);
	BITPIT_UNUSED(enabled);

	return false;
}

/*!
 * Checks if the specified point is inside the patch.
 *
 * \param[in] point is the point to be checked
 * \result Returns true if the point is inside the patch, false otherwise.
 */
bool SurfTriPatch::isPointInside(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

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
long SurfTriPatch::locatePoint(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	return false;
}

/*!
	@}
*/

}
