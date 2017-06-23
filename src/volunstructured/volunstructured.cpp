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

/*!
	Creates a new patch.

	\param id is the id of the patch
	\param dimension is the dimension of the patch
*/
VolUnstructured::VolUnstructured(const int &id, const int &dimension)
	: VolumeKernel(id, dimension, true)
{

}

/*!
	Destroys the patch.
*/
VolUnstructured::~VolUnstructured()
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
void VolUnstructured::setExpert(bool expert)
{
	VolumeKernel::setExpert(expert);
}

/*!
	Evaluates the volume of the specified cell.

	\param id is the id of the cell
	\result The volume of the specified cell.
*/
double VolUnstructured::evalCellVolume(const long &id) const
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the characteristic size of the specified cell.

	\param id is the id of the cell
	\result The characteristic size of the specified cell.
*/
double VolUnstructured::evalCellSize(const long &id) const
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the area of the specified interface.

	\param id is the id of the interface
	\result The area of the specified interface.
*/
double VolUnstructured::evalInterfaceArea(const long &id) const
{
	BITPIT_UNUSED(id);

	return 0;
}

/*!
	Evaluates the normal of the specified interface.

	\param id is the id of the interface
	\result The normal of the specified interface.
*/
std::array<double, 3> VolUnstructured::evalInterfaceNormal(const long &id) const
{
	BITPIT_UNUSED(id);

	return {{0., 0., 0.}};
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int VolUnstructured::_getDumpVersion() const
{
	const int DUMP_VERSION = 1;

	return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void VolUnstructured::_dump(std::ostream &stream)
{
    BITPIT_UNUSED(stream);

	throw std::runtime_error ("Dump of volunstructured is not implemented yet.");
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void VolUnstructured::_restore(std::istream &stream)
{
    BITPIT_UNUSED(stream);

	throw std::runtime_error ("Restore of volunstructured is not implemented yet.");
}

/*!
 * Checks if the specified point is inside the patch.
 *
 * \param[in] point is the point to be checked
 * \result Returns true if the point is inside the patch, false otherwise.
 */
bool VolUnstructured::isPointInside(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	return false;
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the idof the cell
	\param[in] point is the point to be checked
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolUnstructured::isPointInside(const long &id, const std::array<double, 3> &point)
{
	BITPIT_UNUSED(id);
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
long VolUnstructured::locatePoint(const std::array<double, 3> &point)
{
	BITPIT_UNUSED(point);

	return false;
}

}
