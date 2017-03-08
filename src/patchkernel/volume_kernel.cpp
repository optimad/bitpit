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

#include "volume_kernel.hpp"

namespace bitpit {


/*!
	\ingroup volumepatches
	@{
*/

/*!
	\class VolumeKernel

	\brief The VolumeKernel class provides an interface for defining
	volume patches.

	VolumeKernel is the base class for defining voulme patches.
*/

/*!
	Creates a new patch.

	\param expert if true, the expert mode will be enabled
*/
VolumeKernel::VolumeKernel(bool expert)
	: PatchKernel(expert)
{
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
VolumeKernel::VolumeKernel(const int &dimension, bool expert)
	: PatchKernel(dimension, expert)
{
}

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
VolumeKernel::VolumeKernel(const int &id, const int &dimension, bool expert)
	: PatchKernel(id, dimension, expert)
{
}

/*!
	Destroys the patch.
*/
VolumeKernel::~VolumeKernel()
{

}

/*!
	Checks if the specified point is inside the patch.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool VolumeKernel::isPointInside(const double &x, const double &y, const double &z)
{
	return isPointInside({{x, y, z}});
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the index of the cells
	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolumeKernel::isPointInside(const long &id, const double &x, const double &y, const double &z)
{
	return isPointInside(id, {{x, y, z}});
}

/*!
	@}
*/

}
