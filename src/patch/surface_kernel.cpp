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

#include "surface_kernel.hpp"

namespace bitpit {


/*!
	\ingroup PatchKernel
	@{
*/

/*!
	\class SurfaceKernel

	\brief The SurfaceKernel class provides an interface for defining
	surface patches.

	SurfaceKernel is the base class for defining surface patches.
*/

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
SurfaceKernel::SurfaceKernel(const int &id, const int &dimension, const int& space_dim, bool expert)
	: PatchKernel(id, dimension, expert)
{
    m_spaceDim = space_dim;
}

/*!
	Destroys the patch.
*/
SurfaceKernel::~SurfaceKernel()
{

}

/*!
 * Returns the number of dimensions of the working space (set at patch construction)
 */
int SurfaceKernel::getSpaceDimensions(void) const
{
    return(m_spaceDim);
}
/*!
	@}
*/

}
