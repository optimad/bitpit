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

#ifndef __BITPIT_SURFACE_KERNEL_HPP__
#define __BITPIT_SURFACE_KERNEL_HPP__

#include "patch_kernel.hpp"

namespace bitpit {

class SurfaceKernel : public PatchKernel {

public:
	SurfaceKernel(const int &id, const int &dimension, bool epxert);

	virtual ~SurfaceKernel();

	virtual double evalCellVolume(const long &id) = 0;
	virtual double evalCellSize(const long &id) = 0;

	virtual double evalInterfaceArea(const long &id) = 0;

};

}

#endif
