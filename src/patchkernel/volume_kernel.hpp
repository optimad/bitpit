/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#ifndef __BITPIT_VOLUME_KERNEL_HPP__
#define __BITPIT_VOLUME_KERNEL_HPP__

#include "patch_kernel.hpp"

namespace bitpit {

class VolumeKernel : public PatchKernel {

public:
	virtual ~VolumeKernel();

	bool isPointInside(const double &x, const double &y, const double &z);
	virtual bool isPointInside(const std::array<double, 3> &point) = 0;
	bool isPointInside(const long &id, const double &x, const double &y, const double &z);
	virtual bool isPointInside(const long &id, const std::array<double, 3> &point) = 0;

	virtual double evalCellVolume(const long &id)const = 0;
	virtual double evalCellSize(const long &id)const = 0;

	virtual double evalInterfaceArea(const long &id)const = 0;
        virtual std::array<double,3> evalInterfaceNormal(const long &id)const = 0;

protected:
	VolumeKernel(bool epxert);
	VolumeKernel(const int &dimension, bool epxert);
	VolumeKernel(const int &id, const int &dimension, bool epxert);

};

}

#endif
