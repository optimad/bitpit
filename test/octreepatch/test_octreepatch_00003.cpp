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

#include <array>
#if ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_octreepatch.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc,argv);
#endif

	std::cout << "Testing octree patch" << "\n";

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	std::array<double, 3> point;
	std::vector<std::array<double, 3>> pointList;

	pointList.push_back(origin);
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 4;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 2;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 2 + 0.01;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] - length / 2;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] - length / 4;
		pointList.push_back(point);
	}

	std::cout << "  >> 2D octree patch" << "\n";

	Patch *patch_2D = new OctreePatch(0, 2, origin, length, dh);
	patch_2D->setName("octree_uniform_patch_2D");
	patch_2D->update();

	std::cout << "\n  >> 2D location test" << std::endl;
	std::cout << std::endl;

	for (auto testPoint : pointList) {
		std::cout << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_2D->isPointInside(testPoint)) {
			std::cout << " is inside the patch" << std::endl;
		} else {
			std::cout << " is outside the patch" << std::endl;
		}

		std::cout << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_2D->isPointInside(testPoint[0], testPoint[1], testPoint[2])) {
			std::cout << " is inside the patch" << std::endl;
		} else {
			std::cout << " is outside the patch" << std::endl;
		}
	}

	std::cout << std::endl;

	delete patch_2D;

	std::cout << "  >> 3D octree patch" << "\n";

	Patch *patch_3D = new OctreePatch(0, 3, origin, length, dh);
	patch_3D->setName("octree_uniform_patch_3D");
	patch_3D->update();

	std::cout << "\n  >> 3D location test" << std::endl;
	std::cout << std::endl;

	for (auto testPoint : pointList) {
		std::cout << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_3D->isPointInside(testPoint)) {
			std::cout << " is inside the patch" << std::endl;
		} else {
			std::cout << " is outside the patch" << std::endl;
		}


		std::cout << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_3D->isPointInside(testPoint[0], testPoint[1], testPoint[2])) {
			std::cout << " is inside the patch" << std::endl;
		} else {
			std::cout << " is outside the patch" << std::endl;
		}

	}

	std::cout << std::endl;

	delete patch_3D;

#if ENABLE_MPI==1
	MPI::Finalize();
#endif

}
