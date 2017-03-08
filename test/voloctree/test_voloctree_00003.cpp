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

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

	log::manager().initialize(log::COMBINED);
	log::cout() << "Testing octree patch" << "\n";

	std::array<double, 3> origin = {{0., 0., 0.}};
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

	log::cout() << "  >> 2D octree patch" << "\n";

	VolOctree *patch_2D = new VolOctree(0, 2, origin, length, dh);
	patch_2D->getVTK().setName("octree_uniform_patch_2D");
	patch_2D->update();

	log::cout() << "\n  >> 2D location test" << std::endl;
	log::cout() << std::endl;

	for (auto testPoint : pointList) {
		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_2D->isPointInside(testPoint)) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_2D->locatePoint(testPoint) << std::endl;

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_2D->isPointInside(testPoint[0], testPoint[1], testPoint[2])) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_2D->locatePoint(testPoint[0], testPoint[1], testPoint[2]) << std::endl;
	}

	log::cout() << std::endl;

	delete patch_2D;

	log::cout() << "  >> 3D octree patch" << "\n";

	VolOctree *patch_3D = new VolOctree(0, 3, origin, length, dh);
	patch_3D->getVTK().setName("octree_uniform_patch_3D");
	patch_3D->update();

	log::cout() << "\n  >> 3D location test" << std::endl;
	log::cout() << std::endl;

	for (auto testPoint : pointList) {
		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_3D->isPointInside(testPoint)) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_3D->locatePoint(testPoint) << std::endl;

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_3D->isPointInside(testPoint[0], testPoint[1], testPoint[2])) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_3D->locatePoint(testPoint[0], testPoint[1], testPoint[2]) << std::endl;
	}

	log::cout() << std::endl;

	delete patch_3D;

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif

}
