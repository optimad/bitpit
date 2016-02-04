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

#include "bitpit_cartesianpatch.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc,argv);
#endif

	std::cout << "Testing Cartesian patch" << "\n";

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	std::cout << "  >> 2D Cartesian patch" << "\n";

	Patch *patch_2D = new CartesianPatch(0, 2, origin, length, dh);
	patch_2D->setName("cartesian_uniform_patch_2D");
	patch_2D->update();
	patch_2D->writeMesh();

	std::cout << "\n  >> 2D neighbour test" << "\n";

	std::vector<long> neighs_2D;

	long cellId_2D = 7;
	std::cout << std::endl;
	std::cout << "Cell id: " << cellId_2D << std::endl << std::endl;

	std::cout << "Face neighbours (complete list): " << std::endl;
	neighs_2D = patch_2D->extractCellFaceNeighs(cellId_2D);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (complete list): " << std::endl;
	neighs_2D = patch_2D->extractCellVertexNeighs(cellId_2D, true);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (excuding face neighbours): " << std::endl;
	neighs_2D = patch_2D->extractCellVertexNeighs(cellId_2D, false);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		std::cout << " - " << neighs_2D[i] << std::endl;
	}

	std::cout << std::endl;

	delete patch_2D;

	std::cout << "  >> 3D Cartesian patch" << "\n";

	Patch *patch_3D = new CartesianPatch(0, 3, origin, length, dh);
	patch_3D->setName("cartesian_uniform_patch_3D");
	patch_3D->update();
	patch_3D->writeMesh();

	std::cout << "\n  >> 3D neighbour test" << "\n";

	std::vector<long> neighs_3D;

	long cellId_3D = 13;
	std::cout << std::endl;
	std::cout << "Cell id: " << cellId_3D << std::endl << std::endl;

	std::cout << "Face neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->extractCellFaceNeighs(cellId_3D);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Edge neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->extractCellEdgeNeighs(cellId_3D, true);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Edge neighbours (excuding face neighbours): " << std::endl;
	neighs_3D = patch_3D->extractCellEdgeNeighs(cellId_3D, false);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->extractCellVertexNeighs(cellId_3D, true);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << "Vertex neighbours (excuding face and edge neighbours): " << std::endl;
	neighs_3D = patch_3D->extractCellVertexNeighs(cellId_3D, false);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		std::cout << " - " << neighs_3D[i] << std::endl;
	}

	std::cout << std::endl;

	delete patch_3D;

#if ENABLE_MPI==1
	MPI::Finalize();
#endif

}
