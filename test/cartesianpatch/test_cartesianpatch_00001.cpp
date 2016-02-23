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

	std::array<double, 3> minPoint;
	std::array<double, 3> maxPoint;

	std::vector<double> cellData;
	std::vector<double> vertexData;

	log::cout() << "Testing Cartesian patch" << "\n";

	std::array<double, 3> origin = {0., 0., 0.};
	double length = 20;
	double dh = 0.5;

	log::cout() << "  >> 2D Cartesian patch" << "\n";

	CartesianPatch *patch_2D = new CartesianPatch(0, 2, origin, length, dh);
	patch_2D->setName("cartesian_uniform_patch_2D");
	patch_2D->update();
	patch_2D->write();

	log::cout() << "\n  >> 2D data" << "\n";

	cellData.resize(patch_2D->getCellCount());
	for (long i = 0; i < patch_2D->getCellCount(); ++i) {
		cellData[i] = i;
	}

	patch_2D->writeField("cell_data_2D", "cell_data", VTKLocation::CELL, cellData);

	vertexData.resize(patch_2D->getVertexCount());
	for (long i = 0; i < patch_2D->getVertexCount(); ++i) {
		vertexData[i] = i;
	}

	patch_2D->writeField("vertex_data_2D", "vertex_data", VTKLocation::POINT, vertexData);

	log::cout() << std::endl;
	log::cout() << "\n  >> 2D bounding box" << "\n";
	log::cout() << std::endl;

	patch_2D->getBoundingBox(minPoint, maxPoint);

	log::cout() << "  >> Minimum x : " << minPoint[0] << std::endl;
	log::cout() << "  >> Minimum y : " << minPoint[1] << std::endl;
	log::cout() << "  >> Minimum z : " << minPoint[2] << std::endl;

	log::cout() << "  >> Maximum x : " << maxPoint[0] << std::endl;
	log::cout() << "  >> Maximum y : " << maxPoint[1] << std::endl;
	log::cout() << "  >> Maximum z : " << maxPoint[2] << std::endl;

	log::cout() << std::endl;
	log::cout() << "\n  >> 2D neighbour test" << "\n";

	std::vector<long> neighs_2D;

	long cellId_2D = 7;
	log::cout() << std::endl;
	log::cout() << "Cell id: " << cellId_2D << std::endl << std::endl;

	log::cout() << "Face neighbours (complete list): " << std::endl;
	neighs_2D = patch_2D->findCellFaceNeighs(cellId_2D);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		log::cout() << " - " << neighs_2D[i] << std::endl;
	}

	log::cout() << "Vertex neighbours (complete list): " << std::endl;
	neighs_2D = patch_2D->findCellVertexNeighs(cellId_2D, true);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		log::cout() << " - " << neighs_2D[i] << std::endl;
	}

	log::cout() << "Vertex neighbours (excuding face neighbours): " << std::endl;
	neighs_2D = patch_2D->findCellVertexNeighs(cellId_2D, false);
	for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
		log::cout() << " - " << neighs_2D[i] << std::endl;
	}

	log::cout() << std::endl;

	delete patch_2D;

	log::cout() << "  >> 3D Cartesian patch" << "\n";

	CartesianPatch *patch_3D = new CartesianPatch(0, 3, origin, length, dh);
	patch_3D->setName("cartesian_uniform_patch_3D");
	patch_3D->update();
	patch_3D->write();

	log::cout() << "\n  >> 3D data" << "\n";

	cellData.resize(patch_3D->getCellCount());
	for (long i = 0; i < patch_3D->getCellCount(); ++i) {
		cellData[i] = i;
	}

	patch_3D->writeField("cell_data_3D", "cell_data", VTKLocation::CELL, cellData);

	vertexData.resize(patch_3D->getVertexCount());
	for (long i = 0; i < patch_3D->getVertexCount(); ++i) {
		vertexData[i] = i;
	}

	patch_3D->writeField("vertex_data_3D", "vertex_data", VTKLocation::POINT, vertexData);

	log::cout() << std::endl;
	log::cout() << "\n  >> 3D bounding box" << "\n";
	log::cout() << std::endl;

	patch_3D->getBoundingBox(minPoint, maxPoint);

	log::cout() << "  >> Minimum x : " << minPoint[0] << std::endl;
	log::cout() << "  >> Minimum y : " << minPoint[1] << std::endl;
	log::cout() << "  >> Minimum z : " << minPoint[2] << std::endl;

	log::cout() << "  >> Maximum x : " << maxPoint[0] << std::endl;
	log::cout() << "  >> Maximum y : " << maxPoint[1] << std::endl;
	log::cout() << "  >> Maximum z : " << maxPoint[2] << std::endl;

	log::cout() << std::endl;
	log::cout() << "\n  >> 3D neighbour test" << "\n";

	std::vector<long> neighs_3D;

	long cellId_3D = 13;
	log::cout() << std::endl;
	log::cout() << "Cell id: " << cellId_3D << std::endl << std::endl;

	log::cout() << "Face neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->findCellFaceNeighs(cellId_3D);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		log::cout() << " - " << neighs_3D[i] << std::endl;
	}

	log::cout() << "Edge neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->findCellEdgeNeighs(cellId_3D, true);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		log::cout() << " - " << neighs_3D[i] << std::endl;
	}

	log::cout() << "Edge neighbours (excuding face neighbours): " << std::endl;
	neighs_3D = patch_3D->findCellEdgeNeighs(cellId_3D, false);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		log::cout() << " - " << neighs_3D[i] << std::endl;
	}

	log::cout() << "Vertex neighbours (complete list): " << std::endl;
	neighs_3D = patch_3D->findCellVertexNeighs(cellId_3D, true);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		log::cout() << " - " << neighs_3D[i] << std::endl;
	}

	log::cout() << "Vertex neighbours (excuding face and edge neighbours): " << std::endl;
	neighs_3D = patch_3D->findCellVertexNeighs(cellId_3D, false);
	for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
		log::cout() << " - " << neighs_3D[i] << std::endl;
	}

	log::cout() << std::endl;

	delete patch_3D;

#if ENABLE_MPI==1
	MPI::Finalize();
#endif

}
