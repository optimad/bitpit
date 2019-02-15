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

#include <array>
#include <vector>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/**
    Test for vertex neighbours
*/
void test_vertex_neighs(VolOctree *patch, std::vector<int> &vertexNeighData)
{
    // Group cells according to their level
    std::map<int, std::vector<long>> cellLevelGropus;
    for (const Cell &cell : patch->getCells()) {
        long cellId = cell.getId();
        int level = patch->getCellLevel(cellId);
        cellLevelGropus[level].push_back(cellId);
    }

    // Evaluate the vertex ring for some randome vertices
    int nCellVertices = std::pow(2, patch->getDimension());

    int nLevels = cellLevelGropus.size();
    for (const auto &entry : cellLevelGropus) {
        int level = entry.first;
        if (level < nLevels - 2) {
            continue;
        }

        const std::vector<long> &cellList = entry.second;
        for (long cellId : cellList) {
            int vertex = rand() % nCellVertices;
            std::vector<long> vertexRing = patch->findCellVertexNeighs(cellId, vertex);
            vertexRing.push_back(cellId);

            // Consider only cells that are not part of previously evaluated
            // rings
            bool skipCell = false;
            for (long ringId : vertexRing) {
                long ringFlatIndex = patch->getCells().evalFlatIndex(ringId);
                if (vertexNeighData[ringFlatIndex] != 0) {
                    skipCell = true;;
                    break;
                }
            }

            if (skipCell) {
                continue;
            }

            // Mark the cells that are part of the ring
            for (long ringId : vertexRing) {
                long ringFlatIndex = patch->getCells().evalFlatIndex(ringId);
                vertexNeighData[ringFlatIndex] += (vertex + 1);
                if (ringId == cellId) {
                    vertexNeighData[ringFlatIndex] *= - 1;
                }
            }
        }
    }
}

/**
    Test for edge neighbours
*/
void test_edge_neighs(VolOctree *patch, std::vector<int> &edgeNeighData)
{
    // Group cells according to their level
    std::map<int, std::vector<long>> cellLevelGropus;
    for (const Cell &cell : patch->getCells()) {
        long cellId = cell.getId();
        int level = patch->getCellLevel(cellId);
        cellLevelGropus[level].push_back(cellId);
    }

    // Evaluate the vertex ring for some randome vertices
    int nCellVertices = std::pow(2, patch->getDimension());

    int nLevels = cellLevelGropus.size();
    for (const auto &entry : cellLevelGropus) {
        int level = entry.first;
        if (level < nLevels - 2) {
            continue;
        }

        const std::vector<long> &cellList = entry.second;
        for (long cellId : cellList) {
            int edge = rand() % nCellVertices;
            std::vector<long> edgeRing = patch->findCellEdgeNeighs(cellId, edge);
            edgeRing.push_back(cellId);

            // Consider only cells that are not part of previously evaluated
            // rings
            bool skipCell = false;
            for (long ringId : edgeRing) {
                long ringFlatIndex = patch->getCells().evalFlatIndex(ringId);
                if (edgeNeighData[ringFlatIndex] != 0) {
                    skipCell = true;;
                    break;
                }
            }

            if (skipCell) {
                continue;
            }

            // Mark the cells that are part of the ring
            for (long ringId : edgeRing) {
                long ringFlatIndex = patch->getCells().evalFlatIndex(ringId);
                edgeNeighData[ringFlatIndex] += (edge + 1);
                if (ringId == cellId) {
                    edgeNeighData[ringFlatIndex] *= - 1;
                }
            }
        }
    }
}

/*!
* Subtest 001
*
* Testing adaption of a 2D patch.
*
* \param origin is the origin of the domain
* \param length is the length of the domain
* \param dh is the maximum cell size
*/
int subtest_001(const std::array<double, 3> &origin, double length, double dh)
{
	log::cout() << std::endl;
	log::cout() << "  :: 2D adaption test ::" << std::endl;

	log::cout() << std::endl;
	log::cout() << ">> Creating the patch" << std::endl;

	VolOctree *patch = new VolOctree(2, origin, length, dh);
	patch->getVTK().setName("octree_adapted_patch_2D");
	patch->buildInterfaces();
	patch->update();
	patch->write();

	log::cout() << std::endl;
	log::cout() << ">> Adapting patch" << std::endl;

	for (int k = 0; k < 5; ++k) {
		long nCells = patch->getCellCount();
		log::cout() << std::endl;
		log::cout() << ">> Marking the cells to adapt... " << std::endl;

		for (int i = 0; i < 150; ++i) {
			long cellId = rand() % nCells;
			if (!patch->getCells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch->findCellNeighs(cellId)) {
				for (auto coarseId : patch->findCellNeighs(neighId)) {
					patch->markCellForCoarsening(coarseId);
				}
			}
		}

		for (int i = 0; i < 200; ++i) {
			long cellId = rand() % nCells;
			if (!patch->getCells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch->findCellNeighs(cellId)) {
				patch->markCellForRefinement(neighId);
			}
		}

		log::cout() << std::endl;
		log::cout() << ">> Initial number of cells... " << nCells << std::endl;

		patch->update();

		nCells = patch->getCellCount();
		log::cout() << ">> Final number of cells... " << nCells << std::endl;
	}
	patch->write();

	log::cout() << std::endl;
	log::cout() << ">> Test neighbour search" << std::endl;

	long nCells = patch->getCellCount();
	std::vector<int> vertexNeighData(nCells, 0);

	patch->getVTK().addData("vertexNeighData", VTKFieldType::SCALAR, VTKLocation::CELL, vertexNeighData);

	test_vertex_neighs(patch, vertexNeighData);

	patch->write();

	return 0;
}

/*!
* Subtest 002
*
* Testing adaption of a 3D patch.
*
* \param origin is the origin of the domain
* \param length is the length of the domain
* \param dh is the maximum cell size
*/
int subtest_002(const std::array<double, 3> &origin, double length, double dh)
{
	log::cout() << std::endl;
	log::cout() << "  :: 3D adaption test ::" << std::endl;

	log::cout() << std::endl;
	log::cout() << ">> Creating patch" << std::endl;

	VolOctree *patch = new VolOctree(3, origin, length, dh);
	patch->getVTK().setName("octree_adapted_patch_3D");
	patch->buildInterfaces();
	patch->update();
	patch->write();

	log::cout() << std::endl;
	log::cout() << ">> Adapting patch" << std::endl;

	for (int k = 0; k < 5; ++k) {
		long nCells = patch->getCellCount();
		log::cout() << std::endl;
		log::cout() << ">> Marking the cells to adapt... " << std::endl;

		for (int i = 0; i < 150; ++i) {
			long cellId = rand() % nCells;
			if (!patch->getCells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch->findCellNeighs(cellId)) {
				for (auto coarseId : patch->findCellNeighs(neighId)) {
					patch->markCellForCoarsening(coarseId);
				}
			}
		}

		for (int i = 0; i < 200; ++i) {
			long cellId = rand() % nCells;
			if (!patch->getCells().exists(cellId)) {
				continue;
			}

			for (auto neighId : patch->findCellNeighs(cellId)) {
				patch->markCellForRefinement(neighId);
			}
		}

		log::cout() << std::endl;
		log::cout() << ">> Initial number of cells... " << nCells << std::endl;

		patch->update();

		nCells = patch->getCellCount();
		log::cout() << ">> Final number of cells... " << nCells << std::endl;
	}
	patch->write();

	log::cout() << std::endl;
	log::cout() << ">> Test neighbour search" << std::endl;

	long nCells = patch->getCellCount();
	std::vector<int> vertexNeighData(nCells, 0);
	std::vector<int> edgeNeighData(nCells, 0);

	patch->getVTK().addData("vertexNeighData", VTKFieldType::SCALAR, VTKLocation::CELL, vertexNeighData);
	patch->getVTK().addData("edgeNeighData", VTKFieldType::SCALAR, VTKLocation::CELL, edgeNeighData);

	test_vertex_neighs(patch, vertexNeighData);
	test_edge_neighs(patch, edgeNeighData);

	patch->write();

	return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

	// Initialize the logger
	log::manager().initialize(log::COMBINED);

	// Seed the random function
	std::srand(1);

	// Run the subtests
	log::cout() << "Testing adaption of octree patches" << std::endl;

	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 20;
	double dh = 1.0;

	int status;
	try {
		status = subtest_001(origin, length, dh);
		if (status != 0) {
			return status;
		}

		status = subtest_002(origin, length, dh);
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		log::cout() << exception.what();
		exit(1);
	}

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif
}
