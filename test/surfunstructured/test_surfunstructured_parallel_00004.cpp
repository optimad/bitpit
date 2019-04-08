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
#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_surfunstructured.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing dump/restore of a 2D octree patch.
*
* \param rank is the rank of the process
* \param patch_2D is the patch that will be created by the test
* \param patch_2D_restored is the patch that will be restored by the test
*/
int subtest_001(int rank, SurfUnstructured *patch_2D, SurfUnstructured *patch_2D_restored)
{
    int archiveVersion = 1;

    log::cout() << "  >> 2D unstructured surface patch" << std::endl;

    // Create the patch
    log::cout() << "Creating 2D patch..." << std::endl;

    patch_2D = new SurfUnstructured(2, 2);
    patch_2D->setCommunicator(MPI_COMM_WORLD);
    patch_2D->getVTK().setName("surfunstructured_patch_2D");

    // Fill the patch
    if (rank == 0) {
        const std::string fielname_2D = "./data/cube.stl";
        patch_2D->importSTL(fielname_2D);
    }

    // Partition the patch
    std::vector<int> cellRanks;
    if (rank == 0) {
        // Evaluate the baricenter of the patch
        long nCells = patch_2D->getCellCount();

        std::array<double, 3> baricenter = {{0., 0., 0.}};
        for (const Cell &cell : patch_2D->getCells()) {
            baricenter += patch_2D->evalCellCentroid(cell.getId());
        }
        baricenter = baricenter / ((double) nCells);

        // Generate patch partitioning
        int nProcs;
        MPI_Comm_size(patch_2D->getCommunicator(), &nProcs);

        for (const Cell &cell : patch_2D->getCells()) {
            int side_x = (patch_2D->evalCellCentroid(cell.getId())[0] > baricenter[0]) ? 0 : 1;
            int side_y = (patch_2D->evalCellCentroid(cell.getId())[1] > baricenter[1]) ? 0 : 1;
            int side_z = (patch_2D->evalCellCentroid(cell.getId())[2] > baricenter[2]) ? 0 : 1;

            int rank = -1;
            if (side_z == 0 && side_y == 0 && side_x == 0) {
                rank = 0;
            } else if (side_z == 0 && side_y == 0 && side_x == 1) {
                rank = 1;
            } else if (side_z == 0 && side_y == 1 && side_x == 0) {
                rank = 2;
            } else if (side_z == 0 && side_y == 1 && side_x == 1) {
                rank = 3;
            } else if (side_z == 1 && side_y == 0 && side_x == 0) {
                rank = 4;
            } else if (side_z == 1 && side_y == 0 && side_x == 1) {
                rank = 5;
            } else if (side_z == 1 && side_y == 1 && side_x == 0) {
                rank = 6;
            } else if (side_z == 1 && side_y == 1 && side_x == 1) {
                rank = 7;
            }
            rank = rank % nProcs;

            cellRanks.push_back(rank);
        }
    }

    patch_2D->partition(cellRanks, true);

    // Show patch info
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->write();

    // Dump the patch
    log::cout() << "Dumping 2D patch..." << std::endl;

    std::string header2D = "2D unstructured surface patch";
    OBinaryArchive binaryWriter2D("surfunstructured_patch_2D.dat", archiveVersion, header2D, rank);
    patch_2D->dump(binaryWriter2D.getStream());
    binaryWriter2D.close();

    // Delete the patch
    log::cout() << "Deleting 2D patch..." << std::endl;

    // Restore the patch
    log::cout() << "Restoring 2D patch..." << std::endl;

    patch_2D_restored = new SurfUnstructured();
    patch_2D_restored->setCommunicator(MPI_COMM_WORLD);
    IBinaryArchive binaryReader2D("surfunstructured_patch_2D.dat", rank);
    patch_2D_restored->restore(binaryReader2D.getStream());
    binaryReader2D.close();

    log::cout() << "Restored cell count:   " << patch_2D_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count: " << patch_2D_restored->getVertexCount() << std::endl;

    patch_2D_restored->getVTK().setName("surfunstructured_patch_2D_restored");
    patch_2D_restored->write();

    return 0;
}

/*!
* Subtest 002
*
* Testing dump/restore through the patch manager.
*
* \param rank is the rank of the process
*/
int subtest_002(int rank)
{
    int archiveVersion = 1;

    // Dump all the patches
    log::cout() << "Dumping patch manager..." << std::endl;

    std::string headerPM = "2D unstructured surface patch";
    OBinaryArchive binaryWriterPM("surfunstructured_patch_PM.dat", archiveVersion, headerPM, rank);
    patch::manager().dumpAll(binaryWriterPM.getStream());
    binaryWriterPM.close();

    // Restore all the patches
    log::cout() << "Restoring patches through patch manager..." << std::endl;

    SurfUnstructured *patch_2D_PM_restored = static_cast<SurfUnstructured *>(patch::manager().get(0));
    patch_2D_PM_restored->reset();

    IBinaryArchive binaryReaderPM("surfunstructured_patch_PM.dat", rank);
    patch::manager().restoreAll(binaryReaderPM.getStream());
    binaryReaderPM.close();

    log::cout() << "Restored cell count (2D):   " << patch_2D_PM_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count (2D): " << patch_2D_PM_restored->getVertexCount() << std::endl;

    patch_2D_PM_restored->getVTK().setName("surfunstructured_patch_2D_restored_PM");
    patch_2D_PM_restored->write();

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    // Initialize the logger
    int nProcs;
    int    rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Testing dump/restore of unstructured surface patches" << std::endl;

    int status;
    try {
        SurfUnstructured *patch_2D = nullptr;
        SurfUnstructured *patch_2D_restored = nullptr;

        status = subtest_001(rank, patch_2D, patch_2D_restored);
        if (status != 0) {
            return status;
        }

        status = subtest_002(rank);
        if (status != 0) {
            return status;
        }

        delete patch_2D;
        delete patch_2D_restored;
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
