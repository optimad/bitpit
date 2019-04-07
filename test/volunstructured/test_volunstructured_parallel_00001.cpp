/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
#include "bitpit_volunstructured.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing parallel dump/restore of a 2D unstructured patch.
*
* \param rank is the rank of the process
* \param patch_2D is the patch that will be created by the test
* \param patch_2D_restored is the patch that will be restored by the test
*/
int subtest_001(int rank, VolUnstructured *patch_2D, VolUnstructured *patch_2D_restored)
{
    int archiveVersion = 1;

    // Create the patch
    log::cout() << "Creating 2D patch..." << std::endl;

    patch_2D = new VolUnstructured(2);
    patch_2D->setCommunicator(MPI_COMM_WORLD);
    patch_2D->getVTK().setName("unstructured_patch_2D");

    // Fill the patch
    if (rank == 0) {
        patch_2D->addVertex({{0.00000000, 0.00000000, 0.00000000}},  1);
        patch_2D->addVertex({{0.00000000, 1.00000000, 0.00000000}},  2);
        patch_2D->addVertex({{1.00000000, 1.00000000, 0.00000000}},  3);
        patch_2D->addVertex({{1.00000000, 0.00000000, 0.00000000}},  4);
        patch_2D->addVertex({{1.00000000, 0.50000000, 0.00000000}},  5);
        patch_2D->addVertex({{0.25992107, 1.00000000, 0.00000000}},  6);
        patch_2D->addVertex({{0.58740113, 1.00000000, 0.00000000}},  7);
        patch_2D->addVertex({{0.00000000, 0.75000000, 0.00000000}},  8);
        patch_2D->addVertex({{0.00000000, 0.50000000, 0.00000000}},  9);
        patch_2D->addVertex({{0.00000000, 0.25000000, 0.00000000}}, 10);
        patch_2D->addVertex({{0.25992107, 0.00000000, 0.00000000}}, 11);
        patch_2D->addVertex({{0.58740113, 0.00000000, 0.00000000}}, 12);
        patch_2D->addVertex({{0.42807699, 0.41426491, 0.00000000}}, 13);
        patch_2D->addVertex({{0.30507278, 0.69963441, 0.00000000}}, 14);
        patch_2D->addVertex({{0.64032722, 0.68239464, 0.00000000}}, 15);
        patch_2D->addVertex({{0.24229808, 0.24179558, 0.00000000}}, 16);
        patch_2D->addVertex({{0.67991107, 0.28835559, 0.00000000}}, 17);
        patch_2D->addVertex({{0.22034760, 0.48841527, 0.00000000}}, 18);
        patch_2D->addVertex({{0.43952167, 0.18888322, 0.00000000}}, 19);

        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 3,  7, 15}}));
        patch_2D->addCell(ElementType::QUAD,     std::vector<long>({{ 1, 11, 16, 10}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 8,  9, 18}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 8, 18, 14}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 3, 15,  5}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 9, 10, 18}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{10, 16, 18}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 4, 17, 12}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 4,  5, 17}}));
        patch_2D->addCell(ElementType::QUAD,     std::vector<long>({{13, 17, 15, 14}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{11, 12, 19}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 19, 17}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{12, 17, 19}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 14, 18}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 5, 15, 17}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{11, 19, 16}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 16, 19}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 18, 16}}));
        patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{14, 15,  7}}));
        patch_2D->addCell(ElementType::POLYGON,  std::vector<long>({{ 5,  2,  8, 14, 7, 6}}));
    }

    patch_2D->buildAdjacencies();
    patch_2D->buildInterfaces();

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
    log::cout() << "Internal count:   " << patch_2D->getInternalCount() << std::endl;
    log::cout() << "Ghost count:   " << patch_2D->getGhostCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->write();

    // Create data associated to the patch
    PiercedStorage<bool> booleanStorage(1, &patch_2D->getCells());
    booleanStorage.fill(true);
    booleanStorage.rawAt(0) = false;
    booleanStorage.rawAt(1) = false;
    booleanStorage.rawAt(3) = false;

    PiercedStorage<double> doubleStorage(5, &patch_2D->getCells());
    doubleStorage.fill(100);
    doubleStorage.rawAt(0, 0) = 11;
    doubleStorage.rawAt(0, 1) = 12;
    doubleStorage.rawAt(0, 2) = 13;
    doubleStorage.rawAt(0, 3) = 14;
    doubleStorage.rawAt(0, 4) = 15;
    doubleStorage.rawAt(2, 0) = 31;
    doubleStorage.rawAt(2, 1) = 32;
    doubleStorage.rawAt(2, 2) = 33;
    doubleStorage.rawAt(2, 3) = 34;
    doubleStorage.rawAt(2, 4) = 35;

    PiercedStorage<std::array<double, 3>> arrayDoubleStorage(5, &patch_2D->getCells());
    arrayDoubleStorage.fill({{1., 2., 3.}});
    arrayDoubleStorage.rawAt(0, 0) = {{11., 12., 13.}};

    std::string stringValue = "TEST";

    std::vector<std::string> stringVector(3);
    stringVector[0] = "Test V 0";
    stringVector[1] = "Test V 1";
    stringVector[2] = "Test V 2";

    std::array<std::string, 3> stringArray;
    stringArray[0] = "Test A 0";
    stringArray[1] = "Test A 1";
    stringArray[2] = "Test A 2";

    // Dump the patch
    log::cout() << "Dumping 2D patch..." << std::endl;

    std::string header2D = "2D unstructured patch";
    OBinaryArchive binaryWriter2D("unstructured_patch_2D", archiveVersion, header2D, rank);
    patch_2D->dump(binaryWriter2D.getStream());

    // Dump the data
    booleanStorage.dump(binaryWriter2D.getStream());
    doubleStorage.dump(binaryWriter2D.getStream());
    arrayDoubleStorage.dump(binaryWriter2D.getStream());

    utils::binary::write(binaryWriter2D.getStream(), stringValue);
    utils::binary::write(binaryWriter2D.getStream(), stringVector);
    utils::binary::write(binaryWriter2D.getStream(), stringArray);

    binaryWriter2D.close();

    // Reset the data
    booleanStorage.fill(false);
    doubleStorage.fill(0);
    arrayDoubleStorage.fill({{0., 0., 0.}});
    stringValue = "";
    stringVector.clear();
    stringArray.fill("");

    // Delete the patch
    log::cout() << "Deleting 2D patch..." << std::endl;

    delete patch_2D;

    // Restore the patch
    log::cout() << "Restoring 2D patch..." << std::endl;

    patch_2D_restored = new VolUnstructured();
    patch_2D_restored->setCommunicator(MPI_COMM_WORLD);
    IBinaryArchive binaryReader2D("unstructured_patch_2D", rank);
    patch_2D_restored->restore(binaryReader2D.getStream());

    log::cout() << "Restored cell count:   " << patch_2D_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count: " << patch_2D_restored->getVertexCount() << std::endl;

    // Restore the data
    booleanStorage.setStaticKernel(&patch_2D_restored->getCells());
    doubleStorage.setStaticKernel(&patch_2D_restored->getCells());
    arrayDoubleStorage.setStaticKernel(&patch_2D_restored->getCells());

    booleanStorage.restore(binaryReader2D.getStream());
    doubleStorage.restore(binaryReader2D.getStream());
    arrayDoubleStorage.restore(binaryReader2D.getStream());

    utils::binary::read(binaryReader2D.getStream(), stringValue);
    utils::binary::read(binaryReader2D.getStream(), stringVector);
    utils::binary::read(binaryReader2D.getStream(), stringArray);

    binaryReader2D.close();

    std::cout << "booleanStorage[0] = " << booleanStorage.rawAt(0) << std::endl;
    std::cout << "booleanStorage[1] = " << booleanStorage.rawAt(1) << std::endl;
    std::cout << "booleanStorage[2] = " << booleanStorage.rawAt(2) << std::endl;
    std::cout << "booleanStorage[3] = " << booleanStorage.rawAt(3) << std::endl;
    std::cout << "booleanStorage[4] = " << booleanStorage.rawAt(4) << std::endl;
    std::cout << std::endl;
    std::cout << "doubleStorage[0][0] = " << doubleStorage.rawAt(0, 0) << std::endl;
    std::cout << "doubleStorage[0][1] = " << doubleStorage.rawAt(0, 1) << std::endl;
    std::cout << "doubleStorage[0][2] = " << doubleStorage.rawAt(0, 2) << std::endl;
    std::cout << "doubleStorage[0][3] = " << doubleStorage.rawAt(0, 3) << std::endl;
    std::cout << "doubleStorage[0][4] = " << doubleStorage.rawAt(0, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "doubleStorage[1][0] = " << doubleStorage.rawAt(1, 0) << std::endl;
    std::cout << "doubleStorage[1][1] = " << doubleStorage.rawAt(1, 1) << std::endl;
    std::cout << "doubleStorage[1][2] = " << doubleStorage.rawAt(1, 2) << std::endl;
    std::cout << "doubleStorage[1][3] = " << doubleStorage.rawAt(1, 3) << std::endl;
    std::cout << "doubleStorage[1][4] = " << doubleStorage.rawAt(1, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "doubleStorage[2][0] = " << doubleStorage.rawAt(2, 0) << std::endl;
    std::cout << "doubleStorage[2][1] = " << doubleStorage.rawAt(2, 1) << std::endl;
    std::cout << "doubleStorage[2][2] = " << doubleStorage.rawAt(2, 2) << std::endl;
    std::cout << "doubleStorage[2][3] = " << doubleStorage.rawAt(2, 3) << std::endl;
    std::cout << "doubleStorage[2][4] = " << doubleStorage.rawAt(2, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "arrayDoubleStorage[0][0] = " << arrayDoubleStorage.rawAt(0, 0) << std::endl;
    std::cout << "arrayDoubleStorage[0][1] = " << arrayDoubleStorage.rawAt(0, 1) << std::endl;
    std::cout << "arrayDoubleStorage[0][2] = " << arrayDoubleStorage.rawAt(0, 2) << std::endl;
    std::cout << "arrayDoubleStorage[0][3] = " << arrayDoubleStorage.rawAt(0, 3) << std::endl;
    std::cout << "arrayDoubleStorage[0][4] = " << arrayDoubleStorage.rawAt(0, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "string = " << stringValue << std::endl;
    std::cout << std::endl;
    std::cout << "stringVector[0] = " << stringVector[0] << std::endl;
    std::cout << "stringVector[1] = " << stringVector[1] << std::endl;
    std::cout << "stringVector[2] = " << stringVector[2] << std::endl;
    std::cout << std::endl;
    std::cout << "stringArray[0] = " << stringArray[0] << std::endl;
    std::cout << "stringArray[1] = " << stringArray[1] << std::endl;
    std::cout << "stringArray[2] = " << stringArray[2] << std::endl;

    patch_2D_restored->getVTK().setName("unstructured_patch_2D_restored");
    patch_2D_restored->write();

    // Serialize the patch
    cellRanks.clear();
    cellRanks.resize(patch_2D_restored->getInternalCount(), 0);
    patch_2D_restored->partition(cellRanks, true);

    log::cout() << "Restored serialized cell count:   " << patch_2D_restored->getCellCount() << std::endl;
    log::cout() << "Restored serialized vertex count: " << patch_2D_restored->getVertexCount() << std::endl;

    patch_2D_restored->getVTK().setName("unstructured_patch_2D_restored_serialized");
    patch_2D_restored->write();

    return 0;
}

/*!
* Subtest 002
*
* Testing parallel dump/restore of a 3D unstructured patch.
*
* \param rank is the rank of the process
* \param patch_3D is the patch that will be created by the test
* \param patch_3D_restored is the patch that will be restored by the test
*/
int subtest_002(int rank, VolUnstructured *patch_3D, VolUnstructured *patch_3D_restored)
{
    int archiveVersion = 1;

    log::cout() << "  >> 3D unstructured patch" << std::endl;

    // Create the patch
    log::cout() << "\n\n:: 3D unstructured mesh ::\n";

    patch_3D = new VolUnstructured(3);
    patch_3D->setCommunicator(MPI_COMM_WORLD);
    patch_3D->getVTK().setName("unstructured_patch_3D");

    // Fill the patch
    if (rank == 0) {
        patch_3D->addVertex({{0.00000000, 0.00000000,  0.00000000}},  1);
        patch_3D->addVertex({{1.00000000, 0.00000000,  0.00000000}},  2);
        patch_3D->addVertex({{1.00000000, 1.00000000,  0.00000000}},  3);
        patch_3D->addVertex({{0.00000000, 1.00000000, -0.75000000}},  4);
        patch_3D->addVertex({{0.00000000, 0.00000000,  1.00000000}},  5);
        patch_3D->addVertex({{1.00000000, 0.00000000,  1.00000000}},  6);
        patch_3D->addVertex({{1.00000000, 1.00000000,  1.00000000}},  7);
        patch_3D->addVertex({{0.00000000, 1.00000000,  1.00000000}},  8);
        patch_3D->addVertex({{1.00000000, 0.00000000,  0.54678323}},  9);
        patch_3D->addVertex({{0.50000000, 0.00000000,  1.00000000}}, 10);
        patch_3D->addVertex({{0.00000000, 0.00000000,  0.54678323}}, 11);
        patch_3D->addVertex({{1.00000000, 1.00000000,  0.54678323}}, 12);
        patch_3D->addVertex({{0.50000000, 1.00000000,  1.00000000}}, 13);
        patch_3D->addVertex({{0.00000000, 1.00000000,  0.54678323}}, 14);
        patch_3D->addVertex({{1.00000000, 0.50000000,  1.00000000}}, 15);
        patch_3D->addVertex({{0.00000000, 0.50000000,  1.00000000}}, 16);
        patch_3D->addVertex({{0.51053620, 0.00000000,  0.34680184}}, 17);
        patch_3D->addVertex({{0.36278402, 0.00000000,  0.68603230}}, 18);
        patch_3D->addVertex({{0.69618860, 0.00000000,  0.73234294}}, 19);
        patch_3D->addVertex({{0.51053620, 1.00000000,  0.34680184}}, 20);
        patch_3D->addVertex({{0.36278402, 1.00000000,  0.68603230}}, 21);
        patch_3D->addVertex({{0.69618860, 1.00000000,  0.73234294}}, 22);
        patch_3D->addVertex({{1.00000000, 0.51053620,  0.34680184}}, 23);
        patch_3D->addVertex({{1.00000000, 0.36278402,  0.68603230}}, 24);
        patch_3D->addVertex({{1.00000000, 0.69618860,  0.73234294}}, 25);
        patch_3D->addVertex({{0.00000000, 0.51053620,  0.34680184}}, 26);
        patch_3D->addVertex({{0.00000000, 0.36278402,  0.68603230}}, 27);
        patch_3D->addVertex({{0.00000000, 0.69618860,  0.73234294}}, 28);
        patch_3D->addVertex({{0.50000000, 0.50000000,  1.00000000}}, 29);
        patch_3D->addVertex({{0.75000000, 0.25000000,  1.00000000}}, 30);
        patch_3D->addVertex({{0.25000000, 0.25000000,  1.00000000}}, 31);
        patch_3D->addVertex({{0.00000000, 0.00000000, -1.50000000}}, 32);
        patch_3D->addVertex({{1.00000000, 0.00000000, -1.00000000}}, 33);
        patch_3D->addVertex({{1.00000000, 1.00000000, -0.50000000}}, 34);
        patch_3D->addVertex({{0.00000000, 1.00000000, -1.00000000}}, 35);
        patch_3D->addVertex({{0.00000000, 0.00000000, -2.25000000}}, 36);
        patch_3D->addVertex({{1.00000000, 0.00000000, -1.66666666}}, 37);
        patch_3D->addVertex({{1.00000000, 1.00000000, -2.12500000}}, 38);
        patch_3D->addVertex({{0.00000000, 1.00000000, -2.00000000}}, 39);
        patch_3D->addVertex({{0.00000000, 0.00000000, -3.00000000}}, 40);
        patch_3D->addVertex({{1.00000000, 0.00000000, -3.00000000}}, 41);
        patch_3D->addVertex({{1.00000000, 1.00000000, -3.00000000}}, 42);
        patch_3D->addVertex({{0.00000000, 1.00000000, -3.00000000}}, 43);
        patch_3D->addVertex({{0.00000000, 0.00000000, -4.00000000}}, 44);
        patch_3D->addVertex({{1.00000000, 0.00000000, -4.00000000}}, 45);
        patch_3D->addVertex({{1.00000000, 1.00000000, -4.00000000}}, 46);
        patch_3D->addVertex({{0.00000000, 1.00000000, -4.00000000}}, 47);
        patch_3D->addVertex({{0.50000000, 0.00000000, -4.00000000}}, 48);
        patch_3D->addVertex({{1.00000000, 0.50000000, -4.00000000}}, 49);
        patch_3D->addVertex({{0.50000000, 1.00000000, -4.00000000}}, 50);
        patch_3D->addVertex({{0.33333333, 0.66666666, -4.00000000}}, 51);
        patch_3D->addVertex({{0.00000000, 0.50000000, -4.00000000}}, 52);

        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{29, 22, 25, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 27, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 28, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 26, 23, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 29, 25, 23}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{29, 21, 22, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 28, 27, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23, 26, 29, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 23, 29, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 29, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 29, 23, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 29, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{29, 25, 23, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 23,  3, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{31, 18, 27, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{22, 12, 25, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 11, 27}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 9, 19, 24, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 14, 28}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{15, 22,  7, 25}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{15, 29, 22, 25}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{25, 12, 23, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23,  4,  3,  2}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 18, 11, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 9, 24, 23, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{21, 26, 14, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{28, 21, 13, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 18, 19, 10}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 31, 29, 18}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 31, 18, 10}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{21,  8, 28, 13}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 13, 29, 28}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 13, 28,  8}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{13, 15, 22,  7}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{13, 15, 29, 22}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{15, 24, 29, 25}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{21, 13, 29, 22}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{28, 16, 27, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{31, 18,  5, 27}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 30, 15, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 31, 27, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 11,  5, 27}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 19, 24,  6}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 1, 26, 11, 17}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{19,  9, 24,  6}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{14, 21,  8, 28}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 31,  5, 10}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23, 12,  3, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{19, 30, 10,  6}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 7, 22, 12, 25}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 30,  6, 15}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 31,  5, 27}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23,  9, 17,  2}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26,  4, 14, 20}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 17, 24, 19}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 17, 18, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 24, 30, 19}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 24, 18, 29}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 2, 17, 26,  1}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 2, 26, 17, 23}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4,  2, 26,  1}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 26,  2, 23}}));
        patch_3D->addCell(ElementType::WEDGE,      std::vector<long>({{ 2, 4, 1, 33, 35, 32}}));
        patch_3D->addCell(ElementType::WEDGE,      std::vector<long>({{ 4, 2, 3, 35, 33, 34}}));
        patch_3D->addCell(ElementType::PYRAMID,    std::vector<long>({{36, 37, 38, 39, 33}}));
        patch_3D->addCell(ElementType::PYRAMID,    std::vector<long>({{39, 38, 34, 35, 33}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{35, 36, 39, 33}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{32, 35, 33, 36}}));
        patch_3D->addCell(ElementType::HEXAHEDRON, std::vector<long>({{42, 43, 40, 41, 38, 39, 36, 37}}));
        patch_3D->addCell(ElementType::POLYHEDRON, std::vector<long>({{11,
                                                                            4, 42, 43, 40, 41,
                                                                            5, 52, 51, 50, 49, 48,
                                                                            3, 41, 40, 48,
                                                                            3, 42, 41, 49,
                                                                            3, 43, 42, 50,
                                                                            3, 40, 43, 52,
                                                                            3, 43, 51, 52,
                                                                            3, 43, 50, 51,
                                                                            3, 42, 49, 50,
                                                                            3, 41, 48, 49,
                                                                            3, 40, 52, 48
                                                                        }}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{52, 51, 47, 43}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{51, 50, 47, 43}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{44, 48, 52, 40}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{48, 45, 49, 41}}));
        patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{50, 49, 46, 42}}));
    }

    patch_3D->buildAdjacencies();
    patch_3D->buildInterfaces();

    // Partition the patch
    std::vector<int> cellRanks;
    if (rank == 0) {
        // Evaluate the baricenter of the patch
        long nCells = patch_3D->getCellCount();

        std::array<double, 3> baricenter = {{0., 0., 0.}};
        for (const Cell &cell : patch_3D->getCells()) {
            baricenter += patch_3D->evalCellCentroid(cell.getId());
        }
        baricenter = baricenter / ((double) nCells);

        // Evaluate patch partitioning
        int nProcs;
        MPI_Comm_size(patch_3D->getCommunicator(), &nProcs);

        for (const Cell &cell : patch_3D->getCells()) {
            int side_x = (patch_3D->evalCellCentroid(cell.getId())[0] > baricenter[0]) ? 0 : 1;
            int side_y = (patch_3D->evalCellCentroid(cell.getId())[1] > baricenter[1]) ? 0 : 1;
            int side_z = (patch_3D->evalCellCentroid(cell.getId())[2] > baricenter[2]) ? 0 : 1;

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

    patch_3D->partition(cellRanks, true);

    // Show patch info
    log::cout() << "Cell count:   " << patch_3D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_3D->getVertexCount() << std::endl;

    patch_3D->write();

    // Create data associated to the patch
    PiercedStorage<bool> booleanStorage(1, &patch_3D->getCells());
    booleanStorage.fill(true);
    booleanStorage.rawAt(0) = false;
    booleanStorage.rawAt(1) = false;
    booleanStorage.rawAt(3) = false;

    PiercedStorage<double> doubleStorage(5, &patch_3D->getCells());
    doubleStorage.fill(100);
    doubleStorage.rawAt(0, 0) = 11;
    doubleStorage.rawAt(0, 1) = 12;
    doubleStorage.rawAt(0, 2) = 13;
    doubleStorage.rawAt(0, 3) = 14;
    doubleStorage.rawAt(0, 4) = 15;
    doubleStorage.rawAt(2, 0) = 31;
    doubleStorage.rawAt(2, 1) = 32;
    doubleStorage.rawAt(2, 2) = 33;
    doubleStorage.rawAt(2, 3) = 34;
    doubleStorage.rawAt(2, 4) = 35;

    PiercedStorage<std::array<double, 3>> arrayDoubleStorage(5, &patch_3D->getCells());
    arrayDoubleStorage.fill({{1., 2., 3.}});
    arrayDoubleStorage.rawAt(0, 0) = {{11., 12., 13.}};

    std::string stringValue = "TEST";

    std::vector<std::string> stringVector(3);
    stringVector[0] = "Test V 0";
    stringVector[1] = "Test V 1";
    stringVector[2] = "Test V 2";

    std::array<std::string, 3> stringArray;
    stringArray[0] = "Test A 0";
    stringArray[1] = "Test A 1";
    stringArray[2] = "Test A 2";

    // Dump the patch
    log::cout() << "Dumping 3D patch..." << std::endl;

    std::string header3D = "3D unstructured patch";
    OBinaryArchive binaryWriter3D("unstructured_patch_3D", archiveVersion, header3D, rank);
    patch_3D->dump(binaryWriter3D.getStream());

    // Dump the data
    booleanStorage.dump(binaryWriter3D.getStream());
    doubleStorage.dump(binaryWriter3D.getStream());
    arrayDoubleStorage.dump(binaryWriter3D.getStream());

    utils::binary::write(binaryWriter3D.getStream(), stringValue);
    utils::binary::write(binaryWriter3D.getStream(), stringVector);
    utils::binary::write(binaryWriter3D.getStream(), stringArray);

    binaryWriter3D.close();

    // Reset the data
    booleanStorage.fill(false);
    doubleStorage.fill(0);
    arrayDoubleStorage.fill({{0., 0., 0.}});
    stringValue = "";
    stringVector.clear();
    stringArray.fill("");

    // Delete the patch
    log::cout() << "Deleting 3D patch..." << std::endl;

    delete patch_3D;

    // Restore the patch
    log::cout() << "Restoring 3D patch..." << std::endl;

    patch_3D_restored = new VolUnstructured();
    patch_3D_restored->setCommunicator(MPI_COMM_WORLD);
    IBinaryArchive binaryReader3D("unstructured_patch_3D", rank);
    patch_3D_restored->restore(binaryReader3D.getStream());

    log::cout() << "Restored cell count:   " << patch_3D_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count: " << patch_3D_restored->getVertexCount() << std::endl;

    // Restore the data
    booleanStorage.setStaticKernel(&patch_3D_restored->getCells());
    doubleStorage.setStaticKernel(&patch_3D_restored->getCells());
    arrayDoubleStorage.setStaticKernel(&patch_3D_restored->getCells());

    booleanStorage.restore(binaryReader3D.getStream());
    doubleStorage.restore(binaryReader3D.getStream());
    arrayDoubleStorage.restore(binaryReader3D.getStream());

    utils::binary::read(binaryReader3D.getStream(), stringValue);
    utils::binary::read(binaryReader3D.getStream(), stringVector);
    utils::binary::read(binaryReader3D.getStream(), stringArray);

    binaryReader3D.close();

    std::cout << "booleanStorage[0] = " << booleanStorage.rawAt(0) << std::endl;
    std::cout << "booleanStorage[1] = " << booleanStorage.rawAt(1) << std::endl;
    std::cout << "booleanStorage[2] = " << booleanStorage.rawAt(2) << std::endl;
    std::cout << "booleanStorage[3] = " << booleanStorage.rawAt(3) << std::endl;
    std::cout << "booleanStorage[4] = " << booleanStorage.rawAt(4) << std::endl;
    std::cout << std::endl;
    std::cout << "doubleStorage[0][0] = " << doubleStorage.rawAt(0, 0) << std::endl;
    std::cout << "doubleStorage[0][1] = " << doubleStorage.rawAt(0, 1) << std::endl;
    std::cout << "doubleStorage[0][2] = " << doubleStorage.rawAt(0, 2) << std::endl;
    std::cout << "doubleStorage[0][3] = " << doubleStorage.rawAt(0, 3) << std::endl;
    std::cout << "doubleStorage[0][4] = " << doubleStorage.rawAt(0, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "doubleStorage[1][0] = " << doubleStorage.rawAt(1, 0) << std::endl;
    std::cout << "doubleStorage[1][1] = " << doubleStorage.rawAt(1, 1) << std::endl;
    std::cout << "doubleStorage[1][2] = " << doubleStorage.rawAt(1, 2) << std::endl;
    std::cout << "doubleStorage[1][3] = " << doubleStorage.rawAt(1, 3) << std::endl;
    std::cout << "doubleStorage[1][4] = " << doubleStorage.rawAt(1, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "doubleStorage[2][0] = " << doubleStorage.rawAt(2, 0) << std::endl;
    std::cout << "doubleStorage[2][1] = " << doubleStorage.rawAt(2, 1) << std::endl;
    std::cout << "doubleStorage[2][2] = " << doubleStorage.rawAt(2, 2) << std::endl;
    std::cout << "doubleStorage[2][3] = " << doubleStorage.rawAt(2, 3) << std::endl;
    std::cout << "doubleStorage[2][4] = " << doubleStorage.rawAt(2, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "arrayDoubleStorage[0][0] = " << arrayDoubleStorage.rawAt(0, 0) << std::endl;
    std::cout << "arrayDoubleStorage[0][1] = " << arrayDoubleStorage.rawAt(0, 1) << std::endl;
    std::cout << "arrayDoubleStorage[0][2] = " << arrayDoubleStorage.rawAt(0, 2) << std::endl;
    std::cout << "arrayDoubleStorage[0][3] = " << arrayDoubleStorage.rawAt(0, 3) << std::endl;
    std::cout << "arrayDoubleStorage[0][4] = " << arrayDoubleStorage.rawAt(0, 4) << std::endl;
    std::cout << std::endl;
    std::cout << "string = " << stringValue << std::endl;
    std::cout << std::endl;
    std::cout << "stringVector[0] = " << stringVector[0] << std::endl;
    std::cout << "stringVector[1] = " << stringVector[1] << std::endl;
    std::cout << "stringVector[2] = " << stringVector[2] << std::endl;
    std::cout << std::endl;
    std::cout << "stringArray[0] = " << stringArray[0] << std::endl;
    std::cout << "stringArray[1] = " << stringArray[1] << std::endl;
    std::cout << "stringArray[2] = " << stringArray[2] << std::endl;

    patch_3D_restored->getVTK().setName("unstructured_patch_3D_restored");
    patch_3D_restored->write();

    // Serialize the patch
    cellRanks.clear();
    cellRanks.resize(patch_3D_restored->getInternalCount(), 0);
    patch_3D_restored->partition(cellRanks, true);

    log::cout() << "Restored serialized cell count:   " << patch_3D_restored->getCellCount() << std::endl;
    log::cout() << "Restored serialized vertex count: " << patch_3D_restored->getVertexCount() << std::endl;

    patch_3D_restored->getVTK().setName("unstructured_patch_3D_restored_serialized");
    patch_3D_restored->write();

    return 0;
}

/*!
* Subtest 003
*
* Testing parallel dump/restore through the patch manager.
*
* \param rank is the rank of the process
*/
int subtest_003(int rank)
{
    int archiveVersion = 1;

    // Dump all the patches
    log::cout() << "Dumping patch manager..." << std::endl;

    std::string headerPM = "2D and 3D unstructured patch";
    OBinaryArchive binaryWriterPM("unstructured_patch_PM", archiveVersion, headerPM, rank);
    patch::manager().dumpAll(binaryWriterPM.getStream());
    binaryWriterPM.close();

    // Restore all the patches
    log::cout() << "Restoring patches through patch manager..." << std::endl;

    VolUnstructured *patch_2D_PM_restored = static_cast<VolUnstructured *>(patch::manager().get(0));
    patch_2D_PM_restored->reset();

    VolUnstructured *patch_3D_PM_restored = static_cast<VolUnstructured *>(patch::manager().get(1));
    patch_3D_PM_restored->reset();

    IBinaryArchive binaryReaderPM("unstructured_patch_PM", rank);
    patch::manager().restoreAll(binaryReaderPM.getStream());
    binaryReaderPM.close();

    log::cout() << "Restored cell count (2D):   " << patch_2D_PM_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count (2D): " << patch_2D_PM_restored->getVertexCount() << std::endl;
    log::cout() << "Restored cell count (3D):   " << patch_3D_PM_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count (3D): " << patch_3D_PM_restored->getVertexCount() << std::endl;

    patch_2D_PM_restored->getVTK().setName("unstructured_patch_2D_restored_PM");
    patch_2D_PM_restored->write();

    patch_3D_PM_restored->getVTK().setName("unstructured_patch_3D_restored_PM");
    patch_3D_PM_restored->write();

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
    log::cout() << "Testing dump/restore of unstructured patches" << std::endl;

    int status;
    try {
        VolUnstructured *patch_2D = nullptr;
        VolUnstructured *patch_2D_restored = nullptr;
        VolUnstructured *patch_3D = nullptr;
        VolUnstructured *patch_3D_restored = nullptr;

        status = subtest_001(rank, patch_2D, patch_2D_restored);
        if (status != 0) {
            return status;
        }

        status = subtest_002(rank, patch_3D, patch_3D_restored);
        if (status != 0) {
            return status;
        }

        status = subtest_003(rank);
        if (status != 0) {
            return status;
        }

        delete patch_2D;
        delete patch_2D_restored;
        delete patch_3D_restored;
        delete patch_3D;
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
