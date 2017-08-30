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
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing dump and restore.
*/
int subtest_001()
{
    int archiveVersion = 1;

    // Instantation of a 2D para_tree object
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1.5;

    PabloUniform octree(x_0, y_0, z_0, l, 2);

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    // Refine and write the octree
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    std::vector<uint32_t> refineList;
    refineList.push_back(  7);
    refineList.push_back( 13);
    refineList.push_back( 15);
    refineList.push_back( 26);
    refineList.push_back( 27);
    refineList.push_back( 31);
    refineList.push_back( 37);
    refineList.push_back( 39);
    refineList.push_back( 49);
    refineList.push_back( 50);
    refineList.push_back( 51);
    refineList.push_back( 53);
    refineList.push_back( 55);
    refineList.push_back( 63);
    refineList.push_back( 78);
    refineList.push_back(100);
    refineList.push_back(102);
    refineList.push_back(105);
    refineList.push_back(108);
    refineList.push_back(109);
    refineList.push_back(110);
    refineList.push_back(135);
    refineList.push_back(141);
    refineList.push_back(143);
    refineList.push_back(146);
    refineList.push_back(147);
    refineList.push_back(151);
    refineList.push_back(153);
    refineList.push_back(154);
    refineList.push_back(155);
    refineList.push_back(157);
    refineList.push_back(159);
    refineList.push_back(165);
    refineList.push_back(167);
    refineList.push_back(183);
    refineList.push_back(198);
    refineList.push_back(204);
    refineList.push_back(206);
    refineList.push_back(225);
    refineList.push_back(228);
    refineList.push_back(229);
    refineList.push_back(230);

    for (uint32_t id : refineList) {
        octree.setMarker(id, 3);
    }
    octree.adapt(false);

    octree.computeConnectivity();
    octree.write("Pablo_00006_original");

    // Dump the tree
    std::string header = "2D PABLO";
    OBinaryArchive binaryWriter2D("Pablo_00006_dump", archiveVersion, header);
    octree.dump(binaryWriter2D.getStream());
    binaryWriter2D.close();

    // Create an empty octree
    PabloUniform octreeRestored;

    // Restore the tree
    IBinaryArchive binaryReader2D("Pablo_00006_dump");
    octreeRestored.restore(binaryReader2D.getStream());

    std::cout << " Restored Origin : ( " << octreeRestored.getX0() << ", " << octreeRestored.getY0() << ", " << octreeRestored.getZ0() << " )" << std::endl;
    std::cout << " Restored Length : " << octreeRestored.getL() << std::endl;

    // Write the restored octree
    octreeRestored.computeConnectivity();
    octreeRestored.write("Pablo_00006_restored");

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

    int nProcs;
    int rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcs = 1;
    rank   = 0;
#endif

    // Initialize the logger
    log::manager().initialize(log::SEPARATE, false, nProcs, rank);
    log::cout() << fileVerbosity(log::NORMAL);
    log::cout() << consoleVerbosity(log::QUIET);

    // Run the subtests
    log::cout() << "Testing data smoothing" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
