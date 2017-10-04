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

#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing parallel dump/restore of a 2D octree patch.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    int archiveVersion = 1;

    // Instantation of a 3D para_tree object
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1.5;

    PabloUniform octree(x_0, y_0, z_0, l, 3);

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    // Refine and write the octree
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    std::vector<uint32_t> refineList3D;
    refineList3D.push_back(2351);
    refineList3D.push_back(2365);
    refineList3D.push_back(2367);
    refineList3D.push_back(2422);
    refineList3D.push_back(2423);
    refineList3D.push_back(2431);
    refineList3D.push_back(2477);
    refineList3D.push_back(2479);
    refineList3D.push_back(2533);
    refineList3D.push_back(2534);
    refineList3D.push_back(2535);
    refineList3D.push_back(2541);
    refineList3D.push_back(2543);
    refineList3D.push_back(2559);
    refineList3D.push_back(2878);
    refineList3D.push_back(2988);
    refineList3D.push_back(2990);
    refineList3D.push_back(2997);
    refineList3D.push_back(3004);
    refineList3D.push_back(3005);
    refineList3D.push_back(3006);
    refineList3D.push_back(3375);
    refineList3D.push_back(3389);
    refineList3D.push_back(3391);
    refineList3D.push_back(3430);
    refineList3D.push_back(3431);
    refineList3D.push_back(3439);
    refineList3D.push_back(3445);
    refineList3D.push_back(3446);
    refineList3D.push_back(3447);
    refineList3D.push_back(3453);
    refineList3D.push_back(3455);
    refineList3D.push_back(3501);
    refineList3D.push_back(3503);
    refineList3D.push_back(3567);
    refineList3D.push_back(3886);
    refineList3D.push_back(3900);
    refineList3D.push_back(3902);
    refineList3D.push_back(4005);
    refineList3D.push_back(4012);
    refineList3D.push_back(4013);
    refineList3D.push_back(4014);

    for (uint32_t id : refineList3D) {
        octree.setMarker(id, 3);
    }
    octree.adapt(false);

    octree.loadBalance();

    octree.computeConnectivity();
    octree.write("Pablo_parallel_00003_original");

    // Dump the tree
    std::string header = "3D PABLO";
    OBinaryArchive binaryWriter3D("Pablo_parallel_00003_dump", archiveVersion, header, rank);
    octree.dump(binaryWriter3D.getStream());
    binaryWriter3D.close();

    // Create an empty octree
    PabloUniform octreeRestored;

    // Restore the tree
    IBinaryArchive binaryReader3D("Pablo_parallel_00003_dump", rank);
    octreeRestored.restore(binaryReader3D.getStream());

    octreeRestored.loadBalance();

    std::cout << " Restored Origin : ( " << octreeRestored.getX0() << ", " << octreeRestored.getY0() << ", " << octreeRestored.getZ0() << " )" << std::endl;
    std::cout << " Restored Length : " << octreeRestored.getL() << std::endl;

    // Write the restored octree
    octreeRestored.computeConnectivity();
    octreeRestored.write("Pablo_parallel_00003_restored");

    // Done
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
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);

    // Run the subtests
    log::cout() << "Testing parallel octree dump and restore." << std::endl;

    int status;
    try {
        status = subtest_001(rank);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
