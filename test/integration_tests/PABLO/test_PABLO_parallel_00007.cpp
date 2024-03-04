/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"

#include <limits>
#include <mpi.h>

#include <vector>

using namespace bitpit;



/*!
* Subtest 001
*
* Testing point owner detection in two-dimensional tree.
*
* \param rank is the rank of the process
*/
int subtest_001(int rank)
{
    // Instantation of a 2D para_tree object
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 0.;
    double l   = 1;

    PabloUniform octree(x_0, y_0, z_0, l, 2);

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    // Initialize the tree
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    std::vector<uint32_t> refineList2D;
    refineList2D.push_back(  7);
    refineList2D.push_back( 13);
    refineList2D.push_back( 15);
    refineList2D.push_back( 26);
    refineList2D.push_back( 27);
    refineList2D.push_back( 31);
    refineList2D.push_back( 37);
    refineList2D.push_back( 39);
    refineList2D.push_back( 49);
    refineList2D.push_back( 50);
    refineList2D.push_back( 51);
    refineList2D.push_back( 53);
    refineList2D.push_back( 55);
    refineList2D.push_back( 63);
    refineList2D.push_back( 78);
    refineList2D.push_back(100);
    refineList2D.push_back(102);
    refineList2D.push_back(105);
    refineList2D.push_back(108);
    refineList2D.push_back(109);
    refineList2D.push_back(110);
    refineList2D.push_back(135);
    refineList2D.push_back(141);
    refineList2D.push_back(143);
    refineList2D.push_back(146);
    refineList2D.push_back(147);
    refineList2D.push_back(151);
    refineList2D.push_back(153);
    refineList2D.push_back(154);
    refineList2D.push_back(155);
    refineList2D.push_back(157);
    refineList2D.push_back(159);
    refineList2D.push_back(165);
    refineList2D.push_back(167);
    refineList2D.push_back(183);
    refineList2D.push_back(198);
    refineList2D.push_back(204);
    refineList2D.push_back(206);
    refineList2D.push_back(225);
    refineList2D.push_back(228);
    refineList2D.push_back(229);
    refineList2D.push_back(230);

    for (uint32_t id : refineList2D) {
        octree.setMarker(id, 3);
    }
    octree.adapt(false);

    octree.loadBalance();

    octree.computeConnectivity();
    octree.write("Pablo_parallel_00007_2D");

    // Identify point owners
    log::cout() << std::endl;
    log::cout() << " Identifying point owners" << std::endl;
    log::cout() << std::endl;

    std::vector<std::array<double, 3>> points;
    std::vector<uint32_t> expectedOwners;
    std::vector<bool> expectedGhostFlags;
    if (rank == 0) {
        points.push_back({{0., 0., 0.}});
        points.push_back({{10., 20., 0.}});
        points.push_back({{10.5625, 20., 0.}});
        points.push_back({{10.4375, 20.375, 0.}});
        points.push_back({{10.515625, 20.4138, 0.}});
        points.push_back({{10.53125, 20.52, 0.}});
        points.push_back({{10.0, 20.28, 0.}});

        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());
        expectedOwners.push_back(0);
        expectedOwners.push_back(1277);
        expectedOwners.push_back(1180);
        expectedOwners.push_back(41);
        expectedOwners.push_back(69);
        expectedOwners.push_back(569);

        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(true);
        expectedGhostFlags.push_back(true);
        expectedGhostFlags.push_back(false);
    } else if (rank == 1) {
        points.push_back({{0., 0., 0.}});

        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());

        expectedGhostFlags.push_back(false);
    } else if (rank == 2) {
        points.push_back({{0., 0., 0.}});
        points.push_back({{11., 21., 0.}});
        points.push_back({{10.6875, 20.87, 0.}});
        points.push_back({{10.6214, 20.8082, 0.}});
        points.push_back({{10.0625, 20.700, 0.}});
        points.push_back({{10.0804, 20.700, 0.}});

        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());
        expectedOwners.push_back(1366);
        expectedOwners.push_back(1289);
        expectedOwners.push_back(1054);
        expectedOwners.push_back(26);
        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());

        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(true);
        expectedGhostFlags.push_back(true);
    }

    for (std::size_t i = 0; i < points.size(); ++i) {
        log::cout() << " Identifying owner of point " << points[i] << std::endl;

        bool isOwnerGhost;
        uint32_t owner = octree.getPointOwnerIdx(points[i], isOwnerGhost);
        log::cout() << "   Owner: " << owner;
        if (isOwnerGhost) {
            log::cout() << " (Ghost)";
        } else {
            log::cout() << " (Interior)";
        }
        log::cout() << std::endl;
        if (expectedOwners[i] != owner || expectedGhostFlags[i] != isOwnerGhost) {
            log::cout() << "   Expected owner: " << expectedOwners[i];
            if (expectedGhostFlags[i]) {
                log::cout() << " (Ghost)";
            } else {
                log::cout() << " (Interior)";
            }
            log::cout() << std::endl;
            log::cout() << "   Identification failed. " << std::endl;

            throw std::runtime_error("Error in identification of point owner.");
        } else {
            log::cout() << "   Identification completed successfuly. " << std::endl;
        }
    }

    // Done
    return 0;
}

/*!
* Subtest 002
*
* Testing point owner detection in two-dimensional tree.
*
* \param rank is the rank of the process
*/
int subtest_002(int rank)
{
    // Instantation of a 3D para_tree object
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1;

    PabloUniform octree(x_0, y_0, z_0, l, 3);

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    // Initialize the tree
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
    octree.write("Pablo_parallel_00007_3D");

    // Identify point owners
    log::cout() << std::endl;
    log::cout() << " Identifying point owners" << std::endl;
    log::cout() << std::endl;

    std::vector<std::array<double, 3>> points;
    std::vector<uint32_t> expectedOwners;
    std::vector<bool> expectedGhostFlags;
    if (rank == 0) {
        points.push_back({{0., 0., 0.}});
        points.push_back({{10., 20., 30.}});
        points.push_back({{10.5625, 20., 30.5625}});
        points.push_back({{10.4375, 20.375, 31}});
        points.push_back({{10.4375, 20.375, 30.5}});
        points.push_back({{10.3671, 20.368, 30.992}});
        points.push_back({{10.4375, 20.368, 30.992}});

        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());
        expectedOwners.push_back(0);
        expectedOwners.push_back(120);
        expectedOwners.push_back(110);
        expectedOwners.push_back(2265);
        expectedOwners.push_back(9936);
        expectedOwners.push_back(11160);

        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(true);
        expectedGhostFlags.push_back(true);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
    } else if (rank == 1) {
        points.push_back({{0., 0., 0.}});

        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());

        expectedGhostFlags.push_back(false);
    } else if (rank == 2) {
        points.push_back({{0., 0., 0.}});
        points.push_back({{11., 21., 31.}});
        points.push_back({{10.9, 20.75, 30.75}});

        expectedOwners.push_back(std::numeric_limits<uint32_t>::max());
        expectedOwners.push_back(11248);
        expectedOwners.push_back(11109);

        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
        expectedGhostFlags.push_back(false);
    }

    for (std::size_t i = 0; i < points.size(); ++i) {
        log::cout() << " Identifying owner of point " << points[i] << std::endl;

        bool isOwnerGhost;
        uint32_t owner = octree.getPointOwnerIdx(points[i], isOwnerGhost);
        log::cout() << "   Owner: " << owner;
        if (isOwnerGhost) {
            log::cout() << " (Ghost)";
        } else {
            log::cout() << " (Interior)";
        }
        log::cout() << std::endl;
        if (expectedOwners[i] != owner || expectedGhostFlags[i] != isOwnerGhost) {
            log::cout() << "   Expected owner: " << expectedOwners[i];
            if (expectedGhostFlags[i]) {
                log::cout() << " (Ghost)";
            } else {
                log::cout() << " (Interior)";
            }
            log::cout() << std::endl;
            log::cout() << "   Identification failed. " << std::endl;

            throw std::runtime_error("Error in identification of point owner.");
        } else {
            log::cout() << "   Identification completed successfuly. " << std::endl;
        }
    }

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

    log::manager().initialize(log::MODE_COMBINE, true, nProcs, rank);
    log::cout().setDefaultVisibility(log::VISIBILITY_GLOBAL);

    // Run the subtests
    log::cout() << "Testing parallel octree dump and restore." << std::endl;

    int status;
    try {
        status = subtest_001(rank);
        if (status != 0) {
            return status;
        }

        status = subtest_002(rank);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
