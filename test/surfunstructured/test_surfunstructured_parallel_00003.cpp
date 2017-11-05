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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>
# include <iostream>
# include <mpi.h>

// BitPit
# include "bitpit_operators.hpp"
# include "bitpit_IO.hpp"
# include "bitpit_surfunstructured.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// SUBTEST #001 Test surface orientation                                      //
// ========================================================================== //
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// Test mesh orientation                                                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
//               err = 1  --> error at step #1 (import STL)                   //
//               err = 2  --> error at step #2 (mesh partition)               //
//               err = 3  --> error at step #3 (mesh orientation)             //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
string                          in_name_bin = "./data/buddha.stl";
SurfUnstructured                mesh(2, 3);

// Counters
// none

// ========================================================================== //
// INITIALIZE MESH PARAMETERS                                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set name ------- ----------------------------------------------------- //
    mesh.getVTK().setName("surfunstructured_partition");

    // Set communicator ----------------------------------------------------- //
    mesh.setCommunicator(MPI_COMM_WORLD);
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Sub-test #001 - Orientation surface                               **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

// ========================================================================== //
// STEP #1 (IMPORT MESH FROM BINARY STL AND FLIP EVERY 4TH TRIANGLE)          //
// ========================================================================== //
int myRank = mesh.getRank();
if (myRank == 0) {
    // Scope variables ------------------------------------------------------ //
    // None

    // Import mesh from stl format ------------------------------------------ //
    log::cout() << "** Importing mesh from (binary): \"" << in_name_bin << "\"" << endl;
    if (mesh.importSTL(in_name_bin, true) > 0) return 1;
    mesh.collapseCoincidentVertices();
    mesh.buildAdjacencies();

    // Mess-up the orientation ---------------------------------------------- //
    for (const auto &cell : mesh.getCells()) {
        if (!cell.isInterior()) {
            continue;
        }

        long cellId = cell.getId();
        if (cellId % 4 == 0) {
            mesh.flipCellOrientation(cellId);
        }
    }
}

// ========================================================================== //
// STEP #2 (PARTITION THE MESH)                                               //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    long nCells = mesh.getCellCount();

    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // Evaluation of baricenter ----------------------------------------------//
    std::array<double, 3> baricenter = {{0, 0, 0}};
    for (const auto &cell : mesh.getCells()) {
        baricenter += mesh.evalCellCentroid(cell.getId());
    }
    baricenter = baricenter / ((double) nCells);

    // Partitioning ----------------------------------------------------------//
    log::cout() << "** Mesh partitioning" << endl;

    std::vector<int> cellRanks;
    if (myRank == 0) {
        for (const auto &cell : mesh.getCells()) {
            int side_x = (mesh.evalCellCentroid(cell.getId())[0] > baricenter[0]) ? 0 : 1;
            int side_y = (mesh.evalCellCentroid(cell.getId())[1] > baricenter[1]) ? 0 : 1;
            int side_z = (mesh.evalCellCentroid(cell.getId())[2] > baricenter[2]) ? 0 : 1;

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

    mesh.partition(cellRanks, false);
    mesh.buildInterfaces();

}

// ========================================================================== //
// STEP #3 (ADJUST ORIENTATION)                                               //
// ========================================================================== //
{
    // Adjust orientation --------------------------------------------------- //
    log::cout() << "** Mesh adjust orientation" << endl;
    bool orientable = mesh.adjustCellOrientation();
    if (orientable) {
        log::cout() << "   Mesh successfully oriented" << endl;
    } else {
        log::cout() << "   Error during surface orientation" << endl;
    }

    // Write mesh ----------------------------------------------------------- //
    log::cout() << "** Writing mesh" << endl;
    mesh.write();

    // Output message ------------------------------------------------------- //
    log::cout() << "** Mesh partitioning" << endl;
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Sub-test #001 - completed!                                        **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

return 0;

}

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    //======================================================================= //
    MPI_Init(&argc,&argv);

    // ====================================================================== //
	// Initialize the logger
    // ====================================================================== //
	int nProcs;
	int	rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log::manager().initialize(log::COMBINED, true, nProcs, rank);
	log::cout().setVisibility(log::GLOBAL);

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    int                             status = 0;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
    MPI_Finalize();

    return status;
}
