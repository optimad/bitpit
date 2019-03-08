// ========================================================================== //
//           ** BitPit mesh ** Test 003 for class surftri_patch **            //
//                                                                            //
// Test import output routines for class SurfUnstructured                         //
// ========================================================================== //
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
# include "bitpit_operators.hpp"                                              // BitPit operators
# include "bitpit_IO.hpp"                                                     // BitPit IO
# include "bitpit_surfunstructured.hpp"                                           // BitPit surftri patch

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// SUBTEST #001 Test STL import/export ruotines                               //
// ========================================================================== //
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// Test mesh partitioning                                                     //
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
//               err = 2  --> error at step #2 (mesh partitioning)            //
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
    log::cout() << "** Sub-test #001 - Mesh partitioning                                 **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

// ========================================================================== //
// STEP #1 (IMPORT MESH FROM BINARY STL)                                      //
// ========================================================================== //
int myRank = mesh.getRank();
if (myRank == 0) {
    // Scope variables ------------------------------------------------------ //
    long                                nExpected = 283274;

    // Import mesh from stl format ------------------------------------------ //
    log::cout() << "** Importing mesh from (binary): \"" << in_name_bin << "\"" << endl;
    if (mesh.importSTL(in_name_bin, true) > 0) return 1;
    mesh.collapseCoincidentVertices();

    if (mesh.getVertexCount() != 3*nExpected) return 1;
    if (mesh.getInternalCount() != nExpected) return 1;
}

{
    //  Build adjacencies --------------------------------------------------- //
    mesh.buildAdjacencies();
}

// ========================================================================== //
// STEP #2 (PARTITION THE MESH)                                               //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
	long nCells = mesh.getCellCount();

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
			int rank = (mesh.evalCellCentroid(cell.getId())[0] > baricenter[0]) ? 0 : 1;
			cellRanks.push_back(rank);
		}
	}

	mesh.partition(cellRanks, false);

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

return 0; }

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
