// ========================================================================== //
//           ** BitPit mesh ** Test 001 for class surftri_patch **            //
//                                                                            //
// Test construction, modifiers and communicators for surftri_patch.          //
// ========================================================================== //
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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>
# include <iostream>

// BitPit
# include "bitpit_common.hpp"                                                 // Utilities and common definitions
# include "bitpit_operators.hpp"                                              // STL containers operators
# include "bitpit_patch.hpp"                                                  // BitPit base patch
# include "bitpit_surftripatch.hpp"                                           // BitPit surftri patch

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// GENERATE A TEST NON-MANIFOLD SURFACE TRIANGULATION FOR TESTS.              //
// ========================================================================== //
void generateTestTriangulation(
    SurfTriPatch                &mesh
) {

// ========================================================================== //
// void generateTestTriangulation(                                            //
//     SurfTriPatch                &mesh)                                     //
//                                                                            //
// Generate a non-manifold surface triangulation for tests.                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - mesh    : SurfTriPatch, surface mesh patch                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long                    nV = 27;
long                    nS = 33;

// Counters
int                     i, j;

// ========================================================================== //
// INITIALIZE TRIANGULATION                                                   //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //

    // Reserve memory for vertex & cell storage ----------------------------- //
    mesh.reserveVertices(nV);
    mesh.reserveCells(nS);
}

// ========================================================================== //
// GENERATE VERTEX LIST                                                       //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    double                      off = -0.5;
    array<double, 3>            vertex;
    
    // 0-row ---------------------------------------------------------------- //
    vertex.fill(0.0);
    for (i = 0; i < 8; ++i) {
        vertex[0] = double(i);
        mesh.addVertex(vertex);
    } //next i

    // 1-row ---------------------------------------------------------------- //
    vertex.fill(0.0);
    vertex[1] = 1.0;
    for (i = 0; i < 9; ++i) {
        vertex[0] = double(i) + off;
        mesh.addVertex(vertex);
    } //next i

    // 2-row ---------------------------------------------------------------- //
    vertex.fill(0.0);
    vertex[1] = 2.0;
    for (i = 0; i < 8; ++i) {
        vertex[0] = double(i);
        mesh.addVertex(vertex);
    } //next i

    // Orthogonal element(s) ------------------------------------------------ //
    vertex[0] = 0.5*(mesh.getVertex(3)[0] + mesh.getVertex(12)[0]);
    vertex[1] = 0.5*(mesh.getVertex(3)[1] + mesh.getVertex(12)[1]);
    vertex[2] = 0.5 * sqrt(3.0);
    mesh.addVertex(vertex);
    vertex[0] = 0.5*(mesh.getVertex(12)[0] + mesh.getVertex(21)[0]);
    vertex[1] = 0.5*(mesh.getVertex(12)[1] + mesh.getVertex(21)[1]);
    vertex[2] = 0.5 * sqrt(3.0);
    mesh.addVertex(vertex);
}

// ========================================================================== //
// GENERATE CONNECTIVITY                                                      //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    long int                    off;
    vector<long>                connectivity(3);

    // 0-row ---------------------------------------------------------------- //
    off = 8;
    for (i = 0; i < 7; ++i) {
        connectivity[0] = i;
        connectivity[1] = i + 1 + off;
        connectivity[2] = i + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
        connectivity[0] = i;
        connectivity[1] = i + 1;
        connectivity[2] = i + 1 + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    } //next i
    connectivity[0] = i;
    connectivity[1] = i + 1 + off;
    connectivity[2] = i + off;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);

    // 1-row ---------------------------------------------------------------- //
    off = 9;
    for (i = 8; i < 15; ++i) {
        connectivity[0] = i;
        connectivity[1] = i + 1;
        connectivity[2] = i + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
        connectivity[0] = i + 1;
        connectivity[1] = i + 1 + off;
        connectivity[2] = i + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    } //next i
    connectivity[0] = i;
    connectivity[1] = i + 1;
    connectivity[2] = i + off;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);

    // Orthogonal element --------------------------------------------------- //
    connectivity[0] = 3;
    connectivity[1] = 12;
    connectivity[2] = 25;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    connectivity[0] = 12;
    connectivity[1] = 26;
    connectivity[2] = 25;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    connectivity[0] = 12;
    connectivity[1] = 21;
    connectivity[2] = 26;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    
}

return;
}

// ========================================================================== //
// SUBTEST #001 Communications among 3 processes                              //
// ========================================================================== //
void COM_step(
    Patch                       &mesh,
    short                        snd,
    short                        rcv,
    vector<long>                &id_list,
    int                          step_id
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
SurfTriPatch                     envelope(0);

// Counters
// none

// ========================================================================== //
// MOVE CELL                                                                  //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none
    cout << "* rank#" << snd << " send {" << id_list << "} to rank#" << rcv << endl;
    mesh.sendCells(snd, rcv, id_list);

    // Extracting mesh boundaries ------------------------------------------- //
    cout << "* (rank #" << mesh.getRank() << ") extracting external envelope" << endl;
    mesh.extractEnvelope(envelope);

    // Export mesh ---------------------------------------------------------- //
    stringstream                f1, f2;
    f1 << "step" << step_id;
    mesh.write(f1.str());
    cout << "  (rank #" << mesh.getRank() << ", mesh exported to \"" << f1.str() << ".vtu\")" << endl;
    f2 << "P" << utils::zeroPadNumber(6, mesh.getRank()) << "_env_step" << step_id;
    envelope.write(f2.str());
    cout << "  (rank #" << mesh.getRank() << ", mesh external envelope exported to \"" << f2.str() << ".vtu\")" << endl;
}

return; }
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// Sub-test #001 for test #00001 - parallel. Test mesh partitioning on 3      //
// processes.                                                                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
   
// ========================================================================== //
// SCOPE VARIABLES                                                            //
// ========================================================================== //

// Local variables
SurfTriPatch                     mesh(0);

// Counters
// none

// ========================================================================== //
// SET MESH ATTRIBUTES                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set MPI communicator ------------------------------------------------- //
    mesh.setCommunicator(MPI_COMM_WORLD);

}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "** ================================================================= **" << endl;
    cout << "** Test #00001 - sub-test #001 - Testing mesh communications between **" << endl;
    cout << "**                               3 processes.                        **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

// ========================================================================== //
// GENERATE 2d NON-MANIFOLD TRIANGULATION                                     //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    stringstream                name;
    SurfTriPatch                envelope(0);

    // Generate dummy triangulation ----------------------------------------- //
    if (mesh.getRank() == 0) {
        cout << "** (rank #0) Generating 2D non-manifold triangulation" << endl;
        generateTestTriangulation(mesh);
    }
    mesh.buildAdjacencies();

    // Extracting mesh boundaries ------------------------------------------- //
    cout << "* (rank #" << mesh.getRank() << ") extracting external envelope" << endl;
    mesh.extractEnvelope(envelope);

    // Export mesh ---------------------------------------------------------- //
    mesh.write("step0");
    cout << "  (rank #" << mesh.getRank() << ", mesh exported to \"step0.vtu\")" << endl;
    name << "P" << utils::zeroPadNumber(6, mesh.getRank()) <<"_env_step0";
    envelope.write(name.str());
    cout << "  (rank #" << mesh.getRank() << ", mesh external envelope exported to " << name.str() << ".vtu)" << endl;
    cout << endl;
}

// ========================================================================== //
// RANK#0 SENDING CELLS {30,31} TO RANK #2                                    //
// ========================================================================== //
{
    vector<long>                        cell_list{30,31};
    COM_step(mesh, 0, 2, cell_list, 1);
    cout << endl;
}

// ========================================================================== //
// RANK#0 SENDING CELLS {7,8,9,23,22,24} TO RANK #1                           //
// ========================================================================== //
{
    vector<long>                        cell_list{7,8,9,23,22,24};
    COM_step(mesh, 0, 1, cell_list, 2);
    cout << endl;
}

// ========================================================================== //
// RANK#2 SENDING CELLS {0,1} TO RANK #1                                      //
// ========================================================================== //
{
    vector<long>                        cell_list{0,1};
    COM_step(mesh, 2, 1, cell_list, 3);
    cout << endl;
}

// ========================================================================== //
// RANK#1 SENDING CELLS {3,4,5} TO RANK #2                                    //
// ========================================================================== //
{
    vector<long>                        cell_list{3,4,5};
    COM_step(mesh, 1, 2, cell_list, 4);
    cout << endl;
}

// ========================================================================== //
// RANK#1 SENDING CELLS {0,1,2} TO RANK #0                                    //
// ========================================================================== //
{
    vector<long>                        cell_list{0,1,2};
    COM_step(mesh, 1, 0, cell_list, 5);
    cout << endl;
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "** ================================================================= **" << endl;
    cout << "** Test #00001 - sub-test #001 - completed!                          **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

return 0;

}

// ========================================================================== //
// MAIN FOR TEST #00001 - parallel                                            //
// ========================================================================== //
int main(int argc, char* argv[])
{

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variabels
int                             err = 0;

// ========================================================================== //
// INITIALIZE MPI COMMUNICATOR                                                //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Initialize MPI comm world -------------------------------------------- //
    MPI_Init(&argc, &argv);

}

// ========================================================================== //
// RUN SUB-TEST #001                                                          //
// ========================================================================== //
err = subtest_001();
if (err > 0) return(10 + err);

// ========================================================================== //
// FINALIZE MPI COMMUNICATOR                                                  //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Finalize comm world -------------------------------------------------- //
    MPI_Finalize();

}

return err;

}
