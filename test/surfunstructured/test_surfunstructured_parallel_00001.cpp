// ========================================================================== //
//           ** BitPit mesh ** Test 001 for class surftri_patch **            //
//                                                                            //
// Test construction, modifiers and communicators for SurfUnstructured        //
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
# include <ctime>
# include <chrono>

// BitPit
# include "bitpit_common.hpp"                                                 // Utilities and common definitions
# include "bitpit_IO.hpp"                                                     // Input/output
# include "bitpit_operators.hpp"                                              // STL containers operators
# include "bitpit_patchkernel.hpp"                                                  // BitPit base patch
# include "bitpit_surfunstructured.hpp"                                           // BitPit surftri patch

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;
using namespace chrono;

// ========================================================================== //
// GENERATE A TEST NON-MANIFOLD SURFACE TRIANGULATION FOR TESTS.              //
// ========================================================================== //
void generateTestTriangulation(
    SurfUnstructured                &mesh
) {

// ========================================================================== //
// void generateTestTriangulation(                                            //
//     SurfUnstructured                &mesh)                                 //
//                                                                            //
// Generate a non-manifold surface triangulation for tests.                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - mesh    : SurfUnstructured, surface mesh patch                           //
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
int                     i;

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
void generateTestQuadMesh(
    long                        nx,
    long                        ny,
    SurfUnstructured                &mesh
) {

// ========================================================================== //
// void generateTestQuadMesh(                                                 //
//     long                        nx,                                        //
//     long                        ny,                                        //
//     SurfUnstructured                &mesh)                                 //
//                                                                            //
// Generate surface triangulation for sub-test #002.                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - mesh       : SurfUnstructured, mesh container                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
array<double, 3>                        xlim, ylim;

// Counters
long                                    i, j;

// ========================================================================== //
// INITIALIZE MESH                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set mesh parameters -------------------------------------------------- //
    xlim[0] = -1.;      xlim[1] = 1.;
    ylim[0] = -1.;      ylim[1] = 1.;
}

// ========================================================================== //
// GENERATE MESH                                                              //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    double                              dx, dy;
    array<double, 3>                    point;
    vector<long>                        connect(4);

    // Generate list of vertices -------------------------------------------- //
    dx = (xlim[1] - xlim[0])/double(nx);
    dy = (ylim[1] - ylim[0])/double(ny);
    point[2] = 0.;
    for (j = 0; j < ny+1; ++j) {
        point[1] = ylim[0] + ((double) j) * dy;
        for (i = 0; i < nx+1; ++i) {
            point[0] = xlim[0] + ((double) i) * dx;
            mesh.addVertex(point);
        } //next i
    } //next j

    // Generate list of cells ----------------------------------------------- //
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            connect[0] = (nx+1)*j + i;
            connect[1] = (nx+1)*j + i+1;
            connect[2] = (nx+1)*(j+1) + i+1;
            connect[3] = (nx+1)*(j+1) + i;
            mesh.addCell(ElementInfo::QUAD, true, connect);
        } //next i
    } //next j
}
    
return;
}

// ========================================================================== //
// SUBTEST #001 Communications among 3 processes                              //
// ========================================================================== //
void COM_step(
    PatchKernel                 &mesh,
    short                        snd,
    short                        rcv,
    vector<long>                &id_list,
    int                          step_id
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
SurfUnstructured                     envelope(0);

// Counters
// none

// ========================================================================== //
// MOVE CELL                                                                  //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none
    log::cout() << "rank#" << snd << " send {" << id_list << "} to rank#" << rcv << endl;
    mesh.sendCells(snd, rcv, id_list);

    // Extracting mesh boundaries ------------------------------------------- //
    log::cout() << "(rank #" << mesh.getRank() << ") extracting external envelope" << endl;
    mesh.extractEnvelope(envelope);

    // Export mesh ---------------------------------------------------------- //
    stringstream                f1, f2;
    f1 << "test00001_subtest001_step" << step_id;
    mesh.write(f1.str());
    log::cout() << "(rank #" << mesh.getRank() << ", mesh exported to \"" << f1.str() << ".vtu\")" << endl;
    f2 << "P" << utils::string::zeroPadNumber(6, mesh.getRank()) << "_env_step" << step_id;
    envelope.write(f2.str());
    log::cout() << "(rank #" << mesh.getRank() << ", mesh external envelope exported to \"" << f2.str() << ".vtu\")" << endl;
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
SurfUnstructured                     mesh(0);

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
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00001 - sub-test #001 - Testing mesh communications between **" << endl;
    log::cout() << "**                               3 processes.                        **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

// ========================================================================== //
// GENERATE 2d NON-MANIFOLD TRIANGULATION                                     //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    stringstream                name;
    SurfUnstructured                envelope(0);

    // Generate dummy triangulation ----------------------------------------- //
    if (mesh.getRank() == 0) {
        log::cout() << "(rank #0) Generating 2D non-manifold triangulation" << endl;
        generateTestTriangulation(mesh);
    }
    mesh.buildAdjacencies();

    // Extracting mesh boundaries ------------------------------------------- //
    log::cout() << "(rank #" << mesh.getRank() << ") extracting external envelope" << endl;
    mesh.extractEnvelope(envelope);

    // Export mesh ---------------------------------------------------------- //
    mesh.write("test00001_subtest001_step0");
    log::cout() << "(rank #" << mesh.getRank() << ", mesh exported to \"test00001_subtest001_step0.vtu\")" << endl;
    name << "P" << utils::string::zeroPadNumber(6, mesh.getRank()) <<"_env_step0";
    envelope.write(name.str());
    log::cout() << "(rank #" << mesh.getRank() << ", mesh external envelope exported to " << name.str() << ".vtu)" << endl;
    log::cout() << endl;
}

// ========================================================================== //
// RANK#0 SENDING CELLS {30,31} TO RANK #2                                    //
// ========================================================================== //
{
    vector<long>                        cell_list{30,31};
    COM_step(mesh, 0, 2, cell_list, 1);
    log::cout() << endl;
}

// ========================================================================== //
// RANK#0 SENDING CELLS {7,8,9,23,22,24} TO RANK #1                           //
// ========================================================================== //
{
    vector<long>                        cell_list{7,8,9,23,22,24};
    COM_step(mesh, 0, 1, cell_list, 2);
    log::cout() << endl;
}

// ========================================================================== //
// RANK#2 SENDING CELLS {0,1} TO RANK #1                                      //
// ========================================================================== //
{
    vector<long>                        cell_list{0,1};
    COM_step(mesh, 2, 1, cell_list, 3);
    log::cout() << endl;
}

// ========================================================================== //
// RANK#1 SENDING CELLS {3,4,5} TO RANK #2                                    //
// ========================================================================== //
{
    vector<long>                        cell_list{3,4,5};
    COM_step(mesh, 1, 2, cell_list, 4);
    log::cout() << endl;
}

// ========================================================================== //
// RANK#1 SENDING CELLS {0,1,2} TO RANK #0                                    //
// ========================================================================== //
{
    vector<long>                        cell_list{0,1,2};
    COM_step(mesh, 1, 0, cell_list, 5);
    log::cout() << endl;
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00001 - sub-test #001 - completed!                          **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

return 0;

}

// ========================================================================== //
// SUBTEST #002 Mesh partitioning                                             //
// ========================================================================== //
int subtest_002(
    void
) {

// ========================================================================== //
// int subtest_002(                                                           //
//     void)                                                                  //
//                                                                            //
// Mesh partitioning.                                                         //
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
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
SurfUnstructured            mesh(0);

// Counters
// none

// ========================================================================== //
// SET MESH ATTRIBUTES                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set MPI communicator ------------------------------------------------- //
    mesh.setExpert(true);
    mesh.setCommunicator(MPI_COMM_WORLD);

}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //

    // Output message ------------------------------------------------------- //
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00001 - sub-test #002 - Testing mesh partitioning           **" << endl;
    log::cout() << "** ================================================================= **" << endl;
}

// ========================================================================== //
// GENERATE TEST MESH                                                         // 
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;
    stringstream                        out_msg;

    // Load stl geometry ---------------------------------------------------- //
    log::cout() << "** Rank#0, initializing mesh" << endl;
    log::cout() << "   generating simple quad mesh" << endl;
    generateTestQuadMesh(8, 8, mesh);

    // Build adjacency ------------------------------------------------------ //
    log::cout() << "   building adjacencies" << endl;
    t0 = high_resolution_clock::now();
    mesh.buildAdjacencies();
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    // Display Stats -------------------------------------------------------- //
    mesh.displayTopologyStats(out_msg, 3);
    log::cout() << out_msg.str() << endl;
}

{
    // Scope variables ------------------------------------------------------ //
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;

    // Export final mesh ---------------------------------------------------- //
    log::cout() << "   exporting initial mesh" << endl;
    t0 = high_resolution_clock::now();
    mesh.write("test00001_subtest002_init");
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;
}


// ========================================================================== //
// INITIAL PARTIONING                                                         //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    vector<long>                        cell_list1{32,33,34,35,
                                                   40,41,42,43,
                                                   48,49,50,51,
                                                   56,57,58,59};
    vector<long>                        cell_list2{36,37,38,39,
                                                   44,45,46,47,
                                                   52,53,54,55,
                                                   60,61,62,63};
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;

    // Send cells to neighboring processors --------------------------------- //
    log::cout() << "** Rank#" << mesh.getRank() << ", partitioning mesh" << endl;

    log::cout() << "   sending cell: " << cell_list1 << " from 0 to 1" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 1, cell_list1);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    log::cout() << "   sending cell: " << cell_list2 << " from 0 to 2" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 2, cell_list2);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    // Export mesh ---------------------------------------------------------- //
    log::cout() << "   exporting mesh to \"test00001_subtest002_step0.vtu\"" << endl;
    t0 = high_resolution_clock::now();
    mesh.write("test00001_subtest002_step0");
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;
}

// ========================================================================== //
// STEP 2                                                                     //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    vector<long>                        cell_list{18,19,20,21,22,23,
                                                  26,27,28,29,30,31};
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;

    // Send cells to neighboring processors --------------------------------- //
    log::cout() << "** Rank#" << mesh.getRank() << ", partitioning mesh" << endl;

    log::cout() << "   sending cell: " << cell_list << " from 0 to 2" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 2, cell_list);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    // Export mesh ---------------------------------------------------------- //
    log::cout() << "   exporting mesh to \"test00001_subtest002_step1.vtu\"" << endl;
    t0 = high_resolution_clock::now();
    mesh.write("test00001_subtest002_step1");
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;
}


// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //

    // Output message ------------------------------------------------------- //
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00001 - sub-test #002 - completed!                          **" << endl;
    log::cout() << "** ================================================================= **" << endl;
}
return 0;

}

// ========================================================================== //
// SUBTEST #003 Mesh partitioning                                             //
// ========================================================================== //
int subtest_003(
    void
) {

// ========================================================================== //
// int subtest_003(                                                           //
//     void)                                                                  //
//                                                                            //
// Mesh partitioning.                                                         //
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
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
SurfUnstructured            mesh(0);

// Counters
// none

// ========================================================================== //
// SET MESH ATTRIBUTES                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set MPI communicator ------------------------------------------------- //
    mesh.setExpert(true);
    mesh.setCommunicator(MPI_COMM_WORLD);

}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //

    // Output message ------------------------------------------------------- //
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00001 - sub-test #003 - Testing mesh partitioning           **" << endl;
    log::cout() << "** ================================================================= **" << endl;
}

// ========================================================================== //
// GENERATE TEST MESH                                                         // 
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;
    stringstream                        out_msg;

    // Load stl geometry ---------------------------------------------------- //
    log::cout() << "** Rank#0, initializing mesh" << endl;
    log::cout() << "   generating simple quad mesh" << endl;
    generateTestQuadMesh(9, 9, mesh);

    // Build adjacency ------------------------------------------------------ //
    log::cout() << "   building adjacencies" << endl;
    t0 = high_resolution_clock::now();
    mesh.buildAdjacencies();
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    // Display Stats -------------------------------------------------------- //
    mesh.displayTopologyStats(out_msg, 3);
    log::cout() << out_msg.str() << endl;
}

{
    // Scope variables ------------------------------------------------------ //
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;

    // Export final mesh ---------------------------------------------------- //
    log::cout() << "   exporting initial mesh" << endl;
    t0 = high_resolution_clock::now();
    mesh.write("test00001_subtest003_init");
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

}

// ========================================================================== //
// INITIAL PARTIONING                                                         //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    vector<long>                        cell_list1{54,55,56,
                                                   63,64,65,
                                                   72,73,74};
    vector<long>                        cell_list2{57,58,59,
                                                   66,67,68,
                                                   75,76,77};
    vector<long>                        cell_list3{60,61,62,
                                                   69,70,71,
                                                   78,79,80};
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;

    // Send cells to neighboring processors --------------------------------- //
    log::cout() << "** Rank#" << mesh.getRank() << ", partitioning mesh" << endl;

    log::cout() << "   sending cell: " << cell_list1 << " from 0 to 1" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 1, cell_list1);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    log::cout() << "   sending cell: " << cell_list2 << " from 0 to 2" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 2, cell_list2);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    log::cout() << "   sending cell: " << cell_list3 << " from 0 to 3" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 3, cell_list3);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    // Export mesh ---------------------------------------------------------- //
    log::cout() << "   exporting mesh to \"test00001_subtest002_step0.vtu\"" << endl;
    t0 = high_resolution_clock::now();
    mesh.write("test00001_subtest003_step0");
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;
}

// ========================================================================== //
// STEP 2                                                                     //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    vector<long>                        cell_list{28,29,30,31,32,33,34,35,
                                                  37,38,39,40,41,42,43,44,
                                                  46,47,48,49,50,51,52,53};
    high_resolution_clock::time_point   t0, t1;
    duration<double>                    time_span;

    // Send cells to neighboring processors --------------------------------- //
    log::cout() << "** Rank#" << mesh.getRank() << ", partitioning mesh" << endl;

    log::cout() << "   sending cell: " << cell_list << " from 0 to 3" << endl;
    t0 = high_resolution_clock::now();
    mesh.sendCells(0, 3, cell_list);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;

    // Export mesh ---------------------------------------------------------- //
    log::cout() << "   exporting mesh to \"test00001_subtest002_step1.vtu\"" << endl;
    t0 = high_resolution_clock::now();
    mesh.write("test00001_subtest003_step1");
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    log::cout() << "     (" << time_span.count() << " sec.)" << endl;
}


// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
if (mesh.getRank() == 0) {

    // Scope variables ------------------------------------------------------ //

    // Output message ------------------------------------------------------- //
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00001 - sub-test #003 - completed!                          **" << endl;
    log::cout() << "** ================================================================= **" << endl;
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
// INITIALIZE LOGGER                                                          //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int nProcs;
    int rank;

    // Initialize logger
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);
}

// ========================================================================== //
// RUN SUB-TEST #001                                                          //
// ========================================================================== //
err = subtest_001();
if (err > 0) return(10 + err);

// ========================================================================== //
// RUN SUB-TEST #002                                                          //
// ========================================================================== //
err = subtest_002();
if (err > 0) return(20 + err);

// ========================================================================== //
// RUN SUB-TEST #003                                                          //
// ========================================================================== //
err = subtest_003();
if (err > 0) return(30 + err);

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
