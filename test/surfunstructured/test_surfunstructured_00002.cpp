// ========================================================================== //
//           ** BitPit mesh ** Test 002 for class SurfUnstructured **         //
//                                                                            //
// Test routines for geometrical queries for class SurfUnstructured.          //
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
# define _USE_MATH_DEFINES

# include <cmath>
# include <array>
# include <vector>
# include <iostream>
#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

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

// ========================================================================== //
// SUBTEST #001 Test geometrical checks                                       //
// ========================================================================== //
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// Test routines for geometrical computations                                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
//               err = 1  --> error at step #1 (edge length calculations)     //
//               err = 2  --> error at step #2 (angle calculations)           //
//               err = 3  --> error at step #3 (facet normal calculations)    //
//               err = 4  --> error at step #4 (aspect ratio calculations)    //
// ========================================================================== //
    
// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long                            id;
SurfUnstructured                mesh(2, 3);

// Counters
// none

// ========================================================================== //
// INITIALIZE MESH PARAMETERS                                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Enable changes ------------------------------------------------------- //
    mesh.setExpert(true);
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00002 - sub-test #001 - Testing geometrical check           **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

// ========================================================================== //
// INITIALIZE MESH                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    SurfUnstructured::VertexIterator                    vit;
    SurfUnstructured::CellIterator                      cit;
    vector<long>                                        t_connect(3, Vertex::NULL_ID);
    vector<long>                                        q_connect(4, Vertex::NULL_ID);
    
    // Initialize internal cell --------------------------------------------- //
    log::cout() << "** Initializing mesh" << endl;

    // Place a triangle
    vit = mesh.addVertex(array<double, 3>{{0.0, 0.0, 0.0}});
    t_connect[0] = vit->getId();
    q_connect[0] = vit->getId();
    vit = mesh.addVertex(array<double, 3>{{1.0, 0.0, 0.0}});
    t_connect[1] = vit->getId();
    q_connect[3] = vit->getId();
    vit = mesh.addVertex(array<double, 3>{{0.0, 1.0, 0.0}});
    t_connect[2] = vit->getId();
    cit = mesh.addCell(ElementType::TRIANGLE, t_connect);

    // Place a quad
    vit = mesh.addVertex(array<double, 3>{{0.0, -1.0, 0.0}});
    q_connect[1] = vit->getId();
    vit = mesh.addVertex(array<double, 3>{{1.0, -1.0, 0.0}});
    q_connect[2] = vit->getId();
    cit = mesh.addCell(ElementType::QUAD, q_connect);
    

    // Display mesh data ---------------------------------------------------- //
    log::cout() << "   Topology:" << endl;
    mesh.displayTopologyStats(log::cout(), 5);
    log::cout() << "   Vertex list:" << endl;
    mesh.displayVertices(log::cout(), 5);
    log::cout() << "   Cell list:" << endl;
    mesh.displayCells(log::cout(), 5);
    log::cout() << endl;

}

// ========================================================================== //
// STEP #1 TEST (MIN/MAX) EDGE LENGTH                                         //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         nedges, edge_id;
    double                      length;
    vector<vector<double>>      expected(2);
    double                      m_length, M_length;
    int                         i;

    // Output message ------------------------------------------------------- //
    log::cout() << "** Testing routines for edge length calculations" << endl;

    // Ref. value to check against ------------------------------------------ //
    expected[0] = vector<double>{1.0, sqrt(2.0), 1.0};
    expected[1].resize(4, 1.0);

    // Check edge length for each cell witin the mesh ----------------------- //
    SurfUnstructured::CellIterator  cell_, end_ = mesh.cellEnd();
    i = 0;
    for (cell_ = mesh.cellBegin(); cell_ != end_; ++cell_) {

        // Cell data
        nedges = cell_->getFaceCount();
        id = cell_->getId();

        // Initialize ref. value for checks
        minval(expected[i], m_length);
        maxval(expected[i], M_length);

        // Compute cell edge
        log::cout() << "   Edge(s) length for cell " << id << ": " << endl;
        for (int j = 0; j < nedges; ++j) {
            length =  mesh.evalEdgeLength(id, j);
            log::cout() << "     edge loc. id = " << j << ", edge length = " << length << endl;
            if (abs(length - expected[i][j]) > 1.0e-12) return 1;
        } //next i
        length = mesh.evalMinEdgeLength(id, edge_id);
        log::cout() << "     min. edge = " << length << " on edge: " << edge_id << endl;
        if (abs(length - m_length) > 1.0e-12) return 1;
        length = mesh.evalMaxEdgeLength(id, edge_id);
        log::cout() << "     max. edge = " << length << " on edge: " << edge_id << endl;
        if (abs(length - M_length) > 1.0e-12) return 1;
        
        // Update counters
        ++i;

    } //next cell_
    log::cout() << endl;
}

// ========================================================================== //
// STEO #2 TEST (MIN/MAX) ANGLE AT VERTEX                                     //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         nedges, vertex_id;
    double                      angle;
    vector<vector<double>>      expected(2);
    double                      m_angle, M_angle;
    int                         i;

    // Output message ------------------------------------------------------- //
    log::cout() << "** Testing routines for angle calculations" << endl;

    // Ref. value to check against ------------------------------------------ //
    expected[0] = vector<double>{BITPIT_PI_2, BITPIT_PI_4, BITPIT_PI_4};
    expected[1].resize(4, BITPIT_PI_2);

    // Check edge length ---------------------------------------------------- //
    SurfUnstructured::CellIterator  cell_, end_ = mesh.cellEnd();
    i = 0;
    for (cell_ = mesh.cellBegin(); cell_ != end_; ++cell_) {

        // Cell data
        id = cell_->getId();
        nedges = cell_->getVertexCount();

        // Initialize ref. value for check
        minval(expected[i], m_angle);
        maxval(expected[i], M_angle);

        // Compute cell angles
        log::cout() << "   Angle(s) for cell " << id << ": " << endl;
        for (int j = 0; j < nedges; ++j) {
            angle = mesh.evalAngleAtVertex(id, j);
            log::cout() << "     vertex loc. id = " << j << ", angle = " << angle << " [rad]" << endl;
            if (abs(angle - expected[i][j]) > 1.0e-12) return 2;
        } //next i
        angle = mesh.evalMinAngleAtVertex(id, vertex_id);
        log::cout() << "     min. angle = " << angle << " on vertex: " << vertex_id << endl;
        if (abs(angle - m_angle) > 1.0e-12) return 2;
        angle = mesh.evalMaxAngleAtVertex(id, vertex_id);
        log::cout() << "     max. angle = " << angle << " on vertex: " << vertex_id << endl;
        if (abs(angle - M_angle) > 1.0e-12) return 2;

        // Update counters
        ++i;

    } //next cell_
    log::cout() << endl;

}

// ========================================================================== //
// STEP #3 TEST FACET NORMAL COMPUTATION                                      //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    array<double, 3>             normal;
    vector<array<double, 3>>     expected(2, array<double, 3>{{0.0, 0.0, 1.0}});
    int                          i;

    // Output message ------------------------------------------------------- //
    log::cout() << "** Testing routines for facet normal calculations" << endl;

    // Compute face normal -------------------------------------------------- //
    SurfUnstructured::CellIterator  cell_, end_ = mesh.cellEnd();
    i = 0;
    for (cell_ = mesh.cellBegin(); cell_ != end_; ++cell_) {

        // Cell data
        id = cell_->getId();

        // Compute face normal
        log::cout() << "   Facet normal for cell " << id << ": " << endl;
        normal = mesh.evalFacetNormal(id);
        log::cout() << "     normal: " << normal << endl;
        if (norm2(normal - expected[i]) > 1.0e-12) return 3;

    } //next cell_
    log::cout() << endl;

}

// ========================================================================== //
// STEP #4 TEST ASPECT RATIO COMPUTATION                                      //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    double                      ar;
    vector<double>              expected{sqrt(2), 1.0};
    int                         i;
    int                         edge_id;

    // Output message ------------------------------------------------------- //
    log::cout() << "** Testing routines for aspect ratio computation" << endl;

    // Compute face normal -------------------------------------------------- //
    SurfUnstructured::CellIterator  cell_, end_ = mesh.cellEnd();
    i = 0;
    for (cell_ = mesh.cellBegin(); cell_ != end_; ++cell_) {

        // Cell data
        id = cell_->getId();

        // Compute face normal
        log::cout() << "   Aspect ratio for cell " << id << ": " << endl;
        ar = mesh.evalAspectRatio(id, edge_id);
        log::cout() << "     a.r.: " << ar << " (shortest edge has loc. id: " << edge_id << ")" << endl;
        if (abs(ar - expected[i]) > 1.0e-12) return 4;

        // Update counters
        ++i;

    } //next cell_
    log::cout() << endl;

}

// ========================================================================== //
// STEP #5 TEST FACET AREA COMPUTATION                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    vector<double>              expected_area{0.5, 1.0};
    vector<double>              expected_size{sqrt(0.5), 1.0};
    int                         i;
    double                      area, c_size;

    // Output message ------------------------------------------------------- //
    log::cout() << "** Testing routines for facet area computation" << endl;

    // Compute face normal -------------------------------------------------- //
    SurfUnstructured::CellIterator  cell_, end_ = mesh.cellEnd();
    i = 0;
    for (cell_ = mesh.cellBegin(); cell_ != end_; ++cell_) {

        // Cell data
        id = cell_->getId();

        // Compute face normal
        log::cout() << "   Facet area for cell " << id << ": " << endl;
        area = mesh.evalCellArea(id);
        c_size = mesh.evalCellSize(id);
        log::cout() << "     area: " << area << endl;
        log::cout() << "     cell size: " << c_size << endl;
        if (abs(area - expected_area[i]) > 1.0e-12) return 5;
        if (abs(c_size - expected_size[i]) > 1.0e-12) return 5;

        // Update counters
        ++i;

    } //next cell_
    log::cout() << endl;

}

// ========================================================================== //
// STEP #6 CELL CENTER COMPUTATION                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    vector<array<double, 3>>    expected_center(2);
    int                         i;
    array<double, 3>            C;

    // Output message ------------------------------------------------------- //
    log::cout() << "** Testing routines for cell's center eval" << endl;

    // Initialize ref. values ----------------------------------------------- //
    expected_center[0] = array<double, 3>{{1./3., 1./3., 0.}};
    expected_center[1] = array<double, 3>{{0.5, -0.5, 0.}};

    // Compute face normal -------------------------------------------------- //
    SurfUnstructured::CellIterator  cell_, end_ = mesh.cellEnd();
    i = 0;
    for (cell_ = mesh.cellBegin(); cell_ != end_; ++cell_) {

        // Cell data
        id = cell_->getId();

        // Compute cell center
        log::cout() << "   Cell center for cell " << id << ": " << endl;
        C = mesh.evalCellCentroid(id);
        log::cout() << "     C: (" << C << ")" << endl;
        if (norm2(C - expected_center[i]) > 1.0e-12) return 6;

        // Update counters
        ++i;

    } //next cell_
    log::cout() << endl;

}

// ========================================================================== //
// STEP #7 TEST HISTOGRAM CONSTRUCTION                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    vector<double>              hist, bins;

    // Test histogram construction ------------------------------------------ //
    log::cout() << "** Testing histogram construction" << endl;
    mesh.displayQualityStats(log::cout(), 3);
    log::cout() << endl;

}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00002 - sub-test #001 - completed!                          **" << endl;
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
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

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
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
