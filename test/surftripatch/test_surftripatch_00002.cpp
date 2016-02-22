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
// SUBTEST #001 Test geometrical checks                                       //
// ========================================================================== //
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// Test computation of min edge/max edge and edge length                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
//               err = 1  --> error at step #1                                //
//               err = 2  --> error at step #2                                //
// ========================================================================== //
    
// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long                            id;
SurfTriPatch                    mesh(0);
vector<long>                    c_connect(3, Element::NULL_ID);

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
    cout << "** ================================================================= **" << endl;
    cout << "** Test #00001 - sub-test #001 - Testing geometrical check           **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

// ========================================================================== //
// INITIALIZE MESH                                                            //
// ========================================================================== //
{
    // Scope variables
    SurfTriPatch::VertexIterator                        vit;
    SurfTriPatch::CellIterator                          cit;
    
    // Initialize internal cell
    cout << "** Initializing mesh" << endl;
    vit = mesh.addVertex(array<double, 3>{0.0, 0.0, 0.0});
    c_connect[0] = vit->get_id();
    vit = mesh.addVertex(array<double, 3>{1.0, 0.0, 0.0});
    c_connect[1] = vit->get_id();
    vit = mesh.addVertex(array<double, 3>{0.0, 1.0, 0.0});
    c_connect[2] = vit->get_id();
    cit = mesh.addCell(ElementInfo::TRIANGLE, true, c_connect);
    id = cit->get_id();

}

// ========================================================================== //
// TEST (MIN/MAX) EDGE LENGTH                                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    Cell                        *cell_ = &mesh.getCell(id);
    int                         nedges = cell_->getEdgeCount();

    // Display cell --------------------------------------------------------- //
    cout << "** Mesh data" << endl;
    cout << "   Topology:" << endl;
    mesh.displayTopologyStats(cout, 5);
    cout << "   Vertex list:" << endl;
//TODO:     mesh.displayVertex(cout, 5);
    cout << "   Cell list:" << endl;
    mesh.displayCells(cout, 5);
    cout << endl;

    // Check edge length ---------------------------------------------------- //
    cout << "   Edge length for cell " << id << ": " << endl;
    for (int i = 0; i < nedges; ++i) {
        cout << "     edge loc. id = " << i << ", edge length = " << mesh.evalEdgeLength(id, i) << endl;
    } //next i
    cout << "     min. edge = " << mesh.evalMinEdgeLength(id) << endl;
    cout << "     max. edge = " << mesh.evalMaxEdgeLength(id) << endl;
    cout << endl;

}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    cout << "** ================================================================= **" << endl;
    cout << "** Test #00001 - sub-test #001 - completed!                          **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

return 0; }

// ========================================================================== //
// MAIN FOR TEST #00002                                                       //
// ========================================================================== //
int main(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variabels
int                             err = 0;

// ========================================================================== //
// RUN SUB-TEST #001                                                          //
// ========================================================================== //
err = subtest_001();
if (err > 0) return(10 + err);

return err;

}
