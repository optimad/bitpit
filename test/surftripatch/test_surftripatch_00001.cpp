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
// SUBTEST #001 Test for cell insertion and deletion                          //
// ========================================================================== //
int subtest_001(
    void
) {

// ========================================================================== //
// int test_001(                                                              //
//     void)                                                                  //
//                                                                            //
// Test insertion order.                                                      //
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
SurfTriPatch                    mesh(0);
vector<long>                    c_connect{0, 1, 2};
vector<long>                    g_connect{3, 4, 5};
Cell                            cell(0, ElementInfo::TRIANGLE);
Cell                            ghost(0, ElementInfo::TRIANGLE);
vector<long>                    expected;
vector<bool>                    internal;

// Counters
int                             i;

// ========================================================================== //
// INITIALIZE CELL                                                            //
// ========================================================================== //
{
    // Scope variables
    int                         j;
    int                         n;

    // Initialize internal cell
    cout << "** Initializing cell" << endl;
    n = cell.getVertexCount();
    for (j = 0; j < n; ++j) {
        cell.setVertex(j, c_connect[j]);
    } //next j
    n = ghost.getVertexCount();
    for (j = 0; j < n; ++j) {
        ghost.setVertex(j, g_connect[j]);
    } //next j
}

// ========================================================================== //
// INSERT CELLS (STEP #1)                                                     //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    const int                                   N = 5;
    SurfTriPatch::CellIterator                  it, et;

    // Insert internal cells (IDX 0-4) -------------------------------------- //
    // cells:  {0,1,2,3,4}
    // ghosts: {}
    cout << "** Inserting interal cells" << endl;
    for (i = 0; i < N/2; ++i) {
        mesh.addCell(cell);
        expected.push_back(long(i));
        internal.push_back(true);
    } //next i
    for (i = N/2; i < N; ++i) {
        mesh.addCell(ElementInfo::TRIANGLE, true, c_connect);
        expected.push_back(long(i));
        internal.push_back(true);
    } //next i

    // Check cell ordering -------------------------------------------------- //
    i = 0;
    et = mesh.end();
    for (it = mesh.begin(); it != et; ++it) {
        if (it->get_id() != expected[i]) return 1;
        if (it->isInterior() != internal[i]) return 1;
        ++i;
    } //next it

    // Display mesh content ------------------------------------------------- //
    cout << "** After inserting internal cells" << endl;
    et = mesh.end();
    for (it = mesh.begin(); it != et; ++it) {
        cout << "  cell: " << endl;
        it->display(cout, 4);
    } //next it

    // Insert ghost cells (IDX 5-9) ----------------------------------------- //
    // cells:  {0,1,2,3,4}
    // ghosts: {5,6,7,8,9}
    cout << "** Inserting ghost cells" << endl;
    for (i = 0; i < N/2; ++i) {
        mesh.AddGhost(ghost);
        expected.push_back(long(N + i));
        internal.push_back(false);
    } //next i
    for (i = N/2; i < N; ++i) {
        mesh.AddGhost(ElementInfo::TRIANGLE, false, g_connect);
        expected.push_back(long(N + i));
        internal.push_back(false);
    }

    // Check cells ordering ------------------------------------------------- //
    i = 0;
    et = mesh.end();
    for (it = mesh.begin(); it != et; ++it) {
        if (it->get_id() != expected[i]) return 1;
        if (it->isInterior() != internal[i]) return 1;
        ++i;
    } //next it

    // Display mesh content ------------------------------------------------- //
    cout << "** After inserting ghost cells" << endl;
    et = mesh.end();
    for (it = mesh.begin(); it != et; ++it) {
        cout << "  cell: " << endl;
        it->display(cout, 4);
    } //next it

}

// // ========================================================================== //
// // REMOVE INSERT/CELLS (STEP #2)                                              //
// // ========================================================================== //
// {
//     // Scope variables
//     Class_PMesh::cell_iterator                  it, et;
// 
//     // Remove internal cells
//     //cells:  {0,1,-1,3,-1}
//     //ghosts: {-1,-1,7,8,9}
//     mesh.DeleteCell(4);
//     mesh.DeleteCell(2);
//     mesh.DeleteGhost(5);
//     mesh.DeleteGhost(6);
//     expected.erase(expected.begin() + 6);
//     expected.erase(expected.begin() + 5);
//     expected.erase(expected.begin() + 4);
//     expected.erase(expected.begin() + 2);
//     internal.erase(internal.begin() + 6);
//     internal.erase(internal.begin() + 5);
//     internal.erase(internal.begin() + 4);
//     internal.erase(internal.begin() + 2);
// 
//     // Check element order
//     i = 0;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         if (it->get_id() != expected[i]) return 2;
//         if (it->isInterior() != internal[i]) return 2;
//         ++i;
//     } //next it
// 
//     // Display mesh
//     cout << "** After removing internal/ghost cells" << endl;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         cout << "  cell: " << endl;
//         it->display(cout, 4);
//     } //next it
// 
//     // Remove ghost cells
//     //cells:  {0,1,12,3,13}
//     //ghosts: {11,10,7,8,9}
//     mesh.AddGhost(ElementInfo::TRIANGLE, g_connect);
//     mesh.AddGhost(ghost);
//     mesh.AddCell(cell);
//     mesh.AddCell(ElementInfo::TRIANGLE, c_connect);
//     expected.insert(expected.begin() + 3, 10);
//     expected.insert(expected.begin() + 3, 11);
//     expected.insert(expected.begin() + 2, 12);
//     expected.insert(expected.begin() + 4, 13);
//     internal.insert(internal.begin() + 3, false);
//     internal.insert(internal.begin() + 3, false);
//     internal.insert(internal.begin() + 2, true);
//     internal.insert(internal.begin() + 4, true);
// 
//     // Check element order
//     i = 0;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         if (it->get_id() != expected[i]) return 2;
//         ++i;
//     } //next it
// 
//     cout << "** After inserting internal/ghost cells" << endl;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         cout << "  cell: " << endl;
//         it->display(cout, 4);
//     } //next it
// 
//     // Remove all internal cells and add 2 ghost cells
//     //cells:  {}
//     //ghosts: {14,15,11,10,7,8,9}
//     mesh.DeleteCell(13);
//     mesh.DeleteCell(1);
//     mesh.DeleteCell(12);
//     mesh.DeleteCell(0);
//     mesh.DeleteCell(3);
//     mesh.AddGhost(ghost);
//     mesh.AddGhost(ElementInfo::TRIANGLE, g_connect);
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.insert(expected.begin(), 15);
//     expected.insert(expected.begin(), 14);
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.insert(internal.begin(), false);
//     internal.insert(internal.begin(), false);
// 
//     // Check element order
//     i = 0;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         if (it->get_id() != expected[i]) return 2;
//         if (it->isInterior() != internal[i]) return 2;
//         ++i;
//     } //next it
// 
//     cout << "** After erasing all internal cells and inserting 2 new ghosts" << endl;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         cout << "  cell: " << endl;
//         it->display(cout, 4);
//     } //next it
// 
//     // Remove all ghosts add 2 internal cells
//     //cells:  {16,17}
//     //ghosts: {}
//     mesh.DeleteGhost(14);
//     mesh.DeleteGhost(10);
//     mesh.DeleteGhost(11);
//     mesh.DeleteGhost(15);
//     mesh.DeleteGhost(9);
//     mesh.DeleteGhost(7);
//     mesh.DeleteGhost(8);
//     mesh.AddCell(cell);
//     mesh.AddCell(ElementInfo::TRIANGLE, c_connect);
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.erase(expected.begin());
//     expected.insert(expected.begin(),17);
//     expected.insert(expected.begin(),16);
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.erase(internal.begin());
//     internal.insert(internal.begin(),true);
//     internal.insert(internal.begin(),true);
// 
//     // Check element order
//     i = 0;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         if (it->get_id() != expected[i]) return 2;
//         if (it->isInterior() != internal[i]) return 2;
//         ++i;
//     } //next it
// 
//     cout << "** After erasing all ghost cells and inserting 2 new internal cells" << endl;
//     et = mesh.end();
//     for (it = mesh.begin(); it != et; ++it) {
//         cout << "  cell: " << endl;
//         it->display(cout, 4);
//     } //next it
// 
// }

return 0; }

// ========================================================================== //
// MAIN FOR TEST #00001                                                       //
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

// ========================================================================== //
// RUN SUB-TEST #002                                                          //
// ========================================================================== //

// ========================================================================== //
// RUN SUB-TEST #003                                                          //
// ========================================================================== //

return err;

}
