// ========================================================================== //
//           ** BitPit mesh ** Test 004 for class surftri_patch **            //
//                                                                            //
// Test topological queries on SurfTriPatch.                                  //
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
// SUBTEST #001 Test edge network extraction                                  //
// ========================================================================== //
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// Test extraction of edge network from surface mesh.                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err     : int, error code:                                               //
//             err = 0 --> no error(s) encountered                            //
//             err = 1 --> error at step#1                                    //
//             err = 2 --> error at step#2                                    //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
SurfTriPatch                    mesh(0), edges(1);

// Counters
int                             nV, nS, nE;

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    cout << "** ================================================================= **" << endl;
    cout << "** Test #00004 - sub-test #001 - edge extraction                     **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

// ========================================================================== //
// SET MESH ATTRIBUTES                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set mesh attributes -------------------------------------------------- //
    mesh.setExpert(true);
    edges.setExpert(true);
}

// ========================================================================== //
// STEP#1 IMPORT MESH FROM STL FILE                                           //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    string                      stl_name = "./data/cube.stl";
    bool                        stl_type = false;

    // Load stl file -------------------------------------------------------- //
    cout << "** Importing surface mesh from \"" << stl_name << "\"" << endl;
    mesh.importSTL(stl_name, stl_type);
    cout << "** Removing coincident vertices" << endl;
    mesh.deleteCoincidentVertex();
    cout << "** Building adjacencies" << endl;
    mesh.buildAdjacencies();

    // Check mesh status ---------------------------------------------------- //
    nV = mesh.getVertexCount();
    nS = mesh.getInternalCount();
    nE = mesh.countFaces(); 
//     if (nS != 283274)   return 1;
//     if (nV != 141639)   return 1;
//     if (nE != 424911)   return 1;

    // Export surface to vtk file ------------------------------------------- //
    cout << "** Exporting mesh to \"test_00004_subtest_001_surf.vtu\"" << endl;
    mesh.write("test_00004_subtest_001_surf");
    
}

// ========================================================================== //
// STEP#2 IMPORT MESH FROM STL FILE                                           //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    bool                                        check;
    int                                         n_faces, n_adj, n_vert;
    int                                         i, j;
    long                                        id;
    SurfTriPatch::CellIterator                  c_, ce_;
    SurfTriPatch::VertexIterator                v_, ve_;
    long                                        edge_counter = 0, vertex_counter = 0;
    unordered_map<long, long>                   vertex_mapper;
    vector<array<double, 3>>                    normals, vnormals;

    // Extract edges -------------------------------------------------------- //
    cout << "** Extracting edges" << endl;
    mesh.extractEdgeNetwork(edges);
    
    // Check mesh status ---------------------------------------------------- //
    if (edges.getInternalCount() != nE)         return 2;
    if (edges.getVertexCount()   != nV)         return 2;

    // Compute edge normals ------------------------------------------------- //
    cout << "** Computing edge normals" << endl;
    normals.resize(edges.getCellCount());
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        n_faces = c_->getFaceCount();
        id = c_->get_id();
        for (i = 0; i < n_faces; ++i) {
            check = true;
            n_adj = c_->getAdjacencyCount(i);
            for (j = 0; j < n_adj; ++j) {
                check = check && (id > c_->getAdjacency(i, j));
            } //next j
            if (check) {
                normals[edge_counter] = mesh.evalEdgeNormal(id, i);
                ++edge_counter;
            }
        } //next i
    } //next c_

    // Compute vertex normals ----------------------------------------------- //
    cout << "** Computing vertex normals" << endl;
    vnormals.resize(edges.getVertexCount());
    ve_ = mesh.vertexEnd();
    for (v_ = mesh.vertexBegin(); v_ != ve_; ++v_) {
        vertex_mapper[v_->get_id()] = vertex_counter;
        ++vertex_counter;
    } //next v_
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        n_vert = c_->getVertexCount();
        for (i = 0; i < n_vert; ++i) {
            id = c_->getVertex(i);
            vnormals[vertex_mapper[id]] = mesh.evalVertexNormal(c_->get_id(), i);
        } //next i
    } //next c_

    // Export edges to vtu file --------------------------------------------- //
    cout << "** Exporting edges to \"test_00004_subtest_001_edges.vtu\"" << endl;
    edges.addData("normals", VTKFieldType::VECTOR, VTKLocation::CELL, normals);
    edges.addData("vnormals", VTKFieldType::VECTOR, VTKLocation::POINT, vnormals);
    edges.write("test_00004_subtest_001_edges");
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
    cout << "** Test #00004 - sub-test #001 - completed!                          **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

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

return err;

}
