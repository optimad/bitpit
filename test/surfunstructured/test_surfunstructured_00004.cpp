// ========================================================================== //
//           ** BitPit mesh ** Test 004 for class SurfUnstructured **         //
//                                                                            //
// Test topological queries on SurfUnstructured.                              //
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
#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// BitPit
# include "bitpit_common.hpp"                                                 // Utilities and common definitions
# include "bitpit_IO.hpp"                                                     // Input/output
# include "bitpit_operators.hpp"                                              // STL containers operators
# include "bitpit_patchkernel.hpp"                                            // BitPit base patch
# include "bitpit_surfunstructured.hpp"                                           // BitPit surftri patch

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// SUBTEST #001 Test edge network extraction for 3D surface                   //
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
SurfUnstructured                mesh(0), edges(1);

// Counters
int                             nV, nE;

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00004 - sub-test #001 - edge extraction (surf in 3D)        **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
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
    log::cout() << "** Importing surface mesh from \"" << stl_name << "\"" << endl;
    mesh.importSTL(stl_name, stl_type);
    log::cout() << "** Removing coincident vertices" << endl;
    mesh.deleteCoincidentVertices();
    log::cout() << "** Building adjacencies" << endl;
    mesh.buildAdjacencies();

    // Check mesh status ---------------------------------------------------- //
    nV = mesh.getVertexCount();
    nE = mesh.countFaces(); 
    
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
    SurfUnstructured::CellIterator              c_, ce_;
    SurfUnstructured::VertexIterator            v_, ve_;
    long                                        edge_counter = 0, vertex_counter = 0, cell_counter = 0;
    unordered_map<long, long>                   vertex_mapper;
    vector<array<double, 3>>                    normals, enormals, vnormals;

    // Extract edges -------------------------------------------------------- //
    log::cout() << "** Extracting edges" << endl;
    mesh.extractEdgeNetwork(edges);
    
    // Check mesh status ---------------------------------------------------- //
    if (edges.getInternalCount() != nE)         return 2;
    if (edges.getVertexCount()   != nV)         return 2;

    // Compute cell normals ------------------------------------------------- //
    log::cout() << "** Computing cell normals" << endl;
    normals.resize(mesh.getCellCount());
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        normals[cell_counter] = mesh.evalFacetNormal(c_->getId());
        ++cell_counter;
    } //next c_
    
    // Compute edge normals ------------------------------------------------- //
    log::cout() << "** Computing edge normals" << endl;
    enormals.resize(edges.getCellCount());
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        n_faces = c_->getFaceCount();
        id = c_->getId();
        for (i = 0; i < n_faces; ++i) {
            check = true;
            n_adj = c_->getAdjacencyCount(i);
            for (j = 0; j < n_adj; ++j) {
                check = check && (id > c_->getAdjacency(i, j));
            } //next j
            if (check) {
                enormals[edge_counter] = mesh.evalEdgeNormal(id, i);
                ++edge_counter;
            }
        } //next i
    } //next c_

    // Compute vertex normals ----------------------------------------------- //
    log::cout() << "** Computing vertex normals" << endl;
    vnormals.resize(edges.getVertexCount());
    ve_ = mesh.vertexEnd();
    for (v_ = mesh.vertexBegin(); v_ != ve_; ++v_) {
        vertex_mapper[v_->getId()] = vertex_counter;
        ++vertex_counter;
    } //next v_
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        const long *cellConnect = c_->getConnect();
        n_vert = c_->getVertexCount();
        for (i = 0; i < n_vert; ++i) {
            id = cellConnect[i];
            vnormals[vertex_mapper[id]] = mesh.evalVertexNormal(c_->getId(), i);
        } //next i
    } //next c_

    // Export edges to vtu file --------------------------------------------- //
    {
        log::cout() << "** Exporting edges to \"test_00004_subtest_001_edges.vtu\"" << endl;
        VTKUnstructuredGrid &vtk = edges.getVTK() ;
        vtk.addData("normals", VTKFieldType::VECTOR, VTKLocation::CELL, enormals);
        vtk.addData("vnormals", VTKFieldType::VECTOR, VTKLocation::POINT, vnormals);

        edges.write("test_00004_subtest_001_edges");
    }

    {
        log::cout() << "** Exporting mesh to \"test_00004_subtest_001_surf.vtu\"" << endl;
        VTKUnstructuredGrid &vtk = mesh.getVTK() ;
        vtk.addData("normals", VTKFieldType::VECTOR, VTKLocation::CELL, normals);

        mesh.write("test_00004_subtest_001_surf");
        log::cout() << endl;
    }
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00004 - sub-test #001 - completed!                          **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

return 0; } 

// ========================================================================== //
// SUBTEST #002 Test edge network extraction for 2D curve                     //
// ========================================================================== //
void Generate2DSurfMesh(
    SurfUnstructured                    &mesh
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
const double                            _PI_ = 3.14159265358979;
const double                            R = 1.0;
const long                              N = 32;
double                                  theta;
double                                  dtheta = 2. * _PI_/((double) N);

// Counters
long                                    i;

// ========================================================================== //
// GENERATE MESH                                                              //
// ========================================================================== //
{
    // Local variables ------------------------------------------------------ //
    array<double, 3>                    point;
    vector<long>                        connect(2, Element::NULL_ID);

    // Create vertex list --------------------------------------------------- //
    point[2] = 0.0;
    for (i = 0; i < N; ++i) {
        theta = ((double) i) * dtheta;
        point[0] = R * cos( theta );
        point[1] = R * sin( theta );
        mesh.addVertex(point);
    } //next i

    // Create simplex list -------------------------------------------------- //
    for (i = 0; i < N; ++i) {
        connect[0] = i;
        connect[1] = (i+1) % N;
        mesh.addCell(ElementType::LINE, true, connect);
    } //next i
}

return; }

int subtest_002(
    void
) {

// ========================================================================== //
// int subtest_002(                                                           //
//     void)                                                                  //
//                                                                            //
// Test extraction of edge network from surface mesh for 2D curves.           //
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
SurfUnstructured                mesh(0, 1, 2);

// Counters
int                             nV, nS, nE;

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00004 - sub-test #002 - edge extraction (curve in 2D)       **" << endl;
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << endl;
}

// ========================================================================== //
// SET MESH ATTRIBUTES                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Set mesh attributes -------------------------------------------------- //
    mesh.setExpert(true);
}

// ========================================================================== //
// STEP#1 IMPORT MESH FROM STL FILE                                           //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Generate test surface mesh ------------------------------------------- //
    log::cout() << "** Generating surface mesh" << endl;
    Generate2DSurfMesh(mesh);
    log::cout() << "** Building adjacencies" << endl;
    mesh.buildAdjacencies();

    // Check mesh status ---------------------------------------------------- //
    nV = mesh.getVertexCount();
    nS = mesh.getInternalCount();
    nE = mesh.countFaces(); 
    if (nS != 32)   return 1;
    if (nV != 32)   return 1;
    if (nE != 32)   return 1;
    
}

// ========================================================================== //
// STEP#2 COMPUTE EDGE/VERTEX NORMALS                                         //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                                         n_vert;
    int                                         i;
    long                                        id;
    SurfUnstructured::CellIterator              c_, ce_;
    SurfUnstructured::VertexIterator            v_, ve_;
    long                                        cell_counter = 0, edge_counter = 0, vertex_counter = 0;
    unordered_map<long, long>                   vertex_mapper, edge_mapper;
    vector<array<double, 3>>                    normals, enormals, vnormals;

    // Compute cell normals ------------------------------------------------- //
    log::cout() << "** Computing cell normals" << endl;
    normals.resize(mesh.getCellCount());
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        normals[cell_counter] = mesh.evalFacetNormal(c_->getId());
        ++cell_counter;
    } //next c_

    // Compute edge normals ------------------------------------------------- //
    log::cout() << "** Computing edge normals" << endl;
    enormals.resize(mesh.countFaces());
    ce_ = mesh.cellEnd();
    edge_counter = 0;
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        enormals[edge_counter] = mesh.evalEdgeNormal(c_->getId(), 0);
        ++edge_counter;
    } //next c_

    // Compute vertex normals ----------------------------------------------- //
    log::cout() << "** Computing vertex normals" << endl;
    vnormals.resize(mesh.getVertexCount());
    ve_ = mesh.vertexEnd();
    for (v_ = mesh.vertexBegin(); v_ != ve_; ++v_) {
        vertex_mapper[v_->getId()] = vertex_counter;
        ++vertex_counter;
    } //next v_
    ce_ = mesh.cellEnd();
    for (c_ = mesh.cellBegin(); c_ != ce_; ++c_) {
        const long *cellConnect = c_->getConnect();
        n_vert = c_->getVertexCount();
        for (i = 0; i < n_vert; ++i) {
            id = cellConnect[i];
            vnormals[vertex_mapper[id]] = mesh.evalVertexNormal(c_->getId(), i);
        } //next i
    } //next c_

    // Export edges to vtu file --------------------------------------------- //
    {
        log::cout() << "** Exporting mesh to \"test_00004_subtest_002_surf.vtu\"" << endl;
        VTKUnstructuredGrid& vtk = mesh.getVTK() ;
        vtk.addData("normals", VTKFieldType::VECTOR, VTKLocation::CELL, normals);
        vtk.addData("enormals", VTKFieldType::VECTOR, VTKLocation::POINT, enormals);
        vtk.addData("vnormals", VTKFieldType::VECTOR, VTKLocation::POINT, vnormals);
        
        mesh.write("test_00004_subtest_002_curve");
        log::cout() << endl;
    }
}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    log::cout() << "** ================================================================= **" << endl;
    log::cout() << "** Test #00004 - sub-test #002 - completed!                          **" << endl;
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

        status = subtest_002();
        if (status != 0) {
            return (20 + status);
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
