// ========================================================================== //
//           ** BitPit mesh ** Test 003 for class SurfUnstructured **         //
//                                                                            //
// Test import output routines for class SurfUnstructured                     //
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
# include "bitpit_patchkernel.hpp"                                                  // BitPit base patch
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
// Test I/O routines for STL file format.                                     //
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
//               err = 2  --> error at step #2 (export STL)                   //
//               err = 3  --> error at step #3 (export VTK)                   //
// ========================================================================== //
    
// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
string                          in_name_bin = "./data/buddha.stl";
string                          in_name_ASCII = "./data/cube.stl";
string                          out_name_bin = "./buddha_copy.stl";
string                          out_name_ASCII = "./buddha_cube_copy.stl";
SurfUnstructured                mesh(0);

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
    cout << "** Test #00003 - sub-test #001 - I/O functions for STL format        **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

// ========================================================================== //
// STEP #1 (IMPORT MESH FROM BINARY STL)                                      //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    long                                nExpected = 283274;
    
    // Import mesh from stl format ------------------------------------------ //
    cout << "** Importing mesh from (binary): \"" << in_name_bin << "\"" << endl;
    if (mesh.importSTL(in_name_bin, true) > 0) return 1;
    if (mesh.getVertexCount() != 3*nExpected) return 1;
    if (mesh.getInternalCount() != nExpected) return 1;

}

// ========================================================================== //
// STEP #2 (EXPORT MESH TO BINARY STL)                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "** Exporting mesh to (binary): \"" << out_name_bin << "\"" << endl;
    if (mesh.exportSTL(out_name_bin, true) > 0) return 2;

}

// ========================================================================== //
// STEP #3 (IMPORT MESH FROM ASCII STL)                                       //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    long                                nExpected = 283286;
    
    // Import mesh from stl format ------------------------------------------ //
    cout << "** Appending mesh from (ASCII): \"" << in_name_ASCII << "\"" << endl;
    if (mesh.importSTL(in_name_ASCII, false) > 0) return 3;
    if (mesh.getVertexCount() != 3*nExpected) return 3;
    if (mesh.getInternalCount() != nExpected) return 3;

}

// ========================================================================== //
// STEP #4 (EXPORT MESH TO BINARY STL)                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "** Exporting mesh to (ASCII): \"" << out_name_ASCII << "\"" << endl;
    if (mesh.exportSTL(out_name_ASCII, false) > 0) return 4;

}

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables
    // none

    // Output message
    cout << "** ================================================================= **" << endl;
    cout << "** Test #00003 - sub-test #001 - completed!                          **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

return 0; }

// ========================================================================== //
// SUBTEST #001 Test DGF import/export ruotines                               //
// ========================================================================== //
int subtest_002(
    void
) {

// ========================================================================== //
// int subtest_002(                                                           //
//     void)                                                                  //
//                                                                            //
// Test DGF I/O functions.                                                    //
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
string                          in_name = "./data/NACA0012.dgf";
string                          out_name = "./NACA0012_copy.dgf";
string                          out_name_vtu = "./NACA0012_copy";
SurfUnstructured                mesh(0);

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
    cout << "** Test #00003 - sub-test #002 - I/O functions for DGF format        **" << endl;
    cout << "** ================================================================= **" << endl;
    cout << endl;
}

// ========================================================================== //
// STEP #1 (IMPORT MESH FROM DGF FILE)                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    long                                nVe = 998;
    long                                nSe = 998;
    
    // Import mesh from stl format ------------------------------------------ //
    cout << "** Importing mesh from : \"" << in_name << "\"" << endl;
    if (mesh.importDGF(in_name) > 0) return 1;
    if (mesh.getVertexCount() != nSe) return 1;
    if (mesh.getInternalCount() != nVe) return 1;

}

// ========================================================================== //
// STEP #2 (EXPORT MESH TO BINARY DGF)                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "** Exporting mesh to: \"" << out_name << "\"" << endl;
    if (mesh.exportDGF(out_name) > 0) return 2;
    mesh.write(out_name_vtu);

}

// ========================================================================== //
// STEP #3 (RE-IMPORT MESH)                                                   //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    SurfUnstructured                    mesh_copy(1);
    
    // Import mesh from stl format ------------------------------------------ //
    cout << "** Re-importing mesh from : \"" << out_name << "\"" << endl;
    if (mesh_copy.importDGF(in_name) > 0) return 3;
    if (mesh_copy.getVertexCount() != mesh.getVertexCount()) return 3;
    if (mesh_copy.getInternalCount() != mesh.getInternalCount()) return 3;

}


return 0;
}

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

// ========================================================================== //
// RUN SUB-TEST #002                                                          //
// ========================================================================== //
err = subtest_002();
if (err > 0) return(20 + err);


return err;

}
