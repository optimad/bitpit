// ========================================================================== //
//                           - STL IO FUNCTIONS -                             //
//                                                                            //
// Templated functions for STL format.                                        //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v3.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
// none

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS                                                   //
// ========================================================================== //

// -------------------------------------------------------------------------- //
template <typename ... T2>
void STL_obj::load(
    string           sname,
    int             &nV,
    int             &nT,
    dvector2D       &V,
    dvector2D       &N,
    ivector2D       &T,
    T2        & ... others
) {

// ========================================================================== //
// template <typename ... T2>                                                 //
// void STL_obj::load(                                                        //
//     string           sname,                                                //
//     int             &nV,                                                   //
//     int             &nT,                                                   //
//     dvector2D       &V,                                                    //
//     dvector2D       &N,                                                    //
//     ivector2D       &T,                                                    //
//     T2        & ... others)                                                //
//                                                                            //
// Load stl with name specified in 'sname' from ASCII stl file. If no solid   //
// is found, input data structure is unchanged. If solid name is not          //
// specified (i.e. sname = ""), then returns the first solid found by         //
// circular scanning of the stl file from the current cursor position.        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - sname     : string, stl solid name                                       //
// - nV        : int, number of stl vertices                                  //
// - nT        : int, number of stl facets                                    //
// - V         : dvector2D, vertex coordinate list. V[i][0], V[i][1], and     //
//               V[i][2] are the x, y, z coordinates of the i-th vertex in    //
//               stl triangulation                                            //
// - N         : dvector2D, triangles unit normal. N[i][0], N[i][1], and      //
//               N[i][2] are the x, y, z components of the i-th vertex in     //
//               stl triangulation                                            //
// - T         : ivector2D, triangle-vertex connectivity. T[i][0], T[i][1],   //
//               and T[i][2] are global indices of vertices of the i-th       //
//               triangle in the stl triangulation.                           //
// - others    : T2 & ..., other stl solids                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");
if (err == 1) {
    return;
}

// ========================================================================== //
// LOAD STL SOLID                                                             //
// ========================================================================== //

// Load solid data
if (stl_type) { Read_STL_bin(ifile_handle, nV, nT, V, N, T); }
else          { Read_STLsolid_ASCII(ifile_handle, nV, nT, V, N, T, sname); }

// ========================================================================== //
// ITERATIVELY READ STL SOLIDS                                                //
// ========================================================================== //
load(others ...);

return; }

// -------------------------------------------------------------------------- //
template <typename ... T2>
void STL_obj::save(
    string           sname,
    int             &nV,
    int             &nT,
    dvector2D       &V,
    dvector2D       &N,
    ivector2D       &T,
    T2        & ... others
) {

// ========================================================================== //
// template <typename ... T2>                                                 //
// void STL_obj::save(                                                        //
//     string           sname,                                                //
//     int             &nV,                                                   //
//     int             &nT,                                                   //
//     dvector2D       &V,                                                    //
//     dvector2D       &N,                                                    //
//     ivector2D       &T,                                                    //
//     T2        & ... others)                                                //
//                                                                            //
// Export stl solid data to stl file.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - sname             : string (optional) stl solid name                     //
// - nV                : int, number of stl vertices                          //
// - nT                : int, number of stl facets                            //
// - V                 : dvector2D, vertex coordinate list. V[i][0], V[i][1]  //
//                       and V[i][2] are the x, y, z coordinates of the i-th  //
//                       vertex in the stl triangulation                      //
// - N                 : dvector2D, unit normals to each triangle. N[i][0],   //
//                       N[i][1] and N[i][2] are the x, y, z components of    //
//                       the normal unit vector to the i-th triangle.         //
// - T                 : ivector2D, triangle-vertex connectivity. S[i][0]     //
//                       S[i][1] and S[i][2] are the global indices of the    //
//                       i-th triangle in the stl triangulation.              //
// - others            : other stl solids.                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("out");
if (err == 1) { return; }

// ========================================================================== //
// EXPORT STL SOLID                                                           //
// ========================================================================== //

// Export solid
if (stl_type) { Write_STLsolid_bin(ofile_handle, nV, nT, V, N, T, sname); }
else          { Write_STLsolid_ASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
template <typename ... T2>
void STL_obj::append(
    string           sname,
    int             &nV,
    int             &nT,
    dvector2D       &V,
    dvector2D       &N,
    ivector2D       &T,
    T2        & ... others
) {

// ========================================================================== //
// template <typename ... T2>                                                 //
// void STL_obj::append(                                                      //
//     string           sname,                                                //
//     int             &nV,                                                   //
//     int             &nT,                                                   //
//     dvector2D       &V,                                                    //
//     dvector2D       &N,                                                    //
//     ivector2D       &T,                                                    //
//     T2        & ... others)                                                //
//                                                                            //
// Append stl solid data to stl file.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - sname             : string (optional) stl solid name                     //
// - nV                : int, number of stl vertices                          //
// - nT                : int, number of stl facets                            //
// - V                 : dvector2D, vertex coordinate list. V[i][0], V[i][1]  //
//                       and V[i][2] are the x, y, z coordinates of the i-th  //
//                       vertex in the stl triangulation                      //
// - N                 : dvector2D, unit normals to each triangle. N[i][0],   //
//                       N[i][1] and N[i][2] are the x, y, z components of    //
//                       the normal unit vector to the i-th triangle.         //
// - T                 : ivector2D, triangle-vertex connectivity. S[i][0]     //
//                       S[i][1] and S[i][2] are the global indices of the    //
//                       i-th triangle in the stl triangulation.              //
// - T2                : other stl solids.                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("app");
if (err == 1) { return; }

// ========================================================================== //
// EXPORT STL SOLID                                                           //
// ========================================================================== //

// Export solid
if (stl_type) { Write_STLsolid_bin(ofile_handle, nV, nT, V, N, T, sname); }
else          { Write_STLsolid_ASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }


