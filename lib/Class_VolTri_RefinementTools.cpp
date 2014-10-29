// ========================================================================== //
//                         - Class_VolTri -                                   //
//                                                                            //
// Grid manager for unstructured volume meshes.                               //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v2.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "Class_VolTri.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// REFINEMENT TOOLS ========================================================= //

// -------------------------------------------------------------------------- //
void Class_VolTri::SwapFaceTri(
    int         T,
    int         i
) {

// ========================================================================== //
// void Class_VolTri::SwapFaceTri(                                            //
//     int         T,                                                         //
//     int         i)                                                         //
//                                                                            //
// Swap face between two neighboring triangles.                               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T    : int, triangle global index.                                       //
// - i    : int, face local index.                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int              j;
int              V0, V1, V2, P;
int              A, B, C, D, R;

// Counters
// none

// ========================================================================== //
// SWAP FACE                                                                  //
// ========================================================================== //

// Extract connectivity info for triangle T --------------------------------- //

// Vertices
V0 = Simplex[T][i];
V1 = Simplex[T][(i+1) % 3];
V2 = Simplex[T][(i+2) % 3];

// Adjacencies
R = Adjacency[T][i];
A = Adjacency[T][(i+1) % 3];
B = Adjacency[T][(i+2) % 3];


// Extract connectivity info for triangle R --------------------------------- //

// Vertexes
j = face(R, T);
P = Simplex[R][(j+2) % 3];

// Adjacencies
C = Adjacency[R][(j+1) % 3];
D = Adjacency[R][(j+2) % 3];

// Update triangle T -------------------------------------------------------- //

// Simplex-vertex connectivity
Simplex[T][0] = P;
Simplex[T][1] = V2;
Simplex[T][2] = V0;

// Adjacencies
Adjacency[T][0] = R;
Adjacency[T][1] = B;
Adjacency[T][2] = C;

// Update T-neighbors ------------------------------------------------------- //
if (C >= 0) {
    Adjacency[C][face(C, R)] = T;
}

// Update triangle R -------------------------------------------------------- //

// Simplex-vertex connectivity
Simplex[R][0] = P;
Simplex[R][1] = V1;
Simplex[R][2] = V2;

// Adjacencies
Adjacency[R][0] = D;
Adjacency[R][1] = A;
Adjacency[R][2] = T;

// Update R-neighbors ------------------------------------------------------- //
if (A >= 0) {
    Adjacency[A][face(A, T)] = R;
}

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::SplitTri3Tri(
    int         T,
    int         P
) {

// ========================================================================== //
// void Class_VolTri::SplitTri3Tri(                                           //
//     int         T,                                                         //
//     int         P)                                                         //
//                                                                            //
// Split triangle into 3 new triangles at a given point with global index P.  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T        : int, triangle global index                                    //
// - P        : int, point global index                                       //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
//  - none                                                                    //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int             A, B, C, V1, V2, V3;
ivector1D       idummy1D(3, -1);

// Counters
// none

// ========================================================================== //
// SPLIT TRIANGLE                                                             //
// ========================================================================== //

// Extract vertex data
V1 = Simplex[T][0];
V2 = Simplex[T][1];
V3 = Simplex[T][2];

// Extract adiacency data
A = Adjacency[T][0];
B = Adjacency[T][1];
C = Adjacency[T][2];

// 1st triangle (P, V1, V2)
Simplex[T][0] = P;
Simplex[T][1] = V1;
Simplex[T][2] = V2;
Adjacency[T][0] = nSimplex + 1;
Adjacency[T][1] = A;
Adjacency[T][2] = nSimplex;

// 2nd triangle (P, V2, V3)
idummy1D[0] = P;
idummy1D[1] = V2;
idummy1D[2] = V3;
AddSimplex(idummy1D, 5);
idummy1D[0] = T;
idummy1D[1] = B;
idummy1D[2] = nSimplex;
SetAdjacency(nSimplex-1, idummy1D);

// 3rd triangle (P, V3, V1)
idummy1D[0] = P;
idummy1D[1] = V3;
idummy1D[2] = V1;
AddSimplex(idummy1D, 5);
idummy1D[0] = nSimplex - 2;
idummy1D[1] = C;
idummy1D[2] = T;
SetAdjacency(nSimplex-1, idummy1D);

// =================================================================================== //
// UPDATE ADJACENCIES FOR NEIGHBORING TRIANGLES                                        //
// =================================================================================== //
if (B >= 0) {
    Adjacency[B][face(B, T)] = nSimplex - 2;
}
if (C >= 0) {
    Adjacency[C][face(C, T)] = nSimplex - 1;
}

return; };