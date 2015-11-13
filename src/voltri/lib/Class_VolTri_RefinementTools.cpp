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

// -------------------------------------------------------------------------- //
void Class_VolTri::SplitHexa5Tetra(
    int         T,
    short int   config
) {

// ========================================================================== //
// void Class_VolTri::SplitHexa5Tetra(                                        //
//     int         T,                                                         //
//     short int   config)                                                    //
//                                                                            //
// Split hexaedron into 5 tetrahedra.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T       : int, simplex global index                                      //
// - config  : short int, spliting configuration (0, or 1)                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                   update_adjacency = (Adjacency.size() == nSimplex);
int                    S = T;
ivector1D              s_hexa = Simplex[T];
ivector1D              s_tetra(4, -1);

// Counters
// none

// ========================================================================== //
// CHECK ELEMENT TYPE                                                         //
// ========================================================================== //
if (e_type[T] != 12) { return; }

// ========================================================================== //
// SPLIT HEXAEDRON ACCORDING TO CONFIGURATION                                 //
// ========================================================================== //
switch (config) {

    // Configuration #0
    case 0: {

        // tetra #0
        s_tetra[0] = s_hexa[0];
        s_tetra[1] = s_hexa[2];
        s_tetra[2] = s_hexa[3];
        s_tetra[3] = s_hexa[7];
        Simplex[T] = s_tetra;
        e_type[T] = 10;

        // tetra #1
        s_tetra[0] = s_hexa[0];
        s_tetra[1] = s_hexa[5];
        s_tetra[2] = s_hexa[7];
        s_tetra[3] = s_hexa[4];
        AddSimplex(s_tetra, 10);

        // tetra #2
        s_tetra[0] = s_hexa[5];
        s_tetra[1] = s_hexa[2];
        s_tetra[2] = s_hexa[7];
        s_tetra[3] = s_hexa[6];
        AddSimplex(s_tetra, 10);

        // tetra #3
        s_tetra[0] = s_hexa[0];
        s_tetra[1] = s_hexa[1];
        s_tetra[2] = s_hexa[2];
        s_tetra[3] = s_hexa[5];
        AddSimplex(s_tetra, 10);

        // tetra #4
        s_tetra[0] = s_hexa[0];
        s_tetra[1] = s_hexa[2];
        s_tetra[2] = s_hexa[5];
        s_tetra[3] = s_hexa[7];
        AddSimplex(s_tetra, 10);

        // Update adjacencies
        if (update_adjacency) {

            Adjacency[T].resize(infos[e_type[T]].n_faces);
            Adjacency[T][0] = -1;    Adjacency[T][1] = nSimplex-1;    Adjacency[T][2] = -1;    Adjacency[T][3] = -1;

            T = nSimplex-4;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = nSimplex-1;    Adjacency[T][1] = -1;    Adjacency[T][2] = -1;    Adjacency[T][3] = -1;

            T = nSimplex-3;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = nSimplex-1;    Adjacency[T][1] = -1;    Adjacency[T][2] = -1;    Adjacency[T][3] = -1;

            T = nSimplex-2;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = -1;    Adjacency[T][1] = -1;    Adjacency[T][2] = -1;    Adjacency[T][3] = nSimplex-1;

            T = nSimplex-1;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = nSimplex-2;    Adjacency[T][1] = S;    Adjacency[T][2] = nSimplex-3;    Adjacency[T][3] = nSimplex-4;
        }

        break;
    }

    // Configuration #1
    case 1: {

        // tetra #0
        s_tetra[0] = s_hexa[4];
        s_tetra[1] = s_hexa[7];
        s_tetra[2] = s_hexa[6];
        s_tetra[3] = s_hexa[3];
	Simplex[T] = s_tetra;
	e_type[T] = 10;

        // tetra #1
        s_tetra[0] = s_hexa[0];
        s_tetra[1] = s_hexa[1];
        s_tetra[2] = s_hexa[3];
        s_tetra[3] = s_hexa[4];
        AddSimplex(s_tetra, 10);

        // tetra #2
        s_tetra[0] = s_hexa[1];
        s_tetra[1] = s_hexa[2];
        s_tetra[2] = s_hexa[3];
        s_tetra[3] = s_hexa[6];
        AddSimplex(s_tetra, 10);

        // tetra #3
        s_tetra[0] = s_hexa[4];
        s_tetra[1] = s_hexa[6];
        s_tetra[2] = s_hexa[5];
        s_tetra[3] = s_hexa[1];
        AddSimplex(s_tetra, 10);

        // tetra #4
        s_tetra[0] = s_hexa[1];
        s_tetra[1] = s_hexa[4];
        s_tetra[2] = s_hexa[3];
        s_tetra[3] = s_hexa[6];
        AddSimplex(s_tetra, 10);

        // Update adjacencies
        if (update_adjacency) {

            Adjacency[T].resize(infos[e_type[T]].n_faces);
            Adjacency[T][0] = -1;    Adjacency[T][1] = -1;    Adjacency[T][2] = -1;    Adjacency[T][3] = nSimplex-1;

            T = nSimplex-4;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = -1;    Adjacency[T][1] = -1;    Adjacency[T][2] = nSimplex-1;    Adjacency[T][3] = -1;

            T = nSimplex-3;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = -1;    Adjacency[T][1] = -1;    Adjacency[T][2] = -1;    Adjacency[T][3] = nSimplex-1;

            T = nSimplex-2;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = -1;    Adjacency[T][1] = nSimplex-1;    Adjacency[T][2] = -1;    Adjacency[T][3] = -1;

            T = nSimplex-1;
            SetAdjacency(T, s_tetra);
            Adjacency[T][0] = nSimplex-4;    Adjacency[T][1] = nSimplex-2;    Adjacency[T][2] = S;    Adjacency[T][3] = nSimplex-3;
        }

        break;
    }
}

return; }

// -------------------------------------------------------------------------- //
void Class_VolTri::SplitPyra2Tetra(
    int         T,
    short int   config
) {

// ========================================================================== //
// void Class_VolTri::SplitPyra2Tetra(                                        //
//     int         T,                                                         //
//     short int   config)                                                    //
//                                                                            //
// Split pyramid into 2 tetrahedra.                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T       : int, simplex global index                                      //
// - config  : short int, spliting configuration (0, or 1)                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                   update_adjacency = (Adjacency.size() == nSimplex);
ivector1D              adjacencies = Adjacency[T];
ivector1D              s_pyra = Simplex[T];
ivector1D              s_tetra(4, -1);

// Counters
// none

// ========================================================================== //
// CHECK ELEMENT TYPE                                                         //
// ========================================================================== //
if (e_type[T] != 14) { return; }

// ========================================================================== //
// SPLIT HEXAEDRON ACCORDING TO CONFIGURATION                                 //
// ========================================================================== //
switch (config) {

    // Configuration #0
    case 0: {

        // tetra #0
        s_tetra[0] = s_pyra[0];
        s_tetra[1] = s_pyra[1];
        s_tetra[2] = s_pyra[2];
        s_tetra[3] = s_pyra[4];
        Simplex[T] = s_tetra;
        e_type[T] = 10;

        // Udate adjacency for tetra#0
        if (update_adjacency) {
            Adjacency[T].resize(infos[e_type[T]].n_faces);
            Adjacency[T][0] = -1;
            Adjacency[T][1] = adjacencies[1];
            Adjacency[T][2] = adjacencies[2];
            Adjacency[T][3] = nSimplex;
        }

        // tetra #1
        s_tetra[0] = s_pyra[0];
        s_tetra[1] = s_pyra[2];
        s_tetra[2] = s_pyra[3];
        s_tetra[3] = s_pyra[4];
        AddSimplex(s_tetra, 10);

        // Update adjacency for tetra#1
        if (update_adjacency) {
            SetAdjacency(nSimplex-1, s_tetra);
            Adjacency[nSimplex-1][0] = -1;
            Adjacency[nSimplex-1][1] = T;
            Adjacency[nSimplex-1][2] = adjacencies[3];
            Adjacency[nSimplex-1][3] = adjacencies[4];
        }

        break;
    }

    // Configuration #1
    case 1: {

        // tetra #0
        s_tetra[0] = s_pyra[0];
        s_tetra[1] = s_pyra[1];
        s_tetra[2] = s_pyra[3];
        s_tetra[3] = s_pyra[4];
        Simplex[T] = s_tetra;
        e_type[T] = 10;

        // Udate adjacency for tetra#0
        if (update_adjacency) {
            Adjacency[T].resize(infos[e_type[T]].n_faces);
            Adjacency[T][0] = -1;
            Adjacency[T][1] = adjacencies[1];
            Adjacency[T][2] = nSimplex;
            Adjacency[T][3] = adjacencies[4];
        }

        // tetra #1
        s_tetra[0] = s_pyra[1];
        s_tetra[1] = s_pyra[2];
        s_tetra[2] = s_pyra[3];
        s_tetra[3] = s_pyra[4];
        AddSimplex(s_tetra, 10);

        // Update adjacency for tetra#1
        if (update_adjacency) {
            SetAdjacency(nSimplex-1, s_tetra);
            Adjacency[nSimplex-1][0] = -1;
            Adjacency[nSimplex-1][1] = adjacencies[2];
            Adjacency[nSimplex-1][2] = adjacencies[3];
            Adjacency[nSimplex-1][3] = T;
        }

        break;
    }
}

return; }


