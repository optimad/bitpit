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

// TOPOLOGY ================================================================= //

// -------------------------------------------------------------------------- //
bool Class_VolTri::SameFace(
    int             A,
    int             i,
    int             B,
    int             j
) {

// ========================================================================== //
// bool Class_VolTri::SameFace(                                               //
//     int             A,                                                     //
//     int             i,                                                     //
//     int             B,                                                     //
//     int             j)                                                     //
//                                                                            //
// Check for coincident faces.                                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A     : int, 1st simplex global index                                    //
// - i     : int, face local index on simplex A                               //
// - B     : int, 2nd simplex global index                                    //
// - j     : int, face local index on simplex B                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, returns true if face (A, i) and (B, j) are coincident      //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check = false;
ivector1D           face_vlistA, face_vlistB;

// ========================================================================== //
// CHECK FOR COINCIDENT FACES                                                 //
// ========================================================================== //
face_vlistA = FaceVertices(A, i);
face_vlistB = FaceVertices(B, j);
if (face_vlistA.size() == face_vlistB.size()) {
    sort(face_vlistA.begin(), face_vlistA.end());
    sort(face_vlistB.begin(), face_vlistB.end());
    check = (face_vlistA == face_vlistB);
}

return(check); };

// -------------------------------------------------------------------------- //
bool Class_VolTri::SameSimplex(
    int             A,
    int             B
) {

// ========================================================================== //
// bool Class_VolTri::SameSimplex(                                            //
//     int             A,                                                     //
//     int             B)                                                     //
//                                                                            //
// Check if two simplicies are coincident.                                    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A     : int, 1st simplex global index                                    //
// - B     : int, 2nd simplex global index                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, true if A and B are coincident, false otherwise            //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool        check = false;

// Counters
// none

// ========================================================================== //
// CHECK FOR COINCIDENT SIMPLICIES                                            //
// ========================================================================== //
if (e_type[A] == e_type[B]) {

    // Scope variables
    ivector1D       vlist_A(infos[e_type[A]].n_vert);
    ivector1D       vlist_B(infos[e_type[B]].n_vert);

    // Sort vertex list in ascending order
    vlist_A = Simplex[A];
    sort(vlist_A.begin(), vlist_A.end());
    vlist_B = Simplex[B];
    sort(vlist_B.begin(), vlist_B.end());
    check = (vlist_A == vlist_B);
    
}

return(check); };

// -------------------------------------------------------------------------- //
ivector1D Class_VolTri::FaceVertices(
    int             S,
    int             i
) {

// ========================================================================== //
// ivector1D Class_VolTri::FaceVertices(                                      //
//     int             S,                                                     //
//     int             i)                                                     //
//                                                                            //
// Return list of global indices of vertices for the i-th face of simplex S.  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - S     : int, simplex global index                                        //
// - i     : int, face local index                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - vlist : ivector1D, list of global indices of vertices of the i-th face   //
//           of simplex S.                                                    //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                 n = infos[e_type[S]].faces[i].size();
ivector1D           vlist(n, -1);

// Counters
int                 j;

// ========================================================================== //
// RETURN GLOBAL INDEX OF VERTICES                                            //
// ========================================================================== //
for (j = 0; j < n; ++j) {
    vlist[j] = Simplex[S][infos[e_type[S]].faces[i][j]];
} //next i

return(vlist); };

// GEOMETRY ================================================================= //

// -------------------------------------------------------------------------- //
a3vector1D Class_VolTri::Baricenter(
    int             T
) {

// ========================================================================== //
// a3vector1D Class_VolTri::Baricenter(                                       //
//     int             T)                                                     //
//                                                                            //
// Compute simplex baricenter.                                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T       : int, simplex global index                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - P       : a3vector1D, with simplex baricenter coordinates.               //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                     n = infos[e_type[T]].n_vert;
a3vector1D              P;

// Counters
int                     i, j;

// ========================================================================== //
// FIND ISOLATED VERTEXES                                                     //
// ========================================================================== //
P.fill(0.0);
for (i = 0; i < n; ++i) {
    for (j = 0; j < 3; ++j) {
        P[j] += Vertex[Simplex[T][i]][j];
    } //next j
} //next i
for (j = 0; j < 3; ++j) {
    P[j] = P[j]/((double) n);
} //next j

return(P); };

// -------------------------------------------------------------------------- //
a3vector1D Class_VolTri::FaceCenter(
    int             T,
    int             i
) {

// ========================================================================== //
// a3vector1D Class_VolTri::FaceCenter(                                       //
//     int             T,                                                     //
//     int             i)                                                     //
//                                                                            //
// Returns baricenter of the i-th face of the T-th simplex.                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex globa index                                        //
// - i      : int, face local index                                           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - C      : a3vector1D, with face center coordinates                        //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int              n;
ivector1D        f_vert;
a3vector1D       C;

// Counters
int              V;
int              j, k;

// ========================================================================== //
// COMPUTE FACE CENTER                                                        //
// ========================================================================== //

// Get face vertex list
C.fill(0.0);
f_vert = FaceVertices(T, i);
n = f_vert.size();
for (j = 0; j < n; ++j) {
    V = f_vert[j];
    for (k = 0; k < 3; ++k) {
        C[k] += Vertex[V][k];
    } //next k
} //next j
for (j = 0; j < 3; ++j) {
    C[j] = C[j]/((double) n);
} //next j

return(C); };

// -------------------------------------------------------------------------- //
a3vector1D Class_VolTri::CircumCenter(
    int             T
) {

// ========================================================================== //
// a3vector1D Class_VolTri::CircumCenter(                                     //
//     int             T)                                                     //
//                                                                            //
// Compute circum center of simplex T.                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
int     const   iter_max = 10;
double  const   abs_toll = 1.0e-6;

// Local variables
bool            check;
a3vector1D      C;
dvector1D       dC;
dvector1D       d, r;
dvector1D       dG;
dvector2D       dr;
dvector3D       d2r;
dvector2D       d2G;

// Counters
int             i, j, k;
int             n, m;
int             U, V;
int             iter;

// ========================================================================== //
// SET PARAMETERS                                                             //
// ========================================================================== //

// Linear system dimensions
if ((e_type[T] == 5) || (e_type[T] == 9)) {
    n = 2;
}
else {
    n = 3;
}
dG.resize(n, 0.0);
d2G.resize(n, dvector1D(n, 0.0));
dC.resize(n, 0.0);

// number of simplex vertices
m = infos[e_type[T]].n_vert;
d.resize(m, 0.0);
r.resize(m, 0.0);
dr.resize(m, dvector1D(n, 0.0));
d2r.resize(m, dvector2D(n, dvector1D(n, 0.0)));

// ========================================================================== //
// ITERATIVE PROCEDURE                                                        //
// ========================================================================== //

// Initial estimate for circum center --------------------------------------- //
C = Baricenter(T);

// Newton-Raphson iterations ------------------------------------------------ //
iter = 0;
check = true;
U = Simplex[T][0];
while (check && (iter < iter_max)) {

    // Reset r.h.s. and coeffs. matrix
    for (i = 0; i < n; ++i) {
        dG[i] = 0.0;
        for (j = 0; j < n; ++j) {
            d2G[i][j] = 0.0;
        } //next j
    } //next i

    // Update distances
    for (i = 0; i < m; ++i) {
        V = Simplex[T][i];
        d[i] = norm2(C - Vertex[V]);
    } //next i

    // Update residuals
    for (i = 1; i < m; ++i) {
        r[i] = d[i] - d[0];
    } //next i

    // Update gradient of residuals
    for (i = 1; i < m; ++i) {
        V = Simplex[T][i];
        for (j = 0; j < 3; ++j) {
            dr[i][j] = (C[j] - Vertex[V][j])/d[i] - (C[j] - Vertex[U][j])/d[0];
        } //next j
    } //next i

    // Update Hessians of residuals
    for (i = 1; i < m; ++i) {
        V = Simplex[T][i];
        for (j = 0; j < n; ++j) {
            for (k = 0; k < n; ++k) {
                d2r[i][j][k] = -(C[j] - Vertex[V][j]) * (C[k] - Vertex[V][k]) / (d[i]*d[i]*d[i])
                               +(C[j] - Vertex[U][j]) * (C[k] - Vertex[U][k]) / (d[0]*d[0]*d[0]);
            } //next k
        } //next j
        for (j = 0; j < n; ++j) {
            d2r[i][j][j] += 1.0/d[i] - 1.0/d[0];
        } //next j
    }

    // r.h.s
    for (i = 1; i < m; ++i) {
        V = Simplex[T][i];
        dG = dG + 2.0 * r[i] * dr[i];
    } //next i
    dG = -1.0*dG;

    // coeff. matrix
    for (i = 0; i < n; ++i) {    
        for (j = 0; j < n; ++j) {
            for (k = 1; k < m; ++k) {
                d2G[i][j] += 2.0 * (dr[k][i]*dr[k][j] + r[i] * d2r[k][i][j]);
            } //next k
        } //next j
    } //next i

    // Solve linear system
    bitpit::linearalgebra::cramer(d2G, dG, dC);

    // Convergence criterion
    check = (norm2(dC) > abs_toll);
    iter++;

    // Update solution
    for (i = 0; i < n; ++i) {
        C[i] += dC[i];
    } //next i

} //next iter

return(C); };

// -------------------------------------------------------------------------- //
double Class_VolTri::EdgeLength(
    int             T,
    int             i
) {

// ========================================================================== //
// double Class_VolTri::EdgeLength(                                           //
//     int             T,                                                     //
//     int             i)                                                     //
//                                                                            //
// Compute length of a specified edge on a given simplex.                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T       : int, simplex global index                                      //
// - i       : int, edge local index                                          //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - l       : double, edge length                                            //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              l;

// Counters
int                 A, B;
int                 n;

// ========================================================================== //
// COMPUTE EDGE LENGTH                                                        //
// ========================================================================== //
n = infos[e_type[T]].edges[i].size();
if (n <= 1) { return(0.0); }
else {
    A = Simplex[T][infos[e_type[T]].edges[i][0]];
    B = Simplex[T][infos[e_type[T]].edges[i][1]];
    l = norm2(Vertex[A] - Vertex[B]);
}

return(l); };

// -------------------------------------------------------------------------- //
double Class_VolTri::EdgeLength(
    int             T,
    int             i,
    a3vector2D     &V
) {

// ========================================================================== //
// double Class_VolTri::EdgeLength(                                           //
//     int             T,                                                     //
//     int             i,                                                     //
//     a3vector2D     &V)                                                     //
//                                                                            //
// Compute length of a specified edge on a given simplex. Vertex list is      //
// provided externally.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T       : int, simplex global index                                      //
// - i       : int, edge local index                                          //
// - V       : a3vector2D, external vertex list. V[i][0], V[i][1], ... are    //
//             the x, y, ... coordinates of the i-th vertex.                  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - l       : double, edge length                                            //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              l;

// Counters
int                 A, B;
int                 n;

// ========================================================================== //
// COMPUTE EDGE LENGTH                                                        //
// ========================================================================== //
n = infos[e_type[T]].edges[i].size();
if (n <= 1) { return(0.0); }
else {
    A = Simplex[T][infos[e_type[T]].edges[i][0]];
    B = Simplex[T][infos[e_type[T]].edges[i][1]];
    l = norm2(V[A] - V[B]);
}

return(l); };

// -------------------------------------------------------------------------- //
double Class_VolTri::minEdgeLength(
    int             T
) {

// ========================================================================== //
// double Class_VolTri::minEdgeLength(                                        //
//     int             T)                                                     //
//                                                                            //
// Compute edge of minimal length for a given simplex.                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - min_l  : double, minimal edge length                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double          min_l;

// Counters
int             i, n;

// ========================================================================== //
// FIND EDGE OF MINIMIAL LENGTH                                               //
// ========================================================================== //
n = infos[e_type[T]].n_edges;
min_l = 1.0e+18;
for (i = 0; i < n; ++i) {
    min_l = min(min_l, EdgeLength(T, i));
} //next i

return(min_l); };

// -------------------------------------------------------------------------- //
double Class_VolTri::minEdgeLength(
    int             T,
    a3vector2D     &V
) {

// ========================================================================== //
// double Class_VolTri::minEdgeLength(                                        //
//     int             T,                                                     //
//     a3vector2D      &V)                                                    //
//                                                                            //
// Compute edge of minimal length for a given simplex. Vertex list is         //
// provided externally.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// - V      : a3vector2D, external vertex list. V[i][0], V[i][1], ... are     //
//            the x, y, ... coordinates of the i-th vertex.                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - min_l  : double, minimal edge length                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double          min_l;

// Counters
int             i, n;

// ========================================================================== //
// FIND EDGE OF MINIMIAL LENGTH                                               //
// ========================================================================== //
n = infos[e_type[T]].n_edges;
min_l = 1.0e+18;
for (i = 0; i < n; ++i) {
    min_l = min(min_l, EdgeLength(T, i, V));
} //next i

return(min_l); };

// -------------------------------------------------------------------------- //
double Class_VolTri::maxEdgeLength(
    int             T
) {

// ========================================================================== //
// double Class_VolTri::maxEdgeLength(                                        //
//     int             T)                                                     //
//                                                                            //
// Compute edge of maximal length for a given simplex.                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - max_l  : double, maximal edge length                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double          max_l;

// Counters
int             i, n;

// ========================================================================== //
// FIND EDGE OF MINIMIAL LENGTH                                               //
// ========================================================================== //
n = infos[e_type[T]].n_edges;
max_l = 0.0;
for (i = 0; i < n; ++i) {
    max_l = max(max_l, EdgeLength(T, i));
} //next i

return(max_l); };

// -------------------------------------------------------------------------- //
double Class_VolTri::maxEdgeLength(
    int             T,
    a3vector2D     &V
) {

// ========================================================================== //
// double Class_VolTri::maxEdgeLength(                                        //
//     int             T,                                                     //
//     a3vector2D     &V)                                                     //
//                                                                            //
// Compute edge of maximal length for a given simplex. Vertex list is         //
// provided externally.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// - V      : a3vector2D, external vertex list. V[i][0], V[i][1], ... are     //
//            the x, y, ... coordinates of the i-th vertex.                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - max_l  : double, maximal edge length                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double          max_l;

// Counters
int             i, n;

// ========================================================================== //
// FIND EDGE OF MINIMIAL LENGTH                                               //
// ========================================================================== //
n = infos[e_type[T]].n_edges;
max_l = 0.0;
for (i = 0; i < n; ++i) {
    max_l = max(max_l, EdgeLength(T, i));
} //next i

return(max_l); };

// -------------------------------------------------------------------------- //
double Class_VolTri::FaceArea(
    int             T,
    int             i
) {

// ========================================================================== //
// double Class_VolTri::FaceArea(                                             //
//     int             T,                                                     //
//     int             i)                                                     //
//                                                                            //
// Compute area of the i-th face of the T-th simplex.                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T    : int, simplex global index                                         //
// - i    : int, face local index.                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - A    : double, face area                                                 //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              a, A = 0.0;
array<double, 3>    u, v;

// Counters
int                 U, V0, V1;
int                 j, k, l, m;
int                 n;

// ========================================================================== //
// COMPUTE FACE AREA                                                          //
// ========================================================================== //
n = infos[e_type[T]].faces[i].size();
if (n == 2) {
    V0 = Simplex[T][infos[e_type[T]].faces[i][0]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][1]];
    A = norm2(Vertex[V0] - Vertex[V1]);
}
else if (n == 3) {
    U = Simplex[T][infos[e_type[T]].faces[i][0]];
    V0 = Simplex[T][infos[e_type[T]].faces[i][1]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][2]];
    for (k = 0; k < 3; ++k) {
        u[k] = Vertex[V0][k] - Vertex[U][k];
        v[k] = Vertex[V1][k] - Vertex[U][k];
    } //next k
    A = 0.5 * norm2(crossProduct(v, u));
}
else {
    A = 0.0;
    for (j = 0; j < n; ++j) {
        U = Simplex[T][infos[e_type[T]].faces[i][j]];
        a = 0.0;
        k = 0;
        l = (j+1)%n;
        while (k < n-2) {
            V0 = Simplex[T][infos[e_type[T]].faces[i][l]];
            V1 = Simplex[T][infos[e_type[T]].faces[i][(l+1)%n]];
            for (m = 0; m < 3; ++m) {
                u[m] = Vertex[V0][m] - Vertex[U][m];
                v[m] = Vertex[V1][m] - Vertex[U][m];
            } //next m
            a += 0.5 * norm2(crossProduct(v, u));
            k++;
            l = (l+1)%n;
        } //next k
        A += a/((double) (n-2));
    } //next j
    A = A/((double) n);
}

return(A); };

// -------------------------------------------------------------------------- //
double Class_VolTri::FaceArea(
    int             T,
    int             i,
    a3vector2D     &V
) {

// ========================================================================== //
// double Class_VolTri::FaceArea(                                             //
//     int             T,                                                     //
//     int             i,                                                     //
//     a3vector2D     &V)                                                     //
//                                                                            //
// Compute area of the i-th face of the T-th simplex. Vertex list is provided //
// externally.                                                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T    : int, simplex global index                                         //
// - i    : int, face local index.                                            //
// - V    : a3vector2D, external vertex list. V[i][0], V[i][1], ... are       //
//          the x, y, ... coordinates of the i-th vertex.                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - A    : double, face area                                                 //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              a, A = 0.0;
array<double, 3>    u, v;

// Counters
int                 U, V0, V1;
int                 j, k, l, m;
int                 n;

// ========================================================================== //
// COMPUTE FACE AREA                                                          //
// ========================================================================== //
n = infos[e_type[T]].faces[i].size();
if (n == 2) {
    V0 = Simplex[T][infos[e_type[T]].faces[i][0]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][1]];
    A = norm2(V[V0] - V[V1]);
}
else if (n == 3) {
    U = Simplex[T][infos[e_type[T]].faces[i][0]];
    V0 = Simplex[T][infos[e_type[T]].faces[i][1]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][2]];
    for (k = 0; k < 3; ++k) {
        u[k] = V[V0][k] - V[U][k];
        v[k] = V[V1][k] - V[U][k];
    } //next k
    A = 0.5 * norm2(crossProduct(v, u));
}
else {
    A = 0.0;
    for (j = 0; j < n; ++j) {
        U = Simplex[T][infos[e_type[T]].faces[i][j]];
        a = 0.0;
        k = 0;
        l = (j+1)%n;
        while (k < n-2) {
            V0 = Simplex[T][infos[e_type[T]].faces[i][l]];
            V1 = Simplex[T][infos[e_type[T]].faces[i][(l+1)%n]];
            for (m = 0; m < 3; ++m) {
                u[k] = V[V0][m] - V[U][m];
                v[k] = V[V1][m] - V[U][m];
            } //next m
            a += 0.5 * norm2(crossProduct(v, u));
            k++;
            l = (l+1)%n;
        } //next k
        A += a/((double) (n-2));
    } //next j
    A = A/((double) n);
}

return(A); };

// -------------------------------------------------------------------------- //
double Class_VolTri::minFaceArea(
    int             T
) {

// ========================================================================== //
// double Class_VolTri::minFaceArea(                                          //
//     int             T)                                                     //
//                                                                            //
// Compute minimal face area of a given simplex.                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex globa index                                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - A      : double, min face area                                           //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              A;

// Counters
int                 n;
int                 i;

// ========================================================================== //
// FIND MIN FACE AREA                                                         //
// ========================================================================== //
n = infos[e_type[T]].n_faces;
A = 1.0e+18;
for (i = 0; i < n; ++i) {
    A = min(A, FaceArea(T, i));
} //next i

return(A); };

// -------------------------------------------------------------------------- //
double Class_VolTri::minFaceArea(
    int             T,
    a3vector2D     &V
) {

// ========================================================================== //
// double Class_VolTri::minFaceArea(                                          //
//     int             T,                                                     //
//     a3vector2D     &V)                                                     //
//                                                                            //
// Compute minimal face area of a given simplex. Vertex coordinate list is    //
// provided externally.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T    : int, simplex globa index                                          //
// - V    : a3vector2D, external vertex list. V[i][0], V[i][1], ... are       //
//          the x, y, ... coordinates of the i-th vertex.                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - A      : double, min face area                                           //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              A;

// Counters
int                 n;
int                 i;

// ========================================================================== //
// FIND MIN FACE AREA                                                         //
// ========================================================================== //
n = infos[e_type[T]].n_faces;
A = 1.0e+18;
for (i = 0; i < n; ++i) {
    A = min(A, FaceArea(T, i, V));
} //next i

return(A); };

// -------------------------------------------------------------------------- //
double Class_VolTri::maxFaceArea(
    int             T
) {

// ========================================================================== //
// double Class_VolTri::maxFaceArea(                                          //
//     int             T)                                                     //
//                                                                            //
// Compute maximal face area of a given simplex.                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex globa index                                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - A      : double, max face area                                           //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              A;

// Counters
int                 n;
int                 i;

// ========================================================================== //
// FIND MIN FACE AREA                                                         //
// ========================================================================== //
n = infos[e_type[T]].n_faces;
A = 0.0;
for (i = 0; i < n; ++i) {
    A = max(A, FaceArea(T, i));
} //next i

return(A); };

// -------------------------------------------------------------------------- //
double Class_VolTri::maxFaceArea(
    int             T,
    a3vector2D     &V
) {

// ========================================================================== //
// double Class_VolTri::maxFaceArea(                                          //
//     int             T,                                                     //
//     a3vector2D     &V)                                                     //
//                                                                            //
// Compute maximal face area of a given simplex. Vertex coordinate list is    //
// provided externally.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T    : int, simplex globa index                                          //
// - V    : dvector2D, external vertex list. V[i][0], V[i][1], ... are        //
//          the x, y, ... coordinates of the i-th vertex.                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - A      : double, max face area                                           //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double              A;

// Counters
int                 n;
int                 i;

// ========================================================================== //
// FIND MIN FACE AREA                                                         //
// ========================================================================== //
n = infos[e_type[T]].n_faces;
A = 0.0;
for (i = 0; i < n; ++i) {
    A = max(A, FaceArea(T, i, V));
} //next i

return(A); };

// -------------------------------------------------------------------------- //
a3vector1D Class_VolTri::FaceNormal(
    int             T,
    int             i
) {

// ========================================================================== //
// a3vector1D Class_VolTri::FaceNormal(                                       //
//     int             T,                                                     //
//     int             i)                                                     //
//                                                                            //
// Compute face normal.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T     : int, simplex global index                                        //
// - i     : int, face local index                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - N     : a3vector1D, face normal                                          //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
array<double, 3>    u, v, w, t;
a3vector1D          N;

// Counters
int                 U, V0, V1;
int                 j, k, l, m;
int                 n;

// ========================================================================== //
// COMPUTE FACE NORMAL                                                        //
// ========================================================================== //
N.fill(0.0);
n = infos[e_type[T]].faces[i].size();
if (n <= 1) { return(N); }
else if (n == 2) {
    V0 = Simplex[T][infos[e_type[T]].faces[i][0]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][1]];
    u = Vertex[V1] - Vertex[V0];
    v[0] = 0.0;    v[1] = 0.0;    v[2] = 1.0;
    w = crossProduct(u, v);
    N = w/norm2(w);
}
else if (n == 3) {
    U = Simplex[T][infos[e_type[T]].faces[i][0]];
    V0 = Simplex[T][infos[e_type[T]].faces[i][1]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][2]];
    u = Vertex[V0] - Vertex[U];
    v = Vertex[V1] - Vertex[U];
    t = crossProduct(u, v);
    N = t/norm2(t);
}
else {
    for (j = 0; j < n; ++j) {
        U = Simplex[T][infos[e_type[T]].faces[i][j]];
        t.fill(0.0);
        k = 0;
        l = (j+1) % n;
        while (k < n+2) {
            V0 = Simplex[T][infos[e_type[T]].faces[i][l]];
            V1 = Simplex[T][infos[e_type[T]].faces[i][(l+1)%n]];
            u = Vertex[V0] - Vertex[U];
            v = Vertex[V1] - Vertex[U];
            w = crossProduct(u, v);
            t = t + w/norm2(w);
            k++;
        } //next k
        t = t/((double) (n-2));
        for (k = 0; k < 3; ++k) {
            N[k] += t[k];
        } //next k
    } //next j
    N = N/norm2(N);
}

return(N); };

// -------------------------------------------------------------------------- //
a3vector1D Class_VolTri::FaceNormal(
    int             T,
    int             i,
    a3vector2D     &V
) {

// ========================================================================== //
// a3vector1D Class_VolTri::FaceNormal(                                       //
//     int             T,                                                     //
//     int             i,                                                     //
//     a3vector2D     &V)                                                     //
//                                                                            //
// Compute face normal.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T     : int, simplex global index                                        //
// - i     : int, face local index                                            //
// - V     : a3vector2D, external vertex list. V[i][0], V[i][1], ... are      //
//           the x, y, ... coordinates of the i-th vertex.                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - N     : dvector1D, face normal                                           //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
array<double, 3>    u, v, w, t;
a3vector1D          N;

// Counters
int                 U, V0, V1;
int                 j, k, l, m;
int                 n;

// ========================================================================== //
// COMPUTE FACE NORMAL                                                        //
// ========================================================================== //
N.fill(0.0);
n = infos[e_type[T]].faces[i].size();
if (n <= 1) { return(N); }
else if (n == 2) {
    V0 = Simplex[T][infos[e_type[T]].faces[i][0]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][1]];
    u = V[V1] - V[V0];
    v[0] = 0.0;    v[1] = 0.0;    v[2] = 1.0;
    w = crossProduct(u, v);
    N = w/norm2(w);
}
else if (n == 3) {
    U = Simplex[T][infos[e_type[T]].faces[i][0]];
    V0 = Simplex[T][infos[e_type[T]].faces[i][1]];
    V1 = Simplex[T][infos[e_type[T]].faces[i][2]];
    u = V[V0] - V[U];
    v = V[V1] - V[U];
    t = crossProduct(u, v);
    N = t/norm2(t);
}
else {
    for (j = 0; j < n; ++j) {
        U = Simplex[T][infos[e_type[T]].faces[i][j]];
        t.fill(0.0);
        k = 0;
        l = (j+1) % n;
        while (k < n+2) {
            V0 = Simplex[T][infos[e_type[T]].faces[i][l]];
            V1 = Simplex[T][infos[e_type[T]].faces[i][(l+1)%n]];
            u = V[V0] - V[U];
            v = V[V1] - V[U];
            w = crossProduct(u, v);
            t = t + w/norm2(w);
            k++;
        } //next k
        t = t/((double) (n-2));
        for (k = 0; k < 3; ++k) {
            N[k] += t[k];
        } //next k
    } //next j
    N = N/norm2(N);
}

return(N); };


        // ----------------------------------------------------------------------------------- //
        double Class_VolTri::Volume(int T) {

        // =================================================================================== //
        // double Class_VolTri::Volume(int T)                                                  //
        //                                                                                     //
        // Compute volume of a given simplex.                                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T        : int, simplex global index                                              //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - V        : double, simplex volume                                                 //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        double             V;

        // Counters


        // =================================================================================== //
        // COMPUTE SIMPLEX VOLUME                                                              //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria; break;
            case 8  : goto label_quad; break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;
        };

        label_tria: {
//         if (dim == 2) {
//             V = 0.5*abs(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm2(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), dim);
//         }
        }
        label_quad: {
//         if (dim == 2) {
//             V = 0.5*abs(crossProduct(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                       Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]))
//               + 0.5*abs(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm2(crossProduct(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                          Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]), 3)
//               + 0.5*norm2(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), 3);
//         }
        }
        label_tetra: {
//         if (dim == 2) {
//             V = 0.0;
//         }
//         else if (dim == 3) {
//             V = 1.0/6.0 * abs(dotProduct(Vertex[Simplex[T][0]] - Vertex[Simplex[T][3]],
//                                           crossProduct(Vertex[Simplex[T][1]] - Vertex[Simplex[T][3]],
//                                                         Vertex[Simplex[T][2]] - Vertex[Simplex[T][3]]), dim));
//         }
        }
        label_pyram: {
            V = 0.0; // to be done
        }
        label_prism: {
            V = 0.0; // to be done
        }
        label_dhexa: {
            V = 0.0; // to be done
        }
        label_rhexa: {
            V = 0.0; // to be done
        }

        return(V); };

        // ----------------------------------------------------------------------------------- //
        double Class_VolTri::Volume(a3vector2D &X, int T) {

        // =================================================================================== //
        // double Class_VolTri::Volume(int T)                                                  //
        //                                                                                     //
        // Compute volume of a given simplex.                                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T        : int, simplex global index                                              //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - V        : double, simplex volume                                                 //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        double             V;

        // Counters


        // =================================================================================== //
        // COMPUTE SIMPLEX VOLUME                                                              //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria; break;
            case 8  : goto label_quad; break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;
        };

        label_tria: {
//         if (dim == 2) {
//             V = 0.5*abs(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm2(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), dim);
//         }
        }
        label_quad: {
//         if (dim == 2) {
//             V = 0.5*abs(crossProduct(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                       Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]))
//               + 0.5*abs(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm2(crossProduct(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                          Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]), 3)
//               + 0.5*norm2(crossProduct(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), 3);
//         }
        }
        label_tetra: {
//         if (dim == 2) {
//             V = 0.0;
//         }
//         else if (dim == 3) {
//             V = 1.0/6.0 * abs(dotProduct(Vertex[Simplex[T][0]] - Vertex[Simplex[T][3]],
//                                           crossProduct(Vertex[Simplex[T][1]] - Vertex[Simplex[T][3]],
//                                                         Vertex[Simplex[T][2]] - Vertex[Simplex[T][3]]), dim));
//         }
        }
        label_pyram: {
            V = 0.0; // to be done
        }
        label_prism: {
            V = 0.0; // to be done
        }
        label_dhexa: {
            V = 0.0; // to be done
        }
        label_rhexa: {
            V = 0.0; // to be done
        }

        return(V); };
