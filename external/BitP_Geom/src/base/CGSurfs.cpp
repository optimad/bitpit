// ========================================================================== //
//                 - COMPUTATIONAL GEOMETRY PACKAGE -                         //
//                                                                            //
// Routines for computational geometry.                                       //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v1.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "CGBase.hpp"
# include "Operators.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
bool CGPLSurf::IsClosed(
    std::vector<std::vector<std::vector<int>>>        &Adjacency
) {

// ========================================================================== //
// bool CGPLSurf::IsClosed(                                                  //
//     std::vector<std::vector<std::vector<int>>>        &Adjacency)                                           //
//                                                                            //
// Check if a given surf is closed.                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Adjacency  : std::vector<std::vector<std::vector<int>>>, with simplex-simplex adjacency. Adjacency[i][j]  //
//                contains all simplicies adjacenct to the i-th simplex       //
//                at edge with local index j                                  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag       : bool, flag = true for closed surf, flag = false, otherwise  //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool      flag = true;

// Counters
// none

// ========================================================================== //
// CHECK IF CURVE IS CLOSED                                                   //
// ========================================================================== //
flag = !CGPLSurf::IsOpen(Adjacency);

return(flag); };

// -------------------------------------------------------------------------- //
bool CGPLSurf::IsOpen(
    std::vector<std::vector<std::vector<int>>>        &Adjacency
) {

// ========================================================================== //
// bool CGPLSurf::IsOpen(                                                    //
//     std::vector<std::vector<std::vector<int>>>        &Adjacency)                                           //
//                                                                            //
// Check if a given surf is open.                                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Adjacency : std::vector<std::vector<std::vector<int>>>, with simplex-simplex adjacency. Adjacency[i][j]   //
//               contains all simplicies adjacenct to the i-th simplex at     //
//               edge with local index j                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag      : bool, flag = true for open surf, flag = false, otherwise.    //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool      flag = false;
int       n = Adjacency.size();

// Counters
int       i, j, m;

// ========================================================================== //
// CHECK IF CURVE IS CLOSED                                                   //
// ========================================================================== //
i = 0;
while(i < n && !flag) {
    j = 0;
    m = Adjacency[i].size();
    while (j < m && !flag) {
        flag = (Adjacency[i][j][0] < 0);
        j++;
    } //next j
    i++;
} //next i

return(flag); };

// -------------------------------------------------------------------------- //
bool CGPLSurf::IsManifold(
    std::vector<std::vector<std::vector<int>>>        &Adjacency
) {

// ========================================================================== //
// bool CGPLSurf::IsManifold(                                                //
//     std::vector<std::vector<std::vector<int>>>        &Adjacency)                                           //
//                                                                            //
// Check if a given surf is simple or not (no "T" junction).                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Adjacency : std::vector<std::vector<std::vector<int>>>, with simplex-simplex adjacency. Adjacency[i][j]   //
//               contains all simplicies adjacenct to the i-th simplex at     //
//               edge with local index j                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag      : bool, flag = true for simple surf, flag = false, otherwise   //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool      flag = true;
int       n = Adjacency.size();

// Counters
int       i, j, m;

// ========================================================================== //
// CHECK IF CURVE IS SIMPLE                                                   //
// ========================================================================== //
i = 0;
while(i < n && flag) {
    j = 0;
    m = Adjacency[i].size();
    while (j < m && flag) {
        flag = (Adjacency[i][j].size() == 1);
        j++;
    } //next j
    i++;
} //next i

return(flag); };

// -------------------------------------------------------------------------- //
void CGPLSurf::maxCurvature(
    std::vector<std::array<double,3>>       &V,
    std::vector<std::array<double,3>>       &N,
    std::vector<std::vector<int>>       &S,
    std::vector<double>       &curv
) {

// ========================================================================== //
// void CGPLSurf::maxCurvature(                                              //
//     dvector2D       &V,                                                    //
//     dvector2D       &N,                                                    //
//     std::vector<std::vector<int>>       &S,                                                    //
//     std::vector<double>       &curv)                                                 //
//                                                                            //
// Compute maximal discrete curvature of a simple surf.                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V  : std::vector<double>, with vertex coordinate list. V[i][0], V[i][1], V[i][2]   //
//        are the x, y, z coordinates of the i-th vertex onto a PL surf.      //
// - N  : dvector2D, normal unit vector at surface vertices. N[i][0], N[i][1] //
//        and N[i][2] are the x, y, z components of the normal unit vector    //
//        at the i-th mesh vertex.                                            //
// - S  : std::vector<std::vector<int>>, with simplex-vertex connectivity. S[i][0], S[i][1], ...  //
//        are the vertices of the i-th simplex of a PL surf.                  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - curv    : std::vector<double>, with local curvature computed at vertices.          //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                 nS = S.size();
double              d;
std::array<double,3>             v;
v.fill(0.0) ;

// Counters
int                 i, j;
int                 m;
int                 T, U, W;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //

// Curvature
curv.resize(V.size(), 0.0);

// ========================================================================== //
// COMPUTE CURVATURE                                                          //
// ========================================================================== //
for (T = 0; T < nS; ++T) {
    m = S[T].size();
    for (i = 0; i < m; ++i) {
        j = (i+1) % m;
        U = S[T][i];
        W = S[T][j];
        v = V[W] - V[U] ;
        d = std::pow(norm2(v), 2);
        curv[U] = std::max(curv[U], 2.0*abs(dotProduct(N[U], v))/d);
    } //next i
} //next R

return; };
