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
# include "SortAlgorithms.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
bool CGPLCurve::IsClosed(
    std::vector<std::vector<std::vector<int>>>        &Adjacency
) {

// ========================================================================== //
// bool CGPLCurve::IsClosed(                                                 //
//     std::vector<std::vector<std::vector<int>>>        &Adjacency)                                           //
//                                                                            //
// Check if a given curve is closed.                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Adjacency  : std::vector<std::vector<std::vector<int>>>, with simplex-simplex adjacency. Adjacency[i][j]  //
//                contains all curve segments adjacenct to the i-th curve     //
//                segment at vertex with local index j                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag       : bool, flag = true for closed curve, flag = false, otherwise //
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
flag = !CGPLCurve::IsOpen(Adjacency);

return(flag); };

// -------------------------------------------------------------------------- //
bool CGPLCurve::IsOpen(
    std::vector<std::vector<std::vector<int>>>        &Adjacency
) {

// ========================================================================== //
// bool CGPLCurve::IsOpen(                                                   //
//     std::vector<std::vector<std::vector<int>>>        &Adjacency)                                           //
//                                                                            //
// Check if a given curve is open.                                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Adjacency : std::vector<std::vector<std::vector<int>>>, with simplex-simplex adjacency. Adjacency[i][j]   //
//               contains all curve segments adjacenct to the i-th curve      //
//               segment at vertex with local index j                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag      : bool, flag = true for open curve, flag = false, otherwise.   //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool      flag = false;
int       n = Adjacency.size();

// Counters
int       i, j;

// ========================================================================== //
// CHECK IF CURVE IS CLOSED                                                   //
// ========================================================================== //
i = 0;
while(i < n && !flag) {
    j = 0;
    while (j < 2 && !flag) {
        flag = (Adjacency[i][j][0] < 0);
        j++;
    } //next j
    i++;
} //next i

return(flag); };

// -------------------------------------------------------------------------- //
bool CGPLCurve::IsManifold(
    std::vector<std::vector<std::vector<int>>>        &Adjacency
) {

// ========================================================================== //
// bool CGPLCurve::IsManifold(                                               //
//     std::vector<std::vector<std::vector<int>>>        &Adjacency)                                           //
//                                                                            //
// Check if a given curve is simple or not (no "T" junction).                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Adjacency : std::vector<std::vector<std::vector<int>>>, with simplex-simplex adjacency. Adjacency[i][j]   //
//               contains all curve segments adjacenct to the i-th curve      //
//               segment at vertex with local index j                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag      : bool, flag = true for simple curve, flag = false, otherwise  //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool      flag = true;
int       n = Adjacency.size();

// Counters
int       i, j;

// ========================================================================== //
// CHECK IF CURVE IS SIMPLE                                                   //
// ========================================================================== //
i = 0;
while(i < n && flag) {
    j = 0;
    while (j < 2 && flag) {
        flag = (Adjacency[i][j].size() == 1);
        j++;
    } //next j
    i++;
} //next i

return(flag); };

// -------------------------------------------------------------------------- //
bool CGPLCurve::IsClockWise(
    std::vector<std::array<double,3>>       &V,
    std::vector<std::vector<int>>       &S,
    std::vector<std::vector<std::vector<int>>>       &A
) {

// ========================================================================== //
// bool CGPLCurve::IsClockWise(                                              //
//     dvector2D       &V,                                                    //
//     std::vector<std::vector<int>>       &S,                                                    //
//     std::vector<std::vector<std::vector<int>>>       &A)                                                    //
//                                                                            //
// Check if a given curve is marched in clockwise direction.                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V         : dvector2D, vertex coordinate list. V[i][0], V[i][1], ...     //
//               are the x, y, ... coordinates of the i-th vertex.            //
// - S         : std::vector<std::vector<int>>, with segments-vertex connectivity. S[i][0],       //
//               S[i][1] are the global indices of vertices of the i-th       //
//               segment of the input PL curve.                               //
// - A         : std::vector<std::vector<std::vector<int>>>, with segment-segment adjacencies. A[i][j] stores  //
//               the list of global indices of all segments sharing the       //
//               vertex S[i][j] with the i-th segment.                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag      : bool, 'true' if curve is marched in clockwise direction,     //
//              'false' otherwise.                                            //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool        flag = true;
int         nS = S.size();
double      x1, y1, x2, y2, sArea;

// Counters
int         T, counter;

// ========================================================================== //
// CHECK IF CURVE IS MARCHED IN CLOCKWISE DIRECTION                           //
// ========================================================================== //
T = 0;
counter = 0;
sArea = 0.0;
while (counter < nS) {
    x1 = V[S[T][0]][0];
    y1 = V[S[T][0]][1];
    x2 = V[S[T][1]][0];
    y2 = V[S[T][1]][1];
    sArea += (x1 + x2) * (y1 - y2);
    T = A[T][1][0];
    counter++;
} //next T
flag = (sArea >= 0.0);

return(flag); };

// -------------------------------------------------------------------------- //
bool CGPLCurve::IsCounterClockWise(
    std::vector<std::array<double,3>>       &V,
    std::vector<std::vector<int>>       &S,
    std::vector<std::vector<std::vector<int>>>       &A
) {

// ========================================================================== //
// bool CGPLCurve::IsCounterClockWise(                                       //
//     dvector2D       &V,                                                    //
//     std::vector<std::vector<int>>       &S,                                                    //
//     std::vector<std::vector<std::vector<int>>>       &A)                                                    //
//                                                                            //
// Check if a given curve is marched in counterclockwise direction.           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V    : dvector2D, vertex coordinate list. V[i][0], V[i][1], ...          //
//          are the x, y, ... coordinates of the i-th vertex.                 //
// - S    : std::vector<std::vector<int>>, segments-vertex connectivity. S[i][0], S[i][1] are     //
//          the global indices of vertices of the i-th segment  of the input  //
//          PL curve.                                                         //
// - A    : std::vector<std::vector<std::vector<int>>>, with segment-segment adjacencies. A[i][j] stores the   //
//          list of global indices of all segments sharing the vertex S[i][j] //
//          with the i-th segment.                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - flag : bool, 'true' if curve is marched in counterclockwise direction,   //
//          'false' otherwise.                                                //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool        flag = true;

// Counters
// none

// ========================================================================== //
// CHECK IF CURVE IS MARHED IN CLOCKWISE DIRECTION                            //
// ========================================================================== //
flag = !CGPLCurve::IsClockWise(V, S, A);

return(flag); };

// -------------------------------------------------------------------------- //
void CGPLCurve::AdjustOrder(
    std::vector<std::vector<int>>        &S,
    std::vector<std::vector<std::vector<int>>>        &A,
    int               seed
) {

// ========================================================================== //
// void CGPLCurve::AdjustOrder(                                              //
//     std::vector<std::vector<int>>        &S,                                                   //
//     std::vector<std::vector<std::vector<int>>>        &A,                                                   //
//     int               seed);                                               //
//                                                                            //
// Adjust local direction onto curve to match direction of seed segment.      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - S     : std::vector<std::vector<int>>, simplex-vertex connectivity. S[i][j] is               //
//           the global index of the vertex with local index j for segment i  //
// - A     : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacency. A[i][j]                    //
//           contains all curve segments adjacenct to the i-th curve segment  //
//           at vertex S[i][j]                                                //
// - seed  : int, seed segment whose direction has to be matched.             //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                          i_dummy;
std::vector<bool>                    flag(S.size(), true);
std::vector<int>                    i_dummy1D;
bitpit::LIFOStack<std::vector<int>> next_dir(2);
std::vector<int>                    dir(3, -1);

// Counters
int                   i, j, m, n, I_, J;

// ========================================================================== //
// ADJUST LOCAL DIRECTIONS                                                    //
// ========================================================================== //

// Initialize propagation directions ---------------------------------------- //
flag[seed] = false;
n = A[seed].size();
for (i = 0; i < n; i++) {
    m = A[seed][i].size();
    for (j = 0; j < m; j++) {
        if (A[seed][i][j] >= 0) {
            dir[0] = seed;
            dir[1] = i;
            dir[2] = j;
            next_dir.push(dir);
        }
    } //next j
} //next i

// Loop until stack is empty ------------------------------------------------ //
while (next_dir.TOPSTK > 0) {

    // Pop next direction from stack
    dir = next_dir.pop();
    I_ = dir[0];
    i = dir[1];
    j = dir[2];
    J = A[I_][i][j];

    if (flag[J]) {

        // Update flag
        flag[J] = false;

        // Correct local direction
        if (S[I_][i] != S[J][1-i]) {

            // Update simplex-vertex connectivity
            i_dummy = S[J][1-i];
            S[J][1-i] = S[J][i];
            S[J][i] = i_dummy;

            // Update adjacencies
            i_dummy1D.resize(A[J][1-i].size());
            i_dummy1D = A[J][1-i];
            A[J][1-i].resize(A[J][i].size());
            A[J][1-i] = A[J][i];
            A[J][i].resize(i_dummy1D.size());
            A[J][i] = i_dummy1D;
        }
    }

    // Find next direction
    m = A[J][i].size();
    for (j = 0; j < m; j++) {
        if ((A[J][i][j] >= 0)
        && (flag[A[J][i][j]])) {
            dir[0] = J;
            dir[1] = i;
            dir[2] = j;
            next_dir.push(dir);
        }
    } //next j
    
} // next dir

return; };

// -------------------------------------------------------------------------- //
void CGPLCurve::InvertOrder(
    std::vector<std::vector<int>>        &S,
    std::vector<std::vector<std::vector<int>>>        &A
) {

// ========================================================================== //
// void CGPLCurve::InvertOrder(                                              //
//     std::vector<std::vector<int>>        &S,                                                   //
//     std::vector<std::vector<std::vector<int>>>        &A)                                                   //
//                                                                            //
// Invert local direction onto curve.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - S    : std::vector<std::vector<int>>, simplex-vertex connectivity. S[i][j] is the global     //
//          index of the vertex with local index j for segment i              //
// - A    : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacency. A[i][j] stores the global   //
//          indices of curve segments adjacenct to the i-th segment at vertex //
//          S[i][j]                                                           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int             i_dummy;
std::vector<int>       i_dummy1D;

// Counters
int             i, n;

// ========================================================================== //
// INVERT CURVE DIRECTION                                                     //
// ========================================================================== //
n = S.size();
for (i = 0; i < n; i++) {

    // Update simplex-vertex connectivity
    i_dummy = S[i][0];
    S[i][0] = S[i][1];
    S[i][1] = i_dummy;

    // Update simplex-simplex adjacencies
    i_dummy1D.resize(A[i][0].size());
    i_dummy1D = A[i][0];
    A[i][0].resize(A[i][1].size());
    A[i][0] = A[i][1];
    A[i][1].resize(i_dummy1D.size());
    A[i][1] = i_dummy1D;

} //next i

return; };

// -------------------------------------------------------------------------- //
void CGPLCurve::FindCorners(
    double            angle,
    std::vector<std::array<double,3>>        &V,
    std::vector<std::vector<int>>        &S,
    std::vector<std::vector<std::vector<int>>>        &A,
    std::vector<std::array<double,3>>        *t,
    std::vector<std::vector<int>>        &corners
) {

// ========================================================================== //
// void CGPLCurve::FindCorners(                                              //
//     double            angle,                                               //
//     dvector2D        &V,                                                   //
//     std::vector<std::vector<int>>        &S,                                                   //
//     std::vector<std::vector<std::vector<int>>>        &A,                                                   //
//     dvector2D        *t = NULL,                                            //
//     std::vector<std::vector<int>>        corners)                                              //
//                                                                            //
// Find sharp corners in a PL curve.                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - angle   : double, angle [rad] for corner to be marked as sharp corner    //
// - V       : std::vector<double>, vertex coordinate list. V[i][0], V[i][1], ...  are  //
//             the x, y, ... coordinates of the i-th vertex                   //
// - Simplex : std::vector<std::vector<int>>, with simplex-vertex coordinates. Simplex[i][0] and  //
//             S[i][1] are the global indices of vertices of the i-th segment //
// - Adj     : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacencies. A[i][j] store global   //
//             indices of segments adjacent to the i-th segment               //
//             at vertex S[i][j]                                              //
// - t       : *dvector2D, (optional) tangent unit vector to the curve at     //
//              each curve segment. t[i][0], t[i][1], ... are the x, y, ...   //
//              components of tangent unit vector to the i-th segment         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - corners : std::vector<std::vector<int>>, corner list. corners[i][0] and corners[i][1]        //
//             are the global index of segment containing corner and          //
//             the local index of vertex representing corner.                 //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double                          toll = cos(angle);

// Local variables
int                             nS = S.size();
double                          dp;
std::vector<bool>                       flag(nS, true);
std::vector<int>                       corner(2, -1);

// Counters
int                             i, n, I_, J;

// ========================================================================== //
// COMPUTE TANGENT VECTORS                                                    //
// ========================================================================== //
if (t == NULL) {
    //ComputeTangent(Vertex, Simplex, (*t));
}

// ========================================================================== //
// FIND SHARP CORNERS                                                         //
// ========================================================================== //
for (I_ = 0; I_ < nS; I_++) {
    flag[I_] = false;
    corner[0] = I_;
    n = S[I_].size();
    for (i = 0; i < n; i++) {
        corner[1] = i;
        if (A[I_][i].size() == 1) {
            if (A[I_][i][0] >= 0) {
                J = A[I_][i][0];
                if (flag[J]) {
                    dp = dotProduct((*t)[I_], (*t)[J]);
                    if (dp < toll) {
                        corners.push_back(corner);
                    }
                }
            }
        }
    } //next i
} //next I_

return; };

// -------------------------------------------------------------------------- //
void CGPLCurve::FindTriplePoints(
    std::vector<std::vector<int>>        &S,
    std::vector<std::vector<std::vector<int>>>        &A,
    std::vector<std::vector<int>>        &corners
) {

// ========================================================================== //
// void CGPLCurve::FindTriplePoints(                                         //
//     std::vector<std::vector<int>>        &S,                                                   //
//     std::vector<std::vector<std::vector<int>>>        &A,                                                   //
//     std::vector<std::vector<int>>        corners)                                              //
//                                                                            //
// Find triple points in a PL curve.                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V    : std::vector<double>, vertex coordinate list. V[i][0], V[i][1], ... are the  //
//          x, y, ... coordinates of the i-th vertex onto a PL curve.         //
// - A    : std::vector<std::vector<int>>, simplex-vertex coordinates. S[i][0] and  S[i][1] are   //
//          the vertices of the i-th segment of a PL curve.                   //
// - A    : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacencies. A[i][j] store global      //
//          indices of segments adjacent to i-th segment at vertex S[i][j]    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - corners : std::vector<std::vector<int>>, triple points list. corners[i][0] and corners[i][1] //
//             are the global index of segment containing corner and the      //
//             local index of vertex representing corner.                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                            check;
int                             nS = S.size();
std::vector<bool>                       flag(nS, true);
std::vector<int>                       corner(2, -1);

// Counters
int                             i, j, n, m, I_;

// ========================================================================== //
// FIND CORNERS                                                               //
// ========================================================================== //
for (I_ = 0; I_ < nS; I_++) {
    flag[I_] = false;
    corner[0] = I_;
    n = S[I_].size();
    for (i = 0; i < n; i++) {
        corner[1] = i;
        if (A[I_][i].size() > 1) {
            check = true;
            j = 0;
            m = A[I_][i].size();
            while ((j < m) && check) {
                check = check && flag[A[I_][i][j]];
                j++;
            } //next j
            if (check) {
                corners.push_back(corner);
            }
        }
    } //next i
} //next I_

return; };

// -------------------------------------------------------------------------- //
void CGPLCurve::FindEndPoints(
    std::vector<std::vector<int>>        &S,
    std::vector<std::vector<std::vector<int>>>        &A,
    std::vector<std::vector<int>>        &corners
) {

// ========================================================================== //
// void CGPLCurve::FindEndPoints(                                            //
//     std::vector<std::vector<int>>        &S,                                                   //
//     std::vector<std::vector<std::vector<int>>>        &A,                                                   //
//     std::vector<std::vector<int>>        corners)                                              //
//                                                                            //
// Find end points in a PL curve.                                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - S    : std::vector<std::vector<int>>, with simplex-vertex coordinates. S[i][0] and S[i][1]   //
//          are the vertices of the i-th segment of a PL curve.               //
// - A    : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacencies. A[i][j] store global      //
//          indices of segments adjacent to i-th segment at vertex S[i][j]    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - corners : std::vector<std::vector<int>>, endpoints list. corners[i][0] and corners[i][1]     //
//             are the global index of segment containing corner and the      //
//             local index of vertex representing corner.                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                             nS = S.size();
double                          dp;
std::vector<int>                       corner(2, -1);

// Counters
int                             i, n, I_;

// ========================================================================== //
// FIND CORNERS                                                               //
// ========================================================================== //
for (I_ = 0; I_ < nS; I_++) {
    corner[0] = I_;
    n = S[I_].size();
    for (i = 0; i < n; i++) {
        corner[1] = i;
        if (A[I_][i].size() == 1) {
            if (A[I_][i][0] < 0) {
                corners.push_back(corner);
            }
        }
    } //next i
} //next I_

return; };

// -------------------------------------------------------------------------- //
void CGPLCurve::BreakCurveAtPoint(
    int               T,
    int               i,
    std::vector<std::vector<int>>        &S,
    std::vector<std::vector<std::vector<int>>>        &A
) {

// ========================================================================== //
// void BreakCurveAtPoint(                                                    //
//     int               T,                                                   //
//     int               i,                                                   //
//     std::vector<std::vector<int>>        &S,                                                   //
//     std::vector<std::vector<std::vector<int>>>        &A)                                                   //
//                                                                            //
// Break a PL curve at a specified vertex.                                    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T        : int, global index of segments containing vertex               //
// - i        : int, local index of vertex                                    //
// - S        : std::vector<std::vector<int>>, simplex-vertex connectivity. S[i][j] is the global //
//              index of the vertex with local index j for segment i          //
// - A        : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacency. A[i][j] contains all    //
//              curve segments adjacenct to i-th curve segment at             //
//              vertex S[i][j]                                                //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// LOCAL VARIABLES                                                            //
// ========================================================================== //

// Local variables
// none

// Counters
int             j, n, R;

// ========================================================================== //
// OPEN CURVE AT SPECIFIED VERTEX                                             //
// ========================================================================== //

// Update adjacencies for neighoors of T ------------------------------------ //
n = A[T][i].size();
for (j = 0; j < n; j++) {
    R = A[T][i][j];
    if (R >= 0) {
        A[R][1-i].resize(1);
        A[R][1-i][0] = -1;
    }
} //next j

// Update adjacencies for T ------------------------------------------------- //
A[T][i].resize(1);
A[T][i][0] = -1;

return; };

// -------------------------------------------------------------------------- //
double CGPLCurve::Length(
    std::vector<std::array<double,3>>       &V,
    std::vector<std::vector<int>>       &S
) {

// ========================================================================== //
// double CGPLCurve::Length(                                                 //
//     dvector2D       &V,                                                    //
//     std::vector<std::vector<int>>       &S)                                                    //
//                                                                            //
// Compute curve length.                                                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V  : std::vector<double>, vertex coordinate list. V[i][0], V[i][1], ... are the    //
//        x, y, ... coordinates of the i-th vertex onto a PL curve.           //
// - S : std::vector<std::vector<int>>, simplex-vertex coordinates. S[i][0] and S[i][1] are the   //
//       vertices of the i-th segment of a PL curve.                          //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - L        : double, curve length                                          //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int             nS = S.size();
double          L = 0.0;

// Counters
int             T;

// ========================================================================== //
// COMPUTE CURVE LENGTH                                                       //
// ========================================================================== //
for (T = 0; T < nS; T++) {
    L += norm2(V[S[T][1]] - V[S[T][0]]);
} //next T

return(L); };

// -------------------------------------------------------------------------- //
void CGPLCurve::Curvature(
    std::vector<std::array<double,3>>       &V,
    std::vector<std::vector<int>>       &S,
    std::vector<std::vector<std::vector<int>>>       &A,
    std::vector<double>       &curv
) {

// ========================================================================== //
// void CGPLCurve::Curvature(                                                //
//     dvector2D       &V,                                                    //
//     std::vector<std::vector<int>>       &S,                                                    //
//     std::vector<std::vector<std::vector<int>>>       &A,                                                    //
//     std::vector<double>       &curv)                                                 //
//                                                                            //
// Compute local curvature of a simple curves.                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V  : std::vector<double>, with vertex coordinate list. V[i][0], V[i][1], ... are   //
//        the x, y, ... coordinates of the i-th vertex onto a PL curve.       //
// - S  : std::vector<std::vector<int>>, with simplex-vertex connectivity. S[i][0] and S[i][1]    //
//        are the vertices of the i-th segment of a PL curve.                 //
// - A  : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacencies. A[i][j] store global        //
//        indices of segments adjacent to the i-th segment at vertex S[i][j]  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - curv    : std::vector<double>, with local curvature computed at vertices.          //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int             nS = S.size();
double          ds;
std::array<double,3>         dt;
std::vector<std::array<double,3>>       t(nS);

// Counters
int             R, T;

// ========================================================================== //
// CHECK CURVE TOPOLOGY                                                       //
// ========================================================================== //

// Simple curve
if (!IsManifold(A)) {
    return;
}

// ========================================================================== //
// RESIZE INPUT/OUTPUT VARIABLES                                              //
// ========================================================================== //

// Curvature
curv.resize(V.size(), 0.0);

// ========================================================================== //
// COMPUTE CURVATURE                                                          //
// ========================================================================== //

// Compute tangent vector --------------------------------------------------- //
for (T = 0; T < nS; T++) {
    t[T] = V[S[T][1]] - V[S[T][0]];
    t[T] = t[T]/norm2(t[T]);
} //next T

// Compute curvature -------------------------------------------------------- //
for (T = 0; T < nS; T++) {
    if (A[T][0][0] >= 0) {

        // Adjacenct simplex
        R = A[T][0][0];

        // Normal direction
        dt = t[T] - t[R];
        ds = 0.5*norm2(V[S[R][1]] - V[S[R][0]])
           + 0.5*norm2(V[S[T][1]] - V[S[T][0]]);
        curv[S[T][0]] = norm2(dt)/ds;

    }
} //next T

return; };
