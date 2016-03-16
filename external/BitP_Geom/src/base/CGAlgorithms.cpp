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

// bitpit library
# include <bitpit_SA.hpp>


// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
double CGAlgorithms::Grad1DUpdate(
    int                        nSimplex,
    std::vector<std::array<double,3>>                 &Vertex,
    std::vector<std::vector<int>>                 &Simplex,
    std::vector<std::vector<std::vector<int>>>                 &Adjacency,
    std::vector<double>                 &val,
    int                        T,
    int                        i,
    double                     g,
    std::vector<bool>                 &flag
) {

// ========================================================================== //
// double CGAlgorithms::Grad1DUpdate(                                        //
//     int                        nSimplex,                                   //
//     dvector2D                 &Vertex,                                     //
//     std::vector<std::vector<int>>                 &Simplex,                                    //
//     std::vector<std::vector<std::vector<int>>>                 &Adjacency,                                  //
//     std::vector<double>                 &val,                                        //
//     int                        T,                                          //
//     int                        i,                                          //
//     double                     g,                                          //
//     std::vector<bool>                 &flag)                                       //
//                                                                            //
// Compute local solution to the 1D gradient limiting equation on a           //
// 1D manifold in a 2D Euclidean space.                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nSimplex  : int, number of simplicies                                    //
// - Vertex    : std::vector<double>, vertex coordinate list. Vertex[i][0],             //
//               Vertex[i][1] are the x, y, coordinates of the i-th vertex    //
// - Simplex   : std::vector<std::vector<int>>, simplex-vertex connectivity. Simplex[i][0] and    //
//               Simplex[i][1] are the global indices of vertices of the i-th //
//               segment.                                                     //
// - Adjacency : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacency. Adjacency[i][j] stores //
//               the global indices of all simplicies adjacenct to the i-th   //
//               simplex at vertex Simplex[i][j].                             //
// - val       : std::vector<double>, scalar field at each mesh vertex                  //
// - T         : int, simplex global index                                    //
// - i         : int, vertex local index                                      //
// - g         : double, max slope                                            //
// - flag      : std::vector<bool> with flag for dead (flag = false), and alive       //
//               (flag = true) vertexes.                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - value     : solution to the 1D gradient limiting equation                //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double     value = 1.0e+18, value_0 = 1.0e+18;
std::array<double,3>    P0;
P0.fill(0.0) ;

// Counters
int        A, V, W, j, k, m;

// ========================================================================== //
// COMPUTE THE LOCAL SOLUTION TO THE 1D GRADIENT LIMITING EQUATION            //
// ========================================================================== //

// Vertex global index ------------------------------------------------------ //
V = Simplex[T][i];

// Find dead neighboors ----------------------------------------------------- //
W = Simplex[T][(i+1) % Simplex[T].size()];
if (!flag[W]) {
    value_0 = val[W];
    P0 = Vertex[W];
    value = std::min(value, value_0 + g * norm2(Vertex[V] - P0));
}
if (Adjacency[T][i][0] >= 0) {
    m = Adjacency[T][i].size();
    for (k = 0; k < m; k++) {
        A = Adjacency[T][i][k];
        j = 0;
        if (Simplex[A][0] == V) { j = 1; }
        W = Simplex[A][j];

        if (!flag[W]) {
            value_0 = val[W];
            P0 = Vertex[W];
            value = std::min(value, value_0 + g * norm2(Vertex[V] - P0));
        }
    } //next k
}

// Solve 1D grad limiting equation ------------------------------------------ //
value = std::min(val[V], value);

return(value); };

// -------------------------------------------------------------------------- //
void CGAlgorithms::GradLimiting1D(
    int                        nSimplex,
    std::vector<std::array<double,3>>                 &Vertex,
    std::vector<std::vector<int>>                 &Simplex,
    std::vector<std::vector<std::vector<int>>>                 &Adjacency,
    std::vector<double>                 &val,
    double                     g
) {

// ========================================================================== //
// void CGAlgorithms::GradLimiting1D(                                        //
//     int                        nSimplex,                                   //
//     dvector2D                 &Vertex,                                     //
//     std::vector<std::vector<int>>                 &Simplex,                                    //
//     std::vector<std::vector<std::vector<int>>>                 &Adjacency,                                  //
//     std::vector<double>                 &val,                                        //
//     double                     g)                                          //
//                                                                            //
// Solve the 1D gradient limiting equation onto a 1D manifold in a 2D ,       //
// Euclidean space, using fast marching method.                               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nSimplex  : int, number of simplicies                                    //
// - Vertex    : std::vector<double>, vertex coordinate list. Vertex[i][0],             //
//               Vertex[i][1] are the x, y, coordinates of the i-th vertex    //
// - Simplex   : std::vector<std::vector<int>>, simplex-vertex connectivity. Simplex[i][0] and    //
//               Simplex[i][1] are the global indices of vertices of the i-th //
//               segment.                                                     //
// - Adjacency : std::vector<std::vector<std::vector<int>>>, simplex-simplex adjacency. Adjacency[i][j] stores //
//               the global indices of all simplicies adjacenct to the i-th   //
//               simplex at vertex Simplex[i][j].                             //
// - val       : std::vector<double>, with scalar field to be limited.                  //
// - g         : double, max slope value.                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                                   nVertex = Vertex.size();
double                                ddummy;
std::vector<bool>                     inserted(nVertex, false);
std::vector<int>                             idummy1D(2, -1);
// array<int,2>                          idummy1D;
// std::vector<std::vector<int>>         map(nVertex, std::vector<int>(2, -1)), *map_ = &map;
std::vector< std::array<int,2> >      map(nVertex), *map_ = &map;
bitpit::MinPQueue<double, std::vector<int>>  heap(nVertex, true, map_);

// Counters
int                                 i, j, k;
int                                 T, V, A, W;
int                                 m;

for (V = 0; V < nVertex; V++) {
    map[V].fill(0.) ;
};

// ========================================================================== //
// BUILD MIN HEAP                                                             //
// ========================================================================== //

// Initialize data structure for min heap ----------------------------------- //
k = 0;
for (T = 0; T < nSimplex; T++) {
    m = Simplex[T].size();
    for (i = 0; i < m; i++) {
        V = Simplex[T][i];
        if (!inserted[V]) {

            // Update mapper
            map[k][0] = V;
            map[V][1] = k;

            // Insert vertex into min-heap
            idummy1D[0] = T;
            idummy1D[1] = i;
            heap.insert(val[V], idummy1D);

            // Update flag & counters
            k++;
            inserted[V] = true;
        }
    } //next i
} //next T

// ========================================================================== //
// FAST GRADIENT LIMITING                                                     //
// ========================================================================== //
while (heap.heap_size > 0) {

    // Extract root from min heap
    heap.extract(ddummy, idummy1D);

    // Update value for extracted item
    T = idummy1D[0];
    i = idummy1D[1];
    V = Simplex[T][i];
    inserted[V] = false;

    // Update value for neighbors
    j = (i + 1) % Simplex[T].size();
    W = Simplex[T][j];
    if (inserted[W]) {
        val[W] = Grad1DUpdate(nSimplex, Vertex, Simplex, Adjacency, val, T, j, g, inserted);
        idummy1D[0] = T;
        idummy1D[1] = j;
    }
    if (Adjacency[T][i][0] >= 0) {
        m = Adjacency[T][i].size();
        for (k = 0; k < m; k++) {
            A = Adjacency[T][i][k];
            j = 0;
            if (Simplex[A][0] == V) { j = 1; }
            W = Simplex[A][j];
            if (inserted[W]) {
                val[W] = Grad1DUpdate(nSimplex, Vertex, Simplex, Adjacency, val, A, j, g, inserted);
                idummy1D[0] = A;
                idummy1D[1] = j;
                heap.modify(map[W][1], val[W], idummy1D);
            }
        } //next k
    }
} //next item

return; };

// -------------------------------------------------------------------------- //
double CGAlgorithms::Grad2DUpdate(
    int                        nSimplex,
    std::vector<std::array<double,3>>                 &Vertex,
    std::vector<std::vector<int>>                 &Simplex,
    std::vector<std::vector<std::vector<int>>>                 &ring_1,
    std::vector<double>                 &val,
    int                        T,
    int                        i,
    double                     g,
    std::vector<bool>                 &flag
) {

// ========================================================================== //
// double CGAlgorithms::Grad2DUpdate(                                        //
//     int                        nSimplex,                                   //
//     dvector2D                 &Vertex,                                     //
//     std::vector<std::vector<int>>                 &Simplex,                                    //
//     std::vector<std::vector<std::vector<int>>>                 &ring_1,                                     //
//     std::vector<double>                 &val,                                        //
//     int                        T,                                          //
//     int                        i,                                          //
//     double                     g,                                          //
//     std::vector<bool>                 &flag)                                       //
//                                                                            //
// Compute local solution to the 2D gradient limiting equation on a           //
// 2D manifold in a 3D Euclidean space.                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nSimplex  : int, number of simplicies in the surface tasselation         //
// - Vertex    : std::vector<double>, vertex coordinate list. Vertex[i][0],             //
//               Vertex[i][1], ... are the x, y, ... coordinates of the i-th  //
//               vertex                                                       //
// - Simplex   : std::vector<std::vector<int>>, simplex-vertex connectivity. Simplex[i][0] and    //
//               Simplex[i][1] are the global indices of vertices of the i-th //
//               segment.                                                     //
// - ring_1    : std::vector<std::vector<std::vector<int>>>, 1-ring of simplcies for each vertex.              //
// - val       : std::vector<double>, scalar field at each mesh vertex                  //
// - T         : int, simplex global index                                    //
// - i         : int, vertex local index                                      //
// - g         : double, max slope                                            //
// - flag      : std::vector<bool> with flag for dead (flag = false), and alive       //
//               (flag = true) vertexes.                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - value     : solution to the 1D gradient limiting equation                //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double     value = 1.0e+18;

// Counters
int        A, V, W, j, k, m;

// ========================================================================== //
// COMPUTE THE LOCAL SOLUTION TO THE 1D GRADIENT LIMITING EQUATION            //
// ========================================================================== //

// Vertex global index ------------------------------------------------------ //
V = Simplex[T][i];

// Find dead neighboors ----------------------------------------------------- //
for (j = 0; j < ring_1[V].size(); ++j) {
    A = ring_1[V][j][0];
    k = ring_1[V][j][1];
    m = (k - 1 + Simplex[A].size()) % Simplex[A].size();
    W = Simplex[A][m];
    if (!flag[W]) {
        value = std::min(value, val[W] + g * norm2(Vertex[V] - Vertex[W]));
    }
    m = (k + 1) % Simplex[A].size();
    W = Simplex[A][m];
    if (!flag[W]) {
        value = std::min(value, val[W] + g * norm2(Vertex[V] - Vertex[W]));
    }
} //next j

// Solve 1D grad limiting equation ------------------------------------------ //
value = std::min(val[V], value);

return(value); };

// -------------------------------------------------------------------------- //
void CGAlgorithms::GradLimiting2D(
    int                        nSimplex,
    std::vector<std::array<double,3>>                 &Vertex,
    std::vector<std::vector<int>>                 &Simplex,
    std::vector<double>                 &val,
    double                     g
) {

// ========================================================================== //
// void CGAlgorithms::GradLimiting2D(                                        //
//     int                        nSimplex,                                   //
//     dvector2D                 &Vertex,                                     //
//     std::vector<std::vector<int>>                 &Simplex,                                    //
//     std::vector<double>                 &val,                                        //
//     double                     g)                                          //
//                                                                            //
// Solve the 2D grad-limiting equation, using a fast marching method over     //
// a 2D manifold immersed in a 3D Euclidean space.                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nSimplex       : int, number of simplicies in surface tasselation        //
// - Vertex         : dvector2D, vertex coordinate list. Vertex[i][0],        //
//                    Vertex[i][1], ... are the x, y, ... coordinates of the  //
//                    i-th vertex                                             //
// - Simplex        : std::vector<std::vector<int>>, simplex-vertex connectivity                  //
// - val            : std::vector<double>, field to be smoothed                         //
// - g              : double, max gradient.                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// // Local variables
int                                 nVertex = Vertex.size();
double                              ddummy;
std::vector<bool>                   inserted(nVertex, false);
std::vector<int>                           idummy1D(2, -1);
std::vector<std::vector<std::vector<int>>>  Ring1(nVertex);
std::vector< std::array<int,2> >            map(nVertex), *map_ = &map;
bitpit::MinPQueue<double, std::vector<int>>        heap(nVertex, true, map_);

// // Counters
int                                 i, j, k;
int                                 T, V, A, W;
int                                 m;

for (V = 0; V < nVertex; V++) {
        map[V].fill(0.) ;
};

// ========================================================================== //
// BUILD 1-RING OF VERTICES                                                   //
// ========================================================================== //
for (T = 0; T < nSimplex; T++) {
    m = Simplex[T].size();
    for (i = 0; i < m; ++i) {
        V = Simplex[T][i];
        idummy1D[0] = T;
        idummy1D[1] = i;
        Ring1[V].push_back(idummy1D);
    } //next i
} //next T

// ========================================================================== //
// BUILD MIN HEAP                                                             //
// ========================================================================== //

// Initialize data structure for min heap ----------------------------------- //
k = 0;
for (T = 0; T < nSimplex; T++) {
    m = Simplex[T].size();
    for (i = 0; i < m; i++) {
        V = Simplex[T][i];
        if (!inserted[V]) {

            // Update mapper
            map[k][0] = V;
            map[V][1] = k;

            // Insert vertex into min-heap
            idummy1D[0] = T;
            idummy1D[1] = i;
            heap.insert(val[V], idummy1D);

            // Update flag & counters
            k++;
            inserted[V] = true;
        }
    } //next i
} //next T

// ========================================================================== //
// FAST GRADIENT LIMITING                                                     //
// ========================================================================== //
while (heap.heap_size > 0) {

    // Extract root from min heap
    heap.extract(ddummy, idummy1D);

    // Update value for extracted item
    T = idummy1D[0];
    i = idummy1D[1];
    V = Simplex[T][i];
    inserted[V] = false;

    // Update value for neighboring vertices
    for (j = 0; j < Ring1[V].size(); ++j) {
        A = Ring1[V][j][0];
        k = Ring1[V][j][1];
        m = (k - 1 + Simplex[A].size()) % Simplex[A].size();
        W = Simplex[A][m];
        if (inserted[W]) {
            val[W] = Grad2DUpdate(nSimplex, Vertex, Simplex, Ring1, val, A, m, g, inserted);
            idummy1D[0] = A;
            idummy1D[1] = m;
            heap.modify(map[W][1], val[W], idummy1D);
        };
        m = (k + 1) % Simplex[A].size();
        W = Simplex[A][m];
        if (inserted[W]) {
            val[W] = Grad2DUpdate(nSimplex, Vertex, Simplex, Ring1, val, A, m, g, inserted);
            idummy1D[0] = A;
            idummy1D[1] = m;
            heap.modify(map[W][1], val[W], idummy1D);
        }
    } //next j
} //next item

return; };

// -------------------------------------------------------------------------- //
double CGAlgorithms::Grad2DUpdate(
    int                        nSimplex,
    std::vector<std::array<double,3>>                 &Vertex,
    std::vector<std::vector<int>>                 &Simplex,
    std::vector<std::vector<int>>                 &Adjacency,
    std::vector<double>                 &val,
    int                        T,
    double                     g,
    std::vector<bool>                 &flag
) {

// ========================================================================== //
// double CGAlgorithms::Grad2DUpdate(                                        //
//     int                        nSimplex,                                   //
//     dvector2D                 &Vertex,                                     //
//     std::vector<std::vector<int>>                 &Simplex,                                    //
//     std::vector<std::vector<int>>                 &Adjacency,                                  //
//     std::vector<double>                 &val,                                        //
//     int                        T,                                          //
//     double                     g,                                          //
//     std::vector<bool>                 &flag)                                       //
//                                                                            //
// Compute local solution to the 2D gradient limiting equation on a           //
// 2D volume.                                                                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nSimplex  : int, number of simplicies in the surface tasselation         //
// - Vertex    : std::vector<double>, vertex coordinate list. Vertex[i][0],             //
//               Vertex[i][1], ... are the x, y, ... coordinates of the cell  //
//               center of the i-th simplex                                   //
// - Simplex   : std::vector<std::vector<int>>, simplex-vertex connectivity. Simplex[i][0] and    //
//               Simplex[i][1] are the global indices of vertices of the i-th //
//               segment.                                                     //
// - Adjacency : std::vector<std::vector<int>>, simplex-simplex adjacencies                       //
// - val       : std::vector<double>, scalar field at each mesh vertex                  //
// - T         : int, simplex global index                                    //
// - g         : double, max slope                                            //
// - flag      : std::vector<bool> with flag for dead (flag = false), and alive       //
//               (flag = true) vertexes.                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - value     : solution to the 2D gradient limiting equation                //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double     value = 1.0e+18;

// Counters
int        A, j, m;

// ========================================================================== //
// COMPUTE THE LOCAL SOLUTION TO THE 2D GRADIENT LIMITING EQUATION            //
// ========================================================================== //

// Find dead neighboors ----------------------------------------------------- //
m = Adjacency[T].size();
for (j = 0; j < m; ++j) {
    A = Adjacency[T][j];
    if (!flag[A]) {
        value = std::min(value, val[A] + g * norm2(Vertex[T] - Vertex[A]));
    }
} //next j

// Solve 1D grad limiting equation ------------------------------------------ //
value = std::min(val[T], value);

return(value); };

// -------------------------------------------------------------------------- //
void CGAlgorithms::GradLimiting2D(
    int                        nSimplex,
    std::vector<std::array<double,3>>                 &Vertex,
    std::vector<std::vector<int>>                 &Simplex,
    std::vector<std::vector<int>>                 &Adjacency,
    std::vector<double>                 &val,
    double                     g
) {

// ========================================================================== //
// void CGAlgorithms::GradLimiting2D(                                        //
//     int                        nSimplex,                                   //
//     dvector2D                 &Vertex,                                     //
//     std::vector<std::vector<int>>                 &Simplex,                                    //
//     std::vector<std::vector<int>>                 &Adjacency,                                  //
//     std::vector<double>                 &val,                                        //
//     double                     g)                                          //
//                                                                            //
// Solve the 2D grad-limiting equation, using a fast marching method over     //
// a 2D domain.                                                               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nSimplex       : int, number of simplicies in volume tasselation         //
// - Vertex         : dvector2D, vertex coordinate list. Vertex[i][0],        //
//                    Vertex[i][1], ... are the x, y, ... coordinates of the  //
//                    cell center of the i-th simplex                         //
// - Simplex        : std::vector<std::vector<int>>, simplex-vertex connectivity                  //
// - Adjacency      : std::vector<std::vector<int>>, simplex-simplex adjacency                    //
// - val            : std::vector<double>, field to be smoothed                         //
// - g              : double, max gradient.                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                                 idummy;
double                              ddummy;
std::vector<bool>                   inserted(nSimplex, false);
std::vector< std::array<int,2> >    map(nSimplex), *map_ = &map;
bitpit::MinPQueue<double, int>      heap(nSimplex, true, map_);

// Counters
int                                 i, j, k;
int                                 T, A;
int                                 m;

for (T = 0; T < nSimplex; T++) {
        map[T].fill(0.) ;
};

// ========================================================================== //
// BUILD MIN HEAP                                                             //
// ========================================================================== //

// Initialize data structure for min heap ----------------------------------- //
k = 0;
for (T = 0; T < nSimplex; T++) {
    if (!inserted[T]) {

        // Update mapper
        map[k][0] = T;
        map[T][1] = k;

        // Insert vertex into min-heap
        idummy = T;
        heap.insert(val[T], idummy);

        // Update flag & counters
        k++;
        inserted[T] = true;
    }
} //next T

// ========================================================================== //
// FAST GRADIENT LIMITING                                                     //
// ========================================================================== //
while (heap.heap_size > 0) {

    // Extract root from min heap
    heap.extract(ddummy, idummy);

    // Update value for extracted item
    T = idummy;
    inserted[T] = false;

    // Update value for neighboring vertices
    m = Adjacency[T].size();
    for (j = 0; j < m; ++j) {
        A = Adjacency[T][j];
        if (inserted[A]) {
            val[A] = Grad2DUpdate(nSimplex, Vertex, Simplex, Adjacency, val, A, g, inserted);
            idummy = A;
            heap.modify(map[A][1], val[A], idummy);
        };
    } //next j
} //next item

return; };



