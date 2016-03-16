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
# include <bitpit_LA.hpp>

# include "Delaunay.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
bool Delaunay::Delaunay2D::SWAP_Criterion(
    a3vector1D          &x1,
    a3vector1D          &x2,
    a3vector1D          &x3,
    a3vector1D          &xp
) {

// ========================================================================== //
// bool Delaunay::Delaunay2D::SWAP_Criterion(                                 //
//     a3vector1D          &x1,                                               //
//     a3vector1D          &x2,                                               //
//     a3vector1D          &x3,                                               //
//     a3vector1D          &xp)                                               //
//                                                                            //
// Cline - Renka criterion for edge swapping                                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - x1    : a3vector1D 1st vertex of the triangle opposed to P               //
// - x2    : a3vector1D 2nd vertex of the triangle opposed to P               //
// - x3    : a3vector1D 3rd vertex of the triangle opposed to P               //
// - xp    : a3vector1D last point added to the triangulation                 //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, value is .TRUE., if edges must be swapped,                 //
//          .FALSE., otherwise.                                               //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool    check;

// Counters
// none

// ========================================================================== //
// SWAPPING CRITERION                                                         //
// ========================================================================== //

// Cline-Rinka criterion ---------------------------------------------------- //
{

    // Scope variables
    double  cosa, cosb;
    double  sina, sinb;
    double  x13, x23, x1P, x2P;
    double  y13, y23, y1P, y2P;

    // Parameters
    x13 = x1[0] - x3[0];
    x23 = x2[0] - x3[0];
    x1P = x1[0] - xp[0];
    x2P = x2[0] - xp[0];
    y13 = x1[1] - x3[1];
    y23 = x2[1] - x3[1];
    y1P = x1[1] - xp[1];
    y2P = x2[1] - xp[1];
    cosa = x13*x23 + y13*y23;
    cosb = x2P*x1P + y1P*y2P;

    // Check
    if ((cosa >= 0.0) && (cosb >= 0.0)) {
        check = false;
    }
    else if ((cosa < 0.0) && (cosb < 0.0)) {
        check = true;
    }
    else {
        sina = x13*y23 - x23*y13;
        sinb = x2P*y1P - x1P*y2P;
        if ((sina*cosb + sinb*cosa) < 0.0) {
            check = true;
        }
        else {
            check = false;
        }
    }

}

return (check); };

// -------------------------------------------------------------------------- //
bool Delaunay::Delaunay2D::ConvexPolygon(
    VolTriPatch        &Tri,
    int           const &L,
    int           const &R
) {

// ========================================================================== //
// bool Delaunay::Delaunay2D::ConvexPolygon(                               //
//     VolTriPatch        &Tri,                                              //
//     int           const &L,                                                //
//     int           const &R)                                                //
//                                                                            //
// Check wheter neighboring triangle L and R form a convex polygon.           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri   : VolTriPatch, volume triangulation                               //
// - L     : int, 1st triangle global index                                   //
// - R     : int, 2nd triangle global index                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, 'true' if L and R form a convex polygon, 'false' otherwise //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double const        abs_toll = 1.0e-12;

// Local variables
bool                check = true;
array<double, 2>    u, v;

// Counters
int                 i, j, k;
ivector1D           V(4, -1);

// ========================================================================== //
// CHECK IF A AND B FORM A CONVEX POLYGON                                     //
// ========================================================================== //

// Polygon vertices --------------------------------------------------------- //
i = Tri.face(L, R);
V[1] = Tri.Simplex[L][i];
i = (i + 1) % 3;
V[3] = Tri.Simplex[L][i];
i = (i + 1) % 3;
V[0] = Tri.Simplex[L][i];
i = Tri.face(R, L);
i = (i + 2) % 3;
V[2] = Tri.Simplex[R][i];

// Convexity check ---------------------------------------------------------- //
i = 0;
while (check && (i < 4)) {
    j = (i + 1) % 4;
    k = (j + 1) % 4;
    u[0] = Tri.Vertex[V[k]][0] - Tri.Vertex[V[j]][0];
    u[1] = Tri.Vertex[V[k]][1] - Tri.Vertex[V[j]][1];
    u = u/max(1.e-16, norm2(u));
    v[0] = Tri.Vertex[V[i]][0] - Tri.Vertex[V[j]][0];
    v[1] = Tri.Vertex[V[i]][1] - Tri.Vertex[V[j]][1];
    v = v/max(1.e-16, norm2(v));
    check = (crossProduct(u, v) >= abs_toll);
    i++;
} //next i

return(check); }

// -------------------------------------------------------------------------- //
void Delaunay::Delaunay2D::InsertVertices(
    VolTriPatch        &Tri,
    ivector1D           &Vlist
) {

// ========================================================================== //
// void Delaunay::Delaunay2D::InsertVertices(                              //
//     VolTriPatch        &Tri,                                              //
//     ivector1D           &Vlist)                                            //
//                                                                            //
// Insert vertices into Delaunay triangulation. Vertices are inserted in      //
// the same order specified in Vlist. No duplicated vertices are allowed      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri     : VolTriPatch, unstructured volume mesh manager                 //
// - Vlist   : ivector1D, with global index of vertices to be inserted into   //
//             the Delaunay triangulation                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
int    const                nV = Vlist.size();
int    const                MAXSTK = (int) sqrt((double) nV);

// // Local variables
bool                          check;
int                           A, B, C, T, L, R, ERL, ERA, ERB;
int                           V1, V2, V3;
ivector1D                     STACK_data(2, -1);
a3vector1D                    xp;
bitpit::LIFOStack<ivector1D>  STACK(MAXSTK);

// Counters
int                         i, j, k, P;

// ========================================================================== //
// ADD VERTICES TO THE DELAUNAY TRIANGULATION                                 //
// ========================================================================== //
xp.fill(0.0);
T = 0;
for (i = 0; i < nV; i++) {

    // Extract the i-th vertex from vertex list ----------------------------- //
    P = Vlist[i];
    xp[0] = Tri.Vertex[P][0];
    xp[1] = Tri.Vertex[P][1];

    // Locate triangle enclosing point P ------------------------------------ //
    Delaunay::Delaunay2D::ReturnSimplexID(Tri, xp, T);

    // Add new triangles to the triangulation ------------------------------- //

    // Extract adiacency data
    A = Tri.Adjacency[T][0];
    B = Tri.Adjacency[T][1];
    C = Tri.Adjacency[T][2];

    // Add new triangles
    Tri.SplitTri3Tri(T, P);

    // Put edges in a stack ------------------------------------------------- //
    if (A >= 0) {
        STACK_data[0] = T;
        STACK_data[1] = 1;
        STACK.push(STACK_data);
    }
    if (B >= 0) {
        STACK_data[0] = Tri.nSimplex-2;
        STACK_data[1] = 1;
        STACK.push(STACK_data);
    }
    if (C >= 0) {
        STACK_data[0] = Tri.nSimplex-1;
        STACK_data[1] = 1;
        STACK.push(STACK_data);
    }

    // Swap edges ----------------------------------------------------------- //
    while (STACK.TOPSTK > 0) {

        // Pop last element in the stack list
        STACK_data = STACK.pop();
        L = STACK_data[0];
        j = STACK_data[1];

        // Adjacenct triangle
        R = Tri.Adjacency[L][j];

        if (R >= 0) {

            // Vertex opposite to the face candidate for swapping
            k = (j+2) % 3;
            P = Tri.Simplex[L][k];
            xp = Tri.Vertex[P];

            // Vertex in the swapping region
            ERL = Tri.face(R, L);
            ERA = (ERL + 1) % 3;
            ERB = (ERA + 1) % 3;
            V1 = Tri.Simplex[R][ERL];
            V2 = Tri.Simplex[R][ERA];
            V3 = Tri.Simplex[R][ERB];

            // Swapping criterion
            check = SWAP_Criterion(Tri.Vertex[V1],
                                   Tri.Vertex[V2],
                                   Tri.Vertex[V3],
                                   xp);

            if (check) {

                // Triangle adiacent to the swapping region
                A = Tri.Adjacency[R][ERA];
                B = Tri.Adjacency[R][ERB];

                // Swap face
                Tri.SwapFaceTri(L, j);

                // Put new edges in the stack list
                if (A >= 0) {
                    STACK_data[0] = L;
                    STACK_data[1] = (Tri.vertex(L, P) + 1) % 3;
                    STACK.push(STACK_data);
                }
                if (B >= 0) {
                    STACK_data[0] = R;
                    STACK_data[1] = (Tri.vertex(R, P) + 1) % 3;
                    STACK.push(STACK_data);
                }
            }
        };
        
    } // next TOPSTK

} //next i

return; };

// -------------------------------------------------------------------------- //
void Delaunay::Delaunay2D::InsertEdges(
    VolTriPatch        &Tri,
    ivector2D           &edges,
    ivector2D           &ring1
) {

// ========================================================================== //
// void Delaunay::Delaunay2D::InsertEdges(                                 //
//     VolTriPatch        &Tri,                                              //
//     ivector2D           &edges,                                            //
//     ivector2D           &ring1)                                            //
//                                                                            //
// Insert constrained edges into triangulation.                               //
// Assumption:                                                                //
// - 2D conformant triangulation                                              //
// - constrained edges must not (self) intersect and fully included in the    //
//   external envelope of the vertex cloud                                    //
// - constrained edges vertex must belong to the triangulation.               //
// - each constrained edge cannot contain more than two triangulaiton vert.   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri     : VolTriPatch, volume triangulation                             //
// - edges   : ivector2D, edge list. edges[i] stores the global indices of    //
//             vertices in the triangulation                                  //
// - ring1   : ivector2D, vertex-simplex connectivity.                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double const        abs_toll = 1.0e-12;

// Local variables
bool                check, inner_check;
int                 nE = edges.size();
ivector1D           path;
ivector2D           swap_list;
array<double, 2>    u, v, d, x, eu, ev;
ivector1D::iterator it_;

// Counters
int                 i, j, k, e;
int                 n;
int                 A, B;
int                 EU, EV, U, V, E, T, S;

// ========================================================================== //
// LOOP OVER EDGE LIST                                                        //
// ========================================================================== //
for (E = 0; E < nE; ++E) {

    // Edges infos ========================================================== //
    {
        // edge's vertices
        EU = edges[E][0];
        EV = edges[E][1];
        eu[0] = Tri.Vertex[EU][0];    eu[1] = Tri.Vertex[EU][1];
        ev[0] = Tri.Vertex[EV][0];    ev[1] = Tri.Vertex[EV][1];
        
        // edge's direction
        d[0] = Tri.Vertex[EV][0] - Tri.Vertex[EU][0];
        d[1] = Tri.Vertex[EV][1] - Tri.Vertex[EU][1];
        d = d/norm2(d);

    }

    // Iteratively swap triangulation edges ================================= //
///*deb*/int swap_iter = 0;
///*deb*/stringstream nome;
///*deb*/nome << "debug_e" << E << "_swap" << swap_iter << ".vtu";
///*deb*/Tri.Export_vtu(nome.str());
///*deb*/nome.str("");
///*deb*/swap_iter++;
    check = false;
    while (!check) {

        // Initialize scope variables --------------------------------------- //
        path.resize(0);
   
        // Check if edge is already included in the triangulation ----------- //
        {
            check = false;
            n = ring1[EU].size();
            i = 0;
            while (!check && (i < n)) {
                j = Tri.vertex(ring1[EU][i], EV);
                check = (check || (j >= 0));
                i++;
            } //next i

        }

///*deb*/ if (check) {
///*deb*/     cout << "edge " << E << " is included in triangulation" << endl;
///*deb*/ }
///*deb*/ else {
///*deb*/     cout << "edge " << E << " is not included in triangulation" << endl;
///*deb*/ }

        // Force constrained edge into triangulation ------------------------ //
        if (!check) {

            // Find origin simplex
///*deb*/     cout << "  searching start simplex" << endl;
            {
                n = ring1[EU].size();
                i = 0;
                inner_check = false;
                while (!inner_check && (i < n)) {
                    T = ring1[EU][i];
                    j = Tri.vertex(T, EU);
                    j = (j + 1) % 3;
                    k = (j + 1) % 3;
                    U = Tri.Simplex[T][j];
                    V = Tri.Simplex[T][k];
                    u[0] = Tri.Vertex[U][0] - Tri.Vertex[EU][0];
                    u[1] = Tri.Vertex[U][1] - Tri.Vertex[EU][1];
                    u = u/norm2(u);
                    v[0] = Tri.Vertex[V][0] - Tri.Vertex[EU][0];
                    v[1] = Tri.Vertex[V][1] - Tri.Vertex[EU][1];
                    v = v/norm2(v);
                    inner_check = ((crossProduct(u, d) > abs_toll )
                                && (crossProduct(v, d) < -abs_toll));
                    i++;
                } //next simplex
                i--;
                A = ring1[EU][i];
            }
///*deb*/     cout << "  edge starts @simplex: " << A << endl;
        
            // Find ending simplex
///*deb*/     cout << "  searching end simplex" << endl;
            {
                n = ring1[EV].size();
                i = 0;
                inner_check = false;
                while (!inner_check && (i < n)) {
                    T = ring1[EV][i];
                    j = Tri.vertex(T, EV);
                    j = (j + 1) % 3;
                    k = (j + 1) % 3;
                    U = Tri.Simplex[T][j];
                    V = Tri.Simplex[T][k];
                    u[0] = Tri.Vertex[U][0] - Tri.Vertex[EV][0];
                    u[1] = Tri.Vertex[U][1] - Tri.Vertex[EV][1];
                    u = u/norm2(u);
                    v[0] = Tri.Vertex[V][0] - Tri.Vertex[EV][0];
                    v[1] = Tri.Vertex[V][1] - Tri.Vertex[EV][1];
                    v = v/norm2(v);
                    inner_check = ((crossProduct(u, d) < abs_toll )
                                && (crossProduct(v, d) > -abs_toll));
                    i++;
                } //next simplex
                i--;        
                B = ring1[EV][i];
            }
///*deb*/     cout << "  edge ends @simplex: " << B << endl;

            // Find path from A to B along d
///*deb*/     cout << "  finding path from " << A << " to " << B << endl;
            {
                T = A;
                e = Tri.vertex(A, EU);
                S = Tri.Adjacency[A][e];
                while (T != B) {

                    // Add T to path
                    path.push_back(T);

                    // Find face intersected by edge
                    inner_check = false;
                    e = Tri.face(T, S);
                    i = 0;
                    while (!inner_check && (i < 2)) {
                        e = (e + 1) % 3;
                        j = (e + 1) % 3;
                        U = Tri.Simplex[T][e];
                        V = Tri.Simplex[T][j];
                        u[0] = Tri.Vertex[U][0];    u[1] = Tri.Vertex[U][1];
                        v[0] = Tri.Vertex[V][0];    v[1] = Tri.Vertex[V][1];
                        inner_check = CGElem::IntersectSegmentSegment(u, v, eu, ev, x) ;
                        i++;
                    } //next face

                    S = T;
                    T = Tri.Adjacency[T][e];

                } //next T
                path.push_back(B);
            }
///*deb*/     cout << "  path is: " << path << endl;

            // Insert edge into triangulation
            i = 0;
            n = path.size();
            while (i < n - 1) {
                A = path[i];
                B = path[i+1];

///*deb*/         cout << "    next element is: " << A << endl;
///*deb*/         cout << "    its neighbor is: " << B << endl;

                // Triangles A and B form a convex polygon
                if (Delaunay::Delaunay2D::ConvexPolygon(Tri, A, B)) {
///*deb*/             cout << "    " << A << " and " << B << " form convex polygon" << endl;
                    // Determines swap arguments
                    j = Tri.face(B, A);
                    j = (j + 2) % 3;
                    U = Tri.Simplex[B][j];
                    u[0] = Tri.Vertex[U][0] - Tri.Vertex[EU][0];
                    u[1] = Tri.Vertex[U][1] - Tri.Vertex[EU][1];
                    u = u/norm2(u);

                    // Edge crosses below
                    if (crossProduct(d, u) >= abs_toll) {
///*deb*/                 cout << "      edge crosses below" << endl;
                        // Swap face
                        Tri.SwapFaceTri(A, Tri.face(A, B));
                        path[i] = B;
                        path[i+1] = A;
///*deb*/                 nome << "debug_e" << E << "_swap" << swap_iter << ".vtu";
///*deb*/                 Tri.Export_vtu(nome.str());
///*deb*/                 nome.str("");
///*deb*/                 swap_iter++;

                        // Update ring1
///*deb*/                 cout << "      updating 1-rings" << endl;
                        {

                            j = Tri.face(A, B);
                            U = Tri.Simplex[A][j];
                            ring1[U].push_back(A);

                            j = (j + 2) % 3;
                            U = Tri.Simplex[A][j];
                            it_ = find(ring1[U].begin(), ring1[U].end(), B);
                            ring1[U].erase(it_);

                            j = Tri.face(B, A);
                            U = Tri.Simplex[B][j];
                            ring1[U].push_back(B);

                            j = (j + 2) % 3;
                            U = Tri.Simplex[B][j];
                            it_ = find(ring1[U].begin(), ring1[U].end(), A);
                            ring1[U].erase(it_);
                        }
                    
                    }

                    // Edge crosses above
                    else {
///*deb*/                 cout << "      edge crosses above" << endl;
                        // Swap face
                        Tri.SwapFaceTri(A, Tri.face(A, B));
///*deb*/                 nome << "debug_e" << E << "_swap" << swap_iter << ".vtu";
///*deb*/                 Tri.Export_vtu(nome.str());
///*deb*/                 nome.str("");
///*deb*/                 swap_iter++;

                        // Update ring1
///*deb*/                 cout << "      updating 1-rings" << endl;
                        {
                            j = Tri.face(A, B);
                            U = Tri.Simplex[A][j];
                            ring1[U].push_back(A);

                            j = (j + 2) % 3;
                            U = Tri.Simplex[A][j];
                            it_ = find(ring1[U].begin(), ring1[U].end(), B);
                            ring1[U].erase(it_);

                            j = Tri.face(B, A);
                            U = Tri.Simplex[B][j];
                            ring1[U].push_back(B);

                            j = (j + 2) % 3;
                            U = Tri.Simplex[B][j];
                            it_ = find(ring1[U].begin(), ring1[U].end(), A);
                            ring1[U].erase(it_);

                        }
                    }

                }

                // Triangles A and B form a non-convex polygon
///*deb*/         else {
///*deb*/             cout << "    " << A << "and " << B << " form non-convex polygon" << endl;
///*deb*/         }

                // Next simplex on path
                i++;

            } //next simplex
        }
///*deb*/ cout << "check: " << check << endl;
    } //next swap iteration
} //next E

return; }

// -------------------------------------------------------------------------- //
void Delaunay::Delaunay2D::ReturnSimplexID(
    VolTriPatch  const &Tri,
    a3vector1D    const &P,
    int                 &seed
) {

// ========================================================================== //
// void Delaunay::Delaunay2D::ReturnSimplexID(                                //
//     VolTriPatch  const &Tri,                                              //
//     a3vector1D    const &P,                                                //
//     int                 &seed)                                             //
//                                                                            //
// Returns ID of simplex in the volume mesh R enclosing the point P.          //
// Assumptions:                                                               //
// - 2D unstructured conformant mesh                                          //
// - Mesh is composed of strictly convex simplicies (the case of non-convex   //
//   may result in undefined behaviour)                                       //
// - Mesh does not contains any cavity                                        //
// - Point P is enclosed by volume mesh.                                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri    : VolTriPatch, unstructured volume mesh                          //
// - P      : a3vector1D, point coordinates                                   //
// -.seed   : int, ID of simplex used as seed by searching algorithm.         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
static double const abs_tol = 1.0e-12;

// Local variables
bool                check = false, inner_check;
array<double, 2>    u, v;

// Counters
int                 i, j;
int                 n;
int                 U, V, T;

// ========================================================================== //
// SEARCH SIMPLEX ENCLOSING POINT P                                           //
// ========================================================================== //
T = seed;
while (!check) {

    // Update siplex ID ----------------------------------------------------- //
    seed = T;

    // Simplex infos -------------------------------------------------------- //
    n = Tri.infos[Tri.e_type[T]].n_vert;

    // Check if simplex encloses point P using cross-product rule ----------- //
    inner_check = true;
    i = 0;
    while (inner_check && (i < n)) {

        // Define edge
        j = (i + 1) % n;
        U = Tri.Simplex[T][i];
        V = Tri.Simplex[T][j];

        // Define directions
        u[0] = Tri.Vertex[V][0] - Tri.Vertex[U][0];
        u[1] = Tri.Vertex[V][1] - Tri.Vertex[U][1];
        v[0] = P[0] - Tri.Vertex[U][0];
        v[1] = P[1] - Tri.Vertex[U][1];

        // Check
        inner_check = inner_check && (crossProduct(u, v) >= -abs_tol);

        // Update edge counter
        i++;

    } //next edge
    i--;

    // Update flags --------------------------------------------------------- //
    check = inner_check;
    
    // Move to next simplex ------------------------------------------------- //
    T = Tri.Adjacency[T][i];

} //next simplex

return; }

// -------------------------------------------------------------------------- //
ivector1D Delaunay::Delaunay2D::FindPath(
    VolTriPatch        &Tri,
    ivector1D     const &edge
) {

// ========================================================================== //
// ivector1D Delaunay::Delaunay2D::FindPath(                               //
//     VolTriPatch        &Tri,                                              //
//     ivector1D     const &edge)                                             //
//                                                                            //
// Find path of mesh simplicies crossed by a given edge.                      //
// Assumptions:                                                               //
// - 2D unstructured conformant mesh                                          //
// - Mesh is composed of strictly convex simplicies (the case of non-convex   //
//   may result in undefined behaviour)                                       //
// - Mesh does not contains any cavity                                        //
// - edge is enclosed by volume mesh and its extrema are mesh vertices        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri      : VolTriPatch, volume tasselation                              //
// - edge     : ivector1D, edge. edge[i][0], edge[i][1] are the global        //
//              indices of mesh vertices defining the edge.                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - path     : ivector1D, list of triangles crossed by edge                  //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                        check;
ivector1D                   path;
array<double, 2>            x, E1, E2, P1, P2;

// Counters
int                         n;
int                         i, j, k;
int                         U, V, S, T, T0, T1;

// ========================================================================== //
// SET PARAMETERS                                                             //
// ========================================================================== //
E1[0] = Tri.Vertex[edge[0]][0];    E1[1] = Tri.Vertex[edge[0]][1];
E2[0] = Tri.Vertex[edge[1]][0];    E2[1] = Tri.Vertex[edge[1]][1];

// ========================================================================== //
// FIND SIMPLICIES ENCLOSING EDGE'S EXTREME POINTS                            //
// ========================================================================== //
T0 = 0;
ReturnSimplexID(Tri, Tri.Vertex[edge[0]], T0);
ReturnSimplexID(Tri, Tri.Vertex[edge[1]], T0);

// ========================================================================== //
// FIND PATH OF SIMPLEX CROSSED BY THE EDGE                                   //
// ========================================================================== //

// Starting triangles ------------------------------------------------------- //
T = T0;
path.push_back(T);
check = false;
i = 0;
n = Tri.infos[Tri.e_type[T]].n_vert;
while (!check && (i < n)) {
    j = (i + 1) % n;
    U = Tri.Simplex[T][i];
    V = Tri.Simplex[T][j];
    P1[0] = Tri.Vertex[U][0];    P1[1] = Tri.Vertex[U][1];
    P2[0] = Tri.Vertex[V][0];    P2[1] = Tri.Vertex[V][1];
    check = CGElem::IntersectSegmentSegment(P1, P2, E1, E2, x);
    i++;
} //next i
i--;
S = Tri.Adjacency[T][i];

// Find simplicies crossed by edge ------------------------------------------ //
while (S != T1) {

    // Push simplex into path
    path.push_back(S);

    // Look for edge to be crossed
    n = Tri.infos[Tri.e_type[S]].n_vert;
    check = false;
    k = Tri.face(S, T);
    i = (k + 1) % 3;
    while (!check && (i != k)) {
        j = (i + 1) % 3;
        U = Tri.Simplex[T][i];
        V = Tri.Simplex[T][j];
        P1[0] = Tri.Vertex[U][0];    P1[1] = Tri.Vertex[U][1];
        P2[0] = Tri.Vertex[V][0];    P2[1] = Tri.Vertex[V][1];
        check = CGElem::IntersectSegmentSegment(P1, P2, E1, E2, x);
        i = (i + 1) % n;
    }
    i = (n + i - 1) % n;

    // Next simplex
    T = S;
    S = Tri.Adjacency[T][i];

} //next Simplex
path.push_back(T1);

return(path); }

// -------------------------------------------------------------------------- //
unsigned int Delaunay::Delaunay2D::Delaunay2D(
    VolTriPatch        &Tri
) {

// ========================================================================== //
// unsigned int Delaunay::Delaunay2D::Delaunay2D(                          //
//     VolTriPatch        &Tri)                                              //
//                                                                            //
// Delaunay triangulation algorithm. Compute the Delaunay triangulation of a  //
// set of vertices in 2D.                                                     //
// Assumptions:                                                               //
// - 2D space                                                                 //
// - No repeated vertices are allowed.                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri       : VolTriPatch, volume mesh manager                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err        : unsigned int, error flag                                    //
//                err = 0    --> no errors encounterd during meshing process  //
//                err = 1    --> empty vertex list.                           //
//                err = 2    --> an error occured during mesh process         //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
int     const               dim = Tri.Vertex[0].size();
double  const               toll = 1.0e-12;
double  const               c00000 = 0.0;
double  const               c00100 = 10.0;

// Local variables
bool                        check;
double                      dx, dy, d;
ivector1D                   idummy1D(3, -1);
ivector1D                   vorder;
a3vector1D                  xp;
array<double, 2>            xlim, ylim;

// Counters
int                         i, j, T;
int                         W1, W2, W3;

// ========================================================================== //
// CHECK INPUT GEOMETRY                                                       //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Check input ---------------------------------------------------------- //
    if (Tri.nVertex == 0) {
        return(1);
    }

}

// ========================================================================== //
// NORMALIZE INPUT COORDINATES                                                //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Normalize input vertices --------------------------------------------- //

    // Vertex bounding box
    xlim[0] = xlim[1] = Tri.Vertex[0][0];
    ylim[0] = ylim[1] = Tri.Vertex[0][1];
    for (T = 1; T < Tri.nVertex; ++T) {
        xlim[0] = min(xlim[0], Tri.Vertex[T][0]);
        ylim[0] = min(ylim[0], Tri.Vertex[T][1]);
        xlim[1] = max(xlim[1], Tri.Vertex[T][0]);
        ylim[1] = max(ylim[1], Tri.Vertex[T][1]);
    } //next T

    // Bounding box dimensions
    dx = xlim[1] - xlim[0];
    dy = ylim[1] - ylim[0];
    d = max(dx, dy);

    // Check for degenerate cases
    if (d < toll) { return(2); }

    // Normalize input vertices in [0,1]
    Tri.Translate(-xlim[0], -ylim[0], 0.0);
    Tri.Scale(1.0/d, 1.0/d, 1.0);

}

// ========================================================================== //
// SORT VERTICES                                                              //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    int                         n_bins = 128, bin_index, label;
    bitpit::MinPQueue<int, int> vlist(Tri.nVertex, true);

    // Sort vertices on regular bins ---------------------------------------- //
    dx = 1.0/((double) n_bins);
    dy = 1.0/((double) n_bins);
    for (T = 0; T < Tri.nVertex; ++T) {
        i = Tri.Vertex[T][0]/dx;
        j = Tri.Vertex[T][1]/dy;
        if (j%2 == 0) {
            vlist.keys[T] = n_bins * j + i;
        }
        else {
            vlist.keys[T] = n_bins *j + (n_bins - i);
        }
    } //next T
    vlist.heap_size = Tri.nVertex;

    // Initialize vertex labels --------------------------------------------- //
    for (i = 0; i < Tri.nVertex; ++i) {
        vlist.labels[i] = i;
    } //next i

    // Quick sort algorithm ------------------------------------------------- //
    vlist.buildHeap();

    // Sort algorithm ------------------------------------------------------- //
    vorder.resize(Tri.nVertex);
    i = 0;
    while (vlist.heap_size > 0) {
        vlist.extract(bin_index, label);
        vorder[i] = label;
        i++;
    } //next item

}

// ========================================================================== //
// RESIZE DATA STRUCTURE                                                      //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Resize connectivity data structure ----------------------------------- //

    // Estimated numbe of triangles
    Tri.nSimplex = 2*Tri.nVertex + 1;

    // Reshape simplex-vertex connectivity
    Tri.ReshapeSimplex(5);

    // Resize adjacency data structure -------------------------------------- //
    Tri.ReshapeAdjacency();
    Tri.nSimplex = 0;
}

// ========================================================================== //
// PLACE THE SUPER-TRIANGLE                                                   //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Create the outer triangle -------------------------------------------- //

    // 1st vertex
    xp[0] = -c00100;
    xp[1] = -c00100;
    Tri.AddVertex(xp);
    W1 = Tri.nVertex-1;

    // 2nd vertex
    xp[0] = c00100;
    xp[1] = -c00100;
    Tri.AddVertex(xp);
    W2 = Tri.nVertex-1;

    // 3rd vertex
    xp[0] = c00000;
    xp[1] = c00100;
    Tri.AddVertex(xp);
    W3 = Tri.nVertex-1;

    // Vertex-triangle connectivity for the 1st triangle -------------------- //
    idummy1D[0] = W1;
    idummy1D[1] = W2;
    idummy1D[2] = W3;
    Tri.AddSimplex(idummy1D, 5);

    // Adjacencies for the 1st triangle ------------------------------------- //
    Tri.Adjacency[Tri.nSimplex] = ivector1D(3, -1);

}

// ========================================================================== //
// ADD VERTICES TO THE TRIANGULATION                                          //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Add vertices to the triangulation ------------------------------------ //
    Delaunay2D::InsertVertices(Tri, vorder);
}

// ========================================================================== //
// CHECK TRIANGULATION CONSISTENCY                                            //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Consistency check ---------------------------------------------------- //
    if (Tri.nSimplex != 2*(Tri.nVertex -3) + 1) {
        return (2);
    }
}

// ========================================================================== //
// CLEAN TRIANGULATION                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    ivector1D           removable_list;

    // Mark removable simplex ----------------------------------------------- //
    for (T = 0; T < Tri.nSimplex; T++) {
        i = 0;
        check = true;
        while (check && (i < 3)) {
            if ((Tri.Simplex[T][i] == W1)
            || (Tri.Simplex[T][i] == W2)
            || (Tri.Simplex[T][i] == W3)) {
                removable_list.push_back(T);
                check = false;
            }
            i++;
        }
    } //next T

    // Remove simplex ------------------------------------------------------- //
    Tri.RemoveSimplex(removable_list);

    // Remove free vertices ------------------------------------------------- //
    Tri.RemoveIsolatedVertex();

    // // Resize data structure --------------------------------------------- //
    Tri.ResizeVertex();
    Tri.ResizeSimplex();
    Tri.ResizeAdjacency();

}

// ========================================================================== //
// RETURN TO ORIGINAL COORDINATES                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Remap vertices into the original coordinates ------------------------- //
    Tri.Scale(d, d, 1.0);
    Tri.Translate(xlim[0], ylim[0], 0.0);

}

return(0); };

// -------------------------------------------------------------------------- //
unsigned int Delaunay::Delaunay2D::CDelaunay2D(
    VolTriPatch        &Tri,
    ivector2D           *domain,
    vector<ivector2D*>  &holes,
    vector<ivector2D*>  &segments
) {

// ========================================================================== //
// unsigned int Delaunay::Delaunay2D::CDelaunay2D(                         //
//     VolTriPatch        &Tri,                                              //
//     ivector2D           *domain,                                           //
//     vector<ivector2D*>  &holes,                                            //
//     vector<ivector2D*>  &segments)                                         //
//                                                                            //
// Compute the constrained Delaunay triangulation of a set of 2D vertices.    //
// Assumptions:                                                               //
// - 2D space                                                                 //
// - no repeated vertices are allowed                                         //
// - constrained edges must not (self) intersect and fully included in the    //
//   external envelope of the vertex cloud                                    //
// - constrained edges vertex must belong to the triangulation.               //
// - each constrained edge cannot contain more than two triangulaiton vert.   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri      : VolTriPatch, volume mesh                                     //
// - domain   : ivector2D*, pointer to edge-vertex connectivity for domain    //
//              boundaries.                                                   //
// - holes    : vector<ivector2D*>, vector of pointer to edge-vertex          //
//              connectivity for each domain hole                             //
// - segments : vector<ivector2D*>, vector of pointer to edge-vertex          //
//              connectivity for each domain segment                          //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err        : unsigned int, error flag                                    //
//                err = 0    --> no errors encounterd during meshing process  //
//                err = 1    --> empty vertex list.                           //
//                err = 2    --> an error occured during mesh process         //
// ========================================================================== //

// Parameters
int     const               dim = Tri.Vertex[0].size();
double  const               toll = 1.0e-12;
double  const               c00000 = 0.0;
double  const               c00100 = 10.0;

// Local variables
bool                        check;
int                         n_holes = holes.size();
int                         n_segments = segments.size();
double                      dx, dy, d;
ivector1D                   idummy1D(3, -1);
ivector1D                   vorder;
a3vector1D                  xp;
array<double, 2>            xlim, ylim;
ivector2D                   ring1;

// Counters
int                         i, j, T;
int                         W1, W2, W3;

// ========================================================================== //
// CHECK INPUT GEOMETRY                                                       //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Check input ---------------------------------------------------------- //
    if (Tri.nVertex == 0) {
        return(1);
    }

}

// ========================================================================== //
// NORMALIZE INPUT COORDINATES                                                //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Normalize input vertices --------------------------------------------- //

    // Vertex bounding box
    xlim[0] = xlim[1] = Tri.Vertex[0][0];
    ylim[0] = ylim[1] = Tri.Vertex[0][1];
    for (T = 1; T < Tri.nVertex; ++T) {
        xlim[0] = min(xlim[0], Tri.Vertex[T][0]);
        ylim[0] = min(ylim[0], Tri.Vertex[T][1]);
        xlim[1] = max(xlim[1], Tri.Vertex[T][0]);
        ylim[1] = max(ylim[1], Tri.Vertex[T][1]);
    } //next T

    // Bounding box dimensions
    dx = xlim[1] - xlim[0];
    dy = ylim[1] - ylim[0];
    d = max(dx, dy);

    // Check for degenerate cases
    if (d < toll) { return(2); }

    // Normalize input vertices in [0,1]
    Tri.Translate(-xlim[0], -ylim[0], 0.0);
    Tri.Scale(1.0/d, 1.0/d, 1.0);

}

// ========================================================================== //
// SORT VERTICES                                                              //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    int                         n_bins = 128, bin_index, label;
    bitpit::MinPQueue<int, int> vlist(Tri.nVertex, true);

    // Sort vertices on regular bins ---------------------------------------- //
    dx = 1.0/((double) n_bins);
    dy = 1.0/((double) n_bins);
    for (T = 0; T < Tri.nVertex; ++T) {
        i = Tri.Vertex[T][0]/dx;
        j = Tri.Vertex[T][1]/dy;
        if (j%2 == 0) {
            vlist.keys[T] = n_bins * j + i;
        }
        else {
            vlist.keys[T] = n_bins *j + (n_bins - i);
        }
    } //next T
    vlist.heap_size = Tri.nVertex;

    // Initialize vertex labels --------------------------------------------- //
    for (i = 0; i < Tri.nVertex; ++i) {
        vlist.labels[i] = i;
    } //next i

    // Quick sort algorithm ------------------------------------------------- //
    vlist.buildHeap();

    // Sort algorithm ------------------------------------------------------- //
    vorder.resize(Tri.nVertex);
    i = 0;
    while (vlist.heap_size > 0) {
        vlist.extract(bin_index, label);
        vorder[i] = label;
        i++;
    } //next item

}

// ========================================================================== //
// RESIZE DATA STRUCTURE                                                      //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Resize connectivity data structure ----------------------------------- //

    // Estimated numbe of triangles
    Tri.nSimplex = 2*Tri.nVertex + 1;

    // Reshape simplex-vertex connectivity
    Tri.ReshapeSimplex(5);

    // Resize adjacency data structure -------------------------------------- //
    Tri.ReshapeAdjacency();
    Tri.nSimplex = 0;
}

// ========================================================================== //
// PLACE THE SUPER-TRIANGLE                                                   //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Create the outer triangle -------------------------------------------- //

    // 1st vertex
    xp[0] = -c00100;
    xp[1] = -c00100;
    Tri.AddVertex(xp);
    W1 = Tri.nVertex-1;

    // 2nd vertex
    xp[0] = c00100;
    xp[1] = -c00100;
    Tri.AddVertex(xp);
    W2 = Tri.nVertex-1;

    // 3rd vertex
    xp[0] = c00000;
    xp[1] = c00100;
    Tri.AddVertex(xp);
    W3 = Tri.nVertex-1;

    // Vertex-triangle connectivity for the 1st triangle -------------------- //
    idummy1D[0] = W1;
    idummy1D[1] = W2;
    idummy1D[2] = W3;
    Tri.AddSimplex(idummy1D, 5);

    // Adjacencies for the 1st triangle ------------------------------------- //
    Tri.Adjacency[Tri.nSimplex] = ivector1D(3, -1);

}

// ========================================================================== //
// ADD VERTICES TO THE TRIANGULATION                                          //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Add vertices to the triangulation ------------------------------------ //
    Delaunay2D::InsertVertices(Tri, vorder);
}

// ========================================================================== //
// CHECK TRIANGULATION CONSISTENCY                                            //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Consistency check ---------------------------------------------------- //
    if (Tri.nSimplex != 2*(Tri.nVertex -3) + 1) {
        return (2);
    }

    // Build vertex 1-rings ------------------------------------------------- //
    ring1.resize(Tri.nVertex);
    for (T = 0; T < Tri.nSimplex; ++T) {
        for (i = 0; i < 3; ++i) {
            ring1[Tri.Simplex[T][i]].push_back(T);
        } //next i
    } //next T
}

// ========================================================================== //
// FORCE CONSTRAINTS INTO EXISTING TRIANGULATION                              //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Insert edges into triangulation -------------------------------------- //

    // Domain boundaries
    Delaunay2D::InsertEdges(Tri, *domain, ring1);

    // Domain holes
    for (i = 0; i < n_holes; ++i) {
        Delaunay2D::InsertEdges(Tri, *(holes[i]), ring1);
    } //next i

    // Domain segments
    for (i = 0; i < n_segments; ++i) {
        Delaunay2D::InsertEdges(Tri, *(segments[i]), ring1);
    } //next i
    
}

// ========================================================================== //
// CLEAN TRIANGULATION                                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                    n_s;
    int                    J;
    int                    n, m;
    int                    EU, EV, A, B;
    vector<short int>      flag(Tri.nSimplex, 0);
    ivector1D              neigh_s(2, -1);
    ivector1D              removable_list;
    bitpit::LIFOStack<int> stack(Tri.nSimplex);

    // Flag simplicies adjacent to constrained edges ------------------------ //

    // domain boundaries
    n_s = (*domain).size();
    for (T = 0; T < n_s; ++T) {

        // edge infos
        EU = (*domain)[T][0];
        EV = (*domain)[T][1];
        n = ring1[EU].size();

        // find simplicies adjacent to constrained edge
        j = 0;
        for (i = 0; i < n; ++i) {
            A = ring1[EU][i];
            m = Tri.vertex(A, EV);
            if (m >= 0) {
                neigh_s[j] = A;
                j++;
            }
        } // next i

        // flag simplicies
        for (i = 0; i < 2; ++i) {
            A = neigh_s[i];
            m = Tri.vertex(A, EU);
            n = (m + 1) % 3;
            if (Tri.Simplex[A][n] == EV) {
                flag[A] = 1;
            }
            else {
                flag[A] = 2;
            }
        } //next i
    } //next T

    // domain holes
    for (J = 0; J < n_holes; ++J) {
        n_s = (*(holes[J])).size();
        for (T = 0; T < n_s; ++T) {

            // edge infos
            EU = (*(holes[J]))[T][0];
            EV = (*(holes[J]))[T][1];
            n = ring1[EU].size();

            // find simplicies adjacent to constrained edge
            j = 0;
            for (i = 0; i < n; ++i) {
                A = ring1[EU][i];
                m = Tri.vertex(A, EV);
                if (m >= 0) {
                    neigh_s[j] = A;
                    j++;
                }
            } // next i

            // flag simplicies
            for (i = 0; i < 2; ++i) {
                A = neigh_s[i];
                m = Tri.vertex(A, EU);
                n = (m + 1) % 3;
                if (Tri.Simplex[A][n] == EV) {
                    flag[A] = 1;
                }
                else {
                    flag[A] = 2;
                }
            } //next i
        } //next T
    } //next J

    // domain segments
    for (J = 0; J < n_segments; ++J) {
        n_s = (*(segments[J])).size();
        for (T = 0; T < n_s; ++T) {

            // edge infos
            EU = (*(segments[J]))[T][0];
            EV = (*(segments[J]))[T][1];
            n = ring1[EU].size();

            // find simplicies adjacent to constrained edge
            j = 0;
            for (i = 0; i < n; ++i) {
                A = ring1[EU][i];
                m = Tri.vertex(A, EV);
                if (m >= 0) {
                    neigh_s[j] = A;
                    j++;
                }
            } // next i

            // Cut triangulation
            A = neigh_s[0];
            B = neigh_s[1];
            j = Tri.face(A, B);
            Tri.Adjacency[A][j] = -1;
            j = Tri.face(B, A);
            Tri.Adjacency[B][j] = -1;

        } //next T
    } //next J

    // Propagate flags ------------------------------------------------------ //

    // Initialize stack
    for (T = 0; T < Tri.nSimplex; ++T) {
        if (flag[T] != 0) {
            stack.push(T);
        }
    } //next T

    // Propagate flags
    while (stack.TOPSTK > 0) {

        // pop item from stack
        T = stack.pop();

        // Loop over neighbors
        for (i = 0; i < 3; ++i) {
            A = Tri.Adjacency[T][i];
            if (A >= 0) {
                if (flag[A] == 0) {
                    flag[A] = flag[T];
                    stack.push(A);
                }
            }
        } //next i
    } //next item

    // Remove triangles with 2-flag ----------------------------------------- //
    removable_list.resize(count(flag.begin(), flag.end(), 2), -1);
    i = 0;
    for (T = 0; T < Tri.nSimplex; ++T) {
        if (flag[T] == 2) {
            removable_list[i] = T;
            i++;
        }
    } //next T
    Tri.RemoveSimplex(removable_list);

    // Remove isolated vertices --------------------------------------------- //
    Tri.RemoveIsolatedVertex();

    // Resize data structure ------------------------------------------------ //
    Tri.ResizeVertex();
    Tri.ResizeSimplex();
    Tri.ResizeAdjacency();

}

// ========================================================================== //
// RETURN TO ORIGINAL COORDINATES                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Remap vertices into the original coordinates ------------------------- //
    Tri.Scale(d, d, 1.0);
    Tri.Translate(xlim[0], ylim[0], 0.0);

}

return(0); }

// -------------------------------------------------------------------------- //
unsigned int Delaunay::Delaunay2D::MAT_distance(
    a3vector2D          &vertex,
    SurfTriPatch       &domain,
    Svector1D           &holes,
    Svector1D           &segments,
    double               h
) {

// ========================================================================== //
// unsigned int Delaunay::Delaunay2D::MAT_distance(                           //
//     a3vector2D          &vertex,                                           //
//     SurfTriPatch       &domain,                                           //
//     Svector1D           &holes,                                            //
//     Svector1D           &segments,                                         //
//     double               h)                                                //
//                                                                            //
// Compute distance from approximated medial axis for a closed domain.        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - vertex   : a3vector2D, external vertex list. vertex[i][0], vertex[i][1]  //
//              are the x, y coordinates of the i-th vertex.                  //
// - domain   : SurfTriPatch, surface mesh of domain boundaries              //
// - holes    : Svector1D, holes[i] stores the surface mesh for the i-th      //
//              domain hole                                                   //
// - segments : Svector1D, segments[i] store the surface mesh for the i-th    //
//              domain segment                                                //
// - h        : double, refinement parameter.                                 //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err      : unsigned int, error flag:                                     //
//                err = 0    --> no errors encounterd during meshing process  //
//                err = 1    --> empty vertex list.                           //
//                err = 2    --> an error occured during mesh process         //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                             nV, nH = holes.size(), nS = segments.size();
int                             nD_new;
ivector1D                       nH_new(nH, -1), nS_new(nS, -1);
VolTriPatch                    Mesh;

// Tmp variables
SurfTriPatch                   MAT;

// Counters
int                             iter;
int                             n;
int                             T, H, S;

// ========================================================================== //
// BINARY REFINAMENT OF SURFACE MESH                                          //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Domain boundaries ---------------------------------------------------- //
    domain.BinaryRefinement(vertex, h);
    
    // Domain holes --------------------------------------------------------- //
    for (H = 0; H < nH; ++H) {
        holes[H].BinaryRefinement(vertex, h);
    } //next H
    
    // Domain segments ------------------------------------------------------ //
    for (S = 0; S < nS; ++S) {
        segments[S].BinaryRefinement(vertex, h);
    } //next S

    // Debug only (export refined surface mesh into .vtk format) ------------ //
    {

        // Scope variables
        SurfTriPatch       Tmp;

        // Reconstruct domain
        Tmp.AddVertices(vertex);
        Tmp.AddSimplicies(domain.Simplex);
        for (H = 0; H < nH; ++H) {
            Tmp.AddSimplicies(holes[H].Simplex);
        } //next H
        for (S = 0; S < nS; ++S) {
            Tmp.AddSimplicies(segments[S].Simplex);
        } //next S

        // Export domain
        // Tmp.Export_vtu("domain_refined.vtu");
    }

}

// ========================================================================== //
// COMPUTE CONSTRAINED DELAUNAY TRIANGULATION                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    ivector2D*          domain_;
    vector<ivector2D*>  holes_(nH);
    vector<ivector2D*>  segments_(nS);

    // Initialize scope variables ------------------------------------------- //

    // Iniitialize pointers to constrained edges
    domain_ = &(domain.Simplex);
    for (H = 0; H < nH; ++H) {
        holes_[H] = &(holes[H].Simplex);
    } //next H
    for (S = 0; S < nS; ++S) {
        segments_[S] = &(segments[S].Simplex);
    } //next H

    // Initialize vertex list for volume mesh
    Mesh.Vertex.resize(nV);
    Mesh.AddVertices(vertex);

    // Compute constrained Delaunay triangulation --------------------------- //
    Delaunay::Delaunay2D::CDelaunay2D(Mesh, domain_, holes_, segments_);

    // Debug only (Export delaunay triangulation in .vtk format) ------------ //
    // Mesh.Export_vtu("vmesh.vtu");

}

// ========================================================================== //
// COMPUTE DISTANCE FROM MEDIAL AXIS                                          //
// ========================================================================== //
{
    // Scope variables ====================================================== //
    bool                check;
    int                 i, j, k;
    int                 U, V, W;
    int                 n_faces;
    double              s;
    array<double, 2>    u, v;
    ivector1D           idummy1D(2, -1);
    a3vector1D          xC, P;

    // Initialize scope variables =========================================== //
    MAT.nSimplex    = 3*Mesh.nSimplex;
    MAT.nVertex     = 4*Mesh.nSimplex;
    MAT.ResizeVertex();
    MAT.nVertex = 0;
    MAT.ResizeSimplex(2);
    MAT.nSimplex = 0;

    // Loop over simplicies in the volume mesh ============================== //
    for (T = 0; T < Mesh.nSimplex; ++T) {

        // Count the number of connected faces for triangle T --------------- //
        n_faces = 0;
        for (i = 0; i < 3; ++i) {
            if (Mesh.Adjacency[T][i] >= 0) {
                n_faces++;
            }
        } //next i

        // Case 0: 3 connected faces ---------------------------------------- //
        if (n_faces == 3) {
            // Compute triangle's circumcenter
            xC = Mesh.CircumCenter(T);

            // Check if circum center is internal to triangle T
            k = -1;
            i = 0;
            check = true;
            while ((i < 3) && check) {
                j = (i+1) % 3;
                U = Mesh.Simplex[T][i];
                V = Mesh.Simplex[T][j];
                u[0] = Mesh.Vertex[V][0] - Mesh.Vertex[U][0];
                u[1] = Mesh.Vertex[V][1] - Mesh.Vertex[U][1];
                u = u/norm2(u);
                v[0] = xC[0] - Mesh.Vertex[U][0];
                v[1] = xC[1] - Mesh.Vertex[U][1];
                s = crossProduct(u, v);
                check = (check && (s > 0.0));
                if (!check) { k = i; }
                i++;
            } //next i

            // Case 0.1: Circumcenter is inside triangle
            if (check) {

                // Reconstruct medial axis
                MAT.AddVertex(xC);
                V = MAT.nVertex-1;
                for (i = 0; i < 3; ++i) {
                    if (Mesh.Adjacency[T][i] >= 0) {
                        P = Mesh.FaceCenter(T, i);
                        MAT.AddVertex(P);
                        idummy1D[0] = V;
                        idummy1D[1] = MAT.nVertex-1;
                        MAT.AddSimplex(idummy1D);
                    }
                } //next i
            }
            // Case 0.1: Circumcenter is outside triangle
            else {
                cout << "circumc. is outside" << endl;
                cout << "on simplex " << T << endl;
                cout << "circumcenter is: " << xC << endl;
                cout << "external to face: " << k << endl;
                cout << "starting @vertex: " << Mesh.Simplex[T][k] << endl;
                P = Mesh.FaceCenter(T, k);
                MAT.AddVertex(P);
                V = MAT.nVertex-1;
                for (i = 0; i < 3; ++i) {
                    if (i != k) {
                        if (Mesh.Adjacency[T][i] >= 0) { 
                            P = Mesh.FaceCenter(T, i);
                            MAT.AddVertex(P);
                            idummy1D[0] = V;
                            idummy1D[1] = MAT.nVertex-1;
                            MAT.AddSimplex(idummy1D);
                        }
                    }
                } //next i
            }
        }

        // Case 1: 2 connected faces ---------------------------------------- //
        else if (n_faces == 2) {

            // find non-connected face
            i = 0;
            check = true;
            while (check && (i < 3)) {
                check = (Mesh.Adjacency[T][i] >= 0);
                i++;
            } //next i
            i--;
            i = (i+1) % 3;
            P = Mesh.FaceCenter(T, i);
            MAT.AddVertex(P);
            V = MAT.nVertex-1;
            i = (i+1) % 3;
            P = Mesh.FaceCenter(T, i);
            MAT.AddVertex(P);
            idummy1D[0] = V;
            idummy1D[1] = MAT.nVertex-1;
            MAT.AddSimplex(idummy1D);

        }

        // Case 2: 1 connected faces ---------------------------------------- //
        else if (n_faces == 1) {
            xC = Mesh.Baricenter(T);
            MAT.AddVertex(xC);
            V = MAT.nVertex-1;
            check = true;
            i = 0;
            while (check && (i < 3)) {
                check = (Mesh.Adjacency[T][i] < 0);
                i++;
            } //next i
            i--;
            P = Mesh.FaceCenter(T, i);
            MAT.AddVertex(P);
            idummy1D[0] = V;
            idummy1D[1] = MAT.nVertex-1;
            MAT.AddSimplex(idummy1D);
        }

    } //next T
}

/*debug*/
//MAT.Export_vtu("MAT.vtu");

return(0); }

// -------------------------------------------------------------------------- //
void Delaunay::Delaunay2D::CCDT(
    VolTriPatch        &Tri,
    bvector1D           &movable,
    ivector2D           &Rings1
) {

// ========================================================================== //
// void Delaunay::Delaunay2D::CCDT(                                        //
//     VolTriPatch        &Tri,                                              //
//     bvector1D           &movable,                                          //
//     ivector2D           &Rings1)                                           //
//                                                                            //
// Optimize vertex location using the CCDT (Constrained Capacity Delaunay     //
// Triangulation) algorithm.                                                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri          : VolTriPatch, unstructured volume mesh manager            //
// - movable      : bvector1D, flags for movable vertices. If movable[i] =    //
//                  true, than the i-th vertex can be displaced.              //
// - Rings1       : ivector2D, with 1-ring of simplicies for each             //
//                  mesh vertex                                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
// none

// Local variables
int             max_iter = 20;
int             n_ele, n_vert;
double          A, B, C, detAA;
dvector1D       a, b, c;
dvector1D       S, sol(2, 0.0), BB(2, 0.0);
dvector2D       R, Rt, AA(2, dvector1D(2, 0.0));

// Counters
int             i, j, k, iter, U, V, W, T;

// ========================================================================== //
// VERTEX DISPLACEMENT                                                        //
// ========================================================================== //
for (iter = 0; iter < max_iter; iter++) {
    for (U = 0; U < Tri.nVertex; U++) {
        if (movable[U]) {

            // number of simplicies in the 1-ring of vertex V --------------- //
            n_ele = Rings1[U].size();

            // Resize variables --------------------------------------------- //
            a.resize(n_ele, 0.0);
            b.resize(n_ele, 0.0);
            c.resize(n_ele, 0.0);
            R.resize(n_ele, dvector1D(2, 0.0));
            S.resize(n_ele, 0.0);

            // Initialize variables ----------------------------------------- //
            A = B = C = 0.0;

            // Build matrix coeffs ------------------------------------------ //
            for (i = 0; i < n_ele; i++) {

                // i.th simplex in the 1-ring
                T = Rings1[U][i];

                // Find local index of the other vertices in simplex T
                j = Tri.vertex(T, U);
                n_vert = Tri.Simplex[T].size();
                j = (j+1) % n_vert;
                k = (j+1) % n_vert;

                // Vertices global index
                V = Tri.Simplex[T][j];
                W = Tri.Simplex[T][k];

                // Compute coeffs.
                a[i] = Tri.Vertex[V][1] - Tri.Vertex[W][1];
                b[i] = Tri.Vertex[W][0] - Tri.Vertex[V][0];
                c[i] = Tri.Vertex[V][0]*Tri.Vertex[W][1]
                     - Tri.Vertex[W][0]*Tri.Vertex[V][1];
                A += 1.0/((double) n_ele) * a[i];
                B += 1.0/((double) n_ele) * b[i];
                C += 1.0/((double) n_ele) * c[i];
            } //next i

            // Build matrix ------------------------------------------------- //
            for (i = 0; i < n_ele; i++) {
                R[i][0] = a[i] - A;
                R[i][1] = b[i] - B;
                S[i] = c[i] - C;
            } //next i
            bitpit::linearalgebra::transpose(R, Rt);
            bitpit::linearalgebra::matmul(Rt, R, AA);
            bitpit::linearalgebra::matmul(Rt, S, BB);
            BB = -1.0*BB;
            Rt.resize(0);

            // Solve the linear system AA * x = BB -------------------------- //
            detAA = bitpit::linearalgebra::det(AA);
            if (abs(detAA) > 1.0e-12) {
                bitpit::linearalgebra::cramer(AA, BB, sol);
            }
            else {
                sol[0] = Tri.Vertex[U][0];
                sol[1] = Tri.Vertex[U][1];
            }

            // Update vert location ----------------------------------------- //
            Tri.Vertex[U][0] = sol[0];
            Tri.Vertex[U][1] = sol[1];
        }
    } //next U
} //next iter

return; };

// -------------------------------------------------------------------------- //
void Delaunay::Delaunay2D::Smooth(
    VolTriPatch        &VMesh,
    bvector1D           &movable
) {

// ========================================================================== //
// void Delaunay::Delaunay2D::Smooth(                                         //
//     VolTriPatch        &VMesh,                                            //
//     bvector1D           &movable)                                          //
//                                                                            //
// Optimize vertex distribution for a Delaunay triangulation                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - VMesh     : VolTriPatch, volume mesh.                                   //
// - movable   : bvector1D, flag for movalble vertices. If                    //
//               movable[i] = true, then the i-th vertex can be               //
//               displaced.                                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
int             max_iter = 50;

// Local variables
bool            loop_continue = true;
int             err = 0;
int             n_vert, n_ele;
ivector2D       Ring1(VMesh.nVertex);

// Counters
int             i, T, V, iter = 0;

// ========================================================================== //
// OPTIMAZE VERTEX DISTRIBUTION                                               //
// ========================================================================== //
while (loop_continue) {

    // Reset triangulation -------------------------------------------------- //
    VMesh.nSimplex = 0;
    VMesh.ResizeSimplex();
    VMesh.ResizeAdjacency();

    // Generate a new triangulation ----------------------------------------- //
    err = Delaunay::Delaunay2D::Delaunay2D(VMesh);

    // Compute vertices 1-ring ---------------------------------------------- //

    // Reset 1-ring
    for (V = 0; V < VMesh.nVertex; ++V) {
        Ring1[V].resize(0);
    } //next V

    // Compute 1-ring
    for (T = 0; T < VMesh.nSimplex; T++) {

        // Number of vertices in simplex T
        n_vert = VMesh.infos[VMesh.e_type[T]].n_vert;

        // Loop over vertices in simplex T
        for (i = 0; i < n_vert; i++) {

            // Vertex global index
            V = VMesh.Simplex[T][i];
            Ring1[V].push_back(T);

        } //next i
    } //next T

    // Fixed iteration loop for vertex optimization ------------------------- //
    Delaunay::Delaunay2D::CCDT(VMesh, movable, Ring1);

    // Stopping ctiterion --------------------------------------------------- //
    iter++;
    loop_continue = (iter < max_iter);

} //next iter

return; };

// -------------------------------------------------------------------------- //
void Delaunay::Delaunay2D::Smooth(
    VolTriPatch        &VMesh,
    ivector2D           &edges,
    bvector1D           &movable
) {

// ========================================================================== //
// void Delaunay::Delaunay2D::Smooth(                                      //
//     VolTriPatch        &VMesh,                                            //
//     ivector2D           &edges,                                            //
//     bvector1D           &movable)                                          //
//                                                                            //
// Optimize vertex distribution for a constrained Delaunay triangulation      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - VMesh     : VolTriPatch, volume mesh.                                   //
// - edges     : ivector2D, edges-vertex connectivity                         //
// - movable   : bvector1D, flag for movalble vertices. If                    //
//               movable[i] = true, then the i-th vertex can be               //
//               displaced.                                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
int             max_iter = 50;

// Local variables
bool            loop_continue = true;
int             err = 0;
int             n_vert, n_ele;
ivector2D       Ring1(VMesh.nVertex);

// Counters
int             i, T, V, iter = 0;

// ========================================================================== //
// OPTIMAZE VERTEX DISTRIBUTION                                               //
// ========================================================================== //
while (loop_continue) {

    // Reset triangulation -------------------------------------------------- //
    VMesh.nSimplex = 0;
    VMesh.ResizeSimplex();
    VMesh.ResizeAdjacency();

    // Generate a new triangulation ----------------------------------------- //
    //err = Delaunay::Delaunay2D::CDelaunay2D(VMesh, edges);

    // Compute vertices 1-ring ---------------------------------------------- //

    // Reset 1-ring
    for (V = 0; V < VMesh.nVertex; ++V) {
        Ring1[V].resize(0);
    } //next V

    // Compute 1-ring
    for (T = 0; T < VMesh.nSimplex; T++) {

        // Number of vertices in simplex T
        n_vert = VMesh.infos[VMesh.e_type[T]].n_vert;

        // Loop over vertices in simplex T
        for (i = 0; i < n_vert; i++) {

            // Vertex global index
            V = VMesh.Simplex[T][i];
            Ring1[V].push_back(T);

        } //next i
    } //next T

    // Fixed iteration loop for vertex optimization ------------------------- //
    Delaunay::Delaunay2D::CCDT(VMesh, movable, Ring1);

    // Stopping ctiterion --------------------------------------------------- //
    iter++;
    loop_continue = (iter < max_iter);

} //next iter

return; };


