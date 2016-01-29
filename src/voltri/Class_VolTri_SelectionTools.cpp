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

// Selection tools ========================================================== //

// -------------------------------------------------------------------------- //
int Class_VolTri::vertex(
    int          A,
    int          V
) {

// ========================================================================== //
// int Class_VolTri::vertex(                                                  //
//     int          A,                                                        //
//     int          V)                                                        //
//                                                                            //
// Return the local index on simplex A of the vertex with global index V.     //
// If such a vertex is not found, returns -1.                                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A     : int, simplex global index                                        //
// - V     : int, vertex global index                                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - vert  : int, local index on simplex A of vertex with global index V      //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check = true;
int                 n, vert = -1;

// Counters
int                 j;

// ========================================================================== //
// FIND VERTEX WITH SPECIFIED GLOBAL INDEX                                    //
// ========================================================================== //
n = infos[e_type[A]].n_vert;
j = 0;
while ((j < n) && check) {
    if (Simplex[A][j] == V) {
        vert = j;
        check = false;
    }
    j++;
} //next j

return(vert); };

// -------------------------------------------------------------------------- //
int Class_VolTri::face(
    int          A,
    int          B
) {

// ========================================================================== //
// int Class_VolTri::face(                                                    //
//     int          A,                                                        //
//     int          B)                                                        //
//                                                                            //
// Returns the local index on simplex A of the face shared by simplcies A     //
// and B. If not such face is found, returns -1.                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A    : int, 1st simplex global index                                     //
// - B    : int, 2nd simplex global index                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - face : int, face local index                                             //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool        flag;
int         face = -1;

// Counters
int         i;

// ========================================================================== //
// SELECT FACE                                                                //
// ========================================================================== //
flag = true;
i = 0;
while (flag && (i < Adjacency[A].size())) {
    if (Adjacency[A][i] == B) {
        flag = false;
        face = i;
    }
    i++;
} //next i

return(face); };

// -------------------------------------------------------------------------- //
int Class_VolTri::edge(
    int          T,
    ivector1D   &e
) {

// ========================================================================== //
// int Class_VolTri::edge(                                                    //
//     int          T,                                                        //
//     ivector1D   &e)                                                        //
//                                                                            //
// Return the local index of a specified edge on simplex A.                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T        : int, simplex global index                                     //
// - e        : ivector1D, global indices of vertices of edge                 //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - j        : int, edge local index                                         //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool            check = false;
int             index = -1;
ivector1D       f, d;

// Counters
int             i, j;
int             n, m;

// ========================================================================== //
// LOOP OVER EDGE                                                             //
// ========================================================================== //
f = e;
sort(f.begin(), f.end());
n = infos[e_type[T]].n_edges;
j = 0;
while (!check && (j < n)) {
    d = infos[e_type[T]].edges[j];
    m = d.size();
    for (i = 0; i < m; ++i) {
        d[i] = Simplex[T][d[i]];
    } //next i
    sort(d.begin(), d.end());
    check = (d == f);
    j++;
} //next j
if (check) {
    index = j-1;
}

return(index); };

// -------------------------------------------------------------------------- //
ivector1D Class_VolTri::EdgeNeigh(
    int          T, 
    int          i
) {

// ========================================================================== //
// ivector1D Class_VolTri::EdgeNeigh(                                         //
//     int          T,                                                        //
//     int          i)                                                        //
//                                                                            //
// Find neighboring simplicies of simplex T along the i-th edge.              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// - i      : int, edge local index                                           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - list   : ivector1D, list of edge neighbors                               //
// ========================================================================== //

// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //

// Local variables
int                             n_faces = infos[e_type[T]].n_faces;
bitpit::LIFOStack<int>                  stack(4*n_faces), visited(4*n_faces);
ivector1D                       e;
ivector1D                       neigh;
ivector1D::iterator             it_, end_;

// Counters
int                             n;
int                             j;
int                             S, A;

// ========================================================================== //
// BUILD ADJACENCIES IF NOT ALREADY BUILT                                     //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// ITERATIVELY ADD SIMPLICIES TO STACK                                        //
// ========================================================================== //

// Initialize stack --------------------------------------------------------- //
stack.push(T);
n = infos[e_type[T]].edges[i].size();
e.resize(n, -1);
for (j = 0; j < n; ++j) {
    e[j] = Simplex[T][infos[e_type[T]].edges[i][j]];
} //next i

// Iterate until stack is empty --------------------------------------------- //

// Initialize visited list
for (i = 0; i < visited.STACK.size(); ++i) {
    visited.STACK[i] = -1;
} //next i
visited.push(T);

// Loop over items in stack
while (stack.TOPSTK > 0) {

    // Pop item from stack
    S = stack.pop();

    // Loop over S-neighbors
    n = infos[e_type[S]].n_faces;
    for (j = 0; j < n; ++j) {
        A = Adjacency[S][j];
        if (A >= 0) {
            end_ = visited.STACK.begin() + visited.TOPSTK;
            if (find(visited.STACK.begin(), end_, A) == end_) {
                if (edge(A, e) >= 0) {
                    stack.push(A);
                    visited.push(A);
                    neigh.push_back(A);
                }
            }
        }
    } //next j
} //next item

// Remove face adjacent simplicies ------------------------------------------ //
for (j = 0; j < n_faces; ++j) {
    A = Adjacency[T][j];
    if (A >= 0) {
        it_ = find(neigh.begin(), neigh.end(), A);
        if (it_ != neigh.end()) {
            neigh.erase(it_);
        }
    }
} //next j

return(neigh); };

// -------------------------------------------------------------------------- //
ivector1D Class_VolTri::VertNeigh(
    int          T,
    int          i
) {

// ========================================================================== //
// ivector1D Class_VolTri::VertNeigh(                                         //
//     int          T,                                                        //
//     int          i)                                                        //
//                                                                            //
// Returns the indices of neighbors at the i-th vertex of simplex T.          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// - i      : int, vertex local index                                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - neigh  : ivector1D, vertex neighbors                                     //
// ========================================================================== //

// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //

// Local variables
int                             n_faces = infos[e_type[T]].n_faces;
bitpit::LIFOStack<int>                  stack(4*n_faces), visited(4*n_faces);
ivector1D                       e;
ivector1D                       neigh;
ivector1D                       e_neigh;
ivector1D::iterator             it_, end_;

// Counters
int                             n, m;
int                             j, k;
int                             S, A, V;

// ========================================================================== //
// BUILD ADJACENCIES IF NOT ALREADY BUILT                                     //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// ITERATIVELY ADD SIMPLICIES TO STACK                                        //
// ========================================================================== //

// Initialize stack --------------------------------------------------------- //
stack.push(T);
visited.push(T);
V = Simplex[T][i];

// Iterate until stack is empty --------------------------------------------- //
while (stack.TOPSTK > 0) {

    // Pop item from stack
    S = stack.pop();

    // Loop over S-neighbors
    n = infos[e_type[S]].n_faces;
    for (j = 0; j < n; ++j) {
        A = Adjacency[S][j];
        if (A >= 0) {
            end_ = visited.STACK.begin() + visited.TOPSTK;
            if (find(visited.STACK.begin(), end_, A) == end_) {
                if (vertex(A, V) >= 0) {
                    stack.push(A);
                    visited.push(A);
                    neigh.push_back(A);
                }
            }
        }
    } //next j
} //next item

// Remove face adjacent simplicies ------------------------------------------ //
for (j = 0; j < n_faces; ++j) {
    A = Adjacency[T][j];
    if (A >= 0) {
        it_ = find(neigh.begin(), neigh.end(), A);
        if (it_ != neigh.end()) {
            neigh.erase(it_);
        }
    }
} //next j

// Remove edge-adjacent simplicies ------------------------------------------ //

// find edges incident to vertex (T, i)
m = infos[e_type[T]].n_edges;
for (k = 0; k < m; k++) {
    if (find(infos[e_type[T]].edges[k].begin(),
             infos[e_type[T]].edges[k].end(),
             i) != infos[e_type[T]].edges[k].end()) {
        e.push_back(k);
    }
} //next k

// Remove edge-adjacent simplicies from neigh
m = e.size();
for (k = 0; k < m; ++k) {
    e_neigh = EdgeNeigh(T, e[k]);
    n = e_neigh.size();
    for (j = 0; j < n; ++j) {
        it_ = find(neigh.begin(), neigh.end(), e_neigh[j]);
        if (it_ != neigh.end()) {
            neigh.erase(it_);
        }
    } //next j
} //next k

return(neigh); };

// Find tools =============================================================== //

        // Find triangle --------------------------------------------------------------------- //
        int Class_VolTri::ReturnTriangleID(a3vector1D    &P) {

        // =================================================================================== //
        // int Class_VolTri::ReturnTriangleID(a3vector1D    &P)                                //
        //                                                                                     //
        // Find the index of the triangle enclosing a given point xp (Lawson's search          //
        // algorithm, for 2D volume triangulation only).                                       //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - P    : a3vector1D, with point coordinates                                         //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - T     : int, global index of triangle enclosing point P                           //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int        V1, V2;
        int        T;
        double     check;

        // Counters
        int        i;

        // =================================================================================== //
        // FIND TRIANGLE USING THE LAWSON'S SEARCH ALGORITHM                                   //
        // =================================================================================== //

        // First guess ----------------------------------------------------------------------- //
        T = nSimplex - 1;

        // Find triangle enclosing point P --------------------------------------------------- //
label_10:
        if (T < 0) {
            return(-1);
        }
        for (i = 0; i < Simplex[T].size(); i++) {
            V1 = Simplex[T][i];
            V2 = Simplex[T][(i+1) % 3];
            check = (Vertex[V1][1] - P[1])*(Vertex[V2][0] - P[0])
                  - (Vertex[V1][0] - P[0])*(Vertex[V2][1] - P[1]);
            if (check > 0.0) {
                T = Adjacency[T][i];
                goto label_10;
            }
        } //next i

        return (T); };

// -------------------------------------------------------------------------- //
int Class_VolTri::ReturnSimplexID(
    a3vector1D  &P,
    int          S0
) {

// ========================================================================== //
// int Class_VolTri::ReturnSimplexID(                                         //
//     a3vector1D  &P,                                                        //
//     int          S0)                                                       //
//                                                                            //
// Return the ID of the simplex enclosing point P.                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P     : a3vector1D, point coordinates                                    //
// - S0    : int, seed for search procedure                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - S     : int, ID of simplex enclosing point P                             //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                        check = false;
int                         n, m, p;
double                      max_dp, dp;
bvector1D                   visited(nSimplex, false);
ivector2D                   face_vlist;
a3vector1D                  x, y, xF, xS, dir;
a3vector2D                  face_normals;
bitpit::LIFOStack<int>              stack;

// Counters
int                         i, j, k, l, S, A;
/*debug*/int                         n_visited = 0;

// /*debug*/ofstream       log_file;
// /*debug*/log_file.open("SEARCH.log", ifstream::app);

// ========================================================================== //
// INITIALIZE PARAMETERS                                                      //
// ========================================================================== //
// x.fill(0.0);
// y.fill(0.0);

// /*debug*/log_file << "point coordinates: " << P << endl;

// ========================================================================== //
// BUILD ADJACENCY IF NOT ALREADY BUILT                                       //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// LOOP UNTIL SIMPLEX IS FOUND                                                //
// ========================================================================== //
stack.push(S0);
while ((stack.TOPSTK > 0) && (!check)) {

// /*debug*/n_visited++;

    // Pop item from stack -------------------------------------------------- //
    S = stack.pop();
    visited[S] = true;

// /*debug*/log_file << "  traversing simplex: " << S << endl;
// /*debug*/log_file << "  (visited: " << n_visited << " of " << nSimplex << endl;

    // Simplex infos -------------------------------------------------------- //
    n = infos[e_type[S]].n_vert;
    m = infos[e_type[S]].n_faces;

    // Get face vertices ---------------------------------------------------- //
    {
        face_vlist.resize(m);
        for (j = 0; j < m; ++j) {
            face_vlist[j] = FaceVertices(S, j);
        } //next j
    }

    // Simplex baricenter --------------------------------------------------- //
    {
        xS.fill(0.0);
        for (j = 0; j < n; ++j) {
            for (l = 0; l < 3; ++l) {
                xS[l] += Vertex[Simplex[S][j]][l];
            } //next l
        } //next j
        xS = xS/((double) n);
    }

    // Compute face normals ------------------------------------------------- //
    {
        face_normals.resize(m);
        for (j = 0; j < m; ++j) {
            if (face_vlist[j].size() == 2) {
                x = Vertex[face_vlist[j][1]] - Vertex[face_vlist[j][0]];
                x[2] = 0.0;
                y[0] = 0.0;
                y[1] = 0.0;
                y[2] = 1.0;
                face_normals[j] = crossProduct(x, y);
            }
            else {
                x = Vertex[face_vlist[j][2]] - Vertex[face_vlist[j][1]];
                y = Vertex[face_vlist[j][1]] - Vertex[face_vlist[j][0]];
                face_normals[j] = crossProduct(x, y);
            }
            face_normals[j] = face_normals[j]/max(norm2(face_normals[j]), 2.0e-16);
        } //next j
    }
    
    // Check if Simplex S encloses point P ---------------------------------- //
    {
        check = true;
        for (j = 0; j < m; ++j) {
            n = face_vlist[j].size();
    
            // Compute face center
            xF.fill(0.0);
            for (k = 0; k < n; ++k) {
                for (l = 0; l < 3; l++) {
                    xF[l] += Vertex[face_vlist[j][k]][l];
                } //next l
            } //next k
            xF = xF/((double) n);
    
            // Check if simplex encloses point
            dir = P - xF;
            dir = dir/max(norm2(dir), 2.0e-16);
            check = (check && (dotProduct(dir, face_normals[j]) <= 0.0));
        } //next j
    }

// /*debug*/log_file << "    encloses point: " << check << endl;

    // Look for best direction (Euristic search of best path) --------------- //
    // NOTES:                                                                 //
    // - modificare il criterio euristico per la scelta del path migliore     //
    //   * congiungente P-xS attraversa una faccia.                           //
    //   * distanza minima dalle facce.                                       //
    // ---------------------------------------------------------------------- //
    if (!check) {

        // local direction of searching path
        dir = P - xS;
        dir = dir/max(norm2(dir), 2.0e-16);
    
        // Loop over simplex faces
        max_dp = -2.0;
        i = 0;
        j = -1;
        while (i < m) {
            dp = dotProduct(face_normals[i], dir);
            if (dp > max_dp) {
                j = i;
                max_dp = dp;
            }
            i++;
        } //next i
    
        // Alternative directions
        for (i = 0; i < m; ++i) {
            A = Adjacency[S][i];
            if ((i != j) && (A >= 0) && (!visited[A])) {
                stack.push(A);
            }
        } //next i
        A = Adjacency[S][j];
        if ((A >= 0) && (!visited[A])) {
            stack.push(A);
        }
    }

} //next simplex

// /*debug*/log_file << "  (visited: " << n_visited << " of " << nSimplex << ")" << endl;

// Old algorithm ============================================================ //
{
    // while ((n_visited <= nSimplex) && (!check) && (S0 >= 0)) {
    
        // // Update simplex ID ---------------------------------------------------- //
        // // cout << "    on simplex " << S0;
        // S = S0;
        // n_visited++;
        // visited[S] = true;
        // /*debug*/log_file << "  traversing simplex: " << S << endl;
        // /*debug*/log_file << "  (visited: " << n_visited << " of " << nSimplex << endl;
    
        // // Simplex infos -------------------------------------------------------- //
        // m = infos[e_type[S]].n_faces;
        // p = infos[e_type[S]].n_vert;
    
        // // Compute simplex baricenter ------------------------------------------- //
        // {
            // xS.fill(0.0);
            // for (j = 0; j < p; ++j) {
                // for (l = 0; l < dim; ++l) {
                    // xS[l] += Vertex[Simplex[S][j]][l];
                // } //next l
            // } //next j
            // xS = xS/((double) p);
        // }
    
        // // Get face vertices ---------------------------------------------------- //
        // {
            // face_vlist.resize(m);
            // for (j = 0; j < m; ++j) {
                // face_vlist[j] = FaceVertices(S, j);
            // } //next j
        // }
    
        // // Compute face normals ------------------------------------------------- //
        // {
            // face_normals.resize(m);
            // for (j = 0; j < m; ++j) {
                // if (face_vlist[j].size() == 2) {
                    // x[0] = Vertex[face_vlist[j][1]][0] - Vertex[face_vlist[j][0]][0];
                    // x[1] = Vertex[face_vlist[j][1]][1] - Vertex[face_vlist[j][0]][1];
                    // x[2] = 0.0;
                    // y[0] = 0.0;
                    // y[1] = 0.0;
                    // y[2] = 1.0;
                    // face_normals[j] = crossProduct(x, y);
                // }
                // else {
                    // x[0] = Vertex[face_vlist[j][2]][0] - Vertex[face_vlist[j][1]][0];
                    // x[1] = Vertex[face_vlist[j][2]][1] - Vertex[face_vlist[j][1]][1];
                    // x[2] = Vertex[face_vlist[j][2]][2] - Vertex[face_vlist[j][1]][2];
                    // y[0] = Vertex[face_vlist[j][1]][0] - Vertex[face_vlist[j][0]][0];
                    // y[1] = Vertex[face_vlist[j][1]][1] - Vertex[face_vlist[j][0]][1];
                    // y[2] = Vertex[face_vlist[j][1]][2] - Vertex[face_vlist[j][0]][2];
                    // face_normals[j] = crossProduct(x, y);
                // }
                // face_normals[j] = face_normals[j]/max(norm2(face_normals[j]), 2.0e-16);
            // } //next j
        // }
    
        // // Check if P is enclosed in simplex S0 --------------------------------- //
        // {
            // check = true;
            // for (j = 0; j < m; ++j) {
                // n = face_vlist[j].size();
        
                // // Compute face center
                // xF.fill(0.0);
                // for (k = 0; k < n; ++k) {
                    // for (l = 0; l < dim; l++) {
                        // xF[l] += Vertex[face_vlist[j][k]][l];
                    // } //next l
                // } //next k
                // xF = xF/((double) n);
        
                // // Check if simplex encloses point
                // for (l = 0; l < dim; l++) {
                    // dir[l] = P[l] - xF[l];
                // } //next l
                // dir = dir/max(norm2(dir), 2.0e-16);
                // check = (check && (dotProduct(dir, face_normals[j]) <= 0.0));
                // // cout << " f " << j << ", c: " << check;
            // } //next j
        // }
        // // cout << ", check: " << check << endl;
        // /*debug*/log_file << "    encloses point: " << check << endl;
    
        // // Look for best direction (Euristic search) ---------------------------- //
        // if (!check) {
        
            // // Find face to be crossed (Euristic search of best path)
            // {
        
                // // path local direction
                // for (l = 0; l < dim; ++l) {
                    // dir[l] = P[l] - xS[l];
                // } //next l
                // dir = dir/max(norm2(dir), 2.0e-16);
        
                // // Loop over simplex faces
                // // WARNING: infinite loop at corner simplicies!!
                // max_dp = -2.0;
                // i = -1;
                // j = 0;
                // while (j < m) {
                    // dp = dotProduct(face_normals[j], dir);
                    // A = Adjacency[S][j];
                    // if ((dp > max_dp) && (A >= 0)) {//&& (!visited[A])) {
                        // i = j;
                        // max_dp = dp;
                    // }
                    // j++;
                // } //next j
            // }
    
            // // Move to adjacent simplex
            // if (i >= 0) { S0 = Adjacency[S][i]; }
            // else        { S0 = -1;}
        // }
    
        // /*debug*/log_file << "    next simplex: " << S0 << endl;
    // } //next simplex
}
// /*debug*/log_file.close();

return(S); };

// Bounding boxes =========================================================== //

// -------------------------------------------------------------------------- //
void Class_VolTri::BoundingBox(
    array<double, 2> &x_ext,
    array<double, 2> &y_ext,
    array<double, 2> &z_ext
) {

// ========================================================================== //
// void Class_VolTri::BoundingBox(                                            //
//     array<double, 2> &x_ext,                                               //
//     array<double, 2> &y_ext,                                               //
//     array<double, 2> &z_ext)                                               //
//                                                                            //
// Compute limits of mesh bounding box (3D case).                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - x_ext   : array<double, 2>, extent of bounding box in the x direction    //
// - y_ext   : array<double, 2>, extent of bounding box in the y direction    //
// - z_ext   : array<double, 2>, extent of bounding box in the z direction    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         n;

// Counters
int         i, j, I;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Compute bounding box limits ---------------------------------------------- //
x_ext[0] = x_ext[1] = Vertex[Simplex[0][0]][0];
y_ext[0] = y_ext[1] = Vertex[Simplex[0][0]][1];
z_ext[0] = z_ext[1] = Vertex[Simplex[0][0]][2];
for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I = Simplex[i][j];
        x_ext[0] = min(x_ext[0], Vertex[I][0]);
        x_ext[1] = max(x_ext[1], Vertex[I][0]);
        y_ext[0] = min(y_ext[0], Vertex[I][1]);
        y_ext[1] = max(y_ext[1], Vertex[I][1]);
        z_ext[0] = min(z_ext[0], Vertex[I][2]);
        z_ext[1] = max(z_ext[1], Vertex[I][2]);
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::BoundingBox(
    a3vector2D         &V,
    array<double, 2>   &x_ext,
    array<double, 2>   &y_ext,
    array<double, 2>   &z_ext
) {

// ========================================================================== //
// void Class_VolTri::BoundingBox(                                            //
//     a3vector2D         &V,                                                 //
//     array<double, 2>   &x_ext,                                             //
//     array<double, 2>   &y_ext,                                             //
//     array<double, 2>   &z_ext)                                             //
//                                                                            //
// Compute limits of tasselation bounding box (3D case). Vertex coordinate    //
// list is provided externally).                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V       : a3vector2D vertex coordinate list. V[i][0], V[i][1], ... are   //
//             the x, y, ... coordinates of the i-th vertex.                  //
// - x_ext   : array<double, 2>, extent of bounding box in the x direction    //
// - y_ext   : array<double, 2>, extent of bounding box in the y direction    //
// - z_ext   : array<double, 2>, extent of bounding box in the z direction    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         n;

// Counters
int         i, j, I;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Compute bounding box limits ---------------------------------------------- //
x_ext[0] = x_ext[1] = V[Simplex[0][0]][0];
y_ext[0] = y_ext[1] = V[Simplex[0][0]][1];
z_ext[0] = z_ext[1] = V[Simplex[0][0]][2];
for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I = Simplex[i][j];
        x_ext[0] = min(x_ext[0], V[I][0]);
        x_ext[1] = max(x_ext[1], V[I][0]);
        y_ext[0] = min(y_ext[0], V[I][1]);
        y_ext[1] = max(y_ext[1], V[I][1]);
        z_ext[0] = min(z_ext[0], V[I][2]);
        z_ext[1] = max(z_ext[1], V[I][2]);
    } //next j
} //next i

return; };

// Extract boundaries ======================================================= //

// -------------------------------------------------------------------------- //
void Class_VolTri::ExtractBoundaries(
    Class_SurfTri    &Bounds,
    int               BC_flag
) {

// ========================================================================== //
// void Class_VolTri::ExtractBoundaries(                                      //
//     Class_SurfTri    &Bounds,                                              //
//     int               BC_flag)                                             //
//                                                                            //
// Extract boundaries of volume mesh.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Bounds     : Class_SurfTri, surface mesh                                 //
// - BC_flag    : int, boundary condition flag                                //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
ivector1D                simplex;

// Counters
int                      nS;
int                      i, j;
int                      n, m;
int                      T;

// ========================================================================== //
// BUILD ADJACENCY IF NOT ALREADY BUILT                                       //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// RESIZE DATA STRUCTURES                                                     //
// ========================================================================== //
Bounds.AddVertices(Vertex);
nS = Bounds.nSimplex;
Bounds.nSimplex += CountFreeFaces();
Bounds.ResizeSimplex();
Bounds.nSimplex = nS;

// ========================================================================== //
// LOOP OVER SIMPLICIES                                                       //
// ========================================================================== //
for (T = 0; T < nSimplex; ++T) {
    n = infos[e_type[T]].n_faces;
    for (i = 0; i < n; ++i) {
        if (Adjacency[T][i] == BC_flag) {
            m = infos[e_type[T]].faces[i].size();
            simplex.resize(m, -1);
            for (j = 0; j < m; ++j) {
                simplex[j] = Simplex[T][infos[e_type[T]].faces[i][j]];
            } //next j
            Bounds.AddSimplex(simplex);
        }
    } //next i
} //next T

// ========================================================================== //
// CLEAN SURFACE TASSELATION                                                  //
// ========================================================================== //
Bounds.RemoveIsolatedVertex();
Bounds.ResizeVertex();

return; };







