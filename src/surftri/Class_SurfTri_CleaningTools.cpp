// ========================================================================== //
//                         - Class_SurfTri -                                  //
//                                                                            //
// Grid manager for unstructured meshes.                                      //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.1                                                            //
//                                                                            //
// All rights reserved.                                                       //
//									      //	
// notes: true double simplex checks included                                 //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include "Class_SurfTri.hpp"

// ========================================================================== //
// SORTING METHODS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::BinSortV(
    ivector1D   &index,
    int          n
) {

// ========================================================================== //
// void Class_SurfTri::BinSortV(                                              //
//     ivector1D   &index,                                                    //
//     int          n)                                                        //
//                                                                            //
// Sort tassellation vertices on regular bins.                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - index   : ivector2D, map vertex->bins                                    //
// - n       : int (optional), number of bins used for sorting.               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
//ht int           dim = Vertex[0].size();
dvector1D     dx(3, 0.0);
ivector1D     ix(3, -1);
dvector2D     xlim(3, dvector1D(2, 0.0));

// Counters
int           i, j, k, I_, m;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
index.resize(nVertex);

// ========================================================================== //
// ENCLOSE TASSELATION IN A SQUARE BOX                                        //
// ========================================================================== //

// Mesh limits
for (j = 0; j < 3; j++) {
    xlim[j][0] = Vertex[0][j];
    xlim[j][1] = Vertex[0][j];
    for (i = 0; i < nVertex; i++) {
        xlim[j][0] = min(xlim[j][0], Vertex[i][j]);
        xlim[j][1] = max(xlim[j][1], Vertex[i][j]);
    }
    dx[j] = xlim[j][1] - xlim[j][0];
    dx[j] = max(dx[j], 1.0);
    xlim[j][0] = xlim[j][0] - 0.05 * dx[j];
    xlim[j][1] = xlim[j][1] + 0.05 * dx[j];
} //next j

// Mesh spacing
for (j = 0; j < 3; j++) {
    dx[j] = (xlim[j][1] - xlim[j][0])/((double) n);
} //next j

// ========================================================================== //
// ASSOCIATE EACH VERTEX TO A CELL                                            //
// ========================================================================== //

// Loop over vertexes
//ht if (dim == 2) {
//ht     for (i = 0; i < nSimplex; i++) {
//ht         m = Simplex[i].size();
//ht         for (j = 0; j < m; j++) {
//ht             for (k = 0; k < dim; k++) {
//ht                 ix[k] = (int) floor((Vertex[Simplex[i][j]][k] - xlim[k][0])/dx[k]);
//ht             } //next k
//ht             I_ = n * ix[0] + ix[1];
//ht             index[Simplex[i][j]] = I_;
//ht         } //next j
//ht     } //next i
//ht }
//ht else if (dim == 3) {
//ht     for (i = 0; i < nSimplex; i++) {
//ht         m = Simplex[i].size();
//ht         for (j = 0; j < m; j++) {
//ht             for (k = 0; k < dim; k++) {
//ht                 ix[k] = (int) floor((Vertex[Simplex[i][j]][k] - xlim[k][0])/dx[k]);
//ht             } //next k
//ht             I_ = n * n * ix[0] + n * ix[1] + ix[2];
//ht             index[Simplex[i][j]] = I_;
//ht         } //next j
//ht     } //next i
//ht }

for (i = 0; i < nSimplex; i++) {
    m = Simplex[i].size();
    for (j = 0; j < m; j++) {
        for (k = 0; k < 3; k++) {
            ix[k] = (int) floor((Vertex[Simplex[i][j]][k] - xlim[k][0])/dx[k]);
        } //next k
        I_ = n * n * ix[0] + n * ix[1] + ix[2];
        index[Simplex[i][j]] = I_;
    } //next j
} //next i


return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BinSortV(
    dvecarr3E   &X,
    ivector1D   &index,
    int          n
) {

// ========================================================================== //
// void Class_SurfTri::BinSortV(                                              //
//     dvecarr3E   &X,                                                        //
//     ivector1D   &index,                                                    //
//     int          n)                                                        //
//                                                                            //
// Sort tassellation vertices on regular bins using an external vertex list   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - X       : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are   //
//             the x, y, ... coordinates of the i-th vertex.                  //
// - index   : ivector2D, map vertex->bins                                    //
// - n       : int (optional), number of bins used for sorting.               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int           nV = X.size();
//ht int           dim = X[0].size();
dvector1D     dx(3, 0.0);
ivector1D     ix(3, -1);
dvector2D     xlim(3, dvector1D(2, 0.0));

// Counters
int           i, j, k, I_, m;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
index.resize(nV);

// ========================================================================== //
// ENCLOSE TASSELATION IN A SQUARE BOX                                        //
// ========================================================================== //

// Mesh limits
for (j = 0; j < 3; j++) {
    xlim[j][0] = X[0][j];
    xlim[j][1] = X[0][j];
    for (i = 0; i < nV; i++) {
        xlim[j][0] = min(xlim[j][0], X[i][j]);
        xlim[j][1] = max(xlim[j][1], X[i][j]);
    }
    dx[j] = xlim[j][1] - xlim[j][0];
    dx[j] = max(dx[j], 1.0);
    xlim[j][0] = xlim[j][0] - 0.05 * dx[j];
    xlim[j][1] = xlim[j][1] + 0.05 * dx[j];
} //next j

// Mesh spacing
for (j = 0; j < 3; j++) {
    dx[j] = (xlim[j][1] - xlim[j][0])/((double) n);
} //next j

// ========================================================================== //
// ASSOCIATE EACH VERTEX TO A CELL                                            //
// ========================================================================== //

// Loop over vertexes
//ht if (dim == 2) {
//ht     for (i = 0; i < nSimplex; i++) {
//ht         m = Simplex[i].size();
//ht         for (j = 0; j < m; j++) {
//ht             for (k = 0; k < dim; k++) {
//ht                 ix[k] = (int) floor((X[Simplex[i][j]][k] - xlim[k][0])/dx[k]);
//ht             } //next k
//ht             I_ = n * ix[0] + ix[1];
//ht             index[Simplex[i][j]] = I_;
//ht         } //next j
//ht     } //next i
//ht }
//ht else if (dim == 3) {
//ht     for (i = 0; i < nSimplex; i++) {
//ht         m = Simplex[i].size();
//ht         for (j = 0; j < m; j++) {
//ht             for (k = 0; k < dim; k++) {
//ht                 ix[k] = (int) floor((X[Simplex[i][j]][k] - xlim[k][0])/dx[k]);
//ht             } //next k
//ht             I_ = n * n * ix[0] + n * ix[1] + ix[2];
//ht             index[Simplex[i][j]] = I_;
//ht         } //next j
//ht     } //next i
//ht }

for (i = 0; i < nSimplex; i++) {
    m = Simplex[i].size();
    for (j = 0; j < m; j++) {
        for (k = 0; k < 3; k++) {
            ix[k] = (int) floor((X[Simplex[i][j]][k] - xlim[k][0])/dx[k]);
        } //next k
        I_ = n * n * ix[0] + n * ix[1] + ix[2];
        index[Simplex[i][j]] = I_;
    } //next j
} //next i


return; };

// ========================================================================== //
// TASSELATION PARAMETERS                                                     //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::SetTolerance(
    void
) {

// ========================================================================== //
// void Class_SurfTri::SetTolerance(                                          //
//     void)                                                                  //
//                                                                            //
// Set tolerance for distance check, (based on minimal edge length).          //
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
double     length;

// Counters
int        i, j, m, I_;

// ========================================================================== //
// SET TOLLERANCE                                                             //
// ========================================================================== //

// Initialize value
toll =  1.0e-1*norm2(Vertex[Simplex[0][1]] - Vertex[Simplex[0][0]]);

// Loop over simplicies
for (I_ = 0; I_ < nSimplex; I_++) {
    m = Simplex[I_].size();
    for (i = 0; i < m; i++) {
        j = ((i + 1) % Simplex[I_].size());
        length = norm2(Vertex[Simplex[I_][i]] - Vertex[Simplex[I_][j]]);
        toll = min(toll, 1.0e-1 * length);
    } // next i
} // next I_

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::SetTolerance(
    dvecarr3E   &X
) {

// ========================================================================== //
// void Class_SurfTri::SetTolerance(                                          //
//     dvecarr3E   &X)                                                        //
//                                                                            //
// Set tolerance for distance check, (based on minimal edge length). Vertex   //
// coordinate list is provided externally.                                    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - X   : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are the   //
//         x, y, ... coordinates of the i-th vertex.                          //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double     length;

// Counters
int        i, j, I_, m;

// ========================================================================== //
// SET TOLLERANCE                                                             //
// ========================================================================== //

// Initialize value
toll =  1.0e-1*norm2(X[Simplex[0][1]] - X[Simplex[0][0]]);

// Loop over simplicies
for (I_ = 0; I_ < nSimplex; I_++) {
    m = Simplex[I_].size();
    for (i = 0; i < m; i++) {
        j = ((i + 1) % m);
        length = norm2(X[Simplex[I_][i]] - X[Simplex[I_][j]]);
        toll = min(toll, 1.0e-1 * length);
    } // next i
} //next I_

return; };

// ========================================================================== //
// FIXING ALGORITHMS                                                          //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::FixNodeNumb(
    void
) {

// ========================================================================== //
// void Class_SurfTri::FixNodeNumb(                                           //
//     void)                                                                  //
//                                                                            //
// Change node numbering from clocwise <-> counter-clockwise direction on     //
// each simplex to match normal direction (right-hand-rule).                  //
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
bool                flaga;
//ht int                 dim = Vertex[0].size();
int                 idummy;
darray3E            x, y, z;
ivector1D           idummy1D;
ivector2D           idummy2D;

// Counters
int                 i, j, I_, m;

// ========================================================================== //
// FLAG                                                                       //
// ========================================================================== //
flaga = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

// ========================================================================== //
// FIX NODE NUMBERING                                                         //
// ========================================================================== //
for (I_ = 0; I_ < nSimplex; I_++){

    // Simplex type
    m = Simplex[I_].size();

    if (m < 2) {

        // Point ------------------------------------------------------------ //
        // no action taken

    }
    else if (m == 2) {

        // Segment ---------------------------------------------------------- //

        // normal direction based on vertex numbering
        z = Vertex[Simplex[I_][1]] - Vertex[Simplex[I_][0]];
        z = z/norm2(z);

        // Fix node numbering
        if (dotProduct(z, Normal[I_]) < 0.0) {
            invert_loc_num(I_);
        }
    }
    else {

        // Non-degenerate simplex ------------------------------------------- //

        // normal vector based on vertex numbering
        x = Vertex[Simplex[I_][1]] - Vertex[Simplex[I_][0]];
        y = Vertex[Simplex[I_][2]] - Vertex[Simplex[I_][1]];
        x = x/norm2(x);
        y = y/norm2(y);
        z = crossProduct(x, y);
        z = z/norm2(z);

        // Fix node numbering
        if (dotProduct(z, Normal[I_]) < 0.0) {
            invert_loc_num(I_);
        }
    }
} //next I_

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::FixNodeNumb(
    dvecarr3E   &X
) {

// ========================================================================== //
// void Class_SurfTri::FixNodeNumb(                                           //
//     dvecarr3E   &X)                                                        //
//                                                                            //
// Change node numbering from clocwise <-> counter-clockwise direction on     //
// each simplex to match normal direction (right-hand-rule). Vertex           //
// coordinate list is provided externally.                                    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1],  ... are the //
//          x, y, ... coordinates of the i-th vertex.                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                flaga;
//ht int                 dim = X[0].size();
int                 idummy;
darray3E            x, y, z;
ivector1D           idummy1D;
ivector2D           idummy2D;

// Counters
int                 i, j, I_, m;

// ========================================================================== //
// CHECK IF ADJACENCY ARE ALREADY BUILT                                       //
// ========================================================================== //
flaga = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

// ========================================================================== //
// FIX NODE NUMBERING                                                         //
// ========================================================================== //
for (I_ = 0; I_ < nSimplex; I_++){

    // Simplex type
    m = Simplex[I_].size();

    if (m < 2) {

        // Point ------------------------------------------------------------ //
        // no action taken

    }
    else if (m == 2) {

        // Segment ---------------------------------------------------------- //

        // normal direction based on vertex numbering
        z = X[Simplex[I_][1]] - X[Simplex[I_][0]];
        z = z/norm2(z);

        // Adjust node numbering
        if (dotProduct(z, Normal[I_]) < 0.0) {
            invert_loc_num(I_);
        }
    }
    else {

        // Non-degenerate simplex ------------------------------------------- //

        // normal vector based on vertex numbering
        x = X[Simplex[I_][1]] - X[Simplex[I_][0]];
        y = X[Simplex[I_][2]] - X[Simplex[I_][1]];
        x = x/norm2(x);
        y = y/norm2(y);
        z = crossProduct(x, y);
        z = z/norm2(z);

        // Fix node numbering
        if (dotProduct(z, Normal[I_]) < 0.0) {
            invert_loc_num(I_);
        }
    }
} //next I_

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::FixNodeNumb(
    int          seed
) {

// ========================================================================== //
// void Class_SurfTri::FixNodeNumb(                                           //
//     int          seed)                                                     //
//                                                                            //
// Fix local node numbering according to the local numbering of a given       //
// simplex. (manifold surface only)                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - seed    : int, simplex global index                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bvector1D             flag(nSimplex, false);
int                   m, n, p;
bitpit::LIFOStack<int>        stack(sqrt(nSimplex));

// Counters
int                   I_, i, j, e, J;

// ========================================================================== //
// BUILD ADJACENCY IF NOT ALREADY BUILT                                       //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// ADJUST NODE NUMBERING                                                      //
// ========================================================================== //
stack.push(seed);
while (stack.TOPSTK > 0) {

    // Pop item from stack
    I_ = stack.pop();
    m = Simplex[I_].size();
    flag[I_] = true;

    // Loop over neighbor
    for (i = 0; i < m; i++) {
        if (Adjacency[I_][i][0] >= 0) {
            n = Adjacency[I_][i].size();
            for (j = 0; j < n; j++) {
                J = Adjacency[I_][i][j];
                if (!flag[J]) {
                    p = Simplex[J].size();
                    e = edge(I_, J);
                    if ((Simplex[J][e] == Simplex[I_][i])
                     && (Simplex[J][(e+1) % p] == Simplex[I_][(i+1) % m])) {
                        invert_loc_num(J);
                    }
                    stack.push(J);
                }
           } //next j
        }
    } //next i
} //next item

return; }

// -------------------------------------------------------------------------- //
void Class_SurfTri::InvertOrder(
    void
) {

// ========================================================================== //
// void Class_SurfTri::InvertOrder(                                           //
//     void)                                                                  //
//                                                                            //
// Invert local number of simplicies.                                         //
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
bool            flag_a, flag_n;

// Counters
int             i, j;
int             n, m;
int             T;

// ========================================================================== //
// SET PARAMETERS                                                             //
// ========================================================================== //
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));
flag_n = ((Normal.size() > 0) && (Normal.size() >= nSimplex));

// ========================================================================== //
// INVERT LOCAL NUMBERING                                                     //
// ========================================================================== //

// Update simplex-vertex connectivity --------------------------------------- //
for (T = 0; T < nSimplex; ++T) {
    n = Simplex[T].size();
    m = n/2;
    for (i = 0; i < m; ++i) {
        j = n - i - 1;
        swap(Simplex[T][i], Simplex[T][j]);
    } //next i
} //next T

// Update adjacencies ------------------------------------------------------- //
if (flag_a) {
    for (T = 0; T < nSimplex; ++T) {
        n = Simplex[T].size();
        m = n/2;
        for (i = 0; i < m; ++i) {
            j = n - i - 1;
            swap(Adjacency[T][i], Adjacency[T][j]);
        } //next i
    } //next T
}

// Update normals ----------------------------------------------------------- //
if (flag_n) {
    for (T = 0; T < nSimplex; ++T) {
        Normal[T] = -1.0*Normal[T];
    } //next T
}

return; }


// ========================================================================== //
// OBJECT BUILDER                                                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::GenerateNormals(
    void
) {

// ========================================================================== //
// void Class_SurfTri::GenerateNormals(                                       //
//     void)                                                                  //
//                                                                            //
// Generate normals according to local numeration.                            //
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
darray3E         x, y, z, dummy;

// Counter
int              I_, m;

// ========================================================================== //
// GENERATE NORMALS                                                           //
// ========================================================================== //

// Reshape Normal list ------------------------------------------------------ //
if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
    ResizeNormal();
}

// Update normal for each simplex ------------------------------------------- //
for (I_ = 0; I_ < nSimplex; I_++) {

    // Simplex type
    m = Simplex[I_].size();

    // Compute normal
    if (m < 2) {

        // Point ------------------------------------------------------------ //
        // no action taken

    }
    else if (m == 2) {

        // Segment ---------------------------------------------------------- //
        z = Vertex[Simplex[I_][1]] - Vertex[Simplex[I_][0]];

        if( dimensions == 2 ){
            dummy = z ;

            z[0] = -dummy[1];
            z[1] =  dummy[0];
        };

    }
    else if (m == 3) {

        // Triangle --------------------------------------------------------- //
        x = Vertex[Simplex[I_][1]] - Vertex[Simplex[I_][0]];
        y = Vertex[Simplex[I_][2]] - Vertex[Simplex[I_][0]];
        x = x/norm2(x);
        y = y/norm2(y);
        z = crossProduct(x, y);

    }
    else {

        // Non-degenerate simplex ------------------------------------------- //

        // Scope variables
        int                 i, j, k;
        darray3E            t;

        // Reset normal
        z.fill(0.) ;

        // Compute average normal
        for (i = 0; i < m; i++) {
            j = (i+1) % m;
            k = (j+1) % m;
            x = Vertex[Simplex[I_][j]] - Vertex[Simplex[I_][i]];
            y = Vertex[Simplex[I_][k]] - Vertex[Simplex[I_][i]];
            x = x/norm2(x);
            y = y/norm2(y);
            t = crossProduct(x, y);
            t = t/norm2(t);
            z = z + t/((double) m);
        } //next i
    }

    // Set normal
    Normal[I_] = z/norm2(z);

} //next I_

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::GenerateNormals(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // void Class_SurfTri::GenerateNormals(                                       //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Generate normals according to vertex local numeration. Vertex coordinate   //
    // list is provided by an external list.                                      //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X     : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are the //
    //           x, y, ... coordinates of the i-th vertex                         //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    darray3E         x, y, z, dummy;

    // Counter
    int              I_, m;

    // ========================================================================== //
    // GENERATE NORMALS                                                           //
    // ========================================================================== //

    // Reshape normal list ------------------------------------------------------ //
    if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
        ResizeNormal();
    }

    // Loop over simplicies ----------------------------------------------------- //
    for (I_ = 0; I_ < nSimplex; I_++) {

        // Simplex type
        m = Simplex[I_].size();

        // Compute normal
        if (m < 2) {

            // Point ------------------------------------------------------------ //
            // no action taken

        }
        else if (m == 2) {

            // Segment ---------------------------------------------------------- //
            z = X[Simplex[I_][1]] - X[Simplex[I_][0]];

            if( dimensions == 2 ){
                dummy = z ;

                z[0] = -dummy[1];
                z[1] =  dummy[0];
            };
        }
        else if (m == 3) {

            // Triangle --------------------------------------------------------- //
            x = X[Simplex[I_][1]] - X[Simplex[I_][0]];
            y = X[Simplex[I_][2]] - X[Simplex[I_][1]];
            x = x/norm2(x);
            y = y/norm2(y);
            z = crossProduct(x, y);

        }
        else {

            // Non-degenerate simplex ------------------------------------------- //

            // Scope variables
            int                 i, j, k;
            darray3E            t;

            // Reset normal
            z = t;

            // Compute average normal
            for (i = 0; i < m; i++) {
                j = (i+1) % m;
                k = (j+1) % m;
                x = X[Simplex[I_][j]] - X[Simplex[I_][i]];
                y = X[Simplex[I_][k]] - X[Simplex[I_][i]];
                x = x/norm2(x);
                y = y/norm2(y);
                t = crossProduct(x, y);
                t = t/norm2(t);
                z = z + t/((double) m);
            } //next i
        }

        // Set normal
        Normal[I_] = z/norm2(z);

    }

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::GenerateENormals(
        void
        ) {

    // ========================================================================== //
    // void Class_SurfTri::GenerateENormals(
    //     ivector2D   &Edges,                                                    //
    //     ivector2D   &EdgeAdj,                                                  //
    //     dvecarr3E   &Enormals)                                                 //
    //                                                                            //
    // Compute edge normals.                                                      //
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
    darray3E    tmp;

    // Counters
    int         T, i, j, m, nE;

    // ========================================================================== //
    // GENERATE EDGES                                                             //
    // ========================================================================== //
    if ((Edge.size() == 0)
            || (Simplex2Edge.size() == 0)
            || (Edge.size() != Simplex2Edge.size())) {
        BuildEdges();
    }

    // ========================================================================== //
    // GENERATE NORMALS                                                           //
    // ========================================================================== //
    if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
        GenerateNormals();
    }

    // ========================================================================== //
    // RESIZE INPUT VARIABLES                                                     //
    // ========================================================================== //
    nE = Edge.size();
    tmp.fill(0.);
    ENormal.resize(nE, tmp);

    // ========================================================================== //
    // GENERATE EDGE NORMALS                                                      //
    // ========================================================================== //

    // Compute normals
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            ENormal[Simplex2Edge[T][i]] = ENormal[Simplex2Edge[T][i]] + Normal[T];
        } //next i
    } //next T

    // Normalization
    for (T = 0; T < nE; T++) {
        ENormal[T] = ENormal[T]/norm2(ENormal[T]);
    }

    return; }

    // -------------------------------------------------------------------------- //
    void Class_SurfTri::GenerateENormals(
            dvecarr3E   &X
            ) {

        // ========================================================================== //
        // void Class_SurfTri::GenerateENormals(                                      //
        //     dvecarr3E   &X)                                                        //
        //                                                                            //
        // Compute edge normals. Vertex coordinate list is provided externally.       //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - X        : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ...      //
        //              the x, y, ... coordinates of the i-th vertex.                 //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        darray3E    tmp;

        // Counters
        int         T, i, j, m, nE;

        // ========================================================================== //
        // GENERATE EDGES                                                             //
        // ========================================================================== //
        if ((Edge.size() == 0)
                || (Simplex2Edge.size() == 0)
                || (Edge.size() != Simplex2Edge.size())) {
            BuildEdges();
        }

        // ========================================================================== //
        // GENERATE NORMALS                                                           //
        // ========================================================================== //
        if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
            GenerateNormals(X);
        }

        // ========================================================================== //
        // RESIZE INPUT VARIABLES                                                     //
        // ========================================================================== //
        nE = Edge.size();
        tmp.fill(0.0);
        ENormal.resize(nE, tmp);

        // ========================================================================== //
        // GENERATE EDGE NORMALS                                                      //
        // ========================================================================== //

        // Compute normals
        for (T = 0; T < nSimplex; T++) {
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                ENormal[Simplex2Edge[T][i]] = ENormal[Simplex2Edge[T][i]] + Normal[T] ;
            } //next i
        } //next T

        // Normalization
        for (T = 0; T < nE; T++) {
            ENormal[T] = ENormal[T]/norm2(ENormal[T]);
        }

        return; }

        // -------------------------------------------------------------------------- //
        void Class_SurfTri::GenerateVNormals(
                unsigned char flag
                ) {

            // ========================================================================== //
            // void Class_SurfTri::GenerateVNormals(                                      //
            //     unsigned char flag)                                                    //
            //                                                                            //
            // Generate vertex normals using algorithm #1 or #2 depending on the value    //
            // of the input flag:                                                         //
            // algorithm #1 (flag = 0): compute normals as the weighted average of        //
            //                          normals of incident simplicies.                   //
            // algorithm #2 (flag = 1): compute normals as the weighted average of        //
            //                          normals of incident edges.                        //
            // ========================================================================== //
            // INPUT                                                                      //
            // ========================================================================== //
            // - flag    : unsigend char (default 0), flag for algorithm selection        //
            // ========================================================================== //
            // OUTPUT                                                                     //
            // ========================================================================== //
            // - none                                                                     //
            // ========================================================================== //

            // ========================================================================== //
            // VARIABLES DECLARATION                                                      //
            // ========================================================================== //

            // Parameters
            double const    pi = 3.1415926535897932;

            // Local variables
            int             nE;
            darray3E        tmp;

            // Counters
            int             V, T, i, j, m;

            // ========================================================================== //
            // COMPUTE NORMALS IF NOT ALREADY COMPUTED                                    //
            // ========================================================================== //
            if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
                GenerateNormals();
            }

            // ========================================================================== //
            // COMPUTE EDGES NORMALS IF NOT ALREADY COMPUTED                              //
            // ========================================================================== //
            if ((ENormal.size() == 0) || (ENormal.size() < Edge.size())) {
                GenerateENormals();
            }

            // ========================================================================== //
            // RESIZE INPUT VARIABLES                                                     //
            // ========================================================================== //
            tmp.fill(0.) ;
            VNormal.resize(nVertex, tmp);

            // ========================================================================== //
            // COMPUTE VERTEX NORMALS                                                     //
            // ========================================================================== //

            // Compute vertex normals using algoritm #2 --------------------------------- //
            if (flag == 1) {

                // Scope variables
                // none

                // Compute vertex normals
                nE = Edge.size();
                for (T = 0; T < nE; T++) {
                    m = Edge[T].size();
                    for (i = 0; i < m; i++) {
                        V = Edge[T][i];
                        VNormal[V] = VNormal[V] + ENormal[T];
                    } //next i
                } //next T
            }

            // Compute vertex normals using algoritm #1 --------------------------------- //
            else if (flag == 0) {

                // Scope variables
                double          angle;

                // Compute vertex normals
                for (T = 0; T < nSimplex; ++T) {
                    m = Simplex[T].size();
                    for (i = 0; i < m; ++i) {
                        V = Simplex[T][i];
                        Angle(T, angle, i);
                        VNormal[V] = VNormal[V] + 0.5 * angle * Normal[T]/pi;
                    } //next i
                } //next T
            }

            // Normalization ------------------------------------------------------------ //
            for (T = 0; T < nVertex; T++) {
                VNormal[T] = VNormal[T]/max(1.0e-16, norm2(VNormal[T]));
            } //next T

            return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::GenerateVNormals(
        dvecarr3E   &X,
        unsigned char flag
        ) {

    // ========================================================================== //
    // void Class_SurfTri::GenerateVNormals(                                      //
    //     dvecarr3E   &X,                                                        //
    //     unsigned char flag)                                                    //
    //                                                                            //
    // Generate vertex normals using algorithm #1 or #2 depending on the value    //
    // of the input flag:                                                         //
    // algorithm #1 (flag = 0): compute normals as the weighted average of        //
    //                          normals of incident simplicies.                   //
    // algorithm #2 (flag = 1): compute normals as the weighted average of        //
    //                          normals of incident edges.                        //
    // Vertex coordinate list is provided externally                              //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X        : dvecarr3E, external vertex list. X[i][0], X[i][1], ... are    //
    //              x, y, ... coordinates of the i-th vertex.                     //
    // - flag    : unsigend char (default 0), flag for algorithm selection        //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Parameters
    double const    pi = 3.1415926535897932;

    // Local variables
    int             nE, nV = X.size();
    darray3E        tmp;

    // Counters
    int             V, T, i, j, m;

    // ========================================================================== //
    // COMPUTE NORMALS IF NOT ALREADY COMPUTED                                    //
    // ========================================================================== //
    if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
        GenerateNormals(X);
    }

    // ========================================================================== //
    // COMPUTE EDGES NORMALS IF NOT ALREADY COMPUTED                              //
    // ========================================================================== //
    if ((ENormal.size() == 0) || (ENormal.size() < Edge.size())) {
        GenerateENormals(X);
    }

    // ========================================================================== //
    // RESIZE INPUT VARIABLES                                                     //
    // ========================================================================== //
    tmp.fill(0.);
    VNormal.resize(nV, tmp);

    // ========================================================================== //
    // COMPUTE VERTEX NORMALS                                                     //
    // ========================================================================== //

    // Compute vertex normals using algorithm #1 -------------------------------- //
    if (flag == 1) {

        // Scope variables
        // none

        // Compute vertex normals
        nE = Edge.size();
        for (T = 0; T < nE; T++) {
            m = Edge[T].size();
            for (i = 0; i < m; i++) {
                V = Edge[T][i];
                VNormal[V] = VNormal[V] + ENormal[T];
            } //next i
        } //next T
    }

    // Compute vertex normals using algorithm #2 -------------------------------- //
    else if (flag == 0) {

        // Scope variables
        double          angle;

        // Compute vertex normals
        for (T = 0; T < nSimplex; ++T) {
            m = Simplex[T].size();
            for (i = 0; i < m; ++i) {
                V = Simplex[T][i];
                Angle(T, angle, i);
                VNormal[V] = VNormal[V] + 0.5 * angle * Normal[T]/pi;
            } //next i
        } //next T

    }


    // Normalization ------------------------------------------------------------ //
    for (T = 0; T < nV; T++) {
        VNormal[T] = VNormal[T]/norm2(VNormal[T]);
    } //next T

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BuildAdjacency(
        void
        ) {

    // ========================================================================== //
    // void Class_SurfTri::BuildAdjacency(                                        //
    //     void)                                                                  //
    //                                                                            //
    // Build simplex-simplex adjacency.                                           //
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
    bool                    check;
    bvector2D               flag(nSimplex);
    ivector2D               V2S(nVertex);
    ivector1D::iterator     index_;

    // Counters
    int                     V, T, A;
    int                     i, j, k, ii, jj;
    int                     n, m, p;

    // ========================================================================== //
    // INITIALIZE VARIABLES                                                       //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        flag[T].resize(Simplex[T].size(), false);
    } //next T

    // ========================================================================== //
    // BUILD VERTEX-> SIMPLEX MAP                                                 //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            V = Simplex[T][i];
            V2S[V].push_back(T);
        } //next i
    } //next T

    // ========================================================================== //
    // BUILD ADJACENCIES                                                          //
    // ========================================================================== //

    // Resize Adjacency --------------------------------------------------------- //
    ClearAdjacency();
    ReshapeAdjacency();

    // Loop over simplex -------------------------------------------------------- //
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            flag[T][i] = true;
            V = Simplex[T][i];
            j = (i+1) % m;
            n = V2S[V].size();
            for (k = 0; k < n; k++) {
                A = V2S[V][k];
                if (A != T) {
                    p = Simplex[A].size();
                    if (((p > 2) && (m > 2)) || ((p == 2) && (m == 2))) {
                        ii = 0;
                        check = false;
                        while ((!check) && (ii < p)) {
                            jj = (ii+1) % p;
                            if (m == 2) {
                                check = (Simplex[T][i] == Simplex[A][ii]);
                            }
                            else {
                                check = (((Simplex[T][i] == Simplex[A][jj])
                                            && (Simplex[T][j] == Simplex[A][ii]))
                                        || ((Simplex[T][i] == Simplex[A][ii])
                                            && (Simplex[T][j] == Simplex[A][jj])));
                            }
                            ii++;
                        } //next ii
                        ii--;
                        if ((check) && (!flag[A][ii])) {
                            if (Adjacency[T][i][0] == -1) {
                                Adjacency[T][i][0] = A;
                            }
                            else {
                                index_ = find(Adjacency[T][i].begin(), Adjacency[T][i].end(), A);
                                if (index_ == Adjacency[T][i].end()) {
                                    Adjacency[T][i].push_back(A);
                                }
                            }
                            if (Adjacency[A][ii][0] == -1) {
                                Adjacency[A][ii][0] = T;
                            }
                            else {
                                index_ = find(Adjacency[A][ii].begin(), Adjacency[A][ii].end(), T);
                                if (index_ == Adjacency[A][ii].end()) {
                                    Adjacency[A][ii].push_back(T);
                                }
                            }
                        }
                    }
                }
            } //next k

        } //next i
    } //next T


    return; }

    // -------------------------------------------------------------------------- //
    void Class_SurfTri::BuildAdjacency(
            dvecarr3E   &X
            ) {

        // ========================================================================== //
        // void Class_SurfTri::BuildAdjacency(                                        //
        //     dvecarr3E   &X)                                                        //
        //                                                                            //
        // Build simplex-simplex adjacency. Vertex list is provided externally.       //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - X      : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are    //
        //            the x, y, ... coordinates of the i-th vertex.                   //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        bool                    check;
        int                     nV = X.size();
        bvector2D               flag(nSimplex);
        ivector2D               V2S(nV);
        ivector1D::iterator     index_;

        // Counters
        int                     V, T, A;
        int                     i, j, k, ii, jj;
        int                     n, m, p;

        // ========================================================================== //
        // INITIALIZE VARIABLES                                                       //
        // ========================================================================== //
        for (T = 0; T < nSimplex; T++) {
            flag[T].resize(Simplex[T].size(), false);
        } //next T

        // ========================================================================== //
        // BUILD VERTEX-> SIMPLEX MAP                                                 //
        // ========================================================================== //
        for (T = 0; T < nSimplex; T++) {
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                V = Simplex[T][i];
                V2S[V].push_back(T);
            } //next i
        } //next T

        // ========================================================================== //
        // BUILD ADJACENCIES                                                          //
        // ========================================================================== //

        // Resize Adjacency --------------------------------------------------------- //
        ClearAdjacency();
        ReshapeAdjacency();

        // Loop over simplex -------------------------------------------------------- //
        for (T = 0; T < nSimplex; T++) {
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                flag[T][i] = true;
                V = Simplex[T][i];
                j = (i+1) % m;
                n = V2S[V].size();
                for (k = 0; k < n; k++) {
                    A = V2S[V][k];
                    if (A != T) {
                        p = Simplex[A].size();
                        if (((p > 2) && (m > 2)) || ((p == 2) && (m == 2))) {
                            ii = 0;
                            check = false;
                            while ((!check) && (ii < p)) {
                                jj = (ii+1) % p;
                                if (m == 2) {
                                    check = (Simplex[T][i] == Simplex[A][ii]);
                                }
                                else {
                                    check = (((Simplex[T][i] == Simplex[A][jj])
                                                && (Simplex[T][j] == Simplex[A][ii]))
                                            || ((Simplex[T][i] == Simplex[A][ii])
                                                && (Simplex[T][j] == Simplex[A][jj])));
                                }
                                ii++;
                            } //next ii
                            ii--;
                            if ((check) && (!flag[A][ii])) {
                                if (Adjacency[T][i][0] == -1) {
                                    Adjacency[T][i][0] = A;
                                }
                                else {
                                    index_ = find(Adjacency[T][i].begin(), Adjacency[T][i].end(), A);
                                    if (index_ == Adjacency[T][i].end()) {
                                        Adjacency[T][i].push_back(A);
                                    }
                                }
                                if (Adjacency[A][ii][0] == -1) {
                                    Adjacency[A][ii][0] = T;
                                }
                                else {
                                    index_ = find(Adjacency[A][ii].begin(), Adjacency[A][ii].end(), T);
                                    if (index_ == Adjacency[A][ii].end()) {
                                        Adjacency[A][ii].push_back(T);
                                    }
                                }
                            }
                        }
                    }
                } //next k

            } //next i
        } //next T


        return; }

        // -------------------------------------------------------------------------- //
        void Class_SurfTri::UpdateAdjacency(
                ivector1D   &list
                ) {

            // ========================================================================== //
            // void Class_SurfTri::UpdateAdjacency(                                       //
            //     ivector1D   &list)                                                     //
            //                                                                            //
            // Update simplex-simplex adjacency.                                          //
            // ========================================================================== //
            // INPUT                                                                      //
            // ========================================================================== //
            // - list      : ivector1D, list of simplicies to be updated                  //
            // ========================================================================== //
            // OUTPUT                                                                     //
            // ========================================================================== //
            // - none                                                                     //
            // ========================================================================== //

            // ========================================================================== //
            // VARIABLES DECLARATION                                                      //
            // ========================================================================== //

            // Local variables
            int                     nS = list.size();
            bool                    check;
            bvector2D               flag(nSimplex);
            ivector2D               V2S(nVertex);
            ivector1D::iterator     index_;

            // Counters
            int                     V, T, A;
            int                     i, j, k, ii, jj;
            int                     n, m, p;

            // ========================================================================== //
            // INITIALIZE VARIABLES                                                       //
            // ========================================================================== //
            for (T = 0; T < nSimplex; T++) {
                flag[T].resize(Simplex[T].size(), false);
            } //next T

            // ========================================================================== //
            // BUILD VERTEX-> SIMPLEX MAP                                                 //
            // ========================================================================== //
            for (T = 0; T < nSimplex; T++) {
                m = Simplex[T].size();
                for (i = 0; i < m; i++) {
                    V = Simplex[T][i];
                    V2S[V].push_back(T);
                } //next i
            } //next T

            // ========================================================================== //
            // BUILD ADJACENCIES                                                          //
            // ========================================================================== //

            // Resize Adjacency --------------------------------------------------------- //
            for (T = 0; T < nS; T++) {
                m = Simplex[list[T]].size();
                Adjacency[list[T]].resize(m, ivector1D(1, -1));
                for (i = 0; i < m; i++) {
                    Adjacency[list[T]][i][0] = -1;
                } //next i
            } //next T

            // Loop over simplex -------------------------------------------------------- //
            for (T = 0; T < nS; T++) {
                m = Simplex[list[T]].size();
                for (i = 0; i < m; i++) {
                    flag[list[T]][i] = true;
                    V = Simplex[list[T]][i];
                    j = (i+1) % m;
                    n = V2S[V].size();
                    for (k = 0; k < n; k++) {
                        A = V2S[Simplex[list[T]][i]][k];
                        if (A != list[T]) {
                            p = Simplex[A].size();
                            if (((p > 2) && (m > 2)) || ((p == 2) && (m == 2))) {
                                ii = 0;
                                check = false;
                                while ((!check) && (ii < p)) {
                                    jj = (ii+1) % p;
                                    if (m == 2) {
                                        check = (Simplex[list[T]][i] == Simplex[A][ii]);
                                    }
                                    else {
                                        check = (((Simplex[list[T]][i] == Simplex[A][jj])
                                                    && (Simplex[list[T]][j] == Simplex[A][ii]))
                                                || ((Simplex[list[T]][i] == Simplex[A][ii])
                                                    && (Simplex[list[T]][j] == Simplex[A][jj])));
                                    }
                                    ii++;
                                } //next ii
                                ii--;
                                if (check && (!flag[A][ii])) {
                                    if (Adjacency[list[T]][i][0] == -1) {
                                        Adjacency[list[T]][i][0] = A;
                                    }
                                    else {
                                        index_ = find(Adjacency[list[T]][i].begin(), Adjacency[list[T]][i].end(), A);
                                        if (index_ == Adjacency[list[T]][i].end()) {
                                            Adjacency[list[T]][i].push_back(A);
                                        }
                                    }
                                    if (Adjacency[A][ii][0] == -1) {
                                        Adjacency[A][ii][0] = list[T];
                                    }
                                    else {
                                        index_ = find(Adjacency[A][ii].begin(), Adjacency[A][ii].end(), list[T]);
                                        if (index_ == Adjacency[A][ii].end()) {
                                            Adjacency[A][ii].push_back(list[T]);
                                        }
                                    }
                                }
                            }
                        }
                    } //next k

                } //next i
            } //next T


            return; }

            // -------------------------------------------------------------------------- //
            void Class_SurfTri::UpdateAdjacency(
                    dvecarr3E   &X,
                    ivector1D   &list
                    ) {

                // ========================================================================== //
                // void Class_SurfTri::UpdateAdjacency(                                       //
                //     dvecarr3E   &X,                                                        //
                //     ivector1D   &list)                                                     //
                //                                                                            //
                // Update simplex-simplex adjacency. Vertex list is provided externally.      //
                // ========================================================================== //
                // INPUT                                                                      //
                // ========================================================================== //
                // - X         : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are //
                //               the x, y, ... coordinates of the i-th vertex.                //
                // - list      : ivector1D, list of simplicies to be updated                  //
                // ========================================================================== //
                // OUTPUT                                                                     //
                // ========================================================================== //
                // - none                                                                     //
                // ========================================================================== //

                // ========================================================================== //
                // VARIABLES DECLARATION                                                      //
                // ========================================================================== //

                // Local variables
                int                     nV = X.size();
                int                     nS = list.size();
                bool                    check;
                bvector2D               flag(nSimplex);
                ivector2D               V2S(nV);
                ivector1D::iterator     index_;

                // Counters
                int                     V, T, A;
                int                     i, j, k, ii, jj;
                int                     n, m, p;

                // ========================================================================== //
                // INITIALIZE VARIABLES                                                       //
                // ========================================================================== //
                for (T = 0; T < nSimplex; T++) {
                    flag[T].resize(Simplex[T].size(), false);
                } //next T

                // ========================================================================== //
                // BUILD VERTEX-> SIMPLEX MAP                                                 //
                // ========================================================================== //
                for (T = 0; T < nSimplex; T++) {
                    m = Simplex[T].size();
                    for (i = 0; i < m; i++) {
                        V = Simplex[T][i];
                        V2S[V].push_back(T);
                    } //next i
                } //next T

                // ========================================================================== //
                // BUILD ADJACENCIES                                                          //
                // ========================================================================== //

                // Resize Adjacency --------------------------------------------------------- //
                for (T = 0; T < nS; T++) {
                    m = Simplex[list[T]].size();
                    Adjacency[list[T]].resize(m, ivector1D(1, -1));
                    for (i = 0; i < m; i++) {
                        Adjacency[list[T]][i][0] = -1;
                    } //next i
                } //next T

                // Loop over simplex -------------------------------------------------------- //
                for (T = 0; T < nS; T++) {
                    m = Simplex[list[T]].size();
                    for (i = 0; i < m; i++) {
                        flag[list[T]][i] = true;
                        V = Simplex[list[T]][i];
                        j = (i+1) % m;
                        n = V2S[V].size();
                        for (k = 0; k < n; k++) {
                            A = V2S[Simplex[list[T]][i]][k];
                            if (A != list[T]) {
                                p = Simplex[A].size();
                                if (((p > 2) && (m > 2)) || ((p == 2) && (m == 2))) {
                                    ii = 0;
                                    check = false;
                                    while ((!check) && (ii < p)) {
                                        jj = (ii+1) % p;
                                        if (m == 2) {
                                            check = (Simplex[list[T]][i] == Simplex[A][ii]);
                                        }
                                        else {
                                            check = (((Simplex[list[T]][i] == Simplex[A][jj])
                                                        && (Simplex[list[T]][j] == Simplex[A][ii]))
                                                    || ((Simplex[list[T]][i] == Simplex[A][ii])
                                                        && (Simplex[list[T]][j] == Simplex[A][jj])));
                                        }
                                        ii++;
                                    } //next ii
                                    ii--;
                                    if (check && (!flag[A][ii])) {
                                        if (Adjacency[list[T]][i][0] == -1) {
                                            Adjacency[list[T]][i][0] = A;
                                        }
                                        else {
                                            index_ = find(Adjacency[list[T]][i].begin(), Adjacency[list[T]][i].end(), A);
                                            if (index_ == Adjacency[list[T]][i].end()) {
                                                Adjacency[list[T]][i].push_back(A);
                                            }
                                        }
                                        if (Adjacency[A][ii][0] == -1) {
                                            Adjacency[A][ii][0] = list[T];
                                        }
                                        else {
                                            index_ = find(Adjacency[A][ii].begin(), Adjacency[A][ii].end(), list[T]);
                                            if (index_ == Adjacency[A][ii].end()) {
                                                Adjacency[A][ii].push_back(list[T]);
                                            }
                                        }
                                    }
                                }
                            }
                        } //next k

                    } //next i
                } //next T


                return; }

                // -------------------------------------------------------------------------- //
                void Class_SurfTri::BuildEdges(
                        void
                        ) {

                    // ========================================================================== //
                    // void Class_SurfTri::BuildEdges(                                            //
                    //     void)                                                                  //
                    //                                                                            //
                    // Build edge-vertex connectivity.                                            //
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
                    // none

                    // Counters
                    int             A, T, E;
                    int             i, j, e;
                    int             m, n;

                    // ========================================================================== //
                    // BUILD ADJACENCIES IF NOT ALREADY BUILT                                     //
                    // ========================================================================== //
                    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
                        BuildAdjacency();
                    }

                    // ========================================================================== //
                    // RESIZE DATA STRUCTURE FOR EDGE-VERTEX AND SIMPLEX-EDGE CONNECITIVITY       //
                    // ========================================================================== //
                    Edge.resize(CountEdges());
                    ReshapeSimplex2Edge();

                    // ========================================================================== //
                    // BUILD EDGE LIST                                                            //
                    // ========================================================================== //
                    E = 0;
                    for (T = 0; T < nSimplex; T++) {
                        n = Simplex[T].size();
                        for (i = 0; i < n; i++) {
                            if (Simplex2Edge[T][i] == -1) {

                                // Update edge-vertex connectivity
                                if (n == 2) {
                                    Edge[E].resize(1, Simplex[T][i]);
                                }
                                else {
                                    Edge[E].resize(2);
                                    Edge[E][0] = Simplex[T][i];
                                    Edge[E][1] = Simplex[T][(i+1) % n];
                                }

                                // Update simplex-edge connectivity
                                Simplex2Edge[T][i] = E;
                                m = Adjacency[T][i].size();
                                for (j = 0; j < m; j++) {
                                    A = Adjacency[T][i][j];
                                    if (A >= 0) {
                                        e = edge(A, T);
                                        Simplex2Edge[A][e] = E;
                                    }
                                } //next j

                                // Update edge counter
                                E++;
                            }
                        } //next i
                    } //next T

                    return; }

                    // ========================================================================== //
                    // COUNTERS                                                                   //
                    // ========================================================================== //

                    // -------------------------------------------------------------------------- //
                    int Class_SurfTri::CountIsolatedVertex(
                            void
                            ) {

                        // ========================================================================== //
                        // int Class_SurfTri::CountIsolatedVertex(                                    //
                        //     void)                                                                  //
                        //                                                                            //
                        // Count isolated nodes in the tasselation. A node is isolated if and only if //
                        // there are no simplicies in the tassellation having a vertex in that node.  //
                        // ========================================================================== //
                        // INPUT                                                                      //
                        // ========================================================================== //
                        // - none                                                                     //
                        // ========================================================================== //
                        // OUTPUT                                                                     //
                        // ========================================================================== //
                        // - n    : int, number of isolated vertexes                                  //
                        // ========================================================================== //

                        // ========================================================================== //
                        // VARIABLES DECLARATION                                                      //
                        // ========================================================================== //

                        // Local variables
                        int                       n;
                        ivector1D                 List;

                        // Counters
                        // none

                        // ========================================================================== //
                        // COUNT ISOLATED VERTEXES                                                    //
                        // ========================================================================== //
                        List = FindIsolatedVertex();
                        n = List.size();

                        return (n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountIsolatedVertex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountIsolatedVertex(                                    //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Count isolated vertex in the tasselation. A node is isolated if and only   //
    // if there are no simplicies in the tassellation having a vertex in that     //
    // node. Vertex coordinate list is provided externally.                       //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are the  //
    //          x, y, ... coordinates of the i-th node.                           //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of isolated vertexes                                  //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                       n;
    ivector1D                 List;

    // Counters
    int                       i, I_;

    // ========================================================================== //
    // COUNT ISOLATED VERTEXES                                                    //
    // ========================================================================== //
    List = FindIsolatedVertex(X);
    n = List.size();

    return (n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountFreeVertex(
        void
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountFreeVertex(                                        //
    //     void)                                                                  //
    //                                                                            //
    // Count free vertexes in the tasselation. A free vertex is a vertex          //
    // connected to a free edge.                                                  //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of free vertexes                                      //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int            n;
    ivector1D      List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT FREE VERTICES                                                        //
    // ========================================================================== //
    List = FindFreeVertex();
    n = List.size();

    return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountFreeVertex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountFreeVertex(                                        //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Count free vertexes in the tasselation. A free vertex is a vertex          //
    // connected to a free edge. Vertex coordinate list is provided externally    //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are the  //
    //          x, y, ... coordinates of the i-th node.                           //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of free vertexes                                      //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int            n;
    ivector1D      List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT FREE VERTEXES                                                        //
    // ========================================================================== //
    List = FindFreeVertex(X);
    n = List.size();

    return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountDoubleVertex(
        void
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountDoubleVertex(                                      //
    //     void)                                                                  //
    //                                                                            //
    // Count duplicated vertex in the tasselation. A vertex is duplicated if      //
    // there exists a simplex in the tasselation having a vertex with the same    //
    // coordinates (within a prescribed tollerance) of the given vertex.          //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of double vertex                                      //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                             n;
    ivector1D                       List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT DOUBLED VERTEXES                                                     //
    // ========================================================================== //
    List = FindDoubleVertex();
    n = List.size();

    return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountDoubleVertex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountDoubleVertex(                                      //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Count duplicated vertex in the tasselation. A vertex is duplicated if      //
    // there exists a simplex in the tasselation having a vertex with the same    //
    // coordinates (within a  prescribed tollerance) of the given vertex.         //
    // Vertex coordinate list is provided externally.                             //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, with vertex coordinate list. X[i][0], X[i][1], ...     //
    //          are the x, y, ... coordinates of the i-th node.                   //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of double vertex                                      //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                             n;
    ivector1D                       List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT DOUBLED VERTEXES                                                     //
    // ========================================================================== //
    List = FindDoubleVertex(X);
    n = List.size();

    return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountIsolatedSimplex(
        void
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountIsolatedSimplex(                                   //
    //     void)                                                                  //
    //                                                                            //
    // Count isolated simplex in the tasselation. A isolated simplex is a simplex //
    // whose vertex are not shared by any of the other simplex in the tasselation //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n     : int, number of isolated simplicies                               //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int              n;
    ivector1D        List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT ISOLATED SIMPLICIES                                                  //
    // ========================================================================== //
    List = FindIsolatedSimplex();
    n = List.size();

    return(n); }

    // -------------------------------------------------------------------------- //
    int Class_SurfTri::CountIsolatedSimplex(
            dvecarr3E   &X
            ) {

        // ========================================================================== //
        // int Class_SurfTri::CountIsolatedSimplex(                                   //
        //     dvecarr3E   &X)                                                        //
        //                                                                            //
        // Count isolated simplex in the tasselation. A isolated simplex is a simplex //
        // whose vertex are not shared by any of the other simplex in the tasselation //
        // Vertex coordinate list is provided externally                              //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1],  ... are the //
        //          x, y, ... coordinates of the i-th node.                           //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - n     : int, number of isolated simplicies                               //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        int              n;
        ivector1D        List;

        // Counters
        // none

        // ========================================================================== //
        // COUNT ISOLATED SIMPLICIES                                                  //
        // ========================================================================== //
        List = FindIsolatedSimplex(X);
        n = List.size();

        return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountFreeSimplex(
        void
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountFreeSimplex(                                       //
    //     void)                                                                  //
    //                                                                            //
    // Count free simplex in the tasselation. A free simplex is a simplex having  //
    // at least one free edge.                                                    //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n     : int, number of free simplicies                                   //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int              n;
    ivector1D        List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT FREE SIMPLICIES                                                      //
    // ========================================================================== //
    List = FindFreeSimplex();
    n = List.size();

    return(n); }

    // -------------------------------------------------------------------------- //
    int Class_SurfTri::CountDoubleSimplex(
            void
            ) {

        // ========================================================================== //
        // int Class_SurfTri::CountDoubleSimplex(                                     //
        //     void)                                                                  //
        //                                                                            //
        // Count duplicated simplicies in the tasselation. A duplicated simplex is a  //
        // simplex whose vertices have coordinates coincident (within a prescribed    //
        // tolerance) with coordinates of vertices of another simplex in the          //
        // tasselation.                                                               //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - n    : int, number of duplicated simplicies.                             //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        int               n;
        ivector1D         List;

        // Counters
        // none

        // ========================================================================== //
        // COUNT DOUBLE SIMPLICIES                                                    //
        // ========================================================================== //
        List = FindDoubleSimplex();
        n = List.size();

        return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountDoubleSimplex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountDoubleSimplex(                                     //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Count duplicated simplicies in the tasselation. A duplicated simplex is a  //
    // simplex whose vertices have coordinates coincident (within a prescribed    //
    // tolerance) with coordinates of vertexes of another simplex in the          //
    // tasselation. Vertex list is provided externally.                           //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1],  ... are the //
    //          x, y, ... coordinates of the i-th node.                           //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of duplicated simplicies.                             //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int               n;
    ivector1D         List;

    // Counters
    // none

    // ========================================================================== //
    // COUNT DOUBLE SIMPLICIES                                                    //
    // ========================================================================== //
    List = FindDoubleSimplex(X);
    n = List.size();

    return(n); };

// ----------------------------------------------------------------------------------- //
int Class_SurfTri::CountTrueDoubleSimplex(
        void
        ) {

    // =================================================================================== //
    // int Class_SurfTri::CountTrueDoubleSimplex(                                          //
    //     void)                                                                           //
    //                                                                                     //
    // Count duplicated simplicies in the tasselation. A duplicated simplex is a           //
    // simplex whose vertexes have coordinates coincident (within a prescribed tolerance)  //
    // with coordinates of vertexes of another simplex in the tasselation,                 //
    // without any distinction of vertex ordering. Check is meaningful on simplicies of the//
    // same kind, of course.							       //		
    // =================================================================================== //
    // INPUT                                                                               //
    // =================================================================================== //
    // - none                                                                              //
    // =================================================================================== //
    // OUTPUT                                                                              //
    // =================================================================================== //
    // - n    : int, number of true duplicated simplicies.                                 //
    // =================================================================================== //

    // =================================================================================== //
    // VARIABLES DECLARATION                                                               //
    // =================================================================================== //

    // Local variables
    int               n;
    ivector1D         List;

    // Counters
    // none

    // =================================================================================== //
    // COUNT DOUBLE SIMPLICIES                                                             //
    // =================================================================================== //
    List = FindTrueDoubleSimplex();
    n = List.size();

    return(n); };

// ----------------------------------------------------------------------------------- //
int Class_SurfTri::CountTrueDoubleSimplex(
        dvecarr3E &X
        ) {

    // =================================================================================== //
    // int Class_SurfTri::CountTrueDoubleSimplex(                                          //
    //     dvecarr3E &X)                                                                   //
    //                                                                                     //
    // Count duplicated simplicies in the tasselation. A duplicated simplex is a           //
    // simplex whose vertexes have coordinates coincident (within a prescribed tolerance)  //
    // with coordinates of vertexes of another simplex in the tasselation,                 //
    // without any distinction of vertex ordering. Check is meaningful on simplicies of the//
    // same kind, of course.																																																								       //		
    // =================================================================================== //
    // INPUT                                                                               //
    // =================================================================================== //
    // - X    : [nVertex-by-dim] dvecarr3E, with vertex coordinate list. X[i][0], X[i][1], //
    //          ... are the x, y, ... coordinates of the i-th node.                        //
    // =================================================================================== //
    // OUTPUT                                                                              //
    // =================================================================================== //
    // - n    : int, number of true duplicated simplicies.                                 //
    // =================================================================================== //

    // =================================================================================== //
    // VARIABLES DECLARATION                                                               //
    // =================================================================================== //

    // Local variables
    int               n;
    ivector1D         List;

    // Counters
    // none

    // =================================================================================== //
    // COUNT DOUBLE SIMPLICIES                                                             //
    // =================================================================================== //
    List = FindTrueDoubleSimplex(X);
    n = List.size();

    return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountEdges(
        void
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountEdges(                                             //
    //     void)                                                                  //
    //                                                                            //
    // Count edges in a surface tassselation.                                     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n     : int, number of edges                                             //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int    n, m, p;
    double dcounter = 0.0;

    // Counters
    int    i, T;

    // ========================================================================== //
    // BUILD ADJACENCY IF NOT ALREADY BUILT                                       //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency();
    }

    // ========================================================================== //
    // COUNT EDGES                                                                //
    // ========================================================================== //

    // Loop over simplicies
    for (T = 0; T < nSimplex; T++) {
        p = Simplex[T].size();
        for (i = 0; i < p; i++) {
            m = Adjacency[T][i].size();
            if (Adjacency[T][i][0] >= 0) {
                dcounter += 1.0/((double) (m + 1));
            }
            else {
                dcounter += 1.0;
            }
        } //next i
    } //next T
    n = (int) round(dcounter);

    return(n); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::CountFreeEdges(
        void
        ) {

    // ========================================================================== //
    // int Class_SurfTri::CountFreeEdges(                                         //
    //     void)                                                                  //
    //                                                                            //
    // Count free edges in a tasselation.                                         //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - n    : int, number of free edges.                                        //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int   n, m;

    // Counters
    int   i, j;

    // ========================================================================== //
    // BUILD ADJACENCY IF NOT ALREADY BUILT                                       //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency();
    }

    // ========================================================================== //
    // COUNT FREE EDGES                                                           //
    // ========================================================================== //

    // Loop over simplicies
    n = 0;
    for (i = 0; i < nSimplex; i++){
        m = Simplex[i].size();
        for (j = 0; j < m; j++) {
            if (Adjacency[i][j][0] < 0) {
                n++;
            }
        } //next j
    } //next i

    return(n); }

    // ========================================================================== //
    // FIND ALGORITHMS                                                            //
    // ========================================================================== //

    // -------------------------------------------------------------------------- //
    ivector1D Class_SurfTri::FindIsolatedVertex(
            void
            ) {

        // ========================================================================== //
        // ivector1D Class_SurfTri::FindIsolatedVertex(                               //
        //     void)                                                                  //
        //                                                                            //
        // Find isolated vertex in the tasselation. A node is isolated if there exist //
        // no simplex in the tasselation having a vertex in that node                 //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - list   : ivector1D, global indices of isolated vertexes                  //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        bvector1D       flag(nVertex, true);
        ivector1D       list;

        // Counters
        int             i, I_, m;

        // ========================================================================== //
        // FIND ISOLATED VERTEX                                                       //
        // ========================================================================== //

        // Loop over simplex
        for (I_ = 0; I_ < nSimplex; I_++) {
            m = Simplex[I_].size();
            for (i = 0; i < m; i++) {
                flag[Simplex[I_][i]] = false;
            }  //next i
        } //next I_

        // Loop over vertexes
        list.resize(count(flag.begin(), flag.end(), true));
        i = 0;
        for (I_ = 0; I_ < nVertex; I_++) {
            if (flag[I_]) {
                list[i] = I_;
                i++;
            }
        } //next I_

        return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindIsolatedVertex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindIsolatedVertex(                               //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Find isolated vertex in the tasselation. A node is isolated if there exist //
    // no simplex in the tasselation having a vertex in that node. Vertex         //
    // coordinate list is provided externally                                     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X      : dvecarr3E, vertex coordinates. X[i][0], X[i][1], ... are the x, //
    //            y, ... coordinates of the i-th vertex.                          //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, global indices of isolated vertexes                  //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D       flag(X.size(), true);
    ivector1D       list;

    // Counters
    int             i, m, I_;

    // ========================================================================== //
    // FIND ISOLATED VERTEX                                                       //
    // ========================================================================== //

    // Loop over simplex
    for (I_ = 0; I_ < nSimplex; I_++) {
        m = Simplex[I_].size();
        for (i = 0; i < m; i++) {
            flag[Simplex[I_][i]] = false;
        }  //next i
    } //next I_

    // Loop over vertexes
    list.resize(count(flag.begin(), flag.end(), true));
    i = 0;
    for (I_ = 0; I_ < X.size(); I_++) {
        if (flag[I_]) {
            list[i] = I_;
            i++;
        }
    } //next I_

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindFreeVertex(
        void
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindFreeVertex(                                   //
    //     void)                                                                  //
    //                                                                            //
    // Find free vertex in the tasselation. A free vertex is a vertex on          //
    // tasselation boundaries.                                                    //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   :ivector1D, global indices of free vertices                       //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D       flag(nVertex, false);
    ivector1D       list;

    // Counters
    int             i, j, m, I_;

    // ========================================================================== //
    // COMPUTE THE ADJACECNY MATRIX (IF NOT AREADY BUILT)                         //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency();
    }

    // ========================================================================== //
    // COUNT FREE VERTEXES                                                        //
    // ========================================================================== //

    // Loop over simplex
    for (I_ = 0; I_ < nSimplex; I_++) {
        m = Simplex[I_].size();
        for (i = 0; i < m; i++) {
            if (Simplex[I_].size() == 2) {
                if (Adjacency[I_][i][0] < 0) {
                    flag[Simplex[I_][i]] = true;
                }
            }
            else {
                j = (i+1) % Simplex[I_].size();
                if (Adjacency[I_][i][0] < 0) {
                    flag[Simplex[I_][i]] = true;
                    flag[Simplex[I_][j]] = true;
                }
            }
        } //next i
    } //next I_

    // Loop over vertices
    list.resize(count(flag.begin(), flag.end(), true));
    i = 0;
    for (I_ = 0; I_ < nVertex; I_++) {
        if (flag[I_]) {
            list[i] = I_;
            i++;
        }
    }

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindFreeVertex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindFreeVertex(                                   //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Find free vertex in the tasselation. A free vertex is a vertex on          //
    // tasselation boundaries. Vertex coordinate list is provided exteernally     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are the  //
    //          x, y, ... coordinates of the i-th node.                           //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   :[ ivector1D, global indices of free vertices                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D       flag(X.size(), false);
    ivector1D       list;

    // Counters
    int             i, j, m, I_;

    // ========================================================================== //
    // COMPUTE THE ADJACECNY MATRIX (IF NOT AREADY BUILT)                         //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency(X);
    }

    // ========================================================================== //
    // COUNT FREE VERTEXES                                                        //
    // ========================================================================== //

    // Loop over simplicies
    for (I_ = 0; I_ < nSimplex; I_++) {
        m = Simplex[I_].size();
        for (i = 0; i < m; i++) {
            if (Simplex[I_].size() == 2) {
                if (Adjacency[I_][i][0] < 0) {
                    flag[Simplex[I_][i]] = true;
                }
            }
            else {
                j = (i+1) % Simplex[I_].size();
                if (Adjacency[I_][i][0] < 0) {
                    flag[Simplex[I_][i]] = true;
                    flag[Simplex[I_][j]] = true;
                }
            }
        } //next i
    } //next I_

    // Loop over vertices
    list.resize(count(flag.begin(), flag.end(), true));
    i = 0;
    for (I_ = 0; I_ < X.size(); I_++) {
        if (flag[I_]) {
            list[i] = I_;
            i++;
        }
    }

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindDoubleVertex(
        int          n
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindDoubleVertex(                                 //
    //     int          n)                                                        //
    //                                                                            //
    // Find duplicated vertexes in the tasselation. A vertex is duplicated if     //
    // there exist a simplex in the tasselation having a vertex with coordinates  //
    // within a prescribed tolerance from the given vertex.                       //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - n       : int (optional), number of bins for vertex sorting              //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - doublev : ivector1D, global indices of isolated vertices                 //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Variables declaration
    int          m, ncell;
    //ht int          dim = Vertex[0].size();
    bvector1D    flag(nVertex, false);
    ivector1D    index(nVertex, -1), idummy1D(2, -1);
    ivector1D    doublev;
    ivector3D    cell;

    // Counters
    int          I_, C, S, V, W;
    int          i;

    // ========================================================================== //
    // INITIALIZE VARIABLES                                                       //
    // ========================================================================== //

    // Random number generator -------------------------------------------------- //
    //srand(time(NULL));
    srand(32767);

    // Resize variables --------------------------------------------------------- //
    //ht if      (dim == 2) { cell.resize(n*n); }
    //ht else if (dim == 3) { cell.resize(n*n*n); }

    cell.resize(n*n*n);

    // ========================================================================== //
    // SORT VERTICES ON BINS                                                      //
    // ========================================================================== //

    // Sort vertices ------------------------------------------------------------ //
    BinSortV(index, n);

    // Sort simplicies ---------------------------------------------------------- //
    for (I_ = 0; I_ < nSimplex; I_++) {
        m = Simplex[I_].size();
        for (i = 0; i < m; i++) {
            V = Simplex[I_][i];
            C = index[V];
            idummy1D[0] = I_;
            idummy1D[1] = i;
            cell[C].push_back(idummy1D);
        } //next i
    } //next I_

    // ========================================================================== //
    // FIND DOUBLE VERTICES                                                       //
    // ========================================================================== //
    ncell = cell.size();
    //ht if (dim == 2) {
    //ht     for (C = 0; C < ncell; C++) {
    //ht         m = cell[C].size();
    //ht         if (m > 0) {
    //ht 
    //ht             // Scope variables
    //ht             ivector1D    list;
    //ht             bitpit::KdTree<2, array<double>, int>     kd(m);
    //ht 
    //ht             // Randomize vertex insertion
    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
    //ht             for (I_ = 0; I_ < m; I_++) {
    //ht                 S = cell[C][list[I_]][0];
    //ht                 i = cell[C][list[I_]][1];
    //ht                 V = Simplex[S][i];
    //ht                 if (kd.exist(&Vertex[V], W) >= 0) {
    //ht                     if (!flag[V]) {
    //ht                         flag[V] = true;
    //ht                         doublev.push_back(V);
    //ht                     }
    //ht                 }
    //ht                 else {
    //ht                     flag[V] = true;
    //ht                     kd.insert(&Vertex[V], V);
    //ht                 }
    //ht             } //next I_
    //ht         }
    //ht     } //next C
    //ht }
    //ht else if (dim == 3) {
    //ht     for (C = 0; C < ncell; C++) {
    //ht         m = cell[C].size();
    //ht         if (m > 0) {
    //ht 
    //ht             // Scope variables
    //ht             ivector1D    list;
    //ht             bitpit::KdTree<3, double, int>     kd(m);
    //ht 
    //ht             // Randomize vertex insertion
    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
    //ht             for (I_ = 0; I_ < m; I_++) {
    //ht                 S = cell[C][list[I_]][0];
    //ht                 i = cell[C][list[I_]][1];
    //ht                 V = Simplex[S][i];
    //ht                 if (kd.exist(&Vertex[V], W) >= 0) {
    //ht                     if (!flag[V]) {
    //ht                         flag[V] = true;
    //ht                         doublev.push_back(V);
    //ht                     }
    //ht                 }
    //ht                 else {
    //ht                     flag[V] = true;
    //ht                     kd.insert(&Vertex[V], V);
    //ht                 }
    //ht             } //next I_
    //ht         }
    //ht     } //next C
    //ht }

    for (C = 0; C < ncell; C++) {
        m = cell[C].size();
        if (m > 0) {

            // Scope variables
            ivector1D    list;
            bitpit::KdTree<3, array<double,3>, int>     kd(m);

            // Randomize vertex insertion
            bitpit::utils::extractWithoutReplacement(m, m-1, list);
            for (I_ = 0; I_ < m; I_++) {
                S = cell[C][list[I_]][0];
                i = cell[C][list[I_]][1];
                V = Simplex[S][i];
                if (kd.exist(&Vertex[V], W) >= 0) {
                    if (!flag[V]) {
                        flag[V] = true;
                        doublev.push_back(V);
                    }
                }
                else {
                    flag[V] = true;
                    kd.insert(&Vertex[V], V);
                }
            } //next I_
        }
    } //next C


    return (doublev); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindDoubleVertex(
        dvecarr3E   &X,
        int          n
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindDoubleVertex(                                 //
    //     dvecarr3E   &X,                                                        //
    //     int          n)                                                        //
    //                                                                            //
    // Find duplicated vertexes in the tasselation. A vertex is duplicated if     //
    // there exist a simplex in the tasselation having a vertex with coordinates  //
    // within a prescribed tolerance from the given vertex.                       //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X       : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are   //
    //             the x, y, ... coordinates of the i-th vertex                   //
    // - n       : int (optional), number of bins for vertex sorting              //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - doublev : ivector1D, global indices of isolated vertices                 //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Variables declaration
    int          nV = X.size();
    int          m, ncell;
    //ht int          dim = X[0].size();
    bvector1D    flag(nV, false);
    ivector1D    index(nV, -1), idummy1D(2, -1);
    ivector1D    doublev;
    ivector3D    cell;

    // Counters
    int          I_, C, S, V, W;
    int          i;

    // ========================================================================== //
    // INITIALIZE VARIABLES                                                       //
    // ========================================================================== //

    // Random number generator -------------------------------------------------- //
    //srand(time(NULL));
    srand(32767);

    // Resize variables --------------------------------------------------------- //
    //ht if      (dim == 2) { cell.resize(n*n); }
    //ht else if (dim == 3) { cell.resize(n*n*n); }

    cell.resize(n*n*n); 

    // ========================================================================== //
    // SORT VERTICES ON BINS                                                      //
    // ========================================================================== //

    // Sort vertices ------------------------------------------------------------ //
    BinSortV(index, n);

    // Sort simplicies ---------------------------------------------------------- //
    for (I_ = 0; I_ < nSimplex; I_++) {
        m = Simplex[I_].size();
        for (i = 0; i < m; i++) {
            V = Simplex[I_][i];
            C = index[V];
            idummy1D[0] = I_;
            idummy1D[1] = i;
            cell[C].push_back(idummy1D);
        } //next i
    } //next I_

    // ========================================================================== //
    // FIND DOUBLE VERTICES                                                       //
    // ========================================================================== //
    ncell = cell.size();
    //ht if (dim == 2) {
    //ht     for (C = 0; C < ncell; C++) {
    //ht         m = cell[C].size();
    //ht         if (m > 0) {
    //ht 
    //ht             // Scope variables
    //ht             ivector1D    list;
    //ht             bitpit::KdTree<2, double, int>     kd(m);
    //ht 
    //ht             // Randomize vertex insertion
    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
    //ht             for (I_ = 0; I_ < m; I_++) {
    //ht                 S = cell[C][list[I_]][0];
    //ht                 i = cell[C][list[I_]][1];
    //ht                 V = Simplex[S][i];
    //ht                 if (kd.exist(&X[V], W) >= 0) {
    //ht                     if (!flag[V]) {
    //ht                         flag[V] = true;
    //ht                         doublev.push_back(V);
    //ht                     }
    //ht                 }
    //ht                 else {
    //ht                     flag[V] = true;
    //ht                     kd.insert(&X[V], V);
    //ht                 }
    //ht             } //next I_
    //ht         }
    //ht     } //next C
    //ht }
    //ht else if (dim == 3) {
    //ht     for (C = 0; C < ncell; C++) {
    //ht         m = cell[C].size();
    //ht         if (m > 0) {
    //ht 
    //ht             // Scope variables
    //ht             ivector1D    list;
    //ht             bitpit::KdTree<3, double, int>     kd(m);
    //ht 
    //ht             // Randomize vertex insertion
    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
    //ht             for (I_ = 0; I_ < m; I_++) {
    //ht                 S = cell[C][list[I_]][0];
    //ht                 i = cell[C][list[I_]][1];
    //ht                 V = Simplex[S][i];
    //ht                 if (kd.exist(&X[V], W) >= 0) {
    //ht                     if (!flag[V]) {
    //ht                         flag[V] = true;
    //ht                         doublev.push_back(V);
    //ht                     }
    //ht                 }
    //ht                 else {
    //ht                     flag[V] = true;
    //ht                     kd.insert(&X[V], V);
    //ht                 }
    //ht             } //next I_
    //ht         }
    //ht     } //next C
    //ht }

    for (C = 0; C < ncell; C++) {
        m = cell[C].size();
        if (m > 0) {

            // Scope variables
            ivector1D    list;
            bitpit::KdTree<3, array<double,3>, int>     kd(m);

            // Randomize vertex insertion
            bitpit::utils::extractWithoutReplacement(m, m-1, list);
            for (I_ = 0; I_ < m; I_++) {
                S = cell[C][list[I_]][0];
                i = cell[C][list[I_]][1];
                V = Simplex[S][i];
                if (kd.exist(&X[V], W) >= 0) {
                    if (!flag[V]) {
                        flag[V] = true;
                        doublev.push_back(V);
                    }
                }
                else {
                    flag[V] = true;
                    kd.insert(&X[V], V);
                }
            } //next I_
        }
    } //next C

    return (doublev); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindIsolatedSimplex(
        void
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindIsolatedSimplex(                              //
    //     void)                                                                  //
    //                                                                            //
    // Find isolated simplicies in the tasselation. A isolated simplex is a       //
    // simplex whose vertex are not shared by any other simplex in the            //
    // tasselation                                                                //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   :ivector1D, global indices of isolated simplicies                 //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D          flag(nSimplex, false);
    ivector1D          flag_V(nVertex, 0);
    ivector1D          list;
    int                flag_T;

    // Counters
    int                i, j, m;

    // ========================================================================== //
    // COUNT ISOLATED SIMPLICIES                                                  //
    // ========================================================================== //

    // Compute vertex valence --------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        m = Simplex[i].size();
        for (j = 0; j < m; j++) {
            flag_V[Simplex[i][j]] += 1;
        } //next j
    } //next i

    // Loop over simplicies ----------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        flag[i] = true;
        j = 0;
        m = Simplex[i].size();
        while ((flag[i] == true) && (j < m)) {
            flag[i] = (flag[i] && (flag_V[Simplex[i][j]] == 1));
            j++;
        } //next j
    } //next i

    // Simplex list ------------------------------------------------------------- //
    list.resize(count(flag.begin(), flag.end(), true));
    j = 0;
    for (i = 0; i < nSimplex; i++) {
        if (flag[i]) {
            list[j] = i;
            j++;
        }
    } //next i

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindIsolatedSimplex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindIsolatedSimplex(                              //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Find isolated simplicies in the tasselation. A isolated simplex is a       //
    // simplex whose vertex are not shared by any other simplex in the            //
    // tasselation. Vertex coordinate list is provided externally.                //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are the  //
    //          x, y, ... coordinates of the i-th node.                           //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, global indices of isolated simplicies                //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D          flag(nSimplex, false);
    ivector1D          flag_V(X.size(), 0);
    ivector1D          list;
    int                flag_T;

    // Counters
    int                i, j, m;

    // ========================================================================== //
    // COUNT ISOLATED SIMPLICIES                                                  //
    // ========================================================================== //

    // Compute vertex valence --------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        m = Simplex[i].size();
        for (j = 0; j < m; j++) {
            flag_V[Simplex[i][j]] += 1;
        } //next j
    } //next i

    // Loop over simplicies ----------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        flag[i] = true;
        j = 0;
        m = Simplex[i].size();
        while ((flag[i] == true) && (j < m)) {
            flag[i] = (flag[i] && (flag_V[Simplex[i][j]] == 1));
            j++;
        } //next j
    } //next i

    // Simplex list ------------------------------------------------------------- //
    list.resize(count(flag.begin(), flag.end(), true));
    j = 0;
    for (i = 0; i < nSimplex; i++) {
        if (flag[i]) {
            list[j] = i;
            j++;
        }
    } //next i

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindFreeSimplex(
        void
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindFreeSimplex(                                  //
    //     void)                                                                  //
    //                                                                            //
    // Find free simplex in the tasselation. A free simplex is a simplex having   //
    // at least one free edge.                                                    //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, with global indices of free simplicies               //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D          flag(nSimplex, false);
    ivector1D          list;

    // Counters
    int                i, j, m;

    // ========================================================================== //
    // COMPUTE ADJACENCY MATRIX IF NOT ALREADY COMPUTED                           //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency();
    }

    // ========================================================================== //
    // FIND FREE SIMPLICIES                                                       //
    // ========================================================================== //

    // Loop over simplicies  ---------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        m = Simplex[i].size();
        for (j = 0; j < m; j++) {
            flag[i] = (flag[i] || (Adjacency[i][j][0] < 0));
        } //next j
    } //next i

    // Simplex list ------------------------------------------------------------- //
    list.resize(count(flag.begin(), flag.end(), true));
    j = 0;
    for (i = 0; i < nSimplex; i++) {
        if (flag[i]) {
            list[j] = i;
            j++;
        }
    }

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindFreeSimplex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindFreeSimplex(                                  //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Find free simplex in the tasselation. A free simplex is a simplex having   //
    // at least one free edge. Vertex list is provided externally.                //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X      : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are    //
    //            the x, y, ... coordinates of the i-th node.                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, with global indices of free simplicies               //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bvector1D          flag(nSimplex, false);
    ivector1D          list;

    // Counters
    int                i, j, m;

    // ========================================================================== //
    // COMPUTE ADJACENCY MATRIX IF NOT ALREADY COMPUTED                           //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency(X);
    }

    // ========================================================================== //
    // FIND FREE SIMPLICIES                                                       //
    // ========================================================================== //

    // Loop over simplicies  ---------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        m = Simplex[i].size();
        for (j = 0; j < m; j++) {
            flag[i] = (flag[i] || (Adjacency[i][j][0] < 0));
        } //next j
    } //next i

    // Simplex list ------------------------------------------------------------- //
    list.resize(count(flag.begin(), flag.end(), true));
    j = 0;
    for (i = 0; i < nSimplex; i++) {
        if (flag[i]) {
            list[j] = i;
            j++;
        }
    }

    return(list); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindDoubleSimplex(
        void
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindDoubleSimplex(                                //
    //     void)                                                                  //
    //                                                                            //
    // Find douplicated simplex in the tasselation.                               //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, list of duplicated simplicies                        //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bool                    check;
    bvector1D               flag(nSimplex, false);
    ivector1D               idummy1D, doubles;
    ivector2D               index(nVertex);

    // Counters
    int                     V, T, S;
    int                     i, j, k, l, ii;
    int                     m, n, p;

    // ========================================================================== //
    // 1-RING OF VERTICES                                                         //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            V = Simplex[T][i];
            index[V].push_back(T);
        } //next i
    } //next T

    // ========================================================================== //
    // FIND DUPLICATED SIMPLICIES                                                 //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        if (!flag[T]) {
            flag[T] = true;
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                V = Simplex[T][i];
                n = index[V].size();
                for (j = 0; j < n; j++) {
                    S = index[V][j];
                    p = Simplex[S].size();
                    if ((p == m) && (!flag[S])) {
                        check = false;
                        idummy1D.resize(p);
                        k = 0;
                        while ((!check) && (k < p)) {
                            ii = k;
                            for (l = 0; l < p; l++) {
                                idummy1D[l] = Simplex[S][ii];
                                ii = (ii+1) % p;
                            } //next l
                            check = (Simplex[T] == idummy1D);
                            k++;
                        } //next circular shifting
                        if (check) {
                            doubles.push_back(S);
                            flag[S] = true;
                        }
                    }
                } //next j
            } //next i
        }
    } //next T

    return(doubles); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindDoubleSimplex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindDoubleSimplex(                                //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Find douplicated simplex in the tasselation. Vertex list is provided       //
    // externally.                                                                //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X      : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are    //
    //            the x, y, ... coordinates of the i-th vertex.                   //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, list of duplicated simplicies                        //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bool                    check;
    bvector1D               flag(nSimplex, false);
    ivector1D               idummy1D, doubles;
    ivector2D               index(X.size());

    // Counters
    int                     V, S, T;
    int                     i, j, k, l, ii;
    int                     m, n, p;

    // ========================================================================== //
    // 1-RING OF VERTICES                                                         //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            V = Simplex[T][i];
            index[V].push_back(T);
        } //next i
    } //next T

    // ========================================================================== //
    // FIND DUPLICATED SIMPLICIES                                                 //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        if (!flag[T]) {
            flag[T] = true;
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                V = Simplex[T][i];
                n = index[V].size();
                for (j = 0; j < n; j++) {
                    S = index[V][j];
                    p = Simplex[S].size();
                    if ((p == m) && (!flag[S])) {
                        check = false;
                        idummy1D.resize(p);
                        k = 0;
                        while ((!check) && (k < p)) {
                            ii = k;
                            for (l = 0; l < p; l++) {
                                idummy1D[l] = Simplex[S][ii];
                                ii = (ii+1) % p;
                            } //next l
                            check = (Simplex[T] == idummy1D);
                            k++;
                        } //next circular shifting
                        if (check) {
                            doubles.push_back(S);
                            flag[S] = true;
                        }
                    }
                } //next j
            } //next i
        }
    } //next T

    return(doubles); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindTrueDoubleSimplex(
        void
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindTrueDoubleSimplex(                            //
    //     void)                                                                  //
    //                                                                            //
    // Find douplicated simplex in the tasselation, indipendently from their      //
    // clockwise/counter-clockwise ordering                                       //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, list of duplicated simplicies                        //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bool                    check;
    bvector1D               flag(nSimplex, false);
    ivector1D               idummy1D,idummyC1D, doubles;
    ivector2D               index(nVertex);

    // Counters
    int                     V, T, S;
    int                     i, j, k, l, ii, jj;
    int                     m, n, p;

    // ========================================================================== //
    // 1-RING OF VERTICES                                                         //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            V = Simplex[T][i];
            index[V].push_back(T);
        } //next i
    } //next T

    // ========================================================================== //
    // FIND DUPLICATED SIMPLICIES                                                 //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        if (!flag[T]) {
            flag[T] = true;
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                V = Simplex[T][i];
                n = index[V].size();
                for (j = 0; j < n; j++) {
                    S = index[V][j];
                    p = Simplex[S].size();
                    if ((p == m) && (!flag[S])) {
                        check = false;
                        idummy1D.resize(p);
                        idummyC1D.resize(p);
                        k = 0;
                        while ((!check) && (k < p)) {
                            ii = k;
                            jj = k;
                            for (l = 0; l < p; l++) {
                                idummy1D[l]  = Simplex[S][ii];
                                idummyC1D[l] = Simplex[S][jj];
                                ii = (ii+1) % p;
                                jj = (p+jj-1)%p;
                            } //next l
                            check = ((Simplex[T] == idummy1D) || (Simplex[T] == idummyC1D));
                            k++;
                        } //next circular shifting
                        if (check) {
                            doubles.push_back(S);
                            flag[S] = true;
                        }
                    }
                } //next j
            } //next i
        }
    } //next T

    return(doubles); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::FindTrueDoubleSimplex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::FindTrueDoubleSimplex(                            //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Find douplicated simplex in the tasselation. indipendently from their      //
    // clockwise/counter-clockwise ordering.Vrtex list is provided       									//
    // externally.                                                                //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X      : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ... are    //
    //            the x, y, ... coordinates of the i-th vertex.                   //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list   : ivector1D, list of duplicated simplicies                        //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bool                    check;
    bvector1D               flag(nSimplex, false);
    ivector1D               idummy1D,idummyC1D,doubles;
    ivector2D               index(X.size());

    // Counters
    int                     V, S, T;
    int                     i, j, k, l, ii,jj;
    int                     m, n, p;

    // ========================================================================== //
    // 1-RING OF VERTICES                                                         //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        m = Simplex[T].size();
        for (i = 0; i < m; i++) {
            V = Simplex[T][i];
            index[V].push_back(T);
        } //next i
    } //next T

    // ========================================================================== //
    // FIND DUPLICATED SIMPLICIES                                                 //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        if (!flag[T]) {
            flag[T] = true;
            m = Simplex[T].size();
            for (i = 0; i < m; i++) {
                V = Simplex[T][i];
                n = index[V].size();
                for (j = 0; j < n; j++) {
                    S = index[V][j];
                    p = Simplex[S].size();
                    if ((p == m) && (!flag[S])) {
                        check = false;
                        idummy1D.resize(p);
                        idummyC1D.resize(p);	
                        k = 0;
                        while ((!check) && (k < p)) {
                            ii = k;
                            jj = k;
                            for (l = 0; l < p; l++) {
                                idummy1D[l] = Simplex[S][ii];
                                idummyC1D[l] = Simplex[S][jj];	
                                ii = (ii+1) % p;
                                jj = (p +jj-1) % p;	
                            } //next l
                            check = ((Simplex[T] == idummy1D) || (Simplex[T] == idummyC1D));
                            k++;
                        } //next circular shifting
                        if (check) {
                            doubles.push_back(S);
                            flag[S] = true;
                        }
                    }
                } //next j
            } //next i
        }
    } //next T

    return(doubles); };

// -------------------------------------------------------------------------- //
ivector1D Class_SurfTri::Find0AreaSimplex(
        void
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::Find0AreaSimplex(                                 //
    //     void)                                                                  //
    //                                                                            //
    // Find 0-area simplicies in the surface tasselation.                         //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - list     : ivector1D, with global index of 0-area simplicies             //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Parameters
    double const       toll = 1.0e-12;

    // Local variables
    double             A;
    ivector1D          list;

    // Counters
    int                T;

    // ========================================================================== //
    // LOOP OVER SIMPLICIES                                                       //
    // ========================================================================== //
    for (T = 0; T < nSimplex; T++) {
        A = Area(T);
        if (A < toll) {
            list.push_back(T);
        }
    } //next T

    return(list); };

// ------------------------------------------------------------------------- //
ivector1D Class_SurfTri::Find0AreaSimplex(
        dvecarr3E  &X
        ) {

    // ========================================================================= //
    // ivector1D Class_SurfTri::Find0AreaSimplex(                                //
    //      dvecarr3E  &X)                                                        //
    //                                                                           //
    // Find 0-area simplicies in the surface tasselation. Vertex list is         //
    // provided externally.                                                      //
    // ========================================================================= //
    // INPUT                                                                     //
    // ========================================================================= //
    // - X        : dvecarr3E, with vertex coordinates list. X[i][0], X[i][1],   //
    //              ... are the x, y, ... coordinates of the i-th vertex.        //
    // ========================================================================= //
    // OUTPUT                                                                    //
    // ========================================================================= //
    // - list     : ivector1D, with global index of 0-area simplicies            //
    // ========================================================================= //

    // ========================================================================= //
    // VARIABLES DECLARATION                                                     //
    // ========================================================================= //

    // Parameters
    double const       toll = 1.0e-12;

    // Local variables
    double             A;
    ivector1D          list;

    // Counters
    int                T;

    // ========================================================================= //
    // LOOP OVER SIMPLICIES                                                      //
    // ========================================================================= //
    for (T = 0; T < nSimplex; T++) {
        A = Area(T, X);
        if (A < toll) {
            list.push_back(T);
        }
    } //next T

    return(list); };

// ========================================================================== //
// CLEANING TOOLS                                                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveVertex(
        ivector1D   &list
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveVertex(                                          //
    //     ivector1D   &list)                                                     //
    //                                                                            //
    // Remove vertex from the tasselation                                         //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - list   : ivector1D, list with global index of removable node.            //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLE DECLARATION                                                       //
    // ========================================================================== //

    // Local variables
    ivector1D              map(nVertex, 0);

    // Counters
    int                    counter = 0;
    int                    i, j, m;

    // ========================================================================== //
    // REMOVE VERTEX FROM THE VERTEX LIST                                         //
    // ========================================================================== //

    // Initialize variables ----------------------------------------------------- //

    // vertex mapper
    for(i = 0; i < nVertex; i++) {
        map[i] = i;
    } //next i
    m = list.size();
    for (i = 0; i < m; i++) {
        map[list[i]] = -1;
    } //next i


    // Re-map node positions in the vertex list --------------------------------- //

    // Update vertex list
    j = 0;
    for (i = 0; i < nVertex; i++) {
        if (map[i] >= 0) {
            counter++;
            Vertex[j] = Vertex[i];
            map[i] = j;
            j++;
        }
    } //next i

    // Update vertex number
    nVertex = counter;

    // Update simplex-vertex connectivity
    for (i = 0; i < nSimplex; i++) {
        m = Simplex[i].size();
        for (j = 0; j < m; j++) {
            Simplex[i][j] = map[Simplex[i][j]];
        } //next j
    } //next i

    return; }

    // -------------------------------------------------------------------------- //
    void Class_SurfTri::RemoveVertex(
            ivector1D   &list,
            ivector1D   &map
            ) {

        // ========================================================================== //
        // void Class_SurfTri::RemoveVertex(                                          //
        //     ivector1D   &list,                                                     //
        //     ivector1D   &map)                                                      //
        //                                                                            //
        // Remove vertex from the tasselation, renumber tasselation vertexes and      //
        // returns a mapper between old and new numeration                            //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - list   : ivector1D, list with global indices of removable nodes.         //
        // - map    : ivector1D, map old -> new vertices. map[i] stores               //
        //            the global index of the i-th vertex after renumbering.          //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLE DECLARATION                                                       //
        // ========================================================================== //

        // Local variables
        // none

        // Counters
        int                    counter = 0;
        int                    i, j, m;

        // ========================================================================== //
        // RESIZE OUTPUT VARIABLES                                                    //
        // ========================================================================== //
        map.resize(nVertex, 0);

        // ========================================================================== //
        // REMOVE VERTEX FROM THE VERTEX LIST                                         //
        // ========================================================================== //

        // Initialize variables ----------------------------------------------------- //

        // vertex mapper
        for(i = 0; i < nVertex; i++) {
            map[i] = i;
        } //next i
        m = list.size();
        for (i = 0; i < m; i++) {
            map[list[i]] = -1;
        } //next i



        // Re-map node positions in the vertex list --------------------------------- //

        // Update vertex list
        j = 0;
        for (i = 0; i < nVertex; i++) {
            if (map[i] >= false) {
                counter++;
                Vertex[j] = Vertex[i];
                map[i] = j;
                j++;
            }
        } //next i

        // Update vertex number
        nVertex = counter;

        // Update simplex-vertex connectivity
        for (i = 0; i < nSimplex; i++) {
            m = Simplex[i].size();
            for (j = 0; j < m; j++) {
                Simplex[i][j] = map[Simplex[i][j]];
            } //next j
        } //next i

        return; }

        // -------------------------------------------------------------------------- //
        void Class_SurfTri::RemoveVertex(
                bvector1D   &flag
                ) {

            // ========================================================================== //
            // void Class_SurfTri::RemoveVertex(                                          //
            //     bvector1D   &flag)                                                     //
            //                                                                            //
            // Remove vertex from the tasselation                                         //
            // ========================================================================== //
            // INPUT                                                                      //
            // ========================================================================== //
            // - flag   : [nVertex-by-1] bvector1D, with flag for removable vertex.       //
            //            If flag[i] = true, then the i-th vertex is removed              //
            // ========================================================================== //
            // OUTPUT                                                                     //
            // ========================================================================== //
            // - none                                                                     //
            // ========================================================================== //

            // ========================================================================== //
            // VARIABLE DECLARATION                                                       //
            // ========================================================================== //

            // Local variables
            ivector1D              map(nVertex, 0);

            // Counters
            int                    counter = 0;
            int                    i, j, m;

            // ========================================================================== //
            // REMOVE VERTEX FROM THE VERTEX LIST                                         //
            // ========================================================================== //

            // Initialize variables ----------------------------------------------------- //

            // vertex mapper
            for(i = 0; i < nVertex; i++) {
                map[i] = i;
            } //next i


            // Re-map node positions in the vertex list --------------------------------- //

            // Update vertex list
            j = 0;
            for (i = 0; i < nVertex; i++) {
                if (flag[i] == false) {
                    counter++;
                    Vertex[j] = Vertex[i];
                    map[i] = j;
                    j++;
                }
            } //next i

            // Update vertex number
            nVertex = counter;

            // Update simplex-vertex connectivity
            for (i = 0; i < nSimplex; i++) {
                m = Simplex[i].size();
                for (j = 0; j < m; j++) {
                    Simplex[i][j] = map[Simplex[i][j]];
                } //next j
            } //next i

            return; }

            // -------------------------------------------------------------------------- //
            void Class_SurfTri::RemoveVertex(
                    bvector1D   &flag,
                    ivector1D   &map
                    ) {

                // ========================================================================== //
                // void Class_SurfTri::RemoveVertex(                                          //
                //     bvector1D   &flag,                                                     //
                //     ivector1D   &map)                                                      //
                //                                                                            //
                // Remove vertex from the tasselation, renumber tasselation vertexes and      //
                // returns mapper between old and new numeration                              //
                // ========================================================================== //
                // INPUT                                                                      //
                // ========================================================================== //
                // - flag   : [nVertex-by-1] bvector1D, with flag for removable vertex.       //
                //            If flag[i] = true, then the i-th vertex can be removed          //
                // - map    : [nVertex-by-1] ivector1D, map old -> new vertexes. map[i]       //
                //            stores the global index of the i-th vertex after renumbering.   //
                // ========================================================================== //
                // OUTPUT                                                                     //
                // ========================================================================== //
                // - none                                                                     //
                // ========================================================================== //

                // ========================================================================== //
                // VARIABLE DECLARATION                                                       //
                // ========================================================================== //

                // Local variables
                // none

                // Counters
                int                    counter = 0;
                int                    i, j, m;

                // ========================================================================== //
                // RESIZE OUTPUT VARIABLES                                                    //
                // ========================================================================== //
                map.resize(nVertex, 0);

                // ========================================================================== //
                // REMOVE VERTEX FROM THE VERTEX LIST                                         //
                // ========================================================================== //

                // Initialize variables ----------------------------------------------------- //

                // vertex mapper
                for(i = 0; i < nVertex; i++) {
                    map[i] = i;
                } //next i


                // Re-map node positions in the vertex list --------------------------------- //

                // Update vertex list
                j = 0;
                for (i = 0; i < nVertex; i++) {
                    if (flag[i] == false) {
                        counter++;
                        Vertex[j] = Vertex[i];
                        map[i] = j;
                        j++;
                    }
                } //next i

                // Update vertex number
                nVertex = counter;

                // Update simplex-vertex connectivity
                for (i = 0; i < nSimplex; i++) {
                    m = Simplex[i].size();
                    for (j = 0; j < m; j++) {
                        Simplex[i][j] = map[Simplex[i][j]];
                    } //next j
                } //next i

                return; }

                // -------------------------------------------------------------------------- //
                void Class_SurfTri::CollapseDoubleVertex(
                        ivector1D   &doublev,
                        int          n
                        ) {

                    // ========================================================================== //
                    // ivector1D Class_SurfTri::CollapseDoubleVertex(                             //
                    //     ivector1D   &doublev,                                                  //
                    //     int          n)                                                        //
                    //                                                                            //
                    // Collapse double vertices in the tasselation.                               //
                    // ========================================================================== //
                    // INPUT                                                                      //
                    // ========================================================================== //
                    // - doublev      : ivector1D, list with collapsed vertices                   //
                    // - n            : int (optional), number of bins for bin sorting            //
                    // ========================================================================== //
                    // OUTPUT                                                                     //
                    // ========================================================================== //
                    // - none                                                                     //
                    // ========================================================================== //

                    // ========================================================================== //
                    // VARIABLES DECLARATION                                                      //
                    // ========================================================================== //

                    // Variables declaration
                    int          m, ncell;
                    //ht int          dim = Vertex[0].size();
                    bvector1D    flag(nVertex, false);
                    ivector1D    index(nVertex, -1), idummy1D(2, -1);
                    ivector3D    cell;

                    // Counters
                    int          I_, C, S, V, W;
                    int          i;

                    // ========================================================================== //
                    // INITIALIZE VARIABLES                                                       //
                    // ========================================================================== //

                    // Random number generator -------------------------------------------------- //
                    //srand(time(NULL));
                    srand(32767);

                    // Resize variables --------------------------------------------------------- //
                    //ht if      (dim == 2) { cell.resize(n*n); }
                    //ht else if (dim == 3) { cell.resize(n*n*n); }

                    cell.resize(n*n*n); 

                    // List of collapsed vertices ----------------------------------------------- //
                    doublev.resize(0);

                    // ========================================================================== //
                    // SORT VERTICES ON BINS                                                      //
                    // ========================================================================== //

                    // Sort vertices ------------------------------------------------------------ //
                    BinSortV(index, n);

                    // Sort simplicies ---------------------------------------------------------- //
                    for (I_ = 0; I_ < nSimplex; I_++) {
                        m = Simplex[I_].size();
                        for (i = 0; i < m; i++) {
                            V = Simplex[I_][i];
                            C = index[V];
                            idummy1D[0] = I_;
                            idummy1D[1] = i;
                            cell[C].push_back(idummy1D);
                        } //next i
                    } //next I_

                    // ========================================================================== //
                    // COLLAPSE DOUBLE VERTICES                                                   //
                    // ========================================================================== //
                    ncell = cell.size();
                    //ht if (dim == 2) {
                    //ht     for (C = 0; C < ncell; C++) {
                    //ht         m = cell[C].size();
                    //ht         if (m > 0) {
                    //ht 
                    //ht             // Scope variables
                    //ht             ivector1D    list;
                    //ht             bitpit::KdTree<2, double, int>     kd(m);
                    //ht 
                    //ht             // Randomize vertex insertion
                    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
                    //ht             for (I_ = 0; I_ < m; I_++) {
                    //ht                 S = cell[C][list[I_]][0];
                    //ht                 i = cell[C][list[I_]][1];
                    //ht                 V = Simplex[S][i];
                    //ht                 if (kd.exist(&Vertex[V], W) >= 0) {
                    //ht                     Simplex[S][i] = W;
                    //ht                     if (!flag[V]) {
                    //ht                         flag[V] = true;
                    //ht                         doublev.push_back(V);
                    //ht                     }
                    //ht                 }
                    //ht                 else {
                    //ht                     flag[V] = true;
                    //ht                     kd.insert(&Vertex[V], V);
                    //ht                 }
                    //ht             } //next I_
                    //ht         }
                    //ht     } //next C
                    //ht }
                    //ht else if (dim == 3) {
                    //ht     for (C = 0; C < ncell; C++) {
                    //ht         m = cell[C].size();
                    //ht         if (m > 0) {
                    //ht 
                    //ht             // Scope variables
                    //ht             ivector1D    list;
                    //ht             bitpit::KdTree<3, double, int>     kd(m);
                    //ht 
                    //ht             // Randomize vertex insertion
                    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
                    //ht             for (I_ = 0; I_ < m; I_++) {
                    //ht                 S = cell[C][list[I_]][0];
                    //ht                 i = cell[C][list[I_]][1];
                    //ht                 V = Simplex[S][i];
                    //ht                 if (kd.exist(&Vertex[V], W) >= 0) {
                    //ht                     Simplex[S][i] = W;
                    //ht                     if (!flag[V]) {
                    //ht                         flag[V] = true;
                    //ht                         doublev.push_back(V);
                    //ht                     }
                    //ht                 }
                    //ht                 else {
                    //ht                     flag[V] = true;
                    //ht                     kd.insert(&Vertex[V], V);
                    //ht                 }
                    //ht             } //next I_
                    //ht         }
                    //ht     } //next C
                    //ht }

                    for (C = 0; C < ncell; C++) {
                        m = cell[C].size();
                        if (m > 0) {

                            // Scope variables
                            ivector1D    list;
                            bitpit::KdTree<3, array<double,3>, int>     kd(m);

                            // Randomize vertex insertion
                            bitpit::utils::extractWithoutReplacement(m, m-1, list);
                            for (I_ = 0; I_ < m; I_++) {
                                S = cell[C][list[I_]][0];
                                i = cell[C][list[I_]][1];
                                V = Simplex[S][i];
                                if (kd.exist(&Vertex[V], W) >= 0) {
                                    Simplex[S][i] = W;
                                    if (!flag[V]) {
                                        flag[V] = true;
                                        doublev.push_back(V);
                                    }
                                }
                                else {
                                    flag[V] = true;
                                    kd.insert(&Vertex[V], V);
                                }
                            } //next I_
                        }
                    } //next C

                    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::CollapseDoubleVertex(
        dvecarr3E   &X,
        ivector1D   &doublev,
        int          n
        ) {

    // ========================================================================== //
    // ivector1D Class_SurfTri::CollapseDoubleVertex(                             //
    //     dvecarr3E   &X,                                                        //
    //     ivector1D   &doublev,                                                  //
    //     int          n)                                                        //
    //                                                                            //
    // Collapse double vertices in the tasselation, using an external vertex      //
    // list.                                                                      //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X            : dvecarr3E, external vertex list. X[i][0], X[i][1], ...    //
    //                  are the x, y, ... components of the i-th vertex           //
    // - doublev      : ivector1D, list of collapsed vertices                     //
    // - n            : int (optional), number of bins for bin sorting            //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Variables declaration
    int          nV = X.size();
    int          m, ncell;
    //ht int          dim = X[0].size();
    bvector1D    flag(nV, false);
    ivector1D    index(nV, -1), idummy1D(2, -1);
    ivector3D    cell;

    // Counters
    int          I_, C, S, V, W;
    int          i;

    // ========================================================================== //
    // INITIALIZE VARIABLES                                                       //
    // ========================================================================== //

    // Random number generator -------------------------------------------------- //
    //srand(time(NULL));
    srand(32767);

    // Resize variables --------------------------------------------------------- //
    //ht if      (dim == 2) { cell.resize(n*n); }
    //ht else if (dim == 3) { cell.resize(n*n*n); }

    cell.resize(n*n*n);
    // ========================================================================== //
    // SORT VERTICES ON BINS                                                      //
    // ========================================================================== //

    // Sort vertices ------------------------------------------------------------ //
    BinSortV(X, index, n);

    // Sort simplicies ---------------------------------------------------------- //
    for (I_ = 0; I_ < nSimplex; I_++) {
        m = Simplex[I_].size();
        for (i = 0; i < m; i++) {
            V = Simplex[I_][i];
            C = index[V];
            idummy1D[0] = I_;
            idummy1D[1] = i;
            cell[C].push_back(idummy1D);
        } //next i
    } //next I_

    // ========================================================================== //
    // COLLAPSE DOUBLE VERTICES                                                   //
    // ========================================================================== //
    ncell = cell.size();
    //ht if (dim == 2) {
    //ht     for (C = 0; C < ncell; C++) {
    //ht         m = cell[C].size();
    //ht         if (m > 0) {
    //ht 
    //ht             // Scope variables
    //ht             ivector1D    list;
    //ht             bitpit::KdTree<2, double, int>     kd(m);
    //ht 
    //ht             // Randomize vertex insertion
    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
    //ht             for (I_ = 0; I_ < m; I_++) {
    //ht                 S = cell[C][list[I_]][0];
    //ht                 i = cell[C][list[I_]][1];
    //ht                 V = Simplex[S][i];
    //ht                 if (kd.exist(&X[V], W) >= 0) {
    //ht                     Simplex[S][i] = W;
    //ht                     doublev.push_back(V);
    //ht                 }
    //ht                 else {
    //ht                     kd.insert(&X[V], V);
    //ht                 }
    //ht             } //next I_
    //ht         }
    //ht     } //next C
    //ht }
    //ht else if (dim == 3) {
    //ht     for (C = 0; C < ncell; C++) {
    //ht         m = cell[C].size();
    //ht         if (m > 0) {
    //ht 
    //ht             // Scope variables
    //ht             ivector1D    list;
    //ht             bitpit::KdTree<3, double, int>     kd(m);
    //ht 
    //ht             // Randomize vertex insertion
    //ht             bitpit::utils::extractWithoutReplacement(m, m-1, list);
    //ht             for (I_ = 0; I_ < m; I_++) {
    //ht                 S = cell[C][list[I_]][0];
    //ht                 i = cell[C][list[I_]][1];
    //ht                 V = Simplex[S][i];
    //ht                 if (kd.exist(&X[V], W) >= 0) {
    //ht                     Simplex[S][i] = W;
    //ht                     if (!flag[V]) {
    //ht                         flag[V] = true;
    //ht                         doublev.push_back(V);
    //ht                     }
    //ht                 }
    //ht                 else {
    //ht                     flag[V] = true;
    //ht                     kd.insert(&X[V], V);
    //ht                 }
    //ht             } //next I_
    //ht         }
    //ht     } //next C
    //ht }

    for (C = 0; C < ncell; C++) {
        m = cell[C].size();
        if (m > 0) {

            // Scope variables
            ivector1D    list;
            bitpit::KdTree<3, array<double,3>, int>     kd(m);

            // Randomize vertex insertion
            bitpit::utils::extractWithoutReplacement(m, m-1, list);
            for (I_ = 0; I_ < m; I_++) {
                S = cell[C][list[I_]][0];
                i = cell[C][list[I_]][1];
                V = Simplex[S][i];
                if (kd.exist(&X[V], W) >= 0) {
                    Simplex[S][i] = W;
                    if (!flag[V]) {
                        flag[V] = true;
                        doublev.push_back(V);
                    }
                }
                else {
                    flag[V] = true;
                    kd.insert(&X[V], V);
                }
            } //next I_
        }
    } //next C


    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveSimplex(
        ivector1D   &list
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveSimplex(                                         //
    //     ivector1D   &list)                                                     //
    //                                                                            //
    // Remove simplicies from the simplex list.                                   //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - list   : ivector1D, list with global indices of removable simplicies.    //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector1D    map(nSimplex, 0);
    ivector2D    idummy2D;

    // Counter
    int          counter;
    int          I_, J;
    int          i, j, m;

    // Initialize variables ----------------------------------------------------- //
    for (i = 0; i < nSimplex; i++) {
        map[i] = i;
    } //next i
    m = list.size();
    for (i = 0; i < list.size(); i++) {
        map[list[i]] = -1;
    }

    // Loop over simplicies ----------------------------------------------------- //
    J = 0;
    counter = 0;
    for (I_ = 0; I_ < nSimplex; I_++) {
        if (map[I_] >= 0) {

            // Update counter
            counter++;

            // Update simplex-vertex connectivity
            Simplex[J].resize(Simplex[I_].size());
            Simplex[J] = Simplex[I_];

            // Update map
            map[I_] = J;

            // Update counter
            J++;

        }
    } //next I_

    // Update Normals ----------------------------------------------------------- //
    if ((Normal.size() > 0) && (Normal.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (map[I_] >= 0) {

                // Update normals
                Normal[map[I_]] = Normal[I_];
            }
        } //next I_
    }

    // Update Adjacency --------------------------------------------------------- //
    if ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (map[I_] >= 0) {

                // Update simplex-simplex connectivity
                idummy2D.resize(Adjacency[I_].size());
                m = Adjacency[I_].size();
                for (i = 0; i < m; i++) {
                    idummy2D[i].resize(Adjacency[I_][i].size());
                    idummy2D[i] = Adjacency[I_][i];
                } //next i

                Adjacency[map[I_]].resize(Adjacency[I_].size());
                for (i = 0; i < idummy2D.size(); i++) {
                    Adjacency[map[I_]][i].resize(idummy2D[i].size());
                    j = 0;
                    while (j < idummy2D[i].size()) {
                        J = idummy2D[i][j];
                        if (J >= 0) {
                            if (map[J] >= 0) {
                                Adjacency[map[I_]][i][j] = map[J];
                                j++;
                            }
                            else {
                                if (Adjacency[map[I_]][i].size() == 1) {
                                    Adjacency[map[I_]][i][j] = -1;
                                    j++;
                                }
                                else {
                                    Adjacency[map[I_]][i].erase(Adjacency[map[I_]][i].begin() + j);
                                    idummy2D[i].erase(idummy2D[i].begin() + j);
                                }
                            }
                        }
                        else {
                            Adjacency[map[I_]][i][j] = -1;
                            j++;
                        }
                    } //next j
                } //next i
            }
        } // next I_
    }

    // Update number of simplicies ---------------------------------------------- //
    nSimplex = counter;

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveSimplex(
        ivector1D   &list,
        ivector1D   &map
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveSimplex(                                         //
    //     ivector1D   &list,                                                     //
    //     ivector1D   &map)                                                      //
    //                                                                            //
    // Remove simplicies from the simplex list, renumber simplicies and returns   //
    // mapper between old and new numeration                                      //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - list  : ivector1D, list with global index of removable simplicies.       //
    // - map   : ivector1D, with old -> new map. map[i] contains the              //
    //           global index of the i-th simplex after renumbering               //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector2D    idummy2D;

    // Counter
    int          counter;
    int          I_, J;
    int          i, j, m;

    // Initialize variables ----------------------------------------------------- //
    map.resize(nSimplex,-0);
    for (i = 0; i < nSimplex; i++) {
        map[i] = i;
    }
    m = list.size();
    for (i = 0; i < m; i++) {
        map[list[i]] = -1;
    } //next i

    // Loop over simplicies ----------------------------------------------------- //
    J = 0;
    counter = 0;
    for (I_ = 0; I_ < nSimplex; I_++) {
        if (map[I_] >= 0) {

            // Update counter
            counter++;

            // Update simplex-vertex connectivity
            Simplex[J].resize(Simplex[I_].size());
            Simplex[J] = Simplex[I_];

            // Update map
            map[I_] = J;

            // Update counter
            J++;

        }
    } //next I_

    // Update Normals ----------------------------------------------------------- //
    if ((Normal.size() > 0) && (Normal.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (map[I_] >= 0) {

                // Update normals
                Normal[map[I_]] = Normal[I_];
            }
        } //next I_
    }

    // Update Adjacency --------------------------------------------------------- //
    if ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (map[I_] >= 0) {

                // Update simplex-simplex connectivity
                idummy2D.resize(Adjacency[I_].size());
                for (i = 0; i < Adjacency[I_].size(); i++) {
                    idummy2D[i].resize(Adjacency[I_][i].size());
                    idummy2D[i] = Adjacency[I_][i];
                } //next i

                Adjacency[map[I_]].resize(Adjacency[I_].size());
                for (i = 0; i < idummy2D.size(); i++) {
                    Adjacency[map[I_]][i].resize(idummy2D[i].size());
                    j = 0;
                    while (j < idummy2D[i].size()) {
                        J = idummy2D[i][j];
                        if (J >= 0) {
                            if (map[J] >= 0) {
                                Adjacency[map[I_]][i][j] = map[J];
                                j++;
                            }
                            else {
                                if (Adjacency[map[I_]][i].size() == 1) {
                                    Adjacency[map[I_]][i][j] = -1;
                                    j++;
                                }
                                else {
                                    Adjacency[map[I_]][i].erase(Adjacency[map[I_]][i].begin() + j);
                                    idummy2D[i].erase(idummy2D[i].begin() + j);
                                }
                            }
                        }
                        else {
                            Adjacency[map[I_]][i][j] = -1;
                            j++;
                        }
                    } //next j
                } //next i
            }
        } // next I_
    }

    // Update number of simplicies ---------------------------------------------- //
    nSimplex = counter;

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveSimplex(
        bvector1D   &flag
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveSimplex(                                         //
    //     bvector1D   &flag)                                                     //
    //                                                                            //
    // Remove simplicies from the simplex list.                                   //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - flag   : bvector1D, with flag removable simplicies.                      //
    //            If flag[i] = true, then the i-th simplex can be removed         //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector1D    map(nSimplex,-1);
    ivector2D    idummy2D;

    // Counter
    int         counter;
    int         I_, J;
    int         i, j;

    // Loop over simplicies ----------------------------------------------------- //
    J = 0;
    counter = 0;
    for (I_ = 0; I_ < nSimplex; I_++) {
        if (flag[I_] == false) {

            // Update counter
            counter += 1;

            // Update simplex-vertex connectivity
            Simplex[J].resize(Simplex[I_].size());
            Simplex[J] = Simplex[I_];

            // Update map
            map[I_] = J;

            // Update counter
            J++;

        }
    } //next I_

    // Update Normals ----------------------------------------------------------- //
    if ((Normal.size() > 0) && (Normal.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (flag[I_] == false) {

                // Update normals
                Normal[map[I_]] = Normal[I_];
            }
        } //next I_
    }

    // Update Adjacency --------------------------------------------------------- //
    if ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (flag[I_] == false) {

                // Update simplex-simplex connectivity
                idummy2D.resize(Adjacency[I_].size());
                for (i = 0; i < Adjacency[I_].size(); i++) {
                    idummy2D[i].resize(Adjacency[I_][i].size());
                    idummy2D[i] = Adjacency[I_][i];
                } //next i

                Adjacency[map[I_]].resize(Adjacency[I_].size());
                for (i = 0; i < idummy2D.size(); i++) {
                    Adjacency[map[I_]][i].resize(idummy2D[i].size());
                    j = 0;
                    while (j < idummy2D[i].size()) {
                        J = idummy2D[i][j];
                        if (J >= 0) {
                            if (map[J] >= 0) {
                                Adjacency[map[I_]][i][j] = map[J];
                                j++;
                            }
                            else {
                                if (Adjacency[map[I_]][i].size() == 1) {
                                    Adjacency[map[I_]][i][j] = -1;
                                    j++;
                                }
                                else {
                                    Adjacency[map[I_]][i].erase(Adjacency[map[I_]][i].begin() + j);
                                    idummy2D[i].erase(idummy2D[i].begin() + j);
                                }
                            }
                        }
                        else {
                            Adjacency[map[I_]][i][j] = -1;
                            j++;
                        }
                    } //next j
                } //next i
            }
        } // next I_
    }

    // Update number of simplicies ------------------------------------------------------- //
    nSimplex = counter;

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveSimplex(
        bvector1D   &flag,
        ivector1D   &map
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveSimplex(                                         //
    //     bvector1D   &flag,                                                     //
    //     ivector1D   &map)                                                      //
    //                                                                            //
    // Remove simplicies from the simplex list, renumber simplicies and returns   //
    // mapper between old and new numeration.                                     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - flag   : bvector1D, with flag removable simplicies.                      //
    //            If flag[i] = true, then the i-th simplex can be removed         //
    // - map   : ivector1D, with old -> new map. map[i] contains the              //
    //           global index of the i-th vertex after renumbering                //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector2D    idummy2D;

    // Counter
    int         counter;
    int         I_, J;
    int         i, j;

    // Initialize variables ----------------------------------------------------- //
    map.resize(nSimplex,-1);

    // Loop over simplicies ----------------------------------------------------- //
    J = 0;
    counter = 0;
    for (I_ = 0; I_ < nSimplex; I_++) {
        if (flag[I_] == false) {

            // Update counter
            counter += 1;

            // Update simplex-vertex connectivity
            Simplex[J].resize(Simplex[I_].size());
            Simplex[J] = Simplex[I_];

            // Update map
            map[I_] = J;

            // Update counter
            J++;

        }
    } //next I_

    // Update Normals -------------------------------------------------------------------- //
    if ((Normal.size() > 0) && (Normal.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (flag[I_] == false) {

                // Update normals
                Normal[map[I_]] = Normal[I_];
            }
        } //next I_
    }

    // Update Adjacency ------------------------------------------------------------------ //
    if ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex)) {
        for (I_ = 0; I_ < nSimplex; I_++) {
            if (flag[I_] == false) {

                // Update simplex-simplex connectivity
                idummy2D.resize(Adjacency[I_].size());
                for (i = 0; i < Adjacency[I_].size(); i++) {
                    idummy2D[i].resize(Adjacency[I_][i].size());
                    idummy2D[i] = Adjacency[I_][i];
                } //next i

                Adjacency[map[I_]].resize(Adjacency[I_].size());
                for (i = 0; i < idummy2D.size(); i++) {
                    Adjacency[map[I_]][i].resize(idummy2D[i].size());
                    j = 0;
                    while (j < idummy2D[i].size()) {
                        J = idummy2D[i][j];
                        if (J >= 0) {
                            if (map[J] >= 0) {
                                Adjacency[map[I_]][i][j] = map[J];
                                j++;
                            }
                            else {
                                if (Adjacency[map[I_]][i].size() == 1) {
                                    Adjacency[map[I_]][i][j] = -1;
                                    j++;
                                }
                                else {
                                    Adjacency[map[I_]][i].erase(Adjacency[map[I_]][i].begin() + j);
                                    idummy2D[i].erase(idummy2D[i].begin() + j);
                                }
                            }
                        }
                        else {
                            Adjacency[map[I_]][i][j] = -1;
                            j++;
                        }
                    } //next j
                } //next i
            }
        } // next I_
    }

    // Update number of simplicies ------------------------------------------------------- //
    nSimplex = counter;

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri:: RemoveIsolatedVertex(
        void
        ) {

    // ========================================================================== //
    // void Class_SurfTri:: RemoveIsolatedVertex(                                 //
    //     void)                                                                  //
    //                                                                            //
    // Remove isolated vertex from the tasselation. A vertex is isolated if and   //
    // only if thera are no simplicies in the tasselation having a vertex in that //
    // node.                                                                      //
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
    ivector1D      list;

    // Counters
    // none

    // ========================================================================== //
    // FIND ISOLATED VERTEX                                                       //
    // ========================================================================== //
    list = FindIsolatedVertex();

    // ========================================================================== //
    // REMOVE VERTEX                                                              //
    // ========================================================================== //
    RemoveVertex(list);

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri:: RemoveIsolatedVertex(
        ivector1D   &map
        ) {

    // ========================================================================== //
    // void Class_SurfTri:: RemoveIsolatedVertex(                                 //
    //     ivector1D   &map)                                                      //
    //                                                                            //
    // Remove isolated vertex from the tasselation, and returns map between old   //
    // and new numeration. A vertex is isolated if and only if thera are no       //
    // simplicies in the tasselation having a vertex in that node.                //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - map    : ivector1D, map old -> new vertexes. map[i] stores the           //
    //            global index of the i-th vertex after renumbering.              //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector1D      list;

    // Counters
    // none

    // ========================================================================== //
    // FIND ISOLATED VERTEX                                                       //
    // ========================================================================== //
    list = FindIsolatedVertex();

    // ========================================================================== //
    // REMOVE VERTEX                                                              //
    // ========================================================================== //
    RemoveVertex(list, map);

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveDoubleVertex(
        int          n
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveDoubleVertex(                                    //
    //     int          n)                                                        //
    //                                                                            //
    // Remove double vertex from the tasselation.                                 //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - n         : int (optional) number of bins for vertex sorting             //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector1D            list;

    // Counters
    // none

    // ========================================================================== //
    // COLLAPSE DOUBLE VERTEX                                                     //
    // ========================================================================== //
    CollapseDoubleVertex(list, n);

    // ========================================================================== //
    // REMOVE COLLAPSED VERTEX                                                    //
    // ========================================================================== //
    RemoveVertex(list);

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveDoubleVertex(
        ivector1D   &map,
        int          n
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveDoubleVertex(                                    //
    //     ivector1D   &map,                                                      //
    //     int          n)                                                        //
    //                                                                            //
    // Remove double vertex from the tasselation, and returns renumbering map     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - map    : ivector1D, map old -> new vertexes. map[i] stores the           //
    //            global index of the i-th vertex after renumbering.              //
    // - n      : int (optional), number of bins used for vertex sorting          //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector1D            list;

    // Counters
    // none

    // ========================================================================== //
    // COLLAPSE DOUBLE VERTEX                                                     //
    // ========================================================================== //
    CollapseDoubleVertex(list, n);

    // ========================================================================== //
    // REMOVE COLLAPSED VERTEX                                                    //
    // ========================================================================== //
    RemoveVertex(list, map);

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::RemoveIsolatedSimplex(
        void
        ) {

    // ========================================================================== //
    // void Class_SurfTri::RemoveIsolatedSimplex(                                 //
    //     void)                                                                  //
    //                                                                            //
    // Remove isolated simplicies from the tasselation. A isolated simplex is a   //
    // simplex whose vertices are not shared by any other simplex in the          //
    // tasselation.                                                               //
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
    ivector1D   list;

    // Counters
    // none

    // ========================================================================== //
    // REMOVE ISOLATED SIMPLICIES                                                 //
    // ========================================================================== //

    // Find isolated simplicies
    list = FindIsolatedSimplex();

    // Remove simplicies
    RemoveSimplex(list);

    return; }

    // -------------------------------------------------------------------------- //
    void Class_SurfTri::RemoveIsolatedSimplex(
            ivector1D   &map
            ) {

        // ========================================================================== //
        // void Class_SurfTri::RemoveIsolatedSimplex(                                 //
        //     ivector1D   &map)                                                      //
        //                                                                            //
        // Remove isolated simplicies from the tasselation, and returns map between   //
        // old and new numeration. A isolated simplex is a simplex whose vertex are   //
        // not shared by any other simplex in the tasselation.                        //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - map   : ivector1D, with old -> new map. map[i] contains the              //
        //           global index of the i-th vertex after renumbering                //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        ivector1D   list;

        // Counters
        // none

        // ========================================================================== //
        // REMOVE ISOLATED SIMPLICIES                                                 //
        // ========================================================================== //

        // Find isolated simplicies
        list = FindIsolatedSimplex();

        // Remove simplicies
        RemoveSimplex(list, map);

        return; }

        // -------------------------------------------------------------------------- //
        void Class_SurfTri::RemoveIsolatedSimplex(
                dvecarr3E   &X
                ) {

            // ========================================================================== //
            // void Class_SurfTri::RemoveIsolatedSimplex(                                 //
            //     dvecarr3E   &X)                                                        //
            //                                                                            //
            // Remove isolated simplicies from the tasselation. A isolated simplex is a   //
            // simplex whose vertex are not shared by any other simplex in the            //
            // tasselation. Vertex coordinate list is provided externally.                //
            // ========================================================================== //
            // INPUT                                                                      //
            // ========================================================================== //
            // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1],              //
            //          ... are the x, y, ... coordinates of the i-th node.               //
            // ========================================================================== //
            // OUTPUT                                                                     //
            // ========================================================================== //
            // - none                                                                     //
            // ========================================================================== //

            // ========================================================================== //
            // VARIABLES DECLARATION                                                      //
            // ========================================================================== //

            // Local variables
            ivector1D   list;

            // Counters
            // none

            // ========================================================================== //
            // REMOVE ISOLATED SIMPLICIES                                                 //
            // ========================================================================== //
            list = FindIsolatedSimplex(X);
            RemoveSimplex(list);

            return; }

            // -------------------------------------------------------------------------- //
            void Class_SurfTri::RemoveDoubleSimplex(
                    void
                    ) {

                // ========================================================================== //
                // void Class_SurfTri::RemoveDoubleSimplex(                                   //
                //     void)                                                                  //
                //                                                                            //
                // Remove duplicated simplicies from tasselation. A duplicated simplex is a   //
                // simplex having the same vertices of another simplex.                       //
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
                ivector1D    list;

                // Counters
                // none

                // ========================================================================== //
                // REMOVE DUPLICATED SIMPLICIES                                               //
                // ========================================================================== //

                // Find double simplicies
                list = FindDoubleSimplex();

                // Remove simplicies
                RemoveSimplex(list);

                return; }

                // -------------------------------------------------------------------------- //
                void Class_SurfTri::RemoveDoubleSimplex(
                        ivector1D   &map
                        ) {

                    // ========================================================================== //
                    // void Class_SurfTri::RemoveDoubleSimplex(                                   //
                    //     ivector1D   &map)                                                      //
                    //                                                                            //
                    // Remove duplicated simplicies from tasselation, and return renumbering map. //
                    // A duplicated simplex is a simplex having the same vertices of another      //
                    // simplex in the tasselation.                                                //
                    // ========================================================================== //
                    // INPUT                                                                      //
                    // ========================================================================== //
                    // - map   : ivector1D, map between old and new simplex numbering             //
                    // ========================================================================== //
                    // OUTPUT                                                                     //
                    // ========================================================================== //
                    // - none                                                                     //
                    // ========================================================================== //

                    // ========================================================================== //
                    // VARIABLES DECLARATION                                                      //
                    // ========================================================================== //

                    // Local variables
                    ivector1D    list;

                    // Counters
                    // none

                    // ========================================================================== //
                    // REMOVE DUPLICATED SIMPLICIES                                               //
                    // ========================================================================== //

                    // Find double simplicies
                    list = FindDoubleSimplex();

                    // Remove simplicies
                    RemoveSimplex(list, map);

                    return; }

                    // -------------------------------------------------------------------------- //
                    void Class_SurfTri::RemoveDoubleSimplex(
                            dvecarr3E   &X
                            ) {

                        // ========================================================================== //
                        // void Class_SurfTri::RemoveDoubleSimplex(                                   //
                        //     dvecarr3E   &X)                                                        //
                        //                                                                            //
                        // Remove duplicated simplicies from tasselation using an external vertex     //
                        // list. A duplicated simplex is a simplex having the same vertices of        //
                        // another simplex.                                                           //
                        // ========================================================================== //
                        // INPUT                                                                      //
                        // ========================================================================== //
                        // - X    : dvecarr3E, with vertex coordinate list. X[i][0], X[i][1], ... are //
                        //          the x, y, ... coordinates of the i-th node.                       //
                        // ========================================================================== //
                        // OUTPUT                                                                     //
                        // ========================================================================== //
                        // - none                                                                     //
                        // ========================================================================== //

                        // ========================================================================== //
                        // VARIABLES DECLARATION                                                      //
                        // ========================================================================== //

                        // Local variables
                        ivector1D    list;

                        // Counters
                        // none

                        // ========================================================================== //
                        // REMOVE DUPLICATED SIMPLICIES                                               //
                        // ========================================================================== //

                        // Find double simplicies
                        list = FindDoubleSimplex(X);

                        // Remove simplicies
                        RemoveSimplex(list);

                        return; }

                        // ----------------------------------------------------------------------------------- //
                        void Class_SurfTri::RemoveTrueDoubleSimplex(
                                void
                                ) {

                            // =================================================================================== //
                            // void Class_SurfTri::RemoveTrueDoubleSimplex(                                        //
                            //     void)                                                                           //
                            //                                                                                     //
                            // Remove duplicated simplicies from the tasselation. A duplicated simplex is a        //
                            // simplex whose vertex have the same coordinates (within a prescribed tolerance) of   //
                            // vertex of another simplex in the tasselation, without distinctions on vertex        //
                            // ordering.																																																																				       //	
                            // =================================================================================== //
                            // INPUT                                                                               //
                            // =================================================================================== //
                            // - none                                                                              //
                            // =================================================================================== //
                            // OUTPUT                                                                              //
                            // =================================================================================== //
                            // - none                                                                              //
                            // =================================================================================== //

                            // =================================================================================== //
                            // VARIABLES DECLARATION                                                               //
                            // =================================================================================== //

                            // Local variables
                            ivector1D    list;

                            // Counters
                            // none

                            // =================================================================================== //
                            // REMOVE DUPLICATED SIMPLICIES                                                        //
                            // =================================================================================== //

                            // Find double simplicies
                            list = FindTrueDoubleSimplex();

                            // Remove simplicies
                            RemoveSimplex(list);

                            return; }

                            // ----------------------------------------------------------------------------------- //
                            void Class_SurfTri::RemoveTrueDoubleSimplex(
                                    ivector1D &map
                                    ) {

                                // =================================================================================== //
                                // void Class_SurfTri::RemoveTrueDoubleSimplex(                                        //
                                //     ivector1D &map)                                                                 //
                                //                                                                                     //
                                // Remove duplicated simplicies from the tasselation, and return renumbering map.      //
                                // A duplicated simplex is a simplex whose vertex have the same coordinates (within a  //
                                // prescribed tolerance) of vertex of another simplex in the tasselation,              //
                                // without distinctions on vertex ordering				                                         //
                                // =================================================================================== //
                                // INPUT                                                                               //
                                // =================================================================================== //
                                // - map   : [nSimplex-by-1] ivector1D, with old -> new map. map[i] contains the       //
                                //           global index of the i-th vertex after renumbering                         //
                                // =================================================================================== //
                                // OUTPUT                                                                              //
                                // =================================================================================== //
                                // - none                                                                              //
                                // =================================================================================== //

                                // =================================================================================== //
                                // VARIABLES DECLARATION                                                               //
                                // =================================================================================== //

                                // Local variables
                                ivector1D    list;

                                // Counters
                                // none

                                // =================================================================================== //
                                // REMOVE DUPLICATED SIMPLICIES                                                        //
                                // =================================================================================== //

                                // Find double simplicies
                                list = FindTrueDoubleSimplex();

                                // Remove simplicies
                                RemoveSimplex(list, map);

                                return; }

                                // ----------------------------------------------------------------------------------- //
                                void Class_SurfTri::RemoveTrueDoubleSimplex(
                                        dvecarr3E &X
                                        ) {

                                    // =================================================================================== //
                                    // void Class_SurfTri::RemoveTrueDoubleSimplex(                                        //
                                    //     dvecarr3E &X)                                                                   //
                                    //                                                                                     //
                                    // Remove duplicated simplicies from the tasselation. A duplicated simplex is a        //
                                    // simplex whose vertex have the same coordinates (within a prescribed tolerance) of   //
                                    // vertex of another simplex in the tasselation.without distinctions on vertex	       // 
                                    // ordering. Vertex coordinate list is provided externally.                            //
                                    // =================================================================================== //
                                    // INPUT                                                                               //
                                    // =================================================================================== //
                                    // - X    : [nVertex-by-dim] dvecarr3E, with vertex coordinate list. X[i][0], X[i][1], //
                                    //          ... are the x, y, ... coordinates of the i-th node.                        //
                                    // =================================================================================== //
                                    // OUTPUT                                                                              //
                                    // =================================================================================== //
                                    // - none                                                                              //
                                    // =================================================================================== //

                                    // =================================================================================== //
                                    // VARIABLES DECLARATION                                                               //
                                    // =================================================================================== //

                                    // Local variables
                                    ivector1D    list;

                                    // Counters
                                    // none

                                    // =================================================================================== //
                                    // REMOVE DUPLICATED SIMPLICIES                                                        //
                                    // =================================================================================== //

                                    // Find double simplicies
                                    list = FindTrueDoubleSimplex(X);

                                    // Remove simplicies
                                    RemoveSimplex(list);

                                    return; }

                                    // -------------------------------------------------------------------------- //
                                    void Class_SurfTri::Remove0AreaSimplex(
                                            void
                                            ) {

                                        // ========================================================================== //
                                        // void Class_SurfTri::Remove0AreaSimplex(                                    //
                                        //     void)                                                                  //
                                        //                                                                            //
                                        // Remove simplex with 0-area.                                                //
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
                                        ivector1D               list;

                                        // Counters
                                        // none

                                        // ========================================================================== //
                                        // REMOVE SIMPLEX WITH 0-AREA                                                 //
                                        // ========================================================================== //

                                        // Find simplicies with 0-area ---------------------------------------------- //
                                        list = Find0AreaSimplex();

                                        // Remove simplicies with 0-area -------------------------------------------- //
                                        RemoveSimplex(list);

                                        return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Remove0AreaSimplex(
        ivector1D   &map
        ) {

    // ========================================================================== //
    // void Class_SurfTri::Remove0AreaSimplex(                                    //
    //     ivector1D   &map)                                                      //
    //                                                                            //
    // Remove simplex with 0-area and returns mapper between old and new          //
    // simplex numbering.                                                         //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - map      : ivector1D, map between old-->new numeration                   //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    ivector1D               list;

    // Counters
    // none

    // ========================================================================== //
    // REMOVE SIMPLEX WITH 0-AREA                                                 //
    // ========================================================================== //

    // Find simplicies with 0-area ---------------------------------------------- //
    list = Find0AreaSimplex();

    // Remove simplicies with 0-area -------------------------------------------- //
    RemoveSimplex(list, map);

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Remove0AreaSimplex(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // void Class_SurfTri::Remove0AreaSimplex(                                    //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Remove simplex with 0-area. Vertex list is provided externally.            //
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
    ivector1D               list;

    // Counters
    // none

    // ========================================================================== //
    // REMOVE SIMPLEX WITH 0-AREA                                                 //
    // ========================================================================== //

    // Find simplicies with 0-area ---------------------------------------------- //
    list = Find0AreaSimplex(X);

    // Remove simplicies with 0-area -------------------------------------------- //
    RemoveSimplex(list);

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Clean(
        void
        ) {

    // ========================================================================== //
    // void Class_SurfTri::Clean(                                                 //
    //     void)                                                                  //
    //                                                                            //
    // Clean tasselation from topology errors.                                    //
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
    bool            check_normal, check_adjacency;

    // Counters
    // none

    // ========================================================================== //
    // CLEAN TASSELATION FROM REPEATED VERTEX                                     //
    // ========================================================================== //

    // Set tollerance
    SetTolerance();

    // Veriables to be reshaped
    check_normal = ((Normal.size() > 0) && (Normal.size() >= nSimplex));
    check_adjacency = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

    // Clean
    RemoveDoubleVertex();
    Remove0AreaSimplex();
    RemoveIsolatedVertex();
    RemoveDoubleSimplex();

    // Resize
    ResizeVertex();
    ResizeSimplex();
    if (check_normal){
        ResizeNormal();
    }
    if (check_adjacency) {
        ResizeAdjacency();
    }

    return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Clean(
        dvecarr3E   &X
        ) {

    // ========================================================================== //
    // void Class_SurfTri::Clean(                                                 //
    //     dvecarr3E   &X)                                                        //
    //                                                                            //
    // Clean tasselation from repeated vertex and repeated simplex.               //
    // Vertex coordinate list is provided externally.                             //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - X    : dvecarr3E, vertex coordinate list. X[i][0], X[i][1],              //
    //          ... are the x, y, ... coordinates of the i-th node.               //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bool            check_normal, check_adjacency;
    ivector1D       list;

    // Counters
    // none

    // ========================================================================== //
    // CLEAN TASSELATION FROM REPEATED VERTEX                                     //
    // ========================================================================== //

    // Set tollerance
    SetTolerance(X);

    // Veriables to be reshaped
    check_normal = ((Normal.size() > 0) && (Normal.size() >= nSimplex));
    check_adjacency = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

    // Clean
    CollapseDoubleVertex(X, list);
    Remove0AreaSimplex(X);
    RemoveDoubleSimplex(X);

    // Resize
    ResizeVertex();
    ResizeSimplex();
    if (check_normal){
        ResizeNormal();
    }
    if (check_adjacency) {
        ResizeAdjacency();
    }

    return; };

// ========================================================================== //
// TASSELATION STATS                                                          //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::Stats(
        ostream     &out
        ) {

    // ========================================================================== //
    // void Class_SurfTri::Stats(                                                 //
    //     ostream     &out)                                                      //
    //                                                                            //
    // Compute tasselation stats.                                                 //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - out      : ostream, output stream                                        //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int             n;
    stringstream    msg;

    // Counters
    // none

    // ========================================================================== //
    // UPDATE ADJACENCY MATRIX                                                    //
    // ========================================================================== //
    if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
        BuildAdjacency();
    }

    // ========================================================================== //
    // VERTEX STATS                                                               //
    // ========================================================================== //
    out << "  Vertex ----------------------------------"   << endl;
    out << "    # vertex       " << nVertex                << endl;
    out << "    # isol. vertex " << CountIsolatedVertex()  << endl;
    out << "    # free  vertex " << CountFreeVertex()      << endl;
    out << "    # dupl. vertex " << CountDoubleVertex()    << endl;

    // ========================================================================== //
    // EDGE STATS                                                                 //
    // ========================================================================== //
    out << "  Edges -----------------------------------"   << endl;
    out << "    # edges        " << CountEdges()           << endl;
    out << "    # free  edges  " << CountFreeEdges()       << endl;

    // ========================================================================== //
    // SIMPLEX STATS                                                              //
    // ========================================================================== //
    out << "  Simplicies ------------------------------"   << endl;
    out << "    # simplicies   " << nSimplex               << endl;
    out << "    # isol. simpl. " << CountIsolatedSimplex() << endl;
    out << "    # free  simpl. " << CountFreeSimplex()     << endl;
    out << "    # dupl. simpl. " << CountDoubleSimplex()   << endl;

    return; }

    // -------------------------------------------------------------------------- //
    void Class_SurfTri::Stats(
            ostream     &out,
            dvecarr3E   &X
            ) {

        // ========================================================================== //
        // void Class_SurfTri::Stats(                                                 //
        //     ostream     &out,                                                      //
        //     dvecarr3E   &X)                                                        //
        //                                                                            //
        // Compute tasselation stats. Vertex coordinate list is provided externally   //
        // ========================================================================== //
        // INPUT                                                                      //
        // ========================================================================== //
        // - out      : ostream, output stream                                        //
        // - X        : dvecarr3E, with vertex coordinate list. X[i][0], X[i][1], ... //
        //              are the x, y, ... coordinates of the i-th node.               //
        // ========================================================================== //
        // OUTPUT                                                                     //
        // ========================================================================== //
        // - none                                                                     //
        // ========================================================================== //

        // ========================================================================== //
        // VARIABLES DECLARATION                                                      //
        // ========================================================================== //

        // Local variables
        int             n;
        stringstream    msg;

        // Counters
        // none

        // ========================================================================== //
        // UPDATE ADJACENCY MATRIX                                                    //
        // ========================================================================== //
        if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
            BuildAdjacency(X);
        }

        // ========================================================================== //
        // VERTEX STATS                                                               //
        // ========================================================================== //
        out << "  Vertex ----------------------------------"    << endl;
        out << "    # vertex       " << nVertex                 << endl;
        out << "    # isol. vertex " << CountIsolatedVertex(X)  << endl;
        out << "    # free  vertex " << CountFreeVertex(X)      << endl;
        out << "    # dupl. vertex " << CountDoubleVertex(X)    << endl;

        // ========================================================================== //
        // EDGE STATS                                                                 //
        // ========================================================================== //
        out << "  Edges -----------------------------------"    << endl;
        out << "    # edges        " << CountEdges()            << endl;
        out << "    # free  edges  " << CountFreeEdges()        << endl;

        // ========================================================================== //
        // SIMPLEX STATS                                                              //
        // ========================================================================== //
        out << "  Simplicies ------------------------------"    << endl;
        out << "    # simplicies   " << nSimplex                << endl;
        out << "    # isol. simpl. " << CountIsolatedSimplex(X) << endl;
        out << "    # free  simpl. " << CountFreeSimplex()      << endl;
        out << "    # dupl. simpl. " << CountDoubleSimplex(X)   << endl;

        return; }

