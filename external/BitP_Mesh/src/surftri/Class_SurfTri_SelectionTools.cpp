// ========================================================================== //
//                         - Class_SurfTri -                                  //
//                                                                            //
// Grid manager for unstructured meshes.                                      //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "Class_SurfTri.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
int Class_SurfTri::edge(
    int          I_,
    int          J
) {

// ========================================================================== //
// int Class_SurfTri::edge(                                                   //
//     int          I_,                                                        //
//     int          J)                                                        //
//                                                                            //
// Return the local index (on simplex I_) of the edge shared by simplicies I_   //
// and J. If no edge is found (i.e. if simplicies I_ and J do not share any    //
// edge), returns -1.                                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - I_     : int, global index of 1st simplex                                 //
// - J     : int, global index of 2nd simplex                                 //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index : int, local index (on simplex I_) of the edge shared by simplicies //
//           I_ and J. If no edge is found (i.e. simplicies I_ and J do not     //
//           share any edge), index = -1                                      //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool         flag = true;
int          index = -1;

// Counters
int          i, j, m, n;

// ========================================================================== //
// FIND EDGE                                                                  //
// ========================================================================== //
i = 0;
m = Simplex[I_].size();
while (flag && (i < m)) {
    j = 0;
    n = Adjacency[I_][i].size();
    while (flag && (j < n)) {
        if (Adjacency[I_][i][j] == J) {
            index = i;
            flag = false;
        }
        j++;
    } //next j
    i++;
} //next i

return(index); };

// -------------------------------------------------------------------------- //
int Class_SurfTri::edge(
    int          I_,
    int          V1,
    int          V2
) {

// ========================================================================== //
// int Class_SurfTri::edge(                                                   //
//     int          I_,                                                        //
//     int          V1,                                                       //
//     int          V2)                                                       //
//                                                                            //
// Return the local index (on simplex I_) of the edge having vertices V1       //
// and V2. If no edge is found, returns -1,                                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - I_      : int, simplex global index                                       //
// - V1, V2 : int, global indices of vertices of the edge                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index  : int, local index (on simplex) of the edge having vertices V1    //
//            and V2. If no edge is found index = -1                          //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool         flag = true;
int          index = -1;

// Counters
int          i, j, m = Simplex[I_].size();

// ========================================================================== //
// FIND EDGE                                                                  //
// ========================================================================== //
i = 0;
while (flag && (i < m)) {
    j = (i + 1) % m;
    if ((Simplex[I_][i] == V1) && (Simplex[I_][j] == V2)) {
            index = i;
            flag = false;
    }
    else if ((Simplex[I_][i] == V2) && (Simplex[I_][j] == V1)) {
            index = i;
            flag = false;
    }
    i++;
} //next i

return(index); };

// Vertices ================================================================= //

// -------------------------------------------------------------------------- //
int Class_SurfTri::vertex(
    int          I_,
    int          V1
) {

// ========================================================================== //
// int Class_SurfTri::vertex(                                                 //
//     int          I_,                                                        //
//     int          V1)                                                       //
//                                                                            //
// Return the local index (on simplex I_) of the vertex with global index V1.  //
// If no vertex is found, returns -1.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - I_        : int, simplex global index                                     //
// - V1       : int, vertex global index                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index    : int, local index of vertex on simplex I_. If no vertex is      //
//              found index = -1.                                             //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool     flag = false;
int      index = -1;

// Counters
int      i, m = Simplex[I_].size();

// ========================================================================== //
// FIND EDGE                                                                  //
// ========================================================================== //
i = 0;
while (!flag && (i < m)) {
    flag = (Simplex[I_][i] == V1);
    i++;
} //next i
index = i-1;

return(index); };

        // Rings ============================================================================= //

        // ----------------------------------------------------------------------------------- //
        ivector1D Class_SurfTri::Ring_1(
            int I_,
            int j,
            bool &flag
        ) {

        // =================================================================================== //
        // ivector1D Class_SurfTri::Ring_1(                                                    //
        //     int I_,                                                                          //
        //     int j,                                                                          //
        //     bool &flag)                                                                     //
        //                                                                                     //
        // Returns the list of simplicies in the 1-ring of the j-th vertex of simplex I_        //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - I_     : int, global index simplex which vertex belongs to.                        //
        // - j     : int, local index (on simplex I_) of vertex                                 //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - List  : ivector1D, list of simplicies in the 1-ring of vertex I_[i]                //
        // - flag  : bool, true if the 1-ring is closed, false otherwise                       //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        ivector1D     List;

        // Counters
        int           i, J, K;

        // =================================================================================== //
        // CREATE 1-RING LIST                                                                  //
        // =================================================================================== //

        // Find the starting simplex in the 1-ring ------------------------------------------- //

        // Set seed
        J = Adjacency[I_][j][0];
        K = I_;
        while ((J >= 0) && (J != I_)) {
            i = edge(J, K);
            j = (i + 1) % Simplex[J].size();
            K = J;
            J = Adjacency[J][j][0];
        } //next j
        I_ = K;

        // Compute 1-ring
        List.push_back(I_);
        i = (Simplex[I_].size() + j - 1) % Simplex[I_].size();
        J = Adjacency[I_][i][0];
        K = I_;
        while ((J >= 0) && (J != I_)) {
            List.push_back(J);
            j = edge(J, K);
            i = (Simplex[J].size() + j - 1) % Simplex[J].size();
            K = J;
            J = Adjacency[J][i][0];
        } //next J

        // Set flag for open/closed ring
        flag = (J >= 0);

        return(List); };

        // ----------------------------------------------------------------------------------- //
        ivector1D Class_SurfTri::Ring_1(
            int I_,
            int j,
            bool &flag,
            bool &isRing
        ) {

        // =================================================================================== //
        // ivector1D Class_SurfTri::Ring_1(                                                    //
        //     int I_,                                                                         //
        //     int j,                                                                          //
        //     bool &flag,								       //	
	//     bool &isRing)                                                                   //
        //                                                                                     //
        // Returns the list of simplicies in the 1-ring of the j-th vertex of simplex I_       //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - I_     : int, global index simplex which vertex belongs to.                       //
        // - j     : int, local index (on simplex I_) of vertex                                //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - List  : ivector1D, list of simplicies in the 1-ring of vertex I_[i]               //
        // - flag  : bool, true if the 1-ring is closed, false otherwise                       //
	// - isRing: bool, true if the 1-ring comp is successfull, false otherwise             //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        ivector1D     List, check, check_old;
						
        // Counters
        int           i, J, K,seedcheck, counter;

        // =================================================================================== //
        // CREATE 1-RING LIST                                                                  //
        // =================================================================================== //

        isRing = true;
        
        // Find the starting simplex in the 1-ring ------------------------------------------- //

        // Set seed
        J = Adjacency[I_][j][0];
        K = I_;
	seedcheck = J;
	counter=0;
    
	while ((J >= 0) && (J != I_)) {
            i = edge(J, K);
            j = (i + 1) % Simplex[J].size();
            K = J;
            J = Adjacency[J][j][0];
	    counter += (seedcheck == J );
	if(counter > 1){isRing= false; return(List);}											
        } //next j
        I_ = K;

        // Compute 1-ring
        List.push_back(I_);
        i = (Simplex[I_].size() + j - 1) % Simplex[I_].size();
        J = Adjacency[I_][i][0];
        K = I_;

	seedcheck = J;		
	counter=0;
        while ((J >= 0) && (J != I_)) {
            List.push_back(J);
            j = edge(J, K);
            i = (Simplex[J].size() + j - 1) % Simplex[J].size();
            K = J;
            J = Adjacency[J][i][0];

	   counter += ( seedcheck == J );
       	if(counter > 1){isRing= false; return(List);}											
        } //next J

        // Set flag for open/closed ring
        flag = (J >= 0);

        return(List); };



        // ----------------------------------------------------------------------------------- //
        ivector2D Class_SurfTri::VRing_1(
            int I_,
            int j,
            bool &flag
        ) {

        // =================================================================================== //
        // ivector2D Class_SurfTri::VRing_1(                                                   //
        //     int I_,                                                                          //
        //     int j,                                                                          //
        //     bool &flag)                                                                     //
        //                                                                                     //
        // Returns the list of vertices in the 1-ring of the j-th vertex of simplex I_.         //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - I_     : int, global index of simplex which the vertex belongs to.                 //
        // - j     : int, local index (on simplex I_) of vertex                                 //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - List  : ivector1D, list of vertices in the 1-ring of vertex I_[j]                  //
        // - flag  : bool, true if the 1-ring is closed, false otherwise                       //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        ivector1D           SList, idummy(2, -1);
        ivector2D           List;

        // Counters
        int                 i, k, J, V;

        // =================================================================================== //
        // 1-RING OF SIMPLEX                                                                   //
        // =================================================================================== //
        SList = Ring_1(I_, j, flag);

        // =================================================================================== //
        // COLLECT VERTICES IN THE 1-RING                                                      //
        // =================================================================================== //

        // Vertex global index --------------------------------------------------------------- //
        V = Simplex[I_][j];

        // Loop over simplicies -------------------------------------------------------------- //
        for (i = 0; i < SList.size(); i++) {
            k = vertex(SList[i], V);
            k = (k+1) % Simplex[SList[i]].size();
            idummy[0] = SList[i];
            idummy[1] = k;
            List.push_back(idummy);
        } //next i

        // Add last vertex in case of open ring ---------------------------------------------- //
        if (!flag) {
            i = SList.size()-1;
            k = vertex(SList[i], V);
            k = (k+2) % Simplex[SList[i]].size();
            idummy[0] = SList[i];
            idummy[1] = k;
            List.push_back(idummy);
        }

        return(List); };

// Bounding boxes =========================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::BoundingBox(
    dvector1D   &x_ext,
    dvector1D   &y_ext
) {

// ========================================================================== //
// void Class_SurfTri::BoundingBox(                                           //
//     dvector1D   &x_ext,                                                    //
//     dvector1D   &y_ext)                                                    //
//                                                                            //
// Compute limits of tasselation bounding box (2D case).                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - x_ext   : dvector1D, extent of bounding box in the x direction           //
// - y_ext   : dvector1D, extent of bounding box in the y direction           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         dim = Vertex[0].size();
int         n;

// Counters
int         i, j, I_;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Check data structure coherency ------------------------------------------- //
if (dim < 2) {
    return;
}

// Compute bounding box limits ---------------------------------------------- //
x_ext[0] = x_ext[1] = Vertex[Simplex[0][0]][0];
y_ext[0] = y_ext[1] = Vertex[Simplex[0][0]][1];
for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I_ = Simplex[i][j];
        x_ext[0] = min(x_ext[0], Vertex[I_][0]);
        x_ext[1] = max(x_ext[1], Vertex[I_][0]);
        y_ext[0] = min(y_ext[0], Vertex[I_][1]);
        y_ext[1] = max(y_ext[1], Vertex[I_][1]);
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BoundingBox(
    dvector1D   &x_ext,
    dvector1D   &y_ext,
    dvector1D   &z_ext
) {

// ========================================================================== //
// void Class_SurfTri::BoundingBox(                                           //
//     dvector1D   &x_ext,                                                    //
//     dvector1D   &y_ext,                                                    //
//     dvector1D   &z_ext)                                                    //
//                                                                            //
// Compute limits of tasselation bounding box (3D case).                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - x_ext   : dvector1D, extent of bounding box in the x direction           //
// - y_ext   : dvector1D, extent of bounding box in the y direction           //
// - z_ext   : dvector1D, extent of bounding box in the z direction           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         dim = Vertex[0].size();
int         n;

// Counters
int         i, j, I_;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Check data structure coherency ------------------------------------------- //
if (dim < 3) {
    return;
}

// Compute bounding box limits ---------------------------------------------- //
x_ext[0] = x_ext[1] = Vertex[Simplex[0][0]][0];
y_ext[0] = y_ext[1] = Vertex[Simplex[0][0]][1];
z_ext[0] = z_ext[1] = Vertex[Simplex[0][0]][2];
for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I_ = Simplex[i][j];
        x_ext[0] = min(x_ext[0], Vertex[I_][0]);
        x_ext[1] = max(x_ext[1], Vertex[I_][0]);
        y_ext[0] = min(y_ext[0], Vertex[I_][1]);
        y_ext[1] = max(y_ext[1], Vertex[I_][1]);
        z_ext[0] = min(z_ext[0], Vertex[I_][2]);
        z_ext[1] = max(z_ext[1], Vertex[I_][2]);
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BoundingBox(
    dvecarr3E   &V,
    dvector1D   &x_ext,
    dvector1D   &y_ext
) {

// ========================================================================== //
// void Class_SurfTri::BoundingBox(                                           //
//     dvecarr3E   &V,                                                        //
//     dvector1D   &x_ext,                                                    //
//     dvector1D   &y_ext)                                                    //
//                                                                            //
// Compute limits of tasselation bounding box (2D case). Vertex coordinate    //
// list is provided externally).                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V       : dvecarr3E vertex coordinate list. V[i][0], V[i][1], ... are    //
//             the x, y, ... coordinates of the i-th vertex.                  //
// - x_ext   : dvector1D, extent of bounding box in the x direction           //
// - y_ext   : dvector1D, extent of bounding box in the y direction           //
// - z_ext   : dvector1D, extent of bounding box in the z direction           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         dim = V[0].size();
int         n;

// Counters
int         i, j, I_;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Check data structure coherency ------------------------------------------- //
if (dim < 2) {
    return;
}

// Compute bounding box limits ---------------------------------------------- //
x_ext[0] = x_ext[1] = V[Simplex[0][0]][0];
y_ext[0] = y_ext[1] = V[Simplex[0][0]][1];
for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I_ = Simplex[i][j];
        x_ext[0] = min(x_ext[0], V[I_][0]);
        x_ext[1] = max(x_ext[1], V[I_][0]);
        y_ext[0] = min(y_ext[0], V[I_][1]);
        y_ext[1] = max(y_ext[1], V[I_][1]);
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BoundingBox(
    dvecarr3E   &V,
    dvector1D   &x_ext,
    dvector1D   &y_ext,
    dvector1D   &z_ext
) {

// ========================================================================== //
// void Class_SurfTri::BoundingBox(                                           //
//     dvecarr3E   &V,                                                        //
//     dvector1D   &x_ext,                                                    //
//     dvector1D   &y_ext,                                                    //
//     dvector1D   &z_ext)                                                    //
//                                                                            //
// Compute limits of tasselation bounding box (3D case). Vertex coordinate    //
// list is provided externally).                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V       : dvecarr3E vertex coordinate list. V[i][0], V[i][1], ... are    //
//             the x, y, ... coordinates of the i-th vertex.                  //
// - x_ext   : dvector1D, extent of bounding box in the x direction           //
// - y_ext   : dvector1D, extent of bounding box in the y direction           //
// - z_ext   : dvector1D, extent of bounding box in the z direction           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         dim = V[0].size();
int         n;

// Counters
int         i, j, I_;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Check data structure coherency ------------------------------------------- //
if (dim < 3) {
    return;
}

// Compute bounding box limits ---------------------------------------------- //
x_ext[0] = x_ext[1] = V[Simplex[0][0]][0];
y_ext[0] = y_ext[1] = V[Simplex[0][0]][1];
z_ext[0] = z_ext[1] = V[Simplex[0][0]][2];
for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I_ = Simplex[i][j];
        x_ext[0] = min(x_ext[0], V[I_][0]);
        x_ext[1] = max(x_ext[1], V[I_][0]);
        y_ext[0] = min(y_ext[0], V[I_][1]);
        y_ext[1] = max(y_ext[1], V[I_][1]);
        z_ext[0] = min(z_ext[0], V[I_][2]);
        z_ext[1] = max(z_ext[1], V[I_][2]);
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BoundingBox(
    darray3E   &B1,
    darray3E   &B2
) {


// Local variables
int         n;
darray3E    P;

// Counters
int         i, j, d, I_;

// ========================================================================== //
// COMPUTE BOUNDING BOX                                                       //
// ========================================================================== //

// Compute bounding box points ---------------------------------------------- //
B1 = B2 = Vertex[Simplex[0][0]];

for (i = 1; i < nSimplex; i++) {
    n = Simplex[i].size();
    for (j = 0; j < n; j++) {
        I_ = Simplex[i][j];

        P = Vertex[I_] ;

        for( d=0; d<3; ++d){
            B1[d] = min( B1[d], P[d] ) ;
            B2[d] = max( B2[d], P[d] ) ;
        };

    } //next j
} //next i

return; };

// ----------------------------------------------------------------------------------- //
void Class_SurfTri::Boundaries(
    Class_SurfTri      &Bounds
) {

// =================================================================================== //
// void Class_SurfTri::Boundaries(                                                     //
//     Class_SurfTri      &Bounds)                                                     //
//                                                                                     //
// Extract tasselation boundaries. A tessalation boundary is composed of free edges.   //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - Bounds      : Class_SurfTri, with tasselation boundaries.                         //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
ivector1D      simplex0D(1, -1);
ivector1D      simplex1D(2, -1);

// Counters
int            n, i, j, T;

// =================================================================================== //
// EXTRACT FREE EDGES IN THE TASSELATION                                               //
// =================================================================================== //
for (T = 0; T < nSimplex; T++) {
    n = Simplex[T].size();
    for (i = 0; i < n; i++) {
        if (Adjacency[T][i][0] < 0) {
            if (n == 1) {
            }
            else if (n == 2) {
                Bounds.AddVertex(Vertex[Simplex[T][i]]);
                simplex0D[0] = Bounds.nVertex-1;
                Bounds.AddSimplex(simplex0D);
            }
            else {
                j = (i+1) % n;
                Bounds.AddVertex(Vertex[Simplex[T][i]]);
                Bounds.AddVertex(Vertex[Simplex[T][j]]);
                simplex1D[0] = Bounds.nVertex-2;
                simplex1D[1] = Bounds.nVertex-1;
                Bounds.AddSimplex(simplex1D);
            }
        }
    } //next i
} //next T

// =================================================================================== //
// CLEAN BOUNDARIES                                                                    //
// =================================================================================== //
Bounds.RemoveDoubleVertex();

return; };

// ----------------------------------------------------------------------------------- //
void Class_SurfTri::Boundaries(
    dvecarr3E          &V,
    Class_SurfTri      &Bounds
) {

// =================================================================================== //
// void Class_SurfTri::Boundaries(                                                     //
//     dvecarr3E          &V,                                                          //
//     Class_SurfTri      &Bounds)                                                     //
//                                                                                     //
// Extract tasselation boundaries. A tessalation boundary is composed of free edges.   //
// Vertex list is provided externally.                                                 //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - V           : dvecarr3E, with vertex coordinate list. V[i][0], V[i][1], ...       //
//                 are the x, y, ... coordinates of the i-th vertex.                   //
// - Bounds      : Class_SurfTri, with tasselation boundaries.                         //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
ivector1D      simplex0D(1, -1);
ivector1D      simplex1D(2, -1);

// Counters
int            n, i, j, T;

// =================================================================================== //
// EXTRACT FREE EDGES IN THE TASSELATION                                               //
// =================================================================================== //
for (T = 0; T < nSimplex; T++) {
    n = Simplex[T].size();
    for (i = 0; i < n; i++) {
        if (Adjacency[T][i][0] < 0) {
            if (n == 1) {
            }
            else if (n == 2) {
                Bounds.AddVertex(V[Simplex[T][i]]);
                simplex0D[0] = Bounds.nVertex-1;
                Bounds.AddSimplex(simplex0D);
            }
            else {
                j = (i+1) % n;
                Bounds.AddVertex(V[Simplex[T][i]]);
                Bounds.AddVertex(V[Simplex[T][j]]);
                simplex1D[0] = Bounds.nVertex-2;
                simplex1D[1] = Bounds.nVertex-1;
                Bounds.AddSimplex(simplex1D);
            }
        }
    } //next i
} //next T

// =================================================================================== //
// CLEAN BOUNDARIES                                                                    //
// =================================================================================== //
Bounds.RemoveDoubleVertex();

return; };

// Search algorithms ================================================================= //

// ----------------------------------------------------------------------------------- //
// NOTE:                                                                               //
// - da ottimizzare                                                                    //
// - da generalizzare per il caso non manifold ed elementi di codimensione < d-1       //
// ----------------------------------------------------------------------------------- //
int Class_SurfTri::ReturnTriangleID(
    darray3E  &P,
    int        I_
) {

// =================================================================================== //
// int Class_SurfTri::ReturnTriangleID(                                                //
//     dvector1D &P,                                                                   //
//     int        I_)                                                                   //
//                                                                                     //
// Returns the global index of the simplex, which contains point P. The point P MUST   //
// lie on the surface tasselation.                                                     //
// (Lawson's search algorithm, adapted for surface tasselation).                       //
// If no simplex is found returns -1.                                                  //
// ** WARNING ****                                                                     //
// Under construction...                                                               //
// Works only on 3d manifold tasselation.                                              //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - P          : dvector1D, with vertex coordinates                                   //
// - I_          : int, global index of simplex used as seed in the search algorithm    //
//                If I_ < 0, seed is set to 0.                                          //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - index      : int, global index of the closes simplex to the seed containing       //
//                the projection of point P                                            //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Parameters
int                 max_iter = 5;
double              dtoll = 1.0e-4;

// Local variables
bool                check;
int                 n, dim = Vertex[0].size();
int                 index = -1;
double              n1, n2, max_theta;
bvector1D           visited(nSimplex, false);
ivector1D           i_dummy1D(2, -1);
ivector2D           backup;
darray3E            theta;
darray3E            v, v1, v2;

theta.fill(-1.0) ;

// Counters
int                 i, j, k, m, J, K, iter;

// =================================================================================== //
// INTIALIZE SEARCH PARAMETERS                                                         //
// =================================================================================== //

// Set seed -------------------------------------------------------------------------- //
J = max(I_, 0);

// Simplex normals ------------------------------------------------------------------- //
if (Normal.size() < nSimplex) {
    GenerateNormals();
}

// Adjacencies ----------------------------------------------------------------------- //
if (Adjacency.size() < nSimplex) {
   BuildAdjacency();
}

// =================================================================================== //
// SEARCH ALGORITHM                                                                    //
// =================================================================================== //
max_theta = -2.0;
iter = 0;
while (iter < max_iter) {

    // Update visited flag ----------------------------------------------------------- //
    visited[J] = true;

    // Check if P lies inside simplex J ---------------------------------------------- //
    n = Simplex[J].size();
    for (i = 0; i < n; i++) {

        // Check direction for the i-th edge
        j = (i+1) % n;
        // if (dim < 2) {
            // e = dotProduct(xP - Vertex[Simplex[J][i]],
                            // Vertex[Simplex[J][i]] - Vertex[Simplex[J][j]],
                            // dim);
            // check = (e >= 0.0);
        // }
        //else {
            v1 = P - Vertex[Simplex[J][i]];
            n1 = norm2(v1);
            if (n1 < 1.0e-12) {
                return(-1);
            }
            v1 = v1/n1;
            v2 = Vertex[Simplex[J][j]] - Vertex[Simplex[J][i]];
            n2 = norm2(v2);
            if (n2 < 1.0e-12) {
                return(-1);
            }
            v2 = v2/n2;
            v = crossProduct(v1, v2);
            v = v/norm2(v);
            theta[i] = dotProduct(Normal[J], v);

            // Update min theta
            max_theta = max(max_theta, theta[i]);
        //}
    } //next i

    // Simplex found ----------------------------------------------------------------- //
    if (max_theta < -1.0 + dtoll) {

        // Debug only (Export search path)
        // {
            // dvector2D    out(nSimplex, dvector1D(1, 0.0));
            // svector1D    nomi(1, "s_path");
            // for (i = 0; i < nSimplex; i++) {
                // if (visited[i] == true) {
                    // out[i][0] = 1.0;
                // }
            // }
            // ExportCellData_vtu("search.vtu",nomi,out);
        // }
        return(J);
    }

    // Next direction ---------------------------------------------------------------- //
    else {

        // Next candidate direction
        check = false;
        max_theta = -2.0;
        for (i = 0; i < n; i++) {
            K = Adjacency[J][i][0];
            if (K >= 0) {
                if (!visited[K]) {
                    if (theta[i] > max_theta) {
                        check = true;
                        max_theta = theta[i];
                        j = i;
                    }
                }
            }
        } // next i

        // no admissible directions
        if (!check) {
            J = Adjacency[backup[backup.size()-1][0]][backup[backup.size()-1][1]][0];
            backup.pop_back();
        }
        // admissible directon
        else {

            // save next direction
            i = j;

            // look for back-up direction
            k = (j+1) % n;
            m = 0;
            check = false;
            max_theta = -2.0;
            while (m < n-1) {
                K = Adjacency[J][k][0];
                if (K >= 0) {
                    if (!visited[K]) {
                        if (theta[k] > max_theta) {
                            check = true;
                            max_theta = theta[k];
                            j = k;
                        }
                    }
                }
                k = (k+1) % n;
                m++;
            } //next m
            if (check) {
                i_dummy1D[0] = J;
                i_dummy1D[1] = j;
                backup.push_back(i_dummy1D);
            }

            
            J = Adjacency[J][i][0];
        }
    }
        
    // Update iteration counter ------------------------------------------------------- //
    iter++;

} //next simplex;

return(index); };

// Patch decomposition =============================================================== //

 // ----------------------------------------------------------------------------------- //
 void Class_SurfTri::FindPatch_old(ivector1D &P,
                               double toll,
                               int index0,
                               int index,
                               int &seed) {

 // =================================================================================== //
 // void Class_SurfTri::FindPatch(ivector1D &P,                                         //
 //                               double toll,                                          //
 //                               int index0,                                           //
 //                               int index,                                            //
 //                               int &seed)                                            //
 //                                                                                     //
 // Find the surface patch of a surface patch, which a given triangle belongs to.       //
 // sub-patch is bounded either by a sharp edge or by a free edge (or both).            //
 // =================================================================================== //
 // INPUT                                                                               //
 // =================================================================================== //
 // - toll    : double, tolerance for sharp edge deterction                             //
 // - index0  : int, index of patch to be scanned                                       //
 // - index   : int, index to be assigned to the new patch                              //
 // - seed    : int, index of simplex chosen as seed                                    //
 // =================================================================================== //
 // OUTPUT                                                                              //
 // =================================================================================== //
 // - P       : [Simplex-by-1] ivector1D, with patch index. P[i] is the patch index     //
 //             which the i-th simplex belongs to                                       //
 // =================================================================================== //

 // =================================================================================== //
 // VARIABLES DECLARATION                                                               //
 // =================================================================================== //

 // Local variables
 bool         stop, stop1;
 bool         free_edge, visited, sharp_edge, T_edge;
 ivector1D    e(2, -1), dummy(2, -1);
 ivector2D    n(0, ivector1D(2,-1));

 // Counters
 int          iter;
 int          i, j, k;
 int          I_, J;

 // =================================================================================== //
 // BUILD ADJACENCY MATRIX IF NOT ALREADY COMPUTED                                      //
 // =================================================================================== //
 if ((Adjacency.size() == 0) | (Adjacency.size() != nSimplex)) {
     BuildAdjacency();
 }

 // =================================================================================== //
 // GENERATE NORMALS IF NOT AVAILABLE                                                   //
 // =================================================================================== //
 if ((Normal.size() == 0) | (Normal.size() != nSimplex)) {
     GenerateNormals();
 }

 // =================================================================================== //
 // SET SEED                                                                            //
 // =================================================================================== //

 // Advancing direction
 e[0] = seed;
 e[1] = -1;
 stop1 = false;
 i = 0;
 while ((stop1 == false) && (i < Simplex[e[0]].size())) {
     T_edge = (Adjacency[e[0]][i].size() > 1);
     if (!T_edge) {
         J = Adjacency[e[0]][i][0];
         free_edge = (J < 0);
         if (!free_edge) {
             visited = !(P[J] == index0);
             sharp_edge = (norm2(Normal[e[0]] - Normal[J]) > toll);
             if ((!visited) && (!sharp_edge)) {
                 e[1] = i;
                 stop1 = true;
             }
         }
     }
     i++;
 } //next i

 // Set flag for the seed element
 P[e[0]] = index;

 // =================================================================================== //
 // LOOP UNTIL ADVANCING LIST IS EMPTY                                                  //
 // =================================================================================== //
 iter = 0;
 stop = (e[1] < 0);
 while (stop == false) {

     // Alternative direction
     I_ = e[0];
     i = (e[1] + 1) % Simplex[I_].size();
     stop1 = false;
     while (stop1 == false) {
         T_edge = (Adjacency[I_][i].size() > 1);
         if (!T_edge) {
             J = Adjacency[I_][i][0];
             free_edge = (J < 0);
             if (!free_edge) {
                 visited = !(P[J] == index0);
                 sharp_edge = (norm2(Normal[I_] - Normal[J]) > toll);
                 if ((!visited) && (!sharp_edge)) {
                     dummy[0] = I_;
                     dummy[1] = i;
                     n.push_back(dummy);
                     stop1 = true;
                 }
             }
          }
          i = (i+1) % Simplex[I_].size();
          stop1 = (stop1 || (i == e[1]));
     } //next i

     // Move to next element
     I_ = Adjacency[e[0]][e[1]][0];
     P[I_] = index;
     seed = I_;

     // Next propagation direction
     dummy[0] = I_;
     dummy[1] = -1;
     i = 0;
     stop1 = false;
     while (stop1 == false) {
         T_edge = (Adjacency[I_][i].size() > 1);
         if (!T_edge) {
             J = Adjacency[I_][i][0];
             free_edge = (J < 0);
             if (!free_edge) {
                 visited = !(P[J] == index0);
                 sharp_edge = (norm2(Normal[I_] - Normal[J]) > toll);
                 if ((!visited) && (!sharp_edge)) {
                     dummy[1] = i;
                     stop1 = true;
                 }
             }
         }
         i = (i+1) % Simplex[I_].size();
         stop1 = (stop1 || (i == 0));
     } //next i

     // Takes the alternative direction
     if (dummy[1] < 0) {
         if (n.size() > 0) {
             e = n[n.size()-1];
             n.pop_back();
             }
         else {
             e[0] = -1;
             e[1] = -1;
         }
     }
     else {
         e = dummy;
     }

     // Stopping criterion
     stop = ((n.size() == 0) && (e[0] < 0));

     // iteration counter
     iter += 1;

 } // next iter

 return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::FindPatch(
    ivector1D   &P,
    double       toll,
    int          index0,
    int          index,
    int         &T
) {

// ========================================================================== //
// void Class_SurfTri::FindPatch(                                             //
//     ivector1D   &P,                                                        //
//     double       toll,                                                     //
//     int          index0,                                                   //
//     int          index,                                                    //
//     int         &T)                                                        //
//                                                                            //
// Flag all simplicies belonging to the same surface patch starting from      //
// a given seed. A surface patch is bounded either by a sharp edge or by a    //
// T junction or by a free edge.                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P       : ivector1D, with patch index. P[i] is the patch index           //
//             which the i-th simplex belongs to                              //
// - toll    : double, tolerance for sharp edge deterction                    //
// - index0  : int, index of which seed belongs to                            //
// - index   : int, index to be assigned to the new patch                     //
// - T       : int, global index of simplex chosen as seed                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double                  sharp_e;
bitpit::LIFOStack<int>          stack;

// Counters
int                     i, j;
int                     S;
int                     m;

// ========================================================================== //
// BUILD ADJACENCY MATRIX IF NOT ALREADY BUILT                                //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// GENERATE NORMALS IF NOT ALREADY COMPUTED                                   //
// ========================================================================== //
if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
    GenerateNormals();
}

// ========================================================================== //
// SET SEED                                                                   //
// ========================================================================== //
stack.push(T);

// ========================================================================== //
// PROPAGATE THROUGH NEIGHBORS                                                //
// ========================================================================== //
while (stack.TOPSTK > 0) {

    // Next simplex to be processed
    T = stack.pop();
    P[T] = index;

    // Loop over neighbors
    m = Simplex[T].size();
    for (i = 0; i < m; i++) {
        if (Adjacency[T][i][0] >= 0) {
            if (Adjacency[T][i].size() == 1) {
                S = Adjacency[T][i][0];
                if (P[S] == index0) {
                    sharp_e = norm2(Normal[T] - Normal[S]);
                    if (sharp_e < toll) {
                        stack.push(S);
                    }
                }
            }
        }
    } //next i
    
} //next item

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::FindConnectedLoop(
    ivector1D   &P,
    int          index0,
    int          index,
    int         &T
) {

// ========================================================================== //
// void Class_SurfTri::FindConnectedLoop(                                     //
//     ivector1D   &P,                                                        //
//     int          index0,                                                   //
//     int          index,                                                    //
//     int         &T)                                                        //
//                                                                            //
// Flag all simplicies belonging to the same connected loop starting from     //
// a given seed.                                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P       : ivector1D, with patch index. P[i] is the patch index           //
//             which the i-th simplex belongs to                              //
// - index0  : int, index of which seed belongs to                            //
// - index   : int, index to be assigned to the new patch                     //
// - T       : int, global index of simplex chosen as seed                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bitpit::LIFOStack<int>          stack;

// Counters
int                     i, j;
int                     S;
int                     m, n;

// ========================================================================== //
// BUILD ADJACENCY MATRIX IF NOT ALREADY BUILT                                //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() < nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// SET SEED                                                                   //
// ========================================================================== //
stack.push(T);

// ========================================================================== //
// PROPAGATE THROUGH NEIGHBORS                                                //
// ========================================================================== //
while (stack.TOPSTK > 0) {

    // Next simplex to be processed
    T = stack.pop();
    P[T] = index;

    // Loop over neighbors
    m = Simplex[T].size();
    for (i = 0; i < m; i++) {
        if (Adjacency[T][i][0] >= 0) {
            n = Adjacency[T][i].size();
            for (j = 0; j < n; j++) {
                S = Adjacency[T][i][j];
                if (P[S] == index0) {
                    stack.push(S);
                }
            } //next j
        }
    } //next i
    
} //next item

return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::FindClosedLoop(ivector1D &P,
                                           int index0,
                                           int index,
                                           int &seed) {

        // =================================================================================== //
        // void Class_SurfTri::FindClosedLoop(ivector1D &P,                                    //
        //                                    int index0,                                      //
        //                                    int index,                                       //
        //                                    int &seed)                                       //
        //                                                                                     //
        // Find closed loops in a surface tasselation.                                         //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - index0  : int, index of loop to be scanned                                        //
        // - index   : int, index to be assigned to the new loop                               //
        // - seed    : int, index of simplex chosen as seed                                    //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - P       : [Simplex-by-1] ivector1D, with loops index. P[i] is the loop index      //
        //             which the i-th simplex belongs to                                       //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool         stop, stop1;
        bool         free_edge, visited, T_edge;
        ivector1D    e(2, -1), dummy(2, -1);
        ivector2D    n(0, ivector1D(2,-1));

        // Counters
        int          iter;
        int          i, j, k;
        int          I_, J;

        // =================================================================================== //
        // BUILD ADJACENCY MATRIX IF NOT ALREADY COMPUTED                                      //
        // =================================================================================== //
        if (Adjacency.size() < nSimplex) {
            BuildAdjacency();
        }

        // =================================================================================== //
        // SET SEED                                                                            //
        // =================================================================================== //

        // Advancing direction
        e[0] = seed;
        e[1] = -1;
        stop1 = false;
        i = 0;
        while ((stop1 == false) && (i < Simplex[e[0]].size())) {
            T_edge = (Adjacency[e[0]][i].size() > 1);
            if (!T_edge) {
                J = Adjacency[e[0]][i][0];
                free_edge = (J < 0);
                if (!free_edge) {
                    visited = !(P[J] == index0);
                    if (!visited) {
                        e[1] = i;
                        stop1 = true;
                    }
                }
            }
            i++;
        } //next i

        // Set flag for the seed element
        P[e[0]] = index;

        // =================================================================================== //
        // LOOP UNTIL ADVANCING LIST IS EMPTY                                                  //
        // =================================================================================== //
        iter = 0;
        stop = (e[1] < 0);
        while (stop == false) {

            // Alternative direction
            I_ = e[0];
            i = (e[1] + 1) % Simplex[I_].size();
            stop1 = false;
            while (stop1 == false) {
                T_edge = (Adjacency[I_][i].size() > 1);
                if (!T_edge) {
                    J = Adjacency[I_][i][0];
                    free_edge = (J < 0);
                    if (!free_edge) {
                        visited = !(P[J] == index0);
                        if (!visited) {
                            dummy[0] = I_;
                            dummy[1] = i;
                            n.push_back(dummy);
                            stop1 = true;
                        }
                    }
                 }
                 i = (i+1) % Simplex[I_].size();
                 stop1 = (stop1 || (i == e[1]));
            } //next i

            // Move to next element
            I_ = Adjacency[e[0]][e[1]][0];
            P[I_] = index;
            seed = I_;

            // Next propagation direction
            dummy[0] = I_;
            dummy[1] = -1;
            i = 0;
            stop1 = false;
            while (stop1 == false) {
                T_edge = (Adjacency[I_][i].size() > 1);
                if (!T_edge) {
                    J = Adjacency[I_][i][0];
                    free_edge = (J < 0);
                    if (!free_edge) {
                        visited = !(P[J] == index0);
                        if (!visited) {
                            dummy[1] = i;
                            stop1 = true;
                        }
                    }
                }
                i = (i+1) % Simplex[I_].size();
                stop1 = (stop1 || (i == 0));
            } //next i

            // Takes the alternative direction
            if (dummy[1] < 0) {
                if (n.size() > 0) {
                    e = n[n.size()-1];
                    n.pop_back();
                    }
                else {
                    e[0] = -1;
                    e[1] = -1;
                }
            }
            else {
                e = dummy;
            }

            // Stopping criterion
            stop = ((n.size() == 0) && (e[0] < 0));

            // iteration counter
            iter += 1;

        } // next iter

        return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::FindRegion(ivector1D &P,
                                       ivector1D &flag,
                                       int index0,
                                       int index,
                                       int &seed) {

        // =================================================================================== //
        // void Class_SurfTri::FindRegion(ivector1D &P,                                        //
        //                                ivector1D &flag,                                     //
        //                                int index0,                                          //
        //                                int index,                                           //
        //                                int &seed)                                           //
        //                                                                                     //
        // Find delimited region in a surface tasselation                                      //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - P       : [Simplex-by-1] ivector1D, with ID of each sub-region associated to each //
        //             simplex.                                                                //
        //             which the i-th simplex belongs to                                       //
        // - index0  : int, index of loop to be scanned                                        //
        // - index   : int, index to be assigned to the new loop                               //
        // - seed    : int, index of simplex chosen as seed                                    //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - P       : [Simplex-by-1] ivector1D, with loops index. P[i] is the loop index      //
        //             which the i-th simplex belongs to                                       //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool         stop, stop1;
        bool         free_edge, visited, T_edge, region;
        ivector1D    e(2, -1), dummy(2, -1);
        ivector2D    n(0, ivector1D(2,-1));

        // Counters
        int          iter;
        int          i, j, k;
        int          I_, J;

        // =================================================================================== //
        // BUILD ADJACENCY MATRIX IF NOT ALREADY COMPUTED                                      //
        // =================================================================================== //
        if (Adjacency.size() < nSimplex) {
            BuildAdjacency();
        }

        // =================================================================================== //
        // SET SEED                                                                            //
        // =================================================================================== //

        // Advancing direction
        e[0] = seed;
        e[1] = -1;
        stop1 = false;
        i = 0;
        while ((stop1 == false) && (i < Simplex[e[0]].size())) {
            T_edge = (Adjacency[e[0]][i].size() > 1);
            if (!T_edge) {
                J = Adjacency[e[0]][i][0];
                free_edge = (J < 0);
                if (!free_edge) {
                    visited = !(P[J] == index0);
                    if (!visited) {
                        region = (flag[e[0]] != flag[J]);
                        if (!region) {
                            e[1] = i;
                            stop1 = true;
                        }
                    }
                }
            }
            i++;
        } //next i

        // Set flag for the seed element
        P[e[0]] = index;

        // =================================================================================== //
        // LOOP UNTIL ADVANCING LIST IS EMPTY                                                  //
        // =================================================================================== //
        iter = 0;
        stop = (e[1] < 0);
        while (stop == false) {

            // Alternative direction
            I_ = e[0];
            i = (e[1] + 1) % Simplex[I_].size();
            stop1 = false;
            while (stop1 == false) {
                T_edge = (Adjacency[I_][i].size() > 1);
                if (!T_edge) {
                    J = Adjacency[I_][i][0];
                    free_edge = (J < 0);
                    if (!free_edge) {
                        visited = !(P[J] == index0);
                        if (!visited) {
                            region = (flag[I_] != flag[J]);
                            if (!region) {
                                dummy[0] = I_;
                                dummy[1] = i;
                                n.push_back(dummy);
                                stop1 = true;
                            }
                        }
                    }
                 }
                 i = (i+1) % Simplex[I_].size();
                 stop1 = (stop1 || (i == e[1]));
            } //next i

            // Move to next element
            I_ = Adjacency[e[0]][e[1]][0];
            P[I_] = index;
            seed = I_;

            // Next propagation direction
            dummy[0] = I_;
            dummy[1] = -1;
            i = 0;
            stop1 = false;
            while (stop1 == false) {
                T_edge = (Adjacency[I_][i].size() > 1);
                if (!T_edge) {
                    J = Adjacency[I_][i][0];
                    free_edge = (J < 0);
                    if (!free_edge) {
                        visited = !(P[J] == index0);
                        if (!visited) {
                            region = (flag[I_] != flag[J]);
                            if (!region) {
                                dummy[1] = i;
                                stop1 = true;
                            }
                        }
                    }
                }
                i = (i+1) % Simplex[I_].size();
                stop1 = (stop1 || (i == 0));
            } //next i

            // Takes the alternative direction
            if (dummy[1] < 0) {
                if (n.size() > 0) {
                    e = n[n.size()-1];
                    n.pop_back();
                    }
                else {
                    e[0] = -1;
                    e[1] = -1;
                }
            }
            else {
                e = dummy;
            }

            // Stopping criterion
            stop = ((n.size() == 0) && (e[0] < 0));

            // iteration counter
            iter += 1;

        } // next iter

        return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::PatchDecompose(int        &nPatches,
                                           ivector1D  &Patches,
                                           int         index0,
                                           double      toll) {

        // =================================================================================== //
        // void Class_SurfTri::PatchDecompose(int        &nPatches,                            //
        //                                    ivector1D  &Patches,                             //
        //                                    int         index0)                              //
        //                                    double      toll)                                //
        //                                                                                     //
        // Decompose a surface patch into sub-patches. Each patch is bounded either by sharp   //
        // edges, or by free edges.                                                            //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - nPatches  : int, number of patches generated                                      //
        // - Patches   : ivector1D, with patch flag. Patch[i] is the patch ID                  //
        //               which the i-th simplex belongs to.                                    //
        // - index0    : int, patch ID to be decomposed                                        //
        // - toll      : double, tollerance for sharp edge detection                           //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool                   flag, stop;

        // Counters
        int                    counter;
        int                    i, j;
        int                    I_, J;

        // =================================================================================== //
        // CHECK FOR EMPTY TRIANGULATION                                                       //
        // =================================================================================== //
        if (nSimplex == 0) {
            return;
        }
        
        // =================================================================================== //
        // COMPUTE ADJACENCY MATRIX IF NOT ALREADY BUILT                                       //
        // =================================================================================== //
        if ((Adjacency.size() == 0) || (Adjacency.size() != nSimplex)) {
            BuildAdjacency();
        }

        // =================================================================================== //
        // GENERATE NORMALS IF NOT ALREADY COMPUTED                                            //
        // =================================================================================== //
        if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
            GenerateNormals();
        }

        // =================================================================================== //
        // DECOMPOSE TRIANGULATION IN PATCHES                                                  //
        // =================================================================================== //

        // Initialize variables -------------------------------------------------------------- //

            // First seed
            J = -1;
            I_ = (J + 1) % nSimplex;
            stop = true;
            counter = 0;
            while ((stop) && (counter < nSimplex)) {
                if (Patches[I_] == index0) {
                    J = I_;
                    stop = false;
                }
                I_ = (I_ + 1) % nSimplex;
                counter++;
            }

            // Flag
            flag = (J >= 0);

        // Loop until all surface is decomposed in patches ----------------------------------- //
        while (flag) {

            // Find patch for seed I_
            FindPatch(Patches, toll, index0, nPatches, J);
            nPatches += 1;

            // Next seed
            I_ = (J + 1) % nSimplex;
            stop = true;
            counter = 0;
            while ((stop) && (counter < nSimplex)) {
                if (Patches[I_] == index0) {
                    J = I_;
                    stop = false;
                }
                I_ = (I_ + 1) % nSimplex;
                counter++;
            }

            // Stopping criterion
            flag = (counter != nSimplex);

        }

        return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::LoopDecompose(int        &nPatches,
                                          ivector1D  &Patches,
                                          int         index0) {


        // =================================================================================== //
        // void Class_SurfTri::LoopDecompose(int        &nPatches,                             //
        //                                   ivector1D  &Patches,                              //
        //                                   int         index0)                               //
        //                                                                                     //
        // Decompose a surface patch into closed loops. Each patch is bounded either by a free //
        // or by T-junctions.                                                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - nPatches  : int, number of patches generated                                      //
        // - Patches   : ivector1D, with patch flag. Patch[i] is the patch ID                  //
        //               which the i-th simplex belongs to.                                    //
        // - index0    : int, patch ID to be decomposed                                        //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool                   flag, stop;

        // Counters
        int                    counter;
        int                    i, j;
        int                    I_, J;

        // =================================================================================== //
        // COMPUTE ADJACENCY MATRIX IF NOT ALREADY BUILT                                       //
        // =================================================================================== //
        if ((Adjacency.size() == 0) || (Adjacency.size() != nSimplex)) {
            BuildAdjacency();
        }

        // =================================================================================== //
        // DECOMPOSE TRIANGULATION IN PATCHES                                                  //
        // =================================================================================== //

        // Initialize variables -------------------------------------------------------------- //

            // First seed
            J = -1;
            I_ = (J + 1) % nSimplex;
            stop = true;
            counter = 0;
            while ((stop) && (counter < nSimplex)) {
                if (Patches[I_] == index0) {
                    J = I_;
                    stop = false;
                }
                I_ = (I_ + 1) % nSimplex;
                counter++;
            }

            // Flag
            flag = (J >= 0);

        // Loop until all surface is decomposed in patches ----------------------------------- //
        while (flag) {

            // Find patch for seed I_
            FindClosedLoop(Patches, index0, nPatches, J);
            nPatches += 1;

            // Next seed
            I_ = (J + 1) % nSimplex;
            stop = true;
            counter = 0;
            while ((stop) && (counter < nSimplex)) {
                if (Patches[I_] == index0) {
                    J = I_;
                    stop = false;
                }
                I_ = (I_ + 1) % nSimplex;
                counter++;
            }

            // Stopping criterion
            flag = (counter != nSimplex);

        }

        return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::RegionDecompose(int        &nPatches,
                                            ivector1D  &Patches,
                                            ivector1D  &rflag,
                                            int         index0) {


        // =================================================================================== //
        // void Class_SurfTri::RegionDecompose(int        &nPatches,                           //
        //                                     ivector1D  &Patches,                            //
        //                                     ivector1D  &rflag,                              //
        //                                     int         index0)                             //
        //                                                                                     //
        // Decompose a surface patch into labeled regions. Each patch is bounded either by a   //
        // free edge, or by T-junctions, or by a change in the value of flag.                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - nPatches  : int, number of patches generated                                      //
        // - Patches   : ivector1D, with patch flag. Patch[i] is the patch ID                  //
        //               which the i-th simplex belongs to.                                    //
        // - rflag     : ivector1D, with region labels.
        // - index0    : int, patch ID to be decomposed                                        //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool                   flag, stop;

        // Counters
        int                    counter;
        int                    i, j;
        int                    I_, J;

        // =================================================================================== //
        // COMPUTE ADJACENCY MATRIX IF NOT ALREADY BUILT                                       //
        // =================================================================================== //
        if (Adjacency.size() < nSimplex) {
            BuildAdjacency();
        }

        // =================================================================================== //
        // DECOMPOSE TRIANGULATION IN PATCHES                                                  //
        // =================================================================================== //

        // Initialize variables -------------------------------------------------------------- //

            // First seed
            J = -1;
            I_ = (J + 1) % nSimplex;
            stop = true;
            counter = 0;
            while ((stop) && (counter < nSimplex)) {
                if (Patches[I_] == index0) {
                    J = I_;
                    stop = false;
                }
                I_ = (I_ + 1) % nSimplex;
                counter++;
            }

            // Flag
            flag = (J >= 0);

        // Loop until all surface is decomposed in patches ----------------------------------- //
        while (flag) {

            // Find patch for seed I_
            FindRegion(Patches, rflag, index0, nPatches, J);
            nPatches += 1;

            // Next seed
            I_ = (J + 1) % nSimplex;
            stop = true;
            counter = 0;
            while ((stop) && (counter < nSimplex)) {
                if (Patches[I_] == index0) {
                    J = I_;
                    stop = false;
                }
                I_ = (I_ + 1) % nSimplex;
                counter++;
            }

            // Stopping criterion
            flag = (counter != nSimplex);

        }

        return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::ConnectedLoopDecompose(
    int        &nPatches,
    ivector1D  &Patches,
    int         index0
) {

// ========================================================================== //
// void Class_SurfTri::ConnectedLoopDecompose(                                //
//     int        &nPatches,                                                  //
//     ivector1D  &Patches,                                                   //
//     int         index0)                                                    //
//                                                                            //
// Decompose a surface patch into connected loops.                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nPatches  : int, number of patches generated                             //
// - Patches   : ivector1D, with patch flag. Patch[i] is the patch ID         //
//               which the i-th simplex belongs to.                           //
// - index0    : int, patch ID to be decomposed                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                   flag, stop;

// Counters
int                    counter;
int                    i, j;
int                    I_, J;

// ========================================================================== //
// COMPUTE ADJACENCY MATRIX IF NOT ALREADY BUILT                              //
// ========================================================================== //
if ((Adjacency.size() == 0) || (Adjacency.size() != nSimplex)) {
    BuildAdjacency();
}

// ========================================================================== //
// DECOMPOSE TRIANGULATION INTO CONNECTED LOOPS                               //
// ========================================================================== //

// Initialize variables ----------------------------------------------------- //

    // First seed
    J = -1;
    I_ = (J + 1) % nSimplex;
    stop = true;
    counter = 0;
    while ((stop) && (counter < nSimplex)) {
        if (Patches[I_] == index0) {
            J = I_;
            stop = false;
        }
        I_ = (I_ + 1) % nSimplex;
        counter++;
    }

    // Flag
    flag = (J >= 0);

// Loop until all surface is decomposed in patches -------------------------- //
while (flag) {

    // Find patch for seed I_
    FindConnectedLoop(Patches, index0, nPatches, J);
    nPatches += 1;

    // Next seed
    I_ = (J + 1) % nSimplex;
    stop = true;
    counter = 0;
    while ((stop) && (counter < nSimplex)) {
        if (Patches[I_] == index0) {
            J = I_;
            stop = false;
        }
        I_ = (I_ + 1) % nSimplex;
        counter++;
    }

    // Stopping criterion
    flag = (counter != nSimplex);

}

return; };

        // ----modificare con ricerca più efficiente - DA FINIRE ----------------------------- //
        void Class_SurfTri::FindSegment(int                    index0,
                                        ivector1D             &Patches,
                                        double                 toll_sharp,
                                        int                   &nSegments,
                                        vector<Class_SurfTri> &Segmentation) {

        // =================================================================================== //
        // void Class_SurfTri::FindSegment(int                   index0,                       //
        //                                 ivector1D             &Patches,                     //
        //                                 double                 toll_sharp,                  //
        //                                 int                   &nSegments,                   //
        //                                 vector<Class_SurfTri> &Segmentation)                //
        //                                                                                     //
        // Find segments in a given patch. Segments are sharp edges interior to the patch.     //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - index0       : int, patch index                                                   //
        // - Patches      : ivector1D, array with patches flags. Patches[i] is the patch index //
        //                  for the i-th simplex.                                              //
        // - toll         : double, tolerance for sharp edge detection                         //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - nSegments    : int, number of segments                                            //
        // - Segmentation : vector< Class_SurfTri >, segments                                  //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool                      flag, rule;
        double                    d;
        bvector1D                 check(nSimplex, false);
        ivector1D                 dummy(2, 0);

        // Counters
        int                       I_, J;
        int                       i, j, k;

        // =================================================================================== //
        // BUILD ADJACENCY IF NOT ALREADY COMPUTED                                             //
        // =================================================================================== //
        if (!(Adjacency.size() == nSimplex)) {
            BuildAdjacency();
        }

        // =================================================================================== //
        // COMPUTE NORMALS IF NOT ALREADY COMPUTED                                             //
        // =================================================================================== //
        if (!(Normal.size() == nSimplex)) {
            GenerateNormals();
        }

        // =================================================================================== //
        // FIND SEGMENTS                                                                       //
        // =================================================================================== //
        flag = true;
        for (i = 0; i < nSimplex; i++) {
            if (Patches[i] == index0) {
                for (j = 0; j < Simplex[i].size(); j++) {
                    rule = false;
                    rule = (Adjacency[i][j].size() > 1);
                    if (!rule) {

                        // Load adjcency
                        J = Adjacency[i][j][0];

                        // Check for sharp edge
                        if (J >= 0) {
                            if (Patches[J] == index0) {
                                d = norm2(Normal[i] - Normal[J]);
                                rule = (d > toll_sharp);
                            }
                        }
                    }

                    // Add edge to segmentation
                    if ((rule) && ((check[i] && check[J]) == false)) {

                        // Increase Segmentation
                        if (flag) {
                            nSegments += 1;
                            Segmentation.resize(nSegments);
                            I_ = nSegments - 1;
                            flag = false;
                        }

                        Segmentation[I_].AddVertex(Vertex[Simplex[i][j]]);
                        Segmentation[I_].AddVertex(Vertex[Simplex[i][(j+1) % Simplex[i].size()]]);
                        dummy[0] = Segmentation[I_].nVertex-2;
                        dummy[1] = Segmentation[I_].nVertex-1;
                        //dummy[0] = Simplex[i][j];
                        //dummy[1] = Simplex[i][(j+1) % Simplex[i].size()];
                        Segmentation[I_].AddSimplex(dummy);

                        // Update check
                        check[J] = true;
                    }

                    // Update check
                    check[i] = true;

                } //next j
            }
        } //next i

        return; };

