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
#include "Class_VolTri.hpp"

// ========================================================================== //
// CLASS CONSTRUCTORS                                                         //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_VolTri::Initialize_infos(
    void
) {

// ========================================================================== //
// void Class_VolTri::Initialize_infos(                                       //
//     void)                                                                  //
//                                                                            //
// Initialize element infos.                                                  //
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
int                 id;

// Counters
// none

// ========================================================================== //
// INITIALIZE ELEMENT INFOS                                                   //
// ========================================================================== //

// Undef.  (element ID=0) --------------------------------------------------- //
{
    // element type id
    id = 0;
    
    // General infos
    infos[id].n_vert = 0;
    infos[id].n_faces = 0;
    infos[id].n_edges = 0;

}


// Tri     (element ID=5) --------------------------------------------------- //
{
    // element type id
    id = 5;

    // General infos
    infos[id].n_vert = 3;
    infos[id].n_faces = 3;
    infos[id].n_edges = 3;
    
    // Resize face data structure
    infos[id].faces.resize(infos[id].n_faces, ivector1D(2, -1));
    infos[id].edges.resize(infos[id].n_edges, ivector1D(1, -1));
    
    // Vertex numbering on element faces
    infos[id].faces[0][0] = 0;    infos[id].faces[0][1] = 1;
    infos[id].faces[1][0] = 1;    infos[id].faces[1][1] = 2;
    infos[id].faces[2][0] = 2;    infos[id].faces[2][1] = 0;

    // Vertex numbering on element edges
    infos[id].edges[0][0] = 0;
    infos[id].edges[1][0] = 1;
    infos[id].edges[2][0] = 2;
}

// Quad    (element ID=9) --------------------------------------------------- //
{
    // element type id
    id = 9;

    // General infos
    infos[id].n_vert = 4;
    infos[id].n_faces = 4;
    infos[id].n_edges = 4;
    
    // Resize face data structure
    infos[id].faces.resize(infos[id].n_faces, ivector1D(2, -1));
    infos[id].edges.resize(infos[id].n_edges, ivector1D(1, -1));
    
    // Vertex numbering on element faces
    infos[id].faces[0][0] = 0;    infos[id].faces[0][1] = 1;
    infos[id].faces[1][0] = 1;    infos[id].faces[1][1] = 2;
    infos[id].faces[2][0] = 2;    infos[id].faces[2][1] = 3;
    infos[id].faces[3][0] = 3;    infos[id].faces[3][1] = 0;

    // Vertex numbering on element edges
    infos[id].edges[0][0] = 0;
    infos[id].edges[1][0] = 1;
    infos[id].edges[2][0] = 2;
    infos[id].edges[3][0] = 3;
}

// Tetra   (element ID=10) -------------------------------------------------- //
{
    // element type id
    id = 10;

    // General infos
    infos[id].n_vert = 4;
    infos[id].n_faces = 4;
    infos[id].n_edges = 6;
    
    // Resize face data structure
    infos[id].faces.resize(infos[id].n_faces, ivector1D(3, -1));
    infos[id].edges.resize(infos[id].n_edges, ivector1D(2, -1));
    
    // Vertex numbering on element faces
    infos[id].faces[0][0] = 0;    infos[id].faces[0][1] = 2;    infos[id].faces[0][2] = 1;
    infos[id].faces[1][0] = 0;    infos[id].faces[1][1] = 1;    infos[id].faces[1][2] = 3;
    infos[id].faces[2][0] = 1;    infos[id].faces[2][1] = 2;    infos[id].faces[2][2] = 3;
    infos[id].faces[3][0] = 0;    infos[id].faces[3][1] = 3;    infos[id].faces[3][2] = 2;

    // Vertex numbering on element edges
    infos[id].edges[0][0] = 0;    infos[id].edges[0][1] = 2;
    infos[id].edges[1][0] = 2;    infos[id].edges[1][1] = 1;
    infos[id].edges[2][0] = 1;    infos[id].edges[2][1] = 0;
    infos[id].edges[3][0] = 0;    infos[id].edges[3][1] = 3;
    infos[id].edges[4][0] = 1;    infos[id].edges[4][1] = 3;
    infos[id].edges[5][0] = 2;    infos[id].edges[5][1] = 3;

}

// Hexa    (element ID=12) -------------------------------------------------- //
{
    // element type id
    id = 12;

    // General infos
    infos[id].n_vert = 8;
    infos[id].n_faces = 6;
    infos[id].n_edges = 12;
    
    // Resize face data structure
    infos[id].faces.resize(infos[id].n_faces, ivector1D(4, -1));
    infos[id].edges.resize(infos[id].n_edges, ivector1D(2, -1));
    
    // Vertex numbering on element faces
    infos[id].faces[0][0] = 0;    infos[id].faces[0][1] = 3;    infos[id].faces[0][2] = 2;    infos[id].faces[0][3] = 1;
    infos[id].faces[1][0] = 1;    infos[id].faces[1][1] = 2;    infos[id].faces[1][2] = 6;    infos[id].faces[1][3] = 5;
    infos[id].faces[2][0] = 2;    infos[id].faces[2][1] = 3;    infos[id].faces[2][2] = 7;    infos[id].faces[2][3] = 6;
    infos[id].faces[3][0] = 0;    infos[id].faces[3][1] = 4;    infos[id].faces[3][2] = 7;    infos[id].faces[3][3] = 3;
    infos[id].faces[4][0] = 0;    infos[id].faces[4][1] = 1;    infos[id].faces[4][2] = 5;    infos[id].faces[4][3] = 4;
    infos[id].faces[5][0] = 4;    infos[id].faces[5][1] = 5;    infos[id].faces[5][2] = 6;    infos[id].faces[5][3] = 7;

    // Vertex numbering on element edges
    infos[id].edges[0][0] = 0;    infos[id].edges[0][1] = 1;
    infos[id].edges[1][0] = 1;    infos[id].edges[1][1] = 2;
    infos[id].edges[2][0] = 2;    infos[id].edges[2][1] = 3;
    infos[id].edges[3][0] = 3;    infos[id].edges[3][1] = 0;
    infos[id].edges[4][0] = 4;    infos[id].edges[4][1] = 5;
    infos[id].edges[5][0] = 5;    infos[id].edges[5][1] = 6;
    infos[id].edges[6][0] = 6;    infos[id].edges[6][1] = 7;
    infos[id].edges[7][0] = 7;    infos[id].edges[7][1] = 4;
    infos[id].edges[8][0] = 0;    infos[id].edges[8][1] = 4;
    infos[id].edges[9][0] = 1;    infos[id].edges[9][1] = 5;
    infos[id].edges[10][0] = 2;   infos[id].edges[10][1] = 6;
    infos[id].edges[11][0] = 3;   infos[id].edges[11][1] = 7;
}

// Prism   (element ID=13) -------------------------------------------------- //
{
    // element type id
    id = 13;

    // General infos
    infos[id].n_vert = 6;
    infos[id].n_faces = 5;
    infos[id].n_edges = 9;
    
    // Resize face data structure
    infos[id].faces.resize(infos[id].n_faces);
    infos[id].faces[0].resize(3, -1);
    infos[id].faces[1].resize(3, -1);
    for (int i = 2; i < infos[id].n_faces; ++i) {
        infos[id].faces[i].resize(4, -1);
    }
    infos[id].edges.resize(infos[id].n_edges, ivector1D(2, -1));
    
    // Vertex numbering on element faces
    infos[id].faces[0][0] = 0;    infos[id].faces[0][1] = 1;    infos[id].faces[0][2] = 2;
    infos[id].faces[1][0] = 3;    infos[id].faces[1][1] = 5;    infos[id].faces[1][2] = 4;
    infos[id].faces[2][0] = 0;    infos[id].faces[2][1] = 3;    infos[id].faces[2][2] = 4;    infos[id].faces[2][3] = 1;
    infos[id].faces[3][0] = 1;    infos[id].faces[3][1] = 4;    infos[id].faces[3][2] = 5;    infos[id].faces[3][3] = 2;
    infos[id].faces[4][0] = 0;    infos[id].faces[4][1] = 2;    infos[id].faces[4][2] = 5;    infos[id].faces[4][3] = 3;

    // Vertex numbering on element edges
    infos[id].edges[0][0] = 0;    infos[id].edges[0][1] = 1;
    infos[id].edges[1][0] = 1;    infos[id].edges[1][1] = 2;
    infos[id].edges[2][0] = 2;    infos[id].edges[2][1] = 0;
    infos[id].edges[3][0] = 3;    infos[id].edges[3][1] = 5;
    infos[id].edges[4][0] = 5;    infos[id].edges[4][1] = 4;
    infos[id].edges[5][0] = 4;    infos[id].edges[5][1] = 3;
    infos[id].edges[6][0] = 0;    infos[id].edges[6][1] = 3;
    infos[id].edges[7][0] = 1;    infos[id].edges[7][1] = 4;
    infos[id].edges[8][0] = 2;    infos[id].edges[8][1] = 5;

}

// Pyramid (element ID=14) -------------------------------------------------- //
{
    // element type id
    id = 14;
    
    // General infos
    infos[id].n_vert = 5;
    infos[id].n_faces = 5;
    infos[id].n_edges = 8;

    // Resize face data structure
    infos[id].faces.resize(infos[id].n_faces);
    infos[id].faces[0].resize(4, -1);
    for (int i = 1; i < infos[id].n_faces; ++i) {
        infos[id].faces[i].resize(3, -1);
    } //next i
    infos[id].edges.resize(infos[id].n_edges, ivector1D(2, -1));
    
    // Vertex numbering on element faces
    infos[id].faces[0][0] = 0;    infos[id].faces[0][1] = 3;    infos[id].faces[0][2] = 2;    infos[id].faces[0][3] = 1;
    infos[id].faces[1][0] = 0;    infos[id].faces[1][1] = 1;    infos[id].faces[1][2] = 4;
    infos[id].faces[2][0] = 1;    infos[id].faces[2][1] = 2;    infos[id].faces[2][2] = 4;
    infos[id].faces[3][0] = 2;    infos[id].faces[3][1] = 3;    infos[id].faces[3][2] = 4;
    infos[id].faces[4][0] = 3;    infos[id].faces[4][1] = 0;    infos[id].faces[4][2] = 4;

    // Vertex numbering on element edges
    infos[id].edges[0][0] = 0;    infos[id].edges[0][1] = 3;
    infos[id].edges[1][0] = 3;    infos[id].edges[1][1] = 2;
    infos[id].edges[2][0] = 2;    infos[id].edges[2][1] = 1;
    infos[id].edges[3][0] = 1;    infos[id].edges[3][1] = 0;
    infos[id].edges[4][0] = 0;    infos[id].edges[4][1] = 4;
    infos[id].edges[5][0] = 3;    infos[id].edges[5][1] = 4;
    infos[id].edges[6][0] = 2;    infos[id].edges[6][1] = 4;
    infos[id].edges[7][0] = 1;    infos[id].edges[7][1] = 4;
    
}

return; };

// -------------------------------------------------------------------------- //
Class_VolTri::Class_VolTri(
    void
) {

// ========================================================================== //
// Class_VolTri::Class_VolTri(                                                //
//     void)                                                                  //
//                                                                            //
// Standard constructor for Class_VolTri variables                            //
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
// none

// ========================================================================== //
// CONSTRUCTOR                                                                //
// ========================================================================== //

// # vertex
nVertex = 0;

// # edges
nFace = 0;

// # simplicies
nSimplex = 0;

// element info
Initialize_infos();

return; };

// -------------------------------------------------------------------------- //
Class_VolTri::Class_VolTri(
    int           nV,
    int           nS
) {

// ========================================================================== //
// Class_VolTri::Class_VolTri(                                                //
//     int           nV,                                                      //
//     int           nS)                                                      //
//                                                                            //
// Custom constructor #1 for Class_VolTri variables.                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nV   : int, number of vertexes                                           //
// - nS   : int, number of simplicies                                         //
// ========================================================================== //
        // OUTPUT                                                             //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CONSTRUCTOR                                                                //
// ========================================================================== //

// # vertex
nVertex = nV;

// # Faces
nFace = 0;

// # simplicies
nSimplex = nS;

// Resize vertex list
ResizeVertex();

// Resize simplex-vertex connectivity
ResizeSimplex();

// Element type
e_type.resize(nS, -1);

// element info
Initialize_infos();

return; };

// ========================================================================== //
// CLASS DESTRUCTOR                                                           //
// ========================================================================== //

// -------------------------------------------------------------------------- //
Class_VolTri::~Class_VolTri(
    void
) {

// ========================================================================== //
// Class_VolTri::~Class_VolTri(                                               //
//     void)                                                                  //
//                                                                            //
// Standard destructor for Class_VolTri variables                             //
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
// none

// ========================================================================== //
// DESTROY CONTENTS                                                           //
// ========================================================================== //

// Tasselation dimensions
nVertex = 0;
nSimplex = 0;
nFace = 0;

// Destroy variables
DestroyVertex();
DestroySimplex();
DestroyAdjacency();

return; };

// ========================================================================== //
// RESIZE OPERATORS                                                           //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_VolTri::ResizeVertex(
    void
) {

// ========================================================================== //
// void Class_VolTri::ResizeVertex(                                           //
//     void)                                                                  //
//                                                                            //
// Resize vertex list to nVertex rows. Each of the new rows will have d       //
// entries. Resize is not destructive, i.e. data previously stored in Vertex  //
// will not be erased.                                                        //
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
a3vector1D              tmp;

// Counters
// none

// ========================================================================== //
// RESIZE VERTEX LIST                                                         //
// ========================================================================== //
tmp.fill(0.0);
Vertex.resize(nVertex, tmp);

return;};

// -------------------------------------------------------------------------- //
void Class_VolTri::ResizeSimplex(
    int           d
) {

// ========================================================================== //
// void Class_VolTri::ResizeSimplex(                                          //
//     int           d)                                                       //
//                                                                            //
// Resize simplex-vertex connectivity matrix to nSimplex rows. Each of the    //
// new rows in Simplex will have a number of entries corresponding to element //
// type. Data previously stored in Simplex are not erased during resize.      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - d     : int, element type.                                               //
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
// none

// ========================================================================== //
// RESIZE SIMPLEX-VERTEX CONNECTIVITY MATRIX                                  //
// ========================================================================== //
if ((d < 0) || (d > 14)) { return; };
Simplex.resize(nSimplex, ivector1D(infos[d].n_vert, -1));
e_type.resize(nSimplex, d);

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::ResizeAdjacency(
    int           d
) {

// ========================================================================== //
// void Class_VolTri::ResizeAdjacency(                                        //
//     int           d)                                                       //
//                                                                            //
// Resize simplex-simplex adjacencies matrix to nSimplex rows. Each of the    //
// new rows will have a number of entries depending on the element type.      //
// Data previously stored in Adjacency are not erased during resize.          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - d          : int (optional), element type                                //
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
// none

// ========================================================================== //
// RESIZE SIMPLEX-SIMPLEX CONNECTIVITY MATRIX                                 //
// ========================================================================== //
if ((d < 0) || (d > 14)) { return; }
Adjacency.resize(nSimplex, ivector1D(infos[d].n_faces, -1));

return; };

// ========================================================================== //
// RESHAPE OPERATORS                                                          //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_VolTri::ReshapeSimplex(
    int           d
) {

// ========================================================================== //
// void Class_VolTri::ReshapeSimplex(                                         //
//     int           d)                                                       //
//                                                                            //
// Reshape simplex-vertex connectivity matrix to nSimplex rows. Each row in   //
// Simplex will have d entries. Data previously stored in Simplex are not     //
// erased during reshape.                                                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - d    : int (optional), element type                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int       old_size = Simplex.size();

// Counters
int       i;

// ========================================================================== //
// RESIZE SIMPLEX-VERTEX CONNECTIVITY MATRIX                                  //
// ========================================================================== //

// Resize Simplex
ResizeSimplex(d);

// Reshape Simplex
old_size = min(old_size, nSimplex);
for (i = 0; i < old_size; i++) {
    e_type[i] = d;
    Simplex[i].resize(infos[d].n_vert, -1);
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::ReshapeAdjacency(
    void
) {

// ========================================================================== //
// void Class_VolTri::ReshapeAdjacency(                                       //
//     void)                                                                  //
//                                                                            //
// Reshape simplex-simplex adjacencies matrix to nSimplex rows. Each row      //
// will have a number of entries depending on the element type. Data          //
// previously stored in Adjacency are not erased during reshape.              //
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
int     i;

// ========================================================================== //
// RESIZE SIMPLEX-SIMPLEX CONNECTIVITY MATRIX                                 //
// ========================================================================== //

// Resize adjacency
ResizeAdjacency();

// Reshape Adjacency
for (i = 0; i < nSimplex; i++) {
    Adjacency[i].resize(infos[e_type[i]].n_faces, -1);
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::ReshapeAdjacency(
    int           d
) {

// ========================================================================== //
// void Class_VolTri::ReshapeAdjacency(                                       //
//     int           d)                                                       //
//                                                                            //
// Reshape simplex-simplex adjacencies matrix to nSimplex rows. Each row      //
// will have d entries. Data previously stored in Adjacency are not erased    //
// during reshape.                                                            //
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
int     old_size = Adjacency.size();

// Counters
int     i;

// ========================================================================== //
// RESIZE SIMPLEX-SIMPLEX CONNECTIVITY MATRIX                                 //
// ========================================================================== //

// Resize adjacency
ResizeAdjacency(d);

// Reshape Adjacency
old_size = min(old_size, nSimplex);
for (i = 0; i < old_size; i++) {
    Adjacency[i].resize(infos[d].n_faces, -1);
} //next i

return; };

// ========================================================================== //
// DESTRUCTORS                                                                //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_VolTri::DestroyVertex(
    void
) {

// ========================================================================== //
// void Class_VolTri::DestroyVertex(                                          //
//     void)                                                                  //
//                                                                            //
// Destroy vertex coordinate list. All data stored in Vertex are lost.        //
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
// none

// ========================================================================== //
// DESTROY VERTEX                                                             //
// ========================================================================== //
nVertex = 0;
ResizeVertex();

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::DestroySimplex(
    void
) {

// ========================================================================== //
// void Class_VolTri::DestroySimplex(                                         //
//     void)                                                                  //
//                                                                            //
// Destroy Simplex list. All data stored in Simplex are lost.                 //
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
// none

// ========================================================================== //
// DESTROY SIMPLEX                                                            //
// ========================================================================== //
nSimplex = 0;
ResizeSimplex();

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::DestroyAdjacency(
    void
) {

// ========================================================================== //
// void Class_VolTri::DestroyAdjacency(                                       //
//     void)                                                                  //
//                                                                            //
// Destroy Adjacency matrix. All data stored in Adjacency are lost.           //
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
// none

// ========================================================================== //
// DESTROY SIMPLEX-SIMPLEX ADJACENCY                                          //
// ========================================================================== //
Adjacency.resize(0);

return; };

// ========================================================================== //
// CLEAR                                                                      //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_VolTri::ClearVertex(
    void
) {

// ========================================================================== //
// void Class_VolTri::ClearVertex(                                            //
//     void)                                                                  //
//                                                                            //
// Clear contents in vertex coordinate list, without altering shape.          //
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
int          i, j, d;

// ========================================================================== //
// CLEAR CONTENTS IN VERTEX                                                   //
// ========================================================================== //
for (i = 0; i < nVertex; i++) {
    Vertex[i].fill(0.0);
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::ClearSimplex(
    void
) {

// ========================================================================== //
// void Class_VolTri::ClearSimplex(                                           //
//     void)                                                                  //
//                                                                            //
// Clear contents in simplex-vertex connectivity matrix, without altering     //
// shape.                                                                     //
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
int         i, j, n;

// ========================================================================== //
// CLEAR CONTENTS IN VERTEX                                                   //
// ========================================================================== //
for (i = 0; i < nSimplex; i++) {
    n = infos[e_type[i]].n_vert;
    for (j = 0; j < n; j++) {
        Simplex[i][j] = -1;
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
void Class_VolTri::ClearAdjacency(
    void
) {

// ========================================================================== //
// void Class_VolTri::ClearAdjacency()                                        //
//                                                                            //
// Clear contents in adjacency matrix, without altering shape.                //
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
int     i, j, n;

// ========================================================================== //
// CLEAR CONTENTS IN VERTEX                                                   //
// ========================================================================== //
for (i = 0; i < nSimplex; i++) {
    n = infos[e_type[i]].n_faces;
    for (j = 0; j < n; j++) {
        Adjacency[i][j] = -1;
    } //next j
} //next i

return; };

// ========================================================================== //
// ASSIGNAMENT                                                                //
// ========================================================================== //

// -------------------------------------------------------------------------- //
Class_VolTri& Class_VolTri::operator=(
    const Class_VolTri &B
) {

// ========================================================================== //
// Class_VolTri& Class_VolTri::operator=(const Class_VolTri &B)               //
//                                                                            //
// Assignament operator for Class_VolTri variables.                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - B   : Class_VolTri, source variable                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool     flag_a;

// Counters
int      i, j;

// ========================================================================== //
// PARAMETERS                                                                 //
// ========================================================================== //
flag_a = ((B.Adjacency.size() > 0) && (B.Adjacency.size() >= B.nSimplex));

// ========================================================================== //
// COPY SOURCE VARIABLES                                                      //
// ========================================================================== //

// Mesh dimensions ---------------------------------------------------------- //
nVertex = B.nVertex;
nSimplex = B.nSimplex;
nFace = B.nFace;

// Copy data ---------------------------------------------------------------- //

// Copy vertex coordinate list
Vertex = B.Vertex;

// Copy simplex-vertex connectivity matrix
e_type = B.e_type;
Simplex = B.Simplex;

// Copy simplex-simplex adjacency matrix
if (flag_a) {
    Adjacency = B.Adjacency;
}

return(*this); };


