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
#ifndef __CLASS_VOLTRI__HH__
#define __CLASS_VOLTRI__HH__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <vector>
# include <sstream>
# include <fstream>
# include <iostream>
# include <algorithm>

// bitpit
# include <bitpit_operators.hpp>
# include <bitpit_IO.hpp>
# include <bitpit_SA.hpp>
# include <bitpit_LA.hpp>

// CC_lib
# include "Class_SurfTri.hpp"

// ========================================================================== //
// NAME SPACES                                                                //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITION                                                           //
// ========================================================================== //

// boolean vectors
typedef vector< bool >                 bvector1D;
typedef vector< bvector1D >            bvector2D;
typedef vector< bvector2D >            bvector3D;
typedef vector< bvector3D >            bvector4D;

// characters vectors
typedef vector< char >                 cvector1D;
typedef vector< cvector1D >            cvector2D;
typedef vector< cvector2D >            cvector3D;
typedef vector< cvector3D >            cvector4D;

// integer vectors
typedef vector< int >                  ivector1D;
typedef vector< ivector1D >            ivector2D;
typedef vector< ivector2D >            ivector3D;
typedef vector< ivector3D >            ivector4D;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

// array vectors
typedef array<double, 3>               a3vector1D;
typedef vector< a3vector1D >           a3vector2D;
typedef vector< a3vector2D >           a3vector3D;
typedef vector< a3vector3D >           a3vector4D;

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

// ========================================================================== //
// CLASS DEFINITION                                                           //
// ========================================================================== //
class Class_VolTri {

    // Structures =========================================================== //
    class e_info {
        public:
        unsigned int         n_vert;                                          // number of vertices
        unsigned int         n_faces;                                         // number of faces
        unsigned int         n_edges;                                         // number if edges
        ivector2D            faces;                                           // face-vertex connectivity
        ivector2D            edges;                                           // edge-vertex connectivity
    };

    // Private members ====================================================== //
    public:
    array<e_info, 15>        infos;                                           // element infos
    double                   toll;                                            // tolerance for distance checks
    vector<short>            e_type;                                          // element type ID

    // Public members ======================================================= //
    public:
    int                      nVertex;                                         // number of vertexes
    int                      nSimplex;                                        // number of simplicies
    int                      nFace;                                           // number of faces
    a3vector2D               Vertex;                                          // vertex coordinate list
    ivector2D                Simplex;                                         // Simplex-vertex connectivity
    ivector2D                Adjacency;                                       // Simplex-Simplex adjacency

    // Class constructors =================================================== //
    private:
    void Initialize_infos(                                                    // Initialize element infos at class creation
        void                                                                  // (input) none
    );
    public:
    Class_VolTri(                                                             // Default constructor
        void                                                                  // (input) none
    );
    Class_VolTri(                                                             // Custom constructor #1
        int                  ,                                                // (input) number of vertices
        int                                                                   // (input) number of simplicies
    );

    // Class destructors ==================================================== //
    ~Class_VolTri(                                                            // Default destructor
        void                                                                  // (input) none
    );

    // Resize operators ===================================================== //
    public:

    // Resize --------------------------------------------------------------- //
    void ResizeVertex(                                                        // Resize vertex list
        void                                                                  // (input) none
    );
    void ResizeSimplex(                                                       // Resize simplex-vertex connectivity
        int                  dim = 0                                          // (input) element type
    );
    void ResizeAdjacency(                                                     // Resize simplex-simplex adjacency
        int                  dim = 0                                          // (input/optional) element type
    );

    // Reshape -------------------------------------------------------------- //
    void ReshapeSimplex(                                                      // Reshape simplex-vertex connectivity
        int                  dim = 0                                          // (input/optional) element type
    );
    void ReshapeAdjacency(                                                    // Reshape simplex-simplex adjacency
        void                                                                  // (input) none
    );
    void ReshapeAdjacency(                                                    // Reshape simplex-simplex adjacency
        int                                                                   // (input) element type
    );

    // Destructors ========================================================== //
    public:

    // Destroy -------------------------------------------------------------- //
    void DestroyVertex(                                                       // Destroy vertex list
        void                                                                  // (input) none
    );
    void DestroySimplex(                                                      // Destroy simplex-vertex connectivity matrix
        void                                                                  // (input) none
    );
    void DestroyAdjacency(                                                    // Destroy simplex-simplex adjacency matrix
        void                                                                  // (input) none
    );

    // Clear ---------------------------------------------------------------- //
    void ClearVertex(                                                         // Clear vertex list
        void                                                                  // (input) none
    );
    void ClearSimplex(                                                        // Clear simplex-vertex connectivity matrix
        void                                                                  // (input) none
    );
    void ClearAdjacency(                                                      // Clear simplex-simplex adjacency matrix
        void                                                                  // (input) none
    );

    // Assignament operators ================================================ //
    public:

    // Copy ----------------------------------------------------------------- //
    Class_VolTri& operator=(                                                  // Assignament operator
        const Class_VolTri  &                                                 // (input) assignament operator
    );
    
    // Cleaning tools ======================================================= //

    // Tollerances ---------------------------------------------------------- //
    public:
    void SetTolerance(                                                        // Set tollerance for distance check
        void                                                                  // (input) none
    );
    void SetTolerance(                                                        // Set tollerance for distance check
        a3vector2D          &                                                 // (input) external vertex list
    );

    // Sorting algorithms --------------------------------------------------- //
    public:
    void BinSortV(                                                            // Sort tasselation vertices on regular bins
        ivector1D           &,                                                // (input/output) map vertex->bin
        int                  a = 128                                          // (input/optional) number of bins
    );
    void BinSortV(                                                            // Sort tasselation vertices on regular bins using an external vertex list
        a3vector2D          &,                                                // (input) external vertex list
        ivector1D           &,                                                // (input/output) map vertex->bin
        int                  a = 128                                          // (input/optional) number of bins
    );

    // Objects builder ------------------------------------------------------ //
    public:
    void BuildAdjacency(                                                      // Build simplex-simplex adjacency matrix
        void                                                                  // (input) none
    );
    void BuildAdjacency(                                                      // Build simplex-simplex adjacency matrix using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    void UpdateAdjacency(                                                     // Update simplex-simplex adjacency matrix
        ivector1D           &                                                 // (input) list of simplicies to be updated
    );
    void BuildFaces(                                                          // Build faces' data structure
        ivector2D           &,                                                // (input/output) simplex->face connectivity
        ivector2D           &                                                 // (input/output) face->simplex connectivity
    );

    // Counters ------------------------------------------------------------- //
    public:
    int CountIsolatedVertex(                                                  // Count isolated vertex
        void                                                                  // (input) none
    );
    int CountIsolatedVertex(                                                  // Count isolated vertex using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    int CountFreeVertex(                                                      // Count free vertices
        void                                                                  // (input) none
    );
    int CountFreeVertex(                                                      // Count free vertices using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    int CountDoubleVertex(                                                    // Counte duplicated vertices
        void                                                                  // (input) none
    );
    int CountDoubleVertex(                                                    // Counte duplicated vertices using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    int CountIsolatedSimplex(                                                 // Count isolated simplicies
        void                                                                  // (input) none
    );
    int CountIsolatedSimplex(                                                 // Count isolated simplicies using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    int CountFreeSimplex(                                                     // Count free simplicies
        void                                                                  // (input) none
    );
    int CountDoubleSimplex(                                                   // Count duplicated simplicies
        void                                                                  // (input) none
    );
    int CountDoubleSimplex(                                                   // Count duplicated simplicies using an external vertex list
        a3vector2D          &                                                 // (input) external vertex list
    );
    int CountFaces(                                                           // Count faces
        void                                                                  // (input) none
    );
    int CountFreeFaces(                                                       // Count free faces
        void                                                                  // (input) none
    );
    int Count0VolumeSimplex(                                                  // Count 0-volume simplicies in the tasselation
        void                                                                  // (input) none
    );
    int Count0VolumeSimplex(                                                  // Count 0-volume simplicies in the tasselation using an external vertex list
        a3vector2D          &                                                 // (input) vertex coordinate list
    );

    // Finders -------------------------------------------------------------- //
    public:
    ivector1D FindIsolatedVertex(                                             // Find isolated vertices in tasselation
        void                                                                  // (input) none
    );
    ivector1D FindIsolatedVertex(                                             // Find isolated vertices in tasselation using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    ivector1D FindFreeVertex(                                                 // Find free vertices in tasselation
        void                                                                  // (input) nonde
    );
    ivector1D FindFreeVertex(                                                 // Find free vertices in tasselation using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    ivector1D FindDoubleVertex(                                               // Find duplicated vertices
        int                  a = 128                                          // (input/optional) number of bins for vertex sorting
    );
    ivector1D FindDoubleVertex(                                               // Find duplicated vertices using an external vertex list
        a3vector2D          &,                                                // (input) External vertex list
        int                  a = 128                                          // (input/optional) number of bins for vertex sorting
    );
    ivector1D FindIsolatedSimplex(                                            // Find isolated simplicies
        void                                                                  // (input) none
    );
    ivector1D FindIsolatedSimplex(                                            // Find isolated simplicies using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    ivector1D FindFreeSimplex(                                                // Find free simplicies
        void                                                                  // (input) none
    );
    ivector1D FindFreeSimplex(                                                // Find free simplicies using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    ivector1D FindDoubleSimplex(                                              // Find duplicated simplicies
        void                                                                  // (input)
    );
    ivector1D FindDoubleSimplex(                                              // Find duplicated simplicies using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    ivector1D Find0VolumeSimplex(                                             // Find 0-volume simplicies in the tasselation
        void                                                                  // (input) none
    );
    ivector1D Find0VolumeSimplex(                                             // Find 0-volume simplicies in the tasselation using an external vertex list
        a3vector2D          &                                                 // (input) vertex coordinate list
    );

    // Cleaning tools ------------------------------------------------------- //
    public:
    void RemoveVertex(                                                        // Remove vertices from vertex list and update vertex-simplex connectivity
        ivector1D           &                                                 // (input) list of removable vertices
    );
    void RemoveVertex(                                                        // Remove vertices from vertex list, update vertex-simplex connectivity and returns mapper
        ivector1D           &,                                                // (input) list of removable vertices
        ivector1D           &                                                 // (input/output) mapper between old and new numeration
    );
    void CollapseDoubleVertex(                                                // Remove double vertices
        ivector1D           &,                                                // (input/output) list of collapsed vertices
        int                  a = 128                                          // (input/optional) number of bins used for vertex sorting
    );
    void CollapseDoubleVertex(                                                // Remove double vertices using an external vertex list
        a3vector2D          &,                                                // (input) external vertex list
        ivector1D           &,                                                // (input/output) list of collapsed vertices
        int                  a = 128                                          // (input/optional) number of bins used for vertex sorting
    );
    void RemoveSimplex(                                                       // Remove simplicies from simplex list and update adjacencies
        ivector1D           &                                                 // (input) List of removable simplicies
    );
    void RemoveSimplex(                                                       // Remove simplicies from simplex list, update adjacencies and returns mapper
        ivector1D           &,                                                // (input) List of removable simplicies
        ivector1D           &                                                 // (input/output) Mapper between old and new numeration
    );
    void RemoveIsolatedVertex(                                                // Remove isolated vertices
        void                                                                  // (input) none
    );
    void RemoveIsolatedVertex(                                                // Remove isolated vertices and returns mapper
        ivector1D           &                                                 // (input/output) mapper between old and new numeration
    );
    void RemoveDoubleVertex(                                                  // Remove duplicated vertices
        int                  a = 128                                          // (input/optional) number of bins for vertex sorting
    );
    void RemoveDoubleVertex(                                                  // Remove duplicated vertices and returns mapper
        ivector1D           &,                                                // (input/output) Mapper between old and new numeration
        int                  a = 128                                          // (input/optional) number of bins for vertex sorting
    );
    void RemoveIsolatedSimplex(                                               // Remove isolated simplicies
        void                                                                  // (input) none
    );
    void RemoveIsolatedSimplex(                                               // Remove isolated simplicies and returns mapper
        ivector1D           &                                                 // (input/output) mapper between old and new numeration
    );
    void RemoveIsolatedSimplex(                                               // Remove isolated simplicies using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    void RemoveDoubleSimplex(                                                 // Remove duplicated simplicies
        void                                                                  // (input) none
    );
    void RemoveDoubleSimplex(                                                 // Remove duplicated simplicies and returns map between old and new numbering
        ivector1D           &                                                 // (input) map between old and new numeration
    );
    void RemoveDoubleSimplex(                                                 // Remove duplicated simplicies using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );
    void Remove0VolumeSimplex(                                                // Remove 0-volume simplex from tasselation
        void                                                                  // (input) none
    );
    void Remove0VolumeSimplex(                                                // Remove 0-volume simplex from tasselation
        a3vector2D          &                                                 // (input) vertex coordinate list
    );
    void Clean(                                                               // Clean mesh
        void                                                                  // (input) none
    );
    void Clean(                                                               // Clean mesh using an external vertex list
        a3vector2D          &                                                 // (input) External vertex list
    );

    // Stats Tools ---------------------------------------------------------- //
    void Stats(                                                               // Compute tasselation stats
        ostream             &                                                 // (input) output stream
    );
    void Stats(                                                               // Compute tasselation stats using an external vertex list
        ostream             &,                                                // (input) output stream
        a3vector2D          &                                                 // (input) External vertex list
    );

    // Generation methods =================================================== //
    public:

    // Add elements --------------------------------------------------------- //
    void AddVertex(                                                           // Add a single vertex to vertex list
        a3vector1D          &                                                 // (input) vertex coordinates
    );
    void AddVertex(                                                           // Add a single vertex to vertex list
        dvector1D           &                                                 // (input) vertex coordinates
    );
    void AddVertices(                                                         // Add a multiple vertices to vertex list
        a3vector2D          &                                                 // (input) vertex coordinate list
    );
    void AddVertices(                                                         // Add a multiple vertices to vertex list
        dvector2D           &                                                 // (input) vertex coordinate list
    );
    void AddSimplex(                                                          // Add simplex to simplex list
        ivector1D           &,                                                // (input) simplex vertex coordinates
        int                                                                   // (input) simplex type
    );
    void AddSimplex(                                                          // Add simplex to simplex list
        a3vector2D          &,                                                // (input) coordinates of each simplex vertex
        int                                                                   // (input) simplex type id
    );
    void AddSimplicies(                                                       // Add multiple simplicies to simplex list
        ivector2D           &,                                                // (input) simplex-vertex connectivity
        ivector1D           &                                                 // (input) simplex type
    );
    void AddSimplicies(                                                       // Add multiple simplicies to simplex list
        a3vector2D          &,                                                // (input) vertex coordinate list
        ivector2D           &,                                                // (input) simplex-vertex connectivity
        ivector1D           &                                                 // (input) simplex type ids
    );
    void SetAdjacency(                                                        // Set adjacency for a given simplex
        int                  ,                                                // (input) simplex global index
        ivector1D           &                                                 // (input) adjacency list
    );
    void Append(                                                              // Append another mesh
        Class_VolTri        &                                                 // (input) source mesh
    );

    // Transformations ------------------------------------------------------ //
    void Scale(                                                               // Scale tasselation along each axis
        double               ,                                                // (input) scaling factor along x axis
        double               ,                                                // (input) scaling factor along y axis
        double                                                                // (input) scaling factor along z axis
    );
    void Translate(                                                           // Translate tasselation along axis
        double               ,                                                // (input) displacement along x axis
        double               ,                                                // (input) displacement along y axis
        double                                                                // (input) displacement along z axis
    );

    // Selection tools ====================================================== //
    public:

    // Bounding boxes ------------------------------------------------------- //
    void BoundingBox(                                                         // Compute limits of tasselation bounding box (3D case)
        array<double, 2>    &,                                                // (input/output) limits in the x-direction
        array<double, 2>    &,                                                // (input/output) limits in the y-direction
        array<double, 2>    &                                                 // (input/output) limits in the z-direction
    );
    void BoundingBox(                                                         // Compute limits of tasselation bounding box (3D case) using an external vertex list
        a3vector2D          &,                                                // (input) external vertex list
        array<double, 2>    &,                                                // (input/output) limits in the x-direction
        array<double, 2>    &,                                                // (input/output) limits in the y-direction
        array<double, 2>    &                                                 // (input/output) limits in the z-direction
    );

    // Check tools ========================================================== //
    public:

    // Topology check ------------------------------------------------------- //
    bool SameFace(                                                            // Check for coincident faces
        int                  ,                                                // (input) 1st simplex globla index
        int                  ,                                                // (input) face local index on 1st simplex
        int                  ,                                                // (input) 2nd simplex globla index
        int                                                                   // (input) face local index on 2nd simplex
    );
    bool SameSimplex(                                                         // Check for coincident simplicies 
        int                  ,                                                // (input) 1st simplex (global index)
        int                                                                   // (input) 2nd simplex (global index)
    );
    ivector1D FaceVertices(                                                   // Returns the global indices of vertices for the i-th face of a given simplex
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) face local index
    );

    // Geometry ------------------------------------------------------------- //
    a3vector1D Baricenter(                                                    // Compute simplex baricenter
        int                                                                   // (input) simplex global index
    );
    a3vector1D FaceCenter(                                                    // Compute face center
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) face local index
    );
    a3vector1D CircumCenter(                                                  // Compute simplex curcumcenter
        int                                                                   // (input) simplex global index
    );
    double EdgeLength(                                                        // Compute edge length
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) edge local index
    );
    double EdgeLength(                                                        // Compute edge length using an external vertex list
        int                  ,                                                // (input) simplex global index
        int                  ,                                                // (input) edge local index
        a3vector2D          &                                                 // (input) external vertex list
    );
    double minEdgeLength(                                                     // Compute min edge length
        int                                                                   // (input) simplex global index
    );
    double minEdgeLength(                                                     // Compute min edge length using an external vertex list
        int                  ,                                                // (input) simplex global index
        a3vector2D          &                                                 // (input) external vertex list
    );
    double maxEdgeLength(                                                     // Compute max edge length
        int                                                                   // (input) simplex global index
    );
    double maxEdgeLength(                                                     // Compute max edge length using an external vertex list
        int                  ,                                                // (input) simplex global index
        a3vector2D          &                                                 // (input) external vertex list
    );
    double FaceArea(                                                          // Compute simplex face area
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) face local index
    );
    double FaceArea(                                                          // Compute simplex face area using an external vertex list
        int                  ,                                                // (input) simplex global index
        int                  ,                                                // (input) face local index
        a3vector2D          &                                                 // (input) external vertex list
    );
    double minFaceArea(                                                       // Compute min face area
        int                                                                   // (input) simplex global index
    );
    double minFaceArea(                                                       // Compute min face area using an external vertex list
        int                  ,                                                // (input) simplex global index
        a3vector2D          &                                                 // (input) external vertex list
    );
    double maxFaceArea(                                                       // Compute max face area
        int                                                                   // (input) simplex global index
    );
    double maxFaceArea(                                                       // Compute max face area using an external vertex list
        int                  ,                                                // (input) simplex global index
        a3vector2D          &                                                 // (input) external vertex list
    );
    a3vector1D FaceNormal(                                                    // Compute simplex face normal
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) face local index
    );
    a3vector1D FaceNormal(                                                    // Compute simplex face normal using an external vertex list
        int                  ,                                                // (input) simplex global index
        int                  ,                                                // (input) face local index
        a3vector2D          &                                                 // (input) external vertex list
    );

            // Aggiustare
            double Volume(int T);                                 // Simplex volume
            double Volume(a3vector2D &X, int T);                  // Simplex volume

    // Selection tools ====================================================== //
    public:

    // Selection tools ------------------------------------------------------ //
    int face(                                                                 // Returns the local index of a face shared by two simplicies
        int                  ,                                                // (input) 1st simplex global index
        int                                                                   // (input) 2nd simplex global index
    );
    int vertex(                                                               // Returns the local index on a given simplex of a given vertex
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) vertex global index
    );
    int edge(                                                                 // Return the local index of a edge with specified vertices
        int                  ,                                                // (input) simplex global index
        ivector1D           &                                                 // (input) edge vertices
    );
    int ReturnSimplexID(                                                      // Returns the ID of simplex enclosing a given point
        a3vector1D          &,                                                // (input) point coordinates
        int                                                                   // (input) seed for search algorithm
    );
    ivector1D EdgeNeigh(                                                      // Find edge neighbors
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) edge local index
    );
    ivector1D VertNeigh(                                                      // Find vertex neighbors
        int                  ,                                                // (input) simplex global index
        int                                                                   // (input) vertex local index
    );
    void ExtractBoundaries(                                                   // Extract mesh boundaries
        Class_SurfTri       &,                                                // (input/output) surface tasselation
        int                  a = -1                                           // (input/optional) boundary condition flag
    );
            int ReturnTriangleID(a3vector1D    &);                 // Return triangle global index given a point

    // Refinement tools ===================================================== //
    public:
    void SwapFaceTri(                                                         // Swap face between adjacent triangles
        int                  ,                                                // (input) triangle global index
        int                                                                   // (input) face local index
    );
    void SplitTri3Tri(                                                        // Split  a triangle in 3 triangles using a vertex
        int                  ,                                                // (input) triangle global index
        int                                                                   // (input) vertex global index
    );
    void SplitHexa5Tetra(                                                     // Split an hexedron into 5 tetrahedron
        int                  ,                                                // (input) simplex global index
        short int                                                             // (input) splitting configuration
    );
    void SplitPyra2Tetra(                                                     // Split a pyramid into 2 tetrahedron
        int                  ,                                                // (input) simplex global index
        short int                                                             // (input) splitting configuration
    );

    // I/O methods ========================================================== //

    // dgf format ----------------------------------------------------------- //
    void Import_dgf(                                                          // Import volume mesh from dgf file
        string                                                                // (input) .dgf file name
    );
    void Export_dgf(                                                          // Export volume mesh into .dgf file
        string                                                                // (input) .dgf file name
    );

};

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
// none

#endif