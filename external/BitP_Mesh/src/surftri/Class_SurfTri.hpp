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
#ifndef __CLASS_SURFTRI__HH__
#define __CLASS_SURFTRI__HH__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template library
# include <vector>
# include <sstream>
# include <fstream>
# include <iostream>
# include <algorithm>
# include <unordered_map>

// bitpit
# include <bitpit_operators.hpp>
# include <bitpit_IO.hpp>
# include <bitpit_SA.hpp>

// ========================================================================== //
// NAME SPACES                                                                //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
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

// double array
typedef array< double,3 >              darray3E;  
typedef vector< darray3E >             dvecarr3E; 

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

// ========================================================================== //
// CLASS DEFINITION                                                           //
// ========================================================================== //
class Class_SurfTri {

    // Members ============================================================== //
    public:
        int                     nVertex;                                          // # of vertexes
        int                     nSimplex;                                         // # of simplicies
        int                     nEdge;                                            // # of edges
        int                     dimensions;                                       // # of dimensions
        double                  toll;                                             // tolerance for distance check

        dvecarr3E               Vertex;                                           // vertex list
        dvecarr3E               Normal;                                           // normals list
        dvecarr3E               ENormal;                                          // edge normals
        dvecarr3E               VNormal;                                          // vertex normals
        ivector2D               Simplex;                                          // Simplex-vertex connectivity
        ivector2D               Edge;                                             // Edge-vertex connectivity
        ivector2D               Simplex2Edge;                                     // Simplex-edge connectivity
        ivector3D               Adjacency;                                        // Simplex-Simplex adjacency

        // Class constructors =================================================== //
    public:
        Class_SurfTri(                                                            // Default constructor
                int                 dim=3                                         // (input) number of dimensions
                );

        Class_SurfTri(                                                            // Custom constructor #1
                int                  ,                                                // (input) number of vertices
                int                  ,                                                // (input) number of simplicies
                int                  dim=3                                         // (input) number of dimensions
                );

        // Class destructors ==================================================== //
    public:
        ~Class_SurfTri(                                                           // Default destructor
                void                                                                  // (input) none
                );



        // set methods ========================================================== //
    public:
        void SetDimensions(                                                           // Default destructor
                int                                                              // (input) # of dimensions
                );



        // Resize operators ===================================================== //

        // Resize --------------------------------------------------------------- //
    public:
        void ResizeVertex(                                                        // Resize vertex list to the actual number of vertices
                void                                                                  // (input) none
                );
        void ResizeNormal(                                                        // Resize normals' list
                void                                                                  // (input) none
                );
        void ResizeENormal(                                                       // Resize edge normals' list
                void                                                                  // (input) none
                );
        void ResizeVNormal(                                                       // Resize vertex normals' list
                void                                                                  // (input) none
                );
        void ResizeSimplex(                                                       // Resize simplex-vertex connectivity
                void                                                                  // (input) none
                );
        void ResizeSimplex(                                                       // Resize simplex-vertex connectivity
                int                                                                   // (input) simplex type
                );
        void ResizeEdge(                                                          // Resize edge-vertex connectivity
                int       d = 2                                                       // (input/optional) number of entries in the i-th row of edge-vertex connectivity
                );
        void ResizeSimplex2Edge(                                                  // Resize simplex-edge connectivity
                int       d = 3                                                       // (input/optional) number of edges for each of the new row
                );
        void ResizeAdjacency(                                                     // Resize simplex-simplex adjacency
                void                                                                  // (input) none
                );
        void ResizeAdjacency(                                                     // Resize simplex-simplex adjacency
                int                                                                   // (input/optional) simplex type
                );

        // Reshape -------------------------------------------------------------- //
    public:
        void ReshapeSimplex(                                                      // Reshape simplex-vertex connectivity
                int                                                                   // (input/optional) simplex type
                );
        void ReshapeSimplex2Edge(                                                 // Reshape simplex-edge connectivity
                void                                                                  // (input) none
                );
        void ReshapeAdjacency(                                                    // Reshape simplex-simplex adjacency
                void                                                                  // (input) none
                );
        void ReshapeAdjacency(                                                    // Reshape simplex-simplex adjacency
                int                                                                   // (input/optional) simplex type
                );

        // Destructors ========================================================== //

        // Destroy -------------------------------------------------------------- //
    public:
        void DestroyVertex(                                                       // Destroy vertex list
                void                                                                  // (input) none
                );
        void DestroyNormal(                                                       // Destroy normals list
                void                                                                  // (input) none
                );
        void DestroyENormal(                                                      // Destroy edge normals list
                void                                                                  // (input) none
                );
        void DestroyVNormal(                                                      // Destroy vertex normals list
                void                                                                  // (input) none
                );
        void DestroySimplex(                                                      // Destroy vertex-simplex connectivity matrix
                void                                                                  // (input) none
                );
        void DestroyEdge(                                                         // Destroy edge-vertex connectivity matrix
                void                                                                  // (input) none
                );
        void DestroySimplex2Edge(                                                 // Destroy simplex-edge connectivity matrix
                void                                                                  // (input) none
                );
        void DestroyAdjacency(                                                    // Destroy adjacency matrix
                void                                                                  // (input) none
                );

        // Clear ---------------------------------------------------------------- //
    public:
        void ClearVertex(                                                         // Clear vertex list
                void                                                                  // (input) none
                );
        void ClearNormal(                                                         // Clear normals list
                void                                                                  // (input) none
                );
        void ClearENormal(                                                        // Clear edge normals list
                void                                                                  // (input) none
                );
        void ClearVNormal(                                                        // Clear vertex normals list
                void                                                                  // (input) none
                );
        void ClearSimplex(                                                        // Clear vertex-simplex connectivity matrix
                void                                                                  // (input) none
                );
        void ClearEdge(                                                           // Clear edge-vertex connectivity matrix
                void                                                                  // (input) none
                );
        void ClearSimplex2Edge(                                                   // Clear simplex-edge connectivity matrix
                void                                                                  // (input) none
                );
        void ClearAdjacency(                                                      // Clear adjacency matrix
                void                                                                  // (input) none
                );

        // Assignament operators ================================================ //

        // Copy ----------------------------------------------------------------- //
    public:
        Class_SurfTri& operator=(                                                 // Assignament operator
                const Class_SurfTri &                                                 // (input) source variable
                );

        // Cleaning methods ===================================================== //

        // Sorting algorithms --------------------------------------------------- //
    private:
        void BinSortV(                                                            // Sort tasselation vertices on regular bins
                ivector1D           &,                                                // (input/output) map vertex->bin
                int                  a = 128                                          // (input/optional) number of bins
                );
        void BinSortV(                                                            // Sort tasselation vertices on regular bins using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                ivector1D           &,                                                // (input/output) map vertex->bin
                int                  a = 128                                          // (input/optional) number of bins
                );

        // Tollerances ---------------------------------------------------------- //
    public:
        void SetTolerance(                                                        // Set tollerance for distance check
                void                                                                  // (input) none
                );
        void SetTolerance(                                                        // Set tollerance for distance check
                dvecarr3E           &                                                 // (input) external vertex list
                );

        // Objects builder ------------------------------------------------------ //
    public:
        void FixNodeNumb(                                                         // Fix local node numbering according to normals
                void                                                                  // (input) nonde
                );
        void FixNodeNumb(                                                         // Fix local node numbering according to normals using an external vertex list
                dvecarr3E           &                                                 // (input) external vertex list
                );
        void FixNodeNumb(                                                         // Fix local node numbering given a seed
                int                                                                   // (input) global index of seed simplex
                );
        void InvertOrder(                                                         // Invert simplex numbering
                void                                                                  // (input) none
                );
        void GenerateNormals(                                                     // Generate simplex normals
                void                                                                  // (input) none
                );
        void GenerateNormals(                                                     // Generate simplex normals using vertex coordinate list
                dvecarr3E           &                                                 // (input) external vertex list
                );
        void GenerateENormals(                                                    // Generate edge's normals
                void                                                                  // (input) none
                );
        void GenerateENormals(                                                    // Generate edge's normals using an external vertex list
                dvecarr3E           &                                                 // (input) external vertex list
                );
        void GenerateVNormals(                                                    // Generate vertex normals
                unsigned char        flag = 0                                         // (input/optional) algorithm selection flag
                );
        void GenerateVNormals(                                                    // Generate vertex normals using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                unsigned char        flag = 0                                         // (input/optional) algorithm selection flag
                );
        void BuildAdjacency(                                                      // Build simplex-simplex adjacency matrix
                void                                                                  // (input) none
                );
        void BuildAdjacency(                                                      // Build simplex-simplex adjacency matrix using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        void UpdateAdjacency(                                                     // Update simplex-simplex adjacency matrix
                ivector1D           &                                                 // (input) list of simplicies to be updated
                );
        void UpdateAdjacency(                                                     // Update simplex-simplex adjacency matrix using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                ivector1D           &                                                 // (input) list of simplicies to be updated
                );
        void BuildEdges(                                                          // Build edge list
                void                                                                  // (input) none
                );

        // Counters ------------------------------------------------------------- //
    public:
        int CountIsolatedVertex(                                                  // Count isolated vertex in tasselation
                void                                                                  // (input) none
                );
        int CountIsolatedVertex(                                                  // Count isolated vertex in tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        int CountFreeVertex(                                                      // Count free vertices in the tasselation
                void                                                                  // (input) none
                );
        int CountFreeVertex(                                                      // Count free vertices in the tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        int CountDoubleVertex(                                                    // Counte duplicated vertices
                void                                                                  // (input) none
                );
        int CountDoubleVertex(                                                    // Counte duplicated vertices using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        int CountIsolatedSimplex(                                                 // Count isolated simplicies in the tasselation
                void                                                                  // (input) none
                );
        int CountIsolatedSimplex(                                                 // Count isolated simplicies in the tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        int CountFreeSimplex(                                                     // Count isolated simplicies in the tasselation
                void                                                                  // (input) none
                );
        int CountDoubleSimplex(                                                   // Count duplicated simplex
                void                                                                  // (input) none
                );
        int CountDoubleSimplex(                                                   // Count duplicated simplicies using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );

        /*! Count true duplicated simplicies: if two simplicies of the same kind are defined by the same
          vertices but different vertex ordering, they are considered duplicated */ 
        int CountTrueDoubleSimplex(    
                void			   
                );				    

        /*! Count true duplicated simplicies using external list: if two simplicies of the same kind are
          defined by the same vertices but different vertex ordering, they are considered duplicated. 
          An external dvecarr3E vertex list is required as INPUT. */ 
        int CountTrueDoubleSimplex(    
                dvecarr3E &		    
                );				   

        int CountEdges(                                                           // Count tasselation edges
                void                                                                  // (input) none
                );
        int CountFreeEdges(                                                       // Count tasselation free edges
                void                                                                  // (input) none
                );

        // Find algorithms ------------------------------------------------------ //
    public:
        ivector1D FindIsolatedVertex(                                             // Find isolated vertices in tasselation
                void                                                                  // (input) none
                );
        ivector1D FindIsolatedVertex(                                             // Find isolated vertices in tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        ivector1D FindFreeVertex(                                                 // Find free vertices in tasselation
                void                                                                  // (input) nonde
                );
        ivector1D FindFreeVertex(                                                 // Find free vertices in tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        ivector1D FindDoubleVertex(                                               // Find duplicated vertices in tasselation
                int                  a = 128                                          // (input/optional) number of bins for vertex sorting
                );
        ivector1D FindDoubleVertex(                                               // Find duplicated vertices in tasselation using an external vertex list
                dvecarr3E           &,                                                // (inptu) External vertex list
                int                  a = 128                                          // (input/optional) number of bins for vertex sorting
                );
        ivector1D FindIsolatedSimplex(                                            // Find isolated simplicies in the tasselation
                void                                                                  // (input) none
                );
        ivector1D FindIsolatedSimplex(                                            // Find isolated simplicies in the tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        ivector1D FindFreeSimplex(                                                // Find free simplicies in the tasselation
                void                                                                  // (input) none
                );
        ivector1D FindFreeSimplex(                                                // Find free simplicies in the tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        ivector1D FindDoubleSimplex(                                              // Find duplicated simplicies
                void                                                                  // (input)
                );
        ivector1D FindDoubleSimplex(                                              // Find duplicated simplicies using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );

        /*! Find true duplicated simplicies: if two simplicies of the same kind are defined by the same
          vertices but different vertex ordering, they are considered duplicated */ 
        ivector1D FindTrueDoubleSimplex( 
                void			
                );

        /*! Find true duplicated simplicies using external list: if two simplicies of the same kind are
          defined by the same vertices but different vertex ordering, they are considered duplicated. 
          An external dvecarr3E vertex list is required as INPUT. */ 
        ivector1D FindTrueDoubleSimplex( 
                dvecarr3E &		
                );				

        ivector1D Find0AreaSimplex(                                               // Find 0-area simplicies in the surface tasselation
                void                                                                  // (input) none
                );
        ivector1D Find0AreaSimplex(                                               // Find 0-area simplicies in the surface tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
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
        void RemoveVertex(                                                        // Remove vertices from vertex list and update vertex-simplex connectivity
                bvector1D           &                                                 // (input) flags for removble veritces
                );
        void RemoveVertex(                                                        // Remove vertices from vertex list, update vertex-simplex connectivity and returns mapper
                bvector1D           &,                                                // (input) flags for removble veritces
                ivector1D           &                                                 // (input/output) mapper between old and new numeration
                );
        void CollapseDoubleVertex(                                                // Remove double vertices
                ivector1D           &,                                                // (input/output) list of collapsed vertices
                int                  a = 128                                          // (input/optional) number of bins used for vertex sorting
                );
        void CollapseDoubleVertex(                                                // Remove double vertices using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                ivector1D           &,                                                // (input/output) list of collapsed vertices
                int                  a = 128                                          // (input/optional) number of bins used for vertex sorting
                );
        void RemoveSimplex(                                                       // Remove simplicies from simplex list and update adjacencies and normals
                ivector1D           &                                                 // (input) List of removable simplicies
                );
        void RemoveSimplex(                                                       // Remove simplicies from simplex list, update adjacencies and normals and returns mapper
                ivector1D           &,                                                // (input) List of removable simplicies
                ivector1D           &                                                 // (input/output) Mapper between old and new numeration
                );
        void RemoveSimplex(                                                       // Remove simplicies from simplex list and update adjacencies and normals
                bvector1D           &                                                 // (input) Flags for removable simplicies
                );
        void RemoveSimplex(                                                       // Remove simplicies from simplex list, update adjacencies and normals and returns mapper
                bvector1D           &,                                                // (input) Flags for removable simplicies
                ivector1D           &                                                 // (input/output) Mapper between old and new numeration
                );
        void RemoveIsolatedVertex(                                                // Remove isolated vertices from tasselation
                void                                                                  // (input) none
                );
        void RemoveIsolatedVertex(                                                // Remove isolated vertices from tasselation and returns mapper
                ivector1D           &                                                 // (input/output) mapper between old and new numeration
                );
        void RemoveDoubleVertex(                                                  // Remove duplicated vertices from tasselation
                int                  a = 128                                          // (input/optional) number of bins for vertex sorting
                );
        void RemoveDoubleVertex(                                                  // Remove duplicated vertices from tasselation and returns mapper
                ivector1D           &,                                                // (input/output) Mapper between old and new numeration
                int                  a = 128                                          // (input/optional) number of bins for vertex sorting
                );
        void RemoveIsolatedSimplex(                                               // Remove isolated simplicies from the tasselation
                void                                                                  // (input) none
                );
        void RemoveIsolatedSimplex(                                               // Remove isolated simplicies from the tasselation and returns mapper
                ivector1D           &                                                 // (input/output) mapper between old and new numeration
                );
        void RemoveIsolatedSimplex(                                               // Remove isolated simplicies from the tasselation
                dvecarr3E           &                                                 // (input) External vertex list
                );
        void RemoveDoubleSimplex(                                                 // Remove duplicated simplicies from tasselation
                void                                                                  // (input) none
                );
        void RemoveDoubleSimplex(                                                 // Remove duplicated simplicies from tasselation and returns map between old and new numbering
                ivector1D           &                                                 // (input) map between old and new numeration
                );
        void RemoveDoubleSimplex(                                                 // Remove duplicated simplicies from tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        /*! Remove true duplicated simplicies. */	
        void RemoveTrueDoubleSimplex( 
                void			  
                );
        /*! Remove true duplicated simplicies and returns a mapper ivector1D with indices of removed simplicies. */	
        void RemoveTrueDoubleSimplex(				
                ivector1D &      
                );
        /*! Remove true duplicated simplicies, providing an external vertex list dvecarr3E */	
        void RemoveTrueDoubleSimplex(
                dvecarr3E &             
                );
        void Remove0AreaSimplex(                                                  // Remove 0-area simplicies from tasselation
                void                                                                  // (input) none
                );
        void Remove0AreaSimplex(                                                  // Remove 0-area simplicies from tasselation and returns map between old and new numbering
                ivector1D           &                                                 // (input/output) mapper between new and old numbering
                );
        void Remove0AreaSimplex(                                                  // Remove 0-area simplicies from tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        void Clean(                                                               // Clean surface tasselation
                void                                                                  // (input) none
                );
        void Clean(                                                               // Clean surface tasselation using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );

        // Stats Tools ---------------------------------------------------------- //
    public:
        void Stats(                                                               // Compute tasselation stats
                ostream             &                                                 // (input) output stream
                );
        void Stats(                                                               // Compute tasselation stats using an external vertex list
                ostream             &,                                                // (input) output stream
                dvecarr3E           &                                                 // (input) External vertex list
                );

        // Generation methods =================================================== //

        // Add elements --------------------------------------------------------- //
    public:
        void AddVertex(                                                           // Add a single vertex to vertex list
                darray3E            &                                                 // (input) vertex coordinates
                );
        void AddVertices(                                                         // Add a multiple vertices to tasselation
                dvecarr3E           &                                                 // (input) vertex coordinate list
                );
        void AddSimplex(                                                          // Add simplex to tasselation
                ivector1D           &                                                 // (input) simplex vertex coordinates
                );
        void AddSimplex(                                                          // Add simplex to tasselation
                dvecarr3E           &                                                 // (input) coordinates of each simplex vertex
                );
        void AddSimplicies(                                                       // Add multiple simplicies to tasselation
                ivector2D           &                                                 // (input) simplex-vertex connectivity
                );
        void AddSimplicies(                                                       // Add multiple simplicies to tasselation
                dvecarr3E           &,                                                // (input) vertex coordinate list
                ivector2D           &                                                 // (input) simplex-vertex connectivity
                );
        void SetNormal(                                                           // Set normal for a given simplex
                int                  ,                                                // (input) simplex global index
                darray3E            &                                                 // (input) normal vector
                );
        void SetAdjacency(                                                        // Set adjacencies for a given simplex
                int                  ,                                                // (input) simpex global index
                ivector2D           &                                                 // (input) adjacency list
                );
        void Append(                                                              // Append another tasselation
                Class_SurfTri       &                                                 // (input) source tasselation
                );

        // Transformations ------------------------------------------------------ //
    public:
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

        // Edges ---------------------------------------------------------------- //
    public:
        int edge(                                                                 // Returns the local index of the edge shared by two adjacent simplicies
                int                  ,                                                // (input) 1st simplex global index
                int                                                                   // (input) 2nd simplex global index
                );
        int edge(                                                                 // Returns the local index of the edge with specified vertices
                int                  ,                                                // (input) simplex global index
                int                  ,                                                // (input) global index of 1st vertex
                int                                                                   // (input) global index of 2nd vertex
                );

        // Vertices ------------------------------------------------------------- //
    public:
        int vertex(                                                               // Returns the local index of a given vertex
                int                  ,                                                // (input) simplex global index
                int                                                                   // (input) vertex global index
                );

        // Selection rings --------------------------------------------------------------- //
        /*! Returns the 1-ring of simplicies of a given vertex. Requires:
          int I: index of simplex INPUT
          int j: local vertex index on simplex I INPUT
          bool flag: check if ring is closed (true) or open (false). INPUT/OUTPUT
          */	 	
        ivector1D Ring_1(int, int, bool &);

        /*! Returns the 1-ring of simplicies of a given vertex. Requires:
          int I: index of simplex INPUT
          int j: local vertex index on simplex I INPUT
          bool flag: check if ring is closed (true) or open (false). INPUT/OUTPUT
          bool isRing: check if ring computation is successfull(true) or not (false) INPUT/OUTPUT
          */	 	
        ivector1D Ring_1(int, int, bool &, bool &); 
        /*! Returns the 1-ring of vertices of a given vertex. Requires:
          int I: index of simplex INPUT
          int j: local vertex index on simplex I INPUT
          bool flag: check if ring is closed (true) or open (false). INPUT/OUTPUT
          */	 	
        ivector2D VRing_1(int, int, bool &);                                               

        // Bounding boxes ------------------------------------------------------- //
    public:
        void BoundingBox(                                                         // Compute limits of tasselation bounding box (2D case)
                dvector1D           &,                                                // (input/output) limits in the x-direction
                dvector1D           &                                                 // (input/output) limits in the y-direction
                );
        void BoundingBox(                                                         // Compute limits of tasselation bounding box (3D case)
                dvector1D           &,                                                // (input/output) limits in the x-direction
                dvector1D           &,                                                // (input/output) limits in the y-direction
                dvector1D           &                                                 // (input/output) limits in the z-direction
                );
        void BoundingBox(                                                         // Compute limits of tasselation bounding box (2D case) using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                dvector1D           &,                                                // (input/output) limits in the x-direction
                dvector1D           &                                                 // (input/output) limits in the y-direction
                );
        void BoundingBox(                                                         // Compute limits of tasselation bounding box (3D case) using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                dvector1D           &,                                                // (input/output) limits in the x-direction
                dvector1D           &,                                                // (input/output) limits in the y-direction
                dvector1D           &                                                 // (input/output) limits in the z-direction
                );

        void BoundingBox(                                                         // Compute limits of tasselation bounding box (3D case)
                darray3E            &,                                            // (output) lower/left bounding point
                darray3E            &                                             // (output) upper/right bounding point
                );

            // Boundaries -------------------------------------------------------------------- //
        void Boundaries(                                                          // Returns mesh boundaries
                    Class_SurfTri       &
                    );                                              

        void Boundaries(                                                           // Returns mesh boundaries (external vertex list)
                dvecarr3E           &, 
                Class_SurfTri       &
                );                                  

        // Searching algorithms ---------------------------------------------------------- //
        int ReturnTriangleID(          // RETURN ID OF SIMPLEX WHICH A GIVEN POINT BELONGS TO
                darray3E  &,               // Point coordinates                                (INPUT)
                int                        // Global index of simplex used as seed             (INPUT)
                );

        // Patch decomposition -------------------------------------------------- //
        void FindPatch(                                                           // Find surface patch
                ivector1D           &,                                                // (input/output) patch index for each simplex
                double               ,                                                // (input) tolerance for sharp edge criterion
                int                  ,                                                // (input) original patch index
                int                  ,                                                // (input) next patch index
                int                 &                                                 // (input) global index of seed simplex
                );
        void FindConnectedLoop(                                                   // Find connected loop
                ivector1D           &,                                                // (input/output) loop index for each simplex
                int                  ,                                                // (input) original loop index
                int                  ,                                                // (input) next loop index
                int                 &                                                 // (input) global index of seed simplex
                );
        void ConnectedLoopDecompose(                                              // Decompose surface tasselation into connected loops
                int                 &,                                                // (input/output) number of connected loops
                ivector1D           &,                                                // (input/output) loop index associated to each simplex
                int                                                                   // (input) input of loop to be decomposed
                );
        void FindPatch_old(ivector1D &, double,int, int, int &);  // Find surface patch
        void PatchDecompose(int &, ivector1D &, int, double); // Patch decomposition
        void FindClosedLoop(ivector1D &,int, int, int &);     // Find closed loops
        void LoopDecompose(int &, ivector1D &, int);          // Decomposition into closed loop
        void FindRegion(                                      // Find bounded regions
                ivector1D &,
                ivector1D &,
                int,
                int,
                int &);
        void RegionDecompose(                                 // Decomposition into bounded region
                int &,
                ivector1D &,
                ivector1D &,
                int);
        void FindSegment(int                    ,             // Temporary
                ivector1D             &,
                double                 ,
                int                   &,
                vector<Class_SurfTri> &);

        // Interpolation tools =========================================================== //
    public:

        // Cell data --------------------------------------------------------------------- //
        template <class T>
            void InterpolateCellData();                                                        // Interpolate cell data
        template <class T>
            void CellData2PointData(vector<T> &, vector<T> &);                                 // Convert cell data into point data
        template <class T>
            void RemapCellData(ivector1D &, vector<T> &);                                      // Remap cell data

        // Point data -------------------------------------------------------------------- //
        template <class T>
            void InterpolatePointData(darray3E &, int, vector<T> &, T &);                     // Interpolate point data
        template <class T>
            void PointData2CellData(vector<T> &, vector<T> &);                                 // Convert point data to cell data
        template <class T>
            void RemapPointData(ivector1D &, vector<T> &);                                     // Remap point data

        // Check tools ========================================================== //
    public:
        darray3E Baricenter(                                                     // Compute simplex baricenter
                int                                                                   // (input) Simplex global index
                );
        darray3E Baricenter(                                                     // Compute simplex baricenter using an external vertex list
                int                  ,                                                // (input) Simplex global index
                dvecarr3E           &                                                 // (input) Exteranl vertex list
                );
        void minEdge(                                                             // Compute lenght of minimal edge over all tasselation simplicies
                double              &,                                                // (input/output) edge length
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) edge local index
                );
        void minEdge(                                                             // Compute lenght of minimal edge over all tasselation simplicies using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                double              &,                                                // (input/output) edge length
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) edge local index
                );
        void minEdge(                                                             // Compute lenght of minimal edge for a given simplex
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) edge length
                int                 &                                                 // (input/output) edge local index
                );
        void minEdge(                                                             // Compute lenght of minimal edge for a given simplex using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) edge length
                int                 &                                                 // (input/output) edge local index
                );
        void maxEdge(                                                             // Compute length of maximal edge over all tasselation simplicies
                double              &,                                                // (input/output) edge length
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) edge local index
                );
        void maxEdge(                                                             // Compute length of maximal edge over all tasselation simplicies using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex list
                double              &,                                                // (input/output) edge length
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) edge local index
                );
        void maxEdge(                                                             // Compute length of maximal edge for a given simplex
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) edge length
                int                 &                                                 // (input/output) edge local index
                );
        void maxEdge(                                                             // Compute length of maximal edge for a given simplex using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) edge length
                int                 &                                                 // (input/output) edge local index
                );
        double Area(                                                              // Compute tasselation area
                void                                                                  // (input) none
                );
        double Area(                                                              // Compute tasselation area using an external vertex list
                dvecarr3E           &                                                 // (input) External vertex list
                );
        double Area(                                                              // Compute simplex area using an external vertex list
                int                                                                   // (input) Simplex global index
                );
        double Area(                                                              // Compute simplex area using an external vertex list
                int                  ,                                                // (input) Simplex global index
                dvecarr3E           &                                                 // (input) External vertex list
                );
        void Area(                                                                // Compute area of each simplex in the surface tasselation
                dvector1D           &                                                 // (input/output) area of each simplex
                );
        void Area(                                                                // Compute area of each simplex in the surface tasselation using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                dvector1D           &                                                 // (input/output) area of each simplex
                );
        double AR(                                                                // Compute simplex aspect ratio
                int                                                                   // (input) Simplex global index
                );
        double AR(                                                                // Compute simplex aspect ratio using an external vertex list
                int                  ,                                                // (input) Simplex global index
                dvecarr3E           &                                                 // (input) External vertex list
                );
        void AR(                                                                  // Compute aspect ratio of each simplex in the tasselation
                dvector1D           &                                                 // (input/output) aspect ratio of each simplex
               );
        void AR(                                                                  // Compute aspect ratio of each simplex in the tasselation using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                dvector1D           &                                                 // (input/output) aspect ratio of each simplex
               );
        void Angle(                                                               // Compute the angle between edges incident on simplex vertex
                int                  ,                                                // (input) simplex global index
                double              &,                                                // (input/output) angle value
                int                                                                   // (input) vertex local index
                );
        void Angle(                                                               // Compute the angle between edges incident on simplex vertex using an external vertex list
                dvecarr3E           &,                                                // (input) external vertex coordinate list
                int                  ,                                                // (input) simplex global index
                double              &,                                                // (input/output) angle value
                int                                                                   // (input) vertex local index
                );
        void minAngle(                                                            // Compute min-angle over all simplicies in the tasselation
                double              &,                                                // (input/output) min angle value
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) vertex local index
                );
        void minAngle(                                                            // Compute min-angle over all simplicies in the tasselation using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                double              &,                                                // (input/output) min angle value
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) vertex local index
                );
        void minAngle(                                                            // Compute simplex min-angle
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) min angle value
                int                 &                                                 // (input/output) vertex local index
                );
        void minAngle(                                                            // Compute simplex min-angle using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) min angle value
                int                 &                                                 // (input/output) vertex local index
                );
        void minAngle(                                                            // Compute min angle for each simplex in the tasselation
                dvector1D           &                                                 // (input/output) min angle for each simplex
                );
        void minAngle(                                                            // Compute min angle for each simplex in the tasselation
                dvecarr3E           &,                                                // (input) External vertex list
                dvector1D           &                                                 // (input/output) min angle for each simplex
                );
        void maxAngle(                                                            // Compute max angle over all simplicies
                double              &,                                                // (input/output) simplex max angle
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) vertex local index
                );
        void maxAngle(                                                            // Compute max angle over all simplicies using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                double              &,                                                // (input/output) simplex max angle
                int                 &,                                                // (input/output) Simplex global index
                int                 &                                                 // (input/output) vertex local index
                );
        void maxAngle(                                                            // Compute simplex max angle
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) simplex max angle
                int                 &                                                 // (input/output) vertex local index
                );
        void maxAngle(                                                            // Compute simplex max angle using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                int                  ,                                                // (input) Simplex global index
                double              &,                                                // (input/output) simplex max angle
                int                 &                                                 // (input/output) vertex local index
                );
        void maxAngle(                                                            // Compute max angle for each simplex
                dvector1D           &                                                 // (input/output) max angle for each simplex
                );
        void maxAngle(                                                            // Compute max angle for each simplex using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                dvector1D           &                                                 // (input/output) max angle for each simplex
                );
        darray3E  Edge_midPoint(                                                  // Compute midpoint of simplex edge
                int                  ,                                                // (input) Simplex global index
                int                                                                   // (input) edge local index
                );
        darray3E   Edge_midPoint(                                                  // Compute midpoint of simplex edge using an external vertex list
                dvecarr3E           &,                                                // (input) External vertex list
                int                  ,                                                // (input) Simplex global index
                int                                                                   // (input) edge local index
                );
        darray3E  CircumCenter(                                                   // COMPUTE SIMPLEX CIRCUMCENTER
                int                                                                   // Simplex global index                             (INPUT)
                );
        darray3E  CircumCenter(                                                   // COMPUTE SIMPLEX CIRCUMCENTER USING AN EXTERNAL VERTEX LIST
                int,                                                                  // Simplex global index                             (INPUT)
                dvecarr3E &                                                           // vertex list                                      (INPUT)
                );
        int VertexValence(                                                        // Compute vertex valence
                int                  ,                                                // (input) Global index of simplex, which vertex belongs to
                int                                                                   // (input) local index of vertex
                );
        void VertexValence(                                                       // Compute valence of each vertex
                ivector1D           &                                                 // (input/output) Valence score for each vertex
                );

        // Refinement tools ===================================================== //
    public:

        // Simplex refinement --------------------------------------------------- //
        void invert_loc_num(                                                      // Invert local numbering of a given simplex
                int                                                                   // (input) simplex global index
                );
        void BinaryRefinement(                                                    // Binary refinement of tasselation
                double                                                                // (input) threshold for binary refinement
                );
        void BinaryRefinement(                                                    // Binary refinement of tasselation using an external vertex list
                dvecarr3E           &,                                                // (input/output) vertex coordinate list
                double                                                                // (input) threshold for binary refinement
                );
        void split_1segm2segm(                                                    // Split 2-simplex into two 2-simplcies at specified point
                int                  ,                                                // (input) simplex global index
                darray3E            &                                                 // (input) point coordinates
                );
        void split_1segm2segm(                                                    // Split 2-simplex into two 2-simplcies at specified point using an external vertex list
                int                  ,                                                // (input) simplex global index
                int                  ,                                                // (input) global index of vertex used to split segment
                dvecarr3E           &                                                 // (input) external vertex list
                );
        void Collapse_2Simplex(                                                   // Collapse 2 simplex
                int , 
                int
                );                                        
        void Split_2Simplex(                                                      // Split a 2-simplex at mid-point
                int                  T                                                // (input) simplex global index
                );
        void Split_2Simplex(                                                      // Split a 2-simplex at mid-point using an external vertex list
                dvecarr3E           &,                                                // (input/output) vertex coordinate list
                int                  T                                                // (input) simplex global index
                );
        void SplitEdge(                                                           // Split edge at mid-point
                int , 
                int 
                );                 
        void CollapseEdge(                                                        // Collapse edge
                int , 
                int , 
                int
                );      
    private:

        bool PointsOnSameSide(
                array< double, 3 > const &P1,
                array< double, 3 > const &P2,
                array< double, 3 > const &A,
                array< double, 3 > const &B
                ); 

        array<double, 3> IntersectLines(
                array<double, 3> const &n1,
                array<double, 3> const &P1,
                array<double, 3> const &n2,
                array<double, 3> const &P2
                ) ;


    public:

        // Voronoi diagrams -------------------------------------------------------------- //
        void Voronoi(Class_SurfTri &);                // Compute approximated Voronoi tasselation
        void Voronoi(Class_SurfTri &, dvector2D &);

        // IO functions ========================================================== //
    public:

        // dgf ------------------------------------------------------------------- //
        void Import_dgf(                                                           // Import tasselation from a .dgf file
                string                                                                 // (input) dgf file name
                );
        void Export_dgf(                                                           // Export tasselation to a .dgf file
                string                                                                 // (input) dgf file name
                );

        // stl --------------------------------------------------------------------------- //
        void Import_stl(                                                           // Import tasselation from a .stl file
                string ,                                                               // (input) stl file name
                bool                                                                   // (input) stl file type (ASCII/binary)
                );
        void Export_stl(                                                           // Export tasselation tp a .stl file
                string ,                                                               // (input) stl file name
                bool                                                                   // (input) stl file type (ASCII/binary)
                );

};

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "Class_SurfTri.tpp"

#endif
