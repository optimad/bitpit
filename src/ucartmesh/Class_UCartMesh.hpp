// ========================================================================== //
//                 - GRID MANAGER FOR CARTESIAN MESHES -                      //
//                                                                            //
// Grid manager for cartesian meshes                                          //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v2.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __CLASS_UCARTMESH_HPP__
# define __CLASS_UCARTMESH_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard template library
# include <array>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>

// CC_Lib
# include <bitpit_operators.hpp>
# include <bitpit_IO.hpp>

// ========================================================================== //
// NAMESPACES                                                                 //
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

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

// ========================================================================== //
// FUNCTION PROTOTYPES                                                        //
// ========================================================================== //
class Class_UCartMesh2D {

    // Members ============================================================== //
    public:

    dvector1D      xlim;                                                      // Mesh limits in the x direction
    dvector1D      ylim;                                                      // Mesh limits in the y direction
    int            nx, ny;                                                    // Number of cells in the x, y direction
    double         dx, dy;                                                    // Mesh spacing in the x, y direction
    dvector1D      xnode, ynode;                                              // Cell's centers coordinates
    dvector1D      xedge, yedge;                                              // Cell's edges coordinates

    // Constructors ========================================================= //
    public:

    Class_UCartMesh2D(                                                        // Default constructor for Class_UCartMesh2D
        void                                                                  // (input) none
    );
    Class_UCartMesh2D(                                                        // Custom constructor #1 for Class_UCartMesh2D
        dvector1D       &,                                                    // (input) mesh limits in the x direction
        dvector1D       &,                                                    // (input) mesh limits in the y direction
        int              ,                                                    // (input) number of mesh cells in the x direction
        int                                                                   // (input) number of mesh cells in the y direction
    );

    // Destructors ========================================================== //
    public:

    ~Class_UCartMesh2D(                                                       // Default destructor for Class_UCartMesh2D
        void
    );

    // Assignament operators ================================================ //
    public: 

    Class_UCartMesh2D& operator=(                                             // Default assignament operator
        const Class_UCartMesh2D &                                             // Source mesh
    );

    // Methods ============================================================== //

    // Initialization ------------------------------------------------------- //
    public:

    void ClearMesh(                                                           // Destroy current mesh
        void                                                                  // (input) none
    );
    void ResizeMesh(                                                          // Resize mesh data structure according to mesh size
        void                                                                  // (input) none
    );
    void SetMesh(                                                             // Generate mesh
        dvector1D       &,                                                    // (input) mesh limits in the x-direction
        dvector1D       &,                                                    // (input) mesh limits in the y-direction
        int              ,                                                    // (input) number of mesh cells in the x-direction
        int                                                                   // (input) number of mesh cells in the y-direction
    );

    // Transformations ------------------------------------------------------ //
    void Translate(                                                           // Translate mesh along x, y axis
        dvector1D       &                                                     // (input) displacement for each direction
    );
    void Scale(                                                               // Scale mesh along x, y axis
        dvector1D       &                                                     // (input) scaling factor along each dimension
    );

    // Format conversion ---------------------------------------------------- //
    void Cart2SurfMesh(                                                       // Convertes cartesian mesh to unstructured surface mesh
        int             &,                                                    // (input/output) number of mesh vertices
        int             &,                                                    // (input/output) number of simplicies
        dvector2D       &,                                                    // (input/output) vertex coordinate list
        ivector2D       &,                                                    // (input/output) simplex-vertex connectivity
        ivector3D       &                                                     // (input/output) simplex-simplex adjacency
    );

    // Mapping -------------------------------------------------------------- //
    void ReturnCellID(                                                        // Return cell id which encloses a given point
        dvector1D       &,                                                    // (input) point coordinates
        int             &,                                                    // (input/output) 1st cartesian index
        int             &                                                     // (input/output) 2nd cartesian index
    );
    int AccessCellData(                                                       // Returns the index of cell given its cartesian indices
        int              ,                                                    // (input) 1st cartesian indices
        int                                                                   // (input) 2nd cartesian indices
    );
    void AccessCellData(                                                      // Returns cartesian indices of cell given its global index
        int              ,                                                    // (input) cell global index
        int             &,                                                    // (input/output) 1st cartesian indices
        int             &                                                     // (input/output) 2nd cartesian indices
    );
    int AccessPointData(                                                      // Returns the index of vertex given its cartesian indices
        int              ,                                                    // (input) 1st cartesian index
        int                                                                   // (input) 2nd cartesian index
    );
    void AccessPointData(                                                     // Returns cartesian indices of a vertex given its global index
        int              ,                                                    // (input) vertex global index
        int             &,                                                    // (input/output) 1st cartesian indices
        int             &                                                     // (input/output) 2nd cartesian indices
    );

    // Interpolation -------------------------------------------------------- //
    template <class T>
    void CellData2PointData(                                                  // Convertes cell data into point data
        vector< T >     &,                                                    // (input) cell data
        vector< T >     &                                                     // (input/output) point data
    );
    template < class T, typename ... T2 >
    void CellData2PointData(                                                  // Convertes arbitrary cell dataset into point data
        vector< T >     &,                                                    // (input) cell data
        vector< T >     &,                                                    // (input/output) point data
        T2              &...                                                  // (input/optional) other cell data set to be converted
    );
    template <class T>
    void PointData2CellData(                                                  // Convertes point data into cell data
        vector< T >     &,                                                    // (input) point data
        vector< T >     &                                                     // (input/output) cell data
    );
    template < class T, typename ... T2 >
    void PointData2CellData(                                                  // Convertes arbitrary point dataset into cell data
        vector< T >     &,                                                    // (input) cell data
        vector< T >     &,                                                    // (input/output) point data
        T2              &...                                                  // (input/optional) other cell data set to be converted
    );
    template < class T >
    void interpolateCellData(                                                 // Interpolate cell data at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) cell data
        T               &                                                     // (input/output) interpolation result
    );
    template < class T, typename ... T2 >
    void interpolateCellData(                                                 // Interpolate cell datasets at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) cell data
        T               &,                                                    // (input/output) interpolation result
        T2          &...                                                      // (input/optional) other dataset to be used for interpolation
    );
    template <class T>
    void interpolatePointData(                                                // Interpolate point datasets at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) point data
        T               &                                                     // (input/output) interpolation result
    );
    template < class T, typename ... T2 >
    void interpolatePointData(                                                // Interpolate vertex datasets at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) point data
        T               &,                                                    // (input/output) interpolation result
        T2          &...                                                      // (input/optional) other dataset to be used for interpolation
    );

        // I/O methods ====================================================================== //

};

class Class_UCartMesh3D {

    // Members ============================================================== //
    public:

    dvector1D      xlim;                                                      // Mesh limits in the x direction
    dvector1D      ylim;                                                      // Mesh limits in the y direction
    dvector1D      zlim;                                                      // Mesh limits in the z direction
    int            nx, ny, nz;                                                // Number of cells in the x, y, z direction
    double         dx, dy, dz;                                                // Mesh spacing in the x, y, z direction
    dvector1D      xnode, ynode, znode;                                       // Cell's centers coordinates
    dvector1D      xedge, yedge, zedge;                                       // Cell's edges coordinates

    // Constructors ========================================================= //
    public:

    Class_UCartMesh3D(                                                        // Default constructor for Class_UCartMesh2D
        void                                                                  // (input) none
    );
    Class_UCartMesh3D(                                                        // Custom constructor #1 for Class_UCartMesh2D
        dvector1D       &,                                                    // (input) mesh limits in the x direction
        dvector1D       &,                                                    // (input) mesh limits in the y direction
        dvector1D       &,                                                    // (input) mesh limits in the z direction
        int              ,                                                    // (input) number of mesh cells in the x direction
        int              ,                                                    // (input) number of mesh cells in the y direction
        int                                                                   // (input) number of mesh cells in the z direction
    );

    // Destructors ========================================================== //
    public:

    ~Class_UCartMesh3D(                                                       // Default destructor for Class_UCartMesh2D
        void
    );

    // Assignament operators ================================================ //
    public: 

    Class_UCartMesh3D& operator=(                                             // Default assignament operator
        const Class_UCartMesh3D &                                             // Source mesh
    );

    // Methods ============================================================== //

    // Initialization ------------------------------------------------------- //
    public:

    void ClearMesh(                                                           // Destroy current mesh
        void                                                                  // (input) none
    );
    void ResizeMesh(                                                          // Resize mesh data structure according to mesh size
        void                                                                  // (input) none
    );
    void SetMesh(                                                             // Generate mesh
        dvector1D       &,                                                    // (input) mesh limits in the x-direction
        dvector1D       &,                                                    // (input) mesh limits in the y-direction
        dvector1D       &,                                                    // (input) mesh limits in the z-direction
        int              ,                                                    // (input) number of mesh cells in the x-direction
        int              ,                                                    // (input) number of mesh cells in the y-direction
        int                                                                   // (input) number of mesh cells in the z-direction
    );

    // Transformations ------------------------------------------------------ //
    void Translate(                                                           // Translate mesh along x, y axis
        dvector1D       &                                                     // (input) displacement for each direction
    );
    void Scale(                                                               // Scale mesh along x, y axis
        dvector1D       &                                                     // (input) scaling factor along each dimension
    );

    // Format conversion ---------------------------------------------------- //

    // Mapping -------------------------------------------------------------- //
    void ReturnCellID(                                                        // Return cell id which encloses a given point
        dvector1D       &,                                                    // (input) point coordinates
        int             &,                                                    // (input/output) 1st cartesian index
        int             &,                                                    // (input/output) 2nd cartesian index
        int             &                                                     // (input/output) 3rd cartesian index
    );
    int AccessCellData(                                                       // Returns the index of cell given its cartesian indices
        int              ,                                                    // (input) 1st cartesian indices
        int              ,                                                    // (input) 2nd cartesian indices
        int                                                                   // (input) 3rd cartesian indices
    );
    void AccessCellData(                                                      // Returns cartesian indices of a cell given its global index
        int              ,                                                    // (input) cell global index
        int             &,                                                    // (input/output) 1st cartesian indices
        int             &,                                                    // (input/output) 2nd cartesian indices
        int             &                                                     // (input/output) 3rd cartesian indices
    );
    int AccessPointData(                                                      // Returns the index of vertex given its cartesian indices
        int              ,                                                    // (input) 1st cartesian index
        int              ,                                                    // (input) 2nd cartesian index
        int                                                                   // (input) 3rd cartesian index
    );
    void AccessPointData(                                                     // Returns cartesian indices of a vertex given its global index
        int              ,                                                    // (input) vertex global index
        int             &,                                                    // (input/output) 1st cartesian indices
        int             &,                                                    // (input/output) 2nd cartesian indices
        int             &                                                     // (input/output) 3rd cartesian indices
    );

    // Interpolation -------------------------------------------------------- //
    template <class T>
    void CellData2PointData(                                                  // Convertes cell data into point data
        vector< T >     &,                                                    // (input) cell data
        vector< T >     &                                                     // (input/output) point data
    );
    template < class T, typename ... T2 >
    void CellData2PointData(                                                  // Convertes arbitrary cell dataset into point data
        vector< T >     &,                                                    // (input) cell data
        vector< T >     &,                                                    // (input/output) point data
        T2              &...                                                  // (input/optional) other cell data set to be converted
    );
    template <class T>
    void PointData2CellData(                                                  // Convertes point data into cell data
        vector< T >     &,                                                    // (input) point data
        vector< T >     &                                                     // (input/output) cell data
    );
    template < class T, typename ... T2 >
    void PointData2CellData(                                                  // Convertes arbitrary point dataset into cell data
        vector< T >     &,                                                    // (input) cell data
        vector< T >     &,                                                    // (input/output) point data
        T2              &...                                                  // (input/optional) other cell data set to be converted
    );
    template < class T >
    void interpolateCellData(                                                 // Interpolate cell data at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) cell data
        T               &                                                     // (input/output) interpolation result
    );
    template < class T, typename ... T2 >
    void interpolateCellData(                                                 // Interpolate cell datasets at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) cell data
        T               &,                                                    // (input/output) interpolation result
        T2          &...                                                      // (input/optional) other dataset to be used for interpolation
    );
    template <class T>
    void interpolatePointData(                                                // Interpolate point datasets at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) point data
        T               &                                                     // (input/output) interpolation result
    );
    template < class T, typename ... T2 >
    void interpolatePointData(                                                // Interpolate vertex datasets at a specified point
        dvector1D       &,                                                    // (input) point coordinates
        vector< T >     &,                                                    // (input) point data
        T               &,                                                    // (input/output) interpolation result
        T2          &...                                                      // (input/optional) other dataset to be used for interpolation
    );

        // I/O methods ====================================================================== //

};

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
# include "Class_UCartMesh.tpp"

# endif
