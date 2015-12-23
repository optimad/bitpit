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
# ifndef __UCARTMESH_HPP__
# define __UCARTMESH_HPP__

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
# include "Operators.hpp"
# include "DGF_IOFunct.hpp"
# include "VTK_IOFunct.hpp"

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

typedef array<int,3>                   iarray3E;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

typedef array<double,3>                darray3E;
typedef vector<darray3E>               dvecar3E;

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

class UCartMesh{

    // Members ============================================================== //
    private:

        int                 dim;
        int                 nCells ;
        int                 nPoints ;
        int                 CellsInIJPlane;
        int                 PointsInIJPlane;

        darray3E            B0;                                 // left/lower/front limit point
        darray3E            B1;                                 // right/upper/rear limit point

        iarray3E            nc;                                 // Number of cells in each direction
        iarray3E            np;                                 // Number of points in each direction
        darray3E            h;                                  // grid spacing in each direction

        dvector2D           center;                             // Cell's centers coordinates
        dvector2D           edge;                               // Cell's edges coordinates

        // Constructors ========================================================= //
    public:

        UCartMesh(                                              // Default constructor for Class_UCartMesh2D
                void                                            // (input) none
                );

        UCartMesh(                                              // Custom constructor #1 for Class_UCartMesh2D
                darray3E       &,                               // (input) lower/left limit point
                darray3E       &,                               // (input) upper/right limit point
                iarray3E       &,                               // (input) number of mesh cells in each direction
                int                                            // (input) dimensions
                );

        // Destructors ========================================================== //
        ~UCartMesh(                                             // Default destructor for Class_UCartMesh2D
                void
                );

        // Assignament operators ================================================ //
        UCartMesh& operator=(                                   // Default assignament operator
                const UCartMesh &                               // Source mesh
                );

        // Methods ============================================================== //
    private:
        void ResizeMesh(                                        // Resize mesh data structure according to mesh size
                void                                            // (input) none
                );

        // Initialization ------------------------------------------------------- //
    public:
        int  getNCells();
        int  getNCells( int );

        int  getNPoints();
        int  getNPoints( int d);

        void getBoundingBox(
                darray3E    &,
                darray3E    &
                );

        darray3E  getCenter(
                int
                );

        darray3E  getCenter(
                iarray3E 
                );

        darray3E  getCenter(
                int             ,
                int             ,
                int k=0            
                );

        darray3E  getPoint(
                int
                );

        darray3E  getPoint(
                iarray3E 
                );

        darray3E  getPoint(
                int             ,
                int             ,
                int   k=0          
                );

        darray3E  getFirstPointOfCell(
                int
                );

        darray3E  getLastPointOfCell(
                int
                );

        void SetMesh(                                           // Generate mesh
                darray3E       &,                               // (input) lower/left limit point
                darray3E       &,                               // (input) upper/right limit point
                iarray3E       &,                               // (input) number of mesh cells in each direction
                int             
                );

        void ClearMesh(                                         // Destroy current mesh
                void                                            // (input) none
                );


        // Transformations ------------------------------------------------------ //
        void Translate(                                         // Translate mesh along x, y axis
                darray3E       &                                // (input) displacement for each direction
                );

        void Scale(                                             // Scale mesh
                darray3E       &,                               // (input) scaling factor along each dimension
                darray3E                                        // (input) origin 
                );

        void Scale(                                             // Scale mesh with origin in B0
                darray3E       &                                // (input) scaling factor along each dimension
                );


        // Format conversion ---------------------------------------------------- //
        void Cart2Unstr(                                     // Convertes cartesian mesh to unstructured surface mesh
                int             &,                              // (input/output) number of mesh vertices
                int             &,                              // (input/output) number of simplicies
                dvecarr3E       &,                              // (input/output) vertex coordinate list
                ivector2D       &,                              // (input/output) simplex-vertex connectivity
                ivector3D       &                               // (input/output) simplex-simplex adjacency
                );

        // Mapping -------------------------------------------------------------- //
        iarray3E CellCartesianId(                              // Return cartesian cell id which encloses a given point
                const darray3E        &                               // (input) point coordinates
                );

        iarray3E CellCartesianId(                                // Returns cartesian indices of cell given its global index
                int                                             // (input) cell global index
                );

        int CellLinearId(                                    // Return linear cell id which encloses a given point
                const darray3E        &                               // (input) point coordinates
                );

        int CellLinearId(                                         // Returns the index of cell given its cartesian indices
                iarray3E                                        // (input) 1st cartesian indices
                );

        int CellLinearId(                                         // Returns the index of cell given its cartesian indices
                int              ,                              // (input) 1st cartesian indices
                int              ,                              // (input) 1st cartesian indices
                int  k=0                                        // (input) 2nd cartesian indices
                );

        iarray3E PointCartesianId(                              // Return cartesian id of closest point 
                const darray3E        &                               // (input) point coordinates
                );

        iarray3E PointCartesianId(                               // Returns cartesian indices of a vertex given its global index
                int                                             // (input) vertex global index
                );

        int PointLinearId(                                    // Return linear id of closest point
                const darray3E        &                               // (input) point coordinates
                );

        int PointLinearId(                                         // Returns the index of cell given its cartesian indices
                iarray3E                                        // (input) 1st cartesian indices
                );

        int PointLinearId(                                        // Returns the index of vertex given its cartesian indices
                int              ,                              // (input) 1st cartesian index
                int              ,                              // (input) 2nd cartesian index
                int  k=0                                        // (input) 2nd cartesian index
                );

                );

                );

        // Interpolation -------------------------------------------------------- //
        //
        bool PointInGrid(                                                              // Point-in-triangle condition
                array<double, 3> const              &                                   // (input) Point coordinatesertex
                );

        template <class T>
            void CellData2PointData(                            // Convertes cell data into point data
                    vector< T >     &,                          // (input) cell data
                    vector< T >     &                           // (input/output) point data
                    );
        template < class T, typename ... T2 >
            void CellData2PointData(                            // Convertes arbitrary cell dataset into point data
                    vector< T >     &,                          // (input) cell data
                    vector< T >     &,                          // (input/output) point data
                    T2              &...                        // (input/optional) other cell data set to be converted
                    );
        template <class T>
            void PointData2CellData(                            // Convertes point data into cell data
                    vector< T >     &,                          // (input) point data
                    vector< T >     &                           // (input/output) cell data
                    );
        template < class T, typename ... T2 >
            void PointData2CellData(                            // Convertes arbitrary point dataset into cell data
                    vector< T >     &,                          // (input) cell data
                    vector< T >     &,                          // (input/output) point data
                    T2              &...                        // (input/optional) other cell data set to be converted
                    );
        template < class T >
            void interpolateCellData(                           // Interpolate cell data at a specified point
                    darray3E        &,                          // (input) point coordinates
                    vector< T >     &,                          // (input) cell data
                    T               &                           // (input/output) interpolation result
                    );
        template < class T, typename ... T2 >
            void interpolateCellData(                           // Interpolate cell datasets at a specified point
                    darray3E       &,                          // (input) point coordinates
                    vector< T >     &,                          // (input) cell data
                    T               &,                          // (input/output) interpolation result
                    T2          &...                            // (input/optional) other dataset to be used for interpolation
                    );
        template <class T>
            void interpolatePointData(                          // Interpolate point datasets at a specified point
                    darray3E       &,                          // (input) point coordinates
                    vector< T >     &,                          // (input) point data
                    T               &                           // (input/output) interpolation result
                    );

        template < class T, typename ... T2 >
            void interpolatePointData(                          // Interpolate vertex datasets at a specified point
                    darray3E       &,                          // (input) point coordinates
                    vector< T >     &,                          // (input) point data
                    T               &,                          // (input/output) interpolation result
                    T2          &...                            // (input/optional) other dataset to be used for interpolation
                    );

        // I/O methods ====================================================================== //

        // Paraview ------------------------------------------------------------------------- //
        void Export_vtr(string );
        template <class T>
            void Export_CellData_vtr(string , string , vector<T> &);
        template <class T>
            void Export_PointData_vtr(string , string , vector<T> &);


};

#include"UCartMesh.tpp"

#endif
