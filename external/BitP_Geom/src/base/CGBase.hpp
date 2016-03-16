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
# ifndef __CGBASE_HPP__
# define __CGBASE_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <string>
# include <iostream>

// bitpit library
# include <bitpit_operators.hpp>



// Basic computational geometry ============================================= //

namespace CGElem{

typedef std::array<double,3> array3D ;

static const std::array< std::array<int,2>,12> BoxEdgeVertexConnectivity =
{
    std::array<int,2>{ {0,2} },
    std::array<int,2>{ {1,3} },
    std::array<int,2>{ {0,1} },
    std::array<int,2>{ {2,3} },
    std::array<int,2>{ {4,6} },
    std::array<int,2>{ {5,7} },
    std::array<int,2>{ {4,5} },
    std::array<int,2>{ {6,7} },
    std::array<int,2>{ {0,4} },
    std::array<int,2>{ {1,5} },
    std::array<int,2>{ {2,6} },
    std::array<int,2>{ {3,7} }

};

static const std::array< std::array<int,4>, 6> BoxFaceVertexConnectivity =
{
    std::array<int,4>{ {0,2,6,4} },
    std::array<int,4>{ {1,3,7,5} },
    std::array<int,4>{ {0,4,5,1} },
    std::array<int,4>{ {2,6,7,3} },
    std::array<int,4>{ {0,1,3,2} },
    std::array<int,4>{ {4,5,7,6} }
};


double              distancePointLine( array3D const &, array3D const &, array3D const &, array3D & ) ;
double              distancePointPlane( array3D const &, array3D const &, array3D const &, array3D & ) ;
double              distancePointSegment( array3D const &, array3D const &, array3D const &, array3D &, int & ) ;
double              distancePointTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, int & ) ;
double              distancePointSimplex( array3D const &, std::vector<array3D> const &, array3D &, int & );

std::vector<double> distanceCloudTriangle( std::vector<array3D> const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> &, std::vector<int> & );
std::vector<double> distanceCloudSimplex( std::vector<array3D> const &, std::vector<array3D> const &, std::vector<array3D> &, std::vector<int> & );


bool                intersectPointSegment( array3D const &, array3D const &, array3D const & ) ;
bool                intersectPointTriangle( array3D const &, array3D const &, array3D const &, array3D const & ) ;
bool                intersectPointBox( array3D const &, array3D const &, array3D const &, int dim=3 ) ;

bool                intersectLineLine( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool                intersectLinePlane( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool                intersectLineTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool                intersectLineSimplex( array3D const &, array3D const &, std::vector<array3D> const &, array3D & ) ;

bool                intersectSegmentSegment( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool                intersectSegmentPlane( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool                intersectSegmentTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool                intersectSegmentSimplex( array3D const &, array3D const &, std::vector<array3D> const &, array3D & ) ;
bool                intersectSegmentBox( array3D const &, array3D const &, array3D const &, array3D const &, int  dim = 3 ) ;
bool                intersectSegmentBox( array3D const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> &, int  dim = 3 ) ;

bool                intersectPlanePlane( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & ) ;

bool                intersectBoxBox( array3D const &, array3D const &, array3D const &, array3D const &, int  dim = 3 ) ;
bool                intersectBoxBox( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D &, int  dim = 3 ) ;
bool                intersectBoxTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const & ) ;
bool                intersectBoxTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> & ) ;
bool                intersectBoxSimplex( array3D const &, array3D const &, std::vector<array3D> const &, int dim=3 ) ;
bool                intersectBoxSimplex( array3D const &, array3D const &, std::vector<array3D> const &, std::vector<array3D> &, int dim=3 ) ;

void                computeAABBSegment( array3D const &, array3D const &, array3D &, array3D & ) ;
void                computeAABBTriangle( array3D const &, array3D const &, array3D const &, array3D &, array3D & ) ;
void                computeAABBSimplex( std::vector<array3D> const &, array3D &, array3D & ) ;

void                unionAABB( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & );
void                unionAABB( std::vector<array3D>  const &, std::vector<array3D> const &, array3D &, array3D & );

//levelset    bool IntersectLineSurface(
//levelset            array3D  const  &,                                                // (input)  first point on surface
//levelset            array3D  const  &,                                                // (input)  normal of first point on surface
//levelset            array3D  const  &,                                                // (input)  second point on surface
//levelset            array3D  const  &,                                                // (input)  normal of second point on surface
//levelset            array3D  const  &,                                                // (input)  point on line
//levelset            array3D  const  &,                                                // (input)  versor of line
//levelset            array3D         &,                                                // (output)  intersection of line and surface
//levelset            array3D         &                                                 // (output)  normal at intersection
//levelset            );


void RotateVector(
        std::array<double, 3>                    &,
        std::array<double, 3> const              &,
        double                                     
        );


void VertexOfBox(
        int              const              &,     
        std::array<double, 3> const              &,
        std::array<double, 3> const              &,
        std::array<double, 3>                    & 
        );

void EdgeOfBox(
        int              const              &,     
        std::array<double, 3> const              &,
        std::array<double, 3> const              &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    & 
        );

void FaceOfBox(
        int              const              &,     
        std::array<double, 3> const              &,
        std::array<double, 3> const              &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    & 
        );





};


// Polygon ================================================================== //
namespace CGPolygon2D {

class Class_CG_Polygon2D {

    // Members ---------------------------------------------------------- //
    public:
        struct CG_Polygon2D_info{
            int                         status;                               // flag for entity status
            bool                        is_clockwise;                         // flag for clockwise polygon
            bool                        is_convex;                            // flag for convex polygon
            std::array<double, 2>            xlim, ylim;                           // polygon bound. box limits (x direction)                   
        } info;                                                               // data structure for polygon infos
        int                             n_faces;                              // number of faces
        std::vector<std::array<double,3>>                     v_list;                               // vertex list

        // Constructors ----------------------------------------------------- //
    public:
        Class_CG_Polygon2D(                                                   // Default constructor for Class_CG_Polygon
                void                                                              // (input) none
                );
        Class_CG_Polygon2D(                                                   // Custom constructor #1 for Class_CG_Polygon
                std::vector<std::array<double,3>>                   &                                     // (input) vertex list
                );
        Class_CG_Polygon2D(                                                   // Custom constructor #2 for Class_CG_Polygon
                std::vector<std::array<double,3>>                   &,                                    // (input) vertex coordinate list
                std::vector<std::vector<int>>                   &,                                    // (input) edge-vertex connectivity
                std::vector<std::vector<std::vector<int>>>                   &                                     // (input) edge-edge adjacencies
                );

        // Destructors ------------------------------------------------------ //
    public:
        ~Class_CG_Polygon2D(                                                  // Default destructor for Class_CG_Polygon
                void                                                              // (input) none
                );

        // Assignament operator --------------------------------------------- //
    public:
        Class_CG_Polygon2D & operator=(                                       // Assignament operator
                const Class_CG_Polygon2D    &                                     // (input) source variable
                );

        // Methods ---------------------------------------------------------- //
    public:
        void AddVertex(                                                       // Add vertex to polygon vertex list
                std::array<double,3>                     &                                     // (input) vertex coordinates
                );
        void AddVertices(                                                     // Add multiple vertices to polygon vertex list
                std::vector<std::array<double,3>>                   &                                     // (input) vertex list
                );
        void IsClockWise(                                                     // Check wheter the polygon is clockwise or not
                void                                                              // (input) none
                );
        void IsConvex(                                                        // Check wheter the polygon is convex or not
                void                                                              // (input) none
                );
        void BoundingBox(                                                     // Compute polygon bounding box
                void                                                              // (input) none
                );
        bool ConvexIncludePoint(                                              // Check wheter the polygon includes a given point or not (convex case)
                std::array<double,3>                     &                                     // (input) point coordinates
                );
        bool NonConvexIncludePoint(                                           // Check wheter the polygon includes a given point or not (non-convex case)
                std::array<double,3>                     &                                     // (input) point coordinates
                );
        bool IncludePoint(                                                    // Check wheter the polygon includes a given point or not
                std::array<double,3>                     &                                     // (input) point coordinates
                );
        void display(                                                         // Display polygon infos
                std::ostream                     &                                     // (input/output) output stream
                );
};
};

// Surfs ==================================================================== //
namespace CGPLSurf {
bool IsClosed(                                                            // Check if a given surf is closed
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) simplex-simplex adjacency
        );
bool IsOpen(                                                              // Check if a given surf is open
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) simplex-simplex adjacency
        );
bool IsManifold(                                                          // Check if a given surf is manifold
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) simplex-simplex adjacency
        );
void maxCurvature(                                                        // Compute max discrete curvature (simple surfs only)
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex normals
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<double>                       &                                     // (input/output) curvature at curve vertex
        );
};

// Curves =================================================================== //
namespace CGPLCurve {
bool IsClosed(                                                            // Check if a given curve is closed
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) segment-segment adjacency
        );
bool IsOpen(                                                              // Check if a given curve is open
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) segment-segment adjacency
        );
bool IsManifold(                                                          // Check if a given curve is manifold
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) segment-segment adjacency
        );
bool IsClockWise(                                                         // Check if a given curve is marched in clockwise direction
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) segment-segment adjacency
        );
bool IsCounterClockWise(                                                  // Check if a given curve is marched counter-clockwise direction
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input) segment-segment adjacency
        );
void AdjustOrder(                                                         // Adjust marching order to match seed direction
        std::vector<std::vector<int>>                       &,                                    // (input/output) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input/output) segment-segment adjacency
        int                                                                   // (input) seed global index
        );
void InvertOrder(                                                         // Invert marching order
        std::vector<std::vector<int>>                       &,                                    // (input/output) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input/output) segment-segment adjacency
        );
void FindCorners(                                                         // Find corners on a curve
        double                           ,                                    // (input) tolerance for corner detection
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) segment-segment adjacency
        std::vector<std::array<double,3>>                      *t,                                    // (input/optional) tangent unit vector to each curve segment
        std::vector<std::vector<int>>                       &                                     // (input/output) corner list
        );
void FindTriplePoints(                                                    // Find triple points in a PL curve
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) segment-segment adjacency
        std::vector<std::vector<int>>                       &                                     // (input/output) triple points list
        );
void FindEndPoints(                                                       // Find endpoints in a PL curve
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) segment-segment adjacency
        std::vector<std::vector<int>>                       &                                     // (input/output) endpoints list
        );
void BreakCurveAtPoint(                                                   // Break a PLG curve at a specified vertex
        int                              ,                                    // (input) global index of segment which vertex belongs to
        int                              ,                                    // (input) local index of vertex onto segments
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &                                     // (input/output) segment-segment adjacency
        );
double Length(                                                            // Compute curve lenght
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &                                     // (input) simplex-vertex connectivity
        );
void Curvature(                                                           // Compute curvature (simple curves only)
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connect
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) segment-segment adjacency
        std::vector<double>                       &                                     // (input/output) curvature at curve vertex
        );
};

// Algorithms =============================================================== //
namespace CGAlgorithms {
double Grad1DUpdate(                                                      // Update the local solution to the 1D grad limiting equation on a vertex of a 1D manifold
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input) scalar field to be limited
        int                              ,                                    // (input) global index of simplex containing the vertex
        int                              ,                                    // (input) local index of vertex
        double                           ,                                    // (input) max slope
        std::vector<bool>                       &                                     // (input) flag for dead/alive vertices
        );
void GradLimiting1D(                                                      // Solve the grad limiting eq. on a 1D manifold in a 2D Euclidean space
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input) scalar field to be limited
        double                                                                // (input) max slope
        );
double Grad2DUpdate(                                                      // Update the local solution to the 2D grad limiting equation on a vertex of a 2D manifold
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) list of simplicies in the 1-ring of the given vertex
        std::vector<double>                       &,                                    // (input) scalar field to be limited
        int                              ,                                    // (input) global index of simplex containing the vertex
        int                              ,                                    // (input) local index of vertex
        double                           ,                                    // (input) max slope
        std::vector<bool>                       &                                     // (input) flag for dead/alive vertices
        );
void GradLimiting2D(                                                      // Solve the grad limiting eq. on a 2D manifold in a 3D Euclidean space
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<double>                       &,                                    // (input/output) scalar field to be limited
        double                                                                // (input) max slope
        );
double Grad2DUpdate(                                                      // Update the local solution to the 2D grad limiting equation on a cell of a 2D volume
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input) scalar field to be limited
        int                              ,                                    // (input) global index of simplex to be updated
        double                           ,                                    // (input) max slope
        std::vector<bool>                       &                                     // (input) flag for dead/alive vertices
        );
void GradLimiting2D(                                                      // Solve the grad limiting eq. in a 2D volume
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity,
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input/output) scalar field to be limited
        double                                                                // (input) max slope
        );
};


# endif
