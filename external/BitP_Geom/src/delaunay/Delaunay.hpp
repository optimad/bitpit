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
# ifndef __DY_HPP__
# define __DY_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <queue>
# include <array>
# include <vector>
# include <string>
# include <iostream>

// bitpit library
//# include <SurfTriPatch.hpp>
//# include <VolTri.hpp>
# include <bitpit.hpp>

// CC_Lib
# include "CGBase.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

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

// SurfTriPatch vector
typedef vector< SurfTriPatch >        Svector1D;
typedef vector< Svector1D >            Svector2D;
typedef vector< Svector2D >            Svector3D;
typedef vector< Svector3D >            Svector4D;

// ========================================================================== //
// FUNCTION PROTOTYPES                                                        //
// ========================================================================== //

// Delaunay triangulation =================================================== //
namespace Delaunay {
    namespace Delaunay2D {
        bool SWAP_Criterion(                                                  // Cline-Rinka criterion for edge swapping
            a3vector1D                  &,                                    // (input) coordinates of 1st vertex on tested edge
            a3vector1D                  &,                                    // (input) coordinates of 2nd vertex on tested edge
            a3vector1D                  &,                                    // (input) coordinates of the 1st vertex opposite to tested edge
            a3vector1D                  &                                     // (input) coordinates of the 2nd vertex opposite to tested edge
        );
        bool ConvexPolygon(                                                   // Check if two triangles form a closed polygon
            VolTriPatch                &,                                    // (input) volume triangulation
            int                   const &,                                    // (input) 1st triangle
            int                   const &                                     // (input) 2nd triangle
        );
        void InsertVertices(                                                  // Insert vertices into Delaunay triangulation
            VolTriPatch                &,                                    // (input/output) unstructured volume mesh manager
            ivector1D                   &                                     // (input) list of vertices to be inserted
        );
        void InsertEdges(                                                     // Insert edges into Delaunay triangulation
            VolTriPatch                &,                                    // (input/output) unstructured volume triangulation
            ivector2D                   &,                                    // (input) edge list
            ivector2D                   &                                     // (input/output) vertex-simplex connectivity
        );
        void ReturnSimplexID(                                                 // Returns ID of the simplex enclosing a given point
            VolTriPatch          const &,                                    // (input) unstructured volume tasselation
            a3vector1D            const &,                                    // (input) point coordinates
            int                         &                                     // (input/output) seed for search algorithm
        );
        ivector1D FindPath(                                                   // Returns the path of simplicies crossed by a given edge
            VolTriPatch                &,                                    // (input) volume triangulation
            ivector1D             const &                                     // (input) edge
        );
        unsigned int Delaunay2D(                                              // Compute the 2D Delaunay triangulation of a set of vertices
            VolTriPatch                &                                     // (input/output) unstructured volume triangulation
        );
        unsigned int CDelaunay2D(                                             // Constrained Delaunay triangulation
            VolTriPatch                &,                                    // (input/output) unstructured volume triangulation
            ivector2D                   *,                                    // (input) pointer to edge-vertex connect for domain boundaries
            vector<ivector2D*>          &,                                    // (input) pointers to edge-vertex connect for each domain hole
            vector<ivector2D*>          &                                     // (input) pointers to edge-vertex connect for each domain segment
        );
        unsigned int MAT_distance(                                            // Compute distance from approximated medial axis
            a3vector2D                  &,                                    // (input/output) vertex coordinate list
            SurfTriPatch                &,                                    // (input/output) surface mesh of domain boundaries
            Svector1D                   &,                                    // (input/output) surface mesh of domain holes
            Svector1D                   &,                                    // (input/output) surface mesh of domain segments
            double                                                            // (input) refinement parameter
        );
        void CCDT(                                                            // Constrained Capacity Delaunay algorithm
            VolTriPatch                 &,                                    // (input/output) unstructured volume mesh manager
            bvector1D                   &,                                    // (input) flag for movable vertices
            ivector2D                   &                                     // (input) 1-ring of vertices
        );
        void Smooth(                                                          // Optimize vertex locations using a CCDT algorithm for a Delaunay triangulation
            VolTriPatch                &,                                    // (input/output) unstructured volume mesh manager
            bvector1D                   &                                     // (input) flag for movable vertex
        );
        voilass_VolTri
            VolTriPatch                &,                                    // (input/output) unstructured volume mesh manager
            ivector2D                   &,                                    // (input) edge list
            bvector1D                   &                                     // (input) flag for movable vertex
        );
    };
};

// Quad mesh generation ===================================================== //
namespace QDelaunay {
    namespace QDelaunay2D {
        void PlaceVertices(                                                   // Insert vertices into a closed bounded domain for quad mesh generation
            VolTriPatch                                &,                    // (input/output) unstructured volume mesh
            CGPolygon2D::Class_CG_Polygon2D            &,                    // (input) polygonal shape describing domain boundaries
            vector<CGPolygon2D::Class_CG_Polygon2D>    &,                    // (input) polygonal shapes describing each domain hole
            ivector2D                                   *,                    // (input) pointer to segment-vertex connect. for domain boundaries
            vector<ivector2D*>                          &,                    // (input) pointers to segment-vertex connect. for domain holes
            vector<ivector2D*>                          &                     // (input) pointers to segment-vertex connect. for domain segments
        );
    };
};

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //

# endif
