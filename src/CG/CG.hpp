/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

# ifndef __BITPIT_CG_BASE_HPP__
# define __BITPIT_CG_BASE_HPP__

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
# include "bitpit_operators.hpp"

namespace bitpit{

// Basic computational geometry ============================================= //
namespace CGElem{

typedef std::array<double,3> array3D ;

static const std::array< std::array<int,2>,12> boxEdgeVertexConnectivity =
{{
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

}};

static const std::array< std::array<int,4>, 6> boxFaceVertexConnectivity =
{{
    std::array<int,4>{ {0,2,6,4} },
    std::array<int,4>{ {1,3,7,5} },
    std::array<int,4>{ {0,4,5,1} },
    std::array<int,4>{ {2,6,7,3} },
    std::array<int,4>{ {0,1,3,2} },
    std::array<int,4>{ {4,5,7,6} }
}};


double              distancePointLine( array3D const &, array3D const &, array3D const &, array3D & ) ;
double              distancePointPlane( array3D const &, array3D const &, array3D const &, array3D & ) ;
double              distancePointSegment( array3D const &, array3D const &, array3D const &, array3D &, int & ) ;
double              distancePointSegment( array3D const &, array3D const &, array3D const &, array3D &, std::array<double,2> &, int & ) ;
double              distancePointTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, int & ) ;
double              distancePointTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D &, int & ) ;
double              distancePointSimplex( array3D const &, std::vector<array3D> const &, array3D &, int & );


std::vector<double> distanceCloudTriangle( std::vector<array3D> const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> &, std::vector<int> & );
std::vector<double> distanceCloudSimplex( std::vector<array3D> const &, std::vector<array3D> const &, std::vector<array3D> &, std::vector<int> & );

std::vector<double> distanceCloudTriangle( std::vector<array3D> const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> * const, std::vector<array3D> * const );
std::vector<double> distanceCloudSimplex( std::vector<array3D> const &, std::vector<array3D> const &);

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
void                intersectionAABB( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & );
void                subtractionAABB( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & );
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


void rotateVector(
        std::array<double, 3>                    &,
        std::array<double, 3> const              &,
        double                                     
        );


void vertexOfBox(
        int              const              &,     
        std::array<double, 3> const              &,
        std::array<double, 3> const              &,
        std::array<double, 3>                    & 
        );

void edgeOfBox(
        int              const              &,     
        std::array<double, 3> const              &,
        std::array<double, 3> const              &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    & 
        );

void faceOfBox(
        int              const              &,     
        std::array<double, 3> const              &,
        std::array<double, 3> const              &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    &,
        std::array<double, 3>                    & 
        );





}

// Algorithms =============================================================== //
namespace CGAlgorithms {

double grad1DUpdate(                                                      // Update the local solution to the 1D grad limiting equation on a vertex of a 1D manifold
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

void gradLimiting1D(                                                      // Solve the grad limiting eq. on a 1D manifold in a 2D Euclidean space
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<std::vector<int>>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input) scalar field to be limited
        double                                                                // (input) max slope
        );

double grad2DUpdate(                                                      // Update the local solution to the 2D grad limiting equation on a vertex of a 2D manifold
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

void gradLimiting2D(                                                      // Solve the grad limiting eq. on a 2D manifold in a 3D Euclidean space
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<double>                       &,                                    // (input/output) scalar field to be limited
        double                                                                // (input) max slope
        );

double grad2DUpdate(                                                      // Update the local solution to the 2D grad limiting equation on a cell of a 2D volume
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input) scalar field to be limited
        int                              ,                                    // (input) global index of simplex to be updated
        double                           ,                                    // (input) max slope
        std::vector<bool>                       &                                     // (input) flag for dead/alive vertices
        );

void gradLimiting2D(                                                      // Solve the grad limiting eq. in a 2D volume
        int                              ,                                    // (input) number of simplicies
        std::vector<std::array<double,3>>                       &,                                    // (input) vertex coordinate list
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-vertex connectivity,
        std::vector<std::vector<int>>                       &,                                    // (input) simplex-simplex adjacency
        std::vector<double>                       &,                                    // (input/output) scalar field to be limited
        double                                                                // (input) max slope
        );
}

}


# endif
