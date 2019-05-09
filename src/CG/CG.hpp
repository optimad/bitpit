/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

// Standard Template Library
# include <array>
# include <vector>

// bitpit library
# include "bitpit_common.hpp"

namespace bitpit{


namespace CGElem{

typedef std::array<double,3> array3D ;

/*!
 *  \ingroup CGElem
 *  \{
 */

/*!
 * The edge vertex connectivty of a box
 */
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

/*!
 * The face vertex connectivty of a box
 */
static const std::array< std::array<int,4>, 6> boxFaceVertexConnectivity =
{{
    std::array<int,4>{ {0,2,6,4} },
    std::array<int,4>{ {1,3,7,5} },
    std::array<int,4>{ {0,4,5,1} },
    std::array<int,4>{ {2,6,7,3} },
    std::array<int,4>{ {0,1,3,2} },
    std::array<int,4>{ {4,5,7,6} }
}};

/*!
 * \}
 */

bool validSegment( array3D const &, array3D const & );
bool validLine( array3D const &, array3D const & );
bool validPlane( array3D const &, array3D const & );
bool validTriangle( array3D const &, array3D const &, array3D const & );
bool validBarycentric( double const * &, int );

int convertBarycentricToFlagTriangle( std::array<double,3> const &);
int convertBarycentricToFlagSegment( std::array<double,2> const &);
int convertBarycentricToFlagPolygon( std::vector<double> const &);
int convertBarycentricToFlagPolygon( std::size_t, double const *);

void computeGeneralizedBarycentric( array3D const &, std::vector<array3D> const &, std::vector<double> &);
void computeGeneralizedBarycentric( array3D const &, std::size_t, array3D const *, std::vector<double> &);

array3D reconstructPointFromBarycentricSegment( array3D const &, array3D const &, std::array<double,2> const & ); 
array3D reconstructPointFromBarycentricSegment( array3D const &, array3D const &, double const * ); 
array3D reconstructPointFromBarycentricTriangle( array3D const &, array3D const &, array3D const &, std::array<double,3> const & ); 
array3D reconstructPointFromBarycentricTriangle( array3D const &, array3D const &, array3D const &, double const * ); 
array3D reconstructPointFromBarycentricPolygon( std::vector<array3D> const &, std::vector<double> const & ); 
array3D reconstructPointFromBarycentricPolygon( std::size_t, array3D const *, std::vector<double> const & );

array3D projectPointLine( array3D const &, array3D const &, array3D const & );
array3D projectPointPlane( array3D const &, array3D const &, array3D const & );
array3D projectPointSegment( array3D const &, array3D const &, array3D const & );
array3D projectPointSegment( array3D const &, array3D const &, array3D const &, std::array<double,2> & );
array3D projectPointSegment( array3D const &, array3D const &, array3D const &, double* );
array3D projectPointTriangle( array3D const &, array3D const &, array3D const &, array3D const & );
array3D projectPointTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D & );
array3D projectPointPolygon( array3D const &, std::vector<array3D> const & );
array3D projectPointPolygon( array3D const &, std::size_t, array3D const * );
array3D projectPointPolygon( array3D const &, std::vector<array3D> const &, std::vector<double> & );
array3D projectPointPolygon( array3D const &, std::size_t, array3D const *, std::vector<double> & );
array3D projectPointCone( array3D const &, array3D const &, array3D const &, double const &);
std::vector<array3D> projectCloudTriangle( std::vector<array3D> const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> &);

array3D restrictPointTriangle( array3D const &, array3D const &, array3D const &, array3D &);
array3D restrictPointTriangle( array3D const &, array3D const &, array3D const &, double *);

double distancePointLine( array3D const &, array3D const &, array3D const &, array3D & ) ;
double distancePointPlane( array3D const &, array3D const &, array3D const &, array3D & ) ;

double distancePointSegment( array3D const &, array3D const &, array3D const & );
double distancePointSegment( array3D const &, array3D const &, array3D const &, std::array<double,2> & );

double distancePointTriangle( array3D const &, array3D const &, array3D const &, array3D const &);
double distancePointTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D &);

double distancePointPolygon( array3D const &, std::vector<array3D> const &, array3D &, int & );
double distancePointPolygon( array3D const &, std::size_t, array3D const *, array3D &, int & );
double distancePointPolygon( array3D const &, std::vector<array3D> const & );
double distancePointPolygon( array3D const &, std::size_t, array3D const * );
double distancePointPolygon( array3D const &, std::vector<array3D> const &, std::vector<double> & );
double distancePointPolygon( array3D const &, std::size_t, array3D const *, std::vector<double> & );

double distancePointCone( array3D const &, array3D const &, array3D const &, double const &);

std::vector<double> distanceCloudTriangle( std::vector<array3D> const &, array3D const &, array3D const &, array3D const &);
std::vector<double> distanceCloudTriangle( std::vector<array3D> const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> & );

std::vector<double> distanceCloudPolygon( std::vector<array3D> const &, std::vector<array3D> const &, std::vector<array3D> &, std::vector<int> & );
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &, std::size_t, array3D const *, std::vector<array3D> &, std::vector<int> & );
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &, std::vector<array3D> const &);
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &, std::size_t, array3D const *);
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &, std::vector<array3D> const &, std::vector<std::vector<double>> &);
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &, std::size_t, array3D const *, std::vector<std::vector<double>> &);

double distanceLineLine(array3D const &, array3D const &, array3D const &, array3D const &);
double distanceLineLine(array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D &);


bool intersectPointSegment( array3D const &, array3D const &, array3D const & ) ;
bool intersectPointTriangle( array3D const &, array3D const &, array3D const &, array3D const & ) ;
bool intersectPointBox( array3D const &, array3D const &, array3D const &, int dim=3 ) ;

bool intersectLineLine( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool intersectLinePlane( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool intersectLineTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool intersectLinePolygon( array3D const &, array3D const &, std::vector<array3D> const &, array3D & ) ;
bool intersectLinePolygon( array3D const &, array3D const &, std::size_t, array3D const *, array3D & ) ;

bool intersectSegmentSegment( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool intersectSegmentPlane( array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool intersectSegmentTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, array3D & ) ;
bool intersectSegmentPolygon( array3D const &, array3D const &, std::vector<array3D> const &, array3D & ) ;
bool intersectSegmentPolygon( array3D const &, array3D const &, std::size_t, array3D const *, array3D & ) ;

bool intersectSegmentBox( array3D const &, array3D const &, array3D const &, array3D const &, int  dim = 3 ) ;
bool intersectSegmentBox( array3D const &, array3D const &, array3D const &, array3D const &, bool, bool, std::vector<array3D> &, int dim=3 );
bool intersectSegmentBox( array3D const &, array3D const &, array3D const &, array3D const &, bool, bool, std::vector<array3D> &, std::vector<int> &, int dim=3 );

bool intersectPlanePlane( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & ) ;
bool intersectPlaneBox( array3D const &, array3D const &, array3D const &, array3D const &, int dim=3);
bool intersectPlaneBox( array3D const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> &, int dim=3);

bool intersectBoxBox( array3D const &, array3D const &, array3D const &, array3D const &, int  dim = 3 ) ;
bool intersectBoxBox( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D &, int  dim = 3 ) ;


bool intersectBoxTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, int dim=3 ) ;
bool intersectBoxTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, bool, bool, bool, std::vector<array3D> &, int dim=3 ) ;
bool intersectBoxTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, bool, bool, bool, std::vector<array3D> &, std::vector<int> &, int dim=3 ) ;


bool intersectBoxPolygon( array3D const &, array3D const &, std::vector<array3D> const &, int dim=3 );
bool intersectBoxPolygon( array3D const &, array3D const &, std::size_t, array3D const *, int dim=3 );
bool intersectBoxPolygon( array3D const &, array3D const &, std::vector<array3D> const &, bool, bool, bool, std::vector<array3D> &, int dim=3);
bool intersectBoxPolygon( array3D const &, array3D const &, std::size_t, array3D const *, bool, bool, bool, std::vector<array3D> &, int dim=3);
bool intersectBoxPolygon( array3D const &, array3D const &, std::vector<array3D> const &, bool, bool, bool, std::vector<array3D> &, std::vector<int> &, int dim=3);
bool intersectBoxPolygon( array3D const &, array3D const &, std::size_t, array3D const *, bool, bool, bool, std::vector<array3D> &, std::vector<int> &, int dim=3);


void computeAABBSegment( array3D const &, array3D const &, array3D &, array3D & ) ;
void computeAABBTriangle( array3D const &, array3D const &, array3D const &, array3D &, array3D & ) ;
void computeAABBPolygon( std::vector<array3D> const &, array3D &, array3D & );
void computeAABBPolygon( std::size_t, array3D const *, array3D &, array3D & );

void unionAABB( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & );
void intersectionAABB( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & );
void subtractionAABB( array3D const &, array3D const &, array3D const &, array3D const &, array3D &, array3D & );
void unionAABB( std::vector<array3D>  const &, std::vector<array3D> const &, array3D &, array3D & );

array3D rotateVector( array3D const &, array3D const &, double);
double areaTriangle( array3D const &, array3D const &, array3D const &);

void vertexOfSegment( int const &, array3D const &, array3D const &, array3D &);

void vertexOfTriangle( int const &, array3D const &, array3D const &, array3D const &, array3D &);
void edgeOfTriangle( int const &, array3D const &, array3D const &, array3D const &, array3D &, array3D &);

void vertexOfBox( int const &, array3D const &, array3D const &, array3D &);
void edgeOfBox( int const &, array3D const &, array3D const &, array3D &, array3D &);
void faceOfBox( int const &, array3D const &, array3D const &, array3D &, array3D &, array3D &, array3D & );

int polygonEdgesCount(std::vector<array3D> const &);
int polygonEdgesCount(std::size_t, array3D const *);
int polygonSubtriangleCount(std::vector<array3D> const &);
int polygonSubtriangleCount(std::size_t, array3D const *);
void edgeOfPolygon( int const &, std::vector<array3D> const &, array3D &, array3D &);
void edgeOfPolygon( int const &, std::size_t, array3D const *, array3D &, array3D &);
void subtriangleOfPolygon( int const &, std::vector<array3D> const &, array3D &, array3D &, array3D &);
void subtriangleOfPolygon( int const &, std::size_t, array3D const *, array3D &, array3D &, array3D &);

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
