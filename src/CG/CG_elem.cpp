/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

# include <cmath>
# include <string>
# include <iostream>

# include <assert.h>

# include "bitpit_private_lapacke.hpp"

# include "bitpit_operators.hpp"

# include "CG.hpp"
# include "CG_private.hpp"


namespace bitpit{

namespace CGElem{


/*!
 * \ingroup CGElem
 * \{
*/

/*!
 * \private
 * Projects points on a triangle.
 * Projection points are the closest points to the original points within the triangle.
 *
 * \param[in] nPoints number of points
 * \param[in] points pointer to points' coordinates
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \param[out] proj pointer to the projection point; 
 * \param[out] lambda pointer to barycentric coordinates of projection points
 * \return distances
 */
void _projectPointsTriangle( int nPoints, array3D const *point, array3D const &Q0, array3D const &Q1, array3D const &Q2, array3D *proj, double *lambda )
{

    assert( validTriangle(Q0,Q1,Q2) );

    array3D s0 = Q1-Q0;
    array3D s1 = Q2-Q0;

    double A[4] = { dotProduct(s0,s0), 0, dotProduct(s0,s1), dotProduct(s1,s1) }  ; 
    double *B = new double [2*nPoints];

    for( int i=0; i<nPoints; ++i){
        array3D rP = *point -Q0;
        B[2*i]   = dotProduct(s0,rP); 
        B[2*i+1] = dotProduct(s1,rP); 
        ++point;
    }

    int info =  LAPACKE_dposv( LAPACK_COL_MAJOR, 'U', 2, nPoints, A, 2, B, 2 );
    assert( info == 0 );
    BITPIT_UNUSED( info );

    for( int i=0; i<nPoints; ++i){

        double *b = &B[2*i];

        lambda[0] = 1. -b[0] -b[1];
        lambda[1] = b[0];
        lambda[2] = b[1];

        *proj = restrictPointTriangle( Q0, Q1, Q2, lambda);
        lambda +=3;
        proj += 1;
    }

    delete [] B;
}

/*!
 * \private
 * Project points on a plane described by a triangle
 *
 * \param[in] nPoints number of points
 * \param[in] points pointer to points' coordinates
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \param[out] proj pointer to the projection point; 
 * \param[out] lambda pointer to barycentric coordinates of projection points
 * \return distances
 */
void _projectPointsPlane( int nPoints, array3D const *point, array3D const &Q0, array3D const &Q1, array3D const &Q2, array3D *proj, double *lambda )
{

    assert( validTriangle(Q0,Q1,Q2) );

    array3D s0 = Q1-Q0;
    array3D s1 = Q2-Q0;

    double A[4] = { dotProduct(s0,s0), 0, dotProduct(s0,s1), dotProduct(s1,s1) }  ; 
    double *B = new double [2*nPoints];

    for( int i=0; i<nPoints; ++i){
        array3D rP = *point -Q0;
        B[2*i]   = dotProduct(s0,rP); 
        B[2*i+1] = dotProduct(s1,rP); 
        ++point;
    }

    int info =  LAPACKE_dposv( LAPACK_COL_MAJOR, 'U', 2, nPoints, A, 2, B, 2 );
    assert( info == 0 );
    BITPIT_UNUSED( info );

    for( int i=0; i<nPoints; ++i){

        double *b = &B[2*i];

        lambda[0] = 1. -b[0] -b[1];
        lambda[1] = b[0];
        lambda[2] = b[1];

        *proj = reconstructPointFromBarycentricTriangle( Q0, Q1, Q2, lambda);
        lambda +=3;
        proj += 1;
    }

    delete [] B;
}

/*!
 * \private
 * Computes intersection between an axis aligned bounding box and a triangle
 * \param[in] A1 min point of first box
 * \param[in] A2 max point of first box
 * \param[in] V1 first vertex of triangle
 * \param[in] V2 second vertex of triangle
 * \param[in] V3 third vertex of triangle
 * \param[out] intrPtr pointed vector will hold intersection points between triangle edges and box faces
 * \param[out] flagPtr pointed vector will have the same size of P. If the ith flag=0, the intersection is due to interiorTriangleVertice. If the ith flag=1, the intersection is due to triangleEdgeBoxHullIntersection. If the ith flag=2, the intersection is due to triangleBoxEdgeIntersections.
 * \return if intersect
 */
bool _intersectBoxTriangle(array3D const &A0, array3D const &A1, array3D const &V0, array3D const &V1, array3D const &V2, bool interiorTriangleVertice, bool triangleEdgeBoxHullIntersections, bool triangleBoxEdgeIntersection, std::vector<array3D> *intrPtr, std::vector<int> *flagPtr, int dim)
{

    bool intersect(false);
    bool addFlag( flagPtr!=nullptr);
    bool computeIntersection(interiorTriangleVertice||triangleBoxEdgeIntersection||triangleEdgeBoxHullIntersections);

    assert( ! (computeIntersection && (intrPtr==nullptr) ) );

    if(computeIntersection){
        intrPtr->clear();
    }

    if(addFlag){
        flagPtr->clear();
    }

    //check if Triangle Boundig Box and Box overlap -> necessary condition
    array3D B0, B1;
    computeAABBTriangle( V0, V1, V2, B0, B1);

    if( !intersectBoxBox( A0, A1, B0, B1, dim) ) { 
        return false; 
    }

    //check if triangle vertices are within the box
    for( int i=0; i<3; ++i){
        vertexOfTriangle(i, V0, V1, V2, B0);
        if( intersectPointBox(B0, A0, A1, dim) ){
            intersect = true;
            if(!interiorTriangleVertice) break;

            intrPtr->push_back(B0);
            if(addFlag) flagPtr->push_back(0);
        }
    }

    //check if triangle edges and box faces intersect
    if( !intersect || triangleEdgeBoxHullIntersections) {

        for( int edge=0; edge<3; ++edge){

            edgeOfTriangle(edge, V0, V1, V2, B0, B1);

            if(dim==2){
                array3D p;
                array3D faceVertex0, faceVertex1;

                for( int face=0; face<4; ++face){
                    edgeOfBox( face, A0, A1, faceVertex0, faceVertex1);

                    if( intersectSegmentSegment(B0, B1, faceVertex0, faceVertex1, p) ){
                        intersect=true;
                        if(!triangleEdgeBoxHullIntersections) break;

                        intrPtr->push_back(p);
                        if(addFlag) flagPtr->push_back(1);
                    }

                }

            } else if(dim==3){
                array3D p;
                std::vector<array3D> V(4);

                for( int face=0; face<6; ++face){
                    faceOfBox( face, A0, A1, V[0], V[1], V[2], V[3] );

                    if( intersectSegmentPolygon( B0, B1, V, p ) ) {
                        intersect=true;
                        if(!triangleEdgeBoxHullIntersections) break;

                        intrPtr->push_back(p);
                        if(addFlag) flagPtr->push_back(1);
                    }
                }
            }
        }
    }

    //check if triangle and box edges (dim=3) or box vertices (dim=2) intersect
    if( !intersect || triangleBoxEdgeIntersection ) {

        array3D p;

        if(dim==2){
            for( int i=0; i<4; ++i){
                vertexOfBox( i, A0, A1, B0);
                if( intersectPointTriangle(B0,V0,V1,V2)) {
                    intersect = true;
                    if(!triangleBoxEdgeIntersection) break;

                    intrPtr->push_back(B0);
                    if(addFlag) flagPtr->push_back(2);
                }
            }

        } else if(dim==3){
            for( int i=0; i<12; ++i){
                edgeOfBox( i, A0, A1, B0, B1);
                if( intersectSegmentTriangle(B0,B1,V0,V1,V2,p)) {
                    intersect = true;
                    if(!triangleBoxEdgeIntersection) break;

                    intrPtr->push_back(p);
                    if(addFlag) flagPtr->push_back(2);
                }
            }
        }
    }

    return intersect;
}

/*!
 * \private
 * Computes intersection between an axis aligned bounding box and a segment
 * \param[in] V0 start point of segment
 * \param[in] V1 end point of segment
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] interiorSegmentVertice if the segment vertices within the box should be added to the list of intersection points
 * \param[in] segmentBoxHullIntersection if the intersections between the segment and the outer hull of the box should be added to the list of intersection points
 * \param[in,out] intrPtr pointer to the list of intersection points. If no intersetion points should be calculated nullptr can be passed as argument
 * \param[in,out] flagPtr if (!=nullptr), for each intersection point a flag will be inserted indicatingif it belongs to interiorSegmentVertice (flag=0) or segmentBoxHullIntersection (flag=1)
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool _intersectSegmentBox(array3D const &V0, array3D const &V1, array3D const &A0, array3D const &A1, bool interiorSegmentVertice, bool segmentBoxHullIntersection, std::vector<array3D> *intrPtr, std::vector<int> *flagPtr, int dim)
{

    bool intersect(false);

    bool addFlag(flagPtr!=nullptr);
    bool computeIntersection(interiorSegmentVertice||segmentBoxHullIntersection);

    assert( ! (computeIntersection && (intrPtr==nullptr) ) );

    if(computeIntersection){
        intrPtr->clear();
    }

    if(addFlag){
        flagPtr->clear();
    }

    array3D p, B0, B1;

    //check if segment boundig box and Box overlap
    computeAABBSegment( V0, V1, B0, B1);
    if( !intersectBoxBox( A0, A1, B0, B1, dim) ) { 
        return false;
    }

    //check segment points
    for( int i=0; i<2; ++i){
        vertexOfSegment(i, V0, V1, B0);

        if( intersectPointBox(B0,A0,A1) ){
            intersect = true;
            if(!interiorSegmentVertice) break;

            intrPtr->push_back(V0);
            if(addFlag) flagPtr->push_back(0);
        }
    }

    //check if segment intersects outer hull of box
    if( !intersect || segmentBoxHullIntersection){
        if( dim == 2){ //check if box edge and segment intersect

            for( int i=0; i<4; ++i){
                edgeOfBox( i, A0, A1, B0, B1);

                if( intersectSegmentSegment(B0,B1,V0,V1,p)) {
                    intersect = true;
                    if(!segmentBoxHullIntersection) break;

                    intrPtr->push_back(p);
                    if(addFlag) flagPtr->push_back(1);
                }
            }


        } else if( dim==3 ) { //3D check if box face and segment intersect

            std::vector< array3D > E(4);

            for( int i=0; i<6; ++i){
                faceOfBox( i, A0, A1, E[0], E[1], E[2], E[3]);

                if( intersectSegmentPolygon(V0,V1,E,p) ) {
                    intersect = true;
                    if(!segmentBoxHullIntersection) break;

                    intrPtr->push_back(p);
                    if(addFlag) flagPtr->push_back(1);
                }
            }
        }
    }

    return intersect;

}

/*!
 * \private
 * Computes intersection between a plane and an axis aligned bounding box
 * \param[in] P point on plane
 * \param[in] N plane normal
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in,out] intrPtr if (intrPtr!=nullptr) the list of intersection points between the edges of the box and the plane.
 * if (dim==3) the intersections will be ordered in counter-clockwise sense with respect to the plane normal.
 * \param[in] dim number of dimensions
 * \return if intersect
 */
bool _intersectPlaneBox(array3D const &P, array3D const &N, array3D const &A0, array3D const &A1, std::vector<array3D> *intrPtr, int dim)
{

    bool intersect(false);
    bool computeIntersection(intrPtr);

    if(computeIntersection){
        intrPtr->clear();
    } else if( intersectPointBox(P,A0,A1,dim)) {
        return true;
    }

    // check if plane intersects the edges of the box
    // and store eventually intersection points
    int edgeCount = (dim==2) ? 4 : 12;
    array3D E0, E1, V;

    for(int i=0; i<edgeCount; ++i){
        edgeOfBox(i, A0, A1, E0, E1);

        if( intersectSegmentPlane( E0, E1, P, N, V) ){
            intersect = true;

            if(intrPtr){
                intrPtr->push_back(V);
            } else {
                break;
            }
        }
    }


    // sort intersection points in couterclock-wise sense
    if( dim==3 && intrPtr && intersect){

        const array3D origin = intrPtr->at(0);

        std::sort(intrPtr->begin(), intrPtr->end(), [&](const array3D &lhs, const array3D &rhs) -> bool {
            array3D v = crossProduct(lhs-origin,rhs-origin);
            return dotProduct(v, N) < 0;
        } );
    }

    return intersect;
}

/*!
 * \private
 * Computes intersection between an axis aligned bounding box and a convex polygon
 * \param[in] A0 min point of first box
 * \param[in] A1 max point of first box
 * \param[in] VS polygon vertices coordinates
 * \param[in] innerPolygonPoints simplex vertices within the box should be added to the intersection list
 * \param[in] polygonEdgeBoxFaceIntersection intersection between the edges of the polygon and the hull of the box should be added to the intersection list
 * \param[in] polygonBoxEdgeIntersection intersection between the polygon and the edges of the box should be added to the intersection list
 * \param[out] P intersection points 
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool _intersectBoxPolygon(array3D const &A0, array3D const &A1, std::vector<array3D> const &VS, bool innerPolygonPoints, bool polygonEdgeBoxHullIntersection, bool polygonBoxEdgeIntersection, std::vector<array3D> *intrPtr, std::vector<int> *flagPtr, int dim)
{

    bool intersect(false);
    bool addFlag(flagPtr!=nullptr);
    bool computeIntersection(innerPolygonPoints || polygonEdgeBoxHullIntersection || polygonBoxEdgeIntersection);

    assert( ! (computeIntersection && (intrPtr==nullptr) ) );

    if(computeIntersection){
        intrPtr->clear();
    }

    if(addFlag){
        flagPtr->clear();
    }

    array3D V0, V1, V2;

    //check if simplex boundig box and box overlap -> necessary condition
    computeAABBPolygon( VS, V0, V1);
    if( !intersectBoxBox( A0, A1, V0, V1, dim) ) { 
        return false; 
    }
    
    std::vector<array3D> partialIntr;
    std::vector<int> partialFlag;

    // check if triangle vertices lie within box
    // or if triangles intersect edges of box
    computeIntersection = innerPolygonPoints || polygonBoxEdgeIntersection;
    int trianglesCount = polygonSubtriangleCount(VS);
    for (int triangle=0; triangle<trianglesCount; ++triangle) {
        subtriangleOfPolygon(triangle, VS, V0, V1, V2);

        if( _intersectBoxTriangle( A0, A1, V0, V1, V2, innerPolygonPoints, false, polygonBoxEdgeIntersection, &partialIntr, &partialFlag, dim ) ){

            intersect = true;
            if(!computeIntersection) break;

            int intrCount = partialIntr.size();
            for( int i=0; i<intrCount; ++i){

                array3D &candidateCoord = partialIntr[i];
                int candidateFlag = partialFlag[i];

                //prune duplicate points
                auto PItr = intrPtr->begin();
                bool iterate = (PItr!=intrPtr->end());
                while(iterate){

                    iterate = !utils::DoubleFloatingEqual()( norm2( *PItr -candidateCoord ), 0. );
                
                    if(iterate){
                        ++PItr;
                    }
                    iterate &= PItr!=intrPtr->end();
                }

                if(PItr!=intrPtr->end()){
                    continue;
                }

                intrPtr->push_back(candidateCoord);
                if(addFlag){
                    flagPtr->push_back(candidateFlag);
                }
            }
        }
    }

    // check if edges of polygon intersect box face
    computeIntersection = polygonEdgeBoxHullIntersection;
    int edgesCount = polygonEdgesCount(VS);
    if(!intersect || polygonEdgeBoxHullIntersection){
        for (int edge=0; edge<edgesCount; ++edge) {
            edgeOfPolygon(edge, VS, V0, V1);

            if( _intersectSegmentBox( V0, V1, A0, A1, false, polygonEdgeBoxHullIntersection, &partialIntr, &partialFlag, dim) ){

                intersect = true;
                if(!computeIntersection) break;

                int intrCount = partialIntr.size();
                for( int i=0; i<intrCount; ++i){

                    array3D &candidateCoord = partialIntr[i];
                    int candidateFlag = partialFlag[i];

                    //prune duplicate points
                    auto PItr = intrPtr->begin();
                    bool iterate = (PItr!=intrPtr->end());
                    while(iterate){

                        iterate = !utils::DoubleFloatingEqual()( norm2( *PItr -candidateCoord ), 0. );
                    
                        if(iterate){
                            ++PItr;
                        }
                        iterate &= PItr!=intrPtr->end();
                    }

                    if(PItr!=intrPtr->end()){
                        continue;
                    }

                    intrPtr->push_back(candidateCoord);
                    if(addFlag){
                        flagPtr->push_back(candidateFlag);
                    }
                }
            }
        }
    }

    return intersect;
}

/*!
 * Checks if a segment is valid
 * \param[in] P0 start point of segment
 * \param[in] P1 end point of segment
 * \return true if valid
 */
bool validSegment(const array3D &P0, const array3D &P1 ) 
{
    return !utils::DoubleFloatingEqual()( norm2(P1-P0), 0.);
}

/*!
 * Checks if a line is valid
 * \param[in] P point on line;
 * \param[in] n line direction
 * \return true if valid
 */
bool validLine(const array3D &P, const array3D &n ) 
{
    BITPIT_UNUSED(P);
    return utils::DoubleFloatingEqual()( norm2(n), 1.);
}

/*!
 * Checks if a plane is valid
 * \param[in] P point on plane
 * \param[in] n plane normal 
 * \return true if valid
 */
bool validPlane(const array3D &P, const array3D &n ) 
{
    BITPIT_UNUSED(P);
    return utils::DoubleFloatingEqual()( norm2(n), 1.);
}

/*!
 * Checks if a triangle is valid
 * \param[in] P0 first triangle vertex
 * \param[in] P1 second triangle vertex
 * \param[in] P2 third triangle vertex
 * \return true if valid
 */
bool validTriangle(const array3D &P0, const array3D &P1, const array3D &P2 ) 
{

    array3D v0 = P1 -P0;
    if( utils::DoubleFloatingEqual()( norm2(v0), 0.) ){
        return false;
    }

    array3D v1 = P2 -P1;
    if( utils::DoubleFloatingEqual()( norm2(v1), 0.) ){
        return false;
    }

    array3D v2 = P0 -P2;
    if( utils::DoubleFloatingEqual()( norm2(v2), 0.) ){
        return false;
    }

    array3D n = crossProduct( v0, -1.*v2);
    if( utils::DoubleFloatingEqual()( norm2(n), 0. ) ){
        return false;
    }

    return true;
}

/*!
 * Checks if barycentric coordinates is valid
 * \param[in] lambdaPtr pointer to barycentric coordinates
 * \param[in] n number of points
 * \return true if valid
 */
bool validBarycentric(double const *lambdaPtr, int n )
{

    double sum(-1.);
    double maxValue(0.);

    for(int i=0; i<n; ++i){

        double value = lambdaPtr[i];
        if( std::abs(value) > std::abs(maxValue) ){
            sum -= maxValue;
            maxValue = -value;

        } else {
            sum += value;

        }

    }

    double accuracyFactor = std::max( 1., std::abs(maxValue) );
    return utils::DoubleFloatingEqual()( sum, maxValue, accuracyFactor);
}

/*!
 * Converts barycentric coordinates of a point on a segment to a flag that indicates where the point lies.
 * Flag = 0 Point lies within the segment
 * Flag = 1 Point coincides with the first vertex or is positioned befor the line
 * Flag = 2 Point coincides with the second vertex or is positioned after the line
 * \param[in] lambda barycentric coordinates of point
 * \return flag
 */
int convertBarycentricToFlagSegment( std::array<double,2> const &lambda)
{

    assert( validBarycentric(&lambda[0],2) );

    if (lambda[0]>=1.) {
        return 1;

    } else if (lambda[1]>=1.) {
        return 2;

    } 
    
    return 0;
}

/*!
 * Converts barycentric coordinates of a point on a triangle to a flag that indicates where the point lies.
 * Flag = 0 Point lies within the triangle
 * Flag = i Point coincides with the ith vertex of triangle or lies within the area spanned by the edges incident in the ith vertex
 * Flag = -i Point lies on the edge starting from the ith vertex and connecting the following vertex in clockwise direction or in its shaddowed area
 * \param[in] lambda barycentric coordinates of point
 * \return flag
 */
int convertBarycentricToFlagTriangle( array3D const &lambda)
{

    assert( validBarycentric(&lambda[0],3) );

    int count = 0;
    std::array<int,2> zeros = {{0,0}};

    for( int i=0; i<3; ++i){
        if ( lambda[i] <= 0.) {
            zeros[count] = i;
            ++count;
        }
    }

    if( count == 1){
        int vertex0 = (zeros[0] +1) %3;
        count = -(vertex0+1);

    } else if (count == 2) {
        count = 3 -zeros[0] -zeros[1] +1;

    }

    return count;
}

/*!
 * Converts barycentric coordinates of a point on a convex polygon to a flag that indicates where the point lies.
 * Flag = 0 Point lies within the simplex
 * Flag = i Point coincides with the ith vertex of simplex or lies within the area spanned by the edges incident in the ith vertex
 * Flag = -i Point lies on the edge starting from the ith vertex and connecting the following vertex in clockwise direction or in its shaddowed area
 * \param[in] lambda barycentric coordinates of point
 * \return flag
 */
int convertBarycentricToFlagPolygon( std::vector<double> const &lambda)
{

    int N(lambda.size());
    int count(0);
    int firstPositive(N);

    assert( validBarycentric(&lambda[0],N) );

    for( int i=0; i<N; ++i){
        if ( lambda[i] > 0.) {
            firstPositive = std::min( firstPositive, i);
            ++count;
        }
    }

    if( count == 1){
        count = firstPositive +1;

    } else if (count==2) {
        count = -(firstPositive+1);

    } else {
        count = 0;
    }

    return count;
}

/*!
 * Computes Generalized Barycentric Coordinates of a point in convex polygons or polyedra.
 * No check is performed to check convexity.
 * Formula [6] of \ref{igeometry.caltech.edu/pubs/MHBD02.pdf} is implemented.
 * This formula actually refers to the method of Eugene Wachpress in the manuscript A Rational Finite Elment Basis.
 * \param[in] p point
 * \param[in] vertex vertex coordinates of polygon
 * \parm[out] lambda generalized barycentric coordinates of p
 */
void computeGeneralizedBarycentric( array3D const &p, std::vector<array3D> const &vertex, std::vector<double> &lambda)
{
    int vertexCount=vertex.size();

    lambda.resize(vertexCount);

    std::vector<double> area(vertexCount);
    for( int i=0; i<vertexCount; ++i){
        int next = (i +1) %vertexCount;
        area[i] = areaTriangle( vertex[i], vertex[next], p);
    }

    double sumWeight(0);

    for( int i=0; i<vertexCount; ++i){
        int prev = (i +vertexCount -1) %vertexCount;
        int next = (i +1) %vertexCount;
        lambda[i]  = areaTriangle(vertex[prev], vertex[i], vertex[next]);

        for( int j=0; j<vertexCount; ++j){
            if( j==prev || j==i){
                continue;
            }

            lambda[i] *= area[j];
        }

        sumWeight += lambda[i];
    }

    lambda /= sumWeight;
}

/*!
 * Reconstructs a point from barycentric coordinates of a segment
 * \param[in] Q0 first vertex of segment
 * \param[in] Q1 second vertex of segment
 * \param[in] lambda barycentric coordinates
 * \param[out] reconstructed point
 */
array3D reconstructPointFromBarycentricSegment(array3D const &Q0, array3D const &Q1, std::array<double,2> const &lambda)
{
    assert( validBarycentric(&lambda[0],2) );

    return lambda[0]*Q0 +lambda[1]*Q1;
}

/*!
 * Reconstructs a point from barycentric coordinates of a segment
 * \param[in] Q0 first vertex of segment
 * \param[in] Q1 second vertex of segment
 * \param[in] lambda barycentric coordinates
 * \param[out] reconstructed point
 */
array3D reconstructPointFromBarycentricSegment(array3D const &Q0, array3D const &Q1, double const *lambda)
{
    assert( validBarycentric(&lambda[0],2) );

    return lambda[0]*Q0 +lambda[1]*Q1;
}

/*!
 * Reconstructs a point from barycentric coordinates of a triangle
 * \param[in] Q0 first vertex of triangle
 * \param[in] Q1 second vertex of triangle
 * \param[in] Q2 third vertex of triangle
 * \param[in] lambda barycentric coordinates
 * \param[out] reconstructed point
 */
array3D reconstructPointFromBarycentricTriangle(array3D const &Q0, array3D const &Q1, array3D const &Q2, array3D const &lambda)
{
    assert( validBarycentric(&lambda[0],3) );

    return lambda[0]*Q0 +lambda[1]*Q1 +lambda[2]*Q2;
}

/*!
 * Reconstructs a point from barycentric coordinates of a triangle
 * \param[in] Q0 first vertex of triangle
 * \param[in] Q1 second vertex of triangle
 * \param[in] Q2 third vertex of triangle
 * \param[in] lambda barycentric coordinates
 * \param[out] reconstructed point
 */
array3D reconstructPointFromBarycentricTriangle(array3D const &Q0, array3D const &Q1, array3D const &Q2, double const *lambda)
{
    assert( validBarycentric(&lambda[0],3) );

    return lambda[0]*Q0 +lambda[1]*Q1 +lambda[2]*Q2;
}

/*!
 * Reconstructs a point from barycentric coordinates of a polygon
 * \param[in] V vertices of simplex
 * \param[in] lambda barycentric coordinates
 * \return reconstructed Point
 */
array3D reconstructPointFromBarycentricPolygon( std::vector<array3D> const &V, std::vector<double> const &lambda)
{
    int N(V.size());
    assert( validBarycentric(&lambda[0],N) );

    array3D xP = {{0.,0.,0.}};
    for(int i=0; i<N; ++i){
        xP += lambda[i]*V[i];
    }

    return xP;
}

/*!
 * Computes projection of point on line in 3D
 * \param[in] P point coordinates
 * \param[in] Q point on line
 * \param[in] n line direction
 * \return projection point
 */
array3D projectPointLine( array3D const &P, array3D const &Q, array3D const &n )
{
    assert( validLine(Q,n) );
    return Q + dotProduct(P - Q, n) * n;
}

/*!
 * Computes projection of point on line in 3D
 * \param[in] P point coordinates
 * \param[in] Q point on plane
 * \param[in] n plane normal
 * \return projection point
 */
array3D projectPointPlane( array3D const &P, array3D const &Q, array3D const &n )
{
    assert( validPlane(Q,n) );
    return P - dotProduct(P - Q, n) * n;
}

/*!
 * Computes projection of point on line in 3D
 * \param[in] P point coordinates
 * \param[in] Q point on plane
 * \param[in] n plane normal
 * \return projection point
 */
array3D projectPointSegment( array3D const &P, array3D const &Q0, array3D const &Q1)
{

    std::array<double,2> lambda;
    return projectPointSegment( P, Q0, Q1, &lambda[0] );
}

/*!
 * Computes projection of point on line in 3D
 * \param[in] P point coordinates
 * \param[in] Q point on plane
 * \param[in] n plane normal
 * \param[out] lambda barycentric coordinates of projection point 
 * \return projection point
 */
array3D projectPointSegment( array3D const &P, array3D const &Q0, array3D const &Q1, std::array<double,2> &lambda )
{
    return projectPointSegment( P, Q0, Q1, &lambda[0]);
}

/*!
 * Computes projection of point on line in 3D
 * \param[in] P point coordinates
 * \param[in] Q point on plane
 * \param[in] n plane normal
 * \param[out] lambda barycentric coordinates of projection point 
 * \return projection point
 */
array3D projectPointSegment( array3D const &P, array3D const &Q0, array3D const &Q1, double *lambda )
{

    assert( validSegment(Q0,Q1) );

    array3D n = Q1 -Q0;
    double t =  -dotProduct(n,Q0-P) / dotProduct(n,n);

    // Restrict projection onto the segment
    t = std::max( std::min( t, 1.), 0. );

    lambda[0] = 1. - t;
    lambda[1] = t;

    return reconstructPointFromBarycentricSegment( Q0, Q1, lambda);
}

/*!
 * Computes projection of point on triangle
 * \param[in] P point coordinates
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \return coordinates of projection point
 */
array3D projectPointTriangle( array3D const &P, array3D const &Q0, array3D const &Q1, array3D const &Q2)
{
    array3D xP;
    array3D lambda;
    _projectPointsTriangle( 1, &P, Q0, Q1, Q2, &xP, lambda.data() );
    return xP;
}

/*!
 * Computes projection of point on triangle
 * \param[in] P point coordinates
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \param[out] lambda barycentric coordinates of projection point
 * \return coordinates of projection point
 */
array3D projectPointTriangle( array3D const &P, array3D const &Q0, array3D const &Q1, array3D const &Q2, array3D &lambda)
{
    array3D xP;
    _projectPointsTriangle( 1, &P, Q0, Q1, Q2, &xP, lambda.data() );

    return xP;
}

/*!
 * Restricts a point given in barycentric coordinates on the triangle
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \param[in,out] lambda barycentric coordinates before and after restriction
 * \return coordinates of restricted point
 */
array3D restrictPointTriangle( array3D const &Q0, array3D const &Q1, array3D const &Q2, array3D &lambda)
{
    return restrictPointTriangle( Q0, Q1, Q2, &lambda[0] );
}

/*!
 * Restricts a point given in barycentric coordinates on the triangle
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \param[in,out] lambda barycentric coordinates before and after restriction
 * \return coordinates of restricted point
 */
array3D restrictPointTriangle( array3D const &Q0, array3D const &Q1, array3D const &Q2, double *lambda)
{

    assert( validBarycentric(&lambda[0], 3) );

    std::array<const array3D*,3> r = {{&Q0, &Q1, &Q2}};

    int count = 0;
    std::array<int,2> negatives = {{ 0, 0 }};

    for( int i=0; i<3; ++i){
        if( lambda[i] < 0){
            negatives[count] = i;
            ++count;
        }
    }

    if( count == 0){
        return reconstructPointFromBarycentricTriangle( Q0, Q1, Q2, lambda );

    } else if( count == 1){
        std::array<double,2>   lambdaLocal;
        int vertex0 = (negatives[0] +1) %3;
        int vertex1 = (vertex0     +1) %3;
        array3D P = reconstructPointFromBarycentricTriangle( Q0, Q1, Q2, lambda ); 
        array3D xP = projectPointSegment(P, *r[vertex0], *r[vertex1], lambdaLocal);
        lambda[negatives[0]] = 0.;
        lambda[vertex0] = lambdaLocal[0];
        lambda[vertex1] = lambdaLocal[1];
        return xP;

    } else {
        int vertex0 = 3 -negatives[0] -negatives[1];
        lambda[0] = 0.;
        lambda[1] = 0.;
        lambda[2] = 0.;
        lambda[vertex0] = 1.;
        return *r[vertex0];

    }

    BITPIT_UNREACHABLE("CANNOT REACH");

}

/*!
 * Projects a point cloud onto a triangle. 
 * Projection points are the closest points to the original points within the triangle.
 * \param[in] cloud point cloud coordinates
 * \param[in] Q1 first triangle vertex
 * \param[in] Q2 second triangle vertex
 * \param[in] Q3 third triangle vertex
 * \param[inout] lambda barycentric coordinates of projection points
 * \return distances
 */
std::vector<array3D> projectCloudTriangle( std::vector<array3D> const &cloud, array3D const &Q0, array3D const &Q1, array3D const &Q2, std::vector<array3D> &lambda )
{
   
    int cloudCount(cloud.size()); 

    std::vector<array3D> xP(cloudCount);

    lambda.resize(cloudCount);

    _projectPointsTriangle( cloudCount, cloud.data(), Q0, Q1, Q2, xP.data(), &lambda[0][0]);

    return xP;

}

/*!
 * Computes projection of point onto a convex polygon
 * \param[in] P point coordinates
 * \param[in] V polygon vertices coordinates
 * \return coordinates of projection point
 */
array3D projectPointPolygon( array3D const &P, std::vector<array3D> const &V)
{
    std::vector<double> lambda(V.size());
    return projectPointPolygon( P, V, lambda);
}

/*!
 * Computes projection of point onto a convex polygon
 * \param[in] P point coordinates
 * \param[in] V polygon vertices coordinates
 * \param[out] lambda baycentric coordinates of projection point
 * \return coordinates of projection point
 */
array3D projectPointPolygon( array3D const &P, std::vector<array3D> const &V, std::vector<double> &lambda)
{
    array3D xP;

    int vertexCount(V.size());
    lambda.resize(vertexCount);

    double distance, minDistance(std::numeric_limits<double>::max());
    int minTriangle = -1;
    array3D V0, V1, V2;
    array3D localLambda, minLambda;

    // Compute the distance from each triangle in the simplex
    int triangleCount = polygonSubtriangleCount(V);

    for (int triangle=0; triangle < triangleCount; ++triangle) {

        subtriangleOfPolygon( triangle, V, V0, V1, V2);

        distance = distancePointTriangle(P, V0, V1, V2, localLambda);

        if (distance <= minDistance) {

            minDistance = distance;

            minLambda = localLambda;
            minTriangle = triangle;
        }

    } //next triangle

    assert(minTriangle >= 0);
    subtriangleOfPolygon( minTriangle, V, V0, V1, V2);
    xP = reconstructPointFromBarycentricTriangle( V0, V1, V2, minLambda);

    computeGeneralizedBarycentric( xP, V, lambda);

    return xP;

}

/*!
 * Computes projection point on semi-infinite cone surface
 * \param[in] point point coordinates
 * \param[in] apex cone apex
 * \param[in] axis cone axis
 * \param[in] alpha cone half angle
 * \return projection point
 */
array3D projectPointCone( array3D const &point, array3D const &apex, array3D const &axis, double const &alpha)
{


    if( alpha <= M_PI/2. ) { //accute cone angle

        array3D versor = point-apex;
        versor /= norm2(versor);

        double cosPointAxis = dotProduct(versor,axis);
        double cosCriticalAngle = cos(alpha+M_PI/2.);

        if( cosPointAxis <= cosCriticalAngle ){ //point projects on cone apex
            return apex;

        } else { // point projects on cone surface

            array3D planeNormal = crossProduct(axis,versor);
            planeNormal /= norm2(planeNormal);

            array3D direction = rotateVector(axis,planeNormal,alpha);

            return projectPointLine(point,apex,direction);

        } 

    } else { // abtuse cone angle -> project on complement
        return projectPointCone( point, apex, -1.*axis, M_PI-alpha);

    }

}

/*!
 * Computes distance point to line in 3D
 * \param[in] P point coordinates
 * \param[in] Q point on line
 * \param[in] n line direction
 * \param[out] xP closest point on line
 * \return distance
 */
double distancePointLine( array3D const &P, array3D const &Q, array3D const &n, array3D &xP) 
{
    xP = projectPointLine(P,Q,n);
    return norm2( P-xP);
}

/*!
 * Computes distance point to plane
 * \param[in] P point coordinates
 * \param[in] Q point on plane
 * \param[in] n plane normal
 * \param[out] xP closest point on line
 * \return distance
 */
double distancePointPlane( array3D const &P, array3D const &Q, array3D const &n, array3D &xP) 
{
    xP = projectPointPlane(P,Q,n);
    return norm2(P-xP);
}

/*!
 * Computes distance point to segment
 * \param[in] P point coordinates
 * \param[in] Q0 segment starting point
 * \param[in] Q1 segment ending point
 * \return distance
 */
double distancePointSegment( array3D const &P, array3D const &Q0, array3D const &Q1)
{ 
    array3D xP = projectPointSegment( P, Q0, Q1);
    return norm2(P-xP); 
}

/*!
 * Computes distance point to segment in 3D using barycentric coordinates
 * \param[in] P point coordinates
 * \param[in] Q1 segment starting point
 * \param[in] Q2 segment ending point
 * \param[out] xP closest point on line
 * \param[out] lambda barycentric coordinates
 * \return distance
 */
double distancePointSegment( array3D const &P, array3D const &Q0, array3D const &Q1, std::array<double,2> &lambda)
{
    array3D xP = projectPointSegment( P, Q0, Q1, lambda);
    return norm2(P-xP); 
}

/*!
 * Computes distance point to triangle
 * \param[in] P point coordinates
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \return distance
 */
double distancePointTriangle( array3D const &P, array3D const &Q0, array3D const &Q1, array3D const &Q2)
{
    array3D xP = projectPointTriangle(P, Q0, Q1, Q2);
    return norm2(P-xP);
}

/*!
 * Computes distance point to triangle
 * \param[in] P point coordinates
 * \param[in] Q0 first triangle vertex
 * \param[in] Q1 second triangle vertex
 * \param[in] Q2 third triangle vertex
 * \param[out] lambda barycentric coordinates of projection point
 * \return distance
 */
double distancePointTriangle( array3D const &P, array3D const &Q0, array3D const &Q1, array3D const &Q2, array3D &lambda)
{
    array3D xP = projectPointTriangle(P, Q0, Q1, Q2, lambda);
    return norm2(P-xP);
}

/*!
 * Computes distance point to semi-infinite cone surface
 * \param[in] point point coordinates
 * \param[in] apex cone apex
 * \param[in] axis cone axis
 * \param[in] alpha cone half angle
 * \return distance
 */
double distancePointCone( array3D const &point, array3D const &apex, array3D const &axis, double const &alpha)
{
    array3D xP = projectPointCone( point, apex, axis, alpha);
    return norm2(point-xP);
}

/*!
 * Computes distances of point cloud to triangle
 * \param[in] cloud point cloud coordinates
 * \param[in] Q1 first triangle vertex
 * \param[in] Q2 second triangle vertex
 * \param[in] Q3 third triangle vertex
 * \return distances
 */
std::vector<double> distanceCloudTriangle( std::vector<array3D> const &cloud, array3D const &Q0, array3D const &Q1, array3D const &Q2)
{ 
    std::vector<array3D> lambda(cloud.size());
    return distanceCloudTriangle( cloud, Q0, Q1, Q2, lambda);
}

/*!
 * Computes distances of point cloud to triangle
 * \param[in] cloud point cloud coordinates
 * \param[in] Q1 first triangle vertex
 * \param[in] Q2 second triangle vertex
 * \param[in] Q3 third triangle vertex
 * \param[out] lambda barycentric coordinates of projection points
 * \return distances
 */
std::vector<double> distanceCloudTriangle( std::vector<array3D> const &cloud, array3D const &Q0, array3D const &Q1, array3D const &Q2, std::vector<array3D> &lambda )
{ 
    int N(cloud.size());

    lambda.resize(N);

    std::vector<array3D> xP(N);
    _projectPointsTriangle( N, cloud.data(), Q0, Q1, Q2, xP.data(), &lambda[0][0]);

    std::vector<double> d(N);
    std::vector<double>::iterator distance = d.begin();

    std::vector<array3D>::iterator projection = xP.begin();
    for( const auto &point : cloud){
        *distance = norm2( point - *projection);

        ++projection;
        ++distance;
    }

    return d;
}

/*!
 * Computes distances of point to a convex polygon
 * \param[in] P point coordinates
 * \param[in] V polygon vertices coordinates
 * \param[out] xP closest points on polygon
 * \param[out] flag point projecting onto polygon's interior (flag = 0), polygon's vertices (flag = 1, 2, ...) or polygon's edges (flag = -1, -2, -...)
 * \return distance
 */
double distancePointPolygon( array3D const &P, std::vector<array3D> const &V, array3D &xP, int &flag)
{
    std::vector<double> lambda(V.size());
    double distance = distancePointPolygon( P, V, lambda);
    xP = reconstructPointFromBarycentricPolygon( V, lambda );
    flag = convertBarycentricToFlagPolygon( lambda );

    return distance; 
}

/*!
 * Computes distances of point to a convex polygon
 * \param[in] P point coordinates
 * \param[in] V simplex vertices coordinates
 * \return distance
 */
double distancePointPolygon( array3D const &P, std::vector<array3D> const &V)
{
    std::vector<double> lambda(V.size());
    return distancePointPolygon( P, V, lambda);
}

/*!
 * Computes distances of point to a convex polygon
 * \param[in] P point coordinates
 * \param[in] V polygon vertices coordinates
 * \param[out] lambda barycentric coordinates
 * \return distance
 */
double distancePointPolygon( array3D const &P, std::vector<array3D> const &V,std::vector<double> &lambda)
{
    array3D xP = projectPointPolygon( P, V, lambda);
    return norm2(P-xP); 
}

/*!
 * Computes distances of point cloud to a convex polygon
 * \param[in] P point cloud coordinates
 * \param[in] V polygon vertices coordinates
 * \param[out] xP closest points on simplex
 * \param[out] flag point projecting onto polygon's interior (flag = 0), polygon's vertices (flag = 1, 2, ...) or polygon's edges (flag = -1, -2, -...)
 * \return distance
 */
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &cloud, std::vector<array3D> const &V, std::vector<array3D> &xP, std::vector<int> &flag)
{
    int cloudCount( cloud.size() );
    int vertexCount( V.size() );

    std::vector<std::vector<double>> lambda( cloudCount, std::vector<double> (vertexCount));
    std::vector<double> d = distanceCloudPolygon( cloud, V, lambda);

    xP.resize(cloudCount);
    std::vector<array3D>::iterator xPItr = xP.begin();

    flag.resize(cloudCount);
    std::vector<int>::iterator flagItr = flag.begin();

    for( const auto &l : lambda){
        *flagItr = convertBarycentricToFlagPolygon( l );
        *xPItr = reconstructPointFromBarycentricPolygon( V, l ); 

        ++xPItr;
        ++flagItr;
    }

    return d; 
}

/*!
 * Computes distances of point cloud to a convex polygon
 * \param[in] P point cloud coordinates
 * \param[in] V polygon vertices coordinates
 * \return distance
 */
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &P, std::vector<array3D> const &V)
{
    int cloudCount(P.size());

    std::vector<double> d(cloudCount,std::numeric_limits<double>::max());

    int triangleCount = polygonSubtriangleCount(V);
    array3D V0, V1, V2;

    for (int triangle=0; triangle < triangleCount; ++triangle) { // foreach triangle
        subtriangleOfPolygon( triangle, V, V0, V1, V2);
        std::vector<double> dT = distanceCloudTriangle(P, V0, V1, V2);

        d = min(d,dT);
    }

    return d; 
}

/*!
 * Computes distances of point cloud to a convex polygon
 * \param[in] cloud point cloud coordinates
 * \param[in] V polygon vertices coordinates
 * \param[out] lambda barycentric coordinates of the projection points
 * \return distance
 */
std::vector<double> distanceCloudPolygon( std::vector<array3D> const &cloud, std::vector<array3D> const &V, std::vector<std::vector<double>> &lambda)
{
    int cloudCount(cloud.size()), vertexCount(V.size());

    std::vector<double> d(cloudCount,std::numeric_limits<double>::max());
    lambda.resize(cloudCount, std::vector<double>(vertexCount,0) );

    std::vector<double> dTemp(cloudCount);
    std::vector<array3D> lambdaTemp(cloudCount);
    int triangleCount = polygonSubtriangleCount(V);
    array3D V0, V1, V2;

    for (int triangle=0; triangle < triangleCount; ++triangle) { // foreach triangle
        subtriangleOfPolygon( triangle, V, V0, V1, V2);
        dTemp = distanceCloudTriangle(cloud, V0, V1, V2, lambdaTemp);

        for(int i=0; i< cloudCount; ++i){
            if( dTemp[i] < d[i]){
                d[i] = dTemp[i];
                std::copy( lambdaTemp[i].begin(), lambdaTemp[i].end(), lambda[i].begin());
            }
        }
    }

    return d; 
}

/*!
 * Computes distance between two lines in 3D
 * \param[in] n0 direction of first line
 * \param[in] P0 point on first line
 * \param[in] n1 direction of second line
 * \param[in] P1 point on second line
 * \return distance
 */
double distanceLineLine( array3D const &P0, array3D const &n0, array3D const &P1, array3D const &n1)
{
    array3D xP0, xP1;
    return distanceLineLine( P0, n0, P1, n1, xP0, xP1);
}

/*!
 * Computes distance between two lines in 3D
 * \param[in] n0 direction of first line
 * \param[in] P0 point on first line
 * \param[in] n1 direction of second line
 * \param[in] P1 point on second line
 * \param[out] xP0 projection of line1 on line0
 * \param[out] xP1 projection of line0 on line1
 * \return distance
 */
double distanceLineLine( array3D const &P0, array3D const &n0, array3D const &P1, array3D const &n1, array3D &xP0, array3D &xP1)
{
    assert( validLine(P0,n0) );
    assert( validLine(P1,n1) );

    double n01 = dotProduct(n0,n1);
    double det = 1. - n01*n01;

    // check if lines are parallel 
    if( std::abs(det) < 1.e-12){
        double distance = distancePointLine(P0, P1, n1, xP1);
        xP0 = projectPointLine(xP1, P0, n0);
        return distance;
    }


    array3D dP = P1-P0;
    double rhs0 =  dotProduct(dP,n0);
    double rhs1 = -dotProduct(dP,n1);

    double det0 = rhs0 +rhs1*n01;
    double det1 = rhs1 +rhs0*n01;

    double s0 = det0/det;
    double s1 = det1/det;

    xP0 = P0 +s0*n0;
    xP1 = P1 +s1*n1;

    return norm2( xP0 - xP1); 
}

/*!
 * Computes intersection between two lines in 3D
 * \param[in] n1 direction of first line
 * \param[in] P1 point on first line
 * \param[in] n2 direction of second line
 * \param[in] P2 point on second line
 * \param[out] P intersection point if intersect, else unaltered
 * \return if intersect
 */
bool intersectLineLine( array3D const &P1, array3D const &n1, array3D const &P2, array3D const &n2, array3D &P)
{
    double tol = 1.e-12;
    array3D xP1, xP2;
    if( distanceLineLine(P1,n1,P2,n2,xP1,xP2) < tol){
        P = xP1;
        return true;
    }

    return false;
}

/*!
 * Computes intersection between two segments in 3D
 * \param[in] P1 start point of first segment
 * \param[in] P2 end point of first segment
 * \param[in] Q1 start point of second segment
 * \param[in] Q2 end point of second segment
 * \param[out] x intersection point if intersect, else unaltered
 * \return if intersect
 */
bool intersectSegmentSegment( array3D const &P1, array3D const &P2, array3D const &Q1, array3D const &Q2, array3D &x)
{
    assert( validSegment(P1,P2) );
    assert( validSegment(Q1,Q2) );

    array3D nP = P2 - P1;
    nP /= norm2(nP);

    array3D nQ = Q2 - Q1;
    nQ /= norm2(nQ);

    array3D temp;
    if( intersectLineLine(P1, nP, Q1, nQ, temp) && intersectPointSegment( temp, P1, P2) && intersectPointSegment(temp, Q1, Q2) ){
        x = temp;
        return true;
    }

    return false;
}

/*!
 * Computes intersection between line and plane
 * \param[in] P1 point on line
 * \param[in] n1 direction of line
 * \param[in] P2 point on plane
 * \param[in] n2 normal to plane
 * \param[out] P intersection point if intersect, else unaltered
 * \return if intersect
 */
bool intersectLinePlane( array3D const &P1, array3D const &n1, array3D const &P2, array3D const &n2, array3D &P)
{

    assert( validLine(P1,n1) );
    assert( validPlane(P2,n2) );

    double const tol = 1.0e-14;

    // ========================================================================== //
    // CHECK DEGENERATE CASES                                                     //
    // ========================================================================== //
    double s = dotProduct(n1, n2);
    if (std::abs(s) < tol) { 
        return false; 
    }

    // ========================================================================== //
    // FIND INTERSECTION POINTS                                                   //
    // ========================================================================== //
    double xi = -dotProduct(P1 - P2, n2) /s;
    P = P1 + xi *n1;

    return true; 
}

/*!
 * Computes intersection between segment and plane
 * \param[in] Q1 start point of segment
 * \param[in] Q2 end point of segment
 * \param[in] P2 point on plane
 * \param[in] n2 normal to plane
 * \param[out] P intersection point if intersect, else unaltered
 * \return if intersect
 */
bool intersectSegmentPlane( array3D const &Q1, array3D const &Q2, array3D const &P2, array3D const &n2, array3D &P)
{

    assert( validSegment(Q1,Q2) );
    assert( validPlane(P2,n2) );

    array3D n = Q2 - Q1;
    n /= norm2(n);

    array3D xP;
    if ( intersectLinePlane(Q1, n, P2, n2, xP) && intersectPointSegment(xP, Q1, Q2) ) {
        P = xP;
        return true;
    }

    return false; 
}

/*!
 * Computes intersection between two planes
 * \param[in] P1 point on first plane
 * \param[in] n1 normal to first plane
 * \param[in] P2 point on second plane
 * \param[in] n2 normal to second plane
 * \param[out] Pl point on intersection line
 * \param[out] nl direction of intersection line
 * \return if intersect
 */
bool intersectPlanePlane( array3D const &P1, array3D const &n1, array3D const &P2, array3D const &n2, array3D &Pl, array3D &nl)
{

    assert( validPlane(P1,n1) );
    assert( validPlane(P2,n2) );

    double const tol = 1.0e-14;
    double n12 = dotProduct(n1, n2);
    double detCB = 1.0-n12*n12;

    // check degenerate condition
    if( std::abs(detCB) <= tol) { 
        return false; 
    }

    nl = crossProduct(n1,n2);
    nl /= norm2(nl);

    // if planes intersect, determine the closest point
    // to P1 and P2 as anchor point. The augmented functional
    // I = 0.5*[ (Pl-P1)^2 + (Pl-P2)^2] + lambda1[ n1.(Pl-P1) ] +lambda2[ n2.(Pl-P2) ]
    // where lambda1 and lambda2 are Lagrange multipliers.
    // The optimality conditions I,Pl I,lambda1 I,lambda2 are 
    // solved using the Schur complment

    array3D  dP = P2-P1;
    std::array<double,2>  rhs = {{ dotProduct(n1,dP) , -dotProduct(n2,dP) }};

    double det1 = rhs[0] - n12*rhs[1];
    double det2 = rhs[1] - n12*rhs[0];
    double lambda1 = det1 /detCB;
    double lambda2 = det2 /detCB;

    Pl = P1 +P2 -lambda1*n1 -lambda2*n2;
    Pl *= 0.5;

    return true; 

}

/*!
 * Computes if plane a box intersect
 * \param[in] P1 point on first plane
 * \param[in] n1 normal to first plane
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] dim number of dimensions
 * \return if intersect
 */
bool intersectPlaneBox( array3D const &P1, array3D const &n1, array3D const &A0, array3D const &A1, int dim)
{
    assert( validPlane(P1,n1) );
    return _intersectPlaneBox( P1, n1, A0, A1, nullptr, dim);
}

/*!
 * Computes if plane a box intersect and intersection points between plane and edges of box
 * \param[in] P1 point on first plane
 * \param[in] n1 normal to first plane
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[out] intersects intersection points
 * \param[in] dim number of dimensions
 * \return if intersect
 */
bool intersectPlaneBox( array3D const &P1, array3D const &n1, array3D const &A0, array3D const &A1, std::vector<array3D> &intersects, int dim)
{
    assert( validPlane(P1,n1) );
    return _intersectPlaneBox( P1, n1, A0, A1, &intersects, dim);
}

/*!
 * Computes intersection between triangle and a line
 * \param[in] P point on plane
 * \param[in] n normal to plane
 * \param[in] A first vertex of triangle
 * \param[in] B second vertex of triangle
 * \param[in] C third vertex of triangle
 * \param[out] Q intersection point
 * \return if intersect
 */
bool intersectLineTriangle( array3D const &P, array3D const &n, array3D const &A, array3D const &B, array3D const &C, array3D &Q)
{
    assert( validLine(P,n) );
    assert( validTriangle(A,B,C) );

    array3D nT  = crossProduct(B - A, C - A);
    nT /= norm2(nT);

    array3D xP;
    if ( intersectLinePlane(P, n, A, nT, xP) && intersectPointTriangle(xP, A, B, C) ) { 
        Q = xP;
        return true; 
    }

    return false; 
}

/*!
 * Computes intersection between triangle and a segment
 * \param[in] P0 start point of segment
 * \param[in] P1 end point of segment
 * \param[in] A first vertex of triangle
 * \param[in] B second vertex of triangle
 * \param[in] C third vertex of triangle
 * \param[out] Q intersection point
 * \return if intersect
 */
bool intersectSegmentTriangle( array3D const &P0, array3D const &P1, array3D const &A, array3D const &B, array3D const &C, array3D &Q)
{
    assert( validSegment(P0,P1) );
    assert( validTriangle(A,B,C) );

    array3D n = P1 - P0;
    n /= norm2(n);

    array3D xP;
    if ( intersectLineTriangle(P0, n, A, B, C, xP) && intersectPointSegment(xP, P0, P1)  ) { 
        Q = xP;
        return true; 
    }

    return false; 
}

/*!
 * Computes intersection between triangle and a convex polygon
 * \param[in] P point on line
 * \param[in] n direction of line
 * \param[in] V polygon vertices coordinates
 * \param[out] Q intersection point
 * \return if intersect
 */
bool intersectLinePolygon( array3D const &P, array3D const &n, std::vector<array3D > const &V, array3D &Q)
{
    assert( validLine(P,n) );

    int nTriangles = polygonSubtriangleCount(V);
    array3D V0, V1, V2;

    for( int i=0; i< nTriangles; ++i){
        subtriangleOfPolygon(i, V, V0, V1, V2);

        if( intersectLineTriangle(P, n, V0, V1, V2, Q) ) { 
            return true; 
        }
    }

    return false; 
}

/*!
 * Computes intersection between a segment and a polygon
 * \param[in] P0 start point of segment
 * \param[in] P1 end point of segment
 * \param[in] V polygon vertices coordinates
 * \param[out] Q intersection point
 * \return if intersect
 */
bool intersectSegmentPolygon( array3D const &P0, array3D const &P1, std::vector<array3D > const &V, array3D &Q)
{
    assert( validSegment(P0,P1) );

    int nTriangles = polygonSubtriangleCount(V);
    array3D V0, V1, V2;

    for( int i=0; i< nTriangles; ++i){
        subtriangleOfPolygon(i, V, V0, V1, V2);

        if( intersectSegmentTriangle(P0, P1, V0, V1, V2, Q) ) { 
            return true; 
        }
    }

    return false;
}

/*!
 * Computes intersection between two axis aligned bounding boxes
 * \param[in] A1 min point of first box
 * \param[in] A2 max point of first box
 * \param[in] B1 min point of second box
 * \param[in] B2 max point of second box
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectBoxBox(array3D const &A1, array3D const &A2, array3D const &B1, array3D const &B2, int dim)
{
    for( int d=0; d<dim; ++d){
        if( B1[d] > A2[d] || B2[d] < A1[d] ){
            return false;
        }
    }

    return true;
}

/*!
 * Computes intersection between two axis aligned bounding boxes
 * \param[in] A1 min point of first box
 * \param[in] A2 max point of first box
 * \param[in] B1 min point of second box
 * \param[in] B2 max point of second box
 * \param[out] I1 min point of intersection box
 * \param[out] I2 max point of intersection box
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectBoxBox(array3D const &A1, array3D const &A2, array3D const &B1, array3D const &B2, array3D &I1, array3D &I2, int dim )
{
    for( int d=0; d<dim; ++d){

        if( B1[d] > A2[d] || B2[d] < A1[d] ){
            return false;
        }

        else{
            I1[d] = std::max( A1[d], B1[d] );
            I2[d] = std::min( A2[d], B2[d] );
        }

    }

    return true;
}

/*!
 * Checks if an axis aligned bounding box and a triangle intersect
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] V0 first vertex of triangle
 * \param[in] V1 second vertex of triangle
 * \param[in] V2 third vertex of triangle
 * \param[in] dim number of dimensions to be checked
 * \return if intersect
 */
bool intersectBoxTriangle(array3D const &A0, array3D const &A1, array3D const &V0, array3D const &V1, array3D const &V2, int dim)
{
    return _intersectBoxTriangle( A0, A1, V0, V1, V2, false, false, false, nullptr, nullptr, dim);
}

/*!
 * Checks if an axis aligned bounding box and a triangle intersect
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] V0 first vertex of triangle
 * \param[in] V1 second vertex of triangle
 * \param[in] V2 third vertex of triangle
 * \param[in] interiorTriangleVertice if true, all triangle vertices within the box will be added to the intersection points P
 * \param[in] triangleEdgeBoxFaceIntersections if true, the intersections between the edges of the triangle and the hull of the box will be added to the intersection points P
 * \param[in] triangleBoxEdgeIntersection if true, the intersections between the edges (dim=3) or vertices (dim=2) of the box and the triangle will be added to the intersection points P
 * \param[out] P calculated intersection points 
 * \param[in] dim number of dimensions to be checked
 * \return if intersect
 */
bool intersectBoxTriangle(array3D const &A0, array3D const &A1, array3D const &V0, array3D const &V1, array3D const &V2, bool interiorTriangleVertice, bool triangleEdgeBoxFaceIntersections, bool triangleBoxEdgeIntersection, std::vector<array3D> &P, int dim)
{
    return _intersectBoxTriangle( A0, A1, V0, V1, V2, interiorTriangleVertice, triangleEdgeBoxFaceIntersections, triangleBoxEdgeIntersection, &P, nullptr, dim);
}

/*!
 * Checks if an axis aligned bounding box and a triangle intersect
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] V0 first vertex of triangle
 * \param[in] V1 second vertex of triangle
 * \param[in] V2 third vertex of triangle
 * \param[in] interiorTriangleVertice if true, all triangle vertices within the box will be added to the intersection points P
 * \param[in] triangleEdgeBoxHullIntersections if true, the intersections between the edges of the triangle and the outer hull of the box will be added to the intersection points P
 * \param[in] triangleBoxEdgeIntersection if true, the intersections between the edges (dim=3) or vertices (dim=2) of the box and the triangle will be added to the intersection points P
 * \param[out] P calculated intersection points 
 * \param[out] flag has the same size of P. If the ith flag=0, the intersection is due to interiorTriangleVertice. If the ith flag=1, the intersection is due to triangleEdgeBoxHullIntersection. If the ith flag=2, the intersection is due to triangleBoxEdgeIntersections.
 * \param[in] dim number of dimensions to be checked
 * \return if intersect
 */
bool intersectBoxTriangle(array3D const &A0, array3D const &A1, array3D const &V0, array3D const &V1, array3D const &V2, bool interiorTriangleVertice, bool triangleEdgeBoxHullIntersections, bool triangleBoxEdgeIntersection, std::vector<array3D> &P, std::vector<int> &flag, int dim)
{
    return _intersectBoxTriangle( A0, A1, V0, V1, V2, interiorTriangleVertice, triangleEdgeBoxHullIntersections, triangleBoxEdgeIntersection, &P, &flag, dim);
}

/*!
 * Computes intersection between an axis aligned bounding box and a segment
 * \param[in] V0 start point of segment
 * \param[in] V1 end point of segment
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectSegmentBox( array3D const &V0, array3D const &V1, array3D const &A0, array3D const &A1, int dim)
{
    return _intersectSegmentBox( V0, V1, A0, A1, false, false, nullptr, nullptr, dim);
}

/*!
 * Computes intersection between an axis aligned bounding box and a segment
 * \param[in] V0 start point of segment
 * \param[in] V1 end point of segment
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] interiorSegmentVertice if the segment vertices within the box should be added to the list of intersection points
 * \param[in] segmentBoxHullIntersection if the intersections between the segment and the outer hull of the box should be added to the list of intersection points
 * \param[out] P list of intersection points. If no intersetion points should be calculated nullptr can be passed as argument
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectSegmentBox( array3D const &V0, array3D const &V1, array3D const &A0, array3D const &A1, bool interiorSegmentVertice, bool segmentBoxHullIntersection, std::vector<array3D> &P, int dim)
{
    return _intersectSegmentBox( V0, V1, A0, A1, interiorSegmentVertice, segmentBoxHullIntersection, &P, nullptr, dim);
}

/*!
 * Computes intersection between an axis aligned bounding box and a segment
 * \param[in] V0 start point of segment
 * \param[in] V1 end point of segment
 * \param[in] A0 min point of box
 * \param[in] A1 max point of box
 * \param[in] interiorSegmentVertice if the segment vertices within the box should be added to the list of intersection points
 * \param[in] segmentBoxHullIntersection if the intersections between the segment and the outer hull of the box should be added to the list of intersection points
 * \param[out] P list of intersection points. If no intersetion points should be calculated nullptr can be passed as argument
 * \param[out] flag indicates for each intersection if it belongs to interiorSegmentVertice (flag=0) or segmentHullIntersection (flag=1)
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectSegmentBox( array3D const &V0, array3D const &V1, array3D const &A0, array3D const &A1, bool interiorSegmentVertice, bool segmentBoxHullIntersection, std::vector<array3D> &P, std::vector<int> &flag, int dim)
{
    return _intersectSegmentBox( V0, V1, A0, A1, interiorSegmentVertice, segmentBoxHullIntersection, &P, &flag, dim);
}

/*!
 * Computes intersection between an axis aligned bounding box and a simplex
 * \param[in] A0 min point of first box
 * \param[in] A1 max point of first box
 * \param[in] VS simplex vertices coordinates
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectBoxPolygon( array3D const &A0, array3D const &A1, std::vector<array3D> const &VS, int dim )
{
    return _intersectBoxPolygon(A0, A1, VS, false, false, false, nullptr, nullptr, dim);
}

/*!
 * Computes intersection between an axis aligned bounding box and a convex polygon
 * \param[in] A0 min point of first box
 * \param[in] A1 max point of first box
 * \param[in] innerPolygonPoints simplex vertices within the box should be added to the intersection list
 * \param[in] polygonEdgeBoxFaceIntersection intersection between the edges of the polygon and the hull of the box should be added to the intersection list
 * \param[in] polygonBoxEdgeIntersection intersection between the polygon and the edges of the box should be added to the intersection list
 * \param[in] VS polygon vertices coordinates
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectBoxPolygon( array3D const &A0, array3D const &A1, std::vector<array3D> const &VS, bool innerPolygonPoints, bool polygonEdgeBoxFaceIntersection, bool polygonBoxEdgeIntersection, std::vector<array3D> &P, int dim)
{
    return _intersectBoxPolygon(A0, A1, VS, innerPolygonPoints, polygonEdgeBoxFaceIntersection, polygonBoxEdgeIntersection, &P, nullptr, dim);
}

/*!
 * Computes intersection between an axis aligned bounding box and a simplex
 * \param[in] A0 min point of first box
 * \param[in] A1 max point of first box
 * \param[in] VS simplex vertices coordinates
 * \param[in] innerPolygonPoints simplex vertices within the box should be added to the intersection list
 * \param[in] polygonEdgeBoxFaceIntersection intersection between the edges of the polygon and the hull of the box should be added to the intersection list
 * \param[in] polygonBoxEdgeIntersection intersection between the polygon and the edges of the box should be added to the intersection list
 * \param[in] dim number of dimension to be checked
 * \return if intersect
 */
bool intersectBoxPolygon( array3D const &A0, array3D const &A1, std::vector<array3D> const &VS, bool innerPolygonPoints, bool polygonEdgeBoxFaceIntersection, bool polygonBoxEdgeIntersection, std::vector<array3D> &P, std::vector<int> &flag, int dim)
{
    return _intersectBoxPolygon(A0, A1, VS, innerPolygonPoints, polygonEdgeBoxFaceIntersection, polygonBoxEdgeIntersection, &P, &flag, dim);
}


//to levelset // -------------------------------------------------------------------------- //
//to levelset bool intersectLineSurface(
//to levelset         array3D  const  &x1,
//to levelset         array3D  const  &n1,
//to levelset         array3D  const  &x2,
//to levelset         array3D  const  &n2,
//to levelset         array3D  const  &xL,
//to levelset         array3D  const  &nL,
//to levelset         array3D         &xp,
//to levelset         array3D         &np
//to levelset         ) {
//to levelset 
//to levelset     // ========================================================================== //
//to levelset     //                                                                            //
//to levelset     // Reconstruct Surface given by two points and their normals and computes     //
//to levelset     // intersection of reconstruction with a line given by point and versor     ' //
//to levelset     // ========================================================================== //
//to levelset     // INPUT                                                                      //
//to levelset     // ========================================================================== //
//to levelset     // ========================================================================== //
//to levelset     // OUTPUT                                                                     //
//to levelset     // ========================================================================== //
//to levelset     // - none                                                                     //
//to levelset     // ========================================================================== //
//to levelset 
//to levelset     // ========================================================================== //
//to levelset     // VARIABLES DECLARATION                                                      //
//to levelset     // ========================================================================== //
//to levelset 
//to levelset     double                  w1(0), w2(0), w(0);
//to levelset     array3D         c1, c2;
//to levelset     bool                    s1, s2;
//to levelset 
//to levelset     s1      =   intersectLinePlane( xL, nL, x1, n1, c1);
//to levelset     s2      =   intersectLinePlane( xL, nL, x2, n2, c2);
//to levelset 
//to levelset     if( s1 && s2) {
//to levelset         w1      =   norm2( c2-x2);
//to levelset         w2      =   norm2( c1-x1);
//to levelset         w       =   w1+w2;
//to levelset 
//to levelset         w1      =   w1 /w;
//to levelset         w2      =   w2 /w;
//to levelset 
//to levelset         xp      =   w1*c1 +w2*c2;
//to levelset         np      =   w1*n1 +w2*n2;
//to levelset 
//to levelset         np      =   np /norm2(np);
//to levelset 
//to levelset         return (true);
//to levelset     }
//to levelset 
//to levelset     else if( s1 ){
//to levelset         xp      =   c1;
//to levelset         np      =   n1;
//to levelset 
//to levelset         return(true);
//to levelset     }
//to levelset 
//to levelset     else if( s2 ){
//to levelset         xp      =   c2;
//to levelset         np      =   n2;
//to levelset 
//to levelset         return(true);
//to levelset     }
//to levelset 
//to levelset     return (false); 
//to levelset 
//to levelset }

/*!
 * checks if points lies on segment in 3D
 * \param[in] P point coordinates
 * \param[in] P1 start point of segment
 * \param[in] P2 end point of segment
 * \return if point lies on segment
 */
bool intersectPointSegment( array3D const &P, array3D const &P1, array3D const &P2)
{
    assert( validSegment(P1,P2) );


    if( utils::DoubleFloatingEqual()( norm2(P-P1), 0.) ){
        return true;
    }

    if( utils::DoubleFloatingEqual()( norm2(P-P2), 0.) ){
        return true;
    }

    array3D segmentDir = P2 - P1;
    double segmentLength = norm2(segmentDir);
    segmentDir /= segmentLength;

    array3D versor = P - P1;
    double length = norm2(versor);
    versor /= length;

    if( length>segmentLength ){
        return false;
    }
    
    if( !utils::DoubleFloatingEqual()( dotProduct(versor,segmentDir), 1.) ){
        return false;
    }

    return true; 

}

/*!
 * checks if points lies on triangle
 * \param[in] P point coordinates
 * \param[in] A first vertex of triangle
 * \param[in] B second vertex of triangle
 * \param[in] C third vertex of triangle
 * \return if point lies on triangle
 */
bool intersectPointTriangle( array3D const &P, array3D const &A, array3D const &B, array3D const &C)
{

    array3D xP, lambda;

    _projectPointsPlane( 1, &P, A, B, C, &xP, &lambda[0]);

    for( const auto &l : lambda){
        if(l<0.) {
            return false;
        }
    }

    return true;
}

/*!
 * checks if points lies within axis aligned box
 * \param[in] P point coordinates
 * \param[in] B1 min coodinates of box
 * \param[in] B2 max coodinates of box
 * \param[in] dim number of dimensions to be checked
 * \return if point in box
 */
bool intersectPointBox( array3D const &P, array3D const &B1, array3D const &B2, int dim)
{

    for( int d=0; d<dim; ++d){

        if( P[d]< B1[d] || P[d] > B2[d] ){
            return false;
        }
    }

    return true;
}

/*!
 * computes axis aligned boundig box of a segment
 * \param[in] A start point of segment
 * \param[in] B end point of segment
 * \param[out] P0 min point of bounding box
 * \param[out] P1 max point of bounding box
 */
void computeAABBSegment(array3D const &A, array3D const &B, array3D &P0, array3D &P1)
{

    P0 = A;
    P1 = A;

    for(int i=0; i<3; ++i){
        P0[i] = std::min( P0[i], B[i] );
        P1[i] = std::max( P1[i], B[i] );
    }

    return;
}

/*!
 * computes axis aligned boundig box of a triangle
 * \param[in] A first vertex of triangle
 * \param[in] B second vertex of triangle
 * \param[in] C third vertex of triangle
 * \param[out] P0 min point of bounding box
 * \param[out] P1 max point of bounding box
 */
void computeAABBTriangle(array3D const &A, array3D const &B, array3D const &C, array3D &P0, array3D &P1)
{
    P0 = A;
    P1 = A;

    for(int i=0; i<3; ++i){
        P0[i] = std::min( P0[i], B[i] );
        P0[i] = std::min( P0[i], C[i] );

        P1[i] = std::max( P1[i], B[i] );
        P1[i] = std::max( P1[i], C[i] );
    }

    return;
}

/*!
 * computes axis aligned boundig box of a polygon
 * \param[in] VS polygon vertices coordinates
 * \param[out] P0 min point of bounding box
 * \param[out] P1 max point of bounding box
 */
void computeAABBPolygon(std::vector<array3D> const &VS, array3D &P0, array3D &P1)
{
    int  vertexCount(VS.size());

    P0 = VS[0];
    P1 = VS[0];

    for( int j=1; j<vertexCount; ++j){
        for( int i=0; i<3; ++i){
            P0[i] = std::min( P0[i], VS[j][i] );
            P1[i] = std::max( P1[i], VS[j][i] );
        }
    }

    return;
}

/*!
 * computes the bounding box of the union of two axis aligned bounding boxes
 * \param[in] A0 min point of first bounding box
 * \param[in] A1 max point of first bounding box
 * \param[in] B0 min point of second bounding box
 * \param[in] B1 max point of second bounding box
 * \param[out] C0 min point of union bounding box
 * \param[out] C1 max point of union bounding box
 */
void unionAABB( array3D const &A0, array3D const &A1, array3D const &B0, array3D const &B1, array3D &C0, array3D &C1)
{
    for( int i=0; i<3; ++i){
        C0[i] = std::min( A0[i], B0[i]);
        C1[i] = std::max( A1[i], B1[i]);
    }

    return;
}

/*!
 * computes the bounding box of the union of several axis aligned bounding boxes
 * \param[in] A0 min points of bounding boxes
 * \param[in] A1 max points of bounding boxes
 * \param[out] C0 min point of union bounding box
 * \param[out] C1 max point of union bounding box
 */
void unionAABB(std::vector<array3D> const &A0, std::vector<array3D> const &A1, array3D &C0, array3D &C1)
{

    int n( std::min(A0.size(), A1.size() ) );

    if( n > 0 ){
        C0 =  A0[0];
        C1 =  A1[0];

        for( int i=1; i<n; ++i){
            unionAABB( A0[i], A1[i], C0, C1, C0, C1);
        }
    }

    return;
}

/*!
 * computes the bounding box of the intersection of two axis aligned bounding boxes
 * \param[in] A0 min point of first bounding box
 * \param[in] A1 max point of first bounding box
 * \param[in] B0 min point of second bounding box
 * \param[in] B1 max point of second bounding box
 * \param[out] C0 min point of intersection of boxes
 * \param[out] C1 max point of intersection of boxes
 */
void intersectionAABB(array3D const &A0, array3D const &A1, array3D const &B0, array3D const &B1, array3D &C0, array3D  &C1)
{
    intersectBoxBox( A0, A1, B0, B1, C0, C1 );
    return;
}

/*!
 * computes the bounding box of the relative complement two axis aligned bounding boxes
 * \param[in] A0 min point of first bounding box
 * \param[in] A1 max point of first bounding box
 * \param[in] B0 min point of second bounding box
 * \param[in] B1 max point of second bounding box
 * \param[out] C0 min point of relative complement
 * \param[out] C1 max point of relative complement
 */
void subtractionAABB(array3D const &A0, array3D const &A1, array3D const &B0, array3D const &B1, array3D &C0, array3D  &C1)
{
    // X direction
    if( B0[1]<=A0[1] && B0[2]<=A0[2] && B1[1]>=A1[1] && B1[2]>=A1[2] ){
        C0[0] = ( B0[0]<=A0[0] && B1[0]>=A0[0] ) ? B1[0] : A0[0];
        C1[0] = ( B0[0]<=A1[0] && B1[0]>=A1[0] ) ? B0[0] : A1[0];
    }

    // Y direction
    if( B0[0]<=A0[0] && B0[2]<=A0[2] && B1[0]>=A1[0] && B1[2]>=A1[2] ){
        C0[1] = ( B0[1]<=A0[1] && B1[1]>=A0[1] ) ? B1[1] : A0[2];
        C1[1] = ( B0[1]<=A1[1] && B1[1]>=A1[1] ) ? B0[1] : A1[2];
    }

    // Z direction
    if( B0[0]<=A0[0] && B0[1]<=A0[1] && B1[0]>=A1[0] && B1[1]>=A1[1] ){
        C0[2] = ( B0[2]<=A0[2] && B1[2]>=A0[2] ) ? B1[2] : A0[2];
        C1[2] = ( B0[2]<=A1[2] && B1[2]>=A1[2] ) ? B0[2] : A1[2];
    }

    return;
}

/*!
 * returns the vertex coordiantes of a segment given its index
 * \param[in] i edge index
 * \param[in] V0 first vertex of segment
 * \param[in] V1 second vertex of segment
 * \param[out] P vertex
 */
void vertexOfSegment(int const &i, array3D const &V0, array3D const &V1, array3D &P)
{
    switch(i){

        case 0:
            P = V0;
            break;

        case 1:
            P = V1;
            break;

        default:
            assert(false);
            break;
    }

    return;
}

/*!
 * returns the vertex coordiantes of a triangle given its index
 * \param[in] i edge index
 * \param[in] V0 first vertex of triangle
 * \param[in] V1 second vertex of triangle
 * \param[in] V2 third vertex of triangle
 * \param[out] P vertex
 */
void vertexOfTriangle(int const &i, array3D const &V0, array3D const &V1, array3D const &V2, array3D &P)
{
    switch(i){

        case 0:
            P = V0;
            break;

        case 1:
            P = V1;
            break;

        case 2:
            P = V2;
            break;

        default:
            assert(false);
            break;
    }

    return;
}

/*!
 * computes the edge coordiantes of a triangle
 * \param[in] i edge index
 * \param[in] V0 first vertex of triangle
 * \param[in] V1 second vertex of triangle
 * \param[in] V2 third vertex of triangle
 * \param[out] P0 first vertex of edge
 * \param[out] P1 first vertex of edge
 */
void edgeOfTriangle(int const &i, array3D const &V0, array3D const &V1, array3D const &V2, array3D &P0, array3D &P1)
{
    switch(i){

        case 0:
            P0 = V0;
            P1 = V1;
            break;

        case 1:
            P0 = V1;
            P1 = V2;
            break;

        case 2:
            P0 = V2;
            P1 = V0;
            break;

        default:
            assert(false);
            break;
    }

    return;
}

/*!
 * computes the face coordiantes of a box
 * \param[in] i face index
 * \param[in] A0 min point of bounding box
 * \param[in] A1 max point of bounding box
 * \param[out] P0 first vertex of face
 * \param[out] P1 first vertex of face
 * \param[out] P2 first vertex of face
 * \param[out] P3 first vertex of face
 */
void faceOfBox(int const &i, array3D const &A0, array3D const &A1, array3D &P0, array3D &P1, array3D &P2, array3D &P3)
{

    assert(i<6);

    int v0 = boxFaceVertexConnectivity[i][0];
    int v1 = boxFaceVertexConnectivity[i][1];
    int v2 = boxFaceVertexConnectivity[i][2];
    int v3 = boxFaceVertexConnectivity[i][3];

    vertexOfBox(v0, A0, A1, P0 );
    vertexOfBox(v1, A0, A1, P1 );
    vertexOfBox(v2, A0, A1, P2 );
    vertexOfBox(v3, A0, A1, P3 );

    return;
}

/*!
 * computes the edge coordiantes of a box
 * \param[in] i edge index
 * \param[in] A0 min point of bounding box
 * \param[in] A1 max point of bounding box
 * \param[out] P0 first vertex of edge
 * \param[out] P1 first vertex of edge
 */
void edgeOfBox(int const &i, array3D const &A0, array3D const &A1, array3D &P0, array3D &P1)
{
    assert(i<12);

    int v0 = boxEdgeVertexConnectivity[i][0];
    int v1 = boxEdgeVertexConnectivity[i][1];

    vertexOfBox(v0, A0, A1, P0 );
    vertexOfBox(v1, A0, A1, P1 );

    return;
}

/*!
 * computes the vertex coordiantes of a box
 * \param[in] i vertex index
 * \param[in] A0 min point of bounding box
 * \param[in] A1 max point of bounding box
 * \param[out] P vertex coordinates
 */
void vertexOfBox(int  const &i, array3D const &A0, array3D const &A1, array3D &P)
{

    switch(i){

        case 0:
            P[0] = A0[0];
            P[1] = A0[1];
            P[2] = A0[2];
            break;

        case 1:
            P[0] = A1[0];
            P[1] = A0[1];
            P[2] = A0[2];
            break;

        case 2:
            P[0] = A0[0];
            P[1] = A1[1];
            P[2] = A0[2];
            break;

        case 3:
            P[0] = A1[0];
            P[1] = A1[1];
            P[2] = A0[2];
            break;

        case 4:
            P[0] = A0[0];
            P[1] = A0[1];
            P[2] = A1[2];
            break;

        case 5:
            P[0] = A1[0];
            P[1] = A0[1];
            P[2] = A1[2];
            break;

        case 6:
            P[0] = A0[0];
            P[1] = A1[1];
            P[2] = A1[2];
            break;

        case 7:
            P[0] = A1[0];
            P[1] = A1[1];
            P[2] = A1[2];
            break;

        default:
            assert(false);
            break;

    }

    return;
}

/*!
 * rotates a vector in 3D using Rodrigues' formula.
 * \param[in] vector vector to be rotated
 * \param[in] axis rotation axis
 * \param[in] theta rotation angle
 * \return rotated vector
 */
array3D rotateVector( array3D const &vector, array3D const &axis, double theta)
{

    array3D rotated;
    double cosTheta = cos(theta);

    rotated  = cosTheta * vector;
    rotated += sin(theta) * crossProduct(axis,vector);
    rotated += (1.0 - cosTheta) * dotProduct(axis,vector) * axis;

    return rotated; 
}

/*!
 * computes the area of an triangle
 * \param[in] a first vertex coordinates
 * \param[in] b second vertex coordinates
 * \param[in] c third vertex coordinates
 */
double areaTriangle( array3D const &a, array3D const &b, array3D const &c)
{
    return 0.5 *norm2(crossProduct(b-a,c-a));
}

/*
 * Gets the number of edges of a polygon
 * \return number of edges
 */
int polygonEdgesCount( std::vector<array3D> const &V)
{
    return V.size();
}

/*
 * Gets the number of subtriangles of a polygon
 * \return number of subtriangles
 */
int polygonSubtriangleCount( std::vector<array3D> const &V)
{
    return V.size()-2;
}

/*
 * Gets the edge coordinates of a convex polygon
 * \param[in] edge index
 * \param[in] V polgon vertices
 * \param[in] V0 first vertice coordinates of edge
 * \param[in] V1 second vertice coordinates of edge
 */
void edgeOfPolygon( int const &edge, std::vector<array3D> const &V, array3D &V0, array3D &V1)
{
    assert(edge<polygonEdgesCount(V));

    V0 = V[edge];
    V1 = V[(edge+1) %V.size()];
    return;
}

/*
 * Gets the subtriangle vertices' coordinates of a convex polygon
 * \param[in] triangle index of triangle
 * \param[in] V polgon vertices
 * \param[in] V0 first vertice coordinates of triangle
 * \param[in] V1 second vertice coordinates of triangle
 * \param[in] V2 third vertice coordinates of triangle
 */
void subtriangleOfPolygon( int const &triangle, std::vector<array3D> const &V, array3D &V0, array3D &V1, array3D &V2)
{
    assert(triangle<polygonSubtriangleCount(V));

    V0 = V[0];
    V1 = V[triangle+1];
    V2 = V[triangle+2];
    return;
}

/*!
 * \}
*/

}

}
