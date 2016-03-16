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

# include "CG.hpp"

# include "bitpit_operators.hpp"

# include <lapacke.h>

namespace bitpit{

/*!
    \brief Namespace for basic computational geometry functions
*/
namespace CGElem{

/*!
    \ingroup CG
    \{
*/

/*!
 * Computes distance point to line in 3D
 * @param[in] P point coordinates
 * @param[in] Q point on line
 * @param[in] n line direction
 * @param[out] xP closest point on line
 * @return distance
 */
double distancePointLine(
        std::array< double, 3 > const &P,
        std::array< double, 3 > const &Q,
        std::array< double, 3 > const &n,
        std::array< double, 3 > &xP
        ) {

    // Local variables
    double         d;

    // Project P onto the line
    xP = Q + dotProduct(P - Q, n) * n;

    // Distance between P and the line
    d = norm(P - xP, 2);

    return(d); };

/*!
 * Computes distance point to plane
 * @param[in] P point coordinates
 * @param[in] Q point on plane
 * @param[in] n plane normal
 * @param[out] xP closest point on line
 * @return distance
 */
double distancePointPlane(
        std::array< double, 3 > const &P,
        std::array< double, 3 > const &Q,
        std::array< double, 3 > const &n,
        std::array< double, 3 > &xP
        ) {

    // Local variables
    double      d;

    // Project point P onto the plane
    xP = P - dotProduct(P - Q, n) * n;

    // Distance between P and its projection onto the plane
    d = norm(P - xP,2);

    return(d); 
}

/*!
 * Computes distance point to segment in 3D
 * @param[in] P point coordinates
 * @param[in] Q1 segment starting point
 * @param[in] Q2 segment ending point
 * @param[out] xP closest point on line
 * @param[out] flag point mapping onto segment vertices (flag = 1, 2), or onto segment interior (flag = 0)
 * @return distance
 */
double distancePointSegment(
        std::array< double, 3 > const &P,
        std::array< double, 3 > const &Q1,
        std::array< double, 3 > const &Q2,
        std::array< double, 3 >       &xP,
        int                      &flag
        ) {

    // Local variables
    double             t, s, d;
    std::array< double, 3 > n;

    // Project P onto the line supporting the segment
    t = norm(Q2 - Q1, 2);
    n = (Q2 - Q1)/t;
    d = distancePointLine(P, Q1, n, xP);

    // Restrict projection onto the segment
    s = dotProduct(xP - Q1, n)/t;
    s = std::min(std::max(s, 0.0), 1.0);
    xP = Q1 + s*t*n;

    // Compute the distance between the point and the segment
    d = norm(P - xP, 2);

    // Flag
    flag = 0;
    if (s == 0.0) {
        flag = 1;
    }
    if (s == 1.0) {
        flag = 2;
    }

    return(d); };

/*!
 * Computes distance point to triangle
 * @param[in] P point coordinates
 * @param[in] Q1 first triangle vertex
 * @param[in] Q2 second triangle vertex
 * @param[in] Q3 third triangle vertex
 * @param[out] xP closest point on triangle
 * @param[out] flag point projecting onto triangle's interior (flag = 0), triangle's vertices (flag = 1, 2, 3) or triangle's edges (flag = -1, -2, -3)
 * @return distance
 */
double distancePointTriangle(
        std::array< double, 3 > const &P,
        std::array< double, 3 > const &Q1,
        std::array< double, 3 > const &Q2,
        std::array< double, 3 > const &Q3,
        std::array< double, 3 >       &xP,
        int                      &flag
        ) {

    std::array< std::array< double, 3 >, 3 > r = { Q1, Q2, Q3} ;

    double              d;
    std::array< double, 3 >  n;
    int                 k, i, j, vertex0, vertex1 ;

    int                 count, oneNegative, flagLocal ;
    std::array<int,2>        twoNegative ;

    double              *A = new double [6]  ; 
    double              *lambda = new double [ 3 ] ;


    // Project P onto the plane supporting the triangle
    n = crossProduct(Q2 - Q1, Q3 - Q1);
    n = n/norm(n, 2);
    d = distancePointPlane(P,Q1,n,xP);

    k = 0 ;
    for( i=0; i<2; ++i){ //columns
        for( j=0; j<3; ++j){ //rows

            A[k] =  r[i][j] - r[2][j] ;
            ++k ;

        };
    };


    d  = distancePointPlane(P,Q1,n,xP );

    lambda[0] = xP[0] -r[2][0] ;
    lambda[1] = xP[1] -r[2][1];
    lambda[2] = xP[2] -r[2][2];


    LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', 3, 2, 1, A, 3, lambda, 3 ) ;


    lambda[2] = 1. -lambda[0] -lambda[1] ;

    count = 0;
    twoNegative.fill(0.) ;
    for( i=0; i<3; ++i){
        if( lambda[i] < 0){
            oneNegative = i;
            twoNegative[count] = i ;
            ++count ;
        }
    };

    if( count == 0){
        flag = 0 ;
    }

    else{
        if( count == 1){
            vertex0 = (oneNegative +1) %3 ;
            vertex1 = (vertex0     +1) %3 ;
            d       =  distancePointSegment(P,r[vertex0],r[vertex1],xP, flagLocal)  ;

            if( flagLocal ==0 ){
                flag    = -(vertex0+1) ;
            }

            else if( flagLocal == 1){
                flag = vertex0+1 ; 
            }
            
            else if( flagLocal == 2){
                flag = vertex1+1 ; 
            };
            
        }

        else{
            vertex0 = 3 -twoNegative[0] -twoNegative[1] ;
            flag    = (vertex0 +1) ;
            xP      = r[vertex0] ;
            d       = norm2( P - xP)  ;
        };

    };

    delete[] A ;
    delete[] lambda ;

    return d;



};

/*!
 * Computes distances of point cloud to triangle
 * @param[in] P point cloud coordinates
 * @param[in] Q1 first triangle vertex
 * @param[in] Q2 second triangle vertex
 * @param[in] Q3 third triangle vertex
 * @param[out] xP closest points on triangle
 * @param[out] flag point projecting onto triangle's interior (flag = 0), triangle's vertices (flag = 1, 2, 3) or triangle's edges (flag = -1, -2, -3)
 * @return distance
 */
std::vector<double> distanceCloudTriangle(
        std::vector< std::array< double, 3 > >    const &P,
        std::array< double, 3 >              const &Q1,
        std::array< double, 3 >              const &Q2,
        std::array< double, 3 >              const &Q3,
        std::vector< std::array< double, 3 > >    &xP,
        std::vector< int >                   &flag
        ) {

    int                         N( P.size() ) ;

    std::array< std::array<double,3>, 3 > r = { Q1, Q2, Q3} ;
    double                      *A = new double [6]  ; 
    double                      *ptrlambda = new double [ 3*N ] ;

    // Local variables
    std::vector<double>      d(N);
    std::array< double, 3 >  n;
    int                 i, j, k, vertex0, vertex1;

    std::array< double, 3 >  lambda ;
    int                 count, oneNegative, flagLocal ;
    std::array<int,2>        twoNegative ;


    xP.resize(N) ;
    flag.resize(N) ;

    // Project P onto the plane supporting the triangle
    n = crossProduct(Q2 - Q1, Q3 - Q1);
    n = n/norm(n, 2);

    k = 0 ;
    for( i=0; i<2; ++i){ //columns
        for( j=0; j<3; ++j){ //rows

            A[k] =  r[i][j] - r[2][j] ;
            ++k ;

        };
    };


    i=0 ;
    k=0 ;
    for( const auto & p :P){

        d[i]  = distancePointPlane( p, Q1, n, xP[i] );

        ptrlambda[k  ] = xP[i][0] -r[2][0] ;
        ptrlambda[k+1] = xP[i][1] -r[2][1];
        ptrlambda[k+2] = xP[i][2] -r[2][2];

        k += 3 ;
        i++ ;
    };

    LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', 3, 2, N, A, 3, ptrlambda, 3 ) ;

    j =0 ;
    for( k=0; k<N; ++k) {

        lambda[0] = ptrlambda[j] ;
        lambda[1] = ptrlambda[j+1] ;
        lambda[2] = 1. -lambda[0] -lambda[1] ;

        j += 3 ;

        i = 0 ;
        count = 0;
        twoNegative.fill(0.) ;

        for( const auto& value : lambda){
            if( value < 0){
                oneNegative = i;
                twoNegative[count] = i ;
                ++count ;
            }

            ++i;
        };

        if( count == 0){
            flag[k] = 0 ;
        }

        else{
            if( count == 1){
                vertex0 = (oneNegative +1) %3 ;
                vertex1 = (vertex0     +1) %3 ;

                d[k]       =  distancePointSegment( P[k],r[vertex0],r[vertex1],xP[k], flagLocal)  ;

                if( flagLocal ==0 ){
                    flag[k]    = -(vertex0+1) ;
                }

                else if( flagLocal == 1){
                    flag[k] = vertex0+1 ; 
                }

                else if( flagLocal == 2){
                    flag[k] = vertex1+1 ; 
                };
            }

            else{
                vertex0 = 3 -twoNegative[0] -twoNegative[1] ;
                flag[k] = (vertex0 +1) ;
                xP[k]   = r[vertex0] ;
                d[k]    =  norm2( P[k] - r[vertex0] )  ;
            };

        };

    };

    delete[] A ;
    delete[] ptrlambda ;

    return d ;

};

/*!
 * Computes distances of point to generic simplex
 * @param[in] P point coordinates
 * @param[in] V simplex vertices coordinates
 * @param[out] xP closest points on simplex
 * @param[out] flag point projecting onto simplex's interior (flag = 0), simplex's vertices (flag = 1, 2, ...) or triangle's edges (flag = -1, -2, -...)
 * @return distance
 */
double distancePointSimplex(
        std::array< double, 3 >            const &P,
        std::vector< std::array< double, 3 > >  const &V,
        std::array< double, 3 >                  &xP,
        int                                 &flag
        ) {
    // Local variables
    int                     local, n( V.size() );
    double                  d = 1.0e+18, dT;
    std::array< double, 3>       xT;

    // Counters
    int                     i, j, p, m;

    if (n == 2) {

        // Segment -------------------------------------------------------------- //
        d = distancePointSegment(P, V[0], V[1], xP, flag);

    }
    if (n == 3) {

        // Triangle ------------------------------------------------------------- //
        d = distancePointTriangle(P, V[0], V[1], V[2], xP, flag);

    }
    else {

        // Generic convex polygon ----------------------------------------------- //

        // Compute the distance from each triangle in the simplex
        p = n - 2;
        m = 0;
        j = 1;
        while (m < p) {
            i = j;
            j = i+1;
            dT = distancePointTriangle(P, V[0], V[i], V[j], xT, local);
            if (dT <= d) {
                d = dT;
                xP = xT;
                switch (local) {
                    case -1 : { if (m == 0) {flag = -1;} break; }
                    case -2 : { flag = -i; break; }
                    case -3 : { if (m == p-1) {flag = -j;} break; }
                    case 0  : { flag =  0; break; }
                    case 1  : { flag =  1; break; }
                    case 2  : { flag =  i; break; }
                    case 3  : { flag =  j; break; }
                }
            }
            m++;
        } //next i

    }

    return(d); };

/*!
 * Computes distances of point cloud to generic simplex
 * @param[in] P point cloud coordinates
 * @param[in] V simplex vertices coordinates
 * @param[out] xP closest points on simplex
 * @param[out] flag point projecting onto simplex's interior (flag = 0), simplex's vertices (flag = 1, 2, ...) or triangle's edges (flag = -1, -2, -...)
 * @return distance
 */
std::vector<double> distanceCloudSimplex(
        std::vector< std::array< double, 3 > >  const &P,
        std::vector< std::array< double, 3 > >  const &V,
        std::vector< std::array< double, 3 > >        &xP,
        std::vector< int >                       &flag
        ) {


    // Local variables
    int                     N( P.size() ), n( V.size() );
    std::vector<double>          d(N,1.0e+18) ;

    // Counters

    if (n == 2) { //Segment

        int     i ;

        xP.resize(N) ;
        flag.resize(N) ;

        i= 0 ;
        for( const auto & point :P ){
            d[i] = distancePointSegment(point, V[0], V[1], xP[i], flag[i]);
            ++i ;
        }

    }

    else if (n == 3) {  // Triangle ------------------------------------------------------------- //

        d = distanceCloudTriangle(P, V[0], V[1], V[2], xP, flag);

    }

    else { // Generic convex polygon ----------------------------------------------- //

        int                             i, j, p, m, t;

        std::vector< int >                   local;
        std::vector< double >                dT;
        std::vector< std::array< double, 3> >     xT;

        p = n - 2;
        m = 0;
        j = 1;

        while (m < p) { // foreach triangle
            i = j;
            j = i+1;
            dT = distanceCloudTriangle(P, V[0], V[i], V[j], xT, local);

            for( t=0; t<N; ++t ){ //foreach point of cloud

                if (dT[t] <= d[t]) {
                    d[t] = dT[t];
                    xP[t] = xT[t];
                    switch (local[t]) {
                        case -1 : { if (m == 0) {flag[t] = -1;} break; }
                        case -2 : { flag[t] = -i; break; }
                        case -3 : { if (m == p-1) {flag[t] = -j;} break; }
                        case 0  : { flag[t] =  0; break; }
                        case 1  : { flag[t] =  1; break; }
                        case 2  : { flag[t] =  i; break; }
                        case 3  : { flag[t] =  j; break; }
                    }
                }
            } //foreach point in cloud

            m++;

        } // end triangles
    } //end generic polygon

    return(d); 

};

/*!
 * Computes intersection between two lines in 3D
 * @param[in] n1 direction of first line
 * @param[in] P1 point on first line
 * @param[in] n2 direction of second line
 * @param[in] P2 point on second line
 * @param[out] xP intersection point if intersect, else unaltered
 * @return if intersect
 */
bool intersectLineLine(
        std::array<double, 3> const &n1,
        std::array<double, 3> const &P1,
        std::array<double, 3> const &n2,
        std::array<double, 3> const &P2,
        std::array<double, 3>       &P
        ) {

    // Local variables
    double                  xi, delta, den;

    delta = dotProduct(n1, n2);
    if (abs((abs(delta) - 1.0)) < 1.0e-12) {
        return false ;
    }

    else {
        den = (1.0 - pow(delta, 2));
        xi = (dotProduct(P2, n1)
                - dotProduct(P1, n1)
                + delta * (dotProduct(P1, n2) - dotProduct(P2, n2))) / den;
        P = P1 + xi * n1;
    }

    return(true); 
};

/*!
 * Computes intersection between two segments in 3D
 * @param[in] P1 star point of first segment
 * @param[in] P2 end point of first segment
 * @param[in] Q1 star point of second segment
 * @param[in] Q2 end point of second segment
 * @param[out] x intersection point if intersect, else unaltered
 * @return if intersect
 */
bool intersectSegmentSegment(
        std::array<double, 3> const &P1,
        std::array<double, 3> const &P2,
        std::array<double, 3> const &Q1,
        std::array<double, 3> const &Q2,
        std::array<double, 3>       &x
        ) {

    // Parameters
    double const                    abs_tol = 1.0e-14;

    // Local variables
    double                          lP, lQ, lxP, lxQ;
    std::array<double, 3>                nP, nQ;

    // Counters
    // none

    // ========================================================================== //
    // COMPUTE INTERSECTION POINT BETWEEN SUPPORTING LINES                        //
    // ========================================================================== //
    nP = P2 - P1;
    lP = norm2(nP);
    nP = nP/lP;
    nQ = Q2 - Q1;
    lQ = norm2(nQ);
    nQ = nQ/lQ;

    // ========================================================================== //
    // CHECK IF INTERSECTION POINTS LIES INSIDE EACH SEGMENT                      //
    // ========================================================================== //
    if ( intersectLineLine(nP, P1, nQ, Q1,x) ) {
        lxP = dotProduct(x - P1, nP);
        lxQ = dotProduct(x - Q1, nQ);

        if ( lxP <= lP + abs_tol && lxP >= -abs_tol && lxQ  <= lQ + abs_tol && lxQ  >= -abs_tol )  {

            return(true);
        }

    }

    return (false);


};

/*!
 * Computes intersection between line and plane
 * @param[in] P1 point on line
 * @param[in] n1 direction of line
 * @param[in] P2 point on plane
 * @param[in] n2 normal to plane
 * @param[out] P intersection point if intersect, else unaltered
 * @return if intersect
 */
bool intersectLinePlane(
        std::array<double, 3> const &P1,
        std::array<double, 3> const &n1,
        std::array<double, 3> const &P2,
        std::array<double, 3> const &n2,
        std::array<double, 3>       &P
        ) {

    // Parameters
    double const                tol = 1.0e-14;

    // Local variables
    double                      s, xi;

    // ========================================================================== //
    // CHECK DEGENERATE CASES                                                     //
    // ========================================================================== //
    s = dotProduct(n1, n2);
    if (abs(s) < tol) { return(false); }

    // ========================================================================== //
    // FIND INTERSECTION POINTS                                                   //
    // ========================================================================== //
    xi = -dotProduct(P1 - P2, n2)/s;
    P = P1 + xi * n1;

    return(true); 
};

/*!
 * Computes intersection between segment and plane
 * @param[in] Q1 start point of segment
 * @param[in] Q2 end point of segment
 * @param[in] P2 point on plane
 * @param[in] n2 normal to plane
 * @param[out] P intersection point if intersect, else unaltered
 * @return if intersect
 */
bool intersectSegmentPlane(
        std::array<double, 3> const &Q1,
        std::array<double, 3> const &Q2,
        std::array<double, 3> const &P2,
        std::array<double, 3> const &n2,
        std::array<double, 3>       &P
        ) {

    std::array<double, 3>        n;


    // ========================================================================== //
    // FIND INTERSECTION BETWEEN SEGMENT'S SUPPORTING LINE AND PLANE              //
    // ========================================================================== //
    n = Q2 - Q1;
    n = n/norm2(n);
    if (   intersectLinePlane(Q1, n, P2, n2,P) ) {
        if (intersectPointSegment(P, Q1, Q2)) {    
            return(true);
        }
    }

    return(false); };

/*!
 * Computes intersection between two planes
 * @param[in] P1 point on first plane
 * @param[in] n1 normal to first plane
 * @param[in] P2 point on second plane
 * @param[in] n2 normal to second plane
 * @param[out] Pl point on intersection line
 * @param[out] nl direction of intersection line
 * @return if intersect
 */
bool intersectPlanePlane(
        std::array<double, 3> const &P1,
        std::array<double, 3> const &n1,
        std::array<double, 3> const &P2,
        std::array<double, 3> const &n2,
        std::array<double, 3>       &Pl,
        std::array<double, 3>       &nl
        ) {

    // Parameters
    double const                    tol = 1.0e-14;

    // Local variables
    double                          s;
    std::array< double, 3>               v;

    // ========================================================================== //
    // CHECK DEGENERATE CONDITIONS                                                //
    // ========================================================================== //
    if (abs(dotProduct(n1, n2) - 1.0) <= tol) { return(false); }

    // ========================================================================== //
    // FIND INTERSECTION LINE                                                     //
    // ========================================================================== //

    // Line direction
    nl = crossProduct(n1, n2);
    nl = nl/norm2(nl);

    // Point onto the line
    s = dotProduct(P1, n1) - dotProduct(P2, n2);
    v = n1 - n2;
    if (abs(v[2]) > tol) {
        Pl[0] = P1[0];
        Pl[1] = P1[1];
        Pl[2] = (s - Pl[0]*v[0] - Pl[1]*v[1])/v[2];
    }

    else if( abs(v[1]) > tol) {

        Pl[0] = P1[0];
        Pl[2] = P1[2];
        Pl[1] = (s - Pl[0]*v[0] - Pl[2]*v[2])/v[1];
    }

    else {
        Pl[1] = P1[1];
        Pl[2] = P1[2];
        Pl[0] = (s - Pl[1]*v[1] - Pl[2]*v[2])/v[0];
    }

    return(true); 

};

/*!
 * Computes intersection between triangle and a line
 * @param[in] P point on plane
 * @param[in] n normal to plane
 * @param[in] A first vertex of triangle
 * @param[in] B second vertex of triangle
 * @param[in] C third vertex of triangle
 * @param[out] Q intersection point
 * @return if intersect
 */
bool intersectLineTriangle(
        std::array<double, 3> const &P,
        std::array<double, 3> const &n,
        std::array<double, 3> const &A,
        std::array<double, 3> const &B,
        std::array<double, 3> const &C,
        std::array<double, 3>       &Q
        ) {

    std::array<double, 3>        nT;

    // Compute triangle normal
    nT = crossProduct(B - A, C - A);
    nT = nT/norm2(nT);

    // Compute intersection with line and triangle's supporting plane
    if ( !intersectLinePlane(P, n, A, nT, Q) ) { return(false); };

    if (!intersectPointTriangle(Q, A, B, C)) { return(false); }

    return(true); 
};

/*!
 * Computes intersection between triangle and a segment
 * @param[in] P0 start point of segment
 * @param[in] P1 end point of segment
 * @param[in] A first vertex of triangle
 * @param[in] B second vertex of triangle
 * @param[in] C third vertex of triangle
 * @param[out] Q intersection point
 * @return if intersect
 */
bool intersectSegmentTriangle(
        std::array<double, 3> const &P0,
        std::array<double, 3> const &P1,
        std::array<double, 3> const &A,
        std::array<double, 3> const &B,
        std::array<double, 3> const &C,
        std::array<double, 3>       &Q
        ) {

    // Local variables
    std::array<double, 3>        n;

    n = P1 - P0;
    n = n/norm2(n);

    if ( !intersectLineTriangle(P0, n, A, B, C, Q) ) { return(false); }

    if ( !intersectPointSegment(Q, P0, P1) ) { return(false); }

    return(true); 

};

/*!
 * Computes intersection between triangle and a segment
 * @param[in] P point on line
 * @param[in] n direction of line
 * @param[in] V simplex vertices coordinates
 * @param[out] Q intersection point
 * @return if intersect
 */
bool intersectLineSimplex(
        std::array<double, 3> const          &P,
        std::array<double, 3> const          &n,
        std::vector<std::array<double, 3> > const &V,
        std::array<double, 3>                &Q
        ) {

    // Local variables
    int                     m;
    std::array<double, 3>        A, B, C;

    // Counters
    int                     i;

    m = V.size();

    i = 1;
    for( i=1; i< m-1; ++i){

        // Split simplex into triangles
        A = V[0];
        B = V[i];
        C = V[(i+1)%m];

        // Compute intersection
        if( intersectLineTriangle(A, B, C, P, n, Q) ) { return(true); }


    } //next i

    return(false); 
};

/*!
 * Computes intersection between triangle and a segment
 * @param[in] P0 start point of segment
 * @param[in] P1 end point of segment
 * @param[in] V simplex vertices coordinates
 * @param[out] Q intersection point
 * @return if intersect
 */
bool intersectSegmentSimplex(
        std::array<double, 3> const          &P0,
        std::array<double, 3> const          &P1,
        std::vector<std::array<double, 3> > const &V,
        std::array<double, 3>                &Q 
        ) {


    // Local variables
    int                     m(V.size());
    std::array<double, 3>        A, B, C;
    int                     i;


    for( i=1; i<m-1; ++i){

        // Split simplex into triangles
        A = V[0];
        B = V[i];
        C = V[(i+1)%m];

        // Compute intersection
        if( intersectSegmentTriangle(A, B, C, P0, P1, Q) ) { return(true); };

    } //next i

    return(false); };

/*!
 * Computes intersection between two axis aligned bounding boxes
 * @param[in] A1 min point of first box
 * @param[in] A2 max point of first box
 * @param[in] B1 min point of second box
 * @param[in] B2 max point of second box
 * @param[in] dim number of dimension to be checked
 * @return if intersect
 */
bool intersectBoxBox(
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        std::array<double, 3> const              &B1,
        std::array<double, 3> const              &B2,
        int                                 dim
        ){

    int     d;

    for( d=0; d<dim; ++d){

        if( B1[d] > A2[d] || B2[d] < A1[d] ){
            return false;
        };

    };

    return true;

};

/*!
 * Computes intersection between two axis aligned bounding boxes
 * @param[in] A1 min point of first box
 * @param[in] A2 max point of first box
 * @param[in] B1 min point of second box
 * @param[in] B2 max point of second box
 * @param[out] I1 min point of intersection box
 * @param[out] I2 max point of intersection box
 * @param[in] dim number of dimension to be checked
 * @return if intersect
 */
bool intersectBoxBox(
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        std::array<double, 3> const              &B1,
        std::array<double, 3> const              &B2,
        std::array<double, 3>                    &I1,
        std::array<double, 3>                    &I2,
        int                                 dim
        ){

    int     d;

    for( d=0; d<dim; ++d){

        if( B1[d] > A2[d] || B2[d] < A1[d] ){
            return false;
        }

        else{
            I1[d] = std::max( A1[d], B1[d] ) ;
            I2[d] = std::min( A2[d], B2[d] ) ;
        };

    };

    return true;

};

/*!
 * Computes intersection between an axis aligned bounding box and a triangle
 * @param[in] A1 min point of first box
 * @param[in] A2 max point of first box
 * @param[in] V1 first vertex of triangle
 * @param[in] V2 first vertex of triangle
 * @param[in] V3 first vertex of triangle
 * @return if intersect
 */
bool intersectBoxTriangle(
        std::array<double, 3> const              &A1, 
        std::array<double, 3> const              &A2,
        std::array<double, 3> const              &V1,
        std::array<double, 3> const              &V2,
        std::array<double, 3> const              &V3
        ){

    std::array<double,3>     B1, B2, E1, E2, P ;

    P.fill(0.0) ;

    //Check if Triangle Boundig Box and Box overlap -> necessary condition
    computeAABBTriangle( V1, V2, V3, B1, B2) ;

    if( !intersectBoxBox( A1, A2, B1, B2) ) { return(false); }

    //Check if Triangle Box edges intersection 
    for( int i=0; i<12; ++i){
        edgeOfBox( i, A1, A2, E1, E2) ;
        if( intersectSegmentTriangle(E1,E2,V1,V2,V3,P)) { return(true); }
    };

    return (false);

};


/*!
 * Computes intersection between an axis aligned bounding box and a triangle
 * @param[in] A1 min point of first box
 * @param[in] A2 max point of first box
 * @param[in] V1 first vertex of triangle
 * @param[in] V2 first vertex of triangle
 * @param[in] V3 first vertex of triangle
 * @param[out] P intersection points between triangle and box edges
 * @return if intersect
 */
bool intersectBoxTriangle(
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        std::array<double, 3> const              &V1,
        std::array<double, 3> const              &V2,
        std::array<double, 3> const              &V3,
        std::vector< std::array<double, 3> >          &P
        ){

    bool                intersect(false) ;
    std::array<double,3>     B1, B2, E1, E2, p ;

    P.clear() ;
    P.reserve(3) ;

    //Check if Triangle Boundig Box and Box overlap -> necessary condition
    computeAABBTriangle( V1, V2, V3, B1, B2) ;

    if( !intersectBoxBox( A1, A2, B1, B2) ) { return(false); }

    //Check if Triangle Box edges intersection 
    for( int i=0; i<12; ++i){
        edgeOfBox( i, A1, A2, E1, E2) ;
        if( intersectSegmentTriangle(E1,E2,V1,V2,V3,p)) {
            P.push_back( p );
            intersect = true ;
        } ;
    };

    return intersect;

};

/*!
 * Computes intersection between an axis aligned bounding box and a segment
 * @param[in] V1 start point of segment
 * @param[in] V2 end point of segment
 * @param[in] A1 min point of box
 * @param[in] A2 max point of box
 * @param[in] dim number of dimension to be checked
 * @return if intersect
 */
bool intersectSegmentBox(
        std::array<double, 3> const              &V1,
        std::array<double, 3> const              &V2,
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        int                                 dim
        ){

    int                 i ;
    std::array<double,3>     B1, B2, p ;

    //Check if Triangle Boundig Box and Box overlap -> necessary condition
    computeAABBSegment( V1, V2, B1, B2) ;
    if( !intersectBoxBox( A1, A2, B1, B2, dim) ) { return(false); }

    //Check if Triangle Box edges intersection 
    if( dim == 2){

        std::array<double,3> E1, E2 ;

        for( i=0; i<4; ++i){
            edgeOfBox( i, A1, A2, E1, E2) ;

            if( intersectSegmentSegment(E1,E2,V1,V2,p)) { return(true); } ;

        };

    }

    { //3D

        std::vector< std::array<double,3> > E ;
        E.resize(4) ;

        for( i=0; i<12; ++i){
            faceOfBox( i, A1, A2, E[0], E[1], E[2], E[3]) ;

            if( intersectSegmentSimplex(V1,V2,E,p) ) { return(true); } ;

        };

    };

    return (false);

};

/*!
 * Computes intersection between an axis aligned bounding box and a segment
 * @param[in] V1 start point of segment
 * @param[in] V2 end point of segment
 * @param[in] A1 min point of box
 * @param[in] A2 max point of box
 * @param[out] P intersection point segment and box faces
 * @param[in] dim number of dimension to be checked
 * @return if intersect
 */
bool intersectSegmentBox(
        std::array<double, 3> const              &V1,
        std::array<double, 3> const              &V2,
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        std::vector< std::array<double, 3> >          &P,
        int                                 dim
        ){

    bool                intersect(false) ;
    int                 i ;
    std::array<double,3>     p, B1, B2 ;

    P.clear() ;

    //Check if segment Boundig Box and Box overlap -> necessary condition
    computeAABBSegment( V1, V2, B1, B2) ;
    if( !intersectBoxBox( A1, A2, B1, B2, dim) ) { return(false); }

    P.reserve(2) ;

    if( dim == 2){ //check if box edge and segment intersect


        std::array<double,3> E1, E2 ;

        for( i=0; i<4; ++i){
            edgeOfBox( i, A1, A2, E1, E2) ;

            if( intersectSegmentSegment(E1,E2,V1,V2,p)) {
                intersect = true ;
                P.push_back(p) ;
            } ;

        };

        return( intersect ) ;

    }

    else { //3D check if box face and segment intersect


        std::vector< std::array<double,3> > E ;
        E.resize(4) ;

        for( i=0; i<12; ++i){
            faceOfBox( i, A1, A2, E[0], E[1], E[2], E[3]) ;

            if( intersectSegmentSimplex(V1,V2,E,p) ) {
                intersect = true ;
                P.push_back(p) ;
            } ;

        };

        return intersect ;

    };


};


/*!
 * Computes intersection between an axis aligned bounding box and a simplex
 * @param[in] A1 min point of first box
 * @param[in] A2 max point of first box
 * @param[in] VS simplex vertices coordinates
 * @param[in] dim number of dimension to be checked
 * @return if intersect
 */
bool intersectBoxSimplex(
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        std::vector< std::array<double, 3> > const    &VS,
        int                                 dim
        ){

    int                 i, j, n, m, p ;
    std::array<double,3>     B1, B2 ;

    //Check if Triangle Boundig Box and Box overlap -> necessary condition
    computeAABBSimplex( VS, B1, B2) ;

    if( !intersectBoxBox( A1, A2, B1, B2, dim) ) { return(false); }

    n = VS.size() ;

    if( n == 2){
        return( intersectSegmentBox( VS[0], VS[1], A1, A2, dim ) ) ;
    }

    else if( n == 3){
        return(  intersectBoxTriangle( A1, A2, VS[0], VS[1], VS[2] ) ) ;
    }

    else{

        // Generic convex polygon ----------------------------------------------- //

        // Compute the distance from each triangle in the simplex
        p = n - 2;
        m = 0;
        j = 1;
        while (m < p) {
            i = j;
            j = i+1;

            if( intersectBoxTriangle( A1, A2, VS[0], VS[i], VS[j] ) ){
                return (true);
            };

            m++;
        } //next i

    };

    return (false);

};


/*!
 * Computes intersection between an axis aligned bounding box and a simplex
 * @param[in] A1 min point of first box
 * @param[in] A2 max point of first box
 * @param[in] VS simplex vertices coordinates
 * @param[out] P intersection points simplex box edges
 * @param[in] dim number of dimension to be checked
 * @return if intersect
 */
bool intersectBoxSimplex(
        std::array<double, 3> const              &A1,
        std::array<double, 3> const              &A2,
        std::vector< std::array<double, 3> > const    &VS,
        std::vector< std::array<double, 3> >          &P,
        int                                 dim
        ){

    int                         n  ;

    P.clear() ;

    { //Check if Triangle Boundig Box and Box overlap -> necessary condition

        std::array<double,3>             B1, B2 ;
        computeAABBSimplex( VS, B1, B2) ;

        if( !intersectBoxBox( A1, A2, B1, B2, dim) ) { return(false); }
    }


    n = VS.size() ;
    if( n == 2){ //segment
        return( intersectSegmentBox( VS[0], VS[1], A1, A2, P, dim ) ) ;
    }

    else if( n == 3){ //triangle
        return(  intersectBoxTriangle( VS[0], VS[1], VS[2], A1, A2, P ) ) ;
    }

    else{  //generic convex polygon split into triangles

        bool                        intersect(false) ;
        int                         i, j, m, p ;
        std::vector< std::array<double, 3> >  partial ;

        p = n - 2;
        m = 0;
        j = 1;
        while (m < p) {
            i = j;
            j = i+1;

            if( intersectBoxTriangle( A1, A2, VS[0], VS[i], VS[j], partial ) ){
                intersect = true ;
                P.insert( P.begin(), partial.begin(), partial.end() ) ;
            };

            m++;
        } //next i

        return intersect ;

    };


};



//to levelset // -------------------------------------------------------------------------- //
//to levelset bool intersectLineSurface(
//to levelset         std::array<double,3>  const  &x1,
//to levelset         std::array<double,3>  const  &n1,
//to levelset         std::array<double,3>  const  &x2,
//to levelset         std::array<double,3>  const  &n2,
//to levelset         std::array<double,3>  const  &xL,
//to levelset         std::array<double,3>  const  &nL,
//to levelset         std::array<double,3>         &xp,
//to levelset         std::array<double,3>         &np
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
//to levelset     double                  w1(0), w2(0), w(0) ;
//to levelset     std::array<double,3>         c1, c2 ;
//to levelset     bool                    s1, s2 ;
//to levelset 
//to levelset     s1      =   intersectLinePlane( xL, nL, x1, n1, c1) ;
//to levelset     s2      =   intersectLinePlane( xL, nL, x2, n2, c2) ;
//to levelset 
//to levelset     if( s1 && s2) {
//to levelset         w1      =   norm2( c2-x2) ;
//to levelset         w2      =   norm2( c1-x1) ;
//to levelset         w       =   w1+w2 ;
//to levelset 
//to levelset         w1      =   w1 /w ;
//to levelset         w2      =   w2 /w ;
//to levelset 
//to levelset         xp      =   w1*c1 +w2*c2 ;
//to levelset         np      =   w1*n1 +w2*n2 ;
//to levelset 
//to levelset         np      =   np /norm2(np) ;
//to levelset 
//to levelset         return (true) ;
//to levelset     }
//to levelset 
//to levelset     else if( s1 ){
//to levelset         xp      =   c1 ;
//to levelset         np      =   n1 ;
//to levelset 
//to levelset         return(true) ;
//to levelset     }
//to levelset 
//to levelset     else if( s2 ){
//to levelset         xp      =   c2 ;
//to levelset         np      =   n2 ;
//to levelset 
//to levelset         return(true) ;
//to levelset     }
//to levelset 
//to levelset     return (false); 
//to levelset 
//to levelset };

/*!
 * checks if points lies on segment in 3D
 * @param[in] P point coordinates
 * @param[in] P1 start point of segment
 * @param[in] P2 end point of segment
 * @return if point lies on segment
 */
bool intersectPointSegment(
        std::array< double, 3> const &P,
        std::array< double, 3> const &P1,
        std::array< double, 3> const &P2
        ) {

    double const        tol = 1.0e-14;

    bool                check;
    std::array<double, 3>    n1, n2;
    double              d1, d2;

    check = (norm2(P - P2) <= tol);
    if (!check) {
        n1 = P1 - P2;
        d1 = norm2(n1);
        n1 = n1/d1;
        n2 = P - P2;
        d2 = norm2(n2);
        n2 = n2/d2;
        check = ((dotProduct(n1, n2) >= 1.0 - tol) && (d2 <= d1));
    }

    return(check); };

/*!
 * checks if points lies on triangle
 * @param[in] P point coordinates
 * @param[in] A first vertex of triangle
 * @param[in] B second vertex of triangle
 * @param[in] C third vertex of triangle
 * @return if point lies on segment
 */
bool intersectPointTriangle(
        std::array<double, 3> const  &P,
        std::array<double, 3> const  &A,
        std::array<double, 3> const  &B,
        std::array<double, 3> const  &C
        ) {

    bool                    check = false;
    std::array<double, 3>        xlim, ylim, zlim;

    // ========================================================================== //
    // CHECK IF POINT LIES OUTSIDE TRIANGLE'S BOUNDING BOX                        //
    // ========================================================================== //
    xlim[0] = std::min(std::min(A[0], B[0]), C[0]);
    xlim[1] = std::max(std::max(A[0], B[0]), C[0]);
    ylim[0] = std::min(std::min(A[1], B[1]), C[1]);
    ylim[1] = std::max(std::max(A[1], B[1]), C[1]);
    zlim[0] = std::min(std::min(A[2], B[2]), C[2]);
    zlim[1] = std::max(std::max(A[2], B[2]), C[2]);

    // ========================================================================== //
    // CHECK IF POINT LIES INSIDE TRIANGLE                                        //
    // ========================================================================== //
    if (((P[0] >= xlim[0]) && (P[0] <= xlim[1]))
            && ((P[1] >= ylim[0]) && (P[1] <= ylim[1]))
            && ((P[2] >= zlim[0]) && (P[2] <= zlim[1]))) {

        // Scope variables
        double                  den;
        std::array<double, 2>        AA, BB, CC, PP;
        std::array<double, 3>        u, v, w, x, y, z;
        std::array<double, 3>        xi;

        // Initialize scope variables
        u = B - A;
        v = C - A;
        w = P - A;

        // Compute local ref. frame
        x = u/norm2(u);
        z = crossProduct(x, v);
        z = z/norm2(z);
        y = crossProduct(z, x);

        // Map triangle's vertices into local ref. frame
        AA[0] = 0.0;                   AA[1] = 0.0;
        BB[0] = dotProduct(u, x);     BB[1] = dotProduct(u, y);
        CC[0] = dotProduct(v, x);     CC[1] = dotProduct(v, y);
        PP[0] = dotProduct(w, x);     PP[1] = dotProduct(w, y);

        // Compute baricentric coordinates
        den = ((BB[1] - CC[1])*(AA[0] - CC[0]) + (CC[0] - BB[0])*(AA[1] - CC[1]));
        xi[0] = ((BB[1] - CC[1]) * (PP[0] - CC[0]) + (CC[0] - BB[0]) * (PP[1] - CC[1]))/den;
        xi[1] = ((CC[1] - AA[1]) * (PP[0] - CC[0]) + (AA[0] - CC[0]) * (PP[1] - CC[1]))/den;
        xi[2] = 1.0 - xi[0] - xi[1];

        // Point-in-triangle condition
        check = !((xi[0] < 0.0) || (xi[1] < 0.0) || (xi[2] < 0.0));

    }

    return(check); };

/*!
 * checks if points lies within axis aligned box
 * @param[in] P point coordinates
 * @param[in] B1 min coodinates of box
 * @param[in] B2 max coodinates of box
 * @param[in] dim number of dimensions to be checked
 * @return if point in box
 */
bool intersectPointBox(
        std::array<double, 3> const              &P,
        std::array<double, 3> const              &B1,
        std::array<double, 3> const              &B2,
        int                                 dim
        ){


    int     d;

    for( d=0; d<dim; ++d){

        if( P[d]< B1[d] || P[d] > B2[d] ){
            return false;
        };
    };

    return true ;
};

/*!
 * computes axis aligned boundig box of a segment
 * @param[in] A start point of segment
 * @param[in] B end point of segment
 * @param[out] P0 min point of bounding box
 * @param[out] P1 max point of bounding box
 */
void computeAABBSegment(
        std::array<double, 3> const &A,
        std::array<double, 3> const &B,
        std::array<double, 3>       &P0,
        std::array<double, 3>       &P1
        ) {

    int     i;

    P0 = A ;
    P1 = A ;

    for(i=0; i<3; ++i){
        P0[i] = std::min( P0[i], B[i] ) ;
        P1[i] = std::max( P1[i], B[i] ) ;
    };

    return;

};

/*!
 * computes axis aligned boundig box of a triangle
 * @param[in] A first vertex of triangle
 * @param[in] B second vertex of triangle
 * @param[in] C third vertex of triangle
 * @param[out] P0 min point of bounding box
 * @param[out] P1 max point of bounding box
 */
void computeAABBTriangle(
        std::array<double, 3> const &A,
        std::array<double, 3> const &B,
        std::array<double, 3> const &C,
        std::array<double, 3>       &P0,
        std::array<double, 3>       &P1
        ) {

    int     i;

    P0 = A ;
    P1 = A ;

    for(i=0; i<3; ++i){
        P0[i] = std::min( P0[i], B[i] ) ;
        P0[i] = std::min( P0[i], C[i] ) ;

        P1[i] = std::max( P1[i], B[i] ) ;
        P1[i] = std::max( P1[i], C[i] ) ;
    };

    return;

};

/*!
 * computes axis aligned boundig box of a simplex
 * @param[in] VS simplex vertices coordinates
 * @param[out] P0 min point of bounding box
 * @param[out] P1 max point of bounding box
 */
void computeAABBSimplex(
        std::vector< std::array<double, 3> > const &VS,
        std::array<double, 3>       &P0,
        std::array<double, 3>       &P1
        ) {

    int     i, j, n(VS.size());

    P0 = VS[0] ;
    P1 = VS[0] ;

    for( j=1; j<n; ++j){
        for(i=0; i<3; ++i){
            P0[i] = std::min( P0[i], VS[j][i] ) ;
            P1[i] = std::max( P1[i], VS[j][i] ) ;
        };
    };

    return;

};

/*!
 * computes the union of two axis aligned bounding boxes
 * @param[in] A0 min point of first bounding box
 * @param[in] A1 max point of first bounding box
 * @param[in] B0 min point of second bounding box
 * @param[in] B1 max point of second bounding box
 * @param[out] C0 min point of union bounding box
 * @param[out] C1 max point of union bounding box
 */
void unionAABB(
        std::array<double,3>  const              & A0,
        std::array<double,3>  const              & A1,
        std::array<double,3>  const              & B0,
        std::array<double,3>  const              & B1,
        std::array<double,3>                     & C0,
        std::array<double,3>                     & C1
        ){


    C0 = std::min( A0, B0) ;
    C1 = std::max( A1, B1) ;


    return;
};

/*!
 * computes the union of several axis aligned bounding boxes
 * @param[in] A0 min points of bounding boxes
 * @param[in] A1 max points of bounding boxes
 * @param[out] C0 min point of union bounding box
 * @param[out] C1 max point of union bounding box
 */
void unionAABB(
        std::vector<std::array<double,3> >  const     & A0,
        std::vector<std::array<double,3> >  const     & A1,
        std::array<double,3>                     & C0,
        std::array<double,3>                     & C1
        ){


    int     i, n( std::min(A0.size(), A1.size() ) ) ;

    if( n > 0 ){
        C0 =  A0[0] ;
        C1 =  A1[0] ;

        for( i=1; i<n; ++i){
            C0 = std::min( C0, A0[i] ) ;
            C1 = std::max( C1, A1[i] ) ;
        };
    };


    return;
};

/*!
 * computes the face coordiantes of a box
 * @param[in] i face index
 * @param[in] A0 min point of bounding box
 * @param[in] A1 max point of bounding box
 * @param[out] P0 first vertex of face
 * @param[out] P1 first vertex of face
 * @param[out] P2 first vertex of face
 * @param[out] P3 first vertex of face
 */
void faceOfBox(
        int              const &i,
        std::array<double, 3> const &A0,
        std::array<double, 3> const &A1,
        std::array<double, 3>       &P0,
        std::array<double, 3>       &P1,
        std::array<double, 3>       &P2,
        std::array<double, 3>       &P3
        ) {

    int     v0, v1, v2, v3;

    v0 = boxFaceVertexConnectivity[i][0] ;
    v1 = boxFaceVertexConnectivity[i][1] ;
    v2 = boxFaceVertexConnectivity[i][2] ;
    v3 = boxFaceVertexConnectivity[i][3] ;

    vertexOfBox(v0, A0, A1, P0 ) ;
    vertexOfBox(v1, A0, A1, P1 ) ;
    vertexOfBox(v2, A0, A1, P2 ) ;
    vertexOfBox(v3, A0, A1, P3 ) ;

    return;

};


/*!
 * computes the edge coordiantes of a box
 * @param[in] i edge index
 * @param[in] A0 min point of bounding box
 * @param[in] A1 max point of bounding box
 * @param[out] P0 first vertex of edge
 * @param[out] P1 first vertex of edge
 */
void edgeOfBox(
        int              const &i,
        std::array<double, 3> const &A0,
        std::array<double, 3> const &A1,
        std::array<double, 3>       &P0,
        std::array<double, 3>       &P1
        ) {

    int     v0, v1;

    v0 = boxEdgeVertexConnectivity[i][0] ;
    v1 = boxEdgeVertexConnectivity[i][1] ;

    vertexOfBox(v0, A0, A1, P0 ) ;
    vertexOfBox(v1, A0, A1, P1 ) ;

    return;

};

/*!
 * computes the vertex coordiantes of a box
 * @param[in] i vertex index
 * @param[in] A0 min point of bounding box
 * @param[in] A1 max point of bounding box
 * @param[out] P vertex coordinates
 */
void vertexOfBox(
        int              const &i,
        std::array<double, 3> const &A0,
        std::array<double, 3> const &A1,
        std::array<double, 3>       &P
        ) {

    switch(i){

        case 0:
            P[0] = A0[0] ;
            P[1] = A0[1] ;
            P[2] = A0[2] ;
            break;

        case 1:
            P[0] = A1[0] ;
            P[1] = A0[1] ;
            P[2] = A0[2] ;
            break;

        case 2:
            P[0] = A0[0] ;
            P[1] = A1[1] ;
            P[2] = A0[2] ;
            break;

        case 3:
            P[0] = A1[0] ;
            P[1] = A1[1] ;
            P[2] = A0[2] ;
            break;

        case 4:
            P[0] = A0[0] ;
            P[1] = A0[1] ;
            P[2] = A1[2] ;
            break;

        case 5:
            P[0] = A1[0] ;
            P[1] = A0[1] ;
            P[2] = A1[2] ;
            break;

        case 6:
            P[0] = A0[0] ;
            P[1] = A1[1] ;
            P[2] = A1[2] ;
            break;

        case 7:
            P[0] = A1[0] ;
            P[1] = A1[1] ;
            P[2] = A1[2] ;
            break;

    };


    return;

};

/*!
 * rotates a vector in 3D using Rodrigues' formula.
 * @param[in,out] v vector to be rotated
 * @param[in] a rotation axis
 * @param[in] theta rotation angle
 */
void rotateVector(
        std::array<double, 3>       &v,
        std::array<double, 3> const &a,
        double                  theta
        ) {

    v = cos(theta) * v
        + sin(theta) * crossProduct(a, v)
        + (1.0 - cos(theta)) * dotProduct(a, v) * a;

    return; };

/*!
    \}
*/

}

}
