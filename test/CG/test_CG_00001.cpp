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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include <cmath>
# include <array>
# include <vector>
# include <chrono>
# include <string>
# include <iostream>
#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

using namespace std;
using namespace bitpit;

/*!
* Subtest 001
*
* Testing basic computational geometry features.
*/
int subtest_001()
{
    // ========================================================================== //
    // INTERSECTION BETWEEN POINT AND BOX                                         //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A0, A1, B;
        bool                        intersect ;

        // Compute intersection (Test 1) ---------------------------------------- //
        A0.fill(0.) ;
        A1.fill(1.) ;

        B = {{0.25,0.25,0.25}} ;

        intersect = CGElem::intersectPointBox( B, A0, A1) ;

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN POINT AND  AXIS ALIGNED BOX" << endl;
        cout << "  BOX 1 " << endl;
        cout << "     A0: " << A0 << endl;
        cout << "     A1: " << A1 << endl;
        cout << "  Point " << endl;
        cout << "     B: " << B << endl;
        cout << " intersect? " << intersect << endl ;
        cout << endl;

    }

    // ========================================================================== //
    // INTERSECTION BETWEEN LINE AND LINE                                         //
    // ========================================================================== //
    {
        array<double,3> P, P0, P1, n0, n1;
        bool intersect;

        // Compute intersection (Test 1) ---------------------------------------- //
        P0 = {{ 0., 0., 0. }};
        n0 = {{ 1., 1., 0. }};
        P1 = {{ 3., 0., 0. }};
        n1 = {{-1., 2., 0. }};

        n0 /= norm2(n0);
        n1 /= norm2(n1);
        intersect = CGElem::intersectLineLine(P0, n0, P1, n1, P);

        cout << " - INTERSECTION BETWEEN LINE AND LINE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line0: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    line1: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection point is: " << P << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        P0 = {{ 0., 0., 0. }};
        n0 = {{ 1., 1., 0. }};
        P1 = {{ 3., 0., 0. }};
        n1 = {{-1., 2., 1. }};

        n0 /= norm2(n0);
        n1 /= norm2(n1);
        intersect = CGElem::intersectLineLine(P0, n0, P1, n1, P);

        cout << "  Test 2 " << endl;
        cout << "    line0: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    line1: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection point is: " << P << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN LINE AND PLANE                                        //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            P, P0, P1, n0, n1;
        bool intersect;

        // Compute intersection (Test 1) ---------------------------------------- //
        P0.fill(0.0);
        n0.fill(0.0);
        n0[2] = 1.0;
        P1.fill(0.0);
        P1[2] = 1.0;
        n1.fill(0.0);
        n1[2] = 1.0;
        intersect = CGElem::intersectLinePlane(P0, n0, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN LINE AND PLANE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection point is: " << P << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        P0.fill(0.0);
        n0.fill(0.0);
        n0[2] = 1.0;
        P1.fill(0.0);
        n1.fill(0.0);
        n1[1] = 1.0;
        intersect = CGElem::intersectLinePlane(P0, n0, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection point is: " << P << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN LINE AND TRIANGLE                                     //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A, B, C, P, Q, n;
        bool intersect;

        // Compute intersection (Test 1) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;
        P.fill(0.0);
        P[2] = 1.0;
        n.fill(0.0);
        n[2] = 1.0;
        intersect = CGElem::intersectLineTriangle(P, n, A, B, C, Q);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN LINE AND TRIANGLE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        if(intersect){
            cout << "    intersection point is: " << Q << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;
        P.fill(0.0);
        P[0] = 1.0;
        P[1] = 1.0;
        P[2] = 1.0;
        n.fill(0.0);
        n[2] = 1.0;
        intersect = CGElem::intersectLineTriangle( P, n, A, B, C, Q);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        if(intersect){
            cout << "    intersection point is: " << Q << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN LINE AND SIMPLEX                                      //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        vector< array<double, 3> >      V(4);
        array<double, 3>                P, Q, n;
        bool intersect;

        // Compute intersection (Test 1) ---------------------------------------- //
        V[0].fill(0.0);

        V[1].fill(0.0);
        V[1][0] = 1.0;

        V[2].fill(1.0);
        V[2][2] = 0.0;

        V[3].fill(0.0);
        V[3][1] = 1.0;

        P[0] = 0.9;
        P[1] = 0.1;
        P[2] = 1.0;
        n.fill(0.0);
        n[2] = 1.0;
        intersect = CGElem::intersectLineSimplex(P, n, V, Q);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN LINE AND SIMPLEX" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    simplex: " << endl;
        cout << "      V: " << V << endl;
        if(intersect){
            cout << "    intersection point is: " << Q << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        V[0].fill(0.0);

        V[1].fill(0.0);
        V[1][0] = 1.0;

        V[2].fill(1.0);
        V[2][2] = 0.0;

        V[3].fill(0.0);
        V[3][1] = 1.0;

        P[0] =-0.1;
        P[1] = 0.9;
        P[2] = 1.0;
        n.fill(0.0);
        n[2] = 1.0;
        intersect = CGElem::intersectLineSimplex(P, n, V, Q);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    simplex: " << endl;
        cout << "      V: " << V << endl;
        if(intersect){
            cout << "    intersection point is: " << Q << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN SEGMENT AND PLANE                                     //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            P, Q0, Q1, P1, n1;
        bool intersect;

        // Compute intersection (Test 1) ---------------------------------------- //
        Q0.fill(0.0);
        Q1.fill(0.0);
        Q1[2] = 1.0;
        P1.fill(0.0);
        P1[2] = 1.0;
        n1.fill(0.0);
        n1[2] = 1.0;
        intersect = CGElem::intersectSegmentPlane(Q0, Q1, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN SEGMENT AND PLANE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    segment: " << endl;
        cout << "      Q0: " << Q0 << endl;
        cout << "      Q1: " << Q1 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection point is: " << P << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        Q0.fill(0.0);
        Q1.fill(0.0);
        Q1[2] = 0.5;
        P1.fill(0.0);
        n1.fill(0.0);
        n1[1] = 1.0;
        intersect = CGElem::intersectSegmentPlane(Q0, Q1, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      Q0: " << Q0 << endl;
        cout << "      Q1: " << Q1 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection point is: " << P << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN SEGMENT AND TRIANGLE                                  //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A, B, C, P0, P1, Q;
        bool intersect ;

        // Compute intersection (Test 1) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;
        P0.fill(0.0);
        P0[0] = P0[1] = 0.5;
        P0[2] = -1.0;
        P1.fill(0.0);
        P1[0] = P1[1] = 0.5;
        P1[2] = 1.0;
        intersect = CGElem::intersectSegmentTriangle(P0, P1, A, B, C, Q);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN LINE AND TRIANGLE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    segment: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      P1: " << P1 << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        if(intersect){
            cout << "    intersection point is: " << Q << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;
        P0.fill(0.0);
        P0[2] = 1.0;
        P1.fill(0.0);
        P1[2] = 2.0;
        intersect = CGElem::intersectSegmentTriangle(P0, P1, A, B, C, Q);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    segment: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      P1: " << P1 << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        if(intersect){
            cout << "    intersection point is: " << Q << endl;
        } else {
            cout << "    do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN PLANE AND PLANE                                       //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            P0, P1, n0, n1;
        array<double, 3>            nL, PL;
        bool intersect;

        // Compute intersection (Test 1) ---------------------------------------- //
        P0.fill(0.0);
        n0.fill(0.0);
        n0[2] = 1.0;
        P1.fill(1.0);
        n1.fill(0.0);
        n1[0] = 1.0;
        intersect = CGElem::intersectPlanePlane(P0, n0, P1, n1, PL, nL);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN PLANES" << endl;
        cout << "  Test 1 " << endl;
        cout << "    1st plane: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    2nd plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection line is: " << "(P: " << PL << ", n: " << nL << ")" << endl;
        } else {
            cout << "    Planes do not intersect " << endl;
        }
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        P0.fill(0.0);
        n0.fill(0.0);
        n0[2] = 1.0;
        P1.fill(0.0);
        n1.fill(0.0);
        n1[2] = 1.0;
        intersect = CGElem::intersectPlanePlane(P0, n0, P1, n1, PL, nL);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    1st plane: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    2nd plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        if(intersect){
            cout << "    intersection line is: " << "(P: " << PL << ", n: " << nL << ")" << endl;
        } else {
            cout << "    Planes do not intersect " << endl;
        }
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN BOX AND BOX                                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A0, A1, B0, B1, C0, C1;
        bool                        intersect ;

        // Compute intersection (Test 1) ---------------------------------------- //
        A0.fill(0.) ;
        A1.fill(1.) ;

        B0 = {{-0.25,0.25,0.25}} ;
        B1 = {{0.75,1.25,0.99}} ;

        intersect = CGElem::intersectBoxBox( A0, A1, B0, B1, C0, C1) ;

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN TWO AXIS ALIGNED BOXES" << endl;
        cout << "  BOX 1 " << endl;
        cout << "     A0: " << A0 << endl;
        cout << "     A1: " << A1 << endl;
        cout << "  BOX 2 " << endl;
        cout << "     B0: " << B0 << endl;
        cout << "     B1: " << B1 << endl;
        cout << "INTERSECTION?" << intersect << endl;
        cout << "     C0: " << C0 << endl;
        cout << "     C1: " << C1 << endl;
        cout << endl;

    }

    // ========================================================================== //
    // Distance Point Triangle
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A, B, C, P0, Q;
        int                         flag ;
        double                      d ;

        // Compute intersection (Test 1) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;

        P0.fill(0.0);
        P0[0] = P0[1] = 1.;
        P0[2] = -2.0;

        d = CGElem::distancePointTriangle(P0, A, B, C, Q, flag);

        // Output message ------------------------------------------------------- //
        cout << " - DISTANCE BETWEEN POINT AND TRIANGLE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    distance to point is: " << d << endl;
        cout << "    projection point is: " << Q << endl;
        cout << "    projecting on : " << flag << endl;
        cout << endl;

    }

    // ========================================================================== //
    // Distance Point Triangle with Calculation of Barycentric Coordinates
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A, B, C, P0, Q, lambda;
        double                      d ;
        int flag ;

        // Compute intersection (Test 1) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;

        P0.fill(0.0);
        P0[0] = P0[1] = 1.;
        P0[2] = -2.0;

        d = CGElem::distancePointTriangle(P0, A, B, C, Q, lambda, flag);

        // Output message ------------------------------------------------------- //
        cout << " - DISTANCE BETWEEN POINT AND TRIANGLE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    distance to point is: " << d << endl;
        cout << "    projection point is: " << Q << endl;
        cout << "    barycentric coordinates are : " << lambda << endl;
        cout << endl;

    }

    // ========================================================================== //
    // Distance Point Cone
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3> point, apex, axis;
        double alpha, distance;

        // Compute distance (Test 1) ---------------------------------------- //
        point = {{ 1., 1., 0. }};
        apex  = {{ 0., 0., 0. }};
        axis  = {{ 1., 0., 0. }};
       
        alpha = 15. /180. *M_PI; 

        distance = CGElem::distancePointCone( point, apex, axis, alpha);

        cout << " - DISTANCE BETWEEN POINT AND CONE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point: " << endl;
        cout << "      p: " << point << endl;
        cout << "    cone: " << endl;
        cout << "      apex: " << apex << endl;
        cout << "      axis: " << axis << endl;
        cout << "      angle: " << alpha << endl;
        cout << "    distance to point is: " << distance << endl;
        cout << endl;

        // Compute distance (Test 2) ---------------------------------------- //
        point = {{ 1., 1., 0. }};
        apex  = {{ 0., 0., 0. }};
        axis  = {{ 1., 0., 0. }};
       
        alpha = 22.5 /180. *M_PI; 

        distance = CGElem::distancePointCone( point, apex, axis, alpha);

        cout << " - DISTANCE BETWEEN POINT AND CONE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point: " << endl;
        cout << "      p: " << point << endl;
        cout << "    cone: " << endl;
        cout << "      apex: " << apex << endl;
        cout << "      axis: " << axis << endl;
        cout << "      angle: " << alpha << endl;
        cout << "    distance to point is: " << distance << endl;
        cout << endl;

        // Compute distance (Test 3) ---------------------------------------- //
        point = {{ -1., -1., -1. }};
        apex  = {{ 0., 0., 0. }};
        axis  = {{ 1., 0., 0. }};
       
        alpha = 22.5/180. *M_PI; 

        distance = CGElem::distancePointCone( point, apex, axis, alpha);

        cout << " - DISTANCE BETWEEN POINT AND CONE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point: " << endl;
        cout << "      p: " << point << endl;
        cout << "    cone: " << endl;
        cout << "      apex: " << apex << endl;
        cout << "      axis: " << axis << endl;
        cout << "      angle: " << alpha << endl;
        cout << "    distance to point is: " << distance << endl;
        cout << endl;
    }

    // ========================================================================== //
    // Distance Cloud Triangle 1st algorithm
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        vector< array<double, 3> >  P(2), Q(2);
        vector<double>              d(2) ;
        vector<int>                 flag(2) ;

        array<double, 3>            A, B, C;

        // Compute intersection (Test 1) ---------------------------------------- //
        A.fill(0.0);
        B.fill(0.0);
        B[0] = 1.0;
        C.fill(0.0);
        C[1] = 1.0;

        P[0].fill(0.0);
        P[0][0] = P[0][1] = 1.;
        P[0][2] = -2.0;

        P[1].fill(-0.5);
        P[1][2] = 1.0;
        
        d = CGElem::distanceCloudTriangle(P, A, B, C, Q, flag);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN POINT CLOUD AND TRIANGLE 1" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point0: " << endl;
        cout << "      P0: " << P[0] << endl;
        cout << "    point1: " << endl;
        cout << "      P1: " << P[1] << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    distance to point 0 is: " << d[0] << endl;
        cout << "    projection point 0 is: " << Q[0] << endl;
        cout << "    projecting 0 on : " << flag[0] << endl;
        cout << "    distance to point 1 is: " << d[1] << endl;
        cout << "    projection point 1 is: " << Q[1] << endl;
        cout << "    projecting 1 on : " << flag[1] << endl;
        cout << endl;

    }

    // ========================================================================== //
    // Distance Cloud Triangle 2nd algorithm
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        vector< array<double, 3> >  P(2), Q(2);
        vector<double>              d(2) ;
        vector< array<double, 3> >  lambda(2) ;

        array<double, 3>            A, B, C;
        std::vector<array<double,3>> VS ;

        // Compute intersection (Test 1) ---------------------------------------- //
        A.fill(0.0);

        B.fill(0.0);
        B[0] = 1.0;

        C.fill(0.0);
        C[1] = 1.0;

        VS.push_back(A) ;
        VS.push_back(B) ;
        VS.push_back(C) ;

        P[0].fill(0.0);
        P[0][0] = P[0][1] = 1.;
        P[0][2] = -2.0;

        P[1].fill(-0.5);
        P[1][2] = 1.0;
        
        //d = CGElem::distanceCloudTriangle(P, A, B, C, nullptr, nullptr);
        d = CGElem::distanceCloudTriangle(P, A, B, C, lambda);
        //d = CGElem::distanceCloudSimplex(P,VS);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN POINT CLOUD AND TRIANGLE 2" << endl;
        cout << "  Test 1 " << endl;
        cout << "    point0: " << endl;
        cout << "      P0: " << P[0] << endl;
        cout << "    point1: " << endl;
        cout << "      P1: " << P[1] << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    distance to point 0 is: " << d[0] << endl;
        cout << "    projection point 0 is: " << CGElem::reconstructPointFromBarycentricTriangle( A, B, C, lambda[0] ) << endl;
        cout << "    barycentric coordinates of projection point: " << lambda[0] << endl;
        cout << "    distance to point 1 is: " << d[1] << endl;
        cout << "    projection point 1 is: " << CGElem::reconstructPointFromBarycentricTriangle( A, B, C, lambda[1] ) << endl;
        cout << "    barycentric coordinates of projection point: " << lambda[1] << endl;
        cout << endl;

    }

    // ========================================================================== //
    // Distance Line Line
    // ========================================================================== //
    {
        array<double,3> P0, n0, P1, n1, xP0, xP1;

        P0 = {{ 1., 1., 0. }};
        n0 = {{ 1., 1., 0. }};
        n0 /= norm2(n0);

        P1 = {{ 4.,-2., 1. }};
        n1 = {{ 4.,-2., 0. }};
        n1 /= norm2(n1);


        double d = CGElem::distanceLineLine(P0,n0,P1,n1,xP0,xP1);

        // Output message ------------------------------------------------------- //
        cout << " - DISTANCE BETWEEN TWO LINES" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line0: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    line1: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        cout << "    distance is: " << d << endl;
        cout << endl;
    }

    return 0; 
};

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    int                             status = 0;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        std::cout << exception.what();
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
