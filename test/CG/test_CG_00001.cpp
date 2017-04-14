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

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

using namespace std;
using namespace bitpit;

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //
// -------------------------------------------------------------------------- //
int main(
    void
    ) {

    // ========================================================================== //
    // void Test_intersection(                                                    //
    //     void)                                                                  //
    //                                                                            //
    // Demo for intersections routines.                                           //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables

    // Counters

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "|=======================================================|" << endl;
        cout << "|       DEMO: intersection between geom. entities       |" << endl;
        cout << "|                                                       |" << endl;
        cout << "| Description: Demo for basic CG routines.              |" << endl;
        cout << "|=======================================================|" << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN PLANES                                                //
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
    // INTERSECTION BETWEEN PLANE AND LINE                                        //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            P, P0, P1, n0, n1;

        // Compute intersection (Test 1) ---------------------------------------- //
        P0.fill(0.0);
        n0.fill(0.0);
        n0[2] = 1.0;
        P1.fill(0.0);
        P1[2] = 1.0;
        n1.fill(0.0);
        n1[2] = 1.0;
        CGElem::intersectLinePlane(P0, n0, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN LINE AND PLANE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        cout << "    intersection point is: " << P << endl;
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        P0.fill(0.0);
        n0.fill(0.0);
        n0[2] = 1.0;
        P1.fill(0.0);
        n1.fill(0.0);
        n1[1] = 1.0;
        CGElem::intersectLinePlane(P0, n0, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      n0: " << n0 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        cout << "    intersection points is: " << P << endl;
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN PLANE AND SEGMENT                                     //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            P, Q0, Q1, P1, n1;

        // Compute intersection (Test 1) ---------------------------------------- //
        Q0.fill(0.0);
        Q1.fill(0.0);
        Q1[2] = 1.0;
        P1.fill(0.0);
        P1[2] = 1.0;
        n1.fill(0.0);
        n1[2] = 1.0;
        CGElem::intersectSegmentPlane(Q0, Q1, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN SEGMENT AND PLANE" << endl;
        cout << "  Test 1 " << endl;
        cout << "    segment: " << endl;
        cout << "      Q0: " << Q0 << endl;
        cout << "      Q1: " << Q1 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        cout << "    intersection point is: " << P << endl;
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        Q0.fill(0.0);
        Q1.fill(0.0);
        Q1[2] = 0.5;
        P1.fill(0.0);
        n1.fill(0.0);
        n1[1] = 1.0;
        CGElem::intersectSegmentPlane(Q0, Q1, P1, n1, P);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      Q0: " << Q0 << endl;
        cout << "      Q1: " << Q1 << endl;
        cout << "    plane: " << endl;
        cout << "      P1: " << P1 << endl;
        cout << "      n1: " << n1 << endl;
        cout << "    intersection points is: " << P << endl;
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN LINE AND TRIANGLE                                     //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A, B, C, P, Q, n;

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
        CGElem::intersectLineTriangle(P, n, A, B, C, Q);

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
        cout << "    intersection point is: " << Q << endl;
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
        CGElem::intersectLineTriangle( P, n, A, B, C, Q);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    intersection point is: " << Q << endl;
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN SEGMENT AND TRIANGLE                                  //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        array<double, 3>            A, B, C, P0, P1, Q;

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
        CGElem::intersectSegmentTriangle(P0, P1, A, B, C, Q);

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
        cout << "    intersection point is: " << Q << endl;
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
        CGElem::intersectSegmentTriangle(P0, P1, A, B, C, Q);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    segment: " << endl;
        cout << "      P0: " << P0 << endl;
        cout << "      P1: " << P1 << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    intersection point is: " << Q << endl;
        cout << endl;
    }

    // ========================================================================== //
    // INTERSECTION BETWEEN LINE AND SIMPLEX                                      //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        vector< array<double, 3> >      V(4);
        array<double, 3>                P, Q, n;

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
        CGElem::intersectLineSimplex(P, n, V, Q);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN LINE AND SIMPLEX" << endl;
        cout << "  Test 1 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    simplex: " << endl;
        cout << "      V: " << V << endl;
        cout << "    intersection point is: " << Q << endl;
        cout << endl;

        // Compute intersection (Test 2) ---------------------------------------- //
        V[0].fill(0.0);

        V[1].fill(0.0);
        V[1][0] = 1.0;

        V[2].fill(1.0);
        V[2][2] = 0.0;

        V[3].fill(0.0);
        V[3][1] = 1.0;

        P[0] = 0.1;
        P[1] = 0.9;
        P[2] = 1.0;
        n.fill(0.0);
        n[2] = 1.0;
        CGElem::intersectLineSimplex(P, n, V, Q);

        // Output message ------------------------------------------------------- //
        cout << "  Test 2 " << endl;
        cout << "    line: " << endl;
        cout << "      P: " << P << endl;
        cout << "      n: " << n << endl;
        cout << "    simplex: " << endl;
        cout << "      V: " << V << endl;
        cout << "    intersection point is: " << Q << endl;
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
        cout << "    projecting on : " << lambda << endl;
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
        cout << "    segment: " << endl;
        cout << "      P0: " << P[0] << endl;
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
        //d = CGElem::distanceCloudTriangle(P, A, B, C, &Q, &lambda);
        d = CGElem::distanceCloudSimplex(P,VS);

        // Output message ------------------------------------------------------- //
        cout << " - INTERSECTION BETWEEN POINT CLOUD AND TRIANGLE 2" << endl;
        cout << "  Test 1 " << endl;
        cout << "    segment: " << endl;
        cout << "      P0: " << P[0] << endl;
        cout << "      P1: " << P[1] << endl;
        cout << "    triangle: " << endl;
        cout << "      A: " << A << endl;
        cout << "      B: " << B << endl;
        cout << "      C: " << C << endl;
        cout << "    distance to point 0 is: " << d[0] << endl;
        cout << "    projection point 0 is: " << Q[0] << endl;
        cout << "    barycentric coordinates of projection point: " << lambda[0] << endl;
        cout << "    distance to point 1 is: " << d[1] << endl;
        cout << "    projection point 1 is: " << Q[1] << endl;
        cout << "    barycentric coordinates of projection point: " << lambda[1] << endl;
        cout << endl;

    }

    // ========================================================================== //
    // intersection Boxes
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
    // Point in Box
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




    return 0; 
};


