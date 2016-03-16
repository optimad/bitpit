// ========================================================================== //
//                  - CG - Computational Geometry Library -                   //
//                                                                            //
// Demo for CG algorithms.                                                    //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v1.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include <cmath>
# include <array>
# include <vector>
# include <chrono>
# include <string>
# include <iostream>

# include <bitpit_operators.hpp>

# include "CGBase.hpp"

using namespace std ;

// -------------------------------------------------------------------------- //
int main(
    void
) {

    // ========================================================================== //
    // void Test_Polygon2D(                                                       //
    //     void)                                                                  //
    //                                                                            //
    // Demo for Class_CG_Polygon2D.                                               //
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
    CGPolygon2D::Class_CG_Polygon2D        A, B, C;

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "|=======================================================|" << endl;
        cout << "|          DEMO: Class_CG_Polygon2D                     |" << endl;
        cout << "|                                                       |" << endl;
        cout << "| Description: Demo for 2D polygons methods.            |" << endl;
        cout << "|=======================================================|" << endl;
    }

    // ========================================================================== //
    // DEMO                                                                       //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        vector< array<double,3> >           xA(5), xB(5), xC(5);

        // Initialize scope variables ------------------------------------------- //

        // Vertex for polygon A
        xA[0][0] = 0.0;    xA[0][1] = 0.0;    xA[0][2] = 0.0;
        xA[1][0] = 1.0;    xA[1][1] = 0.0;    xA[1][2] = 0.0;
        xA[2][0] = 1.5;    xA[2][1] = 0.5;    xA[2][2] = 0.0;
        xA[3][0] = 1.0;    xA[3][1] = 1.0;    xA[3][2] = 0.0;
        xA[4][0] = 0.0;    xA[4][1] = 1.0;    xA[4][2] = 0.0;

        // Vertex for polygon B
        xB[0][0] = 0.0;    xB[0][1] = 0.0;    xB[0][2] = 0.0;
        xB[1][0] = 0.0;    xB[1][1] = 1.0;    xB[1][2] = 0.0;
        xB[2][0] = 1.0;    xB[2][1] = 1.0;    xB[2][2] = 0.0;
        xB[3][0] = 0.5;    xB[3][1] = 0.5;    xB[3][2] = 0.0;
        xB[4][0] = 1.0;    xB[4][1] = 0.0;    xB[4][2] = 0.0;

        // Vertex for polygon C
        xC[0][0] = 0.0;    xC[0][1] = 0.0;    xC[0][2] = 0.0;
        xC[1][0] = 0.0;    xC[1][1] = 1.0;    xC[1][2] = 0.0;
        xC[2][0] = 1.0;    xC[2][1] = 1.0;    xC[2][2] = 0.0;
        xC[3][0] = 1.5;    xC[3][1] = 0.5;    xC[3][2] = 0.0;
        xC[4][0] = 1.0;    xC[4][1] = 0.0;    xC[4][2] = 0.0;

        // Display status after creation ---------------------------------------- //
        cout << " STATUS AFTER INITIALIZATION" << endl;

        cout << "A:";
        A.display(cout);

        cout << "B:";
        B.display(cout);

        cout << "C:";
        C.display(cout);

        cout << endl;

        // Display status after vertex insertion -------------------------------- //
        cout << " STATUS AFTER VERTEX INSERTION" << endl;

        A.AddVertices(xA);
        cout << "A:";
        A.display(cout);

        B.AddVertices(xB);
        cout << "B:";
        B.display(cout);

        C.AddVertices(xC);
        cout << "C:";
        C.display(cout);

        cout << endl;

        // Display status after copy -------------------------------------------- //
        cout << " STATUS AFTER COPY C <- B" << endl;
        C = B;

        cout << "B:";
        B.display(cout);

        cout << "C:";
        C.display(cout);

        cout << endl;

        // Inclusion check ------------------------------------------------------ //
        cout << " INCLUSION CHECK" << endl;

        array<double, 3>                P0, P1, P2;
        P0[0] = 0.1;    P0[1] = 0.1;   P0[2] = 0.0 ;
        P1[0] = -0.1;   P1[1] = 0.1;   P1[2] = 0.0 ;
        P2[0] = 0.7;    P2[1] = 0.5;   P2[2] = 0.0 ;

        cout << "A does ";
        if (!A.IncludePoint(P0)) { cout << "not ";  }
        cout << "contain point P0 = " << P0 << endl;
        cout << "A does ";
        if (!A.IncludePoint(P1)) { cout << "not ";  }
        cout << "contain point P1 = " << P1 << endl;
        cout << "A does ";
        if (!A.IncludePoint(P2)) { cout << "not ";  }
        cout << "contain point P2 = " << P2 << endl;

        cout << "B does ";
        if (!B.IncludePoint(P0)) { cout << "not ";  }
        cout << "contain point P0 = " << P0 << endl;
        cout << "B does ";
        if (!B.IncludePoint(P1)) { cout << "not ";  }
        cout << "contain point P1 = " << P1 << endl;
        cout << "B does ";
        if (!B.IncludePoint(P2)) { cout << "not ";  }
        cout << "contain point P2 = " << P2 << endl;

        cout << endl;

    }
    return 0; 

};
