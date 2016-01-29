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
 *  as published by by the Free Software Foundation.
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
//                         LINEAR ALGEBRA PACKAGE                             //
//                                                                            //
// Demo for linear algebra routines.                                          //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author            : Alessandro Alaia                                       //
// Data              : Sept 26, 2014                                          //
// Version           : v2.0                                                   //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "test_LA_00001.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Test_000(
    void
) {

// ========================================================================== //
// void Test_000(                                                             //
//     void)                                                                  //
//                                                                            //
// Basic matrix template demo.                                                //
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
// none

// Counters
// none

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "BASIC MATRIX TEMPLATES - DEMO" << endl;
}

// ========================================================================== //
// MATRIX TEMPLATES FOR VECTOR<VECTOR<*>>                                     //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    dvector2D                   A;

    // Output message ------------------------------------------------------- //
    cout << "Matrix as vector<vector<*>>" << endl;

    // Matrix templates ----------------------------------------------------- //

    // Zeros
    cout << "zeros(A, 4, 3) = " << endl;
    bitpit::zeros(A, 4, 3);
    bitpit::display_matrix(cout, A);

    // Ones
    cout << "ones(A, 3, 2) = " << endl;
    bitpit::ones(A, 3, 2);
    bitpit::display_matrix(cout, A);

    // Eye
    cout << "eye(A, 5, 5) = " << endl;
    bitpit::eye(A, 5, 5);
    bitpit::display_matrix(cout, A);

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// MATRIX TEMPLATES FOR ARRAY<ARRAY<*, n>, m>                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    array<array<double, 3>, 4>           A1;
    array<array<double, 2>, 3>           A2;
    array<array<double, 5>, 5>           A3;

    // Output message ------------------------------------------------------- //
    cout << "Matrix as array<array<*, n>, m>" << endl;

    // Matrix templates ----------------------------------------------------- //

    // Zeros
    cout << "zeros(A, 4, 3) = " << endl;
    bitpit::zeros(A1);
    bitpit::display_matrix(cout, A1);

    // Ones
    cout << "ones(A, 3, 2) = " << endl;
    bitpit::ones(A2);
    bitpit::display_matrix(cout, A2);

    // Eye
    cout << "eye(A, 5, 5) = " << endl;
    bitpit::eye(A3);
    bitpit::display_matrix(cout, A3);

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "DEMO: done!!" << endl;
}

return; };

// -------------------------------------------------------------------------- //
void Test_001(
    void
) {

// ========================================================================== //
// void Test_001(                                                             //
//     void)                                                                  //
//                                                                            //
// Basic matrix operations demo.                                              //
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
// none

// Counters
// none

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "BASIC MATRIX OPERATIONS - DEMO" << endl;
}

// ========================================================================== //
// MATRIX OPERATIONS FOR VECTOR<VECTOR<*>>                                    //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         i, j;
    int                         k = 2;
    ivector1D                   a(4, 0), b(3, 0), c;
    ivector2D                   A(4, ivector1D(3, 0));
    ivector2D                   B(3, ivector1D(4, 0));
    ivector2D                   P;

    // Initialize scope variables ------------------------------------------- //
    for (i = 0; i < 4; ++i) {
        a[i] = i;
    } //next i
    for (i = 0; i < 3; ++i) {
        b[i] = i;
    } //next i
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 3; ++j) {
            A[i][j] = i+j;
        } //next j
    } //next i
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 4; ++j) {
            B[i][j] = i+j;
        } //next j
    } //next i

    // Output message ------------------------------------------------------- //
    cout << "Matrix as vector<vector<*>>" << endl;
    cout << "A = " << endl;
    bitpit::display_matrix(cout, A);
    cout << "B = " << endl;
    bitpit::display_matrix(cout, B);
    cout << "a = " << endl << a << endl;
    cout << "b = " << endl << b << endl;
    cout << "k = " << k << endl;
    cout << endl;

    // Matrix product ------------------------------------------------------- //
    cout << "bitpit::matmul(k, A, P); P = " << endl;
    bitpit::matmul(k, A, P);
    bitpit::display_matrix(cout, P);
    cout << "bitpit::matmul(A, k, P); P = " << endl;
    bitpit::matmul(A, k, P);
    bitpit::display_matrix(cout, P);
    cout << endl;
    cout << "bitpit::matmul(a, A, c); c = " << endl;
    bitpit::matmul(a, A, c);
    cout << c << endl;;
    cout << "bitpit::matmul(A, b, c); c = " << endl;
    bitpit::matmul(A, b, c);
    cout << c << endl;;
    cout << endl;
    cout << "bitpit::matmul(A, B, P); P = " << endl;
    bitpit::matmul(A, B, P);
    bitpit::display_matrix(cout, P);

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// MATRIX OPERATIONS FOR ARRAY<ARRAY<*, n>, m>                                //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                                  i, j;
    int                                  k = 2;
    array<array<int, 3>, 4>              A, P;
    array<array<int, 4>, 3>              B;
    array<array<int, 4>, 4>              Q;
    array<int, 4>                        a, d;
    array<int, 3>                        b, c;

    // Initialize scope variables ------------------------------------------- //
    for (i = 0; i < 4; ++i) {
        a[i] = i;
    } //next i
    for (i = 0; i < 3; ++i) {
        b[i] = i;
    } //next i
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 3; ++j) {
            A[i][j] = i+j;
        } //next j
    } //next i
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 4; ++j) {
            B[i][j] = i+j;
        } //next j
    } //next i

    // Output message ------------------------------------------------------- //
    cout << "Matrix as array<array<*, n>, m>" << endl;
    cout << "A = " << endl;
    bitpit::display_matrix(cout, A);
    cout << "B = " << endl;
    bitpit::display_matrix(cout, B);
    cout << "a = " << endl << a << endl;
    cout << "b = " << endl << b << endl;
    cout << "k = " << k << endl;
    cout << endl;

    // Matrix product ------------------------------------------------------- //
    cout << "bitpit::matmul(k, A, P); P = " << endl;
    bitpit::matmul(k, A, P);
    bitpit::display_matrix(cout, P);
    cout << "bitpit::matmul(A, k, P); P = " << endl;
    bitpit::matmul(A, k, P);
    bitpit::display_matrix(cout, P);
    cout << endl;

    cout << "bitpit::matmul(a, A, c); c = " << endl;
    bitpit::matmul(a, A, c);
    cout << c << endl;
    cout << "bitpit::matmul(A, b, d); d = " << endl;
    bitpit::matmul(A, b, d);
    cout << d << endl;
    cout << endl;

    cout << "bitpit::matmul(A, B, Q); Q = " << endl;
    bitpit::matmul(A, B, Q);
    bitpit::display_matrix(cout, Q);

    // Closing header ------------------------------------------------------- //
    cout << endl;

}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "DEMO: done!!" << endl;
}

return; };

// -------------------------------------------------------------------------- //
void Test_002(
    void
) {

// ========================================================================== //
// void Test_002(                                                             //
//     void)                                                                  //
//                                                                            //
// Basic matrix manipulation demo.                                            //
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
// none

// Counters
// none

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "BASIC MATRIX MINPULATION - DEMO" << endl;
}

// ========================================================================== //
// MATRIX MANIPULATIONS FOR VECTOR<VECTOR<*>>                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         i, j;
    ivector2D                   A(4, ivector1D(3, 0));
    ivector2D                   P;

    // Initialize scope variables ------------------------------------------- //
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 3; ++j) {
            A[i][j] = i+j;
        } //next j
    } //next i

    // Output message ------------------------------------------------------- //
    cout << "Matrix as vector<vector<*>>" << endl;
    cout << "A = " << endl;
    bitpit::display_matrix(cout, A);
    cout << endl;

    // Matrix manipulations ------------------------------------------------- //
    cout << "bitpit::transpose(A, P); P = " << endl;
    bitpit::transpose(A, P);
    bitpit::display_matrix(cout, P);
    cout << endl;
    cout << "bitpit::complement(1, 1, A, P); P = " << endl;
    bitpit::complement(1, 1, A, P);
    bitpit::display_matrix(cout, P);
    cout << endl;
    cout << "bitpit::triL(A, P); P = " << endl;
    bitpit::triL(A, P);
    bitpit::display_matrix(cout, P);
    cout << endl;
    cout << "bitpit::triU(A, P); P = " << endl;
    bitpit::triU(A, P);
    bitpit::display_matrix(cout, P);
    cout << endl;

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// MATRIX MANIPULATIONS FOR ARRAY<ARRAY<*, n>, m>                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         i, j;
    array<array<double, 3>, 4>  A;
    array<array<double, 3>, 4>  P;
    array<array<double, 2>, 3>  C;
    array<array<double, 4>, 3>  T;

    // Initialize scope variables ------------------------------------------- //
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 3; ++j) {
            A[i][j] = i+j;
        } //next j
    } //next i

    // Output message ------------------------------------------------------- //
    cout << "Matrix as vector<vector<*>>" << endl;
    cout << "A = " << endl;
    bitpit::display_matrix(cout, A);
    cout << endl;

    // Matrix manipulations ------------------------------------------------- //
    cout << "bitpit::transpose(A, P); P = " << endl;
    bitpit::transpose(A, T);
    bitpit::display_matrix(cout, T);
    cout << endl;
    cout << "bitpit::complement(1, 1, A, P); P = " << endl;
    bitpit::complement(1, 1, A, C);
    bitpit::display_matrix(cout, C);
    cout << endl;
    cout << "bitpit::triL(A, P); P = " << endl;
    bitpit::triL(A, P);
    bitpit::display_matrix(cout, P);
    cout << endl;
    cout << "bitpit::triU(A, P); P = " << endl;
    bitpit::triU(A, P);
    bitpit::display_matrix(cout, P);
    cout << endl;

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "DEMO: done!!" << endl;
}

return; };

// -------------------------------------------------------------------------- //
void Test_003(
    void
) {

// ========================================================================== //
// void Test_003(                                                             //
//     void)                                                                  //
//                                                                            //
// Basic linear system solver demo.                                           //
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
// none

// Counters
// none

// ========================================================================== //
// OUTPUT MESSAGE                                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "BASIC LINEAR SYSTEM SOLVER - DEMO" << endl;
}

// ========================================================================== //
// MATRIX DETERMINANT FOR VECTOR<VECTOR<*>>                                   //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         i;
    dvector2D                   A(3, dvector1D(3, 0.0)), B(3, dvector1D(3, 0.0));
    dvector1D                   b(3, 0.0), c, d;
    dvector2D                   P;

    // Initialize scope variables ------------------------------------------- //
    bitpit::eye(B, 3, 3);
    bitpit::matmul(2.0, B, A);
    for (i = 0; i < 3; ++i) {
        b[i] = ((double) i);
    } //next i

    // Output message ------------------------------------------------------- //
    cout << "Matrix as vector<vector<*>>" << endl;
    cout << "A = " << endl;
    bitpit::display_matrix(cout, A);
    cout << "b = " << b << endl;
    cout << endl;

    // Matrix manipulations ------------------------------------------------- //
    cout << "bitpit::det(A) = " << bitpit::det(A) << endl;;
    cout << endl;
    cout << "bitpit::Cramer(A, b, c); c = " << endl;
    bitpit::Cramer(A, b, c);
    cout << c << endl;
    cout << endl;
    cout << "bitpit::SolveLU(A, b, d); d = " << endl;
    bitpit::SolveLU(A, b, d);
    cout << d << endl;

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// MATRIX DETERMINANT FOR ARRAY<ARRAY<*, n>, m>                               //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                         i;
    array<array<double, 3>, 3>  A, B;
    array<double, 3>            b, c, d;
    dvector2D                   P;

    // Initialize scope variables ------------------------------------------- //
    bitpit::eye(B);
    bitpit::matmul(2.0, B, A);
    for (i = 0; i < 3; ++i) {
        b[i] = ((double) i);
    } //next i

    // Output message ------------------------------------------------------- //
    cout << "Matrix as array<array<*, n>, m>" << endl;
    cout << "A = " << endl;
    bitpit::display_matrix(cout, A);
    cout << "b = " << b << endl;
    cout << endl;

    // Matrix manipulations ------------------------------------------------- //
    cout << "bitpit::det(A) = " << bitpit::det(A) << endl;;
    cout << endl;
    cout << "bitpit::Cramer(A, b, c); c = " << endl;
    bitpit::Cramer(A, b, c);
    cout << c << endl;
    cout << endl;
    cout << "bitpit::SolveLU(A, b, d); d = " << endl;
    bitpit::SolveLU(A, b, d);
    cout << d << endl;

    // Closing header ------------------------------------------------------- //
    cout << endl;
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "DEMO: done!!" << endl;
}

return; };

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main()
{

    Test_000();
    Test_001();
    Test_002();
    Test_003();

return(0);

};
