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
//                         LINEAR ALGEBRA PACKAGE                             //
//                                                                            //
// Functions for basic linear algebra computations.                           //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author            : Alessandro Alaia                                       //
// Data              : Sept 26, 2014                                          //
// Version           : v2.0                                                   //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

namespace bitpit{

namespace linearalgebra{
/*!
 * @ingroup   ladisplay
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Display matrix to output stream in a nicely formatted form.

    \param[in,out] out output stream
    \param[in] A matrix to be displayed
*/
template<class T>
void display(
    std::ostream                                &out,
    std::vector< std::vector< T > >             &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int             m = A.size();

// Counters
int             i;

// ========================================================================== //
// DISPLAY MATRIX CONTENT                                                     //
// ========================================================================== //
for (i = 0; i < m; ++i) {
    out << A[i] << std::endl;
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Display matrix to output stream in a nicely formatted form.
    Overloading of display() function for container array.

    \param[in,out] out output stream
    \param[in] A matrix to be displayed
*/
template<class T, size_t m, size_t n>
void display(
    std::ostream                                &out,
    std::array<std::array<T, n>, m>             &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
size_t          i;

// ========================================================================== //
// DISPLAY MATRIX CONTENT                                                     //
// ========================================================================== //
for (i = 0; i < m; ++i) {
    out << A[i] << std::endl;
} //next i

return; };

/*!
 * @}
 */

/*!
 * @ingroup laauxiliary
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix A, compute the complement of the element A[i][j].

    \param[in] i row index of element
    \param[in] j column index of element
    \param[in] A input matrix
    \param[in,out] B complement of alement A[i][j].
*/
template <class T>
void complement(
    int                                          i,
    int                                          j,
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            m, n;

// Counters
int            l, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; }
n = A[0].size();
if (n == 0) { return; }
i--; j--;
if ((i >= m) || (i < 0)) {
    return;
}
if ((j >= n) || (j < 0)) {
    return;
}

// ========================================================================== //
// EXTRACT COMPLEMENT                                                         //
// ========================================================================== //

// Resize output variables
B.resize(m-1);
for (i = 0; i < m-1; ++i) {
    B[i].resize(n-1, 0.0);
} //next i

// Extract complement
for (l = 0; l < i; l++) {
    for (k = 0; k < j; k++) {
        B[l][k] = A[l][k];
    } //next k
    for (k = j+1; k < n; k++) {
        B[l][k-1] = A[l][k];
    } //next k
} //next l
for (l = i+1; l < m; l++) {
    for (k = 0; k < j; k++) {
        B[l-1][k] = A[l][k];
    } //next k
    for (k = j+1; k < n; k++) {
        B[l-1][k-1] = A[l][k];
    } //next k
} //next l

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix A, compute the complement of the element A[i][j].
    Overloading of complement() function for container array.

    \param[in] i row index of element
    \param[in] j column index of element
    \param[in] A input matrix
    \param[in,out] B complement of alement A[i][j].
*/
template <class T, size_t m, size_t n>
void complement(
    int                                          i,
    int                                          j,
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array<T, n-1>, m-1>        &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int            l, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; }
if (n == 0) { return; }
i--; j--;
if ((i >= (long) m) || (i < 0)) {
    return;
}
if ((j >= (long) n) || (j < 0)) {
    return;
}

// ========================================================================== //
// EXTRACT COMPLEMENT                                                         //
// ========================================================================== //
for (l = 0; l < i; l++) {
    for (k = 0; k < j; k++) {
        B[l][k] = A[l][k];
    } //next k
    for (k = j+1; k < (long) n; k++) {
        B[l][k-1] = A[l][k];
    } //next k
} //next l
for (l = i+1; l < (long) m; l++) {
    for (k = 0; k < j; k++) {
        B[l-1][k] = A[l][k];
    } //next k
    for (k = j+1; k < (long) n; k++) {
        B[l-1][k-1] = A[l][k];
    } //next k
} //next l

return; };

/*!
 * @}
 */

/*!
 * @ingroup laspecialmatrix
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n matrix with 0 entries

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T>
void zeros(
    std::vector< std::vector < T > >            &A,
    int                                          m,
    int                                          n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int      i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
A.resize(m);
for (i = 0; i < m; i++) {
    A[i].resize(n, (T) 0.0);
    for (j = 0; j < n; j++) {
        A[i][j] = (T) 0.0;
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initalize a m-by-n matrix with 0-entries. Overloading of zeros() function
    for container array.

    \param[in,out] A container for matrix storage
    \tparam T data type
    \tparam m number of matrix rows
    \tparam n number of matrix columns
*/
template <class T, size_t m, size_t n>
void zeros(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// template <class T, size_t m, size_t n>                                     //
// void zeros(                                                                //
//     array< array < T, n >, m >  &A)                                        //
//                                                                            //
// Initialize a m-by-n matrix of zeros.                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A      : array< array< T, n >, m >, with m-by-n matrix of zeros          //
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
size_t   i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
for (i = 0; i < m; i++) {
    A[i].fill(0.0);
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n matrix of ones.

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T>
void ones(
    std::vector< std::vector < T > >            &A,
    int                                          m,
    int                                          n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int      i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    std::cout << "ERROR: number of rows (columns) must be > 0!!" << std::endl;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
A.resize(m);
for (i = 0; i < m; i++) {
    A[i].resize(n, (T) 1.0);
    for (j = 0; j < n; j++) {
        A[i][j] = (T) 1.0;
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n matrix of ones. Overloading of ones() function for
    container array.

    \param[in,out] A container for matrix storage
    \tparam T data type
    \tparam m number of matrix rows
    \tparam n number of matrix columns
*/
template <class T, size_t m, size_t n>
void ones(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
size_t   i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
for (i = 0; i < m; i++) {
    A[i].fill(1.0);
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n identity matrix having. All the entries on the main
    diagonal are set to 1.
    \tparam T data type
    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T>
void eye(
    std::vector< std::vector < T > >            &A,
    int                                          m,
    int                                          n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int      s = std::min(m, n);

// Counters
int      i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
zeros(A, m, n);
for (i = 0; i < s; i++) {
    A[i][i] = (T) 1.0;
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n identity matrix having. All the entries on the main
    diagonal are set to 1. Overloading of eye() function for container array.

    \param[in,out] A container for matrix storage
    \tparam T data type
    \tparam m number of matrix rows
    \tparam n number of matrix columns
*/
template <class T, size_t m, size_t n>
void eye(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int      s = std::min(m, n);

// Counters
int      i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
zeros(A);
for (i = 0; i < s; i++) {
    A[i][i] = (T) 1.0;
} //next i

return; };

/*!
 * @}
 */


// Matrix determinant ======================================================= //

/*!
 * @ingroup lainfo 
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    End function for recursive calls to det().
    \tparam T data type
    \param[in] A input matrix
    \result determinant of A.
*/
template <class T>
T det(
    std::array< std::array< T, 1 >, 1 >         &A
) {

return(A[0][0]); };

// -------------------------------------------------------------------------- //
/*!
    Compute determinant of a matrix of small dimenions using Laplace rule.
    \tparam T data type
    \param[in] A input matrix
    \result determinant of A.
*/
template <class T>
T det(
    std::vector< std::vector < T > >            &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                     m, n;
T                       d = (T) 1.0e+18;
std::vector< std::vector < T > >  C;

// Counters
int                     i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return(d); };
n = A[0].size();
if (n == 0) { return(d); };
if (m != n) {
    return(d);
}

// ========================================================================== //
// COMPUTE DETERMINANT                                                        //
// ========================================================================== //
if (m == 1) {
    d = A[0][0];
}
else if (m == 2) {
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0];
}
else if (m == 3) {
    d = A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1]
        + A[0][1]*A[1][2]*A[2][0] - A[0][1]*A[1][0]*A[2][2]
        + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0];
}
else {
    d = (T) 0.0;
    for (i = 0; i < m; i++) {
        complement(1, i+1, A, C);
        d += pow(-1.0, i+2) * A[0][i] * det(C);
    } //next i
}

return(d); };

// -------------------------------------------------------------------------- //
/*!
    Compute determinant of a matrix of small dimenions using Laplace rule.
    Overloading of det() function for container array.
    \tparam T data type
    \tparam m number of matrix rows
    \tparam n number of matrix columns
    \param[in] A input matrix
    \result determinant of A.
*/
template <class T, size_t m, size_t n>
T det(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
T                       d = (T) 0.0;

// Counters
int                     i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return(d); };
if (n == 0) { return(d); };
if (m != n) {
    return(d);
}

// ========================================================================== //
// COMPUTE DETERMINANT                                                        //
// ========================================================================== //
if (m == 1) {
    d = A[0][0];
}
else if (m == 2) {
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0];
}
else if (m == 3) {
    d = A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1]
        + A[0][1]*A[1][2]*A[2][0] - A[0][1]*A[1][0]*A[2][2]
        + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0];
}
else {
    std::array< std::array < double, m-1 >, m-1 >     C;
    for (i = 0; i < m; i++) {
        complement(1, i+1, A, C);
        d += pow(-1.0, i+2) * A[0][i] * det(C);
    } //next i
}

return(d); };


/*!
 * @}
 */
}

}
