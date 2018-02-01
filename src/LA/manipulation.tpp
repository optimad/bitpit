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

namespace bitpit{
namespace linearalgebra{

/*!
 * @ingroup lamanipulation
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose.

    \param[in] A input matrix
    \param[in,out] B transpose of A
*/
template <class T>
void transpose(
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            m, n;

// Counters
int            i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; }
n = A[0].size();
if (n == 0) { return; }

// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //

// Resize output variables
B.resize(n);
for (i = 0; i < n; ++i) {
    B[i].resize(m, (T) 0.0);
} //next i

// Transposition
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose. Overloading of transpose()
    function for container array.

    \param[in] A input matrix
    \param[in,out] B transpose of A
*/
template <class T, size_t m, size_t n>
void transpose(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, m >, n >         &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
size_t         i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; }
if (n == 0) { return; }

// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose. Overloading of transpose()
    function.

    \param[in] A input matrix
    \return B transpose of A
*/
template <class T>
std::vector< std::vector< T > > transpose(
    const std::vector< std::vector< T > >             &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                        m = A.size(); 
int                        n(0);

// Counters
int            i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m != 0) n = A[0].size() ;

// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //

std::vector< std::vector< T > >      B( n, std::vector<T> (m,0.) ) ;

// Transposition
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return B; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose. Overloading of transpose()
    function for container array.

    \param[in] A input matrix
    \return B transpose of A
*/
template <class T, size_t m, size_t n>
std::array< std::array< T, m >, n > transpose(
    const std::array< std::array< T, n >, m >         &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
std::array< std::array< T, m >, n >  B ;

// Counters
size_t         i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return B; };

// -------------------------------------------------------------------------- //
/*!
    Extract the lower triangular part of a matrix.

    \param[in] A input matrix
    \param[in,out] L lower triangular part of A
*/
template <class T>
void triL(
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &L
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int        m, n;

// Counters
int        i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
m = A.size();
if (m == 0) { return; };
n = A[0].size();

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(L, m, n);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j <= i; j++) {
        L[i][j] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the lower triangular part of a matrix. Overloading of triL() function
    for container array.

    \param[in] A input matrix
    \param[in,out] L lower triangular part of A
*/

template <class T, size_t m, size_t n>
void triL(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, n >, m >         &L
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
size_t     i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(L);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j <= i; j++) {
        L[i][j] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the upper triangular part of a matrix.

    \param[in] A input matrix
    \param[in,out] U upper triangular part of A
*/
template <class T>
void triU(
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &U
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int        m, n;

// Counters
int        i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
m = A.size();
if (m == 0) { return; };
n = A[0].size();

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(U, m, n);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = i; j < n; j++) {
        U[i][j] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the upper triangular part of a matrix. Overloading of triU() function
    for container array.

    \param[in] A input matrix
    \param[in,out] U upper triangular part of A
*/
template <class T, size_t m, size_t n>
void triU(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, n >, m >         &U
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
size_t     i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(U);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = i; j < n; j++) {
        U[i][j] = A[i][j];
    } //next j
} //next i

return; };

/*!
 * @}
 */
}
}
