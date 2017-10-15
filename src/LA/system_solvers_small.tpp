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
 * @ingroup system_solver_small
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Solve a linear system of small dimenions using cramer's rule.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in,out] x on output stores the solution to the linear system
*/
template <class T>
void cramer(
    std::vector< std::vector < T > >            &A,
    std::vector< T >                            &B,
    std::vector< T >                            &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                    l, m, n;
T                      dA;
std::vector< std::vector < T > > C;

// Counters
int                    i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; }
n = A[0].size();
if (n == 0) { return; }
if (m != n) {
    return;
}
l = B.size();
if (l == 0) { return; }
if (l != m) {
    return;
}

// =================================================================================== //
// SOLVE LINEAR SYSTEM                                                                 //
// =================================================================================== //

// Solvability condition
dA = det(A);
if (dA < 1.0e-14) {
    return;
}

// Solve linear system
x.resize(n, (T) 0.0);
for (i = 0; i < m; i++) {

    // Build B
    C = A;
    for (j = 0; j < m; j++) {
        C[j][i] = B[j];
    } //next j
    x[i] = det(C)/dA;

} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a linear system of small dimenions using cramer's rule. Overloading
    of cramer() function for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in,out] x on output stores the solution to the linear system
*/
template <class T, size_t m, size_t n>
void cramer(
    std::array< std::array < T, n >, m >        &A,
    std::array< T, m >                          &B,
    std::array< T, n >                          &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
T                           dA;
std::array< std::array < T, n >, m >  C;

// Counters
int                         i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; }
if (n == 0) { return; }
if (m != n) {
    return;
}

// ========================================================================== //
// SOLVE LINEAR SYSTEM                                                        //
// ========================================================================== //

// Solvability condition
dA = det(A);
if (dA < 1.0e-14) {
    return;
}

// Solve linear system
for (i = 0; i < m; i++) {

    // Build B
    C = A;
    for (j = 0; j < m; j++) {
        C[j][i] = B[j];
    } //next j
    x[i] = det(C)/dA;

} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Compute the LU factorization (with row pivoting) for a given non-singular matrix.
    Overloading of LU() function for container array.

    \param[in] A input matrix
    \param[in,out] L lower triangular part (L factor)
    \param[in,out] U upper triangular part (U factor)
    \param[in,out] P permutation matrix

    \result returns an error flag:
        err = 0: no error(s) encountered
        err = 1: matrix is ill conditioned
        err = 2: matrix is singular to working precision
        err = 3: wrong dimensions
*/
template<size_t m>
unsigned int factorizeLU(
    std::array< std::array < double, m >, m >   &A,
    std::array< std::array < double, m >, m >   &L,
    std::array< std::array < double, m >, m >   &U,
    std::array< std::array < double, m >, m >   *P
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double            toll_pivot = 1.0e-8;

// Local variables
int                             info = 0;
int                             pivot_row;
double                          pivot, pivot_trial;
std::array<std::array<double, m>, m>      AA;

// Counter
int                             i, j, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return(3); };

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //

// LU matrices
zeros(L);
zeros(U);

// Backup copy of coeffs. matrix
AA = A;

// Pivoting array
eye(*P);

// ========================================================================== //
// COMPUTE LU FACTORIZATION                                                   //
// ========================================================================== //
for (k = 0; k < m; k++) {
    L[k][k] = 1.0;

    // Pivoting ------------------------------------------------------------- //
    pivot_row = k;
    pivot = std::abs(AA[k][k]);
    for (i = k+1; i < m; i++) {
        pivot_trial = std::abs(AA[i][k]);
        if (pivot_trial > pivot) {
            pivot = pivot_trial;
            pivot_row = i;
        }
    } //next i

    // Perform rows permutation --------------------------------------------- //
    if (pivot_row == k) {
        if (pivot < 1.0e-14) {
            info = 2;
            return(info);
        }
        else if ((pivot >= 1.0e-14) && (pivot < toll_pivot)) {
            info = 1;
        }
    }
    else {
        swap(AA[k], AA[pivot_row]);
        if (P != NULL) {
            swap((*P)[k], (*P)[pivot_row]);
        }
    }

    // Gauss elimination ---------------------------------------------------- //
    for (i = k+1; i < m; i++) {
        L[i][k] = AA[i][k]/AA[k][k] ;
        for (j = k+1; j < m; j++) {
            AA[i][j] = AA[i][j] - L[i][k]*AA[k][j];
        } //next j

    } //next i
    for (j = k; j < m; j++) {
        U[k][j] = AA[k][j];
    } //next j
} //next k

return(info); };

// -------------------------------------------------------------------------- //
/*!
    Solve a lower triangular system using backward substitution. Overloading
    of backwardSubstitution() for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in] x on output stores the solution to the linear system
*/
template<size_t m>
void backwardSubstitution(
    std::array< std::array < double, m >, m >   &A,
    std::array< double, m >                     &B,
    std::array< double, m >                     &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double sum, d;

// Counter
int    i, p;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// CHECK SOLVABILITY CONDITION                                                //
// ========================================================================== //
d = 1.0;
for (i = 0; i < m; i++) {
    d = d*A[i][i];
} //next i
if (std::abs(d) < 1.0e-14) {
    return;
}

// ========================================================================== //
// SOLVE LINEAR SYSTEM WITH BACKWARD SUBSTITUTION                             //
// ========================================================================== //
for (i = m-1; i >= 0; i--) {
    sum = 0.0;
    for(p = m-1; p > i; p--) {
        sum += A[i][p]*x[p];
    } //next p
    x[i] = (B[i] - sum)/A[i][i];
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a upper triangular system using forward substitution. Overloading
    of forwardSubstitution() for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in] x on output stores the solution to the linear system
*/
template<size_t m>
void forwardSubstitution(
    std::array< std::array < double, m >, m >   &A,
    std::array< double, m >                     &B,
    std::array< double, m >                     &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double     d, sum;

// Counters
int        i, p;


// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// CHECK SOLVABILITY CONDITION                                                //
// ========================================================================== //
d = 1.0;
for (i = 0; i < m; i++) {
    d = d*A[i][i];
} //next i
if (std::abs(d) < 1.0e-14) {
    return;
}

// ========================================================================== //
// FORWARD SUBSTITUTION                                                       //
// ========================================================================== //
for(i = 0; i < m; i++) {
    sum = 0.0;
    for(p = 0; p < i; p++) {
        sum += A[i][p] * x[p];
    } //next p
    x[i] = (B[i] - sum)/A[i][i];
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a non-singular linear system using LU factorization. Overloading
    of LU() function for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of linear system
    \param[in,out] x on output stores the solution to the linear system
*/
template<size_t m>
void solveLU(
    std::array< std::array< double, m >, m >    &A,
    std::array< double, m >                     &B,
    std::array< double, m >                     &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int                info;
std::array<std::array<double, m>, m>  L, U, P, *P_ = &P;
std::array<double, m>            z, C;

// Counters
// none

// ========================================================================== //
// COMPUTE LU FACTORIZATION                                                   //
// ========================================================================== //
info = factorizeLU(A, L, U, P_);
if ((info == 2) || (info == 3)) {
    return;
}
matmul(P, B, C);

// ========================================================================== //
//  SOLVE THE LINEAR SYSTEM                                                   //
// ========================================================================== //

// Forward substitution
forwardSubstitution(L, C, z);

// Bacward substitution
backwardSubstitution(U, z, x);

return; };

/*!
 * @}
 */
}
}
