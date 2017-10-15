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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "LinearAlgebra.hpp"

namespace bitpit{

namespace linearalgebra{

/*!
 * @ingroup system_solver_small
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Compute the LU factorization (with partial pivoting) of an input matrix of
    small dimensions.

    \param[in] A input matrix
    \param[in,out] L lower triangular part (factor L)
    \param[in,out] U upper triangular part (factor U)
    \param[in,out] P permutation matrix

    \result error flag:
        err = 0: no error(s) encountered
        err = 1: matrix is ill conditioned
        err = 2: matrix is singular to working precision
        err = 3: wrong dimensions
*/
unsigned int factorizeLU(
    std::vector<std::vector<double>>         &A,
    std::vector<std::vector<double>>         &L,
    std::vector<std::vector<double>>         &U,
    std::vector<std::vector<double>>         *P
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double            toll_pivot = 1.0e-8;

// Local variables
int               info = 0;
int               m, n, pivot_row;
double            pivot, pivot_trial;
std::vector<std::vector<double>> AA;

// Counter
int               i, j, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return(3); };
n = A[0].size();
if (m != n) {
    return (3);
}

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //

// LU matrices
zeros(L, n, n);
zeros(U, n, n);

// Backup copy of coeffs. matrix
AA = A;

// Pivoting array
if (P != NULL) {
    eye(*P, n, n);
}

// ========================================================================== //
// COMPUTE LU FACTORIZATION                                                   //
// ========================================================================== //
for (k = 0; k < n; k++) {
    L[k][k] = 1.0;

    // Pivoting ------------------------------------------------------------- //
    pivot_row = k;
    pivot = std::abs(AA[k][k]);
    for (i = k+1; i < n; i++) {
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
    for (i = k+1; i < n; i++) {
        L[i][k] = AA[i][k]/AA[k][k] ;
        for (j = k+1; j < n; j++) {
            AA[i][j] = AA[i][j] - L[i][k]*AA[k][j];
        } //next j

    } //next i
    for (j = k; j < n; j++) {
        U[k][j] = AA[k][j];
    } //next j
} //next k

return(info); };

// -------------------------------------------------------------------------- //
/*!
    Solve a upper triangular linear system, using backward substitution.

    \param[in] A coeffs. matrix
    \param[in] B r.h.s. of linear system
    \param[in,out] x on output store the solution of the linear system
*/
void backwardSubstitution(
    std::vector<std::vector<double>>     &A,
    std::vector<double>                  &B,
    std::vector<double>                  &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int    m, n, l;
double sum, d;

// Counter
int    i, p;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; };
n = A[0].size();
if (m != n) {
    return;
}
l = B.size();
if (l == 0) { return; };
if (l != n) {
    return;
}

// ========================================================================== //
// CHECK SOLVABILITY CONDITION                                                //
// ========================================================================== //
d = 1.0;
for (i = 0; i < n; i++) {
    d = d*A[i][i];
} //next i
if (std::abs(d) < 1.0e-14) {
    return;
}

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
x.resize(n, 0.0);

// ========================================================================== //
// SOLVE LINEAR SYSTEM WITH BACKWARD SUBSTITUTION                             //
// ========================================================================== //
for (i = n-1; i >= 0; i--) {
    sum = 0.0;
    for(p = n-1; p > i; p--) {
        sum += A[i][p]*x[p];
    } //next p
    x[i] = (B[i] - sum)/A[i][i];
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a lower triangular linear system, using forward substitution.

    \param[in] A coeffs. matrix
    \param[in] B r.h.s. of linear system
    \param[in,out] x on output store the solution of the linear system
*/
void forwardSubstitution(
    std::vector<std::vector<double>>     &A,
    std::vector<double>                  &B,
    std::vector<double>                  &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int        m, n, l;
double     d, sum;

// Counters
int        i, p;


// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; };
n = A[0].size();
if (m != n) {
    return;
}
l = B.size();
if (l == 0) { return; };
if (l != n) {
    return;
}

// ========================================================================== //
// CHECK SOLVABILITY CONDITION                                                //
// ========================================================================== //
d = 1.0;
for (i = 0; i < n; i++) {
    d = d*A[i][i];
} //next i
if (std::abs(d) < 1.0e-14) {
    return;
}

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
x.resize(n, 0.0);

// ========================================================================== //
// FORWARD SUBSTITUTION                                                       //
// ========================================================================== //
for(i = 0; i < n; i++) {
    sum = 0.0;
    for(p = 0; p < i; p++) {
        sum += A[i][p] * x[p];
    } //next p
    x[i] = (B[i] - sum)/A[i][i];
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a linear system of small dimenions (coeffs. matrix must have full rank)
    using LU factorization.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in,out] x on output stores the solution of the linear system
*/
void solveLU(
    std::vector<std::vector<double>>     &A,
    std::vector<double>                  &B,
    std::vector<double>                  &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int    info;
std::vector<std::vector<double>> L, U, P, *P_ = &P;
std::vector<double> z, C;

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
