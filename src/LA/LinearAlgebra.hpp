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
# ifndef __BITPIT_LINEAR_ALGEBRA_HPP__
# define __BITPIT_LINEAR_ALGEBRA_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard template library
# include <vector>
# include <cmath>
# include <array>
# include <iostream>

# include "Operators.hpp"

namespace bitpit{

/*!
 * @brief   collection of functions to create and solve small linear systems
 */
namespace linearalgebra{

// Generic routines --------------------------------------------------------- //
template <class T>
void display(                                                          // Display matrix in nicely formatted output
    std::ostream                                &,                            // (input/output) handle to output stream
    std::vector< std::vector< T > >             &                             // (input) matrix to be displayed
); 

template <class T, size_t m, size_t n>
void display(                                                          // Display matrix in nicely formatted output
    std::ostream                                &,                            // (input/output) handle to output stream
    std::array< std::array< T, n >, m >         &                             // (input) matrix to be displayed
);

// Matrix Basic templates --------------------------------------------------- //
template <class T>
void zeros(                                                                   // Initialize an m-by-n matrix of zeros
    std::vector< std::vector < T > >            &,                            // (input/output) matrix
    int                                          ,                            // (input) number of rows
    int                                                                       // (input) number of columns
);

template <class T, size_t m, size_t n>
void zeros(                                                                   // Initialize an m-by-n matrix of zeros
    std::array< std::array < T, n >, m >        &                             // (input/output) matrix
);

template <class T>
void ones(                                                                    // Initialize an m-by-n matrix of ones
    std::vector< std::vector < T > >            &,                            // (input/output) matrix
    int                                          ,                            // (input) number of rows
    int                                                                       // (input) number of columns
);

template <class T, size_t m, size_t n>
void ones(                                                                    // Initialize an m-by-n matrix of ones
    std::array< std::array < T, n >, m >        &                             // (input/output) matrix
);

template <class T>
void eye(                                                                     // Initialize and m-by-n identity matrix
    std::vector< std::vector < T > >            &,                            // (input/output) matrix
    int                                          ,                            // (input) number of rows
    int                                                                       // (input) number of columns
);

template <class T, size_t m, size_t n>
void eye(                                                                     // Initialize and m-by-n identity matrix
    std::array< std::array < T, n >, m >        &                             // (input/output) matrix
);

// Matrix multiplications --------------------------------------------------- //
template <class T>
void matmul(                                                                  // Matrix multiplication with a constant
    T                            ,                                            // (input) 1st input scalar
    std::vector< std::vector< T > >             &,                            // (input) 2nd input matrix
    std::vector< std::vector< T > >             &                             // (input/output) output matrix
);

template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a constant
    T                            ,                                            // (input) 1st input scalar
    std::array< std::array< T, n >, m >         &,                            // (input) 2nd input matrix
    std::array< std::array< T, n >, m >         &                             // (input/output) output matrix
);

template <class T>
void matmul(                                                                  // Matrix multiplication with a constant
    std::vector< std::vector< T > >             &,                            // (input) 1st input matrix
    T                                            ,                            // (input) 2nd input scalar
    std::vector< std::vector< T > >             &                             // (input/output) output matrix
);

template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a constant
    std::array< std::array< T, n >, m >         &,                            // (input) 1st input matrix
    T                                            ,                            // (input) 2nd input scalar
    std::array< std::array< T, n >, m >         &                             // (input/output) output matrix
);

template <class T>
void matmul(                                                                  // Matrix multiplication with a vector
    std::vector< T >                            &,                            // (input) 1st input vector
    std::vector< std::vector< T > >             &,                            // (input) 2nd input matrix
    std::vector< T >                            &                             // (input/output) output vector
);

template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a vector
    std::array< T, m >                          &,                            // (input) 1st input vector
    std::array< std::array < T, n >, m >        &,                            // (input) 2nd input matrix
    std::array< T, n >                          &                             // (input/output) output vector
);

template <class T>
void matmul(                                                                  // Matrix multiplication with a vector
    std::vector< std::vector< T > >             &,                            // (input) 1st input matrix
    std::vector< T >                            &,                            // (input) 2nd input vector
    std::vector< T >                            &                             // (input/output) output vector
);

template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a vector
    std::array< std::array < T, n >, m >        &,                            // (input) 1st input matrix
    std::array< T, n >                          &,                            // (input) 2nd input vector
    std::array< T, m >                          &                             // (input/output) output vector
);

template <class T>
void matmul(                                                                  // Matrix multiplication
    std::vector< std::vector< T > >             &,                            // (input) 1st input matrix
    std::vector< std::vector< T > >             &,                            // (input) 2nd input matrix
    std::vector< std::vector< T > >             &                             // (input/output) output matrix
);

template <class T, size_t m, size_t n, size_t l> 
void matmul(                                                                  // Matrix multiplication
    std::array< std::array< T, n >, m >         &,                            // (input) 1st input matrix
    std::array< std::array< T, l >, n >         &,                            // (input) 2nd input matrix
    std::array< std::array< T, l >, m >         &                             // (input/output) output matrix
);

template <class T, size_t d1, size_t d2, size_t d3>
std::array< std::array<T, d2> , d1> matmul(
    const std::array< std::array<T, d3>, d1>    &, 
    const std::array< std::array<T, d2>, d3>    &
);

template <class T>
std::vector< std::vector<T> > matmul(
    const std::vector< std::vector<T> >         &, 
    const std::vector< std::vector<T> >         &
);

template <class T>
std::vector< std::vector<T> > matmulDiag(
    const std::vector<T>                        &, 
    const std::vector< std::vector<T> >         &
);

template <class T>
std::vector< std::vector<T> > matmulDiag(
    const std::vector< std::vector<T> >         &, 
    const std::vector<T>                        &
);

template <class T, size_t d1, size_t d2>
std::array< std::array<T, d2> , d1> matmulDiag(
    const std::array<T, d1>                     &, 
    const std::array< std::array<T, d2>, d1>    &
);

template <class T, size_t d1, size_t d2>
std::array< std::array<T, d2> , d1> matmulDiag(
    const std::array< std::array<T, d2>, d1>    & , 
    const std::array<T, d2>                     &
);

// Matrix Vector Multiplication ------------------------------------------------------ //
template <class T>
std::vector<T> matmul(
    const std::vector< std::vector<T>>          &, 
    const std::vector<T>                        &
);

template <class T, size_t d1, size_t d2>
std::array<T, d1> matmul(
    const std::array< std::array<T, d2>, d1>    &, 
    const std::array<T, d2>                     &
);

// Operator Tensor_Product ------------------------------------------------------------ //
template <class T>
std::vector<std::vector<T>> tensorProduct(
    std::vector<T> const                        &,
    std::vector<T> const                        &
);

template <class T, size_t n, size_t m>
std::array<std::array<T,m>, n> tensorProduct(
    std::array<T, n> const                      &,
    std::array<T, m> const                      &
);

// Matrix manipulations ----------------------------------------------------- //
template <class T>
void transpose(                                                               // Matrix transposition
    std::vector< std::vector< T > >             &,                            // (input) Input matrix
    std::vector< std::vector< T > >             &                             // (input/output) Output transpose
);

template <class T, size_t m, size_t n>
void transpose(                                                               // Matrix transposition
    std::array< std::array< T, n >, m >         &,                            // (input) Input matrix
    std::array< std::array< T, m >, n >         &                             // (input/output) Output transpose
);

template <class T>
std::vector< std::vector< T > >  transpose(                                   // Matrix transposition
    const std::vector< std::vector< T > >       &                             // (input) Input matrix
);

template <class T, size_t m, size_t n>
std::array< std::array< T, m >, n > transpose(                                // Matrix transposition
    const std::array< std::array< T, n >, m >   &                             // (input) Input matrix
);

template <class T>
void complement(                                                              // Complement extraction
    int                                          ,                            // (input) complement 1st index
    int                                          ,                            // (input) complement 2nd index
    std::vector< std::vector< T > >             &,                            // (input) input matrix
    std::vector< std::vector< T > >             &                             // (input/output) (i,j) complement
);

template <class T, size_t m, size_t n>
void complement(                                                              // Complement extraction
    int                                          ,                            // (input) complement 1st index
    int                                          ,                            // (input) complement 2nd index
    std::array< std::array< T, n >, m >         &,                            // (input) input matrix
    std::array< std::array<T, n-1>, m-1>        &                             // (input/output) (i,j) complement
);

template <class T>
void triL(                                                                    // Lower triangular part extraction
    std::vector< std::vector< T > >             &,                            // (input) Input matrix
    std::vector< std::vector< T > >             &                             // (input/output) Lower triangular part
);

template <class T, size_t m, size_t n>
void triL(                                                                    // Extract the lower triangular part
    std::array< std::array< T, n >, m >         &,                            // (input) Input matrix
    std::array< std::array< T, n >, m >         &                             // (input/output) Lower triangular part
);

template <class T>
void triU(                                                                    // Extract the upper triangular part
    std::vector< std::vector< T > >             &,                            // (input) input matrix
    std::vector< std::vector< T > >             &                             // (input/output) Upper triangular part
);

template <class T, size_t m, size_t n>
void triU(                                                                    // Extract the upper triangular part
    std::array< std::array< T, n >, m >         &,                            // (input) input matrix
    std::array< std::array< T, n >, m >         &                             // (input/output) Upper triangular part
);

// Matrix determinant ------------------------------------------------------- //
template <class T>
T det(                                                                        // Matrix determinant
    std::vector< std::vector< T > >             &                             // (input) input matrix
);

template <class T>
T det(                                                                        // Dummy function for recursive calls to det()
    std::array< std::array < T, 1 >, 1 >        &                             // (input) input matrix
);

template <class T, size_t m, size_t n>
T det(                                                                        // Matrix determinant
    std::array< std::array < T, n >, m >        &                             // (input) input matrix
);

// Linear system ------------------------------------------------------------ //
template <class T>
void cramer(                                                                  // Solve linear system using Cramer's rule
    std::vector< std::vector < T > >            &,                            // (input) coeff. matrix
    std::vector< T >                            &,                            // (input) source term
    std::vector< T >                            &                             // (input/output) solution
);

template <class T, size_t m, size_t n>
void cramer(                                                                  // Solve linear system using Cramer's rule
    std::array< std::array < T, n >, m >        &,                            // (input) coeff. matrix
    std::array< T, m >                          &,                            // (input) source term
    std::array< T, n >                          &                             // (input/output) solution
);

unsigned int factorizeLU(                                                              // LU decomposition of matrix
    std::vector<std::vector<double>>            &,                            // (input) Input matrix
    std::vector<std::vector<double>>            &,                            // (input/output) L factor
    std::vector<std::vector<double>>            &,                            // (input/output) U factor
    std::vector<std::vector<double>>            *                             // (input/optional) permutation matrix
);

template<size_t m>
unsigned int factorizeLU(                                                              // LU decomposition of matrix
    std::array<std::array<double, m>, m>        &,                            // (input) Input matrix
    std::array<std::array<double, m>, m>        &,                            // (input/output) L factor
    std::array<std::array<double, m>, m>        &,                            // (input/output) U factor
    std::array<std::array<double, m>, m>        *                             // (input/optional) permutation matrix
);

void backwardSubstitution(                                                           // Backward substitution algorithm
    std::vector<std::vector<double>>            &,                            // (input) Coeffs. matrix
    std::vector<double>                         &,                            // (input) Source term
    std::vector<double>                         &                             // (input/output) Solution
);

template<size_t m>
void backwardSubstitution(                                                           // Backward substitution algorithm
    std::array<std::array<double, m>, m>        &,                            // (input) Coeffs. matrix 
    std::array<double, m>                       &,                            // (input) Source term
    std::array<double, m>                       &                             // (input/output) Solution
);

void forwardSubstitution(                                                            // Forward substitution algorithm
    std::vector<std::vector<double>>            &,                            // (input) Coeffs. matrix
    std::vector<double>                         &,                            // (input) Source term
    std::vector<double>                         &                             // (input/output) Solution
);

template<size_t m>
void forwardSubstitution(                                                            // Forward substitution algorithm
    std::array<std::array<double, m>, m>        &,                            // (input) Coeffs. matrix
    std::array<double, m>                       &,                            // (input) Source term
    std::array<double, m>                       &                             // (input/output) Solution
);

void solveLU(                                                                 // Solve linear system using LU factorization
    std::vector<std::vector<double>>            &,                            // (input) Input coeffs. matrix
    std::vector<double>                         &,                            // (input) Source term
    std::vector<double>                         &                             // (input/output) Solution
);

template<size_t m>
void solveLU(                                                                 // Solve linear system using LU factorization
    std::array<std::array<double, m>, m>        &,                            // (input) Input coeffs. matrix
    std::array<double, m>                       &,                            // (input) Source term
    std::array<double, m>                       &                             // (input/output) Solution
);


}

}

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "LinearAlgebra.tpp"
# include "Manipulation.tpp"
# include "Multiplication.tpp"
# include "system_solvers_small.tpp"

# endif
