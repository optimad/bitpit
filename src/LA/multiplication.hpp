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
# ifndef __BITPIT_MULTIPLICATION_HPP__
# define __BITPIT_MULTIPLICATION_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard template library
# include <vector>
# include <cmath>
# include <array>
# include <iostream>

# include "bitpit_operators.hpp"

# include "matrix_utilities.hpp"

namespace bitpit{

/*!
 * @brief   collection of functions to create and solve small linear systems
 */
namespace linearalgebra{

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

}

}

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "multiplication.tpp"

# endif
