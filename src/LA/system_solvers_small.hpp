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
# ifndef __BITPIT_SYSTEM_SOLVERS_SMALL_HPP__
# define __BITPIT_SYSTEM_SOLVERS_SMALL_HPP__

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
# include "system_solvers_small.tpp"

# endif
