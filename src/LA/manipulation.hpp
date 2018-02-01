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
# ifndef __BITPIT_MANIPULATION_HPP__
# define __BITPIT_MANIPULATION_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard template library
# include <vector>
# include <cmath>
# include <array>
# include <iostream>

# include "Operators.hpp"

# include "matrix_utilities.hpp"

namespace bitpit{

/*!
 * @brief   collection of functions to create and solve small linear systems
 */
namespace linearalgebra{

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

}

}

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "manipulation.tpp"

# endif
