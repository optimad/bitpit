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
# ifndef __BITPIT_MATRIX_UTILITIES_HPP__
# define __BITPIT_MATRIX_UTILITIES_HPP__

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

// Auxiliary routines --------------------------------------------------------- //
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

}

}

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "matrix_utilities.tpp"

# endif
