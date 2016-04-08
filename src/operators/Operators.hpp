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

// ================================================================================== //
//                                  OPERATORS                                         //
//                                                                                    //
// Operators for standard template library containers.                                //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //
#ifndef __BITPIT_OPERATORS_HPP__
#define __BITPIT_OPERATORS_HPP__

// ================================================================================== //
// INCLUDES                                                                           //
// ================================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <sstream>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <algorithm>
# include <functional> 

// ================================================================================== //
// FUNCTION PROTOTYPES                                                                //
// ================================================================================== //

// STL vectors ====================================================================== //

// Operator "+" --------------------------------------------------------------------- //
template <class T>
std::vector< T > operator+ (                                                          // ELEMENT-WISE SUM BETWEEN TWO VECTORS
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< T > operator+ (                                                          // ELEMENT-WISE SUM BETWEEN VECTOR AND CONSTANT
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T>
std::vector< T > operator+ (                                                          // ELEMENT-WISE SUM BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator+ (                                           // ELEMENT-WISE SUM BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< std::vector< T > >   &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator+ (                                           // ELEMENT-WISE SUM BETWEEN CONSTANT AND VECTOR
        const std::vector< std::vector< T > >   &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "+=" -------------------------------------------------------------------- //
template <class T>
std::vector< T >& operator+= (                                                        // ELEMENT-WISE INCREMENT OF A VECTOR
        std::vector< T >                        &,                                    // 1st argument (std::vector)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< T >& operator+= (                                                        // ELEMENT-WISE INCREMENT OF A VECTOR
        std::vector< T >                        &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T>
std::vector< std::vector< T > >& operator+= (                                         // ELEMENT-WISE INCREMENT OF A VECTOR
        std::vector< std::vector< T > >         &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "-" --------------------------------------------------------------------- //
template <class T>
std::vector< T > operator- (                                                          // ELEMENT-WISE DIFFERENCE BETWEEN TWO VECTORS
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< T > operator- (                                                          // ELEMENT-WISE DIFFERENCE BETWEEN VECTOR AND CONSTANT
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T>
std::vector< T > operator- (                                                          // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator- (                                           // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< std::vector< T > >   &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator- (                                           // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND VECTOR
        const std::vector< std::vector< T > >   &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "-=" -------------------------------------------------------------------- //
template <class T>
std::vector< T >& operator-= (                                                        // ELEMENT-WISE DECREMENT OF A VECTOR
        std::vector< T >                        &,                                    // 1st argument (std::vector)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< T >& operator-= (                                                        // ELEMENT-WISE DECREMENT OF A VECTOR
        std::vector< T >                        &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T>
std::vector< std::vector< T > >& operator-= (                                         // ELEMENT-WISE DECREMENT OF A VECTOR
        std::vector< std::vector< T > >         &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "*" --------------------------------------------------------------------- //
template <class T>
std::vector< T > operator* (                                                          // ELEMENT-WISE PRODUCT BETWEEN TWO VECTORS
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< T > operator* (                                                          // ELEMENT-WISE PRODUCT BETWEEN VECTOR AND CONSTANT
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T>
std::vector< T > operator* (                                                          // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector < T > > operator* (                                          // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< std::vector< T > >   &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator* (                                           // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND VECTOR
        const std::vector< std::vector< T > >   &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "*=" -------------------------------------------------------------------- //
template <class T>
std::vector<T>& operator*= (
        std::vector<T> &,
        const std::vector<T> &
        );
template <class T>
std::vector<T>& operator*= (
        std::vector<T> &,
        const T &
        );
template <class T>
std::vector< std::vector<T> >& operator*= (
        std::vector< std::vector<T> > &,
        const T &
        );

// Operator "/" --------------------------------------------------------------------- //
template <class T>
std::vector< T > operator/ (                                                          // ELEMENT-WISE DIVISION BETWEEN VECTORS
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< T > operator/ (                                                          // ELEMENT-WISE DIVISION BETWEEN VECTOR AND CONSTANT
        const std::vector< T >                  &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T>
std::vector< T > operator/ (                                                          // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< T >                  &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator/ (                                           // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND VECTOR
        const T                                 &,                                    // 1st argument (constant)
        const std::vector< std::vector< T > >   &                                     // 2nd argument (std::vector)
        );
template <class T>
std::vector< std::vector< T > > operator/ (                                           // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND VECTOR
        const std::vector< std::vector< T > >   &,                                    // 1st argument (std::vector)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "/=" -------------------------------------------------------------------- //
template <class T>
std::vector<T>& operator/= (
        std::vector<T> &,
        const std::vector<T> &
        );
template <class T>
std::vector<T>& operator/= (
        std::vector<T> &,
        const T &
        );
template <class T>
std::vector< std::vector<T> >& operator/= (
        std::vector< std::vector<T> > &,
        const T &
        );
// Output operator ------------------------------------------------------------------ //
template <class T>
std::ostream& operator<< (                                                            // INSERTION OPERATOR
        std::ostream                            &,                                    // (input) output stream
        const std::vector< T >                  &                                     // (input) std::vector to be streamed
        );
template <class T>
std::ofstream& operator<< (                                                           // INSERTION OPERATOR
        std::ofstream                           &,                                    // (input) output file stream
        const std::vector< T >                  &                                     // (input) std::vector to be streamed
        );

// Input operators ------------------------------------------------------------------ //
template <class T>
std::istream& operator>> (                                                            // EXTRACTION OPERATOR
        std::istream                            &,                                    // (input) input stream
        std::vector< T >                        &                                     // (input/output) std::vector to be streamed
        );
template <class T>
std::ifstream& operator>> (                                                           // EXTRACTION OPERATOR
        std::ifstream                           &,                                    // (input) input file stream
        std::vector< T >                        &                                     // (input/output) std::vector to be streamed
        );

// Operator "min" ------------------------------------------------------------------- //
template <class T>
std::vector< T > min(                                                                 // RETURNS THE MINIMUM BETWEEN VECTOR X AND VECTOR Y
        const std::vector< T >                  &,                                    // (input) 1st argument of comparison (std::vector)
        const std::vector< T >                  &                                     // (input) 2nd argument of comparison (std::vector)
        );
template <class T>
std::vector< T > min(                                                                 // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const std::vector< T >                  &,                                    // (input) 1st argument of comparison (std::vector)
        const T                                 &                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T>
std::vector< T > min(                                                                 // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const T                                 &,                                    // (input) 1st argument of comparison (scalar)
        const std::vector< T >                  &                                     // (input) 2nd argument of comparison (std::vector)
        );
template <class T>
std::vector< std::vector< T > > min(                                                  // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const std::vector< std::vector< T > >   &,                                    // (input) 1st argument of comparison (std::vector)
        const T                                 &                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T>
std::vector< std::vector< T > > min(                                                  // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const T                                 &,                                    // (input) 1st argument of comparison (scalar)
        const std::vector< std::vector< T > >   &                                     // (input) 2nd argument of comparison (std::vector)
        );

// Operator "minval" ---------------------------------------------------------------- //
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
void inline minval(                                                                   // DUMMY ROUTINE FOR MINVAL SEARCH
        const T                                 &,                                    // (dummy input) 1st argument
        T                                       &                                     // (dummy input/output) 2nd argument
        );
template <class T, class T1>
void minval(                                                                          // RETURNS THE MINIMUM ELEMENT OF A VECTOR
        const std::vector<T>                    &,                                    // (input) input std::vector
        T1                                      &                                     // (input/output) minimum element
        );

// Operator "max" ------------------------------------------------------------------- //
template <class T>
std::vector< T > max(                                                                 // RETURNS THE MAXIMUM BETWEEN X AND Y
        const std::vector< T >                  &,                                    // (input) 1st argument of comparison
        const std::vector< T >                  &                                     // (input) 2nd argument of comparison
        );
template <class T>
std::vector< T > max(                                                                 // RETURNS THE MAXIMUM BETWEEN X AND Y
        const std::vector< T >                  &,                                    // (input) 1st argument of comparison
        const T                                 &                                     // (input) 2nd argument of comparison
        );
template <class T>
std::vector< T > max(                                                                 // RETURNS THE MAXIMUM BETWEEN X AND Y
        const T                                 &,                                    // (input) 1st argument of comparison
        const std::vector< T >                  &                                     // (input) 2nd argument of comparison
        );
template <class T>
std::vector< std::vector< T > > max(                                                  // RETURNS THE MAXIMUM BETWEEN X AND Y
        const std::vector< std::vector< T > >   &,                                    // (input) 1st argument of comparison
        const T                                 &                                     // (input) 2nd argument of comparison
        );
template <class T>
std::vector< std::vector< T > > max(                                                  // RETURNS THE MAXIMUM BETWEEN X AND Y
        const T                                 &,                                    // (input) 1st argument of comparison
        const std::vector< std::vector< T > >   &                                     // (input) 2nd argument of comparison
        );

// Operator "maxval" ---------------------------------------------------------------- //
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
void inline maxval(                                                                   // DUMMY ROUTINE FOR MAXVAL SEARCH
        const T                                 &,                                    // (dummy input) 1st argument
        T                                       &                                     // (dummy input/output) 2nd argument
        );
template <class T, class T1>
void maxval(                                                                          // RETURNS THE MAXIMUM ELEMENT OF A VECTOR
        const std::vector<T>                    &,                                    // (input) input std::vector
        T1                                      &                                     // (input/output) maximum element
        );

// Operator "sum" ------------------------------------------------------------------- //
template <class T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
void inline sum(                                                                      // DUMMY ROUTINE FOR SUM OPERATOR
        const T                                 &,                                    // (dummy input) 1st argument
        T                                       &                                     // (dummy input/output) 2nd argument
        );
template <class T, class T1>
void sum(                                                                             // RETURNS THE SUM OF ELEMENT OF A VECTOR
        const std::vector< T >                  &,                                    // (input) input std::vector
        T1                                      &                                     // (input/output) sum of std::vector's element
        );

// Operator "abs" ------------------------------------------------------------------- //
template <class T>
std::vector<T> abs(                                                                   // RETURNS THE ABSOLUTE VALUE OF A VECTOR
        const std::vector< T >                  &                                     // (input) input std::vector
        );

// Operator "sign" ------------------------------------------------------------------ //
template <class T>
T sign(                                                                               // SIGN FUNCTION
        const  T                                &                                     // (input) input value
      );

// Operator "pow" ------------------------------------------------------------------- //
template <class T>
std::vector< T > pow(                                                                 // RETURNS THE ELEMENTWISE POWER OF VECTOR
        std::vector< T >                        &,                                    // (input) input std::vector
        double                                                                        // (input) power
        );

// Operator "norm" ------------------------------------------------------------------ //
template <class T>
double norm1(                                                                        // RETURNS THE 1-NORM OF A VECTOR
        const std::vector< T >                  &                                     // (input) input std::vector
        );
template <class T>
double norm2(                                                                        // RETURNS THE 2-NORM OF A VECTOR
        const std::vector< T >                  &                                     // (input) input std::vector
        );
template <class T>
double norm(                                                                          // RETURNS THE P-NORM OF A VECTOR
        const std::vector< T >                  &,                                    // (input) input std::vector
        int                                                                           // (input) norm index
        );
template <class T>
double normInf(                                                                      // RETURNS THE inf-NORM OF A VECTOR
        const std::vector< T >                  &                                     // (input) input std::vector
        );

// Operator "dotProduct" ----------------------------------------------------------- //
template <class T>
T dotProduct(                                                                        // COMPUTE THE DOT PRODUCT OF TWO VECTORS
        const std::vector< T >                  &,                                    // (input) 1st argument of dot product
        const std::vector< T >                  &                                     // (input) 2nd argument of dot product
        );

// Operator crossProduct ----------------------------------------------------------- //
template <class T>
std::vector<T> crossProduct(                                                         // COMPUTE THE CROSS-PRODUCT OF TWO VECTORS
        std::vector<T> const                    &,                                    // (input) 1st argument of cross product
        std::vector<T> const                    &                                     // (input) 2nd argument of cross product
        );

// Various -------------------------------------------------------------------------- //
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
std::ostream& display(                                                                // DUMMY ROUTINE FOR RECURSIVE TEMPLATED FUNCTION display
        std::ostream                            &,                                    // (input/output) output stream
        const T                                 &                                     // (input) std::vector to be displayed
        );
template<class T>
std::ostream& display(                                                                // DISPLAY VECTOR IN A NICELY FORMATTED FORM
        std::ostream                            &,                                    // (input/output) output stream
        const std::vector<T>                    &,                                    // (input) std::vector to be displayed
        unsigned int                             padding = 0                          // (input/optional) number of trailing spaces
        );
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
std::ofstream& display(                                                               // DUMMY ROUTINE FOR RECURSIVE TEMPLATED FUNCTION display
        std::ofstream                           &,                                    // (input/output) output stream
        const T                                 &                                     // (input) std::vector to be displayed
        );
template<class T>
std::ofstream& display(                                                               // DISPLAY VECTOR IN A NICELY FORMATTED FORM
        std::ofstream                           &,                                    // (input/output) output stream
        const std::vector<T>                    &,                                    // (input) std::vector to be displayed
        unsigned int                             padding = 0                          // (input/optional) number of trailing spaces
        );

// C++ v10.0 arrays ================================================================= //

// Operator "+" --------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> operator+ (                                                          // ELEMENT-WISE SUM BETWEEN ARRAYS
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array<T, d> operator+ (                                                          // ELEMENT-WISE SUM BETWEEN ARRAY AND CONSTANT
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d>
std::array<T, d> operator+ (                                                          // ELEMENT-WISE SUM BETWEEN CONSTANT AND ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e> 
std::array<std::array<T, e>, d> operator+ (                                           // ELEMENT-WISE SUM BETWEEN CONSTANT AND 2D ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<std::array<T, e>, d>   &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator+ (                                           // ELEMENT-WISE SUM BETWEEN 2D ARRAY AND CONSTANT
        const std::array<std::array<T, e>, d>   &,                                    // 1st argument (constant)
        const T                                 &                                     // 2nd argument (std::array)
        );

// Operator "+=" -------------------------------------------------------------------- //
template <class T, size_t d>
std::array< T, d >& operator+= (                                                      // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const std::array< T, d >                &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array< T, d >& operator+= (                                                      // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator+= (                                     // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< std::array< T, e >, d >     &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "-" --------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> operator- (                                                          // ELEMENT-WISE DIFFERENCE BETWEEN TWO ARRAYS
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array<T, d> operator- (                                                          // ELEMENT-WISE DIFFERENCE BETWEEN ARRAY AND CONSTANT
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d>
std::array<T, d> operator- (                                                          // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator- (                                           // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND 2D ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<std::array<T, e>, d>   &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator- (                                           // ELEMENT-WISE DIFFERENCE BETWEEN 2D ARRAY AND CONSTANT
        const std::array<std::array<T, e>, d>   &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "-=" -------------------------------------------------------------------- //
template <class T, size_t d>
std::array< T, d >& operator-= (                                                      // ELEMENT-WISE DECREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const std::array< T, d >                &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array< T, d >& operator-= (                                                      // ELEMENT-WISE DECREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator-= (                                     // ELEMENT-WISE DECREMENT OF A ARRAY
        std::array< std::array< T, e >, d >     &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "*" --------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> operator* (                                                          // ELEMENT-WISE PRODUCT BETWEEN TWO ARRAYS
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array<T, d> operator* (                                                          // ELEMENT-WISE PRODUCT BETWEEN ARRAY AND CONSTANT
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d>
std::array<T, d> operator* (                                                          // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator* (                                           // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND 2D ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<std::array<T, e>, d>   &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator* (                                           // ELEMENT-WISE PRODUCT BETWEEN 2D ARRAY AND CONSTANT
        const std::array<std::array<T, e>, d>   &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "*=" -------------------------------------------------------------------- //
template <class T, size_t d>
std::array< T, d >& operator*= (                                                      // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const std::array< T, d >                &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array< T, d >& operator*= (                                                      // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator*= (                                     // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< std::array< T, e >, d >     &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "/" --------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> operator/ (                                                          // ELEMENT-WISE DIVISION BETWEEN ARRAYS
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array<T, d> operator/ (                                                          // ELEMENT-WISE DIVISION BETWEEN ARRAY AND CONSTANT
        const std::array<T, d>                  &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d>
std::array<T, d> operator/ (                                                          // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<T, d>                  &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator/ (                                           // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND 2D ARRAY
        const T                                 &,                                    // 1st argument (constant)
        const std::array<std::array<T, e>, d>   &                                     // 2nd argument (std::array)
        );
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator/ (                                           // ELEMENT-WISE PRODUCT BETWEEN 2D ARRAY AND CONSTANT
        const std::array<std::array<T, e>, d>   &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Operator "/=" -------------------------------------------------------------------- //
template <class T, size_t d>
std::array< T, d >& operator/= (                                                      // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const std::array< T, d >                &                                     // 2nd argument (std::array)
        );
template <class T, size_t d>
std::array< T, d >& operator/= (                                                      // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< T, d >                      &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator/= (                                     // ELEMENT-WISE INCREMENT OF A ARRAY
        std::array< std::array< T, e >, d >     &,                                    // 1st argument (std::array)
        const T                                 &                                     // 2nd argument (constant)
        );

// Output operator ------------------------------------------------------------------ //
template <class T, size_t d>
std::ostream& operator<< (                                                            // INSERTION OPERATOR
        std::ostream                            &,                                    // (input) output stream
        const std::array<T, d>                  &                                     // (input) std::array to be streamed
        );
template <class T, size_t d>
std::ofstream& operator<< (                                                           // INSERTION OPERATOR
        std::ofstream                           &,                                    // (input) output file stream
        const std::array<T, d>                  &                                     // (input) std::array to be streamed
        );

// Input operators ------------------------------------------------------------------ //
template <class T, size_t d>
std::istream& operator>> (                                                            // EXTRACTION OPERATOR
        std::istream                            &,                                    // (input) input stream
        std::array< T, d >                      &                                     // (input/output) std::array to be streamed
        );
template <class T, size_t d>
std::ifstream& operator>> (                                                           // EXTRACTION OPERATOR
        std::ifstream                           &,                                    // (input) input file stream
        std::array<T, d>                        &                                     // (input) std::array to be streamed
        );

// Operator "min" ------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> min(                                                                 // RETURNS THE MINIMUM BETWEEN ARRAY X AND ARRAY Y
        const std::array<T, d>                  &,                                    // (input) 1st argument of comparison (std::array)
        const std::array<T, d>                  &                                     // (input) 2nd argument of comparison (std::array)
        );
template <class T, size_t d>
std::array<T, d> min(                                                                 // RETURNS THE MINIMUM BETWEEN ARRAY X AND SCALAR Y
        const std::array<T, d>                  &,                                    // (input) 1st argument of comparison (std::array)
        const T                                 &                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d>
std::array<T, d> min(                                                                 // RETURNS THE MINIMUM BETWEEN SCALAR X AND SCALAR Y
        const T                                 &,                                    // (input) 1st argument of comparison (scalar)
        const std::array<T, d>                  &                                     // (input) 2nd argument of comparison (std::array)
        );
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> min(                                                  // RETURNS THE MINIMUM BETWEEN 2D ARRAY X AND SCALAR Y
        const std::array<std::array<T, n>, d>   &,                                    // (input) 1st argument of comparison (std::array)
        const T                                 &                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> min(                                                  // RETURNS THE MINIMUM BETWEEN SCALAR X AND 2D ARRAY Y
        const T                                 &,                                    // (input) 1st argument of comparison (scalar)
        const std::array<std::array<T, n>, d>   &                                     // (input) 2nd argument of comparison (std::array)
        );

// Operator "minval" ---------------------------------------------------------------- //
template <class T, size_t d, class T1>
void minval(                                                                          // RETURNS THE MINIMUM ELEMENT OF A ARRAY
        const std::array<T, d>                  &,                                    // (input) input std::array
        T1                                      &                                     // (input/output) minimum element
        );

// Operator "max" ------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> max(                                                                 // RETURNS THE MAXIMUM BETWEEN TWO ARRAYS
        const std::array<T, d>                  &,                                    // (input) 1st argument of comparison (std::array)
        const std::array<T, d>                  &                                     // (input) 2nd argument of comparison (std::array)
        );
template <class T, size_t d>
std::array<T, d> max(                                                                 // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const std::array<T, d>                  &,                                    // (input) 1st argument of comparison (std::array)
        const T                                 &                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d>
std::array<T, d> max(                                                                 // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const T                                 &,                                    // (input) 1st argument of comparison (scalar)
        const std::array<T, d>                  &                                     // (input) 2nd argument of comparison (std::array)
        );
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> max(                                                  // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const std::array<std::array<T, n>, d>   &,                                    // (input) 1st argument of comparison (std::array)
        const T                                 &                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> max(                                                  // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const T                                 &,                                    // (input) 1st argument of comparison (scalar)
        const std::array<std::array<T, n>, d>   &                                     // (input) 2nd argument of comparison (std::array)
        );

// Operator "maxval" ---------------------------------------------------------------- //
template <class T, size_t d, class T1>
void maxval(                                                                          // RETURNS THE MAXIMUM ELEMENT OF A ARRAY
        const std::array<T, d>                  &,                                    // (input) input std::array
        T1                                      &                                     // (input/output) maximum element
        );

// Operator "sum" ------------------------------------------------------------------- //
template <class T, size_t d, class T1>
void sum(                                                                             // RETURNS THE SUM OF ARRAY ELEMENTS
        const std::array<T, d>                  &,                                    // (input) input std::array
        T1                                      &                                     // (input/output) sum of std::array's elements
        );

// Operator "abs" ------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> abs(                                                                 // RETURNS THE ABSOLUTE VALUE OF A ARRAY
        const std::array<T, d>                  &                                     // (input) input std::array
        );

// Operator "pow" ------------------------------------------------------------------- //
template <class T, size_t d>
std::array<T, d> pow(                                                                 // RETURNS THE ELEMENTWISE POWER OF ARRAY
        std::array<T, d>                        &,                                    // (input) input std::array
        double                                                                        // (input) power
        );

// Operator "norm" ------------------------------------------------------------------ //
template <class T, size_t d>
double norm1(                                                                        // RETURNS THE 1-NORM OF ARRAY
        const std::array<T, d>                  &                                     // (input) input std::array
        );
template <class T, size_t d>
double norm2(                                                                        // RETURNS THE 2-NORM OF ARRAY
        const std::array<T, d>                  &                                     // (input) input std::array
        );
template <class T, size_t d>
double norm(                                                                          // RETURNS THE p-NORM OF ARRAY
        const std::array<T, d>                  &,                                    // (input) input std::array
        int                                                                           // (input) norm index
        );
template <class T, size_t d>
double normInf(                                                                      // RETURNS THE inf NORM OF ARRAY
        const std::array<T, d>                  &                                     // (input) input std::array
        );

// Operator "dotProduct" ----------------------------------------------------------- //
template <class T, size_t d>
T dotProduct(                                                                        // COMPUTE THE DOT PRODUCT OF TWO ARRAYS
        const std::array<T, d>                  &,                                    // (input) 1st argument of dot product
        const std::array<T, d>                  &                                     // (input) 2nd argument of dot product
        );

// Operator crossProduct ----------------------------------------------------------- //
template <class T, size_t d>
T crossProduct(                                                                      // COMPUTE THE CROSS-PRODUCT OF TWO ARRAYS
        std::array<T, 2> const                  &,                                    // (input) 1st argument of cross product
        std::array<T, 2> const                  &                                     // (input) 2nd argument of cross product
        );

template <class T, size_t d>
std::array<T, 3> crossProduct(                                                       // COMPUTE THE CROSS-PRODUCT OF TWO ARRAYS
        std::array<T, 3> const                  &,                                    // (input) 1st argument of cross product
        std::array<T, 3> const                  &                                     // (input) 2nd argument of cross product
        );

// Various -------------------------------------------------------------------------- //
template<class T, size_t d>
std::ostream& display(                                                                // DISPLAY ARRAY IN A NICELY FORMATTED FORM
        std::ostream                            &,                                    // (input/output) output stream
        const std::array<T, d>                  &,                                    // (input) std::array to be displayed
        unsigned int                             padding = 0                          // (input/optional) number of trailing spaces
        );
template<class T, size_t d>
std::ofstream& display(                                                               // DISPLAY ARRAY IN A NICELY FORMATTED FORM
        std::ofstream                           &,                                    // (input/output) output stream
        const std::array<T, d>                  &,                                    // (input) std::array to be displayed
        unsigned int                             padding = 0                          // (input/optional) number of trailing spaces
        );


// ================================================================================== //
// TEMPLATES                                                                          //
// ================================================================================== //
# include "Operators.tpp"
# include "Operators_vector.tpp"
# include "MathOperators_vector.tpp"
# include "Operators_array.tpp"
# include "MathOperators_array.tpp"

#endif
