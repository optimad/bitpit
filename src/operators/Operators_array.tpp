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

// ===================================================================================//
//                                  OPERATORS                                         //
//                                                                                    //
// Operators for standard template library arrays.                                    //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
/*!
        \file: Operators_array.tpp
        \brief: Basic operator for (C++11) std::array
        \author: Alessandro Alaia
*/
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// ================================================================================== //
// TEMPLATES IMPLEMENTATION                                                           //
// ================================================================================== //

/*!
   \ingroup Operators
   \{
 */

// Operator "+" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::array

    Perform the element-wise sum between two arrays (x and y) and returns a
    array z s.t.
    z[i] = x[i] + y[i] for all i = 0, ..., d
    where d = x.size() = y.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator+ calls itself
    to perform the element-wise sum between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise sum of x and y.
*/
template <class T, size_t d>
std::array<T, d> operator+ (
  const std::array<T, d>                        &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] + y[i];
};

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::array.

    Perform the element-wise sum between a array (x) and a constant (y), 
    and returns a array (z) s.t.
    z[i] = x[i] + y for all i = 0, ..., d
    where d = x.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator+ calls itself to perform the element-wise sum
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise sum of x and y.
*/
template <class T, size_t d>
std::array<T, d> operator+ (
  const std::array<T, d>                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] + y;
};

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::array.

    Perform the element-wise sum between a constant (x), and a array (y)
    and returns a array (z) s.t.
    z[i] = x + y[i] for all i = 0, ..., d
    where d = y.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if the i-th element of y and
    x are std::array, operator+ calls itself to perform the element-wise sum
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise sum of x and y.
*/
template <class T, size_t d>
std::array<T, d> operator+ (
  const T                                       &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>   z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
z = y + x;

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::array.

    Perform the element-wise sum between a array of arrays (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] + y for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = x.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if x[i][j] and
    y are std::array, operator+ calls itself to perform the element-wise sum
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of x and
    storing the element-wise sum of x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator+ (
  const T                                       &x,
  const std::array<std::array<T, e>, d>         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>   z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x + y[i];
};

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::array.

    Perform the element-wise sum between a constant (y) and a array of arrays (y),
    and returns z s.t.
    z[i][j] = x + y[i][j] for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if y[i][j] and
    y are std::array, operator+ calls itself to perform the element-wise sum
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of y and
    storing the element-wise sum of x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator+ (
  const std::array<std::array<T, e>, d>         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>   z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
z = y + x;

return (z); };

// Operator "+=" ==================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator+= for std::array.

    Increment each element in the input array, using the corresping value
    on array at the r.h.s. as increment, i.e.:
    x[i] += y[i] for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator+= calls itself
    to increment x[i][j] by y[i][j], j = 0, ..., x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator+= (
  std::array< T, d >                            &x,
  const std::array< T, d >                      &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    x[i] += y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+= for std::array.

    Increment each element in the input array, using the value
    on the r.h.s. as increment, i.e.:
    x[i] += y for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator+= calls itself
    to increment x[i][j] by y[j], j = 0, ... , x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator+= (
  std::array< T, d >                            &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    x[i] += y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+= for std::array.

    Increment each element in the input array, using the value
    on the r.h.s. as increment, i.e.:
    x[i][j] += y for all j = 0, ..., e-1, i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if x[i][j] and
    y are std::array, operator+= calls itself
    to increment x[i][j][k] by y[k], k = 0, ... , x.size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator+= (
    std::array< std::array< T, e >, d >         &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    x[i] += y;
} //next i

return (x); };

// Operator "-" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::array.

    Perform the element-wise difference between two arrays (x and y) and returns a
    array z s.t.
    z[i] = x[i] - y[i] for all i = 0, ..., d-1
    where d = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator- calls itself
    to perform the element-wise difference between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result array with n elements, storing the element-wise difference between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator- (
  const std::array<T, d>                        &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] - y[i];
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::array.

    Perform the element-wise difference between a array (x) and a constant (y), 
    and returns a array (z) s.t.
    z[i] = x[i] - y for all i = 0, ..., d-1
    where d = x.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator- calls itself to perform the element-wise difference
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise difference between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator- (
  const std::array<T, d>                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] - y;
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::array.

    Perform the element-wise difference between a constant (x), and a array (y)
    and returns a array (z) s.t.
    z[i] = x - y[i] for all i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if the i-th element of y and
    x are std::array, operator- calls itself to perform the element-wise difference
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise difference between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator- (
  const T                                       &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>   z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
z = y - x;

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::array.

    Perform the element-wise difference between a constant (y) and a array of arrays (y),
    and returns z s.t.
    z[i][j] = x - y[i][j] for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if y[i][j] and
    y are std::array, operator- calls itself to perform the element-wise difference
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of y and
    storing the element-wise difference between x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator- (
  const T                                       &x,
  const std::array<std::array<T, e>, d>         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>   z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = y[i] - x;
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::array.

    Perform the element-wise difference between a array of arrays (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] - y for all j = 0, ..., e-1, i = 0, ..., d-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if x[i][j] and
    y are std::array, operator. calls itself to perform the element-wise difference
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of x and
    storing the element-wise difference between x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator- (
  const std::array<std::array<T, e>, d>         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>   z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
z = y - x;

return(z); };

// Operator "-=" ==================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator-= for std::array.

    Decrement each element in the input array, using the corresping value
    on the array at the r.h.s. as negative increment, i.e.:
    x[i] -= y[i] for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator-= is defined.

    The element-wise decrement is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator-= calls itself
    to decrement x[i][j] by y[i][j], j = 0, ..., x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument decremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator-= (
  std::array< T, d >                            &x,
  const std::array< T, d >                      &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    x[i] -= y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator-= for std::array.

    Decrement each element in the input array, using the value
    on the r.h.s. as negative increment, i.e.:
    x[i] -= y for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator-= is defined.

    The element-wise decrement is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator-= calls itself
    to decrement x[i][j] by y[j], j = 0, ... , x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument decremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator-= (
  std::array< T, d >                            &x,
  const T                                       &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array< T, d >& operator-= (                                                        //
//   array< T, d >                 &x,                                                //
//   const T                       &y)                                                //
//                                                                                    //
// Element-wise decrement. Returns:                                                   //
//        x -= y, s.t x[i] -= y                                                       //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array< T, d >, 1st argument of '-=' operator                               //
// - y   : T          , 2nd argument of '-=' operator                                 //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - x   : array< T, d >&, reference to first argument                                //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    x[i] -= y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator-= for std::array.

    Decrement each element in the input array, using the value
    on the r.h.s. as negative increment, i.e.:
    x[i][j] -= y for all j = 0, ..., e-1, i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator-= is defined.

    The element-wise decrement is performed recursively, i.e. if x[i][j] and
    y are std::array, operator-= calls itself
    to decrement x[i][j][k] by y[k], k = 0, ... , x.size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument decremented with r.h.s. values.
*/
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator-= (
    std::array< std::array< T, e >, d >         &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    x[i] -= y;
} //next i

return (x); };

// Operator "*" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::array.

    Perform the element-wise product between two arrays (x and y) and returns a
    array z s.t.
    z[i] = x[i] * y[i] for all i = 0, ..., d-1
    where d = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator* calls itself
    to perform the element-wise product between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise product between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator* (
  const std::array<T, d>                        &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] * y[i];
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::array.

    Perform the element-wise product between a array (x) and a constant (y), 
    and returns a array (z) s.t.
    z[i] = x[i] * y for all i = 0, ..., d-1
    where d = x.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator* calls itself to perform the element-wise product
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise product between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator* (
  const std::array<T, d>                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] * y;
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::array.

    Perform the element-wise product between a constant (x), and a array (y)
    and returns a array (z) s.t.
    z[i] = x * y[i] for all i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if the i-th element of y and
    x are std::array, operator* calls itself to perform the element-wise product
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise product between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator* (
  const T                                       &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>   z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
z = y * x;

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::array.

    Perform the element-wise product between a constant (y) and a array of arrays (y),
    and returns z s.t.
    z[i][j] = x * y[i][j] for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if y[i][j] and
    y are std::array, operator* calls itself to perform the element-wise product
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of y and
    storing the element-wise product between x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator* (
  const T                                       &x,
  const std::array<std::array<T, e>, d>         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>   z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = y[i] * x;
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::array.

    Perform the element-wise product between a array of arrays (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] * y for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = x.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if x[i][j] and
    y are std::array, operator* calls itself to perform the element-wise product
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of x and
    storing the element-wise product between x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator* (
  const std::array<std::array<T, e>, d>         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>   z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = x[i] * y;
} //next i

return(z); };

// Operator "*=" ===================================================================== //
  
// ---------------------------------------------------------------------------------- //
/*!
    Operator*= for std::array.

    Increment each element in the input array, using the corresping value
    on array at the r.h.s. as increment, i.e.:
    x[i] *= y[i] for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator*= calls itself
    to increment x[i][j] by y[i][j], j = 0, ..., x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator*= (
  std::array< T, d >                            &x,
  const std::array< T, d >                      &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    x[i] *= y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator*= for std::array.

    Increment each element in the input array, using the value
    on the r.h.s. as increment, i.e.:
    x[i] *= y for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator*= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator*= calls itself
    to increment x[i][j] by y[j], j = 0, ... , x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator*= (
  std::array< T, d >                            &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    x[i] *= y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator*= for std::array.

    Increment each element in the input array, using the value
    on the r.h.s. as increment, i.e.:
    x[i][j] *= y for all j = 0, ..., e-1, i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator*= is defined.

    The element-wise increment is performed recursively, i.e. if x[i][j] and
    y are std::array, operator*= calls itself
    to increment x[i][j][k] by y[k], k = 0, ... , x.size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator*= (
    std::array< std::array< T, e >, d >         &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    x[i] *= y;
} //next i

return (x); };

// Operator "/" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::array.

    Perform the element-wise division between two arrays (x and y) and returns a
    array z s.t.
    z[i] = x[i] / y[i] for all i = 0, ..., d-1
    where d = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator/ calls itself
    to perform the element-wise division between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise division between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator/ (
  const std::array<T, d>                        &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// PERFORM DIVISION                                                                   //
// ================================================================================== //
for (size_t i = 0; i < d; i++){
    z[i] = x[i] / y[i];
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::array.

    Perform the element-wise division between a array (x) and a constant (y), 
    and returns a array (z) s.t.
    z[i] = x[i] / y for all i = 0, ..., d-1
    where d = x.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator/ calls itself to perform the element-wise division
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise division between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator/ (
  const std::array<T, d>                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
std::array<T, d>    z;

// Counters
size_t         i;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (i = 0; i < d; i++) {
    z[i] = x[i]/y;
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::array.

    Perform the element-wise division between a constant (x), and a array (y)
    and returns a array (z) s.t.
    z[i] = x / y[i] for all i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if the i-th element of y and
    x are std::array, operator/ calls itself to perform the element-wise division
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array with d elements, storing the element-wise division between x and y.
*/
template <class T, size_t d>
std::array<T, d> operator/ (
  const T                                       &x,
  const std::array<T, d>                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>      z;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = x/y[i];
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::array.

    Perform the element-wise division between a constant (y) and a array of arrays (y),
    and returns z s.t.
    z[i][j] = x / y[i][j] for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = y.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if y[i][j] and
    y are std::array, operator/ calls itself to perform the element-wise division
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of y and
    storing the element-wise division between x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator/ (
  const T                                       &x,
  const std::array<std::array<T, e>, d>         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>      z;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = x/y[i];
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::array.

    Perform the element-wise division between a array of arrays (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] / y for all j = 0, ..., e-1, i = 0, ..., d-1
    where d = x.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if x[i][j] and
    y are std::array, operator/ calls itself to perform the element-wise division
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result array of arrays having the same dimensions of x and
    storing the element-wise division between x and y.
*/
template <class T, size_t d, size_t e>
std::array<std::array<T, e>, d> operator/ (
  const std::array<std::array<T, e>, d>         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, e>, d>      z;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = x[i]/y;
}

return(z); };

// Operator "/=" ===================================================================== //
  
// ---------------------------------------------------------------------------------- //
/*!
    Operator/= for std::array.

    Increment each element in the input array, using the corresping value
    on array at the r.h.s. as increment, i.e.:
    x[i] /= y[i] for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator/= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::array, operator/= calls itself
    to increment x[i][j] by y[i][j], j = 0, ..., x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator/= (
  std::array< T, d >                            &x,
  const std::array< T, d >                      &y
) {

for (size_t i = 0; i < d; i++){
    x[i] /= y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/= for std::array.

    Increment each element in the input array, using the value
    on the r.h.s. as increment, i.e.:
    x[i] /= y for all i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator/= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    y are std::array, operator/= calls itself
    to increment x[i][j] by y[j], j = 0, ... , x[i].size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d>
std::array< T, d >& operator/= (
  std::array< T, d >                            &x,
  const T                                       &y
) {

for (size_t i = 0; i < d; i++){
    x[i] /= y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/= for std::array.

    Increment each element in the input array, using the value
    on the r.h.s. as increment, i.e.:
    x[i][j] /= y for all j = 0, ..., e-1, i = 0, ..., d-1
    Template parameters, T, can by any type for which the operator/= is defined.

    The element-wise increment is performed recursively, i.e. if x[i][j] and
    y are std::array, operator/= calls itself
    to increment x[i][j][k] by y[k], k = 0, ... , x.size()-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T, size_t d, size_t e>
std::array< std::array< T, e >, d >& operator/= (
    std::array< std::array< T, e >, d >         &x,
    const T                                     &y
) {

for (size_t i = 0; i < d; i++) {
    x[i] /= y;
} //next i

return (x); };

// Output operator ================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Insertion operator for std::array.

    Flush the content of std::array to std::ostream.
    The content of the input array is flushed with the following format:
    x[0] x[1] x[2] ... x[d-1] where d = x.size();
    (i.e. array elements are separated by blank spaces).
    Template parameter T can be any type such that operator<< is defined.

    \param[in,out] out output stream
    \param[in] x argument of insertion operator

    \result reference to the stream (allows concatenation)
*/
template <class T, size_t d>
std::ostream& operator<< (
    std::ostream                                &out,
    const std::array<T, d>                      &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// OUTPUT ARRAY CONTENT                                                               //
// ================================================================================== //
if (d == 0) {
    return(out);
}
for (size_t i = 0; i < d-1; i++) {
    out << x[i] << " ";
} //next i
out << x[d-1] ;

return(out); };

// ---------------------------------------------------------------------------------- //
/*!
    Insertion operator for std::array.

    Flush the content of std::array to std::ofstream.
    The content of the input array is flushed with the following format:
    x[0] x[1] x[2] ... x[d-1] where d = x.size();
    (i.e. array elements are separated by blank spaces).
    Template parameter T can be any type such that operator<< is defined.

    \param[in,out] out output file stream
    \param[in] x argument of insertion operator

    \result reference to the stream (allows concatenation)
*/
template <class T, size_t d>
std::ofstream& operator<< (
    std::ofstream                               &out,
    const std::array<T, d>                      &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// OUTPUT ARRAY CONTENT                                                               //
// ================================================================================== //
if (d == 0) {
    return(out);
}
for (size_t i = 0; i < d-1; i++) {
    out << x[i] << " ";
} //next i
out << x[d-1];

return(out); };

// Input operator =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Extraction operator for std::array.

    Extract the content of std::array from std::istream.
    The content of the input array is extracted until end-of-stream condition
    is met or all available position in the array are filled.

    \param[in,out] in input stream
    \param[in,out] x argument of extraction operator

    \result reference to input stream (allows concatenation)
*/
template <class T, size_t d>
std::istream& operator>> (
    std::istream                                &in,
    std::array<T, d>                            &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T                    dummy;

// Counters
size_t               i;

// ================================================================================== //
// EXTRACT STREAM CONTENT INTO ARRAY                                                  //
// ================================================================================== //
i = 0;
while ((in.good()) && (i < d)) {
    if (in >> dummy) { x[i] = dummy; }
    i++;
} //next i

return(in); };

// ---------------------------------------------------------------------------------- //
/*!
    Extraction operator for std::array.

    Extract the content of std::array from std::ifstream.
    The content of the input array is extracted until end-of-stream condition
    is met or all available position in the array are filled.

    \param[in,out] in input file stream
    \param[in,out] x argument of extraction operator

    \result reference to input file stream (allows concatenation)
*/
template <class T, size_t d>
std::ifstream& operator>> (
    std::ifstream                               &in,
    std::array<T, d>                            &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T       dummy;

// Counters
int     i;

// ================================================================================== //
// EXTRACT FILE CONTENT INTO ARRAY                                                    //
// ================================================================================== //
i = 0;
while ((in.good()) && (i < d)) {
    if (in >> dummy) { x[i] = dummy; }
    i++;
} //next i

return(in); };

// Miscellanea ====================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Display array in a nicely formatted to a std::ostream

    \param[in,out] out output stream
    \param[in] x array to be displayed
    \param[in] padding (default = 0) number of trailing spaces

    \result reference to output stream
*/
template<class T, size_t d>
std::ostream& display(
    std::ostream                                &out,
    const std::array<T, d>                      &x,
    unsigned int                                 padding
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Counters
typename std::array<T, d>::const_iterator      i, e = x.cend();

// ================================================================================== //
// DISPLAY VECTOR                                                                     //
// ================================================================================== //
if (x.size() == 0) {
    out << "[ ]";
    return(out);
}
out << std::string(padding, ' ') << "[ ";
--e;
for (i = x.begin(); i != e; ++i) {
    display(out, *i) << ", ";
} //next i
display(out, *e) << " ]";

return(out); }

// ---------------------------------------------------------------------------------- //
/*!
    Display array in a nicely formatted to a std::ofstream

    \param[in,out] out output file stream
    \param[in] x array to be displayed
    \param[in] padding (default = 0) number of trailing spaces

    \result reference to output stream
*/
template<class T, size_t d>
std::ofstream& display(
    std::ofstream                               &out,
    const std::array<T, d>                      &x,
    unsigned int                                 padding
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Counters
typename std::array<T, d>::const_iterator      i, e = x.cend();

// ================================================================================== //
// DISPLAY VECTOR                                                                     //
// ================================================================================== //
if (x.size() == 0) {
    out << "[ ]";
    return(out);
}
out << std::string(padding, ' ') << "[ ";
--e;
for (i = x.begin(); i != e; ++i) {
    display(out, *i) << ", ";
} //next i
display(out, *e) << " ]";

return(out); }

/*!
   \}
 */




