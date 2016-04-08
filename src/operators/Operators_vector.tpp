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
// Operators for standard template library vectors.                                   //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// ================================================================================== //
// TEMPLATE IMPLEMENTATION                                                            //
// ================================================================================== //

/*!
   \ingroup Operators
   \{
 */

// Operator "+" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::vector.

    Perform the element-wise sum between two vectors (x and y) and returns a
    vector z s.t.
    z[i] = x[i] + y[i] for all i = 0, ..., n-1
    where n = std::min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator+ calls itself
    to perform the element-wise sum between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise sum of x and y.
*/
template <class T>
std::vector< T > operator+ (
  const std::vector< T >                        &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t             n = std::min(x.size(), y.size());
std::vector< T >        z(n, T(0));

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] + y[i];
};

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::vector.

    Perform the element-wise sum between a vector (x) and a constant (y), 
    and returns a vector (z) s.t.
    z[i] = x[i] + y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator+ calls itself to perform the element-wise sum
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise sum of x and y.
*/
template <class T>
std::vector< T > operator+ (
  const std::vector< T >                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = x.size();
std::vector< T >         z(n, T(0));

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] + y;
};

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::vector.

    Perform the element-wise sum between a constant (x), and a vector (y)
    and returns a vector (z) s.t.
    z[i] = x + y[i] for all i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if the i-th element of y and
    x are std::vector, operator+ calls itself to perform the element-wise sum
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise sum of x and y.
*/
template <class T>
std::vector< T > operator+ (
  const T                                       &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< T >   z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
z = y + x;

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::vector.

    Perform the element-wise sum between a vector of vectors (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] + y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator+ calls itself to perform the element-wise sum
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of x and
    storing the element-wise sum of x and y.
*/
template <class T>
std::vector< std::vector< T > > operator+ (
    const std::vector< std::vector< T > >       &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = x.size();
std::vector< std::vector< T > >   z(n);

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    z[i] = x[i] + y;
} //next i

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+ for std::vector.

    Perform the element-wise sum between a constant (y) and a vector of vectors (y),
    and returns z s.t.
    z[i][j] = x + y[i][j] for all j = 0, ..., y[i].size()-1, i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator+ is defined.

    The element-wise sum is performed recursively, i.e. if y[i][j] and
    y are std::vector, operator+ calls itself to perform the element-wise sum
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of y and
    storing the element-wise sum of x and y.
*/
template <class T>
std::vector< std::vector< T > > operator+ (
    const T                                     &x,
    const std::vector< std::vector< T > >       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< std::vector< T > >   z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
z = y + x;

return (z); };

// Operator "+=" ==================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator+= for std::vector.

    Increment each element in the input vector, using the corresping value
    on vector at the r.h.s. as increment, i.e.:
    x[i] += y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator+= calls itself
    to increment x[i][j] by y[i][j], j = 0, ..., min(x[i].size(), y[i].size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator+= (
  std::vector< T >                              &x,
  const std::vector< T >                        &y
) {

    // ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = std::min(x.size(), y.size());

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] += y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+= for std::vector.

    Increment each element in the input vector, using the value
    on the r.h.s. as increment, i.e.:
    x[i] += y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator+= calls itself
    to increment x[i][j] by y[j], j = 0, ... , min(x[i].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator+= (
  std::vector< T >                              &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] += y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator+= for std::vector.

    Increment each element in the input vector, using the value
    on the r.h.s. as increment, i.e.:
    x[i][j] += y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator+= is defined.

    The element-wise increment is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator+= calls itself
    to increment x[i][j][k] by y[k], k = 0, ... , min(x[i][j].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< std::vector< T > >& operator+= (
    std::vector< std::vector< T > >             &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    x[i] += y;
} //next i

return (x); };

// Operator "-" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::vector.

    Perform the element-wise difference between two vectors (x and y) and returns a
    vector z s.t.
    z[i] = x[i] - y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator- calls itself
    to perform the element-wise difference between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise difference between x and y.
*/
template <class T>
std::vector< T > operator- (
  const std::vector< T >                        &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = std::min(x.size(), y.size());
std::vector< T >         z(n, T(0));

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] - y[i];
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::vector.

    Perform the element-wise difference between a vector (x) and a constant (y), 
    and returns a vector (z) s.t.
    z[i] = x[i] - y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator- calls itself to perform the element-wise difference
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise difference between x and y.
*/
template <class T>
std::vector< T > operator- (
  const std::vector< T >                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = x.size();
std::vector< T >         z(n, T(0));

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] - y;
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::vector.

    Perform the element-wise difference between a constant (x), and a vector (y)
    and returns a vector (z) s.t.
    z[i] = x - y[i] for all i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if the i-th element of y and
    x are std::vector, operator- calls itself to perform the element-wise difference
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise difference between x and y.
*/
template <class T>
std::vector< T > operator- (
  const T                                       &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< T >   z;

// =================================================================================== //
// PERFORM DIFFERENCE                                                                  //
// =================================================================================== //
z = y - x;

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::vector.

    Perform the element-wise difference between a constant (y) and a vector of vectors (y),
    and returns z s.t.
    z[i][j] = x - y[i][j] for all j = 0, ..., y[i].size()-1, i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if y[i][j] and
    y are std::vector, operator- calls itself to perform the element-wise difference
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of y and
    storing the element-wise difference between x and y.
*/
template <class T>
std::vector< std::vector< T > > operator- (
  const T                                       &x,
  const std::vector< std::vector< T > >         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = y.size();
std::vector< std::vector< T > >   z(n);

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    z[i] = x - y[i];
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator- for std::vector.

    Perform the element-wise difference between a vector of vectors (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] - y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator- is defined.

    The element-wise difference is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator. calls itself to perform the element-wise difference
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of x and
    storing the element-wise difference between x and y.
*/
template <class T>
std::vector< std::vector< T > > operator- (
  const std::vector< std::vector< T > >         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< std::vector< T > >   z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
z = y - x;

return(z); };

// Operator "-=" ==================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator-= for std::vector.

    Decrement each element in the input vector, using the corresping value
    on vector at the r.h.s. as negative increment, i.e.:
    x[i] -= y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator-= is defined.

    The element-wise decrement is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator-= calls itself
    to decrement x[i][j] by y[i][j], j = 0, ..., min(x[i].size(), y[i].size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument decremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator-= (
  std::vector< T >                              &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = std::min(x.size(), y.size());

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] -= y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator-= for std::vector.

    Decrement each element in the input vector, using the value
    on the r.h.s. as negative increment, i.e.:
    x[i] -= y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator-= is defined.

    The element-wise decrement is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator-= calls itself
    to decrement x[i][j] by y[j], j = 0, ... , min(x[i].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument decremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator-= (
  std::vector< T >                              &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] -= y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator-= for std::vector.

    Decrement each element in the input vector, using the value
    on the r.h.s. as negative increment, i.e.:
    x[i][j] -= y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator-= is defined.

    The element-wise decrement is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator-= calls itself
    to decrement x[i][j][k] by y[k], k = 0, ... , min(x[i][j].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument decremented with r.h.s. values.
*/
template <class T>
std::vector< std::vector< T > >& operator-= (
    std::vector< std::vector< T > >             &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    x[i] -= y;
} //next i

return (x); };

// Operator "*" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::vector.

    Perform the element-wise product between two vectors (x and y) and returns a
    vector z s.t.
    z[i] = x[i] * y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator* calls itself
    to perform the element-wise product between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise product between x and y.
*/
template <class T>
std::vector< T > operator* (
  const std::vector< T >                        &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t             n = std::min(x.size(), y.size());
std::vector<T>          z(n, T(0));

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] * y[i];
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::vector.

    Perform the element-wise product between a vector (x) and a constant (y), 
    and returns a vector (z) s.t.
    z[i] = x[i] * y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator* calls itself to perform the element-wise product
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise product between x and y.
*/
template <class T>
std::vector< T > operator* (
  const std::vector< T >                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t             n = x.size();
std::vector< T >        z(n, T(0));

// =================================================================================== //
// PERFORM PRODUCT                                                                     //
// =================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] * y;
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::vector.

    Perform the element-wise product between a constant (x), and a vector (y)
    and returns a vector (z) s.t.
    z[i] = x * y[i] for all i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if the i-th element of y and
    x are std::vector, operator* calls itself to perform the element-wise product
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise product between x and y.
*/
template <class T>
std::vector< T > operator* (
  const T                                       &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< T >   z;                                                                     

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
z = y * x;

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::vector.

    Perform the element-wise product between a constant (y) and a vector of vectors (y),
    and returns z s.t.
    z[i][j] = x * y[i][j] for all j = 0, ..., y[i].size()-1, i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if y[i][j] and
    y are std::vector, operator* calls itself to perform the element-wise product
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of y and
    storing the element-wise product between x and y.
*/
template <class T>
std::vector< std::vector< T > > operator* (
  const T                                       &x,
  const std::vector< std::vector< T > >         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = y.size();
std::vector< std::vector< T > >   z(n);

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    z[i] = x * y[i];
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator* for std::vector.

    Perform the element-wise product between a vector of vectors (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] * y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator* is defined.

    The element-wise product is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator* calls itself to perform the element-wise product
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of x and
    storing the element-wise product between x and y.
*/
template <class T>
std::vector< std::vector< T > > operator* (
  const std::vector< std::vector< T > >         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< std::vector< T > >   z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
z = y*x;

return(z); };

// Operator "*=" ==================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator*= for std::vector.

    Increment each element in the input vector, using the corresping value
    on vector at the r.h.s. as increment, i.e.:
    x[i] *= y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator*= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator*= calls itself
    to increment x[i][j] by y[i][j], j = 0, ..., min(x[i].size(), y[i].size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator*= (
  std::vector< T >                              &x,
  const std::vector< T >                        &y
) {

    // ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = std::min(x.size(), y.size());

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] *= y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator*= for std::vector.

    Increment each element in the input vector, using the value
    on the r.h.s. as increment, i.e.:
    x[i] *= y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator*= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator*= calls itself
    to increment x[i][j] by y[j], j = 0, ... , min(x[i].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator*= (
  std::vector< T >                              &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] *= y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator*= for std::vector.

    Increment each element in the input vector, using the value
    on the r.h.s. as increment, i.e.:
    x[i][j] *= y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator*= is defined.

    The element-wise increment is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator*= calls itself
    to increment x[i][j][k] by y[k], k = 0, ... , min(x[i][j].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< std::vector< T > >& operator*= (
    std::vector< std::vector< T > >             &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    x[i] *= y;
} //next i

return (x); };

// Operator "/" ===================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::vector.

    Perform the element-wise division between two vectors (x and y) and returns a
    vector z s.t.
    z[i] = x[i] / y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator/ calls itself
    to perform the element-wise division between x[i] and y[i].

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise division between x and y.
*/
template <class T>
std::vector< T > operator/ (
  const std::vector< T >                        &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t             n = std::min(x.size(), y.size());
std::vector<T>          z(n, T(0));

// ================================================================================== //
// PERFORM DIVISION                                                                   //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    z[i] = x[i] / y[i];
};

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::vector.

    Perform the element-wise division between a vector (x) and a constant (y), 
    and returns a vector (z) s.t.
    z[i] = x[i] / y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator/ calls itself to perform the element-wise division
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise division between x and y.
*/
template <class T>
std::vector< T > operator/ (
  const std::vector< T >                        &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
size_t        n = x.size();
std::vector< T >   z(n, T(0));

// Counters
size_t  i;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (i = 0; i < n; i++) {
    z[i] = x[i]/y;
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::vector.

    Perform the element-wise division between a constant (x), and a vector (y)
    and returns a vector (z) s.t.
    z[i] = x / y[i] for all i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if the i-th element of y and
    x are std::vector, operator/ calls itself to perform the element-wise division
    between x[i] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector with n elements, storing the element-wise division between x and y.
*/
template <class T>
std::vector< T > operator/ (
  const T                                       &x,
  const std::vector< T >                        &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t         n = y.size();
std::vector< T >    z(n);

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    z[i] = x/y[i];
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::vector.

    Perform the element-wise division between a constant (y) and a vector of vectors (y),
    and returns z s.t.
    z[i][j] = x / y[i][j] for all j = 0, ..., y[i].size()-1, i = 0, ..., n-1
    where n = y.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if y[i][j] and
    y are std::vector, operator/ calls itself to perform the element-wise division
    between y[i][j] and x.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of y and
    storing the element-wise division between x and y.
*/
template <class T>
std::vector< std::vector< T > > operator/ (
  const T                                       &x,
  const std::vector< std::vector< T > >         &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = y.size();
std::vector< std::vector< T > >   z(n);

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    z[i] = x/y[i];
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/ for std::vector.

    Perform the element-wise division between a vector of vectors (x) and a constant (y), 
    and returns z s.t.
    z[i][j] = x[i][j] / y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator/ is defined.

    The element-wise division is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator/ calls itself to perform the element-wise division
    between x[i][j] and y.

    \param[in] x first argument
    \param[in] y second argument

    \result vector of vectors having the same dimensions of x and
    storing the element-wise division between x and y.
*/
template <class T>
std::vector< std::vector< T > > operator/ (
  const std::vector< std::vector< T > >         &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = x.size();
std::vector< std::vector< T > >   z(n);

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    z[i] = x[i]/y;
}

return(z); };

// Operator "/=" ==================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Operator/= for std::vector.

    Increment each element in the input vector, using the corresping value
    on vector at the r.h.s. as increment, i.e.:
    x[i] /= y[i] for all i = 0, ..., n-1
    where n = min(x.size(), y.size().
    Template parameters, T, can by any type for which the operator/= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    the i-th element of y are std::vector, operator/= calls itself
    to increment x[i][j] by y[i][j], j = 0, ..., min(x[i].size(), y[i].size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator/= (
  std::vector< T >                              &x,
  const std::vector< T >                        &y
) {

    // ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = std::min(x.size(), y.size());

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] /= y[i];
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/= for std::vector.

    Increment each element in the input vector, using the value
    on the r.h.s. as increment, i.e.:
    x[i] /= y for all i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator/= is defined.

    The element-wise increment is performed recursively, i.e. if the i-th element of x and
    y are std::vector, operator/= calls itself
    to increment x[i][j] by y[j], j = 0, ... , min(x[i].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< T >& operator/= (
  std::vector< T >                              &x,
  const T                                       &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t              n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++){
    x[i] /= y;
};

return (x); };

// ---------------------------------------------------------------------------------- //
/*!
    Operator/= for std::vector.

    Increment each element in the input vector, using the value
    on the r.h.s. as increment, i.e.:
    x[i][j] /= y for all j = 0, ..., x[i].size()-1, i = 0, ..., n-1
    where n = x.size().
    Template parameters, T, can by any type for which the operator/= is defined.

    The element-wise increment is performed recursively, i.e. if x[i][j] and
    y are std::vector, operator/= calls itself
    to increment x[i][j][k] by y[k], k = 0, ... , min(x[i][j].size(), y.size())-1

    \param[in] x first argument
    \param[in] y second argument

    \result first argument incremented with r.h.s. values.
*/
template <class T>
std::vector< std::vector< T > >& operator/= (
    std::vector< std::vector< T > >             &x,
    const T                                     &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
size_t                  n = x.size();

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (size_t i = 0; i < n; i++) {
    x[i] /= y;
} //next i

return (x); };

// Output operator ================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Insertion operator for std::vector.

    Flush the content of std::vector to std::ostream.
    The content of the input vector is flushed with the following format:
    x[0] x[1] x[2] ... x[n-1] where n = x.size();
    (i.e. vector elements are separated by blank spaces).
    Template parameter T can be any type such that operator<< is defined.

    \param[in,out] out output stream
    \param[in] x argument of insertion operator

    \result reference to the stream (allows concatenation)
*/
template <class T>
std::ostream& operator<< (
    std::ostream                                &out,
    const std::vector< T >                      &x
) {

    // ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
size_t               n = x.size();

// Counters
// none

// ================================================================================== //
// OUTPUT VECTOR CONTENT                                                              //
// ================================================================================== //
if (n == 0) return(out);
for (size_t i = 0; i < n-1; i++) {
    out << x[i] << " ";
} //next i
out << x[n-1];

return(out); };

// ---------------------------------------------------------------------------------- //
/*!
    Insertion operator for std::vector.

    Flush the content of std::vector to std::ofstream.
    The content of the input vector is flushed with the following format:
    x[0] x[1] x[2] ... x[n-1] where n = x.size();
    (i.e. vector elements are separated by blank spaces).
    Template parameter T can be any type such that operator<< is defined.

    \param[in,out] out output file stream
    \param[in] x argument of insertion operator

    \result reference to the stream (allows concatenation)
*/
template <class T>
std::ofstream& operator<< (
    std::ofstream                               &out,
    const std::vector< T >                      &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
size_t               n = x.size();

// Counters
// none

// ================================================================================== //
// OUTPUT VECTOR CONTENT                                                              //
// ================================================================================== //
if (n == 0) {
    return(out);
}
for (size_t i = 0; i < n-1; i++) {
    out << x[i] << " ";
} //next i
out << x[n-1];

return(out); };

// Input operator =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Extraction operator for std::vector.

    Extract the content of std::vector from std::istream.
    The content of the input vector is extracted until end-of-stream condition
    is met. Element extracted from the stream are copyied into vector starting from position
    0. When all the availble position within the vector have been overwritten, further
    elements will be added to the container using push_back, therefore increasing
    container's size.

    \param[in,out] in input stream
    \param[in,out] x argument of extraction operator

    \result reference to input stream (allows concatenation)
*/
template <class T>
std::istream& operator>> (
    std::istream                                &in,
    std::vector< T >                            &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T                    dummy;
size_t               n = x.size();

// Counters
size_t               i;

// ================================================================================== //
// EXTRACT STREAM CONTENT INTO VECTOR                                                 //
// ================================================================================== //
i = 0;
while ((in.good()) && (i < n)) {
    if (in >> dummy) { x[i] = dummy; i++; }
} //next i
while (in.good()) {
    if (in >> dummy) { x.push_back(dummy); }
}

return(in); };

// ---------------------------------------------------------------------------------- //
/*!
    Extraction operator for std::vector.

    Extract the content of std::vector from std::ifstream.
    The content of the input vector is extracted until end-of-file condition
    is met. Element extracted from the stream are copyied into vector starting from position
    0. When all the availble position within the vector have been overwritten, further
    elements will be added to the container using push_back, therefore increasing
    container's size.

    \param[in,out] in input file stream
    \param[in,out] x argument of extraction operator

    \result reference to input file stream (allows concatenation)
*/
template <class T>
std::ifstream& operator>> (
    std::ifstream                               &in,
    std::vector< T >                            &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T       dummy;

// Counters
int     i, n = x.size();

// ================================================================================== //
// EXTRACT FILE CONTENT INTO VECTOR                                                   //
// ================================================================================== //
i = 0;
while ((in.good()) && (i < n)) {
    if (in >> dummy) { x[i] = dummy; i++; }
} //next i
while (in.good()) {
    if (in >> dummy) { x.push_back(dummy); }
}

return(in); };

// Miscellanea ====================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Dummy function for recursive templated routines display

    \param[in,out] out output stream
    \param[in] x variable to be displayed

    \result reference to output stream
*/
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type*>
std::ostream& display(
    std::ostream                                &out,
    const T                                     &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// DISPLAY VARIABLE                                                                   //
// ================================================================================== //
out << x;

return(out); }

// ---------------------------------------------------------------------------------- //
/*!
    Display vector in a nicely formatted to a std::ostream

    \param[in,out] out output stream
    \param[in] x vector to be displayed
    \param[in] padding (default = 0) number of trailing spaces

    \result reference to output stream
*/
template<class T>
std::ostream& display(
    std::ostream                                &out,
    const std::vector<T>                        &x,
    unsigned int                                 padding
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Counters
typename std::vector<T>::const_iterator      i, e = x.cend();

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
    Dummy function for recursive templated routines display

    \param[in,out] out output stream
    \param[in] x variable to be displayed

    \result reference to output stream
*/
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type*>
std::ofstream& display(
    std::ofstream                               &out,
    const T                                     &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// DISPLAY VARIABLE                                                                   //
// ================================================================================== //
out << x;

return(out); }

// ---------------------------------------------------------------------------------- //
/*!
    Display vector in a nicely formatted to a std::ofstream

    \param[in,out] out output file stream
    \param[in] x vector to be displayed
    \param[in] padding (default = 0) number of trailing spaces

    \result reference to output stream
*/
template<class T>
std::ofstream& display(
    std::ofstream                               &out,
    const std::vector<T>                        &x,
    unsigned int                                 padding
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Counters
typename std::vector<T>::const_iterator      i, e = x.cend();

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
