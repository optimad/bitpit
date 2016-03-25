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
//                                 MATH OPERATORS                                     //
//                                                                                    //
// Basic math operators for std::array.                                               //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

/*!
   @ingroup MathFunctions
   @{
 */

// ================================================================================== //
// TEMPLATE IMPLEMENTATIONS                                                           //
// ================================================================================== //

// Operator "min" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between two arrays.
    Given two array x and y, returns z such that:
    z[i] = min(x[i], y[i]), for all i = 0, ..., d-1
    where the d is the size of x and y.
    
    Template parameters T can be any type such that min is defined.
    For instance if T = std::array<T1, e>, min function call itself to return
    the element-wise minimum between x[i] and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of x and y, storing
    the element-wise minimum between x and y.
*/
template <class T, size_t d>
std::array<T, d> min(
    const std::array<T, d>                      &x,
    const std::array<T, d>                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>          z;

// ================================================================================== //
// COMPARE VECTORS                                                                    //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = min(x[i], y[i]);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between array and constant.
    Given a array x and a constant y, returns z such that:
    z[i] = min(x[i], y), for all i = 0, ..., d-1
    where d is the size of x.

    Template parameters T can be any type such that min is defined.
    For instance if T = std::array<T1, e>, min function call itself to return
    the element-wise minimum between x[i] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of x, storing
    the element-wise minimum between x and y.
*/
template <class T, size_t d>
std::array<T, d> min(
    const std::array<T, d>                      &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>    z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = min(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between constant and array.
    Given a constant x and a array y, returns z such that:
    z[i] = min(x, y[i]), for all i = 0, ..., d-1
    where d is the size of y.

    Template parameters T can be any type such that min is defined.
    For instance if T = std::array<T1, e>, min function call itself to return
    the element-wise minimum between x and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of y, storing
    the element-wise minimum between x and y.
*/
template <class T, size_t d>
std::array<T, d> min(
    const T                                     &x,
    const std::array<T, d>                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>    z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = min(y, x);

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between 2D array and constant.
    Given a 2D array x and a constant y, returns z such that:
    z[i][j] = min(x[i][j], y), for all j = 0, ..., n-1, i = 0, ..., d-1
    where d is the size of x.

    Template parameters T can be any type such that min is defined.
    For instance if T = std::array<T1, e>, min function call itself to return
    the element-wise minimum between x[i][j] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of x, storing
    the element-wise minimum between x and y.
*/
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> min(
    const std::array<std::array<T, n>, d>       &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, n>, d>    z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = min(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between constant and 2D array.
    Given a constant x and a 2D array y, returns z such that:
    z[i][j] = min(x, y[i][j]), for all j = 0, ..., n-1, i = 0, ..., d-1
    where d is the size of y.

    Template parameters T can be any type such that min is defined.
    For instance if T = std::array<T1, e>, min function call itself to return
    the element-wise minimum between x and y[i][j].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a array having the same dimensions of y, storing
    the element-wise minimum between x and y.
*/
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> min(
    const T                                     &x,
    const std::array<std::array<T, n>, d>       &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, n>, d>    z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = min(y, x);

return(z); };

// Operator "minval" ================================================================ //

// ---------------------------------------------------------------------------------- //
/*!
    Returns the element with the smallest value within a array, i.e. given an input array, x
    min_value = min(x[i]) over all i = 0, ..., d-1

    Parameters template can be of any type with the following requirements:
    1. minval must be defined for any type T
    2. type T1 must be a scalar type
    (for instance, T = std::array<double, e>, T1 = double)

    \param[in] x input array
    \param[in,out] min_value on output stores the elements with the smallest value.
*/
template <class T, size_t d, class T1>
void minval(
    const std::array<T, d>                      &x,
    T1                                          &min_value
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T1               value;

// Counters
size_t           i;

// ================================================================================== //
// FIND THE MIN-VALUE                                                                 //
// ================================================================================== //
if (d > 0) {
    minval(x[0], min_value);
    for (i = 1; i < d; i++) {
        minval(x[i], value);
        if (value < min_value) {
            min_value = value;
        }
    } //next i
}

return; };

// Operator "max" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between two arrays.
    Given two array x and y, returns z such that:
    z[i] = max(x[i], y[i]), for all i = 0, ..., d-1
    where the d is the size of x and y.
    
    Template parameters T can be any type such that max is defined.
    For instance if T = std::array<T1, e>, max function call itself to return
    the element-wise maximum between x[i] and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of x and y, storing
    the element-wise maximum between x and y.
*/
template <class T, size_t d>
std::array<T, d> max(
    const std::array<T, d>                      &x,
    const std::array<T, d>                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>    z;

// ================================================================================== //
// COMPARE VECTORS                                                                    //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = max(x[i], y[i]);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between array and constant.
    Given a array x and a constant y, returns z such that:
    z[i] = max(x[i], y), for all i = 0, ..., d-1
    where d is the size of x.

    Template parameters T can be any type such that max is defined.
    For instance if T = std::array<T1, e>, max function call itself to return
    the element-wise maximum between x[i] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of x, storing
    the element-wise maximum between x and y.
*/
template <class T, size_t d>
std::array<T, d> max(
    const std::array<T, d>                      &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>      z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = max(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between constant and array.
    Given a constant x and a array y, returns z such that:
    z[i] = max(x, y[i]), for all i = 0, ..., d-1
    where d is the size of y.

    Template parameters T can be any type such that max is defined.
    For instance if T = std::array<T1, e>, max function call itself to return
    the element-wise maximum between x and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of y, storing
    the element-wise maximum between x and y.
*/
template <class T, size_t d>
std::array<T, d> max(
    const T                                     &x,
    const std::array<T, d>                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<T, d>       z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = max(y, x);

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between 2D array and constant.
    Given a 2D array x and a constant y, returns z such that:
    z[i][j] = max(x[i][j], y), for all j = 0, ..., n-1, i = 0, ..., d-1
    where d is the size of x.

    Template parameters T can be any type such that max is defined.
    For instance if T = std::array<T1, e>, max function call itself to return
    the element-wise maximum between x[i][j] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns an array having the same dimensions of x, storing
    the element-wise maximum between x and y.
*/
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> max(
    const std::array<std::array<T, n>, d>       &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, n>, d>      z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (size_t i = 0; i < d; i++) {
    z[i] = max(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between constant and 2D array.
    Given a constant x and a 2D array y, returns z such that:
    z[i][j] = max(x, y[i][j]), for all j = 0, ..., n-1, i = 0, ..., d-1
    where d is the size of y.

    Template parameters T can be any type such that max is defined.
    For instance if T = std::array<T1, e>, max function call itself to return
    the element-wise maximum between x and y[i][j].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a array having the same dimensions of y, storing
    the element-wise maximum between x and y.
*/
template <class T, size_t d, size_t n>
std::array<std::array<T, n>, d> max(
    const T                                     &x,
    const std::array<std::array<T, n>, d>       &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::array<std::array<T, n>, d>       z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = max(y, x);

return(z); };


// Operator "maxval" ================================================================ //

// ---------------------------------------------------------------------------------- //
/*!
    Returns the element with the largest value within a array, i.e. given an input array, x
    max_value = max(x[i]) over all i = 0, ..., n-1
    where n = x.size().

    Parameters template can be of any type with the following requirements:
    1. minval must be defined for any type T
    2. type T1 must be a scalar type
    (for instance, T = std::array<double, e>, T1 = double)

    \param[in] x input array
    \param[in,out] max_value on output stores the elements with the largest value.
*/
template <class T, size_t d, class T1>
void maxval(
    const std::array<T, d>                      &x,
    T1                                          &max_value
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T1               value;

// Counters
size_t           i;

// ================================================================================== //
// FIND THE MIN-VALUE                                                                 //
// ================================================================================== //
if (d > 0) {
    maxval(x[0], max_value);
    for (i = 1; i < d; i++) {
        maxval(x[i], value);
        if (value > max_value) {
            max_value = value;
        }
    } //next i
}

return; };

// Operator "sum" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Given a input array, x, returns the sum of its elements, i.e.:
    s = sum (x[i]) over all i = 0, ..., n-1
    where n = x.size().

    Parameters template can be of any type with the following requirements:
    1. operator += must be defined for type T
    2. type T1 must be a scalar type
    (for instance, T = std::array<double, e>, T1 = double)

    \param[in] x input array
    \param[in,out] s sum of element in x.
*/
template <class T, size_t d, class T1>
void sum(
    const std::array<T, d>                      &x,
    T1                                          &s
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T1            value;

// Counters
size_t        i;

// ================================================================================== //
// PERFORM SUMMATION                                                                  //
// ================================================================================== //
if (d > 0) {
    sum(x[0], s);
    for (i = 1; i < d; i++) {
        sum(x[i], value);
        s += value;
    } //next i
}

return; };

// Operator "abs" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Given a input array, x, returns a array storing the absolute value of elements
    in x, i.e.:
    z[i] = abs(x[i]), for i = 0, ..., d-1
    where d is the size of x.

    Template parameter can be any type such that abs() function is defined.
    For instance if  T = std::array<T1, e>, the abs() function will call itself
    to return the absolute value of elements in x[i].

    \param[in] x input array

    \result array having the same dimensions of x storing the absolute value
    of the elements in x.
*/  
template <class T, size_t d>
std::array<T, d> abs(
    const std::array<T, d>                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
std::array<T, d>             z;

// Counters
size_t                       i;

// ================================================================================== //
// COMPUTE THE ABSOLUTE VALUE OF A VECTOR                                             //
// ================================================================================== //
if (d > 0) {
    for (i = 0; i < d; ++i) {
        z[i] = abs(x[i]);
    } //next i
}

return(z); };

// Operator "pow" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Given a input array, x, returns a array storing the p-th power of its elements, i.e.:
    z[i] = pow(x[i], p), for i = 0, ..., d-1
    where d is the size of x.

    Template parameter can be any type such that pow function is defined.
    For instance if  T = std::array<T1, e>, the pow() function will call itself
    to return the p-th power of the elements in x[i].

    \param[in] x input array
    \param[in] p power index

    \result array having the same dimensions of x storing the p-th power
    of the elements in x.
*/  
template <class T, size_t d>
std::array<T, d> pow(
    std::array<T, d>                            &x,
    double                                       p
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
std::array<T, d>      y;

// Counters
size_t                i;

// ================================================================================== //
// COMPUTE ELEMENT-WISE POWER                                                         //
// ================================================================================== //
for (i = 0; i < d; i++) {
    y[i] = pow(x[i], p);
} //next i

return(y); };

// Operator "norm" ================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Compute the 1-norm of a input array, x, i.e.:
    n = sum(abs(x)).

    Template parameter can be any scalar type

    \param[in] x input array

    \result on output returns the 1-norm of the input array.
*/
template <class T, size_t d>
double norm1(
    const std::array<T, d>                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
double          z = 0.0;

// Counters
size_t          i;

// ================================================================================== //
// COMPUTE THE P-NORM                                                                 //
// ================================================================================== //
if (d > 0) {
    for (i = 0; i < d; i++) {
        z += abs(x[i]);
    } //next i
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the 2-norm of a input array, x, i.e.:
    n = sqrt(sum(pow(x, 2))).

    Template parameter can be any scalar type

    \param[in] x input array

    \result on output returns the 2-norm of the input array.
*/
template <class T, size_t d>
double norm2(
    const std::array<T, d>                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
double          z = 0.0;

// Counters
size_t          i;

// ================================================================================== //
// COMPUTE THE P-NORM                                                                 //
// ================================================================================== //
if (d > 0) {
    for (i = 0; i < d; i++) {
        z += x[i]*x[i];
    } //next i
}

return(sqrt(z)); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the generic norm of a input array, x, i.e.:
    n = pow( sum( pow( abs(x), p ) ), 1/p ).

    Template parameter can be any scalar type

    \param[in] x input array
    \param[in] p norm index

    \result on output returns the p-norm of the input array.
*/
template <class T, size_t d>
double norm(
    const std::array<T, d>                      &x,
    int                                          p
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
double          z = 0.0;
double          t, y;

// Counters
size_t          i;
int             j;

// ================================================================================== //
// COMPUTE THE P-NORM                                                                 //
// ================================================================================== //
if (p == 1) { return(norm1(x)); }
if (p == 2) { return(norm2(x)); }

if (d > 0) {
    for (i = 0; i < d; i++) {
        y = 1.0;
        t = abs(x[i]);
        for (j = 1; j <= p; j++) {
            y = y*t;
        } //next j
        //z += pow(abs(x[i]), p);
        z += y;
    } //next i
}

return(std::exp(std::log(std::max(z, 1.0e-307))/((double) p))); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the infinity-norm of a input array, x, i.e.:
    n = maxval(abs(x)).

    Template parameter can be any scalar type

    \param[in] x input array

    \result on output returns the inf-norm of the input array.
*/
template <class T, size_t d>
double normInf(
    const std::array<T, d>                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
double          z = 0.0, y;

// Counters
size_t          i;

// ================================================================================== //
// COMPUTE THE inf-NORM                                                               //
// ================================================================================== //
if (d > 0) {
    z = abs(x[0]);
    for (i = 1; i < d; i++) {
        y = abs(x[i]);
        if (y > z) {
            z = y;
        }
    } //next i
}

return(z); };

// Operator "dotProduct" =========================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Compute the scalar product of 2 input arrays, x and y, i.e.:
    d = sum(x * y).

    Template parameter can be any scalar type

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result on output returns the scalar product of x and y.
*/
template <class T, size_t d>
T dotProduct(
    const std::array<T, d>                      &x,
    const std::array<T, d>                      &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T                   dp = ((T) 0.0);

// Counters
size_t              i;

// ================================================================================== //
// COMPUTE THE DOT PRODUCT                                                            //
// ================================================================================== //
if (d > 0) {
    for (i = 0; i < d; i++) {
        dp += x[i]*y[i];
    } //next i
}

return(dp); };

// Operator "crossProduct" ========================================================= //

// ---------------------------------------------------------------------------------- //
/*!
    Compute the cross product in R3 of input arrays, x and y.
    Template parameter can be any scalar type

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result on output returns the cross product product of x and y.
*/
template <class T>
std::array<T, 3> crossProduct(
    const std::array<T, 3>                      &x,
    const std::array<T, 3>                      &y
) {

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //
std::array<T, 3>      z;

// =================================================================================== //
// COMPUTE THE EXTERNAL PRODUCT                                                        //
// =================================================================================== //
z[0] = x[1] * y[2] - x[2] * y[1];
z[1] = x[2] * y[0] - x[0] * y[2];
z[2] = x[0] * y[1] - x[1] * y[0];

return (z); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the cross product in R2 of input arrays, x and y.
    Template parameter can be any scalar type

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result on output returns the cross product product of x and y.
*/
template <class T>
T crossProduct(
    const std::array<T, 2>                      &x,
    const std::array<T, 2>                      &y
) {

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //
T      z;

// =================================================================================== //
// COMPUTE THE EXTERNAL PRODUCT                                                        //
// =================================================================================== //
z = x[0] * y[1] - x[1] * y[0];

return (z); };

/*!
   @}
 */
