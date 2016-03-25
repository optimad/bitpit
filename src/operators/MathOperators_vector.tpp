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
// Basic math operators for vectors.                                                  //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// ================================================================================== //
// TEMPLATE IMPLEMENTATIONS                                                           //
// ================================================================================== //
/*!
   @ingroup MathFunctions
   @{
 */

// Operator "min" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between two vectors.
    Given two vectors x and y, returns z such that:
    z[i] = min(x[i], y[i]), for all i = 0, ..., n-1
    where the n = min(x.size(), y.size()).
    
    Template parameters T can be any type such that min is defined.
    For instance if T = std::vector<T1>, min function call itself to return
    the element-wise minimum between x[i] and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of x and y, storing
    the element-wise minimum between x and y.
*/
template <class T>
std::vector< T > min(
    const std::vector< T >                      &x,
    const std::vector< T >                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int            n = std::min(x.size(), y.size());
int            m = std::max(x.size(), y.size());
std::vector< T >    z(m);

// ================================================================================== //
// COMPARE VECTORS                                                                    //
// ================================================================================== //
for (int i = 0; i < n; i++) {
    z[i] = min(x[i], y[i]);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between vector and constant.
    Given a vector x and a constant y, returns z such that:
    z[i] = min(x[i], y), for all i = 0, ..., n-1
    where n = x.size().

    Template parameters T can be any type such that min is defined.
    For instance if T = std::vector<T1>, min function call itself to return
    the element-wise minimum between x[i] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of x, storing
    the element-wise minimum between x and y.
*/
template <class T>
std::vector< T > min(
    const std::vector< T >                      &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int            n = x.size();
std::vector< T >    z(n);

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (int i = 0; i < n; i++) {
    z[i] = min(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between constant and vector.
    Given a constant x and a vector y, returns z such that:
    z[i] = min(x, y[i]), for all i = 0, ..., d-1
    where d is the size of y.

    Template parameters T can be any type such that min is defined.
    For instance if T = std::vector<T1>, min function call itself to return
    the element-wise minimum between x and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of y, storing
    the element-wise minimum between x and y.
*/
template <class T>
std::vector< T > min(
    const T                                     &x,
    const std::vector< T >                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< T >    z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = min(y, x);

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between 2D vector and constant.
    Given a 2D vector x and a constant y, returns z such that:
    z[i][j] = min(x[i][j], y), for all j = 0, ..., x[i].size(), i = 0, ..., n-1
    where n = x.size().

    Template parameters T can be any type such that min is defined.
    For instance if T = std::vector<T1>, min function call itself to return
    the element-wise minimum between x[i][j] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of x, storing
    the element-wise minimum between x and y.
*/
template <class T>
std::vector< std::vector< T > > min(
    const std::vector< std::vector < T > >      &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int                      n = x.size();
std::vector< std::vector< T > >    z(n);

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (int i = 0; i < n; i++) {
    z[i] = min(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise minimum between constant and 2D vector.
    Given a constant x and a 2D vector y, returns z such that:
    z[i][j] = min(x, y[i][j]), for all j = 0, ..., y[i].size(), i = 0, ..., n-1
    where n = y.size().

    Template parameters T can be any type such that min is defined.
    For instance if T = std::vector<T1>, min function call itself to return
    the element-wise minimum between x and y[i][j].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of y, storing
    the element-wise minimum between x and y.
*/
template <class T>
std::vector< std::vector < T > > min(
    const T                                     &x,
    const std::vector< std::vector< T > >       &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< std::vector< T > >    z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = min(y, x);

return(z); };

// Operator "minval" ================================================================ //

// ---------------------------------------------------------------------------------- //
/*!
    Overloading of minval() function for scalar types.
    
    \param[in] x input scalar
    \param[in,out] min_value on output returns the input value
*/
template <typename T, typename std::enable_if< std::is_scalar< T >::value>::type*>
void inline minval(
    const T                                     &x,
    T                                           &min_value
) {

// ================================================================================== //
// template <typename T,                                                              //
//           typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>   //
// void inline minval(                                                                //
//     const T              &x,                                                       //
//     T                    &min_value)                                               //
//                                                                                    //
// Overloading of minval function for scalar types.                                   //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x         : T, input scalar                                                      //
// - min_value : T, on output stores the input value.                                 //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
// none

// Counters
// none

// ================================================================================== //
// COMPUTE THE MIN VALUE                                                              //
// ================================================================================== //
min_value = x;

return; };

// ---------------------------------------------------------------------------------- //
/*!
    Returns the element with the smallest value within a vector, i.e. given an input vector, x
    min_value = min(x[i]) over all i = 0, ..., n-1
    where n = x.size().
    
    Parameters template can be of any type with the following requirements:
    1. minval must be defined for any type T
    2. type T1 must be a scalar type
    (for instance, T = std::vector<double>, T1 = double)

    \param[in] x input vector
    \param[in,out] min_value on output stores the elements with the smallest value.
*/
template <class T, class T1>
void minval(
    const std::vector<T>                        &x,
    T1                                          &min_value
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T1               value;

// Counters
int              i, n;

// ================================================================================== //
// FIND THE MIN-VALUE                                                                 //
// ================================================================================== //
n = x.size();
if (n > 0) {
    minval(x[0], min_value);
    for (i = 1; i < n; i++) {
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
    Element-wise maximum between two vectors.
    Given two vectors x and y, returns z such that:
    z[i] = max(x[i], y[i]), for all i = 0, ..., n-1
    where the n = min(x.size(), y.size()).
    
    Template parameters T can be any type such that max is defined.
    For instance if T = std::vector<T1>, max function call itself to return
    the element-wise maximum between x[i] and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of x and y, storing
    the element-wise maximum between x and y.
*/
template <class T>
std::vector<T> max(
    const std::vector< T >                      &x,
    const std::vector< T >                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int            n = std::min(x.size(), y.size());
int            m = std::max(x.size(), y.size());
std::vector< T >    z(m);

// ================================================================================== //
// COMPARE VECTORS                                                                    //
// ================================================================================== //
for (int i = 0; i < n; i++) {
    z[i] = max(x[i], y[i]);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between vector and constant.
    Given a vector x and a constant y, returns z such that:
    z[i] = max(x[i], y), for all i = 0, ..., n-1
    where n = x.size().

    Template parameters T can be any type such that max is defined.
    For instance if T = std::vector<T1>, max function call itself to return
    the element-wise maximum between x[i] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of x, storing
    the element-wise maximum between x and y.
*/
template <class T>
std::vector< T > max(
    const std::vector< T >                      &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int            n = x.size();
std::vector<T>      z(n);

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (int i = 0; i < n; i++) {
    z[i] = max(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between constant and vector.
    Given a constant x and a vector y, returns z such that:
    z[i] = max(x, y[i]), for all i = 0, ..., d-1
    where d is the size of y.

    Template parameters T can be any type such that max is defined.
    For instance if T = std::vector<T1>, max function call itself to return
    the element-wise maximum between x and y[i].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of y, storing
    the element-wise maximum between x and y.
*/
template <class T>
std::vector< T > max(
    const T                                     &x,
    const std::vector< T >                      &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< T >       z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = max(y, x);

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between 2D vector and constant.
    Given a 2D vector x and a constant y, returns z such that:
    z[i][j] = max(x[i][j], y), for all j = 0, ..., x[i].size(), i = 0, ..., n-1
    where n = x.size().

    Template parameters T can be any type such that max is defined.
    For instance if T = std::vector<T1>, max function call itself to return
    the element-wise maximum between x[i][j] and y.

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of x, storing
    the element-wise maximum between x and y.
*/
template <class T>
std::vector< std::vector< T > > max(
    const std::vector< std::vector< T > >       &x,
    const T                                     &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int                        n = x.size();
std::vector< std::vector< T > >      z(n);

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
for (int i = 0; i < n; i++) {
    z[i] = max(x[i], y);
} //next i

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Element-wise maximum between constant and 2D vector.
    Given a constant x and a 2D vector y, returns z such that:
    z[i][j] = max(x, y[i][j]), for all j = 0, ..., y[i].size(), i = 0, ..., n-1
    where n = y.size().

    Template parameters T can be any type such that max is defined.
    For instance if T = std::vector<T1>, max function call itself to return
    the element-wise maximum between x and y[i][j].

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result returns a vector having the same dimensions of y, storing
    the element-wise maximum between x and y.
*/
template <class T>
std::vector< std::vector< T > > max(
    const T                                     &x,
    const std::vector< std::vector< T > >       &y
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
std::vector< std::vector< T > >       z;

// ================================================================================== //
// COMPARE VECTOR AND SCALAR                                                          //
// ================================================================================== //
z = max(y, x);

return(z); };

// Operator "maxval" ================================================================ //

// ---------------------------------------------------------------------------------- //
/*!
    Overloading of maxval() function for scalar types.
    
    \param[in] x input scalar
    \param[in,out] max_value on output returns the input value
*/
template <typename T, typename std::enable_if< std::is_scalar< T >::value>::type*>
void inline maxval(
    const T                                     &x,
    T                                           &max_value
) {

// ================================================================================== //
// template <typename T,                                                              //
//           typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>   //
// void inline maxval(                                                                //
//     const T              &x,                                                       //
//     T                    &max_value)                                               //
//                                                                                    //
// Overloading of maxval function for scalar types.                                   //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x         : T, input scalar                                                      //
// - max_value : T, on output stores the input value                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
// none

// Counters
// none

// ================================================================================== //
// COMPUTE THE MIN VALUE                                                              //
// ================================================================================== //
max_value = x;

return; };

// ---------------------------------------------------------------------------------- //
/*!
    Returns the element with the largest value within a vector, i.e. given an input vector, x
    max_value = max(x[i]) over all i = 0, ..., n-1
    where n = x.size().
    
    Parameters template can be of any type with the following requirements:
    1. minval must be defined for any type T
    2. type T1 must be a scalar type
    (for instance, T = std::vector<double>, T1 = double)

    \param[in] x input vector
    \param[in,out] max_value on output stores the elements with the largest value.
*/
template <class T, class T1>
void maxval(
    const std::vector<T>                        &x,
    T1                                          &max_value
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T1               value;

// Counters
int              i, n;

// ================================================================================== //
// FIND THE MIN-VALUE                                                                 //
// ================================================================================== //
n = x.size();
if (n > 0) {
    maxval(x[0], max_value);
    for (i = 1; i < n; i++) {
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
    Overloading of sum() function for scalar type.

    \param[in] x input scalar
    \param[in,out] s on output stores the value of the input scalar
*/
template<class T, typename std::enable_if< std::is_scalar< T >::value>::type*>
void inline sum(
    const T                                     &x,
    T                                           &s
) {

// ================================================================================== //
// template<class T,                                                                  //
//          typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>    //
// void inline sum(                                                                   //
//     const T              &x,                                                       //
//     T                    &s)                                                       //
//                                                                                    //
// Overloading of sum() function for scalar type.                                     //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x         : int, dummy input                                                     //
// - s         : int, dummy output                                                    //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
// none

// Counters
// none

// ================================================================================== //
// PERFORM DUMMY SUM                                                                  //
// ================================================================================== //
s = x;

return; };

// ---------------------------------------------------------------------------------- //
/*!
    Given a input vector, x, returns the sum of its elements, i.e.:
    s = sum (x[i]) over all i = 0, ..., n-1
    where n = x.size().

    Parameters template can be of any type with the following requirements:
    1. operator += must be defined for type T
    2. type T1 must be a scalar type
    (for instance, T = std::vector<double>, T1 = double)

    \param[in] x input vector
    \param[in,out] s sum of element in x.
*/
template <class T, class T1>
void sum(
    const std::vector< T >                      &x,
    T1                                          &s
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int           n = x.size();
T1            value;

// Counters
int           i;

// ================================================================================== //
// PERFORM SUMMATION                                                                  //
// ================================================================================== //
if (n > 0) {
    sum(x[0], s);
    for (i = 1; i < n; i++) {
        sum(x[i], value);
        s += value;
    } //next i
}

return; };

// Operator "abs" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Given a ipnut vector, x, returns a vector storing the absolute value of elements
    in x, i.e.:
    z[i] = abs(x[i]), for i = 0, ..., n-1
    where n is the size of x.

    Template parameter can be any type such that abs() function is defined.
    For instance if  T = std::vector<T1>, the abs() function will call itself
    to return the absolute value of elements in x[i].

    \param[in] x input vector

    \result vector having the same dimensions of x storing the absolute value
    of the elements in x.
*/  
template <class T>
std::vector<T> abs(
    const std::vector< T >                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int                     n = x.size();
std::vector< T >             z;

// Counters
int                     i;

// ================================================================================== //
// COMPUTE THE ABSOLUTE VALUE OF A VECTOR                                             //
// ================================================================================== //
if (n > 0) {
    z.resize(n);
    for (i = 0; i < n; i++) {
        z[i] = abs(x[i]);
    } //next i
}

return(z); };

// Operator "pow" =================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Given a input vector, x, returns a vector storing the p-th power of its elements, i.e.:
    z[i] = pow(x[i], p), for i = 0, ..., n-1
    where n is the size of x.

    Template parameter can be any type such that pow function is defined.
    For instance if  T = std::vector<T1>, the pow() function will call itself
    to return the p-th power of the elements in x[i].

    \param[in] x input vector
    \param[in] p power index

    \result vector having the same dimensions of x storing the p-th power
    of the elements in x.
*/  
template <class T>
std::vector< T > pow(
    std::vector< T >                            &x,
    double                                       p
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int              n;
std::vector< T >      y;

// Counters
int              i;

// ================================================================================== //
// COMPUTE ELEMENT-WISE POWER                                                         //
// ================================================================================== //
n = x.size();
y.resize(n);
for (i = 0; i < n; i++) {
    y[i] = pow(x[i], p);
} //next i

return(y); };

// Operator "norm" ================================================================== //

// ---------------------------------------------------------------------------------- //
/*!
    Compute the 1-norm of a input vector, x, i.e.:
    n = sum(abs(x)).

    Template parameter can be any scalar type

    \param[in] x input vector

    \result on output returns the 1-norm of the input vector.
*/
template <class T>
double norm1(
    const std::vector< T >                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int             n = x.size();
double          z = 0.0;

// Counters
int             i;

// ================================================================================== //
// COMPUTE THE P-NORM                                                                 //
// ================================================================================== //
if (n > 0) {
    for (i = 0; i < n; i++) {
        z += abs(x[i]);
    } //next i
}

return(z); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the 2-norm of a input vector, x, i.e.:
    n = sqrt(sum(pow(x, 2))).

    Template parameter can be any scalar type

    \param[in] x input vector

    \result on output returns the 2-norm of the input vector.
*/
template <class T>
double norm2(
    const std::vector< T >                      &x
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int             n = x.size();
double          z = 0.0;

// Counters
int             i;

// ================================================================================== //
// COMPUTE THE P-NORM                                                                 //
// ================================================================================== //
if (n > 0) {
    for (i = 0; i < n; i++) {
        z += x[i]*x[i];
    } //next i
}

return(sqrt(z)); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the generic norm of a input vector, x, i.e.:
    n = pow( sum( pow( abs(x), p ) ), 1/p ).

    Template parameter can be any scalar type

    \param[in] x input vector
    \param[in] p norm index

    \result on output returns the p-norm of the input vector.
*/
template <class T>
double norm(
    const std::vector< T >                      &x,
    int                                          p
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int             n = x.size();
double          z = 0.0;
double          t, y;

// Counters
int             i, j;

// ================================================================================== //
// COMPUTE THE P-NORM                                                                 //
// ================================================================================== //
if (p == 1) { return(norm1(x)); }
if (p == 2) { return(norm2(x)); }

if (n > 0) {
    for (i = 0; i < n; i++) {
        y = 1.0;
        t = x[i];
        for (j = 1; j <= p; j++) {
            y = y*t;
        } //next j
        z += abs(y);
    } //next i
}

return(std::exp(std::log(std::max(z, 1.0e-307))/((double) p))); };

// ---------------------------------------------------------------------------------- //
/*!
    Compute the infinity-norm of a input vector, x, i.e.:
    n = maxval(abs(x)).

    Template parameter can be any scalar type

    \param[in] x input vector

    \result on output returns the inf-norm of the input vector.
*/
template <class T>
double normInf(
    const std::vector< T >                      &x
) {

using namespace std;

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int             n = x.size();
double          z = 0.0, y;

// Counters
int             i;

// ================================================================================== //
// COMPUTE THE inf-NORM                                                               //
// ================================================================================== //
if (n > 0) {
    z = abs(x[0]);
    for (i = 1; i < n; i++) {
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
    Compute the scalar product of 2 input vectors, x and y, i.e.:
    d = sum(x * y).

    Template parameter can be any scalar type

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result on output returns the scalar product of x and y.
*/
template <class T>
T dotProduct(
    const std::vector< T >                      &x,
    const std::vector< T >                      &y
) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
int                 n = x.size(), m = y.size();
T                   dp = ((T) 0.0);

// Counters
int                 i;

// ================================================================================== //
// COMPUTE THE DOT PRODUCT                                                            //
// ================================================================================== //
if ((n > 0) && (n == m)) {
    for (i = 0; i < n; i++) {
        dp += x[i]*y[i];
    } //next i
}

return(dp); };

// Operator "crossProduct" ========================================================= //

// ---------------------------------------------------------------------------------- //
/*!
    Compute the cross product in R3 of input vectors, x and y.
    Template parameter can be any scalar type

    \param[in] x 1st argument
    \param[in] y 2nd argument

    \result on output returns the cross product product of x and y.
*/
template <class T>
std::vector<T> crossProduct(
    const std::vector<T>                        &x,
    const std::vector<T>                        &y
) {

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //
std::vector<T>      z(3,0.0);

// =================================================================================== //
// COMPUTE THE EXTERNAL PRODUCT                                                        //
// =================================================================================== //
z[0] = x[1] * y[2] - x[2] * y[1];
z[1] = x[2] * y[0] - x[0] * y[2];
z[2] = x[0] * y[1] - x[1] * y[0];

return (z);}
/*!
   @}
 */
