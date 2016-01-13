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
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// ================================================================================== //
// TEMPLATES IMPLEMENTATION                                                           //
// ================================================================================== //

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
array<T, d> operator+ (
  const array<T, d> &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator+ (                                                            //
//   const array<T, d> &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise sum of C++ v10.0 arrays. Returns:                                     //
//        z = x + y, s.t. z[i] = x[i] + y[i]                                          //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '+' operator                                  //
// - y   : array<T, d>, 2nd argument of '+' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, sum of x, y                                                   //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator+ (
  const array<T, d> &x,
  const T           &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator+ (                                                            //
//   const array<T, d> &x,                                                            //
//   const T           &y)                                                            //
//                                                                                    //
// Element-wise sum between C++ v10.0 array and constant. Returns:                    //
//        z = x + y, s.t. z[i] = x[i] + y                                             //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '+' operator                                  //
// - y   : T, 2nd argument of '+' operator                                            //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, sum of x, y                                                   //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator+ (
  const T           &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator+ (                                                            //
//   const T           &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise sum between constant and C++ v10.0 array. Retunrs:                    //
//        z = x + y, s.t. z[i] = x + y[i]                                             //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - y   : T, 1st argument of '+' operator                                            //
// - x   : array<T, d>, 2nd argument of '+' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, sum of x, y                                                   //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>   z;

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
array<array<T, e>, d> operator+ (
  const T                       &x,
  const array<array<T, e>, d>   &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator+ (                                                  //
//   const T                       &x,                                                //
//   const array<array<T, e>, d>   &y)                                                //
//                                                                                    //
// Element-wise sum between constant and C++ v10.0 array. Retunrs:                    //
//        z = x + y, s.t. z[i][j] = x + y[i][j]                                       //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : T, 1st argument of '+' operator                                            //
// - y   : array<array<T, e>, d>, 1st argument of '+' operator                        //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, sum of x, y                                                   //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>   z;

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<array<T, e>, d> operator+ (
  const array<array<T, e>, d>   &x,
  const T                       &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator+ (                                                  //
//   const T                       &x,                                                //
//   const array<array<T, e>, d>   &y)                                                //
//                                                                                    //
// Element-wise sum between constant and C++ v10.0 array. Retunrs:                    //
//        z = x + y, s.t. z[i][j] = x[i][j] + y                                       //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<array<T, e>, d>, 1st argument of '+' operator                        //
// - y   : T, 2nd argument of '+' operator                                            //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<array<T, e>, d>, sum of x, y                                         //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>   z;

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
array< T, d >& operator+= (
  array< T, d >                 &x,
  const array< T, d >           &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array< T, d >& operator+= (                                                        //
//   array< T, d >                 &x,                                                //
//   const array< T, d >           &y)                                                //
//                                                                                    //
// Element-wise increment. Returns:                                                   //
//      x += y, s.t. x[i] += y[i]                                                     //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array< T, d >, 1st argument of '+=' operator                               //
// - y   : array< T, d >, 2nd argument of '+=' operator                               //
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
for (int i = 0; i < d; i++){
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
array< T, d >& operator+= (
  array< T, d >                 &x,
  const T                       &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array< T, d >& operator+= (                                                        //
//   array< T, d >                 &x,                                                //
//   const T                       &y)                                                //
//                                                                                    //
// Element-wise increment. Returns:                                                   //
//        x += y, s.t x[i] += y                                                       //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array< T, d >, 1st argument of '+=' operator                               //
// - y   : T          , 2nd argument of '+=' operator                                 //
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
for (int i = 0; i < d; i++){
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
array< array< T, e >, d >& operator+= (
    array< array< T, e >, d >       &x,
    const T                         &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array< array< T, e >, d >& operator+= (                                            //
//     array< array< T, e >, d >       &x,                                            //
//     const T                         &y)                                            //
//                                                                                    //
// Element-wise increment. Returns:                                                   //
//     x += y, s.t. x[i][j] += y                                                      //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array< array< T, e >, d >, 1st argument of '+=' operator                   //
// - y   : T, 2nd argument of '+=' operator                                           //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - x   : array< array< T, e >, d >&, reference to first argument                    //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (int i = 0; i < d; i++) {
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
array<T, d> operator- (
  const array<T, d> &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator- (                                                            //
//   const array<T, d> &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise difference between C++ v10.0 arrays. Returns:                         //
//     z = x - y, s.t. z[i] = x[i] - y[i]                                             //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '-' operator                                  //
// - y   : array<T, d>, 2nd argument of '-' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, difference of x, y                                            //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator- (
  const array<T, d> &x,
  const T           &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator- (                                                            //
//   const array<T, d> &x,                                                            //
//   const T           &y)                                                            //
//                                                                                    //
// Element-wise difference between C++ v10.0 array and constant. Returns:             //
//     z = x - y, s.t. z[i] = x[i] - y                                                //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '-' operator                                  //
// - y   : T, 2nd argument of '-' operator                                            //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, difference of x, y                                            //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator- (
  const T           &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator- (                                                            //
//   const T           &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise difference between constant and C++ v10.0 array. Returns:             //
//     z = x - y, s.t. z[i] = x - y[i]                                                //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : T, 1st argument of '-' operator                                            //
// - y   : array<T, d>, 2nd argument of '-' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, difference of x, y                                            //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>   z;

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
array<array<T, e>, d> operator- (
  const T                       &x,
  const array<array<T, e>, d>   &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator- (                                                  //
//   const T                       &x,                                                //
//   const array<array<T, e>, d>   &y)                                                //
//                                                                                    //
// Element-wise difference between constant and C++ v10.0 array. Returns:             //
//     z = x - y, s.t. z[i][j] = x - y[i][j]                                          //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : T, 1st argument of '-' operator                                            //
// - y   : array<array<T, e>, d>, 2nd argument of '-' operator                        //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<array<T, e>, d>, difference of x, y                                  //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>   z;

// ================================================================================== //
// PERFORM DIFFERENCE                                                                 //
// ================================================================================== //
for (int i = 0; i < d; i++) {
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
array<array<T, e>, d> operator- (
  const array<array<T, e>, d>   &x,
  const T                       &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator- (                                                  //
//   const array<array<T, e>, d>   &x,                                                //
//   const T                       &y)                                                //
//                                                                                    //
// Element-wise difference between constant and C++ v10.0 array. Returns:             //
//     z = x - y, s.t. z[i][j] = x[i][j] - y                                          //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<array<T, e>, d>, 1st argument of '-' operator                        //
// - y   : T, 2nd argument of '-' operator                                            //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<array<T, e>, d>, difference of x, y                                  //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>   z;

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
array< T, d >& operator-= (
  array< T, d >                 &x,
  const array< T, d >           &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array< T, d >& operator-= (                                                        //
//   array< T, d >                 &x,                                                //
//   const array< T, d >           &y)                                                //
//                                                                                    //
// Element-wise decrement. Returns:                                                   //
//      x -= y, s.t. x[i] -= y[i]                                                     //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array< T, d >, 1st argument of '-=' operator                               //
// - y   : array< T, d >, 2nd argument of '-=' operator                               //
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
for (int i = 0; i < d; i++){
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
array< T, d >& operator-= (
  array< T, d >                 &x,
  const T                       &y
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
for (int i = 0; i < d; i++){
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
array< array< T, e >, d >& operator-= (
    array< array< T, e >, d >       &x,
    const T                         &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array< array< T, e >, d >& operator-= (                                            //
//     array< array< T, e >, d >       &x,                                            //
//     const T                         &y)                                            //
//                                                                                    //
// Element-wise decrement. Returns:                                                   //
//     x -= y, s.t. x[i][j] -= y                                                      //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array< array< T, e >, d >, 1st argument of '-=' operator                   //
// - y   : T, 2nd argument of '-=' operator                                           //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - x   : array< array< T, e >, d >&, reference to first argument                    //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
// none

// ================================================================================== //
// PERFORM SUM                                                                        //
// ================================================================================== //
for (int i = 0; i < d; i++) {
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
array<T, d> operator* (
  const array<T, d> &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator* (                                                            //
//   const array<T, d> &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise product between C++ v10.0 arrays. Returns:                            //
//      z = x * y, s.t. z[i] = x[i]*y[i]                                              //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '*' operator                                  //
// - y   : array<T, d>, 2nd argument of '*' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, elementwise product between x, y                              //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator* (
  const array<T, d> &x,
  const T           &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator* (                                                            //
//   const array<T, d> &x,                                                            //
//   const T           &y)                                                            //
//                                                                                    //
// Element-wise product between C++ v10.0 array and constant. Returns:                //
//      z = x * y, s.t. z[i] = x[i]*y                                                 //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '*' operator                                  //
// - y   : T, 2nd argument of '*' operator                                            //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, elementwise product between x, y                              //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator* (
  const T           &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator* (                                                            //
//   const T           &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise product between constant and C++ v10.0 array. Returns:                //
//      z = x * y, s.t. z[i] = x*y[i]                                                 //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : T, 1st argument of '*' operator                                            //
// - y   : array<T, d>, 2nd argument of '*' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, elementwise product between x, y                              //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>   z;

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
array<array<T, e>, d> operator* (
  const T                       &x,
  const array<array<T, e>, d>   &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator* (                                                  //
//   const T                       &x,                                                //
//   const array<array<T, e>, d>   &y)                                                //
//                                                                                    //
// Element-wise product between constant and C++ v10.0 array. Returns:                //
//      z = x * y, s.t. z[i][j] = x*y[i][j]                                           //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : T, 1st argument of '*' operator                                            //
// - y   : array<array<T, e>, d>, 2nd argument of '*' operator                        //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<array<T, e>, d>, elementwise product between x, y                    //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>   z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (int i = 0; i < d; i++) {
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
array<array<T, e>, d> operator* (
  const array<array<T, e>, d>   &x,
  const T                       &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator* (                                                  //
//   const array<array<T, e>, d>   &x,                                                //
//   const T                       &y)                                                //
//                                                                                    //
// Element-wise product between constant and C++ v10.0 array. Returns:                //
//      z = x * y, s.t. z[i][j] = x[i][j]*y                                           //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<array<T, e>, d>, 1st argument of '*' operator                        //
// - y   : T, 2nd argument of '*' operator                                            //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<array<T, e>, d>, elementwise product between x, y                    //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>   z;

// ================================================================================== //
// PERFORM PRODUCT                                                                    //
// ================================================================================== //
for (int i = 0; i < d; i++) {
    z[i] = x[i] * y;
} //next i

return(z); };

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
array<T, d> operator/ (
  const array<T, d> &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator/ (                                                            //
//   const array<T, d> &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise division between C++ v10.0 arrays. Retuns:                            //
//    z = x / y, s.t. z[i] = x[i]/y[i].                                               //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x   : array<T, d>, 1st argument of '/' operator                                  //
// - y   : array<T, d>, 2nd argument of '/' operator                                  //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z   : array<T, d>, elementwise division between x, y                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>          z;

// ================================================================================== //
// PERFORM DIVISION                                                                   //
// ================================================================================== //
for (int i = 0; i < d; i++){
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
array<T, d> operator/ (
  const array<T, d> &x,
  const T           &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator/ (                                                            //
//   const array<T, d> &x,                                                            //
//   const T           &y)                                                            //
//                                                                                    //
// Element-wise division between C++ v10.0 array and constant. Retuns:                //
//    z = x / y, s.t. z[i] = x[i]/y.                                                  //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x    : array<T, d>, 1st argument of "/" operator                                 //
// - y    : T, 2nd argument of "/" operator                                           //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z    : array<T, d>>, result of element-wise division                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
array<T, d>    z;

// Counters
int            i;

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
array<T, d> operator/ (
  const T           &x,
  const array<T, d> &y
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// array<T, d> operator/ (                                                            //
//   const T           &x,                                                            //
//   const array<T, d> &y)                                                            //
//                                                                                    //
// Element-wise division between constant and C++ v10.0 array. Retuns:                //
//    z = x / y, s.t. z[i] = x/y[i].                                                  //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x    : T, 1st argument of "/" operator                                           //
// - y    : array<T, d>, 2nd argument of "/" operator                                 //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z    : array<T, d>, result of element-wise division                              //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<T, d>      z;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (int i = 0; i < d; i++) {
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
array<array<T, e>, d> operator/ (
  const T                       &x,
  const array<array<T, e>, d>   &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator/ (                                                  //
//   const T                       &x,                                                //
//   const array<array<T, e>, d>   &y)                                                //
//                                                                                    //
// Element-wise division between constant and C++ v10.0 array. Retuns:                //
//    z = x / y, s.t. z[i][j] = x/y[i][j].                                            //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x    : T, 1st argument of "/" operator                                           //
// - y    : array<array<T, e>, d>, 2nd argument of "/" operator                       //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z    : array<array<T, e>, d>, result of element-wise division                    //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>      z;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (int i = 0; i < d; i++) {
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
array<array<T, e>, d> operator/ (
  const array<array<T, e>, d>   &x,
  const T                       &y
) {

// ================================================================================== //
// template <class T, size_t d, size_t e>                                             //
// array<array<T, e>, d> operator/ (                                                  //
//   const array<array<T, e>, d>   &x,                                                //
//   const T                       &y)                                                //
//                                                                                    //
// Element-wise division between constant and C++ v10.0 array. Retuns:                //
//    z = x / y, s.t. z[i][j] = x[i][j]/y.                                            //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - x    : T, 1st argument of "/" operator                                           //
// - y    : array<array<T, e>, d>, 2nd argument of "/" operator                       //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - z    : array<array<T, e>, d>, result of element-wise division                    //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
array<array<T, e>, d>      z;

// ================================================================================== //
// DIVIDE X BY Y                                                                      //
// ================================================================================== //
for (int i = 0; i < d; i++) {
    z[i] = x[i]/y;
}

return(z); };

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
ostream& operator<< (
    ostream              &out,
    const array<T, d>    &x
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// ostream& operator<< (                                                              //
//     ostream              &out,                                                     //
//     const array<T, d>    &x)                                                       //
//                                                                                    //
// Output stream for array.                                                           //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - out       : ostream, with output stream                                          //
// - x         : array<T, d>, with array to be streamed                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - out       : ostream, with updated output stream                                  //
// ================================================================================== //

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
for (int i = 0; i < d-1; i++) {
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
ofstream& operator<< (
    ofstream             &out,
    const array<T, d>    &x
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// ofstream& operator<< (                                                             //
//     ofstream              &out,                                                    //
//     const array<T, d>    &x)                                                       //
//                                                                                    //
// Output file stream for array.                                                      //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - out       : ofstream, with output stream                                         //
// - x         : array<T, d>, with array to be streamed                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - out       : ofstream, with updated output stream                                 //
// ================================================================================== //

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
for (int i = 0; i < d-1; i++) {
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

    \param[in,out] out input stream
    \param[in,out] x argument of extraction operator

    \result reference to input stream (allows concatenation)
*/
template <class T, size_t d>
istream& operator>> (
    istream              &in,
    array<T, d>          &x
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// istream& operator>> (                                                              //
//     istream              &in,                                                      //
//     array<T, d>          &x)                                                       //
//                                                                                    //
// Input stream for array.                                                            //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - in        : istream, with input stream                                           //
// - x         : array<T, d>, with array to be streamed                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - x         : array<T, d>, updated streamed array.                                 //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
T                    dummy;

// Counters
int                  i;

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

    \param[in,out] out input file stream
    \param[in,out] x argument of extraction operator

    \result reference to input file stream (allows concatenation)
*/
template <class T, size_t d>
ifstream& operator>> (
    ifstream             &in,
    array<T, d>          &x
) {

// ================================================================================== //
// template <class T, size_t d>                                                       //
// ifstream& operator>> (                                                             //
//     ifstream             &in,                                                      //
//     array<T, d>          &x)                                                       //
//                                                                                    //
// Input file stream for array.                                                       //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - in        : ifstream, with input file stream                                     //
// - x         : array<T, d>, with array to be streamed                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - x         : array<T, d>, with updated streamed array                             //
// ================================================================================== //

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
ostream& display(
    ostream             &out,
    const array<T, d>   &x,
    unsigned int         padding
) {

// ================================================================================== //
// template<class T, size_t d>                                                        //
// ostream& display(                                                                  //
//     ostream             &out,                                                      //
//     const array<T, d>   &x,                                                        //
//     unsigned int         padding)                                                  //
//                                                                                    //
// Display array in a nicely formatted form.                                          //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - out      : ostream, output stream                                                //
// - x        : array<T, d>, array to be displayed                                    //
// - padding  : unsigned int (default = 0), number of trailing spaces                 //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - out      : ostream&, reference to output stream                                  //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Counters
typename array<T, d>::const_iterator      i, e = x.cend();

// ================================================================================== //
// DISPLAY VECTOR                                                                     //
// ================================================================================== //
if (x.size() == 0) {
    out << "[ ]";
    return(out);
}
out << string(padding, ' ') << "[ ";
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
ofstream& display(
    ofstream            &out,
    const array<T, d>   &x,
    unsigned int         padding
) {

// ================================================================================== //
// template<class T, size_t d>                                                        //
// ofstream& display(                                                                 //
//     ofstream            &out,                                                      //
//     const array<T, d>   &x,                                                        //
//     unsigned int         padding)                                                  //
//                                                                                    //
// Display array in a nicely formatted form.                                          //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - out      : ofstream, output stream                                               //
// - x        : array<T, d>, array to be displayed                                    //
// - padding  : unsigned int (default = 0), number of trailing spaces                 //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - out      : ofstream&, reference to output stream                                 //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Counters
typename array<T, d>::const_iterator      i, e = x.cend();

// ================================================================================== //
// DISPLAY VECTOR                                                                     //
// ================================================================================== //
if (x.size() == 0) {
    out << "[ ]";
    return(out);
}
out << string(padding, ' ') << "[ ";
--e;
for (i = x.begin(); i != e; ++i) {
    display(out, *i) << ", ";
} //next i
display(out, *e) << " ]";

return(out); }




