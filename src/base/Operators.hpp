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
#ifndef __OPERATORS_HH__
#define __OPERATORS_HH__

// ================================================================================== //
// INCLUDES                                                                           //
// ================================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <algorithm>
# include <functional> 
// CC_lib
// none

// ================================================================================== //
// NAMESPACES                                                                         //
// ================================================================================== //
using namespace std;

// ================================================================================== //
// FUNCTION PROTOTYPES                                                                //
// ================================================================================== //

// STL vectors ====================================================================== //

// Operator "+" --------------------------------------------------------------------- //
template <class T>
vector< T > operator+ (                                                               // ELEMENT-WISE SUM BETWEEN TWO VECTORS
        const vector< T >             &,                                                    // 1st argument (vector)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< T > operator+ (                                                               // ELEMENT-WISE SUM BETWEEN VECTOR AND CONSTANT
        const vector< T >             &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T>
vector< T > operator+ (                                                               // ELEMENT-WISE SUM BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator+ (                                                     // ELEMENT-WISE SUM BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< vector< T > >   &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator+ (                                                     // ELEMENT-WISE SUM BETWEEN CONSTANT AND VECTOR
        const vector< vector< T > >   &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "+=" -------------------------------------------------------------------- //
template <class T>
vector< T >& operator+= (                                                             // ELEMENT-WISE INCREMENT OF A VECTOR
        vector< T >                   &,                                                    // 1st argument (vector)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< T >& operator+= (                                                             // ELEMENT-WISE INCREMENT OF A VECTOR
        vector< T >                   &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T>
vector< vector< T > >& operator+= (                                                   // ELEMENT-WISE INCREMENT OF A VECTOR
        vector< vector< T > >         &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "-" --------------------------------------------------------------------- //
template <class T>
vector< T > operator- (                                                               // ELEMENT-WISE DIFFERENCE BETWEEN TWO VECTORS
        const vector< T >             &,                                                    // 1st argument (vector)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< T > operator- (                                                               // ELEMENT-WISE DIFFERENCE BETWEEN VECTOR AND CONSTANT
        const vector< T >             &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T>
vector< T > operator- (                                                               // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator- (                                                     // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< vector< T > >   &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator- (                                                     // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND VECTOR
        const vector< vector< T > >   &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "-=" -------------------------------------------------------------------- //
template <class T>
vector< T >& operator-= (                                                             // ELEMENT-WISE DECREMENT OF A VECTOR
        vector< T >                   &,                                                    // 1st argument (vector)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< T >& operator-= (                                                             // ELEMENT-WISE DECREMENT OF A VECTOR
        vector< T >                   &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T>
vector< vector< T > >& operator-= (                                                   // ELEMENT-WISE DECREMENT OF A VECTOR
        vector< vector< T > >         &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "*" --------------------------------------------------------------------- //
template <class T>
vector< T > operator* (                                                               // ELEMENT-WISE PRODUCT BETWEEN TWO VECTORS
        const vector< T >             &,                                                    // 1st argument (vector)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< T > operator* (                                                               // ELEMENT-WISE PRODUCT BETWEEN VECTOR AND CONSTANT
        const vector< T >             &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T>
vector< T > operator* (                                                               // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector < T > > operator* (                                                    // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< vector< T > >   &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator* (                                                     // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND VECTOR
        const vector< vector< T > >   &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "/" --------------------------------------------------------------------- //
template <class T>
vector< T > operator/ (                                                               // ELEMENT-WISE DIVISION BETWEEN VECTORS
        const vector< T >             &,                                                    // 1st argument (vector)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< T > operator/ (                                                               // ELEMENT-WISE DIVISION BETWEEN VECTOR AND CONSTANT
        const vector< T >             &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T>
vector< T > operator/ (                                                               // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< T >             &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator/ (                                                     // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND VECTOR
        const T                       &,                                                    // 1st argument (constant)
        const vector< vector< T > >   &                                                     // 2nd argument (vector)
        );
template <class T>
vector< vector< T > > operator/ (                                                     // ELEMENT-WISE DIVISION BETWEEN CONSTANT AND VECTOR
        const vector< vector< T > >   &,                                                    // 1st argument (vector)
        const T                       &                                                     // 2nd argument (constant)
        );

// Output operator ------------------------------------------------------------------ //
template <class T>
ostream& operator<< (                                                                 // INSERTION OPERATOR
        ostream                     &,                                                    // (input) output stream
        const vector< T >           &                                                     // (input) vector to be streamed
        );
template <class T>
ofstream& operator<< (                                                                // INSERTION OPERATOR
        ofstream                    &,                                                    // (input) output file stream
        const vector< T >           &                                                     // (input) vector to be streamed
        );

// Input operators ------------------------------------------------------------------ //
template <class T>
istream& operator>> (                                                                 // EXTRACTION OPERATOR
        istream                     &,                                                    // (input) input stream
        vector< T >                 &                                                     // (input/output) vector to be streamed
        );
template <class T>
ifstream& operator>> (                                                                // EXTRACTION OPERATOR
        ifstream                    &,                                                    // (input) input file stream
        vector< T >                 &                                                     // (input/output) vector to be streamed
        );

// Operator "min" ------------------------------------------------------------------- //
template <class T>
vector< T > min(                                                                      // RETURNS THE MINIMUM BETWEEN VECTOR X AND VECTOR Y
        const vector< T >           &,                                                    // (input) 1st argument of comparison (vector)
        const vector< T >           &                                                     // (input) 2nd argument of comparison (vector)
        );
template <class T>
vector< T > min(                                                                      // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const vector< T >           &,                                                    // (input) 1st argument of comparison (vector)
        const T                     &                                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T>
vector< T > min(                                                                      // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const T                     &,                                                    // (input) 1st argument of comparison (scalar)
        const vector< T >           &                                                     // (input) 2nd argument of comparison (vector)
        );
template <class T>
vector< vector< T > > min(                                                            // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const vector< vector< T > > &,                                                    // (input) 1st argument of comparison (vector)
        const T                     &                                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T>
vector< vector< T > > min(                                                            // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const T                     &,                                                    // (input) 1st argument of comparison (scalar)
        const vector< vector< T > > &                                                     // (input) 2nd argument of comparison (vector)
        );

// Operator "minval" ---------------------------------------------------------------- //
void inline minval(                                                                   // DUMMY ROUTINE FOR MINVAL SEARCH
        const double                &,                                                    // (dummy input) 1st argument
        double                      &                                                     // (dummy input/output) 2nd argument
        );
void inline minval(                                                                   // DUMMY ROUTINE FOR MINVAL SEARCH
        const int                   &,                                                    // (dummy input) 1st argument
        int                         &                                                     // (dummy input/output) 2nd argument
        );
template <class T, class T1>
void minval(                                                                          // RETURNS THE MINIMUM ELEMENT OF A VECTOR
        const vector<T>             &,                                                    // (input) input vector
        T1                          &                                                     // (input/output) minimum element
        );

// Operator "max" ------------------------------------------------------------------- //
template <class T>
vector< T > max(                                                                      // RETURNS THE MAXIMUM BETWEEN X AND Y
        const vector< T >           &,                                                    // (input) 1st argument of comparison
        const vector< T >           &                                                     // (input) 2nd argument of comparison
        );
template <class T>
vector< T > max(                                                                      // RETURNS THE MAXIMUM BETWEEN X AND Y
        const vector< T >           &,                                                    // (input) 1st argument of comparison
        const T                     &                                                     // (input) 2nd argument of comparison
        );
template <class T>
vector< T > max(                                                                      // RETURNS THE MAXIMUM BETWEEN X AND Y
        const T                     &,                                                    // (input) 1st argument of comparison
        const vector< T >           &                                                     // (input) 2nd argument of comparison
        );
template <class T>
vector< vector< T > > max(                                                            // RETURNS THE MAXIMUM BETWEEN X AND Y
        const vector< vector< T > > &,                                                    // (input) 1st argument of comparison
        const T                     &                                                     // (input) 2nd argument of comparison
        );
template <class T>
vector< vector< T > > max(                                                            // RETURNS THE MAXIMUM BETWEEN X AND Y
        const T                     &,                                                    // (input) 1st argument of comparison
        const vector< vector< T > > &                                                     // (input) 2nd argument of comparison
        );

// Operator "maxval" ---------------------------------------------------------------- //
void inline maxval(                                                                   // DUMMY ROUTINE FOR MAXVAL SEARCH
        const double                &,                                                    // (dummy input) 1st argument
        double                      &                                                     // (dummy input/output) 2nd argument
        );
void inline maxval(                                                                   // DUMMY ROUTINE FOR MAXVAL SEARCH
        const int                   &,                                                    // (dummy input) 1st argument
        int                         &                                                     // (dummy input/output) 2nd argument
        );
template <class T, class T1>
void maxval(                                                                          // RETURNS THE MAXIMUM ELEMENT OF A VECTOR
        const vector<T>             &,                                                    // (input) input vector
        T1                          &                                                     // (input/output) maximum element
        );

// Operator "sum" ------------------------------------------------------------------- //
void inline sum(                                                                      // DUMMY ROUTINE FOR SUM OPERATOR
        const int                   &,                                                    // (dummy input) 1st argument
        int                         &                                                     // (dummy input/output) 2nd argument
        );
void inline sum(                                                                      // DUMMY ROUTINE FOR SUM OPERATOR
        const double                &,                                                    // (dummy input) 1st argument
        double                      &                                                     // (dummy input/output) 2nd argument
        );
template <class T>
void inline sum(                                                                      // DUMMY ROUTINE FOR SUM OPERATOR
        const T                     &,                                                    // (dummy input) 1st argument
        T                           &                                                     // (dummy input/output) 2nd argument
        );
template <class T, class T1>
void sum(                                                                             // RETURNS THE SUM OF ELEMENT OF A VECTOR
        const vector< T >           &,                                                    // (input) input vector
        T1                          &                                                     // (input/output) sum of vector's element
        );

// Operator "abs" ------------------------------------------------------------------- //
template <class T>
vector<T> abs(                                                                        // RETURNS THE ABSOLUTE VALUE OF A VECTOR
        const vector< T >           &                                                     // (input) input vector
        );

// Operator "sign" ------------------------------------------------------------------- //
template <class T>
T sign(                                                                               // RETURNS THE SIGN 
        const  T                    &                                                     // (input) input value
      );

// Operator "pow" ------------------------------------------------------------------- //
template <class T>
vector< T > pow(                                                                      // RETURNS THE ELEMENTWISE POWER OF VECTOR
        vector< T >                 &,                                                    // (input) input vector
        double                                                                            // (input) power
        );

// Operator "norm" ------------------------------------------------------------------ //
template <class T>
double norm_1(                                                                        // RETURNS THE 1-NORM OF A VECTOR
        const vector< T >           &                                                     // (input) input vector
        );
template <class T>
double norm_2(                                                                        // RETURNS THE 2-NORM OF A VECTOR
        const vector< T >           &                                                     // (input) input vector
        );
template <class T>
double norm(                                                                          // RETURNS THE P-NORM OF A VECTOR
        const vector< T >           &,                                                    // (input) input vector
        int                                                                               // (input) norm index
        );
template <class T>
double norm_inf(                                                                      // RETURNS THE inf-NORM OF A VECTOR
        const vector< T >           &                                                     // (input) input vector
        );

// Operator "Dot_Product" ----------------------------------------------------------- //
template <class T>
T Dot_Product(                                                                        // COMPUTE THE DOT PRODUCT OF TWO VECTORS
        const vector< T >           &,                                                    // (input) 1st argument of dot product
        const vector< T >           &                                                     // (input) 2nd argument of dot product
        );

// Operator Cross_Product ----------------------------------------------------------- //
template <class T>
vector<T> Cross_Product(                                                              // COMPUTE THE CROSS-PRODUCT OF TWO VECTORS
        vector<T> const             &,                                                    // (input) 1st argument of cross product
        vector<T> const             &                                                     // (input) 2nd argument of cross product
        );

// Various -------------------------------------------------------------------------- //
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
ostream& display(                                                                     // DUMMY ROUTINE FOR RECURSIVE TEMPLATED FUNCTION display
        ostream                     &,                                                    // (input/output) output stream
        const T                     &                                                     // (input) vector to be displayed
        );
template<class T>
ostream& display(                                                                     // DISPLAY VECTOR IN A NICELY FORMATTED FORM
        ostream                     &,                                                    // (input/output) output stream
        const vector<T>             &,                                                    // (input) vector to be displayed
        unsigned int                 padding = 0                                          // (input/optional) number of trailing spaces
        );
template<typename T, typename std::enable_if< std::is_scalar< T >::value>::type* = nullptr>
ofstream& display(                                                                    // DUMMY ROUTINE FOR RECURSIVE TEMPLATED FUNCTION display
        ofstream                    &,                                                    // (input/output) output stream
        const T                     &                                                     // (input) vector to be displayed
        );
template<class T>
ofstream& display(                                                                    // DISPLAY VECTOR IN A NICELY FORMATTED FORM
        ofstream                    &,                                                    // (input/output) output stream
        const vector<T>             &,                                                    // (input) vector to be displayed
        unsigned int                 padding = 0                                          // (input/optional) number of trailing spaces
        );

// C++ v10.0 arrays ================================================================= //

// Operator "+" --------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> operator+ (                                                               // ELEMENT-WISE SUM BETWEEN ARRAYS
        const array<T, d>             &,                                                    // 1st argument (array)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d>
array<T, d> operator+ (                                                               // ELEMENT-WISE SUM BETWEEN ARRAY AND CONSTANT
        const array<T, d>             &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T, size_t d>
array<T, d> operator+ (                                                               // ELEMENT-WISE SUM BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e> 
array<array<T, e>, d> operator+ (                                                     // ELEMENT-WISE SUM BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<array<T, e>, d>   &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator+ (                                                     // ELEMENT-WISE SUM BETWEEN CONSTANT AND ARRAY
        const array<array<T, e>, d>   &,                                                    // 1st argument (constant)
        const T                       &                                                     // 2nd argument (array)
        );

// Operator "+=" -------------------------------------------------------------------- //
template <class T, size_t d>
array< T, d >& operator+= (                                                           // ELEMENT-WISE INCREMENT OF A ARRAY
        array< T, d >                 &,                                                    // 1st argument (array)
        const array< T, d >           &                                                     // 2nd argument (vector)
        );
template <class T, size_t d>
array< T, d >& operator+= (                                                           // ELEMENT-WISE INCREMENT OF A ARRAY
        array< T, d >                 &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T, size_t d, size_t e>
array< array< T, e >, d >& operator+= (                                               // ELEMENT-WISE INCREMENT OF A ARRAY
        array< array< T, e >, d >     &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "-" --------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> operator- (                                                               // ELEMENT-WISE DIFFERENCE BETWEEN TWO ARRAYS
        const array<T, d>             &,                                                    // 1st argument (array)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d>
array<T, d> operator- (                                                               // ELEMENT-WISE DIFFERENCE BETWEEN ARRAY AND CONSTANT
        const array<T, d>             &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T, size_t d>
array<T, d> operator- (                                                               // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator- (                                                     // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<array<T, e>, d>   &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator- (                                                     // ELEMENT-WISE DIFFERENCE BETWEEN CONSTANT AND ARRAY
        const array<array<T, e>, d>   &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "-=" -------------------------------------------------------------------- //
template <class T, size_t d>
array< T, d >& operator-= (                                                           // ELEMENT-WISE DECREMENT OF A ARRAY
        array< T, d >                 &,                                                    // 1st argument (array)
        const array< T, d >           &                                                     // 2nd argument (vector)
        );
template <class T, size_t d>
array< T, d >& operator-= (                                                           // ELEMENT-WISE DECREMENT OF A ARRAY
        array< T, d >                 &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T, size_t d, size_t e>
array< array< T, e >, d >& operator-= (                                               // ELEMENT-WISE DECREMENT OF A ARRAY
        array< array< T, e >, d >     &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "*" --------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> operator* (                                                               // ELEMENT-WISE PRODUCT BETWEEN TWO ARRAYS
        const array<T, d>             &,                                                    // 1st argument (array)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d>
array<T, d> operator* (                                                               // ELEMENT-WISE PRODUCT BETWEEN ARRAY AND CONSTANT
        const array<T, d>             &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T, size_t d>
array<T, d> operator* (                                                               // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator* (                                                     // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<array<T, e>, d>   &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator* (                                                     // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const array<array<T, e>, d>   &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );

// Operator "/" --------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> operator/ (                                                               // ELEMENT-WISE PRODUCT BETWEEN ARRAYS
        const array<T, d>             &,                                                    // 1st argument (array)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d>
array<T, d> operator/ (                                                               // ELEMENT-WISE PRODUCT BETWEEN ARRAY AND CONSTANT
        const array<T, d>             &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );
template <class T, size_t d>
array<T, d> operator/ (                                                               // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<T, d>             &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator/ (                                                     // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const T                       &,                                                    // 1st argument (constant)
        const array<array<T, e>, d>   &                                                     // 2nd argument (array)
        );
template <class T, size_t d, size_t e>
array<array<T, e>, d> operator/ (                                                     // ELEMENT-WISE PRODUCT BETWEEN CONSTANT AND ARRAY
        const array<array<T, e>, d>   &,                                                    // 1st argument (array)
        const T                       &                                                     // 2nd argument (constant)
        );

// Output operator ------------------------------------------------------------------ //
template <class T, size_t d>
ostream& operator<< (                                                                 // INSERTION OPERATOR
        ostream                     &,                                                    // (input) output stream
        const array<T, d>           &                                                     // (input) vector to be streamed
        );
template <class T, size_t d>
ofstream& operator<< (                                                                // INSERTION OPERATOR
        ofstream                    &,                                                    // (input) output file stream
        const array<T, d>           &                                                     // (input) vector to be streamed
        );

// Input operators ------------------------------------------------------------------ //
template <class T, size_t d>
istream& operator>> (                                                                 // EXTRACTION OPERATOR
        istream                    &,                                                     // (input) input stream
        array< T, d >              &                                                     // (input/output) vector to be streamed
        );
template <class T, size_t d>
ifstream& operator>> (                                                                // EXTRACTION OPERATOR
        ifstream                    &,                                                    // (input) input file stream
        array<T, d>                 &                                                     // (input) vector to be streamed
        );

// Operator "min" ------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> min(                                                                      // RETURNS THE MINIMUM BETWEEN VECTOR X AND VECTOR Y
        const array<T, d>           &,                                                    // (input) 1st argument of comparison (vector)
        const array<T, d>           &                                                     // (input) 2nd argument of comparison (vector)
        );
template <class T, size_t d>
array<T, d> min(                                                                      // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const array<T, d>           &,                                                    // (input) 1st argument of comparison (vector)
        const T                     &                                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d>
array<T, d> min(                                                                      // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const T                     &,                                                    // (input) 1st argument of comparison (scalar)
        const array<T, d>           &                                                     // (input) 2nd argument of comparison (vector)
        );
template <class T, size_t d, size_t n>
array<array<T, n>, d> min(                                                            // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const array<array<T, n>, d> &,                                                    // (input) 1st argument of comparison (vector)
        const T                     &                                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d, size_t n>
array<array<T, n>, d> min(                                                            // RETURNS THE MINIMUM BETWEEN VECTOR X AND SCALAR Y
        const T                     &,                                                    // (input) 1st argument of comparison (scalar)
        const array<array<T, n>, d> &                                                     // (input) 2nd argument of comparison (vector)
        );

// Operator "minval" ---------------------------------------------------------------- //
template <class T, size_t d, class T1>
void minval(                                                                          // RETURNS THE MINIMUM ELEMENT OF A VECTOR
        const array<T, d>           &,                                                    // (input) input vector
        T1                          &                                                     // (input/output) minimum element
        );

// Operator "max" ------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> max(                                                                      // RETURNS THE MAXIMUM BETWEEN TWO ARRAYS
        const array<T, d>           &,                                                    // (input) 1st argument of comparison (array)
        const array<T, d>           &                                                     // (input) 2nd argument of comparison (array)
        );
template <class T, size_t d>
array<T, d> max(                                                                      // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const array<T, d>           &,                                                    // (input) 1st argument of comparison (array)
        const T                     &                                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d>
array<T, d> max(                                                                      // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const T                     &,                                                    // (input) 1st argument of comparison (scalar)
        const array<T, d>           &                                                     // (input) 2nd argument of comparison (array)
        );
template <class T, size_t d, size_t n>
array<array<T, n>, d> max(                                                            // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const array<array<T, n>, d> &,                                                    // (input) 1st argument of comparison (array)
        const T                     &                                                     // (input) 2nd argument of comparison (scalar)
        );
template <class T, size_t d, size_t n>
array<array<T, n>, d> max(                                                            // RETURNS THE MAXIMUM BETWEEN ARRAY AND SCALAR
        const T                     &,                                                    // (input) 1st argument of comparison (scalar)
        const array<array<T, n>, d> &                                                     // (input) 2nd argument of comparison (array)
        );

// Operator "maxval" ---------------------------------------------------------------- //
template <class T, size_t d, class T1>
void maxval(                                                                          // RETURNS THE MAXIMUM ELEMENT OF A ARRAY
        const array<T, d>           &,                                                    // (input) input array
        T1                          &                                                     // (input/output) maximum element
        );

// Operator "sum" ------------------------------------------------------------------- //
template <class T, size_t d, class T1>
void sum(                                                                             // RETURNS THE SUM OF ARRAY ELEMENTS
        const array<T, d>           &,                                                    // (input) input array
        T1                          &                                                     // (input/output) sum of array's elements
        );

// Operator "abs" ------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> abs(                                                                      // RETURNS THE ABSOLUTE VALUE OF A ARRAY
        const array<T, d>           &                                                     // (input) input array
        );

// Operator "pow" ------------------------------------------------------------------- //
template <class T, size_t d>
array<T, d> pow(                                                                      // RETURNS THE ELEMENTWISE POWER OF ARRAY
        array<T, d>                 &,                                                    // (input) input array
        double                                                                            // (input) power
        );

// Operator "norm" ------------------------------------------------------------------ //
template <class T, size_t d>
double norm_1(                                                                        // RETURNS THE 1-NORM OF ARRAY
        const array<T, d>  &x                                                             // (input) input array
        );
template <class T, size_t d>
double norm_2(                                                                        // RETURNS THE 2-NORM OF ARRAY
        const array<T, d>  &x                                                             // (input) input array
        );
template <class T, size_t d>
double norm(                                                                          // RETURNS THE p-NORM OF ARRAY
        const array<T, d>  &x,                                                            // (input) input array
        int                 p                                                             // (input) norm index
        );
template <class T, size_t d>
double norm_inf(                                                                      // RETURNS THE inf NORM OF ARRAY
        const array<T, d>  &x                                                             // (input) input array
        );

// Operator "Dot_Product" ----------------------------------------------------------- //
template <class T, size_t d>
T Dot_Product(                                                                        // COMPUTE THE DOT PRODUCT OF TWO ARRAYS
        const array<T, d>           &,                                                    // (input) 1st argument of dot product
        const array<T, d>           &                                                     // (input) 2nd argument of dot product
        );

// Operator Cross_Product ----------------------------------------------------------- //
template <class T, size_t d>
T Cross_Product(                                                                      // COMPUTE THE CROSS-PRODUCT OF TWO ARRAYS
        array<T, 2> const           &,                                                    // (input) 1st argument of cross product
        array<T, 2> const           &                                                     // (input) 2nd argument of cross product
        );

template <class T, size_t d>
array<T, 3> Cross_Product(                                                            // COMPUTE THE CROSS-PRODUCT OF TWO ARRAYS
        array<T, 3> const           &,                                                    // (input) 1st argument of cross product
        array<T, 3> const           &                                                     // (input) 2nd argument of cross product
        );

// Various -------------------------------------------------------------------------- //
template<class T, size_t d>
ostream& display(                                                                     // DISPLAY ARRAY IN A NICELY FORMATTED FORM
        ostream                     &,                                                    // (input/output) output stream
        const array<T, d>           &,                                                    // (input) array to be displayed
        unsigned int                 padding = 0                                          // (input/optional) number of trailing spaces
        );
template<class T, size_t d>
ofstream& display(                                                                    // DISPLAY ARRAY IN A NICELY FORMATTED FORM
        ofstream                    &,                                                    // (input/output) output stream
        const array<T, d>           &,                                                    // (input) array to be displayed
        unsigned int                 padding = 0                                          // (input/optional) number of trailing spaces
        );

// String =========================================================================== //

// Trimming operators --------------------------------------------------------------- //
static inline string &ltrim(                                                          // STRING LEFT TRIMMING
        string                      &                                                     // (input) string to be trimmed
        );
static inline string &rtrim(                                                          // STRING RIGHT TRIMMING
        string                      &                                                     // (input) string to be trimmed
        );
static inline string &trim(                                                           // STRING TRIMMING
        string                      &                                                     // (input) string to be trimmed
        );

// Padding operators ---------------------------------------------------------------- //
static inline string ZeroPadNumber(                                                   // PERFORMS CONVERSION OF INTEGER INTO STRING
        int                          ,                                                    // (input) number of char in string
        int                                                                               // (input) integer to be padded
        );

// Input stream operator ------------------------------------------------------------ //
bool Get_After_Keyword( 							      // EXTRACT FIELD AFTER SPECIFIC KEYWORD
        string line_, 								      // (input) string
        string key_, 								      // (input) keyword
        char del_,   								      // (input) field delimiter
        string& result_						                      // (output) field found
        );

// returns true if key_ is present in line ------------------------------------------ //
static inline bool Keyword_In_String( 
        string line_, 
        string key_
        ) ;

// converts a string to fundamental data types and vectors or arrays of them -------- //
template <class T>
void  convert_string( 
        string input_, 
        T &output_
        );

template <class T>
void  convert_string( 
        string input_, 
        vector<T> &output_
        );

template <class T, size_t n>
void  convert_string( 
        string input_, 
        array<T,n> &output_
        );

// ================================================================================== //
// TEMPLATES                                                                          //
// ================================================================================== //
# include "Operators.tpp"
# include "Operators_vector.tpp"
# include "MathOperators_vector.tpp"
# include "Operators_array.tpp"
# include "MathOperators_array.tpp"
# include "Operators_string.tpp"

#endif
