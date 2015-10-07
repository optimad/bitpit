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
# ifndef __LINEAR_ALGEBRA_HH__
# define __LINEAR_ALGEBRA_HH__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard template library
# include <vector>
# include <cmath>
# include <array>
# include <iostream>

// CC_lib
# include "Operators.hpp"

// Others
// none

// ========================================================================== //
// NAMESAPCES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //

// boolean vectors
typedef vector< bool >                 bvector1D;
typedef vector< bvector1D >            bvector2D;
typedef vector< bvector2D >            bvector3D;
typedef vector< bvector3D >            bvector4D;

// characters vectors
typedef vector< char >                 cvector1D;
typedef vector< cvector1D >            cvector2D;
typedef vector< cvector2D >            cvector3D;
typedef vector< cvector3D >            cvector4D;

// integer vectors
typedef vector< int >                  ivector1D;
typedef vector< ivector1D >            ivector2D;
typedef vector< ivector2D >            ivector3D;
typedef vector< ivector3D >            ivector4D;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

// ========================================================================== //
// FUNCTIONS PROTOTYPES                                                       //
// ========================================================================== //

// Generic routines --------------------------------------------------------- //
template <class T>
void display_matrix(                                                          // Display matrix in nicely formatted output
    ostream                     &,                                            // (input/output) handle to output stream
    vector< vector< T > >       &                                             // (input) matrix to be displayed
); 
template <class T, size_t m, size_t n>
void display_matrix(                                                          // Display matrix in nicely formatted output
    ostream                     &,                                            // (input/output) handle to output stream
    array< array< T, n >, m >   &                                             // (input) matrix to be displayed
);

// Matrix Basic templates --------------------------------------------------- //
template <class T>
void zeros(                                                                   // Initialize an m-by-n matrix of zeros
    vector< vector < T > >      &,                                            // (input/output) matrix
    int                          ,                                            // (input) number of rows
    int                                                                       // (input) number of columns
);
template <class T, size_t m, size_t n>
void zeros(                                                                   // Initialize an m-by-n matrix of zeros
    array< array < T, n >, m >  &                                             // (input/output) matrix
);
template <class T>
void ones(                                                                    // Initialize an m-by-n matrix of ones
    vector< vector < T > >      &,                                            // (input/output) matrix
    int                          ,                                            // (input) number of rows
    int                                                                       // (input) number of columns
);
template <class T, size_t m, size_t n>
void ones(                                                                    // Initialize an m-by-n matrix of ones
    array< array < T, n >, m >  &                                             // (input/output) matrix
);
template <class T>
void eye(                                                                     // Initialize and m-by-n identity matrix
   vector< vector < T > >       &,                                            // (input/output) matrix
    int                          ,                                            // (input) number of rows
    int                                                                       // (input) number of columns
);
template <class T, size_t m, size_t n>
void eye(                                                                     // Initialize and m-by-n identity matrix
    array< array < T, n >, m >  &                                             // (input/output) matrix
);

// Matrix multiplications --------------------------------------------------- //
template <class T>
void matmul(                                                                  // Matrix multiplication with a constant
    T                            ,                                            // (input) 1st input scalar
    vector< vector< T > >       &,                                            // (input) 2nd input matrix
    vector< vector< T > >       &                                             // (input/output) output matrix
);
template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a constant
    T                            ,                                            // (input) 1st input scalar
    array< array< T, n >, m >   &,                                            // (input) 2nd input matrix
    array< array< T, n >, m >   &                                             // (input/output) output matrix
);
template <class T>
void matmul(                                                                  // Matrix multiplication with a constant
    vector< vector< T > >       &,                                            // (input) 1st input matrix
    T                            ,                                            // (input) 2nd input scalar
    vector< vector< T > >       &                                             // (input/output) output matrix
);
template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a constant
    array< array< T, n >, m >   &,                                            // (input) 1st input matrix
    T                            ,                                            // (input) 2nd input scalar
    array< array< T, n >, m >   &                                             // (input/output) output matrix
);
template <class T>
void matmul(                                                                  // Matrix multiplication with a vector
    vector< T >                 &,                                            // (input) 1st input vector
    vector< vector< T > >       &,                                            // (input) 2nd input matrix
    vector< T >                 &                                             // (input/output) output vector
);
template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a vector
    array< T, m >               &,                                            // (input) 1st input vector
    array< array < T, n >, m >  &,                                            // (input) 2nd input matrix
    array< T, n >               &                                             // (input/output) output vector
);
template <class T>
void matmul(                                                                  // Matrix multiplication with a vector
    vector< vector< T > >       &,                                            // (input) 1st input matrix
    vector< T >                 &,                                            // (input) 2nd input vector
    vector< T >                 &                                             // (input/output) output vector
);
template <class T, size_t m, size_t n>
void matmul(                                                                  // Matrix multiplication with a vector
    array< array < T, n >, m >  &,                                            // (input) 1st input matrix
    array< T, n >               &,                                            // (input) 2nd input vector
    array< T, m >               &                                             // (input/output) output vector
);
template <class T>
void matmul(                                                                  // Matrix multiplication
    vector< vector< T > >       &,                                            // (input) 1st input matrix
    vector< vector< T > >       &,                                            // (input) 2nd input matrix
    vector< vector< T > >       &                                             // (input/output) output matrix
);
template <class T, size_t m, size_t n, size_t l> 
void matmul(                                                                  // Matrix multiplication
    array< array< T, n >, m >   &,                                            // (input) 1st input matrix
    array< array< T, l >, n >   &,                                            // (input) 2nd input matrix
    array< array< T, l >, m >   &                                             // (input/output) output matrix
);

// Matrix manipulations ----------------------------------------------------- //
template <class T>
void transpose(                                                               // Matrix transposition
    vector< vector< T > >       &,                                            // (input) Input matrix
    vector< vector< T > >       &                                             // (input/output) Output transpose
);
template <class T, size_t m, size_t n>
void transpose(                                                               // Matrix transposition
    array< array< T, n >, m >   &,                                            // (input) Input matrix
    array< array< T, m >, n >   &                                             // (input/output) Output transpose
);

template <class T>
vector< vector< T > >  transpose(                                            // Matrix transposition
    vector< vector< T > >       &                                            // (input) Input matrix
);
template <class T, size_t m, size_t n>
array< array< T, m >, n > transpose(                                         // Matrix transposition
    array< array< T, n >, m >   &                                            // (input) Input matrix
);

template <class T>
void complement(                                                              // Complement extraction
    int                          ,                                            // (input) complement 1st index
    int                          ,                                            // (input) complement 2nd index
    vector< vector< T > >       &,                                            // (input) input matrix
    vector< vector< T > >       &                                             // (input/output) (i,j) complement
);
template <class T, size_t m, size_t n>
void complement(                                                              // Complement extraction
    int                          ,                                            // (input) complement 1st index
    int                          ,                                            // (input) complement 2nd index
    array< array< T, n >, m >   &,                                            // (input) input matrix
    array< array<T, n-1>, m-1>  &                                             // (input/output) (i,j) complement
);
template <class T>
void triL(                                                                    // Lower triangular part extraction
    vector< vector< T > >       &,                                            // (input) Input matrix
    vector< vector< T > >       &                                             // (input/output) Lower triangular part
);
template <class T, size_t m, size_t n>
void triL(                                                                    // Extract the lower triangular part
    array< array< T, n >, m >   &,                                            // (input) Input matrix
    array< array< T, n >, m >   &                                             // (input/output) Lower triangular part
);
template <class T>
void triU(                                                                    // Extract the upper triangular part
    vector< vector< T > >       &,                                            // (input) input matrix
    vector< vector< T > >       &                                             // (input/output) Upper triangular part
);
template <class T, size_t m, size_t n>
void triU(                                                                    // Extract the upper triangular part
    array< array< T, n >, m >   &,                                            // (input) input matrix
    array< array< T, n >, m >   &                                             // (input/output) Upper triangular part
);

// Matrix determinant ------------------------------------------------------- //
template <class T>
T det(                                                                        // Matrix determinant
    vector< vector< T > >       &                                             // (input) input matrix
);
template <class T>
T det(                                                                        // Dummy function for self-recursive templated routine "det"
    array< array < T, 1 >, 1 >  &                                             // (input) input matrix
);
template <class T, size_t m, size_t n>
T det(                                                                        // Matrix determinant
    array< array < T, n >, m >  &                                             // (input) input matrix
);

// Linear system ------------------------------------------------------------ //
template <class T>
void Cramer(                                                                  // Solve linear system using Cramer's rule
    vector< vector < T > >      &,                                            // (input) coeff. matrix
    vector< T >                 &,                                            // (input) source term
    vector< T >                 &                                             // (input/output) solution
);
template <class T, size_t m, size_t n>
void Cramer(                                                                  // Solve linear system using Cramer's rule
    array< array < T, n >, m >  &,                                            // (input) coeff. matrix
    array< T, m >               &,                                            // (input) source term
    array< T, n >               &                                             // (input/output) solution
);
unsigned int LU(                                                              // LU decomposition of matrix
    dvector2D                   &,                                            // (input) Input matrix
    dvector2D                   &,                                            // (input/output) L factor
    dvector2D                   &,                                            // (input/output) U factor
    dvector2D                   *                                             // (input/optional) permutation matrix
);
template<size_t m>
unsigned int LU(                                                              // LU decomposition of matrix
    array<array<double, m>, m>  &,                                            // (input) Input matrix
    array<array<double, m>, m>  &,                                            // (input/output) L factor
    array<array<double, m>, m>  &,                                            // (input/output) U factor
    array<array<double, m>, m>  *                                             // (input/optional) permutation matrix
);
void BackwardSubst(                                                           // Backward substitution algorithm
    dvector2D                   &,                                            // (input) Coeffs. matrix
    dvector1D                   &,                                            // (input) Source term
    dvector1D                   &                                             // (input/output) Solution
);
template<size_t m>
void BackwardSubst(                                                           // Backward substitution algorithm
    array<array<double, m>, m>  &,                                            // (input) Coeffs. matrix 
    array<double, m>            &,                                            // (input) Source term
    array<double, m>            &                                             // (input/output) Solution
);
void ForwardSubst(                                                            // Forward substitution algorithm
    dvector2D                   &,                                            // (input) Coeffs. matrix
    dvector1D                   &,                                            // (input) Source term
    dvector1D                   &                                             // (input/output) Solution
);
template<size_t m>
void ForwardSubst(                                                            // Forward substitution algorithm
    array<array<double, m>, m>  &,                                            // (input) Coeffs. matrix
    array<double, m>            &,                                            // (input) Source term
    array<double, m>            &                                             // (input/output) Solution
);
void SolveLU(                                                                 // Solve linear system using LU factorization
    dvector2D                   &,                                            // (input) Input coeffs. matrix
    dvector1D                   &,                                            // (input) Source term
    dvector1D                   &                                             // (input/output) Solution
);
template<size_t m>
void SolveLU(                                                                 // Solve linear system using LU factorization
    array<array<double, m>, m>  &,                                            // (input) Input coeffs. matrix
    array<double, m>            &,                                            // (input) Source term
    array<double, m>            &                                             // (input/output) Solution
);

// =================================================================================== //
// Functions using return                                                              //

// Operator Tensor_Product ------------------------------------------------------------ //
template <class T>
vector<vector<T>> Tensor_Product(vector<T> const       &,
                                 vector<T> const       &);

template <class T, size_t n, size_t m>
array<array<T,m>, n> Tensor_Product(array<T, n> const   &,
                                    array<T, m> const   &);

// Matrix Vector Multiplication ------------------------------------------------------ //
template <class T>
vector<T> Mat_Mul( const vector< vector<T>>  &, 
                   const vector<T>           &);

template <class T, size_t d1, size_t d2>
array<T, d1> Mat_Mul( const array< array<T, d2>, d1> &, 
                      const array<T, d2>             &);


// Matrix Matrix Multiplication ------------------------------------------------------ //
template <class T>
vector< vector<T> > Mat_Mul( const vector< vector<T> > &, 
                             const vector< vector<T> > &);

template <class T>
vector< vector<T> > Dia_Mat_Mul( const vector<T>           &, 
                                 const vector< vector<T> > &);

template <class T>
vector< vector<T> > Dia_Mat_Mul( const vector< vector<T> > &, 
                                 const vector<T>           &);

template <class T, size_t d1, size_t d2, size_t d3>
array< array<T, d2> , d1> Mat_Mul( const array< array<T, d3>, d1>     &, 
                                   const array< array<T, d2>, d3>     &) ;

template <class T, size_t d1, size_t d2>
array< array<T, d2> , d1> Dia_Mat_Mul( const array<T, d1>             & , 
                                       const array< array<T, d2>, d1> & );

template <class T, size_t d1, size_t d2>
array< array<T, d2> , d1> Dia_Mat_Mul( const array< array<T, d2>, d1> & , 
                                       const array<T, d2>             & );


// Matrix transposition -------------------------------------------------------------- //
template <class T>
vector < vector<T>> Transpose( const vector< vector<T> > & );

template <class T, size_t m, size_t n>
array < array<T, n>, m> Transpose( const array< array<T, m>, n> & );


// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "LinearAlgebra.tpp"
# include "LinearAlgebra2.tpp"


# endif
