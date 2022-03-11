/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

namespace bitpit{
namespace linearalgebra{
/*!
 * @ingroup lamultiplication
 * @{
 */
// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between a scalar and a matrix.

    \param[in] A input scalar
    \param[in] B input matrix
    \param[in,out] C product between A and B
*/
template <class T>
void matmul(
    T                                            A,
    const  std::vector< std::vector< T > >      &B,
    std::vector< std::vector< T > >             &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
std::size_t          m, n;

// Counters
std::size_t          i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Matrix
m = B.size();
if (m == 0) {
    return;
}
n = B[0].size();
if (n == 0) {
    return;
}

// ========================================================================== //
// PERFORM PRODUCT                                                            //
// ========================================================================== //

// Resize output variable
C.resize(m);

// Perform product
for (i = 0; i < m; i++) {
    C[i].resize(n, (T) 0.0);
    for (j = 0; j < n; j++) {
        C[i][j] = A * B[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between a scalar and a matrix. Overloading of 
    matmul() function for container array.

    \param[in] A input scalar
    \param[in] B input matrix
    \param[in,out] C product between A and B
*/
template <class T, std::size_t m, std::size_t n>
void matmul(
    T                                            A,
    const std::array< std::array< T, n >, m >   &B,
    std::array< std::array< T, n >, m >         &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
// none

// Counters
std::size_t       i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Matrix
if (m == 0) {
    return;
}
if (n == 0) {
    return;
}

// ========================================================================== //
// PERFORM PRODUCT                                                            //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        C[i][j] = A * B[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between matrix and scalar.

    \param[in] A input matrix
    \param[in] B input scalar
    \param[in,out] C product between A and B
*/
template <class T>
void matmul(
    const std::vector< std::vector< T > >       &B,
    T                                            A,
    std::vector< std::vector< T > >             &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
std::size_t          m, n;

// Counters
// none

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Matrix
m = B.size();
if (m == 0) {
    return;
}
n = B[0].size();
if (n == 0) {
    return;
}

// ========================================================================== //
// PERFORM PRODUCT                                                            //
// ========================================================================== //
matmul(A, B, C);

return; };

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between matrix and scalar. Overloading of matmul()
    function for container array.

    \param[in] A input matrix
    \param[in] B input scalar
    \param[in,out] C product between A and B
*/
template <class T, std::size_t m, std::size_t n>
void matmul(
    const std::array< std::array< T, n >, m >   &B,
    T                                            A,
    std::array< std::array< T, n >, m >         &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Matrix
if (m == 0) {
    return;
}
if (n == 0) {
    return;
}

// ========================================================================== //
// PERFORM PRODUCT                                                            //
// ========================================================================== //
matmul(A, B, C);

return; };

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between vector and matrix.

    \param[in] A input vector
    \param[in] B input matrix
    \param[in,out] C product between A and B
*/
template <class T>
void matmul(
    const std::vector< T >                      &A,
    const std::vector< std::vector < T > >      &B,
    std::vector< T >                            &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
std::size_t           l, m, n;

// Counters
std::size_t           i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Input vector
l = A.size();
if (l == 0) {
    return;
}

// Input matrix
m = B.size();
if (m == 0) {
    return;
}
n = B[0].size();
if (n == 0) {
    return;
}

// Check dimensions coherency
if (l != m) {
    return;
}

// ========================================================================== //
// COMPUTE THE MATRIX PRODUCT                                                 //
// ========================================================================== //

// Resize vector
C.resize(n, 0.0);

// Compute matrix product
for (i = 0; i < n; i++) {
    C[i] = 0.0;
    for (j = 0; j < m; j++) {
        C[i] += A[j]*B[j][i];
    } //next j
} //next i

return; }

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between vector and matrix. Overloading of matmul()
    function for container array.

    \param[in] A input vector
    \param[in] B input matrix
    \param[in,out] C product between A and B
*/
template <class T, std::size_t m, std::size_t n>
void matmul(
    const std::array< T, m >                    &A,
    const std::array< std::array < T, n >, m >  &B,
    std::array< T, n >                          &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
std::size_t        i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Input matrix
if (m == 0) {
    return;
}
if (n == 0) {
    return;
}

// ========================================================================== //
// COMPUTE THE MATRIX PRODUCT                                                 //
// ========================================================================== //
for (i = 0; i < n; i++) {
    C[i] = 0.0;
    for (j = 0; j < m; j++) {
        C[i] += A[j]*B[j][i];
    } //next j
} //next i

return; }

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between matrix and vector.

    \param[in] A input matrix
    \param[in] B input vector
    \param[in,out] C product between A and B
*/
template <class T>
void matmul(
    const std::vector< std::vector < T > >      &A,
    const std::vector< T >                      &B,
    std::vector< T >                            &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
std::size_t           l, m, n;

// Counters
std::size_t           i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Input vector
l = B.size();
if (l == 0) {
    return;
}

// Input matrix
m = A.size();
if (m == 0) {
    return;
}
n = A[0].size();
if (n == 0) {
    return;
}

// Check dimensions coherency
if (l != n) {
    return;
}

// ========================================================================== //
// COMPUTE THE MATRIX PRODUCT                                                 //
// ========================================================================== //

// Resize vector
C.resize(m, 0.0);

// Compute matrix product
for (i = 0; i < m; i++) {
    C[i] = 0.0;
    for (j = 0; j < n; j++) {
        C[i] += B[j]*A[i][j];
    } //next j
} //next i

return; }

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between matrix and vector. Overloading of matmul()
    function for container array.

    \param[in] A input matrix
    \param[in] B input vector
    \param[in,out] C product between A and B
*/
template <class T, std::size_t m, std::size_t n>
void matmul(
    const std::array< std::array < T, n >, m >  &A,
    const std::array< T, n >                    &B,
    std::array< T, m >                          &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
std::size_t        i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// Input matrix
if (m == 0) {
    return;
}
if (n == 0) {
    return;
}

// ========================================================================== //
// COMPUTE THE MATRIX PRODUCT                                                 //
// ========================================================================== //
for (i = 0; i < m; i++) {
    C[i] = 0.0;
    for (j = 0; j < n; j++) {
        C[i] += B[j]*A[i][j];
    } //next j
} //next i

return; }

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between two matrices.

    \param[in] A input matrix
    \param[in] B input matrix
    \param[in,out] C product between A and B
*/
template <class T>
void matmul(
    const std::vector< std::vector< T > >       &A,
    const std::vector< std::vector< T > >       &B,
    std::vector< std::vector< T > >             &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
std::size_t          m1, n1, n2, m2;

// Counters
std::size_t          i, j, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// 1st matrix
m1 = A.size();
if (m1 == 0) {
    return;
}
n1 = A[0].size();
if (n1 == 0) {
    return;
}

// 2nd matrix
m2 = B.size();
if (m2 == 0) {
    return;
}
n2 = B[0].size();
if (n2 == 0) {
    return;
}

// Check dimensions coherency
if (n1 != m2) {
    return;
}

// ========================================================================== //
// PERFORM PRODUCT                                                            //
// ========================================================================== //

// Resiz output variable
C.resize(m1);

for (i = 0; i < m1; i++) {
    C[i].resize(n2, (T) 0.0);
    for (j = 0; j < n2; j++) {
        C[i][j] = (T) 0.0;
        for (k = 0; k < n1; k++) {
            C[i][j] += A[i][k] * B[k][j];
        } //next k
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Matrix multiplication between two matrices. Overloading of matmul() for
    container array.

    \param[in] A input matrix
    \param[in] B input matrix
    \param[in,out] C product between A and B
*/
template <class T, std::size_t m, std::size_t n, std::size_t l>
void matmul(
    const std::array< std::array< T, n >, m >   &A,
    const std::array< std::array< T, l >, n >   &B,
    std::array< std::array< T, l >, m >         &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
// none

// Counters
std::size_t       i, j, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //

// 1st matrix
if (m == 0) {
    return;
}
if (n == 0) {
    return;
}

// 2nd matrix
if (l == 0) {
    return;
}

// ========================================================================== //
// PERFORM PRODUCT                                                            //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j < l; j++) {
        C[i][j] = (T) 0.0;
        for (k = 0; k < n; k++) {
            C[i][j] += A[i][k] * B[k][j];
        } //next k
    } //next j
} //next i

return; };

// ----------------------------------------------------------------------------------- //
/*!
    Matrix product.

    \param[in] M 1st argument of matrix multiplication
    \param[in] N 2nd argument of matrix multiplication

    \result product between M and N.
*/
template <class T>
std::vector< std::vector<T> > matmul(
    const std::vector< std::vector<T> >         &M,
    const std::vector<std::vector<T> >          &N
) {

    std::vector< std::vector<T> > Tr = transpose(N);

    std::size_t d1= M.size();
    std::size_t d2= Tr.size();

    std::vector< std::vector<T> > Q(d1, std::vector<T> (d2, T()) );

    for( std::size_t i=0; i<d1; i++){
        for( std::size_t j=0; j<d2; j++){
            Q[i][j]= dotProduct( M[i], Tr[j] );
        };
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Matrix multiplication. Overloading of matmul() function for container array.
    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result product of M and N.
*/
template <class T, std::size_t d1, std::size_t d2, std::size_t d3>
std::array< std::array<T, d2> , d1> matmul(
    const std::array< std::array<T, d3>, d1>    &M,
    const std::array<std::array<T, d2>, d3>     &N
){
    std::size_t i, j;

    std::array< std::array<T, d2> , d1> Q;
    std::array< std::array<T, d3> , d2> Tr;

    Tr = transpose( N ) ;

    for( i=0; i<d1; i++){
        for( j=0; j<d2; j++){
            Q[i][j]= dotProduct( M[i], Tr[j] );
        };
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplicationx

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T>
std::vector< std::vector<T> > matmulDiag(
    const std::vector<T>                        &M,
    const std::vector<std::vector<T> >          &N
) {

    std::size_t d1= M.size();
    std::vector< std::vector<T> > Q( N );

    for( std::size_t i=0; i<d1; i++){
        Q[i] *= M[i] ;
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplication

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T>
std::vector< std::vector<T> > matmulDiag(
    const std::vector< std::vector<T> >         &M,
    const std::vector<T>                        &N
) {

    std::vector< std::vector<T> > Q( M );

    std::size_t d1= M.size() ;

    for( std::size_t i=0; i<d1; i++ ){
        Q[i] = M[i] * N ;
    };

    return (Q);

};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplication. Overloading of matmulDiag() function for container array.

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T, std::size_t d1, std::size_t d2>
std::array< std::array<T, d2> , d1> matmulDiag(
    const std::array< T, d1>                    &M,
    const std::array<std::array<T, d2>, d1>     &N
){

    std::size_t i;
    std::array< std::array<T, d2> , d1> Q(N);

    for( i=0; i<d1; i++){
        Q[i] *= M[i]  ;
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplication. Overloading of matmulDiag() function for container array.

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T, std::size_t d1, std::size_t d2>
std::array< std::array<T, d2> , d1> matmulDiag(
    const std::array<std::array<T, d2>, d1>     &M,
    const std::array< T, d2>                    &N
) {

    std::size_t i;
    std::array< std::array<T, d2> , d1> Q;

    for( i=0; i<d1; i++){
        Q[i] = M[i] *N ;
    };

    return (Q);
};

// Matrix Vector Multiplication ====================================================== //

// ----------------------------------------------------------------------------------- //
/*!
    Matrix-vector multiplication.

    \param[in] M input matrix
    \param[in] x input vector

    \result product between M and x
*/
template <class T>
std::vector<T> matmul(
    const std::vector< std::vector<T>>          &M,
    const std::vector<T>                        &x
) {

    std::size_t d1 = M.size();

    std::vector<T>      z(d1,0.0);

    for( std::size_t i=0; i<d1; i++){
        z[i]= dotProduct( M[i], x );
    }

    return (z);
};

// ----------------------------------------------------------------------------------- //
/*!
    Matrix-vector multiplication. Overloading of matmul() for container array.

    \param[in] M input matrix
    \param[in] x input vector

    \result product between M and x
*/
template <class T, std::size_t d1, std::size_t d2>
std::array<T, d1> matmul(
    const std::array< std::array<T, d2>, d1>    &M,
    const std::array<T, d2>                     &x
) {

    std::array<T, d1>      z;

    for( std::size_t i=0; i<d1; i++){
        z[i]= dotProduct( M[i], x);
    }

    return (z);
};

// ----------------------------------------------------------------------------------- //
/*!
    Tensor product.

    \param[in] x 1st argument of the tensor product
    \param[in] y 2nd argument of the tensor product

    \result tensor product between x and y
*/
template <class T>
std::vector<std::vector<T>> tensorProduct(
    const std::vector<T>                        &x,
    const std::vector<T>                        &y
) {

    std::size_t  n = x.size();
    std::size_t  m = y.size();
    std::vector<T>      row(m,0.0);
    std::vector<std::vector<T>> z(n,row) ;

    for( std::size_t i=0; i<n; i++){
        for( std::size_t j=0; j<m; j++){
            z[i][j] = x[i] *y[j] ;
        };
    };

return (z);}

// ----------------------------------------------------------------------------------- //
/*!
    Tensor product. Overloading of tensorProduct() for container array.

    \param[in] x 1st argument of the tensor product
    \param[in] y 2nd argument of the tensor product

    \result tensor product between x and y
*/
template <class T, std::size_t n, std::size_t m>
std::array<std::array<T,m>,n> tensorProduct(
    const std::array<T,n>                       &x,
    const std::array<T,m>                       &y
) {
    std::array<std::array<T,m>,n> z ;

    for( std::size_t i=0; i<n; i++){
        for( std::size_t j=0; j<m; j++){
            z[i][j] = x[i] *y[j] ;
        };
    };

    return (z);
}

/*!
 * @}
 */
}
}
