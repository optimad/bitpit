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

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //

/*!
 \ingroup   LinearAlgebra
 \{
 */
// Generic routines ========================================================= //

// -------------------------------------------------------------------------- //
/*!
    Display matrix to output stream in a nicely formatted form.

    \param[in,out] out output stream
    \param[in] A matrix to be displayed
*/
template<class T>
void display_matrix(
    std::ostream                                &out,
    std::vector< std::vector< T > >             &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int             m = A.size();

// Counters
int             i;

// ========================================================================== //
// DISPLAY MATRIX CONTENT                                                     //
// ========================================================================== //
for (i = 0; i < m; ++i) {
    out << A[i] << std::endl;
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Display matrix to output stream in a nicely formatted form.
    Overloading of display_matrix() function for container array.

    \param[in,out] out output stream
    \param[in] A matrix to be displayed
*/
template<class T, size_t m, size_t n>
void display_matrix(
    std::ostream                                &out,
    std::array<std::array<T, n>, m>             &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int             i;

// ========================================================================== //
// DISPLAY MATRIX CONTENT                                                     //
// ========================================================================== //
for (i = 0; i < m; ++i) {
    out << A[i] << std::endl;
} //next i

return; };

// Matrix basic templates =================================================== //

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n matrix with 0 entries

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T>
void zeros(
    std::vector< std::vector < T > >            &A,
    int                                          m,
    int                                          n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int      i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
A.resize(m);
for (i = 0; i < m; i++) {
    A[i].resize(n, (T) 0.0);
    for (j = 0; j < n; j++) {
        A[i][j] = (T) 0.0;
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initalize a m-by-n matrix with 0-entries. Overloading of zeros() function
    for container array.

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T, size_t m, size_t n>
void zeros(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// template <class T, size_t m, size_t n>                                     //
// void zeros(                                                                //
//     array< array < T, n >, m >  &A)                                        //
//                                                                            //
// Initialize a m-by-n matrix of zeros.                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A      : array< array< T, n >, m >, with m-by-n matrix of zeros          //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int      i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
for (i = 0; i < m; i++) {
    A[i].fill(0.0);
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n matrix of ones.

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T>
void ones(
    std::vector< std::vector < T > >            &A,
    int                                          m,
    int                                          n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int      i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    std::cout << "ERROR: number of rows (columns) must be > 0!!" << std::endl;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
A.resize(m);
for (i = 0; i < m; i++) {
    A[i].resize(n, (T) 1.0);
    for (j = 0; j < n; j++) {
        A[i][j] = (T) 1.0;
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n matrix of ones. Overloading of ones() function for
    container array.

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T, size_t m, size_t n>
void ones(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int      i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
for (i = 0; i < m; i++) {
    A[i].fill(1.0);
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n identity matrix having. All the entries on the main
    diagonal are set to 1.

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T>
void eye(
    std::vector< std::vector < T > >            &A,
    int                                          m,
    int                                          n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int      s = std::min(m, n);

// Counters
int      i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
zeros(A, m, n);
for (i = 0; i < s; i++) {
    A[i][i] = (T) 1.0;
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Initialize a m-by-n identity matrix having. All the entries on the main
    diagonal are set to 1. Overloading of eye() function for container array.

    \param[in,out] A container for matrix storage
    \param[in] m number of matrix rows
    \param[in] n number of matrix columns
*/
template <class T, size_t m, size_t n>
void eye(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int      s = std::min(m, n);

// Counters
int      i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if ((m == 0) || (n == 0)) {
    return;
}

// ========================================================================== //
// CREATE MATRIX                                                              //
// ========================================================================== //
zeros(A);
for (i = 0; i < s; i++) {
    A[i][i] = (T) 1.0;
} //next i

return; };

// Matrix multiplications =================================================== //

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
    std::vector< std::vector< T > >             &B,
    std::vector< std::vector< T > >             &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
int          m, n;

// Counters
int          i, j;

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
template <class T, size_t m, size_t n>
void matmul(
    T                                            A,
    std::array< std::array< T, n >, m >         &B,
    std::array< std::array< T, n >, m >         &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
// none

// Counters
int          i, j;

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
    std::vector< std::vector< T > >             &B,
    T                                            A,
    std::vector< std::vector< T > >             &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
int          m, n;

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
template <class T, size_t m, size_t n>
void matmul(
    std::array< std::array< T, n >, m >         &B,
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
    std::vector< T >                            &A,
    std::vector< std::vector < T > >            &B,
    std::vector< T >                            &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int           l, m, n;

// Counters
int           i, j;

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
template <class T, size_t m, size_t n>
void matmul(
    std::array< T, m >                          &A,
    std::array< std::array < T, n >, m >        &B,
    std::array< T, n >                          &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int           i, j;

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
    std::vector< std::vector < T > >            &A,
    std::vector< T >                            &B,
    std::vector< T >                            &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int           l, m, n;

// Counters
int           i, j;

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
template <class T, size_t m, size_t n>
void matmul(
    std::array< std::array < T, n >, m >        &A,
    std::array< T, n >                          &B,
    std::array< T, m >                          &C
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int           i, j;

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
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &B,
    std::vector< std::vector< T > >             &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
int          m1, n1, n2, m2;

// Counters
int          i, j, k;

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
template <class T, size_t m, size_t n, size_t l>
void matmul(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, l >, n >         &B,
    std::array< std::array< T, l >, m >         &C
) {

// ========================================================================== //
// VARIABLE DECLARATION                                                       //
// ========================================================================== //

// Local variables
// none

// Counters
int          i, j, k;

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


// Matrix manipulation ====================================================== //

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose.

    \param[in] A input matrix
    \param[in,out] B transpose of A
*/
template <class T>
void transpose(
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            m, n;

// Counters
int            i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; }
n = A[0].size();
if (n == 0) { return; }

// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //

// Resize output variables
B.resize(n);
for (i = 0; i < n; ++i) {
    B[i].resize(m, (T) 0.0);
} //next i

// Transposition
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose. Overloading of transpose()
    function for container array.

    \param[in] A input matrix
    \param[in,out] B transpose of A
*/
template <class T, size_t m, size_t n>
void transpose(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, m >, n >         &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int            i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; }
if (n == 0) { return; }

// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose. Overloading of transpose()
    function.

    \param[in] A input matrix
    \param[in,out] B transpose of A
*/
template <class T>
std::vector< std::vector< T > > transpose(
    std::vector< std::vector< T > >             &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                        m = A.size(); 
int                        n(0);

// Counters
int            i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m != 0) n = A[0].size() ;

// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //

std::vector< std::vector< T > >      B( n, std::vector<T> (m,0.) ) ;

// Transposition
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix, compute its transpose. Overloading of transpose()
    function for container array.

    \param[in] A input matrix
    \param[in,out] B transpose of A
*/
template <class T, size_t m, size_t n>
std::array< std::array< T, m >, n > transpose(
    std::array< std::array< T, n >, m >         &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
std::array< std::array< T, m >, n >  B ;

// Counters
int            i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
// ========================================================================== //
// MATRIX TRANSPOSITION                                                       //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
        B[j][i] = A[i][j];
    } //next j
} //next i

return B; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix A, compute the complement of the element A[i][j].

    \param[in] i row index of element
    \param[in] j column index of element
    \param[in] A input matrix
    \param[in,out] B complement of alement A[i][j].
*/
template <class T>
void complement(
    int                                          i,
    int                                          j,
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            m, n;

// Counters
int            l, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; }
n = A[0].size();
if (n == 0) { return; }
i--; j--;
if ((i >= m) || (i < 0)) {
    return;
}
if ((j >= n) || (j < 0)) {
    return;
}

// ========================================================================== //
// EXTRACT COMPLEMENT                                                         //
// ========================================================================== //

// Resize output variables
B.resize(m-1);
for (i = 0; i < m-1; ++i) {
    B[i].resize(n-1, 0.0);
} //next i

// Extract complement
for (l = 0; l < i; l++) {
    for (k = 0; k < j; k++) {
        B[l][k] = A[l][k];
    } //next k
    for (k = j+1; k < n; k++) {
        B[l][k-1] = A[l][k];
    } //next k
} //next l
for (l = i+1; l < m; l++) {
    for (k = 0; k < j; k++) {
        B[l-1][k] = A[l][k];
    } //next k
    for (k = j+1; k < n; k++) {
        B[l-1][k-1] = A[l][k];
    } //next k
} //next l

return; };

// -------------------------------------------------------------------------- //
/*!
    Given an input matrix A, compute the complement of the element A[i][j].
    Overloading of complement() function for container array.

    \param[in] i row index of element
    \param[in] j column index of element
    \param[in] A input matrix
    \param[in,out] B complement of alement A[i][j].
*/
template <class T, size_t m, size_t n>
void complement(
    int                                          i,
    int                                          j,
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array<T, n-1>, m-1>        &B
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int            l, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; }
if (n == 0) { return; }
i--; j--;
if ((i >= m) || (i < 0)) {
    return;
}
if ((j >= n) || (j < 0)) {
    return;
}

// ========================================================================== //
// EXTRACT COMPLEMENT                                                         //
// ========================================================================== //
for (l = 0; l < i; l++) {
    for (k = 0; k < j; k++) {
        B[l][k] = A[l][k];
    } //next k
    for (k = j+1; k < n; k++) {
        B[l][k-1] = A[l][k];
    } //next k
} //next l
for (l = i+1; l < m; l++) {
    for (k = 0; k < j; k++) {
        B[l-1][k] = A[l][k];
    } //next k
    for (k = j+1; k < n; k++) {
        B[l-1][k-1] = A[l][k];
    } //next k
} //next l

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the lower triangular part of a matrix.

    \param[in] A input matrix
    \param[in,out] L lower triangular part of A
*/
template <class T>
void triL(
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &L
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int        m, n;

// Counters
int        i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
m = A.size();
if (m == 0) { return; };
n = A[0].size();

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(L, m, n);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j <= i; j++) {
        L[i][j] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the lower triangular part of a matrix. Overloading of triL() function
    for container array.

    \param[in] A input matrix
    \param[in,out] L lower triangular part of A
*/

template <class T, size_t m, size_t n>
void triL(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, n >, m >         &L
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int        i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(L);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = 0; j <= i; j++) {
        L[i][j] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the upper triangular part of a matrix.

    \param[in] A input matrix
    \param[in,out] U upper triangular part of A
*/
template <class T>
void triU(
    std::vector< std::vector< T > >             &A,
    std::vector< std::vector< T > >             &U
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int        m, n;

// Counters
int        i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
m = A.size();
if (m == 0) { return; };
n = A[0].size();

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(U, m, n);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = i; j < n; j++) {
        U[i][j] = A[i][j];
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Extract the upper triangular part of a matrix. Overloading of triU() function
    for container array.

    \param[in] A input matrix
    \param[in,out] U upper triangular part of A
*/
template <class T, size_t m, size_t n>
void triU(
    std::array< std::array< T, n >, m >         &A,
    std::array< std::array< T, n >, m >         &U
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int        i, j;

// ========================================================================== //
// CHECK INPUT COHERENCY                                                      //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
zeros(U);

// ========================================================================== //
// EXTRACT THE LOWER TRIANGULAR PART OF A                                     //
// ========================================================================== //
for (i = 0; i < m; i++) {
    for (j = i; j < n; j++) {
        U[i][j] = A[i][j];
    } //next j
} //next i

return; };

// Matrix determinant ======================================================= //

// -------------------------------------------------------------------------- //
/*!
    End function for recursive calls to det().

    \param[in] A input matrix

    \result determinant of A.
*/
template <class T>
T det(
    std::array< std::array< T, 1 >, 1 >         &A
) {

return(A[0][0]); };

// -------------------------------------------------------------------------- //
/*!
    Compute determinant of a matrix of small dimenions using Laplace rule.

    \param[in] A input matrix

    \result determinant of A.
*/
template <class T>
T det(
    std::vector< std::vector < T > >            &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                     m, n;
T                       d = (T) 1.0e+18;
std::vector< std::vector < T > >  C;

// Counters
int                     i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return(d); };
n = A[0].size();
if (n == 0) { return(d); };
if (m != n) {
    return(d);
}

// ========================================================================== //
// COMPUTE DETERMINANT                                                        //
// ========================================================================== //
if (m == 1) {
    d = A[0][0];
}
else if (m == 2) {
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0];
}
else if (m == 3) {
    d = A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1]
        + A[0][1]*A[1][2]*A[2][0] - A[0][1]*A[1][0]*A[2][2]
        + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0];
}
else {
    d = (T) 0.0;
    for (i = 0; i < m; i++) {
        complement(1, i+1, A, C);
        d += pow(-1.0, i+2) * A[0][i] * det(C);
    } //next i
}

return(d); };

// -------------------------------------------------------------------------- //
/*!
    Compute determinant of a matrix of small dimenions using Laplace rule.
    Overloading of det() function for container array.

    \param[in] A input matrix

    \result determinant of A.
*/
template <class T, size_t m, size_t n>
T det(
    std::array< std::array < T, n >, m >        &A
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
T                       d = (T) 0.0;

// Counters
int                     i;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return(d); };
if (n == 0) { return(d); };
if (m != n) {
    return(d);
}

// ========================================================================== //
// COMPUTE DETERMINANT                                                        //
// ========================================================================== //
if (m == 1) {
    d = A[0][0];
}
else if (m == 2) {
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0];
}
else if (m == 3) {
    d = A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1]
        + A[0][1]*A[1][2]*A[2][0] - A[0][1]*A[1][0]*A[2][2]
        + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0];
}
else {
    std::array< std::array < double, m-1 >, m-1 >     C;
    for (i = 0; i < m; i++) {
        complement(1, i+1, A, C);
        d += pow(-1.0, i+2) * A[0][i] * det(C);
    } //next i
}

return(d); };

// Linear system ============================================================ //

// -------------------------------------------------------------------------- //
/*!
    Solve a linear system of small dimenions using Cramer's rule.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in,out] x on output stores the solution to the linear system
*/
template <class T>
void Cramer(
    std::vector< std::vector < T > >            &A,
    std::vector< T >                            &B,
    std::vector< T >                            &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                    l, m, n;
T                      dA;
std::vector< std::vector < T > > C;

// Counters
int                    i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
m = A.size();
if (m == 0) { return; }
n = A[0].size();
if (n == 0) { return; }
if (m != n) {
    return;
}
l = B.size();
if (l == 0) { return; }
if (l != m) {
    return;
}

// =================================================================================== //
// SOLVE LINEAR SYSTEM                                                                 //
// =================================================================================== //

// Solvability condition
dA = det(A);
if (dA < 1.0e-14) {
    return;
}

// Solve linear system
x.resize(n, (T) 0.0);
for (i = 0; i < m; i++) {

    // Build B
    C = A;
    for (j = 0; j < m; j++) {
        C[j][i] = B[j];
    } //next j
    x[i] = det(C)/dA;

} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a linear system of small dimenions using Cramer's rule. Overloading
    of Cramer() function for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in,out] x on output stores the solution to the linear system
*/
template <class T, size_t m, size_t n>
void Cramer(
    std::array< std::array < T, n >, m >        &A,
    std::array< T, m >                          &B,
    std::array< T, n >                          &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
T                           dA;
std::array< std::array < T, n >, m >  C;

// Counters
int                         i, j;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; }
if (n == 0) { return; }
if (m != n) {
    return;
}

// ========================================================================== //
// SOLVE LINEAR SYSTEM                                                        //
// ========================================================================== //

// Solvability condition
dA = det(A);
if (dA < 1.0e-14) {
    return;
}

// Solve linear system
for (i = 0; i < m; i++) {

    // Build B
    C = A;
    for (j = 0; j < m; j++) {
        C[j][i] = B[j];
    } //next j
    x[i] = det(C)/dA;

} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Compute the LU factorization (with row pivoting) for a given non-singular matrix.
    Overloading of LU() function for container array.

    \param[in] A input matrix
    \param[in,out] L lower triangular part (L factor)
    \param[in,out] U upper triangular part (U factor)
    \param[in,out] P permutation matrix

    \result returns an error flag:
        err = 0: no error(s) encountered
        err = 1: matrix is ill conditioned
        err = 2: matrix is singular to working precision
        err = 3: wrong dimensions
*/
template<size_t m>
unsigned int LU(
    std::array< std::array < double, m >, m >   &A,
    std::array< std::array < double, m >, m >   &L,
    std::array< std::array < double, m >, m >   &U,
    std::array< std::array < double, m >, m >   *P
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double            toll_pivot = 1.0e-8;

// Local variables
int                             info = 0;
int                             pivot_row;
double                          pivot, pivot_trial;
std::array<std::array<double, m>, m>      AA;

// Counter
int                             i, j, k;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return(3); };

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //

// LU matrices
zeros(L);
zeros(U);

// Backup copy of coeffs. matrix
AA = A;

// Pivoting array
eye(*P);

// ========================================================================== //
// COMPUTE LU FACTORIZATION                                                   //
// ========================================================================== //
for (k = 0; k < m; k++) {
    L[k][k] = 1.0;

    // Pivoting ------------------------------------------------------------- //
    pivot_row = k;
    pivot = abs(AA[k][k]);
    for (i = k+1; i < m; i++) {
        pivot_trial = abs(AA[i][k]);
        if (pivot_trial > pivot) {
            pivot = pivot_trial;
            pivot_row = i;
        }
    } //next i

    // Perform rows permutation --------------------------------------------- //
    if (pivot_row == k) {
        if (pivot < 1.0e-14) {
            info = 2;
            return(info);
        }
        else if ((pivot >= 1.0e-14) && (pivot < toll_pivot)) {
            info = 1;
        }
    }
    else {
        swap(AA[k], AA[pivot_row]);
        if (P != NULL) {
            swap((*P)[k], (*P)[pivot_row]);
        }
    }

    // Gauss elimination ---------------------------------------------------- //
    for (i = k+1; i < m; i++) {
        L[i][k] = AA[i][k]/AA[k][k] ;
        for (j = k+1; j < m; j++) {
            AA[i][j] = AA[i][j] - L[i][k]*AA[k][j];
        } //next j

    } //next i
    for (j = k; j < m; j++) {
        U[k][j] = AA[k][j];
    } //next j
} //next k

return(info); };

// -------------------------------------------------------------------------- //
/*!
    Solve a lower triangular system using backward substitution. Overloading
    of BackwardSubst() for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in] x on output stores the solution to the linear system
*/
template<size_t m>
void BackwardSubst(
    std::array< std::array < double, m >, m >   &A,
    std::array< double, m >                     &B,
    std::array< double, m >                     &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double sum, d;

// Counter
int    i, p;

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// CHECK SOLVABILITY CONDITION                                                //
// ========================================================================== //
d = 1.0;
for (i = 0; i < m; i++) {
    d = d*A[i][i];
} //next i
if (abs(d) < 1.0e-14) {
    return;
}

// ========================================================================== //
// SOLVE LINEAR SYSTEM WITH BACKWARD SUBSTITUTION                             //
// ========================================================================== //
for (i = m-1; i >= 0; i--) {
    sum = 0.0;
    for(p = m-1; p > i; p--) {
        sum += A[i][p]*x[p];
    } //next p
    x[i] = (B[i] - sum)/A[i][i];
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a upper triangular system using forward substitution. Overloading
    of ForwardSubst() for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of the linear system
    \param[in] x on output stores the solution to the linear system
*/
template<size_t m>
void ForwardSubst(
    std::array< std::array < double, m >, m >   &A,
    std::array< double, m >                     &B,
    std::array< double, m >                     &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
double     d, sum;

// Counters
int        i, p;


// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
if (m == 0) { return; };

// ========================================================================== //
// CHECK SOLVABILITY CONDITION                                                //
// ========================================================================== //
d = 1.0;
for (i = 0; i < m; i++) {
    d = d*A[i][i];
} //next i
if (abs(d) < 1.0e-14) {
    return;
}

// ========================================================================== //
// FORWARD SUBSTITUTION                                                       //
// ========================================================================== //
for(i = 0; i < m; i++) {
    sum = 0.0;
    for(p = 0; p < i; p++) {
        sum += A[i][p] * x[p];
    } //next p
    x[i] = (B[i] - sum)/A[i][i];
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Solve a non-singular linear system using LU factorization. Overloading
    of LU() function for container array.

    \param[in] A coeffs matrix
    \param[in] B r.h.s. of linear system
    \param[in,out] x on output stores the solution to the linear system
*/
template<size_t m>
void SolveLU(
    std::array< std::array< double, m >, m >    &A,
    std::array< double, m >                     &B,
    std::array< double, m >                     &x
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int                info;
std::array<std::array<double, m>, m>  L, U, P, *P_ = &P;
std::array<double, m>            z, C;

// Counters
// none

// ========================================================================== //
// COMPUTE LU FACTORIZATION                                                   //
// ========================================================================== //
info = LU(A, L, U, P_);
if ((info == 2) || (info == 3)) {
    return;
}
matmul(P, B, C);

// ========================================================================== //
//  SOLVE THE LINEAR SYSTEM                                                   //
// ========================================================================== //

// Forward substitution
ForwardSubst(L, C, z);

// Bacward substitution
BackwardSubst(U, z, x);

return; };

/*!
 \}
 */
