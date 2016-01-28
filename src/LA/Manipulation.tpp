namespace bitpit{

/*!
 * @ingroup Manipulation
 * @{
 */

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


// ----------------------------------------------------------------------------------- //
/*!
    Compute matrix transpose.

    \param[in] M input matrix
    
    \result transpose of M
*/
template <class T>
std::vector < std::vector<T> > Transpose(
    const std::vector< std::vector<T>>          &M
) {

    int n = M.size();
    int m = M[0].size() ;
    std::vector< std::vector<T> >   z(m, std::vector<T>(n, 0.0));

    for( int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            z[i][j] = M[j][i];
        };
    };

    return (z);
};

// ----------------------------------------------------------------------------------- //
/*!
    Compute matrix transpose. Overloading of Transpose() function for container
    array.

    \param[in] M input matrix
    
    \result transpose of M
*/
template <class T, size_t m, size_t n>
std::array< std::array<T, n>, m> Transpose(
    const std::array< std::array<T, m>, n>      &M
) {

    std::array< std::array<T, n> ,m>      z;

    for( int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            z[i][j] = M[j][i];
        };
    };

    return (z);
};

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

/*!
 * @}
 */
}
